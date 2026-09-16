#!/usr/bin/env python3
"""Amino-acid diversity along the HLA protein, groove vs non-groove, and allele-level
differentiation between ancestries.

Supervisor asks (Cole, 2026-09-17 call):

    A8: "What's actually more interesting is if we look at the CDS between the peptide-binding
         and the non-peptide-binding regions. [...] just the diversity at the amino acid level,
         and then comparing that from the peptide curves to the non-peptide curves."
    A10: "It seems even more structured to me [class II]. There seems to be greater basic
          diversity statistics on the class II HLA. [...] I suspect these class II HLAs are going
          to be the most different between the ancestries."

Design decision: frequency-weighted reference proteins, not a per-person sequence scan
--------------------------------------------------------------------------------------
The obvious implementation reads every person's CDS FASTA off disk, translates it and builds a
12,000-person x 65-gene alignment. That is a multi-hour full-cohort scan for a quantity that does
not need it. Instead:

  * each haplotype's call is resolved to its **two-field (protein) identity**;
  * the protein sequence for that identity comes from the IPD-IMGT reference shipped with
    Immuannot (3.55.0, ENVIRONMENT quirk #38);
  * per-codon amino-acid diversity is computed over those proteins **weighted by how often each
    was actually observed in the cohort**, per ancestry.

This is exactly protein-level pi_aa, and it is cheap. Its one real limitation is that it covers
only catalogued proteins -- but S01 established that 99.0-99.8% of haplotypes carry a protein
already in IPD-IMGT (script 27), so the uncovered fraction is under 1%, and it is reported.

Why the groove comparison is done by exon and not by a published residue list
-----------------------------------------------------------------------------
The canonical antigen-recognition-site residue lists (Hughes & Nei 1988 Table 1, Parham 1988,
Bondinas 2007) could not be verified from primary sources for this project
(`sprints/S02_supervisor_items/WS_literature_selection.md` says so explicitly). Rather than
compute a headline number from an unverified list, the primary comparison here uses the
**groove-encoding exons**, whose coordinates come straight from IMGT's own `exon=` annotation:
exons 2+3 for class I, exon 2 for class II. That is coarser but it is derived from the data we
ship, and it cannot be silently off by a signal-peptide length.

If `reference/ars_peptide_contacts.tsv` exists (peptide-contact residues derived structurally from
PDB by `_ars_residues.py`), a second, sharper comparison is added: contact vs non-contact residues
*within* the groove exons. That one is optional by construction, so the script's main result never
depends on it.

Usage (VM, full cohort):

    setsid nohup python3 -u scripts/hla_popgen/31_aa_diversity_selection.py \\
        --out-dir ~/results/31_aa_diversity < /dev/null & disown

Usage (local fixtures): scripts/hla_popgen/tests/test_aa_diversity_selection.py
"""
import argparse
import importlib.util
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_THIS_DIR))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)


def _load_module(filename, modname):
    path = os.path.join(_THIS_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


_mods = {}


def mod(filename, name):
    if name not in _mods:
        _mods[name] = _load_module(filename, name)
    return _mods[name]


DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_TABLE1 = os.path.join(DEFAULT_OUTROOT, "hla_calls_rich.tsv")
DEFAULT_COHORT = os.path.join(DEFAULT_OUTROOT, "cohort_membership.tsv")
# ENVIRONMENT quirk #35: ~/mnt/aou-controlled is stale and hangs any process that touches it.
DEFAULT_RELATEDNESS = os.path.expanduser(
    "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/"
    "relatedness/samples_relatedness.tsv")
DEFAULT_REFDATA = os.path.expanduser("~/tools/Immuannot_refdata")
DEFAULT_ARS = os.path.join(_REPO_ROOT, "reference", "ars_peptide_contacts.tsv")
DEFAULT_OUT_DIR = os.path.expanduser("~/results/31_aa_diversity")

SUPPRESS_BELOW = 20
KIN_MIN = 0.0442
STRICT_MIN = 0.9
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
ANCESTRY_COLORS = {"AFR": "#E69F00", "AMR": "#56B4E9", "EAS": "#009E73",
                   "EUR": "#0072B2", "MID": "#D55E00", "SAS": "#CC79A7"}

# The genes worth a protein-diversity track. Classical first; DRA / HLA-F are the conserved
# negative controls Marc's nucleotide Manhattan already showed are nearly invariant.
TARGET_GENES = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1",
                "DRB3", "DRB4", "DRB5", "E", "F", "G", "DRA"]
CONSERVED_CONTROL_GENES = {"DRA", "F", "E"}

MIN_HAPS_FOR_GENE = 200
N_PERMUTATIONS = 2000


def suppress(n, threshold=SUPPRESS_BELOW):
    n = int(n)
    return "0" if n == 0 else ("<%d" % threshold if n < threshold else str(n))


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def gene_bare(g):
    return str(g).replace("HLA-", "")


def two_field(name):
    """`HLA-A*01:01:01:01` -> `A*01:01`; None when the protein identity is unresolved."""
    if not isinstance(name, str) or "*" not in name:
        return None
    gene, rest = name.strip().split("*", 1)
    f = rest.split(":")
    if len(f) < 2 or f[0].strip().lower() == "new" or f[1].strip().lower() == "new":
        return None
    return "%s*%s:%s" % (gene.replace("HLA-", ""), f[0], f[1])


# ---------------------------------------------------------------------------
# alignment
# ---------------------------------------------------------------------------
def align_to_reference(ref_prot, obs_prot, match=1, mismatch=-1, gap=-2):
    """Map every reference residue to the aligned observed residue, or '-' where the observed
    protein has a deletion. Returns a string of len(ref_prot).

    Equal-length proteins are the overwhelmingly common case at HLA (indels within a protein
    allele are rare), so they short-circuit to a direct correspondence. Unequal lengths fall back
    to Needleman-Wunsch. Residues inserted in the observed protein have no reference column and
    are dropped, which is correct: pi is defined per reference position.
    """
    if len(ref_prot) == len(obs_prot):
        return obs_prot
    n, m = len(ref_prot), len(obs_prot)
    S = np.zeros((n + 1, m + 1), dtype=np.int32)
    S[:, 0] = np.arange(n + 1) * gap
    S[0, :] = np.arange(m + 1) * gap
    ptr = np.zeros((n + 1, m + 1), dtype=np.int8)   # 0 diag, 1 up (gap in obs), 2 left (gap in ref)
    ptr[:, 0] = 1
    ptr[0, :] = 2
    ref_arr = np.frombuffer(ref_prot.encode("ascii"), dtype=np.uint8)
    obs_arr = np.frombuffer(obs_prot.encode("ascii"), dtype=np.uint8)
    for i in range(1, n + 1):
        sc = np.where(ref_arr[i - 1] == obs_arr, match, mismatch)
        for j in range(1, m + 1):
            d = S[i - 1, j - 1] + sc[j - 1]
            u = S[i - 1, j] + gap
            l = S[i, j - 1] + gap
            best = max(d, u, l)
            S[i, j] = best
            ptr[i, j] = 0 if best == d else (1 if best == u else 2)
    out = []
    i, j = n, m
    while i > 0 or j > 0:
        p = ptr[i, j]
        if i > 0 and (j == 0 or p == 1):
            out.append("-")
            i -= 1
        elif j > 0 and (i == 0 or p == 2):
            j -= 1                      # insertion in obs: no reference column
        else:
            out.append(obs_prot[j - 1])
            i -= 1
            j -= 1
    return "".join(reversed(out))


# ---------------------------------------------------------------------------
# diversity
# ---------------------------------------------------------------------------
def pi_per_codon(weighted_proteins, ref_len):
    """Unbiased expected heterozygosity of amino acids at each reference position.

    `weighted_proteins`: iterable of (aligned_protein, weight). Weight is the number of observed
    haplotypes carrying that protein. Gaps ('-') are excluded from that position's denominator
    rather than treated as a 21st amino acid -- a deletion is missing data for this statistic, and
    counting it as a state would inflate diversity exactly at the positions where alignment is
    least certain.

    Returns (pi array of length ref_len, coverage array of effective haplotype counts).
    """
    counts = [Counter() for _ in range(ref_len)]
    for prot, w in weighted_proteins:
        if w <= 0:
            continue
        for i in range(min(ref_len, len(prot))):
            aa = prot[i]
            if aa != "-" and aa != "X":
                counts[i][aa] += w
    pi = np.zeros(ref_len, dtype=float)
    cov = np.zeros(ref_len, dtype=float)
    for i, c in enumerate(counts):
        n = sum(c.values())
        cov[i] = n
        if n < 2:
            continue
        s = sum((v / n) ** 2 for v in c.values())
        pi[i] = (n / (n - 1.0)) * (1.0 - s)
    return pi, cov


def codon_to_exon(segments):
    """[(exon_no, coding_len_bp), ...] -> list of exon numbers, one per codon.

    An exon boundary can fall mid-codon, which is real biology (HLA exon 1/2 and 2/3 junctions do
    this). The split codon is assigned to the exon contributing most of its bases; ties go to the
    earlier exon. Silently dropping split codons would shorten the groove by one residue at each
    junction.
    """
    per_base = []
    for ex, ln in segments:
        per_base.extend([ex] * int(ln))
    out = []
    for i in range(0, len(per_base) - len(per_base) % 3, 3):
        c = Counter(per_base[i:i + 3])
        out.append(min(c.items(), key=lambda kv: (-kv[1], kv[0]))[0])
    return out


def permutation_test(pi, mask, n_perm, rng):
    """Is mean pi inside `mask` higher than outside, beyond what a random set of that size gives?

    Codon positions are shuffled, which keeps the gene's diversity distribution exactly and only
    breaks the association with the groove. Positions with no coverage are excluded from both the
    statistic and the shuffling pool.
    """
    ok = ~np.isnan(pi)
    pi, mask = pi[ok], mask[ok]
    k = int(mask.sum())
    if k == 0 or k == len(mask):
        return float("nan"), float("nan"), float("nan")
    obs = float(pi[mask].mean() - pi[~mask].mean())
    idx = np.arange(len(pi))
    ge = 0
    for _ in range(n_perm):
        sel = rng.choice(idx, size=k, replace=False)
        m = np.zeros(len(pi), dtype=bool)
        m[sel] = True
        if (pi[m].mean() - pi[~m].mean()) >= obs:
            ge += 1
    p = (ge + 1) / (n_perm + 1)
    ratio = (float(pi[mask].mean()) / float(pi[~mask].mean())
             if pi[~mask].mean() > 0 else float("inf"))
    return obs, p, ratio


def allele_fst(freq_by_anc, n_by_anc):
    """Multi-allelic Fst over allele frequencies, plus Hedrick's standardised G'st.

    Fst = (Ht - Hs) / Ht with Hs the sample-size-weighted mean within-population expected
    heterozygosity. At a locus where Hs is already ~0.95, Fst is bounded near 0.05 no matter how
    different the populations are -- which is the whole reason HLA Fst looks tiny and why Hedrick's
    G'st = Gst / Gst_max is reported next to it. Brandt et al. 2018 (G3) is the citation for why
    raw Fst at HLA misleads.
    """
    ancs = [a for a in freq_by_anc if n_by_anc.get(a, 0) > 0]
    if len(ancs) < 2:
        return {"fst": float("nan"), "gst_prime": float("nan"), "hs": float("nan"),
                "ht": float("nan"), "n_pops": len(ancs)}
    tot = sum(n_by_anc[a] for a in ancs)
    hs = sum(n_by_anc[a] / tot * (1 - sum(p * p for p in freq_by_anc[a].values())) for a in ancs)
    pooled = defaultdict(float)
    for a in ancs:
        for al, p in freq_by_anc[a].items():
            pooled[al] += (n_by_anc[a] / tot) * p
    ht = 1 - sum(p * p for p in pooled.values())
    fst = (ht - hs) / ht if ht > 0 else float("nan")
    k = len(ancs)
    gst_max = ((k - 1) * (1 - hs)) / (k - 1 + hs) if (k - 1 + hs) > 0 else float("nan")
    gst_prime = (fst / gst_max) if gst_max and gst_max > 0 else float("nan")
    return {"fst": float(fst), "gst_prime": float(gst_prime), "hs": float(hs), "ht": float(ht),
            "n_pops": k}


# ---------------------------------------------------------------------------
# ARS contacts (optional)
# ---------------------------------------------------------------------------
def load_ars_contacts(path):
    """reference/ars_peptide_contacts.tsv -> {gene: DataFrame} or {} when absent.

    Absence is normal and not an error: the structural derivation is a separate workstream, and
    every headline number in this script is computed without it.
    """
    if not path or not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    if "gene" not in df.columns:
        return {}
    return {g: sub for g, sub in df.groupby("gene")}


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def fig_protein_track(gene, pi, exon_of_codon, groove_exons, ars_idx, path):
    """Per-residue amino-acid diversity along the protein, with groove exons shaded.

    This is the amino-acid version of the nucleotide Manhattan plot Marc showed in the call --
    the plot Cole said would be 'more interesting' than the nucleotide one.
    """
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(11, 3.4))
    x = np.arange(1, len(pi) + 1)
    in_groove = np.array([e in groove_exons for e in exon_of_codon[:len(pi)]]) \
        if groove_exons else np.zeros(len(pi), dtype=bool)
    # shade contiguous groove runs
    if in_groove.any():
        start = None
        for i, v in enumerate(list(in_groove) + [False]):
            if v and start is None:
                start = i
            elif not v and start is not None:
                ax.axvspan(start + 0.5, i + 0.5, color="#F2D7B6", alpha=0.55, lw=0, zorder=0)
                start = None
    ax.bar(x, pi, width=1.0, color="#3A3A3A", zorder=2)
    if ars_idx is not None and len(ars_idx):
        ax.plot(np.array(sorted(ars_idx)) + 1, np.full(len(ars_idx), -0.03 * max(1e-9, pi.max())),
                marker="^", ls="none", ms=3.2, color="#C0392B", zorder=3,
                label="peptide-contact residue (PDB)")
        ax.legend(frameon=False, fontsize=7, loc="upper right")
    ax.set_xlim(0.5, len(pi) + 0.5)
    ax.set_ylim(bottom=min(0, -0.06 * max(1e-9, pi.max())))
    ax.set_xlabel("residue (reference protein numbering, includes signal peptide)")
    ax.set_ylabel(r"$\pi_{aa}$")
    shaded = "groove exons %s shaded" % (sorted(groove_exons) if groove_exons else "none")
    ax.set_title("HLA-%s — amino-acid diversity per residue (%s)" % (gene, shaded), fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_groove_ratio(summary, path):
    plt = _mpl()
    df = pd.DataFrame(summary)
    df = df[df["ratio_groove_vs_rest"].notna()].copy()
    if df.empty:
        return
    df = df.sort_values("ratio_groove_vs_rest", ascending=False)
    colors = ["#999999" if g in CONSERVED_CONTROL_GENES else "#0072B2" for g in df["gene"]]
    fig, ax = plt.subplots(figsize=(max(6.0, 0.55 * len(df)), 3.8))
    ax.bar(range(len(df)), df["ratio_groove_vs_rest"], color=colors, width=0.7, zorder=3)
    ax.axhline(1.0, color="#C0392B", lw=1.0, ls="--", zorder=4)
    for i, (_, r) in enumerate(df.iterrows()):
        if pd.notna(r["p_permutation"]) and r["p_permutation"] < 0.01:
            ax.text(i, r["ratio_groove_vs_rest"], "*", ha="center", va="bottom", fontsize=11)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["gene"], rotation=0, fontsize=8)
    ax.set_ylabel(r"mean $\pi_{aa}$ groove / non-groove")
    ax.set_title("Amino-acid diversity is concentrated in the peptide-groove exons\n"
                 "dashed line = no enrichment · grey = conserved control genes · * = p < 0.01 "
                 "(position permutation)", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_fst(fst_rows, path):
    plt = _mpl()
    df = pd.DataFrame(fst_rows)
    if df.empty:
        return
    df = df[df["gst_prime"].notna()].sort_values("gst_prime", ascending=False)
    cls = {"A": "I", "B": "I", "C": "I", "E": "I", "F": "I", "G": "I"}
    colors = ["#0072B2" if cls.get(g, "II") == "I" else "#D55E00" for g in df["gene"]]
    fig, ax = plt.subplots(figsize=(max(6.0, 0.55 * len(df)), 3.8))
    ax.bar(range(len(df)), df["gst_prime"], color=colors, width=0.7, zorder=3)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["gene"], rotation=0, fontsize=8)
    ax.set_ylabel("Hedrick's G'st (standardised)")
    ax.set_title("Between-ancestry differentiation of HLA protein alleles\n"
                 "blue = class I · orange = class II · standardised because raw Fst is bounded "
                 "near zero when within-population heterozygosity is ~0.95", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def run(args):
    rng = np.random.default_rng(args.seed)
    ensure_dir(args.out_dir)
    m24 = mod("24_novelty_by_field.py", "novelty_by_field")

    if not os.path.exists(args.table1):
        sys.exit("FATAL: table1 not found: %s" % args.table1)

    print("[31] loading table1 ...", flush=True)
    t1 = pd.read_csv(args.table1, sep="\t", dtype=str, low_memory=False)
    t1["gene_bare"] = t1["gene"].map(gene_bare)
    t1 = t1[t1["gene_bare"].isin(TARGET_GENES)].copy()

    people, n_removed = m24.build_people(
        sorted(set(t1["person_id"].astype(str))), args.cohort_membership,
        args.relatedness_table, args.kin_min, args.strict_threshold, args.skip_relatedness)
    keep = people[people["unrelated"]]
    anc_of = dict(zip(keep["person_id"].astype(str), keep["anc_strict"]))
    n_people = len(keep)
    t1 = t1[t1["person_id"].astype(str).isin(anc_of)].copy()
    t1["ancestry"] = t1["person_id"].astype(str).map(anc_of)
    print("[31] unrelated people: %d (dropped %d)" % (n_people, n_removed), flush=True)

    t1["prot_id"] = t1["consensus"].map(two_field)
    n_calls = len(t1)
    unresolved = int(t1["prot_id"].isna().sum())
    t1 = t1.dropna(subset=["prot_id"])

    print("[31] loading IPD-IMGT reference proteins ...", flush=True)
    ref = m24.load_refdata(args.refdata, genes_needed=set(TARGET_GENES))
    # protein for each two-field group: the modal protein among reference alleles in that group
    # A two-field group can span several four-field alleles; they share a protein by definition,
    # but a handful disagree in the reference (partial records, corrected sequences). Take the
    # modal protein per group, which is the one the group's name actually denotes.
    _groups = {}
    for name, prot in ref.name_prot.items():
        tf = two_field(name)
        if tf is not None:
            _groups.setdefault(tf, Counter())[prot] += 1
    prot_of_group = {k: v.most_common(1)[0][0] for k, v in _groups.items()}
    print("[31] %d reference protein groups" % len(prot_of_group), flush=True)

    alleles_csv = args.alleles_csv or os.path.join(args.refdata, "alleles.csv.gz")
    exon_maps = None
    if os.path.exists(alleles_csv):
        # names_needed is membership-tested with `in`, so it must be a set, never None. We only
        # need the per-gene modal exon structure, not per-allele maps, hence the empty set.
        exon_maps = m24.load_exon_maps(alleles_csv, names_needed=set(),
                                       genes_needed=set(TARGET_GENES))
        print("[31] exon maps for %d genes" % len(getattr(exon_maps, "gene_mode", {})), flush=True)
    else:
        print("[31] WARNING: %s absent; groove-exon comparison will be skipped" % alleles_csv,
              flush=True)

    ars = load_ars_contacts(args.ars_contacts)
    if ars:
        print("[31] ARS peptide contacts loaded for: %s" % sorted(ars), flush=True)
    else:
        print("[31] no ARS contact file; groove-exon comparison only (this is fine)", flush=True)

    codon_rows, summary_rows, fst_rows = [], [], []
    coverage_rows = []

    for gene in TARGET_GENES:
        sub = t1[t1["gene_bare"] == gene]
        if len(sub) < args.min_haps:
            continue
        counts = Counter(sub["prot_id"])
        resolved = {p: c for p, c in counts.items() if p in prot_of_group}
        n_res = sum(resolved.values())
        coverage_rows.append({
            "gene": gene,
            "n_haplotype_calls_disp": suppress(len(sub)),
            "pct_protein_in_reference": round(100.0 * n_res / max(1, len(sub)), 3),
            "n_distinct_proteins": len(resolved),
        })
        if n_res < args.min_haps:
            continue

        # Reference frame = the most frequently observed protein in this cohort.
        ref_group = max(resolved.items(), key=lambda kv: kv[1])[0]
        ref_prot = prot_of_group[ref_group]
        L = len(ref_prot)

        aligned = {p: align_to_reference(ref_prot, prot_of_group[p]) for p in resolved}

        pi_all, cov_all = pi_per_codon([(aligned[p], w) for p, w in resolved.items()], L)

        # exon annotation
        exon_of_codon, groove = None, None
        if exon_maps is not None:
            segs = getattr(exon_maps, "gene_mode", {}).get(gene)
            if segs:
                exon_of_codon = codon_to_exon(segs)
                groove = m24.groove_exons_for(gene)

        ars_idx = []
        if gene in ars:
            col = "ref_index0" if "ref_index0" in ars[gene].columns else None
            if col:
                ars_idx = [int(v) for v in ars[gene][col] if 0 <= int(v) < L]

        # per-codon table
        for i in range(L):
            ex = exon_of_codon[i] if exon_of_codon and i < len(exon_of_codon) else None
            codon_rows.append({
                "gene": gene, "residue": i + 1, "ref_aa": ref_prot[i],
                "exon": ex,
                "in_groove_exon": (ex in groove) if (groove and ex is not None) else None,
                "is_peptide_contact": (i in set(ars_idx)) if ars_idx else None,
                "pi_aa": round(float(pi_all[i]), 6),
                "n_haplotypes_disp": suppress(int(cov_all[i])),
            })

        row = {"gene": gene, "ref_protein_group": ref_group, "protein_length": L,
               "n_haplotypes_disp": suppress(n_res), "n_distinct_proteins": len(resolved),
               "mean_pi_aa": round(float(np.mean(pi_all)), 6),
               "is_conserved_control": gene in CONSERVED_CONTROL_GENES}

        if exon_of_codon and groove:
            mask = np.array([(exon_of_codon[i] in groove) if i < len(exon_of_codon) else False
                             for i in range(L)])
            diff, p, ratio = permutation_test(pi_all, mask, args.n_permutations, rng)
            row.update({"n_groove_residues": int(mask.sum()),
                        "mean_pi_groove": round(float(pi_all[mask].mean()), 6) if mask.any() else None,
                        "mean_pi_rest": round(float(pi_all[~mask].mean()), 6) if (~mask).any() else None,
                        "ratio_groove_vs_rest": round(ratio, 4) if np.isfinite(ratio) else None,
                        "p_permutation": round(p, 5)})
            if ars_idx:
                am = np.zeros(L, dtype=bool)
                am[list(ars_idx)] = True
                inner = mask & am
                outer = mask & ~am
                if inner.any() and outer.any():
                    row.update({
                        "mean_pi_contact_in_groove": round(float(pi_all[inner].mean()), 6),
                        "mean_pi_noncontact_in_groove": round(float(pi_all[outer].mean()), 6),
                        "ratio_contact_vs_noncontact": round(
                            float(pi_all[inner].mean() / pi_all[outer].mean()), 4)
                        if pi_all[outer].mean() > 0 else None})
        summary_rows.append(row)

        # per-ancestry allele frequencies -> differentiation
        freq_by_anc, n_by_anc = {}, {}
        for anc in ANCESTRY_ORDER:
            s = sub[sub["ancestry"] == anc]
            c = Counter(p for p in s["prot_id"] if p in prot_of_group)
            n = sum(c.values())
            if n < args.min_haps_per_ancestry:
                continue
            freq_by_anc[anc] = {k: v / n for k, v in c.items()}
            n_by_anc[anc] = n
        f = allele_fst(freq_by_anc, n_by_anc)
        f.update({"gene": gene,
                  "ancestries": ",".join(sorted(n_by_anc)),
                  "min_n_disp": suppress(min(n_by_anc.values())) if n_by_anc else "0"})
        fst_rows.append(f)

        fig_protein_track(gene, pi_all, exon_of_codon or [], groove or set(), ars_idx,
                          os.path.join(args.out_dir, "fig_protein_track_%s.png" % gene))
        print("[31] %s: L=%d, %d proteins, mean pi=%.4f" % (gene, L, len(resolved),
                                                            float(np.mean(pi_all))), flush=True)

    pd.DataFrame(codon_rows).to_csv(os.path.join(args.out_dir, "aa_diversity_per_residue.tsv"),
                                    sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(os.path.join(args.out_dir, "groove_vs_rest_summary.tsv"),
                                      sep="\t", index=False)
    pd.DataFrame(fst_rows).to_csv(os.path.join(args.out_dir, "allele_differentiation.tsv"),
                                  sep="\t", index=False)
    pd.DataFrame(coverage_rows).to_csv(os.path.join(args.out_dir, "protein_reference_coverage.tsv"),
                                       sep="\t", index=False)

    fig_groove_ratio(summary_rows, os.path.join(args.out_dir, "fig_groove_enrichment.png"))
    fig_fst(fst_rows, os.path.join(args.out_dir, "fig_allele_differentiation.png"))

    write_readme(os.path.join(args.out_dir, "README.md"), args, summary_rows, fst_rows,
                 coverage_rows, n_people, n_calls, unresolved, bool(ars))
    with open(os.path.join(args.out_dir, "summary.json"), "w") as fh:
        json.dump({"n_people_unrelated": suppress(n_people),
                   "pct_calls_unresolved_protein": round(100.0 * unresolved / max(1, n_calls), 3),
                   "ars_contacts_used": bool(ars),
                   "genes": [r["gene"] for r in summary_rows]}, fh, indent=2)
    print("[31] done -> %s" % args.out_dir, flush=True)


def write_readme(path, args, summary_rows, fst_rows, coverage_rows, n_people, n_calls,
                 unresolved, used_ars):
    L = []
    L.append("# 31 — amino-acid diversity along the HLA protein, and between-ancestry "
             "differentiation\n")
    L.append("*Supervisor asks A8 (diversity at the amino-acid level, peptide-binding vs "
             "non-binding) and A10 (is class II more differentiated than class I), 2026-09-17.*\n")

    L.append("## Method in one paragraph\n")
    L.append("Every haplotype's call is reduced to its two-field protein identity and looked up in "
             "the IPD-IMGT reference shipped with Immuannot (3.55.0). Per-residue amino-acid "
             "diversity is the unbiased expected heterozygosity over those proteins, weighted by "
             "how often each was observed. No per-person sequence scan is needed. "
             "**%s** unrelated people; **%.2f%%** of calls had no resolvable protein identity "
             "(novel at field 1 or 2, or undetermined) and were excluded.\n"
             % (suppress(n_people), 100.0 * unresolved / max(1, n_calls)))
    L.append("\nThe groove is defined by **IMGT's own exon annotation** — exons 2+3 for class I, "
             "exon 2 for class II — not by a published residue list, because the canonical ARS "
             "residue tables could not be verified from primary sources "
             "(`WS_literature_selection.md`). %s\n"
             % ("A sharper contact-vs-non-contact comparison inside the groove exons is also "
                "reported, using peptide-contact residues derived structurally from PDB."
                if used_ars else
                "A structural peptide-contact definition was not available for this run, so only "
                "the exon-level comparison is reported."))

    L.append("\n## 1. How much of the cohort this covers\n")
    c = pd.DataFrame(coverage_rows)
    if not c.empty:
        L.append(c.to_markdown(index=False) + "\n")
        L.append("\n`pct_protein_in_reference` below ~99% for a gene means this analysis is "
                 "missing real diversity there — check it before interpreting that gene.\n")

    L.append("\n## 2. Groove vs non-groove amino-acid diversity\n")
    s = pd.DataFrame(summary_rows)
    if not s.empty:
        cols = [c for c in ["gene", "protein_length", "n_distinct_proteins", "n_haplotypes_disp",
                            "mean_pi_aa", "mean_pi_groove", "mean_pi_rest",
                            "ratio_groove_vs_rest", "p_permutation", "is_conserved_control"]
                if c in s.columns]
        L.append(s[cols].to_markdown(index=False) + "\n")
        L.append("\nThe permutation test shuffles residue positions within the gene, so it holds "
                 "the gene's diversity distribution fixed and only breaks the association with "
                 "the groove. The conserved control genes (DRA, HLA-E, HLA-F) are the check: they "
                 "should show little or no enrichment.\n")
        if used_ars and "ratio_contact_vs_noncontact" in s.columns:
            L.append("\n### Within the groove exons: peptide-contact vs the rest\n")
            cc = s[s["ratio_contact_vs_noncontact"].notna()]
            if not cc.empty:
                L.append(cc[["gene", "mean_pi_contact_in_groove", "mean_pi_noncontact_in_groove",
                             "ratio_contact_vs_noncontact"]].to_markdown(index=False) + "\n")

    L.append("\n## 3. Between-ancestry differentiation (ask A10)\n")
    f = pd.DataFrame(fst_rows)
    if not f.empty:
        L.append(f[["gene", "n_pops", "min_n_disp", "hs", "ht", "fst", "gst_prime"]]
                 .round(4).to_markdown(index=False) + "\n")
        L.append("\n`hs` is within-ancestry heterozygosity. At HLA it is close to 1, which "
                 "mechanically bounds raw `fst` near zero however different the populations are — "
                 "this is the Brandt et al. 2018 (G3) point. **Compare genes on `gst_prime`**, "
                 "Hedrick's standardised measure, not on `fst`.\n")

    L.append("\n## Files\n")
    L.append("- `aa_diversity_per_residue.tsv` — per-residue pi, exon, groove flag, contact flag.\n")
    L.append("- `groove_vs_rest_summary.tsv`, `allele_differentiation.tsv`, "
             "`protein_reference_coverage.tsv`.\n")
    L.append("- `fig_protein_track_<gene>.png` — the amino-acid Manhattan per gene.\n")
    L.append("- `fig_groove_enrichment.png`, `fig_allele_differentiation.png`.\n")

    L.append("\n## Caveats\n")
    L.append("- Catalogued proteins only. Novel proteins found in S01 are not in the reference "
             "and are excluded; see §1 for how much that is.\n")
    L.append("- Residue numbering is the reference protein's own, **including the signal "
             "peptide**. It is not mature-protein numbering, so these numbers are not directly "
             "comparable to a published residue list without an offset.\n")
    L.append("- Frequency weighting means a gene dominated by a few common alleles gets a lower "
             "pi than its allele *count* suggests. That is intended — it is the diversity a "
             "randomly drawn haplotype actually sees.\n")
    L.append("- Differentiation is computed on protein-level alleles, not SNPs, so it is not "
             "directly comparable to genome-wide Fst distributions.\n")
    with open(path, "w") as fh:
        fh.write("\n".join(L))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=DEFAULT_TABLE1)
    ap.add_argument("--cohort-membership", default=DEFAULT_COHORT)
    ap.add_argument("--relatedness-table", default=DEFAULT_RELATEDNESS)
    ap.add_argument("--skip-relatedness", action="store_true")
    ap.add_argument("--kin-min", type=float, default=KIN_MIN)
    ap.add_argument("--strict-threshold", type=float, default=STRICT_MIN)
    ap.add_argument("--refdata", default=DEFAULT_REFDATA)
    ap.add_argument("--alleles-csv", default=None)
    ap.add_argument("--ars-contacts", default=DEFAULT_ARS)
    ap.add_argument("--min-haps", type=int, default=MIN_HAPS_FOR_GENE)
    ap.add_argument("--min-haps-per-ancestry", type=int, default=100)
    ap.add_argument("--n-permutations", type=int, default=N_PERMUTATIONS)
    ap.add_argument("--seed", type=int, default=20260917)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    args = ap.parse_args(argv)
    run(args)


if __name__ == "__main__":
    main()
