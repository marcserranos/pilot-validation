#!/usr/bin/env python3
"""Linkage disequilibrium between HLA alleles at physically-phased cis gene pairs, within ancestry.

Supervisor ask (Cole, 2026-09-17 call, raised twice and repeated in the closing next-steps):

    "I would just want you to compute the LD between alleles within each ancestry and just kind of
     see what looks different. [...] the R2 measure would be the most normal way to do it. [...]
     we expect the LD patterns, just the correlation between any DQA1 and DQB1 allele, to be very
     different between the ancestries besides just having different frequencies."

That last clause is the whole analysis. Two ancestries can differ in an LD statistic for three
completely different reasons, and only the third is the biology Cole is asking about:

  1. **Different allele frequencies.** r^2 is bounded by the marginal frequencies. A pair of alleles
     at 2% cannot reach the r^2 a pair at 40% can. Comparing raw per-pair r^2 across ancestries
     therefore compares frequency spectra, not linkage.
  2. **Different allelic richness and different sample size.** Every multi-allelic association
     statistic built on a contingency table (Cramer's V, mutual information, chi-square) is biased
     upward when the table is large relative to n. AFR has both more alleles and, here, fewer
     haplotypes than EUR, so the naive comparison is confounded in *both* directions at once.
  3. **Genuinely different haplotype structure** -- which alleles travel together.

So this script reports three tiers, and the README says plainly which is which:

  * per-allele-pair r^2 and D' (what Cole asked for literally), with the marginal frequencies next
    to every value so the bound is visible;
  * multi-allelic summaries per gene pair per ancestry: Hedrick's (1987) multi-allelic D',
    Cramer's V, and normalised mutual information;
  * the same multi-allelic summaries **rarefied to a common haplotype count** across ancestries,
    with bootstrap CIs. This is the only cross-ancestry comparison we will defend in the paper.

Gene pairs. DQA1~DQB1 and DPA1~DPB1 are the ask. DRB1~DQB1 (the classic extended class II
haplotype) and HLA-B~HLA-C (very high LD, ~90 kb apart) are included as **positive controls**: if
the method does not recover strong LD there, the method is broken, not the biology. HLA-A~HLA-B
(~1.4 Mb) is a weak-LD control at the other end.

Phasing. Haplotypes come only from rows sharing `(person_id, hap, contig)` -- Immuannot's physical
assembly phase. A pair whose two genes landed on different contigs is NOT phased and is dropped;
the dropped fraction is itself reported per pair per ancestry, because it measures assembly
fragmentation across that interval (SCHEMA.md Table 2).

Usage (VM, full cohort):

    setsid nohup python3 -u scripts/hla_popgen/29_hla_ld_by_ancestry.py \\
        --out-dir ~/results/29_hla_ld < /dev/null & disown

Usage (local fixtures): scripts/hla_popgen/tests/test_hla_ld_by_ancestry.py
"""
import argparse
import importlib.util
import json
import math
import os
import random
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)


def _load_module(filename, modname):
    path = os.path.join(_THIS_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


# Reuse 24's cohort helpers so "unrelated" and "strict ancestry" mean exactly what they mean there.
_m24 = None


def m24():
    global _m24
    if _m24 is None:
        _m24 = _load_module("24_novelty_by_field.py", "novelty_by_field")
    return _m24


DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_TABLE1 = os.path.join(DEFAULT_OUTROOT, "hla_calls_rich.tsv")
DEFAULT_COHORT = os.path.join(DEFAULT_OUTROOT, "cohort_membership.tsv")
# ENVIRONMENT quirk #35: the old manual ~/mnt/aou-controlled mount is STALE and any process that
# touches it hangs in uninterruptible I/O -- Ctrl-C does nothing and the terminal is lost. The
# controlled bucket is auto-mounted under ~/workspace/ instead. Never point this at ~/mnt.
DEFAULT_RELATEDNESS = os.path.expanduser(
    "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/"
    "relatedness/samples_relatedness.tsv")
DEFAULT_OUT_DIR = os.path.expanduser("~/results/29_hla_ld")

SUPPRESS_BELOW = 20          # AoU disclosure: participant/haplotype counts 1-19 are written "<20"
KIN_MIN = 0.0442             # third degree or closer
STRICT_MIN = 0.9             # ancestry-probability floor for "strict" ancestry assignment

ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
ANCESTRY_COLORS = {"AFR": "#E69F00", "AMR": "#56B4E9", "EAS": "#009E73",
                   "EUR": "#0072B2", "MID": "#D55E00", "SAS": "#CC79A7"}

# (gene_a, gene_b, role). Order within a pair is fixed so allele_a/allele_b are stable.
GENE_PAIRS = [
    ("DQA1", "DQB1", "ask"),
    ("DPA1", "DPB1", "ask"),
    ("DRB1", "DQB1", "control_high"),
    ("B", "C", "control_high"),
    ("A", "B", "control_weak"),
]

MIN_ALLELE_HAPS = 20         # an allele must reach this in an ancestry to enter the pairwise table
MIN_HAPS_PER_ANCESTRY = 100  # below this, the ancestry is reported but flagged low-confidence
N_BOOTSTRAP = 500


# ---------------------------------------------------------------------------
# small utilities
# ---------------------------------------------------------------------------
def suppress(n, threshold=SUPPRESS_BELOW):
    """AoU small-cell rule: 0 is disclosable, 1-19 is not."""
    n = int(n)
    return "0" if n == 0 else ("<%d" % threshold if n < threshold else str(n))


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def two_field(name):
    """Truncate an allele call to `GENE*field1:field2`. Returns None if it cannot be resolved.

    A call carrying `new` at field 1 or 2 has no protein-level identity, so it cannot enter an LD
    table -- two different novel proteins would be silently merged under one label. Field-3/4
    novelty is fine: the protein is known, and that is the resolution we work at.
    """
    if not isinstance(name, str):
        return None
    tok = name.strip()
    if not tok or tok.upper() in {"NA", "NAN", "NONE", "UNDETERMINED", "."}:
        return None
    if "*" not in tok:
        return None
    gene, rest = tok.split("*", 1)
    fields = rest.split(":")
    if len(fields) < 2:
        return None
    if fields[0].strip().lower() == "new" or fields[1].strip().lower() == "new":
        return None
    return "%s*%s:%s" % (gene.replace("HLA-", ""), fields[0], fields[1])


# ---------------------------------------------------------------------------
# haplotype extraction
# ---------------------------------------------------------------------------
def extract_haplotypes(t1, gene_a, gene_b):
    """Physically-phased two-locus haplotypes for one gene pair.

    A haplotype exists only when both genes were annotated on the SAME contig of the SAME
    assembly haplotype, each with exactly one copy. Returns (haps_df, diagnostics dict).

    `haps_df` columns: person_id, hap, contig, allele_a, allele_b.
    Diagnostics counts every way a (person, hap) can fail to yield one, which is the
    assembly-fragmentation number SCHEMA.md Table 2 asks us to report.
    """
    need = {"person_id", "hap", "contig", "gene", "copy_index", "consensus"}
    missing = need - set(t1.columns)
    if missing:
        raise ValueError("table1 is missing columns: %s" % sorted(missing))

    bare = t1["gene"].astype(str).str.replace("^HLA-", "", regex=True)
    sub = t1[bare.isin([gene_a, gene_b])].copy()
    sub["gene_bare"] = bare[sub.index]

    diag = Counter()
    rows = []
    # Group by the cis key. Both genes present on one contig is the only phased case.
    for (pid, hap, contig), grp in sub.groupby(["person_id", "hap", "contig"], sort=False):
        a_rows = grp[grp["gene_bare"] == gene_a]
        b_rows = grp[grp["gene_bare"] == gene_b]
        if a_rows.empty or b_rows.empty:
            diag["contig_carries_only_one_gene"] += 1
            continue
        if len(a_rows) > 1 or len(b_rows) > 1:
            # >1 copy of a gene on one contig: real CNV or a mapping artifact. Either way the
            # pairing is ambiguous, so it cannot enter an LD table. Counted, not guessed at.
            diag["multi_copy_ambiguous"] += 1
            continue
        aa = two_field(a_rows.iloc[0]["consensus"])
        bb = two_field(b_rows.iloc[0]["consensus"])
        if aa is None or bb is None:
            diag["unresolved_2field_call"] += 1
            continue
        rows.append((pid, hap, contig, aa, bb))
        diag["phased_haplotypes"] += 1

    # How many (person, hap) assemblies carried both genes at all, on any contig? The difference
    # between that and `phased_haplotypes` is the fragmentation cost for this interval.
    per_assembly = sub.groupby(["person_id", "hap"])["gene_bare"].agg(lambda s: set(s))
    diag["assemblies_with_both_genes_somewhere"] = int(
        sum(1 for s in per_assembly if {gene_a, gene_b} <= s))

    haps = pd.DataFrame(rows, columns=["person_id", "hap", "contig", "allele_a", "allele_b"])
    return haps, dict(diag)


# ---------------------------------------------------------------------------
# LD statistics
# ---------------------------------------------------------------------------
def contingency(alleles_a, alleles_b):
    """Haplotype counts as a dense matrix plus the ordered allele labels."""
    a_labels = sorted(set(alleles_a))
    b_labels = sorted(set(alleles_b))
    ai = {a: i for i, a in enumerate(a_labels)}
    bi = {b: i for i, b in enumerate(b_labels)}
    M = np.zeros((len(a_labels), len(b_labels)), dtype=float)
    for a, b in zip(alleles_a, alleles_b):
        M[ai[a], bi[b]] += 1.0
    return M, a_labels, b_labels


def pairwise_ld(M, a_labels, b_labels, min_haps=MIN_ALLELE_HAPS):
    """Biallelic-collapsed r^2 and D' for every allele pair, in the standard 'this allele vs all
    others' framing -- the only way r^2, which is defined for two biallelic markers, is defined at
    a multi-allelic locus.

    Returns a list of dicts. `n_hap_ij` is the observed haplotype count; callers suppress it.
    """
    n = M.sum()
    if n <= 0:
        return []
    p_a = M.sum(axis=1) / n
    p_b = M.sum(axis=0) / n
    out = []
    for i, a in enumerate(a_labels):
        if M[i, :].sum() < min_haps:
            continue
        pi = p_a[i]
        for j, b in enumerate(b_labels):
            if M[:, j].sum() < min_haps:
                continue
            qj = p_b[j]
            pij = M[i, j] / n
            D = pij - pi * qj
            denom = pi * (1 - pi) * qj * (1 - qj)
            r2 = (D * D / denom) if denom > 0 else float("nan")
            if D >= 0:
                dmax = min(pi * (1 - qj), (1 - pi) * qj)
            else:
                dmax = min(pi * qj, (1 - pi) * (1 - qj))
            dprime = (D / dmax) if dmax > 0 else float("nan")
            # chi-square on the collapsed 2x2, 1 df -- exactly n * r^2.
            chi2 = n * r2 if not math.isnan(r2) else float("nan")
            out.append({
                "allele_a": a, "allele_b": b,
                "freq_a": float(pi), "freq_b": float(qj),
                "freq_hap": float(pij), "n_hap_ij": int(round(M[i, j])),
                "D": float(D), "Dprime": float(dprime), "r2": float(r2),
                "chi2_1df": float(chi2),
            })
    return out


def multiallelic_stats(M):
    """Three whole-table summaries of how strongly the two loci are associated.

    * `Dprime_multi` -- Hedrick (1987) multi-allelic D': sum_i sum_j p_i q_j |D'_ij|. Frequency
      weighted, in [0, 1], and the statistic HLA haplotype papers usually report.
    * `cramers_v`    -- sqrt(chi2 / (n * min(k-1, l-1))). Upward biased for large tables / small n.
    * `nmi`          -- mutual information normalised by min(H_a, H_b). Same bias direction.
    * `bias_corrected_cramers_v` -- Bergsma's bias correction, which removes most of the
      table-size inflation. Reported because the raw V across ancestries with different allele
      counts is not comparable, and that is precisely the comparison being asked for.
    """
    n = M.sum()
    k, l = M.shape
    if n <= 1 or k < 2 or l < 2:
        return {"n_haplotypes": int(n), "n_alleles_a": k, "n_alleles_b": l,
                "Dprime_multi": float("nan"), "cramers_v": float("nan"),
                "cramers_v_bc": float("nan"), "nmi": float("nan"),
                "chi2": float("nan"), "df": (k - 1) * (l - 1)}
    P = M / n
    p_a = P.sum(axis=1)
    p_b = P.sum(axis=0)
    E = np.outer(p_a, p_b)
    D = P - E

    with np.errstate(divide="ignore", invalid="ignore"):
        dmax_pos = np.minimum(np.outer(p_a, 1 - p_b), np.outer(1 - p_a, p_b))
        dmax_neg = np.minimum(np.outer(p_a, p_b), np.outer(1 - p_a, 1 - p_b))
        dmax = np.where(D >= 0, dmax_pos, dmax_neg)
        dprime = np.where(dmax > 0, np.abs(D) / dmax, 0.0)
    dprime_multi = float(np.sum(E * dprime))

    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = float(np.nansum(np.where(E > 0, (P - E) ** 2 / E, 0.0)) * n)
    phi2 = chi2 / n
    v = math.sqrt(phi2 / min(k - 1, l - 1)) if min(k - 1, l - 1) > 0 else float("nan")

    # Bergsma (2013) bias correction.
    phi2_bc = max(0.0, phi2 - (k - 1) * (l - 1) / (n - 1))
    k_bc = k - (k - 1) ** 2 / (n - 1)
    l_bc = l - (l - 1) ** 2 / (n - 1)
    denom_bc = min(k_bc - 1, l_bc - 1)
    v_bc = math.sqrt(phi2_bc / denom_bc) if denom_bc > 0 else float("nan")

    def H(p):
        p = p[p > 0]
        return float(-np.sum(p * np.log(p)))
    with np.errstate(divide="ignore", invalid="ignore"):
        mi = float(np.nansum(np.where(P > 0, P * np.log(np.where(E > 0, P / E, 1.0)), 0.0)))
    hmin = min(H(p_a), H(p_b))
    nmi = (mi / hmin) if hmin > 0 else float("nan")

    return {"n_haplotypes": int(n), "n_alleles_a": k, "n_alleles_b": l,
            "Dprime_multi": dprime_multi, "cramers_v": float(v), "cramers_v_bc": float(v_bc),
            "nmi": float(nmi), "chi2": chi2, "df": (k - 1) * (l - 1)}


def rarefied_stats(pairs, m, n_rep, rng):
    """Multi-allelic stats on `n_rep` subsamples of size `m` drawn without replacement.

    This is the cross-ancestry comparison. Sampling without replacement (a true subsample of
    haplotypes) rather than bootstrapping keeps the table size comparable to a real cohort of
    size m, which is the quantity being equalised.
    """
    if len(pairs) < m or m < 10:
        return None
    keys = ["Dprime_multi", "cramers_v", "cramers_v_bc", "nmi", "n_alleles_a", "n_alleles_b"]
    acc = defaultdict(list)
    idx = list(range(len(pairs)))
    for _ in range(n_rep):
        pick = rng.sample(idx, m)
        sub = [pairs[i] for i in pick]
        M, a_lab, b_lab = contingency([p[0] for p in sub], [p[1] for p in sub])
        st = multiallelic_stats(M)
        for k in keys:
            acc[k].append(st[k])
    out = {"rarefied_to": m, "n_reps": n_rep}
    for k in keys:
        v = np.array(acc[k], dtype=float)
        v = v[~np.isnan(v)]
        # Bergsma's correction is undefined when the rarefied table is as large as the subsample
        # (k_bc - 1 collapses to 0). That is a real "not estimable at this n", so it stays NaN
        # rather than being filled in with a number that would look like an estimate.
        if v.size == 0:
            out[k + "_mean"] = out[k + "_lo"] = out[k + "_hi"] = float("nan")
            continue
        out[k + "_mean"] = float(np.mean(v))
        out[k + "_lo"] = float(np.percentile(v, 2.5))
        out[k + "_hi"] = float(np.percentile(v, 97.5))
    return out


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def fig_multiallelic(rare_rows, path):
    """Rarefied multi-allelic D' per ancestry, one group of bars per gene pair, with 95% CIs."""
    plt = _mpl()
    df = pd.DataFrame(rare_rows)
    if df.empty:
        return
    pairs = [p for p in df["pair"].unique()]
    ancs = [a for a in ANCESTRY_ORDER if a in set(df["ancestry"])]
    fig, ax = plt.subplots(figsize=(1.6 * len(pairs) + 3.0, 4.2))
    width = 0.8 / max(1, len(ancs))
    for ai, anc in enumerate(ancs):
        xs, ys, lo, hi = [], [], [], []
        for pi, pair in enumerate(pairs):
            r = df[(df["pair"] == pair) & (df["ancestry"] == anc)]
            if r.empty:
                continue
            r = r.iloc[0]
            xs.append(pi - 0.4 + width * (ai + 0.5))
            ys.append(r["Dprime_multi_mean"])
            lo.append(max(0.0, r["Dprime_multi_mean"] - r["Dprime_multi_lo"]))
            hi.append(max(0.0, r["Dprime_multi_hi"] - r["Dprime_multi_mean"]))
        if not xs:
            continue
        ax.bar(xs, ys, width=width * 0.92, color=ANCESTRY_COLORS.get(anc, "#777777"),
               label=anc, yerr=[lo, hi], capsize=2, error_kw={"lw": 0.8})
    ax.set_xticks(range(len(pairs)))
    ax.set_xticklabels(pairs, rotation=0)
    ax.set_ylabel("multi-allelic D' (Hedrick 1987)")
    ax.set_ylim(0, 1.02)
    n_rare = int(df["rarefied_to"].iloc[0]) if "rarefied_to" in df.columns else 0
    ax.set_title("Two-locus LD within ancestry, rarefied to %d phased haplotypes per ancestry\n"
                 "(equal n and equal opportunity for allelic richness; 95%% CI over subsamples)"
                 % n_rare, fontsize=9)
    ax.legend(frameon=False, ncol=len(ancs), fontsize=8, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_r2_heatmaps(pairwise, pair_name, path, top_k=12):
    """r^2 heatmap per ancestry for the top-k alleles of each gene, shared allele ordering."""
    plt = _mpl()
    df = pd.DataFrame(pairwise)
    if df.empty:
        return
    df = df[df["pair"] == pair_name]
    if df.empty:
        return
    ancs = [a for a in ANCESTRY_ORDER if a in set(df["ancestry"])]
    if not ancs:
        return
    # One shared allele ordering across panels, chosen on pooled frequency, so the panels are
    # directly comparable cell by cell.
    a_top = (df.groupby("allele_a")["freq_a"].max().sort_values(ascending=False)
             .head(top_k).index.tolist())
    b_top = (df.groupby("allele_b")["freq_b"].max().sort_values(ascending=False)
             .head(top_k).index.tolist())
    ncol = min(3, len(ancs))
    nrow = int(math.ceil(len(ancs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.8 * nrow), squeeze=False)
    for idx, anc in enumerate(ancs):
        ax = axes[idx // ncol][idx % ncol]
        sub = df[df["ancestry"] == anc]
        M = np.full((len(a_top), len(b_top)), np.nan)
        look = {(r["allele_a"], r["allele_b"]): r["r2"] for _, r in sub.iterrows()}
        for i, a in enumerate(a_top):
            for j, b in enumerate(b_top):
                if (a, b) in look:
                    M[i, j] = look[(a, b)]
        im = ax.imshow(M, vmin=0, vmax=1, cmap="magma_r", aspect="auto")
        ax.set_xticks(range(len(b_top)))
        ax.set_xticklabels(b_top, rotation=90, fontsize=6)
        ax.set_yticks(range(len(a_top)))
        ax.set_yticklabels(a_top, fontsize=6)
        ax.set_title(anc, fontsize=10, color=ANCESTRY_COLORS.get(anc, "#333333"))
    for idx in range(len(ancs), nrow * ncol):
        axes[idx // ncol][idx % ncol].axis("off")
    fig.suptitle("%s: pairwise r$^2$ (allele vs all others), top %d alleles by frequency\n"
                 "blank = allele below the %d-haplotype floor in that ancestry"
                 % (pair_name, top_k, MIN_ALLELE_HAPS), fontsize=10)
    fig.colorbar(im, ax=axes, shrink=0.6, label="r$^2$")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_pairable(diag_rows, path):
    """What fraction of assemblies carrying both genes actually had them on one contig."""
    plt = _mpl()
    df = pd.DataFrame(diag_rows)
    if df.empty or "pct_phased" not in df.columns:
        return
    pairs = list(df["pair"].unique())
    fig, ax = plt.subplots(figsize=(1.4 * len(pairs) + 2.5, 3.4))
    ax.bar(range(len(pairs)), [df[df["pair"] == p]["pct_phased"].iloc[0] for p in pairs],
           color="#4C72B0", width=0.6)
    ax.set_xticks(range(len(pairs)))
    ax.set_xticklabels(pairs)
    ax.set_ylim(0, 102)
    ax.set_ylabel("% of assemblies with both genes\nthat had them on one contig")
    ax.set_title("Physical phasing yield per interval (assembly fragmentation)", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def write_readme(path, args, diag_rows, multi_rows, rare_rows, top_r2, n_people):
    L = []
    L.append("# 29 — LD between HLA alleles at phased cis gene pairs, within ancestry\n")
    L.append("*Supervisor ask A3 (Cole, 2026-09-17): \"compute the LD between alleles within each "
             "ancestry and see what looks different\", r^2, for DQA1–DQB1 and DPA1–DPB1.*\n")
    L.append("## What was computed\n")
    L.append("Haplotypes are **physically phased**: both genes annotated on the same contig of the "
             "same assembly haplotype, one copy each, both resolvable to two fields (protein "
             "level). People are restricted to a greedy unrelated set (kinship < %.4f) and, for "
             "the per-ancestry tables, to strict ancestry (probability >= %.2f).\n"
             % (args.kin_min, args.strict_threshold))
    L.append("- unrelated people used: **%s**\n" % suppress(n_people))
    L.append("\n### Why three tiers of statistic\n")
    L.append("r^2 between two alleles is bounded by their marginal frequencies, and every "
             "multi-allelic association statistic is inflated when there are many alleles relative "
             "to sample size. Ancestries differ in both. So a raw cross-ancestry comparison of "
             "either quantity measures the frequency spectrum and the allele count, not linkage. "
             "The defensible comparison is the **rarefied** one: every ancestry subsampled to the "
             "same number of phased haplotypes, statistics recomputed, 95% interval over "
             "subsamples.\n")

    L.append("\n## 1. Phasing yield per interval\n")
    d = pd.DataFrame(diag_rows)
    if not d.empty:
        show = d[["pair", "assemblies_with_both_genes_somewhere_disp", "phased_haplotypes_disp",
                  "pct_phased", "multi_copy_ambiguous_disp", "unresolved_2field_call_disp"]]
        show = show.rename(columns={"assemblies_with_both_genes_somewhere_disp": "both genes present",
                                    "phased_haplotypes_disp": "phased (same contig)",
                                    "pct_phased": "% phased",
                                    "multi_copy_ambiguous_disp": "dropped: >1 copy",
                                    "unresolved_2field_call_disp": "dropped: unresolved call"})
        L.append(show.to_markdown(index=False) + "\n")
        L.append("\nA low percentage here is assembly fragmentation across that interval, not a "
                 "typing failure. DPA1–DPB1 sit ~10 kb apart and should be near-100%; the class I "
                 "pairs span far more sequence and should be lower.\n")

    L.append("\n## 2. Multi-allelic LD per ancestry (raw, NOT comparable across ancestries)\n")
    m = pd.DataFrame(multi_rows)
    if not m.empty:
        show = m[["pair", "ancestry", "n_haplotypes_disp", "n_alleles_a", "n_alleles_b",
                  "Dprime_multi", "cramers_v", "cramers_v_bc", "nmi"]].copy()
        for c in ["Dprime_multi", "cramers_v", "cramers_v_bc", "nmi"]:
            show[c] = show[c].map(lambda v: "%.3f" % v if pd.notna(v) else "NA")
        L.append(show.to_markdown(index=False) + "\n")

    L.append("\n## 3. Rarefied multi-allelic LD (the cross-ancestry comparison)\n")
    r = pd.DataFrame(rare_rows)
    if not r.empty:
        show = r[["pair", "ancestry", "rarefied_to", "Dprime_multi_mean", "Dprime_multi_lo",
                  "Dprime_multi_hi", "cramers_v_bc_mean", "n_alleles_a_mean", "n_alleles_b_mean"]].copy()
        for c in ["Dprime_multi_mean", "Dprime_multi_lo", "Dprime_multi_hi", "cramers_v_bc_mean"]:
            show[c] = show[c].map(lambda v: "%.3f" % v if pd.notna(v) else "NA")
        for c in ["n_alleles_a_mean", "n_alleles_b_mean"]:
            show[c] = show[c].map(lambda v: "%.1f" % v if pd.notna(v) else "NA")
        L.append(show.to_markdown(index=False) + "\n")
        L.append("\nRead this table as: *at equal sample size*, does the same pair of loci travel "
                 "together more tightly in one ancestry than another? Non-overlapping CIs are the "
                 "claim; overlapping CIs are not.\n")

    L.append("\n## 4. Strongest individual allele pairs (what Cole asked for literally)\n")
    t = pd.DataFrame(top_r2)
    if not t.empty:
        show = t[["pair", "ancestry", "allele_a", "allele_b", "freq_a", "freq_b", "freq_hap",
                  "r2", "Dprime"]].copy()
        for c in ["freq_a", "freq_b", "freq_hap", "r2", "Dprime"]:
            show[c] = show[c].map(lambda v: "%.3f" % v if pd.notna(v) else "NA")
        L.append(show.to_markdown(index=False) + "\n")
        L.append("\nThe marginal frequencies are printed next to every r^2 on purpose: r^2 is "
                 "bounded by them, so a small r^2 at freq 0.02 can be *stronger* linkage than a "
                 "larger r^2 at freq 0.4.\n")

    L.append("\n## Files\n")
    L.append("- `ld_pairwise.tsv` — every allele pair passing the %d-haplotype floor, per "
             "ancestry, with r^2, D', D and marginals.\n" % MIN_ALLELE_HAPS)
    L.append("- `ld_multiallelic.tsv` — whole-table statistics per pair per ancestry (raw).\n")
    L.append("- `ld_rarefied.tsv` — the same, rarefied to a common haplotype count.\n")
    L.append("- `phasing_yield.tsv` — physical-phasing yield per interval.\n")
    L.append("- `haplotype_freqs.tsv` — two-locus haplotype frequencies per ancestry "
             "(counts below %d written `<%d`).\n" % (SUPPRESS_BELOW, SUPPRESS_BELOW))
    L.append("- `fig_ld_multiallelic.png`, `fig_r2_<pair>.png`, `fig_phasing_yield.png`.\n")
    L.append("\n## Caveats\n")
    L.append("- Physical phase is Immuannot's assembly phase. S01's QC (script 26) found zero "
             "phase switches over 3,021 testable transitions with 100% power on injected "
             "switches, which is what licenses treating these as true haplotypes.\n")
    L.append("- Two-field resolution. Calls novel at field 1 or 2 are dropped, not merged.\n")
    L.append("- Strict ancestry shrinks AMR and MID most; check `n_haplotypes` before reading "
             "anything into those rows.\n")
    with open(path, "w") as fh:
        fh.write("\n".join(L))


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def run(args):
    rng = random.Random(args.seed)
    ensure_dir(args.out_dir)

    if not os.path.exists(args.table1):
        sys.exit("FATAL: table1 not found: %s" % args.table1)

    print("[29] loading table1 ...", flush=True)
    t1 = pd.read_csv(args.table1, sep="\t", dtype=str, low_memory=False)
    t1["copy_index"] = pd.to_numeric(t1.get("copy_index", 1), errors="coerce").fillna(1).astype(int)

    # build_people -> (DataFrame[person_id, anc_pred, anc_strict, unrelated], n_removed)
    people, n_removed = m24().build_people(
        sorted(set(t1["person_id"].astype(str))), args.cohort_membership,
        args.relatedness_table, args.kin_min, args.strict_threshold, args.skip_relatedness)
    keep = people[people["unrelated"]]
    anc_col = "anc_strict" if args.ancestry_scheme == "strict" else "anc_pred"
    anc_of = dict(zip(keep["person_id"].astype(str), keep[anc_col]))
    n_people = len(keep)
    print("[29] unrelated people: %d (dropped %d relatives); ancestry column %s"
          % (n_people, n_removed, anc_col), flush=True)

    t1 = t1[t1["person_id"].astype(str).isin(anc_of)]

    diag_rows, multi_rows, rare_rows, pairwise_rows, hapfreq_rows = [], [], [], [], []

    for gene_a, gene_b, role in GENE_PAIRS:
        pair_name = "%s~%s" % (gene_a, gene_b)
        print("[29] %s ..." % pair_name, flush=True)
        haps, diag = extract_haplotypes(t1, gene_a, gene_b)
        both = diag.get("assemblies_with_both_genes_somewhere", 0)
        phased = diag.get("phased_haplotypes", 0)
        diag_rows.append({
            "pair": pair_name, "role": role,
            "assemblies_with_both_genes_somewhere": both,
            "assemblies_with_both_genes_somewhere_disp": suppress(both),
            "phased_haplotypes": phased, "phased_haplotypes_disp": suppress(phased),
            "pct_phased": round(100.0 * phased / both, 2) if both else float("nan"),
            "multi_copy_ambiguous_disp": suppress(diag.get("multi_copy_ambiguous", 0)),
            "unresolved_2field_call_disp": suppress(diag.get("unresolved_2field_call", 0)),
            "contig_carries_only_one_gene_disp": suppress(
                diag.get("contig_carries_only_one_gene", 0)),
        })
        if haps.empty:
            continue
        haps["ancestry"] = haps["person_id"].astype(str).map(anc_of)

        by_anc = {}
        for anc, grp in haps.groupby("ancestry"):
            if not isinstance(anc, str) or anc.upper() not in ANCESTRY_ORDER:
                continue
            by_anc[anc.upper()] = list(zip(grp["allele_a"], grp["allele_b"]))

        for anc, pairs in sorted(by_anc.items()):
            M, a_lab, b_lab = contingency([p[0] for p in pairs], [p[1] for p in pairs])
            st = multiallelic_stats(M)
            st.update({"pair": pair_name, "role": role, "ancestry": anc,
                       "n_haplotypes_disp": suppress(st["n_haplotypes"]),
                       "low_confidence": st["n_haplotypes"] < MIN_HAPS_PER_ANCESTRY})
            multi_rows.append(st)

            for rec in pairwise_ld(M, a_lab, b_lab, args.min_allele_haps):
                rec.update({"pair": pair_name, "ancestry": anc,
                            "n_hap_ij_disp": suppress(rec.pop("n_hap_ij"))})
                pairwise_rows.append(rec)

            n = M.sum()
            for i, a in enumerate(a_lab):
                for j, b in enumerate(b_lab):
                    if M[i, j] <= 0:
                        continue
                    hapfreq_rows.append({"pair": pair_name, "ancestry": anc,
                                         "haplotype": "%s~%s" % (a, b),
                                         "n_haplotypes_disp": suppress(int(M[i, j])),
                                         "freq": round(float(M[i, j] / n), 5)})

        # Rarefaction target: the smallest ancestry that clears the floor, so no ancestry is
        # extrapolated. Ancestries below the target are excluded from the rarefied comparison
        # rather than being compared at a different n.
        sizes = [len(v) for v in by_anc.values() if len(v) >= MIN_HAPS_PER_ANCESTRY]
        if not sizes:
            continue
        target = args.rarefy_to or min(sizes)
        for anc, pairs in sorted(by_anc.items()):
            rs = rarefied_stats(pairs, target, args.n_bootstrap, rng)
            if rs is None:
                continue
            rs.update({"pair": pair_name, "role": role, "ancestry": anc})
            rare_rows.append(rs)

    # ---- write ----
    def w(rows, name, cols=None):
        df = pd.DataFrame(rows)
        if cols:
            df = df[[c for c in cols if c in df.columns]]
        df.to_csv(os.path.join(args.out_dir, name), sep="\t", index=False, na_rep="NA")
        return df

    w(diag_rows, "phasing_yield.tsv")
    # Raw haplotype counts never leave the VM; only the suppressed display column does.
    w([{k: v for k, v in r.items() if k != "n_haplotypes"} for r in multi_rows],
      "ld_multiallelic.tsv")
    w(rare_rows, "ld_rarefied.tsv")
    w(pairwise_rows, "ld_pairwise.tsv")
    w(hapfreq_rows, "haplotype_freqs.tsv")

    top_r2 = []
    if pairwise_rows:
        pw = pd.DataFrame(pairwise_rows)
        for pair_name in pw["pair"].unique():
            sub = pw[pw["pair"] == pair_name].sort_values("r2", ascending=False).head(5)
            top_r2.extend(sub.to_dict("records"))

    fig_multiallelic(rare_rows, os.path.join(args.out_dir, "fig_ld_multiallelic.png"))
    fig_pairable(diag_rows, os.path.join(args.out_dir, "fig_phasing_yield.png"))
    for gene_a, gene_b, _ in GENE_PAIRS:
        pair_name = "%s~%s" % (gene_a, gene_b)
        fig_r2_heatmaps(pairwise_rows, pair_name,
                        os.path.join(args.out_dir, "fig_r2_%s_%s.png" % (gene_a, gene_b)))

    write_readme(os.path.join(args.out_dir, "README.md"), args, diag_rows,
                 [dict(r) for r in multi_rows], rare_rows, top_r2, n_people)

    with open(os.path.join(args.out_dir, "summary.json"), "w") as fh:
        json.dump({"n_people_unrelated": suppress(n_people),
                   "pairs": [f"{a}~{b}" for a, b, _ in GENE_PAIRS],
                   "rarefy_to": args.rarefy_to,
                   "strict_threshold": args.strict_threshold,
                   "kin_min": args.kin_min}, fh, indent=2)
    print("[29] done -> %s" % args.out_dir, flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=DEFAULT_TABLE1)
    ap.add_argument("--cohort-membership", default=DEFAULT_COHORT)
    ap.add_argument("--relatedness-table", default=DEFAULT_RELATEDNESS)
    ap.add_argument("--skip-relatedness", action="store_true")
    ap.add_argument("--kin-min", type=float, default=KIN_MIN)
    ap.add_argument("--strict-threshold", type=float, default=STRICT_MIN,
                    help="ancestry probability floor; the call asked for a stricter value (0.98) "
                         "as a sensitivity check")
    ap.add_argument("--ancestry-scheme", choices=["strict", "pred"], default="strict",
                    help="strict = ancestry probability >= --strict-threshold (the call asked for "
                         "this); pred = AoU's point prediction, admixture ignored")
    ap.add_argument("--min-allele-haps", type=int, default=MIN_ALLELE_HAPS)
    ap.add_argument("--rarefy-to", type=int, default=None,
                    help="common haplotype count for the cross-ancestry comparison; default is "
                         "the smallest ancestry clearing the floor")
    ap.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)
    ap.add_argument("--seed", type=int, default=20260917)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    args = ap.parse_args(argv)
    run(args)


if __name__ == "__main__":
    main()
