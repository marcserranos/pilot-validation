#!/usr/bin/env python3
"""Gene deletions, duplications and copy number across the MHC, from phased long-read assemblies.

Supervisor asks (Cole, 2026-09-17 call):

    A4: "did you look at any deletions or duplications? [...] Let's do that. It could be pretty
         crazy. [...] I know in the Immuannot paper they find duplications and deletions in
         1000 Genomes."
    A5: "did you look at the KIR alleles at all? It should be called if you just ran it straight up."

Why this is not simply "count the missing genes"
------------------------------------------------
Immuannot annotates a gene when it finds homology on a contig. A gene absent from a haplotype's
annotation has three possible causes, and only the first is biology:

  1. the haplotype genuinely lacks the gene (a real deletion, e.g. DRB3/4/5, C4, DRB1 haplotype
     groups, and the known MICA/MICB and HLA-DPA2/DPB2 deletion polymorphisms);
  2. the assembly broke: the gene's position fell off a contig end, so nothing could be annotated;
  3. detection failed: divergent sequence below Immuannot's homology filters.

Cause 2 is not rare here -- S01 and SCHEMA.md Table 2 both show real fragmentation across the
region. So this script never calls a deletion from absence alone. A gene is called **deleted** only
when a single contig carries a gene on each side of it, in canonical physical order, and the gene
itself is not on that contig. That is a bridged absence: the sequence was assembled through the
gene's position and the gene was not there.

Controls, which are the point
-----------------------------
* **Positive control, DRB3/4/5.** Their presence is determined by the DRB1 haplotype group
  (DR52 = DRB1*03/11/12/13/14 -> DRB3; DR53 = DRB1*04/07/09 -> DRB4; DR51 = DRB1*15/16 -> DRB5;
  DR1/DR8/DR10 -> none). This is textbook, it is copy-number variation, and we know the answer in
  advance. If the method does not recover it, no other number in this script can be trusted.
* **Negative control, HLA-A/B/C/DRA.** Essentially never deleted in a viable genome. A nonzero
  bridged-deletion rate at these genes is the method's false-positive rate, and it is reported as
  such rather than interpreted.
* **C4A/C4B.** Genuine, well-characterised copy-number variation (0-4 copies per haplotype) with a
  long/short distinction that Immuannot already records in `c4_size` (SCHEMA.md Table 1). This is
  the analysis where long reads have a real advantage over short reads.
* **KIR.** Expected to be **entirely absent**: the KIR cluster is on chr19, outside the chr6 trim
  window used to build these assemblies (SCHEMA.md gene classification). Any KIR call would mean a
  trim bug and is reported loudly. This answers A5 directly rather than by assumption.

Usage (VM, full cohort):

    setsid nohup python3 -u scripts/hla_popgen/30_hla_structural_variation.py \\
        --out-dir ~/results/30_hla_sv < /dev/null & disown

Usage (local fixtures): scripts/hla_popgen/tests/test_hla_structural_variation.py
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
# ENVIRONMENT quirk #35: the old manual ~/mnt/aou-controlled mount is STALE and any process that
# touches it hangs in uninterruptible I/O -- Ctrl-C does nothing and the terminal is lost. The
# controlled bucket is auto-mounted under ~/workspace/ instead. Never point this at ~/mnt.
DEFAULT_RELATEDNESS = os.path.expanduser(
    "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/"
    "relatedness/samples_relatedness.tsv")
DEFAULT_OUT_DIR = os.path.expanduser("~/results/30_hla_sv")

SUPPRESS_BELOW = 20
KIN_MIN = 0.0442
STRICT_MIN = 0.9
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
ANCESTRY_COLORS = {"AFR": "#E69F00", "AMR": "#56B4E9", "EAS": "#009E73",
                   "EUR": "#0072B2", "MID": "#D55E00", "SAS": "#CC79A7"}

NEGATIVE_CONTROL_GENES = ["A", "B", "C", "DRA", "DQA1", "DQB1", "DPA1", "DPB1", "DRB1"]
POSITIVE_CONTROL_GENES = ["DRB3", "DRB4", "DRB5"]

# Textbook DRB1 group -> expected second DRB locus. Keyed on the DRB1 first field.
DRB1_GROUP_EXPECTATION = {
    "03": "DRB3", "11": "DRB3", "12": "DRB3", "13": "DRB3", "14": "DRB3",
    "04": "DRB4", "07": "DRB4", "09": "DRB4",
    "15": "DRB5", "16": "DRB5",
    "01": None, "08": None, "10": None,
}


def suppress(n, threshold=SUPPRESS_BELOW):
    n = int(n)
    return "0" if n == 0 else ("<%d" % threshold if n < threshold else str(n))


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def gene_bare(g):
    return str(g).replace("HLA-", "")


def first_field(consensus):
    """First allele field, e.g. `HLA-DRB1*15:01:01` -> `15`. None if unresolvable."""
    if not isinstance(consensus, str) or "*" not in consensus:
        return None
    f = consensus.split("*", 1)[1].split(":")[0].strip()
    if not f or f.lower() == "new":
        return None
    return f.zfill(2) if f.isdigit() else f


# ---------------------------------------------------------------------------
# bridged-absence deletion calling
# ---------------------------------------------------------------------------
def bridged_absences(t1, gene_order, genes_of_interest):
    """Call a gene deleted on a haplotype only when one contig spans its position.

    For each (person_id, hap, contig), take the set of annotated genes and their canonical ranks.
    A gene g with rank r is BRIDGED-ABSENT on that contig when the contig carries at least one
    gene with rank < r and at least one gene with rank > r, and g itself is absent. A gene that is
    absent while the contig has no flanking gene on one side is UNBRIDGED -- the assembly simply
    ends there -- and is counted separately, never as a deletion.

    Returns (calls_df, diag). calls_df has one row per (person_id, hap, gene, status) with
    status in {present, deleted_bridged, absent_unbridged}. A gene absent from *every* contig of a
    haplotype but bridged on one of them counts as deleted once for that haplotype.
    """
    diag = Counter()
    per_hap = defaultdict(dict)  # (pid, hap) -> gene -> status, 'present' wins
    seen_haps = set()

    cols = ["person_id", "hap", "contig", "gene_bare"]
    for (pid, hap, contig), grp in t1[cols].groupby(["person_id", "hap", "contig"], sort=False):
        seen_haps.add((pid, hap))
        present = set(grp["gene_bare"])
        ranks = [gene_order[g] for g in present if g in gene_order]
        if not ranks:
            diag["contig_with_no_ordered_gene"] += 1
            continue
        lo, hi = min(ranks), max(ranks)
        for g in present:
            per_hap[(pid, hap)][g] = "present"
        for g in genes_of_interest:
            if g in present or g not in gene_order:
                continue
            r = gene_order[g]
            if lo < r < hi:
                # Bridged: the contig was assembled through this gene's position.
                per_hap[(pid, hap)].setdefault(g, "deleted_bridged")
                diag["bridged_absence"] += 1
            else:
                per_hap[(pid, hap)].setdefault(g, "absent_unbridged")
                diag["unbridged_absence"] += 1

    rows = []
    for (pid, hap), gmap in per_hap.items():
        for g, status in gmap.items():
            if g in genes_of_interest:
                rows.append((pid, hap, g, status))
    calls = pd.DataFrame(rows, columns=["person_id", "hap", "gene", "status"])
    diag["haplotypes_seen"] = len(seen_haps)
    return calls, dict(diag)


def duplications(t1, genes_of_interest):
    """Extra copies of a gene on one contig: `copy_index` > 1 at the same (person, hap, contig).

    SCHEMA.md is explicit that copy_index > 1 means a real second mapping cluster, not an
    artifact -- but a segmental duplication can still produce a spurious second cluster, so this
    is reported as a *candidate* rate per gene, next to the negative-control genes where the rate
    should be ~0.
    """
    g = (t1.groupby(["person_id", "hap", "contig", "gene_bare"])["copy_index"]
           .max().reset_index(name="max_copy"))
    g = g[g["gene_bare"].isin(genes_of_interest)]
    dup = (g.groupby("gene_bare")
             .agg(n_haplotype_contigs=("max_copy", "size"),
                  n_multi_copy=("max_copy", lambda s: int((s > 1).sum())),
                  max_copy_seen=("max_copy", "max"))
             .reset_index().rename(columns={"gene_bare": "gene"}))
    dup["pct_multi_copy"] = 100.0 * dup["n_multi_copy"] / dup["n_haplotype_contigs"].clip(lower=1)
    return dup.sort_values("pct_multi_copy", ascending=False)


def drb_expectation_check(t1, calls):
    """The positive control. For each haplotype with a resolvable DRB1 first field, does the
    presence of DRB3/DRB4/DRB5 match the textbook DR51/52/53 expectation?

    Concordance here is the method's accuracy on a copy-number call whose answer is known
    independently of our data.
    """
    drb1 = t1[t1["gene_bare"] == "DRB1"].copy()
    drb1["ff"] = drb1["consensus"].map(first_field)
    drb1 = drb1.dropna(subset=["ff"])
    # one DRB1 per haplotype; ambiguous multi-copy haplotypes are excluded
    per_hap = drb1.groupby(["person_id", "hap"])["ff"].agg(lambda s: s.iloc[0] if len(set(s)) == 1
                                                           else None).dropna()

    present = defaultdict(set)
    for r in calls[calls["status"] == "present"].itertuples(index=False):
        if r.gene in POSITIVE_CONTROL_GENES:
            present[(r.person_id, r.hap)].add(r.gene)

    rows = []
    for (pid, hap), ff in per_hap.items():
        if ff not in DRB1_GROUP_EXPECTATION:
            continue
        expected = DRB1_GROUP_EXPECTATION[ff]
        obs = present.get((pid, hap), set())
        if expected is None:
            ok = len(obs) == 0
        else:
            ok = expected in obs
        rows.append({"drb1_group": ff, "expected": expected or "none",
                     "observed": ",".join(sorted(obs)) if obs else "none", "concordant": ok})
    df = pd.DataFrame(rows)
    if df.empty:
        return df, df
    summ = (df.groupby(["drb1_group", "expected"])
              .agg(n=("concordant", "size"), n_concordant=("concordant", "sum"))
              .reset_index())
    summ["pct_concordant"] = 100.0 * summ["n_concordant"] / summ["n"]
    summ["n_disp"] = summ["n"].map(suppress)
    summ["n_concordant_disp"] = summ["n_concordant"].map(suppress)
    return df, summ.drop(columns=["n", "n_concordant"])


def c4_copy_number(t1):
    """C4 copy number and long/short composition per haplotype.

    C4A/C4B copy number is the best-characterised CNV in the MHC and the one where long reads
    should most clearly beat short reads. `c4_size` (L/S, the intron-9 HERV insertion) comes
    straight from Table 1.
    """
    c4 = t1[t1["gene_bare"].isin(["C4A", "C4B"])].copy()
    if c4.empty:
        return pd.DataFrame(), pd.DataFrame()
    per_hap = (c4.groupby(["person_id", "hap", "gene_bare"])
                 .size().unstack(fill_value=0).reset_index())
    for g in ["C4A", "C4B"]:
        if g not in per_hap.columns:
            per_hap[g] = 0
    per_hap["total_c4"] = per_hap["C4A"] + per_hap["C4B"]
    dist = (per_hap.groupby(["C4A", "C4B", "total_c4"]).size()
            .reset_index(name="n_haplotypes"))
    dist["n_haplotypes_disp"] = dist["n_haplotypes"].map(suppress)
    dist["pct"] = 100.0 * dist["n_haplotypes"] / dist["n_haplotypes"].sum()

    size = pd.DataFrame()
    if "c4_size" in c4.columns:
        s = c4[c4["c4_size"].isin(["L", "S"])]
        if not s.empty:
            size = (s.groupby(["gene_bare", "c4_size"]).size().reset_index(name="n"))
            size["n_disp"] = size["n"].map(suppress)
            size["pct_within_gene"] = size.groupby("gene_bare")["n"].transform(
                lambda v: 100.0 * v / v.sum())
            size = size.drop(columns=["n"])
    return dist.drop(columns=["n_haplotypes"]), size


def kir_audit(t1_all):
    """A5. KIR is on chr19, outside the chr6 trim window, so the expected answer is zero calls.
    Reported as a measured fact rather than an assumption -- and loudly if nonzero, because a
    nonzero count means the trim window is wrong and every assembly-derived number is suspect."""
    if "gene_class" in t1_all.columns:
        kir = t1_all[t1_all["gene_class"].astype(str).str.lower() == "kir"]
    else:
        kir = t1_all[t1_all["gene_bare"].astype(str).str.upper().str.startswith("KIR")]
    genes = sorted(set(kir["gene_bare"])) if not kir.empty else []
    return {"n_kir_rows": int(len(kir)), "kir_genes_seen": genes,
            "verdict": ("EXPECTED: no KIR calls; the KIR cluster is on chr19, outside the chr6 "
                        "trim window used to build these assemblies. Genotyping KIR would need a "
                        "separate extraction from the original BAMs."
                        if len(kir) == 0 else
                        "UNEXPECTED: KIR calls are present. This means the trim window is not "
                        "what SCHEMA.md documents. Stop and re-check the extraction before "
                        "trusting any assembly-derived result.")}


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def fig_deletion_rates(gene_rates, path, gene_order):
    """Bridged-deletion rate per gene, ordered by physical position, controls colour-coded."""
    plt = _mpl()
    df = gene_rates.copy()
    if df.empty:
        return
    df["rank"] = df["gene"].map(lambda g: gene_order.get(g, 1e9))
    df = df.sort_values("rank")
    colors = []
    for g in df["gene"]:
        if g in POSITIVE_CONTROL_GENES:
            colors.append("#009E73")       # known CNV -- expect a high rate
        elif g in NEGATIVE_CONTROL_GENES:
            colors.append("#999999")       # expect ~0; this is the false-positive rate
        else:
            colors.append("#D55E00")
    fig, ax = plt.subplots(figsize=(max(7.0, 0.32 * len(df)), 4.2))
    ax.bar(range(len(df)), df["pct_deleted_bridged"], color=colors, width=0.75, zorder=3)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["gene"], rotation=90, fontsize=7)
    ax.set_ylabel("% of bridged haplotypes lacking the gene")
    ax.set_title("Gene deletions called only where one contig spans the gene's position\n"
                 "green = known copy-number genes (positive control) · grey = genes that should "
                 "never be deleted (false-positive rate) · orange = everything else", fontsize=9)
    ax.grid(axis="y", lw=0.4, alpha=0.4, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_deletion_by_ancestry(anc_rates, path, top_genes):
    plt = _mpl()
    df = anc_rates[anc_rates["gene"].isin(top_genes)]
    if df.empty:
        return
    genes = [g for g in top_genes if g in set(df["gene"])]
    ancs = [a for a in ANCESTRY_ORDER if a in set(df["ancestry"])]
    fig, ax = plt.subplots(figsize=(1.5 * len(genes) + 3.0, 4.0))
    width = 0.8 / max(1, len(ancs))
    for ai, anc in enumerate(ancs):
        xs, ys = [], []
        for gi, g in enumerate(genes):
            r = df[(df["gene"] == g) & (df["ancestry"] == anc)]
            if r.empty:
                continue
            xs.append(gi - 0.4 + width * (ai + 0.5))
            ys.append(r["pct_deleted_bridged"].iloc[0])
        if xs:
            ax.bar(xs, ys, width=width * 0.9, color=ANCESTRY_COLORS.get(anc, "#777"), label=anc)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes)
    ax.set_ylabel("% of bridged haplotypes lacking the gene")
    ax.set_title("Gene-deletion frequency by ancestry (strict ancestry, unrelated people)",
                 fontsize=9)
    ax.legend(frameon=False, ncol=len(ancs), fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_c4(dist, path):
    plt = _mpl()
    if dist.empty:
        return
    d = dist.sort_values("pct", ascending=False).head(10).copy()
    labels = ["A%d/B%d" % (int(r["C4A"]), int(r["C4B"])) for _, r in d.iterrows()]
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    ax.bar(range(len(d)), d["pct"], color="#4C72B0", width=0.7)
    ax.set_xticks(range(len(d)))
    ax.set_xticklabels(labels, rotation=0, fontsize=8)
    ax.set_ylabel("% of haplotypes")
    ax.set_xlabel("C4A / C4B copies on the haplotype")
    ax.set_title("C4 copy number per haplotype (phased, from long-read assemblies)", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def write_readme(path, args, gene_rates, anc_rates, dup, drb_summ, c4_dist, c4_size, kir, diag,
                 n_people):
    L = []
    L.append("# 30 — structural variation across the MHC: deletions, duplications, copy number\n")
    L.append("*Supervisor asks A4 (deletions/duplications) and A5 (KIR), 2026-09-17 call.*\n")

    L.append("## The method, and why absence is not deletion\n")
    L.append("A gene missing from a haplotype can mean a real deletion, a broken assembly, or a "
             "failed homology search. This script calls a deletion **only** when a single contig "
             "carries a gene on each side of the missing gene in canonical physical order — i.e. "
             "the assembly was built through the gene's position and the gene was not there. "
             "Absences without flanking evidence are counted separately and never interpreted.\n")
    L.append("- unrelated people: **%s**; haplotypes examined: **%s**\n"
             % (suppress(n_people), suppress(diag.get("haplotypes_seen", 0))))
    L.append("- bridged absences: **%s**; unbridged (assembly ends): **%s**\n"
             % (suppress(diag.get("bridged_absence", 0)),
                suppress(diag.get("unbridged_absence", 0))))

    L.append("\n## 1. Positive control — DRB3/DRB4/DRB5 against the DR51/52/53 expectation\n")
    L.append("DRB1 haplotype group determines which second DRB locus is present. The answer is "
             "known from immunogenetics, independently of our data, so this measures whether the "
             "deletion caller works at all.\n")
    if not drb_summ.empty:
        L.append(drb_summ.to_markdown(index=False) + "\n")
        L.append("\n**Read this first.** If concordance is high, the deletion calls below are "
                 "trustworthy. If it is not, nothing else in this report is.\n")

    L.append("\n## 2. Negative control — genes that should never be deleted\n")
    neg = gene_rates[gene_rates["gene"].isin(NEGATIVE_CONTROL_GENES)]
    if not neg.empty:
        L.append(neg[["gene", "n_bridged_disp", "pct_deleted_bridged"]].to_markdown(index=False)
                 + "\n")
        L.append("\nThis is the **false-positive rate** of the method, not biology. Every rate in "
                 "the next section should be read against it.\n")

    L.append("\n## 3. Deletion frequency per gene\n")
    if not gene_rates.empty:
        top = gene_rates.sort_values("pct_deleted_bridged", ascending=False).head(25)
        L.append(top[["gene", "gene_class", "n_bridged_disp", "pct_deleted_bridged"]]
                 .to_markdown(index=False) + "\n")

    L.append("\n## 4. Deletion frequency by ancestry\n")
    if not anc_rates.empty:
        piv = anc_rates.pivot_table(index="gene", columns="ancestry",
                                    values="pct_deleted_bridged").round(2)
        L.append(piv.to_markdown() + "\n")

    L.append("\n## 5. Duplication candidates (>1 copy of a gene on one contig)\n")
    if not dup.empty:
        d = dup.head(20).copy()
        d["n_multi_copy_disp"] = d["n_multi_copy"].map(suppress)
        d["n_haplotype_contigs_disp"] = d["n_haplotype_contigs"].map(suppress)
        L.append(d[["gene", "n_haplotype_contigs_disp", "n_multi_copy_disp", "pct_multi_copy",
                    "max_copy_seen"]].to_markdown(index=False) + "\n")
        L.append("\nThese are **candidates**. A second mapping cluster inside a segmental "
                 "duplication can be an alignment artifact; confirming a duplication needs "
                 "read-depth or the assembly graph, neither of which is used here.\n")

    L.append("\n## 6. C4 copy number\n")
    if not c4_dist.empty:
        L.append(c4_dist.sort_values("pct", ascending=False).head(10)
                 [["C4A", "C4B", "total_c4", "n_haplotypes_disp", "pct"]]
                 .to_markdown(index=False) + "\n")
    if not c4_size.empty:
        L.append("\nLong/short (intron-9 HERV) composition:\n\n")
        L.append(c4_size.to_markdown(index=False) + "\n")

    L.append("\n## 7. KIR (ask A5)\n")
    L.append("- KIR rows in the calls table: **%d**\n" % kir["n_kir_rows"])
    L.append("- %s\n" % kir["verdict"])

    L.append("\n## Caveats\n")
    L.append("- Deletion calls are annotation-level, not sequence-level. A confirmed deletion "
             "would show the breakpoint; we show only that the assembled contig skipped the "
             "gene's position.\n")
    L.append("- Canonical gene order is derived empirically from the cohort's own assemblies "
             "(script 16's `derive_canonical_order`), because `gene_start` is contig-relative, "
             "not hg38 (SCHEMA.md Table 1).\n")
    L.append("- Duplication candidates are not validated by depth.\n")
    L.append("- **A gene at either end of the annotated region can never be bridged**, because "
             "nothing flanks it on one side. Such a gene will always report a 0% deletion rate "
             "here — that is 'not testable by this method', not 'never deleted'. Check a gene's "
             "`n_bridged` before reading its rate.\n")
    with open(path, "w") as fh:
        fh.write("\n".join(L))


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def run(args):
    ensure_dir(args.out_dir)
    if not os.path.exists(args.table1):
        sys.exit("FATAL: table1 not found: %s" % args.table1)

    print("[30] loading table1 ...", flush=True)
    t1 = pd.read_csv(args.table1, sep="\t", dtype=str, low_memory=False)
    t1["copy_index"] = pd.to_numeric(t1.get("copy_index", 1), errors="coerce").fillna(1).astype(int)
    t1["gene_bare"] = t1["gene"].map(gene_bare)

    kir = kir_audit(t1)
    print("[30] KIR: %s" % kir["verdict"].split(".")[0], flush=True)

    people, n_removed = mod("24_novelty_by_field.py", "novelty_by_field").build_people(
        sorted(set(t1["person_id"].astype(str))), args.cohort_membership,
        args.relatedness_table, args.kin_min, args.strict_threshold, args.skip_relatedness)
    keep = people[people["unrelated"]]
    anc_of = dict(zip(keep["person_id"].astype(str), keep["anc_strict"]))
    n_people = len(keep)
    t1 = t1[t1["person_id"].astype(str).isin(anc_of)].copy()
    print("[30] unrelated people: %d (dropped %d)" % (n_people, n_removed), flush=True)

    # derive_canonical_order sorts on gene_start, which is a string here because table1 is read
    # with dtype=str. A string sort would order "1000" before "900" and silently produce a wrong
    # physical order -- which would then produce wrong bridging decisions that look plausible.
    t1["gene_start"] = pd.to_numeric(t1.get("gene_start"), errors="coerce")
    gene_order = mod("16_phasing_mendelian_validation.py", "phasing_val").derive_canonical_order(t1)
    if len(gene_order) < 5:
        sys.exit("FATAL: canonical gene order has only %d genes; gene_start is probably missing "
                 "or unparseable. Bridging cannot be decided without it." % len(gene_order))
    genes_of_interest = sorted(set(t1["gene_bare"]) & set(gene_order))
    print("[30] canonical order over %d genes" % len(gene_order), flush=True)

    calls, diag = bridged_absences(t1, gene_order, genes_of_interest)

    gene_class = (t1.drop_duplicates("gene_bare").set_index("gene_bare")["gene_class"].to_dict()
                  if "gene_class" in t1.columns else {})

    def rate_table(df, extra_keys=()):
        keys = list(extra_keys) + ["gene"]
        g = df.groupby(keys)["status"].agg(
            n_present=lambda s: int((s == "present").sum()),
            n_deleted=lambda s: int((s == "deleted_bridged").sum()),
            n_unbridged=lambda s: int((s == "absent_unbridged").sum())).reset_index()
        g["n_bridged"] = g["n_present"] + g["n_deleted"]
        g["pct_deleted_bridged"] = (100.0 * g["n_deleted"] / g["n_bridged"].clip(lower=1)).round(3)
        g["n_bridged_disp"] = g["n_bridged"].map(suppress)
        g["n_deleted_disp"] = g["n_deleted"].map(suppress)
        g["gene_class"] = g["gene"].map(lambda x: gene_class.get(x, "NA"))
        return g.drop(columns=["n_present", "n_deleted", "n_bridged", "n_unbridged"])

    gene_rates = rate_table(calls)
    calls["ancestry"] = calls["person_id"].astype(str).map(anc_of)
    anc_calls = calls[calls["ancestry"].isin(ANCESTRY_ORDER)]
    anc_rates = rate_table(anc_calls, extra_keys=["ancestry"]) if not anc_calls.empty \
        else pd.DataFrame()

    dup = duplications(t1, genes_of_interest)
    drb_detail, drb_summ = drb_expectation_check(t1, calls)
    c4_dist, c4_size = c4_copy_number(t1)

    # ---- write ----
    gene_rates.to_csv(os.path.join(args.out_dir, "deletion_rates_by_gene.tsv"), sep="\t", na_rep="NA",
                      index=False)
    if not anc_rates.empty:
        anc_rates.to_csv(os.path.join(args.out_dir, "deletion_rates_by_gene_ancestry.tsv"),
                         sep="\t", index=False, na_rep="NA")
    dup_out = dup.copy()
    dup_out["n_multi_copy"] = dup_out["n_multi_copy"].map(suppress)
    dup_out["n_haplotype_contigs"] = dup_out["n_haplotype_contigs"].map(suppress)
    dup_out.to_csv(os.path.join(args.out_dir, "duplication_candidates.tsv"), sep="\t", index=False, na_rep="NA")
    if not drb_summ.empty:
        drb_summ.to_csv(os.path.join(args.out_dir, "drb_positive_control.tsv"), sep="\t", na_rep="NA",
                        index=False)
    if not c4_dist.empty:
        c4_dist.to_csv(os.path.join(args.out_dir, "c4_copy_number.tsv"), sep="\t", index=False, na_rep="NA")
    if not c4_size.empty:
        c4_size.to_csv(os.path.join(args.out_dir, "c4_long_short.tsv"), sep="\t", index=False, na_rep="NA")

    fig_deletion_rates(gene_rates, os.path.join(args.out_dir, "fig_deletion_rates.png"),
                       gene_order)
    top_genes = (gene_rates.sort_values("pct_deleted_bridged", ascending=False)
                 .head(8)["gene"].tolist())
    if not anc_rates.empty:
        fig_deletion_by_ancestry(anc_rates,
                                 os.path.join(args.out_dir, "fig_deletion_by_ancestry.png"),
                                 top_genes)
    fig_c4(c4_dist, os.path.join(args.out_dir, "fig_c4_copy_number.png"))

    write_readme(os.path.join(args.out_dir, "README.md"), args, gene_rates, anc_rates, dup,
                 drb_summ, c4_dist, c4_size, kir, diag, n_people)
    with open(os.path.join(args.out_dir, "summary.json"), "w") as fh:
        json.dump({"kir": kir, "n_people_unrelated": suppress(n_people),
                   "haplotypes_seen": suppress(diag.get("haplotypes_seen", 0)),
                   "bridged_absences": suppress(diag.get("bridged_absence", 0)),
                   "unbridged_absences": suppress(diag.get("unbridged_absence", 0))},
                  fh, indent=2)
    print("[30] done -> %s" % args.out_dir, flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=DEFAULT_TABLE1)
    ap.add_argument("--cohort-membership", default=DEFAULT_COHORT)
    ap.add_argument("--relatedness-table", default=DEFAULT_RELATEDNESS)
    ap.add_argument("--skip-relatedness", action="store_true")
    ap.add_argument("--kin-min", type=float, default=KIN_MIN)
    ap.add_argument("--strict-threshold", type=float, default=STRICT_MIN)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    args = ap.parse_args(argv)
    run(args)


if __name__ == "__main__":
    main()
