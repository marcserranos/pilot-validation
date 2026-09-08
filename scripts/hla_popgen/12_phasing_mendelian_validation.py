#!/usr/bin/env python3
"""Validate long-read (Immuannot) haplotype phasing using REAL relative pairs already in the
biobank -- the "easy" case identified by 11_relatedness_cohort_overlap.py: 827 long-read-cohort
people have >=1 relative who is ALSO in the long-read cohort (569 first-degree pairs + 5 likely
duplicate/MZ-twin pairs), so validation is a direct comparison between two already-phased
assemblies -- no alignment pipeline needed.

## The core idea (Mendelian consistency as a phasing-accuracy probe)

A parent transmits exactly ONE of their two haplotype copies to a child at every locus -- so a true
parent-child pair must share at least one allele at every gene along the region (IBD=1, always).
Full siblings are noisier: at any given locus their sharing state (IBD0/1/2) depends on which
grandparental copies were transmitted, BUT the MHC region is well known to have very low internal
recombination (it behaves almost like a single haplotype block across a single meiosis), so a
sibling pair's IBD state at HLA is close to constant across this whole ~4Mb window, not scattered
gene-by-gene -- exactly the same "0.25/0.25/0.25 kinship, mostly one IBD state throughout the
region" fact that transplant medicine relies on for sibling HLA matching. This script cannot tell a
parent-child pair from a full-sibling pair (11_relatedness_cohort_overlap.py's documented
limitation: the AoU relatedness table carries no IBD0 column) -- so instead of assuming, it
EMPIRICALLY classifies each pair by its own genome-wide-in-this-region sharing rate:

  - high_sharing  (>= --high-thresh, default 0.8): consistent with a real parent-child pair, OR an
    IBD1/IBD2 sibling pair -- both give one unbroken shared-haplotype block. INFORMATIVE for the
    error-rate scan below.
  - low_sharing   (<= --low-thresh, default 0.3): consistent with an IBD0 sibling pair -- the two
    relatives simply didn't inherit the same HLA haplotype from either parent. NOT an error, just
    biology -- excluded from the error-rate scan (but reported).
  - mixed         (in between): the interesting middle case -- could be a genuine recombination
    breakpoint landing inside this specific 4Mb window in the meiosis that produced one of these two
    people (rare, but real, and independently interesting), or a switch error muddying the whole
    picture. Reported separately, inspected via the per-gene ideogram figure, never pooled into the
    headline error rate (its own sharing pattern isn't a clean baseline to test individual genes
    against).

**Within high_sharing pairs**, a specific gene where NEITHER of the two people's alleles is shared
by the other, despite the pair's overwhelming genome-wide-in-region consistency, is a genuine
Mendelian-incompatible observation -- almost certainly a genotyping/assembly error in ONE of the
two people at that specific locus (or, far less likely at HLA-typical allele diversity, two
independent people who happen to share a haplotype but not that one rare allele). Aggregating this
per-gene mismatch rate across all high_sharing pairs is the direct empirical error-rate estimate
Marc asked for -- and comparing the classical HLA genes' rate against non-classical/pseudogene
genes' rate answers "are errors concentrated in the HLA genes specifically."

## The secondary idea (switch/crossover detection)

Within high_sharing pairs, walking genes in their real physical order and tracking WHICH of person
A's two haplotypes is linked to WHICH of person B's two haplotypes at each gene turns a Mendelian
scan into a phase-consistency scan: zero linkage switches across the whole region is the expected,
literature-supported result (HLA as one recombination-cold block); >=1 switch is either a genuine
meiotic crossover in this family or an assembly switch-error in one of the two people's phased
contigs -- exactly the failure mode 11's report flagged as a live risk.

Physical gene order is not read off `gene_start` directly (SCHEMA.md: contig-relative, not hg38,
and can differ hap-to-hap when an assembly is fragmented across contigs) -- instead this script
DERIVES a canonical order empirically from the real cohort's own assemblies (median relative rank
of each gene's `gene_start` within every single-contig hap instance, aggregated across the whole
long-read cohort), the same "let real data settle it" principle already used throughout this
project (SCHEMA.md, 00_recon_vm.py) rather than hardcoding remembered hg38 coordinates.

Usage (fixtures):
    python3 scripts/hla_popgen/12_phasing_mendelian_validation.py \\
        --table1 /tmp/hla_fixtures/hla_calls_rich.sample.tsv \\
        --pair-list /tmp/hla_fixtures/relatedness_lr_overlap_pairs.sample.tsv \\
        --out-dir /tmp/hla_fixtures/reports/12_phasing

Real run (VM):
    python3 scripts/hla_popgen/12_phasing_mendelian_validation.py
"""
import argparse
import os
import sys
from collections import defaultdict, Counter

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _viz_common as vc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

DEFAULT_PAIR_LIST = os.path.expanduser("~/pipeline_outputs/relatedness_lr_overlap_pairs.tsv")
INFORMATIVE_DEGREES = {"first_degree_parentchild_or_sibling", "duplicate_or_MZ_twin"}
GENE_CLASS_ORDER = ["classical_I", "classical_II", "class_II_accessory", "class_II_paralog",
                    "nonclassical_I", "pseudogene_I", "mic_tap", "complement", "other"]
GENE_CLASS_COLORS = {
    "classical_I": "#0072B2", "classical_II": "#D55E00", "class_II_accessory": "#CC79A7",
    "class_II_paralog": "#E69F00", "nonclassical_I": "#009E73", "pseudogene_I": "#999999",
    "mic_tap": "#56B4E9", "complement": "#F0E442", "other": "#000000",
}


# ---------------------------------------------------------------------------
# Data assembly
# ---------------------------------------------------------------------------
def build_allele_lookup(table1, resolution):
    """{(person_id, gene_bare, hap): allele_nfield}, plus a parallel {(person_id, gene_bare, hap):
    contig} lookup. copy_index>1 (segmental duplication / DRB paralog copy number) rows are dropped
    -- comparing which specific extra copy matches across two unrelated assemblies is a different,
    harder question than this script asks.

    The contig lookup exists because of a real bug this script's first real run caught: SCHEMA.md
    is explicit that "the only valid cis key is (person_id, hap, contig)" -- a hap1/hap2 LABEL is
    only guaranteed phase-consistent WITHIN one contig, not across contigs, when a person's own
    assembly is fragmented across the ~4Mb HLA region (routine here, per SCHEMA.md Table 2's own
    pairable-fraction discussion). The switch-scan below must never compare linkage state across a
    contig boundary for either person -- doing so measures arbitrary hap-label reshuffling between
    independently-assembled fragments, not a real crossover or a real assembly switch-error."""
    sub = table1[table1["copy_index"].fillna(1) == 1].copy()
    sub["allele_n"] = sub["consensus"].map(lambda a: vc.to_nfield(a, resolution))
    sub = sub.dropna(subset=["allele_n"])
    sub = sub[sub["hap"].isin(["hap1", "hap2"])]
    lookup = {}
    contig_lookup = {}
    for row in sub.itertuples(index=False):
        lookup[(row.person_id, row.gene_bare, row.hap)] = row.allele_n
        contig_lookup[(row.person_id, row.gene_bare, row.hap)] = row.contig
    gene_class = dict(zip(sub["gene_bare"], sub["gene_class"]))
    return lookup, gene_class, contig_lookup


def derive_canonical_order(table1):
    """Median fractional rank of each gene's gene_start within every single-contig (person, hap,
    contig) group that calls >=2 genes -- see module docstring. Returns {gene_bare: rank_float},
    lower = earlier in the region. Robust to fragmented assemblies: a hap contributes whatever
    partial ordering info it has, and thousands of hap-instances aggregate into a stable order."""
    sub = table1[table1["copy_index"].fillna(1) == 1].dropna(subset=["gene_start"])
    ranks = defaultdict(list)
    for (person, hap, contig), grp in sub.groupby(["person_id", "hap", "contig"]):
        grp = grp.drop_duplicates("gene_bare").sort_values("gene_start")
        n = len(grp)
        if n < 2:
            continue
        for i, gene in enumerate(grp["gene_bare"]):
            ranks[gene].append(i / (n - 1))
    return {gene: float(np.median(vals)) for gene, vals in ranks.items() if vals}


# ---------------------------------------------------------------------------
# Pair-level analysis
# ---------------------------------------------------------------------------
def analyze_pair(person_a, person_b, genes, allele_lookup, gene_order, contig_lookup=None):
    """Returns a dict of per-pair summary stats plus a per-gene detail list, walked in canonical
    gene order. `genes`: iterable of gene_bare names to consider (both-called anywhere in cohort).
    `contig_lookup`, if given, records which contig backs the specific hap chosen for each resolved
    linkage state, so count_switches can refuse to compare state across a contig boundary (see
    build_allele_lookup's docstring for why that comparison would be meaningless)."""
    contig_lookup = contig_lookup or {}
    ordered_genes = sorted(
        (g for g in genes if g in gene_order or True),  # keep genes even if order unknown (tail)
        key=lambda g: gene_order.get(g, float("inf")))

    detail = []
    n_compared, n_shared = 0, 0
    for gene in ordered_genes:
        a1 = allele_lookup.get((person_a, gene, "hap1"))
        a2 = allele_lookup.get((person_a, gene, "hap2"))
        b1 = allele_lookup.get((person_b, gene, "hap1"))
        b2 = allele_lookup.get((person_b, gene, "hap2"))
        a_set = {x for x in (a1, a2) if x is not None}
        b_set = {x for x in (b1, b2) if x is not None}
        if not a_set or not b_set:
            continue
        n_compared += 1
        shared_alleles = a_set & b_set
        shared = bool(shared_alleles)
        n_shared += int(shared)

        # Linkage state: which of A's hap indices carries a shared allele, and which of B's.
        # Left None (uninformative) when either side is homozygous at the shared allele (can't
        # tell which copy) or no sharing at all.
        state = None
        a_contig = b_contig = None
        if shared and len(shared_alleles) == 1:
            allele = next(iter(shared_alleles))
            a_idx = [i for i, a in ((1, a1), (2, a2)) if a == allele]
            b_idx = [i for i, b in ((1, b1), (2, b2)) if b == allele]
            if len(a_idx) == 1 and len(b_idx) == 1:
                state = (a_idx[0], b_idx[0])
                a_contig = contig_lookup.get((person_a, gene, f"hap{a_idx[0]}"))
                b_contig = contig_lookup.get((person_b, gene, f"hap{b_idx[0]}"))
        detail.append({"gene": gene, "shared": shared, "state": state,
                        "a_contig": a_contig, "b_contig": b_contig,
                        "a1": a1, "a2": a2, "b1": b1, "b2": b2})

    sharing_frac = n_shared / n_compared if n_compared else np.nan
    return {"person_a": person_a, "person_b": person_b, "n_compared": n_compared,
            "n_shared": n_shared, "sharing_frac": sharing_frac, "detail": detail}


def classify_pair(sharing_frac, n_compared, min_genes, high_thresh, low_thresh):
    if n_compared < min_genes or pd.isna(sharing_frac):
        return "insufficient_data"
    if sharing_frac >= high_thresh:
        return "high_sharing"
    if sharing_frac <= low_thresh:
        return "low_sharing"
    return "mixed"


def count_switches(detail):
    """Walk the gene-ordered detail list (already in canonical order) and count linkage-state
    changes among genes with a resolved state, skipping unresolved/mismatched genes rather than
    treating them as a state (a Mendelian-incompatible or homozygous-uninformative gene carries no
    phase information, it is not evidence of a switch by itself).

    Critically, a transition is only counted when BOTH people's resolved hap at the two adjacent
    genes sits on the SAME contig as it did at the previous resolved gene. A hap1/hap2 label is only
    guaranteed phase-consistent within one contig (SCHEMA.md); comparing state across a contig
    boundary for either person measures arbitrary hap-label reshuffling between independently
    assembled fragments, not a real crossover or assembly switch-error. A first real run of this
    script (2026-09-08) without this guard reported 518/545 high-sharing pairs with >=1 "switch" --
    implausibly high for real meiotic recombination in a ~4Mb window -- which is exactly the
    contig-reshuffling artifact this guard exists to exclude."""
    resolved = [d for d in detail if d["state"] is not None]
    if len(resolved) < 2:
        return 0, len(resolved)
    switches = 0
    for i in range(1, len(resolved)):
        prev, cur = resolved[i - 1], resolved[i]
        same_contig_run = (cur["a_contig"] == prev["a_contig"] and
                            cur["b_contig"] == prev["b_contig"])
        if same_contig_run and cur["state"] != prev["state"]:
            switches += 1
    return switches, len(resolved)


# ---------------------------------------------------------------------------
# Aggregation across all pairs
# ---------------------------------------------------------------------------
def run_all_pairs(pairs_df, allele_lookup, gene_order, all_genes, min_genes, high_thresh, low_thresh,
                   contig_lookup=None):
    pair_rows = []
    per_gene_mismatch = defaultdict(lambda: {"n_compared": 0, "n_mismatch": 0})
    switch_counts = []
    mixed_examples = []

    for row in pairs_df.itertuples(index=False):
        result = analyze_pair(row.person_i, row.person_j, all_genes, allele_lookup, gene_order,
                               contig_lookup)
        cls = classify_pair(result["sharing_frac"], result["n_compared"], min_genes,
                             high_thresh, low_thresh)
        pair_rows.append({
            "person_i": row.person_i, "person_j": row.person_j, "kin": row.kin,
            "kin_degree": row.kin_degree, "n_compared": result["n_compared"],
            "sharing_frac": result["sharing_frac"], "pair_class": cls,
        })
        if cls == "high_sharing":
            for d in result["detail"]:
                pg = per_gene_mismatch[d["gene"]]
                pg["n_compared"] += 1
                pg["n_mismatch"] += int(not d["shared"])
            n_sw, n_states = count_switches(result["detail"])
            switch_counts.append({"person_i": row.person_i, "person_j": row.person_j,
                                   "n_switches": n_sw, "n_states": n_states})
        elif cls == "mixed" and len(mixed_examples) < 6:
            mixed_examples.append(result)

    pairs_out = pd.DataFrame(pair_rows)
    gene_rows = [{"gene": g, "n_compared": v["n_compared"], "n_mismatch": v["n_mismatch"]}
                 for g, v in per_gene_mismatch.items()]
    gene_out = pd.DataFrame(gene_rows)
    switches_out = pd.DataFrame(switch_counts)
    return pairs_out, gene_out, switches_out, mixed_examples


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_per_gene_error_rate(gene_out, gene_class, gene_order, min_n, out_path):
    df = gene_out[gene_out["n_compared"] >= min_n].copy()
    if df.empty:
        return None
    df["gene_class"] = df["gene"].map(gene_class).fillna("other")
    df["rank"] = df["gene"].map(lambda g: gene_order.get(g, float("inf")))
    df = df.sort_values("rank")
    p, lo, hi = vc.wilson_ci(df["n_mismatch"].values, df["n_compared"].values)
    df["p"], df["lo"], df["hi"] = p, lo, hi

    fig, ax = plt.subplots(figsize=(max(10, 0.35 * len(df)), 5))
    colors = [GENE_CLASS_COLORS.get(c, "#000000") for c in df["gene_class"]]
    x = np.arange(len(df))
    ax.bar(x, df["p"] * 100, color=colors, width=0.7, zorder=3)
    # np.clip guards a floating-point sliver at the Wilson-interval boundary (observed: p
    # infinitesimally < lo for a couple of n_mismatch==0 genes) that otherwise trips
    # matplotlib's "yerr must not contain negative values" check.
    yerr_lo = np.clip((df["p"] - df["lo"]) * 100, 0, None)
    yerr_hi = np.clip((df["hi"] - df["p"]) * 100, 0, None)
    ax.errorbar(x, df["p"] * 100, yerr=[yerr_lo, yerr_hi],
                fmt="none", ecolor="black", elinewidth=0.8, capsize=2, zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels([vc.gene_display(g) for g in df["gene"]], rotation=90, fontsize=7)
    ax.set_ylabel("Mendelian-incompatible rate (%)\namong high-sharing relative pairs")
    ax.set_title("Per-gene phasing/calling error rate, ordered by physical position in the region\n"
                  "(derived canonical order; color = gene class; error bars = 95% Wilson CI)")
    for c in sorted(set(df["gene_class"]), key=lambda c: GENE_CLASS_ORDER.index(c)
                     if c in GENE_CLASS_ORDER else 99):
        ax.bar([], [], color=GENE_CLASS_COLORS.get(c, "#000000"), label=c)
    ax.legend(fontsize=7, ncol=3, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(bottom=0)
    vc.savefig(fig, out_path)
    return df


def plot_sharing_distribution(pairs_out, high_thresh, low_thresh, out_path):
    df = pairs_out.dropna(subset=["sharing_frac"])
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(df["sharing_frac"], bins=np.linspace(0, 1, 41), color="#0072B2",
            edgecolor="white", linewidth=0.3)
    ax.axvline(high_thresh, color="#D55E00", linestyle="--", linewidth=1,
               label=f"high_sharing >= {high_thresh}")
    ax.axvline(low_thresh, color="#009E73", linestyle="--", linewidth=1,
               label=f"low_sharing <= {low_thresh}")
    ax.set_xlabel("Per-gene allele-sharing fraction, this relative pair")
    ax.set_ylabel("Number of pairs")
    ax.set_title("Pair-level IBD-sharing profile across the HLA region\n"
                 "(expected trimodal: IBD0 / IBD1-or-parent-child / IBD2, per HLA-matching biology)")
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    vc.savefig(fig, out_path)


def plot_switch_distribution(switches_out, out_path):
    fig, ax = plt.subplots(figsize=(5, 4))
    counts = switches_out["n_switches"].clip(upper=3)
    labels = ["0\n(clean)", "1", "2", "3+"]
    vals = [int((counts == i).sum()) for i in range(3)] + [int((counts >= 3).sum())]
    ax.bar(labels, vals, color="#0072B2")
    for i, v in enumerate(vals):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Number of high-sharing pairs")
    ax.set_xlabel("Linkage-phase switches across the region")
    ax.set_title("Switch-point count per pair\n"
                 "(0 = single unbroken IBD block, the expected/literature result;\n"
                 ">=1 = candidate crossover OR assembly switch-error)")
    ax.spines[["top", "right"]].set_visible(False)
    vc.savefig(fig, out_path)


def plot_mixed_example(result, gene_order, out_path, example_label):
    """Illustrative per-gene sharing ideogram for one 'mixed' pair -- shows whether the region
    splits into a clean high-sharing block followed by a clean low-sharing block (consistent with
    a real recombination breakpoint within the window) or looks like scattered noise instead.
    `example_label` (e.g. "Example pair A") stands in for the real person_ids -- this figure is
    aggregate-only output under reports/, so it must never display a bare person_id (SCHEMA.md hard
    rule 5)."""
    detail = sorted(result["detail"], key=lambda d: gene_order.get(d["gene"], float("inf")))
    if not detail:
        return
    x = np.arange(len(detail))
    shared = [1 if d["shared"] else 0 for d in detail]
    fig, ax = plt.subplots(figsize=(max(8, 0.3 * len(detail)), 2.2))
    ax.bar(x, shared, color=["#0072B2" if s else "#D55E00" for s in shared], width=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([vc.gene_display(d["gene"]) for d in detail], rotation=90, fontsize=6)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["no share", "share"])
    ax.set_title(f"Mixed-sharing pair {example_label} (sharing={result['sharing_frac']:.2f}) "
                 f"-- candidate recombination breakpoint")
    vc.savefig(fig, out_path)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def write_report(pairs_out, gene_out_annotated, switches_out, gene_class, out_dir, resolution):
    lines = ["# Phasing validation via real relative pairs -- Mendelian consistency + switch scan\n"]
    lines.append(f"Resolution: {resolution}-field. See module docstring for full method.\n")

    lines.append("## Pair classification\n")
    class_counts = pairs_out["pair_class"].value_counts()
    lines.append("| class | n_pairs |")
    lines.append("|---|---|")
    for cls, n in class_counts.items():
        lines.append(f"| {cls} | {n} |")

    high = gene_out_annotated
    if high is not None and not high.empty:
        total_compared = int(high["n_compared"].sum())
        total_mismatch = int(high["n_mismatch"].sum())
        overall_p, overall_lo, overall_hi = vc.wilson_ci(
            np.array([total_mismatch]), np.array([total_compared]))
        lines.append("\n## Headline error rate (high-sharing pairs only)\n")
        lines.append(f"Overall, across all genes: **{overall_p[0]*100:.3f}%** "
                     f"(95% CI {overall_lo[0]*100:.3f}-{overall_hi[0]*100:.3f}%), "
                     f"{total_mismatch}/{total_compared} gene-comparisons.\n")

        lines.append("### By gene class\n")
        lines.append("| gene_class | n_compared | n_mismatch | error_rate_% | 95% CI |")
        lines.append("|---|---|---|---|---|")
        by_class = high.groupby("gene_class")[["n_compared", "n_mismatch"]].sum().reset_index()
        for row in by_class.itertuples(index=False):
            p, lo, hi = vc.wilson_ci(np.array([row.n_mismatch]), np.array([row.n_compared]))
            lines.append(f"| {row.gene_class} | {row.n_compared} | {row.n_mismatch} | "
                         f"{p[0]*100:.3f} | {lo[0]*100:.3f}-{hi[0]*100:.3f} |")

        classical = high[high["gene_class"].isin(["classical_I", "classical_II"])]
        if not classical.empty:
            c_p, c_lo, c_hi = vc.wilson_ci(
                np.array([classical["n_mismatch"].sum()]), np.array([classical["n_compared"].sum()]))
            lines.append(f"\n**Classical HLA genes (A/B/C/DPA1/DPB1/DQA1/DQB1/DRB1) specifically: "
                         f"{c_p[0]*100:.3f}% (95% CI {c_lo[0]*100:.3f}-{c_hi[0]*100:.3f}%)**\n")

        lines.append("\n### Worst individual genes (by point estimate, n>=10 comparisons)\n")
        worst = high[high["n_compared"] >= 10].copy()
        worst["p"] = worst["n_mismatch"] / worst["n_compared"]
        worst = worst.sort_values("p", ascending=False).head(10)
        lines.append("| gene | gene_class | n_compared | n_mismatch | error_rate_% |")
        lines.append("|---|---|---|---|---|")
        for row in worst.itertuples(index=False):
            lines.append(f"| {vc.gene_display(row.gene)} | {row.gene_class} | {row.n_compared} | "
                         f"{row.n_mismatch} | {row.p*100:.3f} |")

    if switches_out is not None and not switches_out.empty:
        lines.append("\n## Phase-switch scan\n")
        n_pairs = len(switches_out)
        n_clean = int((switches_out["n_switches"] == 0).sum())
        lines.append(f"{n_clean}/{n_pairs} high-sharing pairs ({100*n_clean/n_pairs:.1f}%) show "
                     f"ZERO linkage-phase switches across the region -- one unbroken shared "
                     f"haplotype block, the expected result given the MHC's known low internal "
                     f"recombination.\n")
        with_switch = switches_out[switches_out["n_switches"] > 0]
        if not with_switch.empty:
            lines.append(f"{len(with_switch)} pairs show >=1 switch -- candidate real crossovers "
                         f"or assembly switch-errors (cannot be distinguished without a third "
                         f"family member); see mixed-pair example figures for the closest "
                         f"available illustration of this pattern.\n")

    lines.append("\n## Caveats\n")
    lines.append("- Cannot distinguish parent-child from full-sibling pairs (no IBD0 data) -- "
                 "pair classification is by empirical sharing rate, not known pedigree.\n")
    lines.append("- copy_index>1 (segmental duplications, DRB paralog CNV) genes excluded.\n")
    lines.append("- Genes below the minimum-comparison-count floor are excluded from the "
                 "gene-level table/figure but included in the overall headline rate.\n")
    lines.append("- The switch-scan only counts a transition between two genes on the SAME contig "
                 "for BOTH people (see count_switches docstring) -- comparing across a contig "
                 "boundary would measure arbitrary hap-label reshuffling between independently "
                 "assembled fragments, not a real crossover or assembly error. A pair's switch count "
                 "can therefore undercount true switches that happen to coincide with a contig "
                 "break, but will not report a spurious one.\n")
    vc.write_report(os.path.join(out_dir, "phasing_validation_report.md"), lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=os.path.join(vc.DEFAULT_OUTROOT, "hla_calls_rich.tsv"))
    ap.add_argument("--pair-list", default=DEFAULT_PAIR_LIST)
    ap.add_argument("--out-dir", default=os.path.join(vc.DEFAULT_REPORT_ROOT,
                                                       "12_phasing_mendelian_validation"),
                    help="Aggregate-only outputs (figures, headline report, per-gene counts -- no "
                         "person_ids). SCHEMA.md hard rule 5.")
    ap.add_argument("--outroot", default=vc.DEFAULT_OUTROOT,
                    help="Where person-id-level outputs (pair_classification.tsv, "
                         "switch_counts.tsv) are written -- pipeline_outputs, never reports/.")
    ap.add_argument("--resolution", type=int, default=2, choices=[2, 3, 4])
    ap.add_argument("--high-thresh", type=float, default=0.8)
    ap.add_argument("--low-thresh", type=float, default=0.3)
    ap.add_argument("--min-genes-per-pair", type=int, default=5)
    ap.add_argument("--min-gene-n", type=int, default=10,
                    help="Minimum high-sharing gene-comparisons for a gene to appear in the "
                         "per-gene figure/table (thin-N suppression).")
    args = ap.parse_args()

    print("Loading Table 1...", file=sys.stderr)
    table1 = vc.load_table1(args.table1)
    print(f"  {len(table1)} rows, {table1['person_id'].nunique()} people.", file=sys.stderr)

    print("Building allele lookup...", file=sys.stderr)
    allele_lookup, gene_class, contig_lookup = build_allele_lookup(table1, args.resolution)

    print("Deriving canonical gene order from real assemblies...", file=sys.stderr)
    gene_order = derive_canonical_order(table1)
    print(f"  ordered {len(gene_order)} genes.", file=sys.stderr)

    print("Loading relative pair list...", file=sys.stderr)
    pairs_df = pd.read_csv(args.pair_list, sep="\t", dtype={"person_i": str, "person_j": str})
    pairs_df = pairs_df[pairs_df["bucket"] == "both_lr"]
    pairs_df = pairs_df[pairs_df["kin_degree"].isin(INFORMATIVE_DEGREES)]
    print(f"  {len(pairs_df)} first-degree/twin both_lr pairs.", file=sys.stderr)

    all_genes = sorted(set(g for (_, g, _) in allele_lookup.keys()))

    print("Analyzing pairs...", file=sys.stderr)
    pairs_out, gene_out, switches_out, mixed_examples = run_all_pairs(
        pairs_df, allele_lookup, gene_order, all_genes,
        args.min_genes_per_pair, args.high_thresh, args.low_thresh, contig_lookup)
    print(pairs_out["pair_class"].value_counts().to_string(), file=sys.stderr)

    os.makedirs(args.out_dir, exist_ok=True)
    gene_annotated = plot_per_gene_error_rate(
        gene_out, gene_class, gene_order, args.min_gene_n,
        os.path.join(args.out_dir, "per_gene_error_rate.png"))
    plot_sharing_distribution(pairs_out, args.high_thresh, args.low_thresh,
                               os.path.join(args.out_dir, "pair_sharing_distribution.png"))
    if not switches_out.empty:
        plot_switch_distribution(switches_out, os.path.join(args.out_dir, "switch_distribution.png"))
    for i, ex in enumerate(mixed_examples):
        label = f"Example {chr(ord('A') + i)}"
        plot_mixed_example(ex, gene_order,
                            os.path.join(args.out_dir, f"mixed_example_{i}.png"), label)

    write_report(pairs_out, gene_annotated, switches_out, gene_class, args.out_dir, args.resolution)

    # Aggregate-only (no person_ids) -- safe under reports/.
    gene_out.to_csv(os.path.join(args.out_dir, "per_gene_mismatch_counts.tsv"), sep="\t", index=False)

    # Person-id-level (pair_classification/switch_counts carry real person_i/person_j) -- stays in
    # pipeline_outputs, never reports/ (SCHEMA.md hard rule 5: "no bare person_ids ... under reports/").
    pairs_out.to_csv(os.path.join(args.outroot, "phasing_pair_classification.tsv"), sep="\t", index=False)
    if not switches_out.empty:
        switches_out.to_csv(os.path.join(args.outroot, "phasing_switch_counts.tsv"), sep="\t", index=False)
    print(f"Done. Wrote aggregate outputs to {args.out_dir}, person-level outputs to {args.outroot}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
