#!/usr/bin/env python3
"""Allele-frequency spectrum by ancestry, full production cohort -- directly serves this
project's Aim 1 (context/TASK_CONTEXT.md: HLA allele frequency across diverse ancestries), never
built before at full cohort scale (only smaller pilots, e.g. Experiment A's ancestry-clustering
work on the AoU-native callset -- this is Immuannot's own calls instead).

For each of the 8 classical genes, computes standard population-genetics allele frequency (2-field
resolution) within each ancestry group: count of copies of an allele / (2 x N people with a
resolved call at that gene in that group). Picks the top-N alleles by POOLED (whole-cohort)
frequency per gene so the same allele set is compared consistently across ancestry groups, then
plots a grouped bar per gene (8-panel facet grid), one color per ancestry.

Reads: ~/pipeline_outputs/immuannot_calls.tsv, ~/pipeline_outputs/immuannot_cohort_full.tsv
Writes: markdown + 1 PNG (8-panel facet) under --out-dir. Aggregate-only (frequencies/counts, no
per-person data) -- keep on the VM per the standing egress discipline.

Usage (via `pixi run -e spechla`):
  python3 scripts/production_analysis/analyze_allele_frequency_by_ancestry.py
"""
import argparse
import os
import sys
from collections import Counter

import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
ANCESTRY_COLOR = {
    "AFR": "#4C72B0", "AMR": "#DD8452", "EAS": "#55A868",
    "EUR": "#C44E52", "MID": "#8172B2", "SAS": "#937860",
}
NULL = {"", "NA", "nan", "None", ".", "-"}
DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")
TOP_N = 6


def to_2field(allele):
    if allele is None or (hasattr(pd, "isna") and pd.isna(allele)):
        return None
    s = str(allele).strip()
    if s in NULL:
        return None
    if "*" in s:
        s = s.split("*", 1)[1]
    fields = [f for f in s.split(":") if f != ""]
    fields = [f for f in fields if f.strip().lower() != "new"]
    if len(fields) < 2:
        return None
    return f"{fields[0]}:{fields[1]}"


def load_calls(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for c in ["person_id", "gene", "immuannot_1", "immuannot_2"]:
        if c not in df.columns:
            sys.exit(f"FATAL: {path} missing column '{c}'. Actual: {list(df.columns)}")
    df["gene_bare"] = df["gene"].str.replace("^HLA-", "", regex=True)
    filtered = df[df["gene_bare"].isin(GENES)]
    if filtered.empty:
        seen = sorted(df["gene_bare"].unique())[:20]
        sys.exit(f"FATAL: 0 of {len(df)} rows in {path} match any of the 8 classical genes "
                 f"{GENES}. Genes actually present (first 20): {seen}. This is the exact "
                 f"symptom of the 2026-08-10 merge_fragments() dedup bug (see "
                 f"scripts/production_orchestrator/run_production_orchestrator.py's merge_fragments "
                 f"comment) -- run scripts/production_orchestrator/rebuild_immuannot_calls.py first.")
    return filtered


def load_ancestry(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    for c in ["person_id", "ancestry_pred"]:
        if c not in df.columns:
            sys.exit(f"FATAL: {path} missing column '{c}'. Actual: {list(df.columns)}")
    return dict(zip(df["person_id"], df["ancestry_pred"]))


def compute_frequencies(calls, ancestry):
    """{gene: {ancestry: Counter(allele -> copy count)}} and {gene: {ancestry: n_alleles_total}}."""
    counts = {g: {a: Counter() for a in ANCESTRY_ORDER} for g in GENES}
    totals = {g: {a: 0 for a in ANCESTRY_ORDER} for g in GENES}
    for _, r in calls.iterrows():
        anc = ancestry.get(r["person_id"])
        if anc not in ANCESTRY_ORDER:
            continue
        gene = r["gene_bare"]
        for raw in (r["immuannot_1"], r["immuannot_2"]):
            allele = to_2field(raw)
            if allele is None:
                continue
            counts[gene][anc][f"{gene}*{allele}"] += 1
            totals[gene][anc] += 1
    return counts, totals


def plot_spectrum(counts, totals, out_path):
    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    for ax, gene in zip(axes.flat, GENES):
        pooled = Counter()
        for a in ANCESTRY_ORDER:
            pooled.update(counts[gene][a])
        top_alleles = [a for a, _ in pooled.most_common(TOP_N)]
        if not top_alleles:
            ax.set_title(f"{gene} (no data)")
            ax.axis("off")
            continue

        width = 0.13
        x = range(len(top_alleles))
        for i, anc in enumerate(ANCESTRY_ORDER):
            total = totals[gene][anc]
            freqs = [100 * counts[gene][anc].get(a, 0) / total if total else 0 for a in top_alleles]
            offs = [xi + (i - 2.5) * width for xi in x]
            ax.bar(offs, freqs, width=width, color=ANCESTRY_COLOR[anc],
                   label=anc if gene == GENES[0] else None)
        ax.set_xticks(list(x))
        ax.set_xticklabels([a.split("*", 1)[1] for a in top_alleles], rotation=45, ha="right", fontsize=8)
        ax.set_title(gene, fontsize=11)
        ax.set_ylabel("Allele frequency (%)", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=6, bbox_to_anchor=(0.5, 1.02), frameon=False)
    fig.suptitle(f"Top {TOP_N} alleles per gene by pooled frequency, split by ancestry "
                 "(Field 2 / 2-field resolution)", y=1.06, fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", default=os.path.join(DEFAULT_OUTROOT, "immuannot_calls.tsv"))
    ap.add_argument("--cohort", default=os.path.join(DEFAULT_OUTROOT, "immuannot_cohort_full.tsv"))
    ap.add_argument("--out-dir", default=os.path.join(DEFAULT_OUTROOT, "production_analysis", "allele_frequency"))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    calls = load_calls(args.calls)
    ancestry = load_ancestry(args.cohort)
    counts, totals = compute_frequencies(calls, ancestry)

    fig_path = os.path.join(args.out_dir, "allele_frequency_spectrum.png")
    plot_spectrum(counts, totals, fig_path)

    md = ["# Allele frequency spectrum by ancestry -- full production cohort\n",
          f"Top {TOP_N} alleles per gene (by pooled frequency), Field 2 resolution.\n"]
    for gene in GENES:
        pooled = Counter()
        for a in ANCESTRY_ORDER:
            pooled.update(counts[gene][a])
        md.append(f"\n## {gene}\n")
        md.append("| Allele | " + " | ".join(ANCESTRY_ORDER) + " | Pooled n copies |")
        md.append("|---|" + "---|" * (len(ANCESTRY_ORDER) + 1))
        for allele, pooled_n in pooled.most_common(TOP_N):
            cells = []
            for a in ANCESTRY_ORDER:
                total = totals[gene][a]
                pct = 100 * counts[gene][a].get(allele, 0) / total if total else 0
                cells.append(f"{pct:.1f}%")
            md.append(f"| {allele} | " + " | ".join(cells) + f" | {pooled_n} |")
    md.append(f"\nFigure: `{fig_path}`\n")
    md_text = "\n".join(md)
    md_path = os.path.join(args.out_dir, "allele_frequency_report.md")
    with open(md_path, "w") as f:
        f.write(md_text)
    print(md_text)
    print(f"\n(written to {md_path} + 1 PNG in {args.out_dir})", file=sys.stderr)


if __name__ == "__main__":
    main()
