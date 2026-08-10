#!/usr/bin/env python3
"""Full-cohort production run: completeness, coverage-per-ancestry, and operational stats.

Answers three questions for the supervisor report, cheaply (reads a handful of small TSVs,
no raw sequencing data, no 96-core VM needed -- see scripts/monitoring/README.md's "resize the
compute" note):
  1. Did the pipeline actually capture data correctly? (per-person gene-count completeness,
     timing-outlier sanity check per ENVIRONMENT.md quirk #18's "fast + tiny = starved input"
     lesson)
  2. Who got covered, and who didn't? (ancestry + platform coverage, intended cohort vs actual)
  3. Basic operational stats/distributions (timing, trim method mix).

Reads (all on the VM, per RESULTS_LOCATION.md -- this script itself never leaves aggregate
counts/rates, but run it ON the Workbench, not off it):
  ~/pipeline_outputs/immuannot_cohort_full.tsv   (person_id, platform, trim_tier, n_rows, ancestry_pred)
  ~/pipeline_outputs/immuannot_calls.tsv          (person_id, gene [HLA-prefixed], immuannot_1/_2)
  ~/pipeline_outputs/immuannot_timing.tsv         (person_id, hap, trim_method, trimmed_mb,
                                                    hap_total_seconds, ...)

Writes a markdown report + PNG figures under --out-dir (default
~/pipeline_outputs/production_analysis/completeness/). All aggregate-only (counts/rates/
distributions) -- no allele values, no bare person_ids in any output file, consistent with this
project's standing egress discipline (context/DECISIONS.md).

Usage (via `pixi run -e spechla` for pandas+matplotlib):
  python3 scripts/production_analysis/analyze_completeness_and_demographics.py
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")  # no display server on the VM
import matplotlib.pyplot as plt

GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS", "NA"]
ANCESTRY_COLOR = {
    "AFR": "#4C72B0", "AMR": "#DD8452", "EAS": "#55A868", "EUR": "#C44E52",
    "MID": "#8172B2", "SAS": "#937860", "NA": "#B0B0B0",
}
NULL = {"", "NA", "nan", "None", ".", "-"}

DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")


def is_null(v):
    return pd.isna(v) or str(v).strip() in NULL


def load_cohort(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    for c in ["person_id", "platform", "trim_tier"]:
        if c not in df.columns:
            sys.exit(f"FATAL: {path} missing expected column '{c}'. Actual columns: {list(df.columns)}")
    if "ancestry_pred" not in df.columns:
        df["ancestry_pred"] = None
    df["ancestry_pred"] = df["ancestry_pred"].fillna("NA")
    df.loc[~df["ancestry_pred"].isin(ANCESTRY_ORDER), "ancestry_pred"] = "NA"
    return df


def load_calls(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for c in ["person_id", "gene", "immuannot_1", "immuannot_2"]:
        if c not in df.columns:
            sys.exit(f"FATAL: {path} missing expected column '{c}'. Actual columns: {list(df.columns)}")
    df["gene_bare"] = df["gene"].str.replace("^HLA-", "", regex=True)
    return df


def load_timing(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for c in ["hap_total_seconds", "trimmed_mb", "trim_method"]:
        if c not in df.columns:
            sys.exit(f"FATAL: {path} missing expected column '{c}'. Actual columns: {list(df.columns)}")
    for c in ["hap_total_seconds", "trimmed_mb", "n_contigs", "padded_mb"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def build_person_summary(cohort, calls):
    """One row per person in the intended cohort: has this person produced ANY output, and how
    many of the 8 classical genes are 'complete' (a real call on BOTH haplotypes)?"""
    by_person = {pid: {} for pid in cohort["person_id"]}
    for _, r in calls.iterrows():
        pid, gene = r["person_id"], r["gene_bare"]
        by_person.setdefault(pid, {})[gene] = (r["immuannot_1"], r["immuannot_2"])

    n_genes_any_typed = calls.groupby("person_id")["gene_bare"].nunique()

    rows = []
    for _, crow in cohort.iterrows():
        pid = crow["person_id"]
        gene_calls = by_person.get(pid, {})
        has_output = len(gene_calls) > 0
        n_complete = n_partial = 0
        for g in GENES:
            pair = gene_calls.get(g)
            if pair is None:
                continue
            a1, a2 = pair
            null1, null2 = is_null(a1), is_null(a2)
            if not null1 and not null2:
                n_complete += 1
            elif null1 != null2:
                n_partial += 1
        rows.append({
            "person_id": pid, "platform": crow["platform"], "trim_tier": crow["trim_tier"],
            "ancestry_pred": crow["ancestry_pred"], "has_output": has_output,
            "n_classical_complete": n_complete, "n_classical_partial": n_partial,
            "n_genes_any_typed": int(n_genes_any_typed.get(pid, 0)),
        })
    return pd.DataFrame(rows)


def human_count(n):
    if n >= 1000:
        return f"{n / 1000:.1f}K"
    if n >= 100:
        return "hundreds"
    return str(n)


def plot_overview(summary, out_path):
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))

    # 1. Funnel: intended -> any output -> fully complete (8/8 classical genes)
    ax = axes[0, 0]
    n_intended = len(summary)
    n_any = (summary["has_output"]).sum()
    n_full = (summary["n_classical_complete"] == 8).sum()
    labels = ["Intended", "Any output", "8/8 genes\ncomplete"]
    vals = [n_intended, n_any, n_full]
    bars = ax.bar(labels, vals, color=["#B0B0B0", "#4C72B0", "#2C5C8A"])
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v}\n({100 * v / n_intended:.1f}%)",
                ha="center", va="bottom", fontsize=9)
    ax.set_title("Completion funnel")
    ax.spines[["top", "right"]].set_visible(False)

    # 2. Histogram of classical-gene completeness (0-8) among people with any output
    ax = axes[0, 1]
    have_output = summary[summary["has_output"]]
    counts = have_output["n_classical_complete"].value_counts().reindex(range(9), fill_value=0)
    ax.bar(counts.index, counts.values, color="#4C72B0")
    ax.set_xlabel("Classical genes complete (both haplotypes) out of 8")
    ax.set_ylabel("n people")
    ax.set_title("Per-person gene-count completeness")
    ax.spines[["top", "right"]].set_visible(False)

    # 3. Ancestry coverage: intended vs any-output vs fully-complete
    ax = axes[0, 2]
    width = 0.27
    x = range(len(ANCESTRY_ORDER))
    for i, (label, mask) in enumerate([
        ("Intended", pd.Series(True, index=summary.index)),
        ("Any output", summary["has_output"]),
        ("8/8 complete", summary["n_classical_complete"] == 8),
    ]):
        sub = summary[mask]
        vals = [(sub["ancestry_pred"] == a).sum() for a in ANCESTRY_ORDER]
        offs = [xi + (i - 1) * width for xi in x]
        ax.bar(offs, vals, width=width, label=label,
               color=["#CCCCCC", "#88AACC", "#2C5C8A"][i])
    ax.set_xticks(list(x))
    ax.set_xticklabels(ANCESTRY_ORDER)
    ax.set_ylabel("n people")
    ax.set_title("Coverage per ancestry")
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    # 4. Platform mix vs completion rate
    ax = axes[1, 0]
    platforms = sorted(summary["platform"].unique())
    any_rate = [100 * summary.loc[summary["platform"] == p, "has_output"].mean() for p in platforms]
    full_rate = [100 * (summary.loc[summary["platform"] == p, "n_classical_complete"] == 8).mean()
                 for p in platforms]
    xp = range(len(platforms))
    ax.bar([i - 0.18 for i in xp], any_rate, width=0.36, label="Any output", color="#88AACC")
    ax.bar([i + 0.18 for i in xp], full_rate, width=0.36, label="8/8 complete", color="#2C5C8A")
    ax.set_xticks(list(xp))
    ax.set_xticklabels(platforms)
    ax.set_ylabel("% of people")
    ax.set_ylim(0, 108)
    ax.set_title("Completion rate by platform")
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    # 5. trim_tier mix (shows the self_align_needed/Phase-2-disabled gap explicitly)
    ax = axes[1, 1]
    tier_counts = summary["trim_tier"].value_counts()
    ax.pie(tier_counts.values, labels=tier_counts.index, autopct="%1.0f%%",
           colors=["#4C72B0", "#DD8452", "#C44E52", "#55A868"][:len(tier_counts)])
    ax.set_title("Cohort trim_tier mix\n(self_align_needed = Phase 2, disabled this launch)")

    # 6. Breadth beyond the 8 classical genes
    ax = axes[1, 2]
    vals = have_output["n_genes_any_typed"]
    lo, hi = int(vals.min()), int(vals.max())
    ax.hist(vals, bins=np.arange(lo, hi + 2) - 0.5, color="#4C72B0")
    ax.set_xticks(range(lo, hi + 1))
    ax.set_xlabel("Distinct genes typed (any resolution, full Immuannot panel)")
    ax.set_ylabel("n people")
    ax.set_title("Breadth of typing beyond the 8 classical genes")
    ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("Full-cohort production run -- completeness & coverage overview", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def plot_timing(timing, out_path):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    ax = axes[0]
    ax.hist(timing["hap_total_seconds"].dropna() / 60, bins=40, color="#4C72B0")
    ax.set_xlabel("Per-haplotype wall time (min)")
    ax.set_ylabel("n haplotypes")
    ax.set_title("Timing distribution")
    ax.spines[["top", "right"]].set_visible(False)

    # Starved-input sanity check (ENVIRONMENT.md quirk #18): a haplotype that is BOTH unusually
    # fast AND unusually small is suspicious -- flag anything below the 5th percentile on both
    # axes simultaneously, a data-driven threshold rather than an arbitrary constant.
    ax = axes[1]
    sub = timing.dropna(subset=["hap_total_seconds", "trimmed_mb"])
    if len(sub) > 20:
        t_lo = sub["hap_total_seconds"].quantile(0.05)
        mb_lo = sub["trimmed_mb"].quantile(0.05)
        suspicious = (sub["hap_total_seconds"] <= t_lo) & (sub["trimmed_mb"] <= mb_lo)
        ax.scatter(sub.loc[~suspicious, "trimmed_mb"], sub.loc[~suspicious, "hap_total_seconds"] / 60,
                   s=6, alpha=0.35, color="#4C72B0", label="normal")
        ax.scatter(sub.loc[suspicious, "trimmed_mb"], sub.loc[suspicious, "hap_total_seconds"] / 60,
                   s=14, alpha=0.9, color="#C44E52", label=f"flagged ({suspicious.sum()})")
        ax.legend(fontsize=8, frameon=False)
        ax.set_title("Fast+tiny outlier check (quirk #18)")
    else:
        ax.text(0.5, 0.5, "not enough data", ha="center", transform=ax.transAxes)
    ax.set_xlabel("Trimmed input size (MB)")
    ax.set_ylabel("Wall time (min)")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[2]
    tm = timing["trim_method"].fillna("NA").value_counts()
    ax.bar(tm.index, tm.values, color="#4C72B0")
    ax.set_ylabel("n haplotypes")
    ax.set_title("Trim method used")
    ax.tick_params(axis="x", rotation=20)
    ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("Operational timing stats", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", default=os.path.join(DEFAULT_OUTROOT, "immuannot_cohort_full.tsv"))
    ap.add_argument("--calls", default=os.path.join(DEFAULT_OUTROOT, "immuannot_calls.tsv"))
    ap.add_argument("--timing", default=os.path.join(DEFAULT_OUTROOT, "immuannot_timing.tsv"))
    ap.add_argument("--out-dir", default=os.path.join(DEFAULT_OUTROOT, "production_analysis", "completeness"))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    cohort = load_cohort(args.cohort)
    calls = load_calls(args.calls)
    timing = load_timing(args.timing)
    summary = build_person_summary(cohort, calls)

    overview_path = os.path.join(args.out_dir, "completeness_overview.png")
    timing_path = os.path.join(args.out_dir, "timing_stats.png")
    plot_overview(summary, overview_path)
    plot_timing(timing, timing_path)

    n_intended, n_any, n_full = len(summary), summary["has_output"].sum(), (summary["n_classical_complete"] == 8).sum()
    md = [
        "# Production run -- completeness & coverage report\n",
        f"Intended cohort: **{n_intended}** people. "
        f"Any output: **{n_any}** ({100 * n_any / n_intended:.1f}%). "
        f"All 8 classical genes complete: **{n_full}** ({100 * n_full / n_intended:.1f}%).\n",
        "## Coverage per ancestry (intended vs any-output vs 8/8-complete)\n",
        "| Ancestry | Intended | Any output | 8/8 complete |",
        "|---|---|---|---|",
    ]
    for a in ANCESTRY_ORDER:
        n_i = (summary["ancestry_pred"] == a).sum()
        n_a = ((summary["ancestry_pred"] == a) & summary["has_output"]).sum()
        n_f = ((summary["ancestry_pred"] == a) & (summary["n_classical_complete"] == 8)).sum()
        if n_i == 0:
            continue
        md.append(f"| {a} | {n_i} | {n_a} ({100 * n_a / n_i:.0f}%) | {n_f} ({100 * n_f / n_i:.0f}%) |")
    md.append("\n## Platform mix\n")
    md.append("| Platform | n intended | any output | 8/8 complete |")
    md.append("|---|---|---|---|")
    for p in sorted(summary["platform"].unique()):
        sub = summary[summary["platform"] == p]
        md.append(f"| {p} | {len(sub)} | {sub['has_output'].sum()} ({100 * sub['has_output'].mean():.0f}%) "
                   f"| {(sub['n_classical_complete'] == 8).sum()} "
                   f"({100 * (sub['n_classical_complete'] == 8).mean():.0f}%) |")
    md.append(f"\nNote: `trim_tier == self_align_needed` (the 991 sequel2 people) is expected to "
               f"show ~0% completion -- Phase 2 was disabled for this production launch "
               f"(EXPERIMENTS.md, 2026-08-05) after Tier 3 failed testing. Not a pipeline bug.\n")
    md.append(f"\nFigures: `{overview_path}`, `{timing_path}`\n")
    md_text = "\n".join(md)
    md_path = os.path.join(args.out_dir, "completeness_report.md")
    with open(md_path, "w") as f:
        f.write(md_text)

    print(md_text)
    print(f"\n(written to {md_path} + 2 PNGs in {args.out_dir} -- aggregate only, no allele "
          f"values or bare IDs in either file)", file=sys.stderr)


if __name__ == "__main__":
    main()
