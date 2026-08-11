#!/usr/bin/env python3
"""Remaining ancestry-gap confound checks, continuing from check_depth_confound.py
(../results/rnaseq_100person_ancestry_writeup.md -- depth explains ~40% of the AFR-vs-EAS
gap, not all of it). Checks two more mundane explanations before treating the residual
gap as a reference/pipeline-bias finding, plus a significance test the writeup was
missing:

  1. RNA quality (RQS, 0-10 scale, higher = less degraded -- confirmed scale from the
     primary-source PDF, unlike depth's units which needed a live sanity check).
     Degraded RNA recovers less of everything, receptors included.
  2. Cell composition proxy -- T-cell chains (TRA+TRB+TRG+TRD) as a fraction of all
     chains recovered. Different blood cell mixes change how much lymphocyte signal is
     even in the library. Computed directly from data already in hand (batch_summary_
     detail.tsv's per-chain columns) -- no extra join needed for this one.
  3. Kruskal-Wallis test on raw CDR3 counts across the 6 ancestry groups -- the
     descriptive "looks real" claim in the writeup (AFR high / EAS low in 5/5 major
     chains) was never backed by an actual significance test, and with n=16-17 per
     group and large within-group spread, that matters.

Joins, on the gcsfuse mount:
  - RNA-seq metadata   v9/multiomics/rnaseq/rnaseq_metadata.tsv (RQS column)
against:
  - the cohort         (research_id, ancestry) from build_rnaseq_cohort.py
  - the batch results  ~/pipeline_outputs/rnaseq/batch_summary_detail.tsv

[MED] RQS column name inferred from Marc's QC chart panel title (REPORT #4 slide 17:
"rqs"). Checked live, fails loudly with the real column list if wrong -- same
discipline as check_depth_confound.py and build_rnaseq_cohort.py.

Outputs, same privacy split as every other script here:
  1. VM-local per-person joined detail -- not committed.
  2. DE-IDENTIFIED per-ancestry summary (RQS, T-cell fraction, CDR3) -- safe to commit,
     written to ../results/. Kruskal-Wallis result printed, not written into the CSV
     (it's a single test result, not a per-ancestry row).

Usage:
  python3 check_ancestry_confounds.py <cohort.tsv> [--mount ~/mnt/aou-controlled] [--detail PATH]
"""
import argparse
import os
import sys

import pandas as pd
from scipy.stats import kruskal

RNASEQ_METADATA = "v9/multiomics/rnaseq/rnaseq_metadata.tsv"
DETAIL_DEFAULT = os.path.expanduser("~/pipeline_outputs/rnaseq/batch_summary_detail.tsv")

ID_CANDIDATES = ["research_id", "Research ID", "sampleid"]
RQS_CANDIDATES = ["rqs", "RQS", "RNA Quality Score", "rna_quality_score"]
T_CHAINS = ["TRA", "TRB", "TRG", "TRD"]
B_CHAINS = ["IGH", "IGK", "IGL"]
ANCESTRY_GROUPS = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def require_cols(df, cols, label):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        die(f"{label}: expected column(s) {missing} not found. Actual columns: "
            f"{list(df.columns)}. If AoU renamed a column, update the constants at the "
            f"top of this script.")


def find_col(df, candidates, label):
    for c in candidates:
        if c in df.columns:
            return c
    die(f"{label}: none of {candidates} found. Actual columns: {list(df.columns)}. "
        f"Add the real column name to the candidates list at the top of this script.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cohort", help="cohort.tsv from build_rnaseq_cohort.py")
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    ap.add_argument("--detail", default=DETAIL_DEFAULT,
                    help="per-person detail from aggregate_rnaseq_results.py")
    ap.add_argument("--repo-results",
                    default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
    args = ap.parse_args()

    meta_path = os.path.join(args.mount, RNASEQ_METADATA)
    for p in (args.cohort, args.detail, meta_path):
        if not os.path.exists(p):
            die(f"not found: {p}")

    cohort = pd.read_csv(args.cohort, sep="\t", dtype=str)
    require_cols(cohort, ["research_id", "ancestry"], "cohort.tsv")

    detail = pd.read_csv(args.detail, sep="\t", dtype=str)
    require_cols(detail, ["research_id", "status", "n_cdr3"] + [f"n_{c}" for c in T_CHAINS + B_CHAINS],
                 "batch_summary_detail.tsv")
    detail = detail[detail["status"] == "ok"].copy()
    for c in ["n_cdr3"] + [f"n_{ch}" for ch in T_CHAINS + B_CHAINS]:
        detail[c] = detail[c].astype(float)

    t_total = detail[[f"n_{c}" for c in T_CHAINS]].sum(axis=1)
    b_total = detail[[f"n_{c}" for c in B_CHAINS]].sum(axis=1)
    detail["tcell_fraction"] = t_total / (t_total + b_total)

    print(f"Reading RNA-seq metadata: {meta_path}", file=sys.stderr)
    meta = pd.read_csv(meta_path, sep="\t", dtype=str)
    id_col = find_col(meta, ID_CANDIDATES, "rnaseq_metadata.tsv (ID column)")
    rqs_col = find_col(meta, RQS_CANDIDATES, "rnaseq_metadata.tsv (RQS column)")
    meta = meta.rename(columns={id_col: "research_id", rqs_col: "rqs"})
    meta["rqs"] = pd.to_numeric(meta["rqs"], errors="coerce")
    print(f"RQS range: min={meta['rqs'].min():.2f}  max={meta['rqs'].max():.2f}  "
          f"(expected ~0-10 per the primary-source PDF -- sanity check this before trusting below)",
          file=sys.stderr)

    merged = cohort.merge(
        detail[["research_id", "n_cdr3", "tcell_fraction"]], on="research_id", how="inner"
    ).merge(meta[["research_id", "rqs"]], on="research_id", how="inner")
    print(f"\nJoined: {len(merged)} people with cohort + result + RQS", file=sys.stderr)
    if merged.empty:
        die("no overlap between cohort, results, and metadata -- check research_id formats match.")

    os.makedirs(os.path.dirname(args.detail), exist_ok=True)
    joined_path = os.path.join(os.path.dirname(args.detail), "ancestry_confounds_detail.tsv")
    merged.to_csv(joined_path, sep="\t", index=False)
    print(f"Per-person joined detail (VM-local, NOT committed): {joined_path}", file=sys.stderr)

    r_rqs = merged["rqs"].corr(merged["n_cdr3"])
    r_tcell = merged["tcell_fraction"].corr(merged["n_cdr3"])
    print(f"\nCorrelation, RQS vs CDR3 count: r = {r_rqs:.3f}", file=sys.stderr)
    print(f"Correlation, T-cell fraction vs CDR3 count: r = {r_tcell:.3f}", file=sys.stderr)

    summary = merged.groupby("ancestry").agg(
        n_people=("research_id", "count"),
        mean_rqs=("rqs", "mean"),
        mean_tcell_fraction=("tcell_fraction", "mean"),
        mean_cdr3=("n_cdr3", "mean"),
    ).round(3)

    os.makedirs(args.repo_results, exist_ok=True)
    out_path = os.path.join(args.repo_results, "rnaseq_remaining_confounds_check.csv")
    summary.to_csv(out_path)
    print(f"\nDe-identified summary (safe to commit): {out_path}", file=sys.stderr)
    print(summary.to_string(), file=sys.stderr)

    groups = [merged.loc[merged["ancestry"] == g, "n_cdr3"].values
              for g in ANCESTRY_GROUPS if (merged["ancestry"] == g).any()]
    stat, pval = kruskal(*groups)
    print(f"\nKruskal-Wallis across ancestry groups (raw CDR3 counts): "
          f"H = {stat:.2f}, p = {pval:.4g}", file=sys.stderr)
    print("(p < 0.05 => the between-group differences are unlikely to be pure chance, "
          "given the within-group spread. This does NOT by itself say what's causing "
          "them -- reference bias and confounds are equally consistent with a significant p.)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
