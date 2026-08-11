#!/usr/bin/env python3
"""Depth-confound check for the ancestry recovery gap found in the 100-person cohort
(../results/rnaseq_cohort_ancestry_summary.csv -- AFR 7,406 vs EAS 4,538 mean CDR3s,
~1.63x spread, consistent across all 5 major chain types).

Tests the mundane alternative explanation BEFORE trusting a reference-bias
interpretation: does AoU's own RNA-seq sequencing depth vary by ancestry in this
cohort, and does that alone explain the CDR3 recovery gap?

Joins, on the gcsfuse mount:
  - RNA-seq metadata   v9/multiomics/rnaseq/rnaseq_metadata.tsv (per-sample QC incl. depth)
against:
  - the cohort         (research_id, ancestry) from build_rnaseq_cohort.py
  - the batch results  ~/pipeline_outputs/rnaseq/batch_summary_detail.tsv (research_id, n_cdr3)

[MED] Column names for rnaseq_metadata.tsv are inferred from Marc's QC chart panel
titles (REPORT #4 slide 17: "reads_aligned_in_pairs"), not confirmed against the file's
actual header. This script checks live and fails with the real column list if wrong --
same discipline as build_rnaseq_cohort.py's require_cols().

[UNVERIFIED] Whether depth_read_pairs is stored as a raw count or already in millions
is unknown until this runs -- the script prints min/max of the raw column before
computing anything derived from it. Sanity-check that printout (values around
85-400 => already millions, don't divide again; values around 8.5e7-4e8 => raw count,
the /1e6 in cdr3_per_million_reads is correct as written) before trusting the final
normalized numbers.

Outputs TWO things:
  1. VM-local per-person joined detail (has research_ids) -- not committed.
  2. A DE-IDENTIFIED per-ancestry summary: mean depth, mean CDR3s, and CDR3s-per-
     million-read-pairs (depth-normalized) -- the number that actually answers the
     question. If the ancestry gap in raw CDR3 count mostly disappears once
     normalized, the gap is a cohort-depth-composition effect. If it persists, that's
     the finding worth taking seriously. Safe to commit, written to ../results/.

Usage:
  python3 check_depth_confound.py <cohort.tsv> [--mount ~/mnt/aou-controlled] [--detail PATH]
"""
import argparse
import os
import sys

import pandas as pd

RNASEQ_METADATA = "v9/multiomics/rnaseq/rnaseq_metadata.tsv"
DETAIL_DEFAULT = os.path.expanduser("~/pipeline_outputs/rnaseq/batch_summary_detail.tsv")

# Candidate column names, most-likely first -- see [MED] note above.
ID_CANDIDATES = ["research_id", "Research ID", "sampleid"]
DEPTH_CANDIDATES = ["reads_aligned_in_pairs", "aligned_read_pairs",
                     "number_of_aligned_read_pairs", "Aligned Read Pairs"]


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
    require_cols(detail, ["research_id", "status"], "batch_summary_detail.tsv")
    detail = detail[detail["status"] == "ok"].copy()
    detail["n_cdr3"] = detail["n_cdr3"].astype(float)

    print(f"Reading RNA-seq metadata: {meta_path}", file=sys.stderr)
    meta = pd.read_csv(meta_path, sep="\t", dtype=str)
    id_col = find_col(meta, ID_CANDIDATES, "rnaseq_metadata.tsv (ID column)")
    depth_col = find_col(meta, DEPTH_CANDIDATES, "rnaseq_metadata.tsv (depth column)")
    meta = meta.rename(columns={id_col: "research_id", depth_col: "depth_read_pairs"})
    meta["depth_read_pairs"] = pd.to_numeric(meta["depth_read_pairs"], errors="coerce")

    print(f"\nRaw depth_read_pairs range: min={meta['depth_read_pairs'].min():,.1f}  "
          f"max={meta['depth_read_pairs'].max():,.1f}", file=sys.stderr)
    print("SANITY CHECK before trusting anything below: values ~85-400 => column is "
          "ALREADY in millions (remove the /1e6 further down). Values ~8.5e7-4e8 => "
          "raw count, current code is correct as written.\n", file=sys.stderr)

    merged = cohort.merge(detail[["research_id", "n_cdr3"]], on="research_id", how="inner") \
                    .merge(meta[["research_id", "depth_read_pairs"]], on="research_id", how="inner")
    print(f"Joined: {len(merged)} people with cohort + result + depth", file=sys.stderr)
    if merged.empty:
        die("no overlap between cohort, results, and metadata -- check research_id formats match "
            "(e.g. one side zero-padded, leading/trailing whitespace).")

    os.makedirs(os.path.dirname(args.detail), exist_ok=True)
    joined_path = os.path.join(os.path.dirname(args.detail), "depth_confound_detail.tsv")
    merged.to_csv(joined_path, sep="\t", index=False)
    print(f"Per-person joined detail (VM-local, NOT committed): {joined_path}", file=sys.stderr)

    overall_r = merged["depth_read_pairs"].corr(merged["n_cdr3"])
    print(f"\nOverall correlation, depth vs CDR3 count: r = {overall_r:.3f}", file=sys.stderr)

    merged["cdr3_per_million_reads"] = merged["n_cdr3"] / (merged["depth_read_pairs"] / 1e6)

    summary = merged.groupby("ancestry").agg(
        n_people=("research_id", "count"),
        mean_depth_read_pairs=("depth_read_pairs", "mean"),
        mean_cdr3=("n_cdr3", "mean"),
        mean_cdr3_per_million_reads=("cdr3_per_million_reads", "mean"),
    ).round(1)

    os.makedirs(args.repo_results, exist_ok=True)
    out_path = os.path.join(args.repo_results, "rnaseq_depth_confound_check.csv")
    summary.to_csv(out_path)
    print(f"\nDe-identified summary (safe to commit): {out_path}", file=sys.stderr)
    print(summary.to_string(), file=sys.stderr)
    print(f"\nInterpretation: if the ancestry SPREAD in mean_cdr3_per_million_reads is "
          f"much smaller (proportionally) than the spread in mean_cdr3, depth explains "
          f"most of the gap -- a cohort-composition finding. If the spread in the "
          f"normalized column is still ~1.5x+ like the raw one, depth alone doesn't "
          f"explain it, and the reference-bias explanation becomes more credible.",
          file=sys.stderr)


if __name__ == "__main__":
    main()
