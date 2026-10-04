#!/usr/bin/env python3
"""Aggregate per-person TRUST4 results from a batch run into two outputs:

  1. VM-local, per-person detail (research_id, ancestry, n_cdr3, per-chain counts) --
     stays under ~/pipeline_outputs. Same privacy posture as the cohort file itself:
     bare research_ids are the open question flagged in build_experiment_d_cohort.py /
     DECISIONS.md, so this file is never committed to git.

  2. A DE-IDENTIFIED, ancestry-GROUP-level summary (no individual research_ids) -- the
     only thing this script writes into the repo, under ../results/.

Chain is assigned from the first gene TRUST4 actually called -- V, then J, then C -- the
same rule as embed_cdr3s.py's clonotype step (2026-10-04; previously V only, which left
CDR3s without a V call out of every per-chain count). The V-only counts are still computed
and written side by side to ../results/chain_assignment_check.csv so the change is visible.
Group spread is reported as 5th/95th percentiles: min and max are single-person values.

Usage:
  python3 aggregate_rnaseq_results.py <cohort.tsv> [--outdir ~/pipeline_outputs/rnaseq]
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from embed_cdr3s import _chain_of, _top_gene  # noqa: E402

CHAINS = ["TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cohort", help="cohort.tsv from build_rnaseq_cohort.py")
    ap.add_argument("--outdir", default=os.path.expanduser("~/pipeline_outputs/rnaseq"))
    ap.add_argument("--repo-results",
                    default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"),
                    help="where the DE-IDENTIFIED group summary is written (committed to git)")
    args = ap.parse_args()

    cohort = pd.read_csv(args.cohort, sep="\t", dtype=str)
    rows = []
    for _, r in cohort.iterrows():
        rid = r["research_id"]
        report = os.path.join(args.outdir, rid, f"{rid}_report.tsv")
        row = {"research_id": rid, "ancestry": r["ancestry"], "status": "missing"}
        if os.path.exists(report) and os.path.getsize(report) > 0:
            df = pd.read_csv(report, sep="\t", dtype=str)
            row["status"] = "ok"
            row["n_cdr3"] = len(df)
            genes = {g: (df[g].map(_top_gene) if g in df.columns else pd.Series("", index=df.index))
                     for g in ("V", "J", "C")}
            chain = pd.Series([_chain_of(v, j, c) for v, j, c in
                               zip(genes["V"], genes["J"], genes["C"])], index=df.index)
            chain_v = genes["V"].str[:3]
            for c in CHAINS:
                row[f"n_{c}"] = int((chain == c).sum())
                row[f"n_{c}_vonly"] = int((chain_v == c).sum())
        rows.append(row)

    detail = pd.DataFrame(rows)
    os.makedirs(args.outdir, exist_ok=True)
    detail_path = os.path.join(args.outdir, "batch_summary_detail.tsv")
    detail.to_csv(detail_path, sep="\t", index=False)
    print(f"Per-person detail (VM-local, NOT committed): {detail_path}", file=sys.stderr)

    ok = detail[detail["status"] == "ok"].copy()
    n_missing = len(detail) - len(ok)
    if n_missing:
        print(f"  {n_missing}/{len(detail)} people missing a result -- batch may still be "
              f"running, or those failed (check their .batch.log)", file=sys.stderr)

    if ok.empty:
        print("No completed results yet -- nothing to summarize.", file=sys.stderr)
        return

    group = ok.groupby("ancestry").agg(
        n_people=("research_id", "count"),
        cdr3_mean=("n_cdr3", "mean"),
        cdr3_median=("n_cdr3", "median"),
        cdr3_p05=("n_cdr3", lambda x: x.quantile(0.05)),
        cdr3_p95=("n_cdr3", lambda x: x.quantile(0.95)),
    ).round(1)
    for c in CHAINS:
        group[f"{c}_mean"] = ok.groupby("ancestry")[f"n_{c}"].mean().round(1)

    os.makedirs(args.repo_results, exist_ok=True)
    summary_path = os.path.join(args.repo_results, "rnaseq_cohort_ancestry_summary.csv")
    group.to_csv(summary_path)
    print(f"\nDe-identified group summary (safe to commit): {summary_path}", file=sys.stderr)
    print(group.to_string(), file=sys.stderr)

    check = pd.DataFrame({
        "chain": CHAINS,
        "mean_v_only": [ok[f"n_{c}_vonly"].mean() for c in CHAINS],
        "mean_v_j_c": [ok[f"n_{c}"].mean() for c in CHAINS],
    })
    check["pct_gained"] = 100 * (check["mean_v_j_c"] / check["mean_v_only"] - 1)
    check = check.round(2)
    check_path = os.path.join(args.repo_results, "chain_assignment_check.csv")
    check.to_csv(check_path, index=False)
    print(f"\nChain assignment, V-only vs V/J/C (safe to commit): {check_path}", file=sys.stderr)
    print(check.to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
