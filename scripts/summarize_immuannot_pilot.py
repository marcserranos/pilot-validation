#!/usr/bin/env python3
"""Aggregate-only summary of the 60-person Immuannot pilot run (2026-07-21/22). Reads the three
VM-local files produced by run_immuannot_person.py + build_experiment_d_cohort.py -- none of which
leave the VM -- and writes ONLY aggregate stats (completion rate, timing distribution, trim-method
mix, per-gene call-rate) with no bare person_ids and no genotypes. Same discipline as
analyze_experiment_d.py / analyze_experiment_d_field_cascade.py: paste/commit this script's output,
never the raw inputs.

Inputs (all VM-local, never committed):
  ~/pipeline_outputs/experiment_d/cohort.tsv    -- person_id, ancestry (the attempted 60)
  ~/pipeline_outputs/immuannot_calls.tsv        -- person_id, gene, immuannot_1, immuannot_2 (REAL
                                                    genotype calls -- participant data)
  ~/pipeline_outputs/immuannot_timing.tsv       -- person_id, hap, trim_method, sizes, per-stage
                                                    timings (bare person_ids -- operational, not
                                                    genotypes, but still kept VM-local per the
                                                    same "bare ids in a public repo" open question
                                                    as cohort.tsv, DECISIONS.md)

Usage (any pixi env with pandas):
  python3 scripts/summarize_immuannot_pilot.py [--outroot ~/pipeline_outputs] \\
      [--analysis-dir ~/pipeline_outputs/immuannot_pilot/analysis]
"""
import argparse
import os
import sys

import pandas as pd

CLASSICAL_8 = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1"]


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def pct(n, d):
    return f"{n}/{d} ({100*n/d:.0f}%)" if d else "-"


def stats_line(series):
    if series.empty:
        return "n=0"
    return (f"n={len(series)}, mean={series.mean():.1f}, median={series.median():.1f}, "
            f"min={series.min():.1f}, max={series.max():.1f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outroot", default=os.path.expanduser("~/pipeline_outputs"))
    ap.add_argument("--cohort", default=None,
                    help="Default: <outroot>/experiment_d/cohort.tsv")
    ap.add_argument("--analysis-dir", default=None,
                    help="Default: <outroot>/immuannot_pilot/analysis")
    args = ap.parse_args()

    cohort_path = args.cohort or os.path.join(args.outroot, "experiment_d", "cohort.tsv")
    calls_path = os.path.join(args.outroot, "immuannot_calls.tsv")
    timing_path = os.path.join(args.outroot, "immuannot_timing.tsv")
    for p in (cohort_path, calls_path, timing_path):
        if not os.path.exists(p):
            die(f"not found: {p}")

    cohort = pd.read_csv(cohort_path, sep="\t", dtype=str)
    calls = pd.read_csv(calls_path, sep="\t", dtype=str)
    timing = pd.read_csv(timing_path, sep="\t", dtype=str)
    for col in ("n_contigs", "whole_contig_mb", "padded_mb", "trimmed_mb", "resolve_seconds",
                "contig_lookup_seconds", "trim_seconds", "immuannot_seconds", "hap_total_seconds"):
        if col in timing.columns:
            timing[col] = pd.to_numeric(timing[col], errors="coerce")

    n_attempted = cohort["person_id"].nunique()
    completed_ids = set(calls["person_id"].astype(str).unique())
    n_completed = len(completed_ids)
    n_missing = n_attempted - n_completed

    genes_per_person = calls.groupby("person_id").size()

    classical_rate = {}
    for gene in CLASSICAL_8:
        sub = calls[calls["gene"] == gene]
        n_h1 = (sub["immuannot_1"] != "NA").sum()
        n_h2 = (sub["immuannot_2"] != "NA").sum()
        classical_rate[gene] = (n_h1 + n_h2, 2 * n_completed)

    trim_counts = timing["trim_method"].value_counts(dropna=False)
    hap_rows_with_result = timing.dropna(subset=["immuannot_seconds"])

    parts = [
        "# Immuannot 60-person pilot -- aggregate summary",
        "",
        "Aggregate-only: no bare person_ids, no genotype/allele values. Raw calls and per-person "
        "timing stay on the Workbench VM.",
        "",
        f"**Completion: {pct(n_completed, n_attempted)}** of the attempted cohort produced calls "
        f"({n_missing} did not -- see 'assembly data available' caveat below).",
        "",
        "## Genes called per person",
        f"- {stats_line(genes_per_person)}",
        "",
        "## Classical-8 call rate (out of person x haplotype = 2 per completed person)",
        "",
        "| Gene | Called | Rate |",
        "|---|---|---|",
    ]
    for gene, (n, d) in classical_rate.items():
        parts.append(f"| {gene} | {n}/{d} | {100*n/d:.0f}% |")

    parts += [
        "",
        "## Trim method used (per haplotype attempt)",
        "",
        "| Method | Count |",
        "|---|---|",
    ]
    for method, n in trim_counts.items():
        label = method if isinstance(method, str) else "not reached (skipped earlier)"
        parts.append(f"| {label} | {n} |")

    parts += [
        "",
        "## Timing distribution (minutes, per haplotype attempt that reached immuannot.sh)",
        "",
        f"- trim (samtools faidx), seconds: {stats_line(hap_rows_with_result['trim_seconds'])}",
        f"- immuannot.sh, minutes: {stats_line(hap_rows_with_result['immuannot_seconds'] / 60)}",
        f"- hap total, minutes: {stats_line(hap_rows_with_result['hap_total_seconds'] / 60)}",
        "",
        f"Estimated person-level total (sum of both haps' hap_total, minutes): "
        f"{stats_line(timing.groupby('person_id')['hap_total_seconds'].sum().dropna() / 60)}",
        "",
        "## Caveats",
        "- 'Missing' people are not assumed to be a bug -- most likely lack assembly data on their "
        "platform (only revio/sequel2e/sequel2 do, per reports/lr_data_census/README.md). This "
        "summary doesn't break down *why* each missing person was skipped -- that reasoning was "
        "only printed to the run's own stderr, not persisted to a file.",
        "- Real per-person calls and timing remain on the VM "
        "(~/pipeline_outputs/immuannot_calls.tsv, immuannot_timing.tsv) -- never committed or "
        "pasted in full.",
    ]
    md = "\n".join(parts) + "\n"

    adir = args.analysis_dir or os.path.join(args.outroot, "immuannot_pilot", "analysis")
    os.makedirs(adir, exist_ok=True)
    md_path = os.path.join(adir, "immuannot_pilot_summary.md")
    with open(md_path, "w") as f:
        f.write(md)

    print(md)
    print(f"\n(written to {md_path} -- aggregate-only, safe to paste back or bring into the repo)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
