#!/usr/bin/env python3
"""RNA-seq cohort builder -- ancestry-stratified pick of N people with a confirmed-real,
existing RNA-seq BAM, for TRUST4 batch runs. Mirrors ../../../scripts/build_experiment_d_cohort.py's
pattern (same ANCESTRY_GROUPS, same deterministic-sort discipline, same ground-truth
existence check rather than trusting the manifest blindly -- a lesson the long-read
cohort builder paid for the hard way, see that script's own comments).

Joins, on the gcsfuse mount (~/mnt/aou-controlled):
  - RNA-seq manifest      v9/multiomics/rnaseq/manifest.tsv (research_id, markduplicates_bam_file_path, ...)
  - genetic-ancestry TSV  v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv (research_id, ancestry_pred, ...)

Picks people DETERMINISTICALLY (sorted by research_id within each ancestry group, first N)
so re-running yields the exact same cohort. Verifies each candidate's BAM actually exists
on the mount before including it -- lazily, only for the candidates being considered, not
the whole ~9k-row manifest, so this stays fast.

PRIVACY: output contains real research_ids. Same rule as build_experiment_d_cohort.py --
this file is VM-local only (default under ~/pipeline_outputs), never committed to git.

Usage:
  python3 build_rnaseq_cohort.py [--total 100] [--per-group N] [--out PATH] [--force] [--mount ~/mnt/aou-controlled]

--per-group, if given, overrides --total and applies uniformly (N per group, 6N total) --
same semantics as build_experiment_d_cohort.py. Otherwise --total is split as evenly as
possible across the 6 groups (remainder distributed to the first few, alphabetically).
"""
import argparse
import os
import sys

import pandas as pd

BUCKET_PREFIX = "gs://vwb-aou-datasets-controlled/"
ANCESTRY_GROUPS = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]  # TSV excludes "other"

RNASEQ_MANIFEST = "v9/multiomics/rnaseq/manifest.tsv"
ANCESTRY_TSV = "v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def require_cols(df, cols, label):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        die(f"{label}: expected column(s) {missing} not found. Actual columns: "
            f"{list(df.columns)}. If AoU renamed a column, update the constants at the "
            f"top of this script.")


def strip_bucket(uri):
    """gs://vwb-aou-datasets-controlled/pooled/... -> pooled/... (mount-relative)."""
    if not isinstance(uri, str) or not uri:
        return None
    if uri.startswith(BUCKET_PREFIX):
        return uri[len(BUCKET_PREFIX):]
    if uri.startswith("gs://"):
        return None  # some other bucket, can't resolve against our mount
    return uri


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--total", type=int, default=100,
                    help="Total cohort size, split as evenly as possible across 6 "
                         "ancestry groups (default 100).")
    ap.add_argument("--per-group", type=int, default=None,
                    help="Exact people per ancestry group -- overrides --total if given.")
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    ap.add_argument("--out", default=os.path.expanduser("~/pipeline_outputs/rnaseq/cohort.tsv"))
    ap.add_argument("--force", action="store_true",
                    help="Overwrite an existing cohort.tsv. Refused by default so an "
                         "in-flight batch run's cohort can't silently change under it.")
    ap.add_argument("--skip", type=int, default=0,
                    help="Skip this many confirmed-existing candidates per ancestry group "
                         "(same deterministic sort order) before picking. For building a "
                         "second, non-overlapping cohort on a different machine -- e.g. "
                         "--skip 17 to skip past a prior 100-person pick (17/17/17/17/16/16). "
                         "Uses the same value for every group regardless of how many that "
                         "group actually took, so it's always safe to over-skip slightly.")
    args = ap.parse_args()

    if os.path.exists(args.out) and not args.force:
        die(f"{args.out} already exists. Pass --force only if you're sure you want a "
            f"fresh cohort (and haven't started run_rnaseq_batch.sh yet).")

    if args.per_group is not None:
        per_group = {g: args.per_group for g in ANCESTRY_GROUPS}
    else:
        base, remainder = divmod(args.total, len(ANCESTRY_GROUPS))
        per_group = {g: base + (1 if i < remainder else 0)
                     for i, g in enumerate(ANCESTRY_GROUPS)}

    manifest_path = os.path.join(args.mount, RNASEQ_MANIFEST)
    anc_path = os.path.join(args.mount, ANCESTRY_TSV)
    for p in (manifest_path, anc_path):
        if not os.path.exists(p):
            die(f"not found: {p} -- is the gcsfuse mount up? (remount and `ls`-verify first)")

    print(f"Reading RNA-seq manifest: {manifest_path}", file=sys.stderr)
    rna = pd.read_csv(manifest_path, sep="\t", dtype=str)
    require_cols(rna, ["research_id", "markduplicates_bam_file_path"], "RNA-seq manifest")
    rna["bam_rel_path"] = rna["markduplicates_bam_file_path"].map(strip_bucket)
    rna = rna[rna["bam_rel_path"].notna()][["research_id", "bam_rel_path"]].drop_duplicates("research_id")
    print(f"  {len(rna)} people with an RNA-seq BAM path in the manifest", file=sys.stderr)

    print(f"Reading ancestry TSV: {anc_path}", file=sys.stderr)
    anc = pd.read_csv(anc_path, sep="\t", dtype=str,
                       usecols=lambda c: c in ("research_id", "ancestry_pred"))
    require_cols(anc, ["research_id", "ancestry_pred"], "ancestry TSV")
    anc["ancestry"] = anc["ancestry_pred"].astype(str).str.strip().str.upper()

    merged = rna.merge(anc[["research_id", "ancestry"]], on="research_id", how="inner")
    print(f"\nEligible (RNA-seq BAM + ancestry label): {len(merged)} people", file=sys.stderr)

    picks = []
    print("\n=== Per-ancestry availability (target -> confirmed-existing picked) ===",
          file=sys.stderr)
    for grp in ANCESTRY_GROUPS:
        target = per_group[grp]
        pool = merged[merged["ancestry"] == grp].sort_values("research_id")
        confirmed = []
        skipped = 0
        for _, row in pool.iterrows():
            if len(confirmed) >= target:
                break
            full_path = os.path.join(args.mount, row["bam_rel_path"])
            if os.path.exists(full_path):
                if skipped < args.skip:
                    skipped += 1
                    continue
                confirmed.append(row)
        flag = "  <-- SHORT" if len(confirmed) < target else ""
        skip_note = f", skipped {skipped}" if args.skip else ""
        print(f"  {grp}: target {target}, {len(pool)} in manifest, "
              f"{len(confirmed)} confirmed-existing picked{skip_note}{flag}", file=sys.stderr)
        if confirmed:
            picks.append(pd.DataFrame(confirmed))

    if not picks:
        die("no eligible people found in any ancestry group -- check the mount and manifest paths.")

    cohort = pd.concat(picks, ignore_index=True)[
        ["research_id", "ancestry", "bam_rel_path"]
    ].sort_values(["ancestry", "research_id"])

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    cohort.to_csv(args.out, sep="\t", index=False)
    print(f"\n=== Wrote {len(cohort)} people to {args.out} ===", file=sys.stderr)
    print(f"\nTest small first (recommended): head -4 {args.out} > {args.out}.test3.tsv",
          file=sys.stderr)
    print(f"Then compare:  bash run_rnaseq_batch.sh {args.out}.test3.tsv --jobs 1", file=sys.stderr)
    print(f"           vs  bash run_rnaseq_batch.sh {args.out}.test3.tsv --jobs 3", file=sys.stderr)


if __name__ == "__main__":
    main()
