#!/usr/bin/env python3
"""LR x RNA-seq overlap counter -- the gate for the disease study (Task 2).

Answers, purely from manifests already on the gcsfuse mount (NO BigQuery, NO batch run,
NO BAM reads -- this is a lookup, not an experiment):

  1. How many unique people have RNA-seq?              (manifest.tsv, expected ~8,980)
  2. How many unique people have long-read WGS?        (expected 14,521 across 15,424 rows)
  3. How many have BOTH?                               <- THE NUMBER EVERYTHING DEPENDS ON
  4. What is the ancestry breakdown of that overlap?
  5. (--check-bams) How many of the overlap have a real, EXISTING long-read GRCh38 BAM?

Why (5) is separate and optional: ENVIRONMENT.md quirk #13 -- a row in the LR manifest does
NOT guarantee a usable aligned BAM (confirmed floor: only 2,763 of 14,521). Existence must
be checked file-by-file, never by path pattern -- that lesson cost three wrong "general
rules" in the HLA workstream. It costs a few minutes of mount stats, so it is opt-in.

Crucially: the DISEASE study does not need the BAM at all -- it needs the *people*, and
their EHR records. So run without --check-bams first; only the "re-call HLA from long reads"
follow-on needs the BAM-existence number.

PRIVACY: the per-person cohort file contains real research_ids and stays VM-local under
~/pipeline_outputs, never committed -- same rule as build_rnaseq_cohort.py and
build_experiment_d_cohort.py. Only the ancestry-group-level summary (no ids) goes to
results/ and is safe to commit.

Usage:
  python3 check_lr_rnaseq_overlap.py [--mount ~/mnt/aou-controlled] [--check-bams]
      [--out ~/pipeline_outputs/rnaseq/lr_rnaseq_overlap_cohort.tsv]
      [--summary-out ../results/lr_rnaseq_overlap_summary.csv]

Needs only pandas + the mount. Run from any pixi env.
"""
import argparse
import os
import sys

import pandas as pd

BUCKET_PREFIX = "gs://vwb-aou-datasets-controlled/"
ANCESTRY_GROUPS = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]  # the TSV excludes "other"

RNASEQ_MANIFEST = "v9/multiomics/rnaseq/manifest.tsv"
LR_MANIFEST = "v9/wgs/long_read/manifest.tsv"
ANCESTRY_TSV = "v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"

# The two columns known to hold an aligned GRCh38 BAM, checked in this order.
# grch38_bam is the raw aligned BAM; grch38_haplotagged_bam is a phased derivative --
# prefer the raw one. (ENVIRONMENT.md quirk #13.)
LR_BAM_COLS = ["grch38_bam", "grch38_haplotagged_bam"]


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
    return uri.lstrip("/")


def find_existing_bam(row, mount):
    """Return the first LR BAM path in this row that ACTUALLY EXISTS on the mount, else None.

    Deliberately checks existence rather than matching a path/folder pattern -- see the
    module docstring and ENVIRONMENT.md quirk #13. Slower, but correct by construction
    however many more release-folder naming surprises this manifest still holds.
    """
    for col in LR_BAM_COLS:
        rel = strip_bucket(row.get(col))
        if rel and os.path.exists(os.path.join(mount, rel)):
            return rel
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"),
                    help="gcsfuse mount point of vwb-aou-datasets-controlled")
    ap.add_argument("--check-bams", action="store_true",
                    help="Also verify which overlap people have a REAL existing long-read "
                         "GRCh38 BAM. Costs a few minutes of mount stats. Only needed for "
                         "the long-read HLA-recall follow-on, NOT for the disease study.")
    ap.add_argument("--out", default=os.path.expanduser(
                        "~/pipeline_outputs/rnaseq/lr_rnaseq_overlap_cohort.tsv"),
                    help="VM-local per-person output (contains real research_ids)")
    ap.add_argument("--summary-out", default=os.path.join(
                        os.path.dirname(os.path.abspath(__file__)),
                        "..", "results", "lr_rnaseq_overlap_summary.csv"),
                    help="De-identified group-level summary (safe to commit)")
    args = ap.parse_args()

    mount = os.path.expanduser(args.mount)
    if not os.path.isdir(mount):
        die(f"Mount not found: {mount}. The gcsfuse mount does not survive a VM restart "
            f"(ENVIRONMENT.md quirk #14) -- remount, `ls` a known path to warm it, then rerun.")

    # ---------- load the three tables ----------
    rna_path = os.path.join(mount, RNASEQ_MANIFEST)
    lr_path = os.path.join(mount, LR_MANIFEST)
    anc_path = os.path.join(mount, ANCESTRY_TSV)
    for p in (rna_path, lr_path, anc_path):
        if not os.path.exists(p):
            die(f"Missing expected file: {p}")

    rna = pd.read_csv(rna_path, sep="\t", dtype=str)
    lr = pd.read_csv(lr_path, sep="\t", dtype=str)
    anc = pd.read_csv(anc_path, sep="\t", dtype=str)

    require_cols(rna, ["research_id"], RNASEQ_MANIFEST)
    require_cols(lr, ["research_id"], LR_MANIFEST)
    require_cols(anc, ["research_id", "ancestry_pred"], ANCESTRY_TSV)

    rna_ids = set(rna["research_id"].dropna())
    lr_ids = set(lr["research_id"].dropna())
    both = rna_ids & lr_ids

    print("=" * 68)
    print("LR x RNA-seq overlap")
    print("=" * 68)
    print(f"RNA-seq manifest : {len(rna):>7,} rows   {len(rna_ids):>7,} unique people")
    print(f"lrWGS manifest   : {len(lr):>7,} rows   {len(lr_ids):>7,} unique people")
    print(f"OVERLAP (both)   : {len(both):>7,} people")

    # Is the overlap bigger than chance? Under independent sampling from the srWGS base
    # (535,662 per the v9 org PDF p.4), expected overlap = |LR| * |RNA| / base.
    SRWGS_BASE = 535_662
    expected = len(lr_ids) * len(rna_ids) / SRWGS_BASE
    print(f"\nExpected under independent draws from the {SRWGS_BASE:,} srWGS base: "
          f"{expected:,.0f}")
    if expected > 0:
        print(f"Observed / expected: {len(both) / expected:.2f}x")
        print("  >1 means the two sub-cohorts were co-selected (shared biospecimen or")
        print("  recruitment criteria); ~1 means they were drawn independently. Either way")
        print("  this is a fact about AoU's sampling that we should know before")
        print("  interpreting anything about these people.")

    if not both:
        print("\nNo overlap -- Task 2 as specified is not possible. Stop here and re-plan.")

    # ---------- ancestry breakdown ----------
    # Normalize case: the v9 TSV ships lowercase values ("eur", "afr", ...) but every
    # prior results file in this workstream is uppercase, because
    # build_rnaseq_cohort.py:114 upper-cases on read. Match it, or the overlap cohort
    # silently fails to join against the existing ancestry summaries.
    anc = anc[["research_id", "ancestry_pred"]].rename(columns={"ancestry_pred": "ancestry"})
    anc["ancestry"] = anc["ancestry"].astype(str).str.strip().str.upper()
    overlap_df = pd.DataFrame({"research_id": sorted(both)}).merge(anc, on="research_id",
                                                                  how="left")
    overlap_df["ancestry"] = overlap_df["ancestry"].fillna("unlabeled")

    print("\nAncestry breakdown of the overlap "
          "(genetically inferred, from ancestry_preds.tsv):")
    counts = overlap_df["ancestry"].value_counts()
    for grp in ANCESTRY_GROUPS + [g for g in counts.index if g not in ANCESTRY_GROUPS]:
        if grp in counts:
            print(f"  {grp:<10} {counts[grp]:>6,}")

    # ---------- optional: does a real LR BAM exist for these people? ----------
    if args.check_bams:
        print("\nChecking long-read BAM existence for the overlap people "
              "(file-by-file; a few minutes)...")
        lr_sub = lr[lr["research_id"].isin(both)].copy()
        have_cols = [c for c in LR_BAM_COLS if c in lr_sub.columns]
        if not have_cols:
            die(f"None of {LR_BAM_COLS} present in the LR manifest. Columns: "
                f"{list(lr_sub.columns)}")
        resolved = {}
        for _, row in lr_sub.iterrows():
            rid = row["research_id"]
            if rid in resolved:
                continue  # already resolved via another row for this same person
            rel = find_existing_bam(row, mount)
            if rel:
                resolved[rid] = rel
        overlap_df["lr_bam_rel_path"] = overlap_df["research_id"].map(resolved)
        n_bam = overlap_df["lr_bam_rel_path"].notna().sum()
        print(f"  Overlap people with a REAL existing GRCh38 long-read BAM: "
              f"{n_bam:,} / {len(both):,}")
        print("  (The disease study does NOT need these -- only a long-read HLA recall would.)")
    else:
        overlap_df["lr_bam_rel_path"] = pd.NA

    # ---------- write ----------
    out = os.path.expanduser(args.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    overlap_df.to_csv(out, sep="\t", index=False)
    print(f"\nVM-LOCAL per-person cohort (real research_ids -- DO NOT COMMIT):\n  {out}")

    summary = overlap_df.groupby("ancestry").agg(
        n_people=("research_id", "count"),
        n_with_lr_bam=("lr_bam_rel_path", lambda s: int(s.notna().sum())),
    ).reset_index()
    summary["n_rnaseq_total"] = len(rna_ids)
    summary["n_lr_total"] = len(lr_ids)
    summary["n_overlap_total"] = len(both)

    sout = os.path.abspath(os.path.expanduser(args.summary_out))
    os.makedirs(os.path.dirname(sout), exist_ok=True)
    summary.to_csv(sout, index=False)
    print(f"De-identified summary (safe to commit):\n  {sout}")
    print("\nNext: feed the VM-local cohort file to query_overlap_phenotypes.py "
          "(needs a Jupyter notebook -- BigQuery is not reachable from the plain shell).")


if __name__ == "__main__":
    main()
