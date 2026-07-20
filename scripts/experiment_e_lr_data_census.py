#!/usr/bin/env python3
"""Experiment E: exhaustive long-read data-format census.

Motivation (2026-07-20 session, Marc): ENVIRONMENT.md quirk #13 and DECISIONS.md's "Where are
the aligned BAMs for the other ~11,758 long-read people?" only ever checked TWO columns
(grch38_bam, grch38_haplotagged_bam) of the LR manifest for real file existence. DECISIONS.md's
"Assembly-based HLA typing" note found (by directory listing, not exhaustive per-column checking)
that some v9_delta rows resolve to a phased diploid assembly (assembly_hap1_fa/assembly_hap2_fa)
+ assembly-derived VCF instead -- a structurally different data product, not a missing file. This
script answers, exhaustively and per-person: for every one of the ~14,521 people with an LR
manifest row, which of ALL the manifest's columns point at a file that actually exists on the
mount? Single type, multiple types, or none of the types we can detect?

Follows the same rule that already cost real debugging time to learn (quirk #13): never guess
which columns are "the" file columns from their names or from the org PDF -- detect path-like
columns programmatically (do their non-null values look like gs:// URIs with a known genomics
extension?), then verify each candidate's existence directly against the mount. A column report
is printed before any existence-checking so you can sanity-check the auto-detected column list
against what you expect before committing to the (slow, gcsfuse-bound) existence pass.

Two multi-row-per-person subtleties already learned in Experiment D's cohort builder, reused here:
  - ~879 of 14,521 people have more than one manifest row (different center/platform).
  - A person's other row can resolve when their primary-looking one doesn't (person 1008366).
So this script reports BOTH a per-row and a per-person (union across rows) profile distribution.

Usage (any pixi env with pandas; read-only, no writes to the mount):
  python3 experiment_e_lr_data_census.py [--mount ~/mnt/aou-controlled] \
      [--out ~/pipeline_outputs/experiment_e/census.tsv] [--sample-only N]

  --sample-only N   Run column auto-detection + existence-check on only the first N rows, to
                     sanity-check output shape/timing before committing to the full ~15k-row pass.
"""
import argparse
import os
import sys
from collections import Counter

import pandas as pd

BUCKET_PREFIX = "gs://vwb-aou-datasets-controlled/"
LR_MANIFEST = "v9/wgs/long_read/manifest.tsv"

# Recognized genomics file extensions -- used only to auto-DETECT which columns are path-like,
# never to assume a column's semantic role. Add to this list if the manifest has a column whose
# values look like a real file path but use an extension not listed here (the column report below
# will surface any non-null string column that was NOT auto-detected as path-like, so nothing
# should be silently missed).
KNOWN_EXTENSIONS = (
    ".bam", ".bai", ".cram", ".crai", ".fastq", ".fastq.gz", ".fq.gz",
    ".fa", ".fasta", ".fa.gz", ".fasta.gz",
    ".vcf", ".vcf.gz", ".g.vcf.gz", ".gvcf", ".gvcf.gz",
    ".tbi", ".csi",
)


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def strip_bucket(uri):
    """gs://vwb-aou-datasets-controlled/pooled/... -> pooled/... (mount-relative).
    Mirrors build_experiment_d_cohort.py's strip_bucket -- same wrong-bucket trap guard."""
    if not isinstance(uri, str) or not uri:
        return None
    if "fc-aou-datasets-controlled" in uri:
        return None  # legacy Firecloud bucket -- not our mount (ENVIRONMENT.md wrong-bucket trap)
    if uri.startswith(BUCKET_PREFIX):
        return uri[len(BUCKET_PREFIX):]
    if uri.startswith("gs://"):
        return None  # some other bucket -- can't resolve against our mount
    return uri


def looks_like_path(series, sample_n=50):
    """A column is path-like if a sample of its non-null values are strings starting with
    gs:// (or already bucket-relative) AND end in a recognized genomics extension."""
    non_null = series.dropna()
    if non_null.empty:
        return False
    sample = non_null.astype(str).head(sample_n)
    hits = sum(
        1 for v in sample
        if (v.startswith("gs://") or "/" in v) and v.lower().endswith(KNOWN_EXTENSIONS)
    )
    return hits >= max(1, len(sample) // 2)  # majority of sampled values look like paths


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    ap.add_argument("--out", default=os.path.expanduser("~/pipeline_outputs/experiment_e/census.tsv"))
    ap.add_argument("--sample-only", type=int, default=None,
                    help="Only process the first N manifest rows (sanity-check run).")
    args = ap.parse_args()

    lr_path = os.path.join(args.mount, LR_MANIFEST)
    if not os.path.exists(lr_path):
        die(f"manifest not found: {lr_path} -- is the gcsfuse mount up? "
            f"(ENVIRONMENT.md quirk #11/#14: remount and `ls`-verify first.)")

    print(f"Reading LR manifest: {lr_path}", file=sys.stderr)
    lr = pd.read_csv(lr_path, sep="\t", dtype=str)
    if args.sample_only:
        lr = lr.head(args.sample_only)
        print(f"  --sample-only {args.sample_only}: using first {len(lr)} rows only", file=sys.stderr)

    n_rows, n_people = len(lr), lr["research_id"].nunique() if "research_id" in lr.columns else None
    print(f"  {n_rows} rows, {len(lr.columns)} columns, "
          f"{n_people if n_people is not None else '?'} unique research_id" , file=sys.stderr)

    # --- Column report: full column list + which were auto-detected as path-like ---
    print("\n=== All manifest columns (non-null count, auto-detected as path-like?) ===", file=sys.stderr)
    path_cols = []
    for col in lr.columns:
        is_path = looks_like_path(lr[col])
        if is_path:
            path_cols.append(col)
        non_null = lr[col].notna().sum()
        sample_val = lr[col].dropna().astype(str).head(1).tolist()
        sample_val = sample_val[0][:90] if sample_val else ""
        flag = "PATH-LIKE" if is_path else ""
        print(f"  {col:35s} non_null={non_null:6d}  {flag:10s} e.g. {sample_val}", file=sys.stderr)

    if not path_cols:
        die("No path-like columns auto-detected -- KNOWN_EXTENSIONS may need updating for this "
            "manifest's actual column contents (see the column report above for what's really "
            "in each column before assuming this is broken).")

    print(f"\n=== {len(path_cols)} path-like columns detected: {path_cols} ===", file=sys.stderr)
    print("Sanity-check this list against the column report above before trusting the "
          "existence-check pass below. Re-run with --sample-only 50 first if unsure.", file=sys.stderr)

    # --- Existence check: every path-like column, every row. Ground truth, not a pattern. ---
    print(f"\nVerifying real file existence for {n_rows} rows x {len(path_cols)} columns "
          f"against the mount (this is the slow part -- a few minutes per column-thousand-rows, "
          f"not hung).", file=sys.stderr)

    exist_cols = []
    for col in path_cols:
        exist_col = f"_exists_{col}"
        exist_cols.append(exist_col)
        results = []
        for i, v in enumerate(lr[col]):
            if i and i % 5000 == 0:
                print(f"    {col}: ...checked {i}/{n_rows}", file=sys.stderr)
            rel = strip_bucket(v)
            results.append(bool(rel) and os.path.exists(os.path.join(args.mount, rel)))
        lr[exist_col] = results
        n_exist = sum(results)
        print(f"  {col}: {n_exist}/{lr[col].notna().sum()} non-null values resolve to a real file "
              f"({n_exist}/{n_rows} of all rows)", file=sys.stderr)

    # --- Per-row profile: which path-like columns actually resolved, as a sorted tuple ---
    def row_profile(row):
        present = tuple(sorted(c[len("_exists_"):] for c in exist_cols if row[c]))
        return present if present else ("NONE_FOUND",)

    lr["profile"] = lr.apply(row_profile, axis=1)

    print("\n=== Per-ROW profile distribution (before per-person dedup) ===", file=sys.stderr)
    for profile, n in Counter(lr["profile"]).most_common():
        print(f"  {n:6d}  {' + '.join(profile)}", file=sys.stderr)

    if "platform" in lr.columns:
        print("\n=== Per-ROW profile x platform (top combos) ===", file=sys.stderr)
        combo_counts = Counter(zip(lr["platform"], lr["profile"]))
        for (platform, profile), n in sorted(combo_counts.items(), key=lambda kv: -kv[1])[:40]:
            print(f"  {n:6d}  platform={platform:15s} profile={' + '.join(profile)}", file=sys.stderr)

    if "center" in lr.columns:
        print("\n=== Per-ROW profile x center (top combos) ===", file=sys.stderr)
        combo_counts = Counter(zip(lr["center"], lr["profile"]))
        for (center, profile), n in sorted(combo_counts.items(), key=lambda kv: -kv[1])[:40]:
            print(f"  {n:6d}  center={center:8s} profile={' + '.join(profile)}", file=sys.stderr)

    # --- Per-PERSON profile: union of what's available across all of a person's rows ---
    # Same rescue logic as build_experiment_d_cohort.py: a person's *other* row can resolve
    # something their primary-looking row doesn't -- so union, don't just take one row.
    if "research_id" in lr.columns:
        def union_profile(group):
            present = set()
            for p in group["profile"]:
                if p != ("NONE_FOUND",):
                    present.update(p)
            return tuple(sorted(present)) if present else ("NONE_FOUND",)

        per_person = lr.groupby("research_id").apply(union_profile, include_groups=False)
        print(f"\n=== Per-PERSON profile distribution (union across rows, n={len(per_person)} unique people) ===",
              file=sys.stderr)
        for profile, n in Counter(per_person).most_common():
            pct = 100 * n / len(per_person)
            print(f"  {n:6d} ({pct:5.1f}%)  {' + '.join(profile)}", file=sys.stderr)

        n_none = sum(1 for p in per_person if p == ("NONE_FOUND",))
        n_multi_row = (lr.groupby("research_id").size() > 1).sum()
        print(f"\n  People with >1 manifest row: {n_multi_row}", file=sys.stderr)
        print(f"  People with NONE of the {len(path_cols)} detected file types resolving: {n_none} "
              f"-- these are the ones worth escalating (different/supplementary manifest? "
              f"ask Aleix per ENVIRONMENT.md quirk #8; check Workbench Data Dictionary directly) "
              f"rather than assumed to genuinely lack data.", file=sys.stderr)

        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        out_df = per_person.reset_index()
        out_df.columns = ["research_id", "profile"]
        out_df["profile"] = out_df["profile"].apply(lambda p: "+".join(p))
        out_df.to_csv(args.out, sep="\t", index=False)
        print(f"\n=== Wrote per-person profile table to {args.out} ===", file=sys.stderr)
    else:
        print("\nNo 'research_id' column found -- skipping per-person aggregation "
              "(check the column report above for the actual person-id column name).",
              file=sys.stderr)


if __name__ == "__main__":
    main()
