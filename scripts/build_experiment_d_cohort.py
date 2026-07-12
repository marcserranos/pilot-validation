#!/usr/bin/env python3
"""Experiment D cohort builder. Derives the ancestry-stratified 3-way-comparison cohort
directly from the mounted v9 manifests on the VM -- no notebook / BigQuery step needed,
so the whole of Experiment D can run from the shell.

Joins, on the gcsfuse mount (~/mnt/aou-controlled), the three files whose exact paths +
schemas are the confirmed ones in ENVIRONMENT.md ("Confirmed data locations"):
  - srWGS CRAM manifest    v9/wgs/cram/manifest.csv             (person_id, cram_uri, cram_index_uri)
  - lrWGS manifest         v9/wgs/long_read/manifest.tsv        (research_id, ..., grch38_bam, ...)
  - genetic-ancestry TSV   .../aux/ancestry/ancestry_preds.tsv  (research_id, ancestry_pred, ...)

Keeps only people who:
  (a) have a REAL /revio/ long-read BAM -- ENVIRONMENT.md quirk #13: the has-LR flag and even
      a manifest grch38_bam entry are NOT enough; the path must literally contain /revio/ or it
      may point at a pacbio assembly folder with no aligned BAM.
  (b) have a CRAM in the srWGS manifest.
  (c) have a genetic-ancestry label (AFR/AMR/EAS/EUR/MID/SAS -- the TSV excludes "other").

Picks --per-group people per ancestry, DETERMINISTICALLY (sorted by id, first N) so re-running
yields the exact same cohort -- important because run_experiment_d.sh is resumable across
multiple nights and must see a stable cohort. Groups with fewer than N eligible people take all
available and are flagged (MID/SAS are the likely-thin ones in the LR sub-cohort).

Writes a reviewable TSV: person_id, ancestry, cram_rel_path, lr_bam_rel_path -- mount-relative
paths ready for run_experiment_d.sh. Aggregate/operational only (ids + paths + ancestry) -- no
genotypes, so this file is not in the genotype-quarantine class (still, ids in a public repo is
the open question in DECISIONS.md -- this file lives under ~/pipeline_outputs, not in the repo).

Usage (from any pixi env -- only needs pandas + the mount):
  python3 build_experiment_d_cohort.py [--per-group 10] [--out PATH] [--force] \
      [--mount ~/mnt/aou-controlled]
"""
import argparse
import os
import sys

import pandas as pd

BUCKET_PREFIX = "gs://vwb-aou-datasets-controlled/"
ANCESTRY_GROUPS = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]  # TSV excludes "other" (ENVIRONMENT.md)

# Relative-to-mount manifest paths -- the confirmed v9 locations (ENVIRONMENT.md). Resolve
# physical CRAM/BAM paths only via these manifests, never hand-build into pooled/ (quirk in
# ENVIRONMENT.md: physical files are spread across v7/v8/v9 base/delta dirs).
CRAM_MANIFEST = "v9/wgs/cram/manifest.csv"
LR_MANIFEST = "v9/wgs/long_read/manifest.tsv"
ANCESTRY_TSV = "v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def require_cols(df, cols, label):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        die(f"{label}: expected column(s) {missing} not found. Actual columns: "
            f"{list(df.columns)}. If AoU renamed a column, update the constants at the top "
            f"of this script (schemas are the confirmed ones in ENVIRONMENT.md).")


def strip_bucket(uri):
    """gs://vwb-aou-datasets-controlled/pooled/... -> pooled/... (mount-relative)."""
    if not isinstance(uri, str) or not uri:
        return None
    if "fc-aou-datasets-controlled" in uri:
        # Wrong-bucket trap (ENVIRONMENT.md): legacy Firecloud naming, not our mount.
        return None
    if uri.startswith(BUCKET_PREFIX):
        return uri[len(BUCKET_PREFIX):]
    if uri.startswith("gs://"):
        # Some other bucket -- can't resolve against our mount.
        return None
    return uri  # already relative


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-group", type=int, default=10,
                    help="People per ancestry group (default 10). Groups with fewer eligible "
                         "take all available.")
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    ap.add_argument("--out", default=os.path.expanduser("~/pipeline_outputs/experiment_d/cohort.tsv"))
    ap.add_argument("--force", action="store_true",
                    help="Overwrite an existing cohort.tsv. Refused by default so an in-flight "
                         "multi-night run's cohort can't silently change under it.")
    args = ap.parse_args()

    if os.path.exists(args.out) and not args.force:
        die(f"{args.out} already exists. A resumable run may be using it -- re-deriving could "
            f"change which people are in the cohort mid-run. Pass --force only if you're sure "
            f"you want a fresh cohort (and haven't started run_experiment_d.sh yet).")

    cram_path = os.path.join(args.mount, CRAM_MANIFEST)
    lr_path = os.path.join(args.mount, LR_MANIFEST)
    anc_path = os.path.join(args.mount, ANCESTRY_TSV)
    for p in (cram_path, lr_path, anc_path):
        if not os.path.exists(p):
            die(f"manifest not found: {p} -- is the gcsfuse mount up? (ENVIRONMENT.md quirk #11: "
                f"remount and `ls`-verify before running anything that reads the bucket.)")

    # --- long-read manifest first (smallest, ~15k rows; it's the binding constraint) ---
    # BAM resolution (ENVIRONMENT.md quirk #13, full investigation 2026-07-11): after TWO rounds of
    # path-pattern guessing each missed a real source (a `/revio/`-only filter missed `sequel2`'s
    # real BAM under a different column; a `/v8_delta/`-folder filter then missed it AGAIN, because
    # `sequel2`'s real files apparently live under yet another release folder, not literally
    # `v8_delta`) -- this resolver stops pattern-matching release-folder names entirely and instead
    # verifies each candidate file's existence directly against the mount. Ground truth, not a guess.
    # Checks both known BAM-holding columns (`grch38_bam`, then `grch38_haplotagged_bam` -- confirmed
    # 2026-07-11 that at least the `sequel2` platform publishes its real BAM under the latter).
    # ~879 of 14,521 people have multiple manifest rows (different center/platform per AoU's own
    # docs); a person's newest row can be unusable while an older row for the same person resolves
    # fine (confirmed: person 1008366's `revio` row is a dead path, but their separate `ont-r10.4.1`
    # row is real) -- so resolution happens per ROW first, then dedupes to any resolving row per person.
    # This does one os.path.exists() per candidate column per row (up to ~2x15,424 mount stats) --
    # takes a few minutes over gcsfuse, not instant like the old pattern-match version. That's expected.
    print(f"Reading LR manifest: {lr_path}", file=sys.stderr)
    lr = pd.read_csv(lr_path, sep="\t", dtype=str)
    require_cols(lr, ["research_id", "center", "platform", "grch38_bam", "grch38_haplotagged_bam"], "LR manifest")
    n_lr_total, n_people_total = len(lr), lr["research_id"].nunique()
    print(f"  Verifying real file existence for {n_lr_total} rows against the mount "
          f"(ground truth, not a path-pattern guess -- this will take a few minutes, not hung)",
          file=sys.stderr)

    def find_existing_bam(row):
        for col in ("grch38_bam", "grch38_haplotagged_bam"):
            v = row.get(col)
            if not isinstance(v, str) or not v:
                continue
            rel_path = strip_bucket(v)
            if rel_path is None:
                continue
            if os.path.exists(os.path.join(args.mount, rel_path)):
                return v, col
        return None, "no_existing_file_found"

    uris, sources = [], []
    for i, row in enumerate(lr.to_dict("records")):
        if i and i % 2000 == 0:
            print(f"    ...checked {i}/{n_lr_total}", file=sys.stderr)
        uri, src = find_existing_bam(row)
        uris.append(uri)
        sources.append(src)
    lr["_bam_uri"], lr["bam_source"] = uris, sources
    print("  BAM resolution by source (per row, before per-person dedup):", file=sys.stderr)
    for src, n in lr["bam_source"].value_counts().items():
        print(f"    {src}: {n}", file=sys.stderr)
    lr = lr[lr["_bam_uri"].notna()].copy()
    lr = lr.rename(columns={"research_id": "person_id"})
    lr["lr_bam_rel_path"] = lr["_bam_uri"].map(strip_bucket)
    lr = lr[lr["lr_bam_rel_path"].notna()]
    # Deterministic per-person pick: a person can have multiple resolving rows (rare) or one
    # resolving + one non-resolving (the 1008366 case) -- sort for reproducibility, keep first.
    lr = lr.sort_values(["person_id", "center", "platform"])[
        ["person_id", "lr_bam_rel_path", "bam_source"]
    ].drop_duplicates("person_id", keep="first")
    print(f"  {n_lr_total} LR rows ({n_people_total} unique people) -> "
          f"{len(lr)} people with a confirmed-usable aligned BAM after per-person dedup", file=sys.stderr)

    # --- ancestry labels ---
    print(f"Reading ancestry TSV: {anc_path}", file=sys.stderr)
    anc = pd.read_csv(anc_path, sep="\t", dtype=str, usecols=lambda c: c in ("research_id", "ancestry_pred"))
    require_cols(anc, ["research_id", "ancestry_pred"], "ancestry TSV")
    anc = anc.rename(columns={"research_id": "person_id", "ancestry_pred": "ancestry"})
    # The TSV ships lowercase labels (afr/amr/...) -- normalize to the uppercase group names
    # used everywhere downstream (ANCESTRY_GROUPS here, GROUPS in analyze_experiment_d.py).
    anc["ancestry"] = anc["ancestry"].astype(str).str.strip().str.upper()

    # --- CRAM manifest (largest, ~535k rows; read only the two columns we need) ---
    print(f"Reading CRAM manifest: {cram_path}", file=sys.stderr)
    cram = pd.read_csv(cram_path, dtype=str, usecols=lambda c: c in ("person_id", "cram_uri"))
    require_cols(cram, ["person_id", "cram_uri"], "CRAM manifest")
    cram["cram_rel_path"] = cram["cram_uri"].map(strip_bucket)
    cram = cram[cram["cram_rel_path"].notna()][["person_id", "cram_rel_path"]]
    cram = cram.drop_duplicates("person_id")

    # --- inner-join all three: must have a confirmed-usable LR BAM + CRAM + ancestry ---
    merged = lr.merge(cram, on="person_id", how="inner").merge(anc, on="person_id", how="inner")
    print(f"\nEligible (confirmed-usable LR BAM + CRAM + ancestry): {len(merged)} people", file=sys.stderr)
    print("  by bam_source:", dict(merged["bam_source"].value_counts()), file=sys.stderr)

    # --- availability report + deterministic pick ---
    picks = []
    print("\n=== Per-ancestry availability (eligible / picked) ===", file=sys.stderr)
    for grp in ANCESTRY_GROUPS:
        pool = merged[merged["ancestry"] == grp].sort_values("person_id")
        take = pool.head(args.per_group)
        flag = "  <-- SHORT" if len(pool) < args.per_group else ""
        print(f"  {grp}: {len(pool):>6} eligible, picking {len(take)}{flag}", file=sys.stderr)
        picks.append(take)
    other = sorted(set(merged["ancestry"]) - set(ANCESTRY_GROUPS))
    if other:
        print(f"  (ignored non-standard ancestry labels present in data: {other})", file=sys.stderr)

    cohort = pd.concat(picks, ignore_index=True)[
        ["person_id", "ancestry", "cram_rel_path", "lr_bam_rel_path", "bam_source"]
    ].sort_values(["ancestry", "person_id"])

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    cohort.to_csv(args.out, sep="\t", index=False)
    print(f"\n=== Wrote {len(cohort)} people to {args.out} ===", file=sys.stderr)
    print(f"Review it, then launch:  bash run_experiment_d.sh {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
