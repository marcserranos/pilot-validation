#!/usr/bin/env python3
"""Immuannot-specific cohort builder. Filters the v9 lrWGS manifest to people who actually have
what run_immuannot_person.py needs -- existence-verified, not inferred from a platform label.

REWRITTEN 2026-08-04 for the full-cohort production run (production_orchestrator/BRIEF.md), from
the original small-test-batch version (git history has it): that version defaulted -n=40 and
stopped as soon as N verified people were found, filtered to platform=="revio" only, and (like
build_experiment_d_cohort.py before quirk #13's full lesson) took the FIRST manifest row per
research_id without checking whether that person's OTHER row might resolve instead. All three are
fixed here:
  1. **No default cap.** -n now means "stop after N eligible people found" only if explicitly
     passed (e.g. for a quick smoke-test cohort); omitted, it enumerates the entire eligible
     population. Every eligible row/person is existence-checked, not just the first N encountered.
  2. **Platform filter widened to revio + sequel2e + sequel2** (BRIEF.md), not just revio.
     revio/sequel2e carry the full assembly+alignment suite (reports/lr_data_census/README.md);
     sequel2 (991 people, ~95% AFR) has the assembly FASTA but NOT the aln-to-hg38 BAM/PAF that
     Tiers 1/2 need -- tracked separately below via the `trim_tier` column, not silently dropped
     and not silently included at full production-launch weight either (see "sequel2" section).
  3. **Row-level check, then person-level union** -- same discipline as ENVIRONMENT.md quirk #13
     and lr_manifest_format_census.py: ~879/14,521 people have >1 manifest row, and a person's
     newer row can dead-end while an older one resolves. Every candidate ROW is existence-checked
     independently; a person is eligible if ANY of their rows resolves the FASTA requirement, and
     `trim_tier` reflects the BEST tier available across all their rows (paf_region beats
     bam_whole_contig beats self_align_needed).
  4. **Parallelized + checkpointed/resumable**, same pattern as lr_manifest_format_census.py
     (ThreadPoolExecutor -- this is FUSE-I/O-bound, not CPU-bound; chunked writes with fsync so a
     restart resumes instead of re-checking from row 0). At ~13,280 revio+sequel2e+sequel2 rows x
     6 columns each, a serial loop was assumed-not-measured to be too slow at this scale (BRIEF.md
     explicit ask: "verify by timing a real run, don't assume").

The one HARD requirement is assembly_hap{1,2}_fa (Immuannot's real dependency -- the .fa is what
gets typed; aln-to-hg38 BAM/PAF are this project's own trim-step optimization, not a real
Immuannot requirement, per DECISIONS.md "Assembly-based HLA typing"). A person with the FASTA but
missing the alignment files is still eligible, just routed to `trim_tier=self_align_needed`
instead of `paf_region`/`bam_whole_contig` -- run_immuannot_person.py's Tier 3
(--enable-self-align-fallback) is what would actually process them, and that tier is UNTESTED AT
SCALE as of 2026-08-04 (see BRIEF.md / run_immuannot_person.py). This script surfaces the split;
it does not itself decide whether to launch on self_align_needed people.

Usage (from ~/repos/pilot-validation, inside `pixi shell -e specimmune` or `pixi run -e specimmune --`):
  python3 scripts/build_immuannot_cohort.py [-n N] [--out PATH] [--force] [--restart]
      [--mount ~/mnt/aou-controlled] [--platforms revio,sequel2e,sequel2]
"""
import argparse
import concurrent.futures
import os
import sys
from collections import Counter

import pandas as pd

LR_MANIFEST = "v9/wgs/long_read/manifest.tsv"
FA_COLS = ["assembly_hap1_fa", "assembly_hap2_fa"]
ALN_COLS = ["assembly_hap1_aln2_hg38_bam", "assembly_hap2_aln2_hg38_bam"]
PAF_COLS = ["assembly_hap1_aln2_hg38_paf", "assembly_hap2_aln2_hg38_paf"]
ALL_CHECK_COLS = FA_COLS + ALN_COLS + PAF_COLS
BUCKET_PREFIX = "gs://vwb-aou-datasets-controlled/"
DEFAULT_PLATFORMS = ["revio", "sequel2e", "sequel2"]
CHUNK_SIZE = 500
MAX_WORKERS = 16


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def strip_bucket(uri):
    if not isinstance(uri, str) or not uri:
        return None
    if "fc-aou-datasets-controlled" in uri:
        return None
    if uri.startswith(BUCKET_PREFIX):
        return uri[len(BUCKET_PREFIX):]
    if uri.startswith("gs://"):
        return None
    return uri


def resolves(mount, uri):
    rel = strip_bucket(uri)
    return rel is not None and os.path.exists(os.path.join(mount, rel))


def check_row(mount, row):
    """Returns {col: bool} for every column in ALL_CHECK_COLS, for one manifest row."""
    return {c: resolves(mount, row.get(c)) for c in ALL_CHECK_COLS}


def row_trim_tier(flags):
    """Best tier a single row alone would support, or None if it doesn't even have the FASTA."""
    has_fa = all(flags[c] for c in FA_COLS)
    if not has_fa:
        return None
    if all(flags[c] for c in PAF_COLS):
        return "paf_region"
    if all(flags[c] for c in ALN_COLS):
        return "bam_whole_contig"
    return "self_align_needed"


TIER_RANK = {"paf_region": 0, "bam_whole_contig": 1, "self_align_needed": 2}


def best_tier(tiers):
    real = [t for t in tiers if t is not None]
    if not real:
        return None
    return min(real, key=lambda t: TIER_RANK[t])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-n", type=int, default=None,
                    help="Stop after N eligible people found (default: none -- enumerate the "
                         "entire eligible population). Only useful for a quick smoke-test cohort; "
                         "omit for the real production cohort build.")
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    ap.add_argument("--out", default=os.path.expanduser("~/pipeline_outputs/immuannot_cohort_full.tsv"))
    ap.add_argument("--platforms", default=",".join(DEFAULT_PLATFORMS),
                    help=f"Comma-separated platform labels to consider (default "
                         f"{','.join(DEFAULT_PLATFORMS)}). ont-r10.4.1/ont-r9.4.1 have no assembly "
                         f"data at all (reports/lr_data_census/README.md) and are never viable "
                         f"for Immuannot regardless of this flag.")
    ap.add_argument("--force", action="store_true", help="Overwrite an existing --out.")
    ap.add_argument("--restart", action="store_true",
                    help="Ignore any existing checkpoint and start over from row 0.")
    args = ap.parse_args()

    if os.path.exists(args.out) and not args.force:
        die(f"{args.out} already exists. Pass --force to overwrite.")

    lr_path = os.path.join(args.mount, LR_MANIFEST)
    if not os.path.exists(lr_path):
        die(f"manifest not found: {lr_path} -- is the gcsfuse mount up? "
            f"(remount and `ls`-verify before running this.)")

    print(f"Reading LR manifest: {lr_path}", file=sys.stderr)
    lr = pd.read_csv(lr_path, sep="\t", dtype=str)
    for c in ["research_id", "platform"] + ALL_CHECK_COLS:
        if c not in lr.columns:
            die(f"expected column '{c}' not found. Actual columns: {list(lr.columns)}")

    platforms = [p.strip() for p in args.platforms.split(",") if p.strip()]
    candidates = lr[lr["platform"].isin(platforms)].reset_index(drop=True)
    print(f"{len(candidates)} rows across platforms {platforms} out of {len(lr)} total manifest "
          f"rows / {candidates['research_id'].nunique()} unique people -- existence-verifying "
          f"each row against the mount (label alone is not proof, per ENVIRONMENT.md quirk #13).",
          file=sys.stderr)

    out_dir = os.path.dirname(args.out) or "."
    checkpoint_path = os.path.join(out_dir, ".build_immuannot_cohort_full.checkpoint.tsv")
    os.makedirs(out_dir, exist_ok=True)

    done_rows = 0
    if os.path.exists(checkpoint_path) and not args.restart:
        with open(checkpoint_path) as f:
            done_rows = max(0, sum(1 for _ in f) - 1)  # minus header
        if done_rows:
            print(f"Resuming: {done_rows}/{len(candidates)} rows already checkpointed in "
                  f"{checkpoint_path} -- skipping those.", file=sys.stderr)
    elif args.restart and os.path.exists(checkpoint_path):
        os.remove(checkpoint_path)
        print(f"--restart: cleared existing checkpoint {checkpoint_path}", file=sys.stderr)

    header_cols = ["row_idx", "research_id", "platform"] + ALL_CHECK_COLS + ["row_trim_tier"]
    if done_rows == 0:
        with open(checkpoint_path, "w") as f:
            f.write("\t".join(header_cols) + "\n")

    print(f"Checking {len(candidates) - done_rows} remaining rows x {len(ALL_CHECK_COLS)} columns "
          f"({MAX_WORKERS} parallel workers, I/O-bound over the mount, not CPU-bound), "
          f"checkpointed every {CHUNK_SIZE} rows so a restart resumes instead of starting over.",
          file=sys.stderr)

    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        for chunk_start in range(done_rows, len(candidates), CHUNK_SIZE):
            chunk_end = min(chunk_start + CHUNK_SIZE, len(candidates))
            chunk = candidates.iloc[chunk_start:chunk_end]

            futures = {pool.submit(check_row, args.mount, row): local_i
                       for local_i, (_, row) in enumerate(chunk.iterrows())}
            results = {}
            for fut in concurrent.futures.as_completed(futures):
                results[futures[fut]] = fut.result()

            lines = []
            for local_i, (_, row) in enumerate(chunk.iterrows()):
                flags = results[local_i]
                tier = row_trim_tier(flags) or "NONE"
                vals = [str(chunk_start + local_i), str(row["research_id"]), str(row["platform"])] + \
                       [str(int(flags[c])) for c in ALL_CHECK_COLS] + [tier]
                lines.append("\t".join(vals))

            with open(checkpoint_path, "a") as f:
                f.write("\n".join(lines) + "\n")
                f.flush()
                os.fsync(f.fileno())

            print(f"  ...checkpointed rows {chunk_start}-{chunk_end}/{len(candidates)}", file=sys.stderr)

    # --- Per-person aggregation: union across rows, best tier wins (quirk #13 discipline) ---
    ckpt = pd.read_csv(checkpoint_path, sep="\t", dtype=str)
    ckpt["row_trim_tier"] = ckpt["row_trim_tier"].replace("NONE", None)

    print("\n=== Per-ROW trim-tier distribution (before per-person dedup) ===", file=sys.stderr)
    for tier, n in Counter(ckpt["row_trim_tier"].fillna("NONE (no FASTA)")).most_common():
        print(f"  {n:6d}  {tier}", file=sys.stderr)

    per_person = ckpt.groupby("research_id").agg(
        platform=("platform", "first"),
        trim_tier=("row_trim_tier", lambda s: best_tier(list(s))),
        n_rows=("row_trim_tier", "size"),
    ).reset_index()
    per_person = per_person.rename(columns={"research_id": "person_id"})

    eligible = per_person[per_person["trim_tier"].notna()].copy()
    if args.n is not None and len(eligible) > args.n:
        eligible = eligible.head(args.n)
        print(f"\n-n {args.n}: capping to the first {args.n} eligible people found "
              f"(smoke-test mode, not the full cohort).", file=sys.stderr)

    print(f"\n=== Per-PERSON eligibility (n={len(per_person)} unique people across "
          f"{platforms}) ===", file=sys.stderr)
    print(f"  Eligible (has assembly FASTA on >=1 row): {len(eligible)}", file=sys.stderr)
    print(f"  NOT eligible (no FASTA on any row): {len(per_person) - len(eligible)}", file=sys.stderr)
    print("\n=== Per-PERSON trim-tier distribution (eligible only) ===", file=sys.stderr)
    for tier, n in Counter(eligible["trim_tier"]).most_common():
        pct = 100 * n / max(len(eligible), 1)
        print(f"  {n:6d} ({pct:5.1f}%)  {tier}", file=sys.stderr)

    n_needs_fallback = (eligible["trim_tier"] == "self_align_needed").sum()
    if n_needs_fallback:
        print(f"\n  {n_needs_fallback} people need Tier 3 (self-align, UNTESTED AT SCALE as of "
              f"2026-08-04) to be processed -- this is the sequel2 gap. They are INCLUDED in "
              f"{args.out} with trim_tier=self_align_needed so the orchestrator can decide "
              f"whether to launch them (e.g. --skip-trim-tier self_align_needed to hold them back "
              f"for a later batch) -- this script does not make that call.", file=sys.stderr)

    eligible[["person_id", "platform", "trim_tier", "n_rows"]].to_csv(args.out, sep="\t", index=False)
    print(f"\nWrote {len(eligible)} eligible person_ids to {args.out}", file=sys.stderr)
    print(f"Full per-row checkpoint (all columns' individual flags) at {checkpoint_path}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
