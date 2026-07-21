#!/usr/bin/env python3
"""Exhaustive long-read data-format census (not part of the lettered A-D roadmap -- "Experiment E"
is already taken by the disease-sanity-check thread, see EXPERIMENTS.md/reports/disease_sanity_check/).

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

**Rewritten 2026-07-21 after a real 5-hour data-loss incident.** The original version did all
existence-checking (16 columns x ~15,424 rows, serial, one os.path.exists() per cell over gcsfuse)
before writing anything to disk -- so a VM auto-stop (quirk #14) or any interruption mid-run lost
100% of the work, no matter how close to finishing. Fixed three ways:
  1. **Checkpointed, resumable, chunked processing** -- results are appended to a checkpoint file
     every CHUNK_SIZE rows, flushed immediately. A restart resumes from the last completed chunk
     instead of starting over (mirrors run_experiment_d.sh's per-person resumability).
  2. **Only "primary" (data-bearing) columns are existence-checked, not their index/companion
     files** (.bai/.tbi/.csi/.pbi/.gzi are still detected and reported, just not individually
     verified -- they're >99.9% correlated with their primary file and checking them roughly
     doubled runtime for no new "which product type does this person have" information).
  3. **Parallelized** the existence checks with a thread pool -- these are I/O-bound network
     round-trips through FUSE, not CPU-bound, so threading (not multiprocessing) gives a large
     speedup without fighting the GIL.

Two multi-row-per-person subtleties already learned in Experiment D's cohort builder, reused here:
  - ~879 of 14,521 people have more than one manifest row (different center/platform).
  - A person's other row can resolve when their primary-looking one doesn't (person 1008366).
So this script reports BOTH a per-row and a per-person (union across rows) profile distribution.

Usage (any pixi env with pandas; read-only, no writes to the mount):
  python3 lr_manifest_format_census.py [--mount ~/mnt/aou-controlled] \
      [--out ~/pipeline_outputs/lr_data_census/census.tsv] [--sample-only N] [--restart]

  --sample-only N   Only process the first N manifest rows -- ALWAYS writes to its own
                     census_sample<N>.tsv / .checkpoint_sample<N>.tsv paths (never the full-run
                     paths), so a sanity-check run can never clobber a real result again.
  --restart         Ignore any existing checkpoint and start over from row 0 (default: resume).
"""
import argparse
import concurrent.futures
import os
import sys
from collections import Counter

import pandas as pd

BUCKET_PREFIX = "gs://vwb-aou-datasets-controlled/"
LR_MANIFEST = "v9/wgs/long_read/manifest.tsv"
CHUNK_SIZE = 500
MAX_WORKERS = 16

# Data-bearing extensions -- existence-checked directly (the "primary" file per product).
PRIMARY_EXTENSIONS = (
    ".bam", ".cram", ".fastq", ".fastq.gz", ".fq.gz",
    ".fa", ".fasta", ".fa.gz", ".fasta.gz",
    ".vcf", ".vcf.gz", ".g.vcf.gz", ".gvcf", ".gvcf.gz",
    ".paf", ".gfa", ".snf",
)
# Companion/index extensions -- detected and reported, but NOT existence-checked on their own.
INDEX_EXTENSIONS = (".bai", ".crai", ".tbi", ".csi", ".pbi", ".gzi")
KNOWN_EXTENSIONS = PRIMARY_EXTENSIONS + INDEX_EXTENSIONS


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


def classify_column(series, sample_n=50):
    """Returns 'primary', 'index', or None. Primary columns get existence-checked; index columns
    are reported but skipped (see module docstring point 2)."""
    non_null = series.dropna()
    if non_null.empty:
        return None
    sample = non_null.astype(str).head(sample_n)
    is_uri = lambda v: v.startswith("gs://") or "/" in v
    primary_hits = sum(1 for v in sample if is_uri(v) and v.lower().endswith(PRIMARY_EXTENSIONS))
    index_hits = sum(1 for v in sample if is_uri(v) and v.lower().endswith(INDEX_EXTENSIONS))
    half = max(1, len(sample) // 2)
    if primary_hits >= half:
        return "primary"
    if index_hits >= half:
        return "index"
    return None


def check_exists(mount, uri):
    rel = strip_bucket(uri)
    return bool(rel) and os.path.exists(os.path.join(mount, rel))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    ap.add_argument("--out", default=None,
                    help="Default: ~/pipeline_outputs/lr_data_census/census.tsv (or "
                         "census_sample<N>.tsv under --sample-only).")
    ap.add_argument("--sample-only", type=int, default=None,
                    help="Only process the first N manifest rows. Always uses its own output/"
                         "checkpoint paths -- can never clobber a full run's result.")
    ap.add_argument("--restart", action="store_true",
                    help="Ignore any existing checkpoint and start over from row 0.")
    ap.add_argument("--force", action="store_true",
                    help="Overwrite an existing final --out even if one is already there.")
    args = ap.parse_args()

    out_dir = os.path.expanduser("~/pipeline_outputs/lr_data_census")
    if args.out is None:
        args.out = os.path.join(out_dir, f"census_sample{args.sample_only}.tsv"
                                 if args.sample_only else "census.tsv")
    checkpoint_path = os.path.join(
        out_dir, f".checkpoint_sample{args.sample_only}.tsv" if args.sample_only else ".checkpoint.tsv"
    )

    if os.path.exists(args.out) and not args.force:
        die(f"{args.out} already exists -- refusing to overwrite a (possibly full, hard-won) "
            f"prior result. Pass --force if you're sure, or --sample-only N to run a throwaway "
            f"check against its own separate path instead.")

    lr_path = os.path.join(args.mount, LR_MANIFEST)
    if not os.path.exists(lr_path):
        die(f"manifest not found: {lr_path} -- is the gcsfuse mount up? "
            f"(ENVIRONMENT.md quirk #11/#14: remount and `ls`-verify first.)")

    print(f"Reading LR manifest: {lr_path}", file=sys.stderr)
    lr = pd.read_csv(lr_path, sep="\t", dtype=str)
    if args.sample_only:
        lr = lr.head(args.sample_only)
        print(f"  --sample-only {args.sample_only}: using first {len(lr)} rows only "
              f"(writing to {args.out}, separate from any full-run output)", file=sys.stderr)

    n_rows = len(lr)
    n_people = lr["research_id"].nunique() if "research_id" in lr.columns else None
    print(f"  {n_rows} rows, {len(lr.columns)} columns, "
          f"{n_people if n_people is not None else '?'} unique research_id", file=sys.stderr)

    # --- Column report: full column list + classification ---
    print("\n=== All manifest columns (non-null count, classification) ===", file=sys.stderr)
    primary_cols, index_cols = [], []
    for col in lr.columns:
        cls = classify_column(lr[col])
        if cls == "primary":
            primary_cols.append(col)
        elif cls == "index":
            index_cols.append(col)
        non_null = lr[col].notna().sum()
        sample_val = lr[col].dropna().astype(str).head(1).tolist()
        sample_val = sample_val[0][:90] if sample_val else ""
        flag = {"primary": "PRIMARY", "index": "index (skipped)"}.get(cls, "")
        print(f"  {col:35s} non_null={non_null:6d}  {flag:18s} e.g. {sample_val}", file=sys.stderr)

    if not primary_cols:
        die("No primary path-like columns auto-detected -- PRIMARY_EXTENSIONS may need updating "
            "for this manifest's actual column contents (see the column report above).")

    print(f"\n=== {len(primary_cols)} PRIMARY columns will be existence-checked: {primary_cols} ===",
          file=sys.stderr)
    print(f"=== {len(index_cols)} companion/index columns detected but not checked "
          f"(correlated with their primary file): {index_cols} ===", file=sys.stderr)
    print("Sanity-check this against the column report above before trusting the existence-check "
          "pass below. Re-run with --sample-only 50 first if unsure.", file=sys.stderr)

    # --- Resumable, chunked, parallelized existence check ---
    os.makedirs(out_dir, exist_ok=True)
    done_rows = 0
    if os.path.exists(checkpoint_path) and not args.restart:
        with open(checkpoint_path) as f:
            done_rows = max(0, sum(1 for _ in f) - 1)  # minus header
        if done_rows:
            print(f"\nResuming: {done_rows}/{n_rows} rows already checkpointed in "
                  f"{checkpoint_path} -- skipping those, continuing from row {done_rows}.",
                  file=sys.stderr)
    elif args.restart and os.path.exists(checkpoint_path):
        os.remove(checkpoint_path)
        print(f"\n--restart: cleared existing checkpoint {checkpoint_path}", file=sys.stderr)

    header_cols = ["row_idx", "research_id"] + primary_cols + ["profile"]
    if done_rows == 0:
        with open(checkpoint_path, "w") as f:
            f.write("\t".join(header_cols) + "\n")

    print(f"\nVerifying real file existence for {n_rows - done_rows} remaining rows x "
          f"{len(primary_cols)} primary columns against the mount, in chunks of {CHUNK_SIZE} "
          f"({MAX_WORKERS} parallel workers -- I/O-bound over the mount, not CPU-bound). "
          f"Checkpointed after every chunk, so a restart resumes here instead of from scratch.",
          file=sys.stderr)

    research_id_col = lr["research_id"] if "research_id" in lr.columns else None
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        for chunk_start in range(done_rows, n_rows, CHUNK_SIZE):
            chunk_end = min(chunk_start + CHUNK_SIZE, n_rows)
            chunk = lr.iloc[chunk_start:chunk_end]

            # Submit every (row, primary_col) cell in this chunk to the thread pool at once.
            futures = {}
            for local_i, (idx, row) in enumerate(chunk.iterrows()):
                for col in primary_cols:
                    futures[pool.submit(check_exists, args.mount, row[col])] = (local_i, col)
            results = {}  # local_i -> {col: bool}
            for fut in concurrent.futures.as_completed(futures):
                local_i, col = futures[fut]
                results.setdefault(local_i, {})[col] = fut.result()

            lines = []
            for local_i, (idx, row) in enumerate(chunk.iterrows()):
                flags = results.get(local_i, {})
                present = tuple(sorted(c for c in primary_cols if flags.get(c)))
                profile = present if present else ("NONE_FOUND",)
                rid = row["research_id"] if research_id_col is not None else ""
                vals = [str(chunk_start + local_i), str(rid)] + \
                       [str(int(flags.get(c, False))) for c in primary_cols] + \
                       ["+".join(profile)]
                lines.append("\t".join(vals))

            with open(checkpoint_path, "a") as f:
                f.write("\n".join(lines) + "\n")
                f.flush()
                os.fsync(f.fileno())

            print(f"  ...checkpointed rows {chunk_start}-{chunk_end}/{n_rows}", file=sys.stderr)

    # --- Read the full checkpoint back (covers this run + any prior resumed chunks) ---
    ckpt = pd.read_csv(checkpoint_path, sep="\t", dtype=str)
    ckpt["profile_tuple"] = ckpt["profile"].apply(lambda s: tuple(s.split("+")))

    print("\n=== Per-ROW profile distribution (before per-person dedup) ===", file=sys.stderr)
    for profile, n in Counter(ckpt["profile_tuple"]).most_common():
        print(f"  {n:6d}  {' + '.join(profile)}", file=sys.stderr)

    if "platform" in lr.columns:
        lr_platform = lr.reset_index().rename(columns={"index": "row_idx"})[["row_idx", "platform"]]
        lr_platform["row_idx"] = lr_platform["row_idx"].astype(str)
        joined = ckpt.merge(lr_platform, on="row_idx", how="left")
        print("\n=== Per-ROW profile x platform (top combos) ===", file=sys.stderr)
        combo_counts = Counter(zip(joined["platform"], joined["profile_tuple"]))
        for (platform, profile), n in sorted(combo_counts.items(), key=lambda kv: -kv[1])[:40]:
            print(f"  {n:6d}  platform={str(platform):15s} profile={' + '.join(profile)}", file=sys.stderr)

    if "center" in lr.columns:
        lr_center = lr.reset_index().rename(columns={"index": "row_idx"})[["row_idx", "center"]]
        lr_center["row_idx"] = lr_center["row_idx"].astype(str)
        joined = ckpt.merge(lr_center, on="row_idx", how="left")
        print("\n=== Per-ROW profile x center (top combos) ===", file=sys.stderr)
        combo_counts = Counter(zip(joined["center"], joined["profile_tuple"]))
        for (center, profile), n in sorted(combo_counts.items(), key=lambda kv: -kv[1])[:40]:
            print(f"  {n:6d}  center={str(center):8s} profile={' + '.join(profile)}", file=sys.stderr)

    # --- Per-PERSON profile: union of what's available across all of a person's rows ---
    def union_profile(profiles):
        present = set()
        for p in profiles:
            if p != ("NONE_FOUND",):
                present.update(p)
        return tuple(sorted(present)) if present else ("NONE_FOUND",)

    per_person = ckpt.groupby("research_id")["profile_tuple"].apply(union_profile)
    print(f"\n=== Per-PERSON profile distribution (union across rows, n={len(per_person)} unique people) ===",
          file=sys.stderr)
    for profile, n in Counter(per_person).most_common():
        pct = 100 * n / len(per_person)
        print(f"  {n:6d} ({pct:5.1f}%)  {' + '.join(profile)}", file=sys.stderr)

    n_none = sum(1 for p in per_person if p == ("NONE_FOUND",))
    n_multi_row = (ckpt.groupby("research_id").size() > 1).sum()
    print(f"\n  People with >1 manifest row: {n_multi_row}", file=sys.stderr)
    print(f"  People with NONE of the {len(primary_cols)} detected primary file types resolving: "
          f"{n_none} -- these are the ones worth escalating (different/supplementary manifest? "
          f"ask Aleix per ENVIRONMENT.md quirk #8; check Workbench Data Dictionary directly) "
          f"rather than assumed to genuinely lack data.", file=sys.stderr)

    out_df = per_person.reset_index()
    out_df.columns = ["research_id", "profile"]
    out_df["profile"] = out_df["profile"].apply(lambda p: "+".join(p))
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"\n=== Wrote per-person profile table to {args.out} ===", file=sys.stderr)
    print(f"=== Full per-row checkpoint (all columns' individual flags) at {checkpoint_path} ===",
          file=sys.stderr)


if __name__ == "__main__":
    main()
