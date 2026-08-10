#!/usr/bin/env python3
"""One-time repair: rebuild immuannot_calls.tsv from the raw per-person hap{1,2}.gtf.gz files,
correctly this time.

## Why this exists

merge_fragments() in run_production_orchestrator.py deduplicated immuannot_calls.tsv on
`person_id` alone (fixed 2026-08-10 to `[person_id, gene]` -- see that file's comment at the fix
site). Since that file has one row PER GENE per person, the old dedup silently collapsed every
person down to a single surviving row each time it ran (every 500 completions) -- specifically
whichever gene sorted alphabetically last, which for every person turned out to be one of
MICA/MICB/TAP1/TAP2 (no "HLA-" prefix, sorts after every "HLA-*" classical gene name). The
per-worker fragment files this could have been re-merged from were already deleted by the time
this was found (merge_fragments() removes them once folded into the canonical file) -- so this
script goes one level further back, to the raw hap{1,2}.gtf.gz files themselves, which
RESULTS_LOCATION.md confirms were deliberately kept, not pruned.

## What it does

For every person directory under --outroot with a hap1.gtf.gz and/or hap2.gtf.gz, re-parses both
with the EXACT same regex logic as run_immuannot_person.py's parse_gtf() (copied, not imported --
same "runnable standalone" convention as every other script in this project), rebuilds the
(person_id, gene, immuannot_1, immuannot_2) rows the same way process_person() originally did, and
writes a fresh immuannot_calls.tsv -- one full rewrite, not an incremental patch, so there's no way
for a partial/interrupted repair to leave the file in a worse state than before (the old, corrupted
file is renamed aside first, not deleted, so nothing is lost even if this script itself has a bug).

Does NOT touch immuannot_timing.tsv -- that file's dedup key was already correct
(`[person_id, hap]`), it was never affected by this bug.

Usage (via `pixi run -e spechla`; this is I/O-bound gzip parsing on the LOCAL disk, not the
gcsfuse-mounted bucket -- the resized cheap VM is fine for this, a bigger VM is unlikely to help
much beyond a handful of threads):
  python3 scripts/production_orchestrator/rebuild_immuannot_calls.py

**Before running the real repair, measure the real rate first** -- don't trust a guessed estimate
for an operation touching your whole production cohort:
  python3 scripts/production_orchestrator/rebuild_immuannot_calls.py --limit 200
This processes only the first 200 people, writes to a SEPARATE timing-test file (never touches the
real immuannot_calls.tsv), and prints a measured people/second rate plus an extrapolated total time
for the full cohort. Only run the real repair (no --limit) once that extrapolation looks
reasonable.
"""
import argparse
import concurrent.futures
import gzip
import os
import re
import shutil
import sys
import time

import pandas as pd

DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")


def parse_gtf(gtf_gz_path):
    """Copied verbatim (regex-for-regex) from run_immuannot_person.py's parse_gtf() -- must stay
    byte-identical to how the original calls were derived, or the rebuild wouldn't actually match
    what the pipeline originally produced."""
    calls = {}
    with gzip.open(gtf_gz_path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = fields[8]
            gene_m = re.search(r'gene_name "([^"]+)"', attrs) or re.search(r'gene_id "([^"]+)"', attrs)
            allele_m = re.search(r'consensus "([^"]+)"', attrs) or re.search(r'allele "([^"]+)"', attrs)
            if gene_m and allele_m:
                calls[gene_m.group(1)] = allele_m.group(1)
    return calls


def rebuild_person(pid, person_dir):
    """Mirrors run_immuannot_person.py's process_person() gene_rows construction exactly."""
    gene_calls = {}
    for hap in ("hap1", "hap2"):
        gz = os.path.join(person_dir, f"{hap}.gtf.gz")
        if not os.path.exists(gz):
            continue
        try:
            gene_calls[hap] = parse_gtf(gz)
        except (OSError, EOFError, gzip.BadGzipFile) as e:
            print(f"  WARNING: {pid}/{hap}.gtf.gz unreadable ({e}) -- skipping this haplotype", file=sys.stderr)

    if not gene_calls:
        return []

    genes = sorted(set(gene_calls.get("hap1", {})) | set(gene_calls.get("hap2", {})))
    return [{
        "person_id": pid, "gene": gene,
        "immuannot_1": gene_calls.get("hap1", {}).get(gene, "NA"),
        "immuannot_2": gene_calls.get("hap2", {}).get(gene, "NA"),
    } for gene in genes]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT)
    ap.add_argument("--out", default=None, help="Defaults to <outroot>/immuannot_calls.tsv")
    ap.add_argument("--threads", type=int, default=8,
                    help="Parallel workers -- this is I/O-bound gzip reading on local disk, not "
                         "CPU-bound, so threading helps despite the GIL (same reasoning as this "
                         "project's other FUSE/disk-I/O-bound scripts, e.g. lr_manifest_format_"
                         "census.py). Default 8; raise if the VM has more cores and step 1's timing "
                         "test shows it helps.")
    ap.add_argument("--limit", type=int, default=None,
                    help="TIMING TEST MODE: only process the first N people, write to a separate "
                         "file (never touches the real immuannot_calls.tsv), print a measured "
                         "rate and an extrapolated full-cohort time estimate, then exit. Run this "
                         "BEFORE the real repair (no --limit) to know how long the real one will "
                         "take, instead of guessing.")
    ap.add_argument("--checkpoint-every", type=int, default=200,
                    help="Print progress every N people (ENVIRONMENT.md's checkpoint/progress-"
                         "visibility discipline for any job over a few minutes).")
    args = ap.parse_args()

    person_dirs = sorted(
        d for d in os.listdir(args.outroot)
        if os.path.isdir(os.path.join(args.outroot, d, "immuannot_output"))
    )
    total_people = len(person_dirs)
    print(f"Found {total_people} person directories with an immuannot_output/ subfolder under "
          f"{args.outroot}.", file=sys.stderr)
    if not person_dirs:
        sys.exit(f"FATAL: no person directories found under {args.outroot} -- check --outroot.")

    is_timing_test = args.limit is not None
    if is_timing_test:
        person_dirs = person_dirs[:args.limit]
        out_path = os.path.join(args.outroot, "immuannot_calls.timing_test.tsv")
        print(f"--limit {args.limit}: TIMING TEST MODE, processing {len(person_dirs)} people, "
              f"writing to {out_path} (the real immuannot_calls.tsv is NOT touched).", file=sys.stderr)
    else:
        out_path = args.out or os.path.join(args.outroot, "immuannot_calls.tsv")
        if os.path.exists(out_path):
            backup_path = out_path + f".corrupted-backup-{time.strftime('%Y%m%d-%H%M%S')}"
            shutil.move(out_path, backup_path)
            print(f"Existing (corrupted) file moved aside, not deleted: {backup_path}", file=sys.stderr)

    all_rows = []
    n_no_gtf = 0
    t0 = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {
            pool.submit(rebuild_person, pid, os.path.join(args.outroot, pid, "immuannot_output")): pid
            for pid in person_dirs
        }
        for i, fut in enumerate(concurrent.futures.as_completed(futures), 1):
            rows = fut.result()
            if not rows:
                n_no_gtf += 1
            all_rows.extend(rows)
            if i % args.checkpoint_every == 0 or i == len(person_dirs):
                elapsed = time.time() - t0
                rate = i / elapsed if elapsed > 0 else 0
                print(f"  [{i}/{len(person_dirs)}] {len(all_rows)} gene rows so far, "
                      f"{n_no_gtf} people with no readable GTF, {elapsed:.0f}s elapsed "
                      f"({rate:.1f} people/sec)", file=sys.stderr)

    if not all_rows:
        sys.exit("FATAL: rebuilt 0 gene rows across every person directory -- something is "
                 "fundamentally wrong (wrong --outroot? GTFs actually missing, not just "
                 "mis-merged?). Not writing an empty file over any backup.")

    elapsed = time.time() - t0
    df = pd.DataFrame(all_rows)
    df.to_csv(out_path, sep="\t", index=False)

    if is_timing_test:
        rate = len(person_dirs) / elapsed if elapsed > 0 else 0
        est_full_seconds = total_people / rate if rate > 0 else float("inf")
        print(f"\n=== Timing test done: {len(person_dirs)} people in {elapsed:.0f}s "
              f"({rate:.1f} people/sec) ===", file=sys.stderr)
        print(f"Extrapolated to the full cohort ({total_people} people): "
              f"~{est_full_seconds:.0f}s (~{est_full_seconds / 60:.1f} min).", file=sys.stderr)
        print(f"Delete the test file when done looking at it: rm {out_path}", file=sys.stderr)
        print(f"If that estimate looks acceptable, rerun WITHOUT --limit for the real repair.",
              file=sys.stderr)
        return

    genes_seen = sorted(df["gene"].unique())
    classical_genes_present = [g for g in genes_seen if g.startswith("HLA-")]
    print(f"\nWrote {len(df)} rows ({df['person_id'].nunique()} people) to {out_path} in "
          f"{elapsed:.0f}s.", file=sys.stderr)
    print(f"{len(genes_seen)} distinct genes present: {genes_seen}", file=sys.stderr)
    if not classical_genes_present:
        print("\nWARNING: still 0 'HLA-'-prefixed genes present after rebuild -- the classical "
              "gene calls may genuinely not be in the raw GTFs either (a deeper problem than the "
              "merge bug). Do not assume this repair fixed things without checking the gene list "
              "above.", file=sys.stderr)
    else:
        print(f"\n{len(classical_genes_present)} classical HLA genes now present -- looks like "
              f"the repair worked. Spot-check a few rows before trusting fully: "
              f"head {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
