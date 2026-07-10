#!/usr/bin/env python3
"""Helpers for the SpecHLA short-read padding sweep (Experiment B, EXPERIMENTS.md
2026-07-10). Two subcommands:
  insert-stats     -- reads `samtools stats` output from stdin, prints real
                       insert-size percentiles (median/p95/p99/max) computed
                       from the IS histogram lines. No guessing -- same rigor
                       as the LR sweep's measured read-length floor.
  windows --pad N  -- prints the padded, merged multi-region samtools window
                       string for the 4 classical-HLA gene clusters (GRCh38,
                       Ensembl-derived, same coordinates as the LR aligner x
                       padding sweep). Merges overlapping/touching windows
                       generically for any pad value, so coarse pad levels
                       that merge clusters (like the LR sweep's pad200k did
                       for DRB1/DQA1/DQB1 + DPA1/DPB1) are handled automatically.
"""
import argparse
import sys

GENE_CLUSTERS = [
    (29941260, 29949606),   # A
    (31268746, 31367067),   # C, B
    (32577902, 32668383),   # DRB1, DQA1, DQB1
    (33064569, 33091655),   # DPA1, DPB1
]


def cmd_insert_stats(args):
    total = 0
    counts = {}
    for line in sys.stdin:
        if not line.startswith("IS\t"):
            continue
        parts = line.rstrip("\n").split("\t")
        size, pairs = int(parts[1]), int(parts[2])
        counts[size] = counts.get(size, 0) + pairs
        total += pairs
    if total == 0:
        print("ERROR: no IS (insert size) lines found -- is this `samtools stats` "
              "output on a paired-end BAM?", file=sys.stderr)
        sys.exit(1)
    sizes = sorted(counts)
    cum, pct, targets = 0, {}, [50, 95, 99]
    for s in sizes:
        cum += counts[s]
        frac = cum / total * 100
        for t in list(targets):
            if frac >= t:
                pct[t] = s
                targets.remove(t)
    print(f"n_pairs={total}")
    print(f"median={pct.get(50, 'NA')}")
    print(f"p95={pct.get(95, 'NA')}")
    print(f"p99={pct.get(99, 'NA')}")
    print(f"max={sizes[-1]}")


def merge_windows(windows):
    windows = sorted(windows)
    merged = [windows[0]]
    for s, e in windows[1:]:
        ls, le = merged[-1]
        if s <= le + 1:
            merged[-1] = (ls, max(le, e))
        else:
            merged.append((s, e))
    return merged


def cmd_windows(args):
    pad = args.pad
    raw = [(max(1, s - pad), e + pad) for s, e in GENE_CLUSTERS]
    merged = merge_windows(raw)
    print(" ".join(f"chr6:{s}-{e}" for s, e in merged))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("insert-stats").set_defaults(func=cmd_insert_stats)
    w = sub.add_parser("windows")
    w.add_argument("--pad", type=int, required=True)
    w.set_defaults(func=cmd_windows)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
