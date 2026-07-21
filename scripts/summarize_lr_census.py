#!/usr/bin/env python3
"""Re-derive the lr_manifest_format_census.py distribution summary from its saved output,
without re-running the slow existence-check pass.

Motivation (2026-07-21, Marc): the full census run's per-person profile TSV
(~/pipeline_outputs/lr_data_census/census.tsv) survives a VM restart (ENVIRONMENT.md quirk #14 --
~/pipeline_outputs does), but the printed distribution tables (profile counts, per-platform/center
breakdown) only went to the terminal and are gone once that session's scrollback is gone. This
script rebuilds the same summary from the saved TSV -- instant, no mount needed for the base
distribution. Pass --with-manifest to also rejoin the LR manifest (needs the gcsfuse mount) for
the platform/center breakdown that isn't stored in census.tsv itself.

Usage:
  python3 summarize_lr_census.py [--census ~/pipeline_outputs/lr_data_census/census.tsv]
  python3 summarize_lr_census.py --with-manifest [--mount ~/mnt/aou-controlled]
"""
import argparse
import os
import sys
from collections import Counter

import pandas as pd

LR_MANIFEST = "v9/wgs/long_read/manifest.tsv"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--census", default=os.path.expanduser("~/pipeline_outputs/lr_data_census/census.tsv"))
    ap.add_argument("--with-manifest", action="store_true",
                    help="Rejoin the LR manifest for a platform/center breakdown (needs the mount).")
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    args = ap.parse_args()

    if not os.path.exists(args.census):
        print(f"FATAL: {args.census} not found -- did lr_manifest_format_census.py's full run "
              f"(not --sample-only) actually finish and write it? Check "
              f"~/pipeline_outputs/lr_data_census/ directly.", file=sys.stderr)
        sys.exit(1)

    census = pd.read_csv(args.census, sep="\t", dtype=str)
    census["profile_tuple"] = census["profile"].apply(lambda s: tuple(s.split("+")))

    n = len(census)
    print(f"=== Per-PERSON profile distribution (from {args.census}, n={n} unique people) ===")
    for profile, cnt in Counter(census["profile_tuple"]).most_common():
        pct = 100 * cnt / n
        print(f"  {cnt:6d} ({pct:5.1f}%)  {' + '.join(profile)}")

    n_none = (census["profile"] == "NONE_FOUND").sum()
    print(f"\n  People with NONE of the detected file types resolving: {n_none} ({100*n_none/n:.1f}%)")

    if not args.with_manifest:
        print("\n(Run with --with-manifest for the platform/center breakdown -- needs the "
              "gcsfuse mount, since that's not stored in census.tsv itself.)")
        return

    lr_path = os.path.join(args.mount, LR_MANIFEST)
    if not os.path.exists(lr_path):
        print(f"\nFATAL: {lr_path} not found -- is the gcsfuse mount up? (quirk #11/#14: remount "
              f"and `ls`-verify first.)", file=sys.stderr)
        sys.exit(1)

    lr = pd.read_csv(lr_path, sep="\t", dtype=str, usecols=lambda c: c in ("research_id", "center", "platform"))
    lr = lr.drop_duplicates("research_id")  # one platform/center label per person for this breakdown
    merged = census.merge(lr, left_on="research_id", right_on="research_id", how="left")

    print(f"\n=== Per-PERSON profile x platform (top combos) ===")
    combo_counts = Counter(zip(merged["platform"], merged["profile_tuple"]))
    for (platform, profile), cnt in sorted(combo_counts.items(), key=lambda kv: -kv[1])[:40]:
        print(f"  {cnt:6d}  platform={str(platform):15s} profile={' + '.join(profile)}")

    print(f"\n=== Per-PERSON profile x center (top combos) ===")
    combo_counts = Counter(zip(merged["center"], merged["profile_tuple"]))
    for (center, profile), cnt in sorted(combo_counts.items(), key=lambda kv: -kv[1])[:40]:
        print(f"  {cnt:6d}  center={str(center):8s} profile={' + '.join(profile)}")

    print(f"\n=== NONE_FOUND people by platform (the ones worth escalating) ===")
    none_rows = merged[merged["profile"] == "NONE_FOUND"]
    for platform, cnt in Counter(none_rows["platform"]).most_common():
        print(f"  {cnt:6d}  platform={platform}")


if __name__ == "__main__":
    main()
