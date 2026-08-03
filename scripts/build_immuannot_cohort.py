#!/usr/bin/env python3
"""Immuannot-specific cohort builder. Filters the v9 lrWGS manifest to people who actually have
what run_immuannot_person.py needs -- existence-verified, not inferred from a platform label.

WHY THIS EXISTS, SEPARATE FROM build_experiment_d_cohort.py (2026-08-03): that script filters on
a real `/revio/`-pattern short-read-adjacent BAM, which was the right check for the original
SpecHLA/SpecImmune 3-way comparison -- but it is NOT the same requirement Immuannot has, and a
cohort built that way can (and did, in practice: person 1000151) include people who pass that
filter but have no assembly data at all. Immuannot needs, per haplotype: assembly_hap{1,2}_fa
(the actual input) AND assembly_hap{1,2}_aln2_hg38_bam (used by run_immuannot_person.py's trim
step -- both currently required, even though only the .fa is Immuannot's own real dependency; see
DECISIONS.md's open "sequel2 exclusion" note about that gap, not fixed here).

Per reports/lr_data_census/README.md, `platform == revio` always carries the full assembly suite
(11,070/11,070 rows) -- so filtering to revio AND existence-verifying is belt-and-suspenders, not
redundant: the label tells you where to look, existence-checking confirms it's actually there
(same discipline as ENVIRONMENT.md quirk #13, learned the hard way after three rounds of trusting
a label/pattern instead of checking).

Usage (from ~/repos/pilot-validation, inside `pixi shell -e specimmune` or `pixi run -e specimmune --`):
  python3 scripts/build_immuannot_cohort.py [-n 40] [--out PATH] [--force] [--mount ~/mnt/aou-controlled]
"""
import argparse
import os
import sys

import pandas as pd

LR_MANIFEST = "v9/wgs/long_read/manifest.tsv"
FA_COLS = ["assembly_hap1_fa", "assembly_hap2_fa"]
ALN_COLS = ["assembly_hap1_aln2_hg38_bam", "assembly_hap2_aln2_hg38_bam"]
BUCKET_PREFIX = "gs://vwb-aou-datasets-controlled/"


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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-n", type=int, default=40,
                    help="How many verified people to collect (default 40 -- comfortably covers "
                         "a 32-way concurrency test with margin for a few more skips).")
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    ap.add_argument("--out", default=os.path.expanduser("~/pipeline_outputs/immuannot_cohort.tsv"))
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if os.path.exists(args.out) and not args.force:
        die(f"{args.out} already exists. Pass --force to overwrite.")

    lr_path = os.path.join(args.mount, LR_MANIFEST)
    if not os.path.exists(lr_path):
        die(f"manifest not found: {lr_path} -- is the gcsfuse mount up? "
            f"(remount and `ls`-verify before running this.)")

    print(f"Reading LR manifest: {lr_path}", file=sys.stderr)
    lr = pd.read_csv(lr_path, sep="\t", dtype=str)
    for c in ["research_id", "platform"] + FA_COLS + ALN_COLS:
        if c not in lr.columns:
            die(f"expected column '{c}' not found. Actual columns: {list(lr.columns)}")

    revio = lr[lr["platform"] == "revio"].copy()
    print(f"{len(revio)} revio-platform rows out of {len(lr)} total manifest rows -- "
          f"existence-verifying each candidate against the mount (label alone is not proof, "
          f"per ENVIRONMENT.md quirk #13).", file=sys.stderr)

    verified = []
    checked = 0
    for _, row in revio.drop_duplicates("research_id").iterrows():
        checked += 1
        if all(resolves(args.mount, row[c]) for c in FA_COLS + ALN_COLS):
            verified.append(row["research_id"])
        if len(verified) >= args.n:
            break

    print(f"Checked {checked} candidates, verified {len(verified)} with all 4 required files "
          f"(assembly_hap1_fa, assembly_hap2_fa, assembly_hap1_aln2_hg38_bam, "
          f"assembly_hap2_aln2_hg38_bam) actually present on the mount.", file=sys.stderr)

    if len(verified) < args.n:
        print(f"WARNING: only found {len(verified)}/{args.n} requested -- ran out of revio "
              f"candidates in this pass, or an unusual number failed verification. Not fatal, "
              f"but check before assuming this cohort is big enough for what you need it for.",
              file=sys.stderr)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    pd.DataFrame({"person_id": verified}).to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(verified)} verified person_ids to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
