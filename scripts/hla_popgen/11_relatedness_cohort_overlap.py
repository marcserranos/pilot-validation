#!/usr/bin/env python3
"""Cross-reference AoU's own genome-wide relatedness table against the hla_popgen cohorts, to scope
the phasing-validation study proposed in DECISIONS.md / the phasing research thread: can we validate
long-read (Immuannot) haplotype phasing at DQA1~DQB1 / DPA1~DPB1 using REAL relative pairs already
present in the biobank, instead of needing purpose-sequenced trios?

AoU computes a KING-style kinship coefficient for every pair of short-read-WGS participants whose
share exceeds ~0.1, in `v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv`
(columns: i.s, j.s, kin -- confirmed empirically 2026-09-08, NOT documented with exact column names
in AoU's public docs). This script buckets every relative pair by which of our own two HLA cohorts
(SCHEMA.md Table 4, `cohort_membership.tsv`) each side falls into:

  - both_lr:  BOTH people already have a phased long-read (Immuannot) assembly. This is the easy,
    computationally-light validation case per Marc's 2026-09-08 framing -- no alignment pipeline
    needed, just a direct sequence comparison between the child's and parent's already-assembled
    haplotypes at the region of interest. A parent-child pair in this bucket can even reveal a real
    meiotic crossover breakpoint, not just validate absence of an assembly switch-error.
  - lr_sr:    One side has a long-read assembly, the other only has AoU-native short-read data. This
    is the harder case -- validating requires actually aligning the short-read side's reads to the
    long-read side's phased assembly contigs.
  - other:    Neither side touches our lr cohort -- not useful for this study, counted only for
    context (what fraction of AoU relative pairs we simply can't use here).

**Aggregate-only discipline (SCHEMA.md hard rule 5):** this script prints/writes ONLY counts and
distributions to reports/. The person_id-level pair list (needed for the next phase -- actually
pulling specific people's data) is written to ~/pipeline_outputs/, never to reports/, and contains
no genotypes, only person_id pairs + which bucket + kinship value.

Kinship-degree binning follows the standard KING thresholds (Manichaikul et al. 2010): >0.354
duplicate/MZ twin, 0.177-0.354 first-degree (parent-child OR full sibling -- indistinguishable from
kinship alone, see module note below), 0.0884-0.177 second-degree, 0.0442-0.0884 third-degree.
**Important limitation, confirmed empirically:** the real file has only 3 columns (i.s, j.s, kin) --
no IBD0 sharing fraction, which is what would normally separate parent-child from full-sibling pairs
within the first-degree band. Both relationship types are reported together as "first_degree" here;
distinguishing them (if it matters for a later analysis) needs a different data source entirely.

Usage (fixtures):
    python3 scripts/hla_popgen/11_relatedness_cohort_overlap.py \\
        --relatedness-table /tmp/hla_fixtures/samples_relatedness.sample.tsv \\
        --cohort-membership /tmp/hla_fixtures/cohort_membership.sample.tsv \\
        --outroot /tmp/hla_fixtures --out-dir /tmp/hla_fixtures/reports

Real run (VM):
    python3 scripts/hla_popgen/11_relatedness_cohort_overlap.py
"""
import argparse
import os
import sys

import pandas as pd

DEFAULT_RELATEDNESS = os.path.expanduser(
    "~/mnt/aou-controlled/v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv")
DEFAULT_COHORT_MEMBERSHIP = os.path.expanduser("~/pipeline_outputs/cohort_membership.tsv")
DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_OUT_DIR = "reports/hla_popgen/11_relatedness_cohort_overlap"

KIN_BINS = [
    (0.354, float("inf"), "duplicate_or_MZ_twin"),
    (0.177, 0.354, "first_degree_parentchild_or_sibling"),
    (0.0884, 0.177, "second_degree"),
    (0.0442, 0.0884, "third_degree"),
    (-float("inf"), 0.0442, "below_third_degree"),
]


def kin_degree(k):
    for lo, hi, label in KIN_BINS:
        if lo < k <= hi:
            return label
    return "unknown"


def load_relatedness(path):
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Relatedness table not found at {path}. Is the controlled-tier bucket mounted? "
            "See RUNBOOK.md Step 3 (gcsfuse mount) before running this script.")
    df = pd.read_csv(path, sep="\t", dtype={"i.s": str, "j.s": str})
    missing = {"i.s", "j.s", "kin"} - set(df.columns)
    if missing:
        raise ValueError(f"Relatedness table missing expected columns {missing} -- schema changed?")
    return df.rename(columns={"i.s": "person_i", "j.s": "person_j"})


def load_cohort_membership(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"cohort_membership.tsv not found at {path} -- run 02_build_cohorts.py first.")
    df = pd.read_csv(path, sep="\t", dtype={"person_id": str})
    for col in ("in_sr", "in_lr"):
        if df[col].dtype != bool:
            df[col] = df[col].astype(str).str.strip().str.lower().isin({"true", "1"})
    return df


def classify_pairs(rel, cohort):
    lr_ids = set(cohort.loc[cohort["in_lr"], "person_id"])
    sr_ids = set(cohort.loc[cohort["in_sr"], "person_id"])

    def bucket(row):
        i, j = row["person_i"], row["person_j"]
        i_lr, j_lr = i in lr_ids, j in lr_ids
        i_sr, j_sr = i in sr_ids, j in sr_ids
        if i_lr and j_lr:
            return "both_lr"
        if (i_lr and j_sr and not j_lr) or (j_lr and i_sr and not i_lr):
            return "lr_sr"
        if i_lr or j_lr:
            # one side is lr, but the other has neither an sr nor lr flag in cohort_membership --
            # i.e. a real AoU relative who simply isn't in either HLA-called cohort. Still worth
            # counting separately: it's a relative we can't use for THIS study, but a signal of how
            # complete cohort_membership's population coverage actually is.
            return "lr_but_relative_uncalled"
        return "other"

    rel = rel.copy()
    rel["bucket"] = rel.apply(bucket, axis=1)
    rel["kin_degree"] = rel["kin"].astype(float).map(kin_degree)
    return rel


def summarize(rel, lr_ids):
    counts = rel["bucket"].value_counts().rename_axis("bucket").reset_index(name="n_pairs")
    useful = rel[rel["bucket"].isin(["both_lr", "lr_sr"])]
    degree_by_bucket = (useful.groupby(["bucket", "kin_degree"]).size()
                        .rename("n_pairs").reset_index())
    # Count distinct LONG-READ-COHORT people who have >=1 usable relative, per bucket. For "both_lr"
    # both sides qualify (both are lr); for "lr_sr" only the side that is actually IN lr_ids counts --
    # the other side is the short-read-only relative and must never be counted as "lr coverage"
    # (a real bug caught by the fixture test: unioning both columns blindly inflated coverage past
    # 100% by counting sr-only relatives as if they were lr-cohort people).
    lr_people_covered = {}
    for b in ("both_lr", "lr_sr"):
        sub = rel[rel["bucket"] == b]
        ids = set(sub["person_i"]).union(sub["person_j"])
        lr_people_covered[b] = ids & lr_ids
    return counts, degree_by_bucket, lr_people_covered


def _df_to_markdown(df):
    """Plain-Python markdown table writer -- avoids an optional `tabulate` dependency that isn't
    guaranteed present in every pixi env this script runs in (confirmed missing in `spechla` env,
    2026-09-08)."""
    cols = list(df.columns)
    header = "| " + " | ".join(cols) + " |"
    sep = "|" + "|".join(["---"] * len(cols)) + "|"
    rows = ["| " + " | ".join(str(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join([header, sep] + rows)


def write_report(counts, degree_by_bucket, lr_people_covered, cohort, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    n_lr_total = int(cohort["in_lr"].sum())
    both_lr_people = lr_people_covered["both_lr"]
    lr_sr_people = lr_people_covered["lr_sr"]
    any_usable = both_lr_people | lr_sr_people

    lines = ["# Relatedness x long-read cohort overlap -- phasing-validation scoping\n"]
    lines.append(f"Long-read (phased) cohort size: **{n_lr_total}**\n")
    lines.append("## Pair counts by bucket\n")
    lines.append(_df_to_markdown(counts))
    lines.append("\n\n## Kinship-degree distribution, usable buckets only\n")
    lines.append(_df_to_markdown(degree_by_bucket))
    lines.append("\n\n## Coverage of the long-read cohort\n")
    lines.append(f"- People with >=1 relative who ALSO has a long-read assembly (`both_lr`, "
                 f"the easy/direct-comparison case): **{len(both_lr_people)}**\n")
    lines.append(f"- People with >=1 relative who only has short-read data (`lr_sr`, needs "
                 f"read-realignment): **{len(lr_sr_people)}**\n")
    lines.append(f"- People with >=1 usable relative of EITHER kind: **{len(any_usable)}** "
                 f"out of {n_lr_total} long-read-cohort people "
                 f"({100 * len(any_usable) / max(n_lr_total, 1):.1f}%)\n")
    lines.append("\n**Known limitation:** the relatedness table has no IBD0 column, so "
                 "`first_degree_parentchild_or_sibling` pairs cannot be split into parent-child vs. "
                 "full-sibling from this data alone.\n")
    with open(os.path.join(out_dir, "relatedness_overlap_report.md"), "w") as f:
        f.write("\n".join(lines))
    counts.to_csv(os.path.join(out_dir, "bucket_counts.tsv"), sep="\t", index=False)
    degree_by_bucket.to_csv(os.path.join(out_dir, "kin_degree_by_bucket.tsv"), sep="\t", index=False)


def write_pair_list(rel, outroot):
    """Person-id-level output -- stays in pipeline_outputs, never in reports/ (SCHEMA.md hard rule 5)."""
    useful = rel[rel["bucket"].isin(["both_lr", "lr_sr"])][
        ["person_i", "person_j", "kin", "kin_degree", "bucket"]]
    path = os.path.join(outroot, "relatedness_lr_overlap_pairs.tsv")
    useful.to_csv(path, sep="\t", index=False)
    return path, len(useful)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--relatedness-table", default=DEFAULT_RELATEDNESS)
    ap.add_argument("--cohort-membership", default=DEFAULT_COHORT_MEMBERSHIP)
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT,
                    help="Where to write the person-id-level pair list (pipeline_outputs, not reports/).")
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    args = ap.parse_args()

    print("Loading relatedness table...", file=sys.stderr)
    rel = load_relatedness(args.relatedness_table)
    print(f"  {len(rel)} relative pairs loaded.", file=sys.stderr)

    print("Loading cohort membership...", file=sys.stderr)
    cohort = load_cohort_membership(args.cohort_membership)
    print(f"  {len(cohort)} people, {int(cohort['in_lr'].sum())} in lr cohort, "
          f"{int(cohort['in_sr'].sum())} in sr cohort.", file=sys.stderr)

    print("Classifying pairs...", file=sys.stderr)
    rel = classify_pairs(rel, cohort)

    lr_ids = set(cohort.loc[cohort["in_lr"], "person_id"])
    counts, degree_by_bucket, lr_people_covered = summarize(rel, lr_ids)
    print(counts.to_string(index=False), file=sys.stderr)

    write_report(counts, degree_by_bucket, lr_people_covered, cohort, args.out_dir)
    pair_path, n_pairs = write_pair_list(rel, args.outroot)
    print(f"Wrote aggregate report to {args.out_dir}/relatedness_overlap_report.md", file=sys.stderr)
    print(f"Wrote {n_pairs} usable pairs (person_ids, pipeline_outputs only) to {pair_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
