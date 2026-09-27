#!/usr/bin/env python3
"""Stage 2 follow-up (2026-07-25): loads AoU's richer RNA-SeQC2 per-sample metrics file
(rnaseqc2/*.metrics.txt.gz -- small, ~2.8MB, safe to load in full, unlike the multi-GB RSEM/
exon-read matrices in the same rnaseq/ tree, which this script deliberately does NOT touch) and:

  1. Checks whether any of its richer columns (Read Length, Fragment GC Content stats,
     non-globin read fraction) explain the mean_insert_size bimodality found in Stage 1
     (rnaseq_metadata_census.py) that nothing in the smaller rnaseq_metadata.tsv explained.
  2. Surfaces globin-contamination and rRNA-rate as candidate per-sample selection criteria for
     the eventual repertoire-extraction pilot (Stage 4) -- high globin/rRNA fraction eats into
     the read budget available for everything else, including TCR/BCR transcripts.
  3. Cross-checks rRNA Rate here against ribosomal_bases_pct from rnaseq_metadata.tsv -- two
     independently-computed metrics that should agree if both pipelines are self-consistent.

Deliberately does NOT load rsem/*.txt.gz or rnaseqc2/*.gct.gz -- those are cohort-wide gene/
transcript expression matrices (exon_reads.gct.gz alone is 5.39 GiB compressed; likely tens of
GB uncompressed) and don't belong in a plain pd.read_csv on a 25GB VM. If/when a specific gene
(e.g. globin genes for a covariate) is actually needed from those, extract it by streaming
(zgrep/awk on the gene_id column) rather than loading the whole matrix -- not implemented here,
not needed yet.

The metrics filename is date-stamped (confirmed live: aou_rnaseq_20260413.metrics.txt.gz) --
resolved by glob, not hardcoded, in case a future CDR update changes the date.

Usage (needs pandas + matplotlib -- the spechla env has both):
  pixi run -e spechla -- python3 scripts/rnaseq_seqc2_metrics_check.py
"""
import argparse
import glob
import os
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RNASEQ_METADATA = "v9/multiomics/rnaseq/rnaseq_metadata.tsv"
SEQC2_METRICS_GLOB = "v9/multiomics/rnaseq/rnaseqc2/*.metrics.txt.gz"

# Exact column names from the real header (verbose/punctuated -- copied exactly, not guessed).
COLS_OF_INTEREST = [
    "Read Length",
    "Fragment GC Content Mean",
    "Fragment GC Content Std",
    "Fragment GC Content Skewness",
    "Fragment GC Content Kurtosis",
    "Duplicate Rate of Mapped, excluding Globins",
    "Non-Globin Reads",
    "Total Reads",
    "rRNA Rate",
    "Genes Detected",
    "Exonic Rate",
]


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    ap.add_argument("--out-dir", default=os.path.expanduser("~/pipeline_outputs/rnaseq_overlap"))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    meta_path = os.path.join(args.mount, RNASEQ_METADATA)
    if not os.path.exists(meta_path):
        die(f"not found: {meta_path} -- is the gcsfuse mount up? Remount + `ls`-verify first.")

    candidates = glob.glob(os.path.join(args.mount, SEQC2_METRICS_GLOB))
    if len(candidates) != 1:
        die(f"expected exactly one *.metrics.txt.gz under {SEQC2_METRICS_GLOB}, found "
            f"{len(candidates)}: {candidates}. Mount up? Filename changed?")
    seqc2_path = candidates[0]

    print(f"Reading Stage 1 metadata (for mean_insert_size, ribosomal_bases_pct): {meta_path}",
          file=sys.stderr)
    meta = pd.read_csv(meta_path, sep="\t", dtype=str)
    meta["mean_insert_size"] = pd.to_numeric(meta["mean_insert_size"], errors="coerce")
    meta["ribosomal_bases_pct"] = pd.to_numeric(meta["ribosomal_bases_pct"], errors="coerce")

    size_mb = os.path.getsize(seqc2_path) / 1e6
    print(f"Reading RNA-SeQC2 metrics ({size_mb:.1f} MB compressed -- small, loading in full): "
          f"{seqc2_path}", file=sys.stderr)
    seqc2 = pd.read_csv(seqc2_path, sep="\t")
    seqc2.columns = [c.strip() for c in seqc2.columns]
    seqc2["sample_id"] = seqc2["sample_id"].astype(str)
    missing = [c for c in COLS_OF_INTEREST if c not in seqc2.columns]
    if missing:
        die(f"expected column(s) not found in {seqc2_path}: {missing}. Actual columns: "
            f"{list(seqc2.columns)}")
    print(f"  {len(seqc2)} rows, {len(seqc2.columns)} columns", file=sys.stderr)

    # candidate Stage 4 sample-selection criterion -- how much of the read budget survives
    # globin filtering, which matters directly for whether enough non-globin depth remains to
    # recover rare TCR/BCR transcripts
    seqc2["non_globin_frac"] = seqc2["Non-Globin Reads"] / seqc2["Total Reads"]

    merged = meta[["research_id", "mean_insert_size", "ribosomal_bases_pct"]].merge(
        seqc2[["sample_id"] + COLS_OF_INTEREST + ["non_globin_frac"]],
        left_on="research_id", right_on="sample_id", how="inner")
    print(f"\nJoined on research_id == sample_id: {len(merged)} of {len(meta)} Stage-1 rows matched",
          file=sys.stderr)

    check_cols = ["Read Length", "Fragment GC Content Mean", "Fragment GC Content Std",
                  "Fragment GC Content Skewness", "Fragment GC Content Kurtosis",
                  "non_globin_frac", "Genes Detected", "Exonic Rate"]
    print(f"\n=== Correlation of mean_insert_size with RNA-SeQC2's richer metrics "
          f"(does any of these finally explain the bimodal split?) ===", file=sys.stderr)
    corr = merged[["mean_insert_size"] + check_cols].corr()["mean_insert_size"].drop("mean_insert_size")
    print(corr.round(3).sort_values(key=abs, ascending=False).to_string(), file=sys.stderr)

    print(f"\n=== rRNA Rate (RNA-SeQC2) vs ribosomal_bases_pct (rnaseq_metadata.tsv) "
          f"-- do two independently-computed metrics agree? ===", file=sys.stderr)
    print(f"  Pearson r: {merged['rRNA Rate'].corr(merged['ribosomal_bases_pct']):.3f}",
          file=sys.stderr)

    print(f"\n=== Non-globin read fraction (candidate Stage 4 sample-selection criterion) ===",
          file=sys.stderr)
    print(merged["non_globin_frac"].describe().to_string(), file=sys.stderr)

    print(f"\n=== Read Length value counts (is it actually uniform across the cohort?) ===",
          file=sys.stderr)
    print(merged["Read Length"].value_counts().to_string(), file=sys.stderr)

    # plot: mean_insert_size vs. the two candidates with the strongest |correlation|
    top2 = corr.abs().sort_values(ascending=False).index[:2].tolist()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, col in zip(axes, top2):
        ax.scatter(merged[col], merged["mean_insert_size"], s=4, alpha=0.15, color="#4C72B0")
        ax.set_xlabel(col, fontsize=9)
        ax.set_ylabel("mean_insert_size (bp)", fontsize=9)
        ax.set_title(f"r = {corr[col]:.3f}", fontsize=10)
    fig.suptitle("mean_insert_size vs. its two strongest RNA-SeQC2 correlates", fontsize=12)
    fig.tight_layout()
    fig_path = os.path.join(args.out_dir, "insert_size_vs_seqc2_metrics.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    print(f"\nWrote {fig_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
