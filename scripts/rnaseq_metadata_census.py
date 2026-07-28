#!/usr/bin/env python3
"""Stage 1 of the RNA-seq data-characterization plan (2026-07-25): a metadata-only census of
AoU's v9 RNA-seq cohort. Cheap and read-only -- no BAM content is ever touched, only the small
per-sample QC metadata table. This answers the "count / completeness / robustness /
reliability" half of the open questions before anything more expensive (BAM headers, actual
repertoire-extraction tooling) is worth running.

Reads:
  - v9/multiomics/rnaseq/rnaseq_metadata.tsv   (sampleid, research_id, alignment_rate_pct, rqs,
                                                 mrna_bases_pct, ribosomal_bases_pct,
                                                 processing_status, reads_aligned_in_pairs,
                                                 mean_insert_size, pipeline_id)
  - v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv   (research_id, ancestry_pred)
  - ~/pipeline_outputs/rnaseq_overlap/overlap_ids.tsv (from compute_wgs_rnaseq_overlap.py, if
    present -- restricts the ancestry breakdown to the 8,326-person "has all three" subgroup,
    the one this project would actually build a pilot on)

Reports (stderr): row count vs manifest count, QC pass/fail rate, missingness per column,
pipeline_id version spread (batch-effect flag), and ancestry breakdown of the RNA-seq-with-lrWGS
subgroup. Writes a 2x3 grid of QC-metric histograms to
~/pipeline_outputs/rnaseq_overlap/rnaseq_qc_distributions.png.

Usage (needs pandas + matplotlib -- the spechla env has both):
  pixi run -e spechla -- python3 scripts/rnaseq_metadata_census.py
"""
import argparse
import os
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RNASEQ_METADATA = "v9/multiomics/rnaseq/rnaseq_metadata.tsv"
RNASEQ_MANIFEST = "v9/multiomics/rnaseq/manifest.tsv"
ANCESTRY_TSV = "v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"

QC_NUMERIC_COLS = [
    "alignment_rate_pct", "rqs", "mrna_bases_pct", "ribosomal_bases_pct",
    "reads_aligned_in_pairs", "mean_insert_size",
]

# Plot metadata per column: (x-axis label with units, divisor for display scale, AoU's own
# published single-sample QC threshold if one exists -- from the "All of Us Genomics &
# Multi-omics Quality Report" found via web research 2026-07-25: RQS >= 5.5, alignment > 80%,
# mRNA bases > 20%. Drawing these as reference lines is *why* processing_status is 100% "Pass"
# in this delivered manifest -- AoU only ships samples that already cleared this bar, so the
# pass/fail column itself can't discriminate anymore; the real signal is the spread above it.
PLOT_META = {
    "alignment_rate_pct": ("% of reads aligned", 1, 80),
    "rqs": ("RNA Quality Score (RQS, 0-10 scale; higher = less-degraded RNA)", 1, 5.5),
    "mrna_bases_pct": ("% of bases in mRNA / exonic regions", 1, 20),
    "ribosomal_bases_pct": ("% of bases mapping to ribosomal RNA (lower = better depletion)", 1, None),
    "reads_aligned_in_pairs": ("aligned read pairs (millions) -- sequencing depth", 1_000_000, None),
    "mean_insert_size": ("mean insert size (bp)", 1, None),
}


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def require_cols(df, cols, label):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        die(f"{label}: expected column(s) {missing} not found. Actual columns: "
            f"{list(df.columns)}. AoU may have renamed something -- update the constants "
            f"at the top of this script.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    ap.add_argument("--out-dir", default=os.path.expanduser("~/pipeline_outputs/rnaseq_overlap"))
    ap.add_argument("--overlap-ids", default=os.path.expanduser(
        "~/pipeline_outputs/rnaseq_overlap/overlap_ids.tsv"),
        help="Output of compute_wgs_rnaseq_overlap.py. Optional -- skipped with a warning if absent.")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    meta_path = os.path.join(args.mount, RNASEQ_METADATA)
    manifest_path = os.path.join(args.mount, RNASEQ_MANIFEST)
    anc_path = os.path.join(args.mount, ANCESTRY_TSV)
    for p in (meta_path, manifest_path, anc_path):
        if not os.path.exists(p):
            die(f"not found: {p} -- is the gcsfuse mount up? Remount + `ls`-verify first "
                f"(ENVIRONMENT.md quirk #11/#14).")

    print(f"Reading RNA-seq metadata: {meta_path}", file=sys.stderr)
    meta = pd.read_csv(meta_path, sep="\t", dtype=str)
    require_cols(meta, ["research_id", "processing_status"] + QC_NUMERIC_COLS, "RNA-seq metadata")
    for c in QC_NUMERIC_COLS:
        meta[c] = pd.to_numeric(meta[c], errors="coerce")

    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str, usecols=["research_id"])
    manifest_ids = set(manifest["research_id"].dropna().unique())
    meta_ids = set(meta["research_id"].dropna().unique())

    print(f"\n=== Row counts ===", file=sys.stderr)
    print(f"  manifest.tsv:        {len(manifest)} rows, {len(manifest_ids)} unique research_id",
          file=sys.stderr)
    print(f"  rnaseq_metadata.tsv: {len(meta)} rows, {len(meta_ids)} unique research_id",
          file=sys.stderr)
    only_in_manifest = manifest_ids - meta_ids
    only_in_meta = meta_ids - manifest_ids
    print(f"  in manifest but no metadata row: {len(only_in_manifest)}", file=sys.stderr)
    print(f"  in metadata but no manifest row: {len(only_in_meta)}", file=sys.stderr)

    print(f"\n=== QC pass/fail (processing_status) ===", file=sys.stderr)
    print(meta["processing_status"].value_counts(dropna=False).to_string(), file=sys.stderr)

    print(f"\n=== Missingness (NaN count per QC column, of {len(meta)} rows) ===", file=sys.stderr)
    for c in QC_NUMERIC_COLS:
        n_na = meta[c].isna().sum()
        print(f"  {c}: {n_na} missing ({100*n_na/len(meta):.1f}%)", file=sys.stderr)

    print(f"\n=== pipeline_id spread (batch/version-effect flag) ===", file=sys.stderr)
    print(meta["pipeline_id"].value_counts(dropna=False).to_string(), file=sys.stderr)

    print(f"\n=== QC metric summary stats ===", file=sys.stderr)
    print(meta[QC_NUMERIC_COLS].describe().to_string(), file=sys.stderr)

    # --- ancestry breakdown, restricted to the "has all three data types" subgroup if we have it ---
    print(f"\nReading ancestry TSV: {anc_path}", file=sys.stderr)
    anc = pd.read_csv(anc_path, sep="\t", dtype=str, usecols=lambda c: c in ("research_id", "ancestry_pred"))
    require_cols(anc, ["research_id", "ancestry_pred"], "ancestry TSV")
    anc["ancestry_pred"] = anc["ancestry_pred"].astype(str).str.strip().str.upper()

    rna_with_ancestry = meta.merge(anc, left_on="research_id", right_on="research_id", how="left")
    print(f"\n=== Ancestry breakdown, full RNA-seq metadata cohort (n={len(rna_with_ancestry)}) ===",
          file=sys.stderr)
    print(rna_with_ancestry["ancestry_pred"].value_counts(dropna=False).to_string(), file=sys.stderr)

    if os.path.exists(args.overlap_ids):
        overlap = pd.read_csv(args.overlap_ids, sep="\t", dtype=str)
        overlap["has_lrwgs"] = overlap["has_lrwgs"].map({"True": True, "False": False})
        overlap["has_rnaseq"] = overlap["has_rnaseq"].map({"True": True, "False": False})
        all_three_ids = set(overlap.loc[overlap["has_lrwgs"] & overlap["has_rnaseq"], "person_id"])
        subgroup = rna_with_ancestry[rna_with_ancestry["research_id"].isin(all_three_ids)]
        print(f"\n=== Ancestry breakdown, RNA-seq-with-lrWGS subgroup (n={len(subgroup)}, "
              f"the pilot-relevant cohort) ===", file=sys.stderr)
        print(subgroup["ancestry_pred"].value_counts(dropna=False).to_string(), file=sys.stderr)
    else:
        print(f"\n(Skipped RNA-seq-with-lrWGS ancestry breakdown -- {args.overlap_ids} not found. "
              f"Run compute_wgs_rnaseq_overlap.py first if you want it.)", file=sys.stderr)

    # --- QC metric distributions -- labeled axes, real units, AoU's own QC threshold lines
    # where one exists (see PLOT_META comment) so "why does nobody fail QC" is visible, not
    # just stated in the console output ---
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    for ax, col in zip(axes.flat, QC_NUMERIC_COLS):
        xlabel, divisor, threshold = PLOT_META[col]
        data = meta[col].dropna() / divisor
        ax.hist(data, bins=40, color="#4C72B0", edgecolor="white", linewidth=0.3)
        if threshold is not None:
            ax.axvline(threshold / divisor, color="#C44E52", linestyle="--", linewidth=1.5,
                       label=f"AoU's own QC cutoff ({threshold})")
            ax.legend(fontsize=8, loc="upper left")
        ax.set_title(col, fontsize=10, fontweight="bold")
        ax.set_xlabel(xlabel, fontsize=8.5)
        ax.set_ylabel("number of samples")
    fig.suptitle(f"AoU v9 RNA-seq per-sample QC metrics (n={len(meta)} -- everyone shown here "
                 f"already passed AoU's own single-sample QC gate)", fontsize=12.5)
    fig.tight_layout()
    fig_path = os.path.join(args.out_dir, "rnaseq_qc_distributions.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    print(f"\nWrote {fig_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
