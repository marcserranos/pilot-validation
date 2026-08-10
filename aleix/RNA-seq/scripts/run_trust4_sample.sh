#!/bin/bash
# Run TRUST4 on one AoU RNA-seq BAM and report the recovered CDR3 count.
#
# [CONFIRMED, 2026-08-10] Proven end-to-end on person 1000291 -- 2,013 CDR3s recovered,
# both T-cell and B-cell chains, all 7 chain types represented (see
# ../results/1000291_trust4_smoke_test.md for the full writeup). The repertoire direction
# is viable on this data.
#
# --abnormalUnmapFlag is REQUIRED and easy to miss: AoU's STAR run used
# --outSAMunmapped Within (unmapped reads interleaved in-BAM near their mapped mate,
# rather than TRUST4's default expected layout). Without this flag, bam-extractor dies
# immediately with "Two reads from the unaligned fragment are not showing up together."
# Discovered live -- see chat history 2026-08-10, not documented anywhere in TRUST4's own
# README as a gotcha for STAR-aligned BAMs specifically.
#
# TCR/BCR loci sit on ordinary chromosomes (TRB chr7, IGH chr14, etc.), NOT on the
# ALT/HLA contigs AoU's RNA-seq reference excludes -- confirmed live, chr6 in this BAM's
# header is the correct full GRCh38 length (170,805,979 bp), and this run needed no
# ALT-contig workaround of any kind.
#
# Usage:  bash run_trust4_sample.sh <SAMPLE_ID> <path_to_rna_bam>
# The BAM path comes from ../reference or by resolving a research_id against
# gs://vwb-aou-datasets-controlled/v9/multiomics/rnaseq/manifest.tsv (columns: sampleid,
# research_id, markduplicates_bam_file_path, markduplicates_bam_index_path) -- mount the
# bucket first (see ../README.md), then pass the local gcsfuse path.
set -uo pipefail

SAMPLE="${1:?need a sample id}"
BAM="${2:?need a path to the RNA BAM}"

REPO=~/repos/pilot-validation
REF_DIR="$REPO/aleix/RNA-seq/reference"
MANIFEST="$REPO/aleix/RNA-seq/pixi.toml"
OUT=~/pipeline_outputs/rnaseq/$SAMPLE
THREADS=$(nproc)

BCRTCR_FA="$REF_DIR/hg38_bcrtcr.fa"
IMGT_FA="$REF_DIR/human_IMGT+C.fa"

rm -rf "$OUT"
mkdir -p "$OUT"
echo "==== TRUST4 :: $SAMPLE :: $(date) ===="

for f in "$BCRTCR_FA" "$IMGT_FA"; do
  [[ -s "$f" ]] || { echo "MISSING $f -- run setup_trust4_refs.sh first"; exit 1; }
done
[[ -s "$BAM" ]] || { echo "MISSING BAM at $BAM"; exit 1; }

pixi run --manifest-path "$MANIFEST" -- \
  run-trust4 -b "$BAM" -f "$BCRTCR_FA" --ref "$IMGT_FA" \
  --abnormalUnmapFlag \
  -t "$THREADS" -o "$SAMPLE" --od "$OUT" \
  > "$OUT/${SAMPLE}.trust4.log" 2>&1

REPORT="$OUT/${SAMPLE}_report.tsv"
if [[ -s "$REPORT" ]]; then
  N_CDR3=$(tail -n +2 "$REPORT" | wc -l)
  echo ""
  echo "==== $SAMPLE :: RESULT ===="
  echo "CDR3s recovered: $N_CDR3"
  echo "(full report: $REPORT -- columns include chain, CDR3aa, count/abundance per TRUST4's docs)"
else
  echo "!! no report at $REPORT -- see $OUT/${SAMPLE}.trust4.log"
fi
