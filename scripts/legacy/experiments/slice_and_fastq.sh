#!/bin/bash
# Slices srWGS CRAM + lrWGS BAM to the chr6 HLA region and converts both to FASTQ,
# ready for SpecHLA (short-read) / SpecImmune (long-read).
#
# Usage: bash slice_and_fastq.sh <person_id> <cram_path_relative_to_mount> <lr_bam_path_relative_to_mount>
#
# Paths are relative to the gcsfuse mount root (~/mnt/aou-controlled), e.g.:
#   pooled/wgs/cram/v8_base/wgs_1017156.cram
#   pooled/longreads/v8_delta/BI/revio/bam/1017156/GRCh38/1017156.bam
# Get these from the notebook's manifest join (see CURRENT_SPRINT.md) -- print the
# full, untruncated path first, pandas truncates long gs:// strings by default.
#
# Requires: gcsfuse mount already active at ~/mnt/aou-controlled (see CURRENT_SPRINT.md
# quirk #11), and the local reference at ~/ref/Homo_sapiens_assembly38.fasta already
# downloaded. Run from inside either pixi env (spechla or specimmune) -- samtools is
# present in both.

set -euo pipefail

PERSON_ID="$1"
CRAM_PATH="$2"
LR_BAM_PATH="$3"

REGION="chr6:29500000-33500000"
REF=~/ref/Homo_sapiens_assembly38.fasta
MOUNT=~/mnt/aou-controlled
OUTDIR="$HOME/pipeline_outputs/${PERSON_ID}"
mkdir -p "$OUTDIR"

echo "=== [$PERSON_ID] Slicing srWGS CRAM to $REGION ==="
samtools view -b -T "$REF" "$MOUNT/$CRAM_PATH" "$REGION" -o "$OUTDIR/chr6.bam"
echo "  -> $(samtools view "$OUTDIR/chr6.bam" | wc -l) reads"

echo "=== [$PERSON_ID] Slicing lrWGS BAM to $REGION ==="
samtools view -b "$MOUNT/$LR_BAM_PATH" "$REGION" -o "$OUTDIR/chr6_LR.bam"
echo "  -> $(samtools view "$OUTDIR/chr6_LR.bam" | wc -l) reads"

echo "=== [$PERSON_ID] Short-read: name-sort + FASTQ ==="
samtools sort -n -o "$OUTDIR/chr6.namesorted.bam" "$OUTDIR/chr6.bam"
samtools fastq \
  -1 "$OUTDIR/R1.fastq.gz" \
  -2 "$OUTDIR/R2.fastq.gz" \
  -0 /dev/null \
  -s "$OUTDIR/singletons.fastq.gz" \
  -F 0x900 \
  "$OUTDIR/chr6.namesorted.bam"

echo "=== [$PERSON_ID] Long-read: FASTQ ==="
samtools fastq -F 0x900 "$OUTDIR/chr6_LR.bam" > "$OUTDIR/LR.fastq"

echo ""
echo "=== [$PERSON_ID] Done. FASTQ files ready in $OUTDIR ==="
echo ""
echo "Next, run SpecHLA (inside 'pixi shell -e spechla', from ~/tools/SpecHLA):"
echo "  bash script/whole/SpecHLA.sh -n $PERSON_ID -1 $OUTDIR/R1.fastq.gz -2 $OUTDIR/R2.fastq.gz -o $OUTDIR/spechla_output/"
echo ""
echo "Then SpecImmune (inside 'pixi shell -e specimmune', from ~/tools/SpecImmune;"
echo "using the proven bwa config, not the untested-and-flagged minimap2 one -- see COMPUTE_OPTIMIZATION_LOG.md):"
echo "  { time python3 scripts/main.py -n $PERSON_ID -o $OUTDIR/specimmune_output/ -j 4 -y pacbio-hifi -i HLA -r $OUTDIR/LR.fastq --db ./db ; } 2> $OUTDIR/specimmune_timing.txt"
