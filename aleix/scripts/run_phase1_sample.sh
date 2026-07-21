#!/bin/bash
# Phase 1 (Experiment 1) — benchmark one truth genome through both long-read callers.
#
# For one sample: slice the HLA window from its remote *indexed* aligned BAM (HTTP range
# requests, no full download), convert to FASTQ, run SpecImmune AND HLA-Resolve, and print
# each tool's 8-gene calls so they can be scored against the Lai ground truth.
#
# Usage:  bash run_phase1_sample.sh <SAMPLE_ID> <ALIGNED_BAM_URL>
# Example (HG002):
#   bash run_phase1_sample.sh HG002 "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam"
#
# No `set -e`: we want to reach the completeness checks even if a tool misbehaves, because
# both tools can exit 0 while producing nothing (SpecImmune/sniffles, SpecHLA-blank lessons
# in ../context/ENVIRONMENT.md). Success is judged by per-gene completeness, not exit code.
set -uo pipefail

SAMPLE="${1:?need a sample id, e.g. HG002}"
BAM_URL="${2:?need the aligned BAM URL}"

REPO=~/repos/pilot-validation
SPECIMMUNE=~/tools/SpecImmune
HLARESOLVE=~/tools/hla_resolve
OUT=~/phase1/$SAMPLE
REGION="chr6:29500000-33500000"          # full classical-MHC window (generous — accuracy, not speed)
THREADS=$(nproc)

SI_MANIFEST=$REPO/pixi.toml               # Marc's manifest: spechla / specimmune envs
HR_MANIFEST=$REPO/aleix/pixi.toml         # HLA-Resolve env

mkdir -p "$OUT"
echo "==== Phase 1 :: $SAMPLE :: $(date) ===="
echo "region=$REGION threads=$THREADS out=$OUT"

# ---- 1. slice HLA window from the remote indexed BAM (range requests) ----
if [[ ! -s "$OUT/$SAMPLE.chr6.fastq.gz" ]]; then
  echo "[1/4] slicing $REGION from remote BAM ..."
  pixi run --manifest-path "$SI_MANIFEST" -e specimmune -- \
    samtools view -b -F 0x900 -o "$OUT/$SAMPLE.chr6.bam" "$BAM_URL" "$REGION" || { echo "SLICE FAILED"; exit 1; }
  # -F 0x900 drops secondary+supplementary so no read is duplicated in the FASTQ
  echo "[2/4] BAM -> FASTQ ..."
  pixi run --manifest-path "$SI_MANIFEST" -e specimmune -- \
    samtools fastq "$OUT/$SAMPLE.chr6.bam" 2>/dev/null | gzip > "$OUT/$SAMPLE.chr6.fastq.gz"
else
  echo "[1-2/4] FASTQ already present, skipping slice."
fi
NREADS=$(zcat "$OUT/$SAMPLE.chr6.fastq.gz" | wc -l); NREADS=$((NREADS/4))
echo "        reads in HLA window: $NREADS"
[[ "$NREADS" -lt 500 ]] && { echo "TOO FEW READS ($NREADS) — aborting $SAMPLE"; exit 1; }

# ---- 3a. SpecImmune (long-read) — minimap2 binning (Marc's adopted default), viz off ----
SI_OUT=$OUT/specimmune
SI_RESULT=$SI_OUT/$SAMPLE/${SAMPLE}.HLA.final.type.result.formatted.txt
if [[ ! -s "$SI_RESULT" ]]; then
  echo "[3/4] SpecImmune ..."
  ( cd "$SPECIMMUNE" && pixi run --manifest-path "$SI_MANIFEST" -e specimmune -- \
      python3 scripts/main.py -n "$SAMPLE" -o "$SI_OUT" -j "$THREADS" \
      -y pacbio-hifi -i HLA -r "$OUT/$SAMPLE.chr6.fastq.gz" --db ./db \
      --align_method_1 minimap2 --visualization "" ) \
      > "$OUT/${SAMPLE}.specimmune.log" 2>&1
else
  echo "[3/4] SpecImmune result already present, skipping."
fi

# ---- 3b. HLA-Resolve (long-read) — feeds reads, aligns itself; WGS scheme ----
HR_OUT=$OUT/hla_resolve
HR_RESULT=$HR_OUT/$SAMPLE/hla_typing_results/allele_output.csv
if [[ ! -s "$HR_RESULT" ]]; then
  echo "[4/4] HLA-Resolve ..."
  ( cd "$HLARESOLVE" && pixi run --manifest-path "$HR_MANIFEST" -- \
      hla_resolve --input_file "$OUT/$SAMPLE.chr6.fastq.gz" --sample_name "$SAMPLE" \
      --platform pacbio --scheme wgs --output_dir "$HR_OUT" --threads "$THREADS" ) \
      > "$OUT/${SAMPLE}.hlaresolve.log" 2>&1
else
  echo "[4/4] HLA-Resolve result already present, skipping."
fi

# ---- report: per-gene completeness + the raw calls (scored against Lai separately) ----
echo ""
echo "==== $SAMPLE :: RESULTS ===="
echo "--- SpecImmune (8 classical genes) ---"
if [[ -s "$SI_RESULT" ]]; then
  grep -E "^HLA-(A|B|C|DRB1|DQA1|DQB1|DPA1|DPB1)\b" "$SI_RESULT" | awk '{print $1"\t"$3}'
else
  echo "  !! NO SpecImmune result at $SI_RESULT — see ${SAMPLE}.specimmune.log"
fi
echo "--- HLA-Resolve (4-field) ---"
if [[ -s "$HR_RESULT" ]]; then
  cat "$HR_RESULT"
else
  echo "  !! NO HLA-Resolve result at $HR_RESULT — see ${SAMPLE}.hlaresolve.log"
fi
echo ""
echo "Done $SAMPLE. Prune bulky intermediates when satisfied:  rm -f $OUT/*.bam $OUT/*.fastq.gz"
