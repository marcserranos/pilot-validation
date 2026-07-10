#!/bin/bash
# Positive control for the SpecHLA padding sweep -- Experiment B sanity check, 2026-07-10.
#
# The padding sweep (run_spechla_pad_sweep.sh) found ZERO call degradation at any of the 8
# genes, all the way down to 325bp of padding. Before trusting that, this script deliberately
# tries to break the pipeline: instead of shrinking the flanking margin around each gene
# cluster (what the padding sweep did), it truncates INTO each cluster around its midpoint,
# to a total width well below what's needed to span a gene's diagnostic exons or even a
# single short-read fragment (median fragment length 419bp, measured 2026-07-10). If this
# also shows zero degradation, that's a real signal something in the pipeline is broken
# (e.g. silently reusing cached results, or a slicing bug). If it shows real degradation --
# missing calls, garbage identity, or SpecHLA erroring outright -- that confirms the earlier
# "no degradation" result was a genuine finding (short fragments need very little flanking
# margin), not a broken experiment.
#
# This is a blunt "can we make it break at all" check, not a refined floor estimate --
# windows are centered on each cluster's midpoint and may not respect exon boundaries.
#
# Usage: bash run_spechla_truncation_sanity_check.sh <person_id>
# Prereqs: same as run_spechla_pad_sweep.sh -- ~/pipeline_outputs/<person_id>/chr6.bam must
#   exist. Does not need a fresh insert-size measurement (reuses the 2026-07-10 measurement:
#   median 419bp, p99 649bp).

set -uo pipefail

PERSON_ID="${1:?Usage: bash run_spechla_truncation_sanity_check.sh <person_id>}"
OUTDIR="$HOME/pipeline_outputs/$PERSON_ID"
SR_BAM="$OUTDIR/chr6.bam"
SPECHLA_DIR="$HOME/tools/SpecHLA"
PIXI_MANIFEST="$HOME/repos/pilot-validation/pixi.toml"
HELPERS="$HOME/repos/pilot-validation/spechla_pad_helpers.py"
COMPARE="$HOME/repos/pilot-validation/compare_hla_results.py"
SWEEP_DIR="$OUTDIR/spechla_truncation_check"
PROGRESS="$SWEEP_DIR/progress.log"
mkdir -p "$SWEEP_DIR"

log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$PROGRESS"; }

if [ ! -f "$SR_BAM" ]; then
  log "FATAL: $SR_BAM not found -- run slice_and_fastq.sh for $PERSON_ID first."
  exit 1
fi
samtools index -f "$SR_BAM"

# Widths chosen to span "probably still fine" down to "definitely too small":
# 5000bp already cuts into the smaller gene clusters (Gene A alone is ~8.3kb); 650bp is
# about one fragment length (p99); 300bp and 100bp are well below a single fragment.
WIDTHS=(5000 1000 650 300 100)

log "=== truncation sanity check starting for person $PERSON_ID ==="
log "Reusing 2026-07-10 fragment-size measurement: median=419bp p95=586bp p99=649bp"

run_one() {
  local label="$1" width="$2"
  local rundir="$SWEEP_DIR/$label"
  mkdir -p "$rundir"

  local windows
  windows=$(pixi run --manifest-path "$PIXI_MANIFEST" -e spechla -- python3 "$HELPERS" truncated-windows --width "$width")
  log "$label: width=${width}bp -> windows: $windows"

  local bam="$rundir/sliced.bam"
  samtools view -b "$SR_BAM" $windows -o "$bam" 2>>"$PROGRESS"
  local nreads
  nreads=$(samtools view -c "$bam" 2>>"$PROGRESS")
  log "$label: sliced $nreads reads"

  if [ "$nreads" -eq 0 ]; then
    log "$label: zero reads sliced -- skipping SpecHLA (expected at the smallest widths)"
    return
  fi

  samtools sort -n -o "$rundir/namesorted.bam" "$bam" 2>>"$PROGRESS"
  local fastq_log="$rundir/fastq.log"
  samtools fastq -1 "$rundir/R1.fastq.gz" -2 "$rundir/R2.fastq.gz" -0 /dev/null \
    -s "$rundir/singletons.fastq.gz" -F 0x900 "$rundir/namesorted.bam" 2>"$fastq_log"
  cat "$fastq_log" >> "$PROGRESS"
  local discarded
  discarded=$(grep -oP 'discarded \K[0-9]+(?= singletons)' "$fastq_log" || echo 0)
  local dropout_pct
  dropout_pct=$(pixi run --manifest-path "$PIXI_MANIFEST" -e spechla -- python3 -c \
    "print(f'{$discarded / max($nreads, 1) * 100:.3f}')")
  log "$label: mate-dropout rate = ${discarded}/${nreads} = ${dropout_pct}%"

  local specout="$rundir/spechla_output"
  local timing="$rundir/timing.txt"
  local result_file="$specout/$PERSON_ID/hla.result.txt"
  log "$label: SpecHLA starting"
  ( cd "$SPECHLA_DIR" && \
    { time pixi run --manifest-path "$PIXI_MANIFEST" -e spechla -- \
        bash script/whole/SpecHLA.sh -n "$PERSON_ID" \
        -1 "$rundir/R1.fastq.gz" -2 "$rundir/R2.fastq.gz" -o "$specout" ; } \
      2> "$timing" )
  local rc=$?
  if [ $rc -ne 0 ]; then
    log "$label: FAILED (exit $rc) -- see $timing -- this IS an expected possible outcome at small widths, continuing"
    return
  fi
  if [ ! -f "$result_file" ]; then
    log "$label: FAILED -- exit 0 but expected output missing at $result_file -- see $timing -- continuing"
    return
  fi
  log "$label: SpecHLA done ($(grep real "$timing" | tail -1))"

  pixi run --manifest-path "$PIXI_MANIFEST" -e spechla -- python3 "$COMPARE" "$PERSON_ID" \
    --run-label "$label" \
    --spechla-result "$result_file" \
    >> "$SWEEP_DIR/comparisons.md" 2>>"$PROGRESS"
  log "$label: DONE"
}

for width in "${WIDTHS[@]}"; do
  run_one "trunc${width}" "$width"
done

log "=== truncation sanity check complete: see $SWEEP_DIR/comparisons.md ==="
