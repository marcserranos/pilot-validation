#!/bin/bash
# SpecHLA (short-read) padding sweep -- Experiment B, single person, controlled comparison.
# See EXPERIMENTS.md / DECISIONS.md 2026-07-10 entries for full reasoning. Key differences
# from the earlier SpecImmune-LR sweep (run_aligner_pad_sweep.sh), on purpose:
#   (1) Short-read fragment sizes (~150-500bp) are nowhere near long-read lengths (15-20kb) --
#       reusing the LR sweep's padding levels would never reach a real short-read floor. This
#       script measures the REAL insert-size distribution from this person's own data first
#       (samtools stats' IS histogram) and derives the floor/subfloor pad levels from it,
#       instead of guessing.
#   (2) Short reads are paired -- a narrow window can silently drop one mate of a pair. This
#       script checks properly-paired % per config (samtools flagstat), which the LR sweep
#       never needed (long reads aren't paired).
#   (3) No aligner dimension -- SpecHLA hardcodes bwa internally (ENVIRONMENT.md), unlike
#       SpecImmune's --align_method_1/2. Padding is the only variable here.
#   (4) SpecHLA's own README says to clear previous results before rerunning -- this script
#       gives every pad level its own fresh output directory, so that's satisfied by
#       construction, never by reusing a directory.
#
# Usage: bash run_spechla_pad_sweep.sh <person_id>
# Prereqs: ~/pipeline_outputs/<person_id>/chr6.bam must already exist (i.e. slice_and_fastq.sh
#   already ran for this person) and the existing SpecHLA fullpad baseline
#   (~/pipeline_outputs/<person_id>/spechla_output/<person_id>/hla.result.txt) should already
#   exist so it can be reused as the fullpad config, not rerun.
#
# Safe to walk away from: wrap in nohup (exact command at the bottom of this file). Does NOT
# use `set -e` -- one failed config is logged and skipped, the rest continues, and every
# completed run is written to disk immediately (not batched), so a mid-sweep VM sleep
# (ENVIRONMENT.md quirk #14) loses at most the in-flight run, never anything already finished.
# Every SpecHLA/python call is wrapped in `pixi run -e spechla` explicitly, regardless of the
# ambient shell -- same discipline as ENVIRONMENT.md quirk #15, applied here even though
# SpecHLA doesn't have SpecImmune's specific cross-env bug, because being explicit costs
# nothing and the alternative has already caused a real, costly failure once this project.

set -uo pipefail

PERSON_ID="${1:?Usage: bash run_spechla_pad_sweep.sh <person_id>}"
OUTDIR="$HOME/pipeline_outputs/$PERSON_ID"
SR_BAM="$OUTDIR/chr6.bam"
SPECHLA_DIR="$HOME/tools/SpecHLA"
PIXI_MANIFEST="$HOME/repos/pilot-validation/pixi.toml"
HELPERS="$HOME/repos/pilot-validation/spechla_pad_helpers.py"
COMPARE="$HOME/repos/pilot-validation/compare_hla_results.py"
SWEEP_DIR="$OUTDIR/spechla_sweep"
PROGRESS="$SWEEP_DIR/progress.log"
mkdir -p "$SWEEP_DIR"

log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$PROGRESS"; }

if [ ! -f "$SR_BAM" ]; then
  log "FATAL: $SR_BAM not found -- run slice_and_fastq.sh for $PERSON_ID first."
  exit 1
fi
samtools index -f "$SR_BAM"

BASELINE_RESULT="$OUTDIR/spechla_output/$PERSON_ID/hla.result.txt"
if [ ! -f "$BASELINE_RESULT" ]; then
  log "FATAL: no existing fullpad baseline at $BASELINE_RESULT -- run the SpecHLA runbook "
  log "  (ENVIRONMENT.md) for $PERSON_ID once at full padding before running this sweep."
  exit 1
fi

log "=== [$PERSON_ID] measuring real short-read insert-size distribution ==="
STATS=$(pixi run --manifest-path "$PIXI_MANIFEST" -e spechla -- \
  bash -c "samtools stats '$SR_BAM' | python3 '$HELPERS' insert-stats" 2>>"$PROGRESS")
echo "$STATS" | tee -a "$PROGRESS"
N_PAIRS=$(echo "$STATS" | awk -F= '/^n_pairs=/ {print $2}')
MEDIAN=$(echo "$STATS" | awk -F= '/^median=/ {print $2}')
P95=$(echo "$STATS" | awk -F= '/^p95=/ {print $2}')
P99=$(echo "$STATS" | awk -F= '/^p99=/ {print $2}')
MAXSZ=$(echo "$STATS" | awk -F= '/^max=/ {print $2}')

if [ -z "$P99" ] || [ "$P99" == "NA" ]; then
  log "FATAL: could not derive p99 insert size -- see stats above, aborting rather than guessing a floor."
  exit 1
fi

FLOOR_PAD=$(( ( (P99 + 49) / 50 ) * 50 ))     # p99 rounded up to nearest 50bp
SUBFLOOR_PAD=$(( FLOOR_PAD / 2 ))              # deliberately below the real floor -- risky control
log "Real fragment-size distribution: n_pairs=$N_PAIRS median=$MEDIAN p95=$P95 p99=$P99 max=$MAXSZ"
log "Derived floor_pad=${FLOOR_PAD}bp (p99 rounded up), subfloor_pad=${SUBFLOOR_PAD}bp (floor/2, expect real degradation)"

# Fixed structural levels (log-spaced, coarse -> fine) + the two data-derived levels above.
# More intermediate levels than the LR sweep's 5, per the roadmap's explicit ask.
FIXED_PADS=(500000 200000 100000 40000 10000 2000)

run_one() {
  local label="$1" pad="$2"
  local rundir="$SWEEP_DIR/$label"
  mkdir -p "$rundir"

  local windows
  windows=$(pixi run --manifest-path "$PIXI_MANIFEST" -e spechla -- python3 "$HELPERS" windows --pad "$pad")
  log "$label: pad=${pad}bp -> windows: $windows"

  local bam="$rundir/sliced.bam"
  samtools view -b "$SR_BAM" $windows -o "$bam" 2>>"$PROGRESS"
  local nreads pp_pct
  nreads=$(samtools view -c "$bam" 2>>"$PROGRESS")
  pp_pct=$(samtools flagstat "$bam" 2>>"$PROGRESS" | awk '/properly paired/ {print $6}' | tr -d '(%')
  log "$label: sliced $nreads reads, properly-paired=${pp_pct:-NA}%"

  samtools sort -n -o "$rundir/namesorted.bam" "$bam" 2>>"$PROGRESS"
  samtools fastq -1 "$rundir/R1.fastq.gz" -2 "$rundir/R2.fastq.gz" -0 /dev/null \
    -s "$rundir/singletons.fastq.gz" -F 0x900 "$rundir/namesorted.bam" 2>>"$PROGRESS"

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
    log "$label: FAILED (exit $rc) -- see $timing -- continuing to next config"
    return
  fi
  # Never trust exit code alone (ENVIRONMENT.md quirk #15) -- verify the real expected output.
  if [ ! -f "$result_file" ]; then
    log "$label: FAILED -- exit 0 but expected output missing at $result_file -- see $timing -- continuing to next config"
    return
  fi
  log "$label: SpecHLA done ($(grep real "$timing" | tail -1))"

  pixi run --manifest-path "$PIXI_MANIFEST" -e spechla -- python3 "$COMPARE" "$PERSON_ID" \
    --run-label "$label" \
    --spechla-result "$result_file" \
    >> "$SWEEP_DIR/comparisons.md" 2>>"$PROGRESS"
  log "$label: DONE"
}

log "=== sweep starting for person $PERSON_ID ==="

log "fullpad: reusing existing baseline (no rerun)"
pixi run --manifest-path "$PIXI_MANIFEST" -e spechla -- python3 "$COMPARE" "$PERSON_ID" \
  --run-label "fullpad" --spechla-result "$BASELINE_RESULT" >> "$SWEEP_DIR/comparisons.md" 2>>"$PROGRESS"

for pad in "${FIXED_PADS[@]}"; do
  run_one "pad${pad}" "$pad"
done
run_one "pad_floor_${FLOOR_PAD}" "$FLOOR_PAD"
run_one "pad_subfloor_${SUBFLOOR_PAD}" "$SUBFLOOR_PAD"

log "=== sweep complete: see $SWEEP_DIR/comparisons.md and $OUTDIR/comparison_log.csv ==="

# Optional push notification (ntfy.sh -- free, no signup, public topic-name-as-secret).
# Uncomment and pick your own unguessable topic name if you want a phone ping on completion:
# curl -s -d "SpecHLA sweep for $PERSON_ID finished" ntfy.sh/YOUR-UNGUESSABLE-TOPIC-NAME >/dev/null
