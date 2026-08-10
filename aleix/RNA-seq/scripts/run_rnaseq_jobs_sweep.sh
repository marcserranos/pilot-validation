#!/bin/bash
# Efficiency sweep: find where --jobs parallelism stops helping.
#
# Uses NON-OVERLAPPING slices of the real 100-person cohort for each --jobs level tested
# (tracked via a running line-offset into cohort.tsv), so this permanently contributes to
# completing that cohort -- unlike a reset-and-rerun design, nothing here is throwaway.
#
# We already have real data for --jobs 1 (666s/person effective, sequential baseline) and
# --jobs 3 (261s/person, 2.55x speedup, near-ideal) from the 3-person manual test earlier
# -- not re-tested here. This sweep starts at --jobs 4 and probes UPWARD PAST the 8-vCPU
# core count (4, 6, 8, 12, 16) specifically because the --jobs 3 result showed no CPU
# contention even when oversubscribed (24 threads requested on 8 cores, no slowdown) --
# strong evidence extraction is network/IO-bound, not CPU-bound. If speedup keeps
# improving smoothly past 8, that's confirmed. If it plateaus, that's the real ceiling.
#
# Chunk sizes scale with jobs level (~2x, so each config gets ~2 full scheduling waves)
# -- otherwise a high --jobs level with too few people to actually run concurrently would
# just silently test a lower effective concurrency and give a misleading result.
#
# Same "verify the real output file, never trust exit code alone" discipline as
# ../../../scripts/run_spechla_pad_sweep.sh -- continues to the next config on any
# individual failure rather than aborting the whole sweep.
#
# KNOWN LIMITATION: offset tracking is in-memory only. If this script itself gets killed
# and restarted (e.g. VM autostop), already-completed configs will be skipped near-
# instantly by run_rnaseq_batch.sh's own resumability check, producing a near-zero wall
# time for that config in summary.tsv -- an obviously wrong number, not a real result.
# If that happens, just discard that row when reading results, don't fix live -- not
# worth checkpoint logic for a one-off overnight sweep.
#
# Usage: bash run_rnaseq_jobs_sweep.sh [<cohort.tsv>]
set -uo pipefail

COHORT="${1:-$HOME/pipeline_outputs/rnaseq/cohort.tsv}"
REPO=~/repos/pilot-validation
SCRIPT_DIR="$REPO/aleix/RNA-seq/scripts"
SWEEP_DIR=~/pipeline_outputs/rnaseq/jobs_sweep
SUMMARY="$SWEEP_DIR/summary.tsv"
LOG="$SWEEP_DIR/sweep.log"

[[ -s "$COHORT" ]] || { echo "MISSING cohort: $COHORT"; exit 1; }
mkdir -p "$SWEEP_DIR"

BASELINE_SEC_PER_PERSON=666   # --jobs 1 manual test, 3 people, 2026-08-10

declare -A CHUNK=( [4]=8 [6]=12 [8]=16 [12]=24 [16]=32 )
JOBS_ORDER=(4 6 8 12 16)

TOTAL_PEOPLE=$(($(wc -l < "$COHORT") - 1))
OFFSET=1   # line 1 is the header; data starts at line 2

echo -e "jobs\tn_people\twall_seconds\tsec_per_person_effective\tspeedup_vs_seq" > "$SUMMARY"
echo "==== RNA-seq jobs sweep :: $(date) ====" | tee -a "$LOG"

for JOBS in "${JOBS_ORDER[@]}"; do
  N="${CHUNK[$JOBS]}"
  START=$((OFFSET + 1))
  END=$((OFFSET + N))
  if (( END > TOTAL_PEOPLE + 1 )); then
    END=$((TOTAL_PEOPLE + 1))
    N=$((END - START))
  fi
  if (( N <= 0 )); then
    echo "!! no people left in $COHORT for --jobs $JOBS -- skipping" | tee -a "$LOG"
    continue
  fi

  SLICE="$SWEEP_DIR/slice_jobs${JOBS}.tsv"
  { head -1 "$COHORT"; sed -n "${START},${END}p" "$COHORT"; } > "$SLICE"

  echo "" | tee -a "$LOG"
  echo "---- --jobs $JOBS :: $N people (cohort lines $START-$END) ----" | tee -a "$LOG"

  T0=$(date +%s)
  bash "$SCRIPT_DIR/run_rnaseq_batch.sh" "$SLICE" --jobs "$JOBS" >> "$LOG" 2>&1
  T1=$(date +%s)
  WALL=$((T1 - T0))

  # Verify actual output files exist -- never trust exit code alone.
  N_DONE=0
  while IFS= read -r RID; do
    [[ -s "$HOME/pipeline_outputs/rnaseq/$RID/${RID}_report.tsv" ]] && N_DONE=$((N_DONE + 1))
  done < <(tail -n +2 "$SLICE" | cut -f1)

  if (( N_DONE == 0 )); then
    echo "!! --jobs $JOBS :: 0/$N actually completed (wall ${WALL}s) -- see $LOG, continuing" \
      | tee -a "$LOG"
    echo -e "$JOBS\t$N\t$WALL\tFAILED\tFAILED" >> "$SUMMARY"
  else
    EFFECTIVE=$(awk -v w="$WALL" -v n="$N_DONE" 'BEGIN{printf "%.1f", w/n}')
    SPEEDUP=$(awk -v b="$BASELINE_SEC_PER_PERSON" -v e="$EFFECTIVE" 'BEGIN{printf "%.2f", b/e}')
    echo "--jobs $JOBS :: $N_DONE/$N completed, wall ${WALL}s, ${EFFECTIVE}s/person effective, ${SPEEDUP}x speedup" \
      | tee -a "$LOG"
    echo -e "$JOBS\t$N_DONE\t$WALL\t$EFFECTIVE\t$SPEEDUP" >> "$SUMMARY"
  fi

  OFFSET=$END
done

echo "" | tee -a "$LOG"
echo "==== sweep finished :: $(date) ====" | tee -a "$LOG"
echo "Summary (aggregate only, no research_ids -- safe to commit as-is):"
cat "$SUMMARY"
echo ""
REMAINING=$((TOTAL_PEOPLE + 1 - OFFSET))
echo "$REMAINING people in $COHORT not yet touched by this sweep or earlier tests."
echo "Pick the best --jobs from summary.tsv, then finish the rest:"
echo "  tail -n +$((OFFSET + 1)) $COHORT > $SWEEP_DIR/remaining.tsv"
echo "  { head -1 $COHORT; cat $SWEEP_DIR/remaining.tsv; } > $SWEEP_DIR/remaining_with_header.tsv"
echo "  bash $SCRIPT_DIR/run_rnaseq_batch.sh $SWEEP_DIR/remaining_with_header.tsv --jobs <BEST>"
