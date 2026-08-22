#!/bin/bash
# Experiment A -- does pre-copying the BAM to local disk before running TRUST4 beat
# running directly against the gcsfuse-mounted network path?
#
# WHY: every TRUST4 run in this project so far has read the BAM straight through gcsfuse.
# TRUST4's extractor does TWO FULL SEQUENTIAL PASSES over the file for paired-end data
# (see reference/TRUST4_DEEP_DIVE.md section 2.3 -- confirmed from source, not the docs).
# If gcsfuse's per-request network latency dominates over raw bandwidth, copying the BAM
# to local disk ONCE and then reading it twice locally could be much faster than reading
# it twice over the network. Never tested. Cheapest and most fundamental of the throughput
# experiments -- run it first; whichever read method wins becomes the baseline for the
# machine-type experiment (B) that follows.
#
# DESIGN DISCIPLINE: the SAME fixed set of people, same machine, same --jobs, run as two
# SEQUENTIAL arms -- never simultaneously, so they don't compete for the same network pipe
# and contaminate each other's timing. This is the fix for the earlier sweep's real flaw
# (different --jobs settings were compared across different, non-overlapping people, which
# confounded population composition with the setting being tested). Each arm's TRUST4
# output goes to a suffixed sample ID (<research_id>_fuse / <research_id>_local) so nothing
# collides with real batch output under ~/pipeline_outputs/rnaseq/<research_id>/, and
# nothing needs to be deleted between arms.
#
# DISK BUDGET, live-measured 2026-08-22: this VM has 28G free, and AoU RNA-seq BAMs run
# 4.8-6.3GB each (avg ~5.5GB) -- so bulk-copying an entire cohort at once, as an earlier
# version of this script did, would blow the disk budget past about 4 people. Fixed by
# batching the LOCAL arm: copy JOBS people, run TRUST4 on them, delete their local copies,
# move to the next batch of JOBS. This decouples reference-cohort SIZE from local disk
# capacity -- the cohort can be as large as we want; only one batch's worth of BAMs is
# ever on disk at a time. JOBS doubles as the batch size for exactly this reason, and is
# also used for the FUSE arm at the same value, deliberately -- the two arms must use
# identical concurrency or the comparison stops isolating read method and starts also
# comparing parallelism level, which is Experiment B's question, not this one's.
#
# Usage:  bash experiment_A_local_vs_fuse.sh <cohort.tsv> [JOBS]
#   cohort.tsv: research_id, ancestry, bam_rel_path -- from build_rnaseq_cohort.py, unmodified.
#   JOBS: concurrency AND local-arm batch size (default 3 -- sized so 3 * 6.3GB worst-case
#         stays well under this VM's 28G free; re-check `df -h ~` and adjust if disk
#         availability has changed).
set -uo pipefail

COHORT="${1:?need the reference cohort tsv (research_id, ancestry, bam_rel_path)}"
JOBS="${2:-3}"

REPO=~/repos/pilot-validation
SCRIPT_DIR="$REPO/aleix/RNA-seq/scripts"
MOUNT=~/mnt/aou-controlled
EXP_DIR=~/pipeline_outputs/rnaseq/experiments/expA
LOCAL_BAM_DIR="$EXP_DIR/local_bams"
RESULTS="$EXP_DIR/results.tsv"
CHUNK_DIR="$EXP_DIR/chunks"

[[ -s "$COHORT" ]] || { echo "MISSING cohort file: $COHORT"; exit 1; }
rm -rf "$CHUNK_DIR"
mkdir -p "$EXP_DIR" "$LOCAL_BAM_DIR" "$CHUNK_DIR"
printf 'research_id\tarm\tseconds\n' > "$RESULTS"

N_TOTAL=$(($(wc -l < "$COHORT") - 1))
echo "==== Experiment A :: $N_TOTAL people :: jobs=$JOBS :: $(date) ===="

# ---------- ARM 1: direct against the gcsfuse mount (current production method) ----------
run_fuse() {
  local research_id="$1" bam_rel="$2"
  local bam="$MOUNT/$bam_rel"
  if [[ ! -s "$bam" ]]; then
    echo "[$research_id/fuse] !! BAM not found at $bam"
    return 1
  fi
  local start end
  start=$(date +%s)
  bash "$SCRIPT_DIR/run_trust4_sample.sh" "${research_id}_fuse" "$bam" \
    > "$EXP_DIR/${research_id}_fuse.log" 2>&1
  end=$(date +%s)
  printf '%s\tfuse\t%d\n' "$research_id" "$((end-start))" >> "$RESULTS"
  echo "[$research_id/fuse] $((end-start))s"
}
export -f run_fuse
export SCRIPT_DIR MOUNT EXP_DIR RESULTS

echo ""
echo "---- ARM 1: FUSE (direct network-mounted read) ----"
ARM1_START=$(date +%s)
tail -n +2 "$COHORT" | cut -f1,3 | \
  xargs -P "$JOBS" -L1 bash -c 'run_fuse "$1" "$2"' _
ARM1_END=$(date +%s)
echo "ARM 1 (fuse) total wall time: $((ARM1_END-ARM1_START))s"

# ---------- ARM 2: batched local copy + run -- JOBS people on disk at a time, never more ----------
copy_one() {
  local research_id="$1" bam_rel="$2"
  local src="$MOUNT/$bam_rel"
  local dst="$LOCAL_BAM_DIR/${research_id}.bam"
  cp "$src" "$dst"
  cp "${src}.bai" "${dst}.bai" 2>/dev/null || true
}
export -f copy_one
export LOCAL_BAM_DIR

run_local() {
  local research_id="$1"
  local bam="$LOCAL_BAM_DIR/${research_id}.bam"
  if [[ ! -s "$bam" ]]; then
    echo "[$research_id/local] !! local copy missing at $bam"
    return 1
  fi
  local start end
  start=$(date +%s)
  bash "$SCRIPT_DIR/run_trust4_sample.sh" "${research_id}_local" "$bam" \
    > "$EXP_DIR/${research_id}_local.log" 2>&1
  end=$(date +%s)
  printf '%s\tlocal\t%d\n' "$research_id" "$((end-start))" >> "$RESULTS"
  echo "[$research_id/local] $((end-start))s"
}
export -f run_local
export SCRIPT_DIR EXP_DIR RESULTS

# research_id + bam_rel_path only (2 cols), then chopped into JOBS-sized chunk files.
tail -n +2 "$COHORT" | cut -f1,3 | split -l "$JOBS" -d -a 3 - "$CHUNK_DIR/chunk_"
N_CHUNKS=$(ls "$CHUNK_DIR"/chunk_* 2>/dev/null | wc -l)

echo ""
echo "---- ARM 2: LOCAL, batched $JOBS-at-a-time ($N_CHUNKS batches) to respect disk space ----"
COPY_S=0
ARM2_S=0
for chunk in "$CHUNK_DIR"/chunk_*; do
  n_here=$(wc -l < "$chunk")
  echo "  batch $(basename "$chunk"): $n_here people"

  CS=$(date +%s)
  cut -f1,2 "$chunk" | xargs -P "$JOBS" -L1 bash -c 'copy_one "$1" "$2"' _
  CE=$(date +%s)
  COPY_S=$((COPY_S + CE - CS))

  RS=$(date +%s)
  cut -f1 "$chunk" | xargs -P "$JOBS" -L1 bash -c 'run_local "$1"' _
  RE=$(date +%s)
  ARM2_S=$((ARM2_S + RE - RS))

  # free the disk before the next batch -- this is the whole point of batching.
  cut -f1 "$chunk" | while read -r rid; do
    rm -f "$LOCAL_BAM_DIR/${rid}.bam" "$LOCAL_BAM_DIR/${rid}.bam.bai"
  done
done
echo "Local copy stage total wall time (summed across batches): ${COPY_S}s"
echo "Local TRUST4 run total wall time (summed across batches): ${ARM2_S}s"

# ---------- summary ----------
ARM1_S=$((ARM1_END-ARM1_START))
EFFECTIVE_LOCAL=$((COPY_S + ARM2_S))
echo ""
echo "==== SUMMARY :: $(date) ===="
echo "Arm 1 (fuse, direct):                 ${ARM1_S}s total"
echo "Arm 2, copy stage:                    ${COPY_S}s"
echo "Arm 2, TRUST4 against local copies:   ${ARM2_S}s"
echo "Arm 2, effective total (copy + run):  ${EFFECTIVE_LOCAL}s"
echo ""
echo "Per-person results: $RESULTS"
if [[ "$EFFECTIVE_LOCAL" -lt "$ARM1_S" ]]; then
  echo "-> LOCAL-COPY WINS: pre-copying is faster even including the copy cost."
  echo "   This becomes the new baseline read method for Experiment B."
else
  echo "-> FUSE direct is as fast or faster. gcsfuse's sequential-read throughput is"
  echo "   already fine here -- the bottleneck is raw network bandwidth, not per-request"
  echo "   latency. Local copying is not worth adopting; keep FUSE as the baseline for B."
fi
echo ""
echo "Each batch's local BAM copies were deleted immediately after that batch's TRUST4"
echo "run, so nothing is left on disk -- $LOCAL_BAM_DIR should be empty."
