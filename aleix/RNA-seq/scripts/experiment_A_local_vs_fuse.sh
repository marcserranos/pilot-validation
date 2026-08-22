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
# BEFORE RUNNING: check local disk space (`df -h ~`). This copies every person's BAM to
# local disk once -- if BAMs are large and the reference cohort is big, this can fill the
# disk. Size the reference cohort to what your disk actually has room for; this script
# does not manage disk space for you.
#
# Usage:  bash experiment_A_local_vs_fuse.sh <cohort.tsv> [JOBS]
#   cohort.tsv: research_id, ancestry, bam_rel_path -- from build_rnaseq_cohort.py, unmodified.
#   JOBS: parallelism, same value used for both arms (default 8).
set -uo pipefail

COHORT="${1:?need the reference cohort tsv (research_id, ancestry, bam_rel_path)}"
JOBS="${2:-8}"

REPO=~/repos/pilot-validation
SCRIPT_DIR="$REPO/aleix/RNA-seq/scripts"
MOUNT=~/mnt/aou-controlled
EXP_DIR=~/pipeline_outputs/rnaseq/experiments/expA
LOCAL_BAM_DIR="$EXP_DIR/local_bams"
RESULTS="$EXP_DIR/results.tsv"

[[ -s "$COHORT" ]] || { echo "MISSING cohort file: $COHORT"; exit 1; }
mkdir -p "$EXP_DIR" "$LOCAL_BAM_DIR"
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

# ---------- Stage: bulk-copy every BAM to local disk, timed separately from the run ----------
echo ""
echo "---- STAGING: parallel copy to local disk ----"
copy_one() {
  local research_id="$1" bam_rel="$2"
  local src="$MOUNT/$bam_rel"
  local dst="$LOCAL_BAM_DIR/${research_id}.bam"
  cp "$src" "$dst"
  cp "${src}.bai" "${dst}.bai" 2>/dev/null || true
}
export -f copy_one
export LOCAL_BAM_DIR

COPY_START=$(date +%s)
tail -n +2 "$COHORT" | cut -f1,3 | \
  xargs -P "$JOBS" -L1 bash -c 'copy_one "$1" "$2"' _
COPY_END=$(date +%s)
echo "Local copy stage total wall time: $((COPY_END-COPY_START))s"

# ---------- ARM 2: against the local copies ----------
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

echo ""
echo "---- ARM 2: LOCAL (pre-copied to disk) ----"
ARM2_START=$(date +%s)
tail -n +2 "$COHORT" | cut -f1 | \
  xargs -P "$JOBS" -L1 bash -c 'run_local "$1"' _
ARM2_END=$(date +%s)
echo "ARM 2 (local) total wall time: $((ARM2_END-ARM2_START))s"

# ---------- summary ----------
COPY_S=$((COPY_END-COPY_START))
ARM1_S=$((ARM1_END-ARM1_START))
ARM2_S=$((ARM2_END-ARM2_START))
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
echo "Local BAM copies are still on disk at $LOCAL_BAM_DIR -- delete with"
echo "  rm -rf $LOCAL_BAM_DIR"
echo "if you need the space back."
