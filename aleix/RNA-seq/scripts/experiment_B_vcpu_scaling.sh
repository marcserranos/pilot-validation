#!/bin/bash
# Experiment B -- does a higher-vCPU machine actually buy more throughput, given the
# TRUST4_DEEP_DIVE.md finding that --jobs parallelism plateaus because we're saturating
# NETWORK EGRESS BANDWIDTH, not CPU (GCP allocates egress bandwidth per vCPU, roughly
# linearly, up to a per-machine cap)? Experiment A answers "which read method" (fuse vs
# local-copy); THIS answers "does the read method's ceiling move if the machine gets bigger".
#
# DESIGN: same fixed reference cohort as Experiment A (pass the identical cohort.tsv --
# this is the "same subcohort, always" discipline), same winning ARM from Experiment A,
# swept across several --jobs levels, run ONCE PER VM SIZE. This script does not resize
# the VM itself -- GCE machine-type changes require a stop/start cycle, which you do from
# the Workbench console between invocations. What this script DOES do: auto-detect the
# current vCPU count (`nproc`) and tag every result row with it, and APPEND (never
# overwrite) to a persistent results file on the persistent disk -- so running this script
# once on the 8-vCPU machine, then resizing, remounting, and running it again unchanged on
# a 16-vCPU machine, naturally accumulates a single comparable table.
#
# JOBS LEVELS: for a fair comparison, run the SAME jobs_csv on every machine size -- e.g.
# "2,4,8" on both the 8- and 16-vCPU machine -- so matched-concurrency wall times are
# directly comparable (isolates "does the SAME parallelism run faster on a bigger machine",
# i.e. per-job bandwidth). THEN, optionally, add one extra point unique to the bigger
# machine (e.g. jobs=16 on the 16-vCPU box) to see whether it unlocks a plateau the smaller
# machine physically cannot reach. Recommended: "2,4,8" on the 8-vCPU baseline, then
# "2,4,8,16" on the 16-vCPU machine.
#
# DISK: only matters for arm=local (see Experiment A's comment for the batching rationale).
# arm=fuse has no local-disk footprint at all -- prefer it here unless A crowned "local" the
# winner, since it removes one more variable (disk headroom) from this comparison.
#
# Usage:  bash experiment_B_vcpu_scaling.sh <cohort.tsv> <fuse|local> [jobs_csv]
#   cohort.tsv : research_id, ancestry, bam_rel_path -- THE SAME FILE used for Experiment A.
#   arm        : "fuse" or "local" -- whichever Experiment A's SUMMARY declared the winner.
#   jobs_csv   : comma-separated --jobs levels to sweep, e.g. "2,4,8" (default "2,4,8").
set -uo pipefail

COHORT="${1:?need the reference cohort tsv (research_id, ancestry, bam_rel_path) -- reuse Experiment A's}"
ARM="${2:?need fuse or local -- whichever arm Experiment A's summary crowned the winner}"
JOBS_CSV="${3:-2,4,8}"

[[ "$ARM" == "fuse" || "$ARM" == "local" ]] || { echo "arm must be 'fuse' or 'local', got: $ARM"; exit 1; }

REPO=~/repos/pilot-validation
SCRIPT_DIR="$REPO/aleix/RNA-seq/scripts"
MOUNT=~/mnt/aou-controlled
EXP_DIR=~/pipeline_outputs/rnaseq/experiments/expB
LOCAL_BAM_DIR="$EXP_DIR/local_bams"
RESULTS="$EXP_DIR/results.tsv"
CHUNK_DIR="$EXP_DIR/chunks"

[[ -s "$COHORT" ]] || { echo "MISSING cohort file: $COHORT"; exit 1; }
mkdir -p "$EXP_DIR" "$LOCAL_BAM_DIR"

VCPUS=$(nproc)
HOST=$(hostname)
if [[ ! -s "$RESULTS" ]]; then
  printf 'vcpus\thostname\tarm\tjobs\ttotal_seconds\tn_people\ttimestamp\n' > "$RESULTS"
fi

N_TOTAL=$(($(wc -l < "$COHORT") - 1))
echo "==== Experiment B :: vcpus=$VCPUS host=$HOST :: arm=$ARM :: $N_TOTAL people :: $(date) ===="
echo "Jobs levels to sweep: $JOBS_CSV"
echo "Results accumulate (append-only) at: $RESULTS"
echo ""

run_fuse_one() {
  local research_id="$1" bam_rel="$2" tag="$3"
  local bam="$MOUNT/$bam_rel"
  [[ -s "$bam" ]] || { echo "[$research_id] !! BAM not found at $bam"; return 1; }
  bash "$SCRIPT_DIR/run_trust4_sample.sh" "${research_id}_expB_${tag}" "$bam" \
    > "$EXP_DIR/${research_id}_${tag}.log" 2>&1
}
export -f run_fuse_one
export SCRIPT_DIR MOUNT EXP_DIR

copy_one() {
  local research_id="$1" bam_rel="$2"
  cp "$MOUNT/$bam_rel" "$LOCAL_BAM_DIR/${research_id}.bam"
  cp "$MOUNT/${bam_rel}.bai" "$LOCAL_BAM_DIR/${research_id}.bam.bai" 2>/dev/null || true
}
export -f copy_one
export LOCAL_BAM_DIR MOUNT

run_local_one() {
  local research_id="$1" tag="$2"
  local bam="$LOCAL_BAM_DIR/${research_id}.bam"
  [[ -s "$bam" ]] || { echo "[$research_id] !! local copy missing at $bam"; return 1; }
  bash "$SCRIPT_DIR/run_trust4_sample.sh" "${research_id}_expB_${tag}" "$bam" \
    > "$EXP_DIR/${research_id}_${tag}.log" 2>&1
}
export -f run_local_one
export SCRIPT_DIR EXP_DIR

IFS=',' read -ra JOBS_LEVELS <<< "$JOBS_CSV"
for JOBS in "${JOBS_LEVELS[@]}"; do
  TAG="j${JOBS}"
  echo "---- jobs=$JOBS (arm=$ARM, vcpus=$VCPUS) ----"
  START=$(date +%s)

  if [[ "$ARM" == "fuse" ]]; then
    tail -n +2 "$COHORT" | cut -f1,3 | \
      xargs -P "$JOBS" -L1 bash -c "run_fuse_one \"\$1\" \"\$2\" \"$TAG\"" _
  else
    rm -rf "$CHUNK_DIR"; mkdir -p "$CHUNK_DIR"
    tail -n +2 "$COHORT" | cut -f1,3 | split -l "$JOBS" -d -a 3 - "$CHUNK_DIR/chunk_"
    for chunk in "$CHUNK_DIR"/chunk_*; do
      cut -f1,2 "$chunk" | xargs -P "$JOBS" -L1 bash -c 'copy_one "$1" "$2"' _
      cut -f1 "$chunk" | xargs -P "$JOBS" -L1 bash -c "run_local_one \"\$1\" \"$TAG\"" _
      cut -f1 "$chunk" | while read -r rid; do
        rm -f "$LOCAL_BAM_DIR/${rid}.bam" "$LOCAL_BAM_DIR/${rid}.bam.bai"
      done
    done
  fi

  END=$(date +%s)
  ELAPSED=$((END-START))
  printf '%s\t%s\t%s\t%s\t%d\t%d\t%s\n' \
    "$VCPUS" "$HOST" "$ARM" "$JOBS" "$ELAPSED" "$N_TOTAL" "$(date -Iseconds)" >> "$RESULTS"
  echo "  jobs=$JOBS total wall time: ${ELAPSED}s ($N_TOTAL people)"
  echo ""
done

echo "==== DONE :: vcpus=$VCPUS :: $(date) ===="
echo ""
column -t "$RESULTS"
echo ""
echo "To compare against a different machine size: stop this environment in the Workbench"
echo "console, change the CPU/RAM (machine type), start it again, remount gcsfuse, and"
echo "re-run this exact command -- the persistent disk (and this results.tsv) survives the"
echo "resize, so results accumulate in one place instead of being scattered across runs."
