#!/bin/bash
# Batch-run TRUST4 across a cohort TSV (research_id, ancestry, bam_rel_path -- from
# build_rnaseq_cohort.py). Resumable: skips anyone whose report.tsv already exists, so
# re-running after an interruption picks up where it left off instead of redoing everyone.
#
# Sequential by default (--jobs 1) -- matches the one proven data point we have (person
# 1000291, ~9 min, mostly the extraction stage scanning the whole BAM over the network
# mount). Before trusting a --jobs N > 1 estimate for a full run, TEST IT: run the same
# small slice at --jobs 1 vs --jobs 3 and compare wall-clock time. If extraction is really
# network-bound (streaming a huge file through gcsfuse), not CPU-bound, higher --jobs may
# not scale the way more vCPUs would suggest -- untested, don't assume either way.
#
# Usage:  bash run_rnaseq_batch.sh <cohort.tsv> [--jobs N]
set -uo pipefail

COHORT="${1:?need a cohort tsv (research_id, ancestry, bam_rel_path)}"
JOBS=1
if [[ "${2:-}" == "--jobs" ]]; then JOBS="${3:?--jobs needs a number}"; fi

REPO=~/repos/pilot-validation
SCRIPT_DIR="$REPO/aleix/RNA-seq/scripts"
MOUNT=~/mnt/aou-controlled
LOG=~/pipeline_outputs/rnaseq/batch.log

[[ -s "$COHORT" ]] || { echo "MISSING cohort file: $COHORT"; exit 1; }
mkdir -p ~/pipeline_outputs/rnaseq

N_TOTAL=$(($(wc -l < "$COHORT") - 1))
echo "==== RNA-seq batch :: $N_TOTAL people :: jobs=$JOBS :: $(date) ====" | tee -a "$LOG"

run_one() {
  local research_id="$1" bam_rel="$2"
  local out="$HOME/pipeline_outputs/rnaseq/$research_id"
  local report="$out/${research_id}_report.tsv"
  if [[ -s "$report" ]]; then
    echo "[$research_id] already done, skipping"
    return 0
  fi
  local bam="$MOUNT/$bam_rel"
  if [[ ! -s "$bam" ]]; then
    echo "[$research_id] !! BAM not found at $bam -- skipping"
    return 1
  fi
  local start end
  start=$(date +%s)
  bash "$SCRIPT_DIR/run_trust4_sample.sh" "$research_id" "$bam" > "${out}.batch.log" 2>&1
  end=$(date +%s)
  if [[ -s "$report" ]]; then
    local n=$(($(wc -l < "$report") - 1))
    echo "[$research_id] done in $((end-start))s -- $n CDR3s"
  else
    echo "[$research_id] !! FAILED after $((end-start))s -- see ${out}.batch.log"
  fi
}
export -f run_one
export SCRIPT_DIR MOUNT

tail -n +2 "$COHORT" | cut -f1,3 | \
  xargs -P "$JOBS" -L1 bash -c 'run_one "$1" "$2"' _ \
  2>&1 | tee -a "$LOG"

echo "==== batch finished :: $(date) ====" | tee -a "$LOG"
echo "Aggregate with:  pixi run python3 $SCRIPT_DIR/aggregate_rnaseq_results.py $COHORT"
