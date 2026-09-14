#!/bin/bash
# Batch TRUST4 with COPY-LOCAL staging -- the production runner for the full cohort.
#
# For each person: `gcloud storage cp` the BAM to local disk, run TRUST4 against the local
# copy, delete it, move on. Why: the gcsfuse mount streams a BAM at ~11 MB/s, so TRUST4's
# two full sequential passes over a 5-20 GB BAM take ~9-60 min and DON'T parallelize (one
# shared gcsfuse pipe -- Experiment B measured 9 people/hr at jobs=4, ~$925 / 38 days for
# 8,327). Native `gcloud storage cp` pulls the same BAM at ~670 MB/s (measured 2026-09-10,
# 4.8 GB in 15 s, even while gcsfuse was saturated -- it bypasses the fuse daemon). TRUST4
# then reads the local copy at GB/s, extraction becomes CPU-bound, and jobs=8-12 scales.
# Projected: ~$30-70 on-demand / ~$10-20 spot, 1-3 days (hours if sharded).
#
# RESUMABLE: skips anyone whose <research_id>_report.tsv already exists. Safe to Ctrl-C and
# re-run, and safe to shard -- run this on N VMs each with a slice of the cohort (or use
# build_rnaseq_cohort.py --skip), outputs are per-research_id so they never collide.
#
# DISK-BOUNDED: only --jobs BAMs are on local disk at once; each is deleted the instant its
# TRUST4 run finishes. Worst-case AoU RNA-seq BAM ~20 GB, so --jobs 8 wants ~180 GB free.
#
# COMPLIANCE: the BAM copy lands on the VM's own disk / staging dir -- inside the same
# controlled-tier boundary as the bucket. Same as gcsfuse's own file cache, same as the
# intermediate FASTQs TRUST4 already writes there. Never stages outside the workspace.
#
# Usage:
#   bash run_rnaseq_batch_local.sh <cohort.tsv> [--jobs N] [--staging DIR] [--threads-per-job T]
#                                   [--results-bucket gs://...]
#     cohort.tsv        research_id, ancestry, bam_rel_path  (from build_rnaseq_cohort.py)
#     --jobs N          concurrency AND max BAMs on disk at once   (default 8)
#     --staging DIR     where BAMs are copied                       (default ~/pipeline_outputs/rnaseq/_staging)
#     --threads-per-job T   -t passed to each run-trust4           (default: nproc / jobs, min 2)
#     --results-bucket  gs:// prefix each finished person's output dir is synced to, right
#                        after that person succeeds (default below). Safe across N sharded
#                        VMs with zero coordination: research_id is globally unique, so
#                        concurrent syncs from different machines never collide, and syncing
#                        per-person (not just at the end) means a preempted/killed VM never
#                        loses more than the one person it was mid-run on.
set -uo pipefail

COHORT="${1:?need a cohort tsv (research_id, ancestry, bam_rel_path) -- from build_rnaseq_cohort.py}"
shift || true

JOBS=8
STAGING="$HOME/pipeline_outputs/rnaseq/_staging"
TPJ=""
RESULTS_BUCKET="gs://aleix-disease-counts-wb-cordial-leechee-9743/rnaseq_results"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --jobs)            JOBS="${2:?--jobs needs a number}"; shift 2;;
    --staging)         STAGING="${2:?--staging needs a dir}"; shift 2;;
    --threads-per-job) TPJ="${2:?--threads-per-job needs a number}"; shift 2;;
    --results-bucket)  RESULTS_BUCKET="${2:?--results-bucket needs a gs:// path}"; shift 2;;
    *) echo "unknown arg: $1"; exit 1;;
  esac
done

REPO=~/repos/pilot-validation
SCRIPT_DIR="$REPO/aleix/RNA-seq/scripts"
BUCKET="gs://vwb-aou-datasets-controlled"
OUT_BASE="$HOME/pipeline_outputs/rnaseq"
LOG="$OUT_BASE/batch_local.log"
BILLING="${GOOGLE_PROJECT:?GOOGLE_PROJECT not set -- needed for the requester-pays bucket}"

# ---- preflight: fail loud and early on a fresh VM missing a tool, instead of 20 minutes
# into a batch. TRUST4 reference files are the one thing a fresh disk never has -- rebuild
# them automatically rather than erroring, since setup_trust4_refs.sh is idempotent and this
# is exactly what it's for.
command -v gcloud >/dev/null || { echo "FATAL: gcloud not on PATH"; exit 1; }
command -v pixi   >/dev/null || { echo "FATAL: pixi not on PATH -- install it before running this"; exit 1; }
[[ -s "$COHORT" ]] || { echo "FATAL: missing/empty cohort file: $COHORT"; exit 1; }
REF_DIR="$REPO/aleix/RNA-seq/reference"
if [[ ! -s "$REF_DIR/hg38_bcrtcr.fa" || ! -s "$REF_DIR/human_IMGT+C.fa" ]]; then
  echo "TRUST4 reference files missing (fresh disk) -- running setup_trust4_refs.sh ..."
  bash "$SCRIPT_DIR/setup_trust4_refs.sh" || { echo "FATAL: setup_trust4_refs.sh failed"; exit 1; }
fi
mkdir -p "$STAGING" "$OUT_BASE"

NPROC=$(nproc)
if [[ -z "$TPJ" ]]; then
  TPJ=$(( NPROC / JOBS )); (( TPJ < 2 )) && TPJ=2
fi

# disk guard -- worst-case 20 GB/BAM * JOBS, plus headroom for outputs
FREE_GB=$(df -BG --output=avail "$STAGING" 2>/dev/null | tail -1 | tr -dc '0-9')
NEED_GB=$(( JOBS * 22 ))
if [[ -n "$FREE_GB" ]] && (( FREE_GB < NEED_GB )); then
  echo "!! WARNING: ${FREE_GB} GB free at $STAGING, want >= ${NEED_GB} GB for --jobs ${JOBS}"
  echo "   (worst-case 20 GB BAMs). Lower --jobs, use --staging on a bigger disk, or free space."
  echo "   Continuing in 8 s -- Ctrl-C to abort."
  sleep 8
fi

N_TOTAL=$(($(wc -l < "$COHORT") - 1))
echo "==== RNA-seq COPY-LOCAL batch :: ${N_TOTAL} people :: jobs=${JOBS} :: -t ${TPJ}/job :: staging=${STAGING} :: results->${RESULTS_BUCKET} :: $(date) ====" | tee -a "$LOG"

run_one() {
  local research_id="$1" bam_rel="$2"
  local out="$OUT_BASE/$research_id"
  local report="$out/${research_id}_report.tsv"
  if [[ -s "$report" ]]; then
    echo "[$research_id] already done -- skip"
    return 0
  fi

  local lbam="$STAGING/${research_id}.bam"
  local errf="$STAGING/${research_id}.cp.err"
  local t0 t1 t2
  t0=$(date +%s)

  # BAM is required. `gcloud storage cp SRC1 SRC2 DST1 DST2` is NOT "copy each source to
  # its paired dest" -- with 2+ sources the LAST arg must be a single existing destination
  # DIRECTORY, so that four-arg form always fails validation. Two separate single-file
  # copies instead: BAM (required), then .bai (best-effort -- bam-extractor does a linear
  # sweep and does not need the index, so a missing/failed .bai copy is not fatal).
  if ! gcloud storage cp --billing-project "$BILLING" "$BUCKET/${bam_rel}" "$lbam" 2>"$errf"; then
    echo "[$research_id] !! COPY FAILED -- $(tail -1 "$errf" 2>/dev/null)"
    rm -f "$lbam" "${lbam}.bai" "$errf"
    return 1
  fi
  gcloud storage cp --billing-project "$BILLING" "$BUCKET/${bam_rel}.bai" "${lbam}.bai" \
    >/dev/null 2>>"$errf" || true
  t1=$(date +%s)

  TRUST4_THREADS="$TPJ" TRUST4_CLEAN=1 \
    bash "$SCRIPT_DIR/run_trust4_sample.sh" "$research_id" "$lbam" > "${out}.batch.log" 2>&1
  t2=$(date +%s)

  rm -f "$lbam" "${lbam}.bai" "$errf"

  if [[ -s "$report" ]]; then
    local n; n=$(($(wc -l < "$report") - 1))
    if gcloud storage cp --billing-project "$BILLING" -r "$out" \
        "${RESULTS_BUCKET}/" >/dev/null 2>"${out}.sync.err"; then
      rm -f "${out}.sync.err"
      echo "[$research_id] OK  copy $((t1-t0))s  trust4 $((t2-t1))s  total $((t2-t0))s  ${n} CDR3s  synced"
    else
      echo "[$research_id] OK (local only, BUCKET SYNC FAILED -- $(tail -1 "${out}.sync.err" 2>/dev/null))  ${n} CDR3s"
    fi
  else
    echo "[$research_id] !! TRUST4 FAILED after $((t2-t0))s -- see ${out}.batch.log"
    return 1
  fi
}
export -f run_one
export SCRIPT_DIR OUT_BASE BUCKET BILLING STAGING RESULTS_BUCKET

tail -n +2 "$COHORT" | cut -f1,3 | \
  xargs -P "$JOBS" -L1 bash -c 'run_one "$1" "$2"' _ \
  2>&1 | tee -a "$LOG"

DONE=$(tail -n +2 "$COHORT" | cut -f1 | while read -r r; do
         [[ -s "$OUT_BASE/$r/${r}_report.tsv" ]] && echo x; done | wc -l)
echo "==== batch finished :: ${DONE}/${N_TOTAL} have a report.tsv :: $(date) ====" | tee -a "$LOG"
echo "Left on disk (should be empty):"; ls -la "$STAGING" 2>/dev/null | tail -n +2
SYNC_FAILS=$(ls "$OUT_BASE"/*.sync.err 2>/dev/null | wc -l)
if (( SYNC_FAILS > 0 )); then
  echo "!! ${SYNC_FAILS} people finished but failed to sync to ${RESULTS_BUCKET} -- see $OUT_BASE/*.sync.err"
  echo "   Retry with:  for d in \$(ls $OUT_BASE/*.sync.err | sed 's/\\.sync\\.err\$//'); do gcloud storage cp --billing-project $BILLING -r \"\$d\" ${RESULTS_BUCKET}/ && rm \"\$d.sync.err\"; done"
else
  echo "All finished people synced to ${RESULTS_BUCKET}"
fi
echo "Aggregate (this shard only) with:  pixi run python3 $SCRIPT_DIR/aggregate_rnaseq_results.py $COHORT"
