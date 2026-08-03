#!/bin/bash
# Core/concurrency/memory/disk scaling experiment v2 (2026-08-02) -- redesigned to run on a
# purpose-built test VM (not just the original 4-vCPU baseline), and to answer three questions at
# once before committing real money to the ~13,000-14,500 person production run:
#   1. Does wall-clock scale linearly with concurrency, or does something (disk/network I/O,
#      shared-mount contention) bite at higher concurrency levels? (never measured past 4 cores)
#   2. Does running N people concurrently multiply memory use ~linearly, or is there sharing
#      (e.g. OS page-cache reuse of the tiny, identical IPD-IMGT/HLA reference across processes)?
#      This matters specifically because High-CPU machine types (the cost-optimal family per
#      DECISIONS.md/reports/immuannot_pilot) have the THINNEST memory ratio (1GB/vCPU on n2-highcpu)
#      of any family -- if per-process memory turns out to be more than that, high concurrency on a
#      High-CPU box would OOM, and we'd need to fall back to n2-standard (4GB/vCPU) instead.
#   3. What's the real per-person disk footprint, so the production data-disk size can be picked
#      with a safe margin instead of guessed from reading the pipeline's own source code.
#
# THREE phases:
#   Phase 0 (baseline profile, one person, ~7-10 min): wraps a single-person run in `/usr/bin/time
#     -v` for a clean peak-RSS number, and `du -sb` on that person's own output folder before/after
#     -- a clean, uncontended "per-unit" measurement to sanity-check the concurrent numbers against.
#   Phase 1 (thread-vs-concurrency comparison, small batch, cheap): only TWO anchor core budgets
#     (8 and 32 -- a small one and the largest one), each swept across every (concurrency, threads)
#     split, on a small fixed 2-person pool -- same "does the tool's own threading scale as well as
#     running separate people" question the original design asked, now extended past 4 cores for the
#     first time. Deliberately NOT a fine sweep across every budget (2026-08-02, cost-trim pass,
#     Marc): the threading-vs-concurrency trade-off is a smooth trend, not something expected to
#     flip between budget levels, so two anchor points show whether the pattern holds at scale just
#     as well as five points do, for a fraction of the wall-clock. Budget 4's "1 person, 4 threads"
#     point is deliberately skipped too -- Phase 0 already measures exactly that, for free, so
#     re-measuring it here would just be paying twice for the same number. The 2-person batch also
#     means any (concurrency, threads) split needing more than 2 concurrent people gets auto-skipped
#     by run_config's own people-count check -- which is intentional, not a gap: those configs
#     (e.g. "8 people, 1 thread each") measure the exact same thing as a Phase 2 row at that same
#     concurrency level, just with a smaller/less representative batch, so letting them fall away
#     here avoids ever paying for that redundant measurement twice.
#   Phase 2 (horizontal stress test, one row per level, threads=1 always, UNTRIMMED -- this is the
#     phase that actually answers the production-relevant question, so it keeps its full resolution
#     even in the cost-trimmed design above): for increasing concurrency levels up to min(nproc,
#     people supplied), runs that many people FULLY IN PARALLEL (batch size == concurrency, so
#     wall-clock per row stays ~1 person's runtime regardless of level) -- the number that actually
#     answers "does this VM contend on disk/network/memory at high concurrency," which no prior run
#     has ever measured.
#
# Every row in Phase 1/2 also records: min free memory observed during the row (via a background
# /proc/meminfo sampler), and the pipeline_outputs disk-usage delta across the row -- so the output
# TSV is a single source for the timing AND resource-planning numbers.
#
# CONCURRENCY-SAFETY NOTE: each concurrent worker is a genuinely separate OS process (xargs -P),
# and run_immuannot_person.py's normal output files are NOT safe for concurrent writers (see its
# own --out-suffix docstring -- a real race was found and fixed 2026-08-02, before this script was
# ever run past 4 cores). This script always passes a unique --out-suffix per worker (the person_id
# itself) so concurrent runs never collide; the per-worker fragment files are left on disk, not
# auto-merged (merge with `pandas.concat` over `immuannot_calls.*.tsv` if you want the real calls
# later -- this experiment only cares about the timing/resource numbers in core_scaling_experiment.tsv).
#
# Usage (from ~/repos/pilot-validation, inside `pixi shell -e specimmune`):
#   bash scripts/run_core_scaling_experiment.sh <person_id_1> <person_id_2> ... <person_id_N>
# Supply as many distinct person_ids as you can (ideally >= this VM's vCPU count) -- Phase 2's
# highest stress level is capped at whichever is smaller: nproc or the number of ids you pass.
#
# Output: appends to ~/pipeline_outputs/core_scaling_experiment.tsv -- aggregate-only timing +
# resource numbers (no genotypes), safe to paste back off the VM. Baseline profile (Phase 0) goes
# to ~/pipeline_outputs/immuannot_baseline_profile.txt (also aggregate-only).

set -uo pipefail

OUTROOT="$HOME/pipeline_outputs"
LOG="$OUTROOT/core_scaling_experiment.tsv"
BASELINE_LOG="$OUTROOT/immuannot_baseline_profile.txt"
SMALL_BATCH=2
IMMUANNOT_EXTRA_ARGS="${IMMUANNOT_EXTRA_ARGS:---force}"

if [ "$#" -lt 4 ]; then
  echo "FATAL: need at least 4 person_ids, got $#." >&2
  echo "Usage: bash scripts/run_core_scaling_experiment.sh <person_id_1> ... <person_id_N>" >&2
  exit 1
fi

ALL_PEOPLE=("$@")
N_PEOPLE=${#ALL_PEOPLE[@]}
NPROC=$(nproc)
TOTAL_MEM_MB=$(($(grep MemTotal /proc/meminfo | awk '{print $2}') / 1024))
TOP_STRESS_LEVEL=$(( NPROC < N_PEOPLE ? NPROC : N_PEOPLE ))
echo "This VM has $NPROC vCPUs (nproc), ${TOTAL_MEM_MB}MB total memory; $N_PEOPLE person_ids supplied." >&2
echo "Phase 2's top stress level is capped at min(nproc, N_PEOPLE) = $TOP_STRESS_LEVEL." >&2

mkdir -p "$OUTROOT"
if [ ! -f "$LOG" ]; then
  printf "timestamp\tphase\tconcurrency\tthreads_per_person\ttotal_cores_used\tbatch_size\twall_clock_seconds\tseconds_per_person_wallclock\tmem_avail_min_mb\tmem_used_peak_mb\tmem_used_peak_pct\tdisk_delta_mb\n" > "$LOG"
fi

# --- Single-instance lock (added 2026-08-03 after a real incident: three separate invocations of
# this script ended up running concurrently against the same person_ids -- likely from the launch
# command being re-run/re-pasted during troubleshooting -- and silently clobbered each other's
# shared per-person intermediate files (the trimmed FASTA / .gtf.gz under
# ~/pipeline_outputs/<pid>/immuannot_output/, which are NOT isolated by --out-suffix -- only the
# final results files are). One instance ended up stalled with zero child processes and zero
# output for 45+ minutes, costing real debugging time before the duplicate was even found. This
# lock makes that impossible to repeat silently: a second launch now fails loudly and immediately
# instead of quietly contending for the same files. ---
LOCK_FILE="$OUTROOT/core_scaling_experiment.lock"
if [ -f "$LOCK_FILE" ]; then
  OLD_PID=$(cat "$LOCK_FILE" 2>/dev/null)
  if [ -n "$OLD_PID" ] && kill -0 "$OLD_PID" 2>/dev/null; then
    echo "FATAL: another instance of this script is already running (PID $OLD_PID, per $LOCK_FILE)." >&2
    echo "Refusing to start a second one -- two concurrent instances WILL corrupt each other's" >&2
    echo "shared per-person intermediate files (real incident, 2026-08-03). If you're certain the" >&2
    echo "other instance is not actually needed: kill $OLD_PID first, confirm with" >&2
    echo "  ps -ef --forest | grep $OLD_PID" >&2
    echo "that it and all its children are gone, THEN rm $LOCK_FILE, THEN retry." >&2
    exit 1
  else
    echo "NOTE: found a stale lock file (PID $OLD_PID is not running) -- removing it and proceeding." >&2
    rm -f "$LOCK_FILE"
  fi
fi
echo "$$" > "$LOCK_FILE"
trap 'rm -f "$LOCK_FILE"' EXIT

# --- Background memory sampler: appends free-memory (MB) to a temp file every 2s until killed. ---
start_mem_sampler() {
  local outfile="$1"
  : > "$outfile"
  ( while true; do
      awk '/MemAvailable/ {print $2/1024}' /proc/meminfo >> "$outfile"
      sleep 2
    done ) &
  echo $!  # sampler PID, for the caller to kill later
}

stop_mem_sampler_and_get_min() {
  local sampler_pid="$1" outfile="$2"
  kill "$sampler_pid" 2>/dev/null
  wait "$sampler_pid" 2>/dev/null
  if [ -s "$outfile" ]; then
    sort -n "$outfile" | head -1
  else
    echo ""
  fi
}

disk_usage_mb() {
  du -sm "$OUTROOT" 2>/dev/null | awk '{print $1}'
}

# --- Diagnostic checkpoint trail (added 2026-08-03, Marc: "add diagnosis blocks if you need
# more information" -- a prior run produced zero rows in $LOG despite the VM's restart-only-
# after-exit behavior confirming it ran to completion or died, and there was no persistent
# record of how far it actually got, only the two data-bearing files). Writes to its OWN small
# file, separately from $LOG/$BASELINE_LOG, so it survives even if this run is NOT wrapped in
# `tee` -- always-on, not opt-in, because the whole point is to not depend on remembering to
# capture output correctly. ---
CHECKPOINT_LOG="$OUTROOT/core_scaling_checkpoints.log"
checkpoint() { echo "$(date -Iseconds) PID=$$ -- $1" >> "$CHECKPOINT_LOG"; }
checkpoint "script started, argv: $*"

# ============ PHASE 0: single-person baseline profile (clean, uncontended) ============
checkpoint "entering Phase 0"
echo "" >&2
echo ">>> PHASE 0: baseline profile, person=${ALL_PEOPLE[0]}, /usr/bin/time -v + disk delta <<<" >&2
BASELINE_PERSON="${ALL_PEOPLE[0]}"
disk_before=$(disk_usage_mb)
{
  echo "=== Baseline profile: $(date -Iseconds) ==="
  echo "person_id=$BASELINE_PERSON threads=4 (production default) nproc=$NPROC total_mem_mb=$TOTAL_MEM_MB"
} >> "$BASELINE_LOG"
if command -v /usr/bin/time >/dev/null 2>&1; then
  /usr/bin/time -v pixi run -e specimmune -- \
    python3 scripts/run_immuannot_person.py "$BASELINE_PERSON" --threads 4 \
    --out-suffix ".baseline" --force >> "$BASELINE_LOG" 2>&1
  grep -E "Maximum resident set size|Elapsed \(wall clock\)|Percent of CPU" "$BASELINE_LOG" | tail -3 >&2
else
  echo "NOTE: /usr/bin/time not found -- skipping peak-RSS capture, still measuring disk." >&2
  pixi run -e specimmune -- python3 scripts/run_immuannot_person.py "$BASELINE_PERSON" \
    --threads 4 --out-suffix ".baseline" --force >> "$BASELINE_LOG" 2>&1
fi
disk_after=$(disk_usage_mb)
person_dir="$OUTROOT/$BASELINE_PERSON/immuannot_output"
person_dir_mb="-"
[ -d "$person_dir" ] && person_dir_mb=$(du -sm "$person_dir" 2>/dev/null | awk '{print $1}')
{
  echo "pipeline_outputs disk delta (whole tree, includes any concurrent noise): $((disk_after - disk_before)) MB"
  echo "this person's immuannot_output/ folder size: ${person_dir_mb} MB"
  echo ""
} >> "$BASELINE_LOG"
echo "Baseline: this person's immuannot_output/ = ${person_dir_mb} MB. Full log: $BASELINE_LOG" >&2
echo "Per-person disk projection (13,000-14,521 people) at this rate: see analysis step after the run." >&2
checkpoint "Phase 0 complete"

run_config() {
  local phase="$1" concurrency="$2" threads="$3"; shift 3
  local people=("$@")
  checkpoint "run_config ENTER phase=$phase concurrency=$concurrency threads=$threads people=${people[*]}"
  local total_cores=$((concurrency * threads))
  if [ "$total_cores" -gt "$NPROC" ]; then
    echo "SKIP phase=$phase concurrency=$concurrency threads=$threads (total_cores=$total_cores > nproc=$NPROC)" >&2
    checkpoint "run_config SKIP (total_cores > nproc)"
    return
  fi
  if [ "${#people[@]}" -lt "$concurrency" ]; then
    echo "SKIP phase=$phase concurrency=$concurrency threads=$threads (need >= $concurrency people, have ${#people[@]})" >&2
    checkpoint "run_config SKIP (not enough people)"
    return
  fi
  echo "=== phase=$phase concurrency=$concurrency threads=$threads (total_cores=$total_cores), batch=${people[*]} ===" >&2

  local sample_file
  sample_file=$(mktemp)
  local sampler_pid
  sampler_pid=$(start_mem_sampler "$sample_file")
  local disk_before disk_after
  disk_before=$(disk_usage_mb)

  local t0 t1 wall
  t0=$(date +%s)
  checkpoint "run_config: about to launch xargs -P $concurrency (this is the line most likely to hang or die silently)"
  printf "%s\n" "${people[@]}" | \
    xargs -P "$concurrency" -I{} pixi run -e specimmune -- \
      python3 scripts/run_immuannot_person.py {} --threads "$threads" \
      --out-suffix ".worker_{}" $IMMUANNOT_EXTRA_ARGS
  local xargs_exit=$?
  t1=$(date +%s)
  wall=$((t1 - t0))
  checkpoint "run_config: xargs returned exit=$xargs_exit after ${wall}s"

  local mem_avail_min mem_used_peak mem_used_peak_pct
  mem_avail_min=$(stop_mem_sampler_and_get_min "$sampler_pid" "$sample_file")
  rm -f "$sample_file"
  if [ -n "$mem_avail_min" ]; then
    mem_used_peak=$(python3 -c "print(f'{$TOTAL_MEM_MB - $mem_avail_min:.0f}')")
    mem_used_peak_pct=$(python3 -c "print(f'{100*($TOTAL_MEM_MB - $mem_avail_min)/$TOTAL_MEM_MB:.1f}')")
  else
    mem_used_peak="-"; mem_used_peak_pct="-"
  fi
  disk_after=$(disk_usage_mb)
  local disk_delta=$((disk_after - disk_before))

  local per_person
  per_person=$(python3 -c "print(f'{$wall / ${#people[@]}:.1f}')")
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$(date -Iseconds)" "$phase" "$concurrency" "$threads" "$total_cores" "${#people[@]}" "$wall" "$per_person" \
    "${mem_avail_min:--}" "$mem_used_peak" "$mem_used_peak_pct" "$disk_delta" \
    >> "$LOG"
  echo "phase=$phase concurrency=$concurrency threads=$threads: wall=${wall}s (${per_person}s/person), " \
       "peak mem used ~${mem_used_peak}MB (${mem_used_peak_pct}% of ${TOTAL_MEM_MB}MB), disk delta ${disk_delta}MB " \
       "-- logged to $LOG" >&2
  checkpoint "run_config EXIT phase=$phase concurrency=$concurrency threads=$threads -- row appended to \$LOG"
}

# ============ PHASE 1: thread-vs-concurrency comparison, small fixed batch ============
checkpoint "entering Phase 1 loop (budgets: 8 32)"
SMALL=("${ALL_PEOPLE[@]:0:$SMALL_BATCH}")
if [ "${#SMALL[@]}" -lt "$SMALL_BATCH" ]; then
  echo "NOTE: only ${#SMALL[@]} people available for phase 1 (wanted $SMALL_BATCH)." >&2
fi
echo "" >&2; echo ">>> PHASE 1: thread-vs-concurrency, batch=${SMALL[*]} <<<" >&2
for budget in 8 32; do
  checkpoint "Phase 1: budget=$budget iteration starting"
  [ "$budget" -gt "$NPROC" ] && continue
  c=1
  while [ "$c" -le "$budget" ]; do
    if [ $((budget % c)) -eq 0 ]; then
      t=$((budget / c))
      run_config "1_thread_vs_concurrency" "$c" "$t" "${SMALL[@]}"
    fi
    c=$((c * 2))
  done
done
checkpoint "Phase 1 loop complete"

# ============ PHASE 2: horizontal stress test, batch size == concurrency ============
checkpoint "entering Phase 2 loop (TOP_STRESS_LEVEL=$TOP_STRESS_LEVEL)"
echo "" >&2; echo ">>> PHASE 2: horizontal stress test (threads=1, batch size = concurrency) <<<" >&2
level=4
while [ "$level" -le "$TOP_STRESS_LEVEL" ]; do
  checkpoint "Phase 2: level=$level iteration starting"
  BATCH=("${ALL_PEOPLE[@]:0:$level}")
  run_config "2_horizontal_stress" "$level" 1 "${BATCH[@]}"
  level=$((level * 2))
done
checkpoint "Phase 2 loop complete"
# Always include the true top level even if it's not a power of 2 (e.g. nproc=32 with 60 people
# supplied -> level doubling hits 32 exactly, but if nproc were e.g. 24 this catches it).
if [ "$level" -ne $((TOP_STRESS_LEVEL * 2)) ] && [ "$TOP_STRESS_LEVEL" -ge 4 ]; then
  BATCH=("${ALL_PEOPLE[@]:0:$TOP_STRESS_LEVEL}")
  run_config "2_horizontal_stress" "$TOP_STRESS_LEVEL" 1 "${BATCH[@]}"
fi

checkpoint "script reached the end normally (Done)"
echo "" >&2
echo "Done. Raw rows in $LOG, baseline profile in $BASELINE_LOG." >&2
echo "Paste both back -- the cost/curve/memory/disk analysis is a separate step once real numbers exist." >&2
echo "Per-worker output fragments (immuannot_calls.worker_<pid>.tsv etc.) are left in $OUTROOT --" >&2
echo "safe to delete once you've confirmed the run looked healthy (this experiment only needed the" >&2
echo "timing/resource numbers above, not the calls themselves)." >&2
