#!/bin/bash
# Core-scaling experiment: how does wall-clock for a fixed batch of people change as available
# compute is split between "more people at once" vs "more threads per person"? Marc, 2026-07-24
# weekly objective: "distributed 2 cores, 4 cores and 8 cores... see how the curve progresses
# (extrapolate from there)". Uses Immuannot (the fastest current LR caller, median 6.9 min/person
# on this cohort per reports/immuannot_pilot/README.md) so a full sweep is cheap to run.
#
# IMPORTANT SCOPE NOTE: run_immuannot_person.py processes its person_ids list SERIALLY (one
# person at a time -- confirmed by reading its main(), no cross-person parallelism exists in the
# tool itself). This script adds that missing axis via `xargs -P` (N separate OS processes, each
# a single-person run_immuannot_person.py invocation), which is the only way to actually exercise
# more than 1 core concurrently with this pipeline as it exists today.
#
# THIS SCRIPT ONLY MEASURES CONFIGURATIONS THAT FIT ON THE CURRENT VM'S ACTUAL CORE COUNT
# (`nproc`). As of 2026-07-21 the Workbench VM is 4 vCPU (context/ENVIRONMENT.md). To get a REAL
# (not oversubscribed/misleading) measurement at an 8-core total budget, resize the Workbench VM
# to an 8-vCPU machine type first (a Workbench UI action -- Marc's to drive) and rerun this same
# script there; don't try to fake it by running concurrency=8 on a 4-vCPU box, that measures
# contention, not scaling.
#
# For each (concurrency, threads_per_person) pair, concurrency * threads_per_person = the total
# core budget for that row -- deliberately swept both ways (e.g. 4 total cores as
# concurrency=1/threads=4, concurrency=2/threads=2, concurrency=4/threads=1) so the summary can
# show whether this tool's OWN internal threading scales as well as running separate people does.
#
# Usage (from ~/repos/pilot-validation, inside `pixi shell -e specimmune` -- same env
# run_immuannot_person.py itself needs):
#   bash scripts/run_core_scaling_experiment.sh <person_id_1> <person_id_2> ... <person_id_8>
#
# Pick 8 people who are NOT already in immuannot_calls.tsv (or pass --force via IMMUANNOT_EXTRA_ARGS
# below) so every row does real fresh work, not a skip. Needs 8 distinct people total: this sweep
# uses BATCH_SIZE=4 people per row by default (adjust below) so each row's timing reflects a
# consistent batch, not a growing "already done" skip count from a prior row.
#
# Output: appends one line per (concurrency, threads) config to
#   ~/pipeline_outputs/core_scaling_experiment.tsv
# columns: timestamp, concurrency, threads_per_person, total_cores_used, batch_size,
#          wall_clock_seconds, seconds_per_person_wallclock
# Aggregate-only (timing numbers, no genotypes) -- fine to bring this TSV back off the VM.

set -uo pipefail

OUTROOT="$HOME/pipeline_outputs"
LOG="$OUTROOT/core_scaling_experiment.tsv"
BATCH_SIZE=4          # people processed per config row -- keep small, this is a timing probe
IMMUANNOT_EXTRA_ARGS="${IMMUANNOT_EXTRA_ARGS:---force}"   # --force: redo people even if already called,
                                                            # so reruns of this script stay comparable

if [ "$#" -lt "$BATCH_SIZE" ]; then
  echo "FATAL: need at least $BATCH_SIZE person_ids (one batch's worth), got $#." >&2
  echo "Usage: bash scripts/run_core_scaling_experiment.sh <person_id_1> ... <person_id_N>" >&2
  exit 1
fi

NPROC=$(nproc)
echo "This VM has $NPROC vCPUs (nproc). Configs whose concurrency*threads > $NPROC will be SKIPPED" >&2
echo "as oversubscribed/misleading -- resize the VM and rerun for those instead." >&2

mkdir -p "$OUTROOT"
if [ ! -f "$LOG" ]; then
  printf "timestamp\tconcurrency\tthreads_per_person\ttotal_cores_used\tbatch_size\twall_clock_seconds\tseconds_per_person_wallclock\n" > "$LOG"
fi

# people[0..BATCH_SIZE-1] used for every config row (same batch each time for comparability)
PEOPLE=("${@:1:$BATCH_SIZE}")

run_config() {
  local concurrency="$1" threads="$2"
  local total_cores=$((concurrency * threads))
  if [ "$total_cores" -gt "$NPROC" ]; then
    echo "SKIP concurrency=$concurrency threads=$threads (total_cores=$total_cores > nproc=$NPROC)" >&2
    return
  fi
  echo "=== concurrency=$concurrency threads=$threads (total_cores=$total_cores), batch=${PEOPLE[*]} ===" >&2
  local t0 t1 wall
  t0=$(date +%s)
  printf "%s\n" "${PEOPLE[@]}" | \
    xargs -P "$concurrency" -I{} pixi run -e specimmune -- \
      python3 scripts/run_immuannot_person.py {} --threads "$threads" $IMMUANNOT_EXTRA_ARGS
  t1=$(date +%s)
  wall=$((t1 - t0))
  local per_person
  per_person=$(python3 -c "print(f'{$wall / ${#PEOPLE[@]}:.1f}')")
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$(date -Iseconds)" "$concurrency" "$threads" "$total_cores" "${#PEOPLE[@]}" "$wall" "$per_person" \
    >> "$LOG"
  echo "concurrency=$concurrency threads=$threads: wall=${wall}s (${per_person}s/person) -- logged to $LOG" >&2
}

# Sweep: same total-core budgets (2, 4, 8), each split two ways, so the summary can compare
# "N people at once, 1 thread each" against "fewer people at once, more threads each".
run_config 1 2   # 2 cores: 1 person, 2 threads
run_config 2 1   # 2 cores: 2 people at once, 1 thread each

run_config 1 4   # 4 cores: 1 person, 4 threads (today's actual default config)
run_config 2 2   # 4 cores: 2 people at once, 2 threads each
run_config 4 1   # 4 cores: 4 people at once, 1 thread each

run_config 1 8   # 8 cores: 1 person, 8 threads
run_config 2 4   # 8 cores: 2 people at once, 4 threads each
run_config 4 2   # 8 cores: 4 people at once, 2 threads each
run_config 8 1   # 8 cores: 8 people at once, 1 thread each

echo "" >&2
echo "Done. Raw rows in $LOG -- any row with total_cores_used > $NPROC was skipped this run." >&2
echo "If 8-core rows were skipped: resize the Workbench VM to 8 vCPU and rerun this script there" >&2
echo "(same command -- already-done people are safely --force-redone, so results stay comparable)." >&2
echo "Paste the resulting TSV back; the cost/curve analysis is a separate step once real numbers exist." >&2
