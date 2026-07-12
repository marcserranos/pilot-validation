#!/bin/bash
# Experiment C, minimap2 variant -- run this AFTER run_experiment_c_comparison.sh's bwa
# job has fully finished, not alongside it. This is a timing experiment; running both
# concurrently on the same 4-vCPU VM would contend for cores and contaminate the exact
# measurement being taken. Reuses the pad100k_minimap2 baseline (2026-07-09 LR sweep,
# ~14.5min vs bwa's ~21min) so this isolates the exact same single variable (restricted
# DB), just against a different, faster aligner baseline. Separate output directory and
# run-label from the bwa script, so results don't collide, but sequencing (not
# concurrency) is what keeps the timing clean.
#
# Caveat (2026-07-09, EXPERIMENTS.md/DECISIONS.md): minimap2 is NOT yet the confirmed
# production default -- the earlier suspected DPB1 regression (person 1017156, n=1) did
# not reproduce on person 2522883 (n=2 overall), but that's still not conclusive. Treat
# this run's result as directional, a second data point, not the answer -- same caveat
# the bwa run doesn't need since bwa is still the standing default.
#
# Endpoint this is feeding toward, if the data ends up backing it: fully deprecating bwa
# for SpecImmune-LR and standardizing on minimap2 (same accuracy, substantially faster) --
# not decided yet, this run (plus the still-open n=3 confirmation) is what that decision
# would rest on.
#
# Prereqs: same three as run_experiment_c_comparison.sh, plus the pad100k_minimap2
# baseline (not just pad100k_bwa) must already exist for this person, and the bwa
# comparison must have already completed (check its progress.log first).
#
# Usage: bash run_experiment_c_comparison_minimap2.sh <person_id>

set -uo pipefail

PERSON_ID="${1:?Usage: bash run_experiment_c_comparison_minimap2.sh <person_id>}"
OUTDIR="$HOME/pipeline_outputs/$PERSON_ID"
SPECIMMUNE_DIR="$HOME/tools/SpecImmune"
PIXI_MANIFEST="$HOME/repos/pilot-validation/pixi.toml"
COMPARE="$HOME/repos/pilot-validation/scripts/compare_hla_results.py"
RESTRICTED_DB="$SPECIMMUNE_DIR/db_classical8"

BASELINE_FASTQ="$OUTDIR/sweep/pad100k_minimap2/LR.fastq"
BASELINE_RESULT="$OUTDIR/sweep/pad100k_minimap2/specimmune_output/$PERSON_ID/$PERSON_ID.HLA.final.type.result.formatted.txt"

RUNDIR="$OUTDIR/experiment_c_restricted_db_minimap2"
PROGRESS="$RUNDIR/progress.log"
mkdir -p "$RUNDIR"

log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$PROGRESS"; }

if [ ! -f "$BASELINE_FASTQ" ]; then
  log "FATAL: baseline FASTQ not found at $BASELINE_FASTQ -- run the LR padding sweep for $PERSON_ID first."
  exit 1
fi
if [ ! -f "$BASELINE_RESULT" ]; then
  log "FATAL: baseline result not found at $BASELINE_RESULT -- the pad100k_minimap2 config didn't complete for $PERSON_ID."
  exit 1
fi
if [ ! -d "$RESTRICTED_DB/HLA" ]; then
  log "FATAL: restricted DB not found at $RESTRICTED_DB -- run build_restricted_specimmune_db.sh first."
  exit 1
fi

log "=== Experiment C (minimap2 variant): restricted-DB run for $PERSON_ID (same FASTQ as pad100k_minimap2 baseline) ==="

specout="$RUNDIR/specimmune_output"
timing="$RUNDIR/timing.txt"
result_file="$specout/$PERSON_ID/$PERSON_ID.HLA.final.type.result.formatted.txt"

( cd "$SPECIMMUNE_DIR" && \
  { time pixi run --manifest-path "$PIXI_MANIFEST" -e specimmune -- \
      python3 scripts/main.py -n "$PERSON_ID" -o "$specout" -j 4 -y pacbio-hifi \
      -i HLA -r "$BASELINE_FASTQ" --db "$RESTRICTED_DB" --align_method_1 minimap2 --visualization "" ; } \
    2> "$timing" )
rc=$?
if [ $rc -ne 0 ]; then
  log "FAILED (exit $rc) -- see $timing"
  exit 1
fi
if [ ! -f "$result_file" ]; then
  log "FAILED -- exit 0 but expected output missing at $result_file -- see $timing"
  exit 1
fi
log "Restricted-DB run done ($(grep real "$timing" | tail -1))"

# Completeness check (same lesson as Experiment B): count non-empty gene calls, don't just
# trust file existence.
completeness=$(python3 -c "
import re
genes = ['A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1']
seen = {g: 0 for g in genes}
with open('$result_file') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 3:
            continue
        locus = parts[0].replace('HLA-', '')
        if locus in seen and parts[2] not in ('NA', '-', ''):
            seen[locus] += 1
called = sum(1 for g in genes if seen[g] > 0)
print(f'{called}/8 classical genes have at least one haplotype called: {seen}')
")
log "Completeness: $completeness"

log "individual_ref/ genes actually processed:"
find "$specout/$PERSON_ID/individual_ref" -maxdepth 1 -mindepth 1 -type d 2>>"$PROGRESS" | xargs -n1 basename | sort | tee -a "$PROGRESS"

pixi run --manifest-path "$PIXI_MANIFEST" -e specimmune -- python3 "$COMPARE" "$PERSON_ID" \
  --run-label "experiment_c_restricted_db_minimap2" \
  --specimmune-result "$result_file" \
  >> "$RUNDIR/comparison.md" 2>>"$PROGRESS"

log "Baseline (full-panel DB, pad100k_minimap2, same FASTQ) result for comparison: $BASELINE_RESULT"
log "=== Experiment C (minimap2 variant) comparison complete -- see $RUNDIR/comparison.md, $OUTDIR/comparison_log.csv, and $PROGRESS ==="
