#!/bin/bash
# Experiment C: SpecImmune gene-panel restriction, controlled timing + accuracy comparison.
# Isolates ONE variable -- the database -- by reusing the exact same FASTQ as the existing
# pad100k-bwa baseline from the 2026-07-09 LR sweep (already the recommended new default
# padding). See EXPERIMENTS.md Experiment C for the full investigation.
#
# Prereqs:
#   1. patch_specimmune_for_gene_restriction.py already run against ~/tools/SpecImmune
#   2. build_restricted_specimmune_db.sh already run (produces ~/tools/SpecImmune/db_classical8)
#   3. The pad100k_bwa baseline from the LR sweep must already exist for this person
#
# Usage: bash run_experiment_c_comparison.sh <person_id>

set -uo pipefail

PERSON_ID="${1:?Usage: bash run_experiment_c_comparison.sh <person_id>}"
OUTDIR="$HOME/pipeline_outputs/$PERSON_ID"
SPECIMMUNE_DIR="$HOME/tools/SpecImmune"
PIXI_MANIFEST="$HOME/repos/pilot-validation/pixi.toml"
COMPARE="$HOME/repos/pilot-validation/compare_hla_results.py"
RESTRICTED_DB="$SPECIMMUNE_DIR/db_classical8"

BASELINE_FASTQ="$OUTDIR/sweep/pad100k_bwa/LR.fastq"
BASELINE_RESULT="$OUTDIR/sweep/pad100k_bwa/specimmune_output/$PERSON_ID/$PERSON_ID.HLA.final.type.result.formatted.txt"

RUNDIR="$OUTDIR/experiment_c_restricted_db"
PROGRESS="$RUNDIR/progress.log"
mkdir -p "$RUNDIR"

log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$PROGRESS"; }

if [ ! -f "$BASELINE_FASTQ" ]; then
  log "FATAL: baseline FASTQ not found at $BASELINE_FASTQ -- run the LR padding sweep for $PERSON_ID first."
  exit 1
fi
if [ ! -f "$BASELINE_RESULT" ]; then
  log "FATAL: baseline result not found at $BASELINE_RESULT -- the pad100k_bwa config didn't complete for $PERSON_ID."
  exit 1
fi
if [ ! -d "$RESTRICTED_DB/HLA" ]; then
  log "FATAL: restricted DB not found at $RESTRICTED_DB -- run build_restricted_specimmune_db.sh first."
  exit 1
fi

log "=== Experiment C: restricted-DB run for $PERSON_ID (same FASTQ as pad100k_bwa baseline) ==="

specout="$RUNDIR/specimmune_output"
timing="$RUNDIR/timing.txt"
result_file="$specout/$PERSON_ID/$PERSON_ID.HLA.final.type.result.formatted.txt"

( cd "$SPECIMMUNE_DIR" && \
  { time pixi run --manifest-path "$PIXI_MANIFEST" -e specimmune -- \
      python3 scripts/main.py -n "$PERSON_ID" -o "$specout" -j 4 -y pacbio-hifi \
      -i HLA -r "$BASELINE_FASTQ" --db "$RESTRICTED_DB" --align_method_1 bwa --visualization "" ; } \
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

# Also confirm no non-classical gene directories leaked into individual_ref/ -- the
# specific open question from the roadmap.
log "individual_ref/ genes actually processed:"
find "$specout/$PERSON_ID/individual_ref" -maxdepth 1 -mindepth 1 -type d 2>>"$PROGRESS" | xargs -n1 basename | sort | tee -a "$PROGRESS"

pixi run --manifest-path "$PIXI_MANIFEST" -e specimmune -- python3 "$COMPARE" "$PERSON_ID" \
  --run-label "experiment_c_restricted_db" \
  --specimmune-result "$result_file" \
  >> "$RUNDIR/comparison.md" 2>>"$PROGRESS"

log "Baseline (full-panel DB, same FASTQ) result for comparison: $BASELINE_RESULT"
log "=== Experiment C comparison complete -- see $RUNDIR/comparison.md, $OUTDIR/comparison_log.csv, and $PROGRESS ==="
