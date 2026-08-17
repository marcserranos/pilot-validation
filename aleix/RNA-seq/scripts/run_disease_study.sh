#!/usr/bin/env bash
# Unattended driver for the LR x RNA-seq disease study (Task 2 steps 2-3).
#
# WHY THIS EXISTS: query_overlap_phenotypes.py's docstring says "notebook only", because
# ENVIRONMENT.md quirk #5 records that the CDR is not exposed via shell env vars. But that
# only means the dataset id must be passed EXPLICITLY -- BigQuery itself is reachable from
# the shell (quirk #4 notes a real query returns real rows). So with --cdr supplied and a
# python that has the BigQuery client, both steps run headless. That makes the whole study
# a fire-and-forget job instead of an interactive session.
#
# Picks a working python automatically rather than assuming: the Workbench system python
# ships google-cloud-bigquery / pandas-gbq, the pixi env probably does not. Tries each,
# reports which it used, fails loudly with the real import error if neither works.
#
# Usage:
#   bash run_disease_study.sh [CDR_DATASET]
#   nohup bash run_disease_study.sh > ~/pipeline_outputs/rnaseq/disease_study.log 2>&1 &
#
# CDR_DATASET defaults to the value recorded in context/ENVIRONMENT.md. That was verified
# in July 2026 and the id changes with every CDR release -- if the query 404s on the
# dataset, that is the first thing to re-check, not a bug in the SQL.

set -uo pipefail

CDR="${1:-${WORKSPACE_CDR:-wb-silky-artichoke-2408.C2025Q4R6}}"
REPO="$HOME/repos/pilot-validation"
SCRIPTS="$REPO/aleix/RNA-seq/scripts"
COHORT="$HOME/pipeline_outputs/rnaseq/lr_rnaseq_overlap_cohort.tsv"
PHENO_DIR="$HOME/pipeline_outputs/rnaseq/pheno"
PIXI_PY="pixi run --manifest-path $REPO/aleix/RNA-seq/pixi.toml python3"

echo "=== disease study started $(date) ==="
echo "CDR dataset : $CDR"
echo "cohort      : $COHORT"

if [[ ! -f "$COHORT" ]]; then
  echo "FATAL: cohort not found: $COHORT"
  echo "Run check_lr_rnaseq_overlap.py first."
  exit 1
fi
echo "cohort rows : $(( $(wc -l < "$COHORT") - 1 ))"

# ---- pick a python that can talk to BigQuery -------------------------------
PY=""
for cand in "python3" "/opt/conda/bin/python3"; do
  if command -v "${cand%% *}" >/dev/null 2>&1 && \
     $cand -c "import google.cloud.bigquery" >/dev/null 2>&1; then
    PY="$cand"; break
  fi
done
if [[ -z "$PY" ]]; then
  echo "--- no system python with google-cloud-bigquery; trying pixi env ---"
  if $PIXI_PY -c "import google.cloud.bigquery" >/dev/null 2>&1; then
    PY="$PIXI_PY"
  fi
fi
if [[ -z "$PY" ]]; then
  echo "FATAL: no python with the BigQuery client. Real import error follows:"
  python3 -c "import google.cloud.bigquery" || true
  echo
  echo "Fix: run step 1 in a notebook instead, OR add google-cloud-bigquery to"
  echo "aleix/RNA-seq/pixi.toml and re-run 'pixi install'."
  exit 1
fi
echo "python      : $PY"

# ---- step 1: pull phenotypes (BigQuery) ------------------------------------
echo
echo "=== [1/2] pulling demographics + EHR window + conditions $(date) ==="
$PY "$SCRIPTS/query_overlap_phenotypes.py" --cohort "$COHORT" --cdr "$CDR" \
    --outdir "$PHENO_DIR"
rc=$?
if [[ $rc -ne 0 ]]; then
  echo "FATAL: phenotype query failed (exit $rc). Most likely causes, in order:"
  echo "  1. stale CDR dataset id -- re-read it from a Dataset Builder snippet"
  echo "  2. no billing project -- export GOOGLE_PROJECT=<your workspace project>"
  echo "  3. controlled-tier permissions"
  exit $rc
fi

# ---- step 2: analyse (pure pandas, needs scipy-free pixi env is fine) ------
echo
echo "=== [2/2] disease burden + rankings + immune feasibility $(date) ==="
$PIXI_PY "$SCRIPTS/analyze_disease_burden.py" --indir "$PHENO_DIR" --cohort "$COHORT"
rc=$?
if [[ $rc -ne 0 ]]; then
  echo "FATAL: analysis failed (exit $rc)."
  exit $rc
fi

echo
echo "=== disease study finished $(date) ==="
echo "De-identified outputs are in $REPO/aleix/RNA-seq/results/ -- safe to commit."
echo "Per-person clinical data is in $PHENO_DIR -- VM-LOCAL, do NOT commit."
