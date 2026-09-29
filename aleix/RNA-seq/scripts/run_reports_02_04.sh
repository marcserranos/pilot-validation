#!/bin/bash
# Reports 02 -> 03 -> 04 in order; stops at the first failure. Run from anywhere on the VM:
#   nohup bash ~/repos/pilot-validation/aleix/RNA-seq/scripts/run_reports_02_04.sh > ~/reports_02_04.log 2>&1 &
# Then copy the de-identified outputs into the repo (figures + CSVs only) -- see STATUS.md.
set -euo pipefail
cd "$(dirname "$0")/.."

PERSON=~/pipeline_outputs/rnaseq/pheno/person.tsv
[[ -s "$PERSON" ]] || { echo "FATAL: $PERSON missing -- run query_overlap_phenotypes.py in a Workbench notebook first"; exit 1; }
[[ -s ~/pipeline_outputs/rnaseq/embeddings/embeddings_sceptr_cohort_full_vcdr3.npy ]] || { echo "FATAL: SCEPTR embeddings missing -- run embed_cdr3s.py first"; exit 1; }

pixi run pip install --quiet olga
for s in 02_repertoire_atlas 03_publicness_and_robustness 04_antigen_specificity; do
  echo "==== $s :: $(date) ===="
  pixi run python3 -u "scripts/$s.py"
done
echo "==== all reports done :: $(date) ===="
