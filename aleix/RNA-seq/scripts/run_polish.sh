#!/bin/bash
# Methods polish run (2026-10-04): audit, fixes, re-embed, rerun the affected reports, push.
#
#   nohup bash ~/repos/pilot-validation/aleix/RNA-seq/scripts/run_polish.sh > ~/polish.log 2>&1 &
#
# RESUMABLE: each finished step leaves ~/polish_state/<step>.done and is skipped on a rerun,
# so after a failure fix the cause and launch the same command again. `rm -r ~/polish_state`
# to start over. Order matters: the audit (08) must see the old top-500 pool, so it runs
# before the re-embed, and the old pool + embedding are moved, not deleted.
#
# What changes and why: see context/DECISIONS.md (2026-10-04) and the methods artifact.
#   - top-500 ties broken by a seeded per-person key, not alphabetically by V gene -> re-embed
#   - every distinct full-cache (TRBV, CDR3) embedded once -> report 06 uses whole repertoires
#   - 06: block-balanced concatenation, z-scored and depth-residualised identifiability, CIs
#   - 07: person-level axis = mean pooling (+ balanced TRBV), epitope CIs, paired vs b_sceptr
#   - 04: epitope-bootstrap CIs, decoy Pgen-tolerance table;  05: V-stratified (MH) enrichment
#   - aggregation: chain from V/J/C, p5/p95 instead of single-person min/max
# Report 02 is not rerun: it uses the full cache, where tie order does not matter.
set -euo pipefail
REPO=~/repos/pilot-validation
RS=$REPO/aleix/RNA-seq
OUT=~/pipeline_outputs/rnaseq
EMB=$OUT/embeddings
STATE=~/polish_state
TAG=cohort_full_vcdr3
mkdir -p "$STATE" "$OUT/reports/08_methods_audit"
cd "$RS"

step() {
  local name="$1"; shift
  if [[ -f "$STATE/$name.done" ]]; then echo "==== $name :: already done, skipping"; return 0; fi
  echo "==== $name :: start $(date -u +%FT%TZ)"
  local t0=$SECONDS
  "$@"
  touch "$STATE/$name.done"
  echo "==== $name :: done in $(( (SECONDS - t0) / 60 )) min"
}
py() { pixi run python3 -u "$@"; }

# ---- preflight
[[ -s $OUT/pheno/person.tsv ]] || { echo "FATAL: person.tsv missing"; exit 1; }
[[ -s $OUT/cohort_full.tsv ]] || { echo "FATAL: cohort_full.tsv missing"; exit 1; }
[[ -s $OUT/atlas/trb_clonotypes_full.pkl ]] || { echo "FATAL: 02's clonotype cache missing"; exit 1; }
[[ -s $EMB/embeddings_sceptr_$TAG.npy || -s $EMB/v1_alphabetical_ties/embeddings_sceptr_$TAG.npy ]] \
  || { echo "FATAL: no SCEPTR embedding found"; exit 1; }
git -C "$REPO" pull --ff-only

step deps          pixi run pip install --quiet olga scikit-learn umap-learn cnsplots
step freeze_env    bash -c "bash scripts/freeze_env.sh > $OUT/reports/08_methods_audit/environment.txt"
step audit_08      py scripts/08_methods_audit.py
step aggregate     py scripts/aggregate_rnaseq_results.py "$OUT/cohort_full.tsv"

backup_v1() {
  mkdir -p "$EMB/v1_alphabetical_ties"
  for f in embeddings_sceptr_$TAG.npy pool_sceptr_$TAG.tsv; do
    if [[ -e $EMB/$f && ! -e $EMB/v1_alphabetical_ties/$f ]]; then
      mv "$EMB/$f" "$EMB/v1_alphabetical_ties/"
    fi
  done
  ls -la "$EMB/v1_alphabetical_ties"
}
step backup_v1     backup_v1
step embed_pool    py scripts/embed_cdr3s.py "$OUT/cohort_full.tsv" --models sceptr
step embed_full    py scripts/embed_full_cache.py
step backup_v2     gcloud storage cp "$EMB/embeddings_sceptr_$TAG.npy" "$EMB/pool_sceptr_$TAG.tsv" \
                     "$EMB/embeddings_sceptr_fullcache.npy" "$EMB/keys_sceptr_fullcache.tsv" \
                     gs://aleix-rnaseq-wb-cordial-leechee-9743/embeddings/v2_random_ties/

step report_01     py scripts/01_sceptr_embedding_viz.py
step report_03     py scripts/03_publicness_and_robustness.py
step report_04     py scripts/04_antigen_specificity.py
step report_05     py scripts/05_embedding_clusters.py
step report_06     py scripts/06_person_representation.py
step report_07     py scripts/07_embedding_settings.py

# ---- copy aggregate outputs only (top-level figures, CSVs, environment.txt; never _cache/)
publish() {
  for d in 01_sceptr_embedding_viz 03_publicness_and_robustness 04_antigen_specificity \
           05_embedding_clusters 06_person_representation 07_embedding_settings 08_methods_audit; do
    mkdir -p "$RS/reports/$d"
    find "$OUT/reports/$d" -maxdepth 1 -type f \
         \( -name '*.csv' -o -name '*.pdf' -o -name '*.svg' -o -name '*.png' -o -name 'environment.txt' \) \
         -exec cp {} "$RS/reports/$d/" \;
  done
  cd "$REPO"
  git add aleix/RNA-seq/reports aleix/RNA-seq/results
  git commit -q -m "Methods polish: audit (08), random tie-break re-embed, reports 01/03-07 rerun

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" || echo "nothing new to commit"
  git push || { echo "!! push failed -- outputs are committed locally; run: git -C $REPO push"; return 1; }
}
step publish       publish
echo "==== polish run complete :: $(date -u +%FT%TZ)"
