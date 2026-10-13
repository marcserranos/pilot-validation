#!/bin/bash
# W2-W5 in one go: disease labels (11) -> feature sets (12) -> association (13) ->
# classification (14). Resumable: finished steps leave ~/disease_state/<step>.done.
# About 45-90 min, most of it W3 (3-mers and SCEPTR means over 9.4 M clonotypes) and W5.
#
#   nohup bash ~/repos/pilot-validation/aleix/RNA-seq/scripts/run_disease.sh > ~/disease.log 2>&1 &
#   progress:  grep '^====' ~/disease.log
set -euo pipefail
REPO=~/repos/pilot-validation
RS=$REPO/aleix/RNA-seq
OUT=~/pipeline_outputs/rnaseq
STATE=~/disease_state
mkdir -p "$STATE"
git -C "$REPO" pull --ff-only
for f in analysis/person_table.pkl analysis/vj_usage.pkl pheno/conditions.tsv \
         embeddings/embeddings_sceptr_fullcache.npy atlas/trb_clonotypes_full.pkl; do
  [[ -s $OUT/$f ]] || { echo "FATAL: missing $OUT/$f"; exit 1; }
done
cd "$RS"
step() {
  local name="$1"; shift
  if [[ -f "$STATE/$name.done" ]]; then echo "==== $name :: already done"; return 0; fi
  echo "==== $name :: start $(date -u +%T)"; local t0=$SECONDS
  "$@"; touch "$STATE/$name.done"
  echo "==== $name :: done in $(( (SECONDS - t0) / 60 )) min"
}
step deps     pixi run pip install --quiet xgboost
step w2_labels      pixi run python3 -u scripts/11_disease_labels.py
step w3_features    pixi run python3 -u scripts/12_feature_sets.py
step w4_association pixi run python3 -u scripts/13_disease_association.py
step w5_classify    pixi run python3 -u scripts/14_disease_classification.py

publish() {
  for d in 11_disease_labels 13_disease_association 14_disease_classification; do
    mkdir -p "$RS/reports/$d"
    find "$OUT/reports/$d" -maxdepth 1 -type f \( -name '*.csv' -o -name '*.pdf' -o -name '*.svg' -o -name '*.png' \) \
         -exec cp {} "$RS/reports/$d/" \;
  done
  cd "$REPO"
  git add aleix/RNA-seq/reports
  git commit -q -m "Disease labels (11), association (13), classification (14)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" || echo "nothing new to commit"
  git pull --rebase -q && git push -q && echo "pushed: $(git log --oneline -1)"
}
step publish publish
echo "==== disease run complete :: $(date -u +%T)"
