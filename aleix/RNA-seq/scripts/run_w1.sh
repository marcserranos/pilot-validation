#!/bin/bash
# W1: V/J gene usage by ancestry (report 10). Re-reads every cdr3.out to get J genes (the
# clonotype cache keeps only V), ~5-10 min. Commits and pushes the aggregate outputs.
#
#   bash ~/repos/pilot-validation/aleix/RNA-seq/scripts/run_w1.sh 2>&1 | tee ~/w1.log
set -euo pipefail
REPO=~/repos/pilot-validation
RS=$REPO/aleix/RNA-seq
OUT=~/pipeline_outputs/rnaseq
git -C "$REPO" pull --ff-only
[[ -s $OUT/analysis/person_table.pkl ]] || { echo "FATAL: run run_w0.sh first"; exit 1; }
cd "$RS"
pixi run python3 -u scripts/10_vj_usage_ancestry.py 2>&1 | grep -v '^  read '
D=reports/10_vj_usage_ancestry
mkdir -p "$D"
find "$OUT/$D" -maxdepth 1 -type f \( -name '*.csv' -o -name '*.pdf' -o -name '*.svg' -o -name '*.png' \) \
     -exec cp {} "$D/" \;
cd "$REPO"
git add "aleix/RNA-seq/$D"
git commit -q -m "10 V/J usage by ancestry (W1)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" || echo "nothing new to commit"
git pull --rebase -q && git push -q && echo "pushed: $(git log --oneline -1)"
