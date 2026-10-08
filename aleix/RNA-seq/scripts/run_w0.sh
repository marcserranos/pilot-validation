#!/bin/bash
# W0: fetch the three small AoU tables, build the per-person analysis table (report 09),
# commit and push the aggregate outputs. About 10-20 min, almost all of it the fixed-depth
# subsampling over 7,922 repertoires.
#
#   bash ~/repos/pilot-validation/aleix/RNA-seq/scripts/run_w0.sh 2>&1 | tee ~/w0.log
#
# The person table itself stays VM-local (~/pipeline_outputs/rnaseq/analysis/); only the
# aggregate CSVs and the figure are committed.
set -euo pipefail
REPO=~/repos/pilot-validation
RS=$REPO/aleix/RNA-seq
OUT=~/pipeline_outputs/rnaseq
BUCKET=gs://vwb-aou-datasets-controlled/v9
LOCAL=~/aou_local/v9
BP=(--billing-project "$GOOGLE_PROJECT")

git -C "$REPO" pull --ff-only

fetch() {   # fetch <bucket-relative path> -- skipped if already present
  local rel="$1"
  if [[ -s "$LOCAL/$rel" ]]; then echo "have $rel"; return 0; fi
  mkdir -p "$(dirname "$LOCAL/$rel")"
  gcloud storage cp "${BP[@]}" "$BUCKET/$rel" "$LOCAL/$rel"
}
fetch wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv
fetch wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv
fetch multiomics/rnaseq/rnaseq_metadata.tsv
fetch multiomics/rnaseq/manifest.tsv

# header check only (column names, no values)
echo "ancestry_preds columns:   $(head -1 "$LOCAL/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv" | tr '\t' ' ')"
echo "relatedness columns:      $(head -1 "$LOCAL/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv" | tr '\t' ' ')"
echo "rnaseq_metadata columns:  $(head -1 "$LOCAL/multiomics/rnaseq/rnaseq_metadata.tsv" | tr '\t' ' ')"

cd "$RS"
pixi run python3 -u scripts/09_analysis_table.py

D=reports/09_analysis_table
mkdir -p "$D"
find "$OUT/$D" -maxdepth 1 -type f \( -name '*.csv' -o -name '*.pdf' -o -name '*.svg' -o -name '*.png' \) \
     -exec cp {} "$D/" \;
cd "$REPO"
git add "aleix/RNA-seq/$D"
git commit -q -m "09 analysis table (W0): waterfall, column summaries, depth figure

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" || echo "nothing new to commit"
git pull --rebase -q && git push -q && echo "pushed: $(git log --oneline -1)"
