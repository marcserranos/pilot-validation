#!/bin/bash
# Final figure pass: journal font on the VM, then re-render the figures that need
# person-level data (reports 02, 09, 10, 13). Report 14's figure is drawn from its CSV and
# was already re-rendered off the VM. Results are deterministic (seed 0), so the CSVs come
# out identical; only the figures change. ~10 min. Commits and pushes.
#
#   bash ~/repos/pilot-validation/aleix/RNA-seq/scripts/run_figs.sh 2>&1 | tee ~/figs.log
#
# Font: the VM has neither Helvetica nor Arial, so matplotlib fell back to DejaVu Sans.
# Liberation Sans is metric-compatible with Arial (free, SIL OFL) and is in cnsplots' font
# list. Installed via apt if passwordless sudo exists, otherwise from the official release
# (github.com/liberationfonts/liberation-fonts, 2.1.5, TTF archive ~2 MB) into ~/.fonts.
set -euo pipefail
REPO=~/repos/pilot-validation
RS=$REPO/aleix/RNA-seq
OUT=~/pipeline_outputs/rnaseq
URL=https://github.com/liberationfonts/liberation-fonts/files/7261482/liberation-fonts-ttf-2.1.5.tar.gz
git -C "$REPO" pull --ff-only

if fc-list 2>/dev/null | grep -qi "Liberation Sans"; then
  echo "Liberation Sans already installed"
elif sudo -n true 2>/dev/null && sudo -n apt-get install -y -q fonts-liberation >/dev/null 2>&1; then
  echo "Liberation Sans installed via apt"
else
  mkdir -p ~/.fonts /tmp/libfonts
  curl -fsSL -o /tmp/libfonts/lib.tar.gz "$URL"
  tar -xzf /tmp/libfonts/lib.tar.gz -C /tmp/libfonts
  cp /tmp/libfonts/liberation-fonts-ttf-*/LiberationSans-*.ttf ~/.fonts/
  fc-cache -f ~/.fonts >/dev/null 2>&1 || true
  echo "Liberation Sans installed into ~/.fonts"
fi
rm -rf ~/.cache/matplotlib
cd "$RS"
pixi run python3 - <<'PY'
import matplotlib.font_manager as fm
try:
    path = fm.findfont(fm.FontProperties(family="Liberation Sans"), fallback_to_default=False)
except ValueError:
    raise SystemExit("FATAL: matplotlib cannot see Liberation Sans; stopping before re-rendering")
print("matplotlib finds Liberation Sans:", path)
PY
for s in 02_repertoire_atlas 09_analysis_table 10_vj_usage_ancestry 13_disease_association; do
  echo "==== $s :: $(date -u +%T)"
  pixi run python3 -u "scripts/$s.py" 2>&1 | grep -v -i "warn\|findfont" | tail -2
  mkdir -p "reports/$s"
  find "$OUT/reports/$s" -maxdepth 1 -type f \( -name '*.csv' -o -name '*.pdf' -o -name '*.svg' -o -name '*.png' \) \
       -exec cp {} "reports/$s/" \;
done

cd "$REPO"
git add aleix/RNA-seq/reports
git commit -q -m "Figures 02, 09, 10, 13 re-rendered in Liberation Sans (Arial metrics); 10 titles fixed

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" || echo "nothing new to commit"
git pull --rebase -q && git push -q && echo "pushed: $(git log --oneline -1)"
