#!/usr/bin/env bash
# S05 WS-F -- 61_kir_by_ancestry_panel.py on the VM (PUBLIC aggregates only).
# Paste block by block into the JupyterLab terminal (type and Return as separate actions).
# Reads ~/pipeline_outputs_kir/<pid>/immuannot_output/hap{1,2}.gtf.gz, cohort_membership.tsv and
# the AoU relatedness table; writes ~/s05/results/61/. No INTERNAL output.
set -euo pipefail
REPO=${REPO:-~/repos/pilot-validation-main}   # the clone 58 runs from; ~/repos/pilot-validation on older apps
cd "$REPO" && git pull --ff-only
mkdir -p ~/s05/results ~/s05/logs

# 0. Pre-flight: KIR outputs present, automount up for the relatedness table.
ls ~/pipeline_outputs_kir | head -3
ls ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv

# 1. Smoke test (own output path, quirk #22): 200 people, tables + figure.
python3 scripts/hla_popgen/61_kir_by_ancestry_panel.py --mode vm --limit 200 --workers 4 \
    --out-dir ~/s05/results/61_smoke 2>&1 | tail -20

# 2. Full run WITHOUT the coverage filter (reference run; also usable before 59 has run). ~3-5 min at 8 workers (51 took ~3 min).
python3 scripts/hla_popgen/61_kir_by_ancestry_panel.py --mode vm --workers 8 \
    --out-dir ~/s05/results/61_nofilter > ~/s05/logs/61_nofilter.log 2>&1
tail -5 ~/s05/logs/61_nofilter.log

# 3. With the coverage filter: needs 59's INTERNAL exclude_people.tsv (59 must have run; 61 drops
#    everyone with exclude=True, i.e. below 59's region proposal in window_MHC or window_KIR by default).
ls -la ~/s05/internal/59/exclude_people.tsv
python3 scripts/hla_popgen/61_kir_by_ancestry_panel.py --mode vm --workers 8 \
    --exclude-people ~/s05/internal/59/exclude_people.tsv \
    --out-dir ~/s05/results/61 > ~/s05/logs/61.log 2>&1
tail -5 ~/s05/logs/61.log

# If matplotlib is missing in this python, add --no-figure and draw locally after pulling:
#   python3 scripts/hla_popgen/61_kir_by_ancestry_panel.py --mode figure --tables-dir <pulled dir>

# 4. Pull: ~/s05/results/61*/ are PUBLIC (aggregate, masked). Check for '<20' and no person ids:
grep -l -E '[0-9]{7}' ~/s05/results/61_nofilter/*.tsv || echo "no 7-digit ids in tables"
