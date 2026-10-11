#!/usr/bin/env bash
# S05 WS-G -- 62_cole_haplotype_units.py on the VM.
# INTERNAL per-person table -> ~/s05/internal/62/ (stays on the VM / in-perimeter bucket; never pulled).
# PUBLIC aggregates + data dictionary -> ~/s05/results/62/ (safe to pull).
set -euo pipefail
REPO=${REPO:-~/repos/pilot-validation-main}
cd "$REPO" && git pull --ff-only
mkdir -p ~/s05/results ~/s05/internal ~/s05/logs

# 0. Pre-flight.
ls -la ~/pipeline_outputs/hla_calls_rich.tsv ~/pipeline_outputs/cohort_membership.tsv
ls -la ~/s05/internal/59/exclude_people.tsv || echo "59 exclusion file not there yet -> use the 58 fallback (drop --exclude-people)"
ls -la ~/s05/internal/58/per_person_depth.tsv || echo "58 depth table not there yet -> fallback coverage_flag = not_measured"

# 1. Smoke test, own paths (quirk #22).
python3 scripts/hla_popgen/62_cole_haplotype_units.py --limit 200 \
    --internal-dir ~/s05/internal/62_smoke --out-dir ~/s05/results/62_smoke 2>&1 | tail -10
zcat ~/s05/internal/62_smoke/haplotype_units.tsv.gz | head -2   # INTERNAL line + column names only, no data rows

# 2. Full run with the coverage flag set by 59 (needs 58 finished, then 59 run; 59 writes
#    ~/s05/internal/59/exclude_people.tsv). Before 59 exists, drop --exclude-people: 62 then falls back
#    to 58's depth normalised per platform, cutoff 0.6, status provisional_until_59.
python3 scripts/hla_popgen/62_cole_haplotype_units.py \
    --exclude-people ~/s05/internal/59/exclude_people.tsv \
    > ~/s05/logs/62.log 2>&1
tail -5 ~/s05/logs/62.log
cat ~/s05/internal/62/MANIFEST.tsv

# 3. Upload to the in-perimeter collaborator bucket: NOT run by the script or by an agent.
#    Marc reviews the MANIFEST (file names, bytes, rows) and runs these himself:
# gsutil -m rsync -n -r ~/s05/internal/62/ gs://hla-calls-share-wb-cordial-leechee-9743/release_<YYYY-MM-DD>/haplotype_units/
# gsutil -m rsync -r    ~/s05/internal/62/ gs://hla-calls-share-wb-cordial-leechee-9743/release_<YYYY-MM-DD>/haplotype_units/
# gsutil -m rsync -n -r -c ~/s05/internal/62/ gs://hla-calls-share-wb-cordial-leechee-9743/release_<YYYY-MM-DD>/haplotype_units/

# 4. Pull only ~/s05/results/62/ (PUBLIC).
