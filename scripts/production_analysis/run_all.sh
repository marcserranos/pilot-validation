#!/usr/bin/env bash
# Single entrypoint: run every production_analysis script, in the right order, with the right
# preflight checks, printing one clear summary at the end -- Marc, 2026-08-10: "I can just run
# one script... I don't have to be micro-managing the task."
#
# Usage (from ~/repos/pilot-validation):
#   pixi run -e spechla -- bash scripts/production_analysis/run_all.sh
# (`pixi run` activates the env for this one command -- no separate `pixi shell` step needed.)
#
# What this does, in order:
#   1. Checks immuannot_calls.tsv actually has classical-gene rows (the 2026-08-10
#      merge_fragments() dedup bug's exact symptom). If not, runs the repair
#      (rebuild_immuannot_calls.py) automatically -- it's safe to auto-run: the old file is
#      renamed aside with a timestamp, never deleted, so nothing is lost even if run repeatedly.
#   2. Mounts gcsfuse (needed only by analyze_confidence_vs_aou_native.py) if not already up.
#      Non-fatal if this fails -- that one script is skipped with a clear reason, the other 4 run.
#   3. Runs all 5 scripts. A failure in one does not stop the others (except step 1's gate, which
#      is fatal -- there's no point running 5 scripts against data already known to be broken).
#   4. Prints a PASS/FAIL summary table and exactly where every output landed.
set -uo pipefail  # deliberately NOT -e -- see "does not stop the others" above

OUTROOT="${OUTROOT:-$HOME/pipeline_outputs}"
CALLS="$OUTROOT/immuannot_calls.tsv"
COHORT="$OUTROOT/immuannot_cohort_full.tsv"
MOUNT="$HOME/mnt/aou-controlled"
AOU_TSV="$MOUNT/v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv"
# Confirmed 2026-08-10 (Marc, `gcloud config get-value project`) -- override with
# AOU_BILLING_PROJECT=... if this run ever happens in a different workspace.
BILLING_PROJECT="${AOU_BILLING_PROJECT:-wb-cordial-leechee-9743}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== production_analysis: run_all.sh starting ==="
echo "OUTROOT=$OUTROOT"

# --- Step 1: gate on the known merge-bug symptom, auto-repair if needed ---
echo ""
echo "--- Step 1/4: checking immuannot_calls.tsv for classical-gene rows ---"
if [ ! -f "$CALLS" ]; then
  echo "FATAL: $CALLS does not exist. Check OUTROOT / that this is the right VM/disk." >&2
  exit 1
fi
N_CLASSICAL=$(python3 -c "
import pandas as pd
GENES = {'A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1'}
df = pd.read_csv('$CALLS', sep='\t', dtype=str, keep_default_na=False)
bare = df['gene'].str.replace('^HLA-', '', regex=True)
print(bare.isin(GENES).sum())
")
echo "  $N_CLASSICAL rows currently match a classical HLA gene."
if [ "$N_CLASSICAL" -eq 0 ]; then
  echo "  0 classical-gene rows -- this is the known 2026-08-10 merge_fragments() dedup bug."
  echo "  Measuring repair speed on a 200-person sample before committing to the full run..."
  python3 "$SCRIPT_DIR/../production_orchestrator/rebuild_immuannot_calls.py" --outroot "$OUTROOT" --limit 200 2>&1 | tee /tmp/repair_timing_test.log
  EST_MIN=$(grep -oE '~[0-9.]+ min' /tmp/repair_timing_test.log | tail -1 | grep -oE '[0-9.]+')
  MAX_MIN="${MAX_AUTO_REPAIR_MIN:-15}"
  EXCEEDS=$(python3 -c "print(1 if float('${EST_MIN:-0}') > float('$MAX_MIN') else 0)" 2>/dev/null || echo 0)
  if [ -z "$EST_MIN" ]; then
    echo "  WARNING: could not parse the timing estimate -- proceeding with the full repair anyway (couldn't gate on an unreadable estimate)."
  elif [ "$EXCEEDS" = "1" ]; then
    echo ""
    echo "STOPPING: estimated full repair time (~${EST_MIN} min) exceeds the auto-run threshold (${MAX_MIN} min)." >&2
    echo "This is NOT a failure -- it's a deliberate pause so a long operation never runs silently." >&2
    echo "Options: (a) rerun this same command with MAX_AUTO_REPAIR_MIN=<bigger number> to raise the threshold and let it proceed," >&2
    echo "         (b) run the repair manually: python3 $SCRIPT_DIR/../production_orchestrator/rebuild_immuannot_calls.py --outroot $OUTROOT" >&2
    echo "         (c) resize to a bigger VM first if you want it faster (see scripts/production_analysis/README.md)." >&2
    exit 1
  fi
  echo "  Estimate (~${EST_MIN:-unknown} min) is within the ${MAX_MIN}-min auto-run threshold -- running the real repair now..."
  python3 "$SCRIPT_DIR/../production_orchestrator/rebuild_immuannot_calls.py" --outroot "$OUTROOT"
  REPAIR_STATUS=$?
  if [ $REPAIR_STATUS -ne 0 ]; then
    echo "FATAL: repair script itself failed (exit $REPAIR_STATUS) -- see its output above. Not proceeding." >&2
    exit 1
  fi
  N_CLASSICAL=$(python3 -c "
import pandas as pd
GENES = {'A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1'}
df = pd.read_csv('$CALLS', sep='\t', dtype=str, keep_default_na=False)
bare = df['gene'].str.replace('^HLA-', '', regex=True)
print(bare.isin(GENES).sum())
")
  if [ "$N_CLASSICAL" -eq 0 ]; then
    echo "FATAL: still 0 classical-gene rows after repair -- this is a deeper problem than the merge bug (see repair script's own WARNING above). Not proceeding." >&2
    exit 1
  fi
  echo "  Repair produced $N_CLASSICAL classical-gene rows -- proceeding."
else
  echo "  Looks fine -- no repair needed."
fi

# --- Step 2: gcsfuse mount (best-effort, non-fatal) ---
echo ""
echo "--- Step 2/4: checking gcsfuse mount (needed only for the AoU-native comparison script) ---"
RUN_CONFIDENCE=1
if ls "$MOUNT/v9/wgs" >/dev/null 2>&1; then
  echo "  Already mounted and resolving."
else
  echo "  Not mounted -- mounting (billing project: $BILLING_PROJECT)..."
  mkdir -p "$MOUNT"
  gcsfuse --billing-project "$BILLING_PROJECT" --implicit-dirs vwb-aou-datasets-controlled "$MOUNT" >/dev/null 2>&1
  sleep 3
  if ls "$MOUNT/v9/wgs" >/dev/null 2>&1; then
    echo "  Mount succeeded."
  else
    echo "  WARNING: mount did not resolve -- skipping analyze_confidence_vs_aou_native.py (the other 4 scripts don't need it)."
    RUN_CONFIDENCE=0
  fi
fi

# --- Step 3: run everything ---
# Plain accumulated text, not an associative array -- `declare -A` needs bash 4+, and the target
# VM's bash version isn't worth assuming. Each line is "name<TAB>status".
echo ""
echo "--- Step 3/4: running all scripts ---"
SUMMARY=""
run_one() {
  local name="$1"; shift
  echo ""
  echo ">>> $name"
  if "$@"; then
    SUMMARY="${SUMMARY}${name}	PASS
"
  else
    SUMMARY="${SUMMARY}${name}	FAIL (exit $?)
"
  fi
}

run_one "completeness_and_demographics" python3 "$SCRIPT_DIR/analyze_completeness_and_demographics.py" \
  --cohort "$COHORT" --calls "$CALLS" --timing "$OUTROOT/immuannot_timing.tsv" \
  --out-dir "$OUTROOT/production_analysis/completeness"
if [ "$RUN_CONFIDENCE" -eq 1 ]; then
  run_one "confidence_vs_aou_native" python3 "$SCRIPT_DIR/analyze_confidence_vs_aou_native.py" \
    --calls "$CALLS" --aou-tsv "$AOU_TSV" --outroot "$OUTROOT" \
    --out-dir "$OUTROOT/production_analysis/confidence"
else
  SUMMARY="${SUMMARY}confidence_vs_aou_native	SKIPPED (mount unavailable)
"
fi
run_one "cluster_hla_by_ancestry" python3 "$SCRIPT_DIR/cluster_hla_by_ancestry.py" \
  --calls "$CALLS" --cohort "$COHORT" --out-dir "$OUTROOT/production_analysis/clustering"
run_one "allele_frequency_by_ancestry" python3 "$SCRIPT_DIR/analyze_allele_frequency_by_ancestry.py" \
  --calls "$CALLS" --cohort "$COHORT" --out-dir "$OUTROOT/production_analysis/allele_frequency"
run_one "drb1_evidence_capstone" python3 "$SCRIPT_DIR/summarize_drb1_evidence_capstone.py" \
  --out-dir "$OUTROOT/production_analysis/drb1_capstone"

# --- Step 4: summary ---
echo ""
echo "--- Step 4/4: summary ---"
printf "%s\n" "$SUMMARY" | while IFS=$'\t' read -r name status; do
  [ -z "$name" ] && continue
  printf "  %-32s %s\n" "$name" "$status"
done
echo ""
echo "All output under: $OUTROOT/production_analysis/"
echo "  completeness/  confidence/  clustering/  allele_frequency/  drb1_capstone/"
echo "Each holds its PNG figure(s) + a *_report.md. Open PNGs in the Jupyter file browser, or"
echo "paste the *_report.md tables directly into your supervisor report."
echo "=== done ==="
