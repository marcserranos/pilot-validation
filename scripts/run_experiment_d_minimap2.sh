#!/bin/bash
# Experiment D -- minimap2 aligner sweep (n=60 bwa-vs-minimap2 confirmation).
# RUN THIS ONLY AFTER run_experiment_d.sh (the bwa pass) has fully finished for the cohort --
# never concurrently. This is a timing experiment; two SpecImmune jobs sharing the VM's 4 vCPUs
# would contaminate the exact runtime being measured (the same rule Experiment C established and
# Marc insisted on -- EXPERIMENTS.md/DECISIONS.md 2026-07-10).
#
# Why this exists: the bwa-vs-minimap2 choice for SpecImmune read binning (--align_method_1) was
# never backed by evidence -- bwa is only the flag's *default*, while SpecImmune's own README
# examples use minimap2, its typing stage (--align_method_2) already defaults to minimap2, and the
# AoU long-read BAMs are pbmm2 (minimap2)-aligned to begin with. minimap2 is also ~30% faster
# (~14.5 vs ~21 min/person at pad100k). This sweep settles it at n=60 on BOTH accuracy and runtime.
#
# Controlled design: reuses each person's EXISTING pad100k LR FASTQ from the bwa pass
# (~/pipeline_outputs/<pid>/expd_LR.fastq) -- identical input, so --align_method_1 is the single
# variable. Separate output dir (expd_specimmune_output_minimap2) + timing file + run-label so
# nothing collides with the bwa results; the two are compared afterward by analyze_experiment_d.py.
#
# Resumable/idempotent exactly like the bwa orchestrator: per-person isolation, no set -e,
# completeness-gated (not exit-code-trusting, quirk #17), expd.minimap2.done markers, attempts
# backstop, mount auto-repair. Safe to re-run the same command until every person is done.
#
# Usage: bash run_experiment_d_minimap2.sh ~/pipeline_outputs/experiment_d/cohort.tsv

set -uo pipefail

COHORT="${1:?Usage: bash run_experiment_d_minimap2.sh <cohort.tsv>}"
[ -f "$COHORT" ] || { echo "FATAL: cohort file not found: $COHORT" >&2; exit 1; }

MIN_SPECIMMUNE=4
MAX_ATTEMPTS=4
NTFY_TOPIC=""                     # set to the same unguessable topic as the bwa run if you want pings

MOUNT="$HOME/mnt/aou-controlled"
OUTROOT="$HOME/pipeline_outputs"
EXPDIR="$OUTROOT/experiment_d"
PROGRESS="$EXPDIR/progress_minimap2.log"
PIXI_MANIFEST="$HOME/repos/pilot-validation/pixi.toml"
COMPARE="$HOME/repos/pilot-validation/scripts/compare_hla_results.py"
SPECIMMUNE_DIR="$HOME/tools/SpecImmune"
BILLING="wb-glacial-potato-8710"
BUCKET="vwb-aou-datasets-controlled"
mkdir -p "$EXPDIR"

log()    { echo "$(date -u +%FT%TZ) $*" | tee -a "$PROGRESS"; }
notify() { [ -n "$NTFY_TOPIC" ] && curl -s -d "$1" "ntfy.sh/$NTFY_TOPIC" >/dev/null 2>&1 || true; }

ensure_mount() {
  if ls "$MOUNT/v9/wgs" >/dev/null 2>&1; then return 0; fi
  log "mount not resolvable -- (re)mounting $BUCKET"
  gcsfuse --billing-project "$BILLING" --implicit-dirs "$BUCKET" "$MOUNT" >>"$PROGRESS" 2>&1 || true
  for _ in 1 2 3 4 5 6; do
    sleep 3
    if ls "$MOUNT/v9/wgs" >/dev/null 2>&1; then log "mount ready"; return 0; fi
  done
  log "FATAL: mount still not resolvable"; return 1
}

count_specimmune() {
python3 - "$1" <<'PY'
import sys
genes={"A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1"}; hdr=None; seen=set()
for l in open(sys.argv[1]):
    if not l.strip() or l.startswith("#"): continue
    c=l.rstrip("\n").split("\t")
    if hdr is None: hdr=c; continue
    def gv(name):
        return c[hdr.index(name)] if name in hdr and hdr.index(name)<len(c) else ""
    locus=gv("Locus").replace("HLA-","")
    call=gv("One_guess") or ""
    if call in ("","NA","-","nan","."): call=(gv("Genotype") or "").split(";")[0]
    if locus in genes and call and call not in ("NA","-","nan","."): seen.add(locus)
print(len(seen))
PY
}

DONE=0; FAIL=0; SKIP=0; WAIT=0

process_person() {
  local pid="$1" ancestry="$2"
  local pdir="$OUTROOT/$pid"
  mkdir -p "$pdir"

  if [ -f "$pdir/expd.minimap2.done" ]; then
    log "[$pid $ancestry] SKIP -- minimap2 already complete"; SKIP=$((SKIP+1)); return
  fi
  # This pass reuses the bwa pass's sliced LR FASTQ -- if it's not there yet, the bwa run hasn't
  # produced this person. Skip (don't fail): re-run this sweep after the bwa pass finishes.
  local fastq="$pdir/expd_LR.fastq"
  if [ ! -s "$fastq" ]; then
    log "[$pid $ancestry] WAIT -- no expd_LR.fastq yet (bwa pass hasn't done this person); re-run later"
    WAIT=$((WAIT+1)); return
  fi
  local attempts; attempts=$(cat "$pdir/expd.minimap2.attempts" 2>/dev/null || echo 0)
  if [ "$attempts" -ge "$MAX_ATTEMPTS" ]; then
    log "[$pid $ancestry] SKIP -- gave up after $attempts attempts (rm $pdir/expd.minimap2.attempts to retry)"
    FAIL=$((FAIL+1)); return
  fi
  echo $((attempts+1)) > "$pdir/expd.minimap2.attempts"

  local si_out="$pdir/expd_specimmune_output_minimap2"
  local si_res="$si_out/$pid/$pid.HLA.final.type.result.formatted.txt"
  local si_n=0
  [ -f "$si_res" ] && si_n=$(count_specimmune "$si_res")
  if [ ! -f "$si_res" ] || [ "$si_n" -lt "$MIN_SPECIMMUNE" ]; then
    rm -rf "$si_out"; mkdir -p "$si_out"
    log "[$pid $ancestry] SpecImmune (minimap2, pad100k) starting"
    ( cd "$SPECIMMUNE_DIR" && \
      { time pixi run --manifest-path "$PIXI_MANIFEST" -e specimmune -- \
          python3 scripts/main.py -n "$pid" -o "$si_out" -j 4 -y pacbio-hifi \
          -i HLA -r "$fastq" --db ./db \
          --align_method_1 minimap2 --visualization "" ; } \
        2> "$pdir/expd_specimmune_minimap2_timing.txt" )
    if [ ! -f "$si_res" ]; then
      log "[$pid] FAILED SpecImmune-minimap2 -- no result (quirk #17) -- see expd_specimmune_minimap2_timing.txt"
      FAIL=$((FAIL+1)); return; fi
    si_n=$(count_specimmune "$si_res")
    if [ "$si_n" -lt "$MIN_SPECIMMUNE" ]; then
      log "[$pid] FAILED SpecImmune-minimap2 -- only $si_n/8 genes called (< $MIN_SPECIMMUNE)"; FAIL=$((FAIL+1)); return; fi
    log "[$pid] SpecImmune-minimap2 done ($(grep real "$pdir/expd_specimmune_minimap2_timing.txt" | tail -1)) -- $si_n/8 called"
  fi

  # Compare with the SAME short-read SpecHLA result as the bwa pass (identical SR input) so the
  # only thing differing between run_label experiment_d and experiment_d_minimap2 is the LR aligner.
  ensure_mount || { log "[$pid] compare skipped -- mount down"; FAIL=$((FAIL+1)); return; }
  local sh_res="$pdir/expd_spechla_output/$pid/hla.result.txt"
  local sh_arg=(); [ -f "$sh_res" ] && sh_arg=(--spechla-result "$sh_res")
  if pixi run --manifest-path "$PIXI_MANIFEST" -e spechla -- python3 "$COMPARE" "$pid" \
       --run-label experiment_d_minimap2 \
       "${sh_arg[@]}" \
       --specimmune-result "$si_res" \
       >> "$pdir/comparison_experiment_d_minimap2.md" 2>>"$PROGRESS"; then
    touch "$pdir/expd.minimap2.done"; rm -f "$pdir/expd.minimap2.attempts"
    log "[$pid $ancestry] === minimap2 DONE (SpecImmune $si_n/8) ==="
    DONE=$((DONE+1)); notify "Exp D minimap2: $pid ($ancestry) done [$DONE done / $FAIL fail / $WAIT waiting]"
  else
    log "[$pid] FAILED at comparison step"; FAIL=$((FAIL+1))
  fi
}

log "================ Experiment D minimap2 sweep starting (cohort: $COHORT) ================"
ensure_mount || { log "FATAL: cannot mount, aborting"; exit 1; }

while IFS=$'\t' read -r pid ancestry cram_rel lr_rel; do
  [ -z "${pid:-}" ] && continue
  process_person "$pid" "$ancestry"
done < <(tail -n +2 "$COHORT")

log "================ minimap2 sweep pass complete: $DONE done, $FAIL failed, $SKIP skipped, $WAIT waiting-on-bwa ================"
[ "$WAIT" -gt 0 ] && log "NOTE: $WAIT people had no LR FASTQ yet -- re-run this same command after the bwa pass finishes them."
log "Then compare aligners: python3 ~/repos/pilot-validation/scripts/analyze_experiment_d.py $COHORT  (bwa vs minimap2 + runtimes)"
notify "Exp D minimap2 sweep pass done: $DONE done, $FAIL failed, $WAIT waiting"
