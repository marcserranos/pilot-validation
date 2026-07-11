#!/bin/bash
# Experiment D orchestrator -- fused ancestry-stratified 3-way comparison (AoU-native-SR /
# SpecHLA-SR / SpecImmune-LR) across the whole cohort from build_experiment_d_cohort.py.
# See EXPERIMENTS.md / DECISIONS.md "Experiment D" for the science; this is the runner.
#
# ---- Designed to be launched once and forgotten (incl. across multiple nights) ----
# The one thing that can interrupt an unattended run is the VM's ~1h idle auto-stop
# (ENVIRONMENT.md quirk #14), which kills running processes AND drops the gcsfuse mount. It
# is NOT confirmed whether active CPU counts as "activity" and prevents that sleep -- so this
# script does not rely on it. Instead it is fully RESUMABLE and IDEMPOTENT:
#   - Re-running the exact same command skips every person already finished (an `expd.done`
#     marker), and resumes any half-finished person at the first incomplete step.
#   - "Complete" is verified by per-gene COMPLETENESS, never exit code or file existence alone
#     (quirks #17 SpecImmune-exits-0-on-failure, #18 SpecHLA-silently-blank-under-starvation).
#   - Every finished person is written to disk immediately, never batched -- a mid-run sleep
#     loses at most the single in-flight person, never anything already done.
# So whether it runs straight through or sleeps-and-dies-and-you-rerun-it-once in the morning,
# the outcome is identical and no work is lost. Set NTFY_TOPIC below to get a phone ping per
# person + on completion, so you can confirm person 1 succeeded (~25 min in) then truly forget it.
#
# Config: SpecHLA at the Experiment-B-validated optimized padding, SpecImmune at pad100k/bwa
# (Experiment C's validated LR default; gene-panel restriction is NOT used -- rejected in C).
# Each person under this config also serves as the "n>=1 more person" confirmation that B's and
# C's n=1 padding optimizations still need -- watch the concordance for any config-induced drift.
#
# Usage (run from inside `pixi shell -e spechla`, or anywhere -- all tool calls are env-pinned):
#   bash run_experiment_d.sh ~/pipeline_outputs/experiment_d/cohort.tsv
# Prereqs: cohort.tsv built (build_experiment_d_cohort.py), gcsfuse mount up (quirk #11) or it
#   will be auto-mounted, reference at ~/ref/, both pixi envs built, SpecImmune ./db present.

set -uo pipefail

COHORT="${1:?Usage: bash run_experiment_d.sh <cohort.tsv>}"
[ -f "$COHORT" ] || { echo "FATAL: cohort file not found: $COHORT" >&2; exit 1; }

# ---- tunables ----
SPECHLA_PAD=10000                 # Experiment B: zero degradation across pad2000-pad10000; 10k = safety margin
SPECIMMUNE_ALIGNER=minimap2       # 2026-07-11: chosen over bwa for LR — no evidence bwa is more
                                  # accurate (DPB1 "regression" was n=1, didn't reproduce), SpecImmune's
                                  # own README examples + typing-stage default use minimap2, the AoU LR
                                  # BAMs are pbmm2(minimap2)-aligned, and it's ~30% faster. See DECISIONS.md.
MIN_SPECHLA=8                     # of 16 gene-haplotype slots -- below = catastrophic failure, redo/flag
MIN_SPECIMMUNE=4                  # of 8 classical genes with >=1 haplotype
MAX_ATTEMPTS=4                    # per person, across reruns -- backstop against infinite retry of a broken person
NTFY_TOPIC=""                     # set to an unguessable topic (e.g. "marc-hla-9f3k2x") for phone pings

# SpecImmune LR windows = pad100k (Experiment C validated), from run_aligner_pad_sweep.sh:
SPECIMMUNE_WINDOWS="chr6:29841260-30049606 chr6:31168746-31467067 chr6:32477902-32768383 chr6:32964569-33191655"

# ---- paths ----
MOUNT="$HOME/mnt/aou-controlled"
REF="$HOME/ref/Homo_sapiens_assembly38.fasta"
OUTROOT="$HOME/pipeline_outputs"
EXPDIR="$OUTROOT/experiment_d"
PROGRESS="$EXPDIR/progress.log"
PIXI_MANIFEST="$HOME/repos/pilot-validation/pixi.toml"
HELPERS="$HOME/repos/pilot-validation/spechla_pad_helpers.py"
COMPARE="$HOME/repos/pilot-validation/compare_hla_results.py"
SPECHLA_DIR="$HOME/tools/SpecHLA"
SPECIMMUNE_DIR="$HOME/tools/SpecImmune"
BILLING="wb-glacial-potato-8710"
BUCKET="vwb-aou-datasets-controlled"
mkdir -p "$EXPDIR"

log()    { echo "$(date -u +%FT%TZ) $*" | tee -a "$PROGRESS"; }
notify() { [ -n "$NTFY_TOPIC" ] && curl -s -d "$1" "ntfy.sh/$NTFY_TOPIC" >/dev/null 2>&1 || true; }
# All samtools goes through the spechla env (v1.21) so the script works regardless of which
# shell launched it -- pixi shell activation is finicky (quirks #1/#2) and this is unattended.
sam()    { pixi run --manifest-path "$PIXI_MANIFEST" -e spechla -- samtools "$@"; }

ensure_mount() {
  # The mount does not survive a VM sleep (quirk #14) and needs a beat to become ready after
  # remount (FUSE race). Only remount if a known path is not resolvable.
  if ls "$MOUNT/v9/wgs" >/dev/null 2>&1; then return 0; fi
  log "mount not resolvable -- (re)mounting $BUCKET"
  gcsfuse --billing-project "$BILLING" --implicit-dirs "$BUCKET" "$MOUNT" >>"$PROGRESS" 2>&1 || true
  for _ in 1 2 3 4 5 6; do
    sleep 3
    if ls "$MOUNT/v9/wgs" >/dev/null 2>&1; then log "mount ready"; return 0; fi
  done
  log "FATAL: mount still not resolvable after remount attempt"
  return 1
}

# Pure-stdlib completeness counters (no pandas/env dependency): print the count to stdout.
count_spechla() {
python3 - "$1" <<'PY'
import sys
rows=[l.rstrip("\n").split("\t") for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")]
if len(rows)<2: print(0); sys.exit()
d=dict(zip(rows[0],rows[1])); genes=["A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1"]
print(sum(1 for g in genes for h in ("1","2")
          if d.get(f"HLA_{g}_{h}","-").strip() not in ("-","NA","nan","")))
PY
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

DONE=0; FAIL=0; SKIP=0

process_person() {
  local pid="$1" ancestry="$2" cram_rel="$3" lr_rel="$4"
  local pdir="$OUTROOT/$pid"
  mkdir -p "$pdir"

  if [ -f "$pdir/expd.done" ]; then
    log "[$pid $ancestry] SKIP -- already complete"; SKIP=$((SKIP+1)); return
  fi
  local attempts; attempts=$(cat "$pdir/expd.attempts" 2>/dev/null || echo 0)
  if [ "$attempts" -ge "$MAX_ATTEMPTS" ]; then
    log "[$pid $ancestry] SKIP -- gave up after $attempts attempts (rm $pdir/expd.attempts to retry)"
    FAIL=$((FAIL+1)); return
  fi
  echo $((attempts+1)) > "$pdir/expd.attempts"

  ensure_mount || { log "[$pid] abort person -- mount down"; FAIL=$((FAIL+1)); return; }
  log "[$pid $ancestry] === start (attempt $((attempts+1))) ==="

  # ---- 1. SR slice -> FASTQ (SpecHLA input), pad${SPECHLA_PAD} windows ----
  if [ ! -s "$pdir/expd_R1.fastq.gz" ]; then
    sam view -b -T "$REF" "$MOUNT/$cram_rel" $SPECHLA_WINDOWS -o "$pdir/expd_sr_sliced.bam" 2>>"$PROGRESS"
    if [ ! -s "$pdir/expd_sr_sliced.bam" ]; then
      log "[$pid] FAILED SR slice (cram: $cram_rel) -- see progress.log"; FAIL=$((FAIL+1)); return; fi
    sam sort -n -o "$pdir/expd_sr_namesorted.bam" "$pdir/expd_sr_sliced.bam" 2>>"$PROGRESS"
    sam fastq -1 "$pdir/expd_R1.fastq.gz" -2 "$pdir/expd_R2.fastq.gz" -0 /dev/null \
      -s "$pdir/expd_singletons.fastq.gz" -F 0x900 "$pdir/expd_sr_namesorted.bam" 2>>"$PROGRESS"
    log "[$pid] SR fastq ready ($(sam view -c "$pdir/expd_sr_sliced.bam" 2>>"$PROGRESS") reads sliced)"
  fi

  # ---- 2. LR slice -> FASTQ (SpecImmune input), pad100k windows ----
  if [ ! -s "$pdir/expd_LR.fastq" ]; then
    sam view -b "$MOUNT/$lr_rel" $SPECIMMUNE_WINDOWS -o "$pdir/expd_lr_sliced.bam" 2>>"$PROGRESS"
    if [ ! -s "$pdir/expd_lr_sliced.bam" ]; then
      log "[$pid] FAILED LR slice (bam: $lr_rel) -- see progress.log"; FAIL=$((FAIL+1)); return; fi
    sam fastq -F 0x900 "$pdir/expd_lr_sliced.bam" > "$pdir/expd_LR.fastq" 2>>"$PROGRESS"
    log "[$pid] LR fastq ready ($(sam view -c "$pdir/expd_lr_sliced.bam" 2>>"$PROGRESS") reads sliced)"
  fi

  # ---- 3. SpecHLA (short-read) ----
  local sh_out="$pdir/expd_spechla_output"
  local sh_res="$sh_out/$pid/hla.result.txt"
  local sh_n=0
  [ -f "$sh_res" ] && sh_n=$(count_spechla "$sh_res")
  if [ ! -f "$sh_res" ] || [ "$sh_n" -lt "$MIN_SPECHLA" ]; then
    rm -rf "$sh_out"; mkdir -p "$sh_out"   # SpecHLA README: clear before rerun
    log "[$pid] SpecHLA starting"
    ( cd "$SPECHLA_DIR" && \
      { time pixi run --manifest-path "$PIXI_MANIFEST" -e spechla -- \
          bash script/whole/SpecHLA.sh -n "$pid" \
          -1 "$pdir/expd_R1.fastq.gz" -2 "$pdir/expd_R2.fastq.gz" -o "$sh_out" ; } \
        2> "$pdir/expd_spechla_timing.txt" )
    if [ ! -f "$sh_res" ]; then
      log "[$pid] FAILED SpecHLA -- no result file -- see expd_spechla_timing.txt"; FAIL=$((FAIL+1)); return; fi
    sh_n=$(count_spechla "$sh_res")
    if [ "$sh_n" -lt "$MIN_SPECHLA" ]; then
      log "[$pid] FAILED SpecHLA -- only $sh_n/16 gene-haplotypes called (< $MIN_SPECHLA); starved/broke (quirk #18)"
      FAIL=$((FAIL+1)); return; fi
    log "[$pid] SpecHLA done ($(grep real "$pdir/expd_spechla_timing.txt" | tail -1)) -- $sh_n/16 called"
  fi

  # ---- 4. SpecImmune (long-read) -- ALWAYS via `pixi run -e specimmune` (quirk #17) ----
  local si_out="$pdir/expd_specimmune_output"
  local si_res="$si_out/$pid/$pid.HLA.final.type.result.formatted.txt"
  local si_n=0
  [ -f "$si_res" ] && si_n=$(count_specimmune "$si_res")
  if [ ! -f "$si_res" ] || [ "$si_n" -lt "$MIN_SPECIMMUNE" ]; then
    rm -rf "$si_out"; mkdir -p "$si_out"
    log "[$pid] SpecImmune starting ($SPECIMMUNE_ALIGNER, pad100k)"
    ( cd "$SPECIMMUNE_DIR" && \
      { time pixi run --manifest-path "$PIXI_MANIFEST" -e specimmune -- \
          python3 scripts/main.py -n "$pid" -o "$si_out" -j 4 -y pacbio-hifi \
          -i HLA -r "$pdir/expd_LR.fastq" --db ./db \
          --align_method_1 "$SPECIMMUNE_ALIGNER" --visualization "" ; } \
        2> "$pdir/expd_specimmune_timing.txt" )
    if [ ! -f "$si_res" ]; then
      log "[$pid] FAILED SpecImmune -- exit maybe 0 but no result file (quirk #17) -- see expd_specimmune_timing.txt"
      FAIL=$((FAIL+1)); return; fi
    si_n=$(count_specimmune "$si_res")
    if [ "$si_n" -lt "$MIN_SPECIMMUNE" ]; then
      log "[$pid] FAILED SpecImmune -- only $si_n/8 genes called (< $MIN_SPECIMMUNE)"; FAIL=$((FAIL+1)); return; fi
    log "[$pid] SpecImmune done ($(grep real "$pdir/expd_specimmune_timing.txt" | tail -1)) -- $si_n/8 called"
  fi

  # ---- 5. 3-way comparison (AoU-native pulled from mounted TSV by compare) ----
  ensure_mount || { log "[$pid] compare skipped -- mount down"; FAIL=$((FAIL+1)); return; }
  if pixi run --manifest-path "$PIXI_MANIFEST" -e spechla -- python3 "$COMPARE" "$pid" \
       --run-label experiment_d \
       --spechla-result "$sh_res" \
       --specimmune-result "$si_res" \
       >> "$pdir/comparison_experiment_d.md" 2>>"$PROGRESS"; then
    touch "$pdir/expd.done"; rm -f "$pdir/expd.attempts"
    # Disk hygiene (2026-07-11, after a disk-full at ~person 10): the calls are now captured in
    # comparison_log.csv + the KB-sized result files, so drop the regenerable slices/FASTQs and the
    # bulky tool intermediates (>1MB). Without this a person leaves ~3GB behind and 60 need ~180GB;
    # with it, disk stays flat. Re-slicing a pruned person costs ~5s if ever needed.
    rm -f "$pdir"/expd_sr_sliced.bam "$pdir"/expd_sr_namesorted.bam "$pdir"/expd_lr_sliced.bam \
          "$pdir"/expd_R1.fastq.gz "$pdir"/expd_R2.fastq.gz "$pdir"/expd_singletons.fastq.gz "$pdir"/expd_LR.fastq
    find "$pdir/expd_spechla_output" "$pdir/expd_specimmune_output" -type f -size +1M -delete 2>/dev/null
    log "[$pid $ancestry] === DONE (SpecHLA $sh_n/16, SpecImmune $si_n/8) ==="
    DONE=$((DONE+1)); notify "Experiment D: $pid ($ancestry) done [$DONE done / $FAIL failed / $SKIP skip]"
  else
    log "[$pid] FAILED at comparison step -- see progress.log"; FAIL=$((FAIL+1))
  fi
}

log "================ Experiment D run starting (cohort: $COHORT) ================"
ensure_mount || { log "FATAL: cannot mount bucket, aborting"; exit 1; }

# Derive SpecHLA windows once (single source of truth = spechla_pad_helpers.py, same as Exp B).
SPECHLA_WINDOWS=$(pixi run --manifest-path "$PIXI_MANIFEST" -e spechla -- \
  python3 "$HELPERS" windows --pad "$SPECHLA_PAD")
[ -n "$SPECHLA_WINDOWS" ] || { log "FATAL: could not derive SpecHLA windows"; exit 1; }
log "SpecHLA windows (pad${SPECHLA_PAD}): $SPECHLA_WINDOWS"
log "SpecImmune windows (pad100k): $SPECIMMUNE_WINDOWS"

while IFS=$'\t' read -r pid ancestry cram_rel lr_rel bam_source; do
  [ -z "${pid:-}" ] && continue
  process_person "$pid" "$ancestry" "$cram_rel" "$lr_rel"
done < <(tail -n +2 "$COHORT")

log "================ Experiment D run complete: $DONE done, $FAIL failed, $SKIP skipped ================"
log "Aggregate analysis: python3 ~/repos/pilot-validation/analyze_experiment_d.py $COHORT"
notify "Experiment D FINISHED: $DONE done, $FAIL failed, $SKIP skipped"
