#!/bin/bash
# Aligner x padding sweep for SpecImmune-LR, single person, controlled comparison.
# See DECISIONS.md / EXPERIMENTS.md 2026-07-09 entries for the full reasoning:
#   (1) --align_method_1 bwa vs minimap2 trades ~3.25x binning speed for divergence-
#       robustness (minimap2 needs exact minimizer anchors and loses reads where an HLA
#       allele diverges a lot from GRCh38; bwa's full Smith-Waterman extension tolerates
#       that divergence better -- mechanistic explanation for the earlier DPB1 regression).
#   (2) The current 4Mb contiguous slice pads a single block, but the 8 classical genes
#       sit in 4 clusters (A | C+B | DRB1+DQA1+DQB1 | DPA1+DPB1) separated by up to
#       ~1.3Mb of non-classical sequence contributing zero classical-gene reads. Carving
#       those gaps out (multi-region samtools query, one FASTQ out) is a bigger win than
#       the aligner swap and doesn't trade away bwa's accuracy.
#
# Usage: bash run_aligner_pad_sweep.sh <person_id>
# Prereqs: any pixi env active is fine for slicing (samtools is present in spechla).
#   SpecImmune itself is ALWAYS invoked via `pixi run -e specimmune`, regardless of which
#   env launched this script -- SpecImmune's long-read pipeline shells out to `sniffles`
#   (SV caller), which only exists on PATH inside the `specimmune` env. An earlier version
#   of this script assumed the active shell's env didn't matter -- it does; running from
#   `spechla` caused every SpecImmune call to fail deep inside its pipeline (sniffles: command
#   not found) while main.py still exited 0, silently producing zero real output across an
#   entire sweep (see EXPERIMENTS.md 2026-07-09 postmortem). Fixed by never trusting the
#   ambient shell env for this step.
#   ~/pipeline_outputs/<person_id>/chr6_LR.bam and LR.fastq must already exist (i.e.
#   slice_and_fastq.sh + the full-pad bwa SpecImmune baseline already ran for this person).
#
# Gene coordinates: GRCh38, fetched live from Ensembl REST 2026-07-09 (not hand-typed).
# Read-length floor: measured from this person's actual chr6_LR.bam via samtools (median
#   16.5kb, p99 34kb, max 42.6kb, n=3332 reads) -- pad40k below uses p99 rounded up.
#
# Safe to walk away from: wrap the whole script in nohup (see bottom of this file for the
# exact command) -- a dropped SSH/browser connection would otherwise SIGHUP-kill it.
# Does NOT use `set -e`: one failed config is logged and skipped, the rest of the sweep
# continues, and every completed run is written to disk immediately (not batched at the
# end) so a mid-sweep VM sleep (see ENVIRONMENT.md quirk #14) loses at most the in-flight
# run, never anything already finished.

set -uo pipefail

PERSON_ID="${1:?Usage: bash run_aligner_pad_sweep.sh <person_id>}"
OUTDIR="$HOME/pipeline_outputs/$PERSON_ID"
LR_BAM="$OUTDIR/chr6_LR.bam"
SPECIMMUNE_DIR="$HOME/tools/SpecImmune"
PIXI_MANIFEST="$HOME/repos/pilot-validation/pixi.toml"
SWEEP_DIR="$OUTDIR/sweep"
PROGRESS="$SWEEP_DIR/progress.log"
mkdir -p "$SWEEP_DIR"

log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$PROGRESS"; }

if [ ! -f "$LR_BAM" ]; then
  log "FATAL: $LR_BAM not found -- run slice_and_fastq.sh for $PERSON_ID first."
  exit 1
fi
samtools index -f "$LR_BAM"

# ---- window definitions (GRCh38, samtools region syntax) ----
# Raw gene clusters (no padding), from Ensembl (2026-07-09):
#   A:              29,941,260-29,949,606      (~1.32Mb gap to C)
#   C,B:            31,268,746-31,367,067      (~1.21Mb gap to DRB1)
#   DRB1,DQA1,DQB1: 32,577,902-32,668,383      (~396kb gap to DPA1)
#   DPA1,DPB1:      33,064,569-33,091,655
declare -A WINDOWS
WINDOWS[fullpad]="chr6:29500000-33500000"
WINDOWS[pad200k]="chr6:29741260-30149606 chr6:31068746-31567067 chr6:32377902-33291655"
WINDOWS[pad100k]="chr6:29841260-30049606 chr6:31168746-31467067 chr6:32477902-32768383 chr6:32964569-33191655"
WINDOWS[pad40k]="chr6:29901260-29989606 chr6:31228746-31407067 chr6:32537902-32708383 chr6:33024569-33131655"
WINDOWS[pad5k]="chr6:29936260-29954606 chr6:31263746-31372067 chr6:32572902-32673383 chr6:33059569-33096655"
# Total bp per config (documented, not computed at runtime): fullpad=4.00Mb, pad200k=1.82Mb
# (clusters 3+4 merge at this pad), pad100k=1.02Mb, pad40k=0.54Mb ("theoretical minimum",
# ~7.3x smaller than fullpad), pad5k=0.26Mb (deliberately below the read-length floor --
# expect real degradation, this is the control-in-the-other-direction).

ALIGNERS="bwa minimap2"

run_one() {
  local pad_label="$1" aligner="$2"
  local label="${pad_label}_${aligner}"
  local rundir="$SWEEP_DIR/$label"
  mkdir -p "$rundir"

  local fastq
  if [ "$pad_label" == "fullpad" ]; then
    fastq="$OUTDIR/LR.fastq"   # reuse existing full-pad FASTQ, no reslicing needed
  else
    local bam="$rundir/LR.bam"
    samtools view -b "$LR_BAM" ${WINDOWS[$pad_label]} -o "$bam" 2>>"$PROGRESS"
    local nreads
    nreads=$(samtools view -c "$bam" 2>>"$PROGRESS")
    log "$label: sliced $nreads reads (window total documented above)"
    fastq="$rundir/LR.fastq"
    samtools fastq -F 0x900 "$bam" > "$fastq" 2>>"$PROGRESS"
  fi

  local specout="$rundir/specimmune_output"
  local timing="$rundir/timing.txt"
  local result_file="$specout/$PERSON_ID/$PERSON_ID.HLA.final.type.result.formatted.txt"
  log "$label: SpecImmune starting (via pixi run -e specimmune)"
  ( cd "$SPECIMMUNE_DIR" && \
    { time pixi run --manifest-path "$PIXI_MANIFEST" -e specimmune -- \
        python3 scripts/main.py -n "$PERSON_ID" -o "$specout" -j 4 -y pacbio-hifi \
        -i HLA -r "$fastq" --db ./db --align_method_1 "$aligner" --visualization "" ; } \
      2> "$timing" )
  local rc=$?
  if [ $rc -ne 0 ]; then
    log "$label: FAILED (exit $rc) -- see $timing -- continuing to next config"
    return
  fi
  # main.py can exit 0 even when an internal step (e.g. sniffles) failed deep in its
  # pipeline -- never trust exit code alone (see the sniffles/wrong-env postmortem above).
  # Verify the actual expected output file exists before declaring success.
  if [ ! -f "$result_file" ]; then
    log "$label: FAILED -- exit 0 but expected output missing at $result_file -- see $timing -- continuing to next config"
    return
  fi
  log "$label: SpecImmune done ($(grep real "$timing" | tail -1))"

  python3 "$HOME/repos/pilot-validation/compare_hla_results.py" "$PERSON_ID" \
    --run-label "$label" \
    --specimmune-result "$result_file" \
    >> "$SWEEP_DIR/comparisons.md" 2>>"$PROGRESS"
  log "$label: DONE"
}

log "=== sweep starting for person $PERSON_ID ==="

# Reuse the already-completed full-pad bwa baseline -- no rerun, just fold it into the matrix.
log "fullpad_bwa: reusing existing baseline (no rerun)"
python3 "$HOME/repos/pilot-validation/compare_hla_results.py" "$PERSON_ID" \
  --run-label "fullpad_bwa" >> "$SWEEP_DIR/comparisons.md" 2>>"$PROGRESS"

for pad_label in fullpad pad200k pad100k pad40k pad5k; do
  for aligner in $ALIGNERS; do
    if [ "$pad_label" == "fullpad" ] && [ "$aligner" == "bwa" ]; then
      continue  # already logged above
    fi
    run_one "$pad_label" "$aligner"
  done
done

log "=== sweep complete: see $SWEEP_DIR/comparisons.md and $OUTDIR/comparison_log.csv ==="

# Optional push notification (ntfy.sh -- free, no signup, public topic-name-as-secret).
# Uncomment and pick your own unguessable topic name if you want a phone ping on completion:
# curl -s -d "HLA sweep for $PERSON_ID finished" ntfy.sh/YOUR-UNGUESSABLE-TOPIC-NAME >/dev/null
