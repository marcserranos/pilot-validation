#!/bin/bash
# Run TRUST4 on one AoU RNA-seq BAM and report the recovered CDR3 count -- the single
# number that answers whether the repertoire direction is viable on this data (see
# ../README.md). Nothing below has been run against a real AoU BAM yet -- treat this as a
# scaffold to verify step by step, same posture as run_phase1_sample.sh's first runs.
#
# TCR/BCR loci sit on ordinary chromosomes (TRB chr7, IGH chr14, etc.), NOT on the
# ALT/HLA contigs AoU's RNA-seq reference excludes -- so the noALT_noHLA_noDecoy problem
# that blocks HLA work does not apply here. [MED, reasoned not yet confirmed live]
#
# Usage:  bash run_trust4_sample.sh <SAMPLE_ID> <path_to_rna_bam>
# The BAM path is whatever find_rnaseq_manifest.sh resolves for a given research_id --
# either a gcsfuse-mounted path (~/mnt/aou-controlled/...) or a locally copied file.
set -uo pipefail

SAMPLE="${1:?need a sample id}"
BAM="${2:?need a path to the RNA BAM}"

REPO=~/repos/pilot-validation
REF_DIR="$REPO/aleix/RNA-seq/reference"
MANIFEST="$REPO/aleix/RNA-seq/pixi.toml"
OUT=~/pipeline_outputs/rnaseq/$SAMPLE
THREADS=$(nproc)

BCRTCR_FA="$REF_DIR/hg38_bcrtcr.fa"
IMGT_FA="$REF_DIR/human_IMGT+C.fa"

mkdir -p "$OUT"
echo "==== TRUST4 :: $SAMPLE :: $(date) ===="

for f in "$BCRTCR_FA" "$IMGT_FA"; do
  [[ -s "$f" ]] || { echo "MISSING $f -- run setup_trust4_refs.sh first"; exit 1; }
done
[[ -s "$BAM" ]] || { echo "MISSING BAM at $BAM"; exit 1; }

pixi run --manifest-path "$MANIFEST" -- \
  run-trust4 -b "$BAM" -f "$BCRTCR_FA" --ref "$IMGT_FA" \
  -t "$THREADS" -o "$SAMPLE" --od "$OUT" \
  > "$OUT/${SAMPLE}.trust4.log" 2>&1

REPORT="$OUT/${SAMPLE}_report.tsv"
if [[ -s "$REPORT" ]]; then
  N_CDR3=$(tail -n +2 "$REPORT" | wc -l)
  echo ""
  echo "==== $SAMPLE :: RESULT ===="
  echo "CDR3s recovered: $N_CDR3"
  echo "(full report: $REPORT -- columns include chain, CDR3aa, count/abundance per TRUST4's docs)"
else
  echo "!! no report at $REPORT -- see $OUT/${SAMPLE}.trust4.log"
fi
