#!/bin/bash
# Fetch TRUST4's two small reference files into ../reference/. Neither is participant data --
# both are public, generic (not person-specific) receptor-locus / IMGT allele references,
# same category as the Lai truth table or IMGT XML vendored elsewhere in this repo.
#
# [MED] hg38_bcrtcr.fa -- TRUST4 ships this file directly in its own repo for hg38. Path
# below is the standard upstream location; confirm it still exists at this path before
# trusting the copy step.
# [LOW] human_IMGT+C.fa -- built by TRUST4's own BuildImgtAnnot.pl against IMGT's live FTP.
# Untested here: needs perl + internet access to IMGT from the VM. If it fails, TRUST4's
# GitHub releases/wiki may have a pre-built copy -- check there before debugging the build.
#
# Usage:  bash setup_trust4_refs.sh
set -uo pipefail

REPO=~/repos/pilot-validation
OUT="$REPO/aleix/RNA-seq/reference"
TRUST4_SRC=~/tools/TRUST4

mkdir -p "$OUT"

if [[ ! -d "$TRUST4_SRC" ]]; then
  echo "[1/3] cloning TRUST4 ..."
  mkdir -p ~/tools && cd ~/tools
  git clone https://github.com/liulab-dfci/TRUST4.git
  cd TRUST4 && make
else
  echo "[1/3] TRUST4 already present at $TRUST4_SRC, skipping clone."
fi

echo "[2/3] hg38 receptor-locus coordinates ..."
if [[ -s "$TRUST4_SRC/hg38_bcrtcr.fa" ]]; then
  cp "$TRUST4_SRC/hg38_bcrtcr.fa" "$OUT/hg38_bcrtcr.fa"
else
  echo "  !! hg38_bcrtcr.fa not found in $TRUST4_SRC -- check TRUST4's repo layout, it may"
  echo "     live under a subfolder now. See github.com/liulab-dfci/TRUST4"
fi

echo "[3/3] IMGT allele reference (needs perl + internet) ..."
if [[ ! -s "$OUT/human_IMGT+C.fa" ]]; then
  ( cd "$TRUST4_SRC" && perl BuildImgtAnnot.pl Homo_sapien > "$OUT/human_IMGT+C.fa" ) \
    || echo "  !! build failed -- check TRUST4's GitHub releases for a pre-built human_IMGT+C.fa instead"
else
  echo "  already present, skipping."
fi

echo ""
echo "Done. Expect two files in $OUT: hg38_bcrtcr.fa and human_IMGT+C.fa"
ls -la "$OUT"
