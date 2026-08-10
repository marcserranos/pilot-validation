#!/bin/bash
# Fetch TRUST4's two small reference files into ../reference/. Neither is participant data --
# both are public, generic (not person-specific) receptor-locus / IMGT allele references,
# same category as the Lai truth table or IMGT XML vendored elsewhere in this repo.
#
# [CONFIRMED, 2026-08-10] Both files ship pre-built in TRUST4's own repo -- no `make`, no
# perl build step, no compiler needed. `human_IMGT+C.fa` is NOT built live against IMGT's
# FTP as earlier assumed; it's already sitting in the repo. This script no longer builds
# the TRUST4 binary from source either -- the pixi-installed `trust4` package (bioconda)
# already provides a working `run-trust4`; cloning is only to grab these two data files.
#
# Usage:  bash setup_trust4_refs.sh
set -uo pipefail

REPO=~/repos/pilot-validation
OUT="$REPO/aleix/RNA-seq/reference"
TRUST4_SRC=~/tools/TRUST4

mkdir -p "$OUT"

if [[ ! -d "$TRUST4_SRC" ]]; then
  echo "[1/2] cloning TRUST4 (source only, not building) ..."
  mkdir -p ~/tools
  git clone --depth 1 https://github.com/liulab-dfci/TRUST4.git "$TRUST4_SRC"
else
  echo "[1/2] TRUST4 source already present at $TRUST4_SRC, skipping clone."
fi

echo "[2/2] copying the two reference files ..."
for f in hg38_bcrtcr.fa "human_IMGT+C.fa"; do
  if [[ -s "$TRUST4_SRC/$f" ]]; then
    cp "$TRUST4_SRC/$f" "$OUT/$f"
  else
    echo "  !! $f not found in $TRUST4_SRC -- check TRUST4's repo layout has changed."
  fi
done

echo ""
echo "Done. Expect two files in $OUT: hg38_bcrtcr.fa and human_IMGT+C.fa"
ls -la "$OUT"
