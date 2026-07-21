#!/usr/bin/env bash
# One-time Immuannot install on the Workbench VM. Idempotent -- safe to re-run.
#
# Immuannot (github.com/YingZhou001/Immuannot) types HLA/KIR alleles by aligning IPD-IMGT/HLA
# reference gene sequences against a PHASED ASSEMBLY contig (one haplotype's collapsed consensus
# FASTA), not raw reads -- this is the tool for the AoU assembly_hap1_fa/assembly_hap2_fa columns
# (DECISIONS.md "Assembly-based HLA typing"), not a BAM/FASTQ consumer like SpecImmune/HLA-Resolve.
#
# Deps: minimap2 + python3 + bash -- all already satisfied by the existing `specimmune` pixi env
# (SpecImmune itself shells out to minimap2), so this does NOT add a new pixi env, matching the
# project's minimal-footprint convention.
#
# Usage: bash scripts/setup_immuannot.sh
set -euo pipefail

TOOLS_DIR="$HOME/tools/Immuannot"
REF_URL="https://zenodo.org/records/10948964/files/Data-2024Feb02.tar.gz"
REF_TARBALL="$HOME/tools/Data-2024Feb02.tar.gz"
REF_DIR="$HOME/tools/Immuannot_refdata"

echo "== Immuannot setup ==" >&2

if [ -d "$TOOLS_DIR/.git" ]; then
    echo "Already cloned at $TOOLS_DIR -- skipping clone." >&2
else
    git clone https://github.com/YingZhou001/Immuannot.git "$TOOLS_DIR"
fi

if [ -d "$REF_DIR" ] && [ -n "$(ls -A "$REF_DIR" 2>/dev/null)" ]; then
    echo "Reference data already present at $REF_DIR -- skipping download." >&2
else
    # Actual observed size (2026-07-21): ~9MB, not "hundreds of MB" -- corrected from an earlier
    # unverified guess. Makes sense: gene-level allele sequences (IPD-IMGT/HLA + IPD-KIR + one
    # RefSeq gene), not whole genomes.
    echo "Downloading reference data (Zenodo 10.5281/zenodo.10948964, ~9MB)..." >&2
    curl -L -o "$REF_TARBALL" "$REF_URL"
    mkdir -p "$REF_DIR"
    tar -xzf "$REF_TARBALL" -C "$REF_DIR" --strip-components=1
    rm -f "$REF_TARBALL"
fi

echo "" >&2
echo "== Sanity checks ==" >&2
echo "Reference dir contents (spot-check this looks like real reference data, not an empty/failed extract):" >&2
ls -la "$REF_DIR" | head -20 >&2

# Never hardcode the in-repo script path -- confirmed 2026-07-21 the upstream repo ships it under
# a VERSIONED folder name (scripts.pub.v3, not plain scripts/), contradicting its own README and
# our earlier assumption. Same lesson as ENVIRONMENT.md quirk #13: verify by finding the real file,
# don't trust a documented path pattern that could silently drift on a future re-clone/version bump.
echo "" >&2
IMMUANNOT_SCRIPT="$(find "$TOOLS_DIR" -maxdepth 3 -name immuannot.sh 2>/dev/null | head -1)"
if [ -n "$IMMUANNOT_SCRIPT" ]; then
    echo "OK: found immuannot.sh at $IMMUANNOT_SCRIPT" >&2
else
    echo "FATAL: immuannot.sh not found anywhere under $TOOLS_DIR (checked 3 levels deep) -- " >&2
    echo "       the clone may be incomplete, or the upstream layout changed again. Run:" >&2
    echo "       find $TOOLS_DIR -iname '*immuannot*'" >&2
    exit 1
fi

echo "" >&2
echo "Checking minimap2/python3 are on PATH inside the specimmune pixi env (reused, no new env):" >&2
echo "  Run this separately: cd ~/repos/pilot-validation && pixi run -e specimmune -- which minimap2 python3" >&2
echo "" >&2
echo "Done. Next: scripts/run_immuannot_person.py <person_id> [<person_id> ...]" >&2
