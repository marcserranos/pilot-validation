#!/bin/bash
# Builds a SpecImmune HLA database restricted to the 8 classical genes (A, B, C, DRB1,
# DQA1, DQB1, DPA1, DPB1), instead of the ~30-gene default panel. See EXPERIMENTS.md
# Experiment C (2026-07-10) for the full investigation: confirmed via source-reading that
# every per-gene processing stage in main.py's pipeline discovers genes dynamically from
# the --db folder's subdirectories (get_folder_list()), not from a hardcoded list -- so
# restricting the DB genuinely restricts processing.
#
# Prereqs: patch_specimmune_for_gene_restriction.py must already have been run once
# against ~/tools/SpecImmune (fixes a real UnboundLocalError bug in make_db.py's
# --HLA_fa/--HLA_exon_fa branches -- without it this script will crash).
#
# Reuses the IMGT source FASTAs already downloaded when the full-panel db/ was first
# built (~/tools/SpecImmune/db/HLA/hla_gen.fasta, hla_nuc.fasta) -- no re-download needed,
# and guarantees the exact same IMGT release as the existing full-panel db for a clean
# apples-to-apples comparison.

set -euo pipefail

SPECIMMUNE_DIR="$HOME/tools/SpecImmune"
FULL_DB="$SPECIMMUNE_DIR/db"
PIXI_MANIFEST="$HOME/repos/pilot-validation/pixi.toml"

GEN_SRC="$FULL_DB/HLA/hla_gen.fasta"
NUC_SRC="$FULL_DB/HLA_CDS/hla_nuc.fasta"

if [ ! -f "$GEN_SRC" ] || [ ! -f "$NUC_SRC" ]; then
  echo "FATAL: expected source FASTAs not found at:"
  echo "  $GEN_SRC"
  echo "  $NUC_SRC"
  echo "These should already exist from the original 'make scripts/make_db.py -o ./db -i HLA' build."
  echo "Check the actual layout with: find $FULL_DB -maxdepth 2 -iname '*.fasta'"
  exit 1
fi

RESTRICTED_DB="$SPECIMMUNE_DIR/db_classical8"
FILTERED_DIR="$SPECIMMUNE_DIR/db_classical8_source"
mkdir -p "$FILTERED_DIR"

echo "=== Filtering to the 8 classical genes ==="
pixi run --manifest-path "$PIXI_MANIFEST" -e specimmune -- python3 -c "
from Bio import SeqIO

CLASSICAL = {'A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1'}

def filt(inp, outp):
    kept, total = 0, 0
    with open(outp, 'w') as out:
        for rec in SeqIO.parse(inp, 'fasta'):
            total += 1
            parts = rec.description.split()
            gene = parts[1].split('*')[0]
            if gene in CLASSICAL:
                kept += 1
                SeqIO.write(rec, out, 'fasta')
    print(f'{inp}: kept {kept} / {total} records')

filt('$GEN_SRC', '$FILTERED_DIR/hla_gen.classical8.fasta')
filt('$NUC_SRC', '$FILTERED_DIR/hla_nuc.classical8.fasta')
"

echo "=== Building restricted DB at $RESTRICTED_DB ==="
rm -rf "$RESTRICTED_DB"
mkdir -p "$RESTRICTED_DB"
cd "$SPECIMMUNE_DIR/scripts"
time pixi run --manifest-path "$PIXI_MANIFEST" -e specimmune -- python3 make_db.py \
  -o "$RESTRICTED_DB" -i HLA \
  --HLA_fa "$FILTERED_DIR/hla_gen.classical8.fasta" \
  --HLA_exon_fa "$FILTERED_DIR/hla_nuc.classical8.fasta"

echo ""
echo "=== Verifying only 8 gene directories were built (not ~30) ==="
echo "HLA/:"
find "$RESTRICTED_DB/HLA" -maxdepth 1 -type d | sort
echo "HLA_CDS/:"
find "$RESTRICTED_DB/HLA_CDS" -maxdepth 1 -type d | sort

echo ""
echo "Restricted DB ready at: $RESTRICTED_DB"
echo "Run SpecImmune against it with: --db $RESTRICTED_DB (same as any other run, just swap --db)"
