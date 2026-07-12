#!/bin/bash
# Experiment C recon -- read-only, fast (seconds). Needed before designing the
# restricted-DB build step correctly, rather than guessing at make_db.py's real interface.
# See EXPERIMENTS.md roadmap, Experiment C.
set -uo pipefail

OUT=~/pipeline_outputs/experiment_c_recon.txt
{
  echo "=== make_db.py --help ==="
  pixi run --manifest-path ~/repos/pilot-validation/pixi.toml -e specimmune -- \
    python3 ~/tools/SpecImmune/scripts/make_db.py --help

  echo ""
  echo "=== make_db.py full source ==="
  cat ~/tools/SpecImmune/scripts/make_db.py

  echo ""
  echo "=== existing full-panel db/ layout (top 2 levels) ==="
  find ~/tools/SpecImmune/db -maxdepth 2 | sort

  echo ""
  echo "=== individual_ref/ from an already-completed run (person 2522883) -- which genes actually got built ==="
  find ~/pipeline_outputs/2522883/specimmune_output/2522883/individual_ref -maxdepth 1 2>/dev/null | sort

} > "$OUT" 2>&1

echo "Recon written to $OUT"
echo "cat it and paste the contents back to Claude:"
echo "  cat $OUT"
