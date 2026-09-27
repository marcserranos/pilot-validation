#!/bin/bash
# Disk policy, as code. Prunes bulky, regenerable intermediates from per-person pipeline output
# dirs while keeping the distilled must-have outputs. Use it to reclaim space after (or during) a
# run, and to clean up old closed-experiment cruft as the project scales.
#
# TIERS (see the disk-full postmortem, 2026-07-11):
#   MUST-KEEP (distilled science, all <1MB): comparison_log.csv, comparison_*.md, tool result files
#     (hla.result.txt, *.formatted.txt), *_timing.txt, and the expd.done / expd.attempts markers.
#   DELETE (regenerable or bulky intermediate): sliced BAMs, FASTQs (re-slice from the mount in ~5s),
#     and every tool intermediate >1MB (alignments, assemblies, per-sample indices).
#
# The policy reduces to one safe rule -- inside a person dir, delete any file >1MB -- because every
# must-keep file is well under 1MB and everything bulky is regenerable. Explicit small-file types are
# additionally protected below as a belt-and-suspenders guard, so this can never eat a result/log.
#
# Only touches numeric person dirs (~/pipeline_outputs/<digits>/). Leaves analysis dirs
# (experiment_d/, experiment_a_ext/, etc.) alone. Never touches the git repo or the mount.
#
# Usage:
#   bash prune_pipeline_outputs.sh --dry-run     # show what would be freed, delete nothing
#   bash prune_pipeline_outputs.sh               # actually prune
#   bash prune_pipeline_outputs.sh --dry-run 2522883 1253627   # limit to specific person ids

set -uo pipefail

OUTROOT="$HOME/pipeline_outputs"
DRY=0
IDS=()
for a in "$@"; do
  case "$a" in
    --dry-run) DRY=1 ;;
    -*) echo "unknown flag: $a" >&2; exit 1 ;;
    *) IDS+=("$a") ;;
  esac
done

# Small-file types we refuse to delete regardless of size (defense-in-depth; none are >1MB anyway).
PROTECT_REGEX='\.(csv|md|txt|json|tsv|log)$|expd\.(done|attempts|minimap2\.done|minimap2\.attempts)$'

before=$(df -h "$OUTROOT" | awk 'NR==2{print $4" free / "$2" ("$5" used)"}')

# Build the list of person dirs to process.
if [ "${#IDS[@]}" -gt 0 ]; then
  dirs=(); for id in "${IDS[@]}"; do [ -d "$OUTROOT/$id" ] && dirs+=("$OUTROOT/$id"); done
else
  dirs=(); for d in "$OUTROOT"/[0-9]*/; do [ -d "$d" ] && dirs+=("$d"); done
fi

total_bytes=0; total_files=0
for d in "${dirs[@]}"; do
  # Regenerable slices/FASTQs by name, plus any >1MB file, minus the protected small types.
  while IFS= read -r -d '' f; do
    sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
    if echo "$f" | grep -Eq "$PROTECT_REGEX"; then continue; fi
    total_bytes=$((total_bytes + sz)); total_files=$((total_files + 1))
    if [ "$DRY" -eq 1 ]; then
      echo "would delete: $(numfmt --to=iec "$sz" 2>/dev/null || echo "${sz}B")  $f"
    else
      rm -f "$f"
    fi
  done < <(find "$d" -type f \( \
              -name 'expd_*_sliced.bam' -o -name 'expd_sr_namesorted.bam' \
              -o -name 'expd_*.fastq.gz' -o -name 'expd_LR.fastq' \
              -o -size +1M \) -print0)
done

human=$(numfmt --to=iec "$total_bytes" 2>/dev/null || echo "${total_bytes}B")
echo ""
if [ "$DRY" -eq 1 ]; then
  echo "DRY RUN: would delete $total_files files, freeing ~$human across ${#dirs[@]} person dir(s)."
  echo "Disk now: $before"
else
  after=$(df -h "$OUTROOT" | awk 'NR==2{print $4" free / "$2" ("$5" used)"}')
  echo "Pruned $total_files files, freed ~$human across ${#dirs[@]} person dir(s)."
  echo "Disk before: $before"
  echo "Disk after:  $after"
fi
