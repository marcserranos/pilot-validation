# 62 VM run, 2026-10-12: PUBLIC aggregates only

- `62_cole_haplotype_units.py`, run with `--exclude-people ~/s05/internal/59/exclude_people.tsv` (empty; 59 proposed no cutoff).
- Relatedness: cached LR pairs (fallback).
- The INTERNAL per-person file `haplotype_units.tsv.gz` (48,906 rows) and its MANIFEST stay on the VM. They have not been uploaded anywhere; any bucket upload is run by Marc himself (see vm_parts/62.sh).
- The VM README.md is identical to the parent folder's README.md.
- Regenerated 2026-10-12 17:14 UTC at main 5fa84b0 (PR #25: cross-row complementary suppression; dq_g1g2 not_called linked to evidence incomplete).
- Checked: no id columns; no ALL-vs-ancestry single mask; no `unit_qc_summary` partition has exactly one masked cell; no masked `n_haps` shows a frequency.
