# 61 VM run, 2026-10-12 (replaces the DRAFT files in the parent folder)

- `61_kir_by_ancestry_panel.py --mode vm`, run with `--exclude-people ~/s05/internal/59/exclude_people.tsv`. 59 proposed no exclusions, so the result is byte-identical to the no-filter run except for run_summary.
- Cohort: 12,261 people with KIR calls; 11,882 unrelated (cached LR relatedness pairs); 23,755 haplotypes, of which 17,905 are complete.
- Pulled from the VM, sizes checked. No ids. People counts 1–19 are masked by the script, with complementary `<20*` cells where needed.
- Checked:
  - no ancestry partition in genes_per_hap has exactly one masked cell;
  - no gene_freq complement is in 1–19;
  - cA/cB counts are never masked one-sided.
- Note: the KIR3DL1 and KIR2DS4 rows are nearly identical (e.g. AFR 3,852 vs 3,848 haplotypes). They are not duplicates; 3DL1 and 2DS4 almost always occur on the same haplotypes (telomeric cA).
