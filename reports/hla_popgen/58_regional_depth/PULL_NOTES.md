# 58 PUBLIC tables (VM full run, 2026-10-11 21:23 → 2026-10-12 12:42 UTC)

Pulled from `~/s05/results/58/` on the VM, one file at a time. Byte sizes were checked against the VM.

How they were made:
- Script `scripts/hla_popgen/58_regional_depth.py`, run from `~/repos/pilot-validation-main` at 386d81a. It was named 50_regional_depth.py at that commit; the code is identical.
- Python and samtools 1.21 from the existing spechla pixi env.
- `--workers 8` on n2-highmem-4. Regions from GENCODE v46 basic (`regions.bed`).
- 12,233 people measured, 0 failed. 55,138 s wall time.

Disclosure:
- No ids, and no per-person rows.
- The tables were regenerated on 2026-10-12 with `--aggregate-only` at 29e17c8 (PR #23, complementary suppression), from the stored per-person depths. No new measurements.
- People counts 1–19 are masked `<20`. The smallest other cell, where needed, is masked `<20*`: 667 and 241 cells in depth_summary_by_region.
- Rates are blanked on masked CN rows.
- Checked:
  - no CN row and no histogram gene has exactly one masked cell;
  - no depth-tail cell can be recovered via revio + sequel2e = ALL, or via the ancestry sum.
- Residual: in 7 trb_cn_by_ancestry cells, ALL is unmasked and exactly one ancestry is masked. 23 people have no listed ancestry (the ancestries sum to 12,210 of 12,233), so subtraction only bounds those cells within ±23 and does not recover them exactly.

Notes for the report:
- TRBV4-3 / TRBV3-2 / TRBV6-3 / TRBV7-8 are not in the GENCODE primary-assembly gene set. They sit on an ALT-only segment. See `../65_trb_assembly_genotype/trb_probe_mapping_grch38_primary_vs_alt.tsv`.
- The excess depth over TRBV4-2 and TRBV6-2 is a stand-in for the insertion.
- `trb_rel_depth_histogram.tsv` uses 0.05 bins with the last bin open (≥1.45). It therefore cannot show a 1.5× / 2× mode split for TRBV4-2 / TRBV6-2.
