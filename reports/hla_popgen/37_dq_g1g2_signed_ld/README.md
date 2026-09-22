# 37 -- LD supplement table (r2 / D') and DQ/DP heatmap

Source: script `37b_ld_supplement_table.py`, re-rendered from the already-committed
`reports/hla_popgen/29_hla_ld/ld_pairwise.tsv`. No VM access, no re-computation from raw calls --
every number here was already public.

Cole (call #8): "this should be a table too -- a simple CSV with r2 and D-prime", with a distinct
colour for "not observed / below n=20".

## Data-gap check (per the task)

The committed `ld_pairwise.tsv` already carries the **full per-allele-pair signed table**: for
every allele pair (per gene pair, per ancestry) where BOTH alleles individually clear the
20-haplotype floor in their own marginal, it has `r2`, `D`, `Dprime`, `freq_a`, `freq_b`,
`freq_hap`, and the suppressed joint-haplotype count `n_hap_ij_disp`. This is **not** a top-12-only
or multiallelic-summary-only table (compare `ld_multiallelic.tsv`, which IS a whole-table summary
and was not sufficient on its own) -- so this script produces the complete CSV for **all five**
committed gene pairs (A~B, B~C, DPA1~DPB1, DQA1~DQB1, DRB1~DQB1), not only the two the heatmap is
scoped to. **No gap; nothing here is deferred to the VM version.**

## `supp_table_ld_pairs.csv`

Tidy, one row per (gene_pair, ancestry, allele_a, allele_b): `gene_pair`, `ancestry`, `allele_a`,
`allele_b`, `n_hap` (the joint haplotype count, blank/NaN where suppressed), `freq_a`, `freq_b`,
`r2`, `Dprime`, `status`.

`status` is one of:
- **`estimated`** -- the joint haplotype count is disclosable (>=20 or a genuine 0) and `r2`/
  `Dprime` are estimated from real data.
- **`not_observed`** -- the joint count is exactly 0 (this allele pair never co-occurred on a
  phased haplotype in this ancestry). Disclosable; not suppressed.
- **`suppressed_lt20`** -- the joint count is 1-19 under the AoU small-cell rule, written `<20` in
  the source table. **Parsed as censored throughout -- never read as 0 and never as 20.**

## Heatmap panel (DQA1~DQB1, DPA1~DPB1)

One small heatmap per (gene pair, ancestry): rows are the gene pair's first-locus alleles,
columns the second-locus alleles, colour is signed D' on a blue-white-red diverging scale
centered at 0 (`_viz_common.diverging_cmap()`/`diverging_norm()`). Cells with status
`not_observed` or `suppressed_lt20` are drawn as a distinct light-grey hatched cell
(`_viz_common.hatch_suppressed()` / `SUPPRESSED_COLOR`) and are **excluded from the colour scale
entirely** -- they never appear as if they were a real, low-magnitude D' value.

**How to read it:** a solid block of one colour within a gene pair/ancestry panel indicates
strong, consistently-signed linkage; hatched cells mean the underlying allele pair either never
co-occurred in this ancestry's phased haplotypes or its count fell under the 20-haplotype
disclosure floor -- in both cases, no D' is estimated or shown for that cell.

## Caveats

- This is a restyle + tabulation of already-committed, already-suppressed numbers -- no new
  statistic is computed here.
- r2/D' at low allele frequency are noisy; see `29_hla_ld/README.md`'s own caveats on rare-allele
  LD estimates (which this script does not repeat or override).
- `A~B`/`B~C`/`DRB1~DQB1` are in the CSV but not the heatmap panel, which is scoped to
  DQA1~DQB1/DPA1~DPB1 per the task.
