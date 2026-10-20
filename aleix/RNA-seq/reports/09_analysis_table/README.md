# 09 — W0: one analysis table, one row per person

Script: `scripts/09_analysis_table.py` via `scripts/run_w0.sh`. Run 2026-10-08; final run with
sex as a covariate 88d08dd. Figure re-rendered in the journal font 2026-10-13 (b47314e).

## Purpose

Every downstream analysis (10 V/J usage by ancestry, 13 disease association, 14
classification) starts from this table. Inclusion rules, covariates and depth-robust
measures are defined once.

The table itself (7,922 × 72, real research_ids) stays on the VM:
`~/pipeline_outputs/rnaseq/analysis/person_table.{tsv.gz,pkl}`. This folder holds aggregates
only.

## Contents of the table (`data_dictionary.csv`)

- **identity:** ancestry, max ancestry probability, purity flags, AoU PC1–10.
- **demographics:** age (2020 − year of birth), sex at birth.
- **EHR:** observation window, distinct conditions, conditions per EHR-year.
- **technical:** TRB reads and clonotypes, all-chain CDR3 counts, RNA-SeQC2 metrics,
  `rnaseq_metadata.tsv` fields.
- **repertoire:**
  - raw Shannon, inverse Simpson, top-10 share, public fraction;
  - at a fixed 823 reads (mean of 10 seeded subsamples): Hurlbert richness, Shannon, e^H,
    inverse Simpson, public fraction.
- **relatedness:** family and unrelated flag, resolved inside this cohort from AoU's pairwise
  kinship at ≥ 0.0884.
- **flags:** read floor, EHR window, adult, `analysis_main`.

## Results

**Waterfall** (cumulative, `waterfall.csv`):

| Step | All | AFR | AMR | EAS | EUR | MID | SAS |
|---|---|---|---|---|---|---|---|
| cohort | 7,922 | 1,288 | 1,699 | 1,134 | 2,351 | 404 | 1,046 |
| TRB reads ≥ 823 | 7,527 | 1,236 | 1,592 | 1,078 | 2,203 | 392 | 1,026 |
| EHR window > 0 | 6,971 | 1,058 | 1,419 | 1,021 | 2,101 | 375 | 997 |
| adult | 6,968 | 1,058 | 1,418 | 1,019 | 2,101 | 375 | 997 |
| unrelated = **analysis_main** | **6,764** | 1,018 | 1,362 | 1,006 | 2,023 | 371 | 984 |

- The read floor costs 5% evenly. The EHR window costs AFR 14% and AMR 11% vs EUR and EAS 5%.
- Sex is a covariate (Female / Male / Unknown), not a filter: requiring it removed 13% of AFR
  vs 1.5% of EAS.
- Relatedness: 445 pairs, 571 people, 206 families. 275 dropped, keeping the relative with more
  TRB reads.
- **Depth control works.** Spearman with log TRB reads:

  | Measure | Raw | At 823 reads |
  |---|---|---|
  | Shannon | 0.56 | 0.15 |
  | inverse Simpson | 0.13 | 0.02 |
  | richness (clonotype count) | 0.90 | 0.30 |
  | public fraction | ~0 | ~0 |

- Ancestry probabilities: argmax reproduces `ancestry_pred` for 100% of people. Purity is
  uneven: ≥ 0.90 keeps 31% of EUR, 39% of AMR, 76% of EAS and 88% of SAS. Use it for coloured
  figures only.
- `rnaseq_metadata.tsv` has no draw date (`metadata_columns.csv`), so age stays 2020 − year of
  birth.

## Figure

`fig_depth_and_diversity`:
- a: raw Shannon vs log TRB reads.
- b: Shannon at 823 reads vs log TRB reads.
- c: share of people kept at each read floor, by ancestry.

Raising the floor to 2,000 would keep 48–63% and hit EUR hardest.

## Files

`summary.csv`, `waterfall.csv`, `waterfall_by_floor.csv`, `column_summary.csv`,
`data_dictionary.csv`, `metadata_columns.csv`.
