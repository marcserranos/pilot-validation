# Embedding comparison: locked-mapping allele-space views (`08_embedding_compare.py`)

3 columns (full_enriched, collapsed_nearest_ref, known_only) x 4 rows (raw/PCA-denoised UMAP x ancestry/%novel-gradient coloring) -- see module docstring for why `collapsed_nearest_ref` is a genuine `.transform()` projection into the `full_enriched` fit, not a second fit. Full 12-panel grid + individual panels: `matrix/`.


## Cohort
9355 complete 8-classical-gene cases (2877 incomplete, dropped). 8893 (95.1%) carry >=1 novel allele.

`collapsed_nearest_ref` slot resolution: 138661 of 149680 slots (7.4% dropped) had a `nearest_allele` that IS a column in the locked fit and could be placed; the rest had nowhere to project (a nearest reference nobody in this cohort's own known calls happens to carry) and are left at 0 for that slot.

`known_only`: 462 of 9355 people (whole-person drop, per 2026-09 chat clarification).


PCA on `full_enriched`: PC1 0.16% var, PC2 0.16% var (10 components retained for the record; the PCA-denoised UMAP row uses the top 15).


Figure (combined grid): `/home/jupyter/repos/pilot-validation/reports/hla_popgen/08_embedding_compare/lr/matrix/embedding_matrix.png`


## Collapse-radius: LR-2field vs SR-2field, SAME 6668 people
Does SR concentrate on fewer, more common alleles than LR at the SAME matched 2-field resolution? `top10_share` = fraction of allele copies in the 10 most common alleles per gene; `simpson_concentration` = sum(p_i^2).

| gene | LR top10 share | SR top10 share | LR Simpson | SR Simpson | LR n_distinct | SR n_distinct |
|---|---|---|---|---|---|---|
| A | 67.1% | 68.0% | 0.065 | 0.066 | 130 | 93 |
| B | 42.1% | 42.2% | 0.027 | 0.027 | 196 | 179 |
| C | 68.2% | 71.5% | 0.065 | 0.071 | 107 | 74 |
| DPA1 | 99.4% | 99.9% | 0.438 | 0.440 | 35 | 18 |
| DPB1 | 84.1% | 90.3% | 0.123 | 0.135 | 98 | 58 |
| DQA1 | 90.2% | 99.9% | 0.099 | 0.158 | 36 | 18 |
| DQB1 | 84.6% | 91.5% | 0.091 | 0.119 | 57 | 35 |
| DRB1 | 54.1% | 54.2% | 0.042 | 0.043 | 108 | 88 |