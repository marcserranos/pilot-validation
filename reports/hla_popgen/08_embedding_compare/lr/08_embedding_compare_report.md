# Embedding comparison: locked-mapping allele-space views (`08_embedding_compare.py`)

3 columns (full_enriched, collapsed_nearest_ref, known_only) x 4 rows (raw/PCA-denoised UMAP x ancestry/%novel-gradient coloring). Full 12-panel grid + individual panels: `matrix/`.


## Cohort
9355 complete 8-classical-gene cases (2877 incomplete, dropped). 8893 (95.1%) carry >=1 novel allele. Dosage matrix: 9355 people x 5453 allele-identity columns.

`collapsed_nearest_ref` slot resolution: 138661 of 149680 slots (7.4% dropped) had a `nearest_allele` that IS a column in the locked fit.

`known_only`: 462 of 9355 people (whole-person drop).


PCA on `full_enriched` (9355 x 5453, standardized): PC1 0.164% var, PC2 0.155% var, cumulative top-15 1.55% var. See module docstring for why this is expected to be small at this encoding's dimensionality, not evidence of weak signal.


Figure (combined grid): `/home/jupyter/repos/pilot-validation/reports/hla_popgen/08_embedding_compare/lr/matrix/embedding_matrix.png`


## Sanity check: independent fits vs locked projections
Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/08_embedding_compare/lr/sanity_check/sanity_check.png`
Compare visually against `matrix/panel_collapsed_nearest_ref_raw_ancestry.png` and `matrix/panel_known_only_raw_ancestry.png` -- if these independent fits also resemble full_enriched, the locked-mapping similarity is real, not a projection artifact.


## HLA-disease-allele carrier coloring
Curated, exploratory list (12 groups) -- see `DISEASE_ALLELE_GROUPS` in the script for caveats. 6126 people matched at least one group (2393 matched more than one, shown by first match). Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/08_embedding_compare/lr/disease/disease_alleles.png`

| disease/allele group | n carriers |
|---|---|

| B*27 (spondyloarthritis) | 423 |
| B*57:01 (abacavir hypersensitivity) | 350 |
| B*58:01 (allopurinol SJS/TEN) | 473 |
| B*15:02 (carbamazepine SJS/TEN) | 114 |
| DRB1*15:01 (multiple sclerosis) | 648 |
| DRB1*04 (rheumatoid arthritis) | 72 |
| DQB1*06:02 (narcolepsy) | 864 |
| DQB1*03:02 / DQ8 (T1D) | 1229 |
| DQA1*05:01 / DQ2 (celiac) | 798 |
| B*51 (Behcet's disease) | 421 |
| C*06:02 (psoriasis) | 549 |
| A*29:02 (birdshot chorioretinopathy) | 185 |

## Collapse-radius: LR-2field vs SR-2field, SAME 6668 people
`top10_share` = fraction of allele copies in the 10 most common alleles per gene; `simpson_concentration` = sum(p_i^2).

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