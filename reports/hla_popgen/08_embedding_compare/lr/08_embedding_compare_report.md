# Embedding comparison: locked-mapping allele-space views (`08_embedding_compare.py`)

3 columns (full_enriched, collapsed_nearest_ref, known_only) x 4 rows (raw/PCA-denoised UMAP x ancestry/%novel-gradient coloring). Full 12-panel grid + individual panels: `matrix/`.


## Cohort
9355 complete 8-classical-gene cases (2877 incomplete, dropped). 8893 (95.1%) carry >=1 novel allele. Dosage matrix: 9355 people x 5453 allele-identity columns.

`collapsed_nearest_ref` slot resolution: 138661 of 149680 slots (7.4% dropped) had a `nearest_allele` that IS a column in the locked fit.

`known_only`: 462 of 9355 people (whole-person drop).


PCA on `full_enriched` (9355 x 5453, standardized): PC1 0.164% var, PC2 0.155% var, cumulative top-15 1.55% var. See module docstring for why this is expected to be small at this encoding's dimensionality, not evidence of weak signal.


Figure (combined grid): `/home/jupyter/repos/pilot-validation/reports/hla_popgen/08_embedding_compare/lr/matrix/embedding_matrix.png`


## Real EHR-confirmed disease diagnosis coloring
Actual diagnosis, not allele-carrier proxy -- see `scripts/hla_popgen/12_disease_phenotypes.py`'s `HLA_LINKED` list (reused verbatim from `origin/aleix/hla-resolve-phase1`'s `deep_immune_breakdown.py`, ICD-10-3char + SNOMED-substring definitions). 853 of 9355 people matched >=1 disease definition (155 matched more than one, shown by first match). Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/08_embedding_compare/lr/disease/disease_real_diagnosis.png`

| disease | n people |
|---|---|

| Ankylosing spondylitis (HLA-B27) | 28 |
| Celiac disease (HLA-DQ2/DQ8) | 37 |
| Type 1 diabetes (HLA-DR3/DR4) | 199 |
| Psoriasis (HLA-Cw6) | 213 |
| Psoriatic arthritis (HLA-Cw6/B27) | 2 |
| Multiple sclerosis (HLA-DRB1*15:01) | 60 |
| Rheumatoid arthritis (HLA shared epitope) | 201 |
| Systemic lupus erythematosus (HLA-DR2/DR3) | 79 |
| Graves disease (HLA-DR3) | 5 |
| Narcolepsy (HLA-DQB1*06:02) | 24 |
| Behcet disease (HLA-B51) | 5 |