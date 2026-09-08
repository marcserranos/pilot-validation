# Ancestry vs. disease manifold structure (`14_manifold_structure.py`)

See `research/ANCESTRY_VS_DISEASE_MANIFOLD.md` for the full hypothesis and experiment checklist this implements.

Cohort: 9355 people x 5453 allele-identity columns (same `full_enriched` matrix as `08_embedding_compare.py`). 853 (9.1%) matched >=1 of 11 HLA-linked diagnoses.


## 1. Variance / loadings audit
Rank correlation (|PC1 loading| vs. ancestry F-stat, all 5453 columns): **0.395**.

Disease-proxy columns (145): mean |PC1 loading| 0.0092, mean ancestry F-stat 10.01.
All other columns (5308): mean |PC1 loading| 0.0066, mean ancestry F-stat 6.25.

Full per-column table: `/home/jupyter/repos/pilot-validation-main/reports/hla_popgen/14_manifold_structure/variance_loadings.tsv`


## 2. Ancestry-residualized re-embedding
Each column's per-ancestry-group mean subtracted before re-fitting PCA + UMAP from scratch (removes between-ancestry variance, keeps within-ancestry variance). Compare `resid_raw` panel below against `full_enriched`'s `raw_ancestry` panel in `08_embedding_compare`'s matrix figure.


## 3. Supervised UMAP (target = any HLA-linked diagnosis)
`umap-learn` fit with `y=any_diag` -- the sharpest test of whether any combination of the allele-dosage space separates diagnosed from undiagnosed people.


Figure: `/home/jupyter/repos/pilot-validation-main/reports/hla_popgen/14_manifold_structure/manifold_structure_embeddings.png`


## 4. KNN-label-enrichment permutation test (k=15, 500 permutations)
Mean fraction of a point's k nearest embedding-space neighbors sharing its label, vs. a label-permutation null (neighbor graph fixed). p-value is one-sided (observed >= null).

| embedding | label set | observed | null mean +/- sd | p |
|---|---|---|---|---|
| locked_raw | ancestry (positive control) | 0.5291 | 0.2049 +/- 0.0015 | 0.001996 |
| locked_raw | any_HLA_linked_diagnosis | 0.8339 | 0.8342 +/- 0.0010 | 0.6427 |
| ancestry_residualized | ancestry (positive control) | 0.4946 | 0.2047 +/- 0.0016 | 0.001996 |
| ancestry_residualized | any_HLA_linked_diagnosis | 0.8337 | 0.8342 +/- 0.0009 | 0.7385 |
| supervised | ancestry (positive control) | 0.5159 | 0.2049 +/- 0.0014 | 0.001996 |
| supervised | any_HLA_linked_diagnosis | 0.9703 | 0.8343 +/- 0.0009 | 0.001996 |

Full table: `/home/jupyter/repos/pilot-validation-main/reports/hla_popgen/14_manifold_structure/knn_label_enrichment.tsv`
