# Population-structure figures -- cohort `lr (Immuannot long-read, full)`

N people (cohort membership) = 12233.

Allele-dosage matrix: 7038 complete-case people (5194 excluded for incomplete 8-gene calls), 7023 with a usable ancestry label, 775 (gene, 2-field allele) features.

PCA: PC1+PC2 explain 1.5% of variance. Silhouette (ancestry labels, PC1-2): 0.041.

UMAP silhouette (ancestry labels): -0.032.


## Pairwise Hudson Fst (averaged over classical genes)

| | AFR | AMR | EAS | EUR | MID | SAS |
|---|---|---|---|---|---|---|
| AFR | 0.000 | 0.032 | 0.057 | 0.048 | 0.049 | 0.047 |
| AMR | 0.032 | 0.000 | 0.063 | 0.010 | 0.013 | 0.023 |
| EAS | 0.057 | 0.063 | 0.000 | 0.076 | 0.077 | 0.059 |
| EUR | 0.048 | 0.010 | 0.076 | 0.000 | 0.007 | 0.020 |
| MID | 0.049 | 0.013 | 0.077 | 0.007 | 0.000 | 0.021 |
| SAS | 0.047 | 0.023 | 0.059 | 0.020 | 0.021 | 0.000 |

## Continuous-ancestry figures (gene B, allele B*07:02, most common in this cohort)

- 1359 carriers of 12210 people with full admixture data.


## Figures

- `reports/hla_popgen/06_figures_structure/lr_full/pca_pc1_pc2_by_ancestry.png`
- `reports/hla_popgen/06_figures_structure/lr_full/umap_by_ancestry.png`
- `reports/hla_popgen/06_figures_structure/lr_full/fst_dendrogram.png`
- `reports/hla_popgen/06_figures_structure/lr_full/continuum_scatter_B.png`
- `reports/hla_popgen/06_figures_structure/lr_full/ternary_B.png`
- `reports/hla_popgen/06_figures_structure/lr_full/admixture_barcode.png`