# Embedding comparison: allele-space resolution vs population structure (`08_embedding_compare.py`)

See this script's module docstring for why `lr_collapsed`/`sr_collapsed` are separately-fit embeddings (not projections of `lr_full`) while `lr_known_only` IS the same fitted `lr_full` coordinates, subset.


## lr_full
Complete 8-classical-gene cases: 9355 (2877 incomplete, dropped). Novel carriers: 8893 (95.1%).

Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/08_embedding_compare/lr/pca_lr_full.png`

Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/08_embedding_compare/lr/umap_lr_full.png`


## lr_known_only (same lr_full coordinates, subset)
462 of 9355 people carry zero novel calls across these 8 genes -- plotted as the SAME PCA/UMAP points as lr_full above, filtered, never a re-fit.

Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/08_embedding_compare/lr/pca_lr_known_only_highlight.png`


## lr_collapsed
Complete 2-field cases among lr_full's 9355 people: 7059 (2296 incomplete under 2-field truncation). This drop is EXPECTED to be substantial, not a bug: a `protein_altering` novel call (depth-2 novelty, reference/IMMUANNOT_GTF_SPEC.md part B) only has its FIRST field resolved before 'new' is spliced in -- `to_nfield(consensus, 2)` needs 2 resolved fields and returns None for it. Since protein-altering is the dominant novelty class (~91% of real novel calls, 03_novel_alleles.py's report), 2-field collapsing doesn't just blur these calls into a coarser bucket -- it makes the PERSON'S ENTIRE 8-gene case incomplete and drops them from this embedding outright. That is itself a quantitative answer to 'how much does collapsing cost': compare this count directly against lr_full's zero incomplete cases.

Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/08_embedding_compare/lr/pca_lr_collapsed.png`

Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/08_embedding_compare/lr/umap_lr_collapsed.png`


## sr_collapsed
Of 7044 people shared between lr and sr: 7044 have a complete 2-field SR case (0 incomplete/no-call, dropped).

Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/08_embedding_compare/lr/pca_sr_collapsed.png`


## Collapse-radius comparison: LR-collapsed vs SR-collapsed, SAME 7044 people, both at 2-field resolution
Hypothesis (2026-09 chat discussion): SR, having less intrinsic resolution/confidence, calls concentrate more heavily on a few common alleles than LR does even at the SAME nominal 2-field resolution. `top10_share` = fraction of all allele copies covered by the 10 most common alleles for that gene; `simpson_concentration` = sum(p_i^2) (higher = more concentrated on few alleles).

| gene | LR top10 share | SR top10 share | LR Simpson | SR Simpson | LR n_distinct | SR n_distinct |
|---|---|---|---|---|---|---|
| A | 67.1% | 68.0% | 0.064 | 0.065 | 131 | 93 |
| B | 42.0% | 42.1% | 0.027 | 0.027 | 198 | 180 |
| C | 68.2% | 71.5% | 0.064 | 0.071 | 108 | 75 |
| DPA1 | 99.4% | 99.9% | 0.436 | 0.439 | 35 | 18 |
| DPB1 | 84.1% | 90.2% | 0.122 | 0.135 | 98 | 58 |
| DQA1 | 90.2% | 99.9% | 0.099 | 0.158 | 37 | 19 |
| DQB1 | 84.5% | 91.4% | 0.091 | 0.119 | 58 | 35 |
| DRB1 | 54.0% | 54.1% | 0.042 | 0.043 | 109 | 89 |