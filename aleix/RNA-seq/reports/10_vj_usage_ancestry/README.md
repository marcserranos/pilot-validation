# 10 — W1: TRB V and J gene usage by genetic ancestry

Script: `scripts/10_vj_usage_ancestry.py` via `scripts/run_w1.sh`. Final run 887db29
(2026-10-08, with the fixed zero-usage estimator). Figure re-rendered 2026-10-13 (b47314e).
Requested at the 2026-10-06 supervisor meeting: a model-free view of whether repertoire
composition differs by ancestry.

## How

- **Unit:** unique clonotypes per gene, not reads. Clonotypes reflect which genes recombine
  and survive selection; reads add expansion and transcription.
- **Genes:** 46 common TRBV + 13 TRBJ (mean usage ≥ 0.1%, present in ≥ 50% of people).
  13 V and 1 J symbols are dropped as rare.
- **Composition:** CLR per composition (V and J separately), +0.5 pseudocount, Aitchison PCA.
  7,527 people at the read floor.
- **Prediction:** macro one-vs-rest AUROC on unrelated people, CV grouped by family; compared
  with covariates only (age, sex, depth) and permuted labels.
- **Zero usage:**
  - observed share of people with zero clonotypes of a gene vs the share expected by chance
    at its pooled usage rate;
  - genes with ≥ 1% chance zeros are marked not assessable;
  - implied absence-allele frequency = √ of the excess-zero fraction (Hardy–Weinberg).

## Results

- **V/J usage predicts ancestry:** macro AUROC 0.804 (V 0.774, J 0.701; covariates only 0.614;
  permuted 0.504).
  - With covariates: 0.815. Without the indel genes: 0.803, so TRBV4-3 is not what carries
    the prediction.
  - Per group (V + J): AFR 0.88, EAS 0.88, SAS 0.83, EUR 0.79, MID 0.73, AMR 0.70.
- **PC1 (19.4%) is almost one gene, TRBV4-3** (loading 0.90). Ancestry explains 7.5% of PC1 and
  12% of PC2. The two clusters in panel a contain every ancestry. Without TRBV4-3 (panel e)
  they disappear.
- **Germline signatures (panel f):**
  - TRBV4-3 is never used by 27.3% of people, against 0.04% expected by chance: SAS 43.5%,
    AFR 33.3%, EAS 30.1%, MID 25.8%, AMR 25.4%, EUR 16.7%. Implied absence-allele frequency
    0.52.
  - TRBV12-4: 15.6% never use it (EAS 25.7%, EUR 17.6%, AFR 7.4%). It is near-identical to
    TRBV12-3, so this may be call assignment rather than germline.
  - Smaller: TRBV7-3 3.2%, TRBV30 2.9%, TRBV6-4 3.7%.
  - 23 rare genes are not assessable.
  - **Report 15 confirms TRBV4-3 against DNA genotypes** (Marc's script 65): same ranking in
    every group, and TRBV6-2 halves in non-users because the insertion carries a second copy.
- **Genes that differ most** (eta²): J genes first.
  - TRBJ2-1 (0.18): higher in EAS and SAS, lower in AFR.
  - TRBJ1-2, 1-1 and 1-6: higher in AFR.
  - TRBV28 is ~2× lower in AFR.
  - Germline and reference-mapping explanations are not separated.
- **Dosage (panel g):** 2,055 people with zero TRBV4-3, a near-zero shoulder (0–0.2%), then one
  broad mode near 1.2%. Heterozygote and homozygote modes are not separable at this depth.

## Figure

`fig_vj_usage_ancestry`:
- a: V + J PCA, ancestry as 50/80% density contours (purity ≥ 0.90).
- b: variance and ancestry share per PC.
- c: top genes by eta².
- d: ancestry AUROC.
- e: PCA without TRBV4-3.
- f: never-use rates against chance.
- g: TRBV4-3 dosage.

## Caveats

- Usage shares sum to one, so a gene strongly up in one group pushes every other gene down a
  little there. Read the top of the gene list, not the tail.
- TRBV deletion and loss-of-function alleles differ by population, so V-usage differences can
  be genetic without involving HLA. This report describes; it does not attribute.

## Files

`summary.csv`, `ancestry_auroc.csv`, `genes_by_ancestry.csv`, `pca_components.csv`,
`pca_components_without_indel_genes.csv`, `pca_loadings.csv`, `zero_usage_by_ancestry.csv`,
`dosage_histogram.csv`, `codeletion_vs_absence.csv`.
