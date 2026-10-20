# 06 — Best person-level representation of a TRB repertoire

Script: `scripts/06_person_representation.py`. Current outputs from the 2026-10-04 polish run
(5b83486). Vectors come from the full-cache embedding: 98.9% of clonotypes have a vector.

## Question

Which way of turning a person's repertoire into a vector captures stable individual identity
and biology, rather than sequencing depth?

## How

- 4,000 people, 4.72 M clonotypes, 16 representations:
  - TRBV usage; CDR3-length profile; CDR3 3-mers;
  - mean and read-weighted mean SCEPTR vector;
  - SCEPTR metacluster profiles (k = 30, 100, 300), raw and CLR;
  - concatenations, raw and block-balanced.
- Axes:
  1. **Identifiability:** split each person's clonotypes into two random halves; is their own
     second half the nearest of 4,000? Chance 0.025%. Given raw (cosine), z-scored, and
     z-scored with log depth removed, with 95% CIs.
  2. **Signal:** out-of-fold age R² and MAE, sex AUROC.
  3. **Nuisance:** out-of-fold log-depth R² (lower is better).
  4. **Ancestry:** macro AUROC, reported as a property, not a target.
- A label-permuted control runs through the same pipeline.

## Results

| Representation | Identifiability raw / z-scored | Age R² | Depth R² |
|---|---|---|---|
| SCEPTR mean + TRBV | **29.4% / 27.0%** (CI 25.6–28.3 z) | **0.34** | 0.17 |
| TRBV usage | 23.8% / 16.5% | 0.26 | 0.18 |
| SCEPTR clusters k=100 | 13.2% / 12.6% | 0.26 | 0.14 |
| SCEPTR mean | 11.2% / 11.5% | 0.28 | 0.12 |
| SCEPTR clusters k=300 | 7.4% / 15.6% | 0.24 | 0.07 |
| CDR3 3-mer | 4.7% / 2.4% | 0.19 | 0.05 |
| SCEPTR mean, read-weighted | 1.6% / 1.5% | 0.05 | 0.04 |
| Permuted labels | 0.025% | −0.02 | — |

- **Headline (DECISIONS 2026-10-06): z-scored 27.0%.** Removing log depth leaves 26.9%, so
  identity is not depth.
- Mean + TRBV also gives the best sex (0.67) and ancestry (0.79) AUROC.
- Read-weighting and CLR both hurt badly. Expansion dominates the weighted vectors.
- Cluster profiles, predicted to win, came below plain TRBV usage.

## Figure

`fig_person_representation`:
- a: identifiability (bars raw; dots z-scored and depth-removed).
- b: age R².
- c: signal vs nuisance.
- d: identifiability against cluster count.
- e: sex and ancestry AUROC.
- f: all metrics z-scored across representations.

Known layout issue: d's y-axis label overlaps c's legend.

## Files

`summary.csv`, `representation_benchmark.csv` (every metric and CI per representation),
`permuted_control.csv`.
