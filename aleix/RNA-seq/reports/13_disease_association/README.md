# 13 — W4: repertoire diversity vs chronic disease, at fixed depth

Script: `scripts/13_disease_association.py` via `scripts/run_disease.sh`. Run 2026-10-12
(outputs 95f27af); figure re-rendered 2026-10-13 (b47314e). Inputs: report 09's table, 11's
labels, and 12's feature sets (VM-local, no report folder).

## How

- **Model**, for each of the 27 diseases with ≥ 50 cases and each of 6 measures (162 tests):

  `measure (SD) ~ case + age + age² + sex + PC1–10 + log10 TRB reads`

  OLS with HC3 errors. The effect is the case − control difference in SD.
- **Measures:**
  - at 823 reads: Shannon, inverse Simpson, richness, public fraction;
  - top-10 clonal share;
  - raw Shannon, for comparison with the existing entropy-wide association (raw Shannon,
    adjusted for total reads).
- **Multiple testing:** BH and Bonferroni across the main model.
- **Sensitivity model:** adds log(1 + conditions per EHR-year). People seen more often collect
  more diagnoses, and acute illness lowers circulating T cells.

## Results

33 tests have BH q < 0.05 and 12 pass Bonferroni; 12 hold BH after the utilisation
adjustment.

- **HIV (positive control)** is strongest on every measure: richness −0.87 SD, Shannon at depth
  −0.74, top-10 share +0.66. Bonferroni-significant, robust to utilisation.
- **Heart failure (373 cases):** top-10 share +0.32 SD, Shannon at depth −0.32, raw Shannon
  −0.27, richness −0.26, all Bonferroni-significant; utilisation-adjusted q 0.002–0.013. The
  fixed-depth effect is at least as large as the raw one, so it is not depth. Pattern: clonal
  expansion.
- **Hypertension** (top-10 +0.13, Shannon −0.12) and **T2D** (top-10 +0.15) are
  Bonferroni-significant in the main model, q ~0.04–0.07 after utilisation.
- Obesity, CKD and CAD are BH-nominal only and lose significance after utilisation.
- **No autoimmune disease is significant.**
- Direction note (2026-10-20 review, not a result): for inverse Simpson the immune-mediated
  diseases lean positive, i.e. fewer dominant clones. This holds for psoriasis, SLE, atopic
  dermatitis, RA, UC, allergic rhinitis, Sjögren and urticaria, while cardiometabolic diseases
  lean negative. Only psoriasis is nominal (p 0.032, q 0.13).
- Panel b: raw and fixed-depth Shannon effects fall on the identity line. Depth control does
  not change the disease conclusions here.

## Figure

`fig_disease_association`:
- a: effect heatmap, disease × measure (dot = BH q < 0.05, ring = Bonferroni).
- b: raw vs fixed-depth Shannon effects.
- c: inverse Simpson effects with 95% CI, main and utilisation-adjusted.

## Caveats

- Cross-sectional.
- No CMV serostatus. CMV drives clonal expansion and is linked to cardiovascular disease, so it
  is the main candidate confounder for the heart-failure signal.
- Treatment (e.g. immunosuppressants) is not modelled.

## Files

`summary.csv`, `association.csv`: every test, with effect, CI, p, BH q, Bonferroni flag and
the utilisation-adjusted versions.
