# 14 — W5: disease classification v0 from the TRB repertoire

Script: `scripts/14_disease_classification.py` via `scripts/run_disease.sh`. Run 2026-10-13
(9,697 s; outputs 95f27af); figure redrawn from `classification.csv` (`--plot-only`, e4926bb).
Requested at the 2026-10-06 supervisor meeting: TRB only, SCEPTR and k-mer features, XGBoost
and LASSO, a working first version.

## Question

Age alone predicts most chronic disease. What does the repertoire add to age, sex, ancestry
and depth?

## How

- 15 diseases with ≥ 150 cases (report 11).
- **Base features:** age, age², sex, PC1–10, log10 TRB reads.
- **Added families** (from 12): diversity at 823 reads + top-10 share; V/J CLR usage; SCEPTR
  mean + TRBV; CDR3 3-mers (100-d SVD); all.
- **Models:**
  - LASSO: L1 logistic, penalty by inner 3-fold CV;
  - XGBoost: fixed conservative settings, not tuned in v0.
- **Evaluation:** 5-fold stratified CV, same folds for every feature set. Gain = AUROC(base +
  family) − AUROC(base), with a paired bootstrap 95% CI.

## Results

- **Base AUROC** is 0.72–0.84 for cardiometabolic disease (PAD 0.84, AF 0.83, CAD 0.81,
  hypertension 0.81, hyperlipidaemia 0.80, HF 0.79, CKD 0.78, COPD 0.78, T2D 0.73). Obesity
  0.65, urticaria 0.64, asthma 0.60, allergic rhinitis 0.58.
- **The repertoire adds at most ~1 AUROC point.** 4 of 150 gains have a 95% CI above 0:
  - CKD +1.3 (XGBoost, diversity);
  - HF +0.7 (LASSO, diversity);
  - T2D +0.5 (LASSO, diversity);
  - hypertension +0.5 (XGBoost, V/J).

  That is about what chance gives (~3.75 at a one-sided 2.5%). Only gains that agree with
  W4 mean anything: heart failure, weakly T2D and hypertension.
- Diversity is the only family that ever helps. SCEPTR, 3-mers and V/J add nothing or hurt.

## Known design problem (open, DECISIONS.md 2026-10-16)

- LASSO puts one L1 penalty on base and repertoire columns together. With hundreds of
  repertoire columns, the inner CV picks a strong penalty that also shrinks age, sex and the
  PCs.
  - So "+ all" lands below base: asthma 0.602 → 0.528, allergic rhinitis 0.575 → 0.525.
  - XGBoost "+ all" also drops (asthma 0.591 → 0.545).
- The negative gains measure the cost of adding columns, not harm from the repertoire.
  Positive gains may be underestimated.
- Fix before quoting the numbers: fit the base once, then fit the repertoire with the base
  logit as a fixed offset.
- Supervisor asks for v1: add B-cell chains and all features combined, and frame small gains as
  an assay depth limit (~2,000 TRB reads) rather than absence of signal. That framing is
  currently an interpretation, not a tested result.

## Figure

`fig_disease_classification`:
- a: out-of-fold AUROC, base vs base + all (LASSO).
- b: gain over base for every disease × family × model. Circles mark CIs above 0.

## Files

`summary.csv`, `classification.csv`: every disease × model × feature set, with AUROC and CI,
gain and CI.
