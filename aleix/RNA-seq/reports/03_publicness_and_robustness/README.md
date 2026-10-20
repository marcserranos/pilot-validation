# 03 — Why TCRs are public, and how robust the person vector is

Script: `scripts/03_publicness_and_robustness.py`. Current outputs from the 2026-10-04 polish
run (5b83486), after the random tie-break fix (report 08).

## Questions

1. Is sharing between people real convergent recombination? If so, shared clonotypes must
   have higher generation probability (Pgen) and shorter CDR3s (Venturi 2006; Elhanati 2018;
   OLGA, Sethna 2019).
2. Is report 01's person vector (mean SCEPTR vector of the top 500) stable under the cap and
   under read-weighting?
3. How much genetic ancestry does a person's repertoire carry?

## How

- Pgen from OLGA's default human TRB model on ~5,000 clonotypes per sharing level (19,006
  total; 0.03% have Pgen = 0).
- Robustness: 1,500 people. Recompute the person vector from the top k clonotypes
  (k = 25 … 500) and read-weighted. Compare with the top-500 unweighted vector by Spearman
  correlation of all pairwise distances, and by overlap of each person's 10 nearest neighbours.
- Ancestry: 5-fold cross-validated one-vs-rest AUROC, SCEPTR person vector vs TRBV usage vs
  permuted labels.

## Results

- **Sharing follows recombination likelihood.** Spearman ρ = 0.66 between people sharing and
  log10 Pgen. Median log10 Pgen −9.3 for private clonotypes, −6.6 for those in ≥ 20 people.
  Median CDR3 length 15 aa private, 12 aa in 20–99 people.
- **Person vector:**
  - Top 250 vs top 500: distance correlation 0.76, 10-NN overlap 19%. Before the tie fix this
    was 0.35 and 7.9%, so the earlier instability was mostly the alphabetical tie-break.
  - Read-weighting: 0.23 and 8.8%. It changes the representation substantially.
- **Ancestry:** macro AUROC 0.73 (SCEPTR vector) vs 0.74 (TRBV usage) vs 0.50 (permuted). The
  SCEPTR vector adds nothing beyond V-gene usage. Report 10 shows much of that is germline.

## Figure

`fig_publicness_and_robustness`:
- a: share of clonotypes in ≥ 2 people by CDR3 length.
- b: Pgen by sharing level.
- c: CDR3 length by sharing level.
- d: robustness against k.
- e: per-ancestry AUROC.

Known layout issue: the x-axis labels of b and c collide ("5–19", "20–99"). Fixing it needs a
plot-only path or a VM rerun.

## Caveats

- Sharing is computed within this cohort at amino-acid level (TRBV + CDR3 aa).
- Pgen uses OLGA's default model, which is not fitted to this cohort.

## Files

`summary.csv`, `public_fraction_by_length.csv`, `pgen_and_length_by_sharing.csv`,
`person_vector_robustness.csv`, `ancestry_auroc.csv`.
