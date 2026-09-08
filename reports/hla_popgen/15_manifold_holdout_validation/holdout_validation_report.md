# Held-out validation of the supervised-UMAP disease signal (`15_manifold_holdout_validation.py`)

Tackles hypothesis #1 from `research/ANCESTRY_VS_DISEASE_MANIFOLD.md`: is `14_manifold_structure.py`'s supervised-UMAP result (obs=0.97 vs null=0.83 same-label KNN fraction) real generalizing signal, or an artifact of `umap.UMAP(y=...)` being fit and scored on the same people? See this script's module docstring for the full design.


## Split
6548 train / 2807 test (9.1% / 9.1% diagnosed -- stratified split, prevalence preserved in both halves).


## supervised (train-fit, held-out test scored)
Held-out AUROC: **0.482** (bootstrap 95% CI 0.460-0.508, 1000 resamples).
Label-permutation null (test labels shuffled, scores fixed, 2000 perms): null mean 0.500 +/- 0.013, **p=0.931**.


## unsupervised (train-fit, held-out test scored)
Held-out AUROC: **0.528** (bootstrap 95% CI 0.494-0.563, 1000 resamples).
Label-permutation null (test labels shuffled, scores fixed, 2000 perms): null mean 0.500 +/- 0.018, **p=0.06097**.


## Where held-out people land
Small dots are train people (colored by true diagnosis, used to fit/supervise the mapping); triangles are held-out test people `.transform()`-projected into that frozen mapping from their allele dosage alone -- their color was never seen by the fit or the transform. Figure: `/home/jupyter/repos/pilot-validation-main/reports/hla_popgen/15_manifold_holdout_validation/holdout_train_test_embedding.png`


## Verdict
Supervised held-out AUROC 0.482 (CI 0.460-0.508) vs. unsupervised baseline 0.528 (CI 0.494-0.563).

Supervised does NOT clear its own bootstrap CI above 0.5, or does not beat the unsupervised baseline, on held-out data -- the original supervised-UMAP finding does not generalize under this test; treat it as overfit to the labels it was given, not evidence of real disease structure.
