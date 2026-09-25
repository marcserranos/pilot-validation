# 47 — Naive ML / simple tests on allele-carriage matrices

*Status: script + synthetic tests done locally (S04 WS-D, 2026-09-25/26). Not yet run on real
AoU data — needs a VM session. This README is a skeleton to be filled with real numbers once
that run happens; see `STATUS.txt`/the metrics TSVs in this folder for what actually shipped.*

## Question

Do simple, transparent classifiers (L1-logistic regression, a depth-3 decision tree) recover
known structure from HLA + KIR allele-carriage alone — ancestry, sequencing platform, and
HLA-genotype-predicts-KIR-content — using nothing more exotic than 5-fold cross-validated AUROC
and a label-permutation null? These are sanity checks and confound warnings, not disease-risk
claims (see `sprints/S04_kir_recurrence_style_share/WS-D_EXPLAINER_AND_IDEAS.md` for the
plain-language framing and why each expected result matters).

## Method

- **Features:** HLA two-field allele carriage (`GENE*allele`, classical genes only, alleles with
  <20 unrelated carriers dropped before modeling) + KIR gene presence/absence (either haplotype,
  17-gene panel). Unrelated subset only (kin < 0.0442, `41_kir_pilot.py`'s `greedy_unrelated`).
- **Tasks:** ancestry one-vs-rest (6 binary models, AFR/AMR/EAS/EUR/MID/SAS), sequencing platform
  one-vs-rest (`--platform-col`, values expected `{revio, sequel2e, sequel2}` per SCHEMA.md Table
  4 / 41's `ASSEMBLY_PLATFORMS` — **verify on the VM instance actually used**, see caveats), and
  person-level KIR cA/cB content predicted from HLA carriage alone (a near-negative control:
  HLA/chr6 and KIR/chr19 have no known direct linkage).
- **Models:** L1-logistic regression (`liblinear`) and a depth-3 decision tree, 5-fold stratified
  CV, pooled out-of-fold AUROC. A task is skipped (not reported as 0) if either class has <20
  people.
- **Permutation baseline:** >=20 label shuffles (`--n-perms`, default 20 locally / raise for the
  real run), refit, empirical p-value against the shuffled-null AUROC distribution.
- **Top features:** |coefficient| (LR) / feature_importance_ (tree), restricted structurally to
  the already->=20-carrier feature set — never an allele name paired with a small carrier count.
- **PCA:** 2 components on the same carriage matrix. Per-person coordinates never leave the VM —
  only per-(ancestry x platform) centroids (n>=20 groups only) and a 2D binned density grid
  (cells 1-19 -> `<20`, true 0 stays `0`) are exported.

## Result

*(fill in after the VM run: `naive_ml_metrics.tsv`, `naive_ml_top_features.tsv`,
`pca_group_centroids.tsv`, `pca_density_grid.tsv`)*

## How to read each output file

- `naive_ml_metrics.tsv` — one row per task: `lr_auroc`/`tree_auroc` (0.5 = chance, 1.0 =
  perfect), `perm_auroc_mean`/`perm_auroc_p95` (the null distribution — a real AUROC should sit
  above `perm_auroc_p95`), `perm_pvalue_lr`.
- `naive_ml_top_features.tsv` — which alleles/genes drove each model, per task.
- `pca_group_centroids.tsv` — mean PC1/PC2 per ancestry x platform group (n>=20 only); platform
  clustering here would flag a batch effect, not biology.
- `pca_density_grid.tsv` — 2D histogram of PC1/PC2 over all people, binned and count-suppressed.

## Caveats

- **`--platform-col` needs VM verification.** `cohort_membership.tsv`'s `platform` column is
  documented in SCHEMA.md as sourced from `immuannot_cohort_full.tsv`, but the sibling script
  `42_repertoire_baseline.py` recorded a VM instance where an expected HLA column was entirely
  absent. This script exits loudly (printing available columns) rather than guessing if
  `--platform-col` doesn't match.
- Local dry run used `--n-perms 10-20` for speed; the real run should raise this (>=20 minimum
  per the WS-D brief, more if VM time allows).
- The cA/cB task uses a coarse person-level label (any B-content gene on either haplotype); it
  does not distinguish AA/AB/BB genotype dosage.
- Ancestry and platform predictability are expected and are not, on their own, evidence of a
  problem — they are the confound warnings this script exists to surface.

## Distilled

- Not yet run on real data. Script + 13 synthetic unit tests pass locally
  (`scripts/hla_popgen/tests/test_47_naive_ml_tests.py`), including a planted-signal test (AUROC
  clears its own permutation null) and an underpowered-class guard (n<20 skips the task).
- Next step: one VM run with `--platform-col` verified against `cohort_membership.tsv`'s actual
  columns, then fill this README's Result section.
