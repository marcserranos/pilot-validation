# 47 — Naive ML / simple tests on allele-carriage matrices

*Status: run on the full cohort 2026-09-25/26 (S04 WS-D VM session). `--platform-col platform`
verified against `cohort_membership.tsv` (values `revio`/`sequel2e`/`sequel2`, counts 11,042/
1,219/991). See "Result" below for the real numbers from that run. Figure added 2026-09-26:
`fig_naive_ml_auroc.png` (+ `.pdf`), from `47b_naive_ml_figure.py`.*

*2026-09-26 update (local-only change, not yet re-run on the VM): `47_naive_ml_tests.py` now
stores the full permutation-null distribution per task (min/5/25/50/75/95th percentile/max,
default `--n-perms` raised to 200, parallelized with a lighter null model), and adds
`platform_<X>_EUR`/`_AFR`/`_adj_ancestry` and `cB_from_HLA_EUR`/`_adj_ancestry` tasks
(`naive_ml_ancestry_adjusted.tsv`) that directly answer "Is platform confounded by ancestry?"
below. The "Result" numbers below are from the pre-upgrade run and do not yet include these new
columns/tasks — a VM rerun is needed before quoting a settled answer to that question.*

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
- **Permutation baseline:** `--n-perms` label shuffles (default 200), parallelized across
  processes (fork-based, no per-shuffle pickling of the feature matrix) using a *lighter* null
  model than the real fit -- 3-fold CV instead of 5, and a capped solver iteration count -- since
  the null's location near 0.5 doesn't depend on fold count or tight convergence. Every task's
  full null distribution is stored: min, 5th/25th/50th/75th/95th percentile, max, and the actual
  `n_perms` completed (a degenerate single-class fold under permutation is dropped, not padded).
  Empirical p-value against the shuffled-null AUROC distribution is still reported but should not
  be over-read (see "Honest read of the p-values" below).
- **Batch-effect-vs-ancestry-echo check (platform, and cB_from_HLA):** for each platform value and
  for cB_from_HLA, three additional views settle whether the raw AUROC is a genuine effect or
  ancestry leaking through the carriage features: (a) `<task>_EUR` / `<task>_AFR` -- the same
  classification restricted to one ancestry, requiring >=20 people per class **per CV fold** (a
  stricter gate than the plain >=20-total rule used elsewhere); (b) `<task>_adj_ancestry` in
  `naive_ml_ancestry_adjusted.tsv` -- an L1-LR fit on carriage + the 6 ancestry probabilities
  (`p_afr..p_sas`, SCHEMA.md Table 4) vs. the same probabilities alone, on identical CV folds,
  reporting `delta_auroc = auroc(carriage+ancestry) - auroc(ancestry-only)`. A delta near 0 means
  the raw AUROC was mostly an ancestry echo; a delta clearly above 0 means carriage predicts the
  outcome independently of ancestry. These ancestry-probability columns are soft-skipped (not
  fatal) with a log message if a VM instance's `cohort_membership.tsv` lacks them.
- **Top features:** |coefficient| (LR) / feature_importance_ (tree), restricted structurally to
  the already->=20-carrier feature set — never an allele name paired with a small carrier count.
- **PCA:** 2 components on the same carriage matrix. Per-person coordinates never leave the VM —
  only per-(ancestry x platform) centroids (n>=20 groups only) and a 2D binned density grid
  (cells 1-19 -> `<20`, true 0 stays `0`) are exported.

## Result

Full cohort: 11,845 people, 303 features (286 HLA two-field alleles with >=20 carriers, 17 KIR
genes), 620 rare HLA alleles dropped before modeling, 9 tasks run, 20 permutations/task.

| task | n_total | n_pos | LR AUROC | tree AUROC | perm null mean | perm p |
|---|---|---|---|---|---|---|
| ancestry_AFR | 11,845 | 3,009 | 0.950 | 0.782 | 0.499 | 0.048 |
| ancestry_AMR | 11,845 | 2,656 | 0.810 | 0.652 | 0.500 | 0.048 |
| ancestry_EAS | 11,845 | 1,460 | 0.981 | 0.813 | 0.500 | 0.048 |
| ancestry_EUR | 11,845 | 2,974 | 0.897 | 0.712 | 0.499 | 0.048 |
| ancestry_MID | 11,845 | 487 | 0.897 | 0.628 | 0.504 | 0.048 |
| ancestry_SAS | 11,845 | 1,236 | 0.966 | 0.671 | 0.496 | 0.048 |
| platform_revio | 11,845 | 10,710 | 0.573 | 0.507 | 0.501 | 0.048 |
| platform_sequel2e | 11,845 | 1,135 | 0.573 | 0.507 | 0.501 | 0.048 |
| cB_from_HLA | 11,845 | 7,978 | 0.550 | 0.526 | 0.501 | 0.048 |

All permutation p-values sit at the 20-shuffle floor (`1/21 ≈ 0.048`) for every task, including
the near-chance ones — a real effect size, not just "distinguishable from the null", is what
separates ancestry (AUROC 0.81-0.98) from platform/cB_from_HLA (AUROC 0.55-0.57): those two are
well above 0.5 but far from the ancestry tier, consistent with a modest real signal (platform: a
technical/batch correlate riding on the calls; cB_from_HLA: plausibly an ancestry-mediated
correlation rather than direct HLA-KIR linkage, since chr6/chr19 have no known direct linkage) —
not proof of a confound-free result. Top AFR-ancestry-informative alleles (L1 coefficients):
HLA-B\*07:05, HLA-C\*16:02, HLA-B\*38:02 (negative), HLA-DPB1\*18:01, HLA-A\*80:01 (positive) — see
`naive_ml_top_features.tsv` for every task's full list.

`pca_group_centroids.tsv`: per-(ancestry x platform) PC1/PC2 means, all groups n>=20 (one
ancestry=<NA>/revio group at n=23). `pca_density_grid.tsv`: 20x20 binned density grid, 126/400
cells suppressed to `<20`, 0 cells stayed `0` (never conflated).

## How to read each output file

- `naive_ml_metrics.tsv` — one row per task (now including `<task>_EUR`/`_AFR` subgroup variants):
  `lr_auroc`/`tree_auroc` (0.5 = chance, 1.0 = perfect), `n_perms` (permutation draws actually
  completed), `perm_auroc_mean` and the full null distribution
  (`perm_auroc_min`/`_p05`/`_p25`/`_p50`/`_p75`/`_p95`/`_max` — a real AUROC should sit above
  `perm_auroc_p95`), `perm_pvalue_lr`.
- `naive_ml_ancestry_adjusted.tsv` — one row per `<task>_adj_ancestry` (currently platform per
  platform value, and `cB_from_HLA`): `auroc_carriage_plus_ancestry` (carriage features + the 6
  ancestry probabilities), `auroc_ancestry_only` (the 6 probabilities alone, same CV folds), and
  `delta_auroc` = the difference — near 0 means the raw task's AUROC was mostly an ancestry echo.
- `naive_ml_top_features.tsv` — which alleles/genes drove each model, per task.
- `pca_group_centroids.tsv` — mean PC1/PC2 per ancestry x platform group (n>=20 only); platform
  clustering here would flag a batch effect, not biology.
- `pca_density_grid.tsv` — 2D histogram of PC1/PC2 over all people, binned and count-suppressed.

## Figure: `fig_naive_ml_auroc.png` (+ `.pdf`, from `47b_naive_ml_figure.py`)

- **Panel a** — one point per task: the observed cross-validated AUROC, plotted against its
  permutation-null band (grey rectangle) with min–max whiskers. When `naive_ml_metrics.tsv` carries
  the percentile columns (min/5/25/50/75/95/max — the current schema), the band is the **exact**
  5–95th-percentile range of the actual shuffles and the whiskers are the true min/max; the
  figure's legend text says so. Against an older TSV missing those columns, `47b` instead
  reconstructs an approximate band from `perm_auroc_mean`/`perm_auroc_p95` alone (lower edge
  `2*mean - p95`, mirrored around the mean, floored at 0.5) and the legend flags it as an
  approximation with no whiskers drawn. Tasks are grouped and colour-coded (ancestry / platform /
  HLA->KIR cB); x-axis tick labels are coloured to match — no legend. Read it as: does the point
  sit clearly above its own grey band? Ancestry: yes, by a wide margin, for all six ancestries.
  Platform and cB_from_HLA: yes (in the pre-upgrade run), but the point is barely above the band
  and far below the ancestry tier — a real but small effect, not a strong one; the new
  `_EUR`/`_AFR`/`_adj_ancestry` tasks (not yet shown in this panel) are what settle whether that
  small effect is genuine or an ancestry echo (see "Is platform confounded by ancestry?").
- **Panel b** — the top 3 |L1 coefficient| HLA alleles per ancestry, restricted to the >=20-carrier
  feature set the model itself was fit on (no allele below that floor ever appears here, by
  construction — see Method). Bar direction shows which side of the one-vs-rest split the allele
  points toward, not "protective" or "risk" in any clinical sense.

## Honest read of the p-values (do not skip this if quoting a number from this script)

With only 20 permutations, the smallest possible empirical p-value is `1/21 ≈ 0.048` — every task
here, including the near-chance ones (platform, cB_from_HLA), sits at or near that floor. **A
p-value at the 20-shuffle floor says only "distinguishable from a null built from 20 shuffles,"
not "the effect is large."** The right way to read this run is the AUROC-vs-null-band comparison
in panel a, not the p-value column: ancestry's AUROC (0.81-0.98) is dramatically above its own
null band's upper edge, while platform's (0.573) and cB_from_HLA's (0.550) are only marginally
above theirs. If a future rerun affords more permutations (hundreds, not tens), the p-values would
start being informative on their own; until then, treat every `perm_pvalue_lr` in
`naive_ml_metrics.tsv` as "non-trivially above chance" or not, never as a precise significance
level.

## Is platform (AUROC 0.573) confounded by ancestry?

**Follow-up now implemented (pending a VM rerun for real numbers):** `47_naive_ml_tests.py` now
emits `platform_<X>_EUR` / `platform_<X>_AFR` (within-one-ancestry, ≥20/class/fold) and
`platform_<X>_adj_ancestry` (`naive_ml_ancestry_adjusted.tsv`, `delta_auroc` = carriage+ancestry
AUROC minus ancestry-only AUROC on identical folds), plus the same two views for `cB_from_HLA`
(`cB_from_HLA_EUR`, `cB_from_HLA_adj_ancestry`). Read it as: if the within-EUR/within-AFR AUROC
stays near the whole-cohort number and `delta_auroc` is clearly above 0, that is a genuine
platform-driven batch effect independent of ancestry; if the within-ancestry AUROC collapses
toward 0.5 and/or `delta_auroc` is near 0, platform's 0.573 (and cB_from_HLA's 0.550) were
substantially ancestry echoes rather than direct effects. This run's real full-cohort numbers for
these new tasks are not yet in this README — rerun on the VM (see "Running" above) and paste the
actual `naive_ml_ancestry_adjusted.tsv` numbers here before quoting a settled answer.

## Caveats

- **`--platform-col` needs VM verification.** `cohort_membership.tsv`'s `platform` column is
  documented in SCHEMA.md as sourced from `immuannot_cohort_full.tsv`, but the sibling script
  `42_repertoire_baseline.py` recorded a VM instance where an expected HLA column was entirely
  absent. This script exits loudly (printing available columns) rather than guessing if
  `--platform-col` doesn't match.
- Local dry runs use a small `--n-perms` for speed; the real run should use the new default (200,
  parallelized, lighter null model — see "Method" above) or more if VM time allows — see "Honest
  read of the p-values" below for why even 200 shuffles isn't a precision instrument.
- The cA/cB task uses a coarse person-level label (any B-content gene on either haplotype); it
  does not distinguish AA/AB/BB genotype dosage.
- Ancestry and platform predictability are expected and are not, on their own, evidence of a
  problem — they are the confound warnings this script exists to surface.
- The permutation-null band drawn in the figure is now the **exact** 5–95th-percentile range (with
  min–max whiskers) when `naive_ml_metrics.tsv` carries the new percentile columns; it falls back
  to the old mean/p95-based mirrored approximation only for an older TSV that predates this run —
  the figure's own legend text says which case applies.
- The ancestry-probability covariates (`p_afr..p_sas`) are soft-skipped, not fatal, if a VM
  instance's `cohort_membership.tsv` lacks them — the ancestry-adjusted tasks and
  `naive_ml_ancestry_adjusted.tsv` simply won't be produced in that case; check the run log.

## Plain language, for Marc

Three things were tested: can we guess someone's ancestry, their sequencing machine, or their KIR
gene content just from which HLA/KIR alleles they carry? Guessing ancestry works very well (as
expected — this is a warning sign for any future carriage-based disease analysis, not a finding).
Guessing the sequencing machine works only a little better than a coin flip (57%), and we can't
yet tell whether that little bit is a real machine effect or just ancestry leaking through (the
follow-up above would settle it). Guessing KIR content from HLA alone also works only a little
better than chance (55%), which is expected since these two genes are on different chromosomes
and don't have to be inherited together.

## Distilled

- Run on the full cohort (11,845 people). Ancestry is strongly predictable (AUROC 0.81-0.98,
  clearly above its permutation-null band) — a confound warning for future disease-carriage work,
  not a finding. Platform (AUROC 0.57) and cB_from_HLA (AUROC 0.55) show a modest but real lift
  above their own null bands, far below the ancestry tier; with only 20 permutations every task's
  p-value sits at the ~0.048 floor, so the AUROC-vs-null-band comparison (not the p-value) is the
  number to trust.
- Whether platform's 0.573 is itself an ancestry echo (rather than a genuine batch effect) is
  unresolved — no ancestry-adjusted or within-ancestry platform number exists yet; flagged as a
  follow-up.
- Figure: `fig_naive_ml_auroc.png`.
- Next step: figures (`45`-style rendering is out of scope for this pass) and a closer look at
  why platform shows any lift at all, before using HLA/KIR carriage in a downstream disease model
  without an ancestry (and possibly platform) covariate.
