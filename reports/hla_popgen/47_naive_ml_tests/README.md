# 47 — Naive ML / simple tests on allele-carriage matrices

*Status: run on the full cohort 2026-09-25/26 (S04 WS-D VM session). `--platform-col platform`
verified against `cohort_membership.tsv` (values `revio`/`sequel2e`/`sequel2`, counts 11,042/
1,219/991). See "Result" below for the real numbers. Figure added 2026-09-26:
`fig_naive_ml_auroc.png` (+ `.pdf`), from `47b_naive_ml_figure.py`.*

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

- `naive_ml_metrics.tsv` — one row per task: `lr_auroc`/`tree_auroc` (0.5 = chance, 1.0 =
  perfect), `perm_auroc_mean`/`perm_auroc_p95` (the null distribution — a real AUROC should sit
  above `perm_auroc_p95`), `perm_pvalue_lr`.
- `naive_ml_top_features.tsv` — which alleles/genes drove each model, per task.
- `pca_group_centroids.tsv` — mean PC1/PC2 per ancestry x platform group (n>=20 only); platform
  clustering here would flag a batch effect, not biology.
- `pca_density_grid.tsv` — 2D histogram of PC1/PC2 over all people, binned and count-suppressed.

## Figure: `fig_naive_ml_auroc.png` (+ `.pdf`, from `47b_naive_ml_figure.py`)

- **Panel a** — one point per task: the observed cross-validated AUROC, plotted against its
  permutation-null band (grey rectangle). The band is a *reconstruction*, not the raw 20-shuffle
  range: `naive_ml_metrics.tsv` stores only `perm_auroc_mean`/`perm_auroc_p95`, not the individual
  shuffle values, so the band's lower edge is `2*mean - p95` (mirrored around the mean, floored at
  0.5) rather than the true observed minimum. Tasks are grouped and colour-coded (ancestry /
  platform / HLA->KIR cB); x-axis tick labels are coloured to match — no legend. Read it as: does
  the point sit clearly above its own grey band? Ancestry: yes, by a wide margin, for all six
  ancestries. Platform and cB_from_HLA: yes, but the point is barely above the band and far below
  the ancestry tier — a real but small effect, not a strong one.
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

Checked: neither `naive_ml_metrics.tsv` nor any other TSV in this run's output contains an
ancestry-adjusted platform task or a within-one-ancestry (e.g. within-EUR-only) platform AUROC —
the platform tasks (`platform_revio`, `platform_sequel2e`) were fit on the whole cohort, ancestry
mixed in. Since ancestry itself is so strongly encoded in the same carriage features (panel a),
and platform assignment likely correlates with ancestry through recruitment-site/cohort/batch
effects rather than genetics, **platform's 0.573 could be substantially an ancestry echo rather
than a genuine platform-driven batch effect on the calls** — this run cannot distinguish the two
explanations. **Follow-up (not done in this pass):** rerun the platform task stratified within
each ancestry (or add ancestry as a covariate/control feature) and see whether the AUROC drops
toward 0.5; if it does, the confound explanation wins and the calls likely don't carry a real
platform batch effect independent of ancestry.

## Caveats

- **`--platform-col` needs VM verification.** `cohort_membership.tsv`'s `platform` column is
  documented in SCHEMA.md as sourced from `immuannot_cohort_full.tsv`, but the sibling script
  `42_repertoire_baseline.py` recorded a VM instance where an expected HLA column was entirely
  absent. This script exits loudly (printing available columns) rather than guessing if
  `--platform-col` doesn't match.
- Local dry run used `--n-perms 10-20` for speed; the real run should raise this (>=20 minimum
  per the WS-D brief, more if VM time allows) — see "Honest read of the p-values" above for why
  20 is already a hard floor on what the p-values themselves can say.
- The cA/cB task uses a coarse person-level label (any B-content gene on either haplotype); it
  does not distinguish AA/AB/BB genotype dosage.
- Ancestry and platform predictability are expected and are not, on their own, evidence of a
  problem — they are the confound warnings this script exists to surface.
- The permutation-null *band* drawn in the figure is a mean/p95-based approximation of the true
  20-shuffle spread, not the spread itself (see "Figure" above) — do not read the band edges as
  exact percentiles.

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
