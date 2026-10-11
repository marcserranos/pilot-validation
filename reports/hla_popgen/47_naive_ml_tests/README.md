# 47 — Naive ML / simple tests on allele-carriage matrices

*Status: v2, full cohort, 2026-09-27 (S04 WS-D VM session, script commit `8a5ea75`, results pulled
back in `c448a2b`). `--platform-col platform` verified against `cohort_membership.tsv` (values
`revio`/`sequel2e`/`sequel2`, counts 11,042/1,219/991). 14 tasks, 200 permutations/task, full
null-distribution columns, ancestry-adjusted delta-AUROC, and within-ancestry platform/cB subsets
are all present in this run — this is the settled numbers for "is platform confounded by
ancestry?" (see that section below). Figure re-rendered 2026-09-27 from this v2 TSV:
`fig_naive_ml_auroc.png` (+ `.pdf`), from `47b_naive_ml_figure.py` — panel a now shows 11 tasks
(ancestry x6, platform all/EUR-only/AFR-only, cB_from_HLA all/EUR-only); `platform_sequel2e*` rows
are omitted from the plot as numerically identical to `platform_revio*` in a two-valued subset
(same ROC curve, opposite label) — see the script's module-level comment — but remain in
`naive_ml_metrics.tsv` for anyone who wants that row directly.*

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
genes), 620 rare HLA alleles dropped before modeling, 14 tasks run, 200 permutations/task.

| task | n_total | n_pos | LR AUROC | tree AUROC | perm null p95 | perm p |
|---|---|---|---|---|---|---|
| ancestry_AFR | 11,845 | 3,009 | 0.950 | 0.782 | 0.513 | 0.005 |
| ancestry_AMR | 11,845 | 2,656 | 0.810 | 0.652 | 0.515 | 0.005 |
| ancestry_EAS | 11,845 | 1,460 | 0.981 | 0.813 | 0.517 | 0.005 |
| ancestry_EUR | 11,845 | 2,974 | 0.897 | 0.712 | 0.512 | 0.005 |
| ancestry_MID | 11,845 | 487 | 0.897 | 0.628 | 0.527 | 0.005 |
| ancestry_SAS | 11,845 | 1,236 | 0.966 | 0.671 | 0.519 | 0.005 |
| platform_revio | 11,845 | 10,710 | 0.574 | 0.507 | 0.522 | 0.005 |
| platform_sequel2e | 11,845 | 1,135 | 0.573 | 0.507 | 0.522 | 0.005 |
| platform_revio_EUR | 2,974 | 2,753 | 0.558 | 0.514 | 0.547 | 0.035 |
| platform_revio_AFR | 3,009 | 2,770 | 0.489 | 0.507 | 0.541 | 0.627 |
| platform_sequel2e_EUR | 2,974 | 221 | 0.558 | 0.514 | 0.547 | 0.035 |
| platform_sequel2e_AFR | 3,009 | 239 | 0.489 | 0.507 | 0.541 | 0.627 |
| cB_from_HLA | 11,845 | 7,978 | 0.550 | 0.526 | 0.512 | 0.005 |
| cB_from_HLA_EUR | 2,974 | 2,087 | 0.521 | 0.509 | 0.525 | 0.100 |

With 200 permutations, all six ancestry tasks and the whole-cohort `platform_*`/`cB_from_HLA`
tasks sit at the 200-shuffle empirical-p floor (`1/201 ≈ 0.005`); the within-ancestry subsets do
not, and that's the informative part of this table. A real effect size — not just
"distinguishable from a 200-shuffle null" — is what separates ancestry (AUROC 0.81-0.98, far above
its own null p95 of 0.51-0.53) from platform/cB_from_HLA (AUROC 0.55-0.57 whole-cohort, well above
0.5 but far below the ancestry tier). The within-ancestry breakdown resolves the open question this
README carried in the v1 run: platform's whole-cohort 0.574 **drops to 0.558 within EUR (p=0.035,
just above its own null p95 of 0.547) and to 0.489 within AFR (p=0.627, indistinguishable from
chance)** — see "Is platform confounded by ancestry?" below for the full reading. Top
AFR-ancestry-informative alleles (L1 coefficients): HLA-B\*07:05, HLA-C\*16:02, HLA-B\*38:02
(negative), HLA-DPB1\*18:01, HLA-A\*80:01 (positive) — see `naive_ml_top_features.tsv` for every
task's full list.

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
  permutation-null band (grey rectangle) with min–max whiskers. This v2 TSV carries the full
  percentile columns (min/5/25/50/75/95/max), so the band shown is the **exact** 5–95th-percentile
  range of the actual 200 shuffles and the whiskers are the true min/max — the figure's own legend
  text confirms this (`47b`'s `null_band()` falls back to an approximate mirrored band only for an
  older, pre-percentile-column TSV, which this is not). 11 of the 14 tasks are plotted (ancestry x6,
  `platform_revio` all/EUR-only/AFR-only, `cB_from_HLA` all/EUR-only); `platform_sequel2e*` are
  omitted as numerically identical to `platform_revio*` (same ROC curve read from the other
  direction in a two-valued subset) but remain in the TSV. Tasks are grouped and colour-coded
  (ancestry / platform / HLA->KIR cB); x-axis tick labels are coloured to match — no legend. Read it
  as: does the point sit clearly above its own grey band? Ancestry: yes, by a wide margin, for all
  six ancestries. Platform (all) and cB_from_HLA (all): yes, but only just, and far below the
  ancestry tier. Platform (EUR only) sits right at the edge of its own band (p=0.035) and platform
  (AFR only) sits inside it, below 0.5 (p=0.627) — visually confirming the ancestry-echo reading
  in "Is platform confounded by ancestry?" below.
- **Panel b** — the top 3 |L1 coefficient| HLA alleles per ancestry, restricted to the >=20-carrier
  feature set the model itself was fit on (no allele below that floor ever appears here, by
  construction — see Method). Bar direction shows which side of the one-vs-rest split the allele
  points toward, not "protective" or "risk" in any clinical sense.

## Honest read of the p-values (do not skip this if quoting a number from this script)

With 200 permutations (this v2 run), the smallest possible empirical p-value is `1/201 ≈ 0.005` —
every whole-cohort task here (all six ancestries, platform all, cB_from_HLA all) sits at or near
that floor, whether or not the effect is large. **A p-value at the 200-shuffle floor says only
"distinguishable from a null built from 200 shuffles," not "the effect is large."** The right way
to read the whole-cohort tasks is the AUROC-vs-null-band comparison in panel a, not the p-value
column: ancestry's AUROC (0.81-0.98) is dramatically above its own null band's upper edge, while
platform's (0.574) and cB_from_HLA's (0.550) are only marginally above theirs. The within-ancestry
subsets (`platform_revio_EUR` p=0.035, `platform_revio_AFR` p=0.627, `cB_from_HLA_EUR` p=0.100) are
past the floor and are informative on their own — and they're the numbers that settle the
confounding question (see "Is platform confounded by ancestry?" above). 200 shuffles is still not a
precision instrument for the whole-cohort p-values; treat every floored `perm_pvalue_lr` in
`naive_ml_metrics.tsv` as "non-trivially above chance," never as a precise significance level.

## Is platform (AUROC 0.574) confounded by ancestry?

**Settled by this v2 run.** `naive_ml_ancestry_adjusted.tsv`:

| task | n_total | n_pos | AUROC (carriage+ancestry) | AUROC (ancestry only) | delta_auroc |
|---|---|---|---|---|---|
| platform_revio_adj_ancestry | 11,822 | 10,687 | 0.6363 | 0.6324 | +0.0039 |
| platform_sequel2e_adj_ancestry | 11,822 | 1,135 | 0.6363 | 0.6324 | +0.0039 |
| cB_from_HLA_adj_ancestry | 11,822 | 7,965 | 0.5662 | 0.5887 | -0.0225 |

Both lines of evidence agree, and both point the same way: **platform's whole-cohort signal is
mostly an ancestry echo, not a genuine sequencing-batch effect, and HLA carries essentially no
independent information about KIR cA/cB content beyond what ancestry already explains.**

- **Within-ancestry AUROC collapses.** Platform drops from 0.574 (whole cohort) to 0.558 within
  EUR (barely above its own null p95 of 0.547, p=0.035) and to 0.489 within AFR (below 0.5, p=0.627
  — indistinguishable from chance). If platform carriage signal were a real, ancestry-independent
  batch effect, it should have stayed near 0.574 in both subgroups; instead it nearly vanishes.
- **Ancestry-adjusted delta is ~0 for platform.** Adding carriage features on top of the 6 ancestry
  probabilities buys essentially nothing: 0.6363 vs. 0.6324 ancestry-alone, delta +0.0039. (Note
  ancestry-only AUROC for predicting platform is itself 0.63, not 0.5 — platform composition
  differs by ancestry group, e.g. different recruitment sites/kit versions per ancestry — so most of
  platform's whole-cohort predictability was already ancestry-explainable before carriage features
  are added at all.)
- **cB_from_HLA's delta is *negative*.** 0.5662 (carriage+ancestry) vs. 0.5887 (ancestry only),
  delta -0.0225 — adding HLA carriage features on top of ancestry makes the cB prediction slightly
  *worse* than ancestry alone (within noise of a single-fold-seed comparison, but never positive).
  This is the strongest form of "no independent signal": whole-cohort cB_from_HLA's 0.550 AUROC is
  fully accounted for by ancestry, consistent with chr6 (HLA) and chr19 (KIR) having no known
  direct linkage — the modest whole-cohort lift was population-structure confounding, not
  HLA-genotype-predicts-KIR-content.

Read together with "Honest read of the p-values" below: none of platform's or cB_from_HLA's
whole-cohort p-values (both at the 200-shuffle floor, `1/201 ≈ 0.005`) were ever evidence of a
large effect — they only said "distinguishable from a 200-shuffle null." The within-ancestry and
ancestry-adjusted numbers above are what actually answer the confounding question, and they answer
it: no, don't use raw platform or cB_from_HLA carriage-predictability as evidence of a real
biological or batch effect independent of ancestry, in this cohort, with this feature set.

## Caveats

- **`--platform-col` needs VM verification.** `cohort_membership.tsv`'s `platform` column is
  documented in SCHEMA.md as sourced from `immuannot_cohort_full.tsv`, but the sibling script
  `42_repertoire_baseline.py` recorded a VM instance where an expected HLA column was entirely
  absent. This script exits loudly (printing available columns) rather than guessing if
  `--platform-col` doesn't match.
- This run used the default `--n-perms 200` (parallelized, lighter null model — see "Method"
  above); local dry runs use a smaller `--n-perms` for speed. See "Honest read of the p-values"
  above for why even 200 shuffles isn't a precision instrument for the whole-cohort p-values.
- The cA/cB task uses a coarse person-level label (any B-content gene on either haplotype); it
  does not distinguish AA/AB/BB genotype dosage.
- Ancestry and platform predictability are expected and are not, on their own, evidence of a
  problem — they are the confound warnings this script exists to surface.
- The permutation-null band drawn in the figure is the **exact** 5–95th-percentile range (with
  min–max whiskers) for this v2 `naive_ml_metrics.tsv`, which carries the percentile columns; `47b`
  falls back to an old mean/p95-based mirrored approximation only for an older TSV that predates
  this run — the figure's own legend text says which case applies.
- The ancestry-probability covariates (`p_afr..p_sas`) are soft-skipped, not fatal, if a VM
  instance's `cohort_membership.tsv` lacks them — the ancestry-adjusted tasks and
  `naive_ml_ancestry_adjusted.tsv` simply won't be produced in that case; check the run log.

## Plain language, for Marc

Three things were tested: can we guess someone's ancestry, their sequencing machine, or their KIR
gene content just from which HLA/KIR alleles they carry? Guessing ancestry works very well (as
expected — this is a warning sign for any future carriage-based disease analysis, not a finding).
Guessing the sequencing machine looked like it worked a little better than a coin flip (57% overall)
— but when we only look within one ancestry group at a time, that edge mostly disappears: 56%
within European-ancestry people (barely above noise) and 49% within African-ancestry people (worse
than a coin flip). And when we let a second model see both HLA carriage and ancestry side by side,
knowing ancestry already gets 63% and adding HLA carriage buys essentially nothing more. So the
"guess the machine" signal wasn't really about the machine — it was mostly picking up that
different ancestry groups happen to have been run on different machines more or less often, not
that HLA carriage itself carries a machine fingerprint. Guessing KIR content (the cA/cB group) from
HLA alone also only worked a little better than chance (55%), and the same check shows that number,
too, disappears once ancestry is accounted for — HLA and KIR sit on different chromosomes and
aren't expected to travel together, and this run confirms that expectation rather than overturning
it. Bottom line: nothing here suggests HLA or KIR carriage directly encodes sequencing-platform
identity or predicts each other's content — the only strong, real signal in this whole exercise is
ancestry, which is exactly the confound every future carriage-based analysis on this cohort will
need to control for.

## Distilled

- Run on the full cohort (11,845 people, 14 tasks, 200 permutations/task — v2, settled numbers).
  Ancestry is strongly predictable (AUROC 0.81-0.98, clearly above its permutation-null band) — a
  confound warning for future disease-carriage work, not a finding.
- Platform (AUROC 0.57 whole-cohort) and cB_from_HLA (AUROC 0.55 whole-cohort) showed a modest lift
  above their own null bands, but this is **resolved as mostly/fully an ancestry echo, not a
  genuine effect**: platform drops to 0.558 within EUR (p=0.035) and 0.489 within AFR (p=0.627,
  indistinguishable from chance), and its ancestry-adjusted delta-AUROC is only +0.004; cB_from_HLA
  within EUR is 0.521 (p=0.100) and its ancestry-adjusted delta is *negative* (-0.023) — i.e. HLA
  carries no information about KIR cA/cB beyond ancestry, consistent with HLA (chr6) and KIR (chr19)
  having no known direct linkage.
- With 200 permutations, whole-cohort p-values still floor at ~0.005 regardless of effect size — the
  AUROC-vs-null-band comparison (and, for platform/cB_from_HLA specifically, the within-ancestry and
  ancestry-adjusted numbers above) are the numbers to trust, not the whole-cohort p-value column.
- Figure: `fig_naive_ml_auroc.png` (v2, re-rendered 2026-09-27; 11 tasks shown, exact 5-95%
  permutation-null band).
- Next step: use HLA/KIR carriage in any downstream disease model only with an ancestry covariate;
  platform does not need its own separate covariate beyond ancestry based on this result, though
  that conclusion rests on `revio`/`sequel2e` group sizes within EUR/AFR (2,753/221 and 2,770/239)
  and hasn't been checked in the other four ancestry groups.
