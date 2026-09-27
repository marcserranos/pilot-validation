# Critic #47 — fresh-context review of 47_naive_ml_tests v2 (S04, 2026-09-27/28)

Reviewed: `reports/hla_popgen/47_naive_ml_tests/` (README.md, all 5 TSVs, `fig_naive_ml_auroc.png`
read at native resolution) and `scripts/hla_popgen/47_naive_ml_tests.py` / `47b_naive_ml_figure.py`,
against `sprints/S04_kir_recurrence_style_share/AGENT_PREAMBLE.md` and
`sprints/S03_call8_figures_kir_prediction/CRITIC_2.md` (format). Scope: script 47/47b only —
44/45/46 untouched, per instructions (another critic owns those).

## Was the README/figure updated to v2?

**No, not before this pass — this was the main finding.** `naive_ml_metrics.tsv` and
`naive_ml_ancestry_adjusted.tsv` were pulled back as v2 in commit `c448a2b` (200 perms, 14 tasks,
within-EUR/AFR platform and cB subsets, ancestry-adjusted delta-AUROC), but `README.md` and both
figure files still dated from the pre-upgrade 20-permutation/9-task run, and still contained
explicit "not yet re-run on the VM" / "pending a VM rerun" placeholder language in the "Is platform
confounded by ancestry?" section — i.e. the one question this v2 run exists to answer was still
marked unanswered in the shipped doc, with the answer sitting unused in a committed TSV one
directory over.

**Fixed:**
- Rewrote the status header, Result table (9→14 tasks, 20→200 perms, added the 5 new
  platform/cB rows and the null-p95 column), and the "Is platform confounded by ancestry?" section
  with the actual `naive_ml_ancestry_adjusted.tsv` numbers and a full delta_auroc reading.
- Updated "Honest read of the p-values" (1/21≈0.048 floor → 1/201≈0.005 floor; added that the
  within-ancestry subsets are past the floor and are the informative numbers).
- Rewrote "Plain language, for Marc" and "Distilled" to state the settled, honest interpretation
  (see below) instead of leaving it as an open follow-up.
- `47b_naive_ml_figure.py`: `TASK_ORDER`/`TASK_LABEL`/`TASK_GROUP` only listed the original 9
  tasks, so the figure — even after a naive re-render — would have silently omitted the new
  `_EUR`/`_AFR` tasks the whole story now depends on. Added `platform_revio_EUR/AFR` and
  `cB_from_HLA_EUR` to the plotted set (`platform_sequel2e*` deliberately excluded and documented
  as numerically identical to `platform_revio*` in a two-valued one-vs-rest subset — verified this
  is true from the TSV: both pairs share identical AUROC/p95/p-value to 4 decimals). Re-rendered;
  first attempt failed `_viz_common.check_layout()`'s strict mode with 17 violations (text overlap
  from a naive addition of all sequel2e rows too, plus point labels landing on whisker lines) —
  fixed by dropping the redundant sequel2e rows, shortening the cB_from_HLA labels, and moving the
  per-point AUROC labels from directly-above to beside the marker. Final render passes
  `check_layout()` clean; inspected the PNG at native resolution — no residual overlap.

The null band itself was already using the exact 5-95th-percentile columns (`null_band()`'s
preferred path was already correct and the TSV already had the percentile columns) — the "may
still use the mirrored null band" concern in the brief did not materialize as an actual bug; the
real gap was the stale task list, not stale band logic.

## Number traceability

Re-derived every number I put in the README from the TSVs directly (not copied from the brief):
`naive_ml_metrics.tsv` for all 14 rows' AUROC/p95/p-value, `naive_ml_ancestry_adjusted.tsv` for the
3 delta rows, `naive_ml_top_features.tsv` for the AFR top-feature list (unchanged from v1 —
spot-checked the full ancestry_AFR/l1_logreg block, confirmed HLA-B\*07:05/-C\*16:02/-B\*38:02
negative and HLA-DPB1\*18:01/-A\*80:01 positive are exactly the top-2-each-direction by
|coefficient|). All figures matched to 3 decimal places; the brief's supplied numbers (platform EUR
0.558/null p95 0.547, AFR 0.489, ancestry-adjusted +0.004/-0.023, 200 perms) all confirmed exactly
against the committed TSVs.

## Interpretation

Confirmed sound and now stated explicitly rather than left open: platform's whole-cohort AUROC
(0.574) collapses within-ancestry (0.558 EUR barely above its own null p95; 0.489 AFR, below
chance) and its ancestry-adjusted delta is only +0.004 — i.e. mostly an ancestry echo (platform mix
differs by ancestry group), not a genuine sequencing-batch effect. cB_from_HLA's ancestry-adjusted
delta is actually *negative* (-0.023), which is the cleanest statement available that HLA carries
no independent information about KIR cA/cB beyond ancestry — consistent with HLA (chr6) and KIR
(chr19) having no known direct linkage. This was already the right hedge in the v1 prose ("plausibly
an ancestry-mediated correlation... not proof") — v2 turns "plausibly" into a settled, quantified
answer, which the README previously failed to state because it was never rewritten.

## Statistical soundness

- **Nested/leakage check:** `run_binary_task` and `run_delta_task` both build the full feature
  matrix once, then do `StratifiedKFold` + `cross_val_predict` for out-of-fold predictions — no
  fitting on the full matrix before CV split, no leakage found. `run_delta_task` uses the *same*
  `random_state=seed` for both the carriage+ancestry and ancestry-only fits, so the delta_auroc
  comparison is apples-to-apples on identical folds — correct design for isolating the marginal
  contribution of carriage features.
- **Lighter null model:** documented in both the script's module docstring and the README's Method
  section (3-fold vs 5-fold CV, capped `max_iter`, justified as not affecting the null's location
  near 0.5). This is a legitimate speed optimization for a permutation null and is disclosed, not
  hidden.
- **Per-fold power gate:** `platform_<X>_EUR/_AFR` and `cB_from_HLA_EUR` require ≥20 people per
  class *per CV fold* (stricter than the plain ≥20-total rule elsewhere) — appropriate given 5-fold
  CV needs each fold to have enough of the minority class.
- Ran `pytest scripts/hla_popgen/tests -k 47` (both 47 and 47b test files): 29 passed, 0 failed.

## Disclosure

- No count 1-19 printed anywhere in README prose or TSVs; `pca_density_grid.tsv` has 126 cells at
  the literal string `<20` (grep-verified), true 0 cells left as `0`.
- No person IDs found in any committed file (grepped for digit runs ≥6 and found none in the 47
  outputs).
- Every feature in `naive_ml_top_features.tsv` and every n_pos/n_total quoted in the README is
  ≥20 by construction (`MIN_CARRIERS_FOR_FEATURE` gate before modeling; smallest n_pos quoted in
  the README is 221, `platform_sequel2e_EUR`) — no allele name paired with a small carrier count.

## Plain-language paragraph for Marc (now in the README, reproduced here)

Three things were tested: can we guess someone's ancestry, sequencing machine, or KIR gene content
just from HLA/KIR carriage? Ancestry: yes, strongly — expected, and a warning for future
carriage-based disease work, not a finding. Sequencing machine looked like a weak but real signal
overall (57%), but splitting by ancestry group shows that signal mostly evaporates — 56% within
European-ancestry people (barely above noise) and 49% within African-ancestry people (worse than a
coin flip) — and a side-by-side model shows ancestry alone already explains almost all of it. So the
"guess the machine" result was mostly "different ancestry groups happened to run on different
machines more often," not HLA/KIR carriage fingerprinting the machine. KIR content guessed from HLA
alone (55%) shows the same pattern even more clearly — accounting for ancestry actually makes the
guess very slightly *worse*, which is about as clean a "no" as this kind of test can give, and
matches the biological expectation that HLA and KIR (different chromosomes) don't have to travel
together.

## Left open / not fixed

- Nothing structurally left open in scope. One judgment call flagged but not changed: the
  `platform_sequel2e_EUR/_AFR` conclusion (echoed from `platform_revio_EUR/_AFR`) rests on small
  sequel2e-within-ancestry group sizes (n_pos 221 EUR, 239 AFR) and wasn't checked outside
  EUR/AFR — noted in the README's Distilled section as a residual caveat, not resolved here (would
  need a VM rerun with more ancestry groups included in the within-ancestry check, a scope decision
  for the human, not a text fix).

Commit: see below (`git log -1`).
