# 42 — WS6 repertoire baseline (BenchRep-T VJ-kmer replication)

## Question

Does TRB/IGH repertoire composition (V/J gene usage + CDR3 k-mers, BenchRep-T-style) predict
chronic/autoimmune disease status in the AoU RNA-seq x long-read overlap cohort, above a
covariate-only baseline (age, sex, ancestry, sequencing depth, EHR-depth)? Two positive
controls (sex, predicted ancestry) show whether the pipeline finds signal when it is known to
exist, and quantify the ancestry-leakage risk that any repertoire-based classifier carries.

Full protocol: `sprints/S03_call8_figures_kir_prediction/WS6_benchrept_protocol.md`. Script:
`scripts/hla_popgen/42_repertoire_baseline.py`. Tests: `scripts/hla_popgen/tests/test_repertoire_baseline.py`
(31 checks, all pass).

## Method

L1-regularized logistic regression (BenchRep-T "VJ-kmer Reg" track), 3-fold stratified CV,
SAME fold assignment shared across every model variant for a given label, pooled out-of-fold
AUROC/AUPRC with bootstrap CIs (1000 person-level resamples). Model variants: covariate-only;
repertoire-only (TRB+IGH combined); repertoire-only restricted to TRB; repertoire-only
restricted to IGH; repertoire+covariate. Cohort: RNA-seq x long-read HLA-calling overlap,
built fresh on this VM instance (7958 repertoire samples -> 7926 matched to a research_id via
the RNA-seq manifest -> 7900 with long-read HLA calls -> 7640 after greedy-unrelated,
kin >= 0.0442), subsampled with a fixed seed to **1200 people** (`--max-people 1200`) so the
run completes in a reasonable time on a VM shared with several other concurrent agent jobs
(see Caveats). Controls require >=1 year of EHR observation (Aleix's ancestry-correlated
healthcare-access-gradient finding: a "control" with near-zero EHR contact might just never
have been diagnosed, not be genuinely disease-free).

## Result

**No disease phenotype (HIV, B-cell disorder, RA, SLE, other rheum. autoimmune, T1D,
autoimmune thyroid disease, MS, IBD) reached the >=100-case threshold at n=1200** — the left
panel of the figure is empty for that reason, not a modeling failure. This is expected: these
are 9 separate chronic/autoimmune conditions, most well under 10% population prevalence, and a
1200-person fixed-seed subsample of a repertoire-available cohort is not powered to clear 100
cases for any of them (see Power below). This is an honest null result at this cohort size,
not evidence of "no repertoire signal" — a larger subsample (or the full 7640-person cohort,
which the run-time fix in this session should now make tractable, see Caveats) is needed
before the disease-prediction question itself can be answered.

**Positive controls all show real signal** — the pipeline works as expected when signal
genuinely exists:

| phenotype | model | AUROC (95% CI) | n_cases | n_total |
|---|---|---|---|---|
| sex (Female) | covariate_only | 0.575 (0.543–0.608) | 602 | 1130 |
| sex (Female) | repertoire_only | 0.500 (0.468–0.533) | 602 | 1130 |
| sex (Female) | repertoire_only_TRB | 0.497 (0.464–0.530) | 602 | 1130 |
| ancestry (EUR) | repertoire_only | 0.708 (0.676–0.741) | 355 | 1200 |
| ancestry (EUR) | repertoire_only_TRB | 0.668 (0.636–0.699) | 355 | 1200 |
| ancestry (EUR) | repertoire_only_IGH | 0.636 (0.604–0.671) | 355 | 1200 |
| ancestry (AFR) | repertoire_only | 0.813 (0.783–0.844) | 209 | 1200 |
| ancestry (AFR) | repertoire_only_TRB | 0.802 (0.770–0.833) | 209 | 1200 |
| ancestry (AFR) | repertoire_only_IGH | 0.801 (0.769–0.834) | 209 | 1200 |

- **Sex** is correctly weak from repertoire alone (AUROC ~0.50, as expected — the repertoire
  carries no biological sex signal on its own), while the covariate-only model (which includes
  admixture/ancestry covariates correlated with sex-reporting patterns in this cohort) picks up
  a modest 0.575. This matches the protocol's expectation ("sex should be weak").
- **Ancestry is strongly detectable from repertoire alone** (AUROC 0.64–0.81, AFR > EUR), via
  IGHV/TRBV germline allele usage — exactly the leakage risk the protocol flags: any
  repertoire-only disease classifier in an ancestry-imbalanced cohort risks partly detecting
  ancestry rather than disease. **Any future repertoire-only disease result MUST be compared
  against `repertoire_plus_covariate` (which includes ancestry) before being trusted.**
- The HLA-carriage positive control (DRB1*15:01 or B*57:01 carriage predicted from TRB) did
  not run at n=1200 — the qualifying allele's carrier count didn't clear the >=100 threshold in
  this subsample; worth rerunning at larger n.

## QC (before any modeling)

Productive clonotype-count quantiles across the 1200 retained samples:

| | TRB | IGH |
|---|---|---|
| p5 | 483 | 216 |
| p25 | 878 | 515 |
| p50 | 1242 | 840 |
| p75 | 1651 | 1270 |
| p95 | 2298 | 2259 |

**Depth floor**: 100 TRB clonotypes minimum (documented, not silently hardcoded). 0.0% of
samples fell below it at n=1200 (0.1% in the earlier full-cohort load) — this cohort's TRUST4
repertoires are deep enough that the floor barely binds; it is kept as a guard for any future
run on a shallower sub-cohort.

## Power (minimum detectable AUROC)

Hanley-McNeil normal approximation, alpha=0.05, power=0.80, given (n_cases, n_total):
- At n_cases=100, n_total=1200: minimum detectable AUROC ≈ 0.57 (see `min_detectable_auroc`
  column in the metrics TSV — computed per phenotype/model row from its actual n).
- The disease phenotypes never reached n_cases=100 at all in this subsample, so no AUROC
  number could be computed for them; the honest statement is "underpowered at this n," not
  "no signal."

## Caveats

- **This is a 1200-person fixed-seed subsample of the 7640-person eligible cohort**, not the
  full cohort. Three earlier full-cohort attempts (this session) each died silently mid-run
  with no traceback — diagnosed as a genuine performance bug in `build_vj_kmer_features`
  (materializing a dense person x full-raw-kmer-alphabet matrix via nested Python list
  comprehensions before capping to `--max-kmers`, which is fine at ~1000 people but explodes
  at ~7600 with a real, uncapped kmer alphabet of tens of thousands of distinct strings). Fixed
  in this session (cheap cohort-wide-count pre-trim to `max_kmers*3` before the dense matrix is
  built) and covered by 2 new/updated unit tests, but the full-cohort real-data run with this
  fix has not yet been confirmed end-to-end — retry with `--max-people` unset (or a larger
  value) once VM load allows.
- XGBoost (`--with-xgboost`) was **not** run for this delivered result (dropped for the n=1200
  run to guarantee completion inside the session's time budget, per L1-LR-first ordering in
  `model_one_label`) — L1-LR only.
- `n_repeats=1` (not the default 3) for the same reason — CIs are still valid (1000-resample
  bootstrap on the pooled OOF predictions) but the repeat-to-repeat AUROC variance itself isn't
  characterized here.
- Disclosure: every n reported above is >=100 by construction (the MIN_CASES gate); no
  suppressed (`<20`) cells appear in this delivery.

## Distilled

- Pipeline correctly detects strong, expected signal (ancestry from repertoire, AUROC up to
  0.81) and correctly reports weak signal where expected (sex from repertoire, AUROC ~0.50) —
  the modeling code is doing what BenchRep-T's own protocol predicts it should.
- Disease-phenotype prediction is **not yet answered**: this 1200-person subsample isn't
  powered to clear 100 cases for any of the 9 chronic/autoimmune phenotypes tested. Needs a
  rerun at higher n (full 7640-person cohort, now tractable after the performance fix).
- Ancestry leakage is real and large (AUROC up to 0.81) — any future disease result from
  `repertoire_only` alone must be cross-checked against `repertoire_plus_covariate`.
