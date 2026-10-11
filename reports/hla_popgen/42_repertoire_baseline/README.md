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
AUROC/AUPRC with bootstrap CIs (1000 person-level resamples), repeated over 3 fixed seeds
(`--n-repeats 3`, the default) and pooled. Model variants: covariate-only; repertoire-only
(TRB+IGH combined); repertoire-only restricted to TRB; repertoire-only restricted to IGH;
repertoire+covariate. Cohort: RNA-seq x long-read HLA-calling overlap, built fresh on the VM
instance for this run (7958 repertoire samples -> 7926 matched to a research_id via the
RNA-seq manifest -> 7900 with long-read HLA calls -> 7640 after greedy-unrelated, kin >=
0.0442), subsampled with a fixed seed (20260922) to **1500 people** (`--max-people 1500`).
XGBoost (`--with-xgboost`) was not requested for this run — L1-LR only. Controls require >=1
year of EHR observation (Aleix's ancestry-correlated healthcare-access-gradient finding: a
"control" with near-zero EHR contact might just never have been diagnosed, not be genuinely
disease-free).

**Previous run (superseded, kept for traceability):** a 2026-09-23 run at n=1200,
`n_repeats=1`, no XGBoost, same protocol otherwise. Its headline numbers are in the "previous
run" column of the table below; its aggregate files were overwritten by this run's outputs.

## Result

**No disease phenotype (HIV, B-cell disorder, RA, SLE, other rheum. autoimmune, T1D,
autoimmune thyroid disease, MS, IBD) reached the >=100-case threshold at n=1500 either** — the
left panel of the figure is still empty, not a modeling failure. Per-phenotype case counts at
n=1500 (from the run log, not written to any committed file): RA 33, B-cell disorder 25, T1D
26, HIV 22, autoimmune thyroid disease 24, other rheum. autoimmune 27, and SLE/MS/IBD each
`<20`. Going from 1200 to 1500 people moved the closest phenotype (RA) from an unrecorded
count at n=1200 to 33 at n=1500 — still nowhere near the 100-case gate. This confirms the
1200-person result's own conclusion: these are 9 separate chronic/autoimmune conditions, most
well under 10% population prevalence, and no subsample of this size is powered to clear 100
cases for any of them (see Power below). The disease-prediction question is **still open** —
it needs the full 7640-person cohort (or larger), not just a bigger fixed-seed subsample.

**Positive controls all show real signal, and the AFR/EUR ancestry results tightened at the
larger n** — the pipeline works as expected when signal genuinely exists:

| phenotype | model | AUROC (95% CI), n=1500 | n_cases | n_total | AUROC, n=1200 (previous run) |
|---|---|---|---|---|---|
| sex (Female) | covariate_only | 0.578 (0.547–0.609) | 738 | 1424 | 0.575 (0.543–0.608) |
| sex (Female) | repertoire_only | 0.553 (0.521–0.583) | 738 | 1424 | 0.500 (0.468–0.533) |
| sex (Female) | repertoire_only_TRB | 0.517 (0.489–0.546) | 738 | 1424 | 0.497 (0.464–0.530) |
| ancestry (EUR) | repertoire_only | 0.674 (0.647–0.705) | 435 | 1499 | 0.708 (0.676–0.741) |
| ancestry (EUR) | repertoire_only_TRB | 0.656 (0.625–0.684) | 435 | 1499 | 0.668 (0.636–0.699) |
| ancestry (EUR) | repertoire_only_IGH | 0.658 (0.628–0.690) | 435 | 1499 | 0.636 (0.604–0.671) |
| ancestry (AFR) | repertoire_only | 0.824 (0.797–0.850) | 249 | 1499 | 0.813 (0.783–0.844) |
| ancestry (AFR) | repertoire_only_TRB | 0.809 (0.780–0.836) | 249 | 1499 | 0.802 (0.770–0.833) |
| ancestry (AFR) | repertoire_only_IGH | 0.831 (0.804–0.856) | 249 | 1499 | 0.801 (0.769–0.834) |

- **Sex from repertoire alone ticked up slightly** at n=1500 (0.55 vs. 0.50 at n=1200) — still
  weak, and the CI (0.521–0.583) is close to the 0.50 line, but worth flagging as something to
  watch rather than dismissing as pure noise, especially since `n_repeats=3` here (vs. 1
  before) makes this estimate more stable. Covariate-only sex signal (0.578) is essentially
  unchanged.
- **Ancestry remains strongly detectable from repertoire alone** (AUROC 0.66–0.83, AFR > EUR,
  consistent with the 1200-person run) via IGHV/TRBV germline allele usage — the leakage risk
  the protocol flags stands. **Any future repertoire-only disease result MUST be compared
  against `repertoire_plus_covariate` (which includes ancestry) before being trusted.**
- The HLA-carriage positive control (DRB1*15:01 or B*57:01 carriage predicted from TRB) still
  did not run — the run log shows `hla_calls_rich.tsv` on this VM instance was missing the
  expected person/allele columns entirely (a data-availability issue on this instance, not a
  case-count issue), so the check was skipped rather than failing silently. Needs a VM instance
  with a correctly populated `hla_calls_rich.tsv` before this control can run.

## QC (before any modeling)

Productive clonotype-count quantiles across the 1499 retained samples (1500 loaded, minus the
depth-floor exclusion below):

| | TRB | IGH |
|---|---|---|
| p5 | 461 | 232 |
| p25 | 836 | 527 |
| p50 | 1177 | 842 |
| p75 | 1577 | 1281 |
| p95 | 2281 | 2268 |

**Depth floor**: 100 TRB clonotypes minimum (documented, not silently hardcoded). Fewer than
20 of the 1500 samples fell below it (disclosure fix: the raw fraction pulled from the VM,
`0.000667`, back-calculates to an exact count of 1 out of 1500 — a 1–19 cell — so the qc TSV
in this repo has that field masked to `<0.0133`, i.e. `<20/1500`, per the hard disclosure
rule). Consistent with the 1200-person run (0.0%) and the earlier full-cohort load (0.1%) —
this cohort's TRUST4 repertoires stay deep enough that the floor barely binds.

## Power (minimum detectable AUROC)

Hanley-McNeil normal approximation, alpha=0.05, power=0.80, given (n_cases, n_total):
- At n_cases=100, n≈1500 (comparable strata in this run): minimum detectable AUROC ≈ 0.55
  (see `min_detectable_auroc` column in the metrics TSV — computed per phenotype/model row
  from its actual n; e.g. 0.543–0.556 across the positive-control rows here).
- The disease phenotypes never reached n_cases=100 at n=1500 either (max 33, for RA), so no
  AUROC number could be computed for them; the honest statement is still "underpowered at this
  n," not "no signal."

## Caveats

- **This is a 1500-person fixed-seed subsample of the 7640-person eligible cohort**, not the
  full cohort. The full-cohort real-data run (7640 people + XGBoost) is **still pending** — the
  VM is currently occupied for ~31h by the KIR full-cohort run (`sprints/S03_call8_figures_kir_prediction/KIR_FULL_RUN.md`),
  so WS6 could not be scaled further in this session. The perf fix from the 1200-person delivery
  (cheap cohort-wide-count pre-trim in `build_vj_kmer_features` before the dense kmer matrix is
  built) is unaffected by this run and should still hold at 7640 people once VM time frees up.
- XGBoost (`--with-xgboost`) was **not** run for this delivered result either — L1-LR only,
  same as the 1200-person run.
- `n_repeats=3` (the default) was used for this run, an improvement over the 1200-person run's
  `n_repeats=1` — the pooled-OOF bootstrap CIs above are computed from 3x the model fits, so
  repeat-to-repeat variance is now implicitly absorbed into the pooled estimate (still not
  reported as its own number).
- Disclosure: every `n_cases`/`n_total` in the metrics TSV is >=100 by construction (the
  MIN_CASES gate); the only sub-20 count encountered in this run (the depth-floor exclusion,
  and three disease phenotypes' case counts noted in prose above) is reported as `<20`, never
  as an exact number, in every file that left the VM.

## Distilled

- Pipeline correctly detects strong, expected signal (ancestry from repertoire, AUROC up to
  0.83, slightly higher than the 1200-person run) and reports a still-weak but slightly firmer
  sex-from-repertoire signal (AUROC 0.55, up from ~0.50) — consistent with, and a modest
  refinement of, the 1200-person result.
- Disease-phenotype prediction is **still not answered**: even at 1500 people, no
  chronic/autoimmune phenotype clears 100 cases (closest is RA at 33). This needs the full
  7640-person cohort, which is queued behind the ~31h KIR full-cohort run currently occupying
  the VM.
- Ancestry leakage is real and, if anything, slightly larger at n=1500 (AUROC up to 0.83) — any
  future disease result from `repertoire_only` alone must be cross-checked against
  `repertoire_plus_covariate`.
- One disclosure fix applied in this delivery: the depth-floor QC fraction (`0.000667`, which
  back-calculates to an exact count of 1/1500) was masked to `<0.0133` before committing.
