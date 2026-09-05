# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-09-05 — hla_popgen pipeline complete end-to-end on the real ~12,233-person cohort

All seven `scripts/hla_popgen/` steps (00 recon → 01 extraction → 02 cohorts → 03 novel alleles →
04 saturation → 05 frequency → 06 structure → 07 cross-cohort) have now run successfully against
real production data, driven live via Chrome computer-use on the Workbench VM (`18aa4228c0ec`,
`full-cohort-hla-calling` workspace, `AoU_Jupyter_ComputeEngine_20260805_big_run` instance — see
ENVIRONMENT.md quirks #31-33 for how this access works). 06 and 07 had been blocked on a
`pd.NA`-in-list crash (fixed, commit `421484d`) and a missing `--sr-genotypes` mount path (worked
around live: remounted gcsfuse, passed
`~/mnt/aou-controlled/v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv` explicitly —
**RUNBOOK.md does not yet document that 05/07 need the mount; fix this next session**).

### Real headline findings (all reports live at `~/repos/pilot-validation/reports/hla_popgen/` on
the VM, not yet copied anywhere else — see "Next" below)

- **Novel alleles (03):** 21,463 distinct novel-allele clusters with resolved gene identity;
  **1,190 pass every QC gate** (recurrence in ≥2 unrelated people, no disqualifying warning, not a
  homopolymer artifact). **90.8% of resolved clusters are protein-altering** (non-synonymous/
  frameshift) — real evidence of biology under balancing selection, not noise. Ancestry gradient
  in raw novel-call rate: AFR 0.380 > EAS 0.345 > SAS 0.316 > MID 0.307 > AMR 0.303 > EUR 0.267 —
  confirms the central IPD-reference-bias hypothesis directionally.
- **Saturation (04):** Rarefaction curves for "all alleles" and "novel alleles only" show **no
  plateau at all** even at full cohort size (novel-only curve is close to linear) — visually
  striking, strong support for "the HLA allele space is far from discovered." BUT: every pooled
  Chao2 richness estimate is flagged `LOW-CONFIDENCE (CI widened)` (Q2=0, too few exact-duplicate
  novel sequences yet to anchor the estimator), and the **ancestry-stratified "% discovered" table
  does not cleanly replicate the novel-rate gradient** (e.g. EUR is often NOT the highest %
  discovered; MID/SAS swing to extremes) — this is very plausibly an artifact of per-ancestry
  sample size, not a refutation. Real per-ancestry N (from `cohort_membership.tsv`, `cut -f9 |
  sort | uniq -c`): **AFR 4053, EUR 3109, AMR 2807, SAS 1262, EAS 1499, MID 499** (+23 blank/1 stray
  header row, quirk #30's known artifact). MID (n=499) and SAS (n=1262) are genuinely thin for a
  per-gene-per-ancestry Chao2 estimate — **treat the stratified saturation table as directional
  only until N grows or a less variance-hungry per-ancestry metric is used.**
- **Cross-cohort (07), template_distance by ancestry — clean confirmatory result:** % of haplotype
  calls at `template_distance<=0` (exact IPD match): AFR 78.9%, AMR 82.3%, EAS 78.1%, **EUR 84.6%**,
  MID 80.8%, SAS 79.5% — EUR highest, AFR/EAS lowest, exactly the hypothesized reference-bias
  signature, and it's a clean gradient (unlike the 04 stratified table above).
- **Cross-cohort (07), SR-vs-LR "EUR-common-allele" directional bias — UNEXPECTED, needs scrutiny
  before citing.** Every ancestry group's `(sr_freq - lr_freq)` bias on the 8 EUR-common alleles is
  **negative**, INCLUDING EUR itself (AFR -0.214%, AMR -1.389%, EAS -0.884%, **EUR -2.123%**, MID
  -3.111%, SAS -0.750%). The naive hypothesis ("SR over-calls EUR-common alleles specifically in
  non-EUR ancestries," i.e. positive bias for non-EUR, ~zero for EUR) does not hold in this simple
  form — EUR's own bias is more negative than most other groups. Plausible explanations not yet
  investigated: 2-field (SR) vs 4-field-truncated (LR) allele-string matching artifact, or the
  "EUR-common alleles" selection itself. **Do not report this as confirmed reference-bias evidence
  without digging into `07_figures_crosscohort.py`'s `compute_sr_lr_disagreement_by_ancestry()`
  first** — the plain `mean_abs_diff` column (AFR 0.616%, AMR 0.392%, EAS 0.495%, EUR 0.427%, MID
  0.771%, SAS 0.541%) is a cleaner, less confounded metric and roughly tracks ancestry as expected
  (MID/AFR highest disagreement) modulo MID's small N.

### Next

1. **Fix RUNBOOK.md**: step 5 (05/07) needs the gcsfuse mount + explicit `--sr-genotypes
   ~/mnt/aou-controlled/v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv` — currently
   undocumented, cost real time this session to discover live.
2. **Investigate the EUR-common-allele directional-bias sign** (`07_figures_crosscohort.py`,
   `compute_sr_lr_disagreement_by_ancestry()`) before using it in any writeup — check 2-field vs
   4-field truncation alignment between SR and LR allele strings specifically.
3. **Decide how to treat the 04 ancestry-stratified saturation table** given MID/SAS thinness —
   either pool to fewer/larger ancestry groups for this specific metric, gate it behind a minimum-N
   threshold the way `_viz_common.py`'s `MIN_CELL_N_PEOPLE` already does elsewhere, or explicitly
   caveat it in any report rather than presenting per-ancestry % discovered as reliable.
4. **Get real figures/reports off the VM into a form Marc and Cole can review** — nothing has left
   `~/repos/pilot-validation/reports/hla_popgen/` yet. Two paths already exist and are compliant:
   commit the aggregate-only report to the public GitHub repo (precedent: `reports/
   full_immuannot_lr_calling/*.png` already there), or copy to the existing share bucket
   (`gs://hla-calls-share-wb-cordial-leechee-9743/aggregate/`, ENVIRONMENT.md "Our own buckets").
   Neither has been done yet this session — decide which, since committing to the public repo and
   copying to the internal share bucket are different audiences (public vs. Cole-only).
5. Re-run 05 with the `--sr-genotypes` mount path too (not yet confirmed working the same way 07
   was) and sanity-check its output the same way 06/07 were checked here.
