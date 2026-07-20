# HLA x disease sanity check — Celiac and narcolepsy on AoU-native calls

**Summary.** Two known, large-effect HLA-disease associations were tested against AoU-native's
population-scale calls (535,658 people) as a first-principles validation of the whole pipeline:
calling quality, the BigQuery diagnosis join, and the phenotype extraction. Both diseases show
the expected allele in a hand-curated positive-control check and in an untargeted model given no
disease-specific hints. Cross-validated (out-of-fold) performance tells a different story per
disease: Celiac shows real, modest, monotonic predictive signal with a substantial ancestry gap
(European AUROC 0.681 vs African AUROC 0.542, both well-powered); narcolepsy shows essentially
chance-level prediction (pooled AUROC 0.592, non-monotonic decile gradient) despite the model
correctly landing on a haplotype partner of the true causal allele — a concrete demonstration that
finding the right biology and predicting well on unseen people are different questions.

Full raw output: [`cv_results.log`](cv_results.log) (verbatim from the VM run, both diseases).

---

## Method

- **Data:** AoU-native `hla_genotypes.tsv` (all 8 classical genes, 2-field), joined to genetic
  ancestry and to a case/control label pulled from `condition_occurrence` via a single SNOMED
  standard concept per disease (Celiac: 194992; Narcolepsy with cataplexy: 437854 — deliberately
  distinct from "without cataplexy"). Raw ICD/SNOMED-code case definition, not a validated
  phenotype algorithm — a real limitation, flagged rather than fixed, for this first pass.
- **(A) Positive control:** hand-curated carrier flag for the known risk marker (Celiac: DQ2.5 or
  DQ8, unphased presence; Narcolepsy: DQB1\*06:02 dosage) vs. case/control, 2x2 table + Fisher's
  exact, pooled and per ancestry.
- **(B1) Generalization:** every 2-field allele at all 8 genes (dosage-encoded) plus ancestry
  dummies -> L1-penalized logistic regression, 5-fold stratified cross-validation, out-of-fold
  predictions (scaling fit per-fold on the training split only, to avoid leakage). AUROC, AUPRC,
  a risk-decile enrichment table, and PPV/sensitivity/specificity/NPV at top-5%/10% thresholds —
  pooled and per ancestry (gated at ≥10 cases and ≥10 controls per group; below that, flagged
  "insufficient" rather than given a number).
- **(B2) Interpretation:** one final fit on all data, top +/- coefficients — a face-validity check
  ("does the model find real biology"), explicitly not a performance claim.

## Celiac — real, modest signal; a clear and well-powered ancestry gap

- Pooled AUROC 0.740, "fair" discrimination. Risk-decile table is clean and monotonic (2.29% ->
  1.03% -> ... -> 0.13% observed case rate from decile 1 to 10) — the model's risk ranking is
  real, not noise. Flagging the top 10% as high-risk gives 3.8x enrichment and captures 38% of
  true cases (PPV 2.29% vs 0.61% baseline; NPV 99.58%, consistent with celiac genetics' known
  strength at ruling out disease rather than ruling it in).
- **Ancestry gap, well-powered on both sides:** EUR AUROC 0.681 (n=2,776 cases) vs AFR AUROC 0.542
  (n=138 cases) — AFR sits barely above chance. AFR's within-group top-10% sensitivity is 14.5%
  vs EUR's 30.6%. This is a concrete, measured instance of the polygenic-score ancestry-transfer
  problem already flagged as a strategic concern in `DECISIONS.md` — not hypothetical here.
  AMR (0.606, n=285) and EAS (0.595, n=28) sit between. MID's 0.763 (n=10) sits exactly at the
  reporting gate and should be treated as noisy, not as the best-performing group.
- (B2): top features are DQB1\*02:01, DRB1\*03:01, DQA1\*05:01 (DR3-DQ2.5), then A\*03:01, B\*08:01 —
  fragments of the well-known "8.1 ancestral haplotype." Real, textbook biology, found with zero
  hints.

## Narcolepsy — essentially chance-level prediction, despite finding real biology

- Pooled AUROC 0.592, barely above 0.5. AUPRC 0.001 sits at baseline (0.109%) — almost no
  enrichment. The risk-decile table is not monotonic (0.172 -> 0.146 -> 0.140 -> 0.142 -> 0.106 ->
  0.091 -> 0.075 -> 0.088 -> 0.056 -> 0.071%) — the signature of noise, not a real gradient.
- Even EUR, the best-powered group (446 cases), only reaches AUROC 0.521 — essentially
  uninformative. AFR/AMR sit at/below 0.5; EAS's 0.286 (n=13) is likely small-n noise. MID/SAS
  correctly gated as insufficient (2 cases each).
- (B2) still finds DRB1\*15:01 as the top feature — a real haplotype partner of the true causal
  allele, DQB1\*06:02 (DR15-DQ6 haplotype), even though DQB1\*06:02 itself doesn't reach the top 15.
- **Read together, (B1) and (B2) demonstrate the central methodological point directly:** a model
  can independently rediscover the correct biology (B2) while still failing to predict well on
  unseen people (B1). The most likely explanation is the loose, code-only case definition (see
  Method) — narcolepsy-with-cataplexy's raw ICD-based prevalence here (~127/100,000) already ran
  roughly 2x even the "claims data runs hot" literature estimate before this analysis was run,
  suggesting the case group is diluted with non-NT1 or miscoded people. Untested, logged as the
  natural next step if this thread is picked up again: require the diagnosis code on ≥2 separate
  encounter dates before counting someone as a case, and re-run.

## Bottom line

The pipeline (calling -> BigQuery join -> phenotype extraction -> modeling) is demonstrably sound:
both diseases recover known biology with zero hints, and Celiac shows real, validated,
out-of-sample predictive signal. The two findings worth carrying forward are the measured
ancestry-transfer gap for Celiac (a concrete data point for the diversity/PRS-transfer aims) and
the feature-recovery-vs-performance gap for narcolepsy (a concrete illustration of why coefficient
inspection alone was never sufficient — and the reason this check now runs cross-validated by
default, not just as an interpretation pass).
