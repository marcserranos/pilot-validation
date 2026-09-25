# Sprint S03: findings from call #8

*Written as a report to a PI. Each claim points to a committed result folder whose README has the
full method, the numbers and the caveats. Branch `s03-call8-figures-kir-prediction`, not pushed.
The chronological record of plans, runs, errors and fixes is [`LOG.md`](LOG.md).*

## Where each call item landed

| Call item (§ of CALL_SUMMARY) | Result | Status |
|---|---|---|
| §3 Recreate Cole's DQ G1/G2 figure | `37_dq_g1g2_signed_ld/` | done |
| §3 Light-blue cells: phase error or recombination? | same folder, `vm_*` files | done |
| §2 Per-person phasing-confidence file | kept on the VM only (`~/s03/results/37/`); only the schema and an aggregate summary are committed | done |
| §3 LD supplement CSV, with a separate colour for unobserved cells | `37_…/supp_table_ld_pairs.csv`, `supp_heatmap_dq_dp.*` | done |
| §4 Simpler colours for the deletion figure | `38_…/supp_deletions.*` | done |
| §4 Is the HLA-A deletion real? | `38_hla_a_deletion_validation/` | done |
| §6 Per-ancestry saturation curves and % of allele space explored | `39_saturation_by_ancestry/` | done |
| §6 Figure 1, panel by panel | `40_figure1_v5/` | done, all six panels |
| §7 KIR scoping and cost | `41_kir_scoping/`, `43_kir_full_cohort/` | done: full cohort called (12,261 people, 31 h, ~$95) |
| §9–10 BenchRep-T baseline | `42_repertoire_baseline/` | built; 1,500-person subsample run; full cohort pending |

---

## 1. DQ G1/G2: the purge is complete

**Question.** Cole's figure, a DQA1 × DQB1 heatmap of signed D′ sorted by the Petersdorf 2022
G1/G2 rule, is near −1 in the "predicted incompatible" quadrants, but a few cells are lighter. Are
those rare real recombinant haplotypes, or phasing errors?

**What we did.** We computed signed, allele-specific D′ from haplotypes where DQA1 and DQB1 sit on
the same assembled contig. There the phase is physically observed, not inferred.

**Result.**
- **0 of 17,255** physically phased DQA1–DQB1 haplotypes (4-field; 0/22,341 at 2-field) pair a G1
  chain with a G2 chain. This holds in the 11,856 unrelated people and in each of the six ancestries.
- If the two genes were independent we would expect ~170–2,000 such haplotypes per ancestry, so the
  observed/expected ratio is 0.
- Figure: `fig_dq_g1g2_committed_MAIN_POOLED.png`. Both incompatible quadrants are a uniform
  D′ = −1.
- The "bimodality" Cole noticed is explained by the same figure. In the compatible quadrants, true
  partner alleles sit near +1, and non-partner alleles of the same group sit near −1.

**Where the light-blue cells come from (mask-and-rephase experiment).** We took people whose phase
is physically known, hid it, re-inferred it statistically and counted the incompatible pairs the
inference created:
- A standard EM haplotype phaser, trained on other people and applied to held-out ones, created
  **0.24%** spurious incompatible haplotypes (46/19,060).
- A naive phaser that assumes the loci are independent created **~24%**, about 100× more.
- Stored hap1/hap2 labels in people whose DQA1 and DQB1 sit on *different* contigs also produce
  some (<20 of 320 pairs).

**Reading.** The G1/G2 incompatibility is enforced completely in vivo. The residual signal in
Cole's figure most likely reflects how his haplotypes were phased or imputed, not recombination.
This is also a methods point for the paper: physically phased long reads are what make the purge
visible cleanly.

## 2. HLA-A deletion: an artifact

**Question.** The raw HLA-A deletion rate is 1.83% of haplotypes. Cole: *"how is that plausible?"*

**Test.** For each person we called as missing HLA-A on one haplotype, we checked the AoU short-read
genotype. If short reads see **two different** HLA-A alleles, the gene is present on both
chromosomes and the deletion call is wrong.

**Result** (`fig_sr_validation.png`):

| Genes | Share of long-read deletion calls contradicted by short reads |
|---|---|
| HLA-A | **82.8%** (n = 412, Fisher P = 4 × 10⁻⁸) |
| B, DQB1, DRB1, DQA1 (genes that should never be deleted) | 73–83% |
| DRB3/4/5 (known real deletions) | <1% |

- All genes except DPA1 survive a Bonferroni correction across the 12 genes tested.
- About 30% of the "deleted" haplotypes have HLA-A annotated elsewhere in the same assembly, which
  points to fragmentation or misplacement in the assembly.
- **True HLA-A deletion rate ≤ ~0.3%.** The deletion caller should carry a short-read concordance
  filter for genes that are not known to vary in copy number.

## 3. Discovery curves: no ancestry is saturated; African ancestry discovers fastest at equal N

**Method.** Rarefaction over people: 25 random orderings, as in the Pakistan Genome Resource Fig. 3e.
Counts are distinct protein alleles across the 8 classical genes. See `39_…/fig1_main_saturation_panel.png`.

**Findings.**
- Every curve is still rising. The Chao2 estimator says several hundred alleles per ancestry remain
  undetected.
- Comparing ancestries at the same N (N = 1,236), **AFR discovers the most new alleles per extra
  1,000 people**: 101.6, against EAS 99.4, AMR 94.8, SAS 81.3, EUR 80.0.
- In total distinct alleles at equal N, AMR leads (AMR 320 vs AFR 287 at N = 487). This probably
  reflects admixture from three source populations. AMR still leads when people are restricted to
  ≥95% ancestry probability.

**Corrected along the way.** A first version fitted a saturating model and predicted *negative* new
alleles; it was replaced with non-parametric estimators. A later version compared slopes each taken
at that ancestry's own sample size, which is not comparable across ancestries; it now compares at
equal N.

## 4. KIR: called on the full cohort

**What ran.** Immuannot on the chr19 KIR window (54.60–54.92 Mb) of every main-tier long-read
assembly: 12,261 people on an 80-core machine, 20 people × 4 threads at a time, in 31.1 h
(~$95 at $2.98/h). HLA outputs were checksum-verified untouched before and after.
Results: [`43_kir_full_cohort/`](../../reports/hla_popgen/43_kir_full_cohort/README.md).

**Results** (11,882 unrelated people, 23,637 haplotypes with KIR calls):
- 99.5% of haplotypes carry ≥1 KIR call; **9.23 KIR genes per haplotype**. Both pilots replicate
  closely (20 people: 9.32; 170 people: 9.06).
- **58.9% of calls are novel** against IPD-KIR (pilots: 60.8%, 61.0%). As in the pilot, this is
  mostly outside the coding sequence. Protein-level novelty is 4–9% for most genes; KIR2DL5B is the
  outlier at 22.5%.
- Framework genes (3DL3, 2DL4, 3DP1, 3DL2) are present on 93.7–96.5% of haplotypes. 85% of the
  misses have called genes on both sides, so they are not truncated assemblies.
- Gene content differs by ancestry. The share of cA (inhibitory-type) haplotypes runs from 66.7%
  in East Asian ancestry to 41.4% in South Asian ancestry. This is the direction reported in the
  literature.

**Caveats.**
- **991 sequel2 people are excluded.** Their fallback alignment step would align chr19 contigs to
  a cached chr6 reference, which gives wrong calls. So the run covers the main tier only.
- **KIR2DL2 and KIR2DL3 co-occur on 0.8% of haplotypes.** These are alternatives at one locus, so
  the expected rate is ~0 (the pilot had 0/40). Likely a small miscall or duplication class; it
  needs a targeted look.
- **104 people (0.85%) lack a call on at least one haplotype.** The pattern repeats identically on
  retry. Their other haplotype is kept.
- **Novelty is not yet cross-checked by recurrence.** This is the equivalent of the HLA script-03
  check.

**Figure:** `43_kir_full_cohort/fig_kir_full_cohort.png`, six panels: presence per gene; novelty
decomposition; presence and protein novelty by ancestry; cA/cB by ancestry; QC.

## 5. Prediction baseline (BenchRep-T protocol)

**Protocol.** Features are V/J gene usage plus CDR3 4-mers. Models are L1-regularised logistic
regression and XGBoost, evaluated with 3-fold stratified cross-validation on identical folds and
pooled out-of-fold AUROC/AUPRC. Covariates include EHR depth, because Aleix found a
healthcare-access gradient across ancestries.

**Cohort.** 7,640 unrelated people with both repertoire and long-read HLA data.

**Result so far** (1,500-person subsample, 3 repeats of 3-fold CV, L1-LR; b9ff562):

| Prediction target | AUROC (95% CI) | Earlier 1,200-person run |
|---|---|---|
| Ancestry, AFR | 0.82 (0.80–0.85) | 0.81 |
| Ancestry, EUR | 0.67 (0.65–0.71) | 0.71 |
| Sex | 0.55 (0.52–0.58) | 0.50 |

No disease reached 100 cases in the subsample; rheumatoid arthritis came closest, with 33.

**Takeaway.** The repertoire encodes ancestry strongly. Any disease model has to condition on
ancestry, or it will partly be an ancestry classifier.

**Next.** Full-cohort run with XGBoost; the HLA-allele positive control.

## 6. Figure 1 v5

`40_figure1_v5/figure1_v5.png` is 183 × 150 mm. Changes against Cole's feedback:

| Panel | Content | Change | Status |
|---|---|---|---|
| a | Admixture strip | ≥98% annotations removed | done |
| b | HLA-B ternary | 0–100 ticks; ≥0.98 ancestry filter; ≥20 carriers | done |
| c | Ancestry × gene novelty heatmap (the panel Cole preferred) | redrawn cleanly | done |
| d | Protein alleles not in IPD-IMGT/HLA | now a marginal of c | done |
| e | Discovery curves | new | done |
| f | DQ purge (observed vs expected) | new | done |

Update (09-24): panels a and c are now filled from two VM exports (0ec28ff). Panel a shows 589
bins of ≥20 unrelated people (11,833 in total). Panel c shows the any-field novelty rate by gene ×
ancestry; only HLA-A × MID is censored, and it is hatched, not drawn as 0. The highest cell is
DRB1 × EAS at 70.9%. Panel c uses a 0.9 ancestry threshold while panel b uses 0.98; this is
documented in the README. A draft legend is in the README.

## 7. Suggested next steps

1. **The DQ result, plus the phasing experiment, is a candidate main-text panel** for the
   stat-gen half. Physical phasing reveals complete purging; inferred phasing creates the residue.
   Worth showing Cole first, since he ran the same analysis.
2. KIR is done on the full cohort. Next: a recurrence check on novel KIR alleles, a look at the
   2DL2/2DL3 co-occurrence haplotypes, and a fix to the self-align reference bug so the 991 sequel2
   people can be added.
3. Finish WS6 on the full cohort, and add ancestry-matched controls before looking at any disease
   AUROC.
4. Competitive landscape: the AoU long-read flagship preprint (medRxiv 2025.10.02.25336942) and
   FuFiHLA (a long-read HLA typing tool) are the closest work. Worth reading before framing.

## 8. Open items and caveats

- **Disclosure.** Every committed count of 1–19 is written as `<20`, and a critic sweep checked this.
  Small-cell publication policy is still an open question in DECISIONS.md.
- **Cole's reference PDF** is not committed. It is excluded locally, since this is a public repo.
- **References.** The ScienceDirect reference is probably Czarnecki 2010 (unconfirmed). The
  "summer 2026 HLA copy-number" paper was not found; ask Cole for it.
- **Remaining VM steps.** Only WS6 on the full cohort (7,640 people + XGBoost). Everything else in [`VM_HANDOFF.md`](VM_HANDOFF.md) is done.
