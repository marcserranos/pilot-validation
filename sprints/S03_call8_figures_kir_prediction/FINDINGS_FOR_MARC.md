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
| §6 Figure 1, panel by panel | `40_figure1_v5/` | panels b, d, e, f done; a and c need one VM run |
| §7 KIR scoping and cost | `41_kir_scoping/` | pilot done, GO; the full run needs your resize decision |
| §9–10 BenchRep-T baseline | `42_repertoire_baseline/` | built; only a 1,200-person subsample has run |

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

## 4. KIR: feasible now

- **Region.** chr19:54.60–54.92 Mb. Immuannot already calls KIR from the same assemblies; only the
  region flag changes.
- **20-person pilot.** All 20 people (40 haplotypes) succeeded, with ~9.3 KIR genes per haplotype.
- **Sanity checks.** KIR2DL2 and KIR2DL3 never co-occur on a haplotype, as expected. Framework genes
  are present on 85–97.5% of haplotypes, and none of the gaps are due to fragmentation.
- **Novelty.** 60.8% of calls are "novel" at the full-sequence level, but mostly outside the coding
  sequence. Protein-level novelty is 0–18% per gene.
- **Full cohort.** ~15 h / ~$45 on a 96-core machine, the same machine type as the HLA run.
  **Decision for you:** resize the app. Resizing restarts the VM.
- A 170-person ancestry-stratified extension was launched and still needs aggregating.

## 5. Prediction baseline (BenchRep-T protocol)

**Protocol.** Features are V/J gene usage plus CDR3 4-mers. Models are L1-regularised logistic
regression and XGBoost, evaluated with 3-fold stratified cross-validation on identical folds and
pooled out-of-fold AUROC/AUPRC. Covariates include EHR depth, because Aleix found a
healthcare-access gradient across ancestries.

**Cohort.** 7,640 unrelated people with both repertoire and long-read HLA data.

**Result so far** (1,200-person subsample, `42_…_auroc.png`):

| Prediction target | AUROC |
|---|---|
| Ancestry, AFR | 0.81 |
| Ancestry, EUR | 0.71 |
| Sex | 0.50 |

No disease reached 100 cases in the subsample.

**Takeaway.** The repertoire encodes ancestry strongly. Any disease model has to condition on
ancestry, or it will partly be an ancestry classifier.

**Next.** Full-cohort run with XGBoost; the HLA-allele positive control.

## 6. Figure 1 v5

`40_figure1_v5/figure1_v5.png` is 183 × 150 mm. Changes against Cole's feedback:

| Panel | Content | Change | Status |
|---|---|---|---|
| a | Admixture strip | ≥98% annotations removed | needs VM export |
| b | HLA-B ternary | 0–100 ticks; ≥0.98 ancestry filter; ≥20 carriers | done |
| c | Ancestry × gene novelty heatmap (the panel Cole preferred) | redrawn cleanly | needs VM export |
| d | Protein alleles not in IPD-IMGT/HLA | now a marginal of c | done |
| e | Discovery curves | new | done |
| f | DQ purge (observed vs expected) | new | done |

Panels a and c need one VM run. Panel a is exported as bins of ≥20 people, so it is
disclosure-safe. Panel c needs uncensored totals; without them, most cells would show only an upper
bound. Copy-paste commands: [`VM_HANDOFF.md`](VM_HANDOFF.md) (a). A draft legend is in the README.

## 7. Suggested next steps

1. **The DQ result, plus the phasing experiment, is a candidate main-text panel** for the
   stat-gen half. Physical phasing reveals complete purging; inferred phasing creates the residue.
   Worth showing Cole first, since he ran the same analysis.
2. Resize the VM and run KIR on the full cohort. A KIR panel in Figure 1 is then one script away.
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
- **Remaining VM steps.** All are listed in [`VM_HANDOFF.md`](VM_HANDOFF.md).
