# Critic #1 — fresh-context review of 37/38/39/41 (S03, 2026-09-23)

Reviewed: READMEs, figures (PNGs read directly), committed TSV/JSON, and the scripts' stated
statistics for `reports/hla_popgen/{37_dq_g1g2_signed_ld,38_hla_a_deletion_validation,
39_saturation_by_ancestry,41_kir_scoping}`. Cross-checked against
`sprints/CALL_SUMMARY_2026-09-22.md` and `sprints/S03_call8_figures_kir_prediction/LOG.md`.

## Ranked findings

1. **[blocker] Disclosure violation — explicit small counts in committed prose.**
   `reports/hla_popgen/41_kir_scoping/README.md:127` writes "2.8% co-occurrence (**1 haplotype of
   36**...)" — an explicit n=1 printed in a committed file, contradicting this project's own
   stated rule (see 37/38 READMEs: "counts 1-19 written as `<20`... never disclosed"). Line 187
   similarly states "n=11" for KIR2DL5B's protein-novelty denominator. Fix: redact both to `<20`
   (or report the ratio without the raw numerator), and re-audit 41 specifically since it has no
   TSVs to check against — the README prose *is* the data release here.

2. **[blocker] Uncorrected relatedness filter, stated and committed as fact.**
   `reports/hla_popgen/37_dq_g1g2_signed_ld/vm_run_summary.json:2` (`"n_people_unrelated":
   "13252"`) and the README's "VM run" section both assert 13,252 unrelated people. LOG.md's own
   `[~14:00] FLAG` entry says this is `cohort_membership.tsv`'s full row count (all people,
   relatives included) via the same `build_people()` that had a confirmed ancestry-column bug in
   this run. Every other S03 result (38, 39) uses 11,856 unrelated. This number is committed,
   unverified, and drives the headline "0/17,255 incompatible haplotypes" claim's cohort-size
   context. Fix: confirm whether `unrelated_ids` was actually applied before computing the 17,255/
   22,341 haplotype counts (LOG says the filter itself "is unaffected," but that claim is not
   independently checked against `cohort_membership.tsv`'s KING-relatedness column) and correct
   the label, or caveat explicitly.

3. **[major] Circularity risk in the mask-and-rephase EM (37d), not caveated.**
   The EM haplotype-frequency estimator is fit on "every unrelated person of that ancestry with
   two calls at each locus" (10,107 pooled) — which **includes** the 9,967-person physically-
   phased truth set it is then evaluated against (same purged, zero-cross-group population
   supplies both the prior and the test set). Finding "0 spurious incompatible haplotypes" from an
   EM trained on a pool that already contains zero cross-group haplotypes is closer to a
   tautology than an independent test. The README's own LOG entry anticipated this ("the EM learns
   the purged pool") but the published verdict text doesn't carry that caveat forward with the
   force it deserves. Fix (single most important missing analysis for 37): refit EM frequencies
   on a **held-out half** of each ancestry (or leave-one-out) and rephase the other half, or fit
   on a population known to include the light-blue signal (Cole's own SNP-array-phased data) if
   accessible.

4. **[major] Figure defect: overlapping labels, Nature would reject.**
   `reports/hla_popgen/38_hla_a_deletion_validation/fig_sr_validation.png`, panel (a): the "DRB1"
   and "B" point labels collide at the top of the plot (renders as an illegible "∆RB1"/"ARB1"
   glyph cluster). Needs label repulsion/offset (matplotlib `adjustText` or manual dx/dy) before
   this goes in any external report.

5. **[major] No multiple-testing correction stated for the 12-gene SR-calibration panel (38).**
   12 Fisher exact tests are reported (Table in §Q2/Q3) with p-values down to 1e-104, framed as a
   single calibration story — the effect sizes are large enough that correction wouldn't change
   the conclusion, but the README should say so explicitly (e.g. "Bonferroni α=0.05/12 still
   clears every reported test except DPA1") rather than silently reporting raw p only.

6. **[minor] Internal inconsistency: bootstrap p=0.0 still in committed TSV after a stated fix.**
   LOG.md `[~10:00] ERROR` says "p reported as 0.0000 from 200 reps" was flagged and fixed; the
   README now correctly prints "<0.005" in prose, but
   `reports/hla_popgen/39_saturation_by_ancestry/afr_vs_rest_bootstrap_test.tsv` still has literal
   `0.0` for AFR-AMR and AFR-EUR. The fix landed in the narrative layer only, not the underlying
   artifact — anyone quoting the TSV directly reproduces the original bug.

7. **[minor] 41's 20-person ancestry×platform cross-tab is a very thin stratification table.**
   "2 people each for AFR/AMR/EAS × {revio, sequel2e}, 1 each EUR/MID/SAS×sequel2e, 1
   UNASSIGNED/revio" (README §3) discloses cell-level counts of 1-2 for a 20-person pilot. Person
   IDs never left the VM, but per CLAUDE.md's own instruction this is exactly the kind of
   small-cell breakdown that should be flagged rather than decided — recommend collapsing to
   coarser bins (e.g. "roughly even across 6 ancestries and 2 platforms") in any external-facing
   version.

8. **[minor] Git hygiene: Cole's unpublished target figure sits untracked in the working tree.**
   `reference/cole_dq_g1g2_target_figure.pdf` is untracked (correctly, per LOG.md's own note "left
   untracked... Marc to decide") but nothing prevents an accidental `git add -A`/`.` from
   committing a supervisor's unpublished figure to this public repo. No other hygiene issues:
   `git ls-files | grep -iE 'per_person|participant|person_id'` returns nothing.

9. **[minor] 39: MID's Chao2 estimate (f0_hat=351) is presented alongside the five well-powered
   ancestries with only a passing caveat**, despite MID having by far the smallest N (487, its
   full sample) — an extrapolation this far past the observed range deserves a stronger "do not
   quote this number standalone" flag than the current one-line caveat.

10. **[minor] 39 fig1 panel (c) silently drops SAS and MID (no visible curve)** — caveated in text
    ("panel c's sampling units are only carriers of ≥1 novel allele") but not called out in the
    figure itself (e.g. via a "SAS/MID: insufficient novel-allele carriers" note in-panel), which
    a reader skimming only the figure would miss.

## Most important missing analysis per result

- **37**: held-out/leave-one-out EM refit (item 3) — the current mask-and-rephase result cannot
  distinguish "EM never manufactures incompatible pairs" from "EM was fit on data that already
  has none."
- **38**: an SR-genotype source independent of the pipeline already used elsewhere in this
  project, to rule out a shared SR-calling artifact driving both the "deletion" calls and their
  own contradiction test.
- **39**: an allele-frequency-weighted correction for IPD-IMGT's EUR submission bias in the
  %-of-catalogue-explored table — currently reported as a raw ratio with only a prose caveat.
- **41**: the framework-gene miss follow-up is good; still missing a per-haplotype look at what
  gene (if any) occupies the expected coordinates for "flanked" misses, to separate real deletion
  from mislabeling.

## Fixes applied (2026-09-23, sixth S03 agent)

1. **[blocker] KIR small counts (41)** — `reports/hla_popgen/41_kir_scoping/README.md`: masked
   "1 haplotype of 36" to `<20`, masked KIR2DL5B's n=11 to `<20`, masked every 1-19 `n calls`
   value in the novelty-decomposition table to `<20`, and collapsed the ancestry×platform
   cross-tab (§3 selection mix, both instances) to marginals only ("roughly even across 6
   ancestries, split revio(13)/sequel2e(7)"). No TSV/JSON in that folder — README prose was the
   only place counts lived.
2. **[blocker] "13,252 unrelated" mislabel (37)** — confirmed via VM: `37_vm_run.py`'s own
   `build_people()` starts from all of `cohort_membership.tsv` and greedily drops one member per
   related PAIR in isolation, not the maximal-independent-set algorithm scripts 29/38/39 use.
   Reran with the correct definition (`24_novelty_by_field.build_people`/`greedy_unrelated` on
   Table-1's own person_ids) → **11,856 unrelated, exactly matching 38/39**. New script
   `scripts/hla_popgen/37e_unrelated_fix_kfold_em.py`; corrected O/E table
   `oe_purge_table_by_ancestry_UNRELATED_FIXED.tsv` (headline result unchanged: still 0 observed
   cross-group haplotypes against ~10,469/8,066 expected). README's old "13,252" claim annotated
   as superseded in place; new "Critic #1 fix" section added with the corrected numbers.
3. **[major] EM mask-and-rephase circularity (37)** — same new script adds (a) 5-fold held-out
   EM (train on 4 folds, rephase+score the held-out fold only — no truth-set person ever scored
   by a model that saw their own data) and (b) a naive linkage-equilibrium baseline (independent
   marginal frequencies, ties broken by a seeded coin flip — the "zero real LD information"
   bound). Result: held-out EM is NOT perfectly circular-clean any more — 46/19,060 spurious
   incompatible haplotypes pooled (0.24%), vs. the LE-naive baseline's ~24% (~100x worse) —
   revised verdict in README. Unit tests: `tests/test_37e_unrelated_fix_kfold_em.py` (9/9 pass).
4. **[major] 38 fig_sr_validation panel (a) label collision** — replaced the fixed-offset
   annotate() calls with x-proximity clustering + vertical label dodging + leader lines
   (mirrors 39's `_dodge_label_positions` convention); "A"/"DRB1" no longer collide. Figure
   re-rendered locally from committed aggregates (`scripts/hla_popgen/38_hla_a_deletion_validation.py`
   `make_figure()`); panel (c) also gained a `contig_n_genes_hist_df` fallback path so it can be
   rebuilt from the disclosure-safe committed histogram without VM access.
5. **[major] 38 multiple-testing correction** — added explicit Bonferroni statement to README
   (12 genes, α=0.05/12=0.00417 — every test clears it except DPA1, already non-significant
   uncorrected) and fixed `sr_calibration_panel.tsv`'s literal `fisher_p=0.0` (DRB3) to `<1e-300`
   in both the TSV and README table.
6. **[minor] 39 residual `p=0.0`** — `afr_vs_rest_bootstrap_test.tsv`'s two literal `0.0` p-values
   (AFR-AMR, AFR-EUR) changed to `<0.005`, matching what the README prose already said.
7. **[major/blocker-adjacent] 39 equal-N discovery-slope comparison** — added, computed locally
   from the committed `curves.tsv` (no VM needed): slope of the mean rarefaction curve over
   [0.9·N*, N*] at N*=1,236 (5 well-powered ancestries) and N*=487 (all 6, incl. MID) —
   `equal_n_slope_comparison.tsv`. This flips the discovery-rate ranking back toward the call's
   original AFR-least-saturated expectation (AFR highest at N*=1,236: 101.6/1,000) versus the
   old own-N slopes (which are now demoted to a clearly-labelled secondary column). MID's Chao2
   flagged low-power directly in the headline table. Panel (c) of `fig1_main_saturation_panel`
   now carries an in-figure note ("MID: too few novel-allele carriers for a visible curve").
8. Not separately re-addressed: item 9 (MID Chao2) and item 10 (panel c note) — both folded into
   fix #7 above. Items 7/8 (thin KIR cross-tab, Cole's PDF hygiene) — #7 fixed in #1 above; #8
   (Cole's PDF) was already handled by LOG.md's `.git/info/exclude` fix, not re-touched here.
