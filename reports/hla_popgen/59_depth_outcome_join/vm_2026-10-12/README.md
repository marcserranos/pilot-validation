# 59 VM run, 2026-10-12: PUBLIC tables

How it was made:
- `59_depth_outcome_join.py`, rerun at 17:14 UTC at main 5fa84b0 (cross-row complementary suppression, PR #25).
- Input: 58's INTERNAL per-person depths.
- 11,856 unrelated people (cached LR relatedness pairs). statsmodels absent, so the numpy logistic fallback was used.

Result (corrected 2026-10-12 after an independent review):
- No coverage cutoff is worth applying. The largest admissible normalised cutoffs (0.85 KIR, 0.90 MHC, both bound by AFR) would
  remove only 4–9% of the KIR events and about 4% of the HLA-absent events, while removing 3.5–4.7% of people. `exclude_people.tsv`
  excludes 0 people.
- The pre-set rule (a cutoff must remove ≥ 50% of an outcome's events while removing ≤ 10% of people) could not be met by any
  outcome this common: with prevalence p, at most 0.10/p of the events can be removed (KIR incomplete, p ≈ 38%: at most 26%).
  So `depth_does_not_explain` is guaranteed by arithmetic for the common outcomes and is not evidence of "no association".
- Depth IS associated with the KIR and MHC outcomes (`logistic_or.tsv`, OR per 0.1 lower normalised depth): KIR incomplete 1.16
  (p 1e-12), KIR framework gene missing 1.47 (p 2e-58), framework missing on complete haplotypes 2.3 (p 6e-62), HLA absent 1.09
  (p 0.02), MHC fragmented 1.12 (p 0.002); DQ outcomes not associated. The KIR-incomplete rate rises from 34% at normalised depth
  1.0–1.1 to 55% at 0.7–0.8 (`outcome_by_depth_bin.tsv`).
- Part of this is reverse causation: a haplotype that lacks genes (contracted KIR haplotypes, HLA gene deletions) has less sequence
  matching the reference window, so its window depth is lower. The strongest association is for framework genes missing on
  complete haplotypes, where no assembly break is involved, which fits this. Depth is therefore a marker of haplotype structure as
  much as of coverage, and a depth filter would remove real structural haplotypes.

Disclosure:
- No ids. People and event counts 1–19 are masked by the script.
- Checked for recovery by subtraction:
  - hla_gene_absent_by_bin: no single-masked partition.
  - outcome_by_depth_bin: 46 partitions have one masked n_events cell. 39 have no published exact total; in the other 7 the implied value is not a 1–19 count.
  - cutoff_scan: every masked n_events_removed belongs to a key with no exact total.
  - So no masked cell can be recovered.
