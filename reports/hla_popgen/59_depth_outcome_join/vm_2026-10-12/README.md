# 59 VM run, 2026-10-12: PUBLIC tables

How it was made:
- `59_depth_outcome_join.py`, rerun at 17:14 UTC at main 5fa84b0 (cross-row complementary suppression, PR #25).
- Input: 58's INTERNAL per-person depths.
- 11,856 unrelated people (cached LR relatedness pairs). statsmodels absent, so the numpy logistic fallback was used.

Result: low regional depth explains none of the outcomes.
- Every outcome is `depth_does_not_explain`; both region proposals are `no_outcome_explained_by_depth`; `exclude_people.tsv` excludes 0 people.
- The largest admissible normalised cutoffs (0.85 KIR, 0.90 MHC, both bound by AFR) would remove only 4–9% of the KIR events and about 4% of the HLA-absent events.

Disclosure:
- No ids. People and event counts 1–19 are masked by the script.
- Checked for recovery by subtraction:
  - hla_gene_absent_by_bin: no single-masked partition.
  - outcome_by_depth_bin: 46 partitions have one masked n_events cell. 39 have no published exact total; in the other 7 the implied value is not a 1–19 count.
  - cutoff_scan: every masked n_events_removed belongs to a key with no exact total.
  - So no masked cell can be recovered.
