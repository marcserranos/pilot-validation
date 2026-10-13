# 64 PUBLIC tables (VM run 2026-10-12)

Pulled from `~/s05/results/64/` on the VM, one file at a time. Byte sizes were checked against the VM.

How they were made:
- Script `scripts/hla_popgen/64_novel_allele_validation.py` at e9502c2 (adds the pooled key-level verdict).
- 8 classical genes (A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1), novelty depths 2 and 3.
- Unrelated people only (cached LR relatedness pairs, `_s05_common` fallback).
- IPD-IMGT/HLA 3.65.0, compared against 3.55.0. Long-read support from the GRCh38 HiFi BAMs; no short reads.
- Default thresholds. Run 10:02–11:33 UTC with 2 workers: 4,616 calls, 4,299 keys, 3,419 people.

Disclosure:
- No person ids, and no per-person or per-call rows.
- Participant counts (`n_calls`, `n_calls_tested`) are masked `<20` by the script.
- Key counts (`n_keys` and the per-step key columns of `funnel_with_evidence.tsv`) count distinct alleles, not people, and are left unmasked. This follows the looser rule for non-participant counts (DECISIONS 2026-09-28 / 2026-10-11): `<20` is applied at publication.
- `now_known_alleles.tsv` (named alleles) is NOT included here, pending Marc's decision on publishing it.
- 2026-10-13: `failures_summary.tsv` re-masked by `scripts/hla_popgen/66_cross_table_suppression.py` (the reasons summed back to the totals). `summary_by_stratum.tsv` still has recoverable cells (verdicts sum to tested calls, hp_only + other = all, platforms/ancestries/genes sum to ALL); not yet re-masked. Earlier versions remain in git history.

Re-match-only run (no reads), for comparison: `~/s05/results/64_ipd_unrel` on the VM; same funnel counts minus the read columns.
