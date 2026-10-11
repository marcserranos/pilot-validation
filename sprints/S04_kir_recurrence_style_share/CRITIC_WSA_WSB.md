# Critic — fresh-context review of S04 WS-A/WS-B (44/46, v4b), 2026-09-28

Reviewed: `reports/hla_popgen/44_kir_recurrence_saturation/` and
`reports/hla_popgen/46_kir_vs_hla_catalogue/` (READMEs, all TSVs, all PNGs read at full size) and
`scripts/hla_popgen/44_kir_recurrence_saturation.py`/`45_kir_recurrence_figure.py`/
`46_kir_vs_hla_catalogue.py`, against `sprints/S04_kir_recurrence_style_share/AGENT_PREAMBLE.md`,
`reference/FIGURE_STYLE.md`, and the last ~15 `LOG.md` entries (v1→v4b bug history).

## Blockers: 0. Major: 1 (fixed). Minor: 4 (fixed).

## Reconciled headline numbers (pooled-ALL, n=11,856; use these)

| level | HLA S_obs | KIR S_obs |
|---|---:|---:|
| genomic (upper bound) | 17,188 | 38,024 |
| CDS | 1,363 | 2,100 |
| protein | 1,079 | 1,444 |
| protein_novel | 198 | 1,063 |

Verified by summing three independent committed tables (`coverage_chao2.tsv`,
`diagnostics_identity.tsv`, `saturation_curves.tsv` via `pooled_recurrence_from_curves()`) — all
three agree exactly. The operator's LOG-only console numbers (163 / 1,441 / 1,061) are a stale
pre-final-pull printout, not present in any committed artifact; now stated explicitly in both
READMEs so a reader doesn't average or "split the difference."

## Findings

1. **[major, fixed] `fig_catalogue_completeness.png` panel (a) had ~20 gene labels overlapping
   each other and sitting on the wrong data points.** Both species' genomic-level points cluster
   near x=90–100%; free per-point label repulsion couldn't separate them (a text-vs-marker
   collision `check_layout(strict=True)`'s text-vs-text/decoration linter cannot see). Fixed by
   extending `_label_curve_ends` (in `45_kir_recurrence_figure.py`) with an optional `label_x`
   override, then chaining two shared-label columns (HLA then KIR) in `46`'s `_scatter_panel` so
   each species gets its own non-overlapping leader-line column. Re-rendered and visually
   re-verified at full size (crop-checked, no residual overlap). Added 2 clustered synthetic
   fixtures to `test_46_kir_vs_hla_catalogue.py` so the layout tests exercise realistic (not
   arbitrarily wide-spread) data — the wide-spread pre-existing fixture surfaced a real
   long-leader-line-sweep edge case that cannot occur on this pipeline's actual data (verified:
   min real pct_novel_any is 61.9%).
2. **[minor, fixed] Disclosure: `artifact_qc.tsv`, `genomic_artifact_qc.tsv`,
   `diagnostics_identity.tsv`, `kir_cds_match_qc.tsv` had unmasked 1–19 call/allele counts**
   (e.g. `KIR2DS2 frameshift_or_stop=2`). Masked to `<20` in place (cheap, matches the existing
   `recurrence_classes.tsv` convention). Confirmed these tables source no headline number, so
   nothing else needed re-deriving. `coverage_chao2.tsv` was deliberately left as-is (already
   flagged/justified in Caveat 4 — masking it would break the headline sums).
3. **[minor, fixed] Numeric discrepancy vs LOG's operator printout** (198/1,444/1,063 vs
   163/1,441/1,061) traced and documented (see table above) in both READMEs.
4. **[minor, fixed] Simpson's-paradox mechanism, per Marc's specific question.** Tested whether
   the pooled-vs-per-ancestry reversal is just unequal N: N-weighting the 5 ancestries' own
   protein completeness gives HLA 52.9% vs KIR 36.2% — essentially the same ~17pp HLA lead as the
   raw per-ancestry table, NOT the pooled 0.7pp near-tie. So it is **not** an unequal-N artifact;
   it's a genuine Chao2/S_obs nonlinear-pooling effect from population structure. Documented in
   both READMEs with the worked check.
5. **[minor, fixed]** Cross-referenced the masking-consistency rationale between `44` and `46`
   READMEs so the two documents don't contradict each other on why some tables are masked and
   others aren't.

## Checked, no issue found

Chao2 (Chao 1987) and Good-Turing (Chao & Jost 2012) formulas match their cited references;
coverage is correctly labeled "incidence coverage," never "fraction of allele space explored."
Two-proportion z-test is textbook-correct, pooled-variance, two-sided; no multiple-testing
correction is honestly disclosed. Genomic-level upper-bound caveat is prominent and no headline
uses genomic numbers. All PNGs viewed at full size: no other overlaps, honest zero-based axes,
NA/pseudogene genes shown as explicit notes never as 0, recurrence classes 1/2/>2/≥20 visible with
direct labels, per-ancestry panels present, both novelty levels shown side by side per Marc's ask.

## Judgment calls for Marc

- Masked the 4 new QC tables' small cells rather than leaving them — reversible, cheap, no
  downstream number depends on them.
- Did not re-derive the pooling mechanism further (e.g. decompose ancestry-private allele share)
  — flagged as a follow-up, not fabricated.
- `equal_n_slope.tsv` pooled-vs-summed mismatch (pre-existing, predates v3) left open per prior
  critic's scope call.

## Tests

`pytest scripts/hla_popgen/tests -q`: 628 passed, 9 errors (pre-existing, unrelated fixture
errors — same baseline as every prior LOG entry).
