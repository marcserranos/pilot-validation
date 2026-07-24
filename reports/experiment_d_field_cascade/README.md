# Experiment D field-level cascade — SpecImmune-truth and Immuannot-truth

Same re-scoring of Experiment D's n=60 AoU-native/SpecHLA calls (cumulative field 1→4
resolution, per gene), run against two different long-read truth sources:

- `specImmune_bench/` — vs. SpecImmune-LR truth (2026-07-20; script:
  `scripts/analyze_experiment_d_field_cascade.py`).
- `Immuannot_bench/` — vs. Immuannot truth (2026-07-24; script:
  `scripts/analyze_experiment_d_field_cascade_immuannot.py`). Immuannot-truth pool is smaller
  (403/480 slots — Immuannot abstains on the rest, doesn't miscall them).

**Reading the two side by side**: DRB1's Field-2 concordance looks notably better against
Immuannot-truth (89% AoU / 58% SpecHLA) than against SpecImmune-truth (77% / 52%) — likely
**not** a real accuracy difference. Immuannot's own confidence signal (see
`reports/confidence_matched_truth/`) shows DRB1 is by far its least-confident locus too, so
this is plausibly a selection effect in which DRB1 calls each tool manages to resolve, not
evidence AoU/SpecHLA are actually more correct there. DQA1 is the one finding that holds up
under **both** truth sources: AoU is the clear weak method there (59-64%), SpecHLA is not
(94-98%) — a genuine, cross-validated AoU-native issue, not a truth-source artifact.

The confidence-restricted, sample-size-matched follow-up (does the DRB1 gap survive once both
truth sources are held to the same strict confidence bar?) is in
**`reports/confidence_matched_truth/`** — not repeated here.
