# Confidence-matched truth comparison: AoU-native / SpecHLA vs. SpecImmune-LR and Immuannot, both restricted to high-confidence calls

**2026-07-24.** Extends `reports/experiment_d_field_cascade/` (AoU/SpecHLA vs. SpecImmune-LR
truth) and `reports/immuannot_pilot/` (AoU/SpecHLA vs. Immuannot truth, see
`scripts/analyze_experiment_d_field_cascade_immuannot.py`) with a single, objective, very
restrictive confidence bar applied to BOTH long-read truth sources before scoring AoU-native
and SpecHLA against them — Marc's explicit design constraints: **one global threshold per
tool (never tuned per gene), grounded in real sequencing-error math (not an arbitrary
number), and the two truth sets forced to the same retained share so neither comparison is
built on a systematically larger or smaller slice of the data.**

## Method

**Thresholds** — mathematically grounded (see `context/DECISIONS.md` / this session's
research pass), not guessed: sequencing/assembly error-rate math for PacBio HiFi/Revio, and
IPD-IMGT/HLA's real minimum inter-allele spacing (as low as 1 nucleotide between two
genuinely distinct catalogued alleles), argue for the strictest defensible bar on each
signal's own scale:
- **SpecImmune-LR**: worst-haplotype `si_identity >= 0.9995` AND not tied.
- **Immuannot**: worst-haplotype `template_distance = 0` (exact match to the nearest IMGT
  template) AND no `template_warning`.

Both applied identically across all 8 classical genes — the per-gene retention numbers below
are a *result* of that single global bar, not a per-gene knob.

**Sample-size matching**: SpecImmune resolves all 480 possible (person, gene) slots in this
60-person/8-gene cohort; Immuannot resolves only 403 (84.0%) — Immuannot **abstaining** on the
rest, a separate, already-documented behavior (see `reports/immuannot_pilot/README.md`), kept
deliberately distinct from the confidence-threshold effect. Given that, sample-size matching
is done by **equal absolute count**, rank-trimming whichever truth source's threshold-passing
set is larger down to the smaller one's N (keeping its own most-confident calls, not a random
subsample) — measured against the **one shared denominator of 480 possible slots** for both
tools, so the final retained *share* is identical by construction (both land at 167/480 =
34.79% here), not just approximately close.

## Headline results

**DRB1 is decimated by the confidence bar under both truth sources** — SpecImmune-truth drops
from 60 to 4 confident calls, Immuannot-truth from 51 to 3 (matched to 4). This is the sharpest,
most direct confirmation yet that DRB1 sits outside both tools' comfortable confidence zone;
read the resulting Field-2 percentages at DRB1 (n=4-8) as illustrative only, not statistically
load-bearing.

**AoU's DQA1 weakness survives strict truth-filtering** — AoU's Field-2 concordance at DQA1
is 64% (unfiltered vs. SpecImmune) / 59% (unfiltered vs. Immuannot), and stays the clear worst
non-DRB1 gene even after restricting to confident truth (73% vs. SpecImmune-confident, 69% vs.
Immuannot-confident) — both truth sources agree, before and after filtering, that this is a
real AoU-native problem, not an artifact of a noisy truth call. SpecHLA, by contrast, is
excellent at DQA1 under every condition (94-100%).

**Not everything moves the "expected" direction**: SpecHLA's DPB1 Field-2 concordance actually
*drops* when restricted to SpecImmune's most confident DPB1 calls (85% unfiltered -> 78%
confident, n=36) — a reminder that "more confident truth" does not mechanically mean "SR
methods look better against it," and that this deserves a closer look before assuming
confidence-filtering is a pure win.

## Figures

- `figures/field2_merged_confidence_matched.png` — the primary plot. Two subplots (top: vs
  SpecImmune-truth, bottom: vs Immuannot-truth); AoU-native and SpecHLA always sit as adjacent
  bars within the same gene group (never split across separate figures) so they're directly
  comparable; confidence-matched bars carry a bold black outline against the muted unfiltered
  color, so the high-confidence number is the one that visually pops.
- `figures/retention_by_gene.png` — diagnostic: per-gene retention at each of the 3 stages
  (all resolved -> after threshold -> after count-matching), for both truth sources. Purely
  informational — the threshold itself was never tuned per gene; this just shows which genes
  absorbed the cut (DRB1, overwhelmingly).

## Data & reproducibility

- Script: `scripts/analyze_confidence_matched_truth.py`. Run via
  `pixi run -e spechla -- python3 scripts/analyze_confidence_matched_truth.py <cohort.tsv> --si-identity-threshold 0.9995 --immuannot-distance-threshold 0`
- Reads the same per-person `comparison_log.csv`, `immuannot_calls.tsv`, and
  `hap{1,2}.gtf.gz` files as `reports/immuannot_pilot/` and
  `reports/experiment_d_field_cascade/` — no new VM run needed.
- Aggregate-only output (rates, counts, retention tables) — raw per-person calls stay on the
  Workbench VM per the standing egress rule (`context/DECISIONS.md`).
- Full field-level tables (Field 1 and Field 2, unfiltered and confidence-matched, both truth
  sources) are in the script's own markdown output
  (`~/pipeline_outputs/experiment_d/analysis_confidence_matched/confidence_matched_truth.md`
  on the VM) — not duplicated here; this README covers the headline reading, not every cell.
