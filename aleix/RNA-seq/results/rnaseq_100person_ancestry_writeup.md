# RNA-seq repertoire recovery across ancestry — 100-person pilot — 2026-08-11

**Headline: recovery is NOT even across ancestry.** AFR shows the highest mean CDR3
recovery (7,406), EAS the lowest (4,538) — a 1.63x spread. Whether this reflects a real
recovery/reference-bias effect or a cohort sequencing-depth confound is **not yet
determined** — see "Open question" below before treating this as a bias finding.

## Methodology

- 100 people, ancestry-stratified (17/17/17/17/16/16 across AFR/AMR/EAS/EUR/MID/SAS),
  built via `../scripts/build_rnaseq_cohort.py --total 100` — ground-truth-verified
  against the mount, not just trusting the manifest.
- TRUST4 v1.1.5-r573, `--abnormalUnmapFlag` fix applied (see
  `1000291_trust4_smoke_test.md`), run via `../scripts/run_rnaseq_batch.sh`, mostly at
  `--jobs 12` after an efficiency sweep (`../scripts/run_rnaseq_jobs_sweep.sh`) showed
  parallelism gives a solid ~2-4x speedup over sequential but couldn't cleanly rank
  exact `--jobs` values against each other (population noise between non-overlapping
  test slices dominated the signal — a real methodological finding in its own right,
  logged in the sweep script's own comments).
- Aggregated with `../scripts/aggregate_rnaseq_results.py` — full table:
  `rnaseq_cohort_ancestry_summary.csv` (this directory).

## The numbers

| Ancestry | n | Mean CDR3s | Median | Min | Max |
|---|---|---|---|---|---|
| AFR | 17 | 7,406 | 7,000 | 2,310 | 17,015 |
| AMR | 17 | 6,752 | 4,841 | 1,546 | 19,176 |
| SAS | 16 | 6,407 | 6,390 | 2,187 | 10,225 |
| EUR | 17 | 5,725 | 5,378 | 926 | 11,700 |
| MID | 16 | 5,589 | 5,614 | 1,783 | 11,740 |
| EAS | 17 | 4,538 | 4,045 | 697 | 9,570 |

**AMR's mean vs median gap (6,752 vs 4,841) flags a right-skewed distribution** — a few
high-CDR3 people are pulling that group's mean up; the typical AMR person looks more
like the median than the mean.

## Why this doesn't look like pure noise

Chain-level breakdown (TRA/TRB/IGH/IGK/IGL — the 5 chains with substantial counts;
TRG/TRD are too low-count to read much into): **AFR is highest, and EAS is lowest, in
every one of the 5 major chain types.** The same ancestry group winning or losing across
five largely independent measurements is much less likely to be chance than a single
pairwise comparison would suggest.

## Open question — resolve before calling this a bias finding

Two very different explanations produce the same observed pattern:

1. **A recovery/reference-bias effect** — the same shape of problem as the `noHLA`
   reference issue (see the AoU RNA-seq dossier) and the DQA1 HLA discordance already
   found in the main HLA-Resolve workstream. If TRUST4 or the underlying STAR alignment
   genuinely recovers less signal for some ancestries, that's a real finding about the
   pipeline/reference.
2. **A cohort sequencing-depth confound** — if AoU's actual RNA-seq depth (or blood
   cell composition) happens to vary by ancestry group in this specific 100-person
   draw, that alone would produce this exact pattern with nothing to do with reference
   bias.

**Not yet distinguished.** Next check: `../scripts/check_depth_confound.py`, comparing
raw CDR3 recovery against depth-normalized recovery (CDR3s per million read pairs) by
ancestry group. If the gap shrinks substantially once normalized by depth, this is a
cohort-composition effect. If it persists after normalizing, that's the more concerning
finding — and the one worth taking seriously.

## Caveats

- n=16–17 per group; within-group ranges are large (e.g. AFR spans 2,310–17,015 —
  7.4x), so individual variance is substantial relative to the between-group gap. No
  formal significance test run yet.
- Cohort pick is deterministic (sorted by research_id within each ancestry group), not
  randomized — unlikely to matter but not ruled out as a selection effect.
- This is one pipeline (TRUST4) and one library prep (AoU's Watchmaker + Polaris
  Depletion, whole blood) — doesn't generalize to repertoire recovery methods generally.
