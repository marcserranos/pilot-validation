# RNA-seq repertoire recovery across ancestry — 100-person pilot — 2026-08-11

**Headline: recovery is NOT even across ancestry, and depth only partly explains it.**
AFR shows the highest mean CDR3 recovery (7,406), EAS the lowest (4,538) — a 1.63x
spread. Normalizing by each person's actual sequencing depth shrinks that to 1.38x —
depth is a real, partial contributor, but a meaningful gap survives after controlling
for it. See "Depth-confound check" below — this is resolved to first order, not fully
closed.

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

## Depth-confound check — run 2026-08-11, `../scripts/check_depth_confound.py`

Two explanations could produce the same raw pattern: a recovery/reference-bias effect
(the same shape of problem as the `noHLA` reference issue and the DQA1 HLA discordance
already found in the main HLA-Resolve workstream), or a cohort sequencing-depth
confound (AoU's actual RNA-seq depth happening to vary by ancestry in this specific
100-person draw). This check distinguishes them by normalizing CDR3 recovery by each
person's actual depth (from `rnaseq_metadata.tsv`) and re-comparing by ancestry group.

**Result: depth is a real, partial contributor — not the full story.**

| | Raw CDR3 mean | Depth-normalized (CDR3s per million read pairs) |
|---|---|---|
| AFR (highest) | 7,406 | 53.7 |
| EAS (lowest) | 4,538 | 39.0 |
| **Ratio** | **1.63x** | **1.38x** |

- Normalizing shrinks the AFR-vs-EAS ratio from 1.63x to 1.38x — roughly a **40%
  reduction** of the excess gap. Real, but incomplete.
- **Depth itself varies by ancestry in this cohort**: AFR has the highest mean depth
  (146.2M read pairs), EAS the lowest (121.4M) — a 1.20x spread, an independent finding
  about this draw worth knowing on its own.
- **AFR and EAS remain the extremes in both the raw and depth-normalized rankings** —
  if depth were the whole story, controlling for it should have reshuffled the ranking
  more than it did (AMR and SAS swap places; the two extremes don't move).
- **Overall correlation, depth vs. CDR3 count: r = 0.492 (r² ≈ 24%)** — depth explains
  about a quarter of the person-to-person variation in recovery. Real, but most of the
  variance (76%) comes from something else — plausibly genuine immune diversity between
  individuals, consistent with the very wide within-group ranges already noted.

**Honest bottom line:** depth is ruled out as a *full* explanation but confirmed as a
*partial* one. The ~1.38x residual gap, with AFR/EAS still at the extremes, is a more
credible signal than the raw 1.63x — but it is not yet proof of reference/pipeline
bias specifically. Unchecked confounds remain: RNA quality (RQS), blood cell
composition, batch/collection site. Treat this as narrowed-down and still open, not
resolved.

## Caveats

- n=16–17 per group; within-group ranges are large (e.g. AFR spans 2,310–17,015 —
  7.4x), so individual variance is substantial relative to the between-group gap. No
  formal significance test run yet.
- Cohort pick is deterministic (sorted by research_id within each ancestry group), not
  randomized — unlikely to matter but not ruled out as a selection effect.
- This is one pipeline (TRUST4) and one library prep (AoU's Watchmaker + Polaris
  Depletion, whole blood) — doesn't generalize to repertoire recovery methods generally.
