# RNA-seq repertoire recovery across ancestry — 100-person pilot — 2026-08-11

**Headline, corrected 2026-08-11: a visually striking ancestry pattern that does NOT
reach statistical significance at this sample size.** AFR shows the highest mean CDR3
recovery (7,406), EAS the lowest (4,538) — a 1.63x spread, partly (not fully) explained
by sequencing depth. But a formal Kruskal-Wallis test across all 6 groups gives
**p = 0.18** — well above the conventional 0.05 threshold. The earlier version of this
document argued the 5-of-5-chain-types consistency made this "unlikely to be chance" —
that was an intuitive pattern-match, not a statistical test, and the actual test does
not support that claim as strongly. Treat this as suggestive and worth a larger cohort,
not as a confirmed finding.

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

## The pattern that looked compelling — and why it isn't proof

Chain-level breakdown (TRA/TRB/IGH/IGK/IGL — the 5 chains with substantial counts;
TRG/TRD are too low-count to read much into): **AFR is highest, and EAS is lowest, in
every one of the 5 major chain types.** That consistency is genuinely part of why this
looked like more than noise on first read.

**But it isn't a substitute for a real significance test, and the real test doesn't
back it up as strongly** — see "Significance test" below. The 5 chain types aren't
independent measurements the way the informal argument implicitly treated them: they
all come from the same 100 people, driven substantially by the same per-person depth
and immune-diversity variation already documented in "Depth-confound check." Five
correlated measurements agreeing is much weaker evidence than five independent ones
would be. Keeping this section rather than deleting it, because the correction itself
is the useful record here.

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

**Bottom line at this stage:** depth is ruled out as a *full* explanation but confirmed
as a *partial* one. The ~1.38x residual gap, with AFR/EAS still at the extremes, was a
more credible signal than the raw 1.63x — but see the significance test below before
weighting this too heavily.

## Remaining confounds + significance test — run 2026-08-11, `../scripts/check_ancestry_confounds.py`

Two more mundane explanations, checked and **ruled out**:

| | Correlation with CDR3 | Range across ancestry |
|---|---|---|
| RQS (RNA quality, 0–10 scale) | r = 0.061 (none) | 7.82–8.00 (flat) |
| T-cell fraction (chain-based proxy) | r = −0.409 (moderate, real) | 0.409–0.446 (flat) |

RQS doesn't predict recovery at all. T-cell fraction does correlate moderately with
CDR3 count (plausibly: more B-cell-heavy samples recover more distinct sequences,
since somatic hypermutation keeps generating new B-cell variants) — but since neither
metric varies meaningfully by ancestry, neither can be driving the ancestry gap.

**The significance test: Kruskal-Wallis across all 6 ancestry groups, raw CDR3 counts.
H = 7.62, p = 0.18.** Above the conventional 0.05 threshold — **the between-ancestry
differences are not statistically distinguishable from noise at this sample size.**
This is the single most important correction in this document: the "5-of-5 chain
types" pattern described above is real in the data but does not clear the bar for
statistical significance once individual variance is properly accounted for.

**What this means:** not "there is no ancestry effect," but "we cannot yet tell,
given n=16–17 per group and this much within-group spread." The direct next step —
the planned second 100-person cohort on a different machine — is now motivated by
statistical necessity (roughly doubling n per group to ~33–34), not just extra data
for its own sake.

## Caveats

- n=16–17 per group in this analysis; a properly powered answer needs more people per
  group, which is exactly what the second cohort is for.
- Cohort pick is deterministic (sorted by research_id within each ancestry group), not
  randomized — unlikely to matter but not ruled out as a selection effect.
- Depth, RQS, and T-cell fraction are checked; batch/collection site is not.
- This is one pipeline (TRUST4) and one library prep (AoU's Watchmaker + Polaris
  Depletion, whole blood) — doesn't generalize to repertoire recovery methods generally.
