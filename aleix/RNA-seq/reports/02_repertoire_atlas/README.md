# 02 — Population-scale TRB repertoire atlas

Script: `scripts/02_repertoire_atlas.py`. Tables from the 2026-09-29 run (d1bacc8); figure
re-rendered in the journal font on 2026-10-13 (b47314e), numbers unchanged.

## Question

What do 7,922 whole-blood TRB repertoires recovered from bulk RNA-seq look like? And do they
carry real immunology? The test is a textbook one: TCR diversity falls and clonal expansion
rises with age (Britanova et al. 2014 J Immunol; Qi et al. 2014 PNAS).

## How

- **Clonotypes:** the shared definition from `embed_cdr3s.py` (`cdr3.out`, CDR3_score ≥ 0.02,
  canonical junction, ≤ 30 aa, TRB, unique TRBV + CDR3 aa per person, reads summed). This uses
  the full repertoire, not the top-500 embedding sample.
- **Diversity:** rarefied richness (Hurlbert), i.e. the expected number of distinct clonotypes
  in 823 reads. 823 is the depth the analysis table (09) later uses for every fixed-depth
  measure. Rarefaction removes sequencing depth and blood T-cell content as explanations.
- **Expansion:** share of TRB reads in the person's 10 largest clonotypes.
- **Age models:** OLS on age per decade, adjusted for sex and ancestry; 7,105 people.
  Sensitivity model adds log total TRB reads.
- **Age:** 2020 − year of birth (draws span ~2017–2023, so ± ~3 years).

## Results

- 9.44 M TRB clonotypes; median 1,136 per person and 2,076 TRB reads per person.
- **Diversity falls 3.3% per decade** (95% CI 3.1–3.5; p 2e-168). With log reads in the model,
  3.2%. Women have 5% higher richness.
- **Expansion rises 1.91 points per decade** (1.79–2.04; p 4e-186): median top-10 share 10.2% at
  18–29, 18.5% at 70+. With log reads, 1.93. Women have 2.0 points lower expansion.
- **Sharing:** 7.2% of distinct clonotypes are carried by ≥ 2 people in the full repertoire.
  Report 01's 14.3% is on the top-500 sample, which keeps the most expanded clones.
- Medians by ancestry: 1,078 (AMR) to 1,270 (SAS) clonotypes per person. Groups also differ in
  median age (40 SAS to 60 EUR), so ancestry comparisons here are not age-matched.

## Figure

`fig_repertoire_atlas`:
- a: cohort size by ancestry.
- b: TRB clonotypes per person.
- c: CDR3 length by ancestry.
- d: TRBV usage across people (median, IQR, 5–95%; 24 most used genes).
- e: sharing spectrum.
- f: rarefied richness vs age.
- g: top-10 share vs age.

Bins with fewer than 20 people are not drawn (120 people sit in suppressed histogram bins).

**Pending supervisor changes (DECISIONS.md, week of 2026-10-12):**
- drop c;
- f becomes Shannon entropy with total TRB reads regressed out, vs age;
- g to be fixed in the same spirit.

These need the per-person table (VM). From `age_models.csv`, adjusting for reads barely moves
either age slope, so the age trends should hold.

## Caveats

- Cross-sectional; age effects are between people, not within.
- The 12 GB BAM cap drops the deepest 4.9% of samples (report 08). Rarefaction protects f.
- No blood cell counts: T-cell fraction of the sample is not adjusted.

## Files

`summary.csv`, `age_models.csv` (all coefficients, both models), `age_bins.csv` (f/g
quantiles), `cohort_by_ancestry.csv`, `cdr3_length_by_ancestry.csv`, `trbv_usage.csv`,
`sharing_spectrum.csv`.
