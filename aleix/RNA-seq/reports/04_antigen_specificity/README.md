# 04 — Known antigen specificities: benchmark and cohort matching

Script: `scripts/04_antigen_specificity.py`. Current outputs from the 2026-10-04 polish run
(5b83486), with epitope-bootstrap CIs.

## Questions

1. Does SCEPTR place TCRs that recognise the same epitope closer together? (public data)
2. Can exact matching to VDJdb detect past exposure (e.g. CMV) in our repertoires?

## How

- **Benchmark:** 4,193 VDJdb human TRB TCRs over 45 epitopes.
  - Metrics: pairwise AUROC (same vs different epitope, score = −distance) and 5-NN balanced
    accuracy (chance 0.022).
  - Methods compared: SCEPTR on TRBV + CDR3 (production), SCEPTR CDR3 only, CDR3 3-mers, TRBV
    alone.
  - CIs: 500 epitope-bootstrap replicates, with paired differences.
- **Cohort:**
  - Exact matches (same TRBV gene, identical CDR3 aa) of each person's full repertoire to
    VDJdb.
  - Null: OLGA decoys matched per reference TCR on TRBV, CDR3 length and Pgen. 249,616 decoys,
    3 replicates; 97.9% within 0.25 log10 Pgen.
  - Enrichment = observed / decoy-expected carriage.
  - Negative controls: pathogens with US prevalence ≤ ~1% (HBV, HCV, HIV-1, DENV) and human
    self-peptides.

## Results

- **Benchmark (positive):**
  - Pairwise AUROC: SCEPTR TRBV + CDR3 0.580 (0.527–0.641), CDR3 only 0.552, 3-mer 0.554,
    TRBV only 0.553.
  - 5-NN balanced accuracy: 0.380, 0.328, 0.267 and 0.106.
  - CDR3-only and TRBV-only are significantly lower (paired CIs); 3-mer is lower on 5-NN only.
  - SCEPTR is ahead for 76% of epitopes.
- **Cohort (negative, quantified):** naive carriage does not track seroprevalence, and every
  pathogen is enriched over the decoy null, negative controls included:
  - EBV 7.6×, MCPyV 6.2×, CMV 5.5×, HCV 5.1×, HBV 4.4×, influenza A 3.9×, DENV 3.4×;
  - self-peptides 2.9×, HIV-1 2.5×.
  - Excess CMV matches vs age: +0.007 per decade, p 0.12. CMV seroprevalence rises with age,
    so a working assay should show a clear rise.
  - Reading: at ~2,000 TRB reads per person, database matching cannot separate exposure from
    chance recombination. The decoy null controls recombination but not thymic selection.
  - The ≤ 1 aa variant matches most people for every pathogen. It is a failure mode, not a
    result.
- **HLA bridge (panel f):** the most carried CMV epitope is NLVPMVATV (HLA-A\*02:01).
  Carriage depends on HLA, so no ancestry difference in carriage is interpretable without HLA.

## Figure

`fig_antigen_specificity`:
- a: benchmark AUROC.
- b: per epitope, SCEPTR vs 3-mer.
- c: carriage vs seroprevalence.
- d: enrichment over the null.
- e: excess matches vs age.
- f: CMV epitopes carried, with the restricting HLA allele.

Known layout issues: c's x tick labels collide at 0.1–0.5, and f's x-axis label is clipped.

## Caveats

SCEPTR's training repertoires may overlap VDJdb, so absolute benchmark values are optimistic
for every SCEPTR variant equally. The ranking is the result.

## Files

`summary.csv`, `benchmark_methods.csv`, `benchmark_per_epitope.csv`, `carriage_by_pathogen.csv`,
`carriage_by_epitope.csv`, `carriage_by_age.csv`, `carriage_models.csv`, `excess_by_age.csv`,
`excess_vs_age_models.csv`, `decoy_tolerance.csv`, `cmv_matches_by_hla.csv`.
