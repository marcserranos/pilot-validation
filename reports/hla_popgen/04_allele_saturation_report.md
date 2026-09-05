# Allele-space saturation analysis (`04_allele_saturation.py`)

Methodology: scripts/hla_popgen/research/NOVEL_LIT.md section 3.4. Headline estimator: **Chao2** (incidence-based, haplotypes as sampling units), with ACE and 1st/2nd-order jackknife as nonparametric robustness checks. HLA violates the neutral infinite-alleles model that Ewens/Watterson-style estimators assume (balancing selection -- see module docstring and the neutral-model comparison section below); Chao-family estimators make no assumption about the shape of the allele-frequency spectrum and are used as the headline number for exactly that reason.


## Headline: % of allele space discovered, novel alleles only, classical genes, pooled cohort

| gene | n_haplotypes | S_obs | Chao2 (95% CI) | % discovered (95% CI) | |
|---|---|---|---|---|---|
| HLA-A | 1447 | 718 | 13212.42 (646.63-13212.42) | 5.4% (5.4-100.0%) | LOW-CONFIDENCE (CI widened) |
| HLA-B | 3133 | 929 | 9730.45 (836.98-9730.45) | 9.5% (9.5-100.0%) | LOW-CONFIDENCE (CI widened) |
| HLA-C | 2594 | 949 | 18823.22 (865.98-18823.22) | 5.0% (5.0-100.0%) | LOW-CONFIDENCE (CI widened) |
| HLA-DPA1 | 3132 | 424 | 2329.27 (365.2-2329.27) | 18.2% (18.2-100.0%) | LOW-CONFIDENCE (CI widened) |
| HLA-DPB1 | 5125 | 637 | 4016.34 (563.52-4016.34) | 15.9% (15.9-100.0%) | LOW-CONFIDENCE (CI widened) |
| HLA-DQA1 | 3931 | 463 | 4483.91 (408.36-4483.91) | 10.3% (10.3-100.0%) | LOW-CONFIDENCE (CI widened) |
| HLA-DQB1 | 2839 | 580 | 4106.76 (515.59-4106.76) | 14.1% (14.1-100.0%) | LOW-CONFIDENCE (CI widened) |
| HLA-DRB1 | 14670 | 664 | 3903.78 (589.52-3903.78) | 17.0% (17.0-100.0%) | LOW-CONFIDENCE (CI widened) |

**LOW-CONFIDENCE rows** carry one or both of: `Q2=0` (zero doubletons -- Chao2's bias-corrected form is degenerate in this regime, growing quadratically in the singleton count with nothing to anchor it; ACE/jackknife columns in the full TSV are more conservative fallbacks) and `CI widened` (the raw bootstrap percentile interval did not bracket the point estimate on the full sample and was widened to guarantee it does -- see `04_allele_saturation.py`'s `analyze()` for why this happens and is expected precisely in the Q2=0 regime). Neither should be read as a reliable point estimate on its own. This mostly hits small, low-haplotype-count (gene, ancestry) strata even on the current recurrence-bearing fixtures -- the pooled/`known`-category rows, which have larger N and a finite underlying allele pool, saturate cleanly and rarely carry this flag.


## Ancestry-stratified % discovered, novel alleles only, classical genes

| gene | AFR | AMR | EAS | EUR | MID | SAS |
|---|---|---|---|---|---|---|
| HLA-A | 6.1% | 7.8% | 8.8% | 7.8% | 22.0% | 16.5% |
| HLA-B | 8.8% | 18.7% | 24.2% | 11.4% | 24.2% | 15.3% |
| HLA-C | 5.7% | 11.4% | 8.1% | 2.9% | 17.9% | 6.8% |
| HLA-DPA1 | 17.6% | 2.5% | 18.4% | 11.7% | 19.2% | 4.6% |
| HLA-DPB1 | 20.3% | 26.1% | 11.5% | 19.6% | 31.8% | 19.2% |
| HLA-DQA1 | 5.9% | 14.1% | 5.3% | 8.6% | 32.0% | 14.1% |
| HLA-DQB1 | 17.2% | 10.7% | 11.2% | 5.1% | 27.4% | 15.9% |
| HLA-DRB1 | 14.0% | 10.9% | 27.3% | 24.3% | 24.1% | 24.6% |

IPD-IMGT/HLA's European bias predicts EUR should show the highest % discovered (closest to saturation) and AFR the lowest -- this table is where that claim either holds up or doesn't against the data actually loaded.


## Extrapolation: projected additional novel alleles at 2x/5x/10x cohort size (pooled, classical genes)

| gene | 2x | 5x | 10x |
|---|---|---|---|
| HLA-A | 635.3 | 2353.8 | 4682.6 |
| HLA-B | 781.1 | 2732.6 | 4988.2 |
| HLA-C | 865.3 | 3218.0 | 6438.5 |
| HLA-DPA1 | 328.0 | 1010.5 | 1557.4 |
| HLA-DPB1 | 482.0 | 1553.2 | 2533.3 |
| HLA-DQA1 | 390.7 | 1349.4 | 2418.4 |
| HLA-DQB1 | 469.6 | 1535.5 | 2552.2 |
| HLA-DRB1 | 497.4 | 1576.5 | 2516.9 |

## Neutral (Ewens/Watterson) comparison -- SECONDARY, quantifying the selection effect
Watterson's theta fit to the observed distinct-allele count, then the neutral model's own expectation of additional alleles at 2x the current haplotype count -- compared against Chao2's nonparametric estimate of currently-unseen alleles. **Do not read the neutral column as a discovery projection** -- it is here specifically to show how much a neutral-model fit would have missed, per NOVEL_LIT.md section 4's instruction to turn the assumption violation into a result.

| gene | S_obs | Watterson theta | neutral additional @2x | Chao2 unseen (now) |
|---|---|---|---|---|
| HLA-A | 1015 | 215.718 | 148.55 | 6372.9 |
| HLA-B | 1431 | 334.995 | 229.85 | 4890.0 |
| HLA-C | 1303 | 296.787 | 203.87 | 9414.5 |
| HLA-DPA1 | 553 | 101.203 | 69.93 | 1737.2 |
| HLA-DPB1 | 942 | 196.505 | 135.39 | 2488.3 |
| HLA-DQA1 | 603 | 112.685 | 77.84 | 2762.4 |
| HLA-DQB1 | 750 | 147.679 | 101.9 | 2479.4 |
| HLA-DRB1 | 814 | 163.994 | 113.1 | 2755.1 |

Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/discovery_curves_novel.png`

Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/discovery_curves_all.png`
