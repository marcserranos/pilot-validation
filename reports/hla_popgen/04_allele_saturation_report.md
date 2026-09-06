# Allele-space saturation analysis (`04_allele_saturation.py`)

Methodology: scripts/hla_popgen/research/NOVEL_LIT.md section 3.4. Headline estimator: **Chao2** (incidence-based, haplotypes as sampling units), with ACE and 1st/2nd-order jackknife as nonparametric robustness checks. HLA violates the neutral infinite-alleles model that Ewens/Watterson-style estimators assume (balancing selection -- see module docstring and the neutral-model comparison section below); Chao-family estimators make no assumption about the shape of the allele-frequency spectrum and are used as the headline number for exactly that reason.


## Headline: % of allele space discovered, novel alleles only, classical genes, pooled cohort

| gene | n_haplotypes | S_obs | Chao2 (95% CI) | % discovered (95% CI) | |
|---|---|---|---|---|---|
| HLA-A | 1447 | 718 | 13212.42 (642.76-13212.42) | 5.4% (5.4-100.0%) | LOW-CONFIDENCE (CI widened) |
| HLA-B | 3133 | 929 | 9730.45 (850.28-9730.45) | 9.5% (9.5-100.0%) | LOW-CONFIDENCE (CI widened) |
| HLA-C | 2594 | 949 | 18823.22 (855.11-18823.22) | 5.0% (5.0-100.0%) | LOW-CONFIDENCE (CI widened) |
| HLA-DPA1 | 3132 | 424 | 2329.27 (369.63-2329.27) | 18.2% (18.2-100.0%) | LOW-CONFIDENCE (CI widened) |
| HLA-DPB1 | 5125 | 637 | 4016.34 (569.86-4016.34) | 15.9% (15.9-100.0%) | LOW-CONFIDENCE (CI widened) |
| HLA-DQA1 | 3931 | 463 | 4483.91 (403.54-4483.91) | 10.3% (10.3-100.0%) | LOW-CONFIDENCE (CI widened) |
| HLA-DQB1 | 2839 | 580 | 4106.76 (512.92-4106.76) | 14.1% (14.1-100.0%) | LOW-CONFIDENCE (CI widened) |
| HLA-DRB1 | 14670 | 664 | 3903.78 (590.31-3903.78) | 17.0% (17.0-100.0%) | LOW-CONFIDENCE (CI widened) |

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


## Ancestry discovery-rate, singleton-inclusive (novel alleles, classical genes)
`S_obs` already counts every distinct novel allele seen at least once, singleton or not -- what changes under a singleton-inclusive reading (see 03_novel_alleles.py's `confidence_tier`) is only which of those are called 'reportable'. `singleton_share` here is a different, useful number: of everything a given ancestry HAS found, how much is still only seen once (i.e. how provisional/still-emerging that ancestry's discovered set is) -- and `good_turing_p_next_new` (`Q1/N`) is the probability the very next haplotype from that ancestry reveals something never seen before, in that ancestry.

| ancestry | n_haplotypes (summed, 8 genes) | S_obs (summed) | singleton share | p(next haplotype = new) |
|---|---|---|---|---|
| EUR | 7618 | 1484 | 83.6% | 16.3% |
| MID | 1506 | 363 | 66.9% | 16.1% |
| AMR | 7827 | 1428 | 80.2% | 14.6% |
| AFR | 10521 | 1759 | 83.1% | 13.9% |
| EAS | 5210 | 854 | 78.1% | 12.8% |
| SAS | 4126 | 681 | 76.1% | 12.6% |

Per-gene breakdown (not summed) in `allele_richness.tsv` (`singleton_share`/`good_turing_p_next_new` are derivable there as `Q1_uniques/S_obs` and `Q1_uniques/n_haplotypes` respectively, per gene x ancestry).


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

## Discovery-rate convergence: a second, curve-based richness estimate
Chao2 above uses only Q1/Q2 (the two lowest rarity classes) at the CURRENT sample size. This section instead fits the Clench/Michaelis-Menten saturating curve `S(N) = S_max * N / (b + N)` to the ENTIRE observed rarefaction trajectory (novel alleles, classical genes) and reads off its asymptote `S_max` directly -- a genuinely different estimator (uses the curve's shape, not just two rarity counts), reported alongside Chao2 rather than in place of it, since the two can disagree and that disagreement is itself informative. The Clench form is used specifically because its rate of discovery, `dS/dN = S_max*b/(b+N)^2`, integrates in closed form: the number of alleles still undiscovered at sample size N is EXACTLY the area under that rate curve from N to infinity, `S_max - S(N)` -- i.e. 'integrating the discovery rate to get the total space size' is not just an intuition here, it's the identity this fit is built on. `clench_s_max_bracket_lo/hi` is a curve refit on the rarefaction curve's own 2.5/97.5 percentile bands, NOT a calibrated resample-and-refit bootstrap -- read it as a sensitivity range, not a formal CI.

| gene | N sampled | S_obs now | Clench S_max (bracket) | % discovered (Clench) | additional @ infinite sampling |
|---|---|---|---|---|---|
| HLA-A | 1447 | 718 | 3821 (2798-6482) | 18.8% | 3104 |
| HLA-B | 3133 | 929 | 3665 (2864-5210) | 25.3% | 2736 |
| HLA-C | 2594 | 949 | 7006 (4253-21499) | 13.5% | 6057 |
| HLA-DPA1 | 3132 | 424 | 1798 (1305-4396) | 23.6% | 1374 |
| HLA-DPB1 | 5125 | 637 | 1914 (1545-2647) | 33.3% | 1277 |
| HLA-DQA1 | 3931 | 463 | 2104 (1421-4647) | 22.0% | 1641 |
| HLA-DQB1 | 2839 | 580 | 2563 (1789-4823) | 22.6% | 1983 |
| HLA-DRB1 | 14670 | 664 | 1850 (1511-2439) | 35.9% | 1186 |

Compare this table's `% discovered (Clench)` against the Chao2-based headline table above -- material disagreement between the two estimators on the same gene is a signal worth investigating (e.g. a gene where the rarefaction curve hasn't started bending yet will push Clench's S_max toward implausibly large values; a gene with Q2=0 will make Chao2 unstable instead). Full ancestry-stratified fits: `discovery_rate_fits.tsv`.


Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/discovery_curves_novel.png`

Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/discovery_curves_all.png`

Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/discovery_rate_convergence.png`

Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/frequency_spectrum_novel.png`

Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/frequency_spectrum_known.png`

Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/frequency_spectrum_all.png`
