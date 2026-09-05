# Cross-cohort figures -- sr vs `lr (Immuannot long-read, full)`

SR: 13228 people. LR side: 12233 people.


## SR vs LR frequency outliers (>2x, combined n>=10)

| Allele | sr freq | lr freq | ratio |
|---|---|---|---|
| DQA1*05:01 | 22.932% | 7.657% | 0.33x |
| DQB1*02:01 | 18.896% | 7.610% | 0.40x |
| DRB1*14:01 | 2.094% | 0.194% | 0.09x |
| DQA1*03:02 | 0.004% | 1.955% | 517.34x |
| C*18:01 | 1.104% | 0.391% | 0.35x |
| DPB1*39:01 | 0.212% | 0.030% | 0.14x |
| DPA1*01:05 | 0.087% | 0.004% | 0.05x |
| DPB1*34:01 | 0.072% | 0.022% | 0.30x |

## SR-vs-LR agreement by ancestry (headline reference-bias claim, made quantitative)

Bars for an ancestry with fewer than 10 people in the more limited of the sr/lr cohorts are greyed out and hatched in the figure -- treat those as noise, not signal. Note this uses the STRICTER --min-cell-n-people threshold (10), not --min-cell-n (5): each bar here averages ACROSS an ancestry's whole allele set, so it needs more independent people to be trustworthy than a single raw allele-copy count does (see MIN_CELL_N_PEOPLE's docstring in _viz_common.py).

| Ancestry | mean abs freq diff | n alleles | EUR-common-allele bias (sr-lr) | n EUR-common alleles | n people (lr) | n people (sr) | thin? |
|---|---|---|---|---|---|---|---|
| AFR | 0.616% | 298 | -0.214% | 8 | 3108 | 4052 |  |
| AMR | 0.392% | 316 | -1.389% | 8 | 2777 | 2807 |  |
| EAS | 0.495% | 239 | -0.884% | 8 | 1494 | 1499 |  |
| EUR | 0.427% | 250 | -2.123% | 8 | 3078 | 3109 |  |
| MID | 0.771% | 191 | -3.111% | 8 | 495 | 499 |  |
| SAS | 0.541% | 216 | -0.750% | 8 | 1258 | 1262 |  |

## template_distance: % of haplotype calls at or under each cut, by ancestry

| Ancestry | td<=0 | td<=1 | td<=2 | td<=5 | td<=10 | n |
|---|---|---|---|---|---|---|
| AFR | 78.9% | 89.9% | 93.6% | 96.7% | 98.1% | 48658 |
| AMR | 82.3% | 91.3% | 94.1% | 97.1% | 98.4% | 43466 |
| EAS | 78.1% | 89.9% | 93.4% | 97.1% | 98.4% | 23417 |
| EUR | 84.6% | 92.7% | 95.2% | 97.5% | 98.5% | 48132 |
| MID | 80.8% | 91.7% | 95.1% | 97.7% | 98.5% | 7745 |
| SAS | 79.5% | 92.0% | 95.2% | 97.6% | 98.7% | 19758 |

## Resolution cascade: H(4-field | 2-field), nats

Grey/hatched bars in the figure: fewer than 10 people for that ancestry (`n_people` below; `n` is haplotype-call ROWS, which can be much larger even for a thin ancestry -- do not read `n` alone as statistical power).

| Gene | Ancestry | n (rows) | n_people | n 2-field groups | n 4-field subtypes | H(4|2) | thin? |
|---|---|---|---|---|---|---|---|
| A | AFR | 5580 | 3108 | 41 | 93 | 0.341 |  |
| A | AMR | 4988 | 2777 | 46 | 114 | 0.266 |  |
| A | EAS | 2694 | 1494 | 39 | 66 | 0.139 |  |
| A | EUR | 5624 | 3078 | 45 | 93 | 0.203 |  |
| A | MID | 907 | 495 | 35 | 59 | 0.214 |  |
| A | SAS | 2299 | 1258 | 33 | 54 | 0.242 |  |
| B | AFR | 5059 | 3108 | 71 | 189 | 0.548 |  |
| B | AMR | 4564 | 2777 | 89 | 212 | 0.603 |  |
| B | EAS | 2396 | 1494 | 69 | 134 | 0.390 |  |
| B | EUR | 5317 | 3078 | 73 | 196 | 0.458 |  |
| B | MID | 823 | 495 | 53 | 99 | 0.493 |  |
| B | SAS | 2101 | 1258 | 56 | 117 | 0.438 |  |
| C | AFR | 5317 | 3108 | 33 | 123 | 0.737 |  |
| C | AMR | 4838 | 2777 | 41 | 136 | 0.845 |  |
| C | EAS | 2484 | 1494 | 31 | 97 | 0.620 |  |
| C | EUR | 5426 | 3078 | 35 | 128 | 0.744 |  |
| C | MID | 862 | 495 | 31 | 72 | 0.669 |  |
| C | SAS | 2174 | 1258 | 33 | 90 | 0.591 |  |
| DPA1 | AFR | 5084 | 3108 | 9 | 69 | 1.678 |  |
| DPA1 | AMR | 4643 | 2777 | 9 | 69 | 2.002 |  |
| DPA1 | EAS | 2468 | 1494 | 6 | 41 | 0.937 |  |
| DPA1 | EUR | 5286 | 3078 | 10 | 62 | 1.879 |  |
| DPA1 | MID | 842 | 495 | 7 | 41 | 1.901 |  |
| DPA1 | SAS | 2071 | 1258 | 6 | 39 | 1.623 |  |
| DPB1 | AFR | 4290 | 3108 | 34 | 138 | 1.275 |  |
| DPB1 | AMR | 4363 | 2777 | 33 | 154 | 1.066 |  |
| DPB1 | EAS | 1951 | 1494 | 29 | 105 | 1.156 |  |
| DPB1 | EUR | 5094 | 3078 | 30 | 171 | 1.429 |  |
| DPB1 | MID | 748 | 495 | 21 | 81 | 1.446 |  |
| DPB1 | SAS | 1840 | 1258 | 23 | 88 | 1.145 |  |
| DQA1 | AFR | 5002 | 3108 | 16 | 91 | 1.152 |  |
| DQA1 | AMR | 4728 | 2777 | 15 | 91 | 1.042 |  |
| DQA1 | EAS | 2405 | 1494 | 16 | 72 | 0.814 |  |
| DQA1 | EUR | 5129 | 3078 | 17 | 87 | 0.982 |  |
| DQA1 | MID | 795 | 495 | 15 | 60 | 1.130 |  |
| DQA1 | SAS | 1963 | 1258 | 15 | 63 | 0.671 |  |
| DQB1 | AFR | 5168 | 3108 | 19 | 73 | 0.697 |  |
| DQB1 | AMR | 4794 | 2777 | 19 | 80 | 0.775 |  |
| DQB1 | EAS | 2571 | 1494 | 16 | 55 | 0.577 |  |
| DQB1 | EUR | 5349 | 3078 | 18 | 83 | 0.646 |  |
| DQB1 | MID | 831 | 495 | 19 | 47 | 0.725 |  |
| DQB1 | SAS | 2197 | 1258 | 16 | 56 | 0.587 |  |
| DRB1 | AFR | 1944 | 3108 | 33 | 81 | 0.513 |  |
| DRB1 | AMR | 2093 | 2777 | 34 | 85 | 0.647 |  |
| DRB1 | EAS | 786 | 1494 | 24 | 50 | 0.509 |  |
| DRB1 | EUR | 2821 | 3078 | 30 | 78 | 0.687 |  |
| DRB1 | MID | 331 | 495 | 20 | 42 | 0.666 |  |
| DRB1 | SAS | 761 | 1258 | 21 | 49 | 0.461 |  |

## Non-classical gene heterozygosity by ancestry (LR-only; sr has no calls for these genes)

| Gene | Ancestry | n | Heterozygosity |
|---|---|---|---|
| MICA | AFR | 5324 | 0.857 |
| MICA | AMR | 4827 | 0.886 |
| MICA | EAS | 2637 | 0.880 |
| MICA | EUR | 5536 | 0.899 |
| MICA | MID | 892 | 0.905 |
| MICA | SAS | 2298 | 0.901 |
| MICB | AFR | 5495 | 0.743 |
| MICB | AMR | 4953 | 0.679 |
| MICB | EAS | 2706 | 0.729 |
| MICB | EUR | 5482 | 0.784 |
| MICB | MID | 898 | 0.780 |
| MICB | SAS | 2277 | 0.794 |
| TAP1 | AFR | 4750 | 0.380 |
| TAP1 | AMR | 4732 | 0.292 |
| TAP1 | EAS | 2686 | 0.338 |
| TAP1 | EUR | 5419 | 0.273 |
| TAP1 | MID | 876 | 0.332 |
| TAP1 | SAS | 2192 | 0.343 |
| TAP2 | AFR | 5096 | 0.673 |
| TAP2 | AMR | 4896 | 0.714 |
| TAP2 | EAS | 2474 | 0.736 |
| TAP2 | EUR | 5644 | 0.688 |
| TAP2 | MID | 882 | 0.685 |
| TAP2 | SAS | 2196 | 0.663 |
| DRB3 | AFR | 2914 | 0.637 |
| DRB3 | AMR | 2000 | 0.612 |
| DRB3 | EAS | 1077 | 0.588 |
| DRB3 | EUR | 2462 | 0.562 |
| DRB3 | MID | 480 | 0.506 |
| DRB3 | SAS | 941 | 0.451 |
| DRB4 | AFR | 1143 | 0.486 |
| DRB4 | AMR | 1645 | 0.351 |
| DRB4 | EAS | 866 | 0.050 |
| DRB4 | EUR | 1770 | 0.330 |
| DRB4 | MID | 250 | 0.207 |
| DRB4 | SAS | 672 | 0.044 |
| DRB5 | AFR | 920 | 0.225 |
| DRB5 | AMR | 610 | 0.510 |
| DRB5 | EAS | 579 | 0.478 |
| DRB5 | EUR | 745 | 0.444 |
| DRB5 | MID | 91 | 0.593 |
| DRB5 | SAS | 576 | 0.554 |
| DQA2 | AFR | 6012 | 0.563 |
| DQA2 | AMR | 5380 | 0.363 |
| DQA2 | EAS | 2909 | 0.210 |
| DQA2 | EUR | 5981 | 0.302 |
| DQA2 | MID | 970 | 0.387 |
| DQA2 | SAS | 2458 | 0.370 |
| DQB2 | AFR | 5823 | 0.511 |
| DQB2 | AMR | 5304 | 0.543 |
| DQB2 | EAS | 2762 | 0.463 |
| DQB2 | EUR | 5884 | 0.587 |
| DQB2 | MID | 916 | 0.508 |
| DQB2 | SAS | 2284 | 0.651 |
| E | AFR | 5813 | 0.509 |
| E | AMR | 5196 | 0.504 |
| E | EAS | 2795 | 0.518 |
| E | EUR | 5751 | 0.508 |
| E | MID | 925 | 0.502 |
| E | SAS | 2367 | 0.506 |
| F | AFR | 5889 | 0.328 |
| F | AMR | 5242 | 0.223 |
| F | EAS | 2809 | 0.069 |
| F | EUR | 5792 | 0.277 |
| F | MID | 936 | 0.257 |
| F | SAS | 2377 | 0.161 |
| G | AFR | 5898 | 0.579 |
| G | AMR | 5297 | 0.517 |
| G | EAS | 2850 | 0.499 |
| G | EUR | 5886 | 0.426 |
| G | MID | 944 | 0.599 |
| G | SAS | 2418 | 0.563 |

## Cis heterodimer frequency + pairable fraction by ancestry


### DQA1~DQB1

Phased cis-haplotype copies observed by ancestry: {'AFR': 5805, 'AMR': 5249, 'EAS': 2884, 'EUR': 5818, 'MID': 946, 'SAS': 2408}

| Ancestry | pairable fraction | n hap-instances with both genes present |
|---|---|---|
| AFR | 97.4% | 5859 |
| AMR | 98.3% | 5251 |
| EAS | 98.8% | 2864 |
| EUR | 98.3% | 5804 |
| MID | 98.7% | 945 |
| SAS | 98.2% | 2429 |

### DPA1~DPB1

Phased cis-haplotype copies observed by ancestry: {'AFR': 5993, 'AMR': 5337, 'EAS': 2880, 'EUR': 5909, 'MID': 962, 'SAS': 2425}

| Ancestry | pairable fraction | n hap-instances with both genes present |
|---|---|---|
| AFR | 99.8% | 5931 |
| AMR | 99.8% | 5268 |
| EAS | 99.7% | 2855 |
| EUR | 99.9% | 5857 |
| MID | 99.9% | 957 |
| SAS | 99.8% | 2407 |

## Figures

- `/home/jupyter/repos/pilot-validation/reports/hla_popgen/07_figures_crosscohort/lr/sr_vs_lr_frequency_scatter.png`
- `/home/jupyter/repos/pilot-validation/reports/hla_popgen/07_figures_crosscohort/lr/sr_lr_disagreement_by_ancestry.png`
- `/home/jupyter/repos/pilot-validation/reports/hla_popgen/07_figures_crosscohort/lr/template_distance_by_ancestry.png`
- `/home/jupyter/repos/pilot-validation/reports/hla_popgen/07_figures_crosscohort/lr/template_distance_sweep.png`
- `/home/jupyter/repos/pilot-validation/reports/hla_popgen/07_figures_crosscohort/lr/resolution_cascade.png`
- `/home/jupyter/repos/pilot-validation/reports/hla_popgen/07_figures_crosscohort/lr/nonclassical_diversity.png`
- `/home/jupyter/repos/pilot-validation/reports/hla_popgen/07_figures_crosscohort/lr/gene_coverage_by_cohort.png`
- `/home/jupyter/repos/pilot-validation/reports/hla_popgen/07_figures_crosscohort/lr/heterodimer_DQA1_DQB1.png`
- `/home/jupyter/repos/pilot-validation/reports/hla_popgen/07_figures_crosscohort/lr/heterodimer_DPA1_DPB1.png`