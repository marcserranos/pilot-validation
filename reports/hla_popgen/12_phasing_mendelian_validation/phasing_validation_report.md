# Phasing validation via real relative pairs -- Mendelian consistency + switch scan

Resolution: 2-field. See module docstring for full method.

## Pair classification

| class | n_pairs |
|---|---|
| high_sharing | 545 |
| mixed | 29 |

## Headline error rate (high-sharing pairs only)

Overall, across all genes: **4.865%** (95% CI 4.566-5.183%), 908/18664 gene-comparisons.

### By gene class

| gene_class | n_compared | n_mismatch | error_rate_% | 95% CI |
|---|---|---|---|---|
| class_II_accessory | 2712 | 27 | 0.996 | 0.685-1.445 |
| class_II_paralog | 2498 | 99 | 3.963 | 3.266-4.802 |
| classical_I | 1617 | 109 | 6.741 | 5.618-8.068 |
| classical_II | 2704 | 135 | 4.993 | 4.234-5.879 |
| mic_tap | 2097 | 172 | 8.202 | 7.103-9.454 |
| nonclassical_I | 1622 | 41 | 2.528 | 1.869-3.411 |
| other | 49 | 7 | 14.286 | 7.096-26.668 |
| pseudogene_I | 5365 | 318 | 5.927 | 5.326-6.591 |

**Classical HLA genes (A/B/C/DPA1/DPB1/DQA1/DQB1/DRB1) specifically: 5.647% (95% CI 4.997-6.375%)**

### Worst individual genes (by point estimate, n>=10 comparisons)

| gene | gene_class | n_compared | n_mismatch | error_rate_% |
|---|---|---|---|---|
| HLA-HFE | other | 49 | 7 | 14.286 |
| MICA | mic_tap | 533 | 73 | 13.696 |
| HLA-H | pseudogene_I | 520 | 69 | 13.269 |
| HLA-K | pseudogene_I | 343 | 45 | 13.120 |
| HLA-U | pseudogene_I | 427 | 53 | 12.412 |
| HLA-DRB3 | class_II_paralog | 256 | 26 | 10.156 |
| HLA-W | pseudogene_I | 322 | 32 | 9.938 |
| MICB | mic_tap | 520 | 44 | 8.462 |
| HLA-P | pseudogene_I | 488 | 41 | 8.402 |
| TAP2 | mic_tap | 535 | 43 | 8.037 |

## Phase-switch scan

545/545 high-sharing pairs (100.0%) show ZERO linkage-phase switches across the region -- one unbroken shared haplotype block, the expected result given the MHC's known low internal recombination.

## Caveats

- Cannot distinguish parent-child from full-sibling pairs (no IBD0 data) -- pair classification is by empirical sharing rate, not known pedigree.

- copy_index>1 (segmental duplications, DRB paralog CNV) genes excluded.

- Genes below the minimum-comparison-count floor are excluded from the gene-level table/figure but included in the overall headline rate.

- The switch-scan only counts a transition between two genes on the SAME contig for BOTH people -- comparing across a contig boundary would measure arbitrary hap-label reshuffling between independently assembled fragments, not a real crossover or assembly error.

## Run provenance

Real VM run, 2026-09-08, against the 574 first-degree-or-twin `both_lr` relative pairs identified by `11_relatedness_cohort_overlap.py` (827-person pool). 545/574 classified `high_sharing` (empirical per-gene allele-sharing >= 0.8 -- consistent with a real parent-child pair or an IBD1/IBD2 sibling pair); 29 classified `mixed` (intermediate sharing, illustrated in the `mixed_example_*.png` figures, anonymized to "Example A/B/..." rather than real person_ids); none classified `low_sharing` in this specific first-degree/twin pair set. 42 genes received a derived canonical physical order from the cohort's own assemblies (median relative rank of `gene_start` within single-contig hap instances).

**Bug found and fixed during this run:** the first version of the switch-scan (no contig awareness) reported 518/545 pairs with >=1 "switch" -- implausibly high for real meiotic recombination in a ~4Mb window. Root cause: a hap1/hap2 label is only guaranteed phase-consistent *within* one assembled contig, not across contigs, when a person's own HLA-region assembly is fragmented (a known, common occurrence per this project's own SCHEMA.md). The fix restricts switch-counting to consecutive genes sharing the same contig for both people; the corrected result (0/545 switches) matches the literature-expected result cleanly. See `12_phasing_mendelian_validation.py`'s `count_switches` docstring for the full explanation.
