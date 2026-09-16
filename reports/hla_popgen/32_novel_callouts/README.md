# 32 — main-text callout list: common novel HLA alleles absent from IPD-IMGT/HLA

*Supervisor ask A6 (Cole, 2026-09-17): "a really great list of those, because those are things we can literally call out in main text".*

Built from committed aggregates only (scripts 24, 25, 27). No VM, no participant data.


## Tier 1 — nameable in main text (38 alleles)

Clean (no artifact flag, not homopolymer-only), recurrent, and carried by at least 20 unrelated people — which is also exactly the threshold above which we are permitted to publish the count.

### 1a. Novel **proteins** (29)

| gene     | nearest_known_protein   |   n_aa_diffs | in_groove   |   n_persons_unrelated_clean | top_ancestry   |   top_ancestry_share |   carriers_per_1000_in_top_ancestry |   pct_gene_haplotypes_not_in_imgt | ancestries_disclosable      |
|:---------|:------------------------|-------------:|:------------|----------------------------:|:---------------|---------------------:|------------------------------------:|----------------------------------:|:----------------------------|
| TAP1     | TAP1*01:01              |            1 | nan         |                         426 | AFR            |                0.885 |                              139.35 |                             7.981 | AFR:292,AMR:38              |
| TAP1     | TAP1*05:01              |            1 | nan         |                         321 | SAS            |                0.354 |                               64.67 |                             7.981 | AFR:26,AMR:33,EUR:67,SAS:69 |
| TAP2     | TAP2*02:01              |            1 | nan         |                         310 | EAS            |                0.648 |                              160.48 |                             6.847 | EAS:188,SAS:102             |
| TAP2     | TAP2*01:01              |            2 | NA,NA       |                         201 | AFR            |                0.832 |                               55.95 |                             6.847 | AFR:119,AMR:24              |
| TAP2     | TAP2*01:01              |            1 | nan         |                         192 | AFR            |                1     |                               63.94 |                             6.847 | AFR:136                     |
| TAP1     | TAP1*04:01              |            1 | nan         |                         191 | AFR            |                0.693 |                               33.4  |                             7.981 | AFR:70,AMR:31               |
| TAP1     | TAP1*01:01              |            1 | nan         |                         189 | AFR            |                1     |                               65.38 |                             7.981 | AFR:137                     |
| TAP1     | TAP1*05:01              |            2 | NA,NA       |                          98 | AFR            |                1     |                               37.7  |                             7.981 | AFR:79                      |
| TAP2     | TAP2*01:01              |            1 | nan         |                          94 | AFR            |                1     |                               31.97 |                             6.847 | AFR:68                      |
| TAP1     | TAP1*02:01              |            2 | NA,NA       |                          94 | SAS            |                1     |                               53.42 |                             7.981 | SAS:57                      |
| MICB     | MICB*001                |            1 | nan         |                          81 | AFR            |                1     |                               25.96 |                             0.847 | AFR:53                      |
| HLA-DQB2 | HLA-DQB2*01:01          |            1 | False       |                          80 | AFR            |                1     |                               22.13 |                             1.263 | AFR:49                      |
| HLA-DQB2 | HLA-DQB2*01:08          |            1 | False       |                          62 | AFR            |                1     |                               21.22 |                             1.263 | AFR:47                      |
| TAP2     | TAP2*02:01              |            1 | nan         |                          56 | EAS            |                1     |                               35    |                             6.847 | EAS:41                      |
| TAP2     | TAP2*01:01              |            1 | nan         |                          53 | AFR            |                1     |                               19.75 |                             6.847 | AFR:42                      |
| TAP2     | TAP2*02:01              |            2 | NA,NA       |                          50 | AMR            |                1     |                               19.78 |                             6.847 | AMR:29                      |
| TAP1     | TAP1*01:01              |            1 | nan         |                          42 | AFR            |                1     |                               10.5  |                             7.981 | AFR:22                      |
| TAP2     | TAP2*01:02              |            1 | nan         |                          40 | AFR            |                1     |                               12.69 |                             6.847 | AFR:27                      |
| TAP2     | TAP2*01:01              |            1 | nan         |                          36 | AFR            |                1     |                                9.4  |                             6.847 | AFR:20                      |
| HLA-DMA  | HLA-DMA*01:01           |            1 | nan         |                          36 | nan            |              nan     |                              nan    |                             0.403 | none                        |
| TAP1     | TAP1*02:01              |            1 | nan         |                          31 | nan            |              nan     |                              nan    |                             7.981 | none                        |
| TAP2     | TAP2*01:01              |            1 | nan         |                          28 | nan            |              nan     |                              nan    |                             6.847 | none                        |
| TAP1     | TAP1*01:01              |            1 | nan         |                          28 | nan            |              nan     |                              nan    |                             7.981 | none                        |
| HLA-F    | HLA-F*01:01             |            1 | False       |                          26 | nan            |              nan     |                              nan    |                             0.543 | none                        |
| HLA-DQB2 | HLA-DQB2*01:02          |            1 | False       |                          23 | nan            |              nan     |                              nan    |                             1.263 | none                        |
| HLA-DQB2 | HLA-DQB2*01:02          |            1 | False       |                          23 | AFR            |                1     |                                9.03 |                             1.263 | AFR:20                      |
| HLA-DMB  | HLA-DMB*01:07           |            1 | nan         |                          22 | nan            |              nan     |                              nan    |                             0.456 | none                        |
| HLA-DQA2 | HLA-DQA2*01:04          |            1 | False       |                          21 | nan            |              nan     |                              nan    |                             0.431 | none                        |
| MICB     | MICB*005:01             |            1 | nan         |                          21 | AMR            |                1     |                               14.2  |                             0.847 | AMR:20                      |


### 1b. Novel synonymous CDS — same protein, new coding sequence (9)

These are *not* new proteins (`n_aa_diffs` is 0 by construction) and must not be counted in a 'new HLA proteins' sentence. They are still real, uncatalogued coding sequences.

| gene     | nearest_known_protein   |   n_aa_diffs |   in_groove |   n_persons_unrelated_clean | top_ancestry   |   top_ancestry_share |   carriers_per_1000_in_top_ancestry |   pct_gene_haplotypes_not_in_imgt | ancestries_disclosable   |
|:---------|:------------------------|-------------:|------------:|----------------------------:|:---------------|---------------------:|------------------------------------:|----------------------------------:|:-------------------------|
| TAP2     | TAP2*02:01              |            0 |         nan |                         102 | AFR            |                    1 |                               31.97 |                             6.847 | AFR:68                   |
| HLA-G    | HLA-G*01:01             |            0 |         nan |                          83 | EAS            |                    1 |                               62.37 |                             0.323 | EAS:73                   |
| HLA-DQB2 | HLA-DQB2*01:09          |            0 |         nan |                          60 | AFR            |                    1 |                               14.9  |                             1.263 | AFR:33                   |
| HLA-DOA  | HLA-DOA*01:05           |            0 |         nan |                          52 | AFR            |                    1 |                               16.25 |                             0.257 | AFR:36                   |
| HLA-DMA  | HLA-DMA*01:01           |            0 |         nan |                          50 | AFR            |                    1 |                               15.46 |                             0.403 | AFR:34                   |
| TAP2     | TAP2*01:01              |            0 |         nan |                          35 | AFR            |                    1 |                                9.87 |                             6.847 | AFR:21                   |
| TAP1     | TAP1*01:01              |            0 |         nan |                          30 | EAS            |                    1 |                               24.21 |                             7.981 | EAS:28                   |
| TAP2     | TAP2*02:03              |            0 |         nan |                          27 | nan            |                  nan |                              nan    |                             6.847 | none                     |
| MICB     | MICB*005:02             |            0 |         nan |                          20 | nan            |                  nan |                              nan    |                             0.847 | none                     |


**The reference-depth control.** `pct_gene_haplotypes_not_in_imgt` is the share of this cohort's haplotypes at that gene whose protein is already absent from IPD-IMGT/HLA, measured independently in script 27. Read every callout against it: a gene with a large gap was under-catalogued to begin with. The classical genes sit at 0.1-0.2%; TAP1/TAP2 are the extreme at ~7-8%.


`in_groove` is whether any of the amino-acid differences fall in the peptide-binding groove exons — a groove difference is the one most likely to change which peptides the molecule presents, so those are the interesting ones biologically, not just numerically.


## Tier 2 — IMGT submission queue (309 alleles)

Clean and seen in at least two unrelated people, but with fewer than 20 carriers, so their counts stay censored. Too many to name individually; this is the list to submit, not to quote.

| gene     |   n_alleles |
|:---------|------------:|
| TAP2     |          66 |
| TAP1     |          58 |
| HLA-DQB2 |          21 |
| HLA-DMA  |          20 |
| MICA     |          20 |
| HLA-DQA2 |          19 |
| HLA-DMB  |          16 |
| HLA-DOA  |          15 |
| HLA-G    |          14 |
| MICB     |          13 |
| HLA-F    |          13 |
| HLA-DOB  |           9 |
| HLA-HFE  |           8 |
| HLA-E    |           3 |
| HLA-DRA  |           3 |
| HLA-DRB4 |           2 |
| HLA-DRB1 |           2 |
| HLA-DPB1 |           2 |
| HLA-DRB5 |           1 |
| HLA-C    |           1 |
| HLA-DQB1 |           1 |
| HLA-DPA1 |           1 |
| HLA-B    |           1 |


## How to read the ancestry columns

`top_ancestry_share` is the fraction of **disclosable** carriers in the single most common ancestry. `ancestries_censored` lists ancestries where this allele has carriers but fewer than 20, so the true share is at least as concentrated as shown and possibly less. Rates per 1,000 use script 27's per-ancestry haplotype counts divided by two; where an assembly was fragmented the person was still typed on one haplotype, so these rates are slight under-estimates.


## What this list is not

- **Not validated orthogonally.** Every allele here rests on the long-read assembly alone. Cole's own suggestion — realigning the short reads, which every one of these people also has, to the long-read assembly — is the confirmation step, and he explicitly parked it for now.

- **Not an IMGT submission yet.** Submission needs full-length genomic sequence and a documented typing method; what we have is the candidate set.

- **Not a frequency estimate for the general population.** These are All of Us participants with long-read data, which is not a random sample of anything.


## Non-coding novelty

Script 25 resolved the pooled non-coding clusters into 184 distinct genomic signatures. Those are catalogue gaps in IPD-IMGT's *genomic* sequences, not new proteins, and they are reported separately in `reports/hla_popgen/25_noncoding_novelty_paf/`.


## Files

- `main_text_candidates.tsv` — tier 1.

- `imgt_submission_candidates.tsv` — tier 2.

- `all_novel_clusters_annotated.tsv` — every cluster with the ancestry annotation.

- `fig_novel_callouts.png`.
