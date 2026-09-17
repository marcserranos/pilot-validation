# 29 — LD between HLA alleles at phased cis gene pairs, within ancestry

*Supervisor ask A3 (Cole, 2026-09-17): "compute the LD between alleles within each ancestry and see what looks different", r^2, for DQA1–DQB1 and DPA1–DPB1.*

## What was computed

Haplotypes are **physically phased**: both genes annotated on the same contig of the same assembly haplotype, one copy each, both resolvable to two fields (protein level). People are restricted to a greedy unrelated set (kinship < 0.0442) and, for the per-ancestry tables, to strict ancestry (probability >= 0.90).

- unrelated people used: **11856**


### Why three tiers of statistic

r^2 between two alleles is bounded by their marginal frequencies, and every multi-allelic association statistic is inflated when there are many alleles relative to sample size. Ancestries differ in both. So a raw cross-ancestry comparison of either quantity measures the frequency spectrum and the allele count, not linkage. The defensible comparison is the **rarefied** one: every ancestry subsampled to the same number of phased haplotypes, statistics recomputed, 95% interval over subsamples.


## 1. Phasing yield per interval

| pair      |   both genes present |   phased (same contig) |   % phased | dropped: >1 copy   |   dropped: unresolved call |
|:----------|---------------------:|-----------------------:|-----------:|:-------------------|---------------------------:|
| DQA1~DQB1 |                22481 |                  21652 |      96.31 | <20                |                        781 |
| DPA1~DPB1 |                22601 |                  21913 |      96.96 | 0                  |                        911 |
| DRB1~DQB1 |                22055 |                  20229 |      91.72 | <20                |                        822 |
| B~C       |                22416 |                  19474 |      86.88 | <20                |                       1138 |
| A~B       |                22187 |                   6709 |      30.24 | <20                |                         94 |


A low percentage here is assembly fragmentation across that interval, not a typing failure. DPA1–DPB1 sit ~10 kb apart and should be near-100%; the class I pairs span far more sequence and should be lower.


## 2. Multi-allelic LD per ancestry (raw, NOT comparable across ancestries)

| pair      | ancestry   |   n_haplotypes_disp |   n_alleles_a |   n_alleles_b |   Dprime_multi |   cramers_v |   cramers_v_bc |   nmi |
|:----------|:-----------|--------------------:|--------------:|--------------:|---------------:|------------:|---------------:|------:|
| DQA1~DQB1 | AFR        |                4146 |            24 |            31 |          0.905 |       0.511 |          0.505 | 0.741 |
| DQA1~DQB1 | AMR        |                2913 |            21 |            30 |          0.932 |       0.615 |          0.609 | 0.794 |
| DQA1~DQB1 | EAS        |                2377 |            19 |            18 |          0.95  |       0.679 |          0.676 | 0.81  |
| DQA1~DQB1 | EUR        |                2383 |            21 |            22 |          0.972 |       0.624 |          0.619 | 0.86  |
| DQA1~DQB1 | MID        |                 563 |            18 |            24 |          0.96  |       0.676 |          0.655 | 0.842 |
| DQA1~DQB1 | SAS        |                2147 |            18 |            26 |          0.946 |       0.648 |          0.641 | 0.803 |
| DPA1~DPB1 | AFR        |                4275 |            29 |            62 |          0.888 |       0.548 |          0.536 | 0.75  |
| DPA1~DPB1 | AMR        |                2930 |            13 |            55 |          0.944 |       0.717 |          0.706 | 0.832 |
| DPA1~DPB1 | EAS        |                2369 |            12 |            34 |          0.867 |       0.613 |          0.603 | 0.736 |
| DPA1~DPB1 | EUR        |                2362 |             9 |            32 |          0.934 |       0.718 |          0.71  | 0.808 |
| DPA1~DPB1 | MID        |                 570 |             8 |            25 |          0.837 |       0.578 |          0.544 | 0.688 |
| DPA1~DPB1 | SAS        |                2145 |            12 |            36 |          0.843 |       0.567 |          0.554 | 0.701 |
| DRB1~DQB1 | AFR        |                3886 |            54 |            31 |          0.871 |       0.525 |          0.514 | 0.72  |
| DRB1~DQB1 | AMR        |                2729 |            60 |            30 |          0.909 |       0.587 |          0.571 | 0.806 |
| DRB1~DQB1 | EAS        |                2227 |            53 |            18 |          0.913 |       0.753 |          0.74  | 0.807 |
| DRB1~DQB1 | EUR        |                2200 |            47 |            23 |          0.94  |       0.698 |          0.686 | 0.844 |
| DRB1~DQB1 | MID        |                 534 |            35 |            24 |          0.949 |       0.687 |          0.653 | 0.846 |
| DRB1~DQB1 | SAS        |                2029 |            53 |            25 |          0.939 |       0.655 |          0.639 | 0.809 |
| B~C       | AFR        |                3803 |            90 |            55 |          0.857 |       0.519 |          0.5   | 0.7   |
| B~C       | AMR        |                2640 |           113 |            48 |          0.866 |       0.591 |          0.559 | 0.714 |
| B~C       | EAS        |                2087 |            80 |            37 |          0.853 |       0.651 |          0.627 | 0.724 |
| B~C       | EUR        |                2141 |            69 |            37 |          0.893 |       0.675 |          0.657 | 0.763 |
| B~C       | MID        |                 512 |            57 |            33 |          0.921 |       0.715 |          0.655 | 0.804 |
| B~C       | SAS        |                1926 |            64 |            42 |          0.895 |       0.663 |          0.645 | 0.773 |
| A~B       | AFR        |                1066 |            39 |            66 |          0.432 |       0.319 |          0.205 | 0.221 |
| A~B       | AMR        |                1339 |            52 |            90 |          0.488 |       0.317 |          0.188 | 0.305 |
| A~B       | EAS        |                 659 |            28 |            65 |          0.49  |       0.484 |          0.378 | 0.341 |
| A~B       | EUR        |                 761 |            28 |            47 |          0.433 |       0.348 |          0.251 | 0.254 |
| A~B       | MID        |                 175 |            24 |            41 |          0.713 |       0.506 |          0.173 | 0.506 |
| A~B       | SAS        |                 641 |            32 |            50 |          0.477 |       0.397 |          0.292 | 0.287 |


## 3. Rarefied multi-allelic LD (the cross-ancestry comparison)

| pair      | ancestry   |   rarefied_to |   Dprime_multi_mean |   Dprime_multi_lo |   Dprime_multi_hi |   cramers_v_bc_mean |   n_alleles_a_mean |   n_alleles_b_mean |
|:----------|:-----------|--------------:|--------------------:|------------------:|------------------:|--------------------:|-------------------:|-------------------:|
| DQA1~DQB1 | AFR        |           563 |               0.908 |             0.89  |             0.925 |               0.632 |               15.3 |               19.8 |
| DQA1~DQB1 | AMR        |           563 |               0.936 |             0.921 |             0.951 |               0.724 |               15.1 |               19.7 |
| DQA1~DQB1 | EAS        |           563 |               0.951 |             0.941 |             0.961 |               0.714 |               16.6 |               16.4 |
| DQA1~DQB1 | EUR        |           563 |               0.974 |             0.964 |             0.983 |               0.71  |               16.2 |               17.6 |
| DQA1~DQB1 | MID        |           563 |               0.96  |             0.96  |             0.96  |               0.655 |               18   |               24   |
| DQA1~DQB1 | SAS        |           563 |               0.947 |             0.934 |             0.959 |               0.685 |               15.7 |               18.9 |
| DPA1~DPB1 | AFR        |           570 |               0.892 |             0.869 |             0.913 |               0.748 |               13.1 |               33.8 |
| DPA1~DPB1 | AMR        |           570 |               0.946 |             0.92  |             0.972 |               0.814 |                8.4 |               32.8 |
| DPA1~DPB1 | EAS        |           570 |               0.867 |             0.835 |             0.899 |               0.724 |                7.2 |               25.1 |
| DPA1~DPB1 | EUR        |           570 |               0.935 |             0.896 |             0.973 |               0.747 |                6.6 |               23.1 |
| DPA1~DPB1 | MID        |           570 |               0.837 |             0.837 |             0.837 |               0.544 |                8   |               25   |
| DPA1~DPB1 | SAS        |           570 |               0.844 |             0.808 |             0.879 |               0.677 |                8   |               24.6 |
| DRB1~DQB1 | AFR        |           534 |               0.88  |             0.86  |             0.899 |               0.626 |               34.9 |               19.4 |
| DRB1~DQB1 | AMR        |           534 |               0.92  |             0.901 |             0.938 |               0.695 |               45.2 |               19.9 |
| DRB1~DQB1 | EAS        |           534 |               0.916 |             0.898 |             0.933 |               0.767 |               35.8 |               16.3 |
| DRB1~DQB1 | EUR        |           534 |               0.947 |             0.932 |             0.96  |               0.777 |               33.9 |               17.6 |
| DRB1~DQB1 | MID        |           534 |               0.949 |             0.949 |             0.949 |               0.653 |               35   |               24   |
| DRB1~DQB1 | SAS        |           534 |               0.945 |             0.928 |             0.961 |               0.725 |               38.9 |               18.4 |
| B~C       | AFR        |           512 |               0.88  |             0.857 |             0.9   |               0.615 |               53.7 |               30.4 |
| B~C       | AMR        |           512 |               0.896 |             0.877 |             0.915 |               0.588 |               74.4 |               34.8 |
| B~C       | EAS        |           512 |               0.874 |             0.851 |             0.896 |               0.663 |               58.3 |               28.9 |
| B~C       | EUR        |           512 |               0.905 |             0.885 |             0.924 |               0.713 |               45.5 |               25   |
| B~C       | MID        |           512 |               0.921 |             0.921 |             0.921 |               0.655 |               57   |               33   |
| B~C       | SAS        |           512 |               0.911 |             0.893 |             0.928 |               0.72  |               48.8 |               29   |
| A~B       | AFR        |           175 |               0.692 |             0.645 |             0.736 |               0.138 |               26.9 |               38.4 |
| A~B       | AMR        |           175 |               0.739 |             0.685 |             0.79  |               0.245 |               29.9 |               51   |
| A~B       | EAS        |           175 |               0.669 |             0.619 |             0.72  |               0.358 |               21.3 |               45.3 |
| A~B       | EUR        |           175 |               0.583 |             0.522 |             0.645 |               0.316 |               20.6 |               32.3 |
| A~B       | MID        |           175 |               0.713 |             0.713 |             0.713 |               0.173 |               24   |               41   |
| A~B       | SAS        |           175 |               0.641 |             0.591 |             0.697 |               0.27  |               23   |               37.1 |


Read this table as: *at equal sample size*, does the same pair of loci travel together more tightly in one ancestry than another? Non-overlapping CIs are the claim; overlapping CIs are not.


## 4. Strongest individual allele pairs (what Cole asked for literally)

| pair      | ancestry   | allele_a   | allele_b    |   freq_a |   freq_b |   freq_hap |    r2 |   Dprime |
|:----------|:-----------|:-----------|:------------|---------:|---------:|-----------:|------:|---------:|
| DQA1~DQB1 | EAS        | DQA1*05:01 | DQB1*02:01  |    0.053 |    0.053 |      0.053 | 1     |    1     |
| DQA1~DQB1 | MID        | DQA1*05:01 | DQB1*02:01  |    0.107 |    0.107 |      0.107 | 1     |    1     |
| DQA1~DQB1 | AMR        | DQA1*01:04 | DQB1*05:03  |    0.016 |    0.016 |      0.016 | 1     |    1     |
| DQA1~DQB1 | EUR        | DQA1*05:01 | DQB1*02:01  |    0.124 |    0.124 |      0.124 | 0.996 |    1     |
| DQA1~DQB1 | SAS        | DQA1*05:01 | DQB1*02:01  |    0.072 |    0.072 |      0.072 | 0.993 |    1     |
| DPA1~DPB1 | AFR        | DPA1*04:02 | DPB1*665:01 |    0.007 |    0.007 |      0.007 | 1     |    1     |
| DPA1~DPB1 | AFR        | DPA1*02:12 | DPB1*85:01  |    0.015 |    0.014 |      0.014 | 0.936 |    0.983 |
| DPA1~DPB1 | EAS        | DPA1*04:01 | DPB1*107:01 |    0.031 |    0.026 |      0.026 | 0.82  |    1     |
| DPA1~DPB1 | AFR        | DPA1*03:01 | DPB1*105:01 |    0.104 |    0.103 |      0.095 | 0.813 |    0.909 |
| DPA1~DPB1 | AMR        | DPA1*03:01 | DPB1*105:01 |    0.018 |    0.019 |      0.016 | 0.768 |    0.902 |
| DRB1~DQB1 | EAS        | DRB1*03:01 | DQB1*02:01  |    0.054 |    0.054 |      0.054 | 1     |    1     |
| DRB1~DQB1 | EAS        | DRB1*13:01 | DQB1*06:03  |    0.01  |    0.01  |      0.01  | 1     |    1     |
| DRB1~DQB1 | EUR        | DRB1*03:01 | DQB1*02:01  |    0.125 |    0.126 |      0.124 | 0.967 |    0.988 |
| DRB1~DQB1 | MID        | DRB1*15:02 | DQB1*06:01  |    0.043 |    0.045 |      0.043 | 0.956 |    1     |
| DRB1~DQB1 | SAS        | DRB1*03:01 | DQB1*02:01  |    0.075 |    0.075 |      0.073 | 0.937 |    0.972 |
| B~C       | EAS        | B*07:05    | C*15:05     |    0.014 |    0.014 |      0.014 | 0.966 |    1     |
| B~C       | EAS        | B*58:01    | C*03:02     |    0.072 |    0.073 |      0.07  | 0.923 |    0.964 |
| B~C       | MID        | B*14:02    | C*08:02     |    0.045 |    0.051 |      0.045 | 0.879 |    1     |
| B~C       | SAS        | B*07:05    | C*15:05     |    0.012 |    0.015 |      0.012 | 0.855 |    1     |
| B~C       | SAS        | B*44:03    | C*07:06     |    0.093 |    0.083 |      0.082 | 0.854 |    0.986 |
| A~B       | EUR        | A*29:02    | B*44:03     |    0.034 |    0.046 |      0.017 | 0.166 |    0.476 |
| A~B       | EAS        | A*02:07    | B*46:01     |    0.099 |    0.094 |      0.038 | 0.108 |    0.338 |
| A~B       | EAS        | A*33:03    | B*58:01     |    0.067 |    0.038 |      0.018 | 0.108 |    0.443 |
| A~B       | AMR        | A*01:01    | B*08:01     |    0.054 |    0.028 |      0.013 | 0.105 |    0.457 |
| A~B       | EUR        | A*01:01    | B*08:01     |    0.155 |    0.116 |      0.055 | 0.104 |    0.381 |


The marginal frequencies are printed next to every r^2 on purpose: r^2 is bounded by them, so a small r^2 at freq 0.02 can be *stronger* linkage than a larger r^2 at freq 0.4.


## Files

- `ld_pairwise.tsv` — every allele pair passing the 20-haplotype floor, per ancestry, with r^2, D', D and marginals.

- `ld_multiallelic.tsv` — whole-table statistics per pair per ancestry (raw).

- `ld_rarefied.tsv` — the same, rarefied to a common haplotype count.

- `phasing_yield.tsv` — physical-phasing yield per interval.

- `haplotype_freqs.tsv` — two-locus haplotype frequencies per ancestry (counts below 20 written `<20`).

- `fig_ld_multiallelic.png`, `fig_r2_<pair>.png`, `fig_phasing_yield.png`.


## Caveats

- Physical phase is Immuannot's assembly phase. S01's QC (script 26) found zero phase switches over 3,021 testable transitions with 100% power on injected switches, which is what licenses treating these as true haplotypes.

- Two-field resolution. Calls novel at field 1 or 2 are dropped, not merged.

- Strict ancestry shrinks AMR and MID most; check `n_haplotypes` before reading anything into those rows.
