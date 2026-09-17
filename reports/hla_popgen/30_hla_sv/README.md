# 30 — structural variation across the MHC: deletions, duplications, copy number

*Supervisor asks A4 (deletions/duplications) and A5 (KIR), 2026-09-17 call.*

## The method, and why absence is not deletion

A gene missing from a haplotype can mean a real deletion, a broken assembly, or a failed homology search. This script calls a deletion **only** when a single contig carries a gene on each side of the missing gene in canonical physical order — i.e. the assembly was built through the gene's position and the gene was not there. Absences without flanking evidence are counted separately and never interpreted.

- unrelated people: **11856**; haplotypes examined: **23699**

- bridged absences: **935929**; unbridged (assembly ends): **3032899**


## 1. Positive control — DRB3/DRB4/DRB5 against the DR51/52/53 expectation

DRB1 haplotype group determines which second DRB locus is present. The answer is known from immunogenetics, independently of our data, so this measures whether the deletion caller works at all.

|   drb1_group | expected   |   pct_concordant |   n_disp |   n_concordant_disp |
|-------------:|:-----------|-----------------:|---------:|--------------------:|
|           01 | none       |          97.6349 |     1649 |                1610 |
|           03 | DRB3       |          97.6056 |     2130 |                2079 |
|           04 | DRB4       |          92.1966 |     2909 |                2682 |
|           07 | DRB4       |          93.2677 |     2436 |                2272 |
|           08 | none       |          97.4603 |     1260 |                1228 |
|           09 | DRB4       |          94.1558 |      616 |                 580 |
|           10 | none       |          97.1302 |      453 |                 440 |
|           11 | DRB3       |          97.5562 |     2578 |                2515 |
|           12 | DRB3       |          96.973  |      925 |                 897 |
|           13 | DRB3       |          97.7537 |     2671 |                2611 |
|           14 | DRB3       |          97.8846 |     1040 |                1018 |
|           15 | DRB5       |          95.4437 |     2941 |                2807 |
|           16 | DRB5       |          96.6346 |      416 |                 402 |


**Read this first.** If concordance is high, the deletion calls below are trustworthy. If it is not, nothing else in this report is.


## 2. Negative control — genes that should never be deleted

| gene   |   n_bridged_disp |   pct_deleted_bridged |
|:-------|-----------------:|----------------------:|
| A      |            23271 |                 1.831 |
| B      |            23039 |                 0.603 |
| C      |            22883 |                 0.035 |
| DPA1   |            23195 |                 1.091 |
| DPB1   |            22982 |                 0.996 |
| DQA1   |            23088 |                 1.182 |
| DQB1   |            23219 |                 1.998 |
| DRA    |            23080 |                 1.629 |
| DRB1   |            22892 |                 1.647 |


This is the **false-positive rate** of the method, not biology. Every rate in the next section should be read against it.


## 3. Deletion frequency per gene

| gene   | gene_class         |   n_bridged_disp |   pct_deleted_bridged |
|:-------|:-------------------|-----------------:|----------------------:|
| DRB5   | class_II_paralog   |            20588 |                83.107 |
| Y      | pseudogene_I       |            20589 |                81.602 |
| DRB4   | class_II_paralog   |            20446 |                69.642 |
| DRB3   | class_II_paralog   |            19121 |                49.297 |
| C4B    | complement         |            21407 |                19.554 |
| K      | pseudogene_I       |            22061 |                12.42  |
| U      | pseudogene_I       |            22080 |                12.373 |
| T      | pseudogene_I       |            22171 |                12.25  |
| H      | pseudogene_I       |            22165 |                12.15  |
| C4A    | complement         |            21978 |                11.025 |
| MICA   | mic_tap            |            22801 |                 4.202 |
| MICB   | mic_tap            |            22895 |                 2.87  |
| W      | pseudogene_I       |            23326 |                 2.208 |
| DQB1   | classical_II       |            23219 |                 1.998 |
| A      | classical_I        |            23271 |                 1.831 |
| DRB1   | classical_II       |            22892 |                 1.647 |
| DRA    | class_II_accessory |            23080 |                 1.629 |
| S      | pseudogene_I       |            23268 |                 1.354 |
| DQA1   | classical_II       |            23088 |                 1.182 |
| DPA1   | classical_II       |            23195 |                 1.091 |
| TAP1   | mic_tap            |            22941 |                 1.081 |
| DPB1   | classical_II       |            22982 |                 0.996 |
| DPA2   | class_II_paralog   |            23080 |                 0.988 |
| F      | nonclassical_I     |            23040 |                 0.985 |
| P      | pseudogene_I       |            23191 |                 0.979 |


## 4. Deletion frequency by ancestry

| gene   |   AFR |   AMR |   EAS |   EUR |   MID |   SAS |
|:-------|------:|------:|------:|------:|------:|------:|
| A      |  2.03 |  1.64 |  1.72 |  1.46 |  1.83 |  1.58 |
| B      |  0.64 |  0.68 |  0.53 |  0.64 |  0.5  |  0.58 |
| C      |  0.02 |  0.06 |  0    |  0.04 |  0    |  0.09 |
| C4A    |  9.77 |  9.79 |  9.92 | 16.43 | 10.65 |  8.15 |
| C4B    | 21.66 | 16.9  | 18.24 | 21.73 | 18.62 | 19.09 |
| DMA    |  0.84 |  0.65 |  0.44 |  0.52 |  0.67 |  0.62 |
| DMB    |  0.95 |  0.55 |  0.44 |  0.63 |  0.84 |  0.62 |
| DOA    |  0.63 |  0.58 |  0.56 |  0.55 |  0.84 |  0.79 |
| DOB    |  0.51 |  0.32 |  0.44 |  0.48 |  0.17 |  0.7  |
| DPA1   |  1.08 |  0.98 |  1.28 |  0.99 |  1    |  1.24 |
| DPA2   |  0.91 |  0.91 |  1.34 |  1.15 |  1.18 |  1.28 |
| DPB1   |  0.98 |  0.98 |  1.06 |  1.24 |  1.01 |  1.16 |
| DPB2   |  0    |  0    |  0    |  0    |  0    |  0    |
| DQA1   |  1.17 |  0.88 |  0.81 |  1.43 |  1.17 |  1.1  |
| DQA2   |  0.64 |  0.45 |  0.56 |  0.59 |  0.33 |  0.75 |
| DQB1   |  2.09 |  1.68 |  1.56 |  2.33 |  1.17 |  1.41 |
| DQB2   |  0.62 |  0.52 |  0.48 |  0.63 |  0.67 |  0.62 |
| DRA    |  1.68 |  1.2  |  1.41 |  1.7  |  1.18 |  2.17 |
| DRB1   |  1.36 |  1.25 |  1.58 |  1.8  |  1.7  |  2.04 |
| DRB3   | 40.64 | 57.22 | 54.17 | 52.03 | 45.97 | 54.8  |
| DRB4   | 79.55 | 63.38 | 65.8  | 65.31 | 69.87 | 69.42 |
| DRB5   | 81.64 | 88.47 | 77.71 | 83.03 | 89.86 | 74.8  |
| E      |  0.02 |  0.03 |  0    |  0    |  0    |  0.04 |
| F      |  1.02 |  0.49 |  1.34 |  1.12 |  0.68 |  0.62 |
| G      |  0.92 |  0.61 |  1.09 |  0.87 |  1.35 |  0.57 |
| H      | 10.06 | 14.23 | 18.74 |  7.73 | 14.16 | 13.83 |
| HFE    |  0    |  0    |  0    |  0    |  0    |  0    |
| J      |  0.64 |  0.39 |  0.68 |  0.76 |  0.34 |  0.58 |
| K      | 10.27 | 14.54 | 19.31 |  8.14 | 14.06 | 14.13 |
| L      |  0.22 |  0.42 |  0.81 |  0.44 |  0    |  0.35 |
| MICA   |  4.28 |  6.03 |  4.32 |  2.7  |  1.69 |  3.2  |
| MICB   |  3.19 |  2.13 |  2.87 |  2.78 |  2.37 |  3.24 |
| N      |  0.88 |  0.77 |  0.92 |  0.91 |  0.66 |  1.14 |
| P      |  0.97 |  0.71 |  1.05 |  1.07 |  1.18 |  0.71 |
| S      |  1.32 |  1.23 |  1.4  |  1.14 |  0.84 |  1.32 |
| T      | 10.08 | 14.2  | 18.99 |  7.79 | 14.51 | 13.93 |
| TAP1   |  1.16 |  0.92 |  1.17 |  1    |  0.68 |  1.2  |
| TAP2   |  0.91 |  0.59 |  1.05 |  0.84 |  0.68 |  0.8  |
| U      | 10.21 | 14.48 | 19.12 |  8.03 | 14.01 | 14.24 |
| V      |  0.97 |  0.68 |  1.09 |  1.23 |  0.68 |  0.44 |
| W      |  2.9  |  1.41 |  4.73 |  1.42 |  1.16 |  0.79 |
| Y      | 74.84 | 79.69 | 76.93 | 93.88 | 83.88 | 82.15 |


## 5. Duplication candidates (>1 copy of a gene on one contig)

| gene   |   n_haplotype_contigs_disp |   n_multi_copy_disp |   pct_multi_copy |   max_copy_seen |
|:-------|---------------------------:|--------------------:|-----------------:|----------------:|
| C4A    |                      20010 |                2972 |         14.8526  |               5 |
| C4B    |                      17524 |                1818 |         10.3743  |               4 |
| DQA1   |                      23363 |                 550 |          2.35415 |               2 |
| DRB1   |                      23027 |                 513 |          2.22782 |               2 |
| DQB1   |                      23255 |                 503 |          2.16298 |               2 |
| DRA    |                      23155 |                 451 |          1.94774 |               2 |
| F      |                      23233 |                 421 |          1.81208 |               3 |
| DQA2   |                      23446 |                 408 |          1.74017 |               3 |
| DQB2   |                      23418 |                 387 |          1.65257 |               2 |
| DOB    |                      23436 |                 378 |          1.6129  |               2 |
| V      |                      23348 |                 367 |          1.57187 |               3 |
| B      |                      23262 |                 363 |          1.56048 |               3 |
| S      |                      23298 |                 360 |          1.5452  |               3 |
| P      |                      23321 |                 360 |          1.54367 |               2 |
| C      |                      23226 |                 351 |          1.51124 |               2 |
| J      |                      23339 |                 351 |          1.50392 |               2 |
| DOA    |                      23461 |                 351 |          1.4961  |               3 |
| DMA    |                      23326 |                 348 |          1.4919  |               3 |
| W      |                      23141 |                 336 |          1.45197 |               2 |
| MICA   |                      22110 |                 317 |          1.43374 |               2 |


These are **candidates**. A second mapping cluster inside a segmental duplication can be an alignment artifact; confirming a duplication needs read-depth or the assembly graph, neither of which is used here.


## 6. C4 copy number

|   C4A |   C4B |   total_c4 |   n_haplotypes_disp |       pct |
|------:|------:|-----------:|--------------------:|----------:|
|     1 |     1 |          2 |               12046 | 54.0786   |
|     1 |     0 |          1 |                3296 | 14.7969   |
|     0 |     1 |          1 |                2336 | 10.4871   |
|     2 |     0 |          2 |                1633 |  7.33109  |
|     1 |     2 |          3 |                1233 |  5.53535  |
|     2 |     1 |          3 |                 990 |  4.44444  |
|     0 |     2 |          2 |                 362 |  1.62514  |
|     3 |     0 |          3 |                 115 |  0.516274 |
|     2 |     2 |          4 |                  91 |  0.40853  |
|     1 |     3 |          4 |                  56 |  0.251403 |


Long/short (intron-9 HERV) composition:


| gene_bare   | c4_size   |   n_disp |   pct_within_gene |
|:------------|:----------|---------:|------------------:|
| C4A         | L         |    20460 |          90.1679  |
| C4A         | S         |     2231 |           9.83209 |
| C4B         | L         |     7996 |          41.831   |
| C4B         | S         |    11119 |          58.169   |


## 7. KIR (ask A5)

- KIR rows in the calls table: **0**

- EXPECTED: no KIR calls; the KIR cluster is on chr19, outside the chr6 trim window used to build these assemblies. Genotyping KIR would need a separate extraction from the original BAMs.


## Caveats

- Deletion calls are annotation-level, not sequence-level. A confirmed deletion would show the breakpoint; we show only that the assembled contig skipped the gene's position.

- Canonical gene order is derived empirically from the cohort's own assemblies (script 16's `derive_canonical_order`), because `gene_start` is contig-relative, not hg38 (SCHEMA.md Table 1).

- Duplication candidates are not validated by depth.

- **A gene at either end of the annotated region can never be bridged**, because nothing flanks it on one side. Such a gene will always report a 0% deletion rate here — that is 'not testable by this method', not 'never deleted'. Check a gene's `n_bridged` before reading its rate.
