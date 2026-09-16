# Novelty by nomenclature field (script 24)

## What this measures, in plain language

Immuannot writes `new` in place of the first allele field it cannot match. The position of `new` says what is new. Depth 2 (`HLA-A*01:new`) means a new protein. Depth 3 (`HLA-A*01:01:new`) means a new coding sequence that encodes a known protein. Depth 4 (`HLA-A*01:01:01:new`) means the coding sequence is catalogued and only introns or UTRs differ. Depth 1 means even the first field is unresolved. (The WS1 brief's background table shows these strings one field shallower. This script follows SCHEMA.md and 01_extract_rich.py.)

Each haplotype-gene call is counted once (copy 1). Artifact flags are *labels*, not deletions. `homopolymer_indel` means every coding difference is a single-base-run indel, the typical HiFi error. `partial_cds` and `inframe_stop` are Immuannot's reconstruction warnings. A flagged call is not proven wrong, so both views are shown.

For depth-2/3 calls the observed CDS was then compared with the IPD-IMGT snapshot that Immuannot ships. Some 'new' names turn out to be catalogued sequences (`cds_known`) or catalogued proteins (`protein_known`). Frameshifts and premature stops are separated out before any protein is called novel.

Counts of 1-19 are shown as `<20` (AoU policy). Relatives: 377 people were removed greedily for kin >= 0.0442, leaving 11856 unrelated people out of 12233.

## Headline numbers

| metric | value |
|---|---|
| n_haplotype_gene_calls | 875179 |
| n_novel_calls_any_depth | 280695 |
| n_f2_protein_calls | 75442 |
| n_f2_protein_calls_clean | 22122 |
| n_f3_synonymous_calls | 30795 |
| n_f4_noncoding_calls | 172272 |
| n_depth23_calls_with_sequence | 104532 |
| n_cds_known_naming_artifacts | 25282 |
| n_frameshift_or_stop_calls | 72667 |
| n_novel_protein_calls | 5288 |
| n_novel_protein_clusters | 1404 |
| n_novel_protein_clusters_all_clean | 1026 |
| n_novel_protein_clusters_clean_recurrent_unrelated | 231 |
| n_novel_cds_synonymous_clusters | 419 |
| frac_f2_calls_flagged | 0.707 |

## A. Field class per gene (all people, pooled; fraction of calls; flagged share in brackets)

| gene | n_calls | known | f4_noncoding | f3_synonymous | f2_protein | f1_undetermined | uncalled |
|---|---|---|---|---|---|---|---|
| HLA-A | 23571 | 0.940 [0.00] | 0.031 [0.02] | 0.000 [0.00] | 0.028 [0.93] | 0 | 0 |
| HLA-B | 23630 | 0.870 [0.00] | 0.097 [0.00] | 0.000 [0.00] | 0.033 [0.89] | 0.000 [0.00] | 0 |
| HLA-C | 23605 | 0.893 [0.00] | 0.071 [0.01] | 0.000 [0.33] | 0.036 [0.91] | 0 | 0 |
| HLA-DPA1 | 23676 | 0.870 [0.00] | 0.112 [0.00] | 0.000 [0.00] | 0.018 [0.93] | 0 | 0 |
| HLA-DPB1 | 23473 | 0.785 [0.00] | 0.190 [0.01] | 0.000 [0.00] | 0.025 [0.94] | 0.000 [1.00] | 0 |
| HLA-DQA1 | 23544 | 0.839 [0.00] | 0.144 [0.00] | 0.000 [0.00] | 0.017 [0.88] | 0 | 0.000 [1.00] |
| HLA-DQB1 | 23482 | 0.884 [0.00] | 0.093 [0.03] | 0.000 [0.22] | 0.023 [0.93] | 0 | 0 |
| HLA-DRB1 | 23235 | 0.375 [0.00] | 0.603 [0.00] | 0.000 [0.00] | 0.022 [0.93] | 0 | 0 |
| HLA-E | 23431 | 0.911 [0.00] | 0.057 [0.01] | 0.001 [0.08] | 0.032 [0.90] | 0 | 0 |
| HLA-F | 23542 | 0.878 [0.00] | 0.085 [0.01] | 0.001 [0.00] | 0.036 [0.80] | 0 | 0.000 [0.00] |
| HLA-G | 23718 | 0.918 [0.03] | 0.047 [0.00] | 0.005 [0.00] | 0.029 [0.85] | 0 | 0 |
| HLA-DMA | 23718 | 0.907 [0.00] | 0.068 [0.00] | 0.004 [0.00] | 0.021 [0.79] | 0 | 0 |
| HLA-DMB | 23738 | 0.701 [0.00] | 0.276 [0.00] | 0.002 [0.00] | 0.021 [0.75] | 0 | 0 |
| HLA-DOA | 23847 | 0.888 [0.00] | 0.084 [0.00] | 0.002 [0.00] | 0.026 [0.78] | 0 | 0 |
| HLA-DOB | 23800 | 0.820 [0.00] | 0.165 [0.00] | 0.001 [0.00] | 0.014 [0.80] | 0 | 0 |
| HLA-DRA | 23423 | 0.864 [0.00] | 0.112 [0.00] | 0.001 [0.00] | 0.023 [0.88] | 0 | 0 |
| HLA-DPA2 | 23579 | 0.341 [1.00] | 0.240 [1.00] | 0.030 [1.00] | 0.388 [0.99] | 0.000 [1.00] | 0.000 [1.00] |
| HLA-DPB2 | 23223 | 0.132 [1.00] | 0.674 [1.00] | 0.011 [1.00] | 0.182 [0.99] | 0 | 0.000 [0.50] |
| HLA-DQA2 | 23775 | 0.547 [0.00] | 0.434 [0.00] | 0.001 [0.00] | 0.017 [0.71] | 0 | 0.000 [0.67] |
| HLA-DQB2 | 23770 | 0.601 [0.00] | 0.350 [0.00] | 0.002 [0.00] | 0.047 [0.32] | 0 | 0 |
| HLA-DRB3 | 9991 | 0.344 [0.00] | 0.639 [0.00] | 0.000 [0.00] | 0.017 [0.93] | 0 | 0 |
| HLA-DRB4 | 6409 | 0.606 [0.10] | 0.382 [0.14] | 0.001 [0.29] | 0.011 [0.94] | 0 | 0 |
| HLA-DRB5 | 3584 | 0.183 [0.00] | 0.798 [0.02] | 0 | 0.019 [0.93] | 0 | 0 |
| HLA-H | 20083 | 0.920 [0.90] | 0.044 [0.92] | 0.001 [0.74] | 0.035 [0.95] | 0 | 0 |
| HLA-J | 23723 | 0.875 [0.17] | 0.081 [0.13] | 0.003 [0.23] | 0.040 [0.76] | 0 | 0.000 [0.00] |
| HLA-K | 19934 | 0.529 [0.83] | 0.047 [0.89] | 0.024 [0.79] | 0.400 [0.60] | 0 | 0.000 [1.00] |
| HLA-L | 23747 | 0.632 [1.00] | 0.226 [1.00] | 0.017 [1.00] | 0.125 [0.99] | 0 | 0 |
| HLA-N | 23845 | 0.789 [1.00] | 0.203 [1.00] | 0 | 0.008 [0.96] | 0 | 0.000 [0.00] |
| HLA-P | 23702 | 0.568 [0.00] | 0.171 [0.00] | 0.002 [0.00] | 0.169 [0.10] | 0.090 [0.01] | 0.000 [1.00] |
| HLA-S | 23682 | 0.813 [1.00] | 0.038 [1.00] | 0.098 [1.00] | 0.052 [0.99] | 0 | 0.000 [1.00] |
| HLA-T | 20067 | 0.549 [0.00] | 0.402 [0.00] | 0.001 [0.00] | 0.048 [0.28] | 0 | 0.000 [1.00] |
| HLA-U | 19962 | 0.274 [0.90] | 0.470 [1.00] | 0.001 [0.89] | 0.255 [0.91] | 0 | 0 |
| HLA-V | 23712 | 0.676 [1.00] | 0.077 [1.00] | 0.001 [1.00] | 0.245 [0.44] | 0 | 0.000 [1.00] |
| HLA-W | 23533 | 0.296 [0.47] | 0.163 [0.61] | 0.001 [0.38] | 0.540 [0.65] | 0.000 [1.00] | 0.000 [1.00] |
| HLA-Y | 3919 | 0.731 [0.96] | 0.225 [1.00] | 0.001 [1.00] | 0.043 [0.89] | 0 | 0 |
| MICA | 22540 | 0.619 [0.33] | 0.101 [0.14] | 0.233 [0.32] | 0.047 [0.84] | 0.001 [0.75] | 0.000 [1.00] |
| MICB | 22944 | 0.150 [0.00] | 0.014 [0.01] | 0.780 [0.01] | 0.056 [0.83] | 0.000 [0.67] | 0.000 [1.00] |
| TAP1 | 23427 | 0.382 [0.00] | 0.485 [0.00] | 0.006 [0.00] | 0.126 [0.37] | 0 | 0.000 [0.00] |
| TAP2 | 23557 | 0.416 [0.00] | 0.459 [0.00] | 0.015 [0.00] | 0.110 [0.37] | 0.000 [0.80] | 0 |
| C4A | 20174 | 0 | 0 | 0 | 0 | 0 | 1.000 [0.00] |
| C4B | 17796 | 0 | 0 | 0 | 0 | 0 | 1.000 [0.00] |
| HLA-HFE | 3068 | 0.139 [0.00] | 0.122 [0.00] | 0.722 [0.00] | 0.017 [0.34] | 0 | 0 |

## B. Sequence truth check for depth-2/3 calls

Rows are sequence classes per gene class. Counts are haplotype-gene calls.

| gene_class | cds_known | frameshift_or_stop | novel_cds_synonymous | novel_protein | protein_known |
|---|---|---|---|---|---|
| class_II_accessory | <20 | 1971 | 214 | 400 | 59 |
| class_II_paralog | 368 | 15044 | 80 | 463 | 70 |
| classical_I | 20 | 2043 | <20 | 147 | 0 |
| classical_II | <20 | 2196 | 29 | 160 | <20 |
| mic_tap | 22657 | 4035 | 587 | 3762 | 40 |
| nonclassical_I | 21 | 1913 | 163 | 294 | <20 |
| other | 2195 | <20 | 21 | 35 | 0 |
| pseudogene_I | 0 | 45447 | 0 | 27 | 0 |

Depth-3 calls whose protein is **not** catalogued (reclassified to `novel_protein`): <20. Depth-2 calls whose protein **is** catalogued (`protein_known`): 185.

## Clusters

Novel proteins are clustered by exact protein sequence, and synonymous novelty by exact CDS. `clean recurrent` means at least 2 unrelated people carry the cluster in a call with no artifact label.

| gene | cluster_type | n_clusters | n_all_clean | n_clean_recurrent_unrelated |
|---|---|---|---|---|
| HLA-A | novel_cds_synonymous | 8 | 8 | 0 |
| HLA-A | novel_protein | 34 | 24 | 0 |
| HLA-B | novel_cds_synonymous | 5 | 5 | 0 |
| HLA-B | novel_protein | 60 | 42 | 1 |
| HLA-C | novel_cds_synonymous | 3 | 3 | 0 |
| HLA-C | novel_protein | 51 | 30 | 1 |
| HLA-DMA | novel_cds_synonymous | 22 | 22 | 9 |
| HLA-DMA | novel_protein | 45 | 35 | 13 |
| HLA-DMB | novel_cds_synonymous | 18 | 18 | 6 |
| HLA-DMB | novel_protein | 39 | 34 | 11 |
| HLA-DOA | novel_cds_synonymous | 25 | 25 | 8 |
| HLA-DOA | novel_protein | 40 | 29 | 8 |
| HLA-DOB | novel_cds_synonymous | 14 | 14 | 3 |
| HLA-DOB | novel_protein | 32 | 28 | 6 |
| HLA-DPA1 | novel_cds_synonymous | 11 | 11 | 1 |
| HLA-DPA1 | novel_protein | 21 | 17 | 0 |
| HLA-DPB1 | novel_cds_synonymous | 8 | 8 | 0 |
| HLA-DPB1 | novel_protein | 32 | 21 | 2 |
| HLA-DPB2 | novel_protein | 14 | 3 | 0 |
| HLA-DQA1 | novel_cds_synonymous | 4 | 4 | 0 |
| HLA-DQA1 | novel_protein | 38 | 30 | 0 |
| HLA-DQA2 | novel_cds_synonymous | 17 | 17 | 8 |
| HLA-DQA2 | novel_protein | 38 | 33 | 12 |
| HLA-DQB1 | novel_cds_synonymous | 8 | 8 | 1 |
| HLA-DQB1 | novel_protein | 30 | 18 | 0 |
| HLA-DQB2 | novel_cds_synonymous | 23 | 23 | 8 |
| HLA-DQB2 | novel_protein | 57 | 42 | 18 |
| HLA-DRA | novel_cds_synonymous | 10 | 10 | 1 |
| HLA-DRA | novel_protein | 42 | 34 | 2 |
| HLA-DRB1 | novel_cds_synonymous | 1 | 1 | 0 |
| HLA-DRB1 | novel_protein | 35 | 22 | 2 |
| HLA-DRB3 | novel_cds_synonymous | 3 | 3 | 0 |
| HLA-DRB3 | novel_protein | 15 | 6 | 0 |
| HLA-DRB4 | novel_cds_synonymous | 3 | 3 | 2 |
| HLA-DRB4 | novel_protein | 5 | 3 | 0 |
| HLA-DRB5 | novel_protein | 5 | 2 | 1 |
| HLA-E | novel_cds_synonymous | 11 | 11 | 1 |
| HLA-E | novel_protein | 57 | 38 | 2 |
| HLA-F | novel_cds_synonymous | 18 | 18 | 2 |
| HLA-F | novel_protein | 74 | 60 | 12 |
| HLA-G | novel_cds_synonymous | 30 | 30 | 6 |
| HLA-G | novel_protein | 70 | 49 | 9 |
| HLA-HFE | novel_cds_synonymous | 7 | 7 | 3 |
| HLA-HFE | novel_protein | 10 | 10 | 5 |
| HLA-L | novel_protein | 4 | 1 | 0 |
| HLA-N | novel_protein | 5 | 1 | 0 |
| HLA-U | novel_protein | 2 | 1 | 0 |
| HLA-V | novel_protein | 1 | 1 | 0 |
| HLA-Y | novel_protein | 2 | 0 | 0 |
| MICA | novel_cds_synonymous | 20 | 20 | 5 |
| MICA | novel_protein | 95 | 63 | 15 |
| MICB | novel_cds_synonymous | 30 | 30 | 9 |
| MICB | novel_protein | 104 | 64 | 7 |
| TAP1 | novel_cds_synonymous | 44 | 44 | 17 |
| TAP1 | novel_protein | 146 | 123 | 51 |
| TAP2 | novel_cds_synonymous | 76 | 76 | 26 |
| TAP2 | novel_protein | 201 | 162 | 53 |

### Largest clean novel-protein clusters

| cluster_id | gene | n_persons_unrelated_clean | nearest_known_protein | nearest_distance | aa_diffs | in_groove | groove_method |
|---|---|---|---|---|---|---|---|
| TAP1_prot_b531d1b8 | TAP1 | 426 | TAP1*01:01 | 1 | 419:G>C | NA | NA |
| TAP1_prot_9f152af0 | TAP1 | 321 | TAP1*05:01 | 1 | 131:P>L | NA | NA |
| TAP2_prot_e11a5381 | TAP2 | 310 | TAP2*02:01 | 1 | 577:M>V | NA | NA |
| TAP2_prot_17ff0ece | TAP2 | 201 | TAP2*01:01 | 2 | 374:A>T,467:V>I | NA,NA | NA |
| TAP2_prot_1dd29f06 | TAP2 | 192 | TAP2*01:01 | 1 | 686:->Q | NA | NA |
| TAP1_prot_193c898b | TAP1 | 191 | TAP1*04:01 | 1 | 244:V>L | NA | NA |
| TAP1_prot_4292c7c4 | TAP1 | 189 | TAP1*01:01 | 1 | 286:S>F | NA | NA |
| TAP1_prot_ade86b09 | TAP1 | 98 | TAP1*05:01 | 2 | 17:G>R,131:P>L | NA,NA | NA |
| TAP1_prot_134b3404 | TAP1 | 94 | TAP1*02:01 | 2 | 17:G>R,728:Q>K | NA,NA | NA |
| TAP2_prot_82da8490 | TAP2 | 94 | TAP2*01:01 | 1 | 15:V>A | NA | NA |
| MICB_prot_9f73eb8c | MICB | 81 | MICB*001 | 1 | 39:E>G | NA | NA |
| DQB2_prot_4ebace1e | HLA-DQB2 | 80 | HLA-DQB2*01:01 | 1 | 198:R>C | False | refdata_exon |
| DQB2_prot_e2bae253 | HLA-DQB2 | 62 | HLA-DQB2*01:08 | 1 | 1:T>M | False | refdata_exon |
| TAP2_prot_a4570c35 | TAP2 | 56 | TAP2*02:01 | 1 | 328:A>G | NA | NA |
| TAP2_prot_ca41dea3 | TAP2 | 53 | TAP2*01:01 | 1 | 374:A>T | NA | NA |

## Groove annotation

A codon counts as in the groove if its middle base falls in exon 2/3 (class I: A, B, C, E, F, G) or exon 2 (class II chains: DRA, DRB*, DQA*, DQB*, DPA*, DPB*). `refdata_exon` uses the nearest allele's own CDS=/exon= record from `alleles.csv.gz`. `refdata_gene_mode` uses that gene's most common coding-exon pattern across its goodTemplate alleles. `approx_canonical` is a hard-coded HLA-A/B/C map (73/270/276 nt) and is used only when the gene is missing from alleles.csv.gz. `NA` means none applied, with the reason below.

| gene | groove_method | groove_na_reason | n_clusters |
|---|---|---|---|
| HLA-A | refdata_exon |  | 34 |
| HLA-B | refdata_exon |  | 60 |
| HLA-C | refdata_exon |  | 51 |
| HLA-DMA | NA | no peptide-groove exon definition for DMA | 45 |
| HLA-DMB | NA | no peptide-groove exon definition for DMB | 39 |
| HLA-DOA | NA | no peptide-groove exon definition for DOA | 40 |
| HLA-DOB | NA | no peptide-groove exon definition for DOB | 32 |
| HLA-DPA1 | refdata_exon |  | 21 |
| HLA-DPB1 | refdata_exon |  | 32 |
| HLA-DPB2 | NA | refdata_gene_mode: CDS segment total != reference CDS length | 14 |
| HLA-DQA1 | refdata_exon |  | 38 |
| HLA-DQA2 | refdata_exon |  | 38 |
| HLA-DQB1 | refdata_exon |  | 30 |
| HLA-DQB2 | refdata_exon |  | 57 |
| HLA-DRA | refdata_exon |  | 42 |
| HLA-DRB1 | refdata_exon |  | 35 |
| HLA-DRB3 | refdata_exon |  | 15 |
| HLA-DRB4 | refdata_exon |  | 5 |
| HLA-DRB5 | refdata_exon |  | 5 |
| HLA-E | refdata_exon |  | 57 |
| HLA-F | refdata_exon |  | 74 |
| HLA-G | refdata_exon |  | 70 |
| HLA-HFE | NA | no peptide-groove exon definition for HFE | 10 |
| HLA-L | NA | no peptide-groove exon definition for L | 4 |
| HLA-N | NA | no peptide-groove exon definition for N | 5 |
| HLA-U | NA | no peptide-groove exon definition for U | 2 |
| HLA-V | NA | no peptide-groove exon definition for V | 1 |
| HLA-Y | NA | no peptide-groove exon definition for Y | 2 |
| MICA | NA | no peptide-groove exon definition for MICA | 95 |
| MICB | NA | no peptide-groove exon definition for MICB | 104 |
| TAP1 | NA | no peptide-groove exon definition for TAP1 | 146 |
| TAP2 | NA | no peptide-groove exon definition for TAP2 | 201 |

## Release drift

Not run (`--imgt-latest-nuc` not given).

## Matching and dropped rows

| metric | count |
|---|---|
| n_ambiguous_copy | 58 |
| n_copy1_duplicate_keys_dropped | 0 |
| n_join_ambiguous | 3715 |
| n_matched | 108247 |
| n_missing_cds | <20 |
| n_novel_table1_rows | 108306 |
| n_seq_calls_joined | 104532 |

Rows left out of the WS3 allele-count tables (uncalled, depth 1, unmatched, `no_refdata`; `clean_only` also drops flagged calls):

| table | n_calls_dropped |
|---|---|
| cds_all | 41892 |
| cds_clean_only | 40690 |
| protein_all | 41892 |
| protein_clean_only | 40690 |

## Assumptions (verify on the VM)

- Refdata: 40 CDSseq files, 38896 records (21423 ending in a stop codon, 15590 with frame!=1 kept out of the protein index, 0 unparseable headers). Headers were parsed as the first `*`-bearing token, with an optional `HLA-` prefix.
- IPD-IMGT/HLA version (alleles.csv.gz `version=`): IPD-IMGT/HLA-V3.55.0.
- Exon maps from alleles.csv.gz: 286 nearest reference alleles resolved by name; gene-mode patterns for A, B, C, DMA, DMB, DOA, DOB, DPA1, DPB1, DPB2, DQA1, DQA2, DQB1, DQB2, DRA, DRB1, DRB3, DRB4, DRB5, E, F, G, HFE, L, MICA, MICB, N, TAP1, TAP2, U, V, Y.
- Relatedness ids are assumed to equal Table 1 person_id strings.
