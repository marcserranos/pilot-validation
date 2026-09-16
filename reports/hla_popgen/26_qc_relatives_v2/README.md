# Are the long-read HLA calls actually right? (QC v2)

This is the load-bearing QC for the paper's central claim. Version 1 of this analysis (scripts 16/17) answered a question it had partly assumed the answer to; this version fixes each of those problems. Everything below is computed on the **exact assembled DNA sequence** of each gene, not on allele names, except where it says otherwise.

## 1. Which pairs of people we looked at, and why that choice is now fair

We took **every** pair of people in the long-read cohort whom All of Us's own genome-wide relatedness table calls first-degree relatives or closer (kinship >= 0.177). That is **574** pairs, of which **<20** are duplicates or identical twins (kinship > 0.354).

Version 1 instead kept only pairs that already agreed at >=80% of their genes, and then reported how often those same pairs disagreed. That is circular: the pairs were chosen using the very number being reported. Nothing here is chosen using HLA data at all.

## 2. How much HLA the relatives share (a distribution, not a filter)

Two first-degree relatives do not have to share HLA. A parent and child always share one copy, but two siblings can easily inherit different copies from each parent and share none. So the honest picture is the whole spread, not a cutoff:

- Median fraction of genes where the pair shares an exact sequence: **1.000**

- Pairs sharing at almost every gene (>=0.95): **356**; pairs sharing at almost none (<=0.1): **<20**. Two modes here is the expected biology, not a problem.


Because any cutoff is a judgement call, the headline disagreement rate is reported at several cutoffs. If the number moves a lot across this table, the cutoff was doing the work rather than the data:

| threshold | n_pairs | n_gene_comparisons | n_discordant | discordance_rate_% |
|---|---|---|---|---|
| 0.5 | 566 | 21241 | 2607 | 12.273 |
| 0.8 | 532 | 19970 | 2176 | 10.896 |
| 0.9 | 483 | 18150 | 1783 | 9.824 |
| 0.95 | 356 | 13390 | 1126 | 8.409 |


**Age-gap stratification was SKIPPED**: no `--yob-tsv` was supplied, so we cannot separate probable parent-child pairs (which must share a copy at every gene) from siblings (which need not). Supply `--yob-tsv` with columns `person_id`, `year_of_birth` to enable it.


## 3. What the disagreements actually are

Every single gene comparison is placed in exactly one box. Crucially, a gene that one person has and the other does not is now COUNTED (as `missing_one`) instead of being quietly deleted from the denominator, which is what hid allele dropout and assembly gaps in version 1.

- `concordant` - both people have at least one byte-identical copy of this gene.
- `missing_one` / `missing_both` - the gene was not assembled in one / either person. A coverage result, not an accuracy result.
- `point_diff` - up to 3 differing bases in one place. This is what ordinary long-read consensus error looks like.
- `dropout_candidate` - a real difference, and at least one of the two people looks homozygous here, which is what a collapsed (lost) second copy looks like.
- `fragmentation_candidate` - a real difference at a gene that sits alone on its own assembly fragment, where there is nothing to anchor it.
- `other_multiblock` - a real, scattered difference with none of the above excuses. This is the bucket that would genuinely worry us.


| gene_class | n_total | n_comparable | n_discordant | discordance_rate_% | concordant | missing_one | missing_both | point_diff | dropout_candidate | fragmentation_candidate | other_multiblock |
|---|---|---|---|---|---|---|---|---|---|---|---|
| class_II_accessory | 2870 | 2869 | 126 | 4.392 | 2743 | 56 | 1 | 43 | 21 | 0 | 6 |
| class_II_paralog | 4018 | 3286 | 652 | 19.842 | 2634 | 479 | 732 | 72 | 90 | 0 | 11 |
| classical_I | 1722 | 1722 | 214 | 12.427 | 1508 | 65 | 0 | 30 | 63 | 33 | 23 |
| classical_II | 2870 | 2865 | 278 | 9.703 | 2587 | 97 | 5 | 45 | 84 | 6 | 46 |
| complement | 1148 | 0 | 0 | nan | 0 | 0 | 1148 | 0 | 0 | 0 | 0 |
| mic_tap | 2296 | 2295 | 296 | 12.898 | 1999 | 55 | 1 | 102 | 98 | 0 | 41 |
| nonclassical_I | 1722 | 1722 | 156 | 9.059 | 1566 | 56 | 0 | 61 | 34 | 0 | 5 |
| other | 574 | 241 | 197 | 81.743 | 44 | 188 | 333 | 9 | 0 | 0 | 0 |
| pseudogene_I | 6888 | 6545 | 877 | 13.4 | 5668 | 397 | 343 | 129 | 307 | 0 | 44 |


Across all genes, **12.98%** of comparable comparisons disagree (95% CI 12.54-13.43%), and **17.6%** of those disagreements are just 1-3 bases.


## 4. The cleanest possible check: the same person sequenced twice

Duplicate/identical-twin pairs should be byte-identical. Comparing both haplotypes gene by gene gives **121 differing bases in 306033 compared bases**, i.e. a per-base quality of **QV 34.0** (higher is better; 40-50 is the normal range for this kind of assembly). This number involves no relatives, no inheritance assumptions and no allele names - it is the most direct accuracy statement in this report.


## 5. Are we silently losing one copy of a gene? (homozygosity check)

If an assembly ever collapses two different copies of a gene into one, people will look homozygous (two identical copies) more often than genetics allows. We compare how often people ARE homozygous with how often they SHOULD be, given the allele frequencies of their own ancestry group. A ratio near 1.0 is healthy; clearly above 1.0 means copies are going missing.

We do this three ways on the SAME people - long-read exact sequence, long-read allele names, and All of Us's independent short-read allele names - so the comparison is like for like. Only people whose ancestry is unambiguous (own admixture proportion >= 0.9) are used, because the expectation depends on which population's frequencies apply.

Expected homozygosity uses the finite-sample-corrected formula `sum x(x-1)/(n(n-1))`; the naive version is biased upward in small groups, which would have hidden exactly this signal where it matters most.


| measure | gene | ancestry | n_people | obs_hom | exp_hom | obs_exp_ratio | ci_lo | ci_hi | thin_n |
|---|---|---|---|---|---|---|---|---|---|
| LR_2field | A | AFR | 2055 | 0.072 | 0.06 | 1.194 | 0.98 | 1.389 | False |
| LR_2field | A | AMR | 1511 | 0.099 | 0.084 | 1.182 | 1.02 | 1.347 | False |
| LR_2field | A | EAS | 1150 | 0.138 | 0.109 | 1.271 | 1.101 | 1.455 | False |
| LR_2field | A | EUR | 1216 | 0.142 | 0.14 | 1.013 | 0.889 | 1.156 | False |
| LR_2field | A | MID | 277 | 0.097 | 0.074 | 1.309 | 0.864 | 1.762 | False |
| LR_2field | A | SAS | 1067 | 0.108 | 0.09 | 1.191 | 0.972 | 1.403 | False |
| LR_2field | B | AFR | 2077 | 0.052 | 0.048 | 1.1 | 0.92 | 1.287 | False |
| LR_2field | B | AMR | 1483 | 0.046 | 0.03 | 1.533 | 1.201 | 1.854 | False |
| LR_2field | B | EAS | 1123 | 0.063 | 0.042 | 1.52 | 1.174 | 1.819 | False |
| LR_2field | B | EUR | 1217 | 0.063 | 0.062 | 1.018 | 0.815 | 1.213 | False |
| LR_2field | B | MID | 281 | 0.053 | 0.036 | 1.472 | 0.763 | 2.063 | False |
| LR_2field | B | SAS | 1050 | 0.057 | 0.05 | 1.134 | 0.799 | 1.455 | False |
| LR_2field | C | AFR | 2112 | 0.092 | 0.082 | 1.125 | 0.98 | 1.277 | False |
| LR_2field | C | AMR | 1466 | 0.088 | 0.076 | 1.155 | 0.976 | 1.333 | False |
| LR_2field | C | EAS | 1117 | 0.099 | 0.092 | 1.082 | 0.914 | 1.278 | False |
| LR_2field | C | EUR | 1200 | 0.096 | 0.09 | 1.062 | 0.886 | 1.207 | False |
| LR_2field | C | MID | 270 | 0.107 | 0.086 | 1.246 | 0.885 | 1.622 | False |
| LR_2field | C | SAS | 1026 | 0.085 | 0.078 | 1.086 | 0.87 | 1.301 | False |
| LR_2field | DPA1 | AFR | 2168 | 0.286 | 0.277 | 1.033 | 0.961 | 1.105 | False |
| LR_2field | DPA1 | AMR | 1519 | 0.623 | 0.599 | 1.041 | 1.012 | 1.074 | False |
| LR_2field | DPA1 | EAS | 1194 | 0.41 | 0.395 | 1.038 | 0.978 | 1.093 | False |
| LR_2field | DPA1 | EUR | 1235 | 0.684 | 0.683 | 1.002 | 0.979 | 1.022 | False |
| LR_2field | DPA1 | MID | 286 | 0.64 | 0.622 | 1.028 | 0.969 | 1.087 | False |
| LR_2field | DPA1 | SAS | 1071 | 0.503 | 0.497 | 1.012 | 0.964 | 1.052 | False |
| LR_2field | DPB1 | AFR | 2097 | 0.159 | 0.143 | 1.112 | 1.012 | 1.206 | False |
| LR_2field | DPB1 | AMR | 1465 | 0.224 | 0.163 | 1.372 | 1.234 | 1.473 | False |
| LR_2field | DPB1 | EAS | 1149 | 0.185 | 0.156 | 1.183 | 1.039 | 1.3 | False |
| LR_2field | DPB1 | EUR | 1188 | 0.221 | 0.216 | 1.023 | 0.932 | 1.122 | False |
| LR_2field | DPB1 | MID | 279 | 0.201 | 0.178 | 1.127 | 0.953 | 1.354 | False |
| LR_2field | DPB1 | SAS | 1047 | 0.203 | 0.196 | 1.037 | 0.924 | 1.131 | False |
| LR_2field | DQA1 | AFR | 2122 | 0.156 | 0.148 | 1.049 | 0.959 | 1.111 | False |
| LR_2field | DQA1 | AMR | 1531 | 0.135 | 0.108 | 1.245 | 1.099 | 1.385 | False |
| LR_2field | DQA1 | EAS | 1189 | 0.129 | 0.101 | 1.278 | 1.114 | 1.435 | False |
| LR_2field | DQA1 | EUR | 1226 | 0.134 | 0.121 | 1.106 | 0.966 | 1.258 | False |
| LR_2field | DQA1 | MID | 282 | 0.167 | 0.125 | 1.329 | 0.98 | 1.632 | False |
| LR_2field | DQA1 | SAS | 1081 | 0.128 | 0.12 | 1.068 | 0.912 | 1.259 | False |
| LR_2field | DQB1 | AFR | 2079 | 0.129 | 0.115 | 1.123 | 0.981 | 1.254 | False |
| LR_2field | DQB1 | AMR | 1493 | 0.151 | 0.129 | 1.172 | 1.033 | 1.308 | False |
| LR_2field | DQB1 | EAS | 1179 | 0.137 | 0.116 | 1.174 | 1.015 | 1.347 | False |
| LR_2field | DQB1 | EUR | 1223 | 0.128 | 0.111 | 1.144 | 0.962 | 1.29 | False |
| LR_2field | DQB1 | MID | 277 | 0.148 | 0.122 | 1.21 | 0.905 | 1.507 | False |
| LR_2field | DQB1 | SAS | 1069 | 0.121 | 0.099 | 1.22 | 1.012 | 1.424 | False |
| LR_2field | DRB1 | AFR | 2083 | 0.068 | 0.064 | 1.06 | 0.883 | 1.21 | False |
| LR_2field | DRB1 | AMR | 1464 | 0.074 | 0.043 | 1.731 | 1.379 | 2.047 | False |
| LR_2field | DRB1 | EAS | 1140 | 0.101 | 0.072 | 1.401 | 1.19 | 1.63 | False |
| LR_2field | DRB1 | EUR | 1180 | 0.088 | 0.081 | 1.088 | 0.895 | 1.271 | False |
| LR_2field | DRB1 | MID | 271 | 0.111 | 0.065 | 1.712 | 1.166 | 2.131 | False |
| LR_2field | DRB1 | SAS | 1037 | 0.102 | 0.083 | 1.238 | 1.06 | 1.448 | False |
| LR_CDS_exact | A | AFR | 2196 | 0.064 | 0.056 | 1.143 | 0.981 | 1.308 | False |
| LR_CDS_exact | A | AMR | 1566 | 0.092 | 0.08 | 1.154 | 0.979 | 1.328 | False |
| LR_CDS_exact | A | EAS | 1214 | 0.127 | 0.101 | 1.251 | 1.096 | 1.436 | False |
| LR_CDS_exact | A | EUR | 1284 | 0.129 | 0.132 | 0.977 | 0.854 | 1.107 | False |
| LR_CDS_exact | A | MID | 291 | 0.093 | 0.072 | 1.297 | 0.815 | 1.722 | False |
| LR_CDS_exact | A | SAS | 1112 | 0.102 | 0.084 | 1.216 | 1.013 | 1.411 | False |
| LR_CDS_exact | B | AFR | 2231 | 0.045 | 0.044 | 1.016 | 0.851 | 1.183 | False |
| LR_CDS_exact | B | AMR | 1555 | 0.037 | 0.027 | 1.34 | 1.019 | 1.644 | False |
| LR_CDS_exact | B | EAS | 1211 | 0.055 | 0.037 | 1.476 | 1.135 | 1.813 | False |
| LR_CDS_exact | B | EUR | 1274 | 0.053 | 0.059 | 0.899 | 0.709 | 1.095 | False |
| LR_CDS_exact | B | MID | 294 | 0.051 | 0.034 | 1.505 | 0.803 | 2.3 | False |
| LR_CDS_exact | B | SAS | 1101 | 0.047 | 0.045 | 1.043 | 0.775 | 1.317 | False |
| LR_CDS_exact | C | AFR | 2253 | 0.08 | 0.072 | 1.103 | 0.967 | 1.234 | False |
| LR_CDS_exact | C | AMR | 1556 | 0.076 | 0.069 | 1.104 | 0.904 | 1.267 | False |
| LR_CDS_exact | C | EAS | 1205 | 0.09 | 0.084 | 1.072 | 0.883 | 1.219 | False |
| LR_CDS_exact | C | EUR | 1262 | 0.078 | 0.083 | 0.945 | 0.777 | 1.141 | False |
| LR_CDS_exact | C | MID | 282 | 0.089 | 0.076 | 1.166 | 0.784 | 1.551 | False |
| LR_CDS_exact | C | SAS | 1089 | 0.075 | 0.072 | 1.05 | 0.833 | 1.252 | False |
| LR_CDS_exact | DPA1 | AFR | 2260 | 0.213 | 0.21 | 1.014 | 0.939 | 1.087 | False |
| LR_CDS_exact | DPA1 | AMR | 1553 | 0.588 | 0.568 | 1.036 | 1.013 | 1.065 | False |
| LR_CDS_exact | DPA1 | EAS | 1227 | 0.394 | 0.382 | 1.032 | 0.967 | 1.098 | False |
| LR_CDS_exact | DPA1 | EUR | 1289 | 0.642 | 0.642 | 1.0 | 0.98 | 1.019 | False |
| LR_CDS_exact | DPA1 | MID | 295 | 0.607 | 0.591 | 1.027 | 0.963 | 1.08 | False |
| LR_CDS_exact | DPA1 | SAS | 1098 | 0.474 | 0.47 | 1.009 | 0.961 | 1.055 | False |
| LR_CDS_exact | DPB1 | AFR | 2222 | 0.121 | 0.113 | 1.07 | 0.96 | 1.184 | False |
| LR_CDS_exact | DPB1 | AMR | 1531 | 0.206 | 0.153 | 1.349 | 1.239 | 1.454 | False |
| LR_CDS_exact | DPB1 | EAS | 1204 | 0.174 | 0.148 | 1.172 | 1.031 | 1.301 | False |
| LR_CDS_exact | DPB1 | EUR | 1264 | 0.202 | 0.199 | 1.012 | 0.927 | 1.098 | False |
| LR_CDS_exact | DPB1 | MID | 291 | 0.192 | 0.173 | 1.115 | 0.892 | 1.321 | False |
| LR_CDS_exact | DPB1 | SAS | 1081 | 0.194 | 0.187 | 1.038 | 0.94 | 1.148 | False |
| LR_CDS_exact | DQA1 | AFR | 2206 | 0.123 | 0.123 | 1.006 | 0.897 | 1.119 | False |
| LR_CDS_exact | DQA1 | AMR | 1560 | 0.121 | 0.098 | 1.227 | 1.079 | 1.404 | False |
| LR_CDS_exact | DQA1 | EAS | 1224 | 0.115 | 0.088 | 1.315 | 1.11 | 1.517 | False |
| LR_CDS_exact | DQA1 | EUR | 1256 | 0.115 | 0.11 | 1.039 | 0.894 | 1.167 | False |
| LR_CDS_exact | DQA1 | MID | 287 | 0.153 | 0.117 | 1.307 | 1.002 | 1.614 | False |
| LR_CDS_exact | DQA1 | SAS | 1110 | 0.118 | 0.111 | 1.068 | 0.894 | 1.228 | False |
| LR_CDS_exact | DQB1 | AFR | 2198 | 0.112 | 0.104 | 1.071 | 0.955 | 1.222 | False |
| LR_CDS_exact | DQB1 | AMR | 1545 | 0.137 | 0.122 | 1.122 | 0.978 | 1.268 | False |
| LR_CDS_exact | DQB1 | EAS | 1218 | 0.128 | 0.109 | 1.18 | 1.008 | 1.341 | False |
| LR_CDS_exact | DQB1 | EUR | 1260 | 0.118 | 0.108 | 1.098 | 0.954 | 1.257 | False |
| LR_CDS_exact | DQB1 | MID | 288 | 0.139 | 0.112 | 1.24 | 0.889 | 1.486 | False |
| LR_CDS_exact | DQB1 | SAS | 1102 | 0.112 | 0.094 | 1.19 | 1.011 | 1.421 | False |
| LR_CDS_exact | DRB1 | AFR | 2174 | 0.06 | 0.06 | 1.009 | 0.827 | 1.126 | False |
| LR_CDS_exact | DRB1 | AMR | 1515 | 0.066 | 0.04 | 1.641 | 1.362 | 1.891 | False |
| LR_CDS_exact | DRB1 | EAS | 1196 | 0.095 | 0.068 | 1.401 | 1.169 | 1.618 | False |
| LR_CDS_exact | DRB1 | EUR | 1243 | 0.083 | 0.077 | 1.083 | 0.879 | 1.276 | False |
| LR_CDS_exact | DRB1 | MID | 281 | 0.107 | 0.059 | 1.807 | 1.17 | 2.531 | False |
| LR_CDS_exact | DRB1 | SAS | 1066 | 0.09 | 0.076 | 1.184 | 0.945 | 1.365 | False |
| SR_2field | A | AFR | 2401 | 0.066 | 0.064 | 1.021 | 0.899 | 1.183 | False |
| SR_2field | A | AMR | 1677 | 0.089 | 0.085 | 1.047 | 0.887 | 1.189 | False |
| SR_2field | A | EAS | 1303 | 0.128 | 0.11 | 1.17 | 1.036 | 1.321 | False |
| SR_2field | A | EUR | 1377 | 0.125 | 0.14 | 0.895 | 0.777 | 0.999 | False |
| SR_2field | A | MID | 310 | 0.097 | 0.077 | 1.254 | 0.895 | 1.686 | False |
| SR_2field | A | SAS | 1169 | 0.105 | 0.091 | 1.153 | 0.978 | 1.337 | False |
| SR_2field | B | AFR | 2401 | 0.051 | 0.049 | 1.039 | 0.869 | 1.202 | False |
| SR_2field | B | AMR | 1677 | 0.039 | 0.03 | 1.332 | 0.999 | 1.665 | False |
| SR_2field | B | EAS | 1303 | 0.063 | 0.042 | 1.497 | 1.215 | 1.847 | False |
| SR_2field | B | EUR | 1377 | 0.067 | 0.064 | 1.047 | 0.84 | 1.249 | False |
| SR_2field | B | MID | 310 | 0.071 | 0.036 | 1.992 | 1.203 | 2.637 | False |
| SR_2field | B | SAS | 1169 | 0.063 | 0.051 | 1.239 | 0.979 | 1.551 | False |
| SR_2field | C | AFR | 2401 | 0.096 | 0.091 | 1.057 | 0.931 | 1.175 | False |
| SR_2field | C | AMR | 1677 | 0.083 | 0.079 | 1.057 | 0.906 | 1.201 | False |
| SR_2field | C | EAS | 1303 | 0.114 | 0.095 | 1.203 | 1.046 | 1.39 | False |
| SR_2field | C | EUR | 1377 | 0.099 | 0.092 | 1.071 | 0.919 | 1.21 | False |
| SR_2field | C | MID | 310 | 0.152 | 0.091 | 1.66 | 1.244 | 2.012 | False |
| SR_2field | C | SAS | 1169 | 0.103 | 0.085 | 1.209 | 1.005 | 1.39 | False |
| SR_2field | DPA1 | AFR | 2401 | 0.292 | 0.283 | 1.032 | 0.977 | 1.096 | False |
| SR_2field | DPA1 | AMR | 1677 | 0.62 | 0.595 | 1.043 | 1.015 | 1.068 | False |
| SR_2field | DPA1 | EAS | 1303 | 0.402 | 0.395 | 1.018 | 0.96 | 1.072 | False |
| SR_2field | DPA1 | EUR | 1377 | 0.675 | 0.679 | 0.994 | 0.976 | 1.012 | False |
| SR_2field | DPA1 | MID | 310 | 0.645 | 0.63 | 1.023 | 0.976 | 1.072 | False |
| SR_2field | DPA1 | SAS | 1169 | 0.505 | 0.499 | 1.011 | 0.963 | 1.056 | False |
| SR_2field | DPB1 | AFR | 2401 | 0.161 | 0.155 | 1.037 | 0.951 | 1.127 | False |
| SR_2field | DPB1 | AMR | 1677 | 0.23 | 0.181 | 1.272 | 1.17 | 1.367 | False |
| SR_2field | DPB1 | EAS | 1303 | 0.2 | 0.18 | 1.112 | 1.016 | 1.232 | False |
| SR_2field | DPB1 | EUR | 1377 | 0.223 | 0.222 | 1.004 | 0.936 | 1.081 | False |
| SR_2field | DPB1 | MID | 310 | 0.245 | 0.193 | 1.27 | 1.06 | 1.475 | False |
| SR_2field | DPB1 | SAS | 1169 | 0.222 | 0.207 | 1.074 | 0.97 | 1.16 | False |
| SR_2field | DQA1 | AFR | 2401 | 0.178 | 0.181 | 0.986 | 0.916 | 1.066 | False |
| SR_2field | DQA1 | AMR | 1677 | 0.192 | 0.181 | 1.061 | 0.955 | 1.164 | False |
| SR_2field | DQA1 | EAS | 1303 | 0.186 | 0.169 | 1.105 | 0.987 | 1.215 | False |
| SR_2field | DQA1 | EUR | 1377 | 0.182 | 0.178 | 1.017 | 0.91 | 1.129 | False |
| SR_2field | DQA1 | MID | 310 | 0.239 | 0.214 | 1.118 | 0.94 | 1.363 | False |
| SR_2field | DQA1 | SAS | 1169 | 0.176 | 0.163 | 1.081 | 0.94 | 1.21 | False |
| SR_2field | DQB1 | AFR | 2401 | 0.153 | 0.153 | 0.996 | 0.91 | 1.078 | False |
| SR_2field | DQB1 | AMR | 1677 | 0.158 | 0.149 | 1.059 | 0.946 | 1.158 | False |
| SR_2field | DQB1 | EAS | 1303 | 0.126 | 0.118 | 1.063 | 0.887 | 1.179 | False |
| SR_2field | DQB1 | EUR | 1377 | 0.145 | 0.138 | 1.05 | 0.933 | 1.176 | False |
| SR_2field | DQB1 | MID | 310 | 0.19 | 0.152 | 1.251 | 0.99 | 1.523 | False |
| SR_2field | DQB1 | SAS | 1169 | 0.136 | 0.118 | 1.155 | 1.001 | 1.295 | False |
| SR_2field | DRB1 | AFR | 2401 | 0.057 | 0.065 | 0.882 | 0.768 | 1.022 | False |
| SR_2field | DRB1 | AMR | 1677 | 0.063 | 0.044 | 1.423 | 1.164 | 1.671 | False |
| SR_2field | DRB1 | EAS | 1303 | 0.092 | 0.072 | 1.283 | 1.07 | 1.456 | False |
| SR_2field | DRB1 | EUR | 1377 | 0.084 | 0.082 | 1.014 | 0.842 | 1.203 | False |
| SR_2field | DRB1 | MID | 310 | 0.126 | 0.065 | 1.931 | 1.408 | 2.408 | False |
| SR_2field | DRB1 | SAS | 1169 | 0.104 | 0.084 | 1.236 | 1.069 | 1.429 | False |


Groups with fewer than 50 people are shown but flagged `thin_n`; treat them as suggestive only.


Split by how fragmented the person's assembly is (how many separate contigs carry their 8 classical genes):

| gene | fragmentation_class | n_people | obs_hom |
|---|---|---|---|
| A | 1 | 1056 | 0.113 |
| A | 2 | 700 | 0.113 |
| A | >=3 | 5907 | 0.092 |
| B | 1 | 1058 | 0.051 |
| B | 2 | 691 | 0.046 |
| B | >=3 | 5917 | 0.046 |
| C | 1 | 1061 | 0.087 |
| C | 2 | 692 | 0.079 |
| C | >=3 | 5894 | 0.079 |
| DPA1 | 1 | 1067 | 0.516 |
| DPA1 | 2 | 724 | 0.477 |
| DPA1 | >=3 | 5931 | 0.423 |
| DPB1 | 1 | 1063 | 0.193 |
| DPB1 | 2 | 715 | 0.19 |
| DPB1 | >=3 | 5815 | 0.167 |
| DQA1 | 1 | 1063 | 0.135 |
| DQA1 | 2 | 729 | 0.123 |
| DQA1 | >=3 | 5851 | 0.117 |
| DQB1 | 1 | 1062 | 0.144 |
| DQB1 | 2 | 728 | 0.128 |
| DQB1 | >=3 | 5821 | 0.117 |
| DRB1 | 1 | 1064 | 0.094 |
| DRB1 | 2 | 715 | 0.076 |
| DRB1 | >=3 | 5696 | 0.074 |


## 6. Does the phase-switch test have any power? (the missing denominator)

Version 1's headline was '0 of 545 pairs show a phase switch'. That is only meaningful if a switch COULD have been seen. A switch is only visible between two neighbouring genes that (a) both carry unambiguous shared-copy information and (b) sit on the same assembly fragment in both people. We now count those opportunities:

- Median genes carrying phase information per pair: **15.0**

- Median *testable* neighbouring-gene transitions per pair: **5.0**; total across all pairs: **3021**

- Pairs where NO switch could ever have been detected (zero testable transitions): **37** of 569

- Pairs where a switch WAS observed: **0**


Then we measured the test instead of trusting it: in each eligible pair we deliberately broke the phasing of one person at a random point inside one contig (a synthetic switch error, with a fixed random seed so this is reproducible) and asked the same detector to find it. It found **100.0%** of them across 532 eligible pairs.

Read this as the ceiling on what '0 switches' can mean: if the detection rate is low, a clean result is mostly a statement about power, not about the assemblies.


## 7. An independent opinion: All of Us's own short-read calls

Steps 1-6 all use relatives. This step does not. For every person who has both kinds of data we ask how many of their two alleles per gene agree with All of Us's independent short-read HLA calls, at two-field resolution (both sides normalized the same way).

The sharp prediction: a long-read call that is novel only OUTSIDE the coding sequence (`f4_noncoding`) is invisible at two fields and must NOT agree less often. A novel call that changes the protein (`f2_protein`) is a two-field difference by construction and SHOULD agree less often. If `f4` novelty also lowered agreement, our novelty calls would look like noise.


| gene | ancestry | field_class | n_person_genes | mean_alleles_matching | frac_both_match | frac_zero_match | thin_n |
|---|---|---|---|---|---|---|---|
| A | AFR | f2_protein | 132 | 0.985 | 0.0 | 0.015 | False |
| A | AFR | f4_noncoding | 149 | 1.832 | 0.832 | 0.0 | False |
| A | AFR | known | 2094 | 1.862 | 0.866 | 0.003 | False |
| A | AMR | f2_protein | 55 | 0.982 | 0.0 | 0.018 | False |
| A | AMR | f3_synonymous | <20 | 2.0 | 1.0 | 0.0 | True |
| A | AMR | f4_noncoding | 70 | 1.829 | 0.829 | 0.0 | False |
| A | AMR | known | 1541 | 1.907 | 0.91 | 0.004 | False |
| A | EAS | f2_protein | 59 | 1.0 | 0.0 | 0.0 | False |
| A | EAS | f4_noncoding | 71 | 1.873 | 0.887 | 0.014 | False |
| A | EAS | known | 1167 | 1.896 | 0.897 | 0.001 | False |
| A | EUR | f2_protein | 68 | 1.0 | 0.0 | 0.0 | False |
| A | EUR | f3_synonymous | <20 | 2.0 | 1.0 | 0.0 | True |
| A | EUR | f4_noncoding | 71 | 1.873 | 0.873 | 0.0 | False |
| A | EUR | known | 1228 | 1.919 | 0.919 | 0.0 | False |
| A | MID | f2_protein | <20 | 1.0 | 0.0 | 0.0 | True |
| A | MID | f4_noncoding | <20 | 2.0 | 1.0 | 0.0 | True |
| A | MID | known | 278 | 1.906 | 0.906 | 0.0 | False |
| A | SAS | f2_protein | 45 | 0.978 | 0.0 | 0.022 | True |
| A | SAS | f3_synonymous | <20 | 2.0 | 1.0 | 0.0 | True |
| A | SAS | f4_noncoding | 56 | 1.857 | 0.857 | 0.0 | False |
| A | SAS | known | 1058 | 1.928 | 0.933 | 0.005 | False |
| B | AFR | f2_protein | 153 | 0.987 | 0.0 | 0.013 | False |
| B | AFR | f4_noncoding | 472 | 1.915 | 0.915 | 0.0 | False |
| B | AFR | known | 1761 | 1.892 | 0.894 | 0.002 | False |
| B | AMR | f2_protein | 77 | 0.974 | 0.0 | 0.026 | False |
| B | AMR | f3_synonymous | <20 | 1.0 | 0.5 | 0.5 | True |
| B | AMR | f4_noncoding | 299 | 1.936 | 0.943 | 0.007 | False |
| B | AMR | known | 1290 | 1.902 | 0.902 | 0.0 | False |
| B | EAS | f2_protein | 85 | 0.988 | 0.0 | 0.012 | False |
| B | EAS | f4_noncoding | 212 | 1.91 | 0.91 | 0.0 | False |
| B | EAS | known | 1002 | 1.893 | 0.894 | 0.001 | False |
| B | EUR | f2_protein | 60 | 1.0 | 0.0 | 0.0 | False |
| B | EUR | f4_noncoding | 185 | 1.93 | 0.941 | 0.011 | False |
| B | EUR | known | 1122 | 1.902 | 0.903 | 0.001 | False |
| B | MID | f2_protein | <20 | 0.917 | 0.0 | 0.083 | True |
| B | MID | f4_noncoding | 55 | 1.891 | 0.891 | 0.0 | False |
| B | MID | known | 240 | 1.946 | 0.946 | 0.0 | False |
| B | SAS | f2_protein | 53 | 0.981 | 0.0 | 0.019 | False |
| B | SAS | f3_synonymous | <20 | 2.0 | 1.0 | 0.0 | True |
| B | SAS | f4_noncoding | 209 | 1.914 | 0.919 | 0.005 | False |
| B | SAS | known | 899 | 1.907 | 0.907 | 0.0 | False |
| C | AFR | f2_protein | 135 | 0.911 | 0.0 | 0.089 | False |
| C | AFR | f4_noncoding | 371 | 1.787 | 0.792 | 0.005 | False |
| C | AFR | known | 1881 | 1.779 | 0.79 | 0.011 | False |
| C | AMR | f2_protein | 93 | 1.0 | 0.0 | 0.0 | False |
| C | AMR | f3_synonymous | <20 | 2.0 | 1.0 | 0.0 | True |
| C | AMR | f4_noncoding | 178 | 1.826 | 0.831 | 0.006 | False |
| C | AMR | known | 1388 | 1.885 | 0.887 | 0.002 | False |
| C | EAS | f2_protein | 84 | 0.964 | 0.0 | 0.036 | False |
| C | EAS | f4_noncoding | 207 | 1.85 | 0.85 | 0.0 | False |
| C | EAS | known | 1003 | 1.846 | 0.853 | 0.007 | False |
| C | EUR | f2_protein | 65 | 0.985 | 0.0 | 0.015 | False |
| C | EUR | f3_synonymous | <20 | 2.0 | 1.0 | 0.0 | True |
| C | EUR | f4_noncoding | 127 | 1.882 | 0.882 | 0.0 | False |
| C | EUR | known | 1164 | 1.882 | 0.882 | 0.0 | False |
| C | MID | f2_protein | <20 | 1.0 | 0.0 | 0.0 | True |
| C | MID | f4_noncoding | 36 | 1.833 | 0.833 | 0.0 | True |
| C | MID | known | 259 | 1.803 | 0.807 | 0.004 | False |
| C | SAS | f2_protein | 73 | 0.918 | 0.0 | 0.082 | False |
| C | SAS | f3_synonymous | <20 | 2.0 | 1.0 | 0.0 | True |
| C | SAS | f4_noncoding | 167 | 1.772 | 0.79 | 0.018 | False |
| C | SAS | known | 924 | 1.751 | 0.762 | 0.011 | False |
| DPA1 | AFR | f2_protein | 87 | 0.943 | 0.0 | 0.057 | False |
| DPA1 | AFR | f3_synonymous | <20 | 1.667 | 0.667 | 0.0 | True |
| DPA1 | AFR | f4_noncoding | 515 | 1.835 | 0.841 | 0.006 | False |
| DPA1 | AFR | known | 1786 | 1.849 | 0.852 | 0.003 | False |
| DPA1 | AMR | f2_protein | 41 | 1.0 | 0.0 | 0.0 | True |
| DPA1 | AMR | f4_noncoding | 351 | 1.937 | 0.937 | 0.0 | False |
| DPA1 | AMR | known | 1277 | 1.914 | 0.914 | 0.0 | False |
| DPA1 | EAS | f2_protein | 32 | 1.0 | 0.0 | 0.0 | True |
| DPA1 | EAS | f4_noncoding | 321 | 1.953 | 0.953 | 0.0 | False |
| DPA1 | EAS | known | 950 | 1.925 | 0.925 | 0.0 | False |
| DPA1 | EUR | f2_protein | 49 | 1.0 | 0.0 | 0.0 | True |
| DPA1 | EUR | f3_synonymous | <20 | 2.0 | 1.0 | 0.0 | True |
| DPA1 | EUR | f4_noncoding | 182 | 1.934 | 0.934 | 0.0 | False |
| DPA1 | EUR | known | 1141 | 1.909 | 0.911 | 0.002 | False |
| DPA1 | MID | f2_protein | <20 | 1.0 | 0.0 | 0.0 | True |
| DPA1 | MID | f4_noncoding | 62 | 1.952 | 0.968 | 0.016 | False |
| DPA1 | MID | known | 238 | 1.941 | 0.941 | 0.0 | False |
| DPA1 | SAS | f2_protein | 33 | 1.0 | 0.0 | 0.0 | True |
| DPA1 | SAS | f3_synonymous | <20 | 2.0 | 1.0 | 0.0 | True |
| DPA1 | SAS | f4_noncoding | 274 | 1.942 | 0.942 | 0.0 | False |
| DPA1 | SAS | known | 857 | 1.931 | 0.931 | 0.0 | False |
| DPB1 | AFR | f2_protein | 119 | 0.756 | 0.0 | 0.244 | False |
| DPB1 | AFR | f3_synonymous | <20 | 1.0 | 0.0 | 0.0 | True |
| DPB1 | AFR | f4_noncoding | 985 | 1.475 | 0.539 | 0.064 | False |
| DPB1 | AFR | known | 1278 | 1.663 | 0.696 | 0.033 | False |
| DPB1 | AMR | f2_protein | 66 | 0.955 | 0.0 | 0.045 | False |
| DPB1 | AMR | f4_noncoding | 406 | 1.7 | 0.719 | 0.02 | False |
| DPB1 | AMR | known | 1199 | 1.806 | 0.813 | 0.008 | False |
| DPB1 | EAS | f2_protein | 52 | 0.962 | 0.0 | 0.038 | False |
| DPB1 | EAS | f3_synonymous | <20 | 1.0 | 0.0 | 0.0 | True |
| DPB1 | EAS | f4_noncoding | 558 | 1.729 | 0.747 | 0.018 | False |
| DPB1 | EAS | known | 686 | 1.815 | 0.822 | 0.007 | False |
| DPB1 | EUR | f2_protein | 63 | 0.952 | 0.0 | 0.048 | False |
| DPB1 | EUR | f4_noncoding | 263 | 1.856 | 0.859 | 0.004 | False |
| DPB1 | EUR | known | 1040 | 1.859 | 0.862 | 0.003 | False |
| DPB1 | MID | f2_protein | <20 | 0.909 | 0.0 | 0.091 | True |
| DPB1 | MID | f3_synonymous | <20 | 2.0 | 1.0 | 0.0 | True |
| DPB1 | MID | f4_noncoding | 94 | 1.777 | 0.787 | 0.011 | False |
| DPB1 | MID | known | 202 | 1.767 | 0.777 | 0.01 | False |
| DPB1 | SAS | f2_protein | 40 | 0.9 | 0.0 | 0.1 | True |
| DPB1 | SAS | f3_synonymous | <20 | 2.0 | 1.0 | 0.0 | True |
| DPB1 | SAS | f4_noncoding | 414 | 1.87 | 0.872 | 0.002 | False |
| DPB1 | SAS | known | 710 | 1.832 | 0.839 | 0.007 | False |
| DQA1 | AFR | f2_protein | 92 | 0.707 | 0.0 | 0.293 | False |
| DQA1 | AFR | f4_noncoding | 661 | 1.327 | 0.431 | 0.104 | False |
| DQA1 | AFR | known | 1635 | 1.416 | 0.496 | 0.08 | False |
| DQA1 | AMR | f2_protein | 34 | 0.735 | 0.0 | 0.265 | True |
| DQA1 | AMR | f4_noncoding | 287 | 1.129 | 0.296 | 0.167 | False |
| DQA1 | AMR | known | 1351 | 1.375 | 0.486 | 0.112 | False |
| DQA1 | EAS | f2_protein | 36 | 0.611 | 0.0 | 0.389 | True |
| DQA1 | EAS | f4_noncoding | 368 | 1.283 | 0.44 | 0.158 | False |
| DQA1 | EAS | known | 896 | 1.198 | 0.348 | 0.151 | False |
| DQA1 | EUR | f2_protein | 36 | 0.75 | 0.0 | 0.25 | True |
| DQA1 | EUR | f4_noncoding | 288 | 1.365 | 0.462 | 0.097 | False |
| DQA1 | EUR | known | 1046 | 1.446 | 0.514 | 0.068 | False |
| DQA1 | MID | f2_protein | <20 | 0.625 | 0.0 | 0.375 | True |
| DQA1 | MID | f3_synonymous | <20 | 1.0 | 0.0 | 0.0 | True |
| DQA1 | MID | f4_noncoding | 71 | 0.972 | 0.239 | 0.268 | False |
| DQA1 | MID | known | 229 | 1.245 | 0.428 | 0.183 | False |
| DQA1 | SAS | f2_protein | 28 | 0.714 | 0.0 | 0.286 | True |
| DQA1 | SAS | f3_synonymous | <20 | 2.0 | 1.0 | 0.0 | True |
| DQA1 | SAS | f4_noncoding | 395 | 1.152 | 0.299 | 0.147 | False |
| DQA1 | SAS | known | 740 | 1.588 | 0.632 | 0.045 | False |
| DQB1 | AFR | f2_protein | 117 | 0.744 | 0.0 | 0.256 | False |
| DQB1 | AFR | f3_synonymous | <20 | 0.5 | 0.0 | 0.5 | True |
| DQB1 | AFR | f4_noncoding | 461 | 1.286 | 0.377 | 0.091 | False |
| DQB1 | AFR | known | 1807 | 1.453 | 0.529 | 0.076 | False |
| DQB1 | AMR | f2_protein | 57 | 0.912 | 0.0 | 0.088 | False |
| DQB1 | AMR | f3_synonymous | <20 | 1.5 | 0.5 | 0.0 | True |
| DQB1 | AMR | f4_noncoding | 244 | 1.467 | 0.537 | 0.07 | False |
| DQB1 | AMR | known | 1364 | 1.735 | 0.753 | 0.018 | False |
| DQB1 | EAS | f2_protein | 38 | 0.974 | 0.0 | 0.026 | True |
| DQB1 | EAS | f3_synonymous | <20 | 2.0 | 1.0 | 0.0 | True |
| DQB1 | EAS | f4_noncoding | 229 | 1.812 | 0.817 | 0.004 | False |
| DQB1 | EAS | known | 1028 | 1.874 | 0.879 | 0.006 | False |
| DQB1 | EUR | f2_protein | 39 | 0.897 | 0.0 | 0.103 | True |
| DQB1 | EUR | f4_noncoding | 179 | 1.615 | 0.637 | 0.022 | False |
| DQB1 | EUR | known | 1151 | 1.753 | 0.767 | 0.014 | False |
| DQB1 | MID | f2_protein | <20 | 0.923 | 0.0 | 0.077 | True |
| DQB1 | MID | f4_noncoding | 47 | 1.596 | 0.596 | 0.0 | True |
| DQB1 | MID | known | 249 | 1.715 | 0.735 | 0.02 | False |
| DQB1 | SAS | f2_protein | 36 | 0.889 | 0.0 | 0.111 | True |
| DQB1 | SAS | f3_synonymous | <20 | 2.0 | 1.0 | 0.0 | True |
| DQB1 | SAS | f4_noncoding | 184 | 1.582 | 0.592 | 0.011 | False |
| DQB1 | SAS | known | 943 | 1.749 | 0.764 | 0.015 | False |
| DRB1 | AFR | f2_protein | 89 | 0.978 | 0.0 | 0.022 | False |
| DRB1 | AFR | f4_noncoding | 1974 | 1.891 | 0.895 | 0.004 | False |
| DRB1 | AFR | known | 313 | 1.716 | 0.716 | 0.0 | False |
| DRB1 | AMR | f2_protein | 47 | 1.0 | 0.0 | 0.0 | True |
| DRB1 | AMR | f4_noncoding | 1320 | 1.87 | 0.876 | 0.005 | False |
| DRB1 | AMR | known | 293 | 1.765 | 0.765 | 0.0 | False |
| DRB1 | EAS | f2_protein | 49 | 0.98 | 0.0 | 0.02 | True |
| DRB1 | EAS | f4_noncoding | 1114 | 1.864 | 0.866 | 0.003 | False |
| DRB1 | EAS | known | 129 | 1.729 | 0.729 | 0.0 | False |
| DRB1 | EUR | f2_protein | 58 | 0.983 | 0.0 | 0.017 | False |
| DRB1 | EUR | f4_noncoding | 884 | 1.862 | 0.869 | 0.007 | False |
| DRB1 | EUR | known | 422 | 1.834 | 0.834 | 0.0 | False |
| DRB1 | MID | f2_protein | <20 | 1.0 | 0.0 | 0.0 | True |
| DRB1 | MID | f4_noncoding | 253 | 1.929 | 0.929 | 0.0 | False |
| DRB1 | MID | known | 44 | 1.727 | 0.727 | 0.0 | True |
| DRB1 | SAS | f2_protein | 33 | 0.97 | 0.0 | 0.03 | True |
| DRB1 | SAS | f3_synonymous | <20 | 2.0 | 1.0 | 0.0 | True |
| DRB1 | SAS | f4_noncoding | 988 | 1.889 | 0.892 | 0.003 | False |
| DRB1 | SAS | known | 142 | 1.789 | 0.789 | 0.0 | False |


Person-genes excluded because the long-read call's very first field was unresolved (`novelty_depth` 1 / undetermined, where a two-field comparison is not defined): **252**.


## What must still be checked before quoting any of this

- Kinship alone cannot tell a parent-child pair from a full sibling pair (the All of Us relatedness table has no IBD0 column - see script 11). Without `--yob-tsv` the parent-child-specific 'must share a copy at every gene' logic cannot be applied to any pair.

- `dropout_candidate` is a *shape*, not a proof: a genuinely homozygous person who differs from their relative lands in it too. Its value is as an upper bound, read together with the homozygosity obs/exp ratio in step 5.

- `fragmentation_candidate` uses the contig labels in Table 1; a gene alone on a contig is suspicious but not wrong by itself.

- The synthetic-switch power estimate measures sensitivity to ONE clean hap-label swap inside one contig. It does not estimate sensitivity to short switch segments, to errors that also corrupt the sequence, or to switches that coincide with a contig break (which are undetectable by construction).

- Different-length sequences fall back to difflib alignment (inherited from script 17), which can misalign in repetitive stretches; equal-length comparisons are position-anchored and exact.

- copy_index > 1 rows (segmental duplication / DRB copy number) are excluded everywhere here.


## Disclosure discipline

No person identifiers appear in this directory. Any person or pair count between 1 and 19 is written as `<20`. Person-level rows, if written at all, went only to the separate local directory.
