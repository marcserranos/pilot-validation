# Raw-sequence divergence -- nomenclature-bypassed phasing/typing validation

Total gene-comparisons (545 high-sharing pairs x genes both sides have a sequence for): **19700**

Exact sequence match (raw_shared): **18680** (94.822%)

No exact match, closest pairing still differs: **1020** (5.178%, 95% CI 4.877-5.496%)


## How many bases actually differ (among the non-matching comparisons)

| stat | value |
|---|---|
| count | 1020.0 |
| mean | 14.3 |
| std | 43.0 |
| min | 1.0 |
| 25% | 1.0 |
| 50% | 3.0 |
| 75% | 12.0 |
| 90% | 29.0 |
| max | 458.0 |


- Same length (pure substitution-style, no indel): **502/1020** (49.2%)

- 1-3 differing bases specifically (near-miss / point-level): **578/1020** (56.7%)


## Shape: one clustered block vs scattered across the sequence

| n_diff_blocks | n_gene_comparisons | meaning |
|---|---|---|
| 1 | 418 | one localized change (point substitution or small indel) |
| 2 | 121 | two separate differing regions |
| 3 | 73 | scattered across >=3 separate regions |
| 4 | 42 | scattered across >=3 separate regions |
| 5 | 31 | scattered across >=3 separate regions |
| 6 | 33 | scattered across >=3 separate regions |
| 7 | 26 | scattered across >=3 separate regions |
| 8 | 24 | scattered across >=3 separate regions |
| 9 | 21 | scattered across >=3 separate regions |
| 10 | 17 | scattered across >=3 separate regions |
| 11 | 8 | scattered across >=3 separate regions |
| 12 | 10 | scattered across >=3 separate regions |
| 13 | 27 | scattered across >=3 separate regions |
| 14 | 10 | scattered across >=3 separate regions |
| 15 | 14 | scattered across >=3 separate regions |
| 16 | 15 | scattered across >=3 separate regions |
| 17 | 9 | scattered across >=3 separate regions |
| 18 | 17 | scattered across >=3 separate regions |
| 19 | 5 | scattered across >=3 separate regions |
| 20 | 6 | scattered across >=3 separate regions |
| 21 | 7 | scattered across >=3 separate regions |
| 22 | 3 | scattered across >=3 separate regions |
| 23 | 7 | scattered across >=3 separate regions |
| 24 | 8 | scattered across >=3 separate regions |
| 25 | 6 | scattered across >=3 separate regions |
| 26 | 8 | scattered across >=3 separate regions |
| 27 | 3 | scattered across >=3 separate regions |
| 28 | 4 | scattered across >=3 separate regions |
| 29 | 2 | scattered across >=3 separate regions |
| 30 | 3 | scattered across >=3 separate regions |
| 31 | 5 | scattered across >=3 separate regions |
| 32 | 3 | scattered across >=3 separate regions |
| 33 | 2 | scattered across >=3 separate regions |
| 34 | 2 | scattered across >=3 separate regions |
| 35 | 1 | scattered across >=3 separate regions |
| 36 | 4 | scattered across >=3 separate regions |
| 37 | 4 | scattered across >=3 separate regions |
| 38 | 4 | scattered across >=3 separate regions |
| 39 | 1 | scattered across >=3 separate regions |
| 40 | 1 | scattered across >=3 separate regions |
| 41 | 2 | scattered across >=3 separate regions |
| 43 | 1 | scattered across >=3 separate regions |
| 44 | 1 | scattered across >=3 separate regions |
| 51 | 1 | scattered across >=3 separate regions |
| 52 | 1 | scattered across >=3 separate regions |
| 53 | 1 | scattered across >=3 separate regions |
| 55 | 3 | scattered across >=3 separate regions |
| 56 | 1 | scattered across >=3 separate regions |
| 80 | 1 | scattered across >=3 separate regions |
| 85 | 1 | scattered across >=3 separate regions |
| 95 | 1 | scattered across >=3 separate regions |
| 107 | 1 | scattered across >=3 separate regions |

## By gene class

| gene_class | n_compared | n_diverged | divergence_rate_% |
|---|---|---|---|
| other | 51 | 9 | 17.647 |
| mic_tap | 2170 | 184 | 8.479 |
| classical_I | 1631 | 98 | 6.009 |
| pseudogene_I | 6049 | 354 | 5.852 |
| classical_II | 2723 | 130 | 4.774 |
| class_II_paralog | 2720 | 127 | 4.669 |
| nonclassical_I | 1632 | 71 | 4.350 |
| class_II_accessory | 2724 | 47 | 1.725 |