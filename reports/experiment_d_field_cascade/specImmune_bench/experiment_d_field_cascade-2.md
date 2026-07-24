# Experiment D -- field-level cascade (AoU-native / SpecHLA vs SpecImmune-LR truth)

Each field level is cumulative (matching through Field N implies Field 1..N-1 also matched). N/A means the comparison isn't assessable at that depth (one caller didn't report that many fields) -- not counted as a disagreement.

## Pooled cascade -- match rate through each field, all ancestries combined

Match rate = n_true / (n_true + n_false) among *assessable* comparisons; N/A = comparisons where one caller simply doesn't report that many fields (e.g. AoU-native never reports Field 4).

**"structural cap"** means the method never independently reaches that field at all, in anyone, anywhere in this cohort (e.g. AoU-native: max observed depth is 3 fields) -- distinct from an ordinary per-call N/A, which just means this particular call/truth pair happened to be shallow. Cells at or below a method's own max depth are real, measured comparisons even when they show 0% (that 0% is inherited from an earlier field that already disagreed, per the cascade definition above -- not fresh evidence at this depth).

### AoU-native vs SpecImmune-LR truth (max observed depth: 3 fields)

| Gene | Field 1 (allele group) | Field 2 (protein -- non-synonymous) | Field 3 (synonymous, coding) | Field 4 (non-coding) |
|---|---|---|---|---|
| A | 92% (N/A 0%) | 87% (N/A 0%) | 85% (N/A 1%) | N/A -- structural cap |
| B | 93% (N/A 0%) | 88% (N/A 0%) | 82% (N/A 0%) | N/A -- structural cap |
| C | 97% (N/A 0%) | 92% (N/A 0%) | 82% (N/A 0%) | N/A -- structural cap |
| DRB1 | 94% (N/A 0%) | 77% (N/A 0%) | 75% (N/A 0%) | N/A -- structural cap |
| DQA1 | 100% (N/A 0%) | 64% (N/A 0%) | 59% (N/A 0%) | N/A -- structural cap |
| DQB1 | 97% (N/A 0%) | 83% (N/A 0%) | 83% (N/A 0%) | N/A -- structural cap |
| DPA1 | 98% (N/A 0%) | 93% (N/A 0%) | 88% (N/A 0%) | N/A -- structural cap |
| DPB1 | 83% (N/A 0%) | 83% (N/A 0%) | 83% (N/A 0%) | N/A -- structural cap |

### SpecHLA vs SpecImmune-LR truth (max observed depth: 4 fields)

| Gene | Field 1 (allele group) | Field 2 (protein -- non-synonymous) | Field 3 (synonymous, coding) | Field 4 (non-coding) |
|---|---|---|---|---|
| A | 92% (N/A 0%) | 83% (N/A 0%) | 80% (N/A 2%) | 72% (N/A 3%) |
| B | 93% (N/A 0%) | 85% (N/A 0%) | 74% (N/A 0%) | 44% (N/A 1%) |
| C | 97% (N/A 0%) | 92% (N/A 0%) | 92% (N/A 0%) | 58% (N/A 0%) |
| DRB1 | 69% (N/A 0%) | 52% (N/A 0%) | 48% (N/A 0%) | 34% (N/A 2%) |
| DQA1 | 100% (N/A 0%) | 94% (N/A 0%) | 92% (N/A 0%) | 45% (N/A 11%) |
| DQB1 | 97% (N/A 0%) | 95% (N/A 0%) | 95% (N/A 0%) | 75% (N/A 3%) |
| DPA1 | 100% (N/A 0%) | 93% (N/A 0%) | 88% (N/A 0%) | 37% (N/A 0%) |
| DPB1 | 85% (N/A 0%) | 85% (N/A 0%) | 83% (N/A 0%) | 47% (N/A 3%) |

## Field 2 (protein-level) match rate by gene x ancestry

The headline correctness metric -- does the caller get the actual HLA protein right, not just the allele family.

### AoU-native

| Gene | AFR | AMR | EAS | EUR | MID | SAS |
|---|---|---|---|---|---|---|
| A | 95% (n=20) | 90% (n=20) | 75% (n=20) | 80% (n=20) | 90% (n=20) | 90% (n=20) |
| B | 80% (n=20) | 85% (n=20) | 80% (n=20) | 100% (n=20) | 85% (n=20) | 95% (n=20) |
| C | 95% (n=20) | 90% (n=20) | 90% (n=20) | 95% (n=20) | 90% (n=20) | 90% (n=20) |
| DRB1 | 60% (n=20) | 75% (n=20) | 80% (n=20) | 85% (n=20) | 80% (n=20) | 80% (n=20) |
| DQA1 | 75% (n=20) | 75% (n=20) | 45% (n=20) | 55% (n=20) | 55% (n=20) | 80% (n=20) |
| DQB1 | 85% (n=20) | 75% (n=20) | 95% (n=20) | 80% (n=20) | 75% (n=20) | 90% (n=20) |
| DPA1 | 85% (n=20) | 95% (n=20) | 90% (n=20) | 100% (n=20) | 95% (n=20) | 95% (n=20) |
| DPB1 | 85% (n=20) | 85% (n=20) | 90% (n=20) | 90% (n=20) | 65% (n=20) | 85% (n=20) |

### SpecHLA

| Gene | AFR | AMR | EAS | EUR | MID | SAS |
|---|---|---|---|---|---|---|
| A | 85% (n=20) | 85% (n=20) | 75% (n=20) | 80% (n=20) | 90% (n=20) | 85% (n=20) |
| B | 80% (n=20) | 70% (n=20) | 85% (n=20) | 95% (n=20) | 85% (n=20) | 95% (n=20) |
| C | 100% (n=20) | 95% (n=20) | 90% (n=20) | 90% (n=20) | 95% (n=20) | 85% (n=20) |
| DRB1 | 33% (n=18) | 56% (n=16) | 58% (n=12) | 61% (n=18) | 67% (n=12) | 44% (n=18) |
| DQA1 | 90% (n=20) | 100% (n=20) | 100% (n=20) | 95% (n=20) | 90% (n=20) | 90% (n=20) |
| DQB1 | 95% (n=20) | 100% (n=20) | 95% (n=20) | 95% (n=20) | 95% (n=20) | 90% (n=20) |
| DPA1 | 85% (n=20) | 95% (n=20) | 85% (n=20) | 95% (n=20) | 100% (n=20) | 100% (n=20) |
| DPB1 | 95% (n=20) | 85% (n=20) | 80% (n=20) | 95% (n=20) | 60% (n=20) | 95% (n=20) |

