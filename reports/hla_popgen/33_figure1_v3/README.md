# 33 — Figure 1 v3 (panels c and d)

*Supervisor ask A1 (Cole, 2026-09-17). Built from committed aggregates only.*


## What is here and what is not

| panel | Cole's ask | status |
|---|---|---|
| a | admixture bar plot | **not in this file** — needs a rerun of script 06 at the stricter admixture threshold (ask A9) |
| b | class I ternary, HLA-A or HLA-B | **not in this file** — script 10 committed figures but not the per-allele frequency table, so it cannot be recomposed offline |
| c | novel rate by ancestry, broken out by gene | **here**, and rebuilt at S01's corrected novelty definition |
| d | SFS over all alleles, novel vs known in two colours | **here** |


## Panel c — read this before quoting a number from it

The rate presented in the call was a single 'novel' rate per ancestry. S01 showed that number is dominated by non-coding novelty and by artifacts, so this panel splits it by nomenclature field — new protein, new synonymous CDS, new non-coding only — and draws the artifact rate as a flat black bar underneath.


**Suppressed counts.** 3,675 of the 8,784 cells in script 24's table are written `<20` under the AoU small-cell rule. Every rate here is therefore an interval: the tabulated number is the **lower bound** (suppressed cells counted as zero) and the grey whisker in the figure shows where the stack would reach if every suppressed cell sat at 19. This matters most for MID, the smallest strict-ancestry group, where a naive reading of the suppressed cells as zero produces an artifact rate of exactly 0.00% across every gene — which is what the first version of this script reported, and it looked entirely plausible.


The artifact bar is the control. Assembly and annotation artifacts are not ancestry-structured, so a gradient in the coloured stack against a flat black bar is reference incompleteness, not uneven assembly quality.


### Total clean novelty rate (%% of called haplotypes), all fields


| gene     |   AFR |   AMR |   EAS |   EUR |   MID |   SAS |
|:---------|------:|------:|------:|------:|------:|------:|
| HLA-A    |  3.48 |  2.25 |  2.83 |  2.9  |  0    |  2.41 |
| HLA-B    | 11.95 |  9.86 |  9.17 |  7.34 |  9.13 | 10.06 |
| HLA-C    |  8.46 |  5.95 |  8.42 |  5.28 |  5.45 |  7.78 |
| HLA-DPA1 | 12.2  | 11.83 | 13.66 |  7.23 |  9.41 | 12.96 |
| HLA-DPB1 | 25    | 13.86 | 26.26 | 11.42 | 17.59 | 20.91 |
| HLA-DQA1 | 16.07 |  9.48 | 16.63 | 11.76 | 12.37 | 18.66 |
| HLA-DQB1 | 10.72 |  7.87 |  9.36 |  7.01 |  7.5  |  8.2  |
| HLA-DRB1 | 64.87 | 58.88 | 69.58 | 44.68 | 58.78 | 64.96 |


### Artifact rate — the control, expected to be flat across ancestries


| gene     |   AFR |   AMR |   EAS |   EUR |   MID |   SAS |
|:---------|------:|------:|------:|------:|------:|------:|
| HLA-A    |  2.81 |  1.77 |  2    |  2.35 |     0 |  1.64 |
| HLA-B    |  2.83 |  2.1  |  2.95 |  2.23 |     0 |  2.32 |
| HLA-C    |  2.65 |  2.96 |  3.22 |  2.44 |     0 |  2.9  |
| HLA-DPA1 |  1.71 |  1.17 |  1.16 |  1.53 |     0 |  1.53 |
| HLA-DPB1 |  2.4  |  2.02 |  1.92 |  2.52 |     0 |  1.49 |
| HLA-DQA1 |  1.76 |  0.9  |  1.08 |  1.23 |     0 |  1.13 |
| HLA-DQB1 |  2.44 |  1.76 |  1.59 |  1.54 |     0 |  1.56 |
| HLA-DRB1 |  2.13 |  1.67 |  2.18 |  2.36 |     0 |  1.45 |


## Panel d

Spectrum for ALL_CLASSICAL. Counts are divided by bin width: the bins double in size, so raw counts per bin produce a sawtooth that looks like structure and is entirely binning. Flagged-artifact alleles are drawn separately and are *not* counted as novel.


## Caveat on panel d's novelty definition

Panel d reuses script 23's spectrum, which predates S01's corrected split. Its 'novel' series therefore still includes calls whose CDS turned out to be already catalogued. The shape of the spectrum — novel alleles collapsing after a couple of occurrences while known alleles reach hundreds — is unaffected, but the absolute allele counts in the novel series are an over-estimate. Regenerating it at the corrected definition needs the VM.
