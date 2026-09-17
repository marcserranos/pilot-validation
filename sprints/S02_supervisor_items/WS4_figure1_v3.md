# WS4 — Figure 1 v3 (ask A1), and what still needs the VM

Cole's four panels, and where each stands.

| panel | ask | status |
|---|---|---|
| a | admixture bar plot | **needs VM.** Script 06 committed figures but not the underlying per-person admixture table. Must be rerun anyway at the stricter threshold (ask A9). |
| b | class I ternary, HLA-A or HLA-B, "the best one you can make" | **needs VM.** Script 10 (`reports/hla_popgen/10_allele_ancestry_geometry/gene_B/`) has the existing figures; the per-allele per-ancestry frequency table was never committed, so it cannot be recomposed offline. |
| c | novel rate by ancestry × gene | **done**, `reports/hla_popgen/33_figure1_v3/`, rebuilt at S01's corrected novelty definition. |
| d | SFS, all alleles, novel vs known in two colours | **done**, same folder, with the caveat below. |

## What panel c now says

Split by nomenclature field rather than one lumped "novel" rate, with the artifact rate drawn
underneath as a flat negative control. Clean novelty rate (lower bound, % of called haplotypes,
strict ancestry, unrelated):

| gene | AFR | AMR | EAS | EUR | MID | SAS |
|---|---|---|---|---|---|---|
| A | 3.5 | 2.3 | 2.8 | 2.9 | 0.0 | 2.4 |
| B | 12.0 | 9.9 | 9.2 | 7.3 | 9.1 | 10.1 |
| C | 8.5 | 5.9 | 8.4 | 5.3 | 5.5 | 7.8 |
| DPA1 | 12.2 | 11.8 | 13.7 | 7.2 | 9.4 | 13.0 |
| DPB1 | 25.0 | 13.9 | 26.3 | 11.4 | 17.6 | 20.9 |
| DQA1 | 16.1 | 9.5 | 16.6 | 11.8 | 12.4 | 18.7 |
| DQB1 | 10.7 | 7.9 | 9.4 | 7.0 | 7.5 | 8.2 |
| DRB1 | 64.9 | 58.9 | 69.6 | 44.7 | 58.8 | 65.0 |

These are **lower bounds** (suppressed cells counted as zero). MID is the smallest
strict-ancestry group and most of its cells are suppressed, so its lower bound is the least
informative of the six — read its interval, not this column.

The artifact control sits between **0.90% and 3.22%** (lower bound) across the five well-powered
ancestries, with an upper bound never above 5.6%. For MID it is bounded only between 0% and
6.5–13.3%, so MID cannot carry the control argument. Against that flat control the
EUR-vs-everyone-else gradient — clearest at DRB1 (44.7% EUR vs 69.6% EAS) and DPB1 (11.4% vs
26.3%) — is not explained by uneven assembly quality.

## The suppression trap, recorded because it nearly shipped

The first version of script 33 read the source table with `pd.to_numeric(...).fillna(0)`. **3,675
of its 8,784 cells are written `<20`** under the AoU small-cell rule, so every one of them became
a zero. The panel then reported a 0.00% artifact rate for Middle Eastern ancestry across all eight
genes: a flat, clean, plausible-looking row that was entirely fabricated. It was caught only
because a whole column of exact zeros looked too tidy.

Every rate in panel c is now an interval: lower bound counts suppressed cells as 0, upper bound as
19, with the numerator and denominator bounds paired so the interval actually brackets the truth.
`tests/test_figure1_v3_compose.py` pins this.

## Panel d caveat

Panel d reuses script 23's spectrum, which predates S01's corrected novelty split, so its "novel"
series still contains calls whose CDS turned out to be already catalogued. The *shape* — novel
alleles collapsing after two or three occurrences while known alleles reach hundreds of
haplotypes — is unaffected; the absolute novel allele counts are an over-estimate. Regenerating at
the corrected definition needs the VM.

## VM rerun needed for a, b and d

1. Script 06 with a `--strict-threshold 0.98` sweep, committing the admixture table, not just the
   figure.
2. Script 10 for HLA-B (and HLA-A as the alternate), committing the per-allele per-ancestry
   frequency table so the ternary can be recomposed and restyled offline in future.
3. Script 23's SFS regenerated against script 24's corrected cluster definitions.
