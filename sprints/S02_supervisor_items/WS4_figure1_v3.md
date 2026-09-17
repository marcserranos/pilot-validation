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
| A | 3.5 | 2.3 | 2.9 | 3.0 | ~0* | 2.5 |
| B | 12.1 | 10.1 | 9.5 | 7.6 | 10.0 | 10.3 |
| C | 8.6 | 6.1 | 8.6 | 5.4 | 6.2 | 8.1 |
| DPA1 | 12.4 | 12.0 | 13.9 | 7.4 | 10.6 | 13.3 |
| DPB1 | 25.5 | 14.1 | 27.1 | 11.8 | 19.3 | 21.6 |
| DQA1 | 16.3 | 9.7 | 17.0 | 12.0 | 14.0 | 19.3 |
| DQB1 | 11.0 | 8.2 | 9.7 | 7.2 | 8.7 | 8.5 |
| DRB1 | **65.7** | 60.0 | **71.2** | **45.7** | 64.7 | 67.2 |

\* MID is the smallest strict-ancestry group and most of its cells are suppressed; read its
interval, not its lower bound.

The artifact control sits at **1.1–3.3%** (lower bound) in every gene and every ancestry except
MID, where suppression makes the interval wide. So the EUR-vs-everyone-else gradient — clearest at
DRB1 (45.7% EUR vs 71.2% EAS) and DPB1 (11.8% vs 27.1%) — is not explained by uneven assembly
quality.

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
