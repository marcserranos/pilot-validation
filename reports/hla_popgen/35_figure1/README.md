# 35 — Figure 1, assembled

The four panels Cole specified on 2026-09-17, laid out as one sheet.

| Panel | Content | Source |
|---|---|---|
| a | Genetic ancestry of the long-read cohort (admixture barcode) | script 06 |
| b | HLA-B allele ancestry simplex — the ternary Cole said he liked "as it is" | script 10 |
| c | Novelty by gene, ancestry and nomenclature field, with the artifact rate as a flat negative control | script 33 |
| d | Allele-frequency spectrum, novel vs known, density-normalised | script 33 |

*(c and d arrive as one wide sub-sheet from script 33 and occupy the bottom row together.)*

## Draft caption

**Figure 1 | Long-read HLA typing across 11,856 unrelated All of Us participants reveals where the
reference catalogue is incomplete.**
**(a)** Genetic-ancestry composition of the long-read cohort.
**(b)** Per-allele ancestry centroids for HLA-B on the AFR/EUR/AMR simplex; circles are alleles
catalogued in IPD-IMGT/HLA, triangles are alleles first observed here. Alleles sit close to the
vertices, i.e. HLA-B allele frequencies are strongly structured by ancestry.
**(c)** Percentage of called haplotypes carrying sequence absent from IPD-IMGT/HLA, per classical
gene and per ancestry, split by the nomenclature field at which the novelty appears (new protein,
new synonymous CDS, new non-coding only). Black bars show the assembly/annotation artifact rate,
which is flat across ancestries (0.90–3.22% in the five well-powered groups) and therefore acts as
a negative control: the ancestry gradient in the coloured stack is reference incompleteness, not
uneven assembly quality. Grey whiskers show how high the stack could reach if every
disclosure-suppressed count sat at its ceiling.
**(d)** Allele-frequency spectrum over all classical-gene alleles, novel versus already catalogued,
with counts divided by bin width. Novel alleles collapse after a few observations; catalogued
alleles reach hundreds of haplotypes.

## Honest limitations of this version

1. **This is a layout, not a re-analysis.** Panels are composited from rendered PNGs. Scripts 06
   and 10 committed figures but not their underlying per-allele frequency tables, so those two
   panels cannot be re-plotted offline. For a manuscript figure both must be rerun to emit vector
   output plus their tables.
2. **Panels (a) and (b) are at the original admixture threshold**, not the stricter ~98% Cole and
   David asked for (ask A9). Scripts 06 and 10 have no such parameter today — that is a code
   change, not just a rerun.
3. **Panel (d) predates S01's corrected novelty definition.** Its shape is right; the absolute
   novel-allele counts are over-estimates.
