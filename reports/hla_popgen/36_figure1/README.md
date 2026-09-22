# 36 — Figure 1 (current version)

Drawn natively from data as a single figure on the VM. **Supersedes `35_figure1/`**, which pasted
four rendered PNGs into one sheet and had four different type scales, three internal working
titles, and an aspect ratio no journal would take.

## Draft caption

**Figure 1 | Long-read HLA typing across the All of Us cohort shows where the reference catalogue
is incomplete.**
**(a)** Genetic-ancestry composition of 11,833 unrelated long-read participants, ordered within
each predicted-ancestry block from least to most admixed. The figure beneath each block is the
share of that group whose dominant ancestry probability reaches 0.98, the strict threshold used in
(b) and (c).
**(b)** HLA-B alleles placed on the AFR/EUR/AMR simplex at the mean renormalised ancestry
composition of their carriers; point area and opacity scale with carrier count. Alleles sit close
to the vertices, i.e. HLA-B allele frequencies are strongly structured by ancestry.
**(c)** Percentage of called haplotypes carrying sequence absent from IPD-IMGT/HLA, per classical
gene and ancestry, split by the nomenclature field at which the novelty appears. Black bars are
the assembly/annotation artifact rate, which is flat across ancestries and therefore acts as a
negative control: the ancestry gradient in the coloured stack is reference incompleteness, not
uneven assembly quality.
**(d)** Distinct novel protein alleles per gene, split by how many unrelated people carry them.
The discovery is concentrated in TAP1/TAP2, MIC and the class II accessory and paralogous genes,
not in the classical HLA genes.

## What changed from script 35, and why

| Problem in 35 | Fix here |
|---|---|
| Collage of PNGs, four type scales | One figure, one type scale, one palette |
| Panels (a)/(b) at the original admixture threshold | `--strict-threshold`, default **0.98** (ask A9) |
| Ternary plotted alleles with **no minimum carrier count** — a one-carrier point is that person's own admixture | `--min-carriers`, default 20 |
| Spectrum panel silently drew an all-zero novel series | Replaced; see below |
| No underlying tables | Panels b/c/d export aggregate tables, so they can be restyled offline |

## The suppressed-count problem behind panel (d)

Panel (d) was originally an allele-frequency spectrum at the corrected novelty definition. It
could not be built. Script 24 applies the All of Us small-cell rule **when it writes its output**,
so every carrier count from 1 to 19 is stored as `<20` even in the copy on the VM — the uncensored
counts are not recoverable from that file. The first run of this script therefore produced a novel
series of exactly zero in every bin, which looked like a plotting bug but was the data.

Panel (d) now shows the three-class split the data actually supports (seen once / 2–19 / ≥20) per
gene. A true per-count histogram, or a spectrum at the corrected definition, needs script 24 rerun
with an uncensored internal dump — a change to script 24, not just a rerun.

## The cost of the strict threshold

Raising the ancestry-probability floor from 0.90 to 0.98 takes the analysable set from **11,856 to
6,136 people**. It costs most where admixture is highest: only **32%** of European-ancestry and
**46%** of American-ancestry participants clear it, against **93%** of South Asian and **84%** of
East Asian. That trade-off is a result in itself and is worth showing the supervisors, since they
asked for the stricter filter.

## Files

- `figure1.png` (300 dpi) / `figure1.pdf`
- `panel_b_ternary_alleles.tsv` — per-allele simplex coordinates and suppressed carrier counts
- `panel_c_novelty_rates.tsv` — exact rates per gene × ancestry × field, plus the artifact rate
- `panel_d_recurrence_by_gene.tsv` — novel proteins per gene by recurrence class
- `panel_a_group_sizes.tsv` — group sizes and the share clearing the strict threshold

## Still open

- Panels are per-panel legible but the figure has not been through a journal-width check (180 mm).
- `--gene` defaults to HLA-B for the ternary; HLA-A is the alternative Cole mentioned.
