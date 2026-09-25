# 38 -- deletion supplement figure (HLA-A / MHC gene deletions, style pass)

Source: script `38b_deletion_supplement_fig.py`, re-rendered from the already-committed tables
in `reports/hla_popgen/30_hla_sv/` (`deletion_rates_by_gene.tsv`,
`deletion_rates_by_gene_ancestry.tsv`). No VM access, no re-computation from raw calls -- every
number here was already public.

Cole (call #8): "simplify the colour scheme -- drop the orange/green distinction, keep
everything else." Script 30's original `fig_deletion_rates()` colour-coded genes green
(positive control) / grey (negative control) / orange (everything else). This version uses one
neutral colour for every gene; the two control groups are called out with a bracket + label
above their bars instead.

## Panel (a) -- per-gene bridged-deletion rate

Every gene in the committed table, bridged-deletion % (Wilson 95% CI), grouped by SCHEMA.md's
`gene_class` (used here as a documented stand-in for genomic position, since the committed table
carries gene class but not gene_start/gene_end coordinates). Thin grey vertical lines mark class
boundaries. The two control groups from script 30 -- `DRB3, DRB4, DRB5` (expected copy-number genes) and
`A, B, C, DRA, DQA1, DQB1, DPA1, DPB1, DRB1` (negative controls that should essentially never be deleted) -- are marked with a bracket
and label above their bars, not by colour.

**How to read it:** a real, high bridged-deletion rate at DRB3/4/5 is the expected biology (DR51/
52/53 haplotype groups routinely lack the paralogous DRB locus). A near-zero rate at the negative
controls is the method's false-positive floor; anything meaningfully above zero there is measuring
assembly/detection error, not biology.

**Caveat:** 0 gene(s) were omitted from this panel because their bridged-haplotype
denominator itself is small-cell suppressed (`<20`) in the committed table -- see
`deletion_rates_by_gene.tsv`.

## Panel (b) -- CNV-gene deletion frequency by ancestry

The five expected-CNV genes (DRB3, DRB4, DRB5, C4A, C4B) as a dot-plot: bridged-deletion % with Wilson 95% CI, one dot
per (gene, ancestry), coloured by the project's standard ancestry palette
(`ANCESTRY_COLORS` in `_viz_common.py`).

**How to read it:** DRB3/4/5 frequencies track DRB1 haplotype-group frequencies, which differ by
ancestry for well-established population-genetic reasons (not an artifact of this pipeline); C4A/
C4B copy number is independently variable. A cell drawn as "suppressed (n<20)" text instead of a
dot means the (gene, ancestry) bridged-haplotype count itself was `<20` in the committed table --
0 such cell(s) here -- and no rate is shown for it, per the AoU small-cell rule.
`<20` was parsed as censored throughout (never as 0).

## Caveats

- This figure cannot detect anything script 30 did not already detect: it is a restyle of
  committed, already-suppressed tables, not a new analysis.
- `pct_deleted_bridged` is the EXACT percentage script 30 computed before applying small-cell
  suppression to the underlying counts, so it is trusted as-is; the Wilson CI's count term is
  reconstructed from that exact percentage and the (uncensored) bridged-haplotype denominator,
  which discloses nothing beyond what the committed TSV already implies algebraically.
- "Deleted" here means *bridged absence* (script 30's stricter definition: a single contig spans
  the gene's position and the gene is not there) -- not raw absence, which conflates real
  deletions with assembly fragmentation. See `30_hla_structural_variation.py`'s module docstring.
