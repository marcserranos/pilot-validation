# 55 — Does blood RNA expression confirm assembly-based presence/absence calls?

Script: `scripts/hla_popgen/55_structure_expression_validation.py` (VM 2026-10-06). Expression: AoU v9 RNA-SeQC gene TPM
(GRCh38 primary annotation, 8,980 samples). Copy number = number of haplotypes on which the gene is called.
KIR "strict" = both haplotypes complete (51); DRB5 "strict" = DRB1 and DQA1 on one contig on both haplotypes.
Table: `structure_expression.tsv`; figure: `fig_structure_expression.png`.

| gene (set) | median TPM at 0 / 1 / 2 copies | Spearman ρ | n |
|---|---|---|---|
| HLA-DRB5 (strict) | 4.6 / 51.2 / 110.0 | 0.72 | 6,750 |
| KIR2DL3 (strict) | 0.50 / 0.90 / 1.36 | 0.30 | 4,994 |
| KIR2DL1 (strict) | 0.59 / 0.94 / 1.10 | 0.12 | 4,994 |
| KIR2DL4 (strict) | <20 people at 0 / 0.18 / 0.32 | 0.16 | 4,994 |
| KIR3DL1 (strict) | 1.90 / 1.53 / 1.20 | −0.14 | 4,994 |
| KIR2DS4 (strict) | 2.88 / 1.88 / 1.25 | −0.25 | 4,994 |

- DRB5: clean dosage (about 2× per copy), and no carrier below 1 TPM. This is orthogonal confirmation of the DR51
  presence/absence calls, and through them of the DRB1 → paralog rules recovered in 52a. The 0-copy floor (median 4.6 TPM)
  is reads from other DRB genes mapping to the DRB5 reference locus.
- KIR: TPMs are low in bulk blood (NK/T-cell genes). KIR2DL3, 2DL1 and 2DL4 rise with copy number. KIR3DL1 and KIR2DS4 fall,
  because people lacking them carry paralogs that are missing from the primary reference (3DS1, and 2DS1/2DS5 on cB telomeric
  haplotypes), and those paralogs' reads land on the 3DL1/2DS4 loci.
- Conclusion: reference-based short-read quantification can validate presence/absence only for genes whose paralogs are in the
  reference. For KIR, expression has to be quantified against each person's own assembled haplotypes, the same design as 53.
