# 53 — Allele-specific expression (ASE) pilot on personal two-haplotype references

Script: `scripts/hla_popgen/53_hla_ase_pilot.py` (VM 2026-10-06, 80 unrelated people with LR assembly + AoU RNA-seq, 0 failures).
Per person: HLA-A/B/C mRNAs spliced out of both assembled haplotypes (+ all other MHC genes as decoys); RNA read 1 over
chr6:29.6–31.6 Mb aligned with minimap2 -ax sr; a read is allele-informative if it aligns to one haplotype's mRNA with MAPQ ≥ 10.
Mapping bias = the same pipeline on 4,000 reads simulated equally from both haplotypes. Table `ase_by_gene_class.tsv`, figure `fig_ase.png`.

| gene | class | people | median informative reads | median abs log2(h1/h2) observed | simulated bias | % imbalanced > 1.41× (p < 0.01) |
|---|---|---|---|---|---|---|
| A | protein-different | 63 | 27,217 | 0.37 | 0.02 | 31.7 |
| B | protein-different | 62 | 51,452 | 0.33 | 0.01 | 32.3 |
| C | protein-different | 42 | 43,842 | 0.33 | 0.01 | masked (<20 imbalanced) |
| A/B/C | protein-identical, mRNA-different | <20 with ≥30 informative reads each gene | — | — | — | — |

- For protein-different heterozygotes the assay works: tens of thousands of informative reads per gene, simulated mapping
  bias far below the observed imbalance, and about a third of people show significant >1.41× allelic imbalance at HLA-A and HLA-B.
- The protein-identical class fails as designed. Most such pairs differ only in introns, so their mRNAs are identical. Where the
  mRNAs do differ, there are only a few distinguishing sites, and MAPQ-based assignment is dominated by bias (the simulated bias
  is as large as the observed signal). This class needs site-level allele counting at the distinguishing positions.
- Limitations: no STAR-like spliced simulation, PCR duplicates not removed, read 1 only, no per-person checkpointing, n = 80.
