# Mutation topology: where novel differences land along the CDS (`09_mutation_topology.py`)

Methodology: scripts/hla_popgen/research/VIZ_LIT.md addendum (2026-09). Each stem = one codon position; height = number of haplotype-level novel-difference observations at that position (recurrence-weighted, not a per-cluster count -- see module docstring). Position comes from `cds_mut`'s aa-diff field, first tied candidate only.


## Top-5 hottest codons per gene

| gene | codon | n observations | dominant novelty_class |
|---|---|---|---|
| HLA-A | 207 | 69 | protein_altering |
| HLA-A | 14 | 41 | protein_altering |
| HLA-A | 208 | 31 | protein_altering |
| HLA-A | 42 | 27 | protein_altering |
| HLA-A | 6 | 19 | protein_altering |
| HLA-B | 207 | 74 | protein_altering |
| HLA-B | 14 | 58 | protein_altering |
| HLA-B | 119 | 39 | protein_altering |
| HLA-B | 6 | 29 | protein_altering |
| HLA-B | 180 | 28 | protein_altering |
| HLA-C | 44 | 39 | protein_altering |
| HLA-C | 79 | 29 | protein_altering |
| HLA-C | 72 | 28 | protein_altering |
| HLA-C | 86 | 26 | protein_altering |
| HLA-C | 174 | 23 | protein_altering |
| HLA-DPA1 | 117 | 35 | protein_altering |
| HLA-DPA1 | 51 | 23 | protein_altering |
| HLA-DPA1 | 253 | 23 | protein_altering |
| HLA-DPA1 | 257 | 23 | protein_altering |
| HLA-DPA1 | 199 | 15 | protein_altering |
| HLA-DPB1 | 9 | 111 | protein_altering |
| HLA-DPB1 | 71 | 44 | protein_altering |
| HLA-DPB1 | 129 | 44 | protein_altering |
| HLA-DPB1 | 231 | 37 | protein_altering |
| HLA-DPB1 | 191 | 35 | protein_altering |
| HLA-DQA1 | 9 | 19 | protein_altering |
| HLA-DQA1 | 19 | 18 | protein_altering |
| HLA-DQA1 | 179 | 15 | protein_altering |
| HLA-DQA1 | 49 | 14 | protein_altering |
| HLA-DQA1 | 79 | 14 | protein_altering |
| HLA-DQB1 | 177 | 40 | protein_altering |
| HLA-DQB1 | 76 | 29 | protein_altering |
| HLA-DQB1 | 9 | 24 | protein_altering |
| HLA-DQB1 | 52 | 20 | protein_altering |
| HLA-DQB1 | 36 | 18 | protein_altering |
| HLA-DRB1 | 73 | 45 | protein_altering |
| HLA-DRB1 | 233 | 28 | protein_altering |
| HLA-DRB1 | 113 | 22 | protein_altering |
| HLA-DRB1 | 137 | 21 | protein_altering |
| HLA-DRB1 | 24 | 21 | protein_altering |

Figure: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/09_mutation_topology/lr/09_mutation_topology.png`
