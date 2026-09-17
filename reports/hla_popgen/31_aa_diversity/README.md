# 31 — amino-acid diversity along the HLA protein, and between-ancestry differentiation

*Supervisor asks A8 (diversity at the amino-acid level, peptide-binding vs non-binding) and A10 (is class II more differentiated than class I), 2026-09-17.*

## Method in one paragraph

Every haplotype's call is reduced to its two-field protein identity and looked up in the IPD-IMGT reference shipped with Immuannot (3.55.0). Per-residue amino-acid diversity is the unbiased expected heterozygosity over those proteins, weighted by how often each was observed. No per-person sequence scan is needed. **11855** unrelated people; **2.70%** of calls had no resolvable protein identity (novel at field 1 or 2, or undetermined) and were excluded.


The groove is defined by **IMGT's own exon annotation** — exons 2+3 for class I, exon 2 for class II — not by a published residue list, because the canonical ARS residue tables could not be verified from primary sources (`WS_literature_selection.md`). A structural peptide-contact definition was not available for this run, so only the exon-level comparison is reported.


## 1. How much of the cohort this covers

| gene   |   n_haplotype_calls_disp |   pct_protein_in_reference |   n_distinct_proteins |
|:-------|-------------------------:|---------------------------:|----------------------:|
| A      |                    22497 |                     99.849 |                   137 |
| B      |                    22485 |                     99.933 |                   211 |
| C      |                    22360 |                     99.942 |                   116 |
| DPA1   |                    22830 |                     99.93  |                    37 |
| DPB1   |                    22470 |                     99.835 |                    91 |
| DQA1   |                    22942 |                     99.974 |                    39 |
| DQB1   |                    22697 |                     99.833 |                    52 |
| DRB1   |                    22479 |                     99.702 |                    89 |
| DRB3   |                     9598 |                     99.875 |                    10 |
| DRB4   |                     6156 |                     99.854 |                     3 |
| DRB5   |                     3423 |                     97.897 |                     6 |
| E      |                    22189 |                     99.986 |                    18 |
| F      |                    22372 |                     99.933 |                     9 |
| G      |                    22617 |                     96.613 |                    13 |
| DRA    |                    22591 |                    100     |                     7 |


`pct_protein_in_reference` below ~99% for a gene means this analysis is missing real diversity there — check it before interpreting that gene.


## 2. Groove vs non-groove amino-acid diversity

| gene   |   protein_length |   n_distinct_proteins |   n_haplotypes_disp |   mean_pi_aa |   mean_pi_groove |   mean_pi_rest |   ratio_groove_vs_rest |   p_permutation | is_conserved_control   |
|:-------|-----------------:|----------------------:|--------------------:|-------------:|-----------------:|---------------:|-----------------------:|----------------:|:-----------------------|
| A      |              365 |                   137 |               22463 |     0.063194 |         0.087559 |       0.038961 |                 2.2474 |         0.001   | False                  |
| B      |              362 |                   211 |               22470 |     0.058908 |         0.091706 |       0.025746 |                 3.562  |         0.0005  | False                  |
| C      |              366 |                   116 |               22347 |     0.0516   |         0.058103 |       0.045167 |                 1.2864 |         0.17641 | False                  |
| DPA1   |              260 |                    37 |               22814 |     0.014141 |         0.021775 |       0.010624 |                 2.0495 |         0.14943 | False                  |
| DPB1   |              258 |                    91 |               22433 |     0.031946 |         0.075496 |       0.009403 |                 8.0293 |         0.0005  | False                  |
| DQA1   |              255 |                    39 |               22936 |     0.065547 |         0.136377 |       0.031367 |                 4.3478 |         0.0005  | False                  |
| DQB1   |              261 |                    52 |               22659 |     0.080563 |         0.142577 |       0.047924 |                 2.9751 |         0.0005  | False                  |
| DRB1   |              266 |                    89 |               22412 |     0.068755 |         0.13303  |       0.035887 |                 3.7069 |         0.0005  | False                  |
| DRB3   |              266 |                    10 |                9586 |     0.023967 |         0.054001 |       0.008609 |                 6.2724 |         0.001   | False                  |
| DRB4   |              266 |                     3 |                6147 |     0.001211 |         4.7e-05  |       0.001806 |                 0.026  |         0.56322 | False                  |
| DRB5   |              266 |                     6 |                3351 |     0.011906 |         0.026792 |       0.004294 |                 6.239  |         0.0015  | False                  |
| E      |              358 |                    18 |               22186 |     0.001481 |         0.00278  |       0.000139 |                20.0398 |         0.3918  | True                   |
| F      |              346 |                     9 |               22357 |     0.000716 |         5.4e-05  |       0.00145  |                 0.0375 |         0.86457 | True                   |
| G      |              338 |                    13 |               21851 |     0.001639 |         0.002466 |       0.000674 |                 3.6566 |         0.25137 | False                  |
| DRA    |              254 |                     7 |               22591 |     0.001975 |         0        |       0.002917 |                 0      |         0.93453 | True                   |


The permutation test shuffles residue positions within the gene, so it holds the gene's diversity distribution fixed and only breaks the association with the groove. The conserved control genes (DRA, HLA-E, HLA-F) are the check: they should show little or no enrichment.


## 3. Between-ancestry differentiation (ask A10)

| gene   |   n_pops |   min_n_disp |     hs |     ht |    fst |   gst_prime |
|:-------|---------:|-------------:|-------:|-------:|-------:|------------:|
| A      |        6 |          579 | 0.9079 | 0.9384 | 0.0325 |      0.4173 |
| B      |        6 |          579 | 0.954  | 0.9735 | 0.02   |      0.5184 |
| C      |        6 |          571 | 0.916  | 0.9365 | 0.0219 |      0.309  |
| DPA1   |        6 |          590 | 0.5281 | 0.6075 | 0.1307 |      0.3063 |
| DPB1   |        6 |          580 | 0.8288 | 0.8921 | 0.0709 |      0.483  |
| DQA1   |        6 |          588 | 0.8763 | 0.902  | 0.0285 |      0.2704 |
| DQB1   |        6 |          580 | 0.8844 | 0.9125 | 0.0308 |      0.3136 |
| DRB1   |        6 |          574 | 0.9326 | 0.9567 | 0.0252 |      0.4437 |
| DRB3   |        6 |          284 | 0.5871 | 0.6191 | 0.0516 |      0.1397 |
| DRB4   |        6 |          167 | 0.2564 | 0.3004 | 0.1463 |      0.2069 |
| DRB5   |        5 |          306 | 0.3601 | 0.431  | 0.1646 |      0.2804 |
| E      |        6 |          568 | 0.5092 | 0.5143 | 0.01   |      0.0225 |
| F      |        6 |          576 | 0.2255 | 0.234  | 0.0364 |      0.0491 |
| G      |        6 |          554 | 0.484  | 0.4977 | 0.0276 |      0.0587 |
| DRA    |        6 |          587 | 0.4671 | 0.4737 | 0.0141 |      0.0289 |


`hs` is within-ancestry heterozygosity. At HLA it is close to 1, which mechanically bounds raw `fst` near zero however different the populations are — this is the Brandt et al. 2018 (G3) point. **Compare genes on `gst_prime`**, Hedrick's standardised measure, not on `fst`.


## Files

- `aa_diversity_per_residue.tsv` — per-residue pi, exon, groove flag, contact flag.

- `groove_vs_rest_summary.tsv`, `allele_differentiation.tsv`, `protein_reference_coverage.tsv`.

- `fig_protein_track_<gene>.png` — the amino-acid Manhattan per gene.

- `fig_groove_enrichment.png`, `fig_allele_differentiation.png`.


## Caveats

- Catalogued proteins only. Novel proteins found in S01 are not in the reference and are excluded; see §1 for how much that is.

- Residue numbering is the reference protein's own, **including the signal peptide**. It is not mature-protein numbering, so these numbers are not directly comparable to a published residue list without an offset.

- Frequency weighting means a gene dominated by a few common alleles gets a lower pi than its allele *count* suggests. That is intended — it is the diversity a randomly drawn haplotype actually sees.

- Differentiation is computed on protein-level alleles, not SNPs, so it is not directly comparable to genome-wide Fst distributions.
