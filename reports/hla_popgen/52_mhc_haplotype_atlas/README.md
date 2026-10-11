# 52 — MHC haplotype atlas (absent haplotypes), AMR equal-N richness decomposition, ASE census

Script: `scripts/hla_popgen/52_mhc_haplotype_atlas.py` (VM 2026-10-06, ~2 min). Unrelated people, copy-1 calls,
same-contig pairs only.

## a. Absent-haplotype atlas — `absent_cells_field{1,2}.tsv`, `absent_pair_summary_field{1,2}.tsv`, `fig_absent_cells.png`
Cell = (lineage of gene A, lineage of gene B) on one contig. Absent = observed 0 and exp(−E) < 0.05/m (E = independence
expectation, m = cells with E ≥ 1 in that pair × stratum), in ≥2 ancestries and pooled.
- DQA1~DQB1: all 12 G1×G2 lineage cells absent in every ancestry tested (E up to 1,549). One within-G1 cell
  (DQA1*06~DQB1*02, E = 77) is also absent.
- DRB1~paralog (DRB3/4/5 on the same contig as DRB1): the DR51/52/53 rules come out unprompted (22 absent cells, most in all
  6 ancestries).
- Other absences: DPA1~DPB1 (e.g. DPA1*01~DPB1*09, E = 141, 4 ancestries), DR–DQ (e.g. DRB1*11~DQA1*02, E = 281, 6 ancestries),
  B~C (e.g. B*08~C*06). Most likely deep haplotype LD; a pairing constraint is not ruled out for DP.
- Next: separate LD-absence from incompatibility (2-field resolution, cross-ancestry consistency, interface residues).

## b. Why AMR > AFR — `pooled_richness.tsv`, `amr_source_decomposition.tsv`, `fig_equal_n_richness.png`, `fig_amr_sources.png`
Equal N = 2,281 haplotypes per pool (strict p ≥ 0.95), 50 resamples. Two-field richness:

| gene | AFR | EUR | EAS | AMR | synthetic 45% EUR / 10% AFR / 45% EAS |
|---|---|---|---|---|---|
| A | 45.8 | 39.0 | 45.4 | 62.5 | 58.2 |
| B | 80.0 | 65.0 | 80.4 | 107.7 | 103.4 |
| C | 45.9 | 38.0 | 39.1 | 47.1 | 50.3 |
| DRB1 | 44.5 | 45.0 | 53.5 | 55.9 | 61.3 |
| DPB1 | 51.1 | 30.0 | 34.4 | 50.6 | 46.0 |

A pool of source haplotypes reaches 93–96% of AMR richness at A and B and exceeds it at C and DRB1, so admixture pooling explains most of the AMR > AFR excess. Caveats: the 45/10/45 weights are arbitrary (not fitted to AMR p_* means) and EAS is only a proxy; a fitted-weight sweep is the next check.
EAS stands in for the Indigenous-American source; AoU labels give no clean source panel. Share of AMR alleles not found in the
AFR, EUR or EAS pools at equal N: A 38%, B 33%, C 23%, DQB1 21%, DRB1 19%, DPB1 18%.

## c. ASE census — `ase_power_census.tsv`
Unrelated heterozygotes whose two alleles share the 2-field (protein) name but differ in sequence, and how many have AoU RNA-seq
(v9 multiomics manifest, 8,980 samples): A 115 (82), B 160 (98), C 423 (279), DRB1 579 (378), DQA1 663 (428),
DQB1 514 (344), DPB1 1,186 (807), DPA1 4,188 (2,937). Feeds the 53 pilot.
