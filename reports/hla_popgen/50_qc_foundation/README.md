# 50 — QC foundation (HLA): exclusion waterfall, physical-phasing levels, platform effect, novel-protein funnel, DQ ascertainment

Script: `scripts/hla_popgen/50_qc_foundation.py` (run on the VM 2026-10-06, full cohort, ~1 min; PUBLIC disclosure).
Inputs: `hla_calls_rich.tsv` (892,856 rows, 12,233 people), `cohort_membership.tsv`, AoU relatedness table.
Unrelated = greedy removal at kin ≥ 0.0442 (same rule as 24/39) → 11,856.

## a. Exclusion waterfall — `waterfall.tsv`, `fig_waterfall.png`

| step | total | AFR | AMR | EAS | EUR | MID | SAS |
|---|---|---|---|---|---|---|---|
| LR cohort rows | 13,252 | 4,053 | 2,807 | 1,499 | 3,109 | 499 | 1,262 |
| HLA called (Revio 11,015 + Sequel IIe 1,218) | 12,233 | 3,108 | 2,777 | 1,494 | 3,078 | 495 | 1,258 |
| unrelated | 11,856 | 3,020 | 2,656 | 1,460 | 2,974 | 487 | 1,236 |
| + ancestry p ≥ 0.90 | 7,924 | 2,324 | 1,582 | 1,270 | 1,293 | 304 | 1,151 |
| + ancestry p ≥ 0.95 | 7,119 | 2,004 | 1,389 | 1,252 | 1,141 | 184 | 1,149 |
| + ancestry p ≥ 0.98 | 6,136 | 1,486 | 1,214 | 1,226 | 965 | 100 | 1,145 |

- All 991 `sequel2` people have no HLA calls (Tier-3 trim path broken); they are overwhelmingly AFR, so AFR is under-sampled.
- Purity thresholds are not ancestry-neutral: p ≥ 0.95 keeps 93% of SAS and 86% of EAS but 38% of EUR and 38% of MID
  (mean p_eur among EUR-labelled = 0.81). Any "strict ancestry" result changes composition; report N per threshold.
- This table is not yet the agreed exclusion list: it has no long-read sample-QC flag (source undefined) and no repertoire filter.

## b. Physical-phasing level — `same_contig_matrix.tsv`, `mhc_contigs_per_hap.tsv`, `fig_same_contig_matrix.png`

% of haplotypes (both genes called, copy 1) with both genes on one contig: DPA1~DPB1 100, DQA1~DQB1 98, DRB1~DQA1 96,
B~MICA 95, DRB1~DQB1 94, B~C 91; A~B 31, A~C 32; class I ~ class II 19–29. Mean number of contigs carrying MHC calls per
haplotype: 4.3 (Revio), 4.9 (Sequel IIe). Haplotype-level claims are robust within class II and within B–C–MICA, not across
the whole MHC.

## c. Platform effect on novelty — `novelty_by_platform.tsv`, `novelty_platform_mh_or.tsv`, `fig_platform_novelty_or.png`

Mantel–Haenszel OR (Revio vs Sequel IIe, strata = 6 ancestries), classical genes: protein-novel 0.50–0.71, non-coding novel
0.63–1.00. Same ancestry, older platform → more "novel" calls. Platform must be a covariate in every novelty claim.

## d. Protein-novel funnel — `protein_novel_waterfall.tsv`, `fig_protein_novel_funnel.png`

Distinct sequence key = (gene, `cds_mut`), unrelated people.

| | distinct | ≥2 unrelated carriers | not homopolymer-indel-only | no template warning | ≥2 platforms | ≥20 carriers, clean |
|---|---|---|---|---|---|---|
| classical 8 | 4,421 | 174 | 8 | 7 | 1 | 0 |
| TAP1/TAP2/MICA/MICB | 4,066 | 364 | 137 | 127 | 55 | 23 |
| all 40 genes | 18,986 | 1,677 | 661 | 276 | 115 | 37 |

Recurrent protein-novel sequences at classical loci are almost all homopolymer-indel-only (the known HiFi artefact class).
Singletons cannot be separated from assembly error without read-level evidence. A classical-HLA protein-novel headline is not
supported as it stands. Reference used by Immuannot: IPD-IMGT/HLA 3.55.0 (GTF header); a newer release has not been checked.

## e. DQ G1/G2 ascertainment — `dq_ascertainment.tsv`, `fig_dq_ascertainment.png`

- Contig-break rate between DQA1 and DQB1: 2.35% (G1/G1, n = 13,287 haplotypes), 1.31% (G2/G2, n = 9,128).
- Cross-group pairings (G1/G2: 30 haplotypes, G2/G1: 34) are 100% and 94% split across contigs.
- On split haplotypes, cross-group pairing is 12.6% against 42.4% expected under random re-pairing: hap labels are partly but not
  reliably carried across contig breaks.
- Same-contig cross-group haplotypes: <20 of 21,985. A follow-up check (`~/s05/inspect_dq2.py`, VM) found these coincide exactly
  with the only same-contig haplotypes whose DQA1–DQB1 spacing is ≥ 100 kb (98% of the 22,686 same-contig haplotypes lie in
  16.6–24.2 kb) and carry a second DQB1 copy, i.e. structurally anomalous contigs.
- Result: no cross-group DQ haplotype among structurally canonical, physically phased haplotypes. This is an unfiltered
  re-derivation and agrees with 37 (0 / 17,255 under 37's filters).
