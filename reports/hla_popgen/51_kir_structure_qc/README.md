# 51 — KIR structural QC: haplotype completeness, framework genes, allelic-pair co-calls

Script: `scripts/hla_popgen/51_kir_structure_qc.py` (VM 2026-10-06, 12,261 KIR output dirs, 24,513 haplotypes, ~3 min).
Reads only `~/pipeline_outputs_kir/<pid>/immuannot_output/hap{1,2}.gtf.gz`.

## Haplotype classes — `hap_class_counts.tsv`
complete (every KIR call on one contig, 3DL3 and 3DL2 both present) 18,497 (75%); multi-contig 4,189 (17%);
one contig with an open end 1,704 (7%); empty 123.

## Framework genes — `framework_by_completeness.tsv`, `fig_framework_by_completeness.png`
| class | 3DL3 | 2DL4 | 3DP1 | 3DL2 |
|---|---|---|---|---|
| all haplotypes | 96.0 | 93.3 | 93.9 | 93.7 |
| complete | 100 | 96.8 | 96.8 | 100 |
| multi-contig | 94.9 | 86.5 | 89.4 | 93.8 |
| one contig, open end | 62.1 | 78.6 | 80.6 | 31.6 |

Fragmentation explains part of the 93–96% framework presence. On complete haplotypes, a follow-up check
(`~/s05/inspect_kir.py`, VM) showed that most remaining misses are a joint absence of 3DP1 + 2DL4 + 3DL1/S1 on haplotypes
that otherwise carry more than 6 KIR genes. That pattern fits the known contracted KIR haplotypes (central-region deletion),
not random dropout. It also explains why report 43 found most framework misses "flanked" by called genes.

## Allelic-pair co-calls — `allelic_pair_events.tsv`, `fig_allelic_pairs.png`
- 2DL2 + 2DL3 on one haplotype: 0.81% of haplotypes carrying either. About 90% of events have the two genes on different
  contigs, none overlap in coordinates, and on complete haplotypes the event count is <20 of 18,328 → artefact of haplotype
  assignment across contig breaks.
- 3DL1 + 3DS1: 1.82%; 64% on one contig as distinct loci (median spacing 33.8 kb, about one KIR gene), 36% on different
  contigs. On complete haplotypes: 1.22%, all distinct non-overlapping loci → consistent with real 3DL1/S1 tandem-duplication
  haplotypes; 0.9–1.7% of complete haplotypes in every ancestry.

## Gene-content haplotypes — `gene_content_haplotypes.tsv`, `cB_share_by_completeness.tsv`
Content strings (17-gene presence/absence, cen→tel order) with ≥20 haplotypes, frequency by ancestry, on complete haplotypes;
cB share on all vs complete haplotypes. Physical-phasing matrix for KIR: `kir_same_contig_matrix.tsv`.

## Rule proposed
KIR haplotype-level statements only on complete haplotypes; allelic-pair co-occurrence only when same-contig.
