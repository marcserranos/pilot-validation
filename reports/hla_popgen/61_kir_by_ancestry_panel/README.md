# 61 — KIR gene content and alleles by ancestry (Figure 1 panel) — DRAFT

Script: `scripts/hla_popgen/61_kir_by_ancestry_panel.py`. **Everything in this folder is a draft**
built on a laptop from committed PUBLIC aggregates (`--mode local`). No new data was read. The
real panel comes from `--mode vm` (per-person KIR GTFs, complete haplotypes only, after the
coverage filter). VM commands: `sprints/S05_qc_phasing_kir_trb/vm_parts/61.sh`.

Figure: `fig_kir_by_ancestry_panel_DRAFT.png` / `.pdf`. Tables: `draft_*.tsv`. Every row has a
`basis` column that says where its number came from.

## What the panel shows

| panel | content | draft source | exact or estimated |
|---|---|---|---|
| a | % of complete haplotypes carrying each of the 17 KIR genes, per ancestry | 51 `gene_content_haplotypes.tsv` | estimated, with bounds |
| b | top 3 alleles (protein field, e.g. `KIR2DL1*003`) of 6 genes, % of the gene's calls | 43 `kir_allele_freq.tsv` | pooled ancestries, known calls, all haplotypes |
| c | cB share (rest cA) per ancestry, 95% Wilson CI | 51 `cB_share_by_completeness.tsv` | exact |
| d | distribution of distinct KIR genes per complete haplotype | 51 `gene_content_haplotypes.tsv` | estimated (renormalised) |

Complete haplotype (51's rule): every KIR call of the haplotype on one contig, and both ends of
the cluster (KIR3DL3 and KIR3DL2) present. In the full cohort 75% of haplotypes pass.
Ancestry is the predicted genetic ancestry; people are unrelated (kinship < 0.0442).

## How the draft numbers are derived

- **a.** 51 publishes the 27 gene-content strings carried by ≥ 20 complete haplotypes and their
  % in each ancestry. For gene *g* and ancestry *A*: lower bound = sum of the % of the strings that
  carry *g*. Not every haplotype is covered: small strings are pooled (352 haplotypes) and
  ancestry cells with 1–19 haplotypes are blank. The uncovered share is 4.0–8.8% depending on
  ancestry (AFR 4.4, AMR 4.5, EAS 5.5, EUR 4.0, MID 8.8, SAS 6.8), so upper bound = lower bound +
  that share. The plotted value divides the lower bound by the covered share, which assumes the
  uncovered haplotypes have the same gene make-up as the covered ones. KIR3DL3 and KIR3DL2 are
  100% by definition. `all_haps_presence_pct_43` gives 43's value on all haplotypes for
  comparison.
- **b.** 43 publishes full-resolution named alleles (≥ 20 haplotypes) of known (exact IPD-KIR)
  calls only, pooled over ancestries and over all haplotypes. The draft sums them to the protein
  field. Alleles that 43 pooled as `<20 each` cannot be assigned to a protein, so each protein
  share is a lower bound. Known calls are 2–73% of each gene's calls (43), which is a biased
  subset. The VM run reports every call that keeps its protein (known, genomic-only novel,
  synonymous novel) and pools protein-novel calls separately.
- **c.** Copied from 51 (unrelated people, complete haplotypes).
- **d.** Covered strings only, renormalised to 100%. A bin that would rest only on blank cells is
  left blank (a gap in the line).

## Headline numbers (draft)

- Complete haplotypes, unrelated, with ancestry: 17,868 (AFR 4,374; AMR 3,984; EAS 2,310;
  EUR 4,500; MID 752; SAS 1,948).
- cB share on complete haplotypes: EAS 32.0% (30.2–34.0), AFR 38.7%, AMR 42.7%, EUR 44.5%,
  MID 52.0%, SAS 57.1% (54.9–59.3). Restricting to complete haplotypes lowers cB by 0.4–1.6
  points against all haplotypes (51), and the ancestry ranking does not change.
- Gene content (estimates; bounds in `draft_gene_freq_by_ancestry.tsv`):
  - KIR3DS1 / KIR2DL5A: AFR 6.4% against 18–27% in the other groups.
  - KIR2DL5B: EAS 4.0%, AFR 23.5%, SAS 26.9%.
  - KIR2DS2 / KIR2DL2: EAS 13.8% against 26–39%.
  - KIR2DL3: EAS 86.2%, SAS 60.9%.
  - KIR3DL1 / KIR2DS4: AFR 89.7%, SAS 67.7%.
  - KIR2DL4 and KIR3DP1: 94.7–100%. 51 attributes the residual misses to contracted haplotypes,
    not to dropout.
  - Against 43's all-haplotype values, the complete-haplotype estimates are 2–9 points higher for
    the framework and common A-content genes (a fragmented haplotype under-counts genes). For the
    B-content genes they are 0–3 points lower: the estimate assumes the uncovered rare strings look
    like the common ones, and the common ones are mostly cA.
- Genes per complete haplotype: mean 9.35 (EAS) to 9.97 (SAS). The 9-gene class, essentially the
  cA01 content string, makes up 70.7% in EAS, 62.9% in AFR and 47.2% in SAS.
- Top alleles (pooled, known calls): KIR2DL1*003 51.1%, KIR2DS4*001 56.5%, KIR2DL4*001 34.7%,
  KIR3DL2*002 30.8%, KIR3DL1*015 25.4%, KIR3DL3*003 21.2%.

## Caveats

- **Draft.** Panels a and d are estimated from content strings with ≥ 20 haplotypes; panel b is
  pooled over ancestries and is not completeness-restricted. Pairs of genes that are always
  together on common strings (KIR3DL1/KIR2DS4, KIR3DS1/KIR2DL5A, KIR2DP1/KIR2DL1, KIR2DL4/KIR3DP1)
  get identical estimates here. The VM run separates them.
- **Completeness rule.** Only 75% of haplotypes are complete. If the share of complete haplotypes
  differs by ancestry or by content (long cB haplotypes may break more often), the
  complete-haplotype frequencies are biased. The VM run writes `completeness_by_ancestry.tsv` to
  check this. 51's cB share is 0.4–1.6 points lower on complete than on all haplotypes in every
  ancestry, consistent with a small bias of this kind.
- **991 Sequel II people absent.** The Tier-3 `sequel2` people (mostly AFR) were never KIR-called
  (43 Caveats; S05 WS-K). AFR is under-sampled and the AFR numbers come from Revio and Sequel IIe
  only.
- **Coverage filter pending.** Script 59 has not set the regional-depth cutoff yet. The VM run
  takes the 59 exclusion list through `--exclude-people`; without it, low-coverage people stay in.
- **Disclosure.** Counts 1–19 are shown as `<20` and their rates blanked. Alleles are named only
  with ≥ 20 haplotypes. When a masked cell could be recovered as the total minus the shown cells,
  the smallest shown cell is also masked (`<20*`) or folded into `other`. Masked heatmap cells are
  hatched and masked points are gaps. In this draft, no gene-frequency or allele cell is masked.
  Eight genes-per-haplotype bins are blank because they would rest only on masked inputs. All
  inputs are PUBLIC tables already in the repo.
- Ancestry is predicted, not self-reported. MID is the smallest group (752 complete haplotypes).

## VM-mode outputs (`~/s05/results/61/`)

`run_summary.tsv`, `completeness_by_ancestry.tsv`, `gene_freq_by_ancestry.tsv` (haplotype
frequency, plus person-level carrier frequency among people with both haplotypes complete),
`content_by_ancestry.tsv`, `genes_per_hap_by_ancestry.tsv`, `top_alleles_by_ancestry.tsv`
(`protein` and `full_known` resolutions, top 5 per gene × ancestry), and
`fig_kir_by_ancestry_panel.{png,pdf}`. All are PUBLIC aggregates.

## Review note (2026-10-12)

"cB" here is 51's B-content set, which includes telomeric B genes (KIR3DS1, 2DS1, 2DS5, 2DL5A). So the share is "haplotypes with any
B-content gene", not the centromeric-B motif. The S05 report calls it the "B-content share". A strict cB/tB split is a follow-up.
