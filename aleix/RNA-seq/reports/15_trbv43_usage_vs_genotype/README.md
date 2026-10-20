# 15 — TRBV4-3: expressed repertoire vs germline genotype

Script: `scripts/15_trbv43_usage_vs_genotype.py`. It runs off the VM: inputs are committed,
de-identified aggregate tables only (no per-person data). Built 2026-10-20.

## Question

Report 10 found that 27.3% of people never use TRBV4-3 in their expressed TRB repertoire. At
its usage level (1.1% of clonotypes), chance would give ~0.04%, so 10 read this as a germline
absence. Does that agree with a direct DNA genotype?

The DNA side is Marc's script 65 on `main`. TRBV4-3 sits on a ~21.7 kb insertion ("insertion A",
also carrying TRBV3-2 and a second TRBV6-2/6-3 copy) that GRCh38 has only on an ALT contig.
65 counts long reads spanning the breakpoint (chr7:142,346,200) and calls 0, 1 or 2 copies.

## Inputs

| Side | Table | People |
|---|---|---|
| RNA | `reports/10_vj_usage_ancestry/zero_usage_by_ancestry.csv`, `dosage_histogram.csv`, `codeletion_vs_absence.csv` | 7,527 with ≥ 823 TRB reads, relatives included |
| RNA | `reports/09_analysis_table/waterfall.csv` (denominators per ancestry) | same |
| DNA | `inputs/65_genotype_by_stratum_unrelated.tsv` (cut-off genotyper, AF 0.2/0.8) | 11,856 unrelated long-read people, 11,075 called |
| DNA | `inputs/65_af_mixture_fit_unrelated.tsv` (mixture genotyper, summary only) | same, per platform |
| Locus | `inputs/65_trb_probe_mapping.tsv` (insertion content); `inputs/trust4_trb_gene_coords_hg38.tsv` (gene positions, from TRUST4's `hg38_bcrtcr.fa` at 4032f90) | public reference |

DNA tables are copied unchanged from `main` @ e2584e2 (`inputs/PROVENANCE.txt`), masking
included: MID genotype counts are masked there, so MID has no DNA value here.

## Results

| Ancestry | RNA: never used | DNA: 0 copies | RNA / DNA | Absence allele, RNA (implied) | Absence allele, DNA (counted) |
|---|---|---|---|---|---|
| SAS | 43.5% | 55.1% | 0.79 | 0.66 | 0.71 |
| AFR | 33.3% | 38.5% | 0.87 | 0.58 | 0.58 |
| EAS | 30.1% | 37.1% | 0.81 | 0.55 | 0.58 |
| AMR | 25.4% | 30.1% | 0.84 | 0.50 | 0.52 |
| EUR | 16.7% | 19.5% | 0.86 | 0.41 | 0.42 |
| MID | 25.8% | masked | — | 0.51 | — |
| ALL | 27.3% | 33.0% | 0.83 | 0.52 | 0.54 |

- **Same ranking in every group** (5 groups with both values; Pearson r 0.995). RNA is below DNA
  by a near-constant factor, 0.79–0.87.
- **The gap closes at about 0.1% usage (panel d).** 27.3% of people have zero TRBV4-3 clonotypes,
  32.1% have at most 0.10%, and DNA says 33.0% (mixture genotyper, Revio: 33.2%). So most 0-copy
  people have no TRBV4-3 calls, and the rest carry a handful. These are plausibly reads from
  similar TRBV4 genes assigned to TRBV4-3. A ≤ 0.1% rule would match DNA at the group level. It
  is not validated per person.
- **In allele terms the two sides agree more closely (panel c):** RNA 0.52 vs DNA 0.54 overall,
  within 0.02 in AFR, AMR and EUR. The RNA value assumes Hardy–Weinberg
  (√ of the excess-zero fraction). The DNA value counts alleles, which is robust to 65's known
  het deficit (HWE p 1e-31 to 0.002: some heterozygotes called homozygous). Under HWE, DNA's
  allele frequency implies 28.7% 0-copy people, close to RNA's 27.3%. Part of the 0.83 ratio may
  therefore be on the DNA side. The two can't be separated at group level.
- **Panel e, other TRBV genes in non-users vs users** (genes with median usage ≥ 0.3%):
  - **TRBV6-2** halves (0.47). Insertion A carries a second TRBV6-2/6-3 copy, so 0-copy people
    have one copy instead of two. This is a mechanistic check of the genotype reading.
  - **TRBV7-2**, the first gene 3′ of the insertion, drops to 0.59. It is on GRCh38 primary,
    so it is not deleted. Unexplained. Candidates: an allele in linkage with the 0-copy
    haplotype that is expressed less or named differently, or a local effect on recombination.
  - **Genes 5′ of the insertion** (TRBV2 to TRBV6-1) are used 19–33% more (median ratio 1.21).
    TRBV4-2, next to the breakpoint, is +8%. Most genes 3′ of the insertion change by at most a
    few percent, as expected when ~5% of the repertoire is redistributed (median ratio 1.04).
    Exceptions: TRBV10-2 (1.35) and TRBV11-2 (1.15), ~80 kb 3′.
  - The 5′ rise is an observation, not an explained effect. Candidates: alleles in linkage
    with the 0-copy haplotype, or recombination favouring the 5′ cluster when the insertion is
    absent.

## What this does and does not show

- Group-level only. The RNA cohort (7,527) and DNA cohort (11,856 unrelated) overlap (8,327
  people have both assays) but are not identical, and no person is matched here. The
  person-level test (does each RNA non-user have 0 copies?) is Marc's script 60, which needs
  our per-person usage table (VM, queued).
- RNA includes relatives; DNA is unrelated only. Ancestry labels are the same AoU predictions
  on both sides.
- Panel e compares medians of non-users vs users, with no adjustment for age, ancestry or
  depth. It is descriptive.
- Panel d's curve ends below 100% because report 10 masks histogram bins under 20 people and
  stops at 4.15% usage.

## Files

- `fig_trbv43_usage_vs_genotype.{png,pdf,svg}`: panels a–e.
- `comparison_by_ancestry.csv`: the table above, with 95% Wilson intervals (implied-allele
  interval by the same transform).
- `cumulative_usage.csv`: panel d.
- `codeletion_by_position.csv`: panel e, with GRCh38 position and side of the insertion.
- `summary.csv`: headline numbers.

Rebuild: `python3 scripts/15_trbv43_usage_vs_genotype.py` (needs matplotlib, pandas, scipy,
cnsplots; Liberation Sans for the journal font).
