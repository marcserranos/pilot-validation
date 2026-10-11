# reports/hla_popgen/ — index

One folder per numbered script (`scripts/hla_popgen/NN_*.py` → `reports/hla_popgen/NN_*/`). This
page says which folders are current, which are superseded, and where each paper-facing result
lives. The standing narrative is [`COMPREHENSIVE_REPORT.md`](COMPREHENSIVE_REPORT.md); the
research map is [`NEXT_STEPS_AND_RESEARCH_MAP.md`](NEXT_STEPS_AND_RESEARCH_MAP.md).

Added 2026-10-11 as part of the repo harmonization pass. Folders were not moved or renamed, so
every existing link and script output path still works.

## Paper-facing results (current)

| Folder | What it is | Status |
|---|---|---|
| [`40_figure1_v5/`](40_figure1_v5/) | **Figure 1 v5** (6 panels) | Current Figure 1. Awaiting Cole's sign-off on layout. |
| [`37_dq_g1g2_signed_ld/`](37_dq_g1g2_signed_ld/) | DQ G1/G2 purge, signed phased D′, mask-and-rephase test, LD supplement CSV | Strongest main-text candidate |
| [`39_saturation_by_ancestry/`](39_saturation_by_ancestry/) | Discovery / saturation curves by ancestry | In Figure 1e |
| [`38_hla_a_deletion_validation/`](38_hla_a_deletion_validation/) | HLA-A "deletion" is an artifact (short-read cross-check) | Supplement + QC methods |
| [`43_kir_full_cohort/`](43_kir_full_cohort/) | KIR on the full cohort (12,261 people) | Descriptive result; known gaps (sequel2, 2DL2/2DL3) |
| [`44_kir_recurrence_saturation/`](44_kir_recurrence_saturation/), [`46_kir_vs_hla_catalogue/`](46_kir_vs_hla_catalogue/) | KIR recurrence and KIR-vs-HLA catalogue completeness | Use **v4b** numbers only (v1–v3 invalid) |
| [`29_hla_ld/`](29_hla_ld/), [`30_hla_sv/`](30_hla_sv/), [`31_aa_diversity/`](31_aa_diversity/) | LD by ancestry, structural variation / copy number, peptide-groove diversity | Supplement |
| [`26_qc_relatives_v2/`](26_qc_relatives_v2/) | Call accuracy from twins/relatives (~1 error per 2,500 bases) | Supplement / methods |
| [`32_novel_callouts/`](32_novel_callouts/), [`34_novel_recurrence/`](34_novel_recurrence/), [`24_novelty_by_field/`](24_novelty_by_field/), [`25_noncoding_novelty_paf/`](25_noncoding_novelty_paf/) | Novel-allele callouts, recurrence, novelty by field, non-coding novelty | Feed Figure 1 c/d and main-text callouts |
| [`47_naive_ml_tests/`](47_naive_ml_tests/), [`48_kir_hla_ligand_cooccurrence/`](48_kir_hla_ligand_cooccurrence/) | Confound battery and KIR–HLA ligand co-occurrence | Negative / QC results, supplement |
| [`42_repertoire_baseline/`](42_repertoire_baseline/) | WS6 BenchRep-T style repertoire baseline (1,500-person subsample) | Superseded in practice by Aleix's W0–W5 plan on his branch |
| [`49_cole_share_package/`](49_cole_share_package/) | Description of the data package shared with Cole (2026-09-25) | Reference |
| [`41_kir_scoping/`](41_kir_scoping/) | KIR region scoping and pilot | Reference (led to 43) |

## Superseded (kept for provenance, do not cite)

| Folder | Superseded by |
|---|---|
| `23_fig1_draft_panels/` | `28_figure1/` → … → `40_figure1_v5/` |
| `28_figure1/` (Figure 1 v2) | `40_figure1_v5/` |
| `33_figure1_v3/` | `40_figure1_v5/` |
| `35_figure1/` | `36_figure1/` (says so in its README), then `40_figure1_v5/` |
| `36_figure1/` | `40_figure1_v5/` |

## Earlier exploratory work (scripts 05–21)

`05_figures_frequency/`, `06_figures_structure/`, `07_figures_crosscohort/`,
`08_embedding_compare/`, `09_mutation_topology/`, `10_allele_ancestry_geometry/`,
`11_relatedness_cohort_overlap/`, `13_disease_allele_association/`, `14_manifold_structure/`,
`15_manifold_holdout_validation/`, `16_phasing_mendelian_validation/`,
`17_raw_sequence_divergence/`, `21_hla_manhattan/`. Background and methods context; several
panels were later redrawn in the numbered Figure 1 scripts.

## Loose files at this folder's top level

`novel_alleles.tsv`, `allele_richness.tsv`, `allele_extrapolation.tsv`, `discovery_rate_*`,
`frequency_spectrum*`, `neutral_model_comparison.tsv`, `allele_space_explored.*`,
`03_novel_alleles_report.md`, `04_allele_saturation_report.md`, `recon_report.*` are the outputs of
scripts 00, 03, 04 and 19, written here before the one-folder-per-script convention. They stay in
place because those scripts and later scripts (19, 25, 49) read and write these exact paths.

> **Disclosure note (resolved 2026-10-11, see `context/DECISIONS.md`):** `novel_alleles.tsv` includes
> single-person carrier counts. Marc decided to keep it as is; the `<20` rule applies when numbers
> go into figures, the manuscript or shared slides.
