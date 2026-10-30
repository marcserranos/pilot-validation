# `hla_popgen` — HLA population genetics, admixture, and novel-allele discovery

> **Role:** the second-generation analysis of the full-cohort Immuannot long-read run — everything
> the first pass (`scripts/production_analysis/`) left on the table.
> **Read:** `SCHEMA.md` before touching any script here. It is the binding data contract.

## Why this exists

The production run called HLA on ~12,000 All of Us participants from phased long-read assemblies.
The first round of post-processing (`scripts/production_analysis/`, results in
`reports/pre_sprint/full_immuannot_lr_calling/`) answered "did the pipeline work, and does HLA structure track
ancestry." It did that on **8 classical genes, discrete ancestry labels, one cohort, and a single
string per gene**.

A full read of Immuannot's upstream source (archived: `reference/IMMUANNOT_GTF_SPEC.md`) established
that **the raw output on disk is far richer than anything we extracted**, and that nothing here
requires re-running the expensive pipeline. Three specific gaps:

1. **Contig identity is discarded by every existing parser.** A `hap1.gtf.gz` can legitimately
   contain genes from more than one physical contig when the assembly is fragmented across
   chr6:29.5–33.5Mb. Two genes sharing a hap file are **not** necessarily in cis. The only valid
   cis key is `(person_id, hap, contig)` — GTF column 1, which nothing in the repo reads. Any
   DQA1~DQB1 heterodimer built from the current tables would be silently wrong.
2. **57 of 65 genes were never analysed.** Non-classical HLA (E/F/G), pseudogenes, DRB paralogs,
   DQA2/DQB2, MICA/MICB/TAP1/TAP2 and C4 are all present in the raw GTFs and dropped by a hardcoded
   8-gene list downstream.
3. **The novelty signal is richer than a boolean, and the observed sequences are on disk.**
   Immuannot writes each gene copy's actual observed CDS to `hap{N}/cds.fa.gz` — a file its own
   cleanup step does *not* delete — and encodes novelty severity in the field depth at which it
   splices `"new"` into the allele string. `cds_mut` carries the exact codon-level difference.
   None of this was captured. **Novel-allele discovery is therefore fully possible from existing
   data, with no re-run.**

## The three cohorts

Every frequency and structure figure is reproducible across three cohorts, so that differences
*between* them are the finding:

| Cohort | N | Resolution | Phased? | Genes |
|---|---|---|---|---|
| `sr` — AoU-native short-read | ~500,000 | 2-field | no | 8 classical |
| `lr` — Immuannot long-read | ~12,000 | up to 4-field | **yes** | 65 |
| `lr_td<N>` — distance-filtered | sweep | up to 4-field | yes | 65 |

`template_distance` (edit distance to the nearest documented IPD allele) is treated as a **sweep,
never a fixed threshold**. How a statistic moves as the filter tightens is itself the result — and
the distance distribution by ancestry is direct evidence of reference bias in IPD-IMGT/HLA.

## Pipeline

Scripts are numbered in the order the story builds, and the number is unique — never reuse one.
`reports/hla_popgen/<NN>_<name>/` mirrors the script that wrote it. A gap in the numbers is fine;
a duplicate is not (see COMPREHENSIVE_REPORT.md for the narrative these map onto).

**Foundation — extraction and cohorts**
```
00_recon_vm.py        → reports/hla_popgen/recon_report.{md,json}   [RUN THIS FIRST]
00b_warning_census.py → template_warning census; sets 03's --disqualifying-warnings policy
01_extract_rich.py    → hla_calls_rich.tsv (Table 1), hla_cis_pairs.tsv (Table 2)
02_build_cohorts.py   → cohort_membership.tsv (Table 4)
```

**Allele discovery and frequency**
```
03_novel_alleles.py       → novel_alleles.tsv (Table 3), novel_alleles_seqs.fa [VM-only]
04_allele_saturation.py   → discovery curves, Chao2/ACE richness, "% of allele space found"
05_figures_frequency.py   → allele frequency by ancestry, diversity indices
06_figures_structure.py   → PCA/UMAP, Fst, continuous-admixture representations
07_figures_crosscohort.py → SR-vs-LR bias, template_distance by ancestry, cis heterodimers,
                            non-classical gene diversity, resolution cascade
08_embedding_compare.py   → allele-space embedding matrix, disease-coloured overlays
09_mutation_topology.py   → where novel differences land along the protein sequence
10_allele_ancestry_geometry.py → per-allele ancestry centroids (ternary / tetrahedron / PCA)
```

**Disease and manifold structure** (12 feeds 13; 14 feeds 15)
```
12_disease_phenotypes.py           → per-person diagnosis labels for the locked UMAP cohort
13_disease_allele_association.py   → Fisher + ancestry-adjusted CMH per allele × disease
14_manifold_structure.py           → variance/loadings audit, ancestry residualisation, supervised UMAP
15_manifold_holdout_validation.py  → held-out test of 14's supervised-UMAP signal (it does not generalise)
```

**Relatedness and phasing validation** (11 → 16 → 17)
```
11_relatedness_cohort_overlap.py      → AoU kinship table × our cohorts, bucketed
16_phasing_mendelian_validation.py    → switch errors across 545 real relative pairs
17_raw_sequence_divergence.py         → raw cds.fa.gz comparison, bypassing allele nomenclature
```

**Saturation follow-ups** (read only aggregate tables 04 wrote — laptop-runnable)
```
18_discovery_rate.py        → dS/dN, expected new alleles per additional haplotype
19_allele_space_coverage.py → richness vs copy coverage, with the full evidence base
```

**Capstone — per-site diversity along the gene** (20 → 21 → 22)
```
20_gene_diversity_track.py  → first-pass positional track (superseded by 21; kept for provenance)
21_hla_manhattan.py         → per-gene Manhattan of per-site pi against one canonical reference
22_panel_overview.py        → cross-gene CDS vs non-CDS comparison, from 21's panel_summary.tsv
```

**Figure 1 drafts** (23 — laptop-only, reads committed aggregates)
```
23_fig1_draft_panels.py     → novel rate by ancestry x gene + SFS known vs novel, split CDS/non-coding/artifact
```

**S01–S04 (2026-09-16 → 09-28) — novelty/QC, supervisor items, Figure 1, KIR** (sprint folders `sprints/S01_…` to `S04_…`; not listed above before the 2026-10-27 cleanup)
```
24_novelty_by_field.py            → Novelty re-defined by nomenclature field, with per-call artifact labels, a sequence-level truth check against the IPD-IMGT CDS snapshot Immuannot ship
25_noncoding_novelty_paf.py       → Split depth-4 ("beyond_cds") novelty into distinct non-coding sequences, using the cs difference string between each person's contig and the call's te
26_qc_relatives_v2.py             → Non-circular QC of the long-read (Immuannot) HLA calls using real relative pairs -- the replacement for 16/17's circular, denominator-free version, pe
27_allele_space_coverage.py       → How much of the HLA allele space have we seen, and what would it take to see more? Spec: sprints/S01_novelty_qc_coverage/WS3_allele_space_coverage.md
28_figure1_compose.py             → Figure 1 for the long-read HLA paper, composed from committed aggregate results.
29_hla_ld_by_ancestry.py          → Linkage disequilibrium between HLA alleles at physically-phased cis gene pairs, within ancestry.
30_hla_structural_variation.py    → Gene deletions, duplications and copy number across the MHC, from phased long-read assemblies.
31_aa_diversity_selection.py      → Amino-acid diversity along the HLA protein, groove vs non-groove, and allele-level differentiation between ancestries.
32_novel_allele_callouts.py       → The main-text callout list: novel HLA alleles that are common in a population and absent from IPD-IMGT/HLA.
33_figure1_v3_compose.py          → Figure 1 v3 — the panel set the supervisors asked for on 2026-09-17.
34_novel_protein_recurrence.py    → Novel protein-altering alleles by recurrence, broken down by gene.
35_figure1_assemble.py            → SUPERSEDED by 36_figure1_native.py — kept only for provenance.
36_figure1_native.py              → Figure 1, drawn natively from data as a single figure.
37_dq_g1g2_signed_ld.py           → Signed phased D' between DQA1 and DQB1, the Petersdorf (2022) G1/G2 "predicted incompatible" cross-quadrant, and the transplant-genetics question: are
37b_ld_supplement_table.py        → LD supplement table + heatmap -- re-rendered from the COMMITTED 29_hla_ld tables only.
37c_dq_g1g2_from_committed.py     → DQ G1/G2 signed-D' figure, purge statistic (O/E), and bimodality check, computed ENTIRELY from the already-committed `reports/hla_popgen/29_hla_ld/ld_
37d_mask_rephase_em.py            → WS1 deliverable 3 (S03 call #8, DQ G1/G2): the mask-and-rephase experiment.
37e_unrelated_fix_kfold_em.py     → S03 critic #1 fix (blocker #2 + major #3): rerun 37's O/E and mask-and-rephase EM restricted to the CORRECT unrelated cohort (matching scripts 29/38's
38_hla_a_deletion_validation.py   → HLA-A 'deletion' validation (WS2, S03 call #8).
38b_deletion_supplement_fig.py    → Deletion supplement figure -- re-rendered from the COMMITTED 30_hla_sv tables only.
39_saturation_by_ancestry.py      → Per-ancestry allele discovery / saturation curves (Pakistan Genome Resource, Nature 2026, Fig 3e style).
40_figure1_v5.py                  → Figure 1 v5 -- the introductory figure of the paper, redrawn to 2026-09-22 panel-by-panel supervisor feedback and then to the orchestrator's 2026-09-2
40a_admixture_bins.py             → Panel (a) source data for 40_figure1_v5.py -- binned, disclosure-safe admixture strip.
40b_novelty_rate_export.py        → Panel (c) source data for 40_figure1_v5.py -- gene x ancestry novelty TOTALS, exported once, suppressed once. MUST run on the VM (recomputes from Tabl
41_kir_pilot.py                   → WS5 (Sprint S03, call #8 summary section 7) -- KIR scoping + pilot.
42_repertoire_baseline.py         → WS6 — BenchRep-T-style repertoire-only disease-classification baseline (VJ-kmer).
43_kir_full_aggregate.py          → WS5 follow-on (S03 call #8 sect 7) -- aggregate the full-cohort KIR Immuannot run.
43b_kir_full_figure.py            → S03 WS6 follow-on -- Nature-grade multi-panel figure for the full-cohort KIR run.
44_kir_recurrence_saturation.py   → WS-A (Sprint S04) -- KIR recurrence + saturation, on the same footing as HLA (WS-B setup).
45_kir_recurrence_figure.py       → WS-A (Sprint S04) -- figures for `44_kir_recurrence_saturation.py`'s aggregate TSVs.
46_kir_vs_hla_catalogue.py        → WS-B (Sprint S04) -- KIR vs HLA reference-catalogue coverage, script 46.
47_naive_ml_tests.py              → S04 WS-D -- naive ML / simple tests on allele-carriage feature matrices (HLA two-field + KIR gene presence), aggregate-only outputs.
47b_naive_ml_figure.py            → WS-D (Sprint S04) -- figure for `47_naive_ml_tests.py`'s aggregate TSVs.
48_kir_hla_ligand_cooccurrence.py → S04 WS-D -- KIR-HLA ligand co-occurrence by ancestry (seed-list item: ORCHESTRATOR_HANDOFF.md sect2 WS-D, "KIR-HLA ligand co-occurrence by ancestry").
48b_ligand_figure.py              → WS-D (Sprint S04) -- figures for `48_kir_hla_ligand_cooccurrence.py`'s aggregate TSVs.
49_cole_share_package.py          → Cole data-share package builder (S04 WS-E).
58b_s05_result_figures.py         → Figures for the S05 VM results (58 regional depth, 64 novel-allele validation), drawn from the committed PUBLIC tables only (no VM, no participant dat
```

**S05 (2026-10-06) — QC foundation, KIR structure, haplotype atlas, ASE pilot** (VM; read existing per-person outputs only)
```
_s05_common.py              → shared loaders (unrelated rule, ancestry upper-casing, HLA- prefix, <20 masking)
50_qc_foundation.py         → exclusion waterfall, MHC same-contig matrix, platform x novelty, protein-novel funnel, DQ ascertainment
51_kir_structure_qc.py      → KIR haplotype completeness classes, framework presence, allelic-pair event typing, content haplotypes
52_mhc_haplotype_atlas.py   → absent-haplotype atlas (DQ G1/G2 + DR51/52/53 controls), AMR equal-N decomposition, ASE census
53_hla_ase_pilot.py         → allele-specific expression of HLA-A/B/C from RNA-seq on personal two-haplotype references
54_s05_figures.py           → laptop-only; figures from the committed 50-53 PUBLIC tables
55_structure_expression_validation.py → assembled copy number vs RNA-SeQC TPM (DRB5, KIR)
56_phasing_typing_benchmark.py → twins / first-degree relatives: typing + haplotype transmission, confidence-score sweep
56b_relative_pair_refinement.py → de-confounded 56 (MHC-sharing pairs; both relatives same-contig)
56c_benchmark_figures.py    → laptop-only; 56/56b figures
56d_score_vs_template_distance.py → 7-point score vs template distance alone, error at matched retention (Marc 10-13)
57_yak_trioeval_pilot.py    → yak trioeval/qv on kinship-inferred trios (MHC, class I/II, KIR, chr20 windows); needs ~/envs/yak
57b_trio_gene_origin.py     → parent of origin per HLA call: same-contig pair phasing truth, Mendelian typing check
57c_trio_figures.py         → laptop-only; 57/57b figures
57d_trio_parent_data.py     → parent platform and depth vs 57 switch/hamming and 57b origin calls (Marc 10-13)
```

**S05 continued (2026-10-11 →)** (sprint folder `sprints/S05_qc_phasing_kir_trb/`)
```
58_regional_depth.py        → samtools bedcov MAPQ0/20: MHC/KIR/TRB depth vs genome-wide control windows; TRB gene copy number + HWE
59_depth_outcome_join.py    → regional depth (58) vs KIR incompleteness, DQ anomalies, HLA absences, MHC fragmentation; cutoff scan + proposed rule
60_trb_germline_usage.py    → TRBV copy number (58 proxy or 65) vs never-used in Aleix's repertoire usage; contingency, diagnostics, dosage
61_kir_by_ancestry_panel.py → Figure 1 KIR panel: gene carrier freq, cA/cB, genes/haplotype, top alleles by ancestry on complete haplotypes (local draft / VM / figure modes)
62_cole_haplotype_units.py  → per-haplotype DQA1~DQB1 and DRB1~DQ units with physical-evidence level, DQ group/anomaly, coverage flag (INTERNAL for the collaborator bucket)
63_sequel2_recall.py        → HLA + KIR calling for the 991 Sequel II people via a fixed Tier-3 self-align (run_immuannot_person.py); new output dirs only, then merge
64_novel_allele_validation.py → novel HLA calls vs a current IPD-IMGT/HLA release (now_known) and HiFi/short-read support per call; homopolymer run-length readout; 50d funnel with evidence
65_trb_assembly_genotype.py → TRB ~21 kb insertion (TRBV4-3/3-2/6-3/7-8) presence per haplotype from phased assemblies; frequency, HWE, co-presence
65c_65b_summary.py          → 65b assembly-vs-reads validation, rates only: completion, failures, concordance (cuts and mixture), HWE by caller
66_cross_table_suppression.py → laptop (and called by 58/65 at the end of a PUBLIC run): complementary <20* masking across committed PUBLIC tables (cross-table sums, waterfall steps, unlisted-ancestry residuals); `--check` audits
67_phasing_error_vs_coverage.py → typing (relative pairs, trios) and phasing-switch error vs per-gene depth; switch topology along the MHC; prune curve vs random removal
67b_pair_error_surface.py   → high-res smoothed pair error over both relatives' gene depth; scatter VM-only (Marc 10-13)
68_rnaseq_composition.py    → AoU blood RNA-seq composition: DRAGEN + RNA-SeQC 2 QC distributions and TPM share by gene family (cohort quantiles only)
```

**S06 (2026-10-20 →) — call-quality metrics vs trio/pair typing error, HLA + KIR** (sprint folder `sprints/S06_call_quality_metrics/`; final report built by 79)
```
70_call_exon_depth.py       → exon-, CDS- and per-base long-read depth of every HLA gene, per person (VM).
71_call_quality_table.py    → one table of per-call quality metrics, the typing-error truth sets, and the depth distributions (VM).
71b_call_quality_stats.py   → how well does each per-call quality metric predict typing error, and what does pruning the worst calls do? (VM)
71c_contig_truth_matrix.py  → contig co-location versus trio truth, for every pair of MHC genes (VM).
71d_support_tables.py       → two small aggregate tables for the report, from 71's INTERNAL tables (VM only).
71e_dropout_outcome.py      → the incomplete-genotype ("dropout") outcome, over the whole cohort (VM).
71f_call_status_audit.py    → what does "one called haplotype" mean, and where do the trio / pair "errors" come from? A read-only audit (VM, aggregate output).
72_call_quality_figures.py  → paper-grade figures for the call-quality sprint, drawn from the AGGREGATE tables of 70, 71, 71b and 71c.
73_depth_schematic.py       → methodology schematics for the call-quality report (public reference coordinates only, simulated reads).
74_build_report.py          → assemble the call-quality report (one self-contained HTML file) from the aggregate tables and the figures.
75_trio_error_anatomy.py    → what is a trio "typing error", exactly? Anatomy of the Mendelian mismatches and four alternative definitions (VM).
76_call_outcomes.py         → the three call outcomes, kept apart (VM; aggregate output + INTERNAL unit tables for 76b).
76b_pruning_evidence.py     → pruning evidence per outcome, with transfer and robustness (VM; reads 76's INTERNAL unit tables, writes aggregate curves).
77_kir_call_table.py        → per-call KIR table for the call-quality replication (VM; reads only the existing Immuannot GTFs, no re-annotation).
77b_kir_trio_origin.py      → per-KIR-gene parental origin in trio children (KIR analogue of 57b; VM).
77c_kir_trio_records.py     → Mendelian-mismatch records for KIR (trio / pair / twin), the KIR analogue of 71's records (VM).
78_kir_gene_depth.py        → exon / CDS depth for the KIR genes that exist on primary chr19 of GRCh38 (VM; the KIR counterpart of 70).
78b_kir_pruning_evidence.py → pruning evidence for the KIR genes (VM; the KIR counterpart of 76b, which supplies all curve / CV / transfer machinery).
79_s06_final_report.py      → the final, standalone S06 report (HLA + KIR), built from AGGREGATE tables only (laptop).
```

**S08 (2026-10-23) — base-level case browser and novel-allele recheck** (sprint folder `sprints/S08_case_browser/`)
```
90_case_browser.py          → base-level case browser for HLA/KIR call errors (Immuannot alignment, Mendelian, switches).
90_case_browser_extract.py  → extract (VM side): pick example calls and pull their sequences into one anonymised cases.json.
90_case_browser_sr.py       → sr (VM side): do the short reads support the long-read variants in a case?
90_case_browser_stats.py    → stats (VM side, aggregate output): how are Immuannot's novel calls made up, over EVERY call in the cohort?
90_case_browser_viewer.html → viewer page for 90_case_browser
91_novel_recurrence.py      → novel-call recurrence (VM side, aggregate output): which 'novel' substitution-only HLA calls recur, in whom, and are they inherited?
92_novel_burden.py          → novel burden by sequencing platform and ancestry (VM side, aggregate output).
93_shortlist_ipd_release.py → re-test the S08 shortlist against a CURRENT IPD-IMGT/HLA release (local; reads the gitignored case_browser/ files).
94_s08_report.py            → S08 final report (HTML, inline SVG figures; no external assets).
95_signature_context.py     → context tests for recurrent novel signatures (VM side; think-tank experiments that need no new sequencing data).
96_f_signature_ase.py       → RNA-seq allele-specific expression of recurrent HLA-F signatures (VM side; think-tank experiment 'independent technology').
```

**S09 (2026-10-24) — discovery loop** (sprint folder `sprints/S09_discovery_loop/`; aggregate-only, counts <20 masked)
```
100_native_vs_lr.py              → native short-read HLA calls vs long-read truth (S09-01; VM side, aggregate output).
101_rna_panel_extract.py         → extract a gene-panel TPM table from AoU blood RNA-seq (VM side; PRIVATE per-person output stays on the VM).
102_expression_dosage_tests.py   → long-read-only features -> blood RNA expression (S09; VM side, aggregate output only).
103_expression_within_g_group.py → expression differences between alleles identical over the peptide-binding exons (S09-11; VM side, aggregate output only).
104_s09_report.py                → S09 discovery-loop report (HTML, inline SVG; aggregate numbers only, every count shown is >= 20 or '<20').
104a_background.py               → Plain-language background drop-downs (inline SVG drawings) for the S09 report. Imported by 104_s09_report.py.
105_novel_allele_evidence.py     → raw-read evidence package for candidate novel HLA alleles (S09-12; VM side, aggregate output only).
106_nonclassical_landscape.py    → diversity landscape of HLA-E/F/G vs the classical genes, by ancestry (S09-13; VM side, aggregate output only).
107_phewas_engine.py             → hypothesis-free HLA/KIR/C4 x ICD-10 phenome scan with controls (S09-14; VM side; person-level tables stay on the VM).
108_mechanistic_rna_tests.py     → mechanistic hypotheses tested against blood RNA (S09-15, S09-16; VM side, aggregate output only).
109_phewas_followup.py           → Within-ancestry follow-up of PheWAS leads from 107 (VM side).
111_olink_nkg2a_test.py          → Plasma protein test of the HLA-E 107G -> NKG2A (KLRC1) association (VM side; All of Us Olink Explore HT, v9 batch-normalised).
112_hla_kir_c4_pqtl_scan.py      → Immunogenetic loci -> plasma proteome scan (VM side; All of Us Olink Explore HT x long-read HLA / KIR / C4 calls).
113_b1402_pqtl_dissection.py     → Dissect the B*14:02 / C*08:02 -> LRPAP1 / LRP1 / SGSH plasma-protein signal (VM side; Olink Explore HT x long-read HLA calls).
114_b1402_rna_test.py            → Does HLA-B*14:02 shift blood RNA of LRPAP1 / LRP1 / SGSH (VM side)?
130_call_confidence.py           → S06 v2: error among CALLED alleles by cohort (template distance 0 / >0 / S10 novel), by distance and CDS depth, pruning policies, AUC (VM side; trio + pair truth sets).
131_call_confidence_report.py    → (local only, not in git) figure + compact report for 130 -> reports/sprint_reports/S06_REPORT_v2.html.
118_novel_ladder.py              → Full distribution of minimally novel HLA alleles (no indel, not exact in the newest IPD) with short-read LOD confidence, buckets, saturation (VM side; writes exact + masked JSON).
119_novel_ladder_report.py       → Builds the interactive selector report from 118's JSON (`119_report_template.html`); refuses to put exact counts outside a *.private.html.
```

`00` exists because the source-derived spec has genuinely ambiguous points that only real data can
settle — chiefly **whether the intermediate files survived on the production VM**. It is cheap
(~2 minutes) and prints one JSON block to paste back. Nothing downstream should be trusted until it
has run.

**On-disk layout:** the ~12,000 per-person `<person_id>/immuannot_output/...` directories live under
`~/pipeline_outputs/people/`, one level deeper than the aggregate `.tsv`/`.fa` files listed above
(which stay directly at `~/pipeline_outputs/`) — see SCHEMA.md "On-disk layout" and RUNBOOK.md
"Step 1c" for why (Jupyter file-browser lag at ~12,000 top-level entries) and the one-time move.
Every script's `--outroot` default already points at `people/`; `02_build_cohorts.py` is the one
exception, since it never scans person directories.

## Research notes

- `research/NOVEL_LIT.md` — novel-allele definition, QC, naming, and the saturation-estimation
  methodology (Chao2/ACE, rarefaction), including the balancing-selection caveat that makes
  neutral-model estimators inappropriate as a headline number.
- `research/VIZ_LIT.md` — figure catalog with literature precedent, ranked, plus the figures that
  specifically exploit our phasing / template_distance / 4-field / continuous-ancestry advantages.
- `reference/IMMUANNOT_GTF_SPEC.md` (repo root `reference/`) — the authoritative output-format spec.

## Testing

The production VM is behind a VPC-SC perimeter that blocks SSH entirely (ENVIRONMENT.md quirk #28),
so every round-trip to real data costs a human paste-cycle. **Nothing gets handed over until it runs
clean against synthetic fixtures.**

```bash
python3 scripts/hla_popgen/tests/make_fixtures.py --outroot /tmp/hla_fixtures -n 300
```

The fixtures deliberately encode the failure modes and signals that matter: ~18% of haplotypes span
two contigs (the cis trap), ~6% have genuine copy number >1, ~8% of people have no output,
`cds_distance`/`cds_mut`/`template_warning` are conditionally absent, `template_distance` is
unquoted, and the novel-allele rate is ancestry-skewed (afr 0.22 vs eur 0.05) with ancestry-biased
allele pools. Any script that fails to recover those gradients has a bug.

## Standing rules

Beyond `SCHEMA.md`'s hard rules: outputs under `reports/` are **aggregate-only** — no bare
person_ids, no raw allele calls, no novel-allele nucleotide sequences. Per-person points in a
derived embedding (PCA/UMAP) are acceptable and unlabeled, matching existing precedent. Novel-allele
sequences stay in `~/pipeline_outputs/` on the VM.
