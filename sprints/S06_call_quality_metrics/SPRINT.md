# Sprint S06 — per-call quality metrics, exon-level depth, and how they predict typing error

*Opened 2026-10-20 on `main`. Owner: Marc. Source: Marc's brief in the project thread, 2026-10-20 (replaces the
subjective 0–7 confidence score from 56 with measurable quantities). Log: [`LOG.md`](LOG.md) (append-only).*

## Marc's question
Which measurable per-call quantities (Immuannot's own distances and flags, read depth over the exons / CDS of the called
gene, assembly context) predict a typing error, and what does removing the worst x% of calls per quantity do to the
typing correctness of what is left? Replace the hand-built 0–7 score with these.

Scope is exactly this. Aleix's repertoire work, KIR and everything else are out for now.

## Rules for this sprint (Marc, 2026-10-20)
- Granular, high-resolution, paper-grade figures: fine bins, rolling curves, 0.25% pruning steps, not deciles or sparse matrices.
- The "<20 people" masking is **deferred** to a later pass. Until then nothing from this sprint is committed or pushed
  (`reports/hla_popgen/70_*` stays untracked); the masking pass happens before any commit.
- Person-level tables stay inside the Workbench (INTERNAL). Only aggregate fine-grained TSVs come out, and figures are drawn
  locally from those TSVs (no point-per-person plots).
- No VM job over one day; every long job checkpoints and resumes.

## Truth sets (what "typing error" means here)
1. **yak trioeval / Mendel (57b)** — primary. A child call is an error when the parent that yak says transmitted that
   haplotype does not carry the child's 2-field allele. ~179 trios, so a few thousand calls and tens of errors per gene
   class. Caveat: the mismatch can come from the parent's call too (57d), and an origin-assignment error counts as a typing error.
2. **Relative pairs (56b)** — secondary, more power: first-degree pairs sharing MHC, error = no shared 2-field allele
   (typing error or recombination). Pair metric = the worse of the two calls.
3. **Twin / duplicate pairs (56)** — tertiary, cleanest: discordant genotype between two samples of one person.

**Error label (Marc, 10-20):** only the Mendelian mismatch — the parent yak says transmitted the child's haplotype does not carry
the child's 2-field allele. Metrics attributed to the child's call (primary), the parent's call, or the worse of both (all
three analysed). No yak rerun is needed: 57b's per-call origin table (179 trios, ~1 h on 2026-10-06) is read as is.

## What counts as CDS (ground truth for the exon depth; goes in the report)
- Annotation: GENCODE v46 basic, GRCh38 (`~/ref/gencode.v46.basic.annotation.gtf.gz`). CDS = the `CDS` rows of the gene's
  REFERENCE transcript (MANE Select, else Ensembl canonical, else the protein-coding transcript with the longest CDS); `exons`
  = that transcript's exons; per-exon and ARS definitions use the same transcript. A union over isoforms was tried first and
  rejected: it adds alternative-isoform coding exons (HLA-A 1,211 bp vs 1,098 bp) and would dilute the depth. UTR and introns are not CDS. Pure-UTR exons are not measured as exons.
- Why this and not another: it is the public reference annotation of the same coordinates the reads are aligned to. Immuannot
  itself uses IPD-IMGT/HLA 3.55.0 allele annotation (`alleles.csv.gz`, `exon=`/`CDS=` per allele) projected onto each contig;
  that is what defines the exons of a CALL, but it lives in contig coordinates, not GRCh38.
- Cross-check to run on the VM and report: per gene, GENCODE CDS length and exon count vs the IPD CDS length (`CDSseq/<gene>.fa.gz`
  headers) for the allele the hg38 reference carries; flag any gene where they differ by more than a few bp.
- ARS exons (class I exons 2–3, class II exon 2) follow the standard domain-exon map; they are a convention, stated in the report.

## Metrics per call (grain: person × hap × gene)
- Immuannot: `template_distance`, `template_distance_per_kb`, `cds_distance`, `n_aa_changes`, `template_warning`,
  `n_tied` (alleles tied for best CDS match), novelty depth, copy index.
- Read depth (new, script 70, GRCh38 long reads): whole gene, exon union, CDS union, each exon, ARS exons (class I exons 2–3,
  class II exon 2) — mean depth at MAPQ ≥ 0 and ≥ 20, relative to the person's genome control depth; per-base minimum, 10th
  percentile and share of CDS bases under 5/10/20×; MAPQ20/MAPQ0 ratio (paralog ambiguity).
- Assembly context: contigs per haplotype across the MHC, same-contig with phasing partner, platform.
- Person: MHC window depth (abs and relative).

## Queue (first part of the sprint)
| # | Script | What | Where |
|---|---|---|---|
| E1 | `70_call_exon_depth.py` | exon / CDS / per-exon depth per person × gene; per-base CDS summary; population depth profile along each gene | VM |
| E2 | `71_call_quality_table.py` | join Immuannot metrics + depth + truth into one call table; distributions (depth per gene, per-person MHC depth) | VM |
| E3 | `71` | correlations of every metric with typing error (fine bins, AUC, rank correlation, logistic with person-clustered CIs, metric × metric correlation) | VM |
| E4 | `71` | pruning curves: remove the worst 0–50% of calls per metric (0.25% steps), error among the rest, vs random removal; combined model vs the 0–7 score | VM |
| E5 | `72_call_quality_figures.py` | paper-grade figures from the aggregate TSVs | laptop |
| E6 | `73_depth_schematic.py` + report | methodology schematic (gene / exon / CDS / ARS, how depth is counted) and the report | laptop |
| E7 | `71`/`57b` extension | **Marc 10-20:** the MHC gene-pair matrix from yak trioeval truth. Panel 1 (existing figure, retitled "% of haplotypes with both genes on one *contig*"); panel 2: % of haplotypes with both genes from the same parent over ALL haplotypes (truth, any contig); panel 3: the same restricted to pairs on one contig (how often contig co-location is a true phase). Needs DRA and C4A added to 57b's gene list (rerun its per-gene step) | VM |

Naming: the old figure title "Physical phasing level" is wrong — co-location on one contig is not physical phasing. Use
"contig co-location" for the assembly property and "same-parent (trio truth)" for the physical-phase truth.

Optional extension: haplotype-resolved depth from the haplotagged BAM (HP tag), if the manifest has it.
