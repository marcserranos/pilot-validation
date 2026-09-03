# `hla_popgen` — HLA population genetics, admixture, and novel-allele discovery

> **Role:** the second-generation analysis of the full-cohort Immuannot long-read run — everything
> the first pass (`scripts/production_analysis/`) left on the table.
> **Read:** `SCHEMA.md` before touching any script here. It is the binding data contract.

## Why this exists

The production run called HLA on ~12,000 All of Us participants from phased long-read assemblies.
The first round of post-processing (`scripts/production_analysis/`, results in
`reports/full_immuannot_lr_calling/`) answered "did the pipeline work, and does HLA structure track
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

```
00_recon_vm.py        → reports/hla_popgen/recon_report.{md,json}   [RUN THIS FIRST]
01_extract_rich.py    → hla_calls_rich.tsv (Table 1), hla_cis_pairs.tsv (Table 2)
02_build_cohorts.py   → cohort_membership.tsv (Table 4)
03_novel_alleles.py   → novel_alleles.tsv (Table 3), novel_alleles_seqs.fa [VM-only]
04_allele_saturation.py → discovery curves, Chao2/ACE richness, "% of allele space found"
05_figures_frequency.py   → allele frequency by ancestry, diversity indices
06_figures_structure.py   → PCA/UMAP, Fst, continuous-admixture representations
07_figures_crosscohort.py → SR-vs-LR bias, template_distance by ancestry, cis heterodimers,
                            non-classical gene diversity, resolution cascade
```

`00` exists because the source-derived spec has genuinely ambiguous points that only real data can
settle — chiefly **whether the intermediate files survived on the production VM**. It is cheap
(~2 minutes) and prints one JSON block to paste back. Nothing downstream should be trusted until it
has run.

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
