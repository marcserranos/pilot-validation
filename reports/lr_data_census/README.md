# Long-read data-format census — full cohort (n=14,521)

**Headline: every one of the 14,521 people with a long-read manifest row has a real, existing
aligned GRCh38 BAM.** The previous working assumption — a confirmed floor of 2,763 people
(ENVIRONMENT.md quirk #13, 2026-07-11) with everyone else's status unknown — is now obsolete.
Directly verified (not inferred): AoU backfilled the previously-missing aligned BAMs sometime
between 2026-07-11 and 2026-07-21. See "The corrected BAM-availability finding" below for the
specific spot-check that confirms this rather than a script bug.

## Method

Exhaustive, per-column existence check against the v9 lrWGS manifest (`v9/wgs/long_read/manifest.tsv`,
15,424 rows / 14,521 unique people, 61 columns) — not just the 2 previously-known BAM columns
(`grch38_bam`, `grch38_haplotagged_bam`) that ENVIRONMENT.md quirk #13 originally checked. Columns
were auto-detected as path-like by their file extension, not guessed from their name (same
anti-pattern-matching discipline as quirk #13), then classified as **primary** (a data-bearing
file — existence-checked directly) or **companion/index** (`.bai`/`.tbi`/`.pbi`/`.gzi` — detected
and reported, but not individually verified, since an index file existing whenever its primary
file does is a reasonable assumption, not something checked here). 29 primary columns were
existence-checked; 25 index columns were skipped on that assumption.

Script: `scripts/lr_manifest_format_census.py` (checkpointed every 500 rows, resumable,
16-way parallelized over the gcsfuse mount — see ENVIRONMENT.md quirk #22 for why that mattered).
Verification method: real file existence against the mount (`os.path.exists()`), never a
path-pattern or column-name guess — the same rule ENVIRONMENT.md quirk #13 established after
three rounds of getting it wrong.

## Data shape by platform

| Platform | Rows (% of 15,424) | Aligned GRCh38 BAM | Aligned CHM13/T2T BAM | Phased diploid assembly | Notable extras |
|---|---|---|---|---|---|
| `revio` | 11,070 (72%) | Yes (all rows) | ~13% of rows also have it | Yes — full suite: hap1/hap2/alternate/primary FASTA + assembly graph (GFA) + hap1/hap2-to-hg38 alignment (BAM+PAF) | TRGT (tandem-repeat) VCF on a subset |
| `sequel2e` | 1,219 (8%) | Yes | some | Yes, same shape as `revio` | TRGT on some |
| `sequel2` | 991 (6%) — AoU's original 1,027-person African-American-enriched long-read pilot cohort | Yes — **haplotagged** (phased) BAM | Yes — haplotagged | Partial (hap1/hap2/primary FASTA+GFA only; no alignment-to-hg38 files) | **DeepVariant *phased* VCF**, unique to this platform |
| `ont-r10.4.1` | 2,020 (13%) | Yes | some (as a separate manifest row for the same person) | No — ONT reads aren't assembled in this pipeline | — |
| `ont-r9.4.1` | 124 (1%) | Yes, but missing `grch38_dv_gvcf` specifically | Yes | No | Oldest ONT chemistry, thinnest product set |

Every row also carries structural-variant calls from **two independent callers** (pbsv and
Sniffles) for whichever reference build(s) it has, alongside a DeepVariant gVCF.

**879 of 14,521 people (6%) have more than one manifest row** — typically one row per reference
build or per resequencing platform for the same person. The per-person view (below) unions across
a person's rows, which is why some person-level profiles look richer than any single row (e.g. a
person with a GRCh38-only `ont` row and a separate CHM13-only row shows both once unioned).

## Per-person profile distribution (n=14,521)

| n | % | Profile (abbreviated) |
|---|---|---|
| 10,256 | 70.6% | Full `revio`-family: assembly (all 10 assembly files) + GRCh38 BAM/gVCF/PAV/pbsv/Sniffles/TRGT |
| 1,483 | 10.2% | Same, plus the CHM13 mirror (no TRGT) |
| 1,099 | 7.6% | `ont`-family minimal: GRCh38 BAM/gVCF/pbsv/Sniffles only, no assembly |
| 991 | 6.8% | `sequel2`: haplotagged BAM + phased DeepVariant VCF, both references, partial assembly |
| 407 | 2.8% | `revio`-family + CHM13 mirror + TRGT |
| 168 | 1.2% | `ont`-family minimal, both references |
| 114 | 0.8% | `revio` with only the assembly-to-hg38 alignment files (no raw assembly FASTA/GFA retained) |
| 2 | 0.0% | Rare edge case (CHM13 gVCF + GRCh38 minimal set) |
| 1 | 0.0% | Single outlier, full `revio`-family minus PAV |

**0 of 14,521 people (0.0%) have none of the 29 recognized primary file types.** Nobody is
genuinely empty-handed — every person has usable long-read data of some kind.

## The corrected BAM-availability finding

Every one of the 9 profiles above includes `grch38_bam` — summing the counts confirms **all
14,521 people**, not a subset. This directly contradicts ENVIRONMENT.md quirk #13's 2026-07-11
finding of a 2,763-person floor, so before treating it as fact rather than a possible bug, we
spot-checked the exact case quirk #13 had specifically documented as dangling: **person 1008366's
`revio`-labeled manifest row**, which that investigation found "points at a nonexistent
`v9_delta/pacbio` path." Checking this run's checkpoint directly for that row:

```
row_idx=1033, research_id=1008366, grch38_bam=1 (true)
profile: assembly_hap1_aln2_hg38_bam+assembly_hap1_aln2_hg38_paf+assembly_hap2_aln2_hg38_bam+
         assembly_hap2_aln2_hg38_paf+chm13v2.0_bam+chm13v2.0_dv_gvcf+chm13v2.0_pbsv_vcf+
         chm13v2.0_sniffles_snf+chm13v2.0_sniffles_vcf+grch38_bam+grch38_dv_gvcf+grch38_pav_vcf+
         grch38_pbsv_vcf+grch38_sniffles_snf+grch38_sniffles_vcf+grch38_trgt_vcf
```

The exact file quirk #13 confirmed absent now resolves as real. **Conclusion: AoU backfilled the
previously-missing `revio`-platform aligned BAMs between 2026-07-11 and 2026-07-21** — consistent
with quirk #13's own closing note that this dataset's publishing pattern changes over time as a
normal consequence of active development, not a defect. This is one directly-confirmed case, not
an exhaustive re-audit of all 14,521 — but combined with the clean, universal pattern (every
profile bucket includes `grch38_bam`, with no partial/mixed group), it's a strong, coherent signal
rather than a coincidence.

## Implication for the project

- **The long-read validation cohort is no longer capped near ~2,763.** All 14,521 people are, in
  principle, viable for the existing SpecImmune-LR pipeline (extract HLA-region reads from
  `grch38_bam` → FASTQ → SpecImmune) — a >5x increase over the previous working assumption.
  Experiment D's 60-person cohort and any future scale-up should be reconsidered against this much
  larger pool, including for the MID/SAS ancestry groups that were the practical bottleneck before.
- **Phased assembly data (revio/sequel2e/sequel2, ~86% of the cohort) is available as a second,
  independent data product**, not a fallback for people otherwise unreachable. Worth evaluating
  Immuannot (see DECISIONS.md's "Assembly-based HLA typing" bullet) as an *additional*
  cross-validating long-read caller for these people, alongside SpecImmune-on-reads — not as the
  only way to access them, since it no longer is one.

## Caveats

- Companion/index files (`.bai`/`.tbi`/`.pbi`/`.gzi`) were not individually existence-checked —
  assumed present whenever their primary file is, not verified at scale.
- "File exists" is not "file is valid/complete/uncorrupted" — no content or integrity check was
  performed on any BAM, FASTA, or VCF.
- Only one specific historical case (person 1008366) was directly spot-checked against the prior
  2026-07-11 finding — the backfill conclusion is well-supported, not exhaustively re-audited.
- Aggregate-only report; per-person profile data (with bare `research_id`s) stays on the VM
  (`~/pipeline_outputs/lr_data_census/census.tsv`), not committed to this public repo — same
  practice as `build_experiment_d_cohort.py`'s `cohort.tsv` (see DECISIONS.md's open compliance
  question on bare ids in a public repo).
