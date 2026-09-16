# WS1 — Novelty re-definition by nomenclature field (L2 brief)

*Owner: orchestrator. Implementers: scripts 24 and 25. Status: audit done → implementation.*

## Objective (plain language)
Say exactly what "novel" means at each level of HLA naming, count it without artifacts or double
counting, and find the alleles that actually change the protein. Marc's question: "what did we
discover at field 2 (new protein) and field 3 (new synonymous coding sequence), and how much of the
coding/protein space is new?"

## Background: how Immuannot names a novel allele
Immuannot writes `new` in place of the deepest field it can't match (`reference/IMMUANNOT_GTF_SPEC.md:130-163`):

| call looks like | `novelty_depth` | meaning |
|---|---|---|
| `HLA-A*new` | 2 | protein differs from every known allele (missense, indel or frameshift) |
| `HLA-A*01:new` | 3 | protein known, CDS differs (synonymous change) |
| `HLA-A*01:01:new` | 4 | CDS identical to a known allele; difference is intron/UTR only |
| `HLA-*new` (rare) | 1 | tied candidates disagree at field 1 |

## Audit findings (2026-09-16, code audit of 01/03/04 + `novel_alleles.tsv`)
1. **Clustering key = (gene, sha1 of observed CDS)** (`03_novel_alleles.py:330-336`). No genomic
   sequence is stored anywhere in Tables 1–3. So a depth-4 ("beyond_cds") cluster pools *every*
   haplotype whose CDS equals a given known CDS but whose genome differs from IMGT in any way. The
   largest pools 7,713 people / 10,124 haplotypes. **Marc's hypothesis is confirmed**: these are not
   single alleles, and their `n_persons` must not be reported as recurrence of one sequence.
2. **Artifact rules** (`03_novel_alleles.py:391-392`), cluster-level OR over members:
   - Rule A: a member carries the `partial_CDS` or `inframe_stop` warning.
   - Rule B: every CDS difference is a homopolymer indel (the HiFi error mode).

   Counts: 17,545 of 21,538 clusters flagged (13,735 by rule B, 3,810 by rule A only).
   **17,020 of 19,486 "protein_altering" clusters (87%) are flagged.** The earlier "90.8% of novel
   clusters are protein-altering" headline was computed before the artifact gate, so it mostly
   counts homopolymer frameshift artifacts. It must be retracted.
3. **Why the artifacts fell through:**
   - `04_allele_saturation.py` assigns every novel call its cluster ID with **no tier filter**
     (`build_allele_identity_table`), so ~81% of the "novel alleles" in the saturation curves and
     Chao2 estimates are flagged clusters.
   - The per-ancestry novel-call rate (`03 novel_rate_by_ancestry`) also counts all novel calls.
   - **No relatedness filter** exists in 04 (827 people have a long-read relative), which inflates
     doubletons.
4. `n_aa_changes` counts `<` across **all tied candidates' diff strings concatenated**, not one
   allele's changes. It is unusable as "number of amino-acid changes".
5. A rule-B flag can hit a real frameshift allele, and a rule-A flag describes reconstruction
   completeness, not call correctness. Flagged ≠ proven false. Keep flags as labels, not deletions,
   and report both views.

## What we can recover without re-running Immuannot
- `cds.fa.gz` (observed CDS per gene copy) → exact CDS and translated protein for depth-2/3 calls.
- `~/tools/Immuannot_refdata/CDSseq/<gene>.fa.gz` (IPD-IMGT snapshot shipped in
  `Data-2024Feb02`) → known CDS for every allele, so we can check exact CDS/protein membership.
- `mm2.ipd.gen.paf.gz` (every IPD genomic allele aligned to the contig, `asm5 --cs`, 100% present)
  → the **cs difference string between the contig and the template genomic allele**. Depth-4
  novelty can therefore be split into distinct non-coding sequences and classified (SNV, indel,
  homopolymer indel).

## Implementation spec — script 24 `24_novelty_by_field.py` (VM)
Inputs: Table 1 (`~/pipeline_outputs/hla_calls_rich.tsv`), `cohort_membership.tsv`, per-haplotype
`cds.fa.gz` under `--outroot` (default `~/pipeline_outputs/people`), `--refdata`
(default `~/tools/Immuannot_refdata`; glob `**/CDSseq/*.fa.gz`), relatedness TSV
(default as in `11_relatedness_cohort_overlap.py`), optional `--imgt-latest-nuc` (a newer
`hla_nuc.fasta` for a release-drift check).
Reuse, by importlib and without editing, 03's `parse_cds_fasta`, `match_novel_rows`,
`warning_tokens` and `is_homopolymer_indel_only`.

1. **Unrelated set:** greedy removal from `samples_relatedness.tsv` pairs with kin ≥ 0.0442
   (third degree or closer), restricted to long-read people. Drop the member with more relatives
   first; break ties by sorted id. Output a boolean per person (VM-local only).
2. **Ancestry:** `ancestry_pred` plus `strict_ancestry` (the label only if that ancestry's admixture
   proportion is ≥ 0.9, else NA). Check the actual proportion column names in
   `_viz_common.load_cohort_membership`.
3. **Field-class counts (A):** per haplotype-gene call (copy_index==1, classical + all genes), the
   field class is one of:
   `known`, `f4_noncoding`, `f3_synonymous`, `f2_protein`, `f1_undetermined`, `uncalled`
   (consensus undetermined).
   Artifact label per call (not per cluster):
   - `homopolymer_indel` (rule B on this call's `cds_mut`)
   - `partial_cds`, `inframe_stop` (warning tokens)
   - else `clean`
   
   Write `novelty_field_counts.tsv`: gene, gene_class, ancestry_scheme (pred|strict), ancestry
   (incl. POOLED), unrelated_only (bool), field_class, artifact_label, n_haplotypes.
4. **Sequence-level truth check (B), depth-2/3 calls with a matched observed CDS:**
   - `cds_known` = observed CDS exactly equals some known CDS for that gene in refdata (a naming
     artifact: "novel" but actually catalogued).
   - `protein_known` = translation equals the translation of some known CDS.
   - frameshift = CDS length % 3 != 0.
   - premature stop = `*` before the last codon.

   Then classify as `novel_protein`, `novel_cds_synonymous`, `cds_known`, `protein_known`, or
   `frameshift_or_stop`. Cluster by protein sequence (novel_protein) or by CDS
   (novel_cds_synonymous). For each cluster:
   - gene, n_haplotypes, n_persons, n_persons_unrelated, per-ancestry person counts (both schemes)
   - artifact labels present (any / all members)
   - nearest known protein allele (minimum Hamming distance when lengths are equal, otherwise
     Levenshtein on protein) and its distance
   - the amino-acid differences as `pos:ref>alt` (1-based on the precursor protein)
   - for each difference, `in_groove` = codon lies in exon 2 or 3 (class I) or exon 2 (class II),
     using exon lengths derived from the known CDS set if refdata carries exon annotation;
     otherwise NA, stating why

   Write `protein_level_novel_clusters.tsv` (VM-local full version) and a committed version with
   counts 1–19 written as `<20`.
5. **Release drift (E):** if `--imgt-latest-nuc` is given, count how many novel CDS and protein
   clusters are present in the newer release.
6. **Allele count tables for WS3 (D, VM-local only, `allele_counts_by_resolution.tsv`):** per gene
   × resolution × ancestry_scheme × ancestry × unrelated_only, with rows `allele_id, n_haplotypes`.
   - `protein` resolution: known calls → first 2 fields; clean novel_protein → its protein cluster
     id; cds_known/protein_known → mapped to the matching known 2-field name.
   - `cds` resolution: first 3 fields / CDS cluster.
   - `exclude_flagged` variants: `all` vs `clean_only`.
   - Uncalled or unmatched calls are dropped and their count reported.
7. **Report** `novelty_by_field_report.md`: the tables plus a plain-language explanation. Also write
   `summary.json` with headline numbers.
8. **Figures** (matplotlib, `_viz_common` palette, dpi 200):
   - (i) stacked horizontal bars per gene: fraction of haplotypes by field class, with artifact
     hatching;
   - (ii) per-ancestry rate of clean f2/f3/f4 novelty per 100 haplotypes (strict ancestry), one
     row per field class;
   - (iii) funnel/waterfall: novel calls → matched CDS → not artifact → not in IMGT CDS →
     novel protein → recurrent in ≥2 unrelated.

Tests with fixtures: field-class mapping, translation/frameshift, cds_known/protein_known detection,
unrelated greedy selection, count suppression, protein diff positions.

## Implementation spec — script 25 `25_noncoding_novelty_paf.py` (VM)
Target: depth-4 calls. For each person/hap/gene copy with `novelty_depth==4`:
1. Find the PAF row in `mm2.ipd.gen.paf.gz` whose query is the call's **template allele**
   (Table 1 `template_allele`; if absent, the consensus-matching allele with the lowest NM). The
   target contig must equal the call's contig, and the target span must overlap
   `gene_start..gene_end`. Reuse the strand handling and row-selection approach of
   `21_hla_manhattan.py` (read it; don't import from it unless its functions are clean to import).
2. Parse the `cs:Z` tag into difference events in **query (reference allele) coordinates**:
   - SNV `*ab`
   - insertion `+seq`
   - deletion `-seq`
   - a homopolymer flag: the indel sequence is a single repeated base and extends a run of the same
     base in the flanking reference
   - the length of the flanking homopolymer run
   
   Normalize on the minus strand so that signatures agree across strands (this is the key
   correctness point; add a test).
3. Signature = (gene, template_allele, sorted tuple of events). Cluster calls by signature.
   Per cluster: n_haplotypes, n_persons, n_persons_unrelated, ancestry counts, n_events,
   n_snv, n_indel, n_homopolymer_indel, `all_homopolymer` flag, and event positions.
   Whether a position is exonic is unknowable without the template's exon map, so use the
   template's CDS alignment only if it is trivially available; otherwise leave NA.
4. Also per original CDS-hash cluster (`novel_id` from `novel_alleles.tsv`): how many distinct
   signatures it contains. This quantifies the pooling.
5. Outputs:
   - `noncoding_signature_clusters.tsv` (committed version with suppression)
   - `beyond_cds_pooling.tsv`
   - `noncoding_novelty_report.md`
   - figures: (i) distribution of distinct signatures per former cluster; (ii) fraction of depth-4
     haplotypes explained by homopolymer-only signatures, per gene and per ancestry; (iii) SFS of
     non-homopolymer signatures.
6. Threaded I/O like 21 (`--threads`). Provide `--limit` for a pilot.

## Open questions this brief will answer
- What fraction of depth-4 novelty is homopolymer-indel-only (a likely HiFi artifact) vs SNV-bearing?
- How many truly novel proteins recur in ≥2 unrelated people, per gene and per ancestry?
- Is the ancestry gradient present in clean f2 (protein) novelty, or only in f4?

## Distilled (for SPRINT.md)
- Novel clusters were keyed on CDS only, so non-coding "alleles" pool many sequences (up to 7,713
  people). Confirmed.
- 87% of "protein-altering" novel clusters are flagged artifacts (mostly homopolymer indels). The
  "90.8% protein-altering" claim is retracted.
- Saturation (04) counted artifacts as novel alleles and did not remove relatives. It must be
  re-run (WS3).
