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

## Critic review of script 25 (2026-09-16, fresh-context agent)
**Verdict: the cs/strand core is correct (hand-checked on all four operators, both strands); the
headline numbers are not paper-ready until the fixes below land.**

- **CRITICAL: only one PAF record per call.** The PAF is raw, unfiltered minimap2, so a large
  intronic indel splits the alignment and only one record is kept (the variant disappears), and
  secondary (`tp:A:S`) rows can win on length. Fix: keep primary rows only, chain multiple primary
  records, or mark the call `split_alignment` and exclude it.
- **CRITICAL: no alignment-coverage floor.** A partly covered allele yields a signature that is a
  subset of the truth and merges with full-length ones; `no_difference` is then an artifact class.
  Fix: `--min-qcov` (≈0.98) and cross-tabulate `no_difference` against coverage.
- **MAJOR: the pooling headline is confounded by template choice.** Signatures include the template
  allele, so identical contigs with different templates count as distinct sequences, inflating
  "one old cluster pooled N sequences". Fix: restrict to the modal template or re-project onto one.
- **MAJOR: small-denominator ratios are committed unsuppressed** (`frac_homopolymer_only`,
  `median_qcov` with `n = <20`) — that can disclose one person's genotype summary under `--genes all`.
- **MAJOR: homopolymer fraction is inflated when the allele FASTA is missing** (every 1 bp indel is
  called homopolymer); the figure carries no flag. Gate the figure on context availability.
- **MAJOR: resume can mix reference contexts** — the parameter hash ignores FASTA contents and the
  script's own mtime.
- **MAJOR: threads won't help.** The per-line PAF parsing is GIL-bound; use a byte-prefix test and
  processes, not threads, on the 4-vCPU VM.
- MINOR: people missing from cohort_membership are silently treated as "related"; cluster class
  taken from an arbitrary row; one figure plots unsuppressed counts; homozygotes contribute 2 to
  the spectrum.
- Six specific missing tests named (split/secondary rows, partial coverage collapse, template
  splitting, suppression on small genes, resume invalidation, corrupt cds.fa.gz).

### Spill-over finding, outside this sprint's scope
**`21_hla_manhattan.py:walk_cs_canonical` anchors a minus-strand deletion at `qend-1-k`, where
script 25's hand-verified logic gives `qend-k` — an off-by-one.** 21 produced the committed
per-site diversity (π) results and the CDS-vs-non-CDS enrichment figures, the "capstone" result.
Roughly half of all rows are minus-strand, so deletion positions in those may be shifted by one
base. Needs checking before that figure is published: re-run a gene with the corrected offset and
compare. Logged here so it is not lost.

## Results — script 25, full cohort, classical genes (2026-09-16)
32,207 depth-4 (non-coding-novel) haplotypes across 11,507 people; 31,886 parsed
(300 excluded by the 0.98 alignment-coverage floor, 21 with no usable alignment row); 9,805
secondary alignment rows filtered out; 0 split/chained alignments.

**The pooling is now measured, and it is severe.**

| | |
|---|---|
| former CDS-hash clusters (classical genes) | 708 |
| distinct non-coding sequences (signatures) behind them | **15,112** |
| median signatures per former cluster (modal template / raw) | 1.0 / 2.0 |
| max signatures in one former cluster (modal template / raw) | **287 / 788** |
| signatures carried by ≥20 unrelated people | 184 |

So "a novel allele carried by 1,531 people" was never one allele: the largest former cluster
contains 287 genuinely different genomic sequences. Any recurrence statistic computed on the old
clusters is an upper bound on nothing meaningful, and the 2026-09-15 draft callout was wrong.

**What the non-coding differences are** (haplotypes): non-homopolymer indel 10,250; SNV only 8,731;
**homopolymer-only 6,778 (21.3%)**; SNV+indel 5,048; no difference within the aligned region 1,079.
So about a fifth of non-coding novelty is the classic long-read error mode and should be discounted.

**Where they fall:** intronic-only dominates; 8.8% of depth-4 haplotypes "touch CDS", and in 2,783
of those 2,808 cases the genomic template differs from the called allele at field 3 — i.e. it is the
template-vs-call mismatch the design anticipated, not a contradiction.

## Where the novel proteins actually are (script 24, full cohort)
Clean novel proteins recurrent in ≥2 unrelated people, by gene class:

| gene class | recurrent clean | total clusters |
|---|---|---|
| mic_tap (MICA/MICB/TAP1/TAP2) | **183** | 716 |
| class_II_accessory (DM/DO/DQA2/DQB2) | 67 | 287 |
| class_II_paralog (DRB3/4/5) | 49 | 180 |
| nonclassical_I (E/F/G) | 32 | 260 |
| **classical_II (DR/DQ/DP)** | **6** | 188 |
| **classical_I (A/B/C)** | **2** | 161 |
| other | 8 | 17 |
| pseudogene_I | 0 | 14 |

Top genes: TAP2 79, TAP1 68, DQB2 26, DMA 22, DQA2 20, MICA 20, DMB 17, DOA 16, MICB 16, G 15.
In the eight classical genes there are only **8** recurrent clean novel proteins in 12,233 people
(B, C, DPA1, DQB1 one each; DPB1 and DRB1 two each).

**Interpretation, and a possible reframing of the paper.** At the protein level the *classical*
HLA genes are close to exhausted in a cohort this size — which is itself a publishable, quantitative
statement, and it is what the coverage analysis (WS3) should confirm. The genuinely undiscovered
space sits in the genes nobody types at scale: TAP1/TAP2, MICA/MICB, the DM/DO accessory genes and
the DRB paralogs. Long reads are the only way to see them, and this cohort is the first at this
size. That is a stronger and more defensible headline than "3,000 novel HLA alleles".

Caveat to check before leaning on it: MIC/TAP novelty may partly reflect a less complete IPD
reference for those genes (note MICB, where 78% of calls are field-3 novel and 22,657 of the 25,282
"catalogued CDS under a new name" cases are MIC/TAP). Distinguish "new biology" from "thin
reference" per gene before publishing.

## Resolution of the 21_hla_manhattan spill-over (2026-09-17)
**CONFIRMED as a real off-by-one, but harmless to the committed results.** The `-` cs operator
consumes no canonical base, so its anchor is a boundary, not a base position: on the minus strand
that is `qend - q_consumed`, while the script used the consumed-base formula `qend - 1 - q_consumed`.
However the affected value only ever lands in `insertions` → `insertion_counts` in the per-gene
JSON, which nothing reads: per-site π, the Manhattan figures and `panel_summary.tsv` are all built
from `observed`, whose substitution and deletion branches are correct. **The capstone
CDS-vs-non-CDS diversity results do not need regenerating.** Patched anyway, with a regression test,
before anything starts consuming `insertion_counts`.


## Correction (2026-09-17, from a fresh-context review of the findings doc)
Three errors in the first write-up of these results, now fixed everywhere:
1. **"8 recurrent clean novel proteins in the classical genes" was wrong as stated.** 8 is the count
   when synonymous-CDS clusters are included; the protein-only count is **6**. Totals: 231 recurrent
   clean novel proteins, 347 including synonymous clusters.
2. **The gene-family breakdown was labelled with the repo's `gene_class` buckets but described with
   different gene names.** Correct per-gene-group counts (protein-only / including synonymous):
   TAP1+TAP2 **104 / 147**; DM+DO+DQA2+DQB2 **68 / 110**; E+F+G **23 / 32**; MIC **22 / 36**;
   DRB3/4/5 **1 / 3**. (`class_II_accessory` in the repo is DMA/DMB/DOA/DOB/DRA, and
   `class_II_paralog` is DPB2/DQA2/DQB2/DRB3/4/5 — not what those labels suggest.)
3. The qualitative conclusion is unchanged and if anything strengthened: recurrent novel proteins
   are overwhelmingly outside the classical genes, and TAP1/TAP2 alone carry 104 of 231.
