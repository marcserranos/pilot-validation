# Canonical data contract — `hla_popgen`

> **Role:** the single source of truth for every intermediate table this sub-project produces and
> consumes. Every script here reads and writes these exact column names. Change this file first,
> then the scripts — never the other way round.
> **Edit:** append a column; never silently rename or repurpose one.
> **Read:** before writing or modifying any script in this directory.

All claims about Immuannot's output format below are sourced from a full read of the upstream
source (`YingZhou001/Immuannot`, `scripts.pub.v3`), archived as `reference/IMMUANNOT_GTF_SPEC.md`.
Where that spec marks something AMBIGUOUS, `00_recon_vm.py` resolves it empirically before any
downstream script is allowed to trust it.

---

## Why these tables exist (the three gaps in the existing pipeline)

The production run's `immuannot_calls.tsv` captures exactly one string per (person, hap, gene) —
the `consensus` typing call. Everything below is already on disk and was simply never extracted:

1. **Contig identity (GTF column 1) is discarded by every existing parser.** A single
   `hap1.gtf.gz` can legitimately hold genes from more than one physical contig when the assembly
   is fragmented across chr6:29.5–33.5Mb (`run_immuannot_person.py` tracks `n_contigs` explicitly).
   Two genes being in the same hap file is **necessary but not sufficient** for them to be in cis.
   The only valid cis key is `(person_id, hap, contig)`. Without it, DQA1~DQB1 heterodimer
   reconstruction is unsound.
2. **57 of 65 genes were never analysed.** All downstream scripts hardcode the 8 classical genes.
   Non-classical, pseudogene, MIC/TAP and C4 rows are present in the raw GTFs and dropped.
3. **The novelty signal is richer than a boolean.** `"new"` is spliced into the allele string at a
   *variable* field depth that encodes severity, and `cds_mut` carries the exact codon-level diff.
   Neither is captured anywhere today.

---

## Gene classification (used as `gene_class` throughout)

| `gene_class` | Genes |
|---|---|
| `classical_I` | HLA-A, HLA-B, HLA-C |
| `classical_II` | HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DRB1 |
| `class_II_accessory` | HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB, HLA-DRA |
| `class_II_paralog` | HLA-DRB2, HLA-DRB3, HLA-DRB4, HLA-DRB5, HLA-DRB6, HLA-DRB7, HLA-DRB8, HLA-DRB9, HLA-DPA2, HLA-DPB2, HLA-DQA2, HLA-DQB2 |
| `nonclassical_I` | HLA-E, HLA-F, HLA-G |
| `pseudogene_I` | HLA-H, HLA-J, HLA-K, HLA-L, HLA-N, HLA-P, HLA-S, HLA-T, HLA-U, HLA-V, HLA-W, HLA-Y |
| `mic_tap` | MICA, MICB, TAP1, TAP2 |
| `complement` | C4A, C4B |
| `other` | HLA-HFE |
| `kir` | the 17 KIR genes — **expected to be entirely absent** (chr19, outside the trim window). Their presence in real output would indicate a trim bug and must be reported loudly, not silently ignored. |

---

## Table 1 — `hla_calls_rich.tsv`

**Grain: one row per `(person_id, hap, contig, gene, copy_index)`.** This is the finest grain
Immuannot emits and the grain every dedup/merge in this sub-project must key on. (The production
run lost an entire cohort's classical-gene calls to a dedup keyed at the wrong grain —
ENVIRONMENT.md quirk #29. Do not repeat it.)

| Column | Type | Source | Notes |
|---|---|---|---|
| `person_id` | str | directory name | |
| `hap` | str | `hap1`/`hap2` | Which physical haplotype assembly. |
| `contig` | str | GTF col 1 | **The cis key.** Never drop this. |
| `gene` | str | `gene_name` attr | For C4, strip the trailing size letter (`C4AL`→`C4A`); see `c4_size`. |
| `copy_index` | int | `.N` suffix on `gene_id` | 1 when no suffix. >1 means a real second mapping cluster (segmental duplication / DRB copy number), not an artifact. |
| `gene_class` | str | lookup above | |
| `consensus` | str | `consensus` attr, transcript row | **The typing call.** May be `undetermined`. Never use `template_allele` as the call. |
| `n_fields` | int | derived | Colon-delimited field count of `consensus`. |
| `is_novel` | bool | derived | `"new"` appears as a field in `consensus`. |
| `novelty_depth` | int/NA | derived | 1-based index of the field equal to `new` (**1**, 2, 3, or 4). |
| `novelty_class` | str/NA | derived | depth 1 → `undetermined` (even the first field is unresolved — arises when `consensusCall()`'s commonprefix truncation collapses tied candidates disagreeing at field 1; **rare but real**, 3 occurrences in a 50-person production sample); depth 2 → `protein_altering` (non-synonymous or frameshift); depth 3 → `synonymous`; depth 4 → `beyond_cds` (intronic/UTR only). Derived from `nameNewHlaAllele()`'s truncation rule. |
| `template_allele` | str | gene row | Structural template only — **not** the genotype. |
| `template_distance` | int | gene row | **Unquoted in the GTF, unlike every other attribute** — the parser regex must not require quotes. Weighted variant count over the *full gene span incl. introns/UTR*: substitution=1, indel run ≤5bp=1, indel run >5bp=2. Not a Levenshtein distance, not normalized. |
| `gene_start`, `gene_end` | int | GTF cols 4/5 of gene row | 1-based inclusive, **contig-relative, not hg38**. This span is the alignment span for `template_distance`. |
| `template_distance_per_kb` | float | derived | `1000 * template_distance / (gene_end - gene_start + 1)`. The comparable-across-genes version — raw counts are not comparable between a 3kb and a 15kb gene. |
| `cds_distance` | int/NA | transcript row | **Conditional** — absent when the gene-level template was already a perfect match (Immuannot never ran a CDS search). NA here means "not applicable", *not* "missing data". Raw minimap2 NM over CDS only. |
| `cds_mut` | str/NA | transcript row | Conditional on `cds_distance > 0`. Pipe-delimited: `ref_allele \| cs_string \| AAref(codon)<AAobs(codon):...`. |
| `n_aa_changes` | int/NA | derived from `cds_mut` | Count of `<`-separated codon diffs in the third pipe field. |
| `template_warning` | str/NA | transcript row | Conditional. Tokens: `partial_CDS`, `inframe_stop`, `no-start_codon`, `no-stop_codon`. |
| `has_warning` | bool | derived | **Do not use as a QC gate or confidence filter.** See the warning policy below. |

### `template_warning` policy — measured, not assumed

`00_recon_vm.py` on the real production cohort (2026-09-03, 50-person sample) measured
**`template_warning` present on 95.3% of transcript rows.** It is near-ubiquitous.

This matters because it invalidates a filter this project already uses elsewhere. The existing
confidence convention (`context/DECISIONS.md`, "Confidence-matched truth comparison";
`scripts/analyze_confidence_matched_truth.py`) is *"`template_distance == 0` AND no
`template_warning`"* — on the real production data that second clause would **reject roughly 95%
of all calls**, not a small tail. Any blanket `has_warning == False` gate silently empties the
analysis.

The reason is semantic: `template_warning` describes whether the **template's CDS could be cleanly
reconstructed from the gene-level alignment** (`searchTemplate.py`'s `checkCDScompleteness()`), not
whether the typing call is wrong. Pseudogenes (`HLA-H/J/K/L/...`) legitimately have no valid start
or stop codon and will always warn, and `partial_CDS` fires whenever the trimmed contig truncates a
gene's span.

**Rule: every filter on warnings must be TOKEN-AWARE and configurable.** Treat `partial_CDS` as
benign by default; `inframe_stop` is the token that genuinely suggests a broken reconstruction and
is the defensible default disqualifier for novel-allele QC. Never gate on mere presence.
| `alleles` | str | transcript row | Comma-joined tied-best candidate reference alleles. |
| `n_tied` | int | derived | Ambiguity of the call. |
| `strand` | str | GTF col 7 | |
| `c4_size` | str/NA | derived | `L`/`S` for C4 rows only (intron-9 length class). |

**Parser requirements (non-negotiable):**
- Read the `#`-prefixed header lines. They already carry `## gene (copy num = 0/1/>1): ...` and
  `## contigs for <gene>: ...`. Use them as a free cross-check that the body parse agrees; a
  mismatch is a parser bug and must fail loudly.
- `template_distance` is unquoted. `cds_mut` contains `|`, `<`, `:` and `*` — do not split on
  those naively.
- A gene absent from the GTF is absent because it failed detection filters
  (`adjNM<=3 or adjNM/match<0.03`) or had no homology — record it as absent, never as NA-imputed.

---

## Table 2 — `hla_cis_pairs.tsv`

**Grain: one row per `(person_id, hap, contig, pair)`.** Built *only* from rows sharing a contig.

| Column | Notes |
|---|---|
| `person_id`, `hap`, `contig` | |
| `pair` | `DQA1~DQB1`, `DPA1~DPB1`, `DRA~DRB1`, plus any pair requested |
| `allele_a`, `allele_b` | the two `consensus` calls, in cis |
| `haplotype_label` | `"DQA1*05:01~DQB1*02:01"` — the analysis unit Cole wants |
| `both_exact` | both had `template_distance == 0` |
| `either_novel` | |
| `cis_confidence` | `physical` when both genes are on the same contig; the row is not emitted otherwise. A companion count of *unpairable* cases (genes split across contigs) must be reported — that number is itself a result about assembly fragmentation. |

**DPA1/DPB1 are adjacent (~10kb) so physical phasing should be near-certain; DQA1/DQB1 span
further (~15–20kb) so expect more fragmentation.** Report the pairable fraction per pair per
ancestry — it is a real finding, not just QC.

---

## Table 3 — `novel_alleles.tsv`

**Grain: one row per distinct novel sequence cluster (identical observed CDS collapsed).**

| Column | Notes |
|---|---|
| `novel_id` | Provisional stable id: `<gene>_nov_<sha1(cds_seq)[:8]>`. Deterministic — recomputable, never a counter that shifts between runs. |
| `gene`, `gene_class` | |
| `nearest_allele` | closest documented IPD allele |
| `cds_distance`, `n_aa_changes`, `novelty_class` | |
| `cds_seq_sha1`, `cds_len` | sequence identity without storing the sequence in shareable outputs |
| `n_haplotypes`, `n_persons` | observation count — **recurrence across unrelated people is the primary evidence it is real, not an assembly artifact** |
| `n_ancestries`, `ancestry_counts` | JSON dict of ancestry → haplotype count |
| `is_homopolymer_indel_only` | QC flag — pure homopolymer indels are the standard long-read artifact class and are excluded from headline counts |
| `passes_qc` | composite: seen in ≥2 unrelated persons, no `template_warning`, not homopolymer-only |

Actual nucleotide sequences go to a **separate** `novel_alleles_seqs.fa` that stays on the VM and
is never written into any file destined for a report — same compliance discipline as the rest of
the repo.

---

## Table 4 — `cohort_membership.tsv`

**Grain: one row per person.** The three-cohort definition lives here and nowhere else.

| Column | Notes |
|---|---|
| `person_id` | |
| `in_sr` | has AoU-native short-read HLA calls (~500K) |
| `in_lr` | has Immuannot long-read calls (~12K) |
| `max_template_distance`, `mean_template_distance` | across their classical genes |
| `n_genes_called`, `n_genes_exact` | |
| `td_stratum` | `exact` (all classical genes td==0) / `near` (max td<=1) / `loose` (max td<=5) / `distant` | 
| `ancestry_pred` | **lowercase in the real file** (quirk #30) — normalize with `.str.upper()` immediately on load, always |
| `p_afr`,`p_amr`,`p_eas`,`p_eur`,`p_mid`,`p_sas` | continuous admixture proportions, parsed from the `probabilities` array (ordered AFR/AMR/EAS/EUR/MID/SAS) |
| `platform`, `trim_tier` | from `immuannot_cohort_full.tsv` |

**Cohort 3 is a sweep, not a fixed threshold.** Every frequency/structure figure must be
reproducible at each `td_stratum` so that the *movement* of a statistic as the filter tightens is
itself the result. Do not hardcode one threshold anywhere.

---

## Hard rules for every script in this directory

1. **Fail fast and loud.** Check the gcsfuse mount as the very first action if the script needs it
   (quirk #26). Never continue past a missing input.
2. **Never `du`/`df` a broad path.** Scope disk checks to specific small directories (quirk #25).
3. **Write incrementally with checkpoints** for anything iterating over ~12,000 people (quirk #22),
   and make reruns resume.
4. **Separate test and real output paths** — never share a path between `--sample` and a full run
   (quirk #22b).
5. **Aggregate-only outputs.** No bare person_ids, no raw allele calls in anything under
   `reports/`. Per-person rows stay in `~/pipeline_outputs/`.
6. **Every script must run against the synthetic fixtures in `tests/`** before it is handed to
   Marc. The VM is expensive and round-trips are slow; a bug found locally costs minutes, a bug
   found on the VM costs a day.
