# TRUST4 — how it actually works, what we can tune, where it breaks

**Written 2026-08-17.** Sourced from TRUST4's **actual source code** (`BamExtractor.cpp`,
`main.cpp`, `Annotator.cpp`, `run-trust4`, v1.1.10-r639) and README, not from summaries —
several things below are not documented anywhere in the README and were read directly off
the implementation. Cross-checked against our own live runs.

Goal: enough understanding to **design well**, not to know every nuance.

---

## 0. The problem TRUST4 solves, in one picture

Every ordinary RNA-seq tool works by matching reads against a reference genome. Immune
receptor genes break that assumption at the root.

A T cell or B cell builds its receptor gene by **physically cutting and re-joining its own
DNA** — picking one V segment, (one D segment), one J segment from a pool, deleting random
numbers of bases at each junction, and inserting random bases via an enzyme called **TdT**.
This is **V(D)J recombination**. The resulting gene **exists in exactly one cell lineage in
one person and in no reference genome anywhere**.

So an aligner facing a receptor read has two options, both bad: force it onto the germline
V or J segment it partially resembles (losing the junction — which is the only part that
carries the cell's identity), or fail to map it at all. Either way, the biologically
interesting sequence is destroyed before any expression tool ever sees it.

**TRUST4's answer: don't align. Assemble.** Find the reads that *look like* receptor reads,
throw away the coordinates, and rebuild the original transcript from overlapping reads —
the way you'd reconstruct a shredded document from fragments rather than looking it up in a
library. Only afterwards does it compare the rebuilt sequence to a reference, and then only
to *name the parts*, not to find them.

The output unit is the **CDR3** (Complementarity-Determining Region 3) — the hypervariable
junction created by the recombination, and the part of the receptor that physically contacts
the antigen. CDR3 ≈ the receptor's fingerprint.

---

## 1. The four stages

`run-trust4` is a thin Perl wrapper that runs four programs in sequence. Knowing the split
matters because **each stage has different costs and different failure modes**, and you can
restart at any of them with `--stage`.

```
  BAM ──[0]──> candidate FASTQ ──[1]──> contigs ──[2]──> annotated ──[3]──> report.tsv
       bam-extractor          trust4          annotator       perl scripts
       ~95% of runtime        ~seconds        ~seconds        instant
```

| Stage | Binary | What it does | Cost |
|---|---|---|---|
| **0** | `bam-extractor` | Find reads that might be receptor-derived | **~8.5 min of our ~9 min** |
| **1** | `trust4` | De novo assembly of those reads into contigs | seconds |
| **2** | `annotator` | Align contigs to IMGT; call V/D/J/C and CDR1/2/3 | seconds |
| **3** | `trust-simplerep.pl`, `trust-airr.pl` | Format into report tables | instant |

**The intuition to keep:** TRUST4 is not really an assembler-bound tool. It is an
**I/O-bound filter with a cheap assembler bolted on the end.** Almost all the wall-clock is
stage 0 reading the BAM. This single fact explains every performance observation we've made.

---

## 2. Stage 0 — candidate extraction (the part that actually costs us)

This is the stage we've spent all our compute on, so it's worth understanding exactly.
Read directly from `BamExtractor.cpp`.

### 2.1 Three capture routes

TRUST4 sweeps the BAM and keeps a read if **any** of these fire:

1. **Coordinate route.** The read's alignment overlaps a V/D/J/C gene interval. The
   intervals come from the `-f` file (`hg38_bcrtcr.fa`), which stores, per gene:
   `name chrom start end strand`. Implemented as a linear sweep with a moving pointer —
   **which is why the BAM must be coordinate-sorted** (AoU's are).
2. **Unmapped route.** The read pair failed to align at all. These are k-mer screened
   against the receptor reference (`refSet.HasHitInSet`) before being kept. **This route is
   where the most heavily recombined receptors live** — the ones so far from germline that
   STAR gave up. It is also exactly the route that broke on AoU data until
   `--abnormalUnmapFlag`.
3. **Alt-contig route.** The read is on a "weird" contig — the check is literally
   *"does the chromosome name contain `_` or `.`"* (`ValidAlternativeChrom`), which catches
   `chr6_..._alt`, decoys, HLA contigs. Also k-mer screened.

> **AoU-specific consequence: route 3 is dead.** AoU's RNA-seq reference excludes ALT, HLA
> and decoy contigs entirely, so no read can ever be on one. We lose that rescue path by
> construction. Small systematic sensitivity cost; not fixable on our side.

### 2.2 The k-mer screen, and its one tunable-by-accident threshold

For routes 2 and 3, a read must share a long enough exact stretch with some receptor
reference sequence. From the source:

```
hitLenRequired = 21                      // paired-end default
hitLenRequired = 17                      // single-end (more aggressive)
if (readLen / 5 > hitLenRequired) hitLenRequired = readLen / 5
if (hitLenRequired > 101)         hitLenRequired = 101
```

For AoU's **151 bp** reads: `151/5 = 30`, so **30 bp** is the required exact match. Not
exposed as a flag — it is derived from read length. Worth knowing because it means
**recovery sensitivity is a function of read length**, and AoU's read length is fixed. We
cannot tune this without patching source, and we should not.

A **low-complexity filter** also drops reads that are ≥50% one base, or have ≥10% Ns, or
have ≤2 occurrences of two or more bases. Sensible junk removal.

### 2.3 The two-pass problem — the real reason our runs are slow

**Not in the README. This is the single most useful thing in this document.**

For **paired-end** data, `bam-extractor` cannot emit a read the moment it finds it, because
it needs *both mates* and the mate may be megabases away in a coordinate-sorted file. So:

- **Pass 1** sweeps the entire BAM and records candidate read **names** into a map.
- Then `alignments.Rewind()`.
- **Pass 2** sweeps the **entire BAM again**, pulling the actual sequences for names in that
  map, emitting a pair once both mates are in hand.

Single-end data returns after pass 1 (there's an explicit early-exit in the source).

**So: paired-end TRUST4 reads every byte of the BAM twice.** AoU BAMs are ~100M+ read
pairs. This is the entire explanation for:
- why extraction is ~8.5 min while assembly is ~25 s;
- why 24 threads on 8 cores produced no contention (threads only help the k-mer screening,
  not the sequential BAM read);
- why `--jobs` parallelism plateaued — we were saturating **network bandwidth to the
  gcsfuse-mounted bucket**, not CPU;
- why choosing a higher-vCPU machine (more default GCP egress bandwidth per vCPU) was the
  right lever and a GPU would be worthless.

Everything we concluded empirically about the bottleneck is confirmed in the source. Good —
that means the model is right and we can now reason about it instead of measuring it.

---

## 3. Stages 1–3 — assembly, annotation, reporting

### 3.1 Assembly (`trust4`)

Overlap-based de novo assembly of the extracted reads into **consensus contigs**, seeded by
a k-mer index (`-k`, default 9 — this is the *contig indexing* k-mer, distinct from the
extraction screen above). Paired-end mate information is used to **extend** contigs beyond
what one read can span, which is how TRUST4 recovers near-full-length receptors rather than
just the CDR3 window.

Key intuition: **it assembles per-sample, not per-cell.** In bulk data all reads from all
lymphocytes are pooled. A contig is "some clone's receptor," with no cell of origin.

### 3.2 Annotation (`annotator`)

Realigns each contig against the IMGT reference (`human_IMGT+C.fa`) and labels:
V, D, J, C genes (up to 3 candidates each, ranked by similarity) and CDR1, CDR2, CDR3.

**IMGT** = the international ImMunoGeneTics database, the canonical catalogue of germline
receptor gene segments and the source of names like `TRBV7-2*01` (gene `TRBV7-2`,
allele `*01`).

**Two behaviours here that will bite us if unknown:**

- **CDR3 imputation is ON by default, for TCRs only.** If the assembly is partial, the
  annotator will *infer* missing CDR3 nucleotides. The `--noImpute` flag turns it off, and
  it is **not exposed through `run-trust4`** — you'd call `annotator` directly. Imputed
  CDR3s are flagged in the output: `CDR3_score = 0.01` (i.e. raw score 1.00). **Any
  downstream sequence-level modelling should probably filter these out**, because an imputed
  sequence is a guess, and feeding guesses into an embedding model manufactures signal that
  isn't there. → design decision, §7.
- **Scores encode confidence.** In `trust_cdr3.out`, `CDR3_score` of `0.00` = *partial*
  CDR3, `0.01` = *imputed*, other values = motif signal strength. This is a free,
  per-sequence quality filter we are currently not using at all.

### 3.3 Reporting

`trust-simplerep.pl` collapses to the final table; `trust-airr.pl` emits **AIRR** format
(Adaptive Immune Receptor Repertoire — the community standard interchange format, which
every downstream repertoire tool reads). **Using the AIRR file rather than `report.tsv` is
the right move for interoperability** and costs nothing, it's already written.

---

## 4. Inputs

| Input | Flag | What it is |
|---|---|---|
| Aligned reads | `-b` | The BAM. Coordinate-sorted. |
| Coordinate reference | `-f hg38_bcrtcr.fa` | V/D/J/C gene **genomic coordinates + sequence**. Required for BAM input — coordinates are what route 1 uses. |
| Annotation reference | `--ref human_IMGT+C.fa` | IMGT allele database. Technically optional, **effectively mandatory** — without it you get no allele names. |

Both reference files **ship pre-built in TRUST4's repo** — we discovered live that the
`BuildImgtAnnot.pl` / `BuildDatabaseFa.pl` build steps in the README are unnecessary for
human hg38, and removed them from `setup_trust4_refs.sh`.

Alternative input is FASTQ (`-1/-2` or `-u`), in which case the IMGT file serves both roles.
Irrelevant for us — we have BAMs and converting would cost more I/O, not less.

---

## 5. Every parameter, and which ones actually matter to us

Full surface from `run-trust4`'s own usage text, annotated with our judgement.

### 5.1 Ones we use or should consider

| Flag | Default | Verdict for us |
|---|---|---|
| `-t INT` | 1 | Threads. Helps the k-mer screen only — **not** the BAM read. We saw no gain past a point; consistent with §2.3. |
| `--abnormalUnmapFlag` | off | **Mandatory on AoU.** Passes `-u` to `bam-extractor`. Without it: hard failure, `"Two reads from the unaligned fragment are not showing up together"`. Already in our script. |
| `--od DIR`, `-o PREFIX` | inferred | Output location/naming. Already used. |
| `--clean INT` | 0 | **Should turn on.** `1` = delete intermediates, `2` = keep only AIRR. At 100–200 people the intermediate FASTQs and `_raw.out` files are the bulk of our disk use and we never read them. **Recommend `--clean 1`** — keeps `report.tsv`, `annot.fa`, `cdr3.out`, AIRR; drops the junk. |
| `--stage INT` | 0 | Restart point (0 extract, 1 assemble, 2 annotate, 3 report). **Very useful**: if we ever want to re-annotate with different settings, `--stage 2` skips the 8.5-minute extraction entirely. |
| `--contigMinCov INT` | 0 | Drop contigs covered by < N reads. **The main precision/sensitivity dial.** Default 0 keeps singletons — a contig from a single read pair. Raising to 2 would kill most sequencing-error artifacts at the cost of the rarest true clones. **Untested; worth one controlled experiment.** |
| `--outputReadAssignment` | off | Writes which read went into which contig. Useful only for debugging a suspicious result. |

### 5.2 Ones that don't apply to us — and why, briefly

| Flag | Why not |
|---|---|
| `--barcode`, `--barcodeLevel`, `--barcodeWhitelist`, `--barcodeTranslate`, `--UMI`, `--readFormat` | All single-cell / 10x machinery. **Bulk RNA-seq has no barcodes.** This is the family of options that would have given us **paired α/β chains**, and their inapplicability is the structural limit described in §6.1. |
| `--repseq` | For targeted TCR-seq/BCR-seq libraries (amplicon panels). Ours is unselected whole-transcriptome. Sets `--trimLevel 2 --skipMateExtension`. |
| `--skipMateExtension` | For SMART-seq. We *want* mate extension — it's how we get long contigs. |
| `--skipReadRealign` | Compute saving for barcode/UMI repseq. Would cost us abundance accuracy. |
| `--noExtraction` | FASTQ input only. |
| `--assembleWithRef` | Assemble against IMGT rather than the coordinate file. Niche. |
| `-k INT` | Contig indexing k-mer (default 9). Low-level; no reason to touch. |
| `--minHitLen`, `--cgeneEnd`, `--mateIdSuffixLen`, `--imgtAdditionalGap` | Fine-tuning / other-species. Leave alone. |

### 5.3 Hidden options — not reachable through `run-trust4`

Read off the binaries' own usage strings. Documenting because "not in the README" ≠
"doesn't exist", and one of them matters:

- `trust4 --trimLevel INT` (0 none / 1 trim low quality / **2 trim unmatched**; default 1)
- `trust4 --keepNoBarcode`, `trust4 -c <kmer count file>`
- **`annotator --noImpute`** — disables TCR CDR3 imputation. **The one worth caring about**
  (see §3.2). Reachable by re-running stage 2 by hand.
- `annotator --geneAlignment`, `--outputFormat`, `--needReverseComplement`

---

## 6. Outputs — every file, and which we should keep

| File | Contents |
|---|---|
| `*_report.tsv` | **The one we currently use.** Columns: `read_count`, `frequency`, `CDR3_dna`, `CDR3_amino_acids`, `V`, `D`, `J`, `C`, `consensus_id`, `consensus_id_complete_vdj` |
| `*_airr.tsv` | Same content in **AIRR standard format**. Better for downstream tooling. |
| `*_cdr3.out` | Richer per-assembly table: adds `CDR1`, `CDR2`, **`CDR3_score`**, `read_fragment_count`, **`CDR3_germline_similarity`**, **`complete_vdj_assembly`** |
| `*_annot.fa` | Full assembled contig sequences with per-gene similarity annotations |
| `*_final.out`, `*_raw.out` | Contigs + per-base read weights. Intermediate. |
| `*_assembled_reads.fa`, `*_toassemble*.fq` | Intermediates. Large. Deletable via `--clean 1`. |

**Three columns we are currently throwing away that we should not be:**

1. **`CDR3_score`** — partial (0.00) / imputed (0.01) / real motif strength. Free quality
   filter.
2. **`CDR3_germline_similarity`** — how far the sequence has drifted from germline. For
   **B cells this is a direct readout of somatic hypermutation**, i.e. how much affinity
   maturation that lineage has undergone. That is a real biological variable — "how
   experienced is this B-cell response" — and it is sitting unused in a file we already
   generate.
3. **`complete_vdj_assembly`** — whether the full V-D-J was reconstructed vs. just the
   junction. A completeness/quality flag, and relevant because SCEPTR-style embeddings want
   the V gene, not just the CDR3.

**Two important normalization facts about `frequency`:**
- It is `read_count` normalized **separately within BCR(IG) and TCR(TR)**. So IG
  frequencies sum to 1 and TR frequencies sum to 1 — they are *not* on a common scale.
  Anyone computing a "clonality" number across both families without knowing this will get
  a wrong answer.
- In the amino acid column, `_` = stop codon and `?` = ambiguous base. Both should be
  filtered before any sequence modelling.

---

## 7. Limitations — the honest list

### 7.1 Structural (cannot be fixed by tuning)

1. **No chain pairing.** Bulk data has no cell barcodes. We get a bag of alpha chains and a
   bag of beta chains with no way to know which went together. **This is the biggest single
   constraint on everything downstream** — see `POST_TRUST4_OPTIONS.md` §3.
2. **No cell-type context.** A CDR3 with no idea whether it came from a naive, memory,
   regulatory, or exhausted cell. TRUST4 reports sequences, not states.
3. **Abundance ≈ transcript count, not cell count.** An activated cell transcribes far more
   receptor mRNA than a resting one. So `read_count` conflates *how many cells* with *how
   active they are*. Clonality metrics computed from bulk RNA-seq mean something subtly
   different from the same metrics computed from DNA-based repertoire sequencing, and the
   literature is not always careful about this.
4. **Blood only, one timepoint.** Inherited from the dataset, not the tool.
5. **Depth dependence.** More reads → more distinct CDR3s, essentially without saturating,
   because the true repertoire is far larger than any sample of it. **This is why raw CDR3
   counts are not comparable across samples of different depth** — and it is exactly the
   confound that partially explained our ancestry gap (r = 0.492, ~24% of variance).

### 7.2 Tool-specific / AoU-specific

6. **Alt-contig rescue route is dead on AoU data** (§2.1). Small systematic loss.
7. **Two full BAM passes** (§2.3). Cost, not correctness.
8. **TCR CDR3 imputation on by default** (§3.2). Correctness risk if unfiltered.
9. **Rare chains are unreliable.** Our smoke test: TRD = 14 sequences, TRG = 80. Real
   gamma-delta T cells, correctly detected — but at counts where any per-person statistic
   is noise. **Restrict quantitative analysis to TRA/TRB/IGH/IGK/IGL.**
10. **Read length is fixed at 151 bp and sets the k-mer threshold** (§2.2). We cannot tune
    sensitivity without patching source.

### 7.3 Where TRUST4 is genuinely strong

- **~281% more CDR3s recovered than MiXCR** on bulk RNA-seq, with a reported zero
  false-positive rate; better precision (>18%) and sensitivity (>74%) in 5 of 6 bulk
  samples (TRUST4 paper). **[MED — these are the authors' own benchmarks; MiXCR's own
  benchmarking reaches different conclusions under error-free simulated conditions. Treat
  as "TRUST4 is a defensible, mainstream choice for unselected bulk RNA-seq," not as
  "TRUST4 is objectively best."]**
- **Full-length recovery**: 93 vs 39 full-length IGH sequences vs MiXCR on the same data —
  matters for us, because V-gene identity is needed for embeddings.
- **Works at low depth** — recovered BCRs from as few as 5,000 read pairs in one benchmark.
- **γδ T cells detected** — confirmed in our own smoke test.
- **Takes BAM directly.** No FASTQ conversion step. Given that I/O is our bottleneck, this
  is worth more than it sounds.

---

## 8. Concrete design recommendations

Ordered by value, all cheap:

1. **Add `--clean 1`** to `run_trust4_sample.sh`. Pure disk saving, zero risk.
2. **Switch downstream parsing to the AIRR file** and start carrying `CDR3_score`,
   `CDR3_germline_similarity`, `complete_vdj_assembly` through the aggregation.
3. **Filter before any modelling**: drop `_` (stop codon), `?` (ambiguous), and
   `CDR3_score ∈ {0.00, 0.01}` (partial/imputed). Report how many are dropped — the drop
   *rate* is itself a per-sample quality metric.
4. **Restrict quantitative work to TRA/TRB/IGH/IGK/IGL.** Report TRG/TRD as
   presence/absence only.
5. **Stop using raw CDR3 counts as the recovery metric.** Use rarefaction / downsampling to
   a common read depth. This directly addresses the depth confound rather than
   post-hoc-normalizing it, and it is the standard fix in the repertoire literature
   (rarefaction dropped richness–depth correlation from 0.95 to 0.15 in one published
   simulation). **This is arguably the most important methodological change available to
   us.**
6. **One controlled `--contigMinCov` experiment** (0 vs 2, same people) — the only parameter
   with a real precision/sensitivity trade-off for us. This is a legitimate "experiment that
   adds to our wanted results," so it passes this week's VM rule, but it is not urgent.
7. **Do not** buy a GPU, do not tune `-k`, do not switch to FASTQ input, do not enable
   `--repseq`.

---

## Sources
- TRUST4 source, v1.1.10-r639: `BamExtractor.cpp`, `main.cpp`, `Annotator.cpp`,
  `FastqExtractor.cpp`, `run-trust4`, README — <https://github.com/liulab-dfci/TRUST4>
- Song, L., Cohen, D., Ouyang, Z. et al. *TRUST4: immune repertoire reconstruction from
  bulk and single-cell RNA-seq data.* Nat Methods (2021).
  <https://doi.org/10.1038/s41592-021-01142-2>
- MiXCR comparative benchmarking — <https://mixcr.com/mixcr/guides/mixcr-benchmarking/>
- AIRR Rearrangement schema —
  <https://docs.airr-community.org/en/latest/datarep/rearrangements.html>
- Our own live results: `../results/1000291_trust4_smoke_test.md`,
  `../results/rnaseq_100person_ancestry_writeup.md`
