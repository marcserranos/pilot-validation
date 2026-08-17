# All of Us RNA-seq — full data report

**Written 2026-08-17, from primary sources only, no VM access.** Everything marked
**[HIGH]** comes from AoU's own v9 organization PDF (`../../../reference/How the All of Us
Genomic Data are Organized v9 .pdf`, pages 3–5, 40–52) or from this repo's own live-verified
notes. **[MED]** = inferred from credible sources. **[LOW]** = needs live verification —
these are the items in the "what I need from you" pipeline at the end.

Acronyms are expanded on first use, with a two-line explanation, per the standing rule.

---

## 0. The one-paragraph version

AoU ships **8,980 whole-blood RNA-seq samples** as STAR-aligned BAM files, plus a large
set of pre-computed derived products (expression, splicing, QTL maps, QC metrics). The
derived products are *expression-centric* — they answer "how much of each gene is on."
The BAMs contain something the derived products structurally cannot: the **immune receptor
repertoire**, which is destroyed by the alignment step and absent from every expression
file AoU ships. That gap is the entire reason our workstream exists. This report
characterizes what is in the dataset, what quality it is, and what it can and cannot
support.

---

## 1. Scale, and how RNA-seq compares to the rest of AoU

**[HIGH]** From the v9 PDF's "List of All of Us genomic data" (p.4) and the RNA-seq /
proteomics sections (p.40, p.52):

| Data type | Count | Unit |
|---|---|---|
| Genotyping array | 553,949 | samples |
| Short-read WGS (srWGS) | 535,662 | samples |
| srWGS structural variants | 96,405 | samples |
| **Long-read WGS (lrWGS)** | **14,521** | unique people (15,424 manifest rows) |
| Proteomics (Olink) | 10,096 samples = **9,969 participants** (127 replicates) | samples/people |
| **RNA-seq** | **8,980** | whole-blood samples |

**Terms, since they matter later:**
- **WGS (whole genome sequencing)** — reading a person's DNA. Fixed at birth, identical in
  every cell. Short-read = many tiny fragments (~151 bp); long-read = fewer, much longer
  fragments (~15–17 kb).
- **RNA-seq** — reading the RNA, i.e. which genes are actively transcribed *in that tissue,
  at that moment*. Not fixed; it is a snapshot of cell state.
- **CDR (Curated Data Repository)** — AoU's versioned data release. Current: **v9**,
  controlled-tier ID `C2025Q4R6`, participant data cutoff **2025-01-01**.

**The scale intuition:** RNA-seq is the **second-rarest** modality in AoU — 8,980 out of
535,662 srWGS people is **1.7%**. It is rarer than long-read WGS. This cuts both ways:
it is a genuinely scarce, under-mined resource (good — that is the opportunity), and it
is far too small for the kind of population-scale association work AoU's srWGS supports
(a constraint we must design around, not wish away).

**[LOW] Samples vs. participants for RNA-seq is not stated.** The PDF explicitly reconciles
samples↔participants for proteomics (10,096 → 9,969, naming the 127 replicates) but does
*not* do so for RNA-seq. So 8,980 is a **sample** count and the participant count is
unknown-but-≤8,980. One `wc -l` on the manifest settles it. → item **N1** in §8.

---

## 2. How the data was generated — and why every step matters downstream

**[HIGH]** From the v9 PDF Table 25 and surrounding text (p.40):

| Step | What AoU used |
|---|---|
| Tissue | **Whole blood** (all 8,980 — no other tissue) |
| Library prep | **Watchmaker RNA Library Prep kit with Polaris Depletion** |
| Layout | **Paired-end** |
| Aligner | **STAR**, two-pass |
| Reference | **hg38 excluding ALT, HLA, and decoy contigs** |
| Annotation | **GENCODE v48** |
| Delivered as | **BAM + BAI index**, mark-duplicates applied (`.md.bam`) |

### 2.1 Whole blood — the single most consequential choice

Every sample is whole blood. Blood is ~99% red blood cells and platelets by count, but
those have little/no RNA; the RNA signal is dominated by **white blood cells (leukocytes)** —
which *are* the immune system. That is why a repertoire is recoverable at all: T cells and
B cells are physically present in every one of these samples. In a solid-tissue biobank
(liver, muscle) the repertoire signal would be a faint tumor/tissue-infiltrate trace. Here
it is the main event.

The flip side: **we can only ever see circulating immune cells.** Tissue-resident T cells —
gut, lung, skin, the compartments where much of autoimmune and mucosal immunology actually
happens — are invisible. Any biological claim must be scoped to "the circulating repertoire."

### 2.2 "Polaris Depletion" — and why the read budget is better than it looks

**Depletion** = chemically removing overwhelmingly abundant RNA species before sequencing,
so the sequencer's finite read budget is spent on informative transcripts instead. Two
species dominate blood:
- **rRNA (ribosomal RNA)** — the RNA of the protein-making machinery. Up to ~90% of total
  cellular RNA, and biologically uninformative for expression work.
- **Globin mRNA** — haemoglobin transcripts leaking from reticulocytes (immature red
  cells). Can be 50–70% of the *messenger* RNA in whole blood specifically.

**[MED]** Polaris is a combined rRNA + globin depletion chemistry — the standard choice
for whole blood, and the reason AoU's QC file has a dedicated
`Duplicate Rate of Mapped, excluding Globins` metric and a `Non-Globin Reads` count
**[HIGH — these fields are named in Table 30]**. The presence of those fields is itself
evidence that globin was a first-class concern in the pipeline design.

**Why this matters to us concretely:** depletion is what makes ~100M read pairs per sample
enough to see a repertoire. Without globin depletion, more than half those reads would be
haemoglobin and the effective depth for rare receptor transcripts would collapse. Our
observed 2,013 CDR3s in the smoke test is downstream of this choice.

**A caution, though:** depletion is a *chemical* step with per-sample efficiency variation.
If depletion worked less well in some samples, those samples have fewer usable reads at the
same nominal depth. This is a **candidate confound we have not yet checked** and it is
directly measurable from the QC file's `rRNA Rate` and globin fields. → item **N4** in §8.
This is a real lead for the still-open ancestry-recovery question.

### 2.3 STAR, and the reference — the load-bearing technical fact

**STAR (Spliced Transcripts Alignment to a Reference)** is an RNA-specific aligner. Its
job: take each short read and find where in the genome it came from. Its RNA-specific
skill is handling reads that **span exon-intron junctions** — RNA has the introns cut out,
so a read can legitimately match two genomic stretches thousands of bases apart. A DNA
aligner would fail or soft-clip these; STAR models them explicitly.

**The reference is `hg38` with ALT, HLA, and decoy contigs excluded.** Unpacking that:
- **ALT contigs** — alternate versions of highly variable regions, for people whose genome
  differs too much from the single "primary" reference sequence.
- **HLA contigs** — specifically the alternate sequences of the immune-recognition genes,
  the most variable region in the human genome.
- **Decoy contigs** — junk-magnet sequences included so that reads which belong nowhere
  stick there instead of mis-mapping onto real genes.

**Consequences, and this is where precision matters:**

1. **For HLA work, this is a real handicap** — reads from an HLA allele that diverges from
   the primary chr6 sequence have nowhere better to go and will mis-map or drop. This is
   the same class of problem as the `grch38_noalt` issue already documented for lrWGS in
   the HLA workstream, and it is *worse* here because HLA contigs are explicitly named as
   excluded.

2. **For TRUST4 / repertoire work, it is mostly NOT a handicap — and I want to be precise
   because it is easy to overstate.** The receptor loci (IGH/IGK/IGL, TRA/TRB/TRG/TRD)
   live on the **primary** contigs — chromosomes 2, 7, 14, 22 — which are all present.
   TRUST4's coordinate-based read capture therefore works normally. What the missing
   contigs *do* cost us is TRUST4's **third rescue route** (reads parked on ALT/decoy
   contigs, which TRUST4 re-screens by k-mer — see `TRUST4_DEEP_DIVE.md` §2). On AoU data
   that route captures **nothing, by construction**, because those contigs don't exist in
   this BAM. Net effect: a small, systematic sensitivity loss, not a structural failure.

3. **Two-pass STAR + `--outSAMunmapped Within`** — AoU keeps unmapped reads inside the BAM
   rather than discarding them. This is *fortunate* for us (the unmapped pool is where
   heavily-recombined receptor reads land) and is exactly why the `--abnormalUnmapFlag`
   fix was needed. **[HIGH — verified live in our own smoke test.]**

**[MED] One documentation inconsistency worth knowing about:** the PDF's RNA-seq section
says alignment used STAR, but the `RNA metadata metrics` field list (p.47) includes
*"the Dragen software version used for alignment."* DRAGEN is Illumina's hardware-accelerated
aligner — a different tool. Most likely a copy-paste artifact from the srWGS section of the
same document; possibly a real two-stage pipeline. The BAM `@PG` header is ground truth and
we already read one: our smoke test recorded **STAR v2.7.11b two-pass**. So STAR is correct
for the file we actually touched. Worth a one-line re-check on a second sample. → **N2**.

---

## 3. What AoU actually ships — every deliverable, and whether we care

**[HIGH]** Table 25 (p.40), Tables 26–31 (pp.41–52).

### 3.1 Raw data
| File | Content | Our interest |
|---|---|---|
| `*.Aligned.sortedByCoord.out.md.bam` + `.bai` | STAR-aligned reads, duplicates marked, unmapped reads retained in-BAM | **This is our entire input.** Everything else below is someone else's derived view of it. |

### 3.2 Expression / splicing derived products

**QTL (Quantitative Trait Locus)** — a genetic variant statistically associated with a
measurable trait. **eQTL** = variant associated with a gene's *expression level*.
**sQTL** = variant associated with a gene's *splicing* pattern (which exons get included).
These are "genotype → molecular consequence" maps.

| Deliverable | What it is |
|---|---|
| Cis eQTL (`TensorQTL`, permutation mode) | Top variant per gene, with p-value, effect size (slope), FDR. **Run separately per genetic ancestry group AND on the combined cohort** — notable, and unusually good practice. |
| Cis eQTL nominal stats | Every variant×gene test, not just the top hit |
| Fine-mapped cis eQTL (`susieR`) | 95% credible sets + Posterior Inclusion Probability (PIP) — moves from "correlated" to "probably causal" |
| Cis sQTL, fine-mapped sQTL | Same, for splicing |
| Splicing junction sQTL (`LeafCutter`) | Intron excision ratios, as `chr:start:end:clu_ID_strand` phenotypes |
| **RSEM** outputs | Gene-level and transcript-level expression: `expected_count`, `TPM`, and `IsoPct` (what % of a gene's expression each isoform contributes) |

**Our relationship to these: they are the approach your supervisors explicitly rejected.**
Worth stating plainly rather than re-litigating — every one of these products is a
refinement of "how much of gene X is on," and none of them contains a single receptor
sequence. They are, however, **useful as covariates and as sanity checks** (e.g. RSEM TPM
of *CD3E*, *CD19*, *MS4A1* gives an independent estimate of T-cell and B-cell abundance,
which is a much better cell-composition control than the chain-ratio proxy we improvised
in `check_ancestry_confounds.py`). → this is a real, cheap improvement, item **N5**.

### 3.3 Auxiliary / QC files — the ones we should actually be using

**`RNA metadata metrics`** — per sample: Research ID, alignment rate %, **RQS quality
score**, **% mRNA bases**, **% ribosomal bases**, **number of aligned read pairs**, mean
insert size, aligner version. *(We already use `research_id` and the read-pair count as our
depth variable.)*

**`RNA-SeQC 2 metrics`** — the deep QC table, ~70 fields. The ones that matter for us:

| Field | Why we care |
|---|---|
| `Expression Profiling Efficiency` | Fraction of *all* reads that land on exons and are usable. AoU calls it "the ultimate benchmark metric." A single number for "how much of this sample was worth sequencing." |
| `rRNA Rate`, `rRNA Reads` | Direct readout of how well depletion worked. **Untested confound.** |
| `Non-Globin Reads`, `Duplicate Rate ... excluding Globins` | Same, for globin. **Untested confound.** |
| `Estimated Library Complexity` | How many genuinely distinct molecules existed before PCR amplification. **This is arguably a better denominator than raw read count** for repertoire recovery — 100M reads off a low-complexity library is not the same as 100M reads off a rich one. **Untested.** |
| `Mean/Median 3' bias` | Whether RNA was degraded. 1 = perfectly uniform coverage, →0 = severe 3'-end degradation. Degraded RNA truncates transcripts — and **CDR3 sits in the middle of a receptor transcript**, so 3' bias should plausibly hurt repertoire recovery specifically. **Untested, and mechanistically the most interesting of the lot.** |
| `Genes Detected` | Blunt overall-richness proxy |
| `Mapping Rate`, `Unique Rate of Mapped`, `Intronic/Intergenic Rate` | Standard alignment sanity |
| `Fragment GC Content Mean/Std` | PCR/library GC bias |

**This is the single biggest actionable finding of this report.** We diagnosed the ancestry
recovery gap using *one* quality variable (RQS) and *one* depth variable (aligned read
pairs). AoU ships **five more** plausible technical explanations — library complexity, 3'
bias, rRNA rate, globin rate, expression profiling efficiency — every one of which is
already computed, sitting in a file we have already located, costing zero compute to check.
Running that check is item **N4**, and it is the highest value-per-minute experiment
available to us right now.

---

## 4. Quality — what we know, and honestly what we don't

**[HIGH] From our own live data (100-person cohort, `rnaseq_*` files in `../results/`):**
- Depth: most people **95–145M aligned read pairs**, tail past 350M. Mean by ancestry
  ranged **121.4M (EAS) → 146.2M (AFR)**, a 1.20× spread.
- **RQS: 7.82–8.00 across all six ancestry groups** — remarkably flat, and uncorrelated
  with recovery (r = 0.061). RNA quality is *not* a differentiator in this cohort.
- Recovery: **2,013 CDR3s** in the smoke test; cohort means **4,538 (EAS) → 7,406 (AFR)**.

**[HIGH] What is genuinely good about this dataset:**
- Uniform protocol — one tissue, one library kit, one aligner, one reference, one
  annotation version. No cross-site batch heterogeneity in the *molecular* pipeline.
- Depth is high for bulk RNA-seq. 100M+ read pairs is roughly 2–5× a typical expression
  study, which is precisely why a repertoire is recoverable.
- Unmapped reads retained in-BAM — not universal, and essential for us.
- QC is unusually thorough and per-sample.

**[MED/LOW] What we do not know and should stop assuming:**
- **How the 8,980 were selected.** Not documented in the PDF, and I could not find it in
  public sources. This is *not* a minor point: if the RNA-seq sub-cohort was enriched for
  anything (a disease, a site, a recruitment wave), every population-level claim we make
  inherits that enrichment. → **N3**, and it is a prerequisite for Task 2's interpretation.
- **Collection-site / batch effects.** Still the one confound listed as unchecked in
  `rnaseq_100person_ancestry_writeup.md`. Site is likely obtainable via OMOP.
- **Time-of-draw, fasting state, acute illness at draw.** A repertoire snapshot moves with
  recent infection. Almost certainly not recorded — a permanent, irreducible noise floor
  on any repertoire phenotype, and worth stating openly rather than discovering later.

---

## 5. HLA content specifically

Restating precisely, because this has been muddled before:

- **AoU ships official HLA calls for srWGS** (HLA-HD + Polysolver + OptiType ensemble),
  at `v9/wgs/short_read/snpindel/hla_variants/hla_genotypes.tsv`, ~331 MiB, 30 gene-pair
  columns. **[HIGH, live-verified in the HLA workstream.]** HLA calling is solved and is
  Marc's workstream.
- **The RNA-seq reference explicitly excludes HLA contigs**, so RNA-seq is a *poor*
  substrate for HLA typing — worse than srWGS, which at least had them.
- **Therefore: do not call HLA from RNA-seq.** Use AoU's srWGS calls as **labels**. The
  RNA-seq BAM's job is to supply the repertoire; the HLA type comes from elsewhere, already
  computed, free.

That division of labour is the cleanest thing about the current design and is worth
defending if it gets questioned.

---

## 6. What this dataset can and cannot support

**Can support:**
- Descriptive repertoire characterization at n≈9,000 — genuinely novel; nobody has mined
  this.
- Repertoire ↔ HLA association, using AoU's own HLA calls as labels.
- Repertoire ↔ phenotype description in any sub-cohort with enough cases (Task 2).
- Technical/methodological work on recovery bias — which is what we have been doing, and
  which is a legitimate contribution in its own right given no one has characterized
  TRUST4-on-AoU before.

**Cannot support:**
- **Paired α/β chains.** Bulk RNA-seq has no cell barcodes, so there is no way to know
  which alpha chain was in the same cell as which beta chain. This is a hard structural
  limit, not a tuning problem, and it constrains every downstream embedding option
  (see `POST_TRUST4_OPTIONS.md` §3).
- **Tissue-resident immunity.** Blood only.
- **Longitudinal / dynamic immunology.** One draw per person.
- **Well-powered genome-wide association on repertoire traits.** n≈9,000 with high
  within-group variance; the 100-person pilot already returned p = 0.18 on a 1.63× effect.
- **Rare-disease work.** 8,980 people means most specific diseases will have single-digit
  case counts. This is precisely the constraint Task 2 exists to quantify.

---

## 7. Confirmed paths (live-verified, this repo)

```
Bucket (requester-pays; needs -u <billing-project> or gcsfuse --billing-project):
  gs://vwb-aou-datasets-controlled/

RNA-seq manifest      v9/multiomics/rnaseq/manifest.tsv
    columns: sampleid, research_id, markduplicates_bam_file_path, markduplicates_bam_index_path
RNA-seq metadata      v9/multiomics/rnaseq/rnaseq_metadata.tsv     <- RQS, depth, alignment rate
RNA-seq QC            v9/multiomics/rnaseq/rnaseqc2/               <- the ~70-field QC table
RNA-seq expression    v9/multiomics/rnaseq/rsem/                   <- TPM, counts, isoform %
RNA-seq QTLs          v9/multiomics/rnaseq/{eqtl,sqtl}/
Physical BAMs         pooled/multiomics/v9_base/rnaseq/bam/<research_id>.Aligned.sortedByCoord.out.md.bam

lrWGS manifest        v9/wgs/long_read/manifest.tsv                <- 15,424 rows / 14,521 people
Genetic ancestry      v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv
AoU native HLA calls  v9/wgs/short_read/snpindel/hla_variants/hla_genotypes.tsv

BigQuery CDR: project wb-silky-artichoke-2408, dataset C2025Q4R6 (fully-qualified paths only;
              not exposed via shell env vars)
```

---

## 8. What I need from you — numbered, each one small

Every item below is a short, cheap command. None requires a batch run. Ordered by value.

| # | Question | How | Why it matters |
|---|---|---|---|
| **N1** | Is 8,980 samples = 8,980 people? | `wc -l` on `rnaseq/manifest.tsv`; count unique `research_id` | Every denominator in this report and in Task 2 |
| **N2** | Confirm aligner on a 2nd sample | `samtools view -H <bam> \| grep '@PG'` | Resolves the STAR-vs-DRAGEN doc inconsistency |
| **N3** | Is there any doc on how the 8,980 were selected? | Workbench UI → Data Dictionary / release notes; also just ask your supervisors | Determines whether cohort-level claims generalize |
| **N4** | **Pull the RNA-SeQC2 QC table** for our 100-person cohort | `ls v9/multiomics/rnaseq/rnaseqc2/`, then join on `research_id` | **Highest value.** 5 untested confounds for the ancestry gap, zero compute |
| **N5** | Pull RSEM TPM for `CD3E`, `CD19`, `MS4A1`, `PTPRC` | `rsem_genes_tpm.txt.gz`, subset rows | Replaces our improvised cell-composition proxy with a real one |
| **N6** | Row count + overlap of LR and RNA-seq manifests | `scripts/check_lr_rnaseq_overlap.py` (written, ready) | **Gates all of Task 2** |

Items N1, N2, N6 are one-liners. N4/N5 are file locates + a join. Nothing here needs a
long-running job, which matches this week's "VM only if an experiment adds to results" rule —
these are all *lookups of data AoU already computed*, not experiments.

---

## Sources
- AoU, *How the All of Us Genomic and multi-omic data are organized*, v9 (2026) — local
  copy at `../../../reference/How the All of Us Genomic Data are Organized v9 .pdf`,
  pp. 3–5, 40–52. Primary source for every **[HIGH]** claim in §1–3.
- AoU *Genomic Research Data Quality Report* (v7 edition) — lrWGS selection criteria and QC
  thresholds.
- This repo: `../../../reference/AOU_DATA_ACCESS_NOTES.md`, `../../../context/ENVIRONMENT.md`
  (quirk #13), `../results/*` — all live-verified facts.
