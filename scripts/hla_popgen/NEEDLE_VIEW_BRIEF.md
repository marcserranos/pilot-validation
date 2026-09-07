# Brief: the "needle view" (variant position vs. CDS)

Written 2026-09-07 to hand this workstream to the next agent/session. Answers who/where/when/state
for the existing lollipop/needle plot showing whether novel HLA sequence differences fall inside or
outside the CDS, plus the one open gap and one stale citation to fix.

## TL;DR

The needle/CDS-position view is **`09_mutation_topology.py`**, not the two newest untracked files
(`10_allele_ancestry_geometry.py`, `tests/test_allele_geometry.py`) — those are a separate,
unrelated ancestry-centroid-geometry workstream that only shares a CDS-hashing helper (for
allele-cluster ID joins), nothing about variant position.

## Who / where / when

- **Author**: Marc Serrano, sole author on all relevant commits.
- **File**: [09_mutation_topology.py](09_mutation_topology.py) — core logic at
  `parse_codon_positions()` (~line 85) and `plot_lollipops()` (~line 130).
- **Commits** (all 2026-09-06):
  - `719e963` 12:34 — script created ("Singleton-inclusive QC, discovery-rate curve fits, and
    allele-space embedding comparison")
  - `6c32f3d` 17:06 — fixed a per-codon novelty bug in the position parsing
  - `a8de0a7` 17:12 — reverted per-codon synonymous/non-synonymous classification after a
    real-data regression (see "Known gap" below) — this is the current state on `main`
- **Output artifacts**: `reports/hla_popgen/09_mutation_topology/lr/mutation_topology_classical.png`
  and `09_mutation_topology_report.md`, timestamped 2026-09-06 17:19 — so this has been run against
  real cohort data, not just fixtures.
- No commit anywhere in the repo (`git log --all`) mentions "needle" literally; this file is the
  only hit for the concept.

## What it does

A gene-faceted lollipop/needle plot (2×4 grid across classical HLA genes): one vertical stem per
distinct codon position along the CDS, stem height = recurrence-weighted count of haplotypes
carrying a novel difference at that codon, color = the call's `novelty_class`
(`protein_altering` / `synonymous` / `beyond_cds` / `undetermined`). **`beyond_cds` is literally the
"outside the CDS" bucket** — that's the inside-vs-outside-CDS distinction being asked about.

Position derivation: `cds_mut` (Table 1) is `[ref allele]|[minimap2 cs string]|[aa diff]`. The `cs`
string's `:N` (run of N identical bases), `*xy` (substitution), `+seq`/`-seq` (indel) tokens are
summed into a running nucleotide offset, converted to 1-based codon position via `nt_offset // 3 +
1`, deduped per codon.

## Known, documented gap (do not silently redo this without reading the docstring first)

An earlier version tried a **true per-codon synonymous/non-synonymous** call by zipping the aa-diff
field's `REFaa(REFcodon)<OBSaa(OBScodon)` entries against the cs-string's deduped codon positions.
It passed against `tests/make_fixtures.py`'s synthetic data (clean single-letter AA codes, always
paired) but broke on real Immuannot output: multi-letter/garbled AA tokens (`Tre`, `Rrg`),
indel-consolidated multi-residue entries with variable-length codon strings, and single-sided entries
with no diff pair at all. Result: 100% `n_position_aa_mismatch` on ~5,000 real calls (0 usable).
**Reverted in `a8de0a7`.** The script's own docstring (lines ~20–45) has the full real-data example
string and reasoning — read it before attempting this again.

Today's coloring is therefore the coarser row-level `novelty_class` label, not a true per-codon call.
That finer-grained per-codon classification remains the actual open gap if "the needle view" needs
to get more precise about *which* codon-level change is synonymous vs. not.

## Also worth fixing

The docstring (line ~3) attributes the lollipop/needle design idiom to "the 2026-09
visualization-literature review in `research/VIZ_LIT.md`'s addendum." **This addendum does not
exist** — `research/VIZ_LIT.md` (348 lines) has zero hits for "needle", "lollipop", "addendum",
"topology", "cBioPortal", "ProteinPaint", or "trackViewer". Either the addendum was written
elsewhere and never merged in, or it was removed later. Worth reconciling or just fixing the
citation before trusting it as design rationale.

## Test coverage

**None.** `scripts/hla_popgen/tests/` has `test_extraction.py`, `test_figures.py`, `test_novel.py`,
and (new, untracked) `test_allele_geometry.py` — nothing named for `09_mutation_topology`. If
picking this up, `tests/test_mutation_topology.py` doesn't exist yet.

## Suggested next steps, in order

1. Read the full docstring in `09_mutation_topology.py` (lines 1–45) — it's the primary source of
   truth on what was tried and why it failed.
2. Decide whether per-codon synonymous/non-synonymous classification is actually needed, given the
   real-data parsing landmines documented above — if so, this needs a real-format-aware parser
   (handle multi-letter AA tokens, variable-length indel spans, single-sided entries), validated
   against real `cds_mut` strings pulled from the VM, not just fixtures.
3. Write `tests/test_mutation_topology.py` covering `parse_codon_positions()` against both the
   synthetic fixture format and the real-format edge cases quoted in the docstring.
4. Fix or remove the stale `VIZ_LIT.md` "addendum" citation.

## Extension (2026-09-07): gene-wide diversity/richness track, inside vs. outside CDS

Follow-up ask: not just "which codon is this one variant at" (existing `09`), but a per-position
variability/mutation-rate/richness track across the **whole gene body** (introns + UTR + CDS), to
compare inside-CDS vs. outside-CDS diversity as a selection-pressure signal.

### Get the hypothesis direction right first

The obvious framing — "CDS under purifying selection → lower diversity inside CDS" — is **likely
backwards for classical HLA class I/II peptide-binding-groove exons**. These are the textbook case
of **diversifying/balancing selection** (already flagged elsewhere in this codebase:
`research/NOVEL_LIT.md` §4's Ewens–Watterson neutrality-rejection results, and
`04_allele_saturation.py`'s own module docstring). Under balancing selection, diversity should be
**elevated**, not suppressed, at the groove-encoding codons relative to flanking sequence — expect a
**non-monotonic, position-dependent** pattern (high at exon 2/3 groove codons, closer to neutral
elsewhere), not a single CDS-vs-non-CDS mean comparison. Collapsing "CDS" to one bucket would wash
out the exact signal that's biologically the point. Frame this as a diversifying-selection test
(CDS ≥ flanking, concentrated at specific exons), not a purifying-selection test.

### What data plumbing exists vs. needs to be built

Nothing today gives per-position diversity outside the CDS — the finest granularity anywhere in the
pipeline is per-codon (CDS-only, via `cds_mut`) or per-gene aggregate (`template_distance` /
`template_distance_per_kb`, SCHEMA.md lines 92-99 — a single integer per haplotype over the *whole*
gene span, introns+UTR+CDS collapsed). No script parses a `.gtf`/`.gff`/`.bed` for exon/CDS/UTR
boundaries today; `01_extract_rich.py:264-291`'s GTF loop only handles `gene`/`transcript` rows
(exon/CDS/UTR rows in the same raw `hap{N}.gtf.gz` are already there, just unparsed — cheap to add).

To get an actual per-position track across the full gene, the needed raw material is
`mm2.ipd.gen.paf.gz` — the raw minimap2 PAF file with the full gene-level `cs:Z:` alignment tag
(vs. `cds_mut`'s CDS-only `cs` substring). **CONFIRMED present** (updated 2026-09-07, resolving the
"ambiguous" flag below): `reports/hla_popgen/recon_report.json`/`.md` (a real `00_recon_vm.py` run,
50-person sample, Sep 5) shows `hap1/mm2.ipd.gen.paf.gz` and `hap2/mm2.ipd.gen.paf.gz` at **100%
presence**, mean ~780KB/hap, ~19GB extrapolated total across the full ~12,000-person cohort per hap.
`mm2.ipd.cds.paf.gz`, `gene.filtered.paf`, and `cds.fa.gz` are also 100% present in the same report.
**But `01_extract_rich.py` does not read it** — confirmed by full read of the file: its only
filesystem reads are `hap{1,2}.gtf.gz` (line 459/258-259); no script in `scripts/hla_popgen/`
currently parses PAF content into any table (`03_novel_alleles.py`/`09_mutation_topology.py` only
parse the CDS-only `cs` substring already embedded in the GTF's `cds_mut` attribute, not the
standalone PAF file). So: **the data exists, but the extract does not yet reach it** — this is new
extraction work, not a plotting-layer change.

### Recommended metric

**SNP/variant density per window**, as the primary, cheaply-derivable signal — directly reuses the
`cs`-string-walking pattern already proven in `parse_codon_positions()`
(`09_mutation_topology.py:82-107`, `CS_TOKEN_RE` at line 74: `:N`/`*xy`/`+seq`/`-seq` → running nt
offset), just run against the PAF's gene-level `cs:Z:` tag instead of `cds_mut`, and *without* the
÷3 codon collapse (intron/UTR positions aren't codon-structured). Per-site heterozygosity
(`He = 1 - Σpᵢ²`, already sketched per-population in `research/VIZ_LIT.md` §1.10, lines 151-161) is
a good secondary, more rigorous cross-check once per-position allele frequencies are available.

**dN/dS or pN/pS** would be the scientifically sharpest test of selection *direction* (the standard
metric in the HLA-groove-selection literature, going back to Hughes & Nei — not yet cited in
`research/`), but it's blocked on the same per-codon synonymous/non-synonymous classification
problem documented above as broken and reverted on real data (`a8de0a7`). Don't reopen that unless
dN/dS specifically is wanted; default to variant density / heterozygosity instead.

### Other gaps to flag

- Multi-copy genes (DRB3/4/5 etc.) complicate joining a positional track back to a specific gene
  copy (`IMMUANNOT_GTF_SPEC.md` part D "Caveat", part G) — fine for the near-always-single-copy
  classical class I genes, a real complication for class II paralogs.
- This is new VM-side extraction plus a new aggregate table, not a same-session plotting change —
  expect a new Table (e.g. "per-position variant density") in SCHEMA.md, a new extraction script
  parsing PAF `cs:Z:` tags at cohort scale, and only then a plotting script analogous to `09`.

### Reusable pointers

- `09_mutation_topology.py:74,82-107` — cs-string walker to adapt (drop the codon collapse).
- `01_extract_rich.py:264-291` — GTF parse loop to extend with exon/CDS/UTR feature-type branches.
- `reference/IMMUANNOT_GTF_SPEC.md` part A (exon/CDS/UTR row coordinate semantics) and part E
  (PAF file retention, marked ambiguous/unconfirmed).
- `SCHEMA.md:92-99` — canonical `template_distance`/`gene_start`/`gene_end`/`cds_distance`/
  `cds_mut` semantics, must not be violated by any new column.
- `research/VIZ_LIT.md` §1.10 (lines 151-161) — heterozygosity formula, reusable per-position.
- `research/NOVEL_LIT.md` §4 (lines 139-151) — balancing-selection caution, cross-reference rather
  than re-derive.

## Pipeline architecture (2026-09-07): ground-up design, not a bolt-on to `01_extract_rich.py`

Decision: this is heavy enough (~24,000 PAF files, ~19GB/hap-slot, cohort-scale) and structurally
different enough (positional aggregation vs. row-per-allele extraction) to warrant its own
pipeline rather than extending the existing GTF-only extraction script. Below is the design,
informed by what's actually proven to work/fail in this environment so far (see pointers).

### Confirmed source data (per haplotype, `~/pipeline_outputs/people/<person_id>/immuannot_output/hap{1,2}/`)

`mm2.ipd.gen.paf.gz` — raw, unfiltered gene-level PAF from `minimap2 -cx asm5 --cs --end-bonus=10
<contig> <reference_gene.fa>` (confirmed against Immuannot's real `annot.IPD.sh` source, not just
inferred). **Target = contig, query = the IPD/IMGT reference gene FASTA.** This matters: walking the
`cs` string while tracking the *query* offset (not target/contig offset) gives position **relative
to the shared reference gene sequence** — the common coordinate frame needed to aggregate positions
across ~12,000 different people's differently-coordinate contigs. Tracking target offset instead
would be useless for cross-person aggregation. **This directionality is inferred from the source
command's argument order and needs a one-line confirmation against a real PAF's `qname`/`tname`
fields on the VM before writing the parser**, not assumed blind.

### Known unresolved risk: candidate vs. accepted alignments

`mm2.ipd.gen.paf.gz` is described as the **full candidate universe**, not just the winning call —
Immuannot's `searchTemplate.py` picks a best match per gene copy from multiple candidate rows.
Naively aggregating every row in this file would mix in discarded/losing candidate alignments and
corrupt the diversity signal. **Must cross-reference PAF rows against the final accepted
gene/copy_index already computed in Table 1 (`01_extract_rich.py` output)** to select only the
rows backing the accepted call, before counting anything toward the diversity track. This needs
verification against a real file (is there an unambiguous key — e.g. matching `qname`/candidate
rank — connecting a PAF row to the specific gene copy Table 1 already resolved?) before assuming
"just take the top-scoring row per gene" is correct.

### Environment facts that shape the design (from live-pattern research, not guesses)

- **VM specs are not reliably documented** — conflicting historical figures in
  `context/ENVIRONMENT.md` (4 vCPU/25GB/88GB vs. a separately-observed 138GB disk vs. an
  undocumented `n2-highcpu-32` test VM) and zero hardware-detection code anywhere in
  `scripts/hla_popgen/`. **First live step: `nproc`, `free -h`, and a targeted `df -h` on just
  `~/pipeline_outputs`** (never `du -sm`/`df` over a whole FUSE-mounted tree — `ENVIRONMENT.md`
  quirks #25/26 document repeated 15min–2hr hangs from that; always `timeout 5` any disk check).
- **Existing precedent uses `ThreadPoolExecutor`, not multiprocessing**, for the one concurrent
  script that exists (`01_extract_rich.py:616`, reasoned as I/O/gzip-bound). PAF parsing here does
  meaningfully more per-file CPU work (full `cs`-string walk over a larger file) — **don't assume
  threads are right; time a small sample (e.g. the same 50-person recon set) with both
  ThreadPoolExecutor and ProcessPoolExecutor before committing**, since which one wins depends on
  actual per-file parse cost, which is currently unmeasured.
- **Checkpointing precedent** (`01_extract_rich.py:602-603` `already_done_person_ids()` reading
  back its own output; `:581-582` flush every `--checkpoint-every`, default 200) — adapt this same
  idea, but for array accumulators: periodically serialize partial per-gene position-count arrays
  to disk (cheap — arrays are bounded by gene length, not cohort size) and merge-resume, so a run
  interrupted by VM contention with the other concurrent agents loses at most one checkpoint
  interval, not the whole job.
- **PID-file lock** — no `hla_popgen/` script currently takes one, but `ENVIRONMENT.md` quirk #23
  documents this exact pattern for exactly this situation (concurrent invocations corrupting shared
  intermediate state). Given other agents may be active on the same VM right now, the new script
  should refuse to start a second concurrent run against the same output directory.
- **Never copy raw PAFs off the VM** — matches existing convention (`RUNBOOK.md:220-223`: per-person
  raw data, ~95% of which is `mm2.*.paf.gz` scratch, is explicitly never pushed to the GCS share
  bucket; only small aggregate `.tsv`/output files are). The new pipeline should parse PAFs in place
  (streamed, never fully decompressed to local scratch) and only ever write/export the small
  per-gene aggregate arrays and a summary table.
- **PAF column layout has never actually been inspected in this codebase** — no doc states the 12
  mandatory PAF columns or where `cs:Z:` sits among optional tags. Parse by **tag name, not fixed
  column index** (`cs:Z:` position among optional trailing tags isn't guaranteed) — but confirm
  against one real `zcat`'d file on the VM before writing the parser, same as the target/query
  question above.

### Proposed shape (small footprint by design)

Per-gene output is bounded by **gene length** (a few kb, at most ~16kb genomic for DRB1), not
cohort size — so the aggregate accumulator (position → {n_haplotypes_observed, n_variant_events,
in_cds: bool}) is tiny in memory/on disk regardless of processing 12,000 people. The heavy part is
purely the embarrassingly-parallel per-file read+parse, which is exactly what chunked
worker-pool + periodic checkpoint-and-merge is for. No need for pandas/polars/parquet at this
scale — plain per-gene numpy arrays (or even plain dict-of-Counter, given gene-length bound) are
sufficient; reach for parquet only if the summary table itself needs to be shared/joined broadly.

### Recommended validation order (do this before scaling to 12k people)

1. Live VM probe: `nproc`, `free -h`, targeted `df -h ~/pipeline_outputs` (timeout 5).
2. `zcat` one real `mm2.ipd.gen.paf.gz` — confirm column layout, confirm `qname`/`tname` identity
   matches the target=contig/query=reference-gene assumption above, confirm `cs:Z:` tag presence
   and format on real data.
3. Resolve the candidate-vs-accepted-alignment question against Table 1 for one real person.
4. Prototype on the existing 50-person recon sample (`reports/hla_popgen/recon_report.json`'s
   sample) — sanity-check the resulting per-position track for HLA-A against known exon/CDS
   boundaries before trusting it.
5. Time the 50-person prototype under both threading and multiprocessing to pick a worker model,
   and extrapolate runtime for the full ~12,000-person cohort before committing to a full run.
6. Scale up with checkpointing + PID lock in place.

### Gene prioritization for first pass

Start with the classical, single-copy, best-annotated genes with the strongest expected
diversifying-selection signal at the peptide-binding groove, avoiding multi-copy-paralog ambiguity
(see "Other gaps to flag" above re: DRB3/4/5):

- **Class I**: HLA-A, HLA-B, HLA-C — single-copy, highest polymorphism, well-characterized exon 2/3
  groove-encoding region, classic Hughes & Nei dN/dS-selection loci.
- **Class II**: HLA-DRB1 — single-copy, most polymorphic class II gene, exon 2 groove.

Treat these four as the first-pass gene set (matches typical clinical/PRS transplant panels too);
expand to DQA1/DQB1/DPA1/DPB1 (single-copy but more linkage complexity) once the first pass is
validated, and leave DRB3/4/5 (multi-copy) and C4 (separate alignment approach entirely, per
Immuannot's `annot.C4.sh`) for a later, explicitly scoped extension.

## Live VM validation (2026-09-07, via computer-use on `AoU_Jupyter_ComputeEngine_20260805_big_run`)

Every open design question above is now empirically resolved against real data (person `1000234`,
hap1). **Do not re-derive these — build directly on them.**

- **VM specs** (previously undocumented/conflicting): `nproc`=4, `free -h`=31Gi total/24Gi free,
  `df -h`=2.0T total/1.4T available (27% used) — this is the "cheap resized VM" (4 vCPU matches
  `ENVIRONMENT.md`'s figure), disk is NOT the constraint the old 88GB/138GB figures implied.
  **Another agent's session was active in a separate terminal (pts/2) at the time of this check**
  (its own `tee`-based session log, PID 526/545) — plan worker-pool sizing to leave headroom (e.g.
  2-3 workers, not all 4), and this VM instance does have other concurrent users in general.
- **gcsfuse already mounted**, using billing project `wb-cordial-leechee-9743` (confirms
  `RUNBOOK.md`'s value is the one actually in use; `ENVIRONMENT.md`'s `wb-glacial-potato-8710` may
  be stale — not yet reconciled in the docs themselves, but this session's live evidence favors
  `wb-cordial-leechee-9743`).
- **PAF format, target/query identity — confirmed exactly as designed.** Real row:
  ```
  HLA-A*01:01:01:01  3503  0  3503  -  h1tg002663l:1-895500  895500  616434  619951  3388  3524  60
  NM:i:136 ms:i:1109 AS:i:1059 nn:i:0 tp:A:P cm:i:177 s1:i:2073 s2:i:1139 de:f:0.0323 rl:i:0
  cg:Z:877M3I235M1I185M3I487M17D154M4D1558M cs:Z::17*ag:31*ga:52*ga:55*gc:29*ac:77*a...
  ```
  Column 1 (`qname`) = reference IPD allele name, `qlen` = full reference gene length; column 6
  (`tname`) = contig region. **Confirms query = reference gene, target = contig** — walk the `cs`
  string tracking query-side offset for reference-gene-relative position, exactly as designed.
  Standard 12-col PAF + optional tags (`cg:Z:` cigar precedes `cs:Z:` in the tag list observed, but
  still parse by tag name, not fixed position, per the original caution).
- **Candidate-vs-accepted risk — real and now solved.** One hap1 file had **32,616 candidate PAF
  rows** across all HLA/immune genes tested (not just the 4 classical genes — `HLA-A/B/C/DMA/DMB/
  DOA/DOB/DPA1/DPA2/DPB1/DPB2/DQA1/DQA2/DQB1/DQB2/DRA/DRB1/DRB3/DRB4/HLA-E/...` all present), e.g. 5
  different `HLA-A*01:01:01:0x` candidates all mapped to the same contig region. **The join key is
  the GTF's own gene row**: `hap1.gtf.gz`'s `gene` feature carries `template_allele
  "HLA-A*26:01:01:01"` and `template_distance 0`, with contig coords (616435-619951, 1-based)
  matching the PAF target region (616434-619951, 0-based) exactly. Filtering the PAF for
  `qname == template_allele` returns exactly one matching row:
  `HLA-A*26:01:01:01 ... NM:i:0 ... cs:Z::3517` (one giant identical-run token, `NM:i:0` confirming
  the perfect match implied by `template_distance 0`). **So: extend the GTF gene-row parse
  (`01_extract_rich.py:264-291`-style) to also capture `template_allele`, then filter
  `mm2.ipd.gen.paf.gz` to `qname == template_allele` per gene per hap — trivial, unambiguous join,
  no fuzzy matching or rank-picking needed.**

### Updated next steps (supersedes the "Recommended validation order" section above)

Steps 1-3 of that section are now DONE (this VM probe). Remaining before scaling to the cohort:

1. Write the extraction script: extend GTF parse for `template_allele` (and exon/CDS/UTR rows for
   boundary overlay), then per-gene-per-hap filter `mm2.ipd.gen.paf.gz` by that allele name, walk
   `cs:Z:` on the query side, accumulate into the small per-gene position arrays.
2. Prototype on the existing 50-person recon sample; sanity-check the HLA-A track.
3. Time thread vs. process pool on that sample before deciding the worker model; keep worker count
   conservative given other concurrent VM usage observed live.
4. Add the PID-file lock and checkpoint-every pattern before any full-cohort run.
5. Scale to the full cohort.

## Results (2026-09-07): the script is built, ran, and the CDS-vs-flanking result holds up

`scripts/hla_popgen/11_gene_diversity_track.py` exists and is committed to this file's design.
Two real bugs surfaced during iteration on real VM data -- both fixed, documented in the script's
own comments, and worth reading before touching this code again:

### Bug 1: raw positional track isn't cross-person-comparable (the coordinate-frame problem)

Confirmed empirically on a 50-person prototype: HLA-DRB1 haplotypes matched templates ranging
10,850-16,110bp (>5kb spread), purely from intron-length polymorphism between different reference
alleles. Summing raw query-offset positions across haplotypes that matched different-length
templates smears/misattributes position. **Fix**: pivoted the PRIMARY result away from a raw
position track to a **CDS-vs-flanking density comparison** -- per haplotype, count variants inside
vs. outside its own (locally-computed, always-correct) CDS ranges, and the bp denominators (CDS bp
straight from the GTF, no translation needed; total bp = qlen). Summing these (count, bp) pairs
across haplotypes needs no shared coordinate axis at all. See `cds_vs_flanking_density()` and the
"Density accounting" comment in `process_person_gene_hap()`/`aggregate_gene()`. The raw positional
track (`variant_counts`, `cds_query_ranges_consensus`/`_spread`) is still computed and written to
each gene's JSON as a secondary/exploratory output -- useful for QC (e.g. `cds_spread` is itself a
diagnostic of this exact template-length heterogeneity) but NOT to be plotted as a trustworthy
absolute-position track without further work (per-template stratification or a proper reference-
to-reference re-alignment, neither implemented).

### Bug 2: exact binomial test overflowed at DRB1's cohort scale

The first version computed the density comparison's significance test by hand
(`sum(math.comb(n, i) * p**i * (1-p)**(n-i) for i in range(n+1))`). This crashed
(`OverflowError: int too large to convert to float`) on a 500-person HLA-DRB1 run once
`total_variants_all` reached the low thousands -- `math.comb(n, i)` produces an integer with more
digits than a float can hold long before the multiplication normalizes it back down. **Fix**:
replaced with `scipy.stats.binomtest` (works in log-space internally, confirmed available on the
VM, scipy 1.17.1), plus `float(...)` around its `.pvalue` since `json.dump` can't serialize
`numpy.float64` -- a second real bug (JSON serialization) caught in the same pass.

### The actual result, n=500 people (~966-969 haplotypes per gene)

| gene | density in CDS (per kb) | density outside CDS (per kb) | ratio | binomial p |
|---|---|---|---|---|
| HLA-A | 0.192 | 0.112 | 1.71x | 1.6e-08 |
| HLA-B | 0.275 | 0.150 | 1.83x | 1.2e-14 |
| HLA-C | 0.236 | 0.110 | 2.14x | 1.0e-18 |
| HLA-DRB1 | 0.205 | 0.172 | 1.19x | 3.5e-02 |

**All four genes show variant density concentrated INSIDE the CDS, not outside** -- the opposite
of the naive purifying-selection expectation, consistent with diversifying/balancing selection at
the peptide-binding groove exons, exactly the hypothesis this brief laid out up front. Class I
genes (A/B/C) show a strong, highly significant ~1.7-2.1x enrichment. **HLA-DRB1 is a real,
biologically-explicable outlier**: its enrichment is much weaker (1.19x, only nominally
significant) because its CDS is a much smaller fraction of total gene length than the class I
genes' (`total_cds_bp / total_qlen_bp` = 5.9% for DRB1 vs. 26-31% for A/B/C, from the same JSON
output) -- DRB1's introns (especially intron 1, independently documented as hypervariable) are
large enough that the "outside CDS" bucket has much more room to accumulate its own diversity,
diluting the contrast. **Worth flagging as itself an interesting finding, not an artifact**: an
intron-by-intron breakdown of DRB1's outside-CDS variants (is it concentrated in one specific
intron, e.g. intron 1?) would be the natural follow-up, not yet implemented.

**Also worth noting for the next run**: at n=50 (prototype scale), HLA-C looked like a null result
(ratio 1.08, p=0.80, only 18 total variant events) -- purely an underpowered-sample artifact that
resolved into the strongest signal of the four genes (2.14x, p=1e-18) once n reached 500. Don't
trust small-sample density-comparison output at face value; this cuts against a real result, not
for it, but the direction of the correction (null -> strong signal, not the reverse) means it's not
a case of an early false-positive being walked back.

Figure: a grouped bar chart (density in vs. outside CDS per gene, annotated with ratio + p-value)
was generated from this run and sent to Marc directly (not persisted in the repo).

### Genuinely not yet done (honest gaps)

- Full ~12,000-person cohort run -- everything above is n=500, a real but partial sample.
- The PID-file lock and checkpoint-every pattern (design section above) -- not implemented, since
  runs so far have been short enough (~2-3 min at n=500) not to need them, but a full-cohort run
  should have both, especially given other concurrent VM usage confirmed live this session.
- Worker-pool parallelism -- not needed yet either; single-threaded throughput was ~10-20
  people/sec/gene, so a full 12k-person x 4-gene run is ballpark 40-80 minutes single-threaded,
  which may be acceptable as-is rather than adding thread/process-pool complexity.
- The DRB1 intron-level breakdown suggested above.
- Reconciling the raw positional track across templates (still an open design question, not
  blocking the density result above, which doesn't need it).
