# TRUST4 smoke test — person 1000291 — 2026-08-10

**Verdict: feasibility confirmed.** TRUST4 recovers a real, biologically plausible immune
repertoire from AoU whole-blood RNA-seq. This was the single question the whole repertoire
direction rested on (see `../README.md`).

## Setup

- Participant: `1000291` (first row of `manifest.tsv` — arbitrary pick, not curated)
- BAM: `gs://vwb-aou-datasets-controlled/pooled/multiomics/v9_base/rnaseq/bam/1000291.Aligned.sortedByCoord.out.md.bam`
  (paired-end, STAR v2.7.11b two-pass, mark-duplicates applied, unmapped reads kept
  in-BAM, WASP-flagged)
- TRUST4 v1.1.5-r573, references: `hg38_bcrtcr.fa` + `human_IMGT+C.fa` (both bundled in
  TRUST4's own repo, no build step needed)
- Machine: `n1-highmem-8` (8 vCPU / 52 GB RAM), Ubuntu 24.04

## The one real bug, and the fix

Default run failed immediately: `bam-extractor` couldn't pair up unmapped reads, because
AoU's STAR run used `--outSAMunmapped Within` -- unmapped reads are interleaved in-BAM
rather than laid out the way TRUST4 expects by default. **Fix: add `--abnormalUnmapFlag`**
(a boolean switch, no value) to the `run-trust4` call. Folded into
`scripts/run_trust4_sample.sh` so nobody has to rediscover this.

## Timing

| Stage | Duration |
|---|---|
| Candidate extraction (scan whole BAM) | ~8.5 min |
| Assembly + annotation + reporting | ~25 sec |
| **Total** | **~9 min** |

## Recovery

- **76,653** candidate reads survived extraction (out of 100M+ total read pairs in the
  BAM -- expected, receptor transcripts are a small slice of a whole-blood library)
- **41,395** assembled directly + **309** rescued ≈ 41,700 reads actually used
- **2,013** distinct receptor sequences (CDR3s) in the final report

## Chain breakdown

| Chain | Count | Share |
|---|---|---|
| IGK (B-cell, kappa light) | 550 | 27.3% |
| IGH (B-cell, heavy) | 500 | 24.8% |
| TRB (T-cell, beta) | 409 | 20.3% |
| TRA (T-cell, alpha) | 377 | 18.7% |
| IGL (B-cell, lambda light) | 81 | 4.0% |
| TRG (T-cell, gamma -- gamma-delta subset) | 80 | 4.0% |
| TRD (T-cell, delta -- gamma-delta subset) | 14 | 0.7% |
| unassigned | 2 | 0.1% |

B-cell chains (1,131) slightly outnumber T-cell chains (880) -- plausible: somatic
hypermutation keeps generating new distinct sequences within existing B-cell lineages,
inflating the *distinct-sequence* count even where T cells are more numerous by raw cell
count. TRG/TRD recovery (gamma-delta T cells, a rarer subset) matches a capability
TRUST4's own validation paper specifically highlighted.

## Content sanity check

Top rows by frequency were properly formatted (`TRBV7-2*01`, `TRAV1-1*02`, etc., correct
IMGT allele notation), CDR3 amino acid sequences bounded by the expected conserved
cysteine/phenylalanine anchor pattern, frequencies ranked sensibly (9.2% -> 6.0% -> 5.4%
-> 3.2% ...). No stop codons (`_`) or ambiguous bases (`?`) in the top entries.

## Not yet answered

- Whether recovery is *even* across ancestry groups (needs the ~25-person pilot, not one sample)
- Whether depth/recovery correlates with the RNA-SeQC2 QC metrics for this or other samples
- Anything about this person's actual HLA type -- no labels attached at this stage
