# Immuannot output spec (from source code, not inspection)

Source: `YingZhou001/Immuannot` @ shallow clone (scripts.pub.v3), read in full. Reference data:
`Data-2024Feb02.tar.gz` (Zenodo 10948964), downloaded and extracted. No minimap2 available in this
environment, so nothing was run empirically — every claim below is a direct read of the pipeline
code. Line numbers refer to files under
`/private/tmp/.../scratchpad/Immuannot/scripts.pub.v3/`.

Pipeline call chain per haplotype (`immuannot.sh`):
`annot.IPD.sh` (HLA/KIR, via `searchTemplate.py` → `callIPDallele.py` → `annotIPDallele.py`)
+ `annot.C4.sh` (C4, via `callC4Allele.py`) → `annot.combine.sh` (`combine.py` merges + cleans up)
→ final `<outpref>.gtf.gz`.

---

## A. Exact GTF output schema

9-column GTF/GFF2. Column 1 (seqname) = the actual **contig name** from the trimmed input FASTA
(NOT a synthetic sample id) — `annotIPDallele.py:150` `pref = [da['sample'], ...]` where
`da['sample']` is `ctgid` carried through from `searchTemplate.py`'s PAF target name. Column 2
(source) = value from `gene.names.txt.gz` column 2, i.e. `"IPD-IMGT/HLA-V3.55.0"`, `"IPD-KIR-V2.13.0"`,
or `"RefSeq"` for C4 (`annotIPDallele.py:77`, `callC4Allele.py:130`). Coordinates (cols 4/5) are
**1-based, inclusive, in the CONTIG's own coordinate system** (the trimmed FASTA's coordinates,
not hg38) — see part C.

### HLA/KIR features (`annotIPDallele.py`)

**gene** row (`annotIPDallele.py:161-181`):
- `gene_id "IAGxxxxxx"` — always. Suffixed `.N` (N=copy index, 2,3,...) when this gene has >1
  mapping cluster on this contig (`geneInfo[genename]['copyindex']`, line 151-156).
- `template_allele "HLA-A*01:01:01:01"` — always. The single best-matching reference allele used
  for gene-structure mapping (see part C for how "best" is chosen). **Not** the final typing call —
  README explicitly warns against using this field as the genotype (`reference/README_Immuannot.md:138`).
- `template_distance <int>` — always, **unquoted** (no `"..."`) unlike every other attribute here.
  Source: `ts1 = da['template_mapping'] split ';' [0] split '=' [1]` which is the `nm=` field written
  by `searchTemplate.py:333` (`attr += "nm=" + da['adjNM'] + ';'`). This is `adjNM`
  (`paftools.py:32-42`), an edit-distance-like count over the minimap2 `cs` string of the **gene-level
  alignment** (reference allele vs. contig, full gene span including introns/UTR — not CDS-only).
  It is NOT normalized; NOT the same as `cds_distance` (below). A `+`/`-` (indel) run of >5bp counts
  as 2, one of ≤5bp counts as 1, each `*` (substitution) counts as 1.
- `gene_name "HLA-A"` — always, bare gene symbol, no allele info.
- (commented out in source, `annotIPDallele.py:175-176`: `template_warning` was originally meant to
  appear on the gene row too but that line is disabled — it only appears on the transcript row.)

**transcript** row (`annotIPDallele.py:184-223`):
- `gene_id`, `transcript_id "IATxxxxxx"[.N]`, `gene_name` — same suffixing rule as gene row.
- `consensus "HLA-A*01:01:01:01"` or `..."HLA-A*01:01:new"` etc — **the field to use as the actual
  typing call** (README:134,138). Comes from `callIPDallele.py`'s `consensusCall()`/`nameNewAllele()`
  (see part B). If the gene/CDS mapping was ambiguous across genes, can be the literal string
  `"undetermined"` (`callIPDallele.py:140`).
- `alleles "HLA-A*01:01:01:01,HLA-A*01:01:01:14"` — comma-joined list of the **candidate reference
  alleles** that were tied for best CDS match (all alleles with the longest matching allele-field
  depth among those tied on NM/match, `callIPDallele.py:212-225`). This is present whether or not
  a novel allele was called. When `cds_mapping` was `"NA"` (edit distance 0 at the gene level, i.e.
  a perfect template — no separate CDS search was even run), `alleles` == `consensus` == the single
  template allele (`annotIPDallele.py:214-216`).
- `template_warning "no-start_codon,inframe_stop,..."` — **conditional**, only emitted if the raw
  `warning=` field from `searchTemplate.py`'s gene-level mapping (its own `checkCDScompleteness()`,
  `searchTemplate.py:60-77`) is non-empty. Possible tokens: `partial_CDS`, `inframe_stop`,
  `no-start_codon`, `no-stop_codon`. This is warning info about the **template's CDS reconstruction
  from the gene-level minimap2 alignment**, separate from `cds_mut`/novelty below.
- `cds_distance <int>` — **conditional**, only present if a CDS-level search was actually run (i.e.
  `cds_mapping != "NA"`, meaning the gene-level template wasn't already a perfect 0-distance match —
  `annotIPDallele.py:202-213`). This is the `NM` (raw minimap2 edit count, unmodified `adjNM`-style)
  of the CDS-vs-CDS alignment (`callIPDallele.py:229`, field `'NM'`), over the **CDS-only span**
  (concatenated exonic sequence, not full gene) — see part C.
- `cds_mut "HLA-A*01:01:01:01|:301*ag:...|K(AAG)<E(GAG)"` — **conditional**, only present when
  `cds_distance > 0` (`annotIPDallele.py:212-213`). Pipe-delimited: `[best-matching reference allele]
  | [raw minimap2 cs string of the CDS alignment] | [amino-acid-level diff string from
  algntools.findCodingDiff()]`. If multiple candidate alleles tied, one `mut` entry per candidate,
  comma-joined (`callIPDallele.py:230`, `sel_muts`). The diff string format is
  `REFaa(REFcodon)<OBSaa(OBScodon)` per changed codon, colon-joined for multiple changed codons
  (`algntools.py:290-294`).

**UTR** row: `gene_id`/`transcript_id`/`gene_name` only, no extra attrs (`annotIPDallele.py:226-243`).
Rows with negative coordinates (unmapped/out-of-range) are silently dropped, not emitted.

**exon** row: adds `exon_number "N"` (integer index into the reference allele's own exon numbering,
`annotIPDallele.py:245-263`) — this numbering follows the REFERENCE allele's exon count/order, not
a re-derived one; for a `-` strand gene it still uses the reference's 5'→3' numbering.

**start_codon** / **stop_codon** rows: add `codon "ATG"` / `codon "TGA"` (the literal 3 bases from
the reference record, not re-read from the contig) — a row is only emitted if the codon actually
equals a valid start (`ATG`) or stop (`TGA`/`TAG`/`TAA`) codon after coordinate mapping
(`annotIPDallele.py:265-307`); a broken/absent start or stop codon is simply **omitted**, not
flagged with an explicit "missing" attribute (though `template_warning` on the transcript row would
usually already say `no-start_codon`/`no-stop_codon`).

**CDS** rows: `gene_id`/`transcript_id`/`gene_name` only. Column 8 (frame/phase) is computed
properly per-segment (`annotIPDallele.py:342-359`), correctly handling `+`/`-` strand ordering. The
stop codon's 3 bases are deliberately excluded from the terminal CDS segment
(`removeStopCodonFromCDS = True`, `annotIPDallele.py:10,321-326`) — a script re-deriving protein
sequence from these coordinates must NOT append the stop codon (there isn't one to append; it's
already excluded, matching Ensembl-style GTF convention).

### C4 features (`callC4Allele.py`)

Same feature types (gene/transcript/exon/CDS/start_codon/stop_codon/UTR) but a **different,
narrower attribute set** — no `template_allele`, `template_distance`, `consensus`, `alleles`,
`cds_mut`, or novelty tagging of any kind. C4 typing is not reference-allele-alignment based; see
part B.
- `gene_name` value is **not a clean gene symbol** — it is the concatenation of the called allele
  (`C4`, `C4A`, or `C4B`) with a single size-code letter appended with **no delimiter**:
  `C4AL` (long, intron 9 > 5000bp) or `C4AS` (short) etc (`callC4Allele.py:8-14,163,235`).
  `combine.py:70-71` strips the last character back off (`genename[0:(len(genename)-1)]`) to
  recover the bare gene symbol for its own internal gene-list bookkeeping — a downstream parser
  reading gene_name literally must replicate that same strip, or match on `startswith("C4")`.
- transcript row gets `pipetide_key "..."` (sic — misspelled "pipetide" for "peptide" in the source,
  `callC4Allele.py:221`) — the translated key exon-26 hexamer used to distinguish the A/B isotype,
  not an allele-precision field.
- gene_id gets a bare `.N` copy-index suffix with **no quote-then-semicolon on the un-suffixed
  case in the same style** — verified identical `'";'` vs `'.'+str(copyindex)+'";'` pattern as HLA
  (`callC4Allele.py:229-234`), so parsing logic can be shared.

### Header lines (`combine.py:101-113`)
Every final `<outpref>.gtf.gz` begins with `#`-prefixed lines (not part of the 9-column body):
`## format: gtf`, `## date: YYYY-MM-DD`, then up to three lines listing genes by observed copy
number — `## gene (copy num = 0): ...` (genes with zero hits — **note: this line lists genes not
found on THIS FILE'S contig(s)**, useful as a fast "was this gene even attempted" check without
re-scanning the whole body), `## gene (copy num = 1): ...`, `## gene (copy num > 1): ...` — then
one `## contigs for <gene>: <ctg1>,<ctg2>` line per gene that has ≥1 copy, listing which contig(s)
it was found on. **This header is cheap, already-computed ground truth for "is gene X even present
in this hap's GTF" without grepping the body** — worth using directly in any consumer script.
Genes in `<refdir>/gene.black.list` are excluded entirely from all three copy-number lines and
from the body (though in the actual reference data downloaded, `gene.black.list` is **empty** — no
genes are excluded by default).

---

## B. Novelty representation

Novel-allele naming happens in `callIPDallele.py`'s `nameNewAllele()` → `nameNewHlaAllele()`
(lines 46-81) or `nameNewKirAllele()` (84-117) (KIR is moot for this project's chr6-only trim, see
part F).

**Mechanism**: starting from the best-matching reference allele string (e.g. `HLA-A*01:01:01:01`),
the code inspects the `cds_mut` diff string (the `algntools.findCodingDiff()` amino-acid diff
described in part A) and truncates the allele's colon-delimited field list, replacing the deepest
retained field with the literal string `"new"`:
- No CDS diff computed at all (gene-level perfect match, `cdsdiff` is falsy/`"NA"`) → 4th field
  becomes `"new"`: e.g. `HLA-A*01:01:01:new` (synonymous-at-genomic-level / untyped-beyond-CDS case).
- CDS diff found and **synonymous** (translated AA strings equal on both sides) → 3rd field becomes
  `"new"`, first 3 fields kept: e.g. `HLA-A*01:01:new`.
- CDS diff found and **non-synonymous**, no frameshift → 2nd field becomes `"new"`, only first 2
  fields kept: e.g. `HLA-A*01:new`.
- **Frameshift** detected (`len(observed_codon_string) != len(reference_codon_string)` for any
  changed segment) → same truncation as non-synonymous (2nd field → `"new"`), regardless of whether
  the AA translation happened to look synonymous.

**So: "new" is embedded INSIDE the allele string itself**, at a *variable* colon-delimited field
depth (2nd, 3rd, or 4th field) that itself encodes coarse novelty severity (protein-changing vs.
synonymous vs. untyped-beyond-CDS). It is not a separate boolean/flag attribute. `consensusCall()`
(lines 127-159) then further collapses this across multiple tied candidate alleles: if the tied
candidates disagree at the field where "new" would be inserted, it takes `os.path.commonprefix()`
and truncates to `"new"` at THAT shorter common depth instead — so the final field depth at which
"new" appears is a genuine (if coarse) signal of confidence/typing depth, not arbitrary.

**Does the output record WHAT the difference is?** Yes, but only for genes that needed a CDS-level
search at all (`cds_distance > 0`): the `cds_mut` transcript-row attribute (part A) gives the exact
minimap2 `cs` alignment string PLUS a derived amino-acid diff string
(`REFaa(REFcodon)<OBSaa(OBScodon)`, per changed codon). This is **CDS-level only** — intronic/UTR
differences are never surfaced as a "what changed" string anywhere in the GTF (only as the coarse
integer `template_distance` on the gene row, which mixes intron+UTR+CDS differences together).
There is **no per-exon or per-CDS-segment breakdown of novelty** — `cds_mut` is one string for the
whole CDS, with per-codon diffs concatenated by `:`; you cannot tell from the GTF alone which exon
number a given diff codon falls in without independently re-mapping the CDS-relative codon
position back to an exon using the `exon "N:a..b"` coordinates already in the GTF (straightforward:
the CDS segments in `cds` transcript coordinates are known from `CDS=...` in `alleles.csv.gz` for
the matched reference allele, so a position in the concatenated CDS can be mapped to an exon offset
arithmetically — not automatic, but derivable).

For **C4**, there is no novelty concept at all — `callC4Allele.py` only ever emits `C4`, `C4A`, or
`C4B` (plus size letter) via a fixed regex/heuristic on two short marker sequences (exon 26 key
hexamer, intron 9 length); there is no reference-allele database comparison, no distance, no "new"
tag possible.

---

## C. Distance metric semantics

Two distinct, differently-scoped distance numbers exist, both integers, **neither normalized**:

1. **`template_distance`** (gene row) = `adjNM` (`paftools.py:32-42`) of the **gene-level** minimap2
   alignment (`minimap2 -cx asm5 --cs --end-bonus=10`, contig as target, full IPD/IMGT gene sequence
   — including introns and UTR — as query, `annot.IPD.sh:34-36`). It is a **weighted count of
   variant events** along the `cs` string, not a classic Levenshtein edit distance and not a
   percent identity: each substitution (`*`) = 1, each indel (`+`/`-`) run ≤5bp = 1, each indel run
   >5bp = 2. Coordinates for this alignment ARE contig-relative (see below) — the aligned span's
   length in contig bases is directly recoverable: it's `trg` = `tto - tfro` from the same
   `template_mapping` attribute string that carries `nm=` (`searchTemplate.py:337`,
   `trg=<tfro>..<tto>` are the 0-based minimap2 target-block start/end used to build the `gene`
   feature row's own start/end, `annotIPDallele.py:24-27,164-166`) — i.e. **the GTF `gene` row's own
   (start,end) IS that aligned span**, so `template_distance / (gene.end - gene.start + 1)` is a
   valid, already-derivable per-base "distance rate" without needing anything not in the GTF. It is
   filtered upstream at the gene-detection stage to `adjNM <= 3 OR adjNM/match < diffcut(0.03)`
   (`searchTemplate.py:375-377`) — so genes failing both cutoffs never make it into the GTF at all
   (silently absent, not flagged).

2. **`cds_distance`** (transcript row, conditional) = raw `NM` (tag straight from minimap2's PAF
   `NM:i:` column, `paftools.py:20-28`, used unmodified — **not** run through `adjNM`) of the
   **CDS-only** alignment: `minimap2 -c --cs --end-bonus=10` of the assembled/reconstructed CDS
   nucleotide sequence (see part D — `cds.fa.gz`, built from the gene-level alignment) against the
   reference alleles' own extracted CDS sequences (`<refdir>/CDSseq/<gene>.fa.gz`,
   `annot.IPD.sh:56-65`). This is a true edit-distance-like count over CDS bases only (no introns/
   UTR). Its aligned span length is recoverable too, though less directly from the final GTF alone:
   sum of `(to - from + 1)` over every `CDS` feature row for that gene/transcript (adding back the
   3 stop-codon bases that were deliberately excluded, part A) gives the CDS length actually used.

Both are computed on **contig-relative 1-based coordinates** throughout the final GTF (confirmed:
`geneStructMapping()`/`intervalMapQry2Ctg()` in `annotIPDallele.py`/`algntools.py` convert minimap2's
0-based target block coordinates to 1-based only at the point of GTF-row emission, e.g.
`annotIPDallele.py:164-166` `line.append(str(int(fro) + 1))`). **Not hg38 coordinates** — column 1
is the trimmed contig's own name/coordinate space, so cross-referencing to hg38 would require
re-applying whatever alignment (the person's own `assembly_hap{N}_aln2_hg38_paf`/`.bam`, already
used by `run_immuannot_person.py` for trimming) maps that specific contig back to chr6.

---

## D. Sequence recoverability — can we get the actual observed nucleotide sequence?

**Yes, and there is a MUCH cheaper path than reading coordinates back out of the trimmed FASTA: the
per-gene-copy observed CDS sequence is already written out by Immuannot itself, in coding (sense)
orientation, ready to use, in a file that Immuannot's own cleanup step does NOT delete.**

### The direct route: `<outpref>/cds.fa.gz` (survives cleanup)

`searchTemplate.py` (the very first stage of `annot.IPD.sh`) reconstructs, for **every** gene copy
it detects on the contig (not just novel ones — this runs before novelty is even determined), the
actual contig-derived CDS nucleotide sequence, by walking the *gene-level* minimap2 `cs` string with
`algntools.recoverTargetSeqFromQuery()` (searchTemplate.py:270-323). It writes these to
`<outpref>/cds.fa.gz` (`outfile2`, opened at `searchTemplate.py:16,21`) with header
`>{ctgname}_{gene}_{i}` where `i` is a 1-based sequential index over ALL gene-copy records found on
that specific contig (across every gene, not per-gene) — plus intron-boundary dinucleotides as a
free-text suffix on the header line. `recoverTargetSeqFromQuery()`'s own docstring
(`algntools.py:129-131`) states it "always return[s] the strand the same as qseq" — i.e. **the CDS
sequence is emitted in coding/5'→3' orientation regardless of the gene's genomic strand**, so no
separate reverse-complementation step is needed to get a translatable CDS.

This file is explicitly **kept**: `annot.combine.sh`'s cleanup block (lines 16-23) removes
`tmp.ipd.gtf.gz`, `tmp.c4.gtf.gz`, `tmp.newgene.txt`, `tmp.allele.csv`, `tmp.gene.csv`, and
`tmp.ipd.cds.fa.gz` — but **not** `cds.fa.gz` (a differently-named, earlier-stage file — easy to
mis-read as "the same file" from the name alone; it is not). Since
`reference/README_Immuannot.md`/`RESULTS_LOCATION.md` confirm this project keeps the whole
per-haplotype output folder (`hap{1,2}/`) uncleaned ("no per-person pruning", DECISIONS.md
2026-08-04), `cds.fa.gz` should exist on disk for every haplotype that was ever run, going all the
way back to the pilot.

**Caveat — matching a `cds.fa.gz` record to a specific final-GTF gene row is unambiguous only for
single-copy genes.** The header index `i` is assigned in the order `searchTemplatePGPC()` discovers
mapping clusters for that contig (across all genes together), matching the row order of the now-
deleted `tmp.gene.csv`. But `combine.py`'s final merge (lines 81-82) **re-sorts** every GTF row by
`(contig name, then genomic start coordinate)` before writing the final `.gtf.gz` — a different
order than detection order. For a gene with exactly one copy on that contig (the overwhelming
majority case — one classical HLA gene per haplotype-contig), this is a non-issue: match by
`(ctgname, gene_name)` alone, since there is only one candidate on each side. For a gene with
copy number > 1 on the same contig (segmental duplication case, part G), **the join is genuinely
ambiguous from `cds.fa.gz` + final GTF alone** — see "WHAT WE CANNOT GET" below for the concrete
resolution.

### The indirect route: reading the trimmed FASTA directly by GTF coordinates

Also possible, and useful as a cross-check or for getting genomic (not just CDS) sequence (e.g. full
gene span with introns, or UTRs): `hap{1,2}.trimmed.fa` is confirmed kept
(`RESULTS_LOCATION.md:14-15`). Procedure:
1. Column 1 of the GTF row = the exact FASTA record ID in `hap{N}.trimmed.fa` to `samtools faidx`
   (these are real contig names carried through unmodified from the assembly, not resynthesized).
2. Columns 4/5 are **1-based, inclusive**, contig-relative (part C) — directly usable as a
   `samtools faidx contig:start-end` region string with no off-by-one adjustment.
3. Column 7 (strand) — if `-`, reverse-complement the extracted sequence to get coding-sense
   orientation; the GTF's own CDS/exon row coordinates are still given as ascending genomic
   (forward-strand) ranges regardless of gene strand (standard GTF convention, confirmed by the
   sort step in `annotIPDallele.py:366-367` sorting `cdsdf`/`bodyarr` by ascending column 4 in both
   strand cases).
4. To reconstruct a spliced CDS: concatenate all `CDS` rows for the same `(contig, gene_id)` in
   **increasing genomic order for `+` strand, decreasing genomic order for `-` strand** (this
   matches how `annotIPDallele.py:343-346` orders `cdsdf` when computing frame) — then reverse-
   complement the whole concatenation once for `-` strand genes (not per-exon).
5. The stop codon is NOT included in the CDS rows (part A) — append it back manually from the
   `stop_codon` row's own coordinates if a complete ORF including the stop is wanted.

Both routes should agree exactly for genes with `template_distance == 0` and no CDS-level warnings;
for genes with any warning (`incomplete_mapping`, part A/C), the FASTA route is the more trustworthy
one since it doesn't depend on `cs`-string reconstruction succeeding cleanly. **Recommend cross-
checking both routes against each other on a handful of real genes on the VM before trusting either
at scale** — this specific claim (that the two routes agree) is inferred from source logic, not
empirically verified (no minimap2 available in this environment to run the pipeline end-to-end).

---

## E. Intermediate files — what's actually on disk in `<outpref>/`

`<outpref>` for this project is `hap{1,2}` (a folder), sibling to the final `hap{1,2}.gtf.gz`.
Files written by the pipeline into that folder, and their fate after `annot.combine.sh`'s cleanup:

| File | Written by | Content | Deleted by cleanup? |
|---|---|---|---|
| `mm2.ipd.gen.paf.gz` | `annot.IPD.sh:36` | **Raw, unfiltered** minimap2 PAF: every IPD/IMGT+KIR reference gene sequence aligned against the trimmed contig, `asm5 --cs`. This is the full candidate universe before ANY of `searchTemplate.py`'s quality filters are applied. | **No** |
| `gene.filtered.paf` | `searchTemplate.py:18,23,388` (`file4`) | The subset of the above PAF rows that **failed** the gene-detection filters (query≤target length, match/qlen>overlapcut, mapped-length ratio, adjNM≤3 or adjNM/match<diffcut) — i.e. candidate alleles/genes that were considered and rejected. | **No** |
| `cds.fa.gz` | `searchTemplate.py:16,21` (`file2`) | Per-gene-copy reconstructed observed CDS sequence, coding orientation (part D). **The richest file for sequence-level novel-allele work.** | **No** |
| `tmp.newgene.txt` | `searchTemplate.py:17,22` (`file3`) | Bare list of gene symbols that had ≥1 detected copy on this contig (used downstream to know which `CDSseq/<gene>.fa.gz` reference files to load for the CDS-level search). Trivially reconstructible from the final GTF's own header `## gene (copy num...)` lines. | **Yes** (low loss — reconstructible) |
| `tmp.gene.csv` | `searchTemplate.py:15,20` (`file1`) — this is `annot.IPD.sh`'s `calltmpfile` | Per-cluster **winning-template** call table: `ctgname, gene, strand, template, allele, attr` — one row per detected gene-copy, in original detection order (the join key `cds.fa.gz` needs for multi-copy genes, see part D/G). | **Yes** — this is the one genuinely useful, non-trivially-reconstructible loss (see "WHAT WE CANNOT GET") |
| `mm2.ipd.cds.paf.gz` | `annot.IPD.sh:65` | minimap2 alignment of **every** reference allele's CDS (`CDSseq/<gene>.fa.gz`, the full per-gene allele database, not just top candidates) against `cds.fa.gz`. **This is the full ranked candidate list with NM/match per candidate** — i.e. exactly the file needed to answer "what was the 2nd/3rd-best candidate allele and its distance," which the final GTF and `tmp.allele.csv` (deleted) both collapse away to "the tied best set only." | **No** |
| `tmp.ipd.cds.fa.gz` | `annot.IPD.sh:50,55,62` (`qcdsseq`) | Concatenation of the **reference** CDS sequences (from `<refdir>/CDSseq/`) for every gene in `tmp.newgene.txt` — pure reference data, zero information loss (100% reconstructible from `refdir` + `tmp.newgene.txt`/GTF header). | **Yes** (zero loss) |
| `tmp.allele.csv` | `annot.IPD.sh:68,70` (`calltmpfile2`, `callIPDallele.py`'s stdout) | Per-gene-copy `consensus`/`nm`/`mut` call, i.e. the same info that ends up embedded in the final GTF's transcript-row `consensus`/`cds_distance`/`cds_mut` attributes. | **Yes** (near-zero loss vs. final GTF, EXCEPT it retains this info even for genes where `annotIPDallele.py` might drop/warn — not separately verified, treat as low-risk) |
| `tmp.ipd.gtf.gz` | `annot.IPD.sh:74-75` | Pre-merge HLA/KIR-only GTF (before C4 is merged in and before the final coordinate-sort). Strict subset of the final GTF's HLA/KIR rows, just unsorted. | **Yes** (zero loss) |
| `tmp.c4.gtf.gz` | `annot.C4.sh:28-30` | Pre-merge C4-only GTF. Same as above for C4. | **Yes** (zero loss) |
| `mm2.c4.exon.paf.gz` | `annot.C4.sh:23-25` | minimap2 split-alignment (`-C5 --cs -cx splice:hq`) of the C4 exon-marker query file against the contig — the raw alignment C4 typing is derived from. | **No** |

**Net effect for this project's pipeline (which runs the unmodified upstream `immuannot.sh`, per
`run_immuannot_person.py:270` calling it with no `--outpref`-folder-skipping flag, and separately
never deletes the outpref folder itself, `RESULTS_LOCATION.md`):** every haplotype's `hap{1,2}/`
directory should currently contain `mm2.ipd.gen.paf.gz`, `gene.filtered.paf`, `cds.fa.gz`,
`mm2.ipd.cds.paf.gz`, and `mm2.c4.exon.paf.gz` — i.e. **"what was the 2nd-best candidate allele and
its distance" IS answerable without re-running anything**, by re-parsing `mm2.ipd.cds.paf.gz`
(and, for gene-level rather than CDS-level ranking, `mm2.ipd.gen.paf.gz`) the same way
`callIPDallele.py`'s `callNewAllele()` does, just without discarding everything past the tied-best
set. This is a genuinely valuable, currently-untapped resource sitting on the VM already.

---

## F. Full gene list actually emitted (from `Data-2024Feb02/gene.names.txt.gz`, 65 genes total)

**Classical HLA class I:** HLA-A, HLA-B, HLA-C

**Classical HLA class II:** HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB, HLA-DPA1, HLA-DPA2, HLA-DPB1,
HLA-DPB2, HLA-DQA1, HLA-DQA2, HLA-DQB1, HLA-DQB2, HLA-DRA, HLA-DRB1, HLA-DRB2, HLA-DRB3, HLA-DRB4,
HLA-DRB5, HLA-DRB6, HLA-DRB7, HLA-DRB8, HLA-DRB9

**Non-classical HLA / pseudogenes:** HLA-E, HLA-F, HLA-G, HLA-H, HLA-J, HLA-K, HLA-L, HLA-N, HLA-P,
HLA-S, HLA-T, HLA-U, HLA-V, HLA-W, HLA-Y, HLA-HFE (note: HLA-HFE is the hemochromatosis gene, not a
true HLA locus, but bundled in the same reference set)

**MIC / TAP:** MICA, MICB, TAP1, TAP2

**Complement C4:** C4A, C4B, C4X (C4X appears in `gene.names.txt.gz` but the actual calling logic
in `callC4Allele.py` only ever assigns `C4`, `C4A`, or `C4B` — "C4X" is present as a gene-id
placeholder in the reference table but is not a real output category the pipeline emits)

**KIR (chr19):** KIR2DL1, KIR2DL2, KIR2DL3, KIR2DL4, KIR2DL5A, KIR2DL5B, KIR2DP1, KIR2DS1, KIR2DS2,
KIR2DS3, KIR2DS4, KIR2DS5, KIR3DL1, KIR3DL2, KIR3DL3, KIR3DP1, KIR3DS1

**Confirmed absent given our chr6:29.5–33.5Mb trim window:** all 17 KIR genes are real chr19 loci
(README.md's own detection-strategy figure references the KIR locus's known chr19 organization).
Since the trimmed input FASTA (`run_immuannot_person.py`'s trim logic, region default
`chr6:29500000-33500000`) never contains chr19 sequence, minimap2 will produce zero PAF hits for
any KIR reference gene against the trimmed contig — KIR genes will not appear in the "copy num=1"
or "copy num>1" header lines, only ever in the "copy num=0" header line (or simply absent from that
line too if `gene.black.list` were used to exclude them, which it currently is not — the shipped
blacklist file is empty). This is a real absence from lack of matching sequence, not a
`gene.black.list` exclusion.

---

## G. Per-haplotype cross-gene linkage (cis-phasing question)

**Confirmed: genes ARE only validly "on the same physical haplotype" if they share the same
column-1 contig name within that hap's GTF — being in the same `hap{N}.gtf.gz` FILE is necessary
but not sufficient.**

`run_immuannot_person.py`'s trim step (`regions_from_paf`/`trim_assembly`,
`run_immuannot_person.py:160-254`) can produce a trimmed FASTA containing **more than one contig**
when the phased assembly is fragmented across the target region (multiple distinct contigs each
overlapping part of chr6:29.5–33.5Mb) — this is explicitly tracked (`row["n_contigs"]`) and is not
a corner case the code treats as impossible. `immuannot.sh` accepts that multi-contig FASTA as a
single `-c` input and processes all contigs in one run, so a single `hap1.gtf.gz` CAN legitimately
contain genes from two (or more) different, non-contiguous physical contigs. Column 1 (seqname) is
the only field distinguishing them — `run_immuannot_person.py`'s own `parse_gtf()`/
`rebuild_immuannot_calls.py`'s `parse_gtf()` currently **discard column 1 entirely**, keying calls
only by `(person_id, hap, gene)` — so as currently written, our pipeline's aggregate
`immuannot_calls.tsv` cannot itself distinguish "two genes on the truly same contig" from "two genes
on two different contigs that both happen to overlap the trim window and both ended up in hap1's
GTF." **For cis-heterodimer reconstruction (e.g. DQA1~DQB1 pairing), the safe join key is
`(person_id, hap, contig_name)`, read fresh from the still-on-disk raw `hap{N}.gtf.gz` column 1 —
not `(person_id, hap)` alone.** This is a real, currently-live gap, not just a theoretical one.

**Duplicated gene clusters on ONE contig:** `searchTemplatePGPC()`'s cluster-merging logic
(`searchTemplate.py:157-181`, `mergeInterval()`) explicitly supports and expects multiple genuinely
separate mapping clusters for the *same* reference gene on one contig (segmental duplication /
true extra copy) — reported via the `gene_id`/`transcript_id` `.N` suffix (part A), consistent with
README's "copy number ... reported naturally with the number of mapping clusters" (README:51-53).
This is real biological copy-number signal, not a pipeline artifact, but as noted in part D, it
does introduce a genuine ambiguity when trying to join a specific `cds.fa.gz` record (or an
`mm2.ipd.cds.paf.gz` candidate-ranking record) back to a specific numbered copy in the sorted final
GTF, because the row-order key that would disambiguate them (`tmp.gene.csv`) is deleted by cleanup.

---

## Section 3 — What our own pipeline captured vs. discarded

- **`run_immuannot_person.py`'s `parse_gtf()`** (lines 296-316): scans every GTF line for
  `gene_name`/`gene_id` AND `consensus`/`allele` co-occurring on the same line — only the
  **transcript** row satisfies both (gene row has `gene_name` but no `consensus`; the `alleles`
  attribute doesn't match the `allele "..."` regex because of the trailing `s` before the quote).
  Result: captures exactly one string per `(hap, gene)` — the `consensus` value — and **discards
  everything else**: `template_allele`, `template_distance`, `template_warning`, `cds_distance`,
  `cds_mut`, `alleles` (the candidate list), all coordinates, all C4-specific fields.
- **`scripts/production_orchestrator/rebuild_immuannot_calls.py`'s `parse_gtf()`** (lines 62-97,
  written later, 2026-08-10) additionally captures `template_distance` and a boolean
  "had-any-template_warning" flag per gene, cached to `immuannot_confidence.tsv`. Still discards
  `cds_distance`, `cds_mut`, `alleles`, and contig identity (column 1) — see part G's gap above.
- **`scripts/diagnose_immuannot_pilot.py`'s `parse_gtf_rich()`** (lines 228-253) is the most
  complete of the three: it merges **every** attribute key seen across every row for a gene
  (`slot.setdefault`, so first-seen value per key wins, not last) — this WOULD capture
  `cds_mut`/`alleles`/`template_warning`/etc. if invoked, but its own docstring flags itself as
  "UNVERIFIED against a real transcript-row attribute string as of 2026-07-22 (only gene-row
  samples have actually been seen printed from a real run)" — i.e. this project has apparently
  never actually confirmed a real transcript-row string end-to-end through this specific parser.
  Worth a real spot-check (see AMBIGUOUS section below) before trusting `cds_mut`/`alleles`
  extraction from it at scale.
- **Intermediate folder**: confirmed NOT pruned by anything in this project's own code
  (`run_immuannot_person.py` never touches the `hap{N}/` outpref folder after `run_immuannot()`
  returns; `RESULTS_LOCATION.md` explicitly documents `hap{1,2}.trimmed.fa` as "kept, not pruned"
  and lists nothing about deleting the folder). Combined with Immuannot's own internal cleanup
  (part E table), the surviving intermediate files per haplotype should be: `mm2.ipd.gen.paf.gz`,
  `gene.filtered.paf`, `cds.fa.gz`, `mm2.ipd.cds.paf.gz`, `mm2.c4.exon.paf.gz` — **AMBIGUOUS: this
  has not been empirically confirmed on the actual VM for this project's real runs** (only inferred
  from source code + this project's own documentation of its own scripts). One-line check to
  resolve: `ls -la ~/pipeline_outputs/<any_completed_person_id>/immuannot_output/hap1/` on the VM.

---

## WHAT WE CANNOT GET WITHOUT RE-RUNNING

1. **The exact per-cluster detection-order join key (`tmp.gene.csv`) linking a specific `cds.fa.gz`
   record / `mm2.ipd.cds.paf.gz` candidate-ranking record to a specific numbered gene copy, for
   genes with copy number >1 on the same contig.** Deleted by `annot.combine.sh`. Only matters for
   multi-copy genes (rare for the 3 classical class I + handful of class II genes this project
   likely cares about most; more relevant if DRB region copy-number variation matters).
   **Cheapest partial re-run**: `hap{1,2}.trimmed.fa` is kept, so re-run just
   `bash annot.IPD.sh <trimmed.fa> <refdir> <fresh_outpref> <threads> 0.9 0.03` (the first stage
   only, skipping C4/`combine.sh`) on the affected person/hap — this regenerates `tmp.gene.csv`
   fresh without needing to touch the read-only gcsfuse mount or re-derive the trim at all. Minutes,
   not hours.

2. **`template_warning`/CDS-completeness detail for genes whose gene-level template was already a
   perfect (0-distance) match** — `cds_mapping` is never even computed in that case
   (`annotIPDallele.py:202` branch skipped), so there is no `cds_distance`/`cds_mut` for those genes
   by design (not a data-loss problem — Immuannot itself never computed a CDS-level alignment
   because it didn't need to). Nothing to recover; this is not missing data, it's "not applicable."

3. **Column-1 contig identity per gene, if the aggregate `immuannot_calls.tsv`/
   `immuannot_confidence.tsv` are the only files consulted.** Not actually unavailable — the raw
   `hap{1,2}.gtf.gz` files still have it in column 1 and are kept — but every current parsing script
   in this repo throws it away. **No re-run needed at all**, just a new parse of files already on
   disk (see the cis-linkage caveat, part G).

4. **Confirmation that `cds.fa.gz` sequences exactly match what `samtools faidx`-extracting the
   trimmed FASTA by GTF coordinates would give** (part D's "recommend cross-checking" note) —
   AMBIGUOUS, inferred from source logic only, not run end-to-end in this environment (no minimap2
   available here). **One-line check on the VM**, no re-run: pick one gene/person/hap, extract its
   CDS both ways, `diff` them.
   ```
   zcat hap1.gtf.gz | awk -F'\t' '$3=="CDS" && /gene_name "HLA-A"/' | sort -k4,4n   # coordinates
   zcat hap1/cds.fa.gz | grep -A1 "^>.*_HLA-A_"                                      # cds.fa.gz seq
   ```

5. **Whether `mm2.ipd.gen.paf.gz`/`gene.filtered.paf`/`cds.fa.gz`/`mm2.ipd.cds.paf.gz`/
   `mm2.c4.exon.paf.gz` genuinely still exist on the production VM's disk for already-completed
   people** (Section 3's AMBIGUOUS item) — no re-run needed if they do; if disk space pressure ever
   caused a manual `rm` of these that isn't reflected in this repo's scripts, this would need
   confirming with `ls`. One-line check: see above.

6. **Anything about genes that failed detection entirely** (didn't pass `searchTemplate.py`'s gene-
   level filters, `searchTemplate.py:369-377`) beyond knowing THAT they failed (visible via the
   GTF's own "copy num = 0" header line, or via records in `gene.filtered.paf` which do survive).
   `gene.filtered.paf` gives the raw failed alignment (NM, match length, coordinates) so a genuinely
   marginal near-miss can be inspected without re-running anything — but if a gene has ZERO PAF
   records at all against it (true absence of any detectable homology), there is nothing to recover
   short of re-running minimap2 with relaxed parameters (`--overlaprate`/`--diff` flags on
   `immuannot.sh`) against the same kept `hap{N}.trimmed.fa` — again a cheap partial re-run, no
   mount access needed.
