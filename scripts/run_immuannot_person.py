#!/usr/bin/env python3
"""Run Immuannot (github.com/YingZhou001/Immuannot) on one or more AoU people.

Immuannot types HLA alleles from a PHASED ASSEMBLY contig (one haplotype's collapsed consensus
FASTA), not reads/BAM -- see DECISIONS.md "Assembly-based HLA typing". Only people on the
revio/sequel2e/sequel2 platforms have this data (reports/pre_sprint/lr_data_census/README.md); ont-r10.4.1/
ont-r9.4.1 do not and will be skipped with a clear reason, not silently dropped.

TRIM STEP (2026-07-21, Marc: "trim the zone we need, not ridiculous paddings" before the 60-person
pilot run). A whole diploid assembly (both haplotypes, genome-wide) is far more than Immuannot
needs to type 8 classical HLA genes on chr6 -- but assembly contigs are NOT reference-coordinate
strings (they're de novo sequence, arbitrary contig IDs), so we can't naively slice
"chr6:29,500,000-33,500,000" out of them the way slice_and_fastq.sh does for aligned BAMs.

Three-tier design, preferring the tighter one (Tier 3 added 2026-08-04, UNTESTED AT SCALE, opt-in
via --enable-self-align-fallback -- see the module-level comment near MINIMAP2_ASM_PRESET below
for the sequel2 motivation and Marc's ~700-750x-blowup concern that shaped this design):
  1. PREFERRED -- true sub-range extraction via each haplotype's own assembly-to-hg38 .paf file
     (assembly_hap{1,2}_aln2_hg38_paf). A .paf line's qstart/qend/tstart/tend are the REAL aligned
     boundaries of that block, not an estimate, so regions_from_paf() takes the union of
     (min qstart, max qend) across every block overlapping the target region, per contig, padded
     by --pad (default 100,000bp -- reuses the exact number EXPERIMENTS.md's SpecImmune-LR sweep
     already proved safe: zero degradation, 3.9x smaller than full-pad) and clamped to the contig's
     own length. Confirmed 2026-07-21 on person 1001871: cut the trimmed input further than
     whole-contig extraction alone -- see immuannot_timing.tsv's whole_contig_mb/padded_mb columns
     for the actual per-person reduction once run.
  2. FALLBACK -- whole-contig extraction via the assembly-to-hg38 .bam (assembly_hap{1,2}_aln2_hg38_bam),
     used only if this person's .paf is missing or has nothing recognizable overlapping the region.
     No CIGAR/indel math, no risk of clipping a gene -- just less tight than tier 1.
  3. TIER 3 / SELF-ALIGN, opt-in via --enable-self-align-fallback -- used only when BOTH the .paf
     AND the aln-to-hg38 .bam are missing (sequel2's 991 people: raw assembly FASTA present, no
     alignment-to-hg38 files at all, per reports/pre_sprint/lr_data_census/README.md). Self-aligns the hap
     FASTA against a chr6-only reference slice with minimap2 (own synthetic .paf, parsed by the
     same regions_from_paf() as tier 1) instead of feeding Immuannot the whole ~3.1Gb/hap untrimmed
     assembly. UNTESTED AT SCALE as of 2026-08-04 -- smoke-test on real sequel2 people before
     trusting it for a production run.
     2026-08-05: FAILED on 2 real sequel2 people (EXPERIMENTS.md). 2026-10-11: rewritten as a
     two-stage self-align (region-slice prefilter, then whole-primary-GRCh38 confirmation via
     --self-align-confirm-ref) with a streaming extractor that reads plain gzip -- see the
     "Tier 3 v2" comment block below and sprints/S05_qc_phasing_kir_trb/SEQUEL2_DIAGNOSIS.md.
     Still opt-in, still to be smoke-tested on real sequel2 people (63_sequel2_recall.py --smoke).
Tiers 1+2 confirmed present alongside the raw assembly FASTA for revio/sequel2e per the census.

Per person: resolves assembly_hap1_fa/assembly_hap2_fa AND assembly_hap1_aln2_hg38_bam/
assembly_hap2_aln2_hg38_bam mount-relative paths from the v9 lrWGS manifest (existence-checked
directly against the mount -- same discipline as build_experiment_d_cohort.py's
find_existing_bam(), never a path-pattern guess -- Marc's "I'm sure they're all revio" is treated
as a hypothesis to verify, not a fact, per this project's own repeated lesson about that). Then,
per haplotype: finds overlapping contig(s), trims, runs immuannot.sh (the tool's own documented
mode -- it does not accept both haplotypes in one run), verifies the expected *.gtf.gz output
actually exists (never trust exit code alone -- ENVIRONMENT.md quirks #17/#18), times every stage,
and writes a per-person, per-gene 2-allele call table plus a timing log for the scale-up decision.

GTF-PARSING NOTE: verified against the real reference/README_Immuannot.md now in this repo (not
just a web summary) -- "consensus"/"allele" attribute (not "template_allele") is correct, and it
lives on the GTF's "transcript" feature row, which also carries the inherited gene_id/gene_name
attribute on the same line, so the same-line regex match here is correct as designed.

Usage (from ~/repos/pilot-validation, inside the specimmune pixi env -- it already has both
minimap2, for immuannot.sh, and samtools>=1.10, for the trim step -- see pixi.toml):
  pixi run -e specimmune -- python3 scripts/run_immuannot_person.py <person_id> [<person_id> ...] \\
      [--mount ~/mnt/aou-controlled] [--immuannot-dir ~/tools/Immuannot] \\
      [--refdir ~/tools/Immuannot_refdata] [--outroot ~/pipeline_outputs] [--threads 4] \\
      [--region chr6:29500000-33500000] [--time-budget-min 30]

First run: pass exactly ONE person_id and read the printed timing summary before deciding whether
to pass all 60 cohort ids in one invocation.
"""
import argparse
import concurrent.futures
import glob
import gzip
import os
import re
import subprocess
import sys
import threading
import time

import pandas as pd

BUCKET_PREFIX = "gs://vwb-aou-datasets-controlled/"
LR_MANIFEST = "v9/wgs/long_read/manifest.tsv"
FA_COLS = {"hap1": "assembly_hap1_fa", "hap2": "assembly_hap2_fa"}
ALN_COLS = {"hap1": "assembly_hap1_aln2_hg38_bam", "hap2": "assembly_hap2_aln2_hg38_bam"}
PAF_COLS = {"hap1": "assembly_hap1_aln2_hg38_paf", "hap2": "assembly_hap2_aln2_hg38_paf"}
DEFAULT_REGION = "chr6:29500000-33500000"  # the project's standing HLA window (ENVIRONMENT.md)
# Reuses the exact number ENVIRONMENT.md/EXPERIMENTS.md already proved safe for SpecImmune-LR
# (pad100k: zero degradation across 8 genes, 3.9x smaller than full-pad) -- Marc, 2026-07-21:
# "a really safe reduction" for the analogous concern here (not clipping a gene near a coordinate-
# mapping boundary), not a new untested number.
DEFAULT_PAD = 100_000
# Tier 3 (self-align) constants -- untested at scale as of 2026-08-04, see --enable-self-align-
# -fallback below. sequel2 (991 people) has assembly_hap{1,2}_fa but no assembly_hap{1,2}_aln2_hg38_
# {bam,paf} -- reports/pre_sprint/lr_data_census/README.md. Rather than feed Immuannot the whole ~3.1Gb/hap
# untrimmed assembly (which the BRIEF.md fallback description originally implied -- estimated
# ~700-750x bigger than the normal ~4.2Mb trimmed input, a real risk of blowing up time/cost, per
# Marc's 2026-08-04 concern), self-align each hap FASTA against chr6 ONLY (not the whole genome)
# with minimap2, producing our own synthetic .paf in the exact format regions_from_paf() already
# parses -- then trim exactly as Tier 1 does. Turns "whole assembly" into "one extra minimap2 pass
# + a normal-sized trim." Needs a chr6-only reference slice, cached once from the full hg38
# reference already on the VM (ENVIRONMENT.md: ~/ref/Homo_sapiens_assembly38.fasta).
HG38_REF_DEFAULT = os.path.expanduser("~/ref/Homo_sapiens_assembly38.fasta")
CHR6_REF_CACHE_DEFAULT = os.path.expanduser("~/ref/chr6.fasta")
MINIMAP2_ASM_PRESET = "asm5"  # same-species assembly-vs-reference alignment, <5% divergence

# Tier 3 v2 (2026-10-11, S05 WS-K, scripts/hla_popgen/63_sequel2_recall.py; diagnosis in
# sprints/S05_qc_phasing_kir_trb/SEQUEL2_DIAGNOSIS.md). The 2026-08-05 test (EXPERIMENTS.md) failed
# for four reasons, each addressed here:
#   A. trim_assembly() runs `samtools faidx` on the source .fa.gz; sequel2's are (very likely) plain
#      gzip, which faidx cannot index -> every hap SKIPped. Tier 3 now extracts with a streaming
#      reader (stream_extract_fasta) that reads plain gzip, BGZF or plain FASTA, and never indexes.
#   B. self-align against chr6 ONLY gives contigs from other chromosomes no true locus, so repeat hits
#      pile onto chr6 and MAPQ cannot be trusted -> ~1.3-1.9 GB "trimmed" output. Tier 3 is now two
#      stages: (1) a cheap prefilter against a small slice (target region +/- flank) keeps contigs
#      with a long PRIMARY block; (2) only those contigs are re-aligned to the WHOLE primary GRCh38
#      (--self-align-confirm-ref), where MAPQ is meaningful again, and the trim comes from stage 2
#      blocks that are primary, long, MAPQ>=20 and on the target region. A size guard refuses any
#      trimmed input above --self-align-max-mb instead of feeding it to immuannot.sh.
#   C. regions_from_paf() spanned min(qstart)..max(qend) of every overlapping block, so one stray
#      block ballooned a contig to full length. It now takes optional filters and per-contig block
#      clustering (defaults unchanged, so Tier 1 behaves exactly as in production).
#   D. ensure_chr6_ref() hardcoded "chr6" (the KIR sequel2 bug: a chr19 --region was self-aligned to
#      chr6) and shared one <cache>.tmp across concurrent workers. The reference slice is now derived
#      from --region, the cache header is verified against the region's chromosome, and the temp
#      file is unique per process.
SELF_ALIGN_FLANK_DEFAULT = 1_000_000
SELF_ALIGN_MIN_BLOCK_DEFAULT = 10_000
SELF_ALIGN_MIN_MAPQ_DEFAULT = 20
SELF_ALIGN_MIN_IDENT_DEFAULT = 0.5   # PAF col10/col11 without -c: a seed-coverage proxy, not identity
# Stage 1 aligns the whole assembly to a ~6 Mb slice, so repeats chain into long low-col10/col11
# blocks. 63 smoke (2026-10-12, one HLA hap): 17,011 of 17,053 primary blocks >= 10 kb sat below
# 0.8 and pulled in ~10,000 contigs (1.5 GB); real MHC contigs sat at 0.95-0.98 (42 contigs >= 0.8).
SELF_ALIGN_STAGE1_MIN_IDENT_DEFAULT = 0.8
SELF_ALIGN_CLUSTER_GAP_DEFAULT = 250_000
SELF_ALIGN_MAX_MB_DEFAULT = 25.0     # normal trimmed input: ~4 MB/hap (HLA), ~1 MB/hap (KIR)
SELF_ALIGN_MAX_CANDIDATE_MB_DEFAULT = 600.0
SELF_ALIGN_REF_DIR_DEFAULT = os.path.expanduser("~/ref")


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def strip_bucket(uri):
    """Mirrors build_experiment_d_cohort.py's strip_bucket -- same wrong-bucket trap guard."""
    if not isinstance(uri, str) or not uri:
        return None
    if "fc-aou-datasets-controlled" in uri:
        return None
    if uri.startswith(BUCKET_PREFIX):
        return uri[len(BUCKET_PREFIX):]
    if uri.startswith("gs://"):
        return None
    return uri


def resolve_cols(lr, person_id, mount, cols):
    """Returns {'hap1': rel_path_or_None, 'hap2': rel_path_or_None} for the given column map --
    existence-checked directly against the mount, never assumed from a non-null manifest cell.
    Returns None (not a dict) if the person has no manifest row at all."""
    rows = lr[lr["research_id"] == str(person_id)]
    if rows.empty:
        return None
    out = {}
    for hap, col in cols.items():
        if col not in lr.columns:
            out[hap] = None
            continue
        found = None
        for v in rows[col].dropna():
            rel = strip_bucket(v)
            if rel and os.path.exists(os.path.join(mount, rel)):
                found = rel
                break
        out[hap] = found
    return out


def contigs_overlapping_region(bam_path, region):
    """Returns (sorted unique contig names, stderr) via samtools view on the hap's
    assembly-to-hg38 alignment BAM. contigs is None (not empty) on a real samtools failure --
    callers must tell that apart from "ran fine, zero contigs overlap" (unexpected either way,
    but a different problem: a broken index/region string vs. a genuinely empty result)."""
    proc = subprocess.run(["samtools", "view", bam_path, region], capture_output=True, text=True)
    if proc.returncode != 0:
        return None, proc.stderr
    contigs = sorted({line.split("\t")[0] for line in proc.stdout.splitlines() if line.strip()})
    return contigs, proc.stderr


def parse_region(region_str):
    chrom, coords = region_str.split(":")
    start_str, end_str = coords.split("-")
    return chrom, int(start_str), int(end_str)


def _paf_is_primary(fields):
    """tp:A:P = primary. A line without a tp tag is treated as primary (cannot tell)."""
    for tag in fields[12:]:
        if tag.startswith("tp:A:"):
            return tag == "tp:A:P"
    return True


def _cluster_blocks(blocks, max_gap):
    """blocks: [(qstart, qend)]. Merges blocks whose query gap is <= max_gap; returns
    [(start, end, aligned_bp)] clusters."""
    clusters = []
    for s, e in sorted(blocks):
        if clusters and s - clusters[-1][1] <= max_gap:
            cs, ce, bp = clusters[-1]
            clusters[-1] = (cs, max(ce, e), bp + (e - s))
        else:
            clusters.append((s, e, e - s))
    return clusters


def regions_from_paf(paf_path, chrom, region_start, region_end, pad, min_block=0, min_mapq=0,
                     min_ident=0.0, primary_only=False, cluster_gap=None):
    """True sub-range trim (2026-07-21, Marc: explored the padding-reduction margin properly
    rather than dismissing it -- see EXPERIMENTS.md's SpecImmune-LR pad100k sweep for the
    precedent this reuses). A .paf line's own qstart/qend/tstart/tend are the REAL aligned
    boundaries of that block (not an extrapolated estimate), so padding here only has to cover
    minor internal-indel drift and the chance a gene sits just past a block's reported edge --
    not a coordinate system we're guessing at. For each contig with >=1 PAF block overlapping the
    target region, takes the union of (min qstart, max qend) across its overlapping blocks, pads
    by `pad` on each side, and clamps to [0, qlen].

    Returns (regions, qlens, seen_tnames):
      regions: {contig: (start0, end0)} -- padded/clamped sub-range, ready for a faidx region string
      qlens: {contig: qlen} -- full contig length, from the PAF itself, no extra samtools call
      seen_tnames: set of every target name actually seen (for diagnosing a naming mismatch,
                   e.g. "chr6" vs "6", rather than silently returning empty and looking like "no
                   overlap" when it's really "wrong name assumed").

    Optional filters (2026-10-11, Tier 3 v2; all off by default, so the Tier 1 call is unchanged):
      primary_only  keep only tp:A:P blocks;
      min_block     PAF col11 (alignment block length) >= this;
      min_mapq      PAF col12 >= this (only meaningful against a whole-genome reference);
      min_ident     PAF col10/col11 >= this;
      cluster_gap   instead of min(qstart)..max(qend) over all blocks, merge blocks whose query gap
                    is <= cluster_gap and keep only the cluster with the most aligned bases.
    """
    opener = gzip.open if paf_path.endswith(".gz") else open
    accept_chrom = {chrom, chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"}
    filtering = primary_only or min_block > 0 or min_mapq > 0 or min_ident > 0
    regions, qlens, seen_tnames = {}, {}, set()
    blocks = {}
    with opener(paf_path, "rt") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qname, qlen, qstart, qend = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
            tname, tstart, tend = fields[5], int(fields[7]), int(fields[8])
            seen_tnames.add(tname)
            if tname not in accept_chrom:
                continue
            if tend <= region_start or tstart >= region_end:
                continue
            if filtering:
                nmatch, blen, mapq = int(fields[9]), int(fields[10]), int(fields[11])
                if primary_only and not _paf_is_primary(fields):
                    continue
                if blen < min_block or mapq < min_mapq or nmatch / max(blen, 1) < min_ident:
                    continue
            qlens[qname] = qlen
            if cluster_gap is not None:
                blocks.setdefault(qname, []).append((qstart, qend))
                continue
            new_start, new_end = max(0, qstart - pad), min(qlen, qend + pad)
            if qname in regions:
                s, e = regions[qname]
                regions[qname] = (min(s, new_start), max(e, new_end))
            else:
                regions[qname] = (new_start, new_end)
    for qname, bl in blocks.items():
        s, e, _bp = max(_cluster_blocks(bl, cluster_gap), key=lambda c: (c[2], -c[0]))
        regions[qname] = (max(0, s - pad), min(qlens[qname], e + pad))
    return regions, qlens, seen_tnames


def candidate_contigs_from_paf(paf_path, min_block, primary_only=True, min_ident=0.0):
    """Tier 3 stage 1: {contig: qlen} for every query contig with >=1 block of length >= min_block
    and col10/col11 >= min_ident (primary only by default) against ANY target in the PAF. Used with
    a small reference slice (target region +/- flank), so the target names are slice names and
    coordinates are not needed."""
    opener = gzip.open if paf_path.endswith(".gz") else open
    out = {}
    with opener(paf_path, "rt") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            if primary_only and not _paf_is_primary(fields):
                continue
            blen = int(fields[10])
            if blen < min_block or int(fields[9]) < min_ident * blen:
                continue
            out[fields[0]] = int(fields[1])
    return out


def is_bgzf(path):
    """True if `path` starts with a BGZF block header (gzip + FEXTRA subfield 'BC'). samtools faidx
    can index BGZF but not plain gzip ("Cannot index files compressed with gzip, please use bgzip",
    EXPERIMENTS.md 2026-08-05)."""
    try:
        with open(path, "rb") as f:
            h = f.read(18)
    except OSError:
        return False
    return (len(h) >= 18 and h[0:2] == b"\x1f\x8b" and h[2] == 8 and bool(h[3] & 4)
            and h[12:14] == b"BC")


_TARGET_RE = re.compile(r"^(.+):(\d+)-(\d+)$")


def _parse_target(t):
    """'ctg' -> ('ctg', None, None); 'ctg:S-E' (1-based, inclusive, samtools style) -> (ctg, S, E)."""
    m = _TARGET_RE.match(t)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    return t, None, None


def stream_extract_fasta(fa_path, targets, out_fa_path, width=60):
    """Writes each target ('ctg' or 'ctg:S-E') to out_fa_path in the given order, with the same
    header and 60-column wrapping `samtools faidx` produces, by reading the source FASTA once
    (plain gzip, BGZF or uncompressed). Never writes an index, so it works on the read-only mount
    and on plain-gzip sequel2 assemblies. Stops reading once every wanted contig has been seen.
    Returns (ok, message)."""
    wanted = {}
    for i, t in enumerate(targets):
        name, s, e = _parse_target(t)
        wanted.setdefault(name, []).append((i, s, e))
    found = {}
    opener = gzip.open if fa_path.endswith(".gz") else open

    def take(name, chunks):
        seq = "".join(chunks)
        for i, s, e in wanted[name]:
            if s is None:
                found[i] = (name, seq)
            else:
                found[i] = (f"{name}:{s}-{e}", seq[s - 1:min(e, len(seq))])

    try:
        with opener(fa_path, "rt") as f:
            name, chunks = None, []
            for line in f:
                if line.startswith(">"):
                    if name in wanted:
                        take(name, chunks)
                        if len(found) == len(targets):
                            name = None
                            break
                    name, chunks = line[1:].split()[0], []
                elif name in wanted:
                    chunks.append(line.strip())
            if name in wanted:
                take(name, chunks)
    except (OSError, EOFError, ValueError) as e:
        return False, f"read error on {fa_path}: {e!r}"
    missing = [targets[i] for i in range(len(targets)) if i not in found]
    if missing:
        return False, f"{len(missing)} target(s) not found in {fa_path}: {missing[:5]}"
    with open(out_fa_path, "w") as out:
        for i in range(len(targets)):
            header, seq = found[i]
            out.write(f">{header}\n")
            for j in range(0, len(seq), width):
                out.write(seq[j:j + width] + "\n")
    return os.path.getsize(out_fa_path) > 0, ""


def fasta_seq_names(path):
    """Sequence names of a FASTA, from its .fai if present, else by scanning headers."""
    fai = path + ".fai"
    if os.path.exists(fai):
        with open(fai) as f:
            return [line.split("\t")[0] for line in f if line.strip()]
    if path.endswith(".mmi"):
        return []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        return [line[1:].split()[0] for line in f if line.startswith(">")]


def _chrom_of(seq_name):
    """'chr6' -> 'chr6'; 'chr6:28500001-34500000' (a samtools faidx slice header) -> 'chr6'."""
    return _parse_target(seq_name)[0]


def ensure_ref_slice(hg38_ref, cache_path, region_or_chrom):
    """Extracts `region_or_chrom` (e.g. 'chr19' or 'chr6:28500001-34500000') from the full
    reference once and caches it. Idempotent. Unlike the old ensure_chr6_ref(): the extracted name
    comes from the caller (not hardcoded), an existing cache is only reused if its first header is
    the expected name, and the temp file is unique per process and thread (no shared <cache>.tmp
    race between concurrent workers). Returns the cache path, or None."""
    want = region_or_chrom
    if os.path.exists(cache_path) and os.path.getsize(cache_path) > 0:
        with open(cache_path) as f:
            first = f.readline()
        got = first[1:].split()[0] if first.startswith(">") else ""
        if got != want:
            print(f"    WARNING: reference cache {cache_path} holds '{got}', not '{want}' -- refusing "
                  f"to use it for this region (pass a matching cache path).", file=sys.stderr)
            return None
        return cache_path
    if not os.path.exists(hg38_ref):
        return None
    os.makedirs(os.path.dirname(cache_path) or ".", exist_ok=True)
    tmp_path = f"{cache_path}.tmp.{os.getpid()}.{threading.get_ident()}"
    with open(tmp_path, "w") as out:
        proc = subprocess.run(["samtools", "faidx", hg38_ref, want], stdout=out,
                              stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0 or os.path.getsize(tmp_path) == 0:
        print(f"    WARNING: could not extract {want} from {hg38_ref}: {proc.stderr}", file=sys.stderr)
        os.remove(tmp_path)
        return None
    os.replace(tmp_path, cache_path)  # atomic; a concurrent identical build just replaces it
    return cache_path


def ensure_chr6_ref(hg38_ref, cache_path):
    """Extracts chr6 alone from the full hg38 reference, once, and caches it -- so self-align
    only ever indexes ~170Mb (one chromosome), not the full ~3.1Gb genome, keeping Tier 3's cost
    close to a normal alignment step rather than a whole-genome-vs-whole-genome pass. Idempotent:
    a second call with an existing, non-empty cache is a no-op. Kept for callers of the old name;
    now a thin wrapper over ensure_ref_slice() (unique temp file, header check)."""
    return ensure_ref_slice(hg38_ref, cache_path, "chr6")


def self_align_paf(fa_path, chr6_ref, out_paf_path, threads, minimap2="minimap2"):
    """Runs minimap2 in PAF-output mode (no -a) of fa_path (query) against a reference (a slice,
    a whole-genome FASTA or a prebuilt .mmi). Output is the PAF format regions_from_paf() parses.
    Returns (ok, stderr)."""
    with open(out_paf_path, "w") as out:
        proc = subprocess.run(
            [minimap2, "-x", MINIMAP2_ASM_PRESET, "-t", str(threads), chr6_ref, fa_path],
            stdout=out, stderr=subprocess.PIPE, text=True)
    ok = proc.returncode == 0 and os.path.getsize(out_paf_path) > 0
    return ok, proc.stderr


def self_align_stage1_ref(args, chrom, rstart, rend):
    """Stage-1 reference for Tier 3: --self-align-ref if given; else an explicitly passed legacy
    --chr6-ref-cache; else a cached slice of the target region +/- --self-align-flank, derived
    from --region (so a chr19 KIR run gets a chr19 slice). Returns (path or None, reason)."""
    if args.self_align_ref:
        ref = args.self_align_ref if os.path.exists(args.self_align_ref) else None
        return ref, f"--self-align-ref {args.self_align_ref} not found"
    if args.chr6_ref_cache:  # legacy flag, only when passed explicitly
        ref = args.chr6_ref_cache if os.path.exists(args.chr6_ref_cache) else None
        return ref, f"--chr6-ref-cache {args.chr6_ref_cache} not found"
    s = max(1, rstart - args.self_align_flank + 1)
    e = rend + args.self_align_flank
    region = f"{chrom}:{s}-{e}"
    cache = os.path.join(args.self_align_ref_dir, f"selfalign_{chrom}_{s}_{e}.fa")
    ref = ensure_ref_slice(args.hg38_ref, cache, region)
    return ref, f"could not prepare {region} from {args.hg38_ref} into {cache}"


def _remove_quietly(path):
    try:
        os.remove(path)
    except OSError:
        pass


def self_align_targets(fa_path, hap, person_dir, args, threads, row):
    """Tier 3 v2 (see the module constants comment). Returns (targets, trim_method, trim_source)
    or None after printing a SKIP reason. trim_source is the local candidate-contig FASTA the trim
    must read from (stream_extract_fasta), never the mount."""
    chrom, rstart, rend = parse_region(args.region)
    accept = {chrom, chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"}
    if not args.self_align_confirm_ref:
        print(f"  {hap}: SKIP -- Tier 3 needs --self-align-confirm-ref (whole primary GRCh38, FASTA "
              f"or .mmi): a single-chromosome self-align was shown to pull in the whole assembly "
              f"(EXPERIMENTS.md 2026-08-05).", file=sys.stderr)
        return None
    if not os.path.exists(args.self_align_confirm_ref):
        print(f"  {hap}: SKIP -- --self-align-confirm-ref {args.self_align_confirm_ref} not found.",
              file=sys.stderr)
        return None
    ref1, why = self_align_stage1_ref(args, chrom, rstart, rend)
    if ref1 is None:
        print(f"  {hap}: SKIP -- {why}.", file=sys.stderr)
        return None
    names = fasta_seq_names(ref1)
    if not any(_chrom_of(n) in accept for n in names):
        print(f"  {hap}: SKIP -- stage-1 self-align reference {ref1} holds {names[:5]}, not "
              f"{chrom} (the region's chromosome). This is the sequel2 KIR bug: a chr19 region "
              f"self-aligned to a chr6 reference.", file=sys.stderr)
        return None

    # Stage 1: cheap prefilter against the region slice.
    t0 = time.perf_counter()
    paf1 = os.path.join(person_dir, f"{hap}.self_align.paf")
    ok, err = self_align_paf(fa_path, ref1, paf1, threads, args.minimap2)
    row["self_align_seconds"] = round(time.perf_counter() - t0, 2)
    if not ok:
        print(f"  {hap}: SKIP -- stage-1 self-align failed after {row['self_align_seconds']}s: {err}",
              file=sys.stderr)
        return None
    cands = candidate_contigs_from_paf(paf1, args.self_align_min_block, primary_only=True,
                                       min_ident=args.self_align_stage1_min_ident)
    cand_bp = sum(cands.values())
    row["self_align_candidates"] = len(cands)
    row["self_align_candidate_mb"] = round(cand_bp / 1e6, 2)
    print(f"  {hap}: stage 1 ({row['self_align_seconds']}s vs {os.path.basename(ref1)}): "
          f"{len(cands)} candidate contig(s), {cand_bp/1e6:.1f} MB, with a primary block >= "
          f"{args.self_align_min_block:,} bp at col10/col11 >= {args.self_align_stage1_min_ident:g}",
          file=sys.stderr)
    if not cands:
        print(f"  {hap}: SKIP -- no contig has a primary block >= {args.self_align_min_block:,} bp "
              f"against the {chrom} region slice.", file=sys.stderr)
        return None
    if cand_bp > args.self_align_max_candidate_mb * 1e6:
        print(f"  {hap}: SKIP -- candidates total {cand_bp/1e6:.0f} MB > --self-align-max-candidate-mb "
              f"{args.self_align_max_candidate_mb:g}; looks like whole-assembly pull-in, not the "
              f"region.", file=sys.stderr)
        return None

    # Stage 2: candidates only, against the whole primary genome (MAPQ meaningful again).
    t1 = time.perf_counter()
    cand_fa = os.path.join(person_dir, f"{hap}.self_align_candidates.fa")
    ok, err = stream_extract_fasta(fa_path, sorted(cands), cand_fa)
    if not ok:
        print(f"  {hap}: SKIP -- could not extract candidate contigs: {err}", file=sys.stderr)
        return None
    paf2 = os.path.join(person_dir, f"{hap}.self_align_confirm.paf")
    ok, err = self_align_paf(cand_fa, args.self_align_confirm_ref, paf2, threads, args.minimap2)
    row["confirm_seconds"] = round(time.perf_counter() - t1, 2)
    if not ok:
        print(f"  {hap}: SKIP -- stage-2 whole-genome alignment failed: {err}", file=sys.stderr)
        _remove_quietly(cand_fa)
        return None
    regions, qlens, seen = regions_from_paf(
        paf2, chrom, rstart, rend, args.pad, min_block=args.self_align_min_block,
        min_mapq=args.self_align_min_mapq, min_ident=args.self_align_min_ident, primary_only=True,
        cluster_gap=args.self_align_cluster_gap)
    if not regions:
        print(f"  {hap}: SKIP -- no candidate contig has a primary, MAPQ>={args.self_align_min_mapq}, "
              f">= {args.self_align_min_block:,} bp block on {args.region} against the whole genome "
              f"(targets seen: {sorted(seen)[:10]}).", file=sys.stderr)
        _remove_quietly(cand_fa)
        return None
    padded_bp = sum(e - s for s, e in regions.values())
    row["n_contigs"] = len(regions)
    row["whole_contig_mb"] = round(sum(qlens[c] for c in regions) / 1e6, 2)
    row["padded_mb"] = round(padded_bp / 1e6, 2)
    if padded_bp > args.self_align_max_mb * 1e6:
        print(f"  {hap}: SKIP -- confirmed sub-range totals {padded_bp/1e6:.1f} MB > "
              f"--self-align-max-mb {args.self_align_max_mb:g} (expected a few MB); not feeding it "
              f"to immuannot.sh.", file=sys.stderr)
        _remove_quietly(cand_fa)
        return None
    print(f"  {hap}: stage 2 ({row['confirm_seconds']}s, whole genome): {len(regions)} contig(s) "
          f"confirmed on {args.region}, padded sub-range {padded_bp/1e6:.2f} MB", file=sys.stderr)
    targets = [f"{c}:{s + 1}-{e}" for c, (s, e) in regions.items()]
    return targets, "self_align_confirmed", cand_fa


def trim_assembly(fa_path, targets, out_fa_path):
    """samtools faidx extracts each entry in `targets` -- either a bare contig name (whole contig,
    the fallback path) or a "contig:start-end" region string (the preferred sub-range path via
    regions_from_paf) -- samtools faidx accepts both forms interchangeably, so this function
    doesn't need to know which kind it got. fa_path lives on the read-only gcsfuse mount, but
    samtools faidx's default index locations are alongside fa_path itself -- confirmed 2026-07-21
    (person 1001871) this fails with "Permission denied" building the .gzi/.fai next to the
    read-only source. Fix: --fai-idx/--gzi-idx redirect the index files into out_fa_path's own
    (writable) directory instead -- samtools still just READS fa_path off the mount, only the
    small index files themselves move. Confirmed present in this VM's samtools via
    `samtools faidx --help` before relying on it, not assumed."""
    fai_idx = out_fa_path + ".src.fai"
    gzi_idx = out_fa_path + ".src.gzi"
    with open(out_fa_path, "w") as out:
        proc = subprocess.run(
            ["samtools", "faidx", "--fai-idx", fai_idx, "--gzi-idx", gzi_idx, fa_path] + targets,
            stdout=out, stderr=subprocess.PIPE, text=True)
    ok = proc.returncode == 0 and os.path.getsize(out_fa_path) > 0
    return ok, proc.stderr


def find_immuannot_script(immuannot_dir):
    """Never hardcode the in-repo path -- confirmed 2026-07-21 the upstream repo ships immuannot.sh
    under a VERSIONED folder (scripts.pub.v3, not plain scripts/ as its own README shows), and that
    name could drift again on a future re-clone. Same lesson as ENVIRONMENT.md quirk #13: find the
    real file, don't trust a documented path pattern."""
    matches = glob.glob(os.path.join(immuannot_dir, "**", "immuannot.sh"), recursive=True)
    if not matches:
        return None
    matches.sort(key=lambda p: p.count(os.sep))  # prefer the shallowest match
    return matches[0]


def run_immuannot(script, refdir, contig_path, outprefix, threads):
    cmd = ["bash", script, "-c", contig_path, "-r", refdir, "-o", outprefix, "-t", str(threads)]
    print(f"    Running: {' '.join(cmd)}", file=sys.stderr)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    # ALWAYS persist immuannot.sh's own stdout/stderr to disk, not just on failure and not just to
    # the terminal (2026-07-22, Marc: the first 60-person run's failure logs were lost to the
    # terminal scrollback after the VM restarted, because they were only printed, never written).
    # This log survives a restart (persistent home disk) and lets a failed hap be diagnosed -- or
    # cheaply re-run on the already-trimmed FASTA -- later, instead of re-doing the whole pipeline.
    log_path = outprefix + ".immuannot.log"
    try:
        with open(log_path, "w") as lf:
            lf.write(f"# cmd: {' '.join(cmd)}\n# exit_code: {proc.returncode}\n"
                     f"# ---- STDOUT ----\n{proc.stdout}\n# ---- STDERR ----\n{proc.stderr}\n")
    except OSError as e:
        print(f"    WARNING: could not write immuannot log {log_path}: {e}", file=sys.stderr)
    expected = outprefix + ".gtf.gz"
    if not os.path.exists(expected):
        # Never trust exit code alone (ENVIRONMENT.md quirks #17/#18) -- surface the real logs.
        print(f"    WARNING: expected output {expected} missing after run "
              f"(exit code {proc.returncode}). Full log saved to {log_path}. Tail:", file=sys.stderr)
        print(proc.stdout[-1500:], file=sys.stderr)
        print(proc.stderr[-1500:], file=sys.stderr)
        return None
    return expected


def parse_gtf(gtf_gz_path, show_raw_sample=True):
    """Extracts (gene, allele) pairs from Immuannot's gtf.gz. Prints the first record's raw
    attributes so the field spelling can be eyeballed against reference/README_Immuannot.md."""
    calls = {}
    printed_sample = False
    with gzip.open(gtf_gz_path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = fields[8]
            if show_raw_sample and not printed_sample:
                print(f"    [sanity check] first GTF attribute string seen: {attrs}", file=sys.stderr)
                printed_sample = True
            gene_m = re.search(r'gene_name "([^"]+)"', attrs) or re.search(r'gene_id "([^"]+)"', attrs)
            allele_m = re.search(r'consensus "([^"]+)"', attrs) or re.search(r'allele "([^"]+)"', attrs)
            if gene_m and allele_m:
                calls[gene_m.group(1)] = allele_m.group(1)
    return calls


def process_haplotype(pid, hap, fa_rel, paf_rel, aln_rel, person_dir, immuannot_script,
                       args, resolve_seconds, threads):
    """One haplotype's full trim + immuannot.sh + parse, fully self-contained (returns its own
    row/calls rather than mutating shared state) so it's safe to run hap1 and hap2 concurrently
    in separate threads (Marc, 2026-07-21: why sum the two haps instead of running them at once?)
    without any locking. `threads` is passed in explicitly (not read from args) because it differs
    between sequential mode (args.threads, e.g. 4) and --parallel-haps mode (args.threads split
    across the two concurrent haps, e.g. 2 each)."""
    row = {"person_id": pid, "hap": hap, "resolve_seconds": round(resolve_seconds, 2),
           "trim_method": None, "n_contigs": None,
           "whole_contig_mb": None, "padded_mb": None, "trimmed_mb": None,
           "contig_lookup_seconds": None, "trim_seconds": None, "immuannot_seconds": None,
           "hap_total_seconds": None, "self_align_seconds": None}
    hap_t0 = time.perf_counter()
    fa_path = os.path.join(args.mount, fa_rel)

    # --- Stage 1: find trim targets. Prefer exact sub-ranges from the .paf (Marc, 2026-07-21:
    # explored the padding-reduction margin properly -- reuses the proven-safe pad100k from
    # the SpecImmune-LR sweep). Falls back to whole-contig via the .bam if this person's .paf
    # is missing or has nothing recognizable overlapping the region, so nobody gets skipped
    # just because the (optional) .paf isn't there. ---
    targets, trim_method, trim_source = None, None, None
    if paf_rel:
        paf_path = os.path.join(args.mount, paf_rel)
        chrom, rstart, rend = parse_region(args.region)
        regions, qlens, seen_tnames = regions_from_paf(paf_path, chrom, rstart, rend, args.pad)
        if regions:
            targets = [f"{c}:{s + 1}-{e}" for c, (s, e) in regions.items()]
            trim_method = "paf_region"
            whole_bp = sum(qlens[c] for c in regions)
            padded_bp = sum(e - s for s, e in regions.values())
            row["n_contigs"] = len(regions)
            row["whole_contig_mb"] = round(whole_bp / 1e6, 2)
            row["padded_mb"] = round(padded_bp / 1e6, 2)
            pct = 100 * (1 - padded_bp / max(whole_bp, 1))
            print(f"  {hap}: {len(regions)} contig(s) overlap {args.region} via .paf "
                  f"(pad={args.pad:,}bp) -- whole-contig would be {whole_bp/1e6:.1f} MB, "
                  f"padded sub-range is {padded_bp/1e6:.1f} MB ({pct:.0f}% smaller)",
                  file=sys.stderr)
        else:
            print(f"  {hap}: .paf present but 0 contigs overlap {args.region} -- target names "
                  f"seen in this .paf: {sorted(seen_tnames)[:10]} -- falling back to "
                  f"whole-contig via .bam.", file=sys.stderr)
    elif aln_rel:
        print(f"  {hap}: no .paf resolved for this person -- using whole-contig fallback via .bam.",
              file=sys.stderr)
    else:
        print(f"  {hap}: no .paf AND no aln-to-hg38 .bam resolved for this person (Tier 1/Tier 2 "
              f"both unavailable -- the sequel2 gap, ENVIRONMENT.md/lr_data_census).", file=sys.stderr)

    # --- Tier 3 (self-align, 2026-08-04, untested at scale -- see module-level comment): only
    # reached when NEITHER a real .paf NOR a real aln-to-hg38 .bam resolved for this haplotype
    # (the sequel2 gap). Self-aligns the raw assembly FASTA against a chr6-only reference slice
    # with minimap2, producing our own synthetic .paf that regions_from_paf() parses exactly like
    # a real one -- so this still gets a proper padded SUB-RANGE trim, not the whole ~3.1Gb/hap
    # assembly. Opt-in only (--enable-self-align-fallback) until validated on real sequel2 people.
    if targets is None and not aln_rel:
        if not args.enable_self_align_fallback:
            print(f"  {hap}: SKIP -- Tier 1/2 unavailable and --enable-self-align-fallback not "
                  f"passed (this is the untested sequel2 fallback path -- see run_immuannot_person.py "
                  f"module docstring / production_orchestrator/BRIEF.md).", file=sys.stderr)
            return row, {}
        t3 = self_align_targets(fa_path, hap, person_dir, args, threads, row)
        if t3 is None:
            return row, {}
        targets, trim_method, trim_source = t3

    if targets is None:
        aln_path = os.path.join(args.mount, aln_rel)
        contigs, err = contigs_overlapping_region(aln_path, args.region)
        if contigs is None:
            print(f"  {hap}: SKIP -- samtools view failed against the aln-to-hg38 BAM: {err}",
                  file=sys.stderr)
            return row, {}
        if not contigs:
            print(f"  {hap}: SKIP -- 0 contigs overlap {args.region} (checked both .paf and "
                  f".bam). Unexpected -- worth a closer look for this specific person, not "
                  f"assumed benign.", file=sys.stderr)
            return row, {}
        targets, trim_method = contigs, "bam_whole_contig"
        row["n_contigs"] = len(contigs)
        print(f"  {hap}: {len(contigs)} contig(s) overlap {args.region} via .bam "
              f"(whole-contig): {contigs}", file=sys.stderr)

    row["trim_method"] = trim_method
    lookup_t1 = time.perf_counter()
    lookup_seconds = lookup_t1 - hap_t0
    row["contig_lookup_seconds"] = round(lookup_seconds, 2)
    print(f"  {hap} [1/3] find trim targets ({trim_method}): {lookup_seconds:.2f}s", file=sys.stderr)

    # --- Stage 2: trim (samtools faidx -- files-to-ready-for-immuannot step Marc flagged) ---
    trimmed_fa = os.path.join(person_dir, f"{hap}.trimmed.fa")
    if trim_source is not None:
        # Tier 3: read the local candidate-contig FASTA, never index the (plain-gzip) source.
        ok, err = stream_extract_fasta(trim_source, targets, trimmed_fa)
        _remove_quietly(trim_source)  # whole candidate contigs; the trimmed FASTA is what is kept
    else:
        ok, err = trim_assembly(fa_path, targets, trimmed_fa)
    trim_t2 = time.perf_counter()
    trim_seconds = trim_t2 - lookup_t1
    row["trim_seconds"] = round(trim_seconds, 2)
    print(f"  {hap} [2/3] trim (samtools faidx): {trim_seconds:.2f}s", file=sys.stderr)

    if not ok:
        print(f"  {hap}: SKIP -- samtools faidx trim failed: {err}", file=sys.stderr)
        return row, {}
    trimmed_mb = os.path.getsize(trimmed_fa) / 1e6
    row["trimmed_mb"] = round(trimmed_mb, 1)
    print(f"  {hap}: trimmed FASTA is {trimmed_mb:.1f} MB, ready for immuannot.sh", file=sys.stderr)

    # --- Stage 3: immuannot.sh itself ---
    outprefix = os.path.join(person_dir, hap)
    gtf = run_immuannot(immuannot_script, args.refdir, trimmed_fa, outprefix, threads)
    hap_t3 = time.perf_counter()
    immuannot_seconds = hap_t3 - trim_t2
    row["immuannot_seconds"] = round(immuannot_seconds, 2)
    row["hap_total_seconds"] = round(hap_t3 - hap_t0, 2)
    calls = parse_gtf(gtf) if gtf else {}
    print(f"  {hap} [3/3] immuannot.sh (threads={threads}): {immuannot_seconds/60:.1f} min "
          f"({len(calls)} genes called)", file=sys.stderr)
    print(f"  {hap} stage breakdown -- lookup {lookup_seconds:.1f}s / trim {trim_seconds:.1f}s / "
          f"immuannot {immuannot_seconds/60:.1f}min / hap total {row['hap_total_seconds']/60:.1f}min",
          file=sys.stderr)
    return row, calls


def process_person(pid, lr, args, immuannot_script):
    """Returns (gene_rows, timing_rows) for one person. Every stage is timestamped separately
    (manifest resolution, contig lookup, trim, immuannot.sh) and printed as it completes -- Marc's
    2026-07-21 ask, to see exactly which stage the "BAM things" (contig lookup / trim) actually
    cost versus immuannot.sh itself, rather than one bundled number. A timing row is recorded for
    every haplotype attempt, including ones that get skipped partway through, so a slow-then-failing
    stage is still visible in immuannot_timing.tsv, not silently dropped.

    hap1/hap2 run SEQUENTIALLY by default (each gets the full --threads) -- total person time is
    the SUM of both haps' immuannot.sh runtime, not the slower one alone. Pass --parallel-haps to
    run them concurrently instead (each then gets --threads/2, minimum 1, so the two together don't
    oversubscribe the VM's cores) -- Marc, 2026-07-21: tested on person 1001871 before deciding
    whether to use this for the full 60-person run, same "test then decide" discipline as the pad
    sweep. Concurrent stderr from both haps interleaves in the terminal -- expected, not a bug."""
    print(f"\n=== Person {pid} ===", file=sys.stderr)
    person_t0 = time.perf_counter()

    resolve_t0 = time.perf_counter()
    fa_paths = resolve_cols(lr, pid, args.mount, FA_COLS)
    aln_paths = resolve_cols(lr, pid, args.mount, ALN_COLS)
    resolve_seconds = time.perf_counter() - resolve_t0
    print(f"  manifest path resolution: {resolve_seconds:.2f}s", file=sys.stderr)

    paf_paths = resolve_cols(lr, pid, args.mount, PAF_COLS)  # optional -- not gated on below

    if fa_paths is None:
        print(f"  SKIP: no manifest row for {pid} at all.", file=sys.stderr)
        return [], []
    missing_fa = [h for h, p in fa_paths.items() if p is None]
    missing_aln = [h for h, p in aln_paths.items() if p is None]
    # assembly_hap*_fa is the ONE hard requirement (Immuannot's real dependency -- see module
    # docstring). Missing aln-to-hg38 BAM is no longer an automatic skip (2026-08-04) -- it's the
    # sequel2 gap that Tier 3 (self-align, --enable-self-align-fallback) exists to cover. Without
    # that flag, a missing aln still resolves to a per-haplotype SKIP inside process_haplotype
    # (clear reason logged there), not a whole-person skip here -- so a person missing only ONE
    # haplotype's aln still gets the other haplotype processed.
    if missing_fa:
        print(f"  SKIP: missing assembly FASTA {missing_fa} -- this person's platform likely lacks "
              f"assembly data entirely (only revio/sequel2e/sequel2 do -- "
              f"reports/pre_sprint/lr_data_census/README.md). Not treated as a bug -- flag which person_ids "
              f"hit this so the real usable fraction of the cohort is known.", file=sys.stderr)
        return [], []
    if missing_aln:
        print(f"  NOTE: missing aln-to-hg38 BAM for {missing_aln} -- Tier 1/2 unavailable for "
              f"{'that haplotype' if len(missing_aln) == 1 else 'both haplotypes'}; will "
              f"{'attempt Tier 3 self-align' if args.enable_self_align_fallback else 'SKIP (pass --enable-self-align-fallback to attempt Tier 3 instead)'}.",
              file=sys.stderr)

    person_dir = os.path.join(args.outroot, str(pid), "immuannot_output")
    os.makedirs(person_dir, exist_ok=True)

    hap_args = {
        hap: (pid, hap, fa_paths[hap], paf_paths.get(hap), aln_paths[hap], person_dir,
              immuannot_script, args, resolve_seconds)
        for hap in ("hap1", "hap2")
    }
    results = {}
    if args.parallel_haps:
        per_hap_threads = max(1, args.threads // 2)
        print(f"  --parallel-haps: running hap1+hap2 concurrently, {per_hap_threads} threads "
              f"each (of --threads {args.threads} total)", file=sys.stderr)
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            futures = {
                pool.submit(process_haplotype, *hap_args[hap], per_hap_threads): hap
                for hap in ("hap1", "hap2")
            }
            for fut in concurrent.futures.as_completed(futures):
                results[futures[fut]] = fut.result()
    else:
        for hap in ("hap1", "hap2"):
            results[hap] = process_haplotype(*hap_args[hap], args.threads)

    timing_rows = [results[hap][0] for hap in ("hap1", "hap2")]
    gene_calls = {hap: results[hap][1] for hap in ("hap1", "hap2")}

    genes = sorted(set(gene_calls.get("hap1", {})) | set(gene_calls.get("hap2", {})))
    gene_rows = [{
        "person_id": pid, "gene": gene,
        "immuannot_1": gene_calls.get("hap1", {}).get(gene, "NA"),
        "immuannot_2": gene_calls.get("hap2", {}).get(gene, "NA"),
    } for gene in genes]

    person_total = time.perf_counter() - person_t0
    budget = args.time_budget_min * 60
    verdict = "UNDER budget" if person_total <= budget else "OVER budget"
    print(f"  --- Person {pid} total: {person_total/60:.1f} min "
          f"(budget {args.time_budget_min} min) -- {verdict}. {len(genes)} genes called. ---",
          file=sys.stderr)
    return gene_rows, timing_rows


def write_incremental(df_new, out_path, key_cols):
    """Appends/overwrites rows for the people just processed, keeping prior people's rows from
    earlier invocations intact -- lets the eventual 60-person run be done in batches/reruns
    without losing already-computed results."""
    if os.path.exists(out_path):
        old = pd.read_csv(out_path, sep="\t", dtype=str)
        keys = set(df_new[key_cols[0]].astype(str))
        old = old[~old[key_cols[0]].astype(str).isin(keys)]
        df_new = pd.concat([old, df_new], ignore_index=True)
    df_new.to_csv(out_path, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("person_ids", nargs="+", help="One or more research_id values (cohort.tsv).")
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    ap.add_argument("--immuannot-dir", default=os.path.expanduser("~/tools/Immuannot"))
    ap.add_argument("--refdir", default=os.path.expanduser("~/tools/Immuannot_refdata"))
    ap.add_argument("--outroot", default=os.path.expanduser("~/pipeline_outputs"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--parallel-haps", action="store_true",
                    help="Run hap1 and hap2 concurrently instead of sequentially -- total "
                         "person time trends toward the slower haplotype instead of the sum of "
                         "both. Each hap gets --threads/2 (min 1) to avoid oversubscribing the "
                         "VM's cores. Untested at scale as of 2026-07-21 -- verify on 1-2 people "
                         "before using for a full 60-person run.")
    ap.add_argument("--region", default=DEFAULT_REGION,
                    help=f"GRCh38 region used to pick overlapping contigs (default {DEFAULT_REGION}).")
    ap.add_argument("--pad", type=int, default=DEFAULT_PAD,
                    help=f"bp padding on each side of a .paf-mapped block before extraction "
                         f"(default {DEFAULT_PAD:,} -- reuses the SpecImmune-LR pad100k value "
                         f"already proven safe in EXPERIMENTS.md). Only affects the .paf-based "
                         f"sub-range path; the .bam whole-contig fallback ignores this.")
    ap.add_argument("--time-budget-min", type=float, default=30,
                    help="Per-person minutes to compare against when printing the verdict "
                         "(informational only -- does not abort a run).")
    ap.add_argument("--enable-self-align-fallback", action="store_true",
                    help="Attempt Tier 3 (two-stage self-align: region slice, then whole primary "
                         "GRCh38 via --self-align-confirm-ref) for any haplotype missing BOTH .paf "
                         "and aln-to-hg38 .bam -- the "
                         "sequel2 gap (991 people, ~95%% AFR, per reports/pre_sprint/lr_data_census/README.md). "
                         "UNTESTED AT SCALE as of 2026-08-04 -- explicit smoke-test-first flag, not "
                         "the default. Without it, those haplotypes are SKIPped with a clear reason "
                         "instead of silently producing no calls.")
    ap.add_argument("--hg38-ref", default=HG38_REF_DEFAULT,
                    help=f"Full hg38 reference FASTA, used once to carve out the chr6-only slice "
                         f"for Tier 3 self-align (default {HG38_REF_DEFAULT}, per ENVIRONMENT.md).")
    ap.add_argument("--chr6-ref-cache", default=None,
                    help="LEGACY (pre-2026-10-11). If passed, used as the Tier 3 stage-1 reference "
                         "as is (its sequence names must include --region's chromosome, else the "
                         "haplotype is SKIPped). Default: unset -> a slice of --region +/- "
                         "--self-align-flank is cached under --self-align-ref-dir.")
    ap.add_argument("--self-align-ref", default=None,
                    help="Tier 3 stage-1 reference FASTA (e.g. the region slice prepared by "
                         "63_sequel2_recall.py --prepare-refs). Overrides the auto-built slice.")
    ap.add_argument("--self-align-ref-dir", default=SELF_ALIGN_REF_DIR_DEFAULT,
                    help="Where the auto-built stage-1 region slice is cached.")
    ap.add_argument("--self-align-flank", type=int, default=SELF_ALIGN_FLANK_DEFAULT,
                    help="bp added on each side of --region for the stage-1 slice.")
    ap.add_argument("--self-align-confirm-ref", default=None,
                    help="REQUIRED for Tier 3: whole primary GRCh38 (chr1-22,X,Y,M; no ALT/HLA/decoy "
                         "contigs) as FASTA or prebuilt minimap2 -x asm5 .mmi. Stage-1 candidate "
                         "contigs are re-aligned to it so MAPQ is meaningful.")
    ap.add_argument("--self-align-min-block", type=int, default=SELF_ALIGN_MIN_BLOCK_DEFAULT,
                    help="Minimum PAF alignment block length (col11) in both Tier 3 stages.")
    ap.add_argument("--self-align-min-mapq", type=int, default=SELF_ALIGN_MIN_MAPQ_DEFAULT,
                    help="Minimum MAPQ for a stage-2 (whole-genome) block.")
    ap.add_argument("--self-align-min-ident", type=float, default=SELF_ALIGN_MIN_IDENT_DEFAULT,
                    help="Minimum col10/col11 for a stage-2 block (approximate without -c).")
    ap.add_argument("--self-align-stage1-min-ident", type=float,
                    default=SELF_ALIGN_STAGE1_MIN_IDENT_DEFAULT,
                    help="Minimum col10/col11 for a stage-1 (region-slice) block; drops the long "
                         "low-identity repeat chains that pull in the whole assembly.")
    ap.add_argument("--self-align-cluster-gap", type=int, default=SELF_ALIGN_CLUSTER_GAP_DEFAULT,
                    help="Stage-2 blocks on one contig closer than this (query bp) form one "
                         "cluster; only the largest cluster per contig is trimmed.")
    ap.add_argument("--self-align-max-mb", type=float, default=SELF_ALIGN_MAX_MB_DEFAULT,
                    help="Refuse a Tier 3 trimmed input larger than this (MB per haplotype).")
    ap.add_argument("--self-align-max-candidate-mb", type=float,
                    default=SELF_ALIGN_MAX_CANDIDATE_MB_DEFAULT,
                    help="Refuse stage 2 if stage-1 candidate contigs total more than this (MB).")
    ap.add_argument("--minimap2", default="minimap2", help="minimap2 binary for Tier 3.")
    ap.add_argument("--force", action="store_true",
                    help="Re-process a person even if already present in immuannot_calls.tsv "
                         "(default: skip -- makes a 60-person run safe to kill and re-launch "
                         "with the exact same command).")
    ap.add_argument("--out-suffix", default="",
                    help="Appended to the output filenames (immuannot_calls<SUFFIX>.tsv, "
                         "immuannot_timing<SUFFIX>.tsv) instead of writing the canonical shared "
                         "files. REQUIRED whenever multiple invocations of this script run "
                         "CONCURRENTLY as separate OS processes (e.g. via `xargs -P`, as the core-"
                         "scaling/production orchestrators both do) -- write_incremental() does a "
                         "non-atomic read-modify-write on the output TSV (read whole file, drop "
                         "this person's old rows, append, write whole file back), so two processes "
                         "finishing within moments of each other silently lose whichever one wrote "
                         "first. A single invocation processing many person_ids in its own for-loop "
                         "(the normal serial usage) has no race and does not need this -- only "
                         "concurrent SEPARATE invocations do. Give each concurrent worker a unique "
                         "suffix (e.g. the person_id or a worker index); merge the fragments "
                         "afterward (pandas.concat over immuannot_calls.*.tsv works, since schemas "
                         "match). The --force skip-check below still reads the CANONICAL "
                         "(unsuffixed) file, so resumability against a prior real run is unaffected.")
    args = ap.parse_args()

    lr_path = os.path.join(args.mount, LR_MANIFEST)
    if not os.path.exists(lr_path):
        die(f"manifest not found: {lr_path} -- is the gcsfuse mount up? "
            f"(ENVIRONMENT.md quirk #11/#14: remount and `ls`-verify first.)")
    immuannot_script = find_immuannot_script(args.immuannot_dir)
    if immuannot_script is None:
        die(f"immuannot.sh not found anywhere under {args.immuannot_dir} -- run setup_immuannot.sh first.")
    print(f"Using immuannot.sh at: {immuannot_script}", file=sys.stderr)

    print(f"Reading LR manifest: {lr_path}", file=sys.stderr)
    lr = pd.read_csv(lr_path, sep="\t", dtype=str)
    if "research_id" not in lr.columns:
        die(f"LR manifest missing research_id column. Actual columns: {list(lr.columns)}")

    calls_path = os.path.join(args.outroot, "immuannot_calls.tsv")
    already_done = set()
    if os.path.exists(calls_path) and not args.force:
        already_done = set(pd.read_csv(calls_path, sep="\t", dtype=str)["person_id"].astype(str))
        print(f"{len(already_done)} people already have calls in {calls_path} -- will be skipped "
              f"(pass --force to redo them). This is what makes killing and re-launching a 60-"
              f"person run with the exact same command safe.", file=sys.stderr)

    run_t0 = time.perf_counter()
    all_gene_rows, all_timing_rows = [], []
    n_skipped_done = 0
    for pid in args.person_ids:
        if str(pid) in already_done:
            print(f"\n=== Person {pid}: SKIP -- already in {calls_path} (--force to redo) ===",
                  file=sys.stderr)
            n_skipped_done += 1
            continue
        gene_rows, timing_rows = process_person(pid, lr, args, immuannot_script)
        all_gene_rows.extend(gene_rows)
        all_timing_rows.extend(timing_rows)

    n_ok = len({r["person_id"] for r in all_gene_rows})
    n_failed = len(args.person_ids) - n_ok - n_skipped_done
    run_total = time.perf_counter() - run_t0
    print(f"\n=== Run summary: {n_ok} produced calls this invocation, {n_skipped_done} already "
          f"done (skipped), {n_failed} failed/skipped this time, {run_total/60:.1f} min elapsed "
          f"({run_total/60/max(n_ok,1):.1f} min/person average, this invocation only) ===",
          file=sys.stderr)

    if all_timing_rows:
        timing_path = os.path.join(args.outroot, f"immuannot_timing{args.out_suffix}.tsv")
        write_incremental(pd.DataFrame(all_timing_rows), timing_path, ["person_id"])
        print(f"Timing log: {timing_path}", file=sys.stderr)

    if not all_gene_rows and n_skipped_done == 0:
        die("No people produced any calls -- check the SKIP reasons above.")
    if not all_gene_rows:
        print("Nothing new to write this invocation (everyone requested was already done).",
              file=sys.stderr)
        return

    out_path = os.path.join(args.outroot, f"immuannot_calls{args.out_suffix}.tsv")
    write_incremental(pd.DataFrame(all_gene_rows), out_path, ["person_id"])
    print(f"Calls: {out_path}", file=sys.stderr)
    print("Aggregate-only note doesn't apply here the way it does for comparison_log.csv -- these "
          "ARE per-person real allele calls (participant data). Keep on the VM; do not commit or "
          "download raw calls (same rule as SMOKE_TEST_PICKS.local.md / comparison_log.csv).",
          file=sys.stderr)


if __name__ == "__main__":
    main()
