#!/usr/bin/env python3
"""Split depth-4 ("beyond_cds") novelty into distinct non-coding sequences, using the cs difference
string between each person's contig and the call's template genomic allele.

Spec: sprints/S01_novelty_qc_coverage/WS1_novelty_by_field.md, section "script 25".

## Why
03_novel_alleles.py keys novel clusters on (gene, sha1 of observed CDS). A depth-4 call
(`HLA-A*01:01:01:new`) has a CDS identical to a known allele, so every depth-4 haplotype with that
CDS lands in ONE cluster, however different its introns/UTRs are (the largest pools 7,713 people).
Those `n_persons` are not the recurrence of one sequence. This script recovers the actual
non-coding difference of every depth-4 call from `mm2.ipd.gen.paf.gz`, which holds EVERY IPD
genomic allele aligned to the contig (`minimap2 -cx asm5 --cs --end-bonus=10`, query = IPD
genomic allele, target = contig; reference/IMMUANNOT_GTF_SPEC.md part C/E).

## cs convention (derivation)
minimap2's cs tag describes the alignment relative to the REFERENCE, and in minimap2 the reference
is always the TARGET. The minimap2 man page defines: `:N` identical run; `*ab` substitution
"ref a to query b"; `+seq` "insertion to the reference" (bases present in the query only);
`-seq` "deletion from the reference" (bases present in the target only). Here target = person's
contig, query = the known IPD allele. Therefore:

  `:N`    N identical bases                                  consumes allele + contig
  `*ab`   a = CONTIG base, b = ALLELE base                   -> SNV  allele b > contig a
  `+seq`  seq is in the allele, absent from the contig       -> CONTIG_DEL (allele coords exist)
  `-seq`  seq is in the contig, absent from the allele       -> CONTIG_INS (between allele bases)

This is the same orientation 21_hla_manhattan.py uses (`walk_cs_canonical`), which was checked
against real PAFs on the VM. Events are written as differences of the PERSON relative to the KNOWN
allele, in the known allele's own 1-based coordinates:
  SNV         `<pos><allele_base>><contig_base>`       e.g. `812A>G`
  CONTIG_DEL  `<start>_<end>del<allele_bases>`         (`<pos>del<b>` for 1 bp)
  CONTIG_INS  `<p>_<p+1>ins<contig_bases>`             inserted between allele bases p and p+1

## Minus strand
For a `-` row minimap2 aligned the reverse complement of the query. The cs string runs along the
target's forward strand and shows reverse-complemented query bases. So the k-th query base consumed
is original allele position `qend-1-k` (0-based). A `+seq` of length n starting after k consumed
bases covers allele positions `qend-k-n .. qend-k-1` and its sequence is `revcomp(seq)`. A `-seq`
after k consumed bases lies between allele positions `qend-k-1` and `qend-k`. SNV bases are
complemented. After this, the same biological difference gives the same signature on either
strand, with one exception: an indel inside a repeat. minimap2 places gaps in TARGET orientation,
so on the minus strand the same indel can land at a different offset inside a homopolymer or
tandem repeat. When the allele genomic sequence is available (see below), every indel is
LEFT-NORMALIZED in allele coordinates (VCF-style), which removes this. Without the sequence,
indel positions inside repeats are not normalized (`context_mode=content_only`), and the per-cluster
strand counts are reported so the effect can be seen.

## Homopolymer flag
- With the allele sequence (`context_mode=sequence`): an indel is a homopolymer indel when its bases
  are all one base b AND the allele carries a run of b of length >= `--homopolymer-min-run`
  (default 3) directly next to the event (bases left of the event + bases right of it, excluding
  the deleted bases themselves). That flank length is reported per event (`event_flank_runs`).
- Without it (`context_mode=content_only`, `homopolymer_context_checked=False`): the flag only
  requires the indel bases to be one repeated base. Every 1-bp indel qualifies, so this OVER-calls.

The allele genomic sequence is not in any Table. It is read from `<--refdata>/gen.fa.gz`
(default ~/tools/Immuannot_refdata; verified on the VM 2026-09-16: the IPD genomic alleles used as
PAF queries, 22,667 records, header `>HLA-A*01:01:01:01 HLA00001 ... 3503bp`), loading only the
alleles this run needs. If gen.fa.gz is absent, any FASTA outside `CDSseq/` whose first header
looks like an allele name is used; `--allele-fasta` names files explicitly. A sequence is only used when its length equals the PAF `qlen` AND every cs
event agrees with it (`seq_consistent`); otherwise that call falls back to content_only.
Event regions (CDS / UTR / intron) come from `alleles.csv.gz` (`--allele-annotation`, the same
file 21 uses, parsed generically as key=value), using its `CDS=`/`UTR=`/`geneRange=` fields; NA if
the file or those fields are missing (the VM check saw `geneRange=`; CDS/UTR presence unconfirmed). Depth-4 calls should have few CDS events because the observed CDS
equals a known CDS, though not necessarily the template's.

## PAF row selection (per call)
Rows whose query name equals the call's `template_allele` (compared after upper-casing and
dropping any `HLA-` prefix), whose target equals the call's contig, and whose target span
overlaps `gene_start..gene_end` (Table 1: 1-based inclusive; PAF: 0-based half-open). Among
them the longest query span wins (ties: lower NM). If there is no such row, fall back to rows
whose allele shares the call's first 3 fields (`HLA-A*01:01:01:`). Among those, the lowest NM wins,
then the longest span. `selection_path` records `template` / `fallback_3field` / `no_row`. The
signature is keyed on the allele actually used.

## Mapping a call to its old novel_id
novel_id = `<gene>_nov_<sha1(observed CDS)[:8]>` (03). The observed CDS comes from the hap's
`cds.fa.gz`, read with 03's own `parse_cds_fasta` (imported with importlib), using 03's exact
ambiguity rule: skip the call if Table 1 has more than one row at (person, hap, contig, gene), or
if cds.fa.gz has 0 or more than 1 record for `<contig>_<gene>`. 03's `match_novel_rows` is NOT
called. Its output rows drop `contig`, so a gene seen on two contigs of one hap could not be joined
back to the right call, and it would re-read every novel call, not only depth-4. The hash comes out
identical because the same parser and the same upper-casing are used. Alternative: `--novel-matches`
TSV with columns person_id, hap, contig, gene, cds_seq_sha1, which skips reading cds.fa.gz.

## Outputs
Committed, aggregate only (`--out-dir`, default ~/results/25_noncoding_novelty_paf). Every person
or haplotype count from 1 to 19 is written as `<20`:
  noncoding_signature_clusters.tsv  only clusters with >= 20 unrelated carriers
  beyond_cds_pooling.tsv            only former clusters with >= 20 persons
  pooling_distribution.tsv          former clusters binned by number of distinct signatures
  homopolymer_fraction.tsv          gene x ancestry scheme x ancestry x unrelated_only
  nonhomopolymer_sfs.tsv            number of signatures per carrier-count bin
  gene_summary.tsv, summary.json, noncoding_novelty_report.md, fig1-3 PNGs
VM-local (`--local-dir`, default ~/pipeline_outputs/s01_local/), never committed:
  25_noncoding_calls.tsv (person level), 25_noncoding_signature_clusters_full.tsv,
  25_beyond_cds_pooling_full.tsv, 25_unrelated_set.tsv, plus unsuppressed copies of the
  aggregate tables, and resumable checkpoints under 25_checkpoints/.

Usage (VM):
    pixi run -e spechla -- python3 scripts/hla_popgen/25_noncoding_novelty_paf.py --limit 200
    pixi run -e spechla -- python3 scripts/hla_popgen/25_noncoding_novelty_paf.py
    # all 65 genes instead of the classical 8:
    pixi run -e spechla -- python3 scripts/hla_popgen/25_noncoding_novelty_paf.py --genes all
"""
import argparse
import concurrent.futures
import gzip
import glob
import hashlib
import importlib.util
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict
from functools import partial

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import _viz_common as vc  # noqa: E402


def _load_module(filename, modname):
    if modname in sys.modules:
        return sys.modules[modname]
    spec = importlib.util.spec_from_file_location(modname, os.path.join(HERE, filename))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


# Clean pure helpers, reused as-is (both modules only define functions/constants at import).
_m21 = _load_module("21_hla_manhattan.py", "hla_manhattan")
_m03 = _load_module("03_novel_alleles.py", "novel_alleles_03")
_m11 = _load_module("11_relatedness_cohort_overlap.py", "relatedness_11")
get_tag = _m21.get_tag
parse_ranges = _m21.parse_ranges
parse_cds_fasta = _m03.parse_cds_fasta

DATA_ROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_TABLE1 = os.path.join(DATA_ROOT, "hla_calls_rich.tsv")
DEFAULT_COHORT = os.path.join(DATA_ROOT, "cohort_membership.tsv")
DEFAULT_OUTROOT = os.path.join(DATA_ROOT, "people")
DEFAULT_OUT_DIR = os.path.expanduser("~/results/25_noncoding_novelty_paf")
DEFAULT_LOCAL_DIR = os.path.join(DATA_ROOT, "s01_local")
DEFAULT_REFDATA = os.path.expanduser("~/tools/Immuannot_refdata")
DEFAULT_NOVEL_ALLELES = os.path.join(vc.DEFAULT_REPORT_ROOT, "novel_alleles.tsv")
# Controlled data is auto-mounted here (verified on the VM 2026-09-16). 11's default
# (~/mnt/aou-controlled) is stale and hangs processes, so it is deliberately NOT reused.
DEFAULT_RELATEDNESS = os.path.expanduser(
    "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/"
    "samples_relatedness.tsv")

KIN_THRESHOLD = 0.0442          # third degree or closer (KING; same as 11 / script 24)
STRICT_ANCESTRY_MIN = 0.9
MIN_COMMIT = 20                 # counts 1..19 are written as "<20" in committed outputs
HOMOPOLYMER_MIN_RUN = 3
BATCH_PERSONS = 500
MIN_QCOV = 0.98                 # calls covering less of the allele than this are excluded
                                 # from clustering (a partial alignment yields a signature that
                                 # is a subset of the truth and falsely merges with full ones)
CONTEXT_SEQUENCE_MIN_FRAC = 0.99  # below this, the homopolymer flag is not trustworthy enough
                                   # to publish (see --allow-content-only)

KINDS = ("SNV", "CONTIG_DEL", "CONTIG_INS")
SIG_CLASSES = ["homopolymer_only", "indel_nonhomopolymer", "snv_plus_indel", "snv_only",
               "no_difference"]
SIG_COLORS = {  # Okabe-Ito
    "homopolymer_only": "#E69F00", "indel_nonhomopolymer": "#CC79A7",
    "snv_plus_indel": "#009E73", "snv_only": "#0072B2", "no_difference": "#999999",
}
REGION_CLASSES = ["intron_only", "involves_utr", "touches_cds", "unannotated", "other",
                  "no_difference"]
SFS_BINS = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 10), (11, 20), (21, 50), (51, 100),
            (101, 1000), (1001, 10 ** 9)]

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CS_RE = re.compile(r":(\d+)|\*([A-Za-z])([A-Za-z])|\+([A-Za-z]+)|-([A-Za-z]+)|=([A-Za-z]+)"
                   r"|(~[A-Za-z]{2}\d+[A-Za-z]{2})")


class CsError(ValueError):
    pass


def revcomp(s):
    return s.translate(_COMP)[::-1]


def norm_allele(name):
    """'HLA-A*01:01:01:01' / 'hla-A*01:01:01:01' / 'A*01:01:01:01' -> 'A*01:01:01:01'."""
    if not isinstance(name, str):
        return None
    s = name.strip().upper()
    return s[4:] if s.startswith("HLA-") else s


def three_field_prefix(consensus):
    """'HLA-A*01:01:01:new' -> 'A*01:01:01:' (normalized), or None if < 3 resolved fields."""
    if not isinstance(consensus, str) or "*" not in consensus:
        return None
    gene = consensus.split("*", 1)[0]
    fields = vc.parse_allele_fields(consensus)
    if len(fields) < 3:
        return None
    return norm_allele(f"{gene}*{':'.join(fields[:3])}:")


def key_prefix3(key):
    """Normalized allele key -> its 3-field prefix with trailing ':' (or None)."""
    if not key or "*" not in key:
        return None
    g, rest = key.split("*", 1)
    parts = rest.split(":")
    if len(parts) < 4:
        return None
    return f"{g}*{':'.join(parts[:3])}:"


# ---------------------------------------------------------------------------
# cs parsing -> events in allele (query) coordinates, strand-normalized
# ---------------------------------------------------------------------------
def parse_cs_events(cs, qstart, qend, strand):
    """Returns a list of events (pos0, kind, allele_seq, contig_seq), sorted, in the allele's
    forward orientation. See module docstring for the derivation.
      SNV:        pos0 = allele position; allele_seq = allele base; contig_seq = contig base
      CONTIG_DEL: pos0 = first deleted allele position; allele_seq = deleted bases
      CONTIG_INS: pos0 = allele position the insertion sits BEFORE; contig_seq = inserted bases
    Raises CsError on unparseable strings or when consumed query length != qend - qstart."""
    if strand not in ("+", "-"):
        raise CsError(f"bad strand {strand!r}")
    events = []
    k = 0
    pos = 0
    minus = strand == "-"
    for m in CS_RE.finditer(cs):
        if m.start() != pos:
            raise CsError(f"unparseable cs at offset {pos}")
        pos = m.end()
        n_match, t_base, q_base, ins_q, del_t, eq_seq, intron = m.groups()
        if n_match is not None:
            k += int(n_match)
        elif eq_seq is not None:
            k += len(eq_seq)
        elif t_base is not None:
            a, c = q_base.upper(), t_base.upper()
            if minus:
                p, a, c = qend - 1 - k, a.translate(_COMP), c.translate(_COMP)
            else:
                p = qstart + k
            events.append((p, "SNV", a, c))
            k += 1
        elif ins_q is not None:          # '+': bases in query (allele) only -> contig deletion
            s = ins_q.upper()
            n = len(s)
            if minus:
                p, s = qend - k - n, revcomp(s)
            else:
                p = qstart + k
            events.append((p, "CONTIG_DEL", s, ""))
            k += n
        elif del_t is not None:          # '-': bases in target (contig) only -> contig insertion
            s = del_t.upper()
            if minus:
                p, s = qend - k, revcomp(s)
            else:
                p = qstart + k
            events.append((p, "CONTIG_INS", "", s))
        else:
            raise CsError("intron op '~' not expected in asm5 PAF")
    if pos != len(cs):
        raise CsError(f"unparseable cs tail at offset {pos}")
    if k != qend - qstart:
        raise CsError(f"cs consumes {k} query bases but qend-qstart={qend - qstart}")
    return sorted(events, key=lambda e: (e[0], KINDS.index(e[1]), e[2], e[3]))


def seq_consistent(events, allele):
    for p, kind, a, _c in events:
        if kind == "SNV" and allele[p:p + 1] != a:
            return False
        if kind == "CONTIG_DEL" and allele[p:p + len(a)] != a:
            return False
        if kind == "CONTIG_INS" and not (0 <= p <= len(allele)):
            return False
    return True


def left_normalize(events, allele):
    """VCF-style left shift of every indel in allele coordinates."""
    out = []
    for p, kind, a, c in events:
        if kind == "CONTIG_DEL":
            n = len(a)
            while p > 0 and allele[p - 1] == allele[p + n - 1]:
                p -= 1
            a = allele[p:p + n]
        elif kind == "CONTIG_INS":
            while p > 0 and allele[p - 1] == c[-1]:
                c = c[-1] + c[:-1]
                p -= 1
        out.append((p, kind, a, c))
    return sorted(out, key=lambda e: (e[0], KINDS.index(e[1]), e[2], e[3]))


def _run_left(seq, i, b):
    n = 0
    j = i - 1
    while j >= 0 and seq[j] == b:
        n += 1
        j -= 1
    return n


def _run_right(seq, i, b):
    n = 0
    j = i
    while j < len(seq) and seq[j] == b:
        n += 1
        j += 1
    return n


def homopolymer_info(event, allele=None, min_run=HOMOPOLYMER_MIN_RUN):
    """-> (is_homopolymer_indel, flank_run_len or None). SNVs -> (False, None)."""
    p, kind, a, c = event
    if kind == "SNV":
        return False, None
    s = a if kind == "CONTIG_DEL" else c
    single = len(s) >= 1 and len(set(s)) == 1
    if allele is None:
        return single, None
    b = s[0]
    end = p + len(s) if kind == "CONTIG_DEL" else p
    flank = _run_left(allele, p, b) + _run_right(allele, end, b)
    return bool(single and flank >= min_run), flank


def event_token(event):
    p, kind, a, c = event
    if kind == "SNV":
        return f"{p + 1}{a}>{c}"
    if kind == "CONTIG_DEL":
        return f"{p + 1}del{a}" if len(a) == 1 else f"{p + 1}_{p + len(a)}del{a}"
    return f"{p}_{p + 1}ins{c}"


def parse_exons(spec):
    """'1:1..373,2:504..773' -> [(1, 1, 373), (2, 504, 773)]."""
    out = []
    for part in (spec or "").split(","):
        num, sep, rng = part.strip().partition(":")
        if not sep:
            continue
        r = parse_ranges(rng)
        try:
            n = int(num)
        except ValueError:
            continue
        out.extend((n, a, b) for a, b in r)
    return sorted(out, key=lambda t: t[1])


def region_of(pos0, ann):
    """Allele-coordinate position -> UTR5 / UTR3 / exon<N> (CDS) / exon<N>_noncds / intron<N> /
    outside / unannotated. Uses alleles.csv.gz CDS=, UTR=, exon=, geneRange= (1-based inclusive,
    allele coordinates; verified layout on the VM 2026-09-16)."""
    if not ann or not ann["cds"]:
        return "unannotated"
    x = pos0 + 1
    cds_lo = min(a for a, _ in ann["cds"])
    cds_hi = max(b for _, b in ann["cds"])
    for a, b in ann["utr"]:
        if a <= x <= b:
            return "UTR5" if x < cds_lo else ("UTR3" if x > cds_hi else "UTR")
    exons = ann["exons"]
    for n, a, b in exons:
        if a <= x <= b:
            in_cds = any(ca <= x <= cb for ca, cb in ann["cds"])
            return f"exon{n}" if in_cds else f"exon{n}_noncds"
    if any(ca <= x <= cb for ca, cb in ann["cds"]):
        return "exon?"
    for (n1, _a1, b1), (_n2, a2, _b2) in zip(exons, exons[1:]):
        if b1 < x < a2:
            return f"intron{n1}"
    gr = ann["gene_range"]
    if gr and any(a <= x <= b for a, b in gr):
        # inside the gene but no exon map: intron if between CDS blocks, else UTR side
        return "intron?" if cds_lo < x < cds_hi else ("UTR5" if x < cds_lo else "UTR3")
    return "outside"


def region_bucket(label):
    if label.startswith("exon"):
        return "cds"
    if label.startswith("UTR"):
        return "utr"
    if label.startswith("intron"):
        return "intron"
    return "other"   # outside / unannotated


def classify_regions(labels):
    """Per-call region class from its event labels."""
    if not labels:
        return "no_difference"
    buckets = {region_bucket(lab) for lab in labels}
    if "unannotated" in labels:
        return "unannotated"
    if "cds" in buckets:
        return "touches_cds"
    if "utr" in buckets:
        return "involves_utr"
    if buckets == {"intron"}:
        return "intron_only"
    return "other"


def classify_signature(n_snv, n_indel, n_hp):
    if n_snv == 0 and n_indel == 0:
        return "no_difference"
    if n_snv == 0 and n_hp == n_indel:
        return "homopolymer_only"
    if n_snv == 0:
        return "indel_nonhomopolymer"
    if n_indel == 0:
        return "snv_only"
    return "snv_plus_indel"


def merge_events(event_lists):
    """Merge several already-parsed (allele-coordinate) event lists from chained PAF records into
    one, in the same sort order parse_cs_events produces."""
    merged = [e for lst in event_lists for e in lst]
    return sorted(merged, key=lambda e: (e[0], KINDS.index(e[1]), e[2], e[3]))


def build_signature(gene, used_allele, cs, qstart, qend, strand, allele_seq=None, ann=None,
                    min_run=HOMOPOLYMER_MIN_RUN):
    """Full per-call signature record for a single PAF record (raises CsError)."""
    events = parse_cs_events(cs, qstart, qend, strand)
    return build_signature_from_events(gene, used_allele, events, allele_seq, ann, min_run)


def build_signature_from_events(gene, used_allele, events, allele_seq=None, ann=None,
                                min_run=HOMOPOLYMER_MIN_RUN):
    """Same as build_signature, but takes already-parsed events (allele coordinates) so that
    multiple chained PAF records (see select_paf_group/check_collinear) can be merged into one
    signature before this is called."""
    ctx = None
    consistent = None
    if allele_seq is not None:
        consistent = seq_consistent(events, allele_seq)
        if consistent:
            ctx = allele_seq
            events = left_normalize(events, allele_seq)
    tokens, flanks, regions = [], [], []
    n_snv = n_indel = n_hp = 0
    for ev in events:
        tokens.append(event_token(ev))
        regions.append(region_of(ev[0], ann))
        if ev[1] == "SNV":
            n_snv += 1
            flanks.append("")
        else:
            n_indel += 1
            hp, flank = homopolymer_info(ev, ctx, min_run)
            n_hp += int(hp)
            flanks.append("" if flank is None else str(flank))
    body = ";".join(tokens) if tokens else "NO_DIFF"
    sig = f"{used_allele}|{body}"
    return {
        "events": events,
        "signature": sig,
        "signature_id": f"{gene}_sig_{hashlib.sha1(sig.encode()).hexdigest()[:10]}",
        "event_tokens": body,
        "event_flank_runs": ",".join(flanks),
        "event_regions": ",".join(regions),
        "region_class": classify_regions(regions),
        "n_events": len(events), "n_snv": n_snv, "n_indel": n_indel,
        "n_homopolymer_indel": n_hp,
        "n_events_intron": sum(region_bucket(r) == "intron" for r in regions),
        "n_events_utr": sum(region_bucket(r) == "utr" for r in regions),
        "n_events_cds": sum(region_bucket(r) == "cds" for r in regions),
        "n_events_unannotated": sum(region_bucket(r) == "other" for r in regions),
        "signature_class": classify_signature(n_snv, n_indel, n_hp),
        "context_mode": "sequence" if ctx is not None else "content_only",
        "homopolymer_context_checked": ctx is not None,
        "seq_consistent": consistent,
    }


# ---------------------------------------------------------------------------
# PAF reading and row selection
# ---------------------------------------------------------------------------
def parse_paf_line(line):
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        return None
    nm = get_tag(f, "NM")
    return {
        "qname": f[0], "qkey": norm_allele(f[0]), "qlen": int(f[1]), "qstart": int(f[2]),
        "qend": int(f[3]), "strand": f[4], "tname": f[5], "tstart": int(f[7]),
        "tend": int(f[8]), "nm": int(nm) if nm is not None else None, "cs": get_tag(f, "cs"),
        "tp": get_tag(f, "tp"),
    }


def read_paf_candidates(paf_path, exact_keys, prefix_keys, contigs):
    """One pass; keep PRIMARY rows (tp:A:P, or no tp tag at all -- treated as primary) whose
    query is a wanted allele (exact or 3-field prefix) on a wanted contig. Secondary rows
    (tp:A:S) that would otherwise have matched are dropped and counted, since the raw minimap2
    PAF can let a shorter secondary alignment win on length if it is not filtered out first.
    -> (rows, n_secondary_filtered)."""
    rows = []
    n_secondary = 0
    genes = {k.split("*", 1)[0] for k in exact_keys | prefix_keys}
    # Cheap reject on the raw (unnormalized) query name before paying for .strip()/.upper()/etc:
    # PAF query names are not guaranteed to carry the 'HLA-' prefix or a fixed case, so we accept
    # any of the plausible raw spellings of "<gene>*" here; normalization still runs on survivors.
    raw_prefixes = tuple(sorted(
        {f"{g}*" for g in genes} | {f"HLA-{g}*" for g in genes} | {f"hla-{g}*" for g in genes}
        | {f"Hla-{g}*" for g in genes}))
    with gzip.open(paf_path, "rt") as fh:
        for line in fh:
            tab = line.find("\t")
            if tab < 0:
                continue
            head = line[:tab]
            if not head.startswith(raw_prefixes):
                continue
            key = norm_allele(head)
            if key.split("*", 1)[0] not in genes:
                continue
            if key not in exact_keys and key_prefix3(key) not in prefix_keys:
                continue
            r = parse_paf_line(line)
            if r is None or r["tname"] not in contigs:
                continue
            if r["tp"] == "S":
                n_secondary += 1
                continue
            rows.append(r)
    return rows, n_secondary


def check_collinear(group):
    """group: >=2 primary PAF records sharing (qname, tname, strand). True when they are
    non-overlapping in query coordinates and ordered consistently with strand in target
    coordinates -- i.e. minimap2 split one true alignment into pieces (typically across a large
    intronic indel) rather than reporting two genuinely different/conflicting hits."""
    segs = sorted(group, key=lambda r: r["qstart"])
    for a, b in zip(segs, segs[1:]):
        if a["qend"] > b["qstart"]:
            return False
        if a["strand"] == "+":
            if a["tend"] > b["tstart"]:
                return False
        else:
            if b["tend"] > a["tstart"]:
                return False
    return True


def select_paf_group(rows, contig, gene_start, gene_end, template_key, prefix_key):
    """rows: primary-only PAF records (see read_paf_candidates). -> (group, selection_path).
    `group` is every primary record sharing (qname, tname, strand) with the winning candidate:
    when minimap2 splits one true alignment into several primary records (e.g. across a large
    indel), they all share that triple and the caller (process_hap) must chain or flag them
    split_alignment. group is None (path 'no_row') if nothing matched."""
    lo, hi = int(gene_start) - 1, int(gene_end)  # 1-based inclusive -> 0-based half-open

    def on_target(r):
        return r["tname"] == contig and r["tstart"] < hi and r["tend"] > lo

    def grouped(cands):
        groups = defaultdict(list)
        for r in cands:
            groups[(r["qname"], r["tname"], r["strand"])].append(r)
        return list(groups.values())

    big = 10 ** 9
    cands = [r for r in rows if r["qkey"] == template_key and on_target(r)]
    if cands:
        def score(grp):
            span = sum(r["qend"] - r["qstart"] for r in grp)
            nms = [r["nm"] for r in grp if r["nm"] is not None]
            nm = sum(nms) if len(nms) == len(grp) else None
            return span, -(nm if nm is not None else big)
        best = max(grouped(cands), key=score)
        return best, "template"
    if prefix_key:
        cands = [r for r in rows if key_prefix3(r["qkey"]) == prefix_key and on_target(r)]
        if cands:
            def score2(grp):
                span = sum(r["qend"] - r["qstart"] for r in grp)
                nms = [r["nm"] for r in grp if r["nm"] is not None]
                nm = sum(nms) if len(nms) == len(grp) else None
                return (nm if nm is not None else big, -span, grp[0]["qkey"])
            best = min(grouped(cands), key=score2)
            return best, "fallback_3field"
    return None, "no_row"


def select_paf_row(rows, contig, gene_start, gene_end, template_key, prefix_key):
    """Back-compat convenience wrapper around select_paf_group returning a single representative
    row instead of the whole group (used by callers/tests that don't need chaining). New code
    that must detect multi-record chains/splits should call select_paf_group directly."""
    group, path = select_paf_group(rows, contig, gene_start, gene_end, template_key, prefix_key)
    return (group[0] if group else None), path


# ---------------------------------------------------------------------------
# Reference data (optional)
# ---------------------------------------------------------------------------
def _fasta_first_header(path):
    try:
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    return line[1:].split()[0]
                if line.strip():
                    return None
    except (OSError, UnicodeDecodeError, EOFError):
        return None
    return None


def find_allele_fastas(refdata):
    """`<refdata>/gen.fa.gz` if present (verified on the VM: the IPD genomic alleles used as PAF
    queries, header `>HLA-A*01:01:01:01 HLA00001 ... 3503bp`); otherwise glob for FASTA-like
    files outside CDSseq/ whose first header looks like an allele name."""
    if not refdata or not os.path.isdir(refdata):
        return []
    gen = os.path.join(refdata, "gen.fa.gz")
    if os.path.exists(gen):
        return [gen]
    pats = ("*.fa", "*.fa.gz", "*.fasta", "*.fasta.gz", "*.fna", "*.fna.gz")
    out = set()
    for p in pats:
        out.update(glob.glob(os.path.join(refdata, "**", p), recursive=True))
    keep = []
    for p in sorted(out):
        if os.sep + "CDSseq" + os.sep in p:
            continue
        h = _fasta_first_header(p)
        if h and "*" in h:
            keep.append(p)
    return keep


def load_allele_sequences(paths, wanted_exact, wanted_prefix):
    """{normalized allele key: SEQ} for wanted alleles only. The first file that has an allele wins."""
    seqs = {}

    def want(key):
        return key in wanted_exact or key_prefix3(key) in wanted_prefix

    for path in paths:
        op = gzip.open if path.endswith(".gz") else open
        key, chunks = None, []
        with op(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    if key is not None and key not in seqs:
                        seqs[key] = "".join(chunks).upper()
                    k = norm_allele(line[1:].split()[0])
                    key = k if want(k) else None
                    chunks = []
                elif key is not None:
                    chunks.append(line.strip())
        if key is not None and key not in seqs:
            seqs[key] = "".join(chunks).upper()
    return seqs


def load_annotations(path, wanted_exact, wanted_prefix):
    """alleles.csv.gz (tab-separated key=value; see 21.load_allele_annotation) -> {key: ann}."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with gzip.open(path, "rt") as fh:
        for line in fh:
            fields = {}
            for tok in line.rstrip("\n").split("\t"):
                k, _, v = tok.partition("=")
                fields[k.strip()] = v.strip()
            key = norm_allele(fields.get("allele"))
            if not key or not (key in wanted_exact or key_prefix3(key) in wanted_prefix):
                continue
            out[key] = {"cds": parse_ranges(fields.get("CDS", "")),
                        "utr": parse_ranges(fields.get("UTR", "")),
                        "exons": parse_exons(fields.get("exon", "")),
                        "gene_range": parse_ranges(fields.get("geneRange", ""))}
    return out


# ---------------------------------------------------------------------------
# Cohort: unrelated set, ancestry
# ---------------------------------------------------------------------------
def greedy_unrelated(pairs, people, kin_threshold=KIN_THRESHOLD):
    """Same semantics as script 24: among `people`, drop members of related pairs
    (kin >= threshold, both sides in `people`) until no pair remains, always removing the member
    with the most remaining relatives, ties broken by sorted id. Returns the set kept."""
    people = set(people)
    adj = defaultdict(set)
    for i, j, k in pairs:
        if i == j or k < kin_threshold or i not in people or j not in people:
            continue
        adj[i].add(j)
        adj[j].add(i)
    removed = set()
    while True:
        live = [(len(nb), n) for n, nb in adj.items() if nb]
        if not live:
            break
        maxdeg = max(d for d, _ in live)
        victim = sorted(n for d, n in live if d == maxdeg)[0]
        for nb in adj[victim]:
            adj[nb].discard(victim)
        adj[victim] = set()
        removed.add(victim)
    return people - removed


def strict_ancestry(cohort, threshold=STRICT_ANCESTRY_MIN):
    """Label if that ancestry's admixture proportion (p_<anc>) >= threshold, else NA."""
    out = pd.Series(pd.NA, index=cohort.index, dtype=object)
    for anc in vc.ANCESTRY_ORDER:
        col = f"p_{anc.lower()}"
        if col not in cohort.columns:
            continue
        hit = (cohort["ancestry_pred"] == anc) & (pd.to_numeric(cohort[col], errors="coerce")
                                                   >= threshold)
        out[hit.fillna(False).astype(bool)] = anc
    return out


# ---------------------------------------------------------------------------
# Per-haplotype worker
# ---------------------------------------------------------------------------
def process_hap(pid, hap, calls, cfg):
    """calls: list of dicts (Table 1 depth-4 rows of this person/hap). -> list of result dicts.

    PAF record selection (fix for the CRITICAL "one PAF record per call" issue): secondary rows
    are dropped in read_paf_candidates; the remaining primary rows are grouped by
    (qname, tname, strand) in select_paf_group. A group of size 1 is the common case. A group of
    >1 means minimap2 emitted several primary records for what should be one alignment (typically
    split across a large indel): if they are collinear and non-overlapping in query coordinates
    (check_collinear) they are CHAINED -- their cs strings are parsed independently and the
    resulting events merged, and qcov is the union of their covered query length over qlen; if
    they conflict, the call is marked status='split_alignment' and excluded from clustering (but
    counted). A --min-qcov floor then excludes any surviving call whose covered fraction is too
    low (status='low_qcov'), since a partial alignment yields a signature that is a strict subset
    of the truth and would otherwise silently merge with full-length calls."""
    hap_dir = os.path.join(cfg["outroot"], pid, "immuannot_output", hap)
    paf = os.path.join(hap_dir, "mm2.ipd.gen.paf.gz")
    exact = {c["template_key"] for c in calls if c["template_key"]}
    prefix = {c["prefix_key"] for c in calls if c["prefix_key"]}
    contigs = {c["contig"] for c in calls}
    rows, paf_status, n_secondary_filtered = [], "ok", 0
    if not os.path.exists(paf):
        paf_status = "paf_missing"
    else:
        try:
            rows, n_secondary_filtered = read_paf_candidates(paf, exact, prefix, contigs)
        except (OSError, EOFError, ValueError) as e:
            paf_status = f"paf_error:{type(e).__name__}"

    cds_index = None
    out = []
    for c in calls:
        res = {k: c[k] for k in ("person_id", "hap", "contig", "gene", "copy_index",
                                 "consensus", "template_allele", "gene_start", "gene_end")}
        res["n_secondary_filtered_hap"] = n_secondary_filtered
        # --- old novel_id via observed CDS hash ---
        sha1, cds_status = None, "ok"
        if cfg["novel_matches"] is not None:
            sha1 = cfg["novel_matches"].get((pid, hap, c["contig"], c["gene"]))
            cds_status = "ok" if sha1 else "not_in_novel_matches"
        elif cfg["group_sizes"].get((pid, hap, c["contig"], c["gene"]), 1) > 1:
            cds_status = "ambiguous_copy"
        else:
            if cds_index is None:
                try:
                    cds_index = parse_cds_fasta(os.path.join(hap_dir, "cds.fa.gz"))
                except (OSError, EOFError, gzip.BadGzipFile, UnicodeDecodeError):
                    cds_index = {}
                    cds_status = "cds_fasta_corrupt"
            if cds_status == "ok":
                cands = cds_index.get(f"{c['contig']}_{c['gene']}", [])
                if len(cands) == 0 or not cands[0][1]:
                    cds_status = "missing_cds"
                elif len(cands) > 1:
                    cds_status = "ambiguous_copy"
                else:
                    sha1 = hashlib.sha1(cands[0][1].encode("ascii")).hexdigest()
        res["cds_status"] = cds_status
        res["cds_seq_sha1"] = sha1
        res["novel_id"] = f"{c['gene']}_nov_{sha1[:8]}" if sha1 else None

        # --- PAF row(s) + signature ---
        res.update(used_allele=None, selection_path="no_row", chain_status="n/a",
                   paf_strand=None, qlen=np.nan, qstart=np.nan, qend=np.nan, qcov=np.nan,
                   nm=np.nan)
        if paf_status != "ok":
            res["status"] = paf_status
            out.append(res)
            continue
        group, path = select_paf_group(rows, c["contig"], c["gene_start"], c["gene_end"],
                                       c["template_key"], c["prefix_key"])
        res["selection_path"] = path
        if group is None:
            res["status"] = "no_row"
            out.append(res)
            continue
        qname, strand, qlen = group[0]["qname"], group[0]["strand"], group[0]["qlen"]
        covered = sum(r["qend"] - r["qstart"] for r in group)
        res.update(used_allele=qname, paf_strand=strand, qlen=qlen,
                   qstart=min(r["qstart"] for r in group), qend=max(r["qend"] for r in group),
                   qcov=(covered / qlen) if qlen else np.nan)
        nms = [r["nm"] for r in group if r["nm"] is not None]
        res["nm"] = sum(nms) if len(nms) == len(group) else np.nan
        if len(group) > 1 and not check_collinear(group):
            res["status"] = "split_alignment"
            res["chain_status"] = "split_alignment"
            out.append(res)
            continue
        res["chain_status"] = "single" if len(group) == 1 else "multi_record_chained"
        if any(r["cs"] is None for r in group):
            res["status"] = "cs_missing"
            out.append(res)
            continue
        try:
            events = merge_events([parse_cs_events(r["cs"], r["qstart"], r["qend"], r["strand"])
                                   for r in group])
        except CsError as e:
            res["status"] = f"cs_error:{e}"
            out.append(res)
            continue
        qkey = group[0]["qkey"]
        seq = cfg["allele_seqs"].get(qkey)
        if seq is not None and len(seq) != qlen:
            seq = None
            res["allele_seq_len_mismatch"] = True
        sig = build_signature_from_events(c["gene"], qname, events, seq,
                                          cfg["annotations"].get(qkey), cfg["min_run"])
        res.update(sig)
        # The floor excludes the call from clustering (its signature is a subset of the truth
        # and would otherwise merge with full-length calls), but the signature is still computed
        # and kept so low-coverage 'no_difference' calls can be cross-tabulated against qcov
        # (see no_difference_by_qcov_decile) -- that cross-tab is exactly what motivates the
        # floor.
        res["status"] = "low_qcov" if res["qcov"] < cfg["min_qcov"] else "ok"
        out.append(res)
    return out


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------
def suppress(v, min_n=MIN_COMMIT):
    if v is None or v is pd.NA or (isinstance(v, float) and np.isnan(v)):
        return "NA"
    try:
        iv = int(v)
    except (TypeError, ValueError):
        return v
    return f"<{min_n}" if 0 < iv < min_n else iv


def suppress_frame(df, cols, min_n=MIN_COMMIT):
    df = df.copy()
    for c in cols:
        if c in df.columns:
            df[c] = df[c].map(lambda v: suppress(v, min_n)).astype(object)
    return df


def _majority(series):
    """-> (mode value, n_disagreeing) for a categorical Series that is expected -- but not
    guaranteed -- to be constant within a signature cluster."""
    counts = series.value_counts()
    mode = counts.idxmax()
    return mode, int(len(series) - counts.max())


def build_clusters(ok):
    """ok: per-call rows with status=='ok'. -> per-signature cluster table (full, unsuppressed).
    signature_class/region_class are supposed to be pure functions of the signature, so every
    member of a cluster should agree; they are still taken by majority vote with a recorded
    disagreement count rather than an arbitrary row's value, as a safety net."""
    recs = []
    for sid, g in ok.groupby("signature_id", sort=False):
        first = g.iloc[0]
        persons = set(g["person_id"])
        unrel = g[g["unrelated"].fillna(False) == True]  # noqa: E712
        sig_class, n_sig_disagree = _majority(g["signature_class"])
        region_class, n_region_disagree = _majority(g["region_class"])
        rec = {
            "signature_id": sid, "gene": first["gene"], "template_allele": first["used_allele"],
            "signature_class": sig_class, "n_signature_class_disagree": n_sig_disagree,
            "n_events": int(first["n_events"]), "n_snv": int(first["n_snv"]),
            "n_indel": int(first["n_indel"]),
            "n_homopolymer_indel": int(first["n_homopolymer_indel"]),
            "all_homopolymer": sig_class == "homopolymer_only",
            "homopolymer_context_checked": bool(g["homopolymer_context_checked"].all()),
            "region_class": region_class, "n_region_class_disagree": n_region_disagree,
            "n_events_intron": int(first["n_events_intron"]),
            "n_events_utr": int(first["n_events_utr"]),
            "n_events_cds": int(first["n_events_cds"]),
            "n_events_unannotated": int(first["n_events_unannotated"]),
            "events": first["event_tokens"], "event_regions": first["event_regions"],
            "event_flank_runs": first["event_flank_runs"],
            "n_haplotypes": len(g), "n_persons": len(persons),
            "n_persons_unrelated": unrel["person_id"].nunique(),
            "n_haplotypes_unrelated": len(unrel),
            "n_plus_strand": int((g["paf_strand"] == "+").sum()),
            "n_minus_strand": int((g["paf_strand"] == "-").sum()),
            "n_fallback_rows": int((g["selection_path"] == "fallback_3field").sum()),
            "n_partial_alignment": int((g["qcov"] < 1.0).sum()),
            "n_multi_record_chained": int((g["chain_status"] == "multi_record_chained").sum()),
            "n_former_novel_ids": g["novel_id"].dropna().nunique(),
        }
        per_person = g.drop_duplicates("person_id")
        for scheme, col in (("pred", "ancestry_pred"), ("strict", "strict_ancestry")):
            vcounts = per_person[col].value_counts()
            for anc in vc.ANCESTRY_ORDER:
                rec[f"n_persons_{scheme}_{anc}"] = int(vcounts.get(anc, 0))
        recs.append(rec)
    df = pd.DataFrame(recs)
    if len(df):
        df = df.sort_values(["gene", "n_persons_unrelated", "n_persons"],
                            ascending=[True, False, False]).reset_index(drop=True)
    return df


def build_pooling(ok, novel_tsv):
    """Per former CDS-hash cluster (novel_id), how many distinct non-coding signatures it pooled.
    n_distinct_signatures counts every signature seen, but the signature includes the template
    allele (see module docstring), so two identical contigs aligned against different templates
    count as two signatures even though the sequence pooling is confounded by template choice,
    not by real biology. n_distinct_signatures_modal_template restricts to calls that used the
    cluster's single most common template, which is the honest headline number for "how many
    distinct sequences did this one old cluster actually pool" -- report and figure both use it,
    with the raw number kept alongside for transparency."""
    recs = []
    sub = ok[ok["novel_id"].notna()]
    for nid, g in sub.groupby("novel_id"):
        sig_counts = g["signature_id"].value_counts()
        nonhp = g[~g["signature_class"].isin(["homopolymer_only", "no_difference"])]
        modal_template = g["used_allele"].mode()
        modal_template = modal_template.iloc[0] if len(modal_template) else None
        g_modal = g[g["used_allele"] == modal_template] if modal_template is not None else g.iloc[0:0]
        rec = {
            "novel_id": nid, "gene": g["gene"].iloc[0],
            "n_haplotypes": len(g), "n_persons": g["person_id"].nunique(),
            "n_persons_unrelated":
                g.loc[g["unrelated"].fillna(False) == True, "person_id"].nunique(),  # noqa
            "modal_template": modal_template,
            "n_distinct_signatures_modal_template": int(g_modal["signature_id"].nunique()),
            "n_distinct_signatures": int(len(sig_counts)),
            "n_distinct_nonhomopolymer_signatures": int(nonhp["signature_id"].nunique()),
            "n_distinct_templates": int(g["used_allele"].nunique()),
            "top_signature_hap_fraction": float(sig_counts.iloc[0] / len(g)),
            "n_haplotypes_homopolymer_only": int((g["signature_class"] == "homopolymer_only").sum()),
            "in_novel_alleles_tsv": nid in novel_tsv.index if novel_tsv is not None else "NA",
        }
        if novel_tsv is not None and nid in novel_tsv.index:
            t = novel_tsv.loc[nid]
            rec["old_n_persons"] = t.get("n_persons")
            rec["old_n_haplotypes"] = t.get("n_haplotypes")
            rec["old_novelty_class"] = t.get("novelty_class")
            rec["old_confidence_tier"] = t.get("confidence_tier")
        recs.append(rec)
    df = pd.DataFrame(recs)
    if len(df):
        df = df.sort_values("n_persons", ascending=False).reset_index(drop=True)
    return df


def pooling_distribution(pool, col="n_distinct_signatures_modal_template"):
    """Distribution of former clusters by number of distinct signatures. Uses the modal-template
    count by default (the headline, honest number -- see build_pooling); pass
    col='n_distinct_signatures' for the raw, template-confounded distribution."""
    bins = [(1, 1), (2, 2), (3, 5), (6, 10), (11, 50), (51, 100), (101, 1000), (1001, 10 ** 9)]
    recs = []
    for lo, hi in bins:
        label = str(lo) if lo == hi else (f">{lo - 1}" if hi >= 10 ** 9 else f"{lo}-{hi}")
        n = int(((pool[col] >= lo) & (pool[col] <= hi)).sum()) if len(pool) else 0
        recs.append({"n_distinct_signatures_bin": label, "n_former_clusters": n})
    return pd.DataFrame(recs)


def no_difference_by_qcov_decile(res):
    """Cross-tab of the no_difference signature class against qcov deciles, over every call that
    reached a signature (status 'ok' or 'low_qcov') -- i.e. before the --min-qcov floor is
    applied to clustering. Motivates the floor: no_difference should not simply track low qcov."""
    d = res[res["status"].isin(["ok", "low_qcov"]) & res["qcov"].notna()].copy()
    edges = np.linspace(0.0, 1.0, 11)
    cols = ["qcov_decile", "n_calls", "n_no_difference", "frac_no_difference"]
    if not len(d):
        return pd.DataFrame(columns=cols)
    d["qcov_decile"] = pd.cut(d["qcov"], bins=edges, include_lowest=True)
    recs = []
    for decile, g in d.groupby("qcov_decile", observed=False):
        n = len(g)
        nd = int((g["signature_class"] == "no_difference").sum())
        recs.append({"qcov_decile": str(decile), "n_calls": n, "n_no_difference": nd,
                     "frac_no_difference": (nd / n) if n else np.nan})
    return pd.DataFrame(recs, columns=cols)


def homopolymer_fractions(ok):
    recs = []
    genes = sorted(ok["gene"].unique()) + ["ALL"]
    for scheme, col in (("pred", "ancestry_pred"), ("strict", "strict_ancestry")):
        for unrel_only in (False, True):
            d = ok[ok["unrelated"].fillna(False) == True] if unrel_only else ok  # noqa: E712
            for gene in genes:
                dg = d if gene == "ALL" else d[d["gene"] == gene]
                for anc in ["POOLED"] + vc.ANCESTRY_ORDER:
                    da = dg if anc == "POOLED" else dg[dg[col] == anc]
                    n = len(da)
                    counts = da["signature_class"].value_counts()
                    rec = {"gene": gene, "ancestry_scheme": scheme, "ancestry": anc,
                           "unrelated_only": unrel_only, "n_haplotypes": n,
                           "n_persons": da["person_id"].nunique()}
                    for cls in SIG_CLASSES:
                        rec[f"n_{cls}"] = int(counts.get(cls, 0))
                    rec["frac_homopolymer_only"] = (counts.get("homopolymer_only", 0) / n) if n else np.nan
                    recs.append(rec)
    return pd.DataFrame(recs)


def nonhp_sfs(ok):
    recs = []
    nonhp = ok[~ok["signature_class"].isin(["homopolymer_only", "no_difference"])]
    for unrel_only in (False, True):
        d = nonhp[nonhp["unrelated"].fillna(False) == True] if unrel_only else nonhp  # noqa: E712
        for gene in sorted(ok["gene"].unique()) + ["ALL"]:
            dg = d if gene == "ALL" else d[d["gene"] == gene]
            counts = dg.groupby("signature_id").size()
            for lo, hi in SFS_BINS:
                label = str(lo) if lo == hi else (f">{lo - 1}" if hi >= 10 ** 9 else f"{lo}-{hi}")
                recs.append({"gene": gene, "unrelated_only": unrel_only,
                             "carrier_haplotypes_bin": label,
                             "n_signatures": int(((counts >= lo) & (counts <= hi)).sum())})
    return pd.DataFrame(recs)


def touches_cds_diag(ok):
    """Depth-4 haplotypes whose differences touch the template's CDS. Expected only when the
    template (gene-level best) is not the allele whose CDS the call matched, otherwise a pipeline
    inconsistency. Split by whether the used allele shares the call's first 3 fields."""
    t = ok[ok["region_class"] == "touches_cds"]
    if not len(t):
        return {"n": 0}
    same3 = [key_prefix3(norm_allele(u)) == three_field_prefix(c)
             for u, c in zip(t["used_allele"], t["consensus"])]
    same3 = pd.Series(same3, index=t.index)
    return {"n": suppress(len(t)),
            "n_template_shares_call_3fields": suppress(int(same3.sum())),
            "n_template_differs_at_3fields": suppress(int((~same3).sum())),
            "by_selection_path": {k: suppress(v) for k, v in
                                  t["selection_path"].value_counts().to_dict().items()},
            "n_signatures": suppress(t["signature_id"].nunique())}


def gene_summary(calls, ok, clusters):
    recs = []
    for gene in sorted(calls["gene"].unique()):
        c = calls[calls["gene"] == gene]
        o = ok[ok["gene"] == gene]
        cl = clusters[clusters["gene"] == gene] if len(clusters) else clusters
        n_parsed = len(o)
        # small-denominator ratios are blanked, not just the small counts they're built from --
        # a lone float like frac_homopolymer_only=1.0 or median_qcov=0.87 can disclose one
        # person's summary under --genes all just as surely as an unsuppressed count would.
        frac_hp = float((o["signature_class"] == "homopolymer_only").mean()) if n_parsed else np.nan
        med_qcov = float(o["qcov"].median()) if n_parsed else np.nan
        if n_parsed < MIN_COMMIT:
            frac_hp = np.nan
            med_qcov = np.nan
        rec = {"gene": gene, "n_depth4_haplotypes": len(c), "n_parsed": n_parsed,
               "n_no_row": int((c["status"] == "no_row").sum()),
               "n_cs_error": int(c["status"].astype(str).str.startswith("cs_error").sum()),
               "n_split_alignment": int((c["status"] == "split_alignment").sum()),
               "n_low_qcov": int((c["status"] == "low_qcov").sum()),
               "n_multi_record_chained": int((c["chain_status"] == "multi_record_chained").sum()),
               "n_fallback_3field": int((o["selection_path"] == "fallback_3field").sum()),
               "n_context_sequence": int((o["context_mode"] == "sequence").sum()),
               "n_minus_strand": int((o["paf_strand"] == "-").sum()),
               "n_signatures": len(cl),
               "n_signatures_ge20_unrelated": int((cl["n_persons_unrelated"] >= MIN_COMMIT).sum())
               if len(cl) else 0,
               "n_former_novel_ids": o["novel_id"].dropna().nunique(),
               "frac_homopolymer_only": frac_hp,
               "median_qcov": med_qcov}
        for rc_ in REGION_CLASSES:
            rec[f"n_{rc_}"] = int((o["region_class"] == rc_).sum())
        recs.append(rec)
    return pd.DataFrame(recs)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def fig_pooling(pool_committed, dist_committed, path, dpi=200):
    """pool_committed: beyond_cds_pooling rows already filtered to n_persons >= MIN_COMMIT (as
    committed). dist_committed: pooling_distribution(pool) already passed through suppress_frame
    -- bins with n_former_clusters in 1..19 carry the string '<20' rather than the real count, and
    are drawn as a bar of height 0 with that label, so the figure never shows a small-denominator
    number the TSV would have blanked."""
    plt = _plt()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    heights = [v if isinstance(v, (int, float)) and not pd.isna(v) else 0
              for v in dist_committed["n_former_clusters"]]
    ax.bar(dist_committed["n_distinct_signatures_bin"], heights, color="#0072B2")
    for i, v in enumerate(dist_committed["n_former_clusters"]):
        if not (isinstance(v, (int, float)) and not pd.isna(v)):
            ax.text(i, 0.3, str(v), ha="center", fontsize=7, rotation=90, color="#555")
    ax.set_xlabel("Distinct non-coding signatures inside one former cluster (modal template)")
    ax.set_ylabel("Former depth-4 clusters")
    ax.set_title("How many sequences each 'beyond_cds' cluster pooled", fontsize=10)
    ax.tick_params(axis="x", rotation=30)
    ax = axes[1]
    big = pool_committed
    if len(big):
        ax.scatter(big["n_haplotypes"], big["n_distinct_signatures_modal_template"], s=18,
                   color="#D55E00", alpha=0.8, edgecolor="none")
        ax.plot([1, big["n_haplotypes"].max()], [1, big["n_haplotypes"].max()], ls="--",
                color="#999999", lw=0.8, label="one signature per haplotype")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend(frameon=False, fontsize=8)
    else:
        ax.text(0.5, 0.5, f"no former cluster with >= {MIN_COMMIT} persons", ha="center",
                transform=ax.transAxes)
    ax.set_xlabel("Haplotypes in former cluster")
    ax.set_ylabel("Distinct signatures (modal template)")
    ax.set_title(f"Former clusters with >= {MIN_COMMIT} persons", fontsize=10)
    for a in axes:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _stacked(ax, frame, labels, title):
    y = np.arange(len(labels))
    left = np.zeros(len(labels))
    for cls in SIG_CLASSES:
        vals = []
        for lab in labels:
            r = frame[lab]
            vals.append(r[f"n_{cls}"] / r["n_haplotypes"] if r["n_haplotypes"] >= MIN_COMMIT else 0)
        vals = np.array(vals)
        ax.barh(y, vals, left=left, color=SIG_COLORS[cls], label=cls.replace("_", " "))
        left += vals
    for i, lab in enumerate(labels):
        if frame[lab]["n_haplotypes"] < MIN_COMMIT:
            ax.text(0.01, i, f"n<{MIN_COMMIT}, not shown", va="center", fontsize=7, color="#555")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of depth-4 haplotypes")
    ax.set_title(title, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)


def fig_homopolymer(frac, path, dpi=200):
    plt = _plt()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    base = frac[(frac["ancestry_scheme"] == "pred") & (frac["unrelated_only"] == False)  # noqa
                & (frac["ancestry"] == "POOLED")]
    genes = [g for g in base["gene"] if g != "ALL"] + ["ALL"]
    _stacked(axes[0], {r["gene"]: r for _, r in base.iterrows()}, genes,
             "Per gene (all haplotypes)")
    anc = frac[(frac["ancestry_scheme"] == "strict") & (frac["unrelated_only"] == True)  # noqa
               & (frac["gene"] == "ALL")]
    labels = ["POOLED"] + vc.ANCESTRY_ORDER
    _stacked(axes[1], {r["ancestry"]: r for _, r in anc.iterrows()}, labels,
             "Per ancestry (strict >= 0.9, unrelated, all genes)")
    axes[1].legend(frameon=False, fontsize=8, loc="lower right", bbox_to_anchor=(1.0, 1.08),
                   ncol=3)
    fig.suptitle("What kind of difference explains depth-4 ('beyond CDS') novelty", y=1.08)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fig_sfs(sfs, path, dpi=200):
    plt = _plt()
    fig, ax = plt.subplots(figsize=(8, 4.2))
    d = sfs[sfs["gene"] == "ALL"]
    labels = [(str(lo) if lo == hi else (f">{lo - 1}" if hi >= 10 ** 9 else f"{lo}-{hi}"))
              for lo, hi in SFS_BINS]
    x = np.arange(len(labels))
    w = 0.4
    for off, unrel, color, name in ((-w / 2, False, "#999999", "all haplotypes"),
                                    (w / 2, True, "#0072B2", "unrelated only")):
        dd = d[d["unrelated_only"] == unrel].set_index("carrier_haplotypes_bin")
        vals = [dd["n_signatures"].get(lab, 0) for lab in labels]
        ax.bar(x + off, vals, width=w, color=color, label=name)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30)
    ax.set_yscale("symlog")
    ax.set_xlabel("Carrier haplotypes per signature")
    ax.set_ylabel("Non-homopolymer signatures")
    ax.set_title("Frequency spectrum of non-homopolymer non-coding signatures (all genes)",
                 fontsize=10)
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def _md(df):
    if df is None or not len(df):
        return "_(empty)_"
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for row in df.itertuples(index=False):
        lines.append("| " + " | ".join(
            (f"{v:.3f}" if isinstance(v, float) else str(v)) for v in row) + " |")
    return "\n".join(lines)


def write_report(path, summary, gsum_c, frac_c, pool_dist, qcov_decile, args):
    s = summary
    L = ["# Depth-4 ('beyond CDS') novelty, split by the real non-coding sequence", ""]
    L += ["## In plain words", "",
          "A depth-4 call such as `HLA-A*01:01:01:new` means the person's coding sequence is "
          "already known, but something outside it (an intron or UTR) differs from every "
          "catalogued allele. The earlier novel-allele table grouped these calls only by their "
          "coding sequence, so one 'novel allele' could hold thousands of different non-coding "
          "sequences. This report compares each call's contig with its template allele and "
          "writes the exact differences as a signature. Calls with identical signatures are "
          "grouped together.", "",
          f"- Depth-4 haplotypes examined: **{s['n_depth4_haplotypes']}**; with a readable "
          f"difference string: **{s['n_parsed']}**. Excluded before clustering: "
          f"**{s['status_counts'].get('split_alignment', 0)}** split alignments, "
          f"**{s['status_counts'].get('low_qcov', 0)}** below the --min-qcov floor "
          f"({args.min_qcov}).",
          f"- PAF selection: **{s['n_secondary_filtered']}** secondary (`tp:A:S`) records were "
          f"discarded before selection; **{s['n_multi_record_chained']}** calls needed multiple "
          f"primary records chained into one alignment; **{s['n_split_alignment']}** could not "
          "be chained (conflicting/overlapping primary records) and were excluded.",
          f"- Distinct signatures: **{s['n_signatures']}**. Former clusters split: "
          f"**{s['n_former_clusters']}** former clusters contain a median of "
          f"**{s['median_signatures_per_former_cluster']}** distinct signatures using each "
          "cluster's single most common (modal) template "
          f"(max {s['max_signatures_per_former_cluster']}); using the raw, template-confounded "
          f"count instead the median is **{s['median_signatures_per_former_cluster_raw']}** "
          f"(max {s['max_signatures_per_former_cluster_raw']}) -- the modal-template number is "
          "the honest headline (see Caveats).",
          f"- Fraction of depth-4 haplotypes whose only differences are homopolymer indels "
          f"(the typical long-read error pattern): **{s['frac_homopolymer_only']}**"
          + (f" ({s['homopolymer_gate_reason']})" if s.get("homopolymer_gate_reason") else "")
          + ".",
          f"- Homopolymer context was checked against the real allele sequence for "
          f"**{s['frac_context_sequence']}** of haplotypes. For the rest, the flag only means "
          "'the inserted or deleted bases are all one letter', which also catches every "
          "1-bp indel, so it over-counts.", "",
          f"Counts from 1 to {MIN_COMMIT - 1} are written as `<{MIN_COMMIT}`. Signatures are "
          f"listed only when they have at least {MIN_COMMIT} unrelated carriers. Any ratio whose "
          f"denominator count is below {MIN_COMMIT} is blanked (`NA`) rather than shown, in "
          "every committed table, `summary.json`, and figure.", "",
          "## How the difference string is read", "",
          "minimap2 aligned every known genomic allele (query) to the person's contig (target). "
          "In its cs string, `*ab` means contig base a and allele base b, `+seq` means bases the "
          "contig lacks, and `-seq` means extra bases in the contig. Events are written on the "
          "known allele's own 1-based coordinates as `posA>G`, `posdelX` or `p_p+1insX`. For "
          "alignments on the reverse strand, positions are counted back from the alignment end "
          "and the bases are reverse-complemented, so the same difference gives the same "
          "signature on either strand.", "",
          "## Per gene", "", _md(gsum_c), "",
          "## Former clusters by number of distinct signatures (modal template, the headline)",
          "", _md(pool_dist), "",
          "## `no_difference` calls by alignment-coverage decile", "",
          "Motivates the `--min-qcov` floor: a call that only partly covers the allele can "
          "wrongly look identical to it. Computed before the floor excludes low-coverage calls "
          "from clustering.", "", _md(qcov_decile), "",
          "## Homopolymer-only fraction (all genes, strict ancestry, unrelated only)", ""]
    if s.get("homopolymer_gate_reason"):
        L += [f"_Not shown: {s['homopolymer_gate_reason']}_", ""]
    else:
        L += [_md(frac_c[(frac_c["gene"] == "ALL") & (frac_c["ancestry_scheme"] == "strict")
                         & (frac_c["unrelated_only"] == True)]  # noqa: E712
                  [["ancestry", "n_haplotypes", "n_homopolymer_only", "n_snv_only",
                    "n_snv_plus_indel", "n_indel_nonhomopolymer", "frac_homopolymer_only"]]), ""]
    L += ["## Where the differences are", "",
          "Each event is placed on the template allele's own map from `alleles.csv.gz` (UTR5, "
          "exonN, intronN, UTR3). A depth-4 call has a known coding sequence, so its differences "
          "should sit in introns or UTRs. Haplotypes whose differences touch the coding sequence "
          "point either to a template that is not the allele the call matched, or to a "
          "pipeline inconsistency.", "",
          f"- Share of parsed haplotypes by region class: `{json.dumps(s['frac_region_class'])}`",
          f"- Touching CDS: `{json.dumps(s['touches_cds'])}`", "",
          "## Cohort coverage", "",
          f"- **{s['n_missing_from_cohort_membership']}** haplotypes belong to a person absent "
          "from `cohort_membership.tsv`; their relatedness status is recorded as `NA` (a third "
          "state, not silently treated as related or unrelated).",
          f"- Cluster-level `signature_class`/`region_class` are taken by majority vote across "
          "member calls; disagreements are recorded per cluster "
          "(`n_signature_class_disagree`, `n_region_class_disagree`) -- they are expected to be "
          "zero, since both are pure functions of the signature.", "",
          "## Caveats", "",
          "- The signature is keyed on the template allele. Two identical sequences with "
          "different templates count as two raw signatures; the modal-template number used as "
          "the headline above avoids this confound by restricting to each cluster's single most "
          "common template.",
          "- Differences outside the aligned part of the allele are not seen "
          "(see `median_qcov`); calls covering less than "
          f"{args.min_qcov} of the allele are excluded from clustering entirely.",
          "- Without the allele sequence, an indel inside a repeat may be placed differently on "
          "the two strands and split one signature into two.",
          "- Former-cluster mapping skips genes with more than one copy on a contig (same rule "
          "as 03).",
          "- Multi-record chaining assumes non-overlapping, order-consistent primary alignments; "
          "anything else is excluded as `split_alignment` rather than guessed at.", "",
          f"Run: genes={args.genes}, limit={args.limit}, homopolymer_min_run="
          f"{args.homopolymer_min_run}, min_qcov={args.min_qcov}, threads={args.threads}."]
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def load_novel_matches(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    need = {"person_id", "hap", "contig", "gene", "cds_seq_sha1"}
    miss = need - set(df.columns)
    if miss:
        sys.exit(f"FATAL: --novel-matches {path} lacks columns {sorted(miss)}")
    return {(r.person_id, r.hap, r.contig, r.gene): r.cds_seq_sha1
            for r in df.itertuples(index=False)}


def _file_stat_tuple(path):
    """(abspath, size, mtime) or (path, None, None) if it can't be stat'd -- used so the resume
    checkpoint is invalidated whenever a FASTA/annotation file's content could plausibly have
    changed, not just when Table 1 or the CLI args change."""
    try:
        st = os.stat(path)
        return [os.path.abspath(path), st.st_size, int(st.st_mtime)]
    except OSError:
        return [path, None, None]


def params_hash(args, extra):
    st = os.stat(args.table1)
    script_st = os.stat(os.path.abspath(__file__))
    blob = json.dumps({"t1": [os.path.abspath(args.table1), st.st_size, int(st.st_mtime)],
                       "genes": args.genes, "limit": args.limit, "min_run": args.homopolymer_min_run,
                       "min_qcov": args.min_qcov,
                       "novel_matches": args.novel_matches, "outroot": os.path.abspath(args.outroot),
                       "script": [os.path.abspath(__file__), script_st.st_size,
                                 int(script_st.st_mtime)],
                       **extra}, sort_keys=True)
    return hashlib.sha1(blob.encode()).hexdigest()[:12]


def _process_hap_task(task, cfg):
    pid, hap, calls = task
    return process_hap(pid, hap, calls, cfg)


def _process_hap_probe():
    """Trivial module-level function used to test whether a fresh worker process can actually
    import/pickle code from this module before committing to ProcessPoolExecutor. A dynamically
    loaded copy of this module (e.g. under a test harness using importlib, or any environment
    where this file isn't reachable as a normal import) will fail this probe even though a
    builtin like `int` would round-trip fine, which is why we probe with our own function."""
    return True


def _make_executor(threads):
    """-> (executor, kind). Prefers ProcessPoolExecutor: PAF-line parsing (parse_cs_events, the
    per-line regex scan) is pure Python and GIL-bound, so on the VM's 4 vCPUs threads barely help.
    Falls back to ThreadPoolExecutor if worker processes can't actually run our code (e.g. a
    sandboxed/fork-restricted environment, or this module loaded in a way children can't import)."""
    if threads <= 1:
        return None, "none"
    try:
        ex = concurrent.futures.ProcessPoolExecutor(max_workers=threads)
        if not ex.submit(_process_hap_probe).result(timeout=60):
            raise RuntimeError("probe returned an unexpected result")
        return ex, "process"
    except Exception as e:
        try:
            ex.shutdown(wait=False, cancel_futures=True)
        except Exception:
            pass
        print(f"WARNING: ProcessPoolExecutor unavailable ({type(e).__name__}: {e}); falling "
              "back to ThreadPoolExecutor. PAF parsing is GIL-bound so this will be slower on a "
              "multi-core VM, but threads still overlap file I/O.", file=sys.stderr)
        return concurrent.futures.ThreadPoolExecutor(max_workers=threads), "thread"


def run_calls(tasks, cfg, threads, ckpt_dir, resume=True):
    """tasks: list of (pid, hap, calls) sorted by pid. Checkpoints per batch of persons."""
    os.makedirs(ckpt_dir, exist_ok=True)
    pids = sorted({t[0] for t in tasks})
    by_pid = defaultdict(list)
    for t in tasks:
        by_pid[t[0]].append(t)
    frames = []
    t0 = time.time()
    worker = partial(_process_hap_task, cfg=cfg)
    executor, kind = _make_executor(threads)
    try:
        for b in range(0, len(pids), BATCH_PERSONS):
            fn = os.path.join(ckpt_dir, f"batch_{b // BATCH_PERSONS:05d}.pkl")
            if resume and os.path.exists(fn):
                frames.append(pd.read_pickle(fn))
                continue
            batch = [t for pid in pids[b:b + BATCH_PERSONS] for t in by_pid[pid]]
            results = list(executor.map(worker, batch)) if executor else [worker(t) for t in batch]
            rows = [r for res in results for r in res]
            df = pd.DataFrame(rows)
            tmp = fn + ".tmp"
            df.to_pickle(tmp)
            os.replace(tmp, fn)
            frames.append(df)
            done = min(b + BATCH_PERSONS, len(pids))
            el = time.time() - t0
            print(f"  {done}/{len(pids)} persons, {el:.0f}s ({kind})", file=sys.stderr)
    finally:
        if executor:
            executor.shutdown()
    frames = [f for f in frames if len(f)]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=DEFAULT_TABLE1)
    ap.add_argument("--cohort-membership", default=DEFAULT_COHORT)
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT)
    ap.add_argument("--relatedness-table", default=DEFAULT_RELATEDNESS)
    ap.add_argument("--allow-missing-relatedness", action="store_true",
                    help="Proceed with unrelated=NA if the relatedness table is absent.")
    ap.add_argument("--novel-alleles", default=DEFAULT_NOVEL_ALLELES)
    ap.add_argument("--novel-matches", default=None,
                    help="Optional TSV person_id,hap,contig,gene,cds_seq_sha1 (skips cds.fa.gz).")
    ap.add_argument("--refdata", default=DEFAULT_REFDATA)
    ap.add_argument("--allele-fasta", nargs="*", default=None,
                    help="Genomic IPD allele FASTA(s); default: auto-discover under --refdata.")
    ap.add_argument("--allele-annotation", default=None,
                    help="alleles.csv.gz; default <refdata>/alleles.csv.gz if present.")
    ap.add_argument("--genes", nargs="+", default=list(vc.CLASSICAL_GENES),
                    help="Gene names, or 'all' for every gene in Table 1.")
    ap.add_argument("--limit", type=int, default=None, help="Pilot: first N persons (sorted).")
    ap.add_argument("--threads", type=int, default=4,
                    help="Worker processes (VM: 4 vCPU / 31 GB). PAF-line parsing is pure-Python "
                    "and GIL-bound, so this uses ProcessPoolExecutor by default, falling back to "
                    "threads only if worker processes fail to start.")
    ap.add_argument("--homopolymer-min-run", type=int, default=HOMOPOLYMER_MIN_RUN)
    ap.add_argument("--min-qcov", type=float, default=MIN_QCOV,
                    help="Calls covering less of the allele than this (chained query coverage / "
                    "qlen) are excluded from clustering and marked status=low_qcov.")
    ap.add_argument("--allow-content-only", action="store_true",
                    help="Emit the homopolymer figure/fraction even when the allele sequence "
                    "was unavailable for most calls (frac_context_sequence < "
                    f"{CONTEXT_SEQUENCE_MIN_FRAC}), overriding the honesty gate. For debugging "
                    "only: without the sequence, every 1-bp indel is flagged as a homopolymer.")
    ap.add_argument("--kin-threshold", type=float, default=KIN_THRESHOLD)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    ap.add_argument("--local-dir", default=DEFAULT_LOCAL_DIR)
    ap.add_argument("--no-resume", action="store_true")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)

    out_dir = os.path.expanduser(args.out_dir)
    local_dir = os.path.expanduser(args.local_dir)
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(local_dir, exist_ok=True)
    if not os.path.isdir(args.outroot):
        sys.exit(f"FATAL: --outroot {args.outroot} does not exist")

    # --- Table 1 ---
    t1 = vc.load_table1(args.table1)
    for col in ("template_allele", "gene_start", "gene_end", "novelty_depth", "contig"):
        if col not in t1.columns:
            sys.exit(f"FATAL: Table 1 lacks column {col!r}")
    t1["novelty_depth"] = pd.to_numeric(t1["novelty_depth"], errors="coerce")
    group_sizes = t1.groupby(["person_id", "hap", "contig", "gene"]).size().to_dict()
    genes_all = "all" in [g.lower() for g in args.genes]
    d4 = t1[t1["novelty_depth"] == 4]
    if not genes_all:
        d4 = d4[d4["gene"].isin(args.genes)]
    d4 = d4[d4["gene_start"].notna() & d4["gene_end"].notna()]
    persons = sorted(d4["person_id"].unique())
    if args.limit:
        persons = persons[:args.limit]
        d4 = d4[d4["person_id"].isin(set(persons))]
    print(f"{len(d4)} depth-4 calls across {len(persons)} persons", file=sys.stderr)
    if d4.empty:
        sys.exit("FATAL: no depth-4 calls selected")

    calls = []
    for r in d4.itertuples(index=False):
        calls.append({
            "person_id": r.person_id, "hap": r.hap, "contig": r.contig, "gene": r.gene,
            "copy_index": getattr(r, "copy_index", 1), "consensus": r.consensus,
            "template_allele": r.template_allele if isinstance(r.template_allele, str) else None,
            "template_key": norm_allele(r.template_allele),
            "prefix_key": three_field_prefix(r.consensus),
            "gene_start": int(r.gene_start), "gene_end": int(r.gene_end)})
    wanted_exact = {c["template_key"] for c in calls if c["template_key"]}
    wanted_prefix = {c["prefix_key"] for c in calls if c["prefix_key"]}

    # --- reference data ---
    fastas = args.allele_fasta if args.allele_fasta is not None else find_allele_fastas(args.refdata)
    allele_seqs = load_allele_sequences(fastas, wanted_exact, wanted_prefix) if fastas else {}
    ann_path = args.allele_annotation or os.path.join(args.refdata, "alleles.csv.gz")
    annotations = load_annotations(ann_path, wanted_exact, wanted_prefix)
    print(f"allele FASTA(s): {fastas or 'none'} -> {len(allele_seqs)} sequences; "
          f"annotations: {len(annotations)}", file=sys.stderr)

    novel_matches = load_novel_matches(args.novel_matches) if args.novel_matches else None

    # --- cohort ---
    cohort = vc.load_cohort_membership(args.cohort_membership)
    cohort["strict_ancestry"] = strict_ancestry(cohort)
    lr = set(cohort.loc[cohort["in_lr"], "person_id"]) if "in_lr" in cohort.columns \
        else set(t1["person_id"])
    try:
        rel = _m11.load_relatedness(args.relatedness_table)
        pairs = zip(rel["person_i"], rel["person_j"], rel["kin"].astype(float))
        unrelated = greedy_unrelated(pairs, lr, args.kin_threshold)
        have_rel = True
    except FileNotFoundError:
        if not args.allow_missing_relatedness:
            raise
        unrelated, have_rel = set(), False
    pd.DataFrame({"person_id": sorted(lr), "unrelated": [p in unrelated for p in sorted(lr)]}) \
        .to_csv(os.path.join(local_dir, "25_unrelated_set.tsv"), sep="\t", index=False)

    # --- per-hap processing ---
    by_hap = defaultdict(list)
    for c in calls:
        by_hap[(c["person_id"], c["hap"])].append(c)
    tasks = [(pid, hap, cs) for (pid, hap), cs in sorted(by_hap.items())]
    cfg = {"outroot": args.outroot, "group_sizes": group_sizes, "novel_matches": novel_matches,
           "allele_seqs": allele_seqs, "annotations": annotations,
           "min_run": args.homopolymer_min_run, "min_qcov": args.min_qcov}
    ann_stat = _file_stat_tuple(ann_path) if os.path.exists(ann_path) else None
    h = params_hash(args, {"fastas": [_file_stat_tuple(f) for f in fastas], "ann": ann_stat})
    ckpt = os.path.join(local_dir, "25_checkpoints", h)
    res = run_calls(tasks, cfg, args.threads, ckpt, resume=not args.no_resume)

    anc = cohort.set_index("person_id")
    res["ancestry_pred"] = res["person_id"].map(anc["ancestry_pred"]) if "ancestry_pred" in anc \
        else pd.NA
    res["strict_ancestry"] = res["person_id"].map(anc["strict_ancestry"])
    # People absent from cohort_membership get a third state (NA), not silently "related"
    # (previously False) or "unrelated" -- their relatedness status is simply unknown. Without a
    # relatedness table every unrelated_only view is empty too (summary flags it via
    # relatedness_available).
    cohort_ids = set(cohort["person_id"])
    n_missing_from_cohort = int((~res["person_id"].isin(cohort_ids)).sum())

    def _unrelated_state(p):
        if p not in cohort_ids or not have_rel:
            return pd.NA
        return p in unrelated

    res["unrelated"] = res["person_id"].map(_unrelated_state)
    for col in ("signature_id", "signature_class", "context_mode", "event_tokens",
                "homopolymer_context_checked", "n_events_cds", "event_regions", "region_class",
                "n_events_intron", "n_events_utr", "n_events_unannotated",
                "event_flank_runs", "seq_consistent", "n_events", "n_snv", "n_indel",
                "n_homopolymer_indel", "signature", "allele_seq_len_mismatch"):
        if col not in res.columns:
            res[col] = pd.NA
    res.to_csv(os.path.join(local_dir, "25_noncoding_calls.tsv"), sep="\t", index=False)

    ok = res[res["status"] == "ok"].copy()
    novel_tsv = None
    if args.novel_alleles and os.path.exists(args.novel_alleles):
        novel_tsv = pd.read_csv(args.novel_alleles, sep="\t").set_index("novel_id")
    else:
        print(f"WARNING: {args.novel_alleles} not found; novel_ids derived from hashes only",
              file=sys.stderr)

    clusters = build_clusters(ok)
    pool = build_pooling(ok, novel_tsv)
    frac = homopolymer_fractions(ok)
    sfs = nonhp_sfs(ok)
    gsum = gene_summary(res, ok, clusters)
    dist = pooling_distribution(pool)                              # modal-template (headline)
    dist_raw = pooling_distribution(pool, col="n_distinct_signatures")  # raw, template-confounded
    qcov_decile = no_difference_by_qcov_decile(res)

    # VM-local full versions
    clusters.to_csv(os.path.join(local_dir, "25_noncoding_signature_clusters_full.tsv"), sep="\t",
                    index=False)
    pool.to_csv(os.path.join(local_dir, "25_beyond_cds_pooling_full.tsv"), sep="\t", index=False)
    frac.to_csv(os.path.join(local_dir, "25_homopolymer_fraction_full.tsv"), sep="\t", index=False)
    sfs.to_csv(os.path.join(local_dir, "25_nonhomopolymer_sfs_full.tsv"), sep="\t", index=False)
    gsum.to_csv(os.path.join(local_dir, "25_gene_summary_full.tsv"), sep="\t", index=False)
    dist_raw.to_csv(os.path.join(local_dir, "25_pooling_distribution_raw_full.tsv"), sep="\t",
                    index=False)
    qcov_decile.to_csv(os.path.join(local_dir, "25_no_difference_by_qcov_decile_full.tsv"),
                       sep="\t", index=False)

    # Committed, suppressed versions
    count_cols_cl = ["n_haplotypes", "n_persons", "n_persons_unrelated", "n_haplotypes_unrelated",
                     "n_plus_strand", "n_minus_strand", "n_fallback_rows", "n_partial_alignment",
                     "n_multi_record_chained", "n_signature_class_disagree",
                     "n_region_class_disagree"] + \
        [c for c in clusters.columns if c.startswith("n_persons_pred_") or
         c.startswith("n_persons_strict_")]
    cl_c = clusters[clusters["n_persons_unrelated"] >= MIN_COMMIT] if len(clusters) else clusters
    suppress_frame(cl_c, count_cols_cl).to_csv(
        os.path.join(out_dir, "noncoding_signature_clusters.tsv"), sep="\t", index=False)
    pool_c = pool[pool["n_persons"] >= MIN_COMMIT] if len(pool) else pool
    suppress_frame(pool_c, ["n_haplotypes", "n_persons", "n_persons_unrelated",
                            "n_distinct_templates", "n_haplotypes_homopolymer_only",
                            "old_n_persons", "old_n_haplotypes"]).to_csv(
        os.path.join(out_dir, "beyond_cds_pooling.tsv"), sep="\t", index=False)
    dist_c = suppress_frame(dist, ["n_former_clusters"])
    dist_c.to_csv(os.path.join(out_dir, "pooling_distribution.tsv"), sep="\t", index=False)
    frac_c = frac.copy()
    frac_c.loc[frac_c["n_haplotypes"] < MIN_COMMIT, "frac_homopolymer_only"] = np.nan
    frac_c = suppress_frame(frac_c, ["n_haplotypes", "n_persons"] +
                            [f"n_{c}" for c in SIG_CLASSES])
    n_context_ok = int((ok["context_mode"] == "sequence").sum()) if len(ok) else 0
    frac_context_sequence = round(n_context_ok / len(ok), 4) if len(ok) else None
    homopolymer_gate = (frac_context_sequence is not None
                        and frac_context_sequence < CONTEXT_SEQUENCE_MIN_FRAC
                        and not args.allow_content_only)
    homopolymer_gate_reason = None
    if homopolymer_gate:
        homopolymer_gate_reason = (
            f"frac_context_sequence={frac_context_sequence} < {CONTEXT_SEQUENCE_MIN_FRAC}: the "
            "allele sequence needed to verify homopolymer context was unavailable for most "
            "calls, so the content-only flag (every 1-bp indel counts) would overstate the "
            "homopolymer fraction. Re-run with --allow-content-only to override for debugging.")
        frac_c["frac_homopolymer_only"] = np.nan
    frac_c.to_csv(os.path.join(out_dir, "homopolymer_fraction.tsv"), sep="\t", index=False)
    suppress_frame(sfs, ["n_signatures"]).to_csv(
        os.path.join(out_dir, "nonhomopolymer_sfs.tsv"), sep="\t", index=False)
    gsum_c = suppress_frame(gsum, ["n_depth4_haplotypes", "n_parsed", "n_no_row", "n_cs_error",
                                   "n_split_alignment", "n_low_qcov", "n_multi_record_chained",
                                   "n_fallback_3field", "n_context_sequence", "n_minus_strand",
                                   "n_signatures", "n_signatures_ge20_unrelated",
                                   "n_former_novel_ids"] + [f"n_{c}" for c in REGION_CLASSES])
    if homopolymer_gate:
        gsum_c["frac_homopolymer_only"] = np.nan
    gsum_c.to_csv(os.path.join(out_dir, "gene_summary.tsv"), sep="\t", index=False)
    qcov_decile_c = suppress_frame(qcov_decile, ["n_calls", "n_no_difference"])
    qcov_decile_c.loc[qcov_decile["n_calls"] < MIN_COMMIT, "frac_no_difference"] = np.nan
    qcov_decile_c.to_csv(os.path.join(out_dir, "no_difference_by_qcov_decile.tsv"), sep="\t",
                         index=False)

    status_counts = res["status"].astype(str).str.split(":").str[0].value_counts().to_dict()
    n_secondary_filtered = int(
        res.drop_duplicates(["person_id", "hap"])["n_secondary_filtered_hap"].sum())
    median_modal = pool["n_distinct_signatures_modal_template"].median() if len(pool) else None
    max_modal = pool["n_distinct_signatures_modal_template"].max() if len(pool) else None
    median_raw = pool["n_distinct_signatures"].median() if len(pool) else None
    max_raw = pool["n_distinct_signatures"].max() if len(pool) else None
    frac_hp_overall = float((ok["signature_class"] == "homopolymer_only").mean()) \
        if len(ok) else None
    if frac_hp_overall is not None:
        frac_hp_overall = round(frac_hp_overall, 4) if len(ok) >= MIN_COMMIT else None
    if homopolymer_gate:
        frac_hp_overall = None
    median_qcov_overall = float(ok["qcov"].median()) if len(ok) >= MIN_COMMIT else None
    summary = {
        "n_persons": suppress(len(persons)),
        "n_depth4_haplotypes": suppress(len(res)),
        "n_parsed": suppress(len(ok)),
        "status_counts": {k: suppress(v) for k, v in status_counts.items()},
        "n_secondary_filtered": suppress(n_secondary_filtered),
        "n_multi_record_chained": suppress(int((res["chain_status"] ==
                                                "multi_record_chained").sum())),
        "n_split_alignment": suppress(int((res["status"] == "split_alignment").sum())),
        "n_low_qcov": suppress(int((res["status"] == "low_qcov").sum())),
        "min_qcov": args.min_qcov,
        "selection_path_counts": {k: suppress(v) for k, v in
                                  ok["selection_path"].value_counts().to_dict().items()},
        "cds_status_counts": {k: suppress(v) for k, v in
                              res["cds_status"].value_counts().to_dict().items()},
        "n_signatures": suppress(len(clusters)),
        "n_signatures_ge20_unrelated": suppress(len(cl_c)),
        "n_former_clusters": suppress(len(pool)),
        "median_signatures_per_former_cluster": float(median_modal) if median_modal is not None
        else None,
        "max_signatures_per_former_cluster": int(max_modal) if max_modal is not None else None,
        "median_signatures_per_former_cluster_raw": float(median_raw) if median_raw is not None
        else None,
        "max_signatures_per_former_cluster_raw": int(max_raw) if max_raw is not None else None,
        "frac_homopolymer_only": frac_hp_overall,
        "homopolymer_gate_reason": homopolymer_gate_reason,
        "signature_class_haplotypes": {k: suppress(v) for k, v in
                                       ok["signature_class"].value_counts().to_dict().items()},
        "region_class_haplotypes": {k: suppress(v) for k, v in
                                    ok["region_class"].value_counts().to_dict().items()},
        "frac_region_class": {k: round(float(v), 4) for k, v in
                              ok["region_class"].value_counts(normalize=True).to_dict().items()}
        if len(ok) >= MIN_COMMIT else {},
        "touches_cds": touches_cds_diag(ok),
        "no_difference_by_qcov_decile": json.loads(
            qcov_decile_c.to_json(orient="records")) if len(qcov_decile_c) else [],
        "frac_context_sequence": frac_context_sequence,
        "n_seq_inconsistent": suppress(int((ok["seq_consistent"] == False).sum())),  # noqa: E712
        "strand_counts": {k: suppress(v) for k, v in ok["paf_strand"].value_counts().to_dict().items()},
        "median_qcov": median_qcov_overall,
        "relatedness_available": have_rel,
        "n_unrelated_lr": suppress(len(unrelated)),
        "n_missing_from_cohort_membership": suppress(n_missing_from_cohort),
        "allele_fastas": [os.path.basename(f) for f in fastas],
        "n_allele_sequences_loaded": len(allele_seqs),
        "n_annotations_loaded": len(annotations),
        "genes": "all" if genes_all else args.genes,
        "limit": args.limit,
        "min_commit_count": MIN_COMMIT,
    }
    with open(os.path.join(out_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    write_report(os.path.join(out_dir, "noncoding_novelty_report.md"), summary, gsum_c, frac_c,
                 dist_c, qcov_decile_c, args)
    if not args.no_figures:
        fig_pooling(pool_c, dist_c, os.path.join(out_dir, "fig1_signatures_per_former_cluster.png"))
        if not homopolymer_gate:
            # frac (not frac_c): _stacked() already re-derives each bar from n_{cls}/n_haplotypes
            # and zeroes it out below MIN_COMMIT itself, so it needs the numeric columns intact.
            fig_homopolymer(frac, os.path.join(out_dir, "fig2_homopolymer_fraction.png"))
        fig_sfs(sfs, os.path.join(out_dir, "fig3_nonhomopolymer_sfs.png"))
    print(json.dumps(summary, indent=2, default=str), file=sys.stderr)
    print(f"committed outputs -> {out_dir}\nlocal outputs -> {local_dir}", file=sys.stderr)
    return summary


if __name__ == "__main__":
    main()
