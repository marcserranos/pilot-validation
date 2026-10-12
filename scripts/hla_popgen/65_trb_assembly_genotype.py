#!/usr/bin/env python3
"""TRB insertion/deletion genotype (TRBV4-3, TRBV3-2, TRBV6-3, TRBV7-8): long reads (default) or assemblies.

Sprint S05, WS-B. The TRB locus carries a common ~21 kb insertion/deletion polymorphism holding
TRBV4-3, TRBV3-2, TRBV7-8 and TRBV6-3 (a near-identical duplicate of TRBV6-2). GRCh38's primary chr7
carries the haplotype WITHOUT the segment; the segment exists only on the ALT contig
chr7_KI270803v1_alt (placed at chr7:142,038,120-143,088,503; confirmed on the VM 10-11).
58_regional_depth.py reads depth on noalt BAMs and cannot see TRBV4-3 at all. Lead it serves: Aleix's
repertoire data, TRBV4-3 never used by ~27% of people (SAS 43.5%, EUR 16.7%). Prior art: READING.md B.

Two genotyping routes (`--mode`), same probes, same output tables:

  reads (default, fast)   per person, `samtools view -F 0x904` of the GRCh38 long-read BAM
      (`--bam-column`, default grch38_bam; grch38_haplotagged_bam carries HP tags) over the primary
      TRB span +/- 50 kb (`--reads-region`, chr7:142,249,011-142,863,287). Two kinds of informative read:
        reference-allele read: aligned contiguously across the insertion site, covering
          [breakpoint_left - flank, breakpoint_right + flank] (`--flank` 1000) with no insertion,
          deletion or skip >= `--max-indel` (50 bp) inside that interval. Only a haplotype without
          the segment can produce one (the segment is 21 kb; such a read would carry it as a huge I
          or be clipped at the breakpoint).
        insertion-allele read: any read in which a probe of TRBV4-3 / TRBV3-2
          (`--ins-genes`) passes --min-ident/--min-cov under best-hit assignment (below). Reads
          spanning a breakpoint into the segment are placed on the flank with the segment
          soft-clipped, so they are in the fetched window; reads entirely inside the segment are
          unmapped or mis-placed and mostly missed (see caveats: this biases n_ins down for genes far
          from both breakpoints).
      Per gene and for 'insertion_any': AF = n_ins / (n_ins + n_ref); copies 0 if AF < 0.2,
      1 if < 0.8, else 2 (`--af-cuts`); status 'ok' if n_ins + n_ref >= `--min-reads` (6) else
      'low_reads'. If reads carry HP tags (haplotagged BAM), per-haplotype calls too: a haplotype is
      'complete' with >= `--min-hap-reads` (3) informative reads, and carries the gene if its own AF
      >= 0.5. TRBV6-3 cannot be called from reads (identical to 6-2, which is present on both alleles).
      The breakpoint comes from the probe build (probe_meta.tsv next to probes.fa) or `--breakpoint`.

  assembly (validation)   per person, hap1/hap2 hifiasm assemblies (.fa.gz) streamed from the mount
      (20-112 s per haplotype on the VM: too slow for 12k people, used on a stratified subset via
      `--people-file`). TRB contigs = contigs whose primary PAF alignments cover >= 1 kb of
      `--window`, cut to the projected window +/- `--pad`. Presence = called on ANY of the
      haplotype's TRB contigs; complete = every anchor (TRBV2 5', TRBC1 3') called among those
      contigs; the old strict rule (both anchors on ONE contig) is kept as `complete_single_contig`.
      copies = hap1 + hap2 presence, status 'ok' only if both haplotypes complete.
      TRBV6-2/6-3 is one merged class; '<class>[copy2]' (>= 2 loci on the haplotype) is its proxy.

Probe build (`--build-probes`, once, VM with internet). Route A: IMGT/GENE-DB V-REGION nucleotides
for all human TRBV alleles (GENElect query 7.2; 240 alleles on 10-11), gap dots removed, all TRBV kept
as background probes for best-hit assignment; TRBC1 genomic from Ensembl REST (its IMGT C-REGION is
spliced). Route B: `--probe-fasta`. Writes probes.fa, probe_classes.tsv and, with `--alt-fasta` (the
ALT contig) and `--ref-fasta` (+ `--ref-region` when it is a whole chromosome, e.g. UCSC chr7.fa.gz
with chr7:142038120-143088503), reference_check.tsv and probe_meta.tsv with the breakpoint:
  1 kb ALT chunks every 250 bp over 40 kb on each side of the insertion genes (as located on the ALT)
  are aligned to the primary span (asm5); chunks with one forward hit (>= 98% identity, >= 90%
  covered) give ref = alt + offset. The left and right offsets come from the innermost run of >= 3
  consecutive chunks sharing one offset on each side (the ALT carries other indels further out); the breakpoint interval on primary is [min, max] of (last
  left-colinear hit end, first right-colinear hit start). An interval rather than a point covers
  microhomology or a duplicated unit at the junction. ins_len_estimate = left offset - right offset.

Alignment (both modes): probes are the query, contigs or reads the target; minimap2 preset `sr`,
base-level (`-c`), secondaries kept (`--secondary=yes -N <n> -p 0.5`). Why not asm5: 300 bp probes
with >= 5 mismatches (98.3%) fall below asm5's minimum score and weaker secondaries are dropped; the
paralog hit must be seen to be rejected. Identity = matches / alignment block; coverage = aligned
probe span / length. Loci = overlapping hits (>= 50% of the shorter); each locus goes to the class of
its best probe (>= --min-cov first, then identity); called if identity >= --min-ident (0.97) and
coverage >= --min-cov (0.9); a tie within --tie-eps between classes is 'ambiguous'. Genes whose
alleles are >= --merge-ident (0.99) identical are merged into one class. Aligners: mappy (in-process),
the minimap2 binary (`--aligner minimap2 --minimap2 PATH`), or a pure-Python k-mer fallback.

Other actions
  --sample-stratified N   unrelated (greedy, kin >= 0.0442), ancestry-balanced (equal shares, spare
                          quota passed to larger groups), seeded (--seed) sample of people with both
                          assembly paths -> --people-file (INTERNAL path). Then run --mode assembly
                          --people-file on it.
  --compare READS ASM     65b: per-person copies from two INTERNAL per_person_gene.tsv tables (status ok in
                          both) -> PUBLIC 65b_method_concordance.tsv (per feature x stratum: n, %
                          concordant + Wilson CI, 3x3 confusion; masked).

Outputs
  INTERNAL (default ~/s05/internal/65_<mode>; refuses the repo and reports/; mode recorded in
  65_mode.txt, a mismatching mode or schema refuses to append):
    per_person_gene.tsv   person x feature: hap1/2_present, hap1/2_n_loci, copies, hap1/2_complete,
                          status, best_ident, method, n_ins, n_ref, allele_fraction,
                          hap1/2_n_ins, hap1/2_n_ref                       (checkpoint, append+fsync)
    per_haplotype.tsv     person x hap: contigs/bp (assembly) or reads (reads), complete,
                          complete_single_contig, reason, timing            (the 'done' marker)
    failures.tsv          person, reason
  PUBLIC (default ~/s05/results/65_<mode>; counts 1-19 -> '<20', dependent rates/frequencies/HWE blank):
    probe_classes.tsv, genotype_by_stratum.tsv (n measured, n called = status ok, genotype counts,
    insertion-allele frequency, % 0 copies, HWE exact p from 58), haplotype_completeness.tsv,
    copresence.tsv + haplotype_patterns.tsv (genes on the same complete haplotype),
    allele_fraction_histogram.tsv (reads mode), failures_summary.tsv, run_summary.tsv

Robustness (as 58): mount + manifest hard check; PID lock; per-person append + fsync; resume skips done
people (`--retry-failed`); `--timeout` per person; `--limit`; `--workers` threads.

VM (exact commands: sprints/S05_qc_phasing_kir_trb/vm_parts/65.sh):
    PY=~/repos/pilot-validation/.pixi/envs/spechla/bin/python3
    export PYTHONPATH=~/envs/mappy38/lib/python3.8/site-packages     # mappy for python 3.8
    $PY 65_trb_assembly_genotype.py --mode reads --limit 20 --workers 4 \\
        --out-dir ~/s05/results/65_reads_smoke --internal-dir ~/s05/internal/65_reads_smoke

Tests: scripts/hla_popgen/tests/test_65_trb_assembly_genotype.py (synthetic assemblies and reads).
Python 3.8 compatible.
"""
import argparse
import collections
import json
import gzip
import importlib.util
import os
import random
import re
import subprocess
import sys
import tempfile
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import _disclosure as _disc  # noqa: E402

try:
    import mappy  # noqa: E402
except ImportError:  # pragma: no cover - depends on the env
    mappy = None


def _load_58():
    """58_regional_depth.py (hwe_exact_p) -- module name starts with a digit, so load by path."""
    spec = importlib.util.spec_from_file_location("_rd58", os.path.join(_THIS_DIR, "58_regional_depth.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


hwe_exact_p = _load_58().hwe_exact_p

AOU_ROOT = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9")
LR_MANIFEST_REL = "v9/wgs/long_read/manifest.tsv"
DEFAULT_COHORT = os.path.expanduser("~/pipeline_outputs/cohort_membership.tsv")
DEFAULT_RELATEDNESS = os.path.join(AOU_ROOT, "v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv")
DEFAULT_PROBE_DIR = os.path.expanduser("~/s05/results/65")
DEFAULT_PROBES = os.path.join(DEFAULT_PROBE_DIR, "probes.fa")
RESULTS_ROOT = os.path.expanduser("~/s05/results")
INTERNAL_ROOT = os.path.expanduser("~/s05/internal")
DEFAULT_WORK_DIR = os.path.expanduser("~/s05/work65")

ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
INSERTION_GENES = ["TRBV4-3", "TRBV3-2", "TRBV6-3", "TRBV7-8"]
CONTROL_GENES = ["TRBV4-2", "TRBV6-2", "TRBV3-1", "TRBV2", "TRBV12-4", "TRBC1"]
DEFAULT_GENES = INSERTION_GENES + CONTROL_GENES
DEFAULT_ANCHORS = ["TRBV2", "TRBC1"]
DEFAULT_COPRESENCE = ["TRBV4-3", "TRBV3-2", "TRBV6-3"]
DEFAULT_WINDOW = "chr7:142250000-142860000"
IMGT_URL = "https://www.imgt.org/genedb/GENElect?query=7.2+TRBV&species=Homo+sapiens"
ENSEMBL = "https://rest.ensembl.org"
ANCHOR_PAD = 700_000
COPY2 = "[copy2]"
DEFAULT_READS_REGION = "chr7:142249011-142863287"  # IMGT TRB locus chr7:142,299,011-142,813,287 +/- 50 kb
# unique genes of the ~21.6 kb insertion (TRBV6-3 is merged with 6-2). TRBV7-8 is NOT here: on GRCh38's
# ALT it sits on a separate ALT-only segment ~180 kb downstream, so it needs its own breakpoint.
READS_INS_GENES = ["TRBV4-3", "TRBV3-2"]
BLOCK_MAX_GAP = 50_000
INS_ANY = "insertion_any"

Hit = collections.namedtuple("Hit", "probe contig t_st t_en strand q_st q_en q_len ident")


def log(msg):
    print(f"[65 {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(_RC)[::-1]


# ---------------------------------------------------------------------------
# FASTA / probes
# ---------------------------------------------------------------------------
def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def iter_fasta(path):
    """(name, seq) from .fa / .fa.gz. mappy.fastx_read (C) when available."""
    if mappy is not None:
        for name, seq, _ in mappy.fastx_read(path):
            yield name, seq.upper()
        return
    name, buf = None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf).upper()
                fields = line[1:].split()
                name, buf = (fields[0] if fields else ""), []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf).upper()


def gene_of(header):
    """'X12345|TRBV4-3*01|Homo sapiens|F|...' or 'TRBV4-3*01 ...' or 'TRBV4-3' -> ('TRBV4-3*01', 'TRBV4-3')."""
    h = header.strip().lstrip(">")
    fields = h.split("|")
    tok = fields[1] if len(fields) > 1 else h.split()[0]
    tok = tok.strip()
    return tok, tok.split("*")[0]


def read_probes(path):
    """-> list of (allele_name, gene, seq). Raw headers kept for gene parsing (IMGT '|' style)."""
    out, name, buf = [], None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:].rstrip("\n"), []
            else:
                buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    probes, seen = [], set()
    for hdr, seq in out:
        allele, gene = gene_of(hdr)
        seq = re.sub(r"[^ACGTN]", "", seq.upper())  # drops IMGT gap dots
        if not seq or allele in seen:
            continue
        seen.add(allele)
        probes.append((allele, gene, seq))
    return probes


def write_probes(probes, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        for allele, _gene, seq in probes:
            fh.write(f">{allele}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def _http(url, accept=None, timeout=120):
    req = urllib.request.Request(url, headers={"Accept": accept} if accept else {})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def parse_imgt_page(html):
    """GENElect pages wrap the FASTA in <pre>; also accept plain FASTA."""
    blocks = re.findall(r"<pre>(.*?)</pre>", html, flags=re.S | re.I) or [html]
    txt = "\n".join(re.sub(r"<[^>]+>", "", b) for b in blocks)
    probes, hdr, buf = [], None, []
    for line in txt.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if hdr:
                probes.append((hdr, "".join(buf)))
            hdr, buf = line[1:], []
        elif hdr:
            buf.append(line)
    if hdr:
        probes.append((hdr, "".join(buf)))
    out = []
    for h, s in probes:
        allele, gene = gene_of(h)
        s = re.sub(r"[^ACGTN]", "", s.upper())
        if s:
            out.append((allele, gene, s))
    return out


def fetch_probes_online(genes):
    v = parse_imgt_page(_http(IMGT_URL))
    v = [p for p in v if p[1].startswith("TRBV")]
    if not v:
        raise SystemExit("FATAL: IMGT returned no TRBV records -- check the page by hand or use --probe-fasta")
    out = list(v)
    for g in [g for g in genes if g.startswith("TRBC")]:
        gid = json.loads(_http(f"{ENSEMBL}/lookup/symbol/homo_sapiens/{g}?content-type=application/json",
                               accept="application/json"))["id"]
        fa = _http(f"{ENSEMBL}/sequence/id/{gid}?type=genomic;content-type=text/x-fasta")
        seq = "".join(l.strip() for l in fa.splitlines() if not l.startswith(">"))
        out.append((f"{g}*genomic_{gid}", g, seq.upper()))
    return out


# ---------------------------------------------------------------------------
# Alignment primitives
# ---------------------------------------------------------------------------
def semiglobal_ed(q, t):
    """Edit distance aligning all of q against any substring of t (free end gaps in t)."""
    if not q:
        return 0
    if not t:
        return len(q)
    qa = np.frombuffer(q.encode(), dtype=np.uint8)
    ta = np.frombuffer(t.encode(), dtype=np.uint8)
    j = np.arange(len(ta) + 1)
    prev = np.zeros(len(ta) + 1, dtype=np.int64)
    for i in range(1, len(qa) + 1):
        cur = np.empty_like(prev)
        cur[0] = i
        cur[1:] = np.minimum(prev[:-1] + (ta != qa[i - 1]), prev[1:] + 1)
        cur = np.minimum.accumulate(cur - j) + j
        prev = cur
    return int(prev.min())


def pair_identity(a, b):
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    return 1.0 - semiglobal_ed(q, t) / max(1, len(q))


def _kmers(s, k):
    return {s[i:i + k] for i in range(0, len(s) - k + 1)}


def gene_classes(probes, merge_ident=0.99, k=15):
    """Union genes whose closest cross-gene allele pair is >= merge_ident identical.
    Returns (gene -> class name, DataFrame gene, class, max_cross_identity, closest_gene)."""
    genes = sorted({g for _, g, _ in probes}, key=_natkey)
    alle = collections.defaultdict(list)
    for _, g, s in probes:
        alle[g].append(s)
    km = {g: [_kmers(s, k) | _kmers(revcomp(s), k) for s in alle[g]] for g in genes}
    parent = {g: g for g in genes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    best = {g: (0.0, "") for g in genes}
    for i, g1 in enumerate(genes):
        for g2 in genes[i + 1:]:
            mx = 0.0
            for s1, k1 in zip(alle[g1], km[g1]):
                for s2, k2 in zip(alle[g2], km[g2]):
                    sh = len(k1 & k2) / max(1, min(len(k1), len(k2)))
                    if sh < 0.3:  # cannot reach ~0.99 identity
                        continue
                    idn = max(pair_identity(s1, s2), pair_identity(s1, revcomp(s2)))
                    mx = max(mx, idn)
            if mx > best[g1][0]:
                best[g1] = (mx, g2)
            if mx > best[g2][0]:
                best[g2] = (mx, g1)
            if mx >= merge_ident:
                parent[find(g1)] = find(g2)
    members = collections.defaultdict(list)
    for g in genes:
        members[find(g)].append(g)
    cls = {}
    for root, ms in members.items():
        name = "/".join(sorted(ms, key=_natkey))
        for g in ms:
            cls[g] = name
    tab = pd.DataFrame([dict(gene=g, gene_class=cls[g], max_cross_gene_identity=round(best[g][0], 4),
                             closest_gene=best[g][1]) for g in genes])
    return cls, tab


def _natkey(g):
    return [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", g)]


# ---- aligners: fn(probes[(name, seq)], targets{contig: seq}, opts) -> [Hit] ------------------
def align_mappy(probes, targets, opts):
    """Probes (query) vs targets. Few targets: one in-memory index each. Many targets (reads mode,
    opts['single_index']): one index over a temp FASTA in the work dir, removed once loaded."""
    if mappy is None:
        raise RuntimeError("mappy not importable")
    preset, best_n = opts.get("preset", "sr"), opts.get("best_n", 50)
    hits = []
    if opts.get("single_index") and len(targets) > 1:
        work = opts.get("work_dir") or None
        if work:
            os.makedirs(work, exist_ok=True)
        fd, path = tempfile.mkstemp(suffix=".fa", dir=work)
        try:
            with os.fdopen(fd, "w") as fh:
                for c, s in targets.items():
                    fh.write(f">{c}\n{s}\n")
            a = mappy.Aligner(fn_idx_in=path, preset=preset, best_n=best_n)
        finally:
            os.remove(path)
        if not a:
            raise RuntimeError("mappy index failed")
        for name, seq in probes:
            for h in a.map(seq):
                hits.append(Hit(name, h.ctg, h.r_st, h.r_en, h.strand, h.q_st, h.q_en, len(seq),
                                h.mlen / max(1, h.blen)))
        return hits
    for ctg, tseq in targets.items():
        a = mappy.Aligner(seq=tseq, preset=preset, best_n=best_n)
        if not a:
            raise RuntimeError("mappy index failed")
        for name, seq in probes:
            for h in a.map(seq):
                hits.append(Hit(name, ctg, h.r_st, h.r_en, h.strand, h.q_st, h.q_en, len(seq),
                                h.mlen / max(1, h.blen)))
    return hits


def parse_paf_hits(text):
    hits = []
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        hits.append(Hit(f[0], f[5], int(f[7]), int(f[8]), 1 if f[4] == "+" else -1, int(f[2]), int(f[3]),
                        int(f[1]), int(f[9]) / max(1, int(f[10]))))
    return hits


def align_minimap2(probes, targets, opts):
    work = opts.get("work_dir") or None
    if work:
        os.makedirs(work, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=work) as td:
        tfa, qfa = os.path.join(td, "t.fa"), os.path.join(td, "q.fa")
        with open(tfa, "w") as fh:
            for c, s in targets.items():
                fh.write(f">{c}\n{s}\n")
        with open(qfa, "w") as fh:
            for n, s in probes:
                fh.write(f">{n}\n{s}\n")
        cmd = [opts.get("minimap2", "minimap2"), "-x", opts.get("preset", "sr"), "-c", "--secondary=yes",
               "-N", str(opts.get("best_n", 50)), "-p", "0.5", "-t", "1", tfa, qfa]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=opts.get("timeout"), check=True)
    return parse_paf_hits(r.stdout)


def align_kmer(probes, targets, opts, k=15, min_seeds=4):
    """Pure-Python fallback: exact k-mer seeds -> diagonal clusters -> identity."""
    idx = collections.defaultdict(list)
    for i, (_, seq) in enumerate(probes):
        for strand, s in ((1, seq), (-1, revcomp(seq))):
            for o in range(len(s) - k + 1):
                idx[s[o:o + k]].append((i, strand, o))
    hits = []
    for ctg, t in targets.items():
        diags = collections.defaultdict(list)
        get = idx.get
        for p in range(len(t) - k + 1):
            lst = get(t[p:p + k])
            if lst:
                for i, strand, o in lst:
                    diags[(i, strand)].append(p - o)
        for (i, strand), ds in diags.items():
            ds.sort()
            groups, cur = [], [ds[0]]
            for d in ds[1:]:
                if d - cur[-1] <= 50:
                    cur.append(d)
                else:
                    groups.append(cur)
                    cur = [d]
            groups.append(cur)
            name, seq = probes[i]
            s = seq if strand == 1 else revcomp(seq)
            L = len(s)
            for g in groups:
                if len(g) < min_seeds:
                    continue
                d = int(np.median(g))
                a, b = max(0, d), min(len(t), d + L)
                if b - a <= k:
                    continue
                qa, qb = a - d, b - d  # probe sub-range inside the contig
                sub = s[qa:qb]
                win = t[a:b]
                ham = sum(1 for x, y in zip(sub, win) if x != y)
                ident = 1.0 - ham / max(1, len(sub))
                if ident < 0.99:
                    lo, hi = max(0, a - 30), min(len(t), b + 30)
                    ident = max(ident, 1.0 - semiglobal_ed(sub, t[lo:hi]) / max(1, len(sub)))
                if strand == 1:
                    q_st, q_en = qa, qb
                else:
                    q_st, q_en = L - qb, L - qa
                hits.append(Hit(name, ctg, a, b, strand, q_st, q_en, L, ident))
    return hits


ALIGNERS = {"mappy": align_mappy, "minimap2": align_minimap2, "kmer": align_kmer}


def resolve_aligner(name, minimap2_path):
    if name == "auto":
        if mappy is not None:
            name = "mappy"
        elif minimap2_path and _which(minimap2_path):
            name = "minimap2"
        else:
            name = "kmer"
    if name == "mappy" and mappy is None:
        raise SystemExit("FATAL: --aligner mappy but mappy is not importable (pip install mappy)")
    if name == "minimap2" and not _which(minimap2_path):
        raise SystemExit(f"FATAL: minimap2 not found at {minimap2_path!r}")
    return name


def _which(p):
    if not p:
        return None
    if os.path.sep in p:
        return p if os.access(p, os.X_OK) else None
    for d in os.environ.get("PATH", "").split(os.pathsep):
        c = os.path.join(d, p)
        if os.access(c, os.X_OK):
            return c
    return None


# ---------------------------------------------------------------------------
# Loci, assignment, haplotype genotype
# ---------------------------------------------------------------------------
def cluster_loci(hits):
    """Group hits into loci: same contig, overlap >= 50% of the shorter interval."""
    loci = []
    by_ctg = collections.defaultdict(list)
    for h in hits:
        by_ctg[h.contig].append(h)
    for ctg, hs in by_ctg.items():
        hs.sort(key=lambda h: (h.t_st, h.t_en))
        cur, cs, ce = [], None, None
        for h in hs:
            if cur:
                ov = min(ce, h.t_en) - max(cs, h.t_st)
                if ov >= 0.5 * min(ce - cs, h.t_en - h.t_st):
                    cur.append(h)
                    ce = max(ce, h.t_en)
                    continue
                loci.append((ctg, cs, ce, cur))
            cur, cs, ce = [h], h.t_st, h.t_en
        if cur:
            loci.append((ctg, cs, ce, cur))
    return loci


def assign_loci(loci, probe_gene, cls, min_ident, min_cov, tie_eps):
    """-> list of dict(contig, start, end, gene_class, ident, cov, called, ambiguous)."""
    out = []
    for ctg, s, e, hs in loci:
        best = {}
        for h in hs:
            c = cls[probe_gene[h.probe]]
            cov = (h.q_en - h.q_st) / max(1, h.q_len)
            key = (cov >= min_cov, round(h.ident, 6), cov)
            if c not in best or key > best[c]:
                best[c] = key
        ranked = sorted(best.items(), key=lambda kv: kv[1], reverse=True)
        top_c, (ok_cov, idn, cov) = ranked[0]
        amb = len(ranked) > 1 and ranked[1][1][0] == ok_cov and abs(ranked[1][1][1] - idn) <= tie_eps
        called = bool(ok_cov and idn >= min_ident and not amb)
        out.append(dict(contig=ctg, start=s, end=e, gene_class=top_c, ident=idn, cov=cov,
                        called=called, ambiguous=bool(amb and ok_cov and idn >= min_ident),
                        runner_up=ranked[1][0] if len(ranked) > 1 else ""))
    return out


def summarise_haplotype(assigned, classes_reported, anchor_classes):
    """-> (calls, complete, complete_single_contig, reason). Presence on any TRB contig; complete =
    every anchor called somewhere among them; single-contig = all anchors on one contig."""
    calls = {}
    for c in classes_reported:
        loc = [a for a in assigned if a["called"] and a["gene_class"] == c]
        allhits = [a for a in assigned if a["gene_class"] == c]
        calls[c] = dict(present=len(loc) > 0, n_loci=len(loc),
                        best_ident=max([a["ident"] for a in allhits], default=float("nan")),
                        contigs=sorted({a["contig"] for a in loc}))
    anchor_ctgs = [set(a["contig"] for a in assigned if a["called"] and a["gene_class"] == c)
                   for c in anchor_classes]
    complete = all(anchor_ctgs)
    single = bool(complete and set.intersection(*anchor_ctgs))
    return calls, complete, single, "" if complete else "anchor_missing"


# ---------------------------------------------------------------------------
# Contig selection
# ---------------------------------------------------------------------------
def parse_window(w):
    m = re.match(r"^([^:]+):([\d,]+)-([\d,]+)$", w)
    if not m:
        raise ValueError(f"bad region {w!r}")
    return m.group(1), int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))


def paf_trb_spans(paf, window, pad, minlen=1000):
    """contig -> (start0, end) on the contig covering the projected window +/- pad."""
    chrom, ws, we = window
    spans = {}
    with _open(paf) as fh:
        for ln in fh:
            x = ln.rstrip("\n").split("\t")
            if len(x) < 12 or x[5] != chrom or "tp:A:P" not in ln:
                continue
            ql, qs, qe, strand, ts, te = int(x[1]), int(x[2]), int(x[3]), x[4], int(x[7]), int(x[8])
            a, b = max(ts, ws), min(te, we)
            if b - a < minlen:
                continue
            f0, f1 = (a - ts) / max(1, te - ts), (b - ts) / max(1, te - ts)
            if strand == "+":
                q0, q1 = qs + f0 * (qe - qs), qs + f1 * (qe - qs)
            else:
                q0, q1 = qe - f1 * (qe - qs), qe - f0 * (qe - qs)
            lo, hi = max(0, int(q0) - pad), min(ql, int(q1) + pad)
            if x[0] in spans:
                lo, hi = min(lo, spans[x[0]][0]), max(hi, spans[x[0]][1])
            spans[x[0]] = (lo, hi)
    return spans


def anchor_kmers(probes, genes, k=25, per_probe=3):
    out = set()
    first = {}
    for allele, g, s in probes:
        if g in genes and g not in first and len(s) >= k:
            first[g] = s
    for s in first.values():
        for i in np.linspace(0, len(s) - k, per_probe).astype(int):
            km = s[i:i + k]
            out.add(km)
            out.add(revcomp(km))
    return sorted(out)


def select_contigs(fa, paf, mode, window, pad, anchors, deadline=None):
    """-> (targets {name: seq}, mode_used). Streams the .fa.gz once."""
    paf = paf if isinstance(paf, str) and paf else None  # NaN from a DataFrame row
    if mode == "auto":
        mode = "paf" if paf and os.path.exists(paf) else "anchor"
    targets = {}
    if mode == "paf":
        spans = paf_trb_spans(paf, window, pad)
        if not spans:
            return {}, mode
        need = set(spans)
        for name, seq in iter_fasta(fa):
            if name in need:
                a, b = spans[name]
                targets[name] = seq[a:b]
                need.discard(name)
                if not need:
                    break
            _check(deadline)
        return targets, mode
    for name, seq in iter_fasta(fa):
        _check(deadline)
        if mode == "all":
            targets[name] = seq
            continue
        pos = [p for p in (seq.find(a) for a in anchors) if p >= 0]
        if len(pos) >= 2:
            targets[name] = seq[max(0, min(pos) - ANCHOR_PAD): max(pos) + ANCHOR_PAD]
    return targets, mode


class _Timeout(Exception):
    pass


def _check(deadline):
    if deadline is not None and time.time() > deadline:
        raise _Timeout()


def _opts(args_d, deadline=None, **kw):
    o = dict(preset=args_d["preset"], minimap2=args_d["minimap2"], work_dir=args_d.get("work_dir"),
             timeout=max(1, (deadline - time.time())) if deadline else None)
    o.update(kw)
    return o


def genotype_targets(targets, probes, cls, reported_classes, anchor_classes, args_d, deadline=None):
    """Probe-map a set of sequences treated as ONE haplotype -> dict(calls, complete, ..., assigned)."""
    if not targets:
        calls = {c: dict(present=False, n_loci=0, best_ident=float("nan"), contigs=[]) for c in reported_classes}
        return dict(calls=calls, complete=False, complete_single=False, reason="no_trb_contig",
                    n_called=0, n_ambiguous=0, assigned=[])
    hits = ALIGNERS[args_d["aligner"]]([(a, s) for a, _, s in probes], targets, _opts(args_d, deadline))
    probe_gene = {a: g for a, g, _ in probes}
    assigned = assign_loci(cluster_loci(hits), probe_gene, cls, args_d["min_ident"], args_d["min_cov"],
                           args_d["tie_eps"])
    calls, complete, single, reason = summarise_haplotype(assigned, reported_classes, anchor_classes)
    return dict(calls=calls, complete=complete, complete_single=single, reason=reason, assigned=assigned,
                n_called=sum(a["called"] for a in assigned), n_ambiguous=sum(a["ambiguous"] for a in assigned))


def genotype_haplotype(fa, paf, probes, cls, reported_classes, anchor_classes, args_d, deadline=None):
    targets, mode = select_contigs(fa, paf, args_d["select"], args_d["window"], args_d["pad"],
                                   args_d["anchors_kmers"], deadline)
    _check(deadline)
    r = genotype_targets(targets, probes, cls, reported_classes, anchor_classes, args_d, deadline)
    r.update(n_contigs=len(targets), bp=sum(len(s) for s in targets.values()), select=mode)
    return r


def features_for(reported_classes):
    """Reported features: every class, plus '<class>[copy2]' for merged classes."""
    feats = []
    for c in reported_classes:
        feats.append(c)
        if "/" in c:
            feats.append(c + COPY2)
    return feats


def feature_of_gene(gene, cls):
    c = cls.get(gene, gene)
    return c + COPY2 if "/" in c else c


def hap_feature_present(calls, feat):
    if feat.endswith(COPY2):
        return calls[feat[:-len(COPY2)]]["n_loci"] >= 2
    return calls[feat]["present"]


def hap_feature_nloci(calls, feat):
    c = feat[:-len(COPY2)] if feat.endswith(COPY2) else feat
    return calls[c]["n_loci"]


GENE_COLS = ["person_id", "feature", "hap1_present", "hap2_present", "hap1_n_loci", "hap2_n_loci", "copies",
             "hap1_complete", "hap2_complete", "status", "best_ident", "method", "n_ins", "n_ref",
             "allele_fraction", "hap1_n_ins", "hap1_n_ref", "hap2_n_ins", "hap2_n_ref"]
HAP_COLS = ["person_id", "hap", "select", "n_trb_contigs", "trb_bp", "complete", "complete_single_contig",
            "incomplete_reason", "n_loci_called", "n_loci_ambiguous", "n_reads", "n_ins", "n_ref", "method",
            "seconds"]
FAIL_COLS = ["person_id", "reason"]


def _row(cols, **kw):
    return {c: kw.get(c, "") for c in cols}


def genotype_person(person, probes, cls, reported_classes, anchor_classes, args_d, timeout=None):
    """Assembly mode. -> (gene_rows, hap_rows, None) or (None, None, reason)."""
    t0 = time.time()
    deadline = t0 + timeout if timeout else None
    haps = {}
    for h in (1, 2):
        fa = person.get(f"hap{h}_fa")
        if not isinstance(fa, str) or not fa or not os.path.exists(fa):
            return None, None, "assembly_missing"
    method = f"{args_d['aligner']}:{args_d['preset']}" if args_d["aligner"] != "kmer" else "kmer"
    for h in (1, 2):
        th = time.time()
        try:
            r = genotype_haplotype(person[f"hap{h}_fa"], person.get(f"hap{h}_paf"), probes, cls,
                                   reported_classes, anchor_classes, args_d, deadline)
        except (_Timeout, subprocess.TimeoutExpired):
            return None, None, "timeout"
        except subprocess.CalledProcessError:
            return None, None, "aligner_error"
        except (OSError, EOFError, ValueError, RuntimeError) as ex:
            return None, None, f"read_error:{type(ex).__name__}"
        r["seconds"] = round(time.time() - th, 1)
        haps[h] = r
    pid = person["person_id"]
    hap_rows = [_row(HAP_COLS, person_id=pid, hap=h, select=r["select"], n_trb_contigs=r["n_contigs"],
                     trb_bp=r["bp"], complete=r["complete"], complete_single_contig=r["complete_single"],
                     incomplete_reason=r["reason"] or "", n_loci_called=r["n_called"],
                     n_loci_ambiguous=r["n_ambiguous"], method="assembly:" + method, seconds=r["seconds"])
                for h, r in haps.items()]
    gene_rows = []
    for feat in features_for(reported_classes):
        p1, p2 = (hap_feature_present(haps[h]["calls"], feat) for h in (1, 2))
        c = feat[:-len(COPY2)] if feat.endswith(COPY2) else feat
        idents = [haps[h]["calls"][c]["best_ident"] for h in (1, 2)]
        idents = [x for x in idents if np.isfinite(x)]
        gene_rows.append(_row(GENE_COLS, person_id=pid, feature=feat, hap1_present=int(p1), hap2_present=int(p2),
                              hap1_n_loci=hap_feature_nloci(haps[1]["calls"], feat),
                              hap2_n_loci=hap_feature_nloci(haps[2]["calls"], feat),
                              copies=int(p1) + int(p2), hap1_complete=int(haps[1]["complete"]),
                              hap2_complete=int(haps[2]["complete"]),
                              status="ok" if haps[1]["complete"] and haps[2]["complete"] else "incomplete",
                              best_ident=round(max(idents), 4) if idents else "", method="assembly:" + method))
    return gene_rows, hap_rows, None


# ---------------------------------------------------------------------------
# Reads mode
# ---------------------------------------------------------------------------
def run_samtools_view(bam, region, timeout, samtools="samtools"):
    """SAM text of primary, mapped reads over region (patchable in tests, as 58's run_bedcov)."""
    cmd = [samtools, "view", "-F", "0x904", bam, region]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=True).stdout


_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_sam(text):
    reads = []
    for line in text.splitlines():
        if not line or line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 11:
            continue
        hp = 0
        for t in f[11:]:
            if t.startswith("HP:i:"):
                hp = int(t[5:])
        reads.append(dict(name=f[0], flag=int(f[1]), pos=int(f[3]), cigar=f[5], seq=f[9].upper(), hp=hp))
    return reads


def spans_contiguously(pos, cigar, lo, hi, max_indel):
    """True if the alignment (1-based pos) covers ref [lo, hi] with no I/D/N >= max_indel inside it."""
    if cigar == "*" or not cigar:
        return False
    r = pos
    for n, op in _CIG.findall(cigar):
        n = int(n)
        if op in "M=X":
            r += n
        elif op in "DN":
            if n >= max_indel and r <= hi and r + n - 1 >= lo:
                return False
            r += n
        elif op == "I":
            if n >= max_indel and lo <= r <= hi:
                return False
    return pos <= lo and r - 1 >= hi


def call_copies(af, cuts):
    return 0 if af < cuts[0] else (1 if af < cuts[1] else 2)


def rederive_reads_calls(gene, min_reads, cuts):
    """Reads mode: recompute status/copies from the stored n_ins/n_ref with the current --min-reads and
    --af-cuts, so the PUBLIC tables can be rebuilt with --aggregate-only at another floor."""
    g = gene.copy()
    ni, nr = pd.to_numeric(g["n_ins"], errors="coerce"), pd.to_numeric(g["n_ref"], errors="coerce")
    tot = ni + nr
    ok = tot >= min_reads
    af = ni / tot.where(tot > 0)
    g["status"] = np.where(ok, "ok", "low_reads")
    g["copies"] = [call_copies(a, cuts) if o else "" for a, o in zip(af, ok)]
    return g


def genotype_person_reads(person, probes, cls, args_d, timeout=None):
    """Reads mode. -> (gene_rows, hap_rows, None) or (None, None, reason)."""
    t0 = time.time()
    bam = person.get("bam")
    if not isinstance(bam, str) or not bam or not os.path.exists(bam):
        return None, None, "bam_missing"
    try:
        text = run_samtools_view(bam, args_d["reads_region"], timeout, args_d["samtools"])
    except subprocess.TimeoutExpired:
        return None, None, "timeout"
    except subprocess.CalledProcessError:
        return None, None, "samtools_error"
    except OSError as ex:
        return None, None, f"samtools_error:{type(ex).__name__}"
    seen, reads = set(), []
    for r in parse_sam(text):
        if r["flag"] & 0x904 or r["seq"] == "*" or r["name"] in seen:
            continue
        seen.add(r["name"])
        reads.append(r)
    _chrom, bl, br = args_d["breakpoint"]
    lo, hi = bl - args_d["flank"], br + args_d["flank"]
    ref_names = {r["name"] for r in reads if spans_contiguously(r["pos"], r["cigar"], lo, hi, args_d["max_indel"])}
    targets = {r["name"]: r["seq"] for r in reads if r["name"] not in ref_names}
    deadline = t0 + timeout if timeout else None
    try:
        hits = ALIGNERS[args_d["aligner"]]([(a, s) for a, _, s in probes], targets,
                                           _opts(args_d, deadline, single_index=True, best_n=5000)) if targets else []
    except subprocess.TimeoutExpired:
        return None, None, "timeout"
    except (subprocess.CalledProcessError, RuntimeError):
        return None, None, "aligner_error"
    probe_gene = {a: g for a, g, _ in probes}
    assigned = assign_loci(cluster_loci(hits), probe_gene, cls, args_d["min_ident"], args_d["min_cov"],
                           args_d["tie_eps"])
    ins_classes = args_d["ins_classes"]
    read_genes, best_id = collections.defaultdict(set), {}
    for a in assigned:
        if a["called"] and a["gene_class"] in ins_classes:
            read_genes[a["contig"]].add(a["gene_class"])
            best_id[a["gene_class"]] = max(best_id.get(a["gene_class"], 0.0), a["ident"])
    hp = {r["name"]: r["hp"] for r in reads}
    has_hp = any(v in (1, 2) for v in hp.values())

    def count(feat, h=None):
        ins = sum(1 for n, gs in read_genes.items() if (feat == INS_ANY or feat in gs) and (h is None or hp[n] == h))
        ref = sum(1 for n in ref_names if h is None or hp[n] == h)
        return ins, ref

    cuts, method = args_d["af_cuts"], f"reads:{args_d['aligner']}"
    hap_ok = {}
    hap_rows = []
    for h in (1, 2):
        ni, nr = count(INS_ANY, h)
        ok = has_hp and ni + nr >= args_d["min_hap_reads"]
        hap_ok[h] = ok
        hap_rows.append(_row(HAP_COLS, person_id=person["person_id"], hap=h, select="reads", complete=ok,
                             incomplete_reason="" if ok else ("low_reads" if has_hp else "no_hp_tag"),
                             n_reads=sum(1 for v in hp.values() if v == h), n_ins=ni, n_ref=nr, method=method,
                             seconds=round(time.time() - t0, 1)))
    gene_rows = []
    for feat in list(ins_classes) + [INS_ANY]:
        ni, nr = count(feat)
        tot = ni + nr
        af = ni / tot if tot else float("nan")
        hc = {}
        pres = {}
        for h in (1, 2):
            hi_, hr_ = count(feat, h)
            hc[h] = (hi_, hr_)
            pres[h] = int(hap_ok[h] and hi_ + hr_ > 0 and hi_ / (hi_ + hr_) >= 0.5)
        ok = tot >= args_d["min_reads"]
        bi = max([best_id.get(c, 0.0) for c in ins_classes] if feat == INS_ANY else [best_id.get(feat, 0.0)])
        gene_rows.append(_row(GENE_COLS, person_id=person["person_id"], feature=feat,
                              hap1_present=pres[1], hap2_present=pres[2],
                              copies=call_copies(af, cuts) if ok else "",
                              hap1_complete=int(hap_ok[1]), hap2_complete=int(hap_ok[2]),
                              status="ok" if ok else "low_reads", best_ident=round(bi, 4) if bi else "",
                              method=method, n_ins=ni, n_ref=nr,
                              allele_fraction=round(af, 4) if tot else "",
                              hap1_n_ins=hc[1][0], hap1_n_ref=hc[1][1], hap2_n_ins=hc[2][0], hap2_n_ref=hc[2][1]))
    return gene_rows, hap_rows, None


# ---------------------------------------------------------------------------
# Breakpoint (probe build)
# ---------------------------------------------------------------------------
def load_ref(path, region=None):
    """-> (chrom, offset0, seq). With region 'chr7:S-E' (1-based) cut it from contig chrom; else the
    first record, a samtools-faidx header 'chr7:S-E' giving the offset."""
    if region:
        chrom, s, e = parse_window(region)
        for name, seq in iter_fasta(path):
            if name == chrom:
                return chrom, s - 1, seq[s - 1:e]
        raise SystemExit(f"FATAL: {chrom} not found in {path}")
    for name, seq in iter_fasta(path):
        m = re.match(r"^(.+):([\d,]+)-([\d,]+)$", name)
        if m:
            return m.group(1), int(m.group(2).replace(",", "")) - 1, seq
        return name, 0, seq
    raise SystemExit(f"FATAL: no sequence in {path}")


def find_breakpoint(alt, ins_lo, ins_hi, ref, aligner, opts, span=40000, chunk=1000, step=250, tol=50):
    """ALT flank chunks -> colinear offsets on ref -> breakpoint interval (0-based ref coords)."""
    def chunks(lo, hi, tag):
        lo = max(0, lo)
        return [(f"{tag}{s}", alt[s:s + chunk]) for s in range(lo, max(lo, hi - chunk) + 1, step)
                if len(alt[s:s + chunk]) == chunk]
    q = chunks(ins_lo - span, ins_lo, "L") + chunks(ins_hi, min(len(alt), ins_hi + span), "R")
    hits = ALIGNERS[aligner](q, {"ref": ref}, opts)
    by = collections.defaultdict(list)
    for h in hits:
        if h.q_en - h.q_st >= 0.9 * h.q_len and h.ident >= 0.98:
            by[h.probe].append(h)
    pts = {"L": [], "R": []}
    for name, hs in by.items():
        if len(hs) == 1 and hs[0].strand == 1:
            h, s = hs[0], int(name[1:])
            pts[name[0]].append((s, h.t_st - (s + h.q_st), h.t_st, h.t_en))
    if len(pts["L"]) < 3 or len(pts["R"]) < 3:
        return None
    # Offsets come from the colinear block NEXT TO the insertion genes (innermost chunks), not from the
    # outer flank ends: the GRCh38 ALT differs from primary by further indels within the flank span,
    # so the outermost chunks sit on other offset blocks and give a too-wide breakpoint.
    inner_l = sorted(pts["L"], reverse=True)   # nearest the insertion first
    inner_r = sorted(pts["R"])
    off_l, blk_l = _inner_block(inner_l, tol)
    off_r, blk_r = _inner_block(inner_r, tol)
    if off_l is None or off_r is None:
        return None
    bp_l = max(te for _, _, _, te in blk_l)
    bp_r = min(ts for _, _, ts, _ in blk_r)
    return dict(lo0=min(bp_l, bp_r), hi0=max(bp_l, bp_r), ins_len=int(off_l - off_r),
                n_left=len(blk_l), n_right=len(blk_r))


def _inner_block(pts, tol, min_chunks=3):
    """pts ordered from the insertion outward -> (offset, chunks) of the first run of >= min_chunks
    consecutive chunks sharing one offset (within tol); (None, []) if there is none."""
    i = 0
    while i < len(pts):
        o = pts[i][1]
        j = i
        while j < len(pts) and abs(pts[j][1] - o) <= tol:
            j += 1
        if j - i >= min_chunks:
            blk = pts[i:j]
            return float(np.median([x[1] for x in blk])), blk
        i = j if j > i else i + 1
    return None, []


def breakpoint_from_reference(alt_seqs, ref_chrom, ref_off, ref_seq, assigned_alt, ins_classes, aligner, opts):
    """Locate the insertion genes on the ALT, then find_breakpoint (retrying with the ALT reversed)."""
    loc = [a for a in assigned_alt if a["called"] and a["gene_class"] in ins_classes]
    if not loc:
        return None
    ctg = collections.Counter(a["contig"] for a in loc).most_common(1)[0][0]
    alt = alt_seqs[ctg]
    block, elsewhere = insertion_block(loc, ctg, ins_classes)
    lo = min(a["start"] for a in block)
    hi = max(a["end"] for a in block)
    bp = find_breakpoint(alt, lo, hi, ref_seq, aligner, opts)
    if bp is None:
        bp = find_breakpoint(revcomp(alt), len(alt) - hi, len(alt) - lo, ref_seq, aligner, opts)
    if bp is None:
        return None
    # 1-based: last base of the left flank .. first base of the right flank
    bp.update(chrom=ref_chrom, left=ref_off + bp["lo0"], right=ref_off + bp["hi0"] + 1, alt_contig=ctg,
              block_classes=sorted({a["gene_class"] for a in block}, key=_natkey),
              classes_elsewhere=sorted(elsewhere, key=_natkey))
    return bp


def insertion_block(loc, ctg, ins_classes, max_gap=BLOCK_MAX_GAP):
    """Insertion-gene loci on the ALT contig, grouped into blocks (a gap > max_gap starts a new one).
    GRCh38's KI270803v1 carries TWO ALT-only segments relative to primary: the ~21.6 kb insertion
    (TRBV3-2, TRBV4-3, second TRBV6-2/6-3 copy, between TRBV4-2 and TRBV6-2) and a separate ~40 kb
    segment ~180 kb downstream holding TRBV7-8. One breakpoint can only describe one of them, so the
    block holding the lead class (the first of ins_classes, TRBV4-3) is used.
    -> (loci of that block, classes found only outside it)."""
    on = sorted((a for a in loc if a["contig"] == ctg), key=lambda a: a["start"])
    blocks, cur = [], [on[0]]
    for a in on[1:]:
        if a["start"] - max(x["end"] for x in cur) > max_gap:
            blocks.append(cur)
            cur = []
        cur.append(a)
    blocks.append(cur)
    lead = ins_classes[0]
    block = next((b for b in blocks if any(a["gene_class"] == lead for a in b)), max(blocks, key=len))
    inside = {a["gene_class"] for a in block}
    elsewhere = {a["gene_class"] for b in blocks if b is not block for a in b} - inside
    return block, elsewhere


def read_probe_meta(path):
    if not path or not os.path.exists(path):
        return {}
    d = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    return dict(zip(d["key"], d["value"]))


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def uri_to_local(uri, root=AOU_ROOT):
    """gs://<bucket>/path -> root/path (57/58's rule); a relative path -> root/path; else None."""
    if not isinstance(uri, str) or not uri or uri == "NA":
        return None
    if uri.startswith("gs://"):
        parts = uri.split("/", 3)
        return os.path.join(root, parts[3]) if len(parts) > 3 else None
    return uri if os.path.isabs(uri) else os.path.join(root, uri)


ASM_COLS = (("hap1_fa", "assembly_hap1_fa"), ("hap2_fa", "assembly_hap2_fa"),
            ("hap1_paf", "assembly_hap1_aln2_hg38_paf"), ("hap2_paf", "assembly_hap2_aln2_hg38_paf"))


def load_people(cohort_path, manifest_path, root=AOU_ROOT, bam_column="grch38_bam"):
    """One row per long-read person: person_id, ancestry, platform, bam, hap1/2_fa, hap1/2_paf.
    Manifest rows matching the cohort platform are preferred. No existence checks (slow on the mount)."""
    coh = pd.read_csv(cohort_path, sep="\t", dtype={"person_id": str})
    if "in_lr" in coh.columns:
        coh = coh[coh["in_lr"].astype(str).str.lower().isin(["true", "1"])]
    coh["ancestry"] = coh.get("ancestry_pred", pd.Series(index=coh.index, dtype=str)).astype(str).str.upper()
    man = pd.read_csv(manifest_path, sep="\t", dtype=str).rename(columns={"research_id": "person_id"})
    for c in [c for _, c in ASM_COLS if c not in man.columns] + ([bam_column] if bam_column not in man.columns else []):
        man[c] = None
    by_pid = {pid: g for pid, g in man.groupby("person_id")}
    rows = []
    for _, r in coh.iterrows():
        g = by_pid.get(r["person_id"])
        rec = dict(person_id=r["person_id"], ancestry=r["ancestry"],
                   platform=str(r.get("platform", "")).lower(), bam=None)
        if g is not None:
            cand = pd.concat([g[g["platform"].astype(str).str.lower() == rec["platform"]], g])
            for _, m in cand.iterrows():
                if isinstance(m[bam_column], str):
                    rec["bam"] = uri_to_local(m[bam_column], root)
                    break
            for _, m in cand.iterrows():
                if isinstance(m["assembly_hap1_fa"], str) and isinstance(m["assembly_hap2_fa"], str):
                    for k, c in ASM_COLS:
                        rec[k] = uri_to_local(m[c], root)
                    break
            if rec["platform"] in ("", "nan", "none"):
                rec["platform"] = str(g["platform"].iloc[0]).lower()
        for k, _ in ASM_COLS:
            rec.setdefault(k, None)
        rows.append(rec)
    return pd.DataFrame(rows)


def read_people_file(path):
    return set(pd.read_csv(path, sep="\t", dtype={"person_id": str}, comment="#")["person_id"])


def sample_stratified(people, n, seed, unrelated_ids, need="assembly"):
    """Unrelated, ancestry-balanced, seeded sample. Equal shares per ancestry; a group smaller than
    its share gives all it has and the remainder is re-split among the others."""
    p = people[people["person_id"].isin(unrelated_ids)]
    if need == "assembly":
        p = p[p["hap1_fa"].notna() & p["hap2_fa"].notna()]
    rng = random.Random(seed)
    groups = {}
    for a in ANCESTRY_ORDER:
        ids = sorted(p.loc[p["ancestry"] == a, "person_id"])
        rng.shuffle(ids)
        if ids:
            groups[a] = ids
    alloc = {a: 0 for a in groups}
    remaining = n
    while remaining > 0:
        open_ = [a for a in groups if alloc[a] < len(groups[a])]
        if not open_:
            break
        share = max(1, remaining // len(open_))
        for a in open_:
            take = min(share, len(groups[a]) - alloc[a], remaining)
            alloc[a] += take
            remaining -= take
            if remaining == 0:
                break
    pick = [pid for a in groups for pid in groups[a][:alloc[a]]]
    out = p[p["person_id"].isin(set(pick))][["person_id", "ancestry", "platform"]]
    return out.sort_values(["ancestry", "person_id"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Checkpoint + lock
# ---------------------------------------------------------------------------
def acquire_lock(internal_dir):
    os.makedirs(internal_dir, exist_ok=True)
    lock = os.path.join(internal_dir, "65.lock")
    if os.path.exists(lock):
        try:
            pid = int(open(lock).read().strip())
            os.kill(pid, 0)
            raise SystemExit(f"FATAL: another 65_trb_assembly_genotype.py (pid {pid}) holds {lock}")
        except (ValueError, ProcessLookupError):
            log(f"removing stale lock {lock}")
    with open(lock, "w") as fh:
        fh.write(str(os.getpid()))
    return lock


def check_internal_dir(internal_dir, mode):
    """Refuse to append to a checkpoint written by another mode or an older column layout."""
    _disc.assert_internal_path_allowed(internal_dir)
    os.makedirs(internal_dir, exist_ok=True)
    mpath = os.path.join(internal_dir, "65_mode.txt")
    if os.path.exists(mpath):
        old = open(mpath).read().strip()
        if old != mode:
            raise SystemExit(f"FATAL: {internal_dir} holds a --mode {old} run; use another --internal-dir")
    for name, cols in (("per_person_gene.tsv", GENE_COLS), ("per_haplotype.tsv", HAP_COLS)):
        p = os.path.join(internal_dir, name)
        if os.path.exists(p):
            with open(p) as fh:
                hdr = [ln for ln in (fh.readline(), fh.readline()) if not ln.startswith("#")]
            if hdr and hdr[0].rstrip("\n").split("\t") != cols:
                raise SystemExit(f"FATAL: {p} has an older column layout; use a new --internal-dir")
    with open(mpath, "w") as fh:
        fh.write(mode + "\n")


def append_rows(path, rows, cols):
    new = not os.path.exists(path)
    _disc.assert_internal_path_allowed(path)
    with open(path, "a") as fh:
        if new:
            fh.write(_disc.INTERNAL_HEADER_COMMENT + "\n")
            fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def _read(path, cols):
    if not os.path.exists(path):
        return pd.DataFrame(columns=cols)
    return pd.read_csv(path, sep="\t", dtype={"person_id": str}, comment="#", keep_default_na=False)


def done_people(internal_dir, retry_failed=False):
    done = set(_read(os.path.join(internal_dir, "per_haplotype.tsv"), HAP_COLS)["person_id"])
    if not retry_failed:
        done |= set(_read(os.path.join(internal_dir, "failures.tsv"), FAIL_COLS)["person_id"])
    return done


# ---------------------------------------------------------------------------
# Aggregation (PUBLIC)
# ---------------------------------------------------------------------------
def _m(n, mode):
    return _disc.mask(int(n), _disc.PARTICIPANT, mode).display


def _small(mode, *ks):
    return mode == _disc.PUBLIC and any(0 < int(k) < _disc.SUPPRESS_BELOW for k in ks)


def _true(s):
    return s.astype(str).isin(["True", "1", "true"])


def strata(df, platform=True):
    out = [("ancestry", a, df[df["ancestry"] == a]) for a in ANCESTRY_ORDER]
    out.append(("ancestry", "ALL", df))
    if platform:
        out += [("platform", p, df[df["platform"] == p]) for p in sorted(df["platform"].dropna().unique())]
    return [(k, v, g) for k, v, g in out if len(g)]


def genotype_table(gene, people, mode):
    d = gene.merge(people[["person_id", "ancestry", "platform"]], on="person_id", how="left")
    rows = []
    for feat in d["feature"].unique():
        f = d[d["feature"] == feat]
        for kind, val, g in strata(f):
            n = len(g)
            ok = g[g["status"] == "ok"]
            m = len(ok)
            cp = pd.to_numeric(ok["copies"], errors="coerce")
            c0, c1, c2 = (int((cp == k).sum()) for k in (0, 1, 2))
            rate = _disc.mask_rate(m, n, mode, pct_decimals=1)
            row = dict(feature=feat, stratum_type=kind, stratum=val, n_people_measured=_m(n, mode),
                       n_called=_m(m, mode), pct_called=rate["pct"],
                       n_0copies=_m(c0, mode), n_1copy=_m(c1, mode), n_2copies=_m(c2, mode))
            blank = m == 0 or _small(mode, m, c0, c1, c2)
            row["insertion_allele_freq"] = "" if blank else round((c1 + 2 * c2) / (2 * m), 4)
            row["pct_0copies"] = "" if blank else round(100 * c0 / m, 2)
            row["hwe_p"] = "" if blank else float(f"{hwe_exact_p(c0, c1, c2):.3g}")
            rows.append(row)
    return pd.DataFrame(rows)


INCOMPLETE_REASONS = ["no_trb_contig", "anchor_missing", "low_reads", "no_hp_tag"]


def completeness_table(hap, people, mode):
    d = hap.merge(people[["person_id", "ancestry", "platform"]], on="person_id", how="left")
    d["complete"] = _true(d["complete"])
    d["single"] = _true(d["complete_single_contig"])
    rows = []
    for kind, val, g in strata(d):
        n, k = len(g), int(g["complete"].sum())
        rate = _disc.mask_rate(k, n, mode, pct_decimals=1)
        row = dict(stratum_type=kind, stratum=val, n_haplotypes=_m(n, mode), n_complete=_m(k, mode),
                   pct_complete=rate["pct"], n_complete_single_contig=_m(int(g["single"].sum()), mode))
        for r in INCOMPLETE_REASONS:
            row[f"n_{r}"] = _m(int((g["incomplete_reason"] == r).sum()), mode)
        rows.append(row)
    return pd.DataFrame(rows)


def hap_presence(gene, hap, feats):
    """person, hap, one 0/1 column per feature, complete haplotypes only."""
    h = hap.copy()
    h["complete"] = _true(h["complete"])
    keep = h[h["complete"]][["person_id", "hap"]]
    g = gene[gene["feature"].isin(feats)]
    if g.empty:
        return pd.DataFrame()
    rows = []
    for hp in (1, 2):
        x = g.pivot_table(index="person_id", columns="feature", values=f"hap{hp}_present", aggfunc="first")
        x = x.reset_index()
        x["hap"] = hp
        rows.append(x)
    t = pd.concat(rows, ignore_index=True)
    t["hap"] = t["hap"].astype(int)
    keep = keep.assign(hap=keep["hap"].astype(int))
    return t.merge(keep, on=["person_id", "hap"])


def copresence_table(gene, hap, people, feats, mode):
    t = hap_presence(gene, hap, feats)
    if t.empty:
        return pd.DataFrame(), pd.DataFrame()
    t = t.merge(people[["person_id", "ancestry", "platform"]], on="person_id", how="left")
    feats = [f for f in feats if f in t.columns]
    rows, pats = [], []
    for kind, val, g in strata(t, platform=False):
        for i, a in enumerate(feats):
            for b in feats[i + 1:]:
                x, y = g[a].astype(int), g[b].astype(int)
                n11, n10 = int(((x == 1) & (y == 1)).sum()), int(((x == 1) & (y == 0)).sum())
                n01, n00 = int(((x == 0) & (y == 1)).sum()), int(((x == 0) & (y == 0)).sum())
                n = n11 + n10 + n01 + n00
                row = dict(feature_a=a, feature_b=b, stratum_type=kind, stratum=val, n_haplotypes=_m(n, mode),
                           n_both=_m(n11, mode), n_a_only=_m(n10, mode), n_b_only=_m(n01, mode),
                           n_neither=_m(n00, mode))
                blank = n == 0 or _small(mode, n11, n10, n01, n00, n)
                pa, pb = (n11 + n10) / n if n else 0, (n11 + n01) / n if n else 0
                D = n11 / n - pa * pb if n else 0
                den = pa * (1 - pa) * pb * (1 - pb)
                dmax = (min(pa * (1 - pb), (1 - pa) * pb) if D > 0 else min(pa * pb, (1 - pa) * (1 - pb)))
                row["concordance"] = "" if blank else round((n11 + n00) / n, 4)
                row["r2"] = "" if blank or den == 0 else round(D * D / den, 4)
                row["d_prime"] = "" if blank or dmax == 0 else round(D / dmax, 4)
                rows.append(row)
        pat = g[feats].astype(int).astype(str).apply(
            lambda r: ",".join(f"{f}{'+' if v == '1' else '-'}" for f, v in zip(feats, r)), axis=1)
        for p, k in pat.value_counts().sort_index().items():
            pats.append(dict(stratum_type=kind, stratum=val, pattern=p, n_haplotypes=_m(k, mode)))
    return pd.DataFrame(rows), pd.DataFrame(pats)


def af_histogram(gene, people, mode, width=0.05):
    """Reads mode: per feature x stratum, people per allele-fraction bin (status ok only)."""
    d = gene[gene["status"] == "ok"].merge(people[["person_id", "ancestry", "platform"]], on="person_id", how="left")
    d["af"] = pd.to_numeric(d["allele_fraction"], errors="coerce")
    d = d[d["af"].notna()]
    edges = np.round(np.arange(0, 1 + width / 2, width), 4)
    rows = []
    for feat in sorted(d["feature"].unique()):
        f = d[d["feature"] == feat]
        for kind, val, g in strata(f):
            counts, _ = np.histogram(g["af"].clip(0, 1), bins=edges)
            for k, lo in zip(counts, edges[:-1]):
                rows.append(dict(feature=feat, stratum_type=kind, stratum=val, af_lo=float(lo),
                                 af_hi=round(float(lo) + width, 4), n_people=_m(int(k), mode)))
    return pd.DataFrame(rows, columns=["feature", "stratum_type", "stratum", "af_lo", "af_hi", "n_people"])


def compare_methods(reads, asm, people, mode):
    """65b: concordance of per-person copies between two per_person_gene tables (status ok in both)."""
    def ok(d):
        d = d[d["status"] == "ok"][["person_id", "feature", "copies"]].copy()
        d["copies"] = pd.to_numeric(d["copies"], errors="coerce")
        return d[d["copies"].isin([0, 1, 2])]
    m = ok(reads).merge(ok(asm), on=["person_id", "feature"], suffixes=("_reads", "_asm"))
    m = m.merge(people[["person_id", "ancestry", "platform"]], on="person_id", how="left")
    rows = []
    for feat in sorted(m["feature"].unique()):
        f = m[m["feature"] == feat]
        for kind, val, g in strata(f):
            n = len(g)
            k = int((g["copies_reads"] == g["copies_asm"]).sum())
            rate = _disc.mask_rate(k, n, mode, pct_decimals=1)
            row = dict(feature=feat, stratum_type=kind, stratum=val, n_both_called=_m(n, mode),
                       n_concordant=_m(k, mode), pct_concordant=rate["pct"], lo=rate["lo"], hi=rate["hi"])
            for i in (0, 1, 2):
                for j in (0, 1, 2):
                    row[f"n_reads{i}_asm{j}"] = _m(int(((g["copies_reads"] == i) & (g["copies_asm"] == j)).sum()), mode)
            rows.append(row)
    return pd.DataFrame(rows)


def _unrelated_ids(people, relatedness):
    import _s05_common as C
    pairs = C.load_relatedness_pairs(relatedness)
    return C.greedy_unrelated(people["person_id"], pairs)


def _cohort_people(path):
    d = pd.read_csv(path, sep="\t", dtype={"person_id": str})
    d["ancestry"] = d["ancestry_pred"].astype(str).str.upper()
    d["platform"] = d["platform"].astype(str).str.lower()
    return d[["person_id", "ancestry", "platform"]]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def build_probes(args):
    if args.probe_fasta:
        probes = read_probes(args.probe_fasta)
    else:
        probes = fetch_probes_online(args.genes)
    probes = [p for p in probes if p[1].startswith("TRBV") or p[1] in args.genes]
    missing = [g for g in args.genes if g not in {p[1] for p in probes}]
    if missing:
        log(f"WARNING: no probe for {missing} (nomenclature? see the ledger caveats)")
    out = os.path.join(args.out_dir, "probes.fa")
    write_probes(probes, out)
    cls, tab = gene_classes(probes, args.merge_ident)
    _disc.write_tsv(tab, os.path.join(args.out_dir, "probe_classes.tsv"), _disc.PUBLIC)
    log(f"wrote {out}: {len(probes)} alleles, {tab['gene'].nunique()} genes; merged classes: "
        f"{sorted({c for c in cls.values() if '/' in c})}")
    if not (args.alt_fasta or args.ref_fasta):
        return 0
    ad = _args_dict(args, probes, cls)
    rep = sorted({cls[g] for g in args.genes if g in cls}, key=_natkey)
    anc = [cls[g] for g in args.anchors if g in cls]
    rows, res = [], {}
    alt_seqs = dict(iter_fasta(args.alt_fasta)) if args.alt_fasta else {}
    ref = load_ref(args.ref_fasta, args.ref_region) if args.ref_fasta else None
    srcs = []
    if alt_seqs:
        srcs.append(("alt", alt_seqs))
    if ref:
        srcs.append(("primary", {f"{ref[0]}:{ref[1] + 1}-{ref[1] + len(ref[2])}": ref[2]}))
    for label, targets in srcs:
        r = genotype_targets(targets, probes, cls, rep, anc, ad)
        res[label] = r
        for c in rep:
            x = r["calls"][c]
            rows.append(dict(source=label, gene_class=c, present=x["present"], n_loci=x["n_loci"],
                             best_ident=round(x["best_ident"], 4) if np.isfinite(x["best_ident"]) else "",
                             contigs=",".join(x["contigs"]), locus_complete=r["complete"]))
    _disc.write_tsv(pd.DataFrame(rows), os.path.join(args.out_dir, "reference_check.tsv"), _disc.PUBLIC)
    log("reference_check.tsv written; insertion genes must be present on 'alt', absent on 'primary'")
    if alt_seqs and ref:
        ins = [cls[g] for g in args.ins_genes if g in cls and "/" not in cls[g]]
        bp = breakpoint_from_reference(alt_seqs, ref[0], ref[1], ref[2], res["alt"]["assigned"], ins,
                                       args.aligner_resolved, _opts(ad, preset="asm5"))
        if bp is None:
            log("WARNING: breakpoint not found (ALT flanks did not align colinearly); use --breakpoint")
        else:
            meta = [("breakpoint", f"{bp['chrom']}:{bp['left']}-{bp['right']}"),
                    ("ins_len_estimate", bp["ins_len"]), ("alt_contig", bp["alt_contig"]),
                    ("n_left_chunks", bp["n_left"]), ("n_right_chunks", bp["n_right"]),
                    ("block_classes", ",".join(bp["block_classes"])),
                    ("classes_elsewhere", ",".join(bp["classes_elsewhere"])),
                    ("ref_span", f"{ref[0]}:{ref[1] + 1}-{ref[1] + len(ref[2])}")]
            _disc.write_tsv(pd.DataFrame(meta, columns=["key", "value"]),
                            os.path.join(args.out_dir, "probe_meta.tsv"), _disc.PUBLIC)
            log(f"breakpoint {meta[0][1]} (insertion ~{bp['ins_len']} bp, genes {bp['block_classes']}) -> probe_meta.tsv")
            if bp["classes_elsewhere"]:
                log(f"NOTE: {bp['classes_elsewhere']} lie on a separate ALT-only segment; reads mode cannot "
                    "genotype them against this breakpoint and drops them")
    return 0


def _args_dict(args, probes, cls):
    return dict(select=args.select, window=parse_window(args.window), pad=args.pad,
                anchors_kmers=anchor_kmers(probes, set(args.genes)), aligner=args.aligner_resolved,
                preset=args.preset, minimap2=args.minimap2, work_dir=args.work_dir, min_ident=args.min_ident,
                min_cov=args.min_cov, tie_eps=args.tie_eps)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", default="reads", choices=["reads", "assembly"])
    ap.add_argument("--build-probes", action="store_true", help="write probes.fa (+ reference check), then exit")
    ap.add_argument("--probe-fasta", help="build route B: user FASTA instead of IMGT/Ensembl download")
    ap.add_argument("--alt-fasta", help="GRCh38 ALT contig carrying the insertion (reference check, breakpoint)")
    ap.add_argument("--ref-fasta", help="GRCh38 primary sequence (a window, or a chromosome + --ref-region)")
    ap.add_argument("--ref-region", help="e.g. chr7:142038120-143088503 when --ref-fasta is a whole chromosome")
    ap.add_argument("--probes", default=DEFAULT_PROBES, help="probes.fa written by --build-probes")
    ap.add_argument("--probe-meta", help="default: probe_meta.tsv next to --probes")
    ap.add_argument("--breakpoint", help="override, chr7:LEFT-RIGHT (1-based; last left-flank base, first right)")
    ap.add_argument("--genes", default=",".join(DEFAULT_GENES))
    ap.add_argument("--ins-genes", default=",".join(READS_INS_GENES))
    ap.add_argument("--anchors", default=",".join(DEFAULT_ANCHORS))
    ap.add_argument("--copresence", help="default: TRBV4-3,TRBV3-2,TRBV6-3 (assembly); --ins-genes (reads)")
    ap.add_argument("--cohort", default=DEFAULT_COHORT)
    ap.add_argument("--aou-root", default=AOU_ROOT)
    ap.add_argument("--out-dir", help="default ~/s05/results/65_<mode> (65 for --build-probes, 65b for --compare)")
    ap.add_argument("--internal-dir", help="default ~/s05/internal/65_<mode>")
    ap.add_argument("--work-dir", default=DEFAULT_WORK_DIR, help="temp FASTAs (mappy read index, minimap2 binary)")
    ap.add_argument("--people-file", help="TSV with person_id: restrict the run to these people")
    ap.add_argument("--sample-stratified", type=int, metavar="N",
                    help="write N unrelated ancestry-balanced people with assemblies to --people-file, exit")
    ap.add_argument("--seed", type=int, default=65)
    ap.add_argument("--compare", nargs=2, metavar=("READS_TSV", "ASSEMBLY_TSV"),
                    help="65b: concordance of two INTERNAL per_person_gene.tsv tables -> PUBLIC table, exit")
    ap.add_argument("--aligner", default="auto", choices=["auto", "mappy", "minimap2", "kmer"])
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--preset", default="sr")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--bam-column", default="grch38_bam")
    ap.add_argument("--reads-region", default=DEFAULT_READS_REGION)
    ap.add_argument("--flank", type=int, default=1000)
    ap.add_argument("--max-indel", type=int, default=50)
    ap.add_argument("--min-reads", type=int, default=6)
    ap.add_argument("--min-hap-reads", type=int, default=3)
    ap.add_argument("--af-cuts", default="0.2,0.8")
    ap.add_argument("--select", default="auto", choices=["auto", "paf", "anchor", "all"])
    ap.add_argument("--window", default=DEFAULT_WINDOW)
    ap.add_argument("--pad", type=int, default=50_000)
    ap.add_argument("--min-ident", type=float, default=0.97)
    ap.add_argument("--min-cov", type=float, default=0.9)
    ap.add_argument("--merge-ident", type=float, default=0.99)
    ap.add_argument("--tie-eps", type=float, default=0.001)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=1800, help="seconds per person")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--retry-failed", action="store_true")
    ap.add_argument("--unrelated", action="store_true",
                    help="PUBLIC tables on the greedy-unrelated subset (kin >= 0.0442, _s05_common rule)")
    ap.add_argument("--relatedness", default=DEFAULT_RELATEDNESS)
    ap.add_argument("--aggregate-only", action="store_true")
    _disc.add_disclosure_arg(ap)
    args = ap.parse_args(argv)
    args.genes = [g for g in args.genes.split(",") if g]
    args.ins_genes = [g for g in args.ins_genes.split(",") if g]
    args.anchors = [g for g in args.anchors.split(",") if g]
    mode = args.disclosure

    if args.compare:
        out = args.out_dir or os.path.join(RESULTS_ROOT, "65b")
        reads, asm = (_read(p, GENE_COLS) for p in args.compare)
        tab = compare_methods(reads, asm, _cohort_people(args.cohort), mode)
        _disc.write_tsv(tab, os.path.join(out, "65b_method_concordance.tsv"), mode)
        log(f"wrote {out}/65b_method_concordance.tsv ({len(tab)} rows)")
        return 0

    args.aligner_resolved = resolve_aligner(args.aligner, args.minimap2)
    if args.work_dir:
        _disc.assert_internal_path_allowed(args.work_dir)  # holds participant sequence

    if args.build_probes:
        args.out_dir = args.out_dir or DEFAULT_PROBE_DIR
        return build_probes(args)

    args.out_dir = args.out_dir or os.path.join(RESULTS_ROOT, f"65_{args.mode}")
    args.internal_dir = args.internal_dir or os.path.join(INTERNAL_ROOT, f"65_{args.mode}")
    manifest = os.path.join(args.aou_root, LR_MANIFEST_REL)

    if args.sample_stratified:
        if not args.people_file:
            ap.error("--sample-stratified needs --people-file (output path, INTERNAL)")
        _disc.assert_internal_path_allowed(args.people_file)
        people = load_people(args.cohort, manifest, args.aou_root, args.bam_column)
        sel = sample_stratified(people, args.sample_stratified, args.seed, _unrelated_ids(people, args.relatedness))
        _disc.write_tsv(sel, args.people_file, _disc.INTERNAL)
        log(f"wrote {args.people_file}: {len(sel)} people; per ancestry {sel['ancestry'].value_counts().to_dict()}")
        return 0

    if not os.path.exists(args.probes):
        ap.error(f"--probes {args.probes} not found (build it first with --build-probes)")
    probes = read_probes(args.probes)
    cls, cls_tab = gene_classes(probes, args.merge_ident)
    reported = sorted({cls[g] for g in args.genes if g in cls}, key=_natkey)
    anchor_classes = [cls[g] for g in args.anchors if g in cls]
    if len(anchor_classes) != len(args.anchors):
        raise SystemExit(f"FATAL: anchors {args.anchors} not all in the probes")
    ins_classes = []
    for g in args.ins_genes:
        if g not in cls:
            log(f"WARNING: --ins-genes {g} has no probe; skipped")
        elif "/" in cls[g]:
            log(f"WARNING: {g} is merged into {cls[g]} (present on both alleles); not usable from reads, skipped")
        else:
            ins_classes.append(cls[g])
    cop_genes = args.copresence.split(",") if args.copresence else (
        DEFAULT_COPRESENCE if args.mode == "assembly" else args.ins_genes)
    cop = [feature_of_gene(g, cls) for g in cop_genes if g in cls]
    if args.mode == "reads":
        cop = [c for c in cop if c in ins_classes]

    ad = _args_dict(args, probes, cls)
    if args.mode == "reads":
        bps = args.breakpoint or read_probe_meta(
            args.probe_meta or os.path.join(os.path.dirname(os.path.abspath(args.probes)), "probe_meta.tsv")).get("breakpoint")
        if not bps:
            raise SystemExit("FATAL: no breakpoint: run --build-probes with --alt-fasta and --ref-fasta, or give --breakpoint")
        meta = read_probe_meta(args.probe_meta or os.path.join(os.path.dirname(os.path.abspath(args.probes)),
                                                               "probe_meta.tsv"))
        away = [c for c in str(meta.get("classes_elsewhere", "")).split(",") if c and c in ins_classes]
        if away and not args.breakpoint:
            log(f"WARNING: {away} lie outside the insertion block of breakpoint {bps}; dropped from reads mode")
            ins_classes = [c for c in ins_classes if c not in away]
            cop = [c for c in cop if c not in away]
        ad.update(breakpoint=parse_window(bps), flank=args.flank, max_indel=args.max_indel,
                  min_reads=args.min_reads, min_hap_reads=args.min_hap_reads, samtools=args.samtools,
                  af_cuts=tuple(float(x) for x in args.af_cuts.split(",")), reads_region=args.reads_region,
                  ins_classes=ins_classes)
        log(f"reads mode: breakpoint {bps}, insertion classes {ins_classes}, region {args.reads_region}")
    log(f"{len(probes)} probe alleles; reported classes {reported}; aligner {args.aligner_resolved}")

    if not args.aggregate_only and not os.path.exists(manifest):
        raise SystemExit(f"FATAL: {manifest} not found -- is the controlled-data automount up? "
                         f"(ENVIRONMENT.md quirks #26/#35). Refusing to start.")
    if os.path.exists(manifest):
        people = load_people(args.cohort, manifest, args.aou_root, args.bam_column)
    else:
        people = _cohort_people(args.cohort)
    if args.people_file:
        keep = read_people_file(args.people_file)
        people = people[people["person_id"].isin(keep)]
    if args.limit:
        people = people.head(args.limit)
    log(f"{len(people)} long-read people; platforms {people['platform'].value_counts().to_dict()}")

    gene_path = os.path.join(args.internal_dir, "per_person_gene.tsv")
    hap_path = os.path.join(args.internal_dir, "per_haplotype.tsv")
    fail_path = os.path.join(args.internal_dir, "failures.tsv")
    t0 = time.time()
    if not args.aggregate_only:
        check_internal_dir(args.internal_dir, args.mode)
        lock = acquire_lock(args.internal_dir)
        try:
            done = done_people(args.internal_dir, args.retry_failed)
            todo = [r for r in people.to_dict("records") if r["person_id"] not in done]
            log(f"{len(done)} already done (resume), {len(todo)} to genotype, {args.workers} workers")
            n_ok = n_fail = 0
            with ThreadPoolExecutor(max_workers=args.workers) as ex:
                if args.mode == "reads":
                    futs = {ex.submit(genotype_person_reads, p, probes, cls, ad, args.timeout): p for p in todo}
                else:
                    futs = {ex.submit(genotype_person, p, probes, cls, reported, anchor_classes, ad, args.timeout): p
                            for p in todo}
                for i, fut in enumerate(as_completed(futs), 1):
                    p = futs[fut]
                    try:
                        grows, hrows, reason = fut.result()
                    except Exception as e:  # one bad person never stops the run
                        log(f"person {p['person_id']}: {type(e).__name__}: {e}")
                        grows, hrows, reason = None, None, f"error:{type(e).__name__}"
                    if grows:
                        append_rows(gene_path, grows, GENE_COLS)
                        append_rows(hap_path, hrows, HAP_COLS)  # written last = the 'done' marker
                        n_ok += 1
                    else:
                        append_rows(fail_path, [dict(person_id=p["person_id"], reason=reason)], FAIL_COLS)
                        n_fail += 1
                    if i % 100 == 0 or i == len(todo):
                        log(f"{i}/{len(todo)} done ({n_ok} ok, {n_fail} failed), "
                            f"{i / max(1e-9, time.time() - t0):.2f} people/s")
        finally:
            os.remove(lock)

    gene = _read(gene_path, GENE_COLS).drop_duplicates(["person_id", "feature"], keep="last")
    hap = _read(hap_path, HAP_COLS).drop_duplicates(["person_id", "hap"], keep="last")
    fails = _read(fail_path, FAIL_COLS)
    fails = fails[~fails["person_id"].isin(hap["person_id"])].drop_duplicates("person_id", keep="last")
    ids = set(people["person_id"])
    if args.unrelated:
        ids &= _unrelated_ids(people, args.relatedness)
    gene, hap = gene[gene["person_id"].isin(ids)], hap[hap["person_id"].isin(ids)]
    if args.mode == "reads" and len(gene):
        gene = rederive_reads_calls(gene, args.min_reads, ad["af_cuts"])
    os.makedirs(args.out_dir, exist_ok=True)
    _disc.write_tsv(cls_tab, os.path.join(args.out_dir, "probe_classes.tsv"), mode)
    _disc.write_tsv(genotype_table(gene, people, mode), os.path.join(args.out_dir, "genotype_by_stratum.tsv"), mode)
    _disc.write_tsv(completeness_table(hap, people, mode),
                    os.path.join(args.out_dir, "haplotype_completeness.tsv"), mode)
    cp, pats = copresence_table(gene, hap, people, cop, mode)
    _disc.write_tsv(cp, os.path.join(args.out_dir, "copresence.tsv"), mode)
    _disc.write_tsv(pats, os.path.join(args.out_dir, "haplotype_patterns.tsv"), mode)
    if args.mode == "reads":
        _disc.write_tsv(af_histogram(gene, people, mode), os.path.join(args.out_dir, "allele_fraction_histogram.tsv"), mode)
    fr = fails[fails["person_id"].isin(people["person_id"])]["reason"].value_counts()
    _disc.write_tsv(pd.DataFrame([dict(reason=k, n_people=_m(v, mode)) for k, v in sorted(fr.items())],
                                 columns=["reason", "n_people"]),
                    os.path.join(args.out_dir, "failures_summary.tsv"), mode)
    measured = hap["person_id"].nunique()
    called = gene[gene["status"] == "ok"]["person_id"].nunique()
    run = [dict(metric="mode", value=args.mode),
           dict(metric="people_in_cohort", value=_m(len(people), mode)),
           dict(metric="people_genotyped", value=_m(measured, mode)),
           dict(metric="people_called_any_feature", value=_m(called, mode)),
           dict(metric="people_failed", value=_m(int(fr.sum()), mode)),
           dict(metric="unrelated_only", value=args.unrelated),
           dict(metric="aligner", value=args.aligner_resolved), dict(metric="preset", value=args.preset),
           dict(metric="min_ident", value=args.min_ident), dict(metric="min_cov", value=args.min_cov),
           dict(metric="merge_ident", value=args.merge_ident),
           dict(metric="reported_classes", value=",".join(reported)),
           dict(metric="probe_alleles", value=len(probes))]
    if args.mode == "reads":
        anyf = gene[gene["feature"] == INS_ANY]
        tot = pd.to_numeric(anyf["n_ins"], errors="coerce") + pd.to_numeric(anyf["n_ref"], errors="coerce")
        run += [dict(metric="breakpoint", value=f"{ad['breakpoint'][0]}:{ad['breakpoint'][1]}-{ad['breakpoint'][2]}"),
                dict(metric="reads_region", value=args.reads_region), dict(metric="flank", value=args.flank),
                dict(metric="max_indel", value=args.max_indel), dict(metric="min_reads", value=args.min_reads),
                dict(metric="af_cuts", value=args.af_cuts), dict(metric="bam_column", value=args.bam_column),
                dict(metric="insertion_classes", value=",".join(ins_classes)),
                dict(metric="informative_reads_median", value=float(tot.median()) if len(tot) else "")]
    else:
        run += [dict(metric="select", value=args.select), dict(metric="window", value=args.window),
                dict(metric="anchors", value=",".join(args.anchors))]
    run.append(dict(metric="wall_seconds_this_invocation", value=int(time.time() - t0)))
    _disc.write_tsv(pd.DataFrame(run), os.path.join(args.out_dir, "run_summary.tsv"), mode)
    log(f"done: {measured} people genotyped ({called} called); tables in {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
