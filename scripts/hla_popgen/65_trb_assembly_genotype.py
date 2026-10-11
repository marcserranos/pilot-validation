#!/usr/bin/env python3
"""TRB insertion/deletion genotype from each person's phased assemblies (hap1/hap2), by probe mapping.

Sprint S05, WS-B. Why a second TRB script: the TRB locus carries a common ~21 kb insertion/deletion
polymorphism holding TRBV4-3, TRBV3-2 and TRBV6-3 (a near-identical duplicate of TRBV6-2; TRBV7-8 is
listed too because some nomenclatures place it in or next to the polymorphic segment). GRCh38's
primary chr7 carries the haplotype WITHOUT the segment; the inserted sequence exists only on a GRCh38
ALT contig. 58_regional_depth.py reads depth from noalt-GRCh38 BAMs, so reads from the inserted genes
have nowhere to map: 58 cannot genotype TRBV4-3 at all (it is not in the GENCODE primary-assembly
gene set it builds regions from). The phased hifiasm assemblies do not depend on the reference: if a
haplotype carries the insertion, its contig contains the genes. This script reads them directly.
Lead it serves: Aleix's repertoire data, TRBV4-3 never used by ~27% of people (SAS 43.5%, EUR 16.7%).
Prior art (deletion known since the 1990s, ~21 kb, TRBV4-3/3-2/6-2-or-6-3): READING.md section B.

------------------------------------------------------------------------------------------------
Step 1. Probe build (`--build-probes`, run once on the VM, which has open internet)
------------------------------------------------------------------------------------------------
Probes = one FASTA record per allele, header `>ALLELE` or IMGT-style `>acc|GENE*NN|Homo sapiens|...`;
the gene is the text before `*`. Two routes, both write `<out-dir>/probes.fa` + `probe_classes.tsv`:

  Route A (default): IMGT/GENE-DB V-REGION nucleotide FASTA, human, all TRBV alleles (F + ORF + P, so
    pseudogene-classified segments such as TRBV3-2 in some releases are kept), fetched from
        https://www.imgt.org/genedb/GENElect?query=7.2+TRBV&species=Homo+sapiens
    (query code 7.2 = "F+ORF+all P" V-REGION nucleotides; the code is an IMGT convention to CONFIRM
    on the VM by eyeballing the first header lines; if IMGT changes the page, use route B). IMGT gap
    dots are removed. Every TRBV gene is kept as a BACKGROUND probe (used for best-hit assignment
    only, so a locus is never credited to TRBV4-3 when TRBV4-1/4-2 explains it better); only the
    REPORTED genes (below) appear in the outputs. TRBC1 is not a V-REGION and its IMGT C-REGION is
    spliced (it would never align over >=90% of its length to unspliced genomic sequence), so its
    probe is the GENOMIC gene sequence from Ensembl REST:
        https://rest.ensembl.org/lookup/symbol/homo_sapiens/TRBC1  -> gene id
        https://rest.ensembl.org/sequence/id/<id>?type=genomic     (exons + introns, GRCh38)
  Route B: `--probe-fasta FILE` given by the user (e.g. a manual IMGT download), filtered to TRBV/TRBC
    records with a `GENE*allele` or `GENE` header.

  Optional reference check (`--alt-fasta`, `--ref-fasta`, each treated as one pseudo-haplotype):
    the ALT contig that carries the insertion in GRCh38 -- expected to be chr7_KI270803v1_alt
    (GenBank KI270803.1), TO BE CONFIRMED ON THE VM from UCSC's altLocations table:
        curl -s https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/altLocations.txt.gz \\
          | zcat | awk '$2=="chr7" && $3<142900000 && $4>142200000'
    then fetch it with
        curl -sO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/<name>.fa.gz
    (NCBI equivalent: nuccore KI270803.1). `--ref-fasta` = the primary TRB window, e.g.
        samtools faidx ~/ref/<GRCh38>.fa chr7:142250000-142860000 > trb_primary.fa
    Expected: TRBV4-3/3-2 present on the ALT and absent on the primary window; TRBV6-2/6-3 class
    with 2 loci on the ALT and 1 on the primary. Written to `<out-dir>/reference_check.tsv`
    (public reference sequence only; no participant data). If the expectation fails, stop: either
    the probes or the nomenclature are not what this script assumes.

Reported genes (`--genes`): insertion candidates TRBV4-3, TRBV3-2, TRBV6-3, TRBV7-8; controls present
on both haplotypes TRBV4-2, TRBV6-2, TRBV3-1, TRBV2, TRBV12-4, TRBC1.

------------------------------------------------------------------------------------------------
Step 2. Per person (one worker per person; hap1 and hap2 assemblies, hifiasm dual assemblies)
------------------------------------------------------------------------------------------------
  a. Assemblies: `assembly_hap{1,2}_fa` (.fa.gz) and `assembly_hap{1,2}_aln2_hg38_paf` from
     v9/wgs/long_read/manifest.tsv, resolved on the controlled-data automount exactly as 57 does
     (gs://bucket/path -> ~/workspace/vwb-aou-datasets-controlled-v9/path). The mount is READ-ONLY:
     `samtools faidx` there fails trying to write .fai/.gzi (57 copied each assembly locally to avoid
     it). This script never indexes on the mount: it streams the .fa.gz once (mappy.fastx_read, C
     kseq; pure-Python gzip fallback) and keeps only the TRB contigs, stopping as soon as it has them.
  b. TRB contigs (`--select`, default auto = paf, else anchor):
       paf     contigs with primary (tp:A:P) alignments >= 1 kb to the GRCh38 TRB window
               (`--window`, default chr7:142,250,000-142,860,000; the IMGT TRB locus is
               chr7:142,299,011-142,813,287) in the person's own hap->hg38 PAF. The contig is cut to
               [min, max] of the projected window +/- `--pad` (50 kb). min/max, not merged pieces,
               so an inserted 21 kb segment that the PAF leaves unaligned stays inside the cut.
       anchor  no PAF: contigs containing >= 2 of a few exact 25-mers (both strands) taken from the
               reported genes' probes; cut to anchor span +/- 700 kb.
       all     map against every contig (slow: indexes the whole haplotype; use few workers).
  c. Aligner (`--aligner`, default auto = mappy if importable, else the minimap2 binary, else kmer):
       mappy / minimap2  probes (query) vs TRB contig pieces (target), preset `sr` with
                 secondaries kept (`minimap2 -x sr -c --secondary=yes -N 50 -p 0.5`). Why sr and
                 not asm5: V-REGION probes are ~290 bp; asm5 scores -A1 -B19 with a 200 minimum, so a
                 300 bp probe with >=5 mismatches (98.3%) is never reported and lower-identity
                 copies never appear as secondaries; best-hit assignment needs to SEE the paralog hit
                 to reject it. sr (k21/w11 dense minimizers, -A2 -B8, -p 0.5) reports 95% copies as
                 secondaries (checked with mappy 2.31 on synthetic data). `-c` gives base-level
                 alignments: identity = matching bases / alignment block length (PAF col10/col11,
                 mappy mlen/blen); coverage = aligned probe span / probe length.
       kmer      pure-Python fallback (no minimap2 at all): 15-mer seeds on a probe-diagonal,
                 then ungapped identity, refined by a semi-global edit distance when < 99%.
  d. Loci = overlapping hits on a contig (>= 50% of the shorter). Best-hit assignment: each locus
     goes to the gene class whose best probe has the highest identity (probes covering >=
     `--min-cov` ranked first). A locus is CALLED for that class if identity >= `--min-ident`
     (0.97) and coverage >= `--min-cov` (0.9). Ties within `--tie-eps` between two classes ->
     'ambiguous', counted, called for neither. So a 95%-identical paralog copy is never called
     TRBV4-3 (fails 0.97), and a TRBV4-2 locus is never TRBV4-3 even at a lower threshold.
  e. Gene classes: genes whose closest alleles are >= `--merge-ident` (0.99) identical cannot be
     told apart by this method and are merged into one class (TRBV6-2*01 and TRBV6-3*01 share the
     V-REGION sequence in IMGT, so expect class 'TRBV6-2/TRBV6-3'). For a merged class the per-haplotype
     number of called loci is reported (deletion haplotype 1, insertion haplotype 2), plus a derived
     feature '<class>[copy2]' = ">= 2 loci on this haplotype", which is the class's insertion proxy.
  f. Haplotype completeness (QC flag): COMPLETE if one contig carries called loci of every anchor
     (`--anchors`, default TRBV2 = 5' end, TRBC1 = 3' end; the polymorphic segment lies between).
     Otherwise incomplete: no_trb_contig, anchor_missing, or anchors_split (fragmented locus).
  g. copies = hap1_present + hap2_present; status 'ok' if both haplotypes complete, else
     'incomplete' (copies is then a lower bound and excluded from every PUBLIC genotype count).

Outputs
  INTERNAL (`--internal-dir`, default ~/s05/internal/65; refuses any path in the repo or reports/):
    per_person_gene.tsv   person x feature: hap1/2_present, hap1/2_n_loci, copies, hap1/2_complete,
                          status, best_ident, method                     (checkpoint, append+fsync)
    per_haplotype.tsv     person x hap: contigs, bp, complete, reason, loci called/ambiguous, timing
    failures.tsv          person, reason (assembly_missing, timeout, aligner_error, read_error)
  PUBLIC (`--out-dir`, default ~/s05/results/65; aggregate only; participant counts 1-19 -> '<20'
  and any rate/frequency/HWE depending on such a cell blanked, via _disclosure):
    probe_classes.tsv            gene -> class, max cross-gene identity (reference data)
    genotype_by_stratum.tsv      feature x {ancestry, ALL, platform}: n measured, n complete-both,
                                 completeness %, genotype counts 0/1/2 copies (complete-both people),
                                 insertion(presence)-allele frequency, % with 0 copies, HWE exact p
    haplotype_completeness.tsv   stratum: haplotypes, complete, and each incomplete reason
    copresence.tsv               pairs of TRBV4-3 / TRBV3-2 / TRBV6-3 (or the 6-2/6-3[copy2] proxy) on
                                 the SAME complete haplotype: 2x2 counts, concordance, r2, D'
    haplotype_patterns.tsv       joint +/- pattern of those genes per haplotype, counts by stratum
    failures_summary.tsv, run_summary.tsv

Robustness (as 58): mount + manifest hard check; PID lock in the internal dir; per-person rows
appended + fsynced as each person finishes; reruns skip done people (`--retry-failed` retries the
failures); `--timeout` per person (soft deadline between steps, hard on the minimap2 subprocess);
`--limit`; `--workers` threads (mappy and gzip release the GIL).

VM commands (python env ~/repos/pilot-validation/.pixi/envs/spechla/bin/python3, repo
~/repos/pilot-validation-main; full list in sprints/S05_qc_phasing_kir_trb/vm_parts/65.sh):
    PY=~/repos/pilot-validation/.pixi/envs/spechla/bin/python3
    cd ~/repos/pilot-validation-main/scripts/hla_popgen
    $PY -m pip install mappy        # else: --minimap2 ~/repos/pilot-validation/.pixi/envs/spechla/bin/minimap2
    $PY 65_trb_assembly_genotype.py --build-probes --out-dir ~/s05/results/65 \\
        --alt-fasta ~/ref/trb_alt.fa --ref-fasta ~/ref/trb_primary.fa
    $PY 65_trb_assembly_genotype.py --probes ~/s05/results/65/probes.fa --limit 20 --workers 4 \\
        --out-dir ~/s05/results/65_smoke --internal-dir ~/s05/internal/65_smoke     # smoke, own paths
    setsid nohup $PY -u 65_trb_assembly_genotype.py --probes ~/s05/results/65/probes.fa \\
        --workers 12 > ~/s05/logs/65.log 2>&1 < /dev/null & disown

Tests: scripts/hla_popgen/tests/test_65_trb_assembly_genotype.py (synthetic assemblies, kmer + mappy).
Python 3.8 compatible (the spechla env pins python 3.8).
"""
import argparse
import collections
import gzip
import importlib.util
import os
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
DEFAULT_OUT_DIR = os.path.expanduser("~/s05/results/65")
DEFAULT_INTERNAL_DIR = os.path.expanduser("~/s05/internal/65")
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
        import json
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
    if mappy is None:
        raise RuntimeError("mappy not importable")
    hits = []
    for ctg, tseq in targets.items():
        a = mappy.Aligner(seq=tseq, preset=opts.get("preset", "sr"), best_n=50)
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
               "-N", "50", "-p", "0.5", "-t", "1", tfa, qfa]
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
    calls = {}
    for c in classes_reported:
        loc = [a for a in assigned if a["called"] and a["gene_class"] == c]
        allhits = [a for a in assigned if a["gene_class"] == c]
        calls[c] = dict(present=len(loc) > 0, n_loci=len(loc),
                        best_ident=max([a["ident"] for a in allhits], default=float("nan")),
                        contigs=sorted({a["contig"] for a in loc}))
    anchor_ctgs = [set(a["contig"] for a in assigned if a["called"] and a["gene_class"] == c)
                   for c in anchor_classes]
    if any(not s for s in anchor_ctgs):
        complete, reason = False, "anchor_missing"
    elif set.intersection(*anchor_ctgs):
        complete, reason = True, ""
    else:
        complete, reason = False, "anchors_split"
    return calls, complete, reason


# ---------------------------------------------------------------------------
# Contig selection
# ---------------------------------------------------------------------------
def parse_window(w):
    m = re.match(r"^([^:]+):([\d,]+)-([\d,]+)$", w)
    if not m:
        raise ValueError(f"bad --window {w!r}")
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


def genotype_haplotype(fa, paf, probes, cls, reported_classes, anchor_classes, args_d, deadline=None):
    """One haplotype -> dict(calls, complete, reason, n_contigs, bp, n_called, n_ambiguous, mode)."""
    targets, mode = select_contigs(fa, paf, args_d["select"], args_d["window"], args_d["pad"],
                                   args_d["anchors_kmers"], deadline)
    base = dict(n_contigs=len(targets), bp=sum(len(s) for s in targets.values()), select=mode)
    if not targets:
        calls = {c: dict(present=False, n_loci=0, best_ident=float("nan"), contigs=[]) for c in reported_classes}
        return dict(base, calls=calls, complete=False, reason="no_trb_contig", n_called=0, n_ambiguous=0)
    _check(deadline)
    opts = dict(preset=args_d["preset"], minimap2=args_d["minimap2"], work_dir=args_d.get("work_dir"),
                timeout=max(1, (deadline - time.time())) if deadline else None)
    hits = ALIGNERS[args_d["aligner"]]([(a, s) for a, _, s in probes], targets, opts)
    probe_gene = {a: g for a, g, _ in probes}
    assigned = assign_loci(cluster_loci(hits), probe_gene, cls, args_d["min_ident"], args_d["min_cov"],
                           args_d["tie_eps"])
    calls, complete, reason = summarise_haplotype(assigned, reported_classes, anchor_classes)
    return dict(base, calls=calls, complete=complete, reason=reason,
                n_called=sum(a["called"] for a in assigned), n_ambiguous=sum(a["ambiguous"] for a in assigned))


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


def genotype_person(person, probes, cls, reported_classes, anchor_classes, args_d, timeout=None):
    """-> (gene_rows, hap_rows, None) or (None, None, reason)."""
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
    hap_rows = [dict(person_id=pid, hap=h, select=r["select"], n_trb_contigs=r["n_contigs"], trb_bp=r["bp"],
                     complete=r["complete"], incomplete_reason=r["reason"] or "", n_loci_called=r["n_called"],
                     n_loci_ambiguous=r["n_ambiguous"], method=method, seconds=r["seconds"])
                for h, r in haps.items()]
    gene_rows = []
    for feat in features_for(reported_classes):
        p1, p2 = (hap_feature_present(haps[h]["calls"], feat) for h in (1, 2))
        c = feat[:-len(COPY2)] if feat.endswith(COPY2) else feat
        idents = [haps[h]["calls"][c]["best_ident"] for h in (1, 2)]
        idents = [x for x in idents if np.isfinite(x)]
        bi = max(idents) if idents else float("nan")
        gene_rows.append(dict(person_id=pid, feature=feat, hap1_present=int(p1), hap2_present=int(p2),
                              hap1_n_loci=hap_feature_nloci(haps[1]["calls"], feat),
                              hap2_n_loci=hap_feature_nloci(haps[2]["calls"], feat),
                              copies=int(p1) + int(p2), hap1_complete=int(haps[1]["complete"]),
                              hap2_complete=int(haps[2]["complete"]),
                              status="ok" if haps[1]["complete"] and haps[2]["complete"] else "incomplete",
                              best_ident=round(float(bi), 4) if np.isfinite(bi) else "", method=method))
    return gene_rows, hap_rows, None


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def uri_to_local(uri, root=AOU_ROOT):
    """gs://<bucket>/path -> root/path (57's rule); a relative path -> root/path; else None."""
    if not isinstance(uri, str) or not uri or uri == "NA":
        return None
    if uri.startswith("gs://"):
        parts = uri.split("/", 3)
        return os.path.join(root, parts[3]) if len(parts) > 3 else None
    return uri if os.path.isabs(uri) else os.path.join(root, uri)


ASM_COLS = (("hap1_fa", "assembly_hap1_fa"), ("hap2_fa", "assembly_hap2_fa"),
            ("hap1_paf", "assembly_hap1_aln2_hg38_paf"), ("hap2_paf", "assembly_hap2_aln2_hg38_paf"))


def load_people(cohort_path, manifest_path, root=AOU_ROOT):
    """One row per long-read person: person_id, ancestry, platform, hap1/2_fa, hap1/2_paf.
    Manifest row preference: cohort platform first, then any row with both assembly paths.
    No existence check here (mount stats are slow); the worker checks."""
    coh = pd.read_csv(cohort_path, sep="\t", dtype={"person_id": str})
    if "in_lr" in coh.columns:
        coh = coh[coh["in_lr"].astype(str).str.lower().isin(["true", "1"])]
    coh["ancestry"] = coh.get("ancestry_pred", pd.Series(index=coh.index, dtype=str)).astype(str).str.upper()
    man = pd.read_csv(manifest_path, sep="\t", dtype=str).rename(columns={"research_id": "person_id"})
    for c in [c for _, c in ASM_COLS if c not in man.columns]:
        man[c] = None
    by_pid = {pid: g for pid, g in man.groupby("person_id")}
    rows = []
    for _, r in coh.iterrows():
        g = by_pid.get(r["person_id"])
        rec = dict(person_id=r["person_id"], ancestry=r["ancestry"],
                   platform=str(r.get("platform", "")).lower())
        if g is not None:
            cand = pd.concat([g[g["platform"].astype(str).str.lower() == rec["platform"]], g])
            for _, m in cand.iterrows():
                if isinstance(m["assembly_hap1_fa"], str) and isinstance(m["assembly_hap2_fa"], str):
                    for k, c in ASM_COLS:
                        rec[k] = uri_to_local(m[c], root)
                    if rec["platform"] in ("", "nan", "none"):
                        rec["platform"] = str(m["platform"]).lower()
                    break
        for k, _ in ASM_COLS:
            rec.setdefault(k, None)
        rows.append(rec)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Checkpoint + lock
# ---------------------------------------------------------------------------
GENE_COLS = ["person_id", "feature", "hap1_present", "hap2_present", "hap1_n_loci", "hap2_n_loci", "copies",
             "hap1_complete", "hap2_complete", "status", "best_ident", "method"]
HAP_COLS = ["person_id", "hap", "select", "n_trb_contigs", "trb_bp", "complete", "incomplete_reason",
            "n_loci_called", "n_loci_ambiguous", "method", "seconds"]
FAIL_COLS = ["person_id", "reason"]


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


def strata(df):
    out = [("ancestry", a, df[df["ancestry"] == a]) for a in ANCESTRY_ORDER]
    out.append(("ancestry", "ALL", df))
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
            c0, c1, c2 = (int((ok["copies"].astype(int) == k).sum()) for k in (0, 1, 2))
            rate = _disc.mask_rate(m, n, mode, pct_decimals=1)
            row = dict(feature=feat, stratum_type=kind, stratum=val, n_people_measured=_m(n, mode),
                       n_complete_both=_m(m, mode), pct_complete_both=rate["pct"],
                       n_0copies=_m(c0, mode), n_1copy=_m(c1, mode), n_2copies=_m(c2, mode))
            blank = m == 0 or _small(mode, m, c0, c1, c2)
            row["insertion_allele_freq"] = "" if blank else round((c1 + 2 * c2) / (2 * m), 4)
            row["pct_0copies"] = "" if blank else round(100 * c0 / m, 2)
            row["hwe_p"] = "" if blank else float(f"{hwe_exact_p(c0, c1, c2):.3g}")
            rows.append(row)
    return pd.DataFrame(rows)


def completeness_table(hap, people, mode):
    d = hap.merge(people[["person_id", "ancestry", "platform"]], on="person_id", how="left")
    d["complete"] = d["complete"].astype(str).isin(["True", "1", "true"])
    reasons = ["no_trb_contig", "anchor_missing", "anchors_split"]
    rows = []
    for kind, val, g in strata(d):
        n, k = len(g), int(g["complete"].sum())
        rate = _disc.mask_rate(k, n, mode, pct_decimals=1)
        row = dict(stratum_type=kind, stratum=val, n_haplotypes=_m(n, mode), n_complete=_m(k, mode),
                   pct_complete=rate["pct"])
        for r in reasons:
            row[f"n_{r}"] = _m(int((g["incomplete_reason"] == r).sum()), mode)
        rows.append(row)
    return pd.DataFrame(rows)


def hap_presence(gene, hap, feats):
    """Long table person, hap, ancestry-free: one column per feature (0/1), complete haps only."""
    h = hap.copy()
    h["complete"] = h["complete"].astype(str).isin(["True", "1", "true"])
    keep = h[h["complete"]][["person_id", "hap"]]
    rows = []
    g = gene[gene["feature"].isin(feats)]
    for hp in (1, 2):
        x = g.pivot_table(index="person_id", columns="feature", values=f"hap{hp}_present", aggfunc="first")
        x = x.reset_index()
        x["hap"] = hp
        rows.append(x)
    if not rows:
        return pd.DataFrame()
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
    for kind, val, g in strata(t):
        if kind == "platform":
            continue
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


def _unrelated_ids(people, relatedness):
    import _s05_common as C
    pairs = C.load_relatedness_pairs(relatedness)
    return C.greedy_unrelated(people["person_id"], pairs)


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
    refs = [("alt", args.alt_fasta), ("primary", args.ref_fasta)]
    refs = [(k, p) for k, p in refs if p]
    if refs:
        ad = _args_dict(args, probes, cls)
        ad["select"] = "all"
        rep = sorted({cls[g] for g in args.genes if g in cls}, key=_natkey)
        anc = [cls[g] for g in args.anchors if g in cls]
        rows = []
        for label, path in refs:
            r = genotype_haplotype(path, None, probes, cls, rep, anc, ad)
            for c in rep:
                x = r["calls"][c]
                rows.append(dict(source=label, gene_class=c, present=x["present"], n_loci=x["n_loci"],
                                 best_ident=round(x["best_ident"], 4) if np.isfinite(x["best_ident"]) else "",
                                 contigs=",".join(x["contigs"]), locus_complete=r["complete"]))
        _disc.write_tsv(pd.DataFrame(rows), os.path.join(args.out_dir, "reference_check.tsv"), _disc.PUBLIC)
        log("reference_check.tsv written; insertion genes must be present on 'alt', absent on 'primary'")
    return 0


def _args_dict(args, probes, cls):
    return dict(select=args.select, window=parse_window(args.window), pad=args.pad,
                anchors_kmers=anchor_kmers(probes, set(args.genes)), aligner=args.aligner_resolved,
                preset=args.preset, minimap2=args.minimap2, work_dir=args.work_dir, min_ident=args.min_ident,
                min_cov=args.min_cov, tie_eps=args.tie_eps)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-probes", action="store_true", help="write probes.fa (+ reference check), then exit")
    ap.add_argument("--probe-fasta", help="build route B: user FASTA instead of IMGT/Ensembl download")
    ap.add_argument("--alt-fasta", help="GRCh38 ALT contig carrying the insertion (reference check)")
    ap.add_argument("--ref-fasta", help="GRCh38 primary TRB window (reference check)")
    ap.add_argument("--probes", help="probes.fa written by --build-probes (genotyping run)")
    ap.add_argument("--genes", default=",".join(DEFAULT_GENES))
    ap.add_argument("--anchors", default=",".join(DEFAULT_ANCHORS))
    ap.add_argument("--copresence", default=",".join(DEFAULT_COPRESENCE))
    ap.add_argument("--cohort", default=DEFAULT_COHORT)
    ap.add_argument("--aou-root", default=AOU_ROOT)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    ap.add_argument("--internal-dir", default=DEFAULT_INTERNAL_DIR)
    ap.add_argument("--work-dir", default=DEFAULT_WORK_DIR, help="temp FASTAs for the minimap2 binary only")
    ap.add_argument("--aligner", default="auto", choices=["auto", "mappy", "minimap2", "kmer"])
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--preset", default="sr")
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
    args.anchors = [g for g in args.anchors.split(",") if g]
    args.aligner_resolved = resolve_aligner(args.aligner, args.minimap2)
    if args.work_dir:
        _disc.assert_internal_path_allowed(args.work_dir)  # holds participant contig sequence

    if args.build_probes:
        return build_probes(args)

    if not args.probes:
        ap.error("--probes is required (build it first with --build-probes)")
    probes = read_probes(args.probes)
    cls, cls_tab = gene_classes(probes, args.merge_ident)
    reported = sorted({cls[g] for g in args.genes if g in cls}, key=_natkey)
    anchor_classes = [cls[g] for g in args.anchors if g in cls]
    if len(anchor_classes) != len(args.anchors):
        raise SystemExit(f"FATAL: anchors {args.anchors} not all in the probes")
    cop = [feature_of_gene(g, cls) for g in args.copresence.split(",") if g in cls]
    log(f"{len(probes)} probe alleles; reported classes {reported}; aligner {args.aligner_resolved}")

    manifest = os.path.join(args.aou_root, LR_MANIFEST_REL)
    if not args.aggregate_only and not os.path.exists(manifest):
        raise SystemExit(f"FATAL: {manifest} not found -- is the controlled-data automount up? "
                         f"(ENVIRONMENT.md quirks #26/#35). Refusing to start.")
    people = load_people(args.cohort, manifest, args.aou_root) if os.path.exists(manifest) else \
        pd.read_csv(args.cohort, sep="\t", dtype={"person_id": str}).assign(
            ancestry=lambda d: d["ancestry_pred"].astype(str).str.upper(),
            platform=lambda d: d["platform"].astype(str).str.lower())
    if args.limit:
        people = people.head(args.limit)
    log(f"{len(people)} long-read people; platforms {people['platform'].value_counts().to_dict()}")

    gene_path = os.path.join(args.internal_dir, "per_person_gene.tsv")
    hap_path = os.path.join(args.internal_dir, "per_haplotype.tsv")
    fail_path = os.path.join(args.internal_dir, "failures.tsv")
    ad = _args_dict(args, probes, cls)
    t0 = time.time()
    if not args.aggregate_only:
        lock = acquire_lock(args.internal_dir)
        try:
            done = done_people(args.internal_dir, args.retry_failed)
            todo = [r for r in people.to_dict("records") if r["person_id"] not in done]
            log(f"{len(done)} already done (resume), {len(todo)} to genotype, {args.workers} workers")
            n_ok = n_fail = 0
            with ThreadPoolExecutor(max_workers=args.workers) as ex:
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
    mode = args.disclosure
    os.makedirs(args.out_dir, exist_ok=True)
    _disc.write_tsv(cls_tab, os.path.join(args.out_dir, "probe_classes.tsv"), mode)
    _disc.write_tsv(genotype_table(gene, people, mode), os.path.join(args.out_dir, "genotype_by_stratum.tsv"), mode)
    _disc.write_tsv(completeness_table(hap, people, mode),
                    os.path.join(args.out_dir, "haplotype_completeness.tsv"), mode)
    cp, pats = copresence_table(gene, hap, people, cop, mode)
    _disc.write_tsv(cp, os.path.join(args.out_dir, "copresence.tsv"), mode)
    _disc.write_tsv(pats, os.path.join(args.out_dir, "haplotype_patterns.tsv"), mode)
    fr = fails[fails["person_id"].isin(people["person_id"])]["reason"].value_counts()
    _disc.write_tsv(pd.DataFrame([dict(reason=k, n_people=_m(v, mode)) for k, v in sorted(fr.items())],
                                 columns=["reason", "n_people"]),
                    os.path.join(args.out_dir, "failures_summary.tsv"), mode)
    measured = hap["person_id"].nunique()
    both = gene[gene["status"] == "ok"]["person_id"].nunique()
    run = [dict(metric="people_in_cohort", value=_m(len(people), mode)),
           dict(metric="people_genotyped", value=_m(measured, mode)),
           dict(metric="people_complete_both_haplotypes", value=_m(both, mode)),
           dict(metric="people_failed", value=_m(int(fr.sum()), mode)),
           dict(metric="unrelated_only", value=args.unrelated),
           dict(metric="aligner", value=args.aligner_resolved), dict(metric="preset", value=args.preset),
           dict(metric="select", value=args.select), dict(metric="window", value=args.window),
           dict(metric="min_ident", value=args.min_ident), dict(metric="min_cov", value=args.min_cov),
           dict(metric="merge_ident", value=args.merge_ident), dict(metric="anchors", value=",".join(args.anchors)),
           dict(metric="reported_classes", value=",".join(reported)),
           dict(metric="probe_alleles", value=len(probes)),
           dict(metric="wall_seconds_this_invocation", value=int(time.time() - t0))]
    _disc.write_tsv(pd.DataFrame(run), os.path.join(args.out_dir, "run_summary.tsv"), mode)
    log(f"done: {measured} people genotyped ({both} complete on both haplotypes); tables in {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
