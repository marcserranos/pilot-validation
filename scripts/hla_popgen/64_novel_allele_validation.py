#!/usr/bin/env python3
"""Are novel HLA calls real alleles or assembly / basecalling artifacts?

Sprint S05 (`SPRINT.md` row L, numbered 64). 50d showed that the classical-gene protein-novel
funnel collapses (8 genes: 4,421 distinct protein-novel sequences -> 174 with >= 2 unrelated
carriers -> 8 not homopolymer-indel-only -> 0 clean at >= 20 carriers), that Sequel IIe calls more
novel alleles than Revio within ancestry, and that Immuannot typed against IPD-IMGT/HLA 3.55.0.
DECISIONS.md (S05, open) proposes demoting classical protein-novel counts until (1) a current
IPD-IMGT/HLA release and (2) read-level support have been checked. This script does both.

------------------------------------------------------------------------------------------------
Which calls
------------------------------------------------------------------------------------------------
Table 1 rows (`hla_calls_rich.tsv`) with `novelty_depth` in `--depths` (default 2,3: protein and
synonymous novelty; 4 = beyond-CDS is optional, see below), `copy_index == 1`, consensus not
`undetermined`, carried by the greedy-unrelated set (kin >= 0.0442, `_s05_common` rule).
Sequence key (the same key as 50d): `(gene, cds_mut)`. For depth 4 `cds_mut` is empty (the CDS
equals a known allele), so the key is `(gene, sha1 of the genomic gene span)` and is only known for
tested calls. Per key at most `--max-per-key` carrier haplotypes are tested, picked round-robin over
platform x ancestry strata in a fixed hash order (reruns pick the same calls).

------------------------------------------------------------------------------------------------
(a) Re-match against a current IPD-IMGT/HLA release
------------------------------------------------------------------------------------------------
`--fetch-ipd` (VM only, open internet) downloads from the ANHIG/IMGTHLA GitHub repository, at git
ref `--ipd-release` (default `Latest`, the repository's default branch = current release):
    https://raw.githubusercontent.com/ANHIG/IMGTHLA/<ref>/hla_nuc.fasta        (CDS per allele)
    https://raw.githubusercontent.com/ANHIG/IMGTHLA/<ref>/hla_gen.fasta        (genomic per allele)
    https://raw.githubusercontent.com/ANHIG/IMGTHLA/<ref>/Allelelist.txt       (version line)
    https://raw.githubusercontent.com/ANHIG/IMGTHLA/<ref>/Allelelist_history.txt (allele per release)
If a file comes back as a Git-LFS pointer, it is re-fetched from
    https://media.githubusercontent.com/media/ANHIG/IMGTHLA/<ref>/<file>
into `--ipd-dir` (default ~/s05/ref/ipd_<ref>). `--ipd-dir` alone (no fetch) uses files already
there, e.g. a manual download.
"Added since 3.55.0" = the allele's accession (HLA00001 ...) has no name in the `3550` column of
`Allelelist_history.txt` (`--baseline-col`). Fallback when that file is absent: accessions found in
Immuannot's own refdata (`--refdata`/CDSseq/*.fa.gz, the 3.55.0 snapshot Immuannot typed against).

Observed sequence per tested call: the called CDS is rebuilt from the haplotype's own Immuannot
files: GTF `CDS` rows of `(contig, gene, copy)` in `<outroot>/<pid>/immuannot_output/<hap>.gtf.gz`
sliced from `<hap>.trimmed.fa` (reference/IMMUANNOT_GTF_SPEC.md "indirect route"), cross-checked
against `<hap>/cds.fa.gz` (03's direct route; `cds_check` column). Stop codons are stripped on
both sides before comparing (24's `cds_core`).
Classes per key (first that applies):
    known_in_baseline          exact CDS (depth 4: genomic, 44's containment rule) of an allele
                               already in 3.55.0 -- the "novel" label itself was wrong
    protein_known_in_baseline  depth 2 only: the protein equals a 3.55.0 allele's protein
    now_known                  exact CDS (depth 4: genomic) of an allele added since 3.55.0
    now_known_protein          depth 2 only: protein equals that of an allele added since 3.55.0
    now_known_partial          a partial-sequence allele added since 3.55.0 is contained in the
                               observed CDS and NOT in the nearest known CDS (so it covers at
                               least one of the novel differences)
    catalogued_age_unknown     exact match, but no baseline list was available
    still_novel                none of the above
    no_sequence                no tested call produced a sequence

------------------------------------------------------------------------------------------------
(b) Read-level support from the person's own HiFi reads
------------------------------------------------------------------------------------------------
Per tested call:
 1. Called haplotype window Wc: the call's contig from `<hap>.trimmed.fa`, gene span +-
    `--window-pad` (2 kb), in coding orientation (reverse-complemented on the `-` strand).
 2. Nearest known allele: the reference allele in `cds_mut` (Immuannot's own nearest), its CDS from
    the current release (fallback: Immuannot refdata). Called CDS vs nearest CDS are aligned
    (banded global alignment, unit costs) and every difference is placed on Wc through the GTF exon
    coordinates. Tn = Wc with those differences reverted = "the same haplotype carrying the nearest
    known allele". Depth 4: the nearest known GENOMIC allele (`template_allele`) is aligned to Wc
    with minimap2 asm5 (mappy) and its differences inside the aligned span are used instead.
 3. Wo = the other haplotype's window for the same gene (Table 1 copy 1 of the other hap), if any.
 4. Reads: `samtools view -F 0x904 <grch38_bam> <region>` (primary alignments only), region = the
    gene in 58's `regions.bed` (GENCODE GRCh38) +- `--region-pad`; DRB paralogs missing from the
    primary assembly use the span of the DR genes present. At most `--max-reads` reads (fixed hash
    order). BAM path: `grch38_bam` in v9/wgs/long_read/manifest.tsv via the controlled-data automount.
 5. Each read is aligned to Wc (mappy, `map-hifi`). Differences are grouped into clusters (events
    closer than 2 x `--flank` merged; homopolymer runs included whole). For every cluster window the
    read fully spans, edit distance read-segment vs Wc-segment and vs Tn-segment; the read votes
    `called` if the summed distance is lower to Wc, `nearest` if lower to Tn, else no vote.
 6. Haplotype of the read (when Wo exists): the read is also aligned to Wo; over the query span both
    alignments cover, the edits OUTSIDE the variant windows are counted against Wc and against Wo
    (the tested sites never decide the haplotype). Difference >= `--hap-margin` -> `this` or
    `other`; else `unknown`. `other` reads are dropped.
 7. Counting: three-target mode (`3t`) when >= `--min-support` informative reads are assigned to
    this haplotype: expected fraction ~1 for a real allele, ~0 for an artifact. Otherwise two-target
    mode (`2t`) on this + unknown reads: ~0.5 (heterozygous) or ~1 (homozygous) for a real allele,
    ~0 for an artifact.
    n_reads_support_called, n_reads_support_nearest, frac_called = called / (called + nearest).
    supported     called >= --min-support and frac >= --min-frac (3t, 0.7) / --min-frac-2t (0.3)
    contradicted  nearest >= --min-support and frac <= --max-frac-contradicted (0.15)
    ambiguous     otherwise (lr_reason: low_reads / mixed / no_window_spanned ...)
 8. Homopolymer indels (indel of one repeated base next to a run, longer side >= 3, the 25 rule on
    sequence context): per spanning read, the run length is read out between 6-bp exact flanks and
    counted as called length / nearest length / other length (hp_n_* columns). A call is
    `hp_only_seq` when every difference is such an indel. The 50d proxy on `cds_mut`
    (`is_homopolymer_indel_only`) is kept as `hp_only_proxy` and is the split used in the PUBLIC
    tables (comparable with 50d).

(c) Optional short-read check (`--sr`): the same procedure with srWGS reads from the person's CRAM
(`cram_uri` in v9/wgs/cram/manifest.csv, decoded with `-T --sr-ref`), mappy preset `sr`, a read
votes only on cluster windows it spans; `--sr-min-reads` reads must vote. Columns sr_*.

------------------------------------------------------------------------------------------------
Outputs
------------------------------------------------------------------------------------------------
PUBLIC (`--out-dir`, default ~/s05/results/64; aggregate only; via _disclosure, PUBLIC default:
call / carrier / read counts 1-19 -> "<20" and a rate over such a cell blanked; distinct-key counts
exact, as allele properties):
    ipd_release.tsv              release ref + version, records, alleles added since baseline/gene
    summary_by_stratum.tsv       gene (+ classical_8, ALL) x novelty class x platform x ancestry
                                 (with ALL margins) x homopolymer split (hp_only / other / all):
                                 n keys, n calls, n tested, % keys now known (each class),
                                 % supported / contradicted / ambiguous among tested calls
    funnel_with_evidence.tsv     the 50d funnel steps per gene and novelty class, each step with
                                 n keys and how many are known_in_baseline / now_known* / read-
                                 supported / contradicted / mixed / ambiguous / untested
    now_known_alleles.tsv        public IPD allele names matched by a novel key, carriers masked,
                                 review_flag where carriers < 20 (publication is a human call)
    hp_runlength_by_gene.tsv     homopolymer run-length readout totals (reads) per gene and class
    run_summary.tsv, failures_summary.tsv
INTERNAL (`--internal-dir`, default ~/s05/internal/64; refuses the repo and any reports/ path):
    per_call_evidence.tsv        one row per tested call (person_id), append + fsync checkpoint
    key_evidence.tsv             one row per key: carriers, IPD class, key-level read verdict

Robustness (as 58/65): mount + manifest hard check before reads are touched; PID lock in the
internal dir; per-person rows appended + fsynced as each person finishes; reruns skip done people
(`--retry-failed` redoes people whose calls all errored); `--timeout` on every samtools call; one
failing person never stops the run; `--limit` (people); `--workers` threads (mappy and samtools
release the GIL; the Python edit-distance path only runs on reads that match neither window).
Runners are module-level and patchable in tests: `fetch_reads`, `fetch_url`.

VM commands (python ~/repos/pilot-validation/.pixi/envs/spechla/bin/python3; samtools in the same
env; repo ~/repos/pilot-validation-main; full list in sprints/S05_qc_phasing_kir_trb/vm_parts/64.sh):
    PY=~/repos/pilot-validation/.pixi/envs/spechla/bin/python3
    cd ~/repos/pilot-validation-main/scripts/hla_popgen
    $PY 64_novel_allele_validation.py --fetch-ipd --ipd-release Latest --ipd-dir ~/s05/ref/ipd_latest
    $PY -u 64_novel_allele_validation.py --ipd-dir ~/s05/ref/ipd_latest --genes A,B,C,DRB1,DQA1,DQB1,DPA1,DPB1 \\
        --limit 30 --workers 4 --samtools ~/repos/pilot-validation/.pixi/envs/spechla/bin/samtools \\
        --out-dir ~/s05/results/64_smoke --internal-dir ~/s05/internal/64_smoke        # smoke, own paths
    setsid nohup $PY -u 64_novel_allele_validation.py --ipd-dir ~/s05/ref/ipd_latest --workers 12 \\
        --samtools ~/repos/pilot-validation/.pixi/envs/spechla/bin/samtools \\
        > ~/s05/logs/64.log 2>&1 < /dev/null & disown

Tests: scripts/hla_popgen/tests/test_64_novel_allele_validation.py (synthetic contigs, GTFs and
read sets; fake IPD release; no samtools, no network). Python 3.8 compatible.
"""
import argparse
import bisect
import gzip
import hashlib
import importlib.util
import os
import re
import subprocess
import sys
import time
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import _disclosure as _disc  # noqa: E402
import _s05_common as C  # noqa: E402

try:
    import mappy  # noqa: E402
except ImportError:  # pragma: no cover - depends on the env
    mappy = None


def _load_50():
    """50_qc_foundation.py (is_homopolymer_indel_only, the 50d proxy) -- loaded by path."""
    spec = importlib.util.spec_from_file_location("_qc50", os.path.join(_THIS_DIR, "50_qc_foundation.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


is_homopolymer_indel_only = _load_50().is_homopolymer_indel_only

AOU_ROOT = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9")
LR_MANIFEST_REL = "v9/wgs/long_read/manifest.tsv"
SR_MANIFEST_REL = "v9/wgs/cram/manifest.csv"
DEFAULT_SR_REF = os.path.expanduser("~/ref/Homo_sapiens_assembly38.fasta")
DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs/people")
DEFAULT_REGIONS = os.path.expanduser("~/s05/results/58/regions.bed")
DEFAULT_REFDATA = os.path.expanduser("~/tools/Immuannot_refdata")
DEFAULT_OUT_DIR = os.path.expanduser("~/s05/results/64")
DEFAULT_INTERNAL_DIR = os.path.expanduser("~/s05/internal/64")

IPD_REPO = "ANHIG/IMGTHLA"
IPD_RAW = "https://raw.githubusercontent.com/" + IPD_REPO + "/{ref}/{file}"
IPD_MEDIA = "https://media.githubusercontent.com/media/" + IPD_REPO + "/{ref}/{file}"
IPD_FILES = ["hla_nuc.fasta", "hla_gen.fasta", "Allelelist.txt", "Allelelist_history.txt"]
BASELINE_COL = "3550"

ANCESTRY_ORDER = C.ANCESTRY_ORDER
CLASSICAL8 = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
DEPTH_CLASS = {2: "protein", 3: "synonymous", 4: "noncoding"}
IPD_CLASSES = ["known_in_baseline", "protein_known_in_baseline", "now_known", "now_known_protein",
               "now_known_partial", "catalogued_age_unknown", "still_novel", "no_sequence"]
NOW_KNOWN = ("now_known", "now_known_protein", "now_known_partial")
VERDICTS = ["supported", "contradicted", "ambiguous"]
HP_FLANK = 6


def log(msg):
    print(f"[64 {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------
_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AA[16 * i + 4 * j + k]
         for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}
STOPS = {"TAA", "TAG", "TGA"}


def revcomp(s):
    return s.translate(_RC)[::-1]


def sha1(s):
    return hashlib.sha1(s.encode("ascii")).hexdigest()


def cds_core(seq):
    """Strip one terminal stop codon. Unlike 24's cds_core this does not require the length to be
    a multiple of 3: Immuannot's GTF CDS rows exclude the stop while IPD's hla_nuc.fasta includes
    it, and a frameshifted called CDS must still be compared stop-free on both sides."""
    seq = seq.upper()
    if len(seq) >= 3 and seq[-3:] in STOPS:
        return seq[:-3]
    return seq


def protein_of(core):
    """Translation of a stop-stripped CDS, or None when it is out of frame / has a stop."""
    if not core or len(core) % 3:
        return None
    p = "".join(CODON.get(core[i:i + 3], "X") for i in range(0, len(core), 3))
    return None if "*" in p else p


def _open_text(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


def iter_fasta(path):
    """(header_line_without_>, seq) from a plain or gzipped FASTA."""
    head, buf = None, []
    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if head is not None:
                    yield head, "".join(buf).upper()
                head, buf = line[1:], []
            else:
                buf.append(line.strip())
    if head is not None:
        yield head, "".join(buf).upper()


def norm_allele(name):
    """'HLA-A*01:01:01:01' / 'A*01:01:01:01' -> 'A*01:01:01:01'."""
    name = str(name).strip()
    return name[4:] if name.startswith("HLA-") else name


_ACC_RE = re.compile(r"^(?:HLA:)?(HLA\d{4,})$")


def parse_allele_header(head):
    """IPD '>HLA:HLA00001 A*01:01:01:01 1098 bp' or Immuannot refdata
    '>HLA-A*01:01:01:01 HLA00001 frame=1 1098bp' -> (name, gene, accession or None)."""
    toks = head.split()
    name = next((t for t in toks if "*" in t), None)
    if name is None:
        return None, None, None
    name = norm_allele(name)
    acc = None
    for t in toks:
        m = _ACC_RE.match(t)
        if m:
            acc = m.group(1)
            break
    return name, C.norm_gene(name.split("*", 1)[0]), acc


# ---------------------------------------------------------------------------
# IPD-IMGT/HLA release
# ---------------------------------------------------------------------------
def fetch_url(url, timeout=300):
    """Patchable downloader -> bytes."""
    req = urllib.request.Request(url, headers={"User-Agent": "pilot-validation/64"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def fetch_ipd(ref, dest, files=IPD_FILES):
    """Downloads the IPD-IMGT/HLA files at git ref `ref` into `dest`. LFS pointer -> media URL."""
    os.makedirs(dest, exist_ok=True)
    got = {}
    for f in files:
        url = IPD_RAW.format(ref=ref, file=f)
        try:
            data = fetch_url(url)
        except Exception as e:  # noqa: BLE001 - report and continue; the history file is optional
            log(f"fetch failed {url}: {type(e).__name__}: {e}")
            got[f] = "failed"
            continue
        if data[:40].startswith(b"version https://git-lfs"):
            url = IPD_MEDIA.format(ref=ref, file=f)
            log(f"{f} is a Git-LFS pointer; fetching {url}")
            data = fetch_url(url)
        tmp = os.path.join(dest, f + ".part")
        with open(tmp, "wb") as fh:
            fh.write(data)
        os.replace(tmp, os.path.join(dest, f))
        got[f] = url
        log(f"wrote {os.path.join(dest, f)} ({len(data)} bytes)")
    return got


def read_version(ipd_dir):
    p = os.path.join(ipd_dir, "Allelelist.txt")
    if not os.path.exists(p):
        return ""
    with open(p) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if "version" in line.lower():
                return line.split(":", 1)[-1].strip()
    return ""


def load_baseline_history(path, col=BASELINE_COL):
    """Accessions that had a name in release column `col` of Allelelist_history.txt."""
    if not path or not os.path.exists(path):
        return None
    lines = [ln.rstrip("\n\r") for ln in open(path) if ln.strip() and not ln.startswith("#")]
    if not lines:
        return None
    delim = "\t" if "\t" in lines[0] else ","
    header = [h.strip() for h in lines[0].split(delim)]
    if col not in header:
        log(f"WARNING: column {col} not in {path}; baseline unavailable from history")
        return None
    j, accs = header.index(col), set()
    for ln in lines[1:]:
        f = ln.split(delim)
        if len(f) > j and f[j].strip() and f[j].strip().upper() not in ("NA", "NAN", "NONE", "-"):
            accs.add(f[0].strip())
    return accs


def load_baseline_refdata(refdata):
    d = os.path.join(refdata or "", "CDSseq")
    if not refdata or not os.path.isdir(d):
        return None
    accs = set()
    for fn in sorted(os.listdir(d)):
        if fn.endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz")):
            for head, _ in iter_fasta(os.path.join(d, fn)):
                _, _, acc = parse_allele_header(head)
                if acc:
                    accs.add(acc)
    return accs or None


class IpdRelease:
    """Current-release allele sequences + the 3.55.0 baseline accession set."""

    def __init__(self, ipd_dir, baseline_accs, baseline_source, refdata=None, load_gen=False):
        nuc = os.path.join(ipd_dir, "hla_nuc.fasta")
        if not os.path.exists(nuc):
            raise SystemExit(f"FATAL: {nuc} not found -- run --fetch-ipd first or point --ipd-dir at "
                             f"a directory holding hla_nuc.fasta (and hla_gen.fasta for depth 4).")
        self.version = read_version(ipd_dir)
        self.baseline = baseline_accs
        self.baseline_source = baseline_source
        self.refdata = refdata
        self.recs = defaultdict(list)          # gene -> [(name, acc, core)]
        self.exact = defaultdict(lambda: defaultdict(list))
        self.by_name = {}
        self.acc_of = {}
        for head, seq in iter_fasta(nuc):
            name, gene, acc = parse_allele_header(head)
            if name is None:
                continue
            core = cds_core(seq)
            self.recs[gene].append((name, acc, core))
            self.exact[gene][core].append(name)
            self.by_name[name] = core
            self.acc_of[name] = acc
        self.n_nuc = sum(len(v) for v in self.recs.values())
        self._prot = None
        self.gen = defaultdict(list)
        self.gen_by_name = {}
        self.n_gen = 0
        gen = os.path.join(ipd_dir, "hla_gen.fasta")
        if load_gen and os.path.exists(gen):
            for head, seq in iter_fasta(gen):
                name, gene, acc = parse_allele_header(head)
                if name is None:
                    continue
                self.gen[gene].append((name, acc, seq))
                self.gen_by_name[name] = seq
                self.acc_of.setdefault(name, acc)
            self.n_gen = sum(len(v) for v in self.gen.values())
        self._ref_cache = {}
        self._cache = {}

    def in_baseline(self, name):
        if self.baseline is None:
            return None
        acc = self.acc_of.get(name)
        return None if acc is None else acc in self.baseline

    def prot_index(self):
        if self._prot is None:
            idx = defaultdict(lambda: defaultdict(list))
            for gene, recs in self.recs.items():
                for name, _acc, core in recs:
                    p = protein_of(core)
                    if p and p.startswith("M"):
                        idx[gene][p].append(name)
            self._prot = idx
        return self._prot

    def new_counts(self):
        rows = []
        for gene in sorted(self.recs):
            names = [n for n, _a, _c in self.recs[gene]]
            new = [n for n in names if self.in_baseline(n) is False]
            rows.append(dict(gene=gene, n_alleles_nuc=len(names),
                             n_added_since_baseline=len(new) if self.baseline is not None else ""))
        return rows

    # -- nearest known sequences ------------------------------------------------------------
    def _refdata_seqs(self, gene, kind):
        key = (gene, kind)
        if key not in self._ref_cache:
            seqs = {}
            if self.refdata:
                if kind == "cds":
                    cands = [os.path.join(self.refdata, "CDSseq", f"{p}{gene}.fa{z}")
                             for p in ("HLA-", "") for z in (".gz", "")]
                else:
                    cands = [os.path.join(self.refdata, "gen.fa.gz"), os.path.join(self.refdata, "gen.fa")]
                for path in cands:
                    if os.path.exists(path):
                        for head, seq in iter_fasta(path):
                            name, g, _acc = parse_allele_header(head)
                            if name and g == gene:
                                seqs[name] = cds_core(seq) if kind == "cds" else seq
                        break
            self._ref_cache[key] = seqs
        return self._ref_cache[key]

    def nearest_cds(self, gene, name, core=None):
        """-> (name, core, source) for the nearest known allele; name first, then same-length
        Hamming over the current release; (None, None, reason) when nothing usable."""
        if name:
            name = norm_allele(name)
            if name in self.by_name:
                return name, self.by_name[name], "ipd_current"
            ref = self._refdata_seqs(gene, "cds").get(name)
            if ref:
                return name, ref, "immuannot_refdata"
        if core:
            same = [(n, c) for n, _a, c in self.recs.get(gene, []) if len(c) == len(core) and c != core]
            if same:
                a = np.frombuffer(core.encode(), dtype=np.uint8)
                best = min(same, key=lambda nc: int((np.frombuffer(nc[1].encode(), dtype=np.uint8) != a).sum()))
                return best[0], best[1], "hamming_fallback"
        return None, None, "no_nearest"

    def nearest_gen(self, gene, name):
        name = norm_allele(name) if name else None
        if name and name in self.gen_by_name:
            return name, self.gen_by_name[name], "ipd_current"
        if name:
            ref = self._refdata_seqs(gene, "gen").get(name)
            if ref:
                return name, ref, "immuannot_refdata"
        return None, None, "no_nearest_genomic"

    # -- classification --------------------------------------------------------------------
    def classify(self, gene, depth, seq, nearest_core=None):
        """-> (ipd_class, matched names). seq = called CDS core (depth 2/3) or genomic span (4).
        Cached by (gene, depth, sha1(seq)) -- the same sequence gets the same answer."""
        if not seq:
            return "no_sequence", []
        ck = (gene, depth, sha1(seq), sha1(nearest_core or ""))
        if ck in self._cache:
            return self._cache[ck]
        res = self._classify(gene, depth, seq, nearest_core)
        self._cache[ck] = res
        return res

    def _bucket(self, names):
        base = [n for n in names if self.in_baseline(n) is True]
        new = [n for n in names if self.in_baseline(n) is False]
        unk = [n for n in names if self.in_baseline(n) is None]
        return base, new, unk

    def _classify(self, gene, depth, seq, nearest_core):
        if depth == 4:
            rel = [n for n, _a, g in self.gen.get(gene, []) if g == seq or g in seq or seq in g]
            base, new, unk = self._bucket(rel)
            if base:
                return "known_in_baseline", sorted(base)
            if new:
                return "now_known", sorted(new)
            if unk:
                return "catalogued_age_unknown", sorted(unk)
            return "still_novel", []
        base, new, unk = self._bucket(self.exact.get(gene, {}).get(seq, []))
        if base:
            return "known_in_baseline", sorted(base)
        prot = protein_of(seq) if depth == 2 else None
        pbase, pnew = [], []
        if prot:
            pbase, pnew, _ = self._bucket(self.prot_index().get(gene, {}).get(prot, []))
            if pbase:
                return "protein_known_in_baseline", sorted(pbase)
        if new:
            return "now_known", sorted(new)
        if unk:
            return "catalogued_age_unknown", sorted(unk)
        if pnew:
            return "now_known_protein", sorted(pnew)
        part = [n for n, _a, c in self.recs.get(gene, [])
                if self.in_baseline(n) is False and 60 <= len(c) < len(seq) and c in seq
                and not (nearest_core and c in nearest_core)]
        if part:
            return "now_known_partial", sorted(part)[:5]
        return "still_novel", []


# ---------------------------------------------------------------------------
# Alignment primitives
# ---------------------------------------------------------------------------
def lev(a, b, band=None):
    """Banded Levenshtein distance (unit costs)."""
    if a == b:
        return 0
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return n + m
    if band is None:
        band = abs(n - m) + 12
    inf = n + m + 1
    prev = [j if j <= band else inf for j in range(m + 1)]
    for i in range(1, n + 1):
        cur = [inf] * (m + 1)
        if i <= band:
            cur[0] = i
        ai = a[i - 1]
        for j in range(max(1, i - band), min(m, i + band) + 1):
            v = prev[j - 1] + (ai != b[j - 1])
            x = prev[j] + 1
            if x < v:
                v = x
            x = cur[j - 1] + 1
            if x < v:
                v = x
            cur[j] = v
        prev = cur
    return prev[m]


def align_events(a, b, extra_band=30):
    """Banded global alignment of `a` (called) to `b` (nearest known). Returns events in `a`
    coordinates, each (pos, n_deleted_from_a, inserted_bases_from_b): applying them to `a` gives
    `b`. SNV = (i, 1, 'X'); bases only in `a` = (i, k, ''); bases only in `b` = (i, 0, 'SEQ')."""
    n, m = len(a), len(b)
    band = abs(n - m) + extra_band
    width = 2 * band + 1
    inf = 10 ** 9
    D, P = [], []
    row, ptr = [inf] * width, [3] * width
    for j in range(0, min(m, band) + 1):
        row[j + band] = j
        ptr[j + band] = 2
    ptr[band] = 3
    D.append(row)
    P.append(ptr)
    for i in range(1, n + 1):
        row, ptr = [inf] * width, [3] * width
        prow = D[i - 1]
        ai = a[i - 1]
        for j in range(max(0, i - band), min(m, i + band) + 1):
            k = j - i + band
            if j == 0:
                best, bp = i, 1
            else:
                best, bp = prow[k] + (ai != b[j - 1]), 0
                if k - 1 >= 0 and row[k - 1] + 1 < best:
                    best, bp = row[k - 1] + 1, 2
            if k + 1 < width and prow[k + 1] + 1 < best:
                best, bp = prow[k + 1] + 1, 1
            row[k], ptr[k] = best, bp
        D.append(row)
        P.append(ptr)
    ops = []
    i, j = n, m
    while i > 0 or j > 0:
        p = P[i][j - i + band]
        if p == 0:
            if a[i - 1] != b[j - 1]:
                ops.append(("X", i - 1, b[j - 1]))
            i, j = i - 1, j - 1
        elif p == 1:
            ops.append(("D", i - 1, ""))
            i -= 1
        elif p == 2:
            ops.append(("I", i, b[j - 1]))
            j -= 1
        else:  # pragma: no cover - band too narrow
            raise ValueError("align_events: traceback left the band")
    ops.reverse()
    events = []
    for op, pos, base in ops:
        if events:
            lp, ld, ls = events[-1]
            if op == "D" and ls == "" and ld > 0 and lp + ld == pos:
                events[-1] = (lp, ld + 1, "")
                continue
            if op == "I" and ld == 0 and lp == pos:
                events[-1] = (lp, 0, ls + base)
                continue
        events.append((pos, 1, base) if op == "X" else ((pos, 1, "") if op == "D" else (pos, 0, base)))
    return events


def apply_events(seq, events):
    out, cur = [], 0
    for p, d, s in sorted(events, key=lambda e: (e[0], 0 if e[1] == 0 else 1)):
        out.append(seq[cur:p])
        out.append(s)
        cur = p + d
    out.append(seq[cur:])
    return "".join(out)


def shift_before(events, x):
    """Length change applied to positions before boundary x (x not inside an event)."""
    return sum(len(s) - d for p, d, s in events if p < x)


def event_hp(seq, ev, min_run=3):
    """Homopolymer context of an indel event on `seq` (the called window). None if not a
    homopolymer indel; else dict(rs, re, base, called_len, nearest_len)."""
    p, d, s = ev
    if d > 0 and s == "":
        bases = seq[p:p + d]
        lo, hi, delta = p, p + d, -d
    elif d == 0 and s:
        bases = s
        lo, hi, delta = p, p, len(s)
    else:
        return None
    if len(set(bases)) != 1:
        return None
    b = bases[0]
    while lo > 0 and seq[lo - 1] == b:
        lo -= 1
    while hi < len(seq) and seq[hi] == b:
        hi += 1
    called_len = hi - lo
    nearest_len = called_len + delta
    if max(called_len, nearest_len) < min_run or min(called_len, nearest_len) < 1:
        return None
    return dict(rs=lo, re=hi, base=b, called_len=called_len, nearest_len=nearest_len)


class Cluster(object):
    __slots__ = ("a", "b", "wc", "tn", "d_cn", "hp", "n_events")


def build_clusters(wc, tn, events, hps, flank):
    spans = []
    for ev, hp in zip(events, hps):
        lo, hi = ev[0], ev[0] + ev[1]
        if hp:
            lo, hi = min(lo, hp["rs"] - HP_FLANK), max(hi, hp["re"] + HP_FLANK)
        spans.append([max(0, lo - flank), min(len(wc), hi + flank), [ev], [hp] if hp else []])
    spans.sort(key=lambda s: s[0])
    merged = []
    for s in spans:
        if merged and s[0] <= merged[-1][1]:
            m = merged[-1]
            m[1] = max(m[1], s[1])
            m[2] += s[2]
            m[3] += s[3]
        else:
            merged.append(s)
    out = []
    for a, b, evs, hp in merged:
        c = Cluster()
        c.a, c.b = a, b
        c.wc = wc[a:b]
        c.tn = tn[a + shift_before(events, a): b + shift_before(events, b)]
        c.d_cn = lev(c.wc, c.tn)
        c.hp = hp
        c.n_events = len(evs)
        out.append(c)
    return out


CS_RE = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+)", re.I)


class CsWalk(object):
    """Target<->query bookkeeping for one mappy hit (short cs; forward strand)."""

    def __init__(self, cs, r_st, q_st):
        self.ts, self.qs, self.tl, self.kind = [], [], [], []
        self.eq, self.en = [], []
        t, q = r_st, q_st
        for tok in CS_RE.findall(cs):
            c = tok[0]
            if c == ":":
                n = int(tok[1:])
                self._seg(t, q, n, "m")
                t += n
                q += n
            elif c == "*":
                self._seg(t, q, 1, "m")
                self._edit(q, 1)
                t += 1
                q += 1
            elif c == "+":
                n = len(tok) - 1
                self._edit(q, n)
                q += n
            elif c == "-":
                n = len(tok) - 1
                self._seg(t, q, n, "d")
                self._edit(q, n)
                t += n
        self.r_st, self.r_en, self.q_st, self.q_en = r_st, t, q_st, q

    def _seg(self, t, q, n, k):
        self.ts.append(t)
        self.qs.append(q)
        self.tl.append(n)
        self.kind.append(k)

    def _edit(self, q, n):
        self.eq.append(q)
        self.en.append(n)

    def q_at(self, t):
        if t == self.r_en:
            return self.q_en
        i = bisect.bisect_right(self.ts, t) - 1
        if i < 0 or t >= self.ts[i] + self.tl[i]:
            return None
        return self.qs[i] if self.kind[i] == "d" else self.qs[i] + (t - self.ts[i])

    def edits_in(self, qa, qb):
        lo = bisect.bisect_left(self.eq, qa)
        hi = bisect.bisect_left(self.eq, qb)
        return sum(self.en[lo:hi])


def genomic_events(wc, gen_seq, aligner_cls=None):
    """Differences between the called window and the nearest known genomic allele (depth 4),
    from minimap2 asm5 (target = window, query = allele). Events in window coordinates."""
    if mappy is None:
        return None, "no_mappy"
    a = (aligner_cls or mappy.Aligner)(seq=wc, preset="asm5")
    hits = [h for h in a.map(gen_seq, cs=True)]
    if not hits:
        return None, "genomic_unaligned"
    h = max(hits, key=lambda x: x.mlen)
    if h.strand != 1:
        return None, "genomic_strand_mismatch"
    events, t = [], h.r_st
    for tok in CS_RE.findall(h.cs):
        c = tok[0]
        if c == ":":
            t += int(tok[1:])
        elif c == "*":
            events.append((t, 1, tok[2].upper()))
            t += 1
        elif c == "+":
            events.append((t, 0, tok[1:].upper()))
        elif c == "-":
            events.append((t, len(tok) - 1, ""))
            t += len(tok) - 1
    return events, None


# ---------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------
def parse_sam_reads(text):
    out = []
    for line in text.splitlines():
        if not line or line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) > 9 and f[9] != "*":
            out.append((f[0], f[9].upper()))
    return out


def fetch_reads(path, region, timeout, samtools="samtools", ref=None):
    """Patchable runner: primary alignments over `region` -> [(name, seq)]. Raises on failure."""
    cmd = [samtools, "view", "-F", "0x904"]
    if ref:
        cmd += ["-T", ref]
    cmd += [path, region]
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=True).stdout
    return parse_sam_reads(out)


def _aligner(seq, preset):
    """mappy index of one window; `map-hifi` needs minimap2 >= 2.19, else fall back to `map-pb`."""
    a = mappy.Aligner(seq=seq, preset=preset)
    if not a and preset == "map-hifi":
        a = mappy.Aligner(seq=seq, preset="map-pb")
    if not a:
        raise RuntimeError(f"mappy index failed (preset {preset})")
    return a


def _best_hit(aligner, seq):
    best = None
    for h in aligner.map(seq, cs=True):
        if best is None or (h.is_primary and not best.is_primary) or h.mlen > best.mlen:
            best = h
    return best


def _oriented_hit(aligner, seq):
    h = _best_hit(aligner, seq)
    if h is None:
        return None, seq
    if h.strand == -1:
        seq = revcomp(seq)
        h = _best_hit(aligner, seq)
        if h is None or h.strand != 1:
            return None, seq
    return h, seq


def hp_readout(seg, a, hp, wc_full):
    """Run length of the homopolymer in a read segment (segment = window [a, ...) of Wc)."""
    lf = wc_full[hp["rs"] - HP_FLANK: hp["rs"]]
    rf = wc_full[hp["re"]: hp["re"] + HP_FLANK]
    if len(lf) < HP_FLANK or len(rf) < HP_FLANK:
        return None
    i = seg.find(lf)
    if i < 0:
        return None
    j = seg.find(rf, i + HP_FLANK)
    if j < 0:
        return None
    mid = seg[i + HP_FLANK: j]
    if mid and set(mid) != {hp["base"]}:
        return "other"
    n = len(mid)
    if n == hp["called_len"]:
        return "called"
    if n == hp["nearest_len"]:
        return "nearest"
    return "other"


def read_support(reads, wc, clusters, wo, preset, p):
    """Votes of each read; returns the aggregate dict (see module docstring step 5-7)."""
    res = dict(n_fetched=len(reads), n_aligned=0, n_spanning=0, n_this=0, n_other=0, n_unknown=0,
               this_called=0, this_nearest=0, pool_called=0, pool_nearest=0,
               hp_called=0, hp_nearest=0, hp_other=0)
    if not reads:
        return res
    ac = _aligner(wc, preset)
    ao = _aligner(wo, preset) if wo else None
    for _name, seq in reads:
        h, seq = _oriented_hit(ac, seq)
        if h is None:
            continue
        res["n_aligned"] += 1
        walk = CsWalk(h.cs, h.r_st, h.q_st)
        diff, spanned, var_q = 0, 0, []
        hp_votes = []
        for cl in clusters:
            if cl.a < h.r_st or cl.b > h.r_en:
                continue
            qa, qb = walk.q_at(cl.a), walk.q_at(cl.b)
            if qa is None or qb is None or qb < qa:
                continue
            seg = seq[qa:qb]
            if seg == cl.wc:
                ec, en = 0, cl.d_cn
            elif seg == cl.tn:
                ec, en = cl.d_cn, 0
            else:
                ec, en = lev(seg, cl.wc), lev(seg, cl.tn)
            diff += en - ec
            spanned += 1
            var_q.append((qa, qb))
            for hp in cl.hp:
                v = hp_readout(seg, cl.a, hp, wc)
                if v:
                    hp_votes.append(v)
        if not spanned:
            continue
        res["n_spanning"] += 1
        hap = "unknown"
        if ao is not None:
            ho = _best_hit(ao, seq)
            if ho is not None and ho.strand == 1:
                wo_walk = CsWalk(ho.cs, ho.r_st, ho.q_st)
                qa, qb = max(h.q_st, ho.q_st), min(h.q_en, ho.q_en)
                if qb - qa >= 0.5 * (h.q_en - h.q_st):
                    e_this = walk.edits_in(qa, qb) - sum(walk.edits_in(max(x, qa), min(y, qb))
                                                         for x, y in var_q if y > qa and x < qb)
                    e_oth = wo_walk.edits_in(qa, qb) - sum(wo_walk.edits_in(max(x, qa), min(y, qb))
                                                           for x, y in var_q if y > qa and x < qb)
                    if e_oth - e_this >= p["hap_margin"]:
                        hap = "this"
                    elif e_this - e_oth >= p["hap_margin"]:
                        hap = "other"
        res["n_" + hap] += 1
        if hap == "other":
            continue
        vote = "called" if diff > 0 else ("nearest" if diff < 0 else None)
        if vote:
            res["pool_" + vote] += 1
            if hap == "this":
                res["this_" + vote] += 1
        for v in hp_votes:
            res["hp_" + v] += 1
    return res


def verdict(res, has_other, min_support, min_frac, min_frac_2t, max_frac_contra):
    """-> (mode, n_called, n_nearest, frac, verdict, reason)."""
    if has_other and res["this_called"] + res["this_nearest"] >= min_support:
        mode, nc, nn, thr = "3t", res["this_called"], res["this_nearest"], min_frac
    else:
        mode, nc, nn, thr = "2t", res["pool_called"], res["pool_nearest"], min_frac_2t
    tot = nc + nn
    frac = nc / tot if tot else float("nan")
    if res.get("n_spanning", 0) == 0:
        return mode, nc, nn, frac, "ambiguous", "no_window_spanned"
    if tot < min_support:
        return mode, nc, nn, frac, "ambiguous", "low_reads"
    if nc >= min_support and frac >= thr:
        return mode, nc, nn, frac, "supported", ""
    if nn >= min_support and frac <= max_frac_contra:
        return mode, nc, nn, frac, "contradicted", ""
    return mode, nc, nn, frac, "ambiguous", "mixed"


# ---------------------------------------------------------------------------
# Per-haplotype Immuannot files
# ---------------------------------------------------------------------------
_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


class HapAssets(object):
    """`<outroot>/<pid>/immuannot_output/<hap>.trimmed.fa`, `<hap>.gtf.gz`, `<hap>/cds.fa.gz`."""

    def __init__(self, outroot, pid, hap):
        base = os.path.join(outroot, str(pid), "immuannot_output")
        self.fa_path = os.path.join(base, f"{hap}.trimmed.fa")
        self.gtf_path = os.path.join(base, f"{hap}.gtf.gz")
        self.cdsfa_path = os.path.join(base, hap, "cds.fa.gz")
        self._contigs = None
        self._cds = None
        self._cdsfa = None

    def contig(self, name):
        if self._contigs is None:
            self._contigs = {}
            if os.path.exists(self.fa_path):
                for head, seq in iter_fasta(self.fa_path):
                    self._contigs[head.split()[0]] = seq
        return self._contigs.get(name)

    def cds_rows(self):
        """{(contig, gene_name, copy): (strand, [(start, end) 1-based])} from GTF CDS rows."""
        if self._cds is None:
            self._cds = {}
            if os.path.exists(self.gtf_path):
                acc = defaultdict(list)
                strand = {}
                with _open_text(self.gtf_path) as fh:
                    for line in fh:
                        if line.startswith("#"):
                            continue
                        f = line.rstrip("\n").split("\t")
                        if len(f) < 9 or f[2] != "CDS":
                            continue
                        at = dict(_ATTR_RE.findall(f[8]))
                        gid = at.get("gene_id", "")
                        copy = int(gid.rsplit(".", 1)[1]) if re.search(r"\.\d+$", gid) else 1
                        k = (f[0], at.get("gene_name", ""), copy)
                        acc[k].append((int(f[3]), int(f[4])))
                        strand[k] = f[6]
                self._cds = {k: (strand[k], sorted(v)) for k, v in acc.items()}
        return self._cds

    def cdsfa(self, contig, gene_raw):
        if self._cdsfa is None:
            self._cdsfa = defaultdict(list)
            if os.path.exists(self.cdsfa_path):
                for head, seq in iter_fasta(self.cdsfa_path):
                    key = head.split()[0]
                    m = re.match(r"^(.*)_(\d+)$", key)
                    if m:
                        self._cdsfa[m.group(1)].append(seq)
        return self._cdsfa.get(f"{contig}_{gene_raw}", [])


def coding_window(ctg, gene_start, gene_end, strand, exons, pad):
    """-> (Wc, cds_idx, genomic_span_coding) with Wc in coding orientation and cds_idx the Wc index
    of every CDS base in coding order."""
    lo = min([gene_start] + [s for s, _ in exons])
    hi = max([gene_end] + [e for _, e in exons])
    ws0, we0 = max(0, lo - 1 - pad), min(len(ctg), hi + pad)
    w = ctg[ws0:we0]
    span = ctg[gene_start - 1: gene_end]
    if strand == "-":
        idx = [(we0 - 1) - (p - 1) for s, e in sorted(exons, reverse=True) for p in range(e, s - 1, -1)]
        return revcomp(w), idx, revcomp(span)
    idx = [(p - 1) - ws0 for s, e in sorted(exons) for p in range(s, e + 1)]
    return w, idx, span


def map_events_to_window(events, cds_idx, core_len):
    """CDS-coordinate events -> window-coordinate events; None if an indel spans an exon junction."""
    out = []
    for p, d, s in events:
        if d > 0:
            ix = cds_idx[p:p + d]
            if len(ix) != d or ix[-1] - ix[0] != d - 1:
                return None
            out.append((ix[0], d, s))
        else:
            if p >= core_len:
                w = cds_idx[core_len - 1] + 1
            elif p > 0 and cds_idx[p] != cds_idx[p - 1] + 1:
                w = cds_idx[p - 1] + 1
            else:
                w = cds_idx[p]
            out.append((w, 0, s))
    return out


# ---------------------------------------------------------------------------
# One call
# ---------------------------------------------------------------------------
PER_CALL_COLS = [
    "person_id", "hap", "contig", "gene", "copy_index", "novelty_depth", "novelty_class", "key",
    "platform", "ancestry", "hp_only_proxy", "status", "seq_sha1", "cds_check", "ipd_class",
    "ipd_match", "nearest_allele", "nearest_source", "n_events", "n_snv", "n_indel", "n_hp_indel",
    "hp_only_seq", "n_clusters", "other_hap_target",
    "lr_n_fetched", "lr_n_aligned", "lr_n_spanning", "lr_n_this", "lr_n_other", "lr_n_unknown",
    "lr_mode", "n_reads_support_called", "n_reads_support_nearest", "frac_called", "lr_verdict",
    "lr_reason", "hp_n_called_len", "hp_n_nearest_len", "hp_n_other_len",
    "sr_n_fetched", "sr_n_spanning", "sr_n_support_called", "sr_n_support_nearest", "sr_frac_called",
    "sr_verdict", "sr_reason", "seconds"]


def _blank_row(call):
    r = {c: "" for c in PER_CALL_COLS}
    for c in ("person_id", "hap", "contig", "gene", "copy_index", "novelty_depth", "novelty_class",
              "key", "platform", "ancestry", "hp_only_proxy"):
        r[c] = call.get(c, "")
    r["lr_verdict"] = "not_run"
    r["sr_verdict"] = "not_run"
    return r


def gene_region(regions, gene, pad):
    g = regions.get(gene)
    if g is None and gene.startswith("DRB"):
        dr = [v for k, v in regions.items() if k.startswith("DR")]
        if dr:
            g = (dr[0][0], min(v[1] for v in dr), max(v[2] for v in dr))
    if g is None:
        return None
    chrom, s, e = g
    return f"{chrom}:{max(1, s + 1 - pad)}-{e + pad}"


def _subsample(reads, n):
    if n and len(reads) > n:
        reads = sorted(reads, key=lambda r: sha1(r[0] + "|64"))[:n]
    return reads


def analyse_call(call, assets, other, ctx):
    t0 = time.time()
    p = ctx["params"]
    row = _blank_row(call)
    gene, depth = call["gene"], int(call["novelty_depth"])
    a = assets
    ctg = a.contig(call["contig"])
    if ctg is None:
        row["status"] = "no_trimmed_contig"
        return row
    k = (call["contig"], call["gene_raw"], int(call["copy_index"]))
    cds = a.cds_rows().get(k)
    if cds is None:
        row["status"] = "no_cds_rows"
        return row
    strand, exons = cds
    wc, cds_idx, gspan = coding_window(ctg, int(call["gene_start"]), int(call["gene_end"]), strand,
                                       exons, p["window_pad"])
    called_cds = "".join(wc[i] for i in cds_idx)
    core = cds_core(called_cds)
    cds_idx = cds_idx[:len(core)]
    fa = [cds_core(s) for s in a.cdsfa(call["contig"], call["gene_raw"])]
    row["cds_check"] = ("no_cds_fa" if not fa else ("agree" if len(fa) == 1 and fa[0] == core else
                                                    ("ambiguous" if len(fa) > 1 else "disagree")))
    ipd = ctx["ipd"]
    seq = gspan if depth == 4 else core
    row["seq_sha1"] = sha1(seq)[:12]
    if depth == 4:
        row["key"] = f"{gene}|gen:{row['seq_sha1']}"
    nname = None
    if isinstance(call.get("cds_mut"), str) and call["cds_mut"]:
        nname = call["cds_mut"].split(",")[0].split("|")[0]
    if not nname:
        nname = call.get("template_allele") if depth == 4 else (str(call.get("alleles") or "").split(",")[0] or None)
    if depth == 4:
        nn, nseq, nsrc = ipd.nearest_gen(gene, call.get("template_allele"))
        ncore = None
    else:
        nn, nseq, nsrc = ipd.nearest_cds(gene, nname, core)
        ncore = nseq
    row["nearest_allele"], row["nearest_source"] = nn or "", nsrc
    cls, names = ipd.classify(gene, depth, seq, ncore)
    row["ipd_class"], row["ipd_match"] = cls, ",".join(names[:5])
    if row["cds_check"] == "disagree" and not p.get("ignore_cds_check"):
        row["status"] = "cds_disagree"
        return row
    if nseq is None:
        row["status"] = nsrc
        return row
    if depth == 4:
        events, why = genomic_events(wc, nseq)
        if events is None:
            row["status"] = why
            return row
    else:
        cev = align_events(core, nseq)
        if len(cev) > p["max_events"]:
            row["status"] = "nearest_too_distant"
            row["n_events"] = len(cev)
            return row
        events = map_events_to_window(cev, cds_idx, len(core))
        if events is None:
            row["status"] = "indel_spans_exon_junction"
            return row
    if not events:
        row["status"] = "no_difference"
        return row
    hps = [event_hp(wc, ev, p["hp_min_run"]) for ev in events]
    row["n_events"] = len(events)
    row["n_snv"] = sum(1 for _p, d, s in events if d == 1 and len(s) == 1)
    row["n_indel"] = len(events) - row["n_snv"]
    row["n_hp_indel"] = sum(1 for h in hps if h)
    row["hp_only_seq"] = all(hps)
    tn = apply_events(wc, events)
    clusters = build_clusters(wc, tn, events, hps, p["flank"])
    row["n_clusters"] = len(clusters)
    wo = None
    if other is not None:
        oa, orow = other
        octg = oa.contig(orow["contig"])
        ok = oa.cds_rows().get((orow["contig"], orow["gene_raw"], 1))
        if octg is not None and ok is not None:
            wo, _ix, _sp = coding_window(octg, int(orow["gene_start"]), int(orow["gene_end"]), ok[0],
                                         ok[1], p["window_pad"])
    row["other_hap_target"] = wo is not None
    row["status"] = "ok"
    if p["no_reads"]:
        row["lr_verdict"] = "not_run"
        row["lr_reason"] = "no_reads_mode"
    elif mappy is None:
        row["lr_reason"] = "no_mappy"
    else:
        bam = call.get("bam")
        region = gene_region(ctx["regions"], gene, p["region_pad"])
        if not bam or not os.path.exists(bam):
            row["lr_reason"] = "bam_missing"
        elif region is None:
            row["lr_reason"] = "no_grch38_region"
        else:
            try:
                reads = _subsample(fetch_reads(
                    bam, region, p["timeout"], p["samtools"], None), p["max_reads"])
                res = read_support(reads, wc, clusters, wo, "map-hifi", p)
                mode, nc, nn_, frac, v, why = verdict(res, wo is not None, p["min_support"], p["min_frac"],
                                                      p["min_frac_2t"], p["max_frac_contra"])
                row.update(lr_n_fetched=res["n_fetched"], lr_n_aligned=res["n_aligned"],
                           lr_n_spanning=res["n_spanning"], lr_n_this=res["n_this"],
                           lr_n_other=res["n_other"], lr_n_unknown=res["n_unknown"], lr_mode=mode,
                           n_reads_support_called=nc, n_reads_support_nearest=nn_,
                           frac_called="" if np.isnan(frac) else round(frac, 3),
                           lr_verdict=v, lr_reason=why, hp_n_called_len=res["hp_called"],
                           hp_n_nearest_len=res["hp_nearest"], hp_n_other_len=res["hp_other"])
            except subprocess.TimeoutExpired:
                row["lr_reason"] = "samtools_timeout"
            except subprocess.CalledProcessError:
                row["lr_reason"] = "samtools_error"
    if p["sr"] and not p["no_reads"] and mappy is not None:
        cram = call.get("cram")
        region = gene_region(ctx["regions"], gene, p["region_pad"])
        if not cram or not os.path.exists(cram):
            row["sr_reason"] = "cram_missing"
        elif region is None:
            row["sr_reason"] = "no_grch38_region"
        else:
            try:
                reads = _subsample(fetch_reads(
                    cram, region, p["timeout"], p["samtools"], p["sr_ref"]), p["sr_max_reads"])
                res = read_support(reads, wc, clusters, None, "sr", p)
                _m, nc, nn_, frac, v, why = verdict(res, False, p["sr_min_reads"], p["min_frac_2t"],
                                                    p["min_frac_2t"], p["max_frac_contra"])
                row.update(sr_n_fetched=res["n_fetched"], sr_n_spanning=res["n_spanning"],
                           sr_n_support_called=nc, sr_n_support_nearest=nn_,
                           sr_frac_called="" if np.isnan(frac) else round(frac, 3), sr_verdict=v,
                           sr_reason=why)
            except subprocess.TimeoutExpired:
                row["sr_reason"] = "samtools_timeout"
            except subprocess.CalledProcessError:
                row["sr_reason"] = "samtools_error"
    row["seconds"] = round(time.time() - t0, 2)
    return row


def process_person(pid, calls, ctx):
    """All selected calls of one person -> rows. Never raises (errors become rows)."""
    rows = []
    assets = {}
    for call in calls:
        try:
            hap = call["hap"]
            if hap not in assets:
                assets[hap] = HapAssets(ctx["outroot"], pid, hap)
            other = None
            ohap = "hap2" if hap == "hap1" else "hap1"
            orow = ctx["t1_idx"].get((pid, ohap, call["gene"]))
            if orow is not None:
                if ohap not in assets:
                    assets[ohap] = HapAssets(ctx["outroot"], pid, ohap)
                other = (assets[ohap], orow)
            rows.append(analyse_call(call, assets[hap], other, ctx))
        except Exception as e:  # noqa: BLE001 - one bad call never stops the person
            r = _blank_row(call)
            r["status"] = f"error:{type(e).__name__}"
            rows.append(r)
            log(f"person {pid} {call.get('hap')} {call.get('gene')}: {type(e).__name__}: {e}")
    return rows


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
T1_COLS = ["person_id", "hap", "contig", "gene", "copy_index", "consensus", "novelty_depth",
           "gene_start", "gene_end", "strand", "cds_mut", "template_warning", "template_allele",
           "alleles"]


def load_t1(path):
    df = pd.read_csv(path, sep="\t", usecols=lambda c: c in T1_COLS, dtype=str)
    df["gene_raw"] = df["gene"]
    df["gene"] = df["gene"].map(C.norm_gene)
    df["copy_index"] = pd.to_numeric(df["copy_index"], errors="coerce").fillna(1).astype(int)
    df["novelty_depth"] = pd.to_numeric(df["novelty_depth"], errors="coerce")
    for c in ("gene_start", "gene_end"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def uri_to_local(uri, root=AOU_ROOT):
    if not isinstance(uri, str) or not uri or uri == "NA":
        return None
    if uri.startswith("gs://"):
        parts = uri.split("/", 3)
        return os.path.join(root, parts[3]) if len(parts) > 3 else None
    return uri if os.path.isabs(uri) else os.path.join(root, uri)


def load_bams(manifest_path, root, platform_of):
    if not manifest_path or not os.path.exists(manifest_path):
        return {}
    man = pd.read_csv(manifest_path, sep="\t", dtype=str).rename(columns={"research_id": "person_id"})
    if "grch38_bam" not in man.columns:
        raise SystemExit(f"FATAL: {manifest_path} has no grch38_bam column")
    out = {}
    for pid, g in man.groupby("person_id"):
        plat = platform_of.get(pid)
        cand = pd.concat([g[g["platform"].astype(str).str.lower() == str(plat)], g]) if "platform" in g else g
        for u in cand["grch38_bam"]:
            loc = uri_to_local(u, root)
            if loc:
                out[pid] = loc
                break
    return out


def load_crams(path, root):
    if not path or not os.path.exists(path):
        return {}
    man = pd.read_csv(path, dtype=str, usecols=lambda c: c in ("person_id", "cram_uri"))
    return {pid: uri_to_local(u, root) for pid, u in zip(man["person_id"], man["cram_uri"])}


def read_regions(path):
    if not path or not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t", header=None, comment="#",
                     names=["chrom", "start", "end", "name", "region_class"])
    df = df[df["region_class"].astype(str).isin(["gene_HLA", "gene"])]
    return {C.norm_gene(n): (c, int(s), int(e)) for c, s, e, n in
            zip(df["chrom"], df["start"], df["end"], df["name"])}


def novel_calls(t1, depths, genes):
    d = t1[t1["novelty_depth"].isin(depths) & (t1["copy_index"] == 1)]
    d = d[d["consensus"].notna() & (d["consensus"] != "undetermined")]
    if genes:
        d = d[d["gene"].isin(genes)]
    d = d.copy()
    d["novelty_class"] = d["novelty_depth"].astype(int).map(DEPTH_CLASS)
    has = d["cds_mut"].notna() & (d["cds_mut"].astype(str) != "")
    d["key"] = np.where(has, d["gene"] + "|" + d["cds_mut"].astype(str),
                        d["gene"] + "|pre:" + d["template_allele"].astype(str) + "|" + d["consensus"].astype(str))
    d["hp_only_proxy"] = d["cds_mut"].map(is_homopolymer_indel_only)
    return d


def select_calls(nov, max_per_key):
    """Round-robin over platform x ancestry strata within each key, fixed hash order."""
    nov = nov.copy()
    nov["_h"] = [sha1(f"{p}|{h}|{k}|64") for p, h, k in zip(nov["person_id"], nov["hap"], nov["key"])]
    picked = []
    for _key, g in nov.groupby("key", sort=True):
        g = g.sort_values("_h")
        strata = [list(s.index) for _, s in g.groupby(["platform", "ancestry"], sort=True)]
        strata.sort(key=lambda ix: g.loc[ix[0], "_h"])
        take, i = [], 0
        while len(take) < max_per_key and any(strata):
            s = strata[i % len(strata)]
            if s:
                take.append(s.pop(0))
            i += 1
        picked.extend(take)
    return nov.loc[sorted(picked)].drop(columns="_h")


# ---------------------------------------------------------------------------
# Checkpoint + lock
# ---------------------------------------------------------------------------
def acquire_lock(internal_dir):
    os.makedirs(internal_dir, exist_ok=True)
    lock = os.path.join(internal_dir, "64.lock")
    if os.path.exists(lock):
        try:
            pid = int(open(lock).read().strip())
            os.kill(pid, 0)
            raise SystemExit(f"FATAL: another 64_novel_allele_validation.py (pid {pid}) holds {lock}")
        except (ValueError, ProcessLookupError):
            log(f"removing stale lock {lock}")
    with open(lock, "w") as fh:
        fh.write(str(os.getpid()))
    return lock


def append_rows(path, rows, cols):
    _disc.assert_internal_path_allowed(path)
    new = not os.path.exists(path)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "a") as fh:
        if new:
            fh.write(_disc.INTERNAL_HEADER_COMMENT + "\n")
            fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in cols) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def read_per_call(path):
    if not os.path.exists(path):
        return pd.DataFrame(columns=PER_CALL_COLS)
    d = pd.read_csv(path, sep="\t", dtype=str, comment="#", keep_default_na=False)
    return d.drop_duplicates(["person_id", "hap", "gene", "copy_index"], keep="last")


def done_people(path, retry_failed=False):
    d = read_per_call(path)
    if d.empty:
        return set()
    if not retry_failed:
        return set(d["person_id"])
    ok = d[~d["status"].str.startswith("error")]
    return set(ok["person_id"])


# ---------------------------------------------------------------------------
# Aggregation (PUBLIC)
# ---------------------------------------------------------------------------
def _mp(n, mode):
    return _disc.mask(int(n), _disc.PARTICIPANT, mode).display


def _mk(n, mode):
    return _disc.mask(int(n), _disc.ALLELE_DISTINCT, mode).display


def _gene_groups(genes):
    out = [(g, [g]) for g in sorted(genes)]
    cl = [g for g in CLASSICAL8 if g in genes]
    if cl:
        out.append(("classical_8", cl))
    out.append(("ALL", sorted(genes)))
    return out


def key_read_class(verdicts):
    v = [x for x in verdicts if x in VERDICTS]
    if not v:
        return "untested"
    s, c = v.count("supported"), v.count("contradicted")
    if s and c:
        return "mixed"
    if s:
        return "supported"
    if c:
        return "contradicted"
    return "ambiguous"


def key_table(nov_unrel, per):
    """One row per key: carriers among unrelated, platforms, ancestries, warning, hp proxy, IPD class
    and key-level read class from its tested calls. Depth-4 keys come from tested calls only."""
    t = per[per["status"] != ""].copy()
    rows = []
    base = nov_unrel[nov_unrel["novelty_depth"] != 4]
    for key, g in base.groupby("key"):
        tg = t[t["key"] == key]
        rows.append(_key_row(key, g, tg))
    t4 = t[t["novelty_depth"].astype(str).isin(["4", "4.0"])]
    for key, tg in t4.groupby("key"):
        g = nov_unrel[nov_unrel.set_index(["person_id", "hap", "gene"]).index.isin(
            list(zip(tg["person_id"], tg["hap"], tg["gene"])))]
        rows.append(_key_row(key, g, tg))
    return pd.DataFrame(rows)


def _key_row(key, g, tg):
    warn = g["template_warning"].notna().any() and (g["template_warning"].astype(str) != "nan").any()
    cls = tg.loc[tg["ipd_class"] != "", "ipd_class"]
    shas = tg.loc[tg["seq_sha1"] != "", "seq_sha1"].nunique()
    depth = int(g["novelty_depth"].iloc[0]) if len(g) else int(float(tg["novelty_depth"].iloc[0]))
    hp_proxy = bool(g["hp_only_proxy"].iloc[0]) if len(g) else False
    if depth == 4 and len(tg):
        hp_proxy = str(tg["hp_only_seq"].iloc[0]) == "True"
    return dict(key=key, gene=(g["gene"].iloc[0] if len(g) else tg["gene"].iloc[0]),
                novelty_depth=depth, novelty_class=DEPTH_CLASS[depth],
                carriers=g["person_id"].nunique(), platforms=g["platform"].nunique(),
                ancestries=g["ancestry"].nunique(), warn=bool(warn), hp_only=hp_proxy,
                n_tested=len(tg), n_seq_variants=shas,
                ipd_class=cls.mode().iloc[0] if len(cls) else "no_sequence",
                ipd_match=",".join(sorted({m for s in tg["ipd_match"] for m in str(s).split(",") if m}))[:200],
                read_class=key_read_class(list(tg["lr_verdict"])),
                n_supported=int((tg["lr_verdict"] == "supported").sum()),
                n_contradicted=int((tg["lr_verdict"] == "contradicted").sum()))


FUNNEL_STEPS = ["distinct", "ge2_unrelated_carriers", "not_homopolymer_only", "no_template_warning",
                "ge2_platforms", "ge2_ancestries", "ge20_carriers_clean"]


def funnel(keys, mode):
    """50d steps (same rules) per gene group and novelty class, with evidence columns."""
    rows = []
    if keys.empty:
        return pd.DataFrame()
    for ncls, kc in keys.groupby("novelty_class"):
        for label, genes in _gene_groups(set(kc["gene"])):
            k0 = kc[kc["gene"].isin(genes)]
            f1 = k0[k0["carriers"] >= 2]
            f2 = f1[~f1["hp_only"]]
            f3 = f2[~f2["warn"]]
            f4 = f3[f3["platforms"] >= 2]
            f5 = f4[f4["ancestries"] >= 2]
            f6 = f3[f3["carriers"] >= 20]
            for step, f in zip(FUNNEL_STEPS, (k0, f1, f2, f3, f4, f5, f6)):
                row = dict(novelty_class=ncls, gene=label, step=step, n_keys=_mk(len(f), mode))
                for c in ("known_in_baseline", "protein_known_in_baseline") + NOW_KNOWN:
                    row["n_" + c] = _mk(int((f["ipd_class"] == c).sum()), mode)
                for c in ("supported", "contradicted", "mixed", "ambiguous", "untested"):
                    row["n_read_" + c] = _mk(int((f["read_class"] == c).sum()), mode)
                rows.append(row)
    return pd.DataFrame(rows)


def stratum_summary(nov_unrel, per, keys, mode):
    """gene group x novelty class x platform x ancestry (with ALL) x hp split."""
    key_ipd = dict(zip(keys["key"], keys["ipd_class"])) if len(keys) else {}
    t = per[per["status"] != ""].copy()
    t["novelty_depth"] = pd.to_numeric(t["novelty_depth"], errors="coerce")
    tested = t[t["lr_verdict"].isin(VERDICTS)]
    calls = nov_unrel[["person_id", "hap", "gene", "key", "novelty_class", "platform", "ancestry",
                       "hp_only_proxy"]].copy()
    # depth-4 calls: key and hp flag from the tested row (sequence only known there)
    d4 = t[t["novelty_depth"] == 4].set_index(["person_id", "hap", "gene"])
    if len(d4):
        ix = list(zip(calls["person_id"], calls["hap"], calls["gene"]))
        k4 = d4["key"].to_dict()
        h4 = d4["hp_only_seq"].to_dict()
        calls["key"] = [k4.get(i, k) for i, k in zip(ix, calls["key"])]
        calls["hp_only_proxy"] = [(str(h4[i]) == "True") if i in h4 else h for i, h in zip(ix, calls["hp_only_proxy"])]
    calls = calls.merge(tested[["person_id", "hap", "gene", "lr_verdict"]], on=["person_id", "hap", "gene"],
                        how="left")
    plats = sorted(calls["platform"].dropna().unique())
    strata = [(pl, an) for pl in plats + ["ALL"] for an in ANCESTRY_ORDER + ["ALL"]]
    rows = []
    for ncls, cc in calls.groupby("novelty_class"):
        for label, genes in _gene_groups(set(cc["gene"])):
            cg = cc[cc["gene"].isin(genes)]
            for pl, an in strata:
                s = cg
                if pl != "ALL":
                    s = s[s["platform"] == pl]
                if an != "ALL":
                    s = s[s["ancestry"] == an]
                if s.empty:
                    continue
                for split in ("all", "hp_only", "other"):
                    x = s if split == "all" else s[s["hp_only_proxy"] == (split == "hp_only")]
                    if x.empty:
                        continue
                    ks = [k for k in x["key"].unique() if not str(k).split("|")[1].startswith("pre:")]
                    kcls = [key_ipd.get(k, "no_sequence") for k in ks]
                    nk = len(ks)
                    nt = int(x["lr_verdict"].notna().sum())
                    row = dict(novelty_class=ncls, gene=label, platform=pl, ancestry=an, hp_split=split,
                               n_keys=_mk(nk, mode), n_calls=_mp(len(x), mode), n_calls_tested=_mp(nt, mode))
                    for c in ("known_in_baseline", "protein_known_in_baseline") + NOW_KNOWN:
                        n = kcls.count(c)
                        row["pct_keys_" + c] = round(100 * n / nk, 2) if nk else ""
                    for v in VERDICTS:
                        n = int((x["lr_verdict"] == v).sum())
                        r = _disc.mask_rate(n, nt, mode, pct_decimals=1)
                        row["n_" + v] = r["n"]
                        row["pct_" + v] = r["pct"]
                    rows.append(row)
    return pd.DataFrame(rows)


def now_known_table(keys, mode):
    k = keys[keys["ipd_class"].isin(NOW_KNOWN)]
    rows = []
    for (gene, ncls, cls, match), g in k.groupby(["gene", "novelty_class", "ipd_class", "ipd_match"]):
        car = int(g["carriers"].sum())
        rows.append(dict(gene=gene, novelty_class=ncls, ipd_class=cls, ipd_alleles=match,
                         n_keys=_mk(len(g), mode), n_unrelated_carriers=_mp(car, mode),
                         review_flag=bool(car < _disc.SUPPRESS_BELOW),
                         hp_only_proxy=",".join(sorted({str(b) for b in g["hp_only"]})),
                         read_class=",".join(sorted(set(g["read_class"])))))
    return pd.DataFrame(rows, columns=["gene", "novelty_class", "ipd_class", "ipd_alleles", "n_keys",
                                       "n_unrelated_carriers", "review_flag", "hp_only_proxy", "read_class"])


def hp_runlength_table(per, mode):
    t = per[per["status"] == "ok"].copy()
    rows = []
    if t.empty:
        return pd.DataFrame()
    for c in ("hp_n_called_len", "hp_n_nearest_len", "hp_n_other_len"):
        t[c] = pd.to_numeric(t[c], errors="coerce").fillna(0).astype(int)
    for (gene, ncls), g in t.groupby(["gene", "novelty_class"]):
        tot = int(g[["hp_n_called_len", "hp_n_nearest_len", "hp_n_other_len"]].sum().sum())
        if tot == 0:
            continue
        row = dict(gene=gene, novelty_class=ncls, n_calls_with_hp=_mp(int((g["n_hp_indel"].astype(str) != "0").sum()), mode),
                   n_read_observations=_mp(tot, mode))
        for c, lab in (("hp_n_called_len", "called_len"), ("hp_n_nearest_len", "nearest_len"),
                       ("hp_n_other_len", "other_len")):
            r = _disc.mask_rate(int(g[c].sum()), tot, mode, pct_decimals=1)
            row["pct_" + lab] = r["pct"]
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fetch-ipd", action="store_true", help="download the IPD release into --ipd-dir, then exit")
    ap.add_argument("--ipd-release", default="Latest", help="git ref of ANHIG/IMGTHLA (branch or tag)")
    ap.add_argument("--ipd-dir", help="directory with hla_nuc.fasta [hla_gen.fasta, Allelelist*.txt]")
    ap.add_argument("--baseline-col", default=BASELINE_COL, help="Allelelist_history.txt release column")
    ap.add_argument("--refdata", default=DEFAULT_REFDATA, help="Immuannot refdata (3.55.0 fallback)")
    ap.add_argument("--table1", default=C.DEFAULT_TABLE1)
    ap.add_argument("--cohort", default=C.DEFAULT_COHORT)
    ap.add_argument("--relatedness", default=C.DEFAULT_RELATEDNESS)
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT)
    ap.add_argument("--regions-bed", default=DEFAULT_REGIONS, help="58's regions.bed (GRCh38 gene spans)")
    ap.add_argument("--aou-root", default=AOU_ROOT)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    ap.add_argument("--internal-dir", default=DEFAULT_INTERNAL_DIR)
    ap.add_argument("--depths", default="2,3")
    ap.add_argument("--genes", default="", help="comma list of bare gene names (default: all)")
    ap.add_argument("--max-per-key", type=int, default=5)
    ap.add_argument("--window-pad", type=int, default=2000)
    ap.add_argument("--region-pad", type=int, default=1000)
    ap.add_argument("--flank", type=int, default=25)
    ap.add_argument("--hp-min-run", type=int, default=3)
    ap.add_argument("--max-events", type=int, default=40)
    ap.add_argument("--max-reads", type=int, default=80)
    ap.add_argument("--min-support", type=int, default=3)
    ap.add_argument("--min-frac", type=float, default=0.7)
    ap.add_argument("--min-frac-2t", type=float, default=0.3)
    ap.add_argument("--max-frac-contradicted", type=float, default=0.15)
    ap.add_argument("--hap-margin", type=int, default=1)
    ap.add_argument("--no-reads", action="store_true", help="IPD re-match only (no BAM access)")
    ap.add_argument("--ignore-cds-check", action="store_true",
                    help="test calls whose GTF-rebuilt CDS differs from cds.fa.gz (default: skip them)")
    ap.add_argument("--sr", action="store_true", help="also check srWGS CRAM reads")
    ap.add_argument("--sr-ref", default=DEFAULT_SR_REF)
    ap.add_argument("--sr-min-reads", type=int, default=3)
    ap.add_argument("--sr-max-reads", type=int, default=200)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--timeout", type=int, default=600, help="seconds per samtools call")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, help="first N people with selected calls")
    ap.add_argument("--retry-failed", action="store_true")
    ap.add_argument("--aggregate-only", action="store_true")
    _disc.add_disclosure_arg(ap)
    args = ap.parse_args(argv)
    ipd_dir = args.ipd_dir or os.path.expanduser(f"~/s05/ref/ipd_{args.ipd_release.lower()}")

    if args.fetch_ipd:
        got = fetch_ipd(args.ipd_release, ipd_dir)
        log(f"IPD release {args.ipd_release}: {read_version(ipd_dir) or 'version line not found'}; {got}")
        return 0

    _disc.assert_internal_path_allowed(args.internal_dir)
    depths = [int(x) for x in args.depths.split(",") if x]
    genes = [C.norm_gene(g) for g in args.genes.split(",") if g]

    hist = os.path.join(ipd_dir, "Allelelist_history.txt")
    base = load_baseline_history(hist, args.baseline_col)
    bsrc = f"Allelelist_history.txt column {args.baseline_col}"
    if base is None:
        base, bsrc = load_baseline_refdata(args.refdata), "immuannot_refdata_CDSseq"
    if base is None:
        bsrc = "none"
        log("WARNING: no 3.55.0 baseline list; exact matches will be 'catalogued_age_unknown'")
    ipd = IpdRelease(ipd_dir, base, bsrc, refdata=args.refdata, load_gen=4 in depths)
    log(f"IPD {args.ipd_release} ({ipd.version or 'version n/a'}): {ipd.n_nuc} CDS records, "
        f"{ipd.n_gen} genomic; baseline from {bsrc}")

    cohort = C.load_cohort(args.cohort)
    platform_of = dict(zip(cohort["person_id"], cohort["platform"]))
    ancestry_of = dict(zip(cohort["person_id"], cohort["ancestry_pred"]))
    t1 = load_t1(args.table1)
    unrelated = C.greedy_unrelated(t1["person_id"].unique(), C.load_relatedness_pairs(args.relatedness))
    nov = novel_calls(t1, depths, genes)
    nov["platform"] = nov["person_id"].map(platform_of).fillna("unknown")
    nov["ancestry"] = nov["person_id"].map(ancestry_of).fillna("UNKNOWN")
    nov_unrel = nov[nov["person_id"].isin(unrelated)].copy()
    sel = select_calls(nov_unrel, args.max_per_key)
    people = sorted(sel["person_id"].unique())
    if args.limit:
        people = people[:args.limit]
        sel = sel[sel["person_id"].isin(set(people))]
        nov_unrel = nov_unrel[nov_unrel["person_id"].isin(set(people))]
    log(f"{len(nov)} novel calls (depths {depths}); {len(nov_unrel)} in unrelated people; "
        f"{nov_unrel['key'].nunique()} keys; {len(sel)} calls selected in {len(people)} people")

    per_path = os.path.join(args.internal_dir, "per_call_evidence.tsv")
    t0 = time.time()
    if not args.aggregate_only:
        manifest = os.path.join(args.aou_root, LR_MANIFEST_REL)
        if not args.no_reads and not os.path.exists(manifest):
            raise SystemExit(f"FATAL: {manifest} not found -- is the controlled-data automount up? "
                             f"(ENVIRONMENT.md quirks #26/#35). Refusing to start (or use --no-reads).")
        regions = read_regions(args.regions_bed)
        if not args.no_reads and not regions:
            raise SystemExit(f"FATAL: no gene regions in {args.regions_bed} (build with 58 --build-regions)")
        bams = {} if args.no_reads else load_bams(manifest, args.aou_root, platform_of)
        crams = {}
        if args.sr:
            cm = os.path.join(args.aou_root, SR_MANIFEST_REL)
            if not os.path.exists(cm) or not os.path.exists(args.sr_ref):
                raise SystemExit(f"FATAL: --sr needs {cm} and --sr-ref {args.sr_ref}")
            crams = load_crams(cm, args.aou_root)
        t1_idx = {}
        sub = t1[t1["person_id"].isin(set(people)) & (t1["copy_index"] == 1)]
        for r in sub.to_dict("records"):
            t1_idx.setdefault((r["person_id"], r["hap"], r["gene"]), r)
        params = dict(window_pad=args.window_pad, region_pad=args.region_pad, flank=args.flank,
                      hp_min_run=args.hp_min_run, max_events=args.max_events, max_reads=args.max_reads,
                      min_support=args.min_support, min_frac=args.min_frac, min_frac_2t=args.min_frac_2t,
                      max_frac_contra=args.max_frac_contradicted, hap_margin=args.hap_margin,
                      no_reads=args.no_reads, ignore_cds_check=args.ignore_cds_check, sr=args.sr, sr_ref=args.sr_ref, sr_min_reads=args.sr_min_reads,
                      sr_max_reads=args.sr_max_reads, samtools=args.samtools, timeout=args.timeout)
        ctx = dict(outroot=args.outroot, ipd=ipd, regions=regions, t1_idx=t1_idx, params=params)
        lock = acquire_lock(args.internal_dir)
        try:
            done = done_people(per_path, args.retry_failed)
            todo = [p for p in people if p not in done]
            log(f"{len(done)} people already done (resume), {len(todo)} to test, {args.workers} workers")
            by_pid = {pid: g.to_dict("records") for pid, g in sel[sel["person_id"].isin(set(todo))].groupby("person_id")}
            for pid, calls in by_pid.items():
                for c in calls:
                    c["bam"] = bams.get(pid)
                    c["cram"] = crams.get(pid)
            n_calls = 0
            with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
                futs = {ex.submit(process_person, pid, by_pid.get(pid, []), ctx): pid for pid in todo}
                for i, fut in enumerate(as_completed(futs), 1):
                    pid = futs[fut]
                    try:
                        rows = fut.result()
                    except Exception as e:  # noqa: BLE001 - belt and braces; process_person never raises
                        log(f"person {pid}: {type(e).__name__}: {e}")
                        rows = []
                    if rows:
                        append_rows(per_path, rows, PER_CALL_COLS)
                        n_calls += len(rows)
                    if i % 100 == 0 or i == len(todo):
                        log(f"{i}/{len(todo)} people done, {n_calls} calls, "
                            f"{i / max(1e-9, time.time() - t0):.2f} people/s")
        finally:
            os.remove(lock)

    per = read_per_call(per_path)
    per = per[per["person_id"].isin(set(people))]
    mode = args.disclosure
    os.makedirs(args.out_dir, exist_ok=True)
    keys = key_table(nov_unrel, per)
    _disc.write_tsv(keys, os.path.join(args.internal_dir, "key_evidence.tsv"), _disc.INTERNAL)
    rel = [dict(item="ipd_repository", value=IPD_REPO), dict(item="ipd_release_ref", value=args.ipd_release),
           dict(item="ipd_version", value=ipd.version), dict(item="baseline", value=f"3.55.0 via {bsrc}"),
           dict(item="n_cds_records", value=ipd.n_nuc), dict(item="n_genomic_records", value=ipd.n_gen)]
    rel += [dict(item=f"gene_{r['gene']}", value=f"{r['n_alleles_nuc']} alleles; "
                 f"{r['n_added_since_baseline']} added since baseline") for r in ipd.new_counts()]
    _disc.write_tsv(pd.DataFrame(rel), os.path.join(args.out_dir, "ipd_release.tsv"), mode)
    _disc.write_tsv(stratum_summary(nov_unrel, per, keys, mode),
                    os.path.join(args.out_dir, "summary_by_stratum.tsv"), mode)
    _disc.write_tsv(funnel(keys, mode), os.path.join(args.out_dir, "funnel_with_evidence.tsv"), mode)
    _disc.write_tsv(now_known_table(keys, mode), os.path.join(args.out_dir, "now_known_alleles.tsv"), mode)
    _disc.write_tsv(hp_runlength_table(per, mode), os.path.join(args.out_dir, "hp_runlength_by_gene.tsv"), mode)
    st = per["status"].where(~per["status"].str.startswith("error"), "error").value_counts()
    lr = per.loc[per["status"] == "ok", "lr_reason"].replace("", "ran").value_counts()
    fail = [dict(stage="call_status", reason=k, n_calls=_mp(v, mode)) for k, v in sorted(st.items())]
    fail += [dict(stage="long_read", reason=k, n_calls=_mp(v, mode)) for k, v in sorted(lr.items())]
    _disc.write_tsv(pd.DataFrame(fail, columns=["stage", "reason", "n_calls"]),
                    os.path.join(args.out_dir, "failures_summary.tsv"), mode)
    run = [dict(metric="people_with_selected_calls", value=_mp(len(people), mode)),
           dict(metric="people_tested", value=_mp(per["person_id"].nunique(), mode)),
           dict(metric="novel_calls_unrelated", value=_mp(len(nov_unrel), mode)),
           dict(metric="calls_selected", value=_mp(len(sel), mode)),
           dict(metric="calls_with_evidence_rows", value=_mp(len(per), mode)),
           dict(metric="keys", value=_mk(len(keys), mode)),
           dict(metric="keys_with_inconsistent_sequence",
                value=_mk(int((keys["n_seq_variants"] > 1).sum()) if len(keys) else 0, mode)),
           dict(metric="cds_check_disagree", value=_mp(int((per["cds_check"] == "disagree").sum()), mode)),
           dict(metric="depths", value=args.depths), dict(metric="genes", value=args.genes or "all"),
           dict(metric="max_per_key", value=args.max_per_key), dict(metric="min_support", value=args.min_support),
           dict(metric="min_frac", value=args.min_frac), dict(metric="min_frac_2t", value=args.min_frac_2t),
           dict(metric="max_frac_contradicted", value=args.max_frac_contradicted),
           dict(metric="flank", value=args.flank), dict(metric="window_pad", value=args.window_pad),
           dict(metric="hp_min_run", value=args.hp_min_run), dict(metric="sr", value=args.sr),
           dict(metric="no_reads", value=args.no_reads), dict(metric="limit", value=args.limit or ""),
           dict(metric="mappy", value=getattr(mappy, "__version__", "absent") if mappy else "absent"),
           dict(metric="wall_seconds_this_invocation", value=int(time.time() - t0))]
    _disc.write_tsv(pd.DataFrame(run), os.path.join(args.out_dir, "run_summary.tsv"), mode)
    log(f"done: {len(per)} call rows, {len(keys)} keys; tables in {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
