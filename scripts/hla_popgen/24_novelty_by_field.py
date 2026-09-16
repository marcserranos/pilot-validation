#!/usr/bin/env python3
"""Novelty re-defined by nomenclature field, with per-call artifact labels, a sequence-level truth
check against the IPD-IMGT CDS snapshot Immuannot ships, protein-level clustering, an unrelated
subset and strict ancestry. Spec: sprints/S01_novelty_qc_coverage/WS1_novelty_by_field.md
("script 24").

## Why this script exists

The WS1 audit found three problems with the earlier novelty numbers:
- 03 clusters novel calls by the observed-CDS hash, so a depth-4 (non-coding) "cluster" pools every
  haplotype whose CDS matches a known allele. Its person count is not the recurrence of one allele.
- 87% of the "protein_altering" clusters carry an artifact flag (mostly homopolymer indels, the
  HiFi error mode). The "90.8% protein-altering" headline was computed before that gate.
- 04's saturation curves counted flagged clusters and did not remove relatives.

This script counts novelty per haplotype-gene call rather than per cluster, labels artifacts per
call (flags are labels, not deletions), checks depth-2/3 calls against the actual sequences, and
writes allele-count tables at protein and CDS resolution for WS3.

## Method

1. **Field class per call** (Table 1, `copy_index == 1`, all genes). The field class comes from
   `novelty_depth`: 1 -> `f1_undetermined`, 2 -> `f2_protein`, 3 -> `f3_synonymous`,
   4 -> `f4_noncoding`. A non-novel call is `known`, and `consensus == undetermined` (or empty) is
   `uncalled`. When `novelty_depth` is missing, the depth is re-derived from `consensus` with the
   same rule as `01_extract_rich.parse_consensus` (1-based index of `new` among the fields after
   `*`). So `HLA-A*01:new` has depth 2 (new protein) and `HLA-A*new` has depth 1.
   **Note:** the example strings in the WS1 brief's background table are shifted by one field
   relative to 01's parser and SCHEMA.md. This script follows SCHEMA.md and 01.
2. **Artifact label per call**, checked in this order: `homopolymer_indel` (03's rule B,
   `is_homopolymer_indel_only`, applied to this call's own `cds_mut`), then `partial_cds` and
   `inframe_stop` (03's `warning_tokens`, which read the literal `NA` as clean), else `clean`.
3. **Unrelated set:** relative pairs with kin >= 0.0442 (third degree or closer) where both people
   are in the long-read cohort. The script greedily removes the person with the most remaining
   relatives, breaking ties by the smallest id (string sort).
4. **Ancestry:** `ancestry_pred`, plus `strict` (the same label, kept only if that ancestry's
   `p_<anc>` is >= 0.9; otherwise UNASSIGNED).
5. **Sequence truth check (depth 2/3):** 03's `match_novel_rows` is run once per person (threaded)
   and reads each `cds.fa.gz` once. The observed CDS is compared with every known CDS for the
   gene. Both sides first have a terminal stop codon removed, if present (`cds_core`). The
   translations are compared the same way. Each call gets one class:
   `cds_known` > `frameshift_or_stop` > (`protein_known` if depth 2 | `novel_cds_synonymous` if
   depth 3) > `novel_protein`. A depth-3 call whose protein is not catalogued becomes
   `novel_protein`. That disagreement is counted and reported. A gene with no CDS file in the
   refdata gets `no_refdata`.
6. **Clusters:** `novel_protein` calls are clustered by protein sequence. Calls with a new CDS but
   a known protein (`novel_cds_synonymous` plus the depth-2 `protein_known` calls) are clustered by
   `cds_core`. The nearest known protein uses Hamming distance against equal-length known
   proteins (numpy). If no known protein has the same length, it uses Levenshtein distance
   against a shortlist: the call's own tied candidates plus the 5 known proteins closest in
   length. Differences are written as `pos:ref>alt`, 1-based on the nearest known precursor
   protein. `-` marks an indel, and an insertion is anchored to the preceding reference residue.
7. **Groove annotation** (`annotate_groove`, isolated): a codon counts as in the groove if its
   middle base falls in exon 2 or 3 for class I (A, B, C, E, F, G) or in exon 2 for the class II
   alpha/beta chains (DRA, DRB*, DQA*, DQB*, DPA*, DPB*). The exon structure is derived from
   `<refdata>/alleles.csv.gz`: each `CDS=` range is assigned to the exon whose `exon=` range
   contains its start, giving per-exon coding lengths. The nearest known allele's own record is
   used when present (`refdata_exon`), else the gene's most common pattern among its
   `goodTemplate` alleles (`refdata_gene_mode`). A hard-coded HLA-A/B/C map
   (`APPROX_CANONICAL_CDS`, `approx_canonical`) is used only if the gene is absent from
   alleles.csv.gz entirely, and only when both CDS have the canonical length. Every other case is
   NA with the reason in `groove_na_reason`. The IPD-IMGT/HLA `version=` string from
   alleles.csv.gz is reported.
8. **Release drift (optional):** `--imgt-latest-nuc hla_nuc.fasta` counts how many novel CDS and
   protein clusters appear in a newer IPD-IMGT/HLA release.
9. **WS3 allele-count tables** (VM-local): allele counts per gene x resolution (protein|cds) x
   exclude_flagged (all|clean_only) x ancestry scheme x ancestry x unrelated_only. Rows that are
   uncalled, depth 1, unmatched or `no_refdata` are dropped, and their numbers are reported.
   `frameshift_or_stop` calls are kept in `all` under a `<gene>_nonfunc_<sha8>` id and removed
   from `clean_only`.

## Assumptions about the real files (confirmed by the coordinator on the VM, 2026-09-16)

- `<refdata>/CDSseq/<gene>.fa.gz` (62 files, e.g. `HLA-A.fa.gz`): header
  `>HLA-A*01:01:01:01 HLA00001 frame=1 1098bp`. The parser takes the first whitespace token that
  contains `*` (falling back to the first token and the file stem), so an absent `HLA-` prefix is
  fine. Records with `frame=` other than 1 are kept in the CDS index but left out of the protein
  index, since they are not translatable from base 1.
- `<refdata>/alleles.csv.gz` (41,128 lines): tab-separated `key=value` with `gene`, `version`,
  `allele`, `goodTemplate`, `CDS=`, `exon=` in 1-based genomic-allele coordinates.
- Relatedness: `~/workspace/vwb-aou-datasets-controlled-v9/...`/`samples_relatedness.tsv` with
  `i.s`, `j.s`, `kin` (per 11). **The old `~/mnt/aou-controlled` path is stale and hangs; this
  script never touches it.** Ids are assumed to match Table 1 `person_id` as strings.
- `cohort_membership.tsv` and Table 1 carry exactly the SCHEMA.md columns.
- `cds.fa.gz` headers are `{contig}_{gene}_{i} GT_AG:...` (03/17).

**Still unverified and worth a one-minute check on the VM:** whether the CDSseq sequences include
the terminal stop codon (this script compares both sides stop-stripped, so it works either way, and
the count of refdata records ending in a stop is reported), and whether any CDSseq record is a
partial (exon-2/3-only) sequence, which would only ever fail to match, never match wrongly.

## Known limitations

- Depth-4 novelty cannot be resolved here (Table 1 has no genomic sequence). See script 25.
- `(person, hap, gene)` keys with more than one Table 1 row (copy number > 1 or a gene split
  across contigs) are excluded from the sequence check and counted. They can't be joined
  unambiguously.
- Rule B is 03's proxy: it sees only the indel run, not the flanking homopolymer context.
- A flag does not prove an error. Both the `all` and `clean` views are reported.

## Outputs

Committed (`--out-dir`, aggregate only, person or haplotype counts from 1 to 19 written as `<20`):
`novelty_field_counts.tsv`, `sequence_class_counts.tsv`, `protein_level_novel_clusters.tsv`,
`novelty_by_field_report.md`, `summary.json`, `fig_field_class_by_gene.png`,
`fig_clean_novelty_rate_by_ancestry.png`, `fig_novelty_funnel.png`.
VM-local (`--local-dir`, never commit): `unrelated_set.tsv` (the only file with person ids),
`novelty_field_counts.full.tsv`, `protein_level_novel_clusters.full.tsv`,
`allele_counts_by_resolution.tsv`, `novel_protein_clusters.faa`.
With `--limit N`, both directories get a `pilot_limitN/` subdirectory, so a pilot never
overwrites a full run (SCHEMA.md hard rule 4).

Usage:
    # full run (VM)
    pixi run -e spechla -- python3 scripts/hla_popgen/24_novelty_by_field.py \\
        --out-dir ~/results/24_novelty_by_field
    # pilot: same command restricted to the first N people (own output subdirectory)
    pixi run -e spechla -- python3 scripts/hla_popgen/24_novelty_by_field.py \\
        --out-dir ~/results/24_novelty_by_field --limit 200
    # tests (local)
    cd scripts/hla_popgen && python3 tests/test_novelty_by_field.py
"""
import argparse
import concurrent.futures
import glob
import gzip
import hashlib
import importlib.util
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _viz_common as vc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_module(filename, modname):
    spec = importlib.util.spec_from_file_location(modname, os.path.join(_HERE, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


novel03 = _load_module("03_novel_alleles.py", "hla_popgen_novel_alleles_for24")
parse_cds_fasta = novel03.parse_cds_fasta  # re-exported (used via match_novel_rows)
match_novel_rows = novel03.match_novel_rows
warning_tokens = novel03.warning_tokens
is_homopolymer_indel_only = novel03.is_homopolymer_indel_only

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DATA_ROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_TABLE1 = os.path.join(DATA_ROOT, "hla_calls_rich.tsv")
DEFAULT_COHORT = os.path.join(DATA_ROOT, "cohort_membership.tsv")
DEFAULT_OUTROOT = os.path.join(DATA_ROOT, "people")
DEFAULT_LOCAL_DIR = os.path.join(DATA_ROOT, "s01_local")
DEFAULT_REFDATA = os.path.expanduser("~/tools/Immuannot_refdata")
DEFAULT_RELATEDNESS = os.path.expanduser(
    "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/"
    "samples_relatedness.tsv")  # the old ~/mnt/aou-controlled mount is stale and hangs: never use
DEFAULT_OUT_DIR = os.path.join(vc.DEFAULT_REPORT_ROOT, "24_novelty_by_field")

KIN_MIN = 0.0442  # third degree or closer (KING; 11_relatedness_cohort_overlap.KIN_BINS)
STRICT_MIN = 0.9
SUPPRESS_BELOW = 20

FIELD_CLASSES = ["known", "f4_noncoding", "f3_synonymous", "f2_protein", "f1_undetermined",
                 "uncalled"]
DEPTH_TO_CLASS = {1: "f1_undetermined", 2: "f2_protein", 3: "f3_synonymous", 4: "f4_noncoding"}
ARTIFACT_LABELS = ["clean", "homopolymer_indel", "partial_cds", "inframe_stop"]
SEQ_CLASSES = ["novel_protein", "novel_cds_synonymous", "protein_known", "cds_known",
               "frameshift_or_stop", "no_refdata"]
GENE_CLASS_ORDER = ["classical_I", "classical_II", "nonclassical_I", "class_II_accessory",
                    "class_II_paralog", "pseudogene_I", "mic_tap", "complement", "other", "kir"]

# APPROXIMATE fallback exon map, used only when alleles.csv.gz has no entry at all for the gene. Standard IPD-IMGT/HLA reference CDS for HLA-A/B/C: exon 1 = 73 nt (leader), exon 2 = 270,
# exon 3 = 276. The total is the CDS length WITHOUT the stop codon (A*01:01 1098 nt incl. stop,
# B*07:02 1089, C*01:02 1101). Applied only when BOTH the observed and the reference cds_core
# have exactly this length. Class II is deliberately left out: exon-1 CDS length varies by locus
# and is not known well enough here to hard-code.
APPROX_CANONICAL_CDS = {
    "A": {"core_len": 1095, "exon_cds_lengths": [73, 270, 276]},
    "B": {"core_len": 1086, "exon_cds_lengths": [73, 270, 276]},
    "C": {"core_len": 1098, "exon_cds_lengths": [73, 270, 276]},
}

NO_PREFIX_GENES = {"MICA", "MICB", "TAP1", "TAP2", "C4A", "C4B"}

TABLE1_COLS = ["person_id", "hap", "contig", "gene", "copy_index", "gene_class", "consensus",
               "is_novel", "novelty_depth", "novelty_class", "template_warning", "cds_mut",
               "alleles", "template_allele", "cds_distance", "n_aa_changes"]

# Okabe-Ito-derived field-class palette (distinct from nothing else on the same figure).
FIELD_COLORS = {
    "known": "#DDDDDD", "f4_noncoding": "#56B4E9", "f3_synonymous": "#009E73",
    "f2_protein": "#D55E00", "f1_undetermined": "#CC79A7", "uncalled": "#555555",
}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def suppress(n, threshold=SUPPRESS_BELOW):
    """AoU small-cell rule for committed tables: a count from 1 to threshold-1 becomes '<20'.
    Zero stays '0' and NA stays 'NA'. Rates and fractions are not passed through this."""
    if n is None or (isinstance(n, float) and np.isnan(n)) or n is pd.NA:
        return "NA"
    n = int(n)
    if 1 <= n < threshold:
        return f"<{threshold}"
    return str(n)


def suppress_df(df, cols, threshold=SUPPRESS_BELOW):
    out = df.copy()
    for c in cols:
        if c in out.columns:
            out[c] = out[c].map(lambda v: suppress(v, threshold)).astype(object)
    return out


def gene_bare(g):
    if not isinstance(g, str):
        return g
    return g[4:] if g.startswith("HLA-") else g


def gene_display(bare):
    if bare in NO_PREFIX_GENES or bare.startswith("KIR"):
        return bare
    return f"HLA-{bare}"


def normalize_allele_name(token, fallback_gene_bare=None):
    """'HLA-A*01:01:01:01' / 'A*01:01:01:01' / '01:01:01:01' (+file-stem gene) ->
    ('A', 'HLA-A*01:01:01:01'). Returns (None, None) if no gene can be determined."""
    if token is None:
        return None, None
    tok = str(token).strip().lstrip(">")
    if not tok:
        return None, None
    if "*" in tok:
        g, _, fields = tok.partition("*")
        g = gene_bare(g.split("|")[-1].split(":")[-1].strip()) or fallback_gene_bare
    else:
        g, fields = fallback_gene_bare, tok
    if not g:
        return None, None
    return g, f"{gene_display(g)}*{fields}"


def n_field_name(name, n):
    """First n colon fields of a normalized allele name (fewer if the name has fewer)."""
    if not isinstance(name, str) or "*" not in name:
        return None
    g, _, f = name.partition("*")
    fields = [x for x in f.split(":") if x]
    if not fields or "new" in [x.lower() for x in fields[:n]]:
        return None
    return f"{g}*{':'.join(fields[:n])}"


def sha8(s):
    return hashlib.sha1(s.encode("ascii")).hexdigest()[:8]


def open_maybe_gz(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


_UNSUPPRESSED_KEYS = {"limit", "runtime_s", "kin_min", "strict_threshold", "n_files",
                      "n_exon_maps"}


def suppress_tree(obj, key=None):
    """Recursively apply suppress() to every int in a JSON-able structure (bools, floats and the
    parameter keys in _UNSUPPRESSED_KEYS are left alone). Used for the committed summary.json."""
    if isinstance(obj, dict):
        return {k: suppress_tree(v, k) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [suppress_tree(v, key) for v in obj]
    if isinstance(obj, (int, np.integer)) and not isinstance(obj, (bool, np.bool_)):
        return int(obj) if key in _UNSUPPRESSED_KEYS else suppress(obj)
    return obj


def df_to_md(df):
    if df is None or len(df) == 0:
        return "_(no rows)_"
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for row in df.itertuples(index=False):
        cells = []
        for v in row:
            if isinstance(v, float):
                cells.append("NA" if np.isnan(v) else f"{v:.3g}")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Translation
# ---------------------------------------------------------------------------
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: _AA[16 * i + 4 * j + k]
               for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}
STOP_CODONS = {"TAA", "TAG", "TGA"}


def translate(seq):
    seq = seq.upper()
    n = len(seq) - len(seq) % 3
    return "".join(CODON_TABLE.get(seq[i:i + 3], "X") for i in range(0, n, 3))


def cds_core(seq):
    """Strip one terminal stop codon if the sequence is in frame and ends with one."""
    seq = seq.upper()
    if len(seq) >= 3 and len(seq) % 3 == 0 and seq[-3:] in STOP_CODONS:
        return seq[:-3]
    return seq


def protein_info(seq):
    """-> dict(protein=translation minus one terminal '*', frameshift, premature_stop).
    frameshift: CDS length % 3 != 0. premature_stop: a '*' before the last translated codon."""
    aa = translate(seq)
    core = aa[:-1] if aa.endswith("*") else aa
    return {"protein": core, "frameshift": len(seq) % 3 != 0, "premature_stop": "*" in core}


# ---------------------------------------------------------------------------
# Known-allele reference (Immuannot refdata CDSseq)
# ---------------------------------------------------------------------------
def parse_ref_header(line, fallback_gene_bare):
    toks = line[1:].split() if line.startswith(">") else line.split()
    if not toks:
        return None, None
    star = [t for t in toks if "*" in t]
    return normalize_allele_name(star[0] if star else toks[0], fallback_gene_bare)


_FRAME_RE = re.compile(r"\bframe=(\d+)")


def header_frame(line):
    """`frame=N` in a CDSseq header (real format: '>HLA-A*01:01:01:01 HLA00001 frame=1 1098bp').
    Anything other than 1 means the record does not start at codon position 1."""
    m = _FRAME_RE.search(line)
    return int(m.group(1)) if m else 1


def iter_fasta(path, fallback_gene_bare=None):
    """Yield (gene_bare, allele_name, seq, frame) from a (gzipped) FASTA of known alleles."""
    g = name = None
    frame = 1
    chunks = []
    with open_maybe_gz(path) as f:
        for line in f:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield g, name, "".join(chunks).upper(), frame
                g, name = parse_ref_header(line, fallback_gene_bare)
                frame = header_frame(line)
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield g, name, "".join(chunks).upper(), frame


class RefIndex:
    """Known CDS and protein sets per bare gene name.

    cds[g][cds_core] -> sorted names; prot[g][protein] -> sorted names (only in-frame, no
    premature stop); name_cds[name] = cds_core; len_mats[g][L] = (names, proteins, uint8 matrix)."""

    def __init__(self):
        self.cds = defaultdict(lambda: defaultdict(list))
        self.prot = defaultdict(lambda: defaultdict(list))
        self.name_cds = {}
        self.name_prot = {}
        self.n_records = 0
        self.n_unparsed_headers = 0
        self.n_with_terminal_stop = 0
        self.n_nonstandard_frame = 0
        self.files = []
        self._len_mats = {}

    def add(self, g, name, seq, frame=1):
        if g is None or not seq:
            self.n_unparsed_headers += 1
            return
        self.n_records += 1
        core = cds_core(seq)
        if core != seq:
            self.n_with_terminal_stop += 1
        self.cds[g][core].append(name)
        self.name_cds[name] = core
        if frame != 1:
            # Not translatable from base 1 (partial record). Usable for exact CDS identity only.
            self.n_nonstandard_frame += 1
            return
        pi = protein_info(seq)
        if not pi["frameshift"] and not pi["premature_stop"] and pi["protein"]:
            self.prot[g][pi["protein"]].append(name)
            self.name_prot[name] = pi["protein"]

    def finalize(self):
        for d in (self.cds, self.prot):
            for g in d:
                for k in d[g]:
                    d[g][k] = sorted(set(d[g][k]))

    def genes(self):
        return set(self.cds)

    def len_mats(self, g):
        if g not in self._len_mats:
            by_len = defaultdict(list)
            for p, names in self.prot.get(g, {}).items():
                by_len[len(p)].append((names[0], p))
            mats = {}
            for L, items in by_len.items():
                items.sort()
                names = [n for n, _ in items]
                prots = [p for _, p in items]
                mat = np.frombuffer("".join(prots).encode("ascii"), dtype=np.uint8).reshape(
                    len(prots), L)
                mats[L] = (names, prots, mat)
            self._len_mats[g] = mats
        return self._len_mats[g]


def load_refdata(refdata_dir, genes_needed=None):
    files = sorted(set(glob.glob(os.path.join(refdata_dir, "**", "CDSseq", "*.fa.gz"),
                                 recursive=True)))
    if not files:
        sys.exit(f"FATAL: no CDSseq/*.fa.gz under {refdata_dir!r} (glob **/CDSseq/*.fa.gz).")
    ref = RefIndex()
    for path in files:
        stem = os.path.basename(path)
        for suf in (".fa.gz", ".fasta.gz"):
            if stem.endswith(suf):
                stem = stem[: -len(suf)]
        stem_bare = gene_bare(stem)
        if genes_needed is not None and stem_bare not in genes_needed:
            # Header genes normally equal the stem. Skip unneeded files to save time.
            continue
        ref.files.append(path)
        for g, name, seq, frame in iter_fasta(path, stem_bare):
            ref.add(g, name, seq, frame)
    ref.finalize()
    return ref


def two_field_consensus(names):
    """Most common 2-field name among names (ties -> sorted first)."""
    c = Counter(n_field_name(n, 2) for n in names if n_field_name(n, 2))
    if not c:
        return None
    top = max(c.values())
    return sorted(k for k, v in c.items() if v == top)[0]


def three_field_consensus(names):
    c = Counter(n_field_name(n, 3) for n in names if n_field_name(n, 3))
    if not c:
        return None
    top = max(c.values())
    return sorted(k for k, v in c.items() if v == top)[0]


# ---------------------------------------------------------------------------
# Sequence classification
# ---------------------------------------------------------------------------
def classify_sequence(seq, gene_b, depth, ref):
    """-> dict with cds_core, protein, frameshift, premature_stop, cds_known, protein_known,
    seq_class, mapped_prot2 (2-field known name when the protein is catalogued), mapped_cds3."""
    core = cds_core(seq)
    pi = protein_info(seq)
    out = {"cds_core": core, "protein": pi["protein"], "frameshift": pi["frameshift"],
           "premature_stop": pi["premature_stop"], "cds_known": None, "protein_known": None,
           "mapped_prot2": None, "mapped_cds3": None}
    if gene_b not in ref.genes():
        out["seq_class"] = "no_refdata"
        return out
    cds_names = ref.cds[gene_b].get(core)
    out["cds_known"] = bool(cds_names)
    nonfunctional = pi["frameshift"] or pi["premature_stop"]
    prot_names = None if nonfunctional else ref.prot[gene_b].get(pi["protein"])
    out["protein_known"] = bool(prot_names)
    if cds_names:
        out["seq_class"] = "cds_known"
        out["mapped_prot2"] = two_field_consensus(cds_names)
        out["mapped_cds3"] = three_field_consensus(cds_names)
    elif nonfunctional:
        out["seq_class"] = "frameshift_or_stop"
    elif prot_names:
        out["seq_class"] = "protein_known" if depth == 2 else "novel_cds_synonymous"
        out["mapped_prot2"] = two_field_consensus(prot_names)
    else:
        out["seq_class"] = "novel_protein"
    return out


def hamming_diffs(ref_p, obs_p):
    return [(i + 1, a, b) for i, (a, b) in enumerate(zip(ref_p, obs_p)) if a != b]


def levenshtein_align(ref_p, obs_p):
    """-> (distance, diffs) with diffs [(ref_pos_1based, ref_aa|'-', obs_aa|'-')]."""
    n, m = len(ref_p), len(obs_p)
    D = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        D[i][0] = i
    for j in range(m + 1):
        D[0][j] = j
    for i in range(1, n + 1):
        Di, Dp, r = D[i], D[i - 1], ref_p[i - 1]
        for j in range(1, m + 1):
            c = Dp[j - 1] + (r != obs_p[j - 1])
            d = Dp[j] + 1
            e = Di[j - 1] + 1
            Di[j] = c if c <= d and c <= e else (d if d <= e else e)
    diffs = []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and D[i][j] == D[i - 1][j - 1] + (ref_p[i - 1] != obs_p[j - 1]):
            if ref_p[i - 1] != obs_p[j - 1]:
                diffs.append((i, ref_p[i - 1], obs_p[j - 1]))
            i, j = i - 1, j - 1
        elif i > 0 and D[i][j] == D[i - 1][j] + 1:
            diffs.append((i, ref_p[i - 1], "-"))
            i -= 1
        else:
            diffs.append((i, "-", obs_p[j - 1]))  # insertion after ref residue i
            j -= 1
    diffs.reverse()
    return D[n][m], diffs


_NN_CACHE = {}


def nearest_known_protein(protein, gene_b, ref, hint_names=(), shortlist_n=5):
    key = (id(ref), protein, gene_b, tuple(hint_names), shortlist_n)
    if key not in _NN_CACHE:
        _NN_CACHE[key] = _nearest_known_protein(protein, gene_b, ref, hint_names, shortlist_n)
    return _NN_CACHE[key]


def _nearest_known_protein(protein, gene_b, ref, hint_names=(), shortlist_n=5):
    """-> dict(name, distance, metric, diffs). Hamming vs equal-length known proteins when any
    exist, else Levenshtein vs a shortlist (hint names + the shortlist_n closest in length)."""
    mats = ref.len_mats(gene_b)
    L = len(protein)
    if L in mats:
        names, prots, mat = mats[L]
        q = np.frombuffer(protein.encode("ascii"), dtype=np.uint8)
        d = (mat != q).sum(axis=1)
        k = int(np.argmin(d))  # rows are sorted by name, so ties go to the smallest name
        return {"name": names[k], "distance": int(d[k]), "metric": "hamming",
                "diffs": hamming_diffs(prots[k], protein)}
    cands = {}
    for h in hint_names:
        _, nm = normalize_allele_name(h, gene_b)
        if nm in ref.name_prot:
            cands[nm] = ref.name_prot[nm]
    # all proteins of the closest lengths
    flat = []
    for Lk in sorted(mats, key=lambda x: (abs(x - L), x)):
        names, prots, _m = mats[Lk]
        flat.extend(zip(names, prots))
        if len(flat) >= shortlist_n:
            break
    for nm, p in flat[:shortlist_n]:
        cands.setdefault(nm, p)
    if not cands:
        return {"name": None, "distance": None, "metric": None, "diffs": []}
    best = None
    for nm in sorted(cands):
        dist, diffs = levenshtein_align(cands[nm], protein)
        if best is None or dist < best[0]:
            best = (dist, nm, diffs)
    return {"name": best[1], "distance": best[0], "metric": "levenshtein", "diffs": best[2]}


def format_diffs(diffs):
    return ",".join(f"{p}:{r}>{a}" for p, r, a in diffs)


# ---------------------------------------------------------------------------
# Exon / groove annotation (isolated)
# ---------------------------------------------------------------------------
_RANGE_RE = re.compile(r"(\d+)\.\.(\d+)")
_TRUTHY = {"1", "true", "t", "yes", "y"}


def groove_exons_for(gene_b):
    """Peptide-groove exons: 2+3 for class I (A,B,C,E,F,G), 2 for class II alpha/beta chains
    (DRA, DRB*, DQA*, DQB*, DPA*, DPB*). None for every other gene."""
    if gene_b in {"A", "B", "C", "E", "F", "G"}:
        return {2, 3}
    if gene_b == "DRA" or re.match(r"^(DRB|DQB|DPB|DQA|DPA)\d", str(gene_b)):
        return {2}
    return None


def cds_segments(cds_field, exon_field):
    """alleles.csv.gz CDS= / exon= -> [(exon_no, coding_len), ...] in CDS order.
    Each CDS range is assigned the exon whose range contains its start. If exon= is missing or
    does not cover a range, the k-th CDS range is taken as exon k."""
    cds = sorted((int(a), int(b)) for a, b in _RANGE_RE.findall(cds_field or ""))
    exons = []
    for part in (exon_field or "").split(","):
        idx, sep, rng = part.partition(":")
        m = _RANGE_RE.search(rng)
        if sep and m and idx.strip().isdigit():
            exons.append((int(idx), int(m.group(1)), int(m.group(2))))
    segs = []
    for k, (a, b) in enumerate(cds, 1):
        ex = next((n for n, lo, hi in exons if lo <= a <= hi), k)
        segs.append((ex, b - a + 1))
    return tuple(segs)


class ExonMaps:
    """Exon structure from Immuannot's alleles.csv.gz (one tab-separated key=value line per
    genomic allele: gene, version, ID, allele, ..., goodTemplate, CDS=, exon=, ...).
    by_name[name] = segments; gene_mode[gene_b] = most common segment pattern among the gene's
    goodTemplate alleles (all alleles if none is marked truthy); genes = genes present."""

    def __init__(self):
        self.by_name = {}
        self.gene_mode = {}
        self.genes = set()
        self.version = None
        self.path = None

    def __len__(self):
        return len(self.by_name)


def load_exon_maps(path, names_needed, genes_needed=None):
    em = ExonMaps()
    em.path = path
    if not path or not os.path.exists(path):
        return em
    good, allp = defaultdict(Counter), defaultdict(Counter)
    with open_maybe_gz(path) as f:
        for line in f:
            if "allele=" not in line:
                continue
            fields = {}
            for tok in line.rstrip("\n").split("\t"):
                k, _, v = tok.partition("=")
                fields[k.strip()] = v.strip()
            if em.version is None and fields.get("version"):
                em.version = fields["version"]
            g, nm = normalize_allele_name(fields.get("allele"),
                                          gene_bare(fields.get("gene", "")) or None)
            if g is None:
                continue
            em.genes.add(g)
            if genes_needed is not None and g not in genes_needed and nm not in names_needed:
                continue
            segs = cds_segments(fields.get("CDS"), fields.get("exon"))
            if not segs:
                continue
            if nm in names_needed:
                em.by_name[nm] = segs
            allp[g][segs] += 1
            if str(fields.get("goodTemplate", "")).strip().lower() in _TRUTHY:
                good[g][segs] += 1
    for g, c in allp.items():
        src = good[g] if good[g] else c
        top = max(src.values())
        em.gene_mode[g] = sorted(k for k, v in src.items() if v == top)[0]
    return em


def _segments_to_array(segs):
    return np.concatenate([np.full(n, ex) for ex, n in segs]) if segs else np.array([], int)


def annotate_groove(diffs, gene_b, ref_name, ref_core_len, obs_core_len, exon_maps):
    """-> (flags aligned with diffs [True/False/None], method, reason).

    A codon counts as in the groove if its middle base (1-based CDS position 3p-1) lies in a
    groove exon (groove_exons_for). Methods, in order:
      refdata_exon       the reference allele's own CDS=/exon= record
      refdata_gene_mode  the gene's most common CDS segment pattern (goodTemplate alleles)
      approx_canonical   APPROX_CANONICAL_CDS (HLA-A/B/C only), used only when alleles.csv.gz
                         has no entry at all for the gene AND both CDS have the canonical length
    A segment map is used only when its total equals the reference cds_core length (+3 if it
    includes the stop codon). Otherwise every flag is None and the reason is returned."""
    groove = groove_exons_for(gene_b)
    none = [None] * len(diffs)
    if groove is None:
        return none, "NA", f"no peptide-groove exon definition for {gene_b}"
    exmap, method, reason = None, None, ""
    if exon_maps is None:
        exon_maps = ExonMaps()
    for m, segs in (("refdata_exon", exon_maps.by_name.get(ref_name)),
                    ("refdata_gene_mode", exon_maps.gene_mode.get(gene_b))):
        if segs is None:
            continue
        if sum(n for _, n in segs) in (ref_core_len, ref_core_len + 3):
            exmap, method = _segments_to_array(segs), m
            break
        reason = f"{m}: CDS segment total != reference CDS length"
    if exmap is None:
        if gene_b in exon_maps.genes:
            return none, "NA", reason or "no usable exon record"
        canon = APPROX_CANONICAL_CDS.get(gene_b)
        if canon is None:
            return none, "NA", f"gene absent from alleles.csv.gz; no canonical map for {gene_b}"
        if not (ref_core_len == canon["core_len"] == obs_core_len):
            return none, "NA", "gene absent from alleles.csv.gz; CDS length != canonical length"
        lens = canon["exon_cds_lengths"]
        segs = [(k + 1, n) for k, n in enumerate(lens)] + [(99, canon["core_len"] - sum(lens))]
        exmap, method = _segments_to_array(segs), "approx_canonical"
    flags = []
    for p, _r, _a in diffs:
        mid = 3 * p - 2  # 0-based index of the codon's middle base
        flags.append(bool(exmap[mid] in groove) if 0 <= mid < len(exmap) else None)
    return flags, method, ""


# ---------------------------------------------------------------------------
# Field class / artifact label
# ---------------------------------------------------------------------------
def depth_from_consensus(consensus):
    if not isinstance(consensus, str):
        return None
    s = consensus.strip()
    if not s or s.lower() == "undetermined":
        return None
    part = s.split("*", 1)[1] if "*" in s else s
    fields = part.split(":")
    return fields.index("new") + 1 if "new" in fields else None


def field_class_of(consensus, novelty_depth=None):
    if not isinstance(consensus, str) or consensus.strip().upper() in vc.NULL_TOKENS:
        return "uncalled"
    d = novelty_depth
    try:
        d = None if d is None or pd.isna(d) else int(float(d))
    except (TypeError, ValueError):
        d = None
    if d is None:
        d = depth_from_consensus(consensus)
    if d is None:
        return "known"
    return DEPTH_TO_CLASS.get(d, "f1_undetermined")


def artifact_label_of(template_warning, cds_mut):
    if isinstance(cds_mut, str) and is_homopolymer_indel_only(cds_mut):
        return "homopolymer_indel"
    toks = warning_tokens(template_warning)
    if "partial_CDS" in toks:
        return "partial_cds"
    if "inframe_stop" in toks:
        return "inframe_stop"
    return "clean"


def add_call_labels(df):
    """Vectorized over unique values: field_class, depth, artifact_label."""
    depth_col = df["novelty_depth"] if "novelty_depth" in df.columns else pd.Series(
        np.nan, index=df.index)
    depth_num = pd.to_numeric(depth_col, errors="coerce")
    cons = df["consensus"].astype(object)
    miss = depth_num.isna()
    if miss.any():
        derived = cons[miss].map(depth_from_consensus)
        depth_num = depth_num.copy()
        depth_num[miss] = pd.to_numeric(derived, errors="coerce")
    df["depth"] = depth_num
    uncalled = cons.map(lambda c: not isinstance(c, str) or c.strip().upper() in vc.NULL_TOKENS)
    fc = depth_num.map(lambda d: DEPTH_TO_CLASS.get(int(d), "f1_undetermined")
                       if pd.notna(d) else "known")
    fc[uncalled.astype(bool)] = "uncalled"
    df["field_class"] = fc.astype(str)
    tw = df["template_warning"] if "template_warning" in df.columns else pd.Series(
        None, index=df.index, dtype=object)
    cm = df["cds_mut"] if "cds_mut" in df.columns else pd.Series(None, index=df.index,
                                                                 dtype=object)
    tw_label = {}
    for v in pd.unique(tw.astype(object)):
        toks = warning_tokens(v)
        tw_label[v if isinstance(v, str) else None] = (
            "partial_cds" if "partial_CDS" in toks else
            "inframe_stop" if "inframe_stop" in toks else "clean")
    hp_vals = {v: is_homopolymer_indel_only(v) for v in pd.unique(cm.dropna().astype(object))}
    hp = cm.astype(object).map(lambda v: hp_vals.get(v, False) if isinstance(v, str) else False)
    lab = tw.astype(object).map(lambda v: tw_label.get(v if isinstance(v, str) else None,
                                                       "clean"))
    lab[hp.astype(bool)] = "homopolymer_indel"
    df["artifact_label"] = lab.astype(str)
    return df


# ---------------------------------------------------------------------------
# People: ancestry + unrelated set
# ---------------------------------------------------------------------------
def strict_ancestry(people, threshold=STRICT_MIN):
    """ancestry_pred if p_<anc> >= threshold, else NA."""
    out = []
    for anc, row in zip(people["ancestry_pred"], people.to_dict("records")):
        if not isinstance(anc, str) or anc in ("", "NA"):
            out.append(None)
            continue
        p = row.get(f"p_{anc.lower()}")
        try:
            ok = p is not None and not pd.isna(p) and float(p) >= threshold
        except (TypeError, ValueError):
            ok = False
        out.append(anc if ok else None)
    return pd.Series(out, index=people.index, dtype=object)


def greedy_unrelated(person_ids, pairs, kin_min=KIN_MIN):
    """Remove people until no pair with kin >= kin_min remains, among person_ids only.
    Each step drops the person with the most remaining relatives; ties go to the smallest id
    (string sort). pairs: iterable of (i, j, kin). -> (kept set, removed list in order)."""
    ids = set(map(str, person_ids))
    adj = defaultdict(set)
    for i, j, k in pairs:
        i, j = str(i), str(j)
        try:
            k = float(k)
        except (TypeError, ValueError):
            continue
        if i == j or k < kin_min or i not in ids or j not in ids:
            continue
        adj[i].add(j)
        adj[j].add(i)
    removed = []
    while True:
        live = [(len(v), n) for n, v in adj.items() if v]
        if not live:
            break
        maxdeg = max(d for d, _ in live)
        victim = min(n for d, n in live if d == maxdeg)
        for nb in adj[victim]:
            adj[nb].discard(victim)
        adj[victim] = set()
        removed.append(victim)
    return ids - set(removed), removed


def load_relatedness_pairs(path):
    if not os.path.exists(path):
        sys.exit(f"FATAL: relatedness table not found at {path!r}. Is the controlled-tier bucket "
                 f"mounted (RUNBOOK.md)? Use --skip-relatedness only for a local dry run.")
    df = pd.read_csv(path, sep="\t", dtype={"i.s": str, "j.s": str})
    missing = {"i.s", "j.s", "kin"} - set(df.columns)
    if missing:
        sys.exit(f"FATAL: relatedness table missing {missing}; columns: {list(df.columns)}")
    return list(zip(df["i.s"], df["j.s"], df["kin"]))


def build_people(person_ids, cohort_path, relatedness_path, kin_min, strict_min, skip_rel):
    cohort = vc.load_cohort_membership(cohort_path)
    pcols = [f"p_{a.lower()}" for a in vc.ANCESTRY_ORDER]
    missing = [c for c in ["ancestry_pred"] + pcols if c not in cohort.columns]
    if missing:
        sys.exit(f"FATAL: {cohort_path} lacks {missing} (needed for strict ancestry). "
                 f"Columns: {list(cohort.columns)}")
    people = pd.DataFrame({"person_id": sorted(set(person_ids))})
    people = people.merge(cohort[["person_id", "ancestry_pred"] + pcols], on="person_id",
                          how="left")
    people["anc_pred"] = people["ancestry_pred"].astype(object).where(
        people["ancestry_pred"].notna(), None)
    people["anc_strict"] = strict_ancestry(people, strict_min)
    for c in ("anc_pred", "anc_strict"):
        people[c] = people[c].map(lambda v: v if isinstance(v, str) and v else "UNASSIGNED")
    if skip_rel:
        people["unrelated"] = True
        n_removed = 0
    else:
        lr_ids = set(people["person_id"])  # Table 1 people ARE the long-read cohort
        kept, removed = greedy_unrelated(lr_ids, load_relatedness_pairs(relatedness_path),
                                         kin_min)
        people["unrelated"] = people["person_id"].isin(kept)
        n_removed = len(removed)
    return people[["person_id", "anc_pred", "anc_strict", "unrelated"]], n_removed


# ---------------------------------------------------------------------------
# Table 1 loading
# ---------------------------------------------------------------------------
def load_table1(path, limit=None, chunksize=400_000):
    """Read Table 1 with usecols + category dtypes, in chunks (the real file is ~1.6M rows on a
    31 GB VM). `alleles`/`template_allele` are only needed on novel rows, so they are dropped
    elsewhere: they are long comma-joined strings on every row otherwise."""
    if not os.path.exists(path):
        sys.exit(f"FATAL: Table 1 not found at {path!r}. Run 01_extract_rich.py first.")
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    use = [c for c in TABLE1_COLS if c in header]
    req = {"person_id", "hap", "contig", "gene", "consensus", "is_novel"}
    if req - set(use):
        sys.exit(f"FATAL: Table 1 lacks required columns {req - set(use)}")
    dtypes = {c: object for c in use}
    for c in ("copy_index", "novelty_depth", "cds_distance", "n_aa_changes"):
        if c in use:
            dtypes[c] = "float64"
    for c in ("gene_class", "novelty_class", "template_warning"):
        if c in use:
            dtypes[c] = "category"
    keep = None
    parts = []
    for chunk in pd.read_csv(path, sep="\t", usecols=use, dtype=dtypes, chunksize=chunksize):
        chunk["person_id"] = chunk["person_id"].astype(str)
        if limit:
            if keep is None:
                keep = set()
            for pid in chunk["person_id"].unique():
                if len(keep) < limit or pid in keep:
                    keep.add(pid)
            keep = set(sorted(keep)[:limit])
            chunk = chunk[chunk["person_id"].isin(keep)]
            if chunk.empty:
                continue
        novel = chunk["is_novel"].astype(str).str.lower().isin(["true", "1"])
        for c in ("alleles", "template_allele"):
            if c in chunk.columns:
                chunk.loc[~novel, c] = None
        parts.append(chunk)
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=use)
    if limit:
        df = df[df["person_id"].isin(set(sorted(df["person_id"].unique())[:limit]))]
        df = df.reset_index(drop=True)
    for c in ("contig", "gene", "hap"):
        df[c] = df[c].fillna("").astype(str)
    if "copy_index" not in df.columns:
        df["copy_index"] = 1.0
    df["copy_index"] = df["copy_index"].fillna(1.0)
    if "gene_class" not in df.columns:
        df["gene_class"] = "other"
    df["gene_class"] = df["gene_class"].astype(object).fillna("other").astype(str)
    return df


# ---------------------------------------------------------------------------
# Sequence matching (threaded, per person)
# ---------------------------------------------------------------------------
def _match_person(sub, outroot):
    sub = sub.copy()
    sub["is_novel"] = sub["depth"].isin([2, 3]).map({True: "True", False: "False"})
    rows, seqs, stats = match_novel_rows(sub, outroot)
    keys_all = Counter(zip(sub["person_id"], sub["hap"], sub["gene"]))
    depth_by_key = {}
    for k, d, ci in zip(zip(sub["person_id"], sub["hap"], sub["gene"]), sub["depth"],
                        sub["copy_index"]):
        depth_by_key[k] = (d, ci)
    out = []
    n_join_ambig = 0
    for r in rows:
        k = (r["person_id"], r["hap"], r["gene"])
        if keys_all[k] != 1:
            n_join_ambig += 1
            continue
        d, ci = depth_by_key[k]
        if ci != 1:
            n_join_ambig += 1
            continue
        out.append({"person_id": k[0], "hap": k[1], "gene": k[2], "depth": int(d),
                    "seq": seqs[r["cds_seq_sha1"]], "nearest_allele": r["nearest_allele"]})
    stats["n_join_ambiguous"] = n_join_ambig
    return out, stats


def match_sequences(t1, outroot, threads):
    target = t1.loc[t1["depth"].isin([2, 3]), "person_id"].unique()
    target_set = set(target)
    groups = {pid: g for pid, g in t1[t1["person_id"].isin(target_set)].groupby("person_id")}
    results, agg = [], Counter()
    t0 = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        futs = [ex.submit(_match_person, g, outroot) for g in groups.values()]
        for n, fut in enumerate(concurrent.futures.as_completed(futs), 1):
            rows, stats = fut.result()
            results.extend(rows)
            agg.update(stats)
            if n % 500 == 0 or n == len(futs):
                print(f"  matched {n}/{len(futs)} people ({time.time() - t0:.0f}s)",
                      file=sys.stderr)
    return results, dict(agg)


# ---------------------------------------------------------------------------
# Stratified counting
# ---------------------------------------------------------------------------
def stratified_counts(df, keys, value_name="n_haplotypes"):
    """Count rows by keys x (ancestry_scheme, ancestry incl. POOLED) x unrelated_only."""
    out = []
    for unrel in (False, True):
        d = df[df["unrelated"]] if unrel else df
        for scheme, col in (("pred", "anc_pred"), ("strict", "anc_strict")):
            g = d.groupby(keys + [col], observed=True).size().reset_index(name=value_name)
            g = g.rename(columns={col: "ancestry"})
            p = d.groupby(keys, observed=True).size().reset_index(name=value_name)
            p["ancestry"] = "POOLED"
            for x in (g, p):
                x["ancestry_scheme"] = scheme
                x["unrelated_only"] = unrel
            out.extend([g, p])
    res = pd.concat(out, ignore_index=True)
    return res[keys[:1] + [c for c in ["ancestry_scheme", "ancestry", "unrelated_only"]]
               + keys[1:] + [value_name]]


# ---------------------------------------------------------------------------
# Clusters
# ---------------------------------------------------------------------------
def build_clusters(seq_calls, ref, exon_maps, gene_class_by_gene, drift=None):
    """seq_calls: call-level frame (novel_protein / novel_cds_synonymous / protein_known) with
    person_id, hap, gene, gene_b, anc_pred, anc_strict, unrelated, artifact_label, seq_class,
    cds_core, protein, depth, nearest_allele. -> DataFrame of clusters."""
    rows = []
    specs = [
        ("novel_protein", seq_calls[seq_calls["seq_class"] == "novel_protein"], "protein"),
        ("novel_cds_synonymous",
         seq_calls[seq_calls["seq_class"].isin(["novel_cds_synonymous", "protein_known"])],
         "cds_core"),
    ]
    ancs = vc.ANCESTRY_ORDER + ["UNASSIGNED"]
    for ctype, sub, keycol in specs:
        for (gb, key), m in sub.groupby(["gene_b", keycol], sort=True):
            gene = m["gene"].iloc[0]
            gclass = gene_class_by_gene.get(gene, "other")
            prot = m["protein"].iloc[0]
            core = m["cds_core"].iloc[0]
            cid = (f"{gb}_prot_{sha8(key)}" if ctype == "novel_protein"
                   else f"{gb}_cds_{sha8(key)}")
            labels = set(m["artifact_label"])
            clean = m[m["artifact_label"] == "clean"]
            rec = {
                "cluster_id": cid, "cluster_type": ctype, "gene": gene, "gene_class": gclass,
                "protein_len": len(prot),
                "cds_len": int(m["cds_core"].str.len().mode().iloc[0]),
                "n_cds_variants": int(m["cds_core"].nunique()),
                "n_haplotypes": int(len(m)),
                "n_persons": int(m["person_id"].nunique()),
                "n_persons_unrelated": int(m.loc[m["unrelated"], "person_id"].nunique()),
                "n_persons_clean": int(clean["person_id"].nunique()),
                "n_persons_unrelated_clean": int(
                    clean.loc[clean["unrelated"], "person_id"].nunique()),
                "n_haplotypes_named_protein_new": int((m["depth"] == 2).sum()),
                "artifact_labels_any": ",".join(sorted(labels)),
                "all_members_flagged": bool("clean" not in labels),
                "all_members_clean": bool(labels == {"clean"}),
            }
            for scheme in ("pred", "strict"):
                col = f"anc_{scheme}"
                cnt = m.groupby(col)["person_id"].nunique()
                for a in ancs:
                    rec[f"n_persons_{scheme}_{a}"] = int(cnt.get(a, 0))
            if ctype == "novel_protein":
                hints = [h for h in m["nearest_allele"].dropna().astype(str)]
                hint = Counter(hints).most_common(1)[0][0] if hints else None
                nn = nearest_known_protein(prot, gb, ref, [hint] if hint else [])
                diffs = nn["diffs"]
                ref_core = ref.name_cds.get(nn["name"], "")
                flags, gmethod, greason = annotate_groove(
                    diffs, gb, nn["name"], len(ref_core), len(core), exon_maps)
                rec.update({
                    "nearest_known_allele": nn["name"],
                    "nearest_known_protein": n_field_name(nn["name"], 2) if nn["name"] else None,
                    "nearest_distance": nn["distance"], "distance_metric": nn["metric"],
                    "n_aa_diffs": len(diffs), "aa_diffs": format_diffs(diffs),
                    "in_groove": ",".join("NA" if f is None else str(f) for f in flags),
                    "n_groove_diffs": (int(sum(1 for f in flags if f)) if gmethod != "NA"
                                       else None),
                    "groove_method": gmethod, "groove_na_reason": greason,
                })
            else:
                names = ref.prot[gb].get(prot, [])
                rec.update({
                    "nearest_known_allele": names[0] if names else None,
                    "nearest_known_protein": two_field_consensus(names) if names else None,
                    "nearest_distance": 0, "distance_metric": "identical_protein",
                    "n_aa_diffs": 0, "aa_diffs": "", "in_groove": "", "n_groove_diffs": 0,
                    "groove_method": "not_applicable", "groove_na_reason": "",
                })
            if drift is not None:
                dc, dp = drift.get(gb, (set(), set()))
                rec["in_imgt_latest_cds"] = bool(set(m["cds_core"]) & dc)
                rec["in_imgt_latest_protein"] = prot in dp
            rec["clean_recurrent_unrelated"] = rec["n_persons_unrelated_clean"] >= 2
            rows.append(rec)
    return pd.DataFrame(rows)


CLUSTER_COUNT_COLS_PREFIX = ("n_haplotypes", "n_persons")


def load_drift(path):
    drift = defaultdict(lambda: (set(), set()))
    n = 0
    for g, _name, seq, frame in iter_fasta(path, None):
        if g is None or frame != 1:
            continue
        n += 1
        cds_set, prot_set = drift[g]
        cds_set.add(cds_core(seq))
        pi = protein_info(seq)
        if not pi["frameshift"] and not pi["premature_stop"]:
            prot_set.add(pi["protein"])
    return dict(drift), n


# ---------------------------------------------------------------------------
# Allele-count tables (WS3)
# ---------------------------------------------------------------------------
def allele_ids(calls):
    """Adds prot_id / cds_id / keep_clean to the call frame (NaN = dropped)."""
    fc = calls["field_class"]
    sc = calls["seq_class"]
    cons_norm = calls["consensus"].map(lambda c: normalize_allele_name(c)[1]
                                       if isinstance(c, str) else None)
    prot = pd.Series(None, index=calls.index, dtype=object)
    cds = pd.Series(None, index=calls.index, dtype=object)
    known = fc.isin(["known", "f4_noncoding"])
    prot[known] = cons_norm[known].map(lambda n: n_field_name(n, 2))
    cds[known] = cons_norm[known].map(lambda n: n_field_name(n, 3))
    gb = calls["gene_b"]
    core_id = calls["cds_core"].map(lambda s: sha8(s) if isinstance(s, str) else None)
    prot_sha = calls["protein"].map(lambda s: sha8(s) if isinstance(s, str) else None)
    m = sc == "cds_known"
    prot[m] = calls.loc[m, "mapped_prot2"]
    cds[m] = calls.loc[m, "mapped_cds3"]
    m = sc.isin(["protein_known", "novel_cds_synonymous"])
    prot[m] = calls.loc[m, "mapped_prot2"]
    cds[m] = gb[m] + "_cds_" + core_id[m]
    m = sc == "novel_protein"
    prot[m] = gb[m] + "_prot_" + prot_sha[m]
    cds[m] = gb[m] + "_cds_" + core_id[m]
    m = sc == "frameshift_or_stop"
    prot[m] = gb[m] + "_nonfunc_" + core_id[m]
    cds[m] = gb[m] + "_nonfunc_" + core_id[m]
    calls["prot_id"] = prot
    calls["cds_id"] = cds
    calls["keep_clean"] = (calls["artifact_label"] == "clean") & (sc != "frameshift_or_stop")
    return calls


def allele_count_table(calls):
    parts, dropped = [], {}
    for res, col in (("protein", "prot_id"), ("cds", "cds_id")):
        for excl in ("all", "clean_only"):
            d = calls if excl == "all" else calls[calls["keep_clean"]]
            ok = d[col].notna()
            dropped[f"{res}_{excl}"] = int((~ok).sum())
            d = d.loc[ok, ["gene", col, "anc_pred", "anc_strict", "unrelated"]].rename(
                columns={col: "allele_id"})
            t = stratified_counts(d, ["gene", "allele_id"])
            t.insert(1, "resolution", res)
            t.insert(2, "exclude_flagged", excl)
            parts.append(t)
    return pd.concat(parts, ignore_index=True), dropped


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _gene_order(df):
    g = df[["gene", "gene_class"]].drop_duplicates()
    g["o"] = g["gene_class"].map(lambda c: GENE_CLASS_ORDER.index(c)
                                 if c in GENE_CLASS_ORDER else 99)
    return list(g.sort_values(["o", "gene"])["gene"])


def fig_field_class_by_gene(calls, path):
    genes = _gene_order(calls)
    tab = calls.groupby(["gene", "field_class", "artifact_label"]).size()
    tot = calls.groupby("gene").size()
    fig, ax = plt.subplots(figsize=(7.2, max(3.0, 0.22 * len(genes) + 1.2)))
    y = np.arange(len(genes))[::-1]
    left = np.zeros(len(genes))
    for fcl in FIELD_CLASSES:
        for flagged in (False, True):
            vals = []
            for g in genes:
                n = 0
                for lab in ARTIFACT_LABELS:
                    if (lab != "clean") == flagged:
                        n += tab.get((g, fcl, lab), 0)
                vals.append(n / tot[g] if tot[g] else 0)
            vals = np.array(vals)
            ax.barh(y, vals, left=left, color=FIELD_COLORS[fcl],
                    hatch="////" if flagged else None, edgecolor="white" if not flagged
                    else "black", linewidth=0.0 if not flagged else 0.3, height=0.8)
            left += vals
    ax.set_yticks(y)
    ax.set_yticklabels([f"{g} (n={vc_n(tot[g])})" for g in genes], fontsize=7)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of haplotype-gene calls (copy 1)", fontsize=8)
    ax.tick_params(axis="x", labelsize=7)
    handles = [Patch(facecolor=FIELD_COLORS[f], label=f) for f in FIELD_CLASSES]
    handles.append(Patch(facecolor="white", edgecolor="black", hatch="////",
                         label="artifact-flagged share"))
    ax.legend(handles=handles, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0),
              frameon=False)
    ax.set_title("Novelty field class per gene (all people, pred ancestry pooled)", fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    vc.savefig(fig, path, dpi=200)


def vc_n(n):
    return suppress(n)


def fig_rate_by_ancestry(calls, path, min_denom=20):
    d = calls[(calls["field_class"] != "uncalled") & calls["unrelated"]
              & calls["gene_class"].isin(["classical_I", "classical_II"])]
    genes = [g for g in vc.CLASSICAL_GENES if g in set(d["gene"])]
    groups = genes + ["classical (all 8)"]
    ancs = [a for a in vc.ANCESTRY_ORDER if a in set(d["anc_strict"])]
    classes = ["f2_protein", "f3_synonymous", "f4_noncoding"]
    fig, axes = plt.subplots(len(classes), 1, figsize=(7.2, 6.4), sharex=True)
    width = 0.8 / max(len(ancs), 1)
    for ax, fcl in zip(np.atleast_1d(axes), classes):
        ymax = 0
        for ai, a in enumerate(ancs):
            da = d[d["anc_strict"] == a]
            xs, ps, los, his = [], [], [], []
            for gi, g in enumerate(groups):
                dg = da if g.startswith("classical (") else da[da["gene"] == g]
                n = len(dg)
                k = int(((dg["field_class"] == fcl) & (dg["artifact_label"] == "clean")).sum())
                if n < min_denom:
                    continue
                p, lo, hi = vc.wilson_ci(k, n)
                xs.append(gi + (ai - (len(ancs) - 1) / 2) * width)
                ps.append(100 * float(p))
                # Wilson is centred on (p + z^2/2n)/denom, NOT on p, so for small k the lower
                # bound can sit above p; matplotlib rejects negative yerr. Clamp at 0.
                los.append(max(0.0, 100 * float(p - lo)))
                his.append(max(0.0, 100 * float(hi - p)))
            if xs:
                ax.bar(xs, ps, width=width, color=vc.ANCESTRY_COLORS[a],
                       label=a if fcl == classes[0] else None)
                ax.errorbar(xs, ps, yerr=[los, his], fmt="none", ecolor="black",
                            elinewidth=0.5, capsize=1)
                ymax = max(ymax, max(np.array(ps) + np.array(his)))
        ax.set_ylabel(f"clean {fcl}\nper 100 hap.", fontsize=7)
        ax.tick_params(labelsize=7)
        ax.set_ylim(0, ymax * 1.1 if ymax > 0 else 1)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[-1].set_xticks(range(len(groups)))
    axes[-1].set_xticklabels(groups, rotation=30, ha="right", fontsize=7)
    axes[0].legend(fontsize=7, ncol=len(ancs), frameon=False, loc="upper left")
    axes[0].set_title(f"Clean novelty rate by strict ancestry (p>=0.9), unrelated people; "
                      f"cells with <{min_denom} haplotypes omitted; 95% Wilson CI", fontsize=8)
    vc.savefig(fig, path, dpi=200)


def funnel_counts(calls, clusters):
    out = {}
    for label, d in (("all_genes", calls), ("classical",
                                           calls[calls["gene_class"].isin(
                                               ["classical_I", "classical_II"])])):
        nov = d[d["field_class"].isin(["f2_protein", "f3_synonymous"])]
        matched = nov[nov["seq_class"].notna() & (nov["seq_class"] != "no_refdata")]
        clean = matched[matched["artifact_label"] == "clean"]
        not_imgt = clean[clean["seq_class"] != "cds_known"]
        nprot = not_imgt[not_imgt["seq_class"] == "novel_protein"]
        rec_ids = set()
        if clusters is not None and len(clusters):
            cc = clusters[(clusters["cluster_type"] == "novel_protein")
                          & clusters["clean_recurrent_unrelated"]]
            rec_ids = set(cc["cluster_id"])
        rec = nprot[nprot["prot_id"].isin(rec_ids)]
        out[label] = [
            ("depth-2/3 novel calls", len(nov)),
            ("matched observed CDS", len(matched)),
            ("no artifact label", len(clean)),
            ("CDS not in IMGT snapshot", len(not_imgt)),
            ("novel protein", len(nprot)),
            ("in clean protein cluster\nrecurrent in >=2 unrelated", len(rec)),
        ]
    return out


def fig_funnel(fc, path):
    """Percentages of the first step, not raw counts: a committed figure must not disclose a
    small cell through its axis (the printed labels go through suppress())."""
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)
    for ax, (label, steps) in zip(axes, fc.items()):
        names = [s for s, _ in steps][::-1]
        vals = [v for _, v in steps][::-1]
        top = max(vals[-1], 1)
        pct = [100 * v / top for v in vals]
        ax.barh(range(len(vals)), pct, color="#0072B2", height=0.7)
        for i, (v, q) in enumerate(zip(vals, pct)):
            ax.text(q + 1.5, i, f"{suppress(v)} ({q:.1f}%)", va="center", fontsize=7)
        ax.set_yticks(range(len(vals)))
        ax.set_yticklabels(names, fontsize=7)
        ax.set_xlim(0, 145)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_title(label.replace("_", " "), fontsize=8)
        ax.tick_params(axis="x", labelsize=7)
        ax.set_xlabel("% of depth-2/3 novel calls", fontsize=7)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.suptitle("From 'new' in the name to a recurrent novel protein", fontsize=9)
    fig.tight_layout()
    vc.savefig(fig, path, dpi=200)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def write_report(path, S, field_pooled, seq_tab, clus_summary, groove_tab, drift_line,
                 top_clusters, limit):
    L = ["# Novelty by nomenclature field (script 24)", ""]
    if limit:
        L += [f"**PILOT RUN: first {limit} people only. Not a result.**", ""]
    L += [
        "## What this measures, in plain language",
        "",
        "Immuannot writes `new` in place of the first allele field it cannot match. The position "
        "of `new` says what is new. Depth 2 (`HLA-A*01:new`) means a new protein. Depth 3 "
        "(`HLA-A*01:01:new`) means a new coding sequence that encodes a known protein. Depth 4 "
        "(`HLA-A*01:01:01:new`) means the coding sequence is catalogued and only introns or UTRs "
        "differ. Depth 1 means even the first field is unresolved. (The WS1 brief's background "
        "table shows these strings one field shallower. This script follows SCHEMA.md and "
        "01_extract_rich.py.)",
        "",
        "Each haplotype-gene call is counted once (copy 1). Artifact flags are *labels*, not "
        "deletions. `homopolymer_indel` means every coding difference is a single-base-run "
        "indel, the typical HiFi error. `partial_cds` and `inframe_stop` are Immuannot's "
        "reconstruction warnings. A flagged call is not proven wrong, so both views are shown.",
        "",
        "For depth-2/3 calls the observed CDS was then compared with the IPD-IMGT snapshot that "
        "Immuannot ships. Some 'new' names turn out to be catalogued sequences (`cds_known`) or "
        "catalogued proteins (`protein_known`). Frameshifts and premature stops are separated "
        "out before any protein is called novel.",
        "",
        "Counts of 1-19 are shown as `<20` (AoU policy). Relatives: "
        f"{suppress(S['n_people_removed_for_relatedness'])} people were removed greedily for "
        f"kin >= {S['kin_min']}, leaving {suppress(S['n_people_unrelated'])} unrelated people "
        f"out of {suppress(S['n_people'])}."
        + (" **Relatedness was skipped (`--skip-relatedness`): everyone is treated as "
           "unrelated.**" if S["relatedness_skipped"] else ""),
        "",
        "## Headline numbers",
        "",
        df_to_md(pd.DataFrame(
            [(k, suppress(v) if isinstance(v, (int, np.integer)) and not isinstance(v, bool)
              else v) for k, v in S["headline"].items()], columns=["metric", "value"])),
        "",
        "## A. Field class per gene (all people, pooled; fraction of calls; flagged share in "
        "brackets)",
        "",
        df_to_md(field_pooled),
        "",
        "## B. Sequence truth check for depth-2/3 calls",
        "",
        "Rows are sequence classes per gene class. Counts are haplotype-gene calls.",
        "",
        df_to_md(seq_tab),
        "",
        f"Depth-3 calls whose protein is **not** catalogued (reclassified to `novel_protein`): "
        f"{suppress(S['n_depth3_protein_not_catalogued'])}. "
        f"Depth-2 calls whose protein **is** catalogued (`protein_known`): "
        f"{suppress(S['n_depth2_protein_catalogued'])}.",
        "",
        "## Clusters",
        "",
        "Novel proteins are clustered by exact protein sequence, and synonymous novelty by exact "
        "CDS. `clean recurrent` means at least 2 unrelated people carry the cluster in a call "
        "with no artifact label.",
        "",
        df_to_md(clus_summary),
        "",
        "### Largest clean novel-protein clusters",
        "",
        df_to_md(top_clusters),
        "",
        "## Groove annotation",
        "",
        "A codon counts as in the groove if its middle base falls in exon 2/3 (class I: A, B, "
        "C, E, F, G) or exon 2 (class II chains: DRA, DRB*, DQA*, DQB*, DPA*, DPB*). "
        "`refdata_exon` uses the nearest allele's own CDS=/exon= record from `alleles.csv.gz`. "
        "`refdata_gene_mode` uses that gene's most common coding-exon pattern across its "
        "goodTemplate alleles. `approx_canonical` is a hard-coded HLA-A/B/C map (73/270/276 nt) "
        "and is used only when the gene is missing from alleles.csv.gz. `NA` means none applied, "
        "with the reason below.",
        "",
        df_to_md(groove_tab),
        "",
        "## Release drift",
        "",
        drift_line,
        "",
        "## Matching and dropped rows",
        "",
        df_to_md(pd.DataFrame([(k, suppress(v)) for k, v in sorted(S["match_stats"].items())],
                              columns=["metric", "count"])),
        "",
        "Rows left out of the WS3 allele-count tables (uncalled, depth 1, unmatched, "
        "`no_refdata`; `clean_only` also drops flagged calls):",
        "",
        df_to_md(pd.DataFrame([(k, suppress(v)) for k, v in
                               sorted(S["allele_table_dropped"].items())],
                              columns=["table", "n_calls_dropped"])),
        "",
        "## Assumptions (verify on the VM)",
        "",
        f"- Refdata: {S['refdata']['n_files']} CDSseq files, {S['refdata']['n_records']} records "
        f"({S['refdata']['n_with_terminal_stop']} ending in a stop codon, "
        f"{S['refdata']['n_nonstandard_frame']} with frame!=1 kept out of the protein index, "
        f"{S['refdata']['n_unparsed_headers']} unparseable headers). Headers were parsed as the "
        "first `*`-bearing token, with an optional `HLA-` prefix.",
        f"- IPD-IMGT/HLA version (alleles.csv.gz `version=`): "
        f"{S.get('imgt_version_alleles_csv') or 'unknown'}.",
        f"- Exon maps from alleles.csv.gz: {S['n_exon_maps']} nearest reference alleles "
        f"resolved by name; gene-mode patterns for "
        f"{', '.join(S.get('exon_map_genes_with_mode') or []) or 'none'}.",
        "- Relatedness ids are assumed to equal Table 1 person_id strings.",
        "",
    ]
    vc.write_report(path, L)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=DEFAULT_TABLE1)
    ap.add_argument("--cohort-membership", default=DEFAULT_COHORT)
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT,
                    help="Root of <person_id>/immuannot_output/hap{1,2}/cds.fa.gz")
    ap.add_argument("--refdata", default=DEFAULT_REFDATA,
                    help="Immuannot refdata dir (glob **/CDSseq/*.fa.gz)")
    ap.add_argument("--alleles-csv", default=None,
                    help="alleles.csv.gz with exon maps (default: <refdata>/alleles.csv.gz, "
                         "else the first **/alleles.csv.gz)")
    ap.add_argument("--relatedness-table", default=DEFAULT_RELATEDNESS)
    ap.add_argument("--skip-relatedness", action="store_true",
                    help="Treat everyone as unrelated (local dry runs only; flagged in report)")
    ap.add_argument("--kin-min", type=float, default=KIN_MIN)
    ap.add_argument("--strict-threshold", type=float, default=STRICT_MIN)
    ap.add_argument("--imgt-latest-nuc", default=None,
                    help="Optional newer hla_nuc.fasta for the release-drift check")
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR,
                    help="Committed, aggregate-only outputs")
    ap.add_argument("--local-dir", default=DEFAULT_LOCAL_DIR,
                    help="VM-local outputs (never commit)")
    ap.add_argument("--threads", type=int, default=8,
                    help="Per-person I/O threads (VM has 4 vCPUs; I/O-bound, so 8 is fine)")
    ap.add_argument("--limit", type=int, default=None,
                    help="Pilot: first N people (sorted person_id); outputs go to pilot_limitN/")
    ap.add_argument("--min-denominator", type=int, default=20,
                    help="Minimum haplotypes for a rate bar in the ancestry figure")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    t0 = time.time()
    out_dir = os.path.expanduser(args.out_dir)
    local_dir = os.path.expanduser(args.local_dir)
    if args.limit:
        out_dir = os.path.join(out_dir, f"pilot_limit{args.limit}")
        local_dir = os.path.join(local_dir, f"pilot_limit{args.limit}")
    vc.ensure_dir(out_dir)
    vc.ensure_dir(local_dir)

    print(f"Loading Table 1 {args.table1} ...", file=sys.stderr)
    t1 = load_table1(args.table1, args.limit)
    t1 = add_call_labels(t1)
    t1["gene_b"] = t1["gene"].map(gene_bare)
    print(f"  {len(t1)} rows, {t1['person_id'].nunique()} people", file=sys.stderr)

    people, n_removed = build_people(t1["person_id"].unique(), args.cohort_membership,
                                     args.relatedness_table, args.kin_min,
                                     args.strict_threshold, args.skip_relatedness)
    people.to_csv(os.path.join(local_dir, "unrelated_set.tsv"), sep="\t", index=False)

    # Calls: copy 1, one row per (person, hap, gene)
    c1 = t1[t1["copy_index"] == 1].sort_values(["person_id", "hap", "gene", "contig"])
    n_dup = int(c1.duplicated(["person_id", "hap", "gene"]).sum())
    calls = c1.drop_duplicates(["person_id", "hap", "gene"]).merge(people, on="person_id",
                                                                   how="left")
    calls["unrelated"] = calls["unrelated"].fillna(False).astype(bool)
    for c in ("anc_pred", "anc_strict"):
        calls[c] = calls[c].fillna("UNASSIGNED")

    # Sequence check
    genes_needed = set(t1.loc[t1["depth"].isin([2, 3]), "gene_b"])
    print(f"Loading refdata ({len(genes_needed)} genes) ...", file=sys.stderr)
    ref = load_refdata(os.path.expanduser(args.refdata), genes_needed)
    print(f"  {ref.n_records} known CDS records from {len(ref.files)} files", file=sys.stderr)
    print("Matching depth-2/3 calls to cds.fa.gz ...", file=sys.stderr)
    matched, mstats = match_sequences(t1, os.path.expanduser(args.outroot), args.threads)
    mstats["n_copy1_duplicate_keys_dropped"] = n_dup
    cache = {}
    seq_rows = []
    for r in matched:
        gb = gene_bare(r["gene"])
        key = (gb, r["seq"], r["depth"])
        if key not in cache:
            cache[key] = classify_sequence(r["seq"], gb, r["depth"], ref)
        info = cache[key]
        seq_rows.append({"person_id": r["person_id"], "hap": r["hap"], "gene": r["gene"],
                         "nearest_allele": r["nearest_allele"], **info})
    seq_df = pd.DataFrame(seq_rows, columns=["person_id", "hap", "gene", "nearest_allele",
                                             "cds_core", "protein", "frameshift",
                                             "premature_stop", "cds_known", "protein_known",
                                             "mapped_prot2", "mapped_cds3", "seq_class"])
    calls = calls.merge(seq_df, on=["person_id", "hap", "gene"], how="left")
    mstats["n_seq_calls_joined"] = int(calls["seq_class"].notna().sum())

    # Field-class counts (A)
    fcounts = stratified_counts(
        calls[["gene", "gene_class", "field_class", "artifact_label", "anc_pred", "anc_strict",
               "unrelated"]], ["gene", "gene_class", "field_class", "artifact_label"])
    fcounts = fcounts[["gene", "gene_class", "ancestry_scheme", "ancestry", "unrelated_only",
                       "field_class", "artifact_label", "n_haplotypes"]]
    fcounts = fcounts.sort_values(["gene", "ancestry_scheme", "ancestry", "unrelated_only",
                                   "field_class", "artifact_label"])
    fcounts.to_csv(os.path.join(local_dir, "novelty_field_counts.full.tsv"), sep="\t",
                   index=False)
    suppress_df(fcounts, ["n_haplotypes"]).to_csv(
        os.path.join(out_dir, "novelty_field_counts.tsv"), sep="\t", index=False)

    # Sequence class counts (committed)
    sc = calls[calls["seq_class"].notna()]
    seq_counts = (sc.groupby(["gene", "gene_class", "field_class", "seq_class", "artifact_label"])
                  .size().reset_index(name="n_haplotypes"))
    suppress_df(seq_counts, ["n_haplotypes"]).to_csv(
        os.path.join(out_dir, "sequence_class_counts.tsv"), sep="\t", index=False)

    # Clusters
    drift = None
    drift_line = "Not run (`--imgt-latest-nuc` not given)."
    if args.imgt_latest_nuc:
        drift, n_drift = load_drift(args.imgt_latest_nuc)
    exon_path = args.alleles_csv or os.path.join(os.path.expanduser(args.refdata),
                                                 "alleles.csv.gz")
    if not os.path.exists(exon_path):
        found = sorted(glob.glob(os.path.join(os.path.expanduser(args.refdata), "**",
                                              "alleles.csv.gz"), recursive=True))
        exon_path = found[0] if found else None
    gene_class_by_gene = dict(zip(calls["gene"], calls["gene_class"]))
    calls["gene_b"] = calls["gene"].map(gene_bare)
    seq_calls = calls[calls["seq_class"].isin(["novel_protein", "novel_cds_synonymous",
                                               "protein_known"])].copy()
    # exon maps for candidate nearest alleles
    names_needed = set()
    for (gb, prot), m in seq_calls[seq_calls["seq_class"] == "novel_protein"].groupby(
            ["gene_b", "protein"]):
        hints = list(m["nearest_allele"].dropna().astype(str))
        hint = Counter(hints).most_common(1)[0][0] if hints else None
        nn = nearest_known_protein(prot, gb, ref, [hint] if hint else [])
        if nn["name"]:
            names_needed.add(nn["name"])
    exon_maps = load_exon_maps(exon_path, names_needed,
                               set(seq_calls.loc[seq_calls["seq_class"] == "novel_protein",
                                                 "gene_b"]))
    clusters = build_clusters(seq_calls, ref, exon_maps, gene_class_by_gene, drift)
    calls = allele_ids(calls)
    if len(clusters):
        clusters = clusters.sort_values(["cluster_type", "gene", "n_haplotypes"],
                                        ascending=[True, True, False])
    clusters.to_csv(os.path.join(local_dir, "protein_level_novel_clusters.full.tsv"), sep="\t",
                    index=False)
    count_cols = [c for c in clusters.columns if c.startswith(CLUSTER_COUNT_COLS_PREFIX)]
    suppress_df(clusters, count_cols).to_csv(
        os.path.join(out_dir, "protein_level_novel_clusters.tsv"), sep="\t", index=False)
    with open(os.path.join(local_dir, "novel_protein_clusters.faa"), "w") as f:
        for r in seq_calls[seq_calls["seq_class"] == "novel_protein"].drop_duplicates(
                ["gene_b", "protein"]).itertuples():
            f.write(f">{r.gene_b}_prot_{sha8(r.protein)}\n{r.protein}\n")
    if drift is not None and len(clusters):
        np_ = clusters[clusters["cluster_type"] == "novel_protein"]
        ns_ = clusters[clusters["cluster_type"] == "novel_cds_synonymous"]
        drift_line = (
            f"`{args.imgt_latest_nuc}` ({n_drift} records): "
            f"{int(np_['in_imgt_latest_protein'].sum())}/{len(np_)} novel-protein clusters and "
            f"{int(np_['in_imgt_latest_cds'].sum() + ns_['in_imgt_latest_cds'].sum())}/"
            f"{len(clusters)} novel-CDS clusters are present in the newer release.")

    # Allele count tables (D)
    act, dropped = allele_count_table(calls)
    act.to_csv(os.path.join(local_dir, "allele_counts_by_resolution.tsv"), sep="\t",
               index=False)

    # Summaries
    npc = clusters[clusters["cluster_type"] == "novel_protein"] if len(clusters) else clusters
    headline = {
        "n_haplotype_gene_calls": int(len(calls)),
        "n_novel_calls_any_depth": int(calls["field_class"].str.startswith("f").sum()),
        "n_f2_protein_calls": int((calls["field_class"] == "f2_protein").sum()),
        "n_f2_protein_calls_clean": int(((calls["field_class"] == "f2_protein")
                                         & (calls["artifact_label"] == "clean")).sum()),
        "n_f3_synonymous_calls": int((calls["field_class"] == "f3_synonymous").sum()),
        "n_f4_noncoding_calls": int((calls["field_class"] == "f4_noncoding").sum()),
        "n_depth23_calls_with_sequence": int(calls["seq_class"].notna().sum()),
        "n_cds_known_naming_artifacts": int((calls["seq_class"] == "cds_known").sum()),
        "n_frameshift_or_stop_calls": int((calls["seq_class"] == "frameshift_or_stop").sum()),
        "n_novel_protein_calls": int((calls["seq_class"] == "novel_protein").sum()),
        "n_novel_protein_clusters": int(len(npc)),
        "n_novel_protein_clusters_all_clean": int(npc["all_members_clean"].sum())
        if len(npc) else 0,
        "n_novel_protein_clusters_clean_recurrent_unrelated":
            int(npc["clean_recurrent_unrelated"].sum()) if len(npc) else 0,
        "n_novel_cds_synonymous_clusters": int((clusters["cluster_type"]
                                                == "novel_cds_synonymous").sum())
        if len(clusters) else 0,
        "frac_f2_calls_flagged": float(
            (calls.loc[calls["field_class"] == "f2_protein", "artifact_label"] != "clean").mean())
        if (calls["field_class"] == "f2_protein").any() else None,
    }
    S = {
        "limit": args.limit, "n_people": int(len(people)),
        "n_people_unrelated": int(people["unrelated"].sum()),
        "n_people_removed_for_relatedness": int(n_removed),
        "relatedness_skipped": bool(args.skip_relatedness),
        "kin_min": args.kin_min, "strict_threshold": args.strict_threshold,
        "headline": headline, "match_stats": {k: int(v) for k, v in mstats.items()},
        "allele_table_dropped": dropped,
        "n_depth3_protein_not_catalogued": int(((calls["depth"] == 3)
                                                & (calls["seq_class"] == "novel_protein")).sum()),
        "n_depth2_protein_catalogued": int((calls["seq_class"] == "protein_known").sum()),
        "refdata": {"n_files": len(ref.files), "n_records": ref.n_records,
                    "n_with_terminal_stop": ref.n_with_terminal_stop,
                    "n_nonstandard_frame": ref.n_nonstandard_frame,
                    "n_unparsed_headers": ref.n_unparsed_headers},
        "alleles_csv": exon_path, "n_exon_maps": len(exon_maps),
        "imgt_version_alleles_csv": exon_maps.version,
        "exon_map_genes_with_mode": sorted(exon_maps.gene_mode),
        "groove_method_counts": (npc["groove_method"].value_counts().to_dict()
                                 if len(npc) else {}),
        "funnel": {},
        "runtime_s": None,
    }
    fc = funnel_counts(calls, clusters)
    S["funnel"] = {k: [[s.replace("\n", " "), int(v)] for s, v in steps]
                   for k, steps in fc.items()}

    # Report tables
    tot = calls.groupby("gene").size()
    fp = calls.groupby(["gene", "field_class"]).size().unstack(fill_value=0)
    ff = calls[calls["artifact_label"] != "clean"].groupby(
        ["gene", "field_class"]).size().unstack(fill_value=0)
    rows = []
    for g in _gene_order(calls):
        rec = {"gene": g, "n_calls": suppress(tot[g])}
        for fcl in FIELD_CLASSES:
            n = int(fp.loc[g, fcl]) if fcl in fp.columns else 0
            nf = int(ff.loc[g, fcl]) if (g in ff.index and fcl in ff.columns) else 0
            rec[fcl] = f"{n / tot[g]:.3f} [{nf / n:.2f}]" if n else "0"
        rows.append(rec)
    field_pooled = pd.DataFrame(rows)
    seq_tab = (sc.groupby(["gene_class", "seq_class"]).size().unstack(fill_value=0)
               if len(sc) else pd.DataFrame())
    if len(seq_tab):
        seq_tab = seq_tab.reset_index()
        for c in seq_tab.columns[1:]:
            seq_tab[c] = seq_tab[c].map(suppress)
    if len(clusters):
        cs = clusters.groupby(["gene", "cluster_type"]).agg(
            n_clusters=("cluster_id", "size"),
            n_all_clean=("all_members_clean", "sum"),
            n_clean_recurrent_unrelated=("clean_recurrent_unrelated", "sum")).reset_index()
    else:
        cs = pd.DataFrame()
    if len(npc):
        top = npc[npc["all_members_clean"]].sort_values("n_persons_unrelated_clean",
                                                        ascending=False).head(15)
        top = top[["cluster_id", "gene", "n_persons_unrelated_clean", "nearest_known_protein",
                   "nearest_distance", "aa_diffs", "in_groove", "groove_method"]].copy()
        top["n_persons_unrelated_clean"] = top["n_persons_unrelated_clean"].map(suppress)
        gt = npc.groupby(["gene", "groove_method", "groove_na_reason"]).size().reset_index(
            name="n_clusters")
    else:
        top = pd.DataFrame()
        gt = pd.DataFrame()

    print("Figures ...", file=sys.stderr)
    fig_field_class_by_gene(calls, os.path.join(out_dir, "fig_field_class_by_gene.png"))
    fig_rate_by_ancestry(calls, os.path.join(out_dir, "fig_clean_novelty_rate_by_ancestry.png"),
                         args.min_denominator)
    fig_funnel(fc, os.path.join(out_dir, "fig_novelty_funnel.png"))

    S["runtime_s"] = round(time.time() - t0, 1)
    write_report(os.path.join(out_dir, "novelty_by_field_report.md"), S, field_pooled, seq_tab,
                 cs, gt, drift_line, top, args.limit)
    with open(os.path.join(out_dir, "summary.json"), "w") as f:
        json.dump(suppress_tree(S), f, indent=2, default=str)
    print(json.dumps(headline, indent=1), file=sys.stderr)
    print(f"Done in {S['runtime_s']}s. Committed: {out_dir}  VM-local: {local_dir}",
          file=sys.stderr)
    return S


if __name__ == "__main__":
    main()
