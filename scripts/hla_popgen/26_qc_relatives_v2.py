#!/usr/bin/env python3
"""Non-circular QC of the long-read (Immuannot) HLA calls using real relative pairs -- the
replacement for 16/17's circular, denominator-free version, per
`sprints/S01_novelty_qc_coverage/WS2_qc_audit.md` ("script 26").

## What the audit found wrong with v1 (16/17), and what this script does instead

1. **v1 selected pairs on the statistic it reported.** A pair was "high_sharing" if >=80% of its
   genes shared an allele, and the headline error rate was then the fraction of non-sharing genes
   *in those same pairs*. Here the pair universe is defined ONLY by AoU's own genome-wide kinship
   (kin >= 0.177, i.e. first degree, plus duplicates/MZ twins at kin > 0.354). Nothing is selected
   on HLA sharing. The sharing fraction is then *reported* as a distribution, and every downstream
   rate is reported at a SWEEP of sharing thresholds (0.5/0.8/0.9/0.95) so the reader sees how
   sensitive the number is to that choice instead of seeing one hidden cut.
2. **v1 dropped genes missing in either person from the denominator**, so allele dropout and
   assembly fragmentation were invisible. Here every gene-comparison lands in exactly one
   mutually-exclusive category, including `missing_one` and `missing_both`.
3. **v1's "0/545 phase switches" had no denominator and no power estimate.** Here the number of
   resolved genes and of *testable same-contig transitions* per pair is reported as a distribution
   (that denominator is the key missing number), and the detector's sensitivity is measured
   empirically by inserting a synthetic hap-label swap into one person of each eligible pair and
   asking how often the same detector sees it.
4. **v1 compared allele NAMES.** Immuannot's naming step collapses tied candidates, so a 1-base
   difference can flip a name and two different sequences can share one. Every sharing /
   discordance decision here is made on the exact observed CDS **sequence** multiset read from
   `cds.fa.gz` (same parser as 17). Names are used in exactly one place: the orthogonal SR-vs-LR
   2-field concordance step, where the short-read side only ever has a name.
5. **v1 had no orthogonal (non-relative) check.** Step 7 compares each person's long-read calls to
   AoU's own short-read HLA genotypes for the same person at 2 fields, split by whether the
   long-read call is novel and at what depth -- a novelty that is noncoding-only (depth 4) must
   NOT reduce 2-field concordance, while a protein-altering one (depth 2) should. That prediction
   is falsifiable and does not depend on relatives at all.

## The eight steps (each one a separately testable pure function; I/O only in the wrappers)

1. `build_pair_universe`      -- kinship-only pair selection, no HLA-based filtering.
2. `pair_sharing` / `sharing_threshold_sweep` -- sequence-multiset sharing per gene, pair-level
   fraction over classical + class-II-accessory genes, reported as a distribution + a threshold
   sweep (and, with `--yob-tsv`, split by age gap >= 15y as a likely parent-child proxy).
3. `classify_gene_comparison` -- the mutually exclusive decomposition of every gene-comparison.
4. `replicate_gene_diff` / `replicate_qv` -- duplicate/MZ pairs get full zygosity-aware
   hap-to-hap equality and a per-base QV.
5. `homozygosity_obs_exp`     -- observed vs expected homozygosity per gene x strict ancestry, for
   LR-CDS-exact / LR-2field / SR-2field on the SAME people, plus a split by assembly
   fragmentation.
6. `resolve_states`, `testable_transitions`, `count_switches`, `apply_synthetic_switch`,
   `switch_power` -- the switch test WITH its denominator and an empirical power estimate.
7. `sr_lr_concordance_rows`   -- per person-gene 0/1/2 alleles matching at 2 fields, by gene x
   strict ancestry x LR field class.
8. figures + `qc_relatives_v2_report.md` + `summary.json`.

## Definitions this script commits to (read before interpreting any number)

- **Strict ancestry**: a person is counted under ancestry label L only if their own admixture
  proportion `p_<l>` (Table 4) is >= `--strict-ancestry-min` (default 0.9). `ancestry_pred` alone
  is a hard assignment of an admixed person to one bin, which is exactly wrong for a
  frequency-based expectation like HWE homozygosity.
- **Expected homozygosity** is `sum_i x_i(x_i-1) / (n(n-1))` over that ancestry's own observed
  allele-copy counts `x_i` at the same resolution (`n` = total copies). That is the standard
  finite-sample-unbiased estimator of `sum_i p_i^2`; the naive `sum (x_i/n)^2` is biased UPWARD
  in small samples, which would make the obs/exp ratio look artificially close to 1 exactly in
  the small ancestry groups where dropout matters most.
- **A "resolved" gene for the switch test** is one where the two people share exactly one distinct
  CDS sequence and that sequence sits on exactly one haplotype in each person (so the hap-to-hap
  linkage is unambiguous). A transition between two adjacent resolved genes is **testable** only
  when both people's resolved hap stayed on the same contig -- a hap1/hap2 label is only
  phase-consistent within one contig (SCHEMA.md).
- **Aggregate-only outputs.** Nothing person-level is written to `--out-dir`; person counts 1-19
  are written as `<20`; person ids appear nowhere in `--out-dir`. Person-level pair rows are
  written only if `--local-dir` is given/kept (default `~/pipeline_outputs/s01_local/`).

Usage (VM, real run):
    pixi run -e spechla -- python3 scripts/hla_popgen/26_qc_relatives_v2.py \\
        --sr-genotypes ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv

Pilot on the first 20 pairs (do this first; it reads two gzipped FASTAs per person):
    pixi run -e spechla -- python3 scripts/hla_popgen/26_qc_relatives_v2.py \\
        --sr-genotypes ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv \\
        --limit-pairs 20 --out-dir ~/results/26_qc_relatives_v2_pilot

Controlled-tier data is auto-mounted at `~/workspace/vwb-aou-datasets-controlled-v9/`; the older
`~/mnt/aou-controlled` gcsfuse path is stale and puts processes into uninterruptible I/O, so it is
referenced nowhere here, not even as a fallback, and no path is stat'ed at import time.

Usage (local fixtures): see scripts/hla_popgen/tests/test_qc_relatives_v2.py, which builds a
complete tiny dataset (Table 1 + cds.fa.gz + relatedness + cohort_membership + SR genotypes) and
runs `main_from_args` end to end.
"""
import argparse
import hashlib
import importlib
import json
import math
import os
import sys
import zlib
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _viz_common as vc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

# Reuse the already-tested pure helpers rather than re-deriving them (WS2 brief): 11's KING
# kinship binning, 16's canonical gene order + contig-aware switch counter, 17's cds.fa.gz parser
# and best-pairing sequence diff.
mod11 = importlib.import_module("11_relatedness_cohort_overlap")
mod16 = importlib.import_module("16_phasing_mendelian_validation")
mod17 = importlib.import_module("17_raw_sequence_divergence")

kin_degree = mod11.kin_degree
derive_canonical_order = mod16.derive_canonical_order
count_switches = mod16.count_switches           # identical logic, reused not re-implemented
best_alignment_diff = mod17.best_alignment_diff
compare_gene_seqs = mod17.compare_gene
parse_cds_fasta = mod17.parse_cds_fasta

AOU_CONTROLLED_ROOT = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9")
DEFAULT_RELATEDNESS = os.path.join(
    AOU_CONTROLLED_ROOT, "v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv")
DEFAULT_SR_GENOTYPES = os.path.join(
    AOU_CONTROLLED_ROOT, "v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv")
DEFAULT_PEOPLE_ROOT = os.path.expanduser("~/pipeline_outputs/people")
DEFAULT_OUT_DIR = os.path.expanduser("~/results/26_qc_relatives_v2")
DEFAULT_LOCAL_DIR = os.path.expanduser("~/pipeline_outputs/s01_local")

FIRST_DEGREE_MIN_KIN = 0.177          # KING first-degree floor (Manichaikul 2010), via 11
DUPLICATE_MIN_KIN = 0.354             # above this: duplicate / MZ twin
SHARING_THRESHOLDS = (0.5, 0.8, 0.9, 0.95)
AGE_GAP_PARENT_CHILD = 15             # years; a proxy, not a pedigree
MIN_ANCESTRY_PEOPLE = 50              # below this a group is reported but flagged as thin
SMALL_COUNT_CEILING = 20              # person/pair counts 1..19 are written as "<20"

CLASSICAL = list(vc.CLASSICAL_GENES_BARE)
# Pair-level sharing fraction is computed over classical + class-II-accessory genes only (WS2
# step 2): pseudogenes and DRB paralogs have both the highest paralog-misassignment rate (12-14%
# at HLA-H/K/U per the audit) and copy-number ambiguity, so including them would inject exactly
# the noise the IBD-state estimate must be robust to.
SHARING_GENE_CLASSES = ("classical_I", "classical_II", "class_II_accessory")

DISCORDANCE_CATEGORIES = ["concordant", "missing_one", "missing_both", "point_diff",
                           "dropout_candidate", "fragmentation_candidate", "other_multiblock"]
# Colorblind-safe (Okabe-Ito), ordered to read as "fine -> missing -> technical -> real".
CATEGORY_COLORS = {
    "concordant": "#0072B2", "missing_one": "#999999", "missing_both": "#666666",
    "point_diff": "#56B4E9", "dropout_candidate": "#E69F00",
    "fragmentation_candidate": "#CC79A7", "other_multiblock": "#D55E00",
}
POINT_DIFF_MAX_BASES = 3              # audit: 57% of divergences are <=3 bases (HiFi-consensus-like)
FIELD_CLASS_ORDER = ["known", "f4_noncoding", "f3_synonymous", "f2_protein"]
FIELD_CLASS_COLORS = {"known": "#0072B2", "f4_noncoding": "#009E73",
                      "f3_synonymous": "#E69F00", "f2_protein": "#D55E00"}


# ===========================================================================
# Small formatting / safety helpers
# ===========================================================================
def safe_count(n):
    """Person- and pair-level counts 1..19 are disclosed as the string '<20' (WS2 step 8 / SCHEMA
    hard rule 5). 0 is not re-identifying and is reported as 0 so a reader can tell "none" from
    "a few"."""
    if n is None or (isinstance(n, float) and math.isnan(n)):
        return "NA"
    n = int(n)
    if n == 0:
        return "0"
    if n < SMALL_COUNT_CEILING:
        return f"<{SMALL_COUNT_CEILING}"
    return str(n)


def _md_table(df):
    """Local markdown writer -- `tabulate` is not present in the spechla pixi env (see 11)."""
    cols = [str(c) for c in df.columns]
    rows = ["| " + " | ".join(str(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(["| " + " | ".join(cols) + " |",
                      "|" + "|".join(["---"] * len(cols)) + "|"] + rows)


def _wilson(count, nobs):
    p, lo, hi = vc.wilson_ci(np.array([count], dtype=float), np.array([nobs], dtype=float))
    return float(p[0]), float(lo[0]), float(hi[0])


# ===========================================================================
# STEP 1 -- pair universe, selected on kinship ONLY
# ===========================================================================
def build_pair_universe(rel, cohort, min_kin=FIRST_DEGREE_MIN_KIN):
    """All pairs where BOTH people are in the long-read cohort and AoU's own genome-wide kinship
    is >= `min_kin` (first degree), duplicates/MZ twins included. Returns a DataFrame with
    person_i, person_j, kin, kin_degree, is_replicate.

    This is the fix for audit finding 1: nothing about HLA sharing enters the selection, so the
    sharing fraction and every discordance rate computed downstream are free to be whatever they
    are instead of being conditioned on being high."""
    lr_ids = set(cohort.loc[cohort["in_lr"].astype(bool), "person_id"])
    df = rel.copy()
    df["kin"] = pd.to_numeric(df["kin"], errors="coerce")
    df = df.dropna(subset=["kin"])
    df = df[df["person_i"].isin(lr_ids) & df["person_j"].isin(lr_ids)]
    df = df[df["kin"] >= min_kin]
    df["kin_degree"] = df["kin"].map(kin_degree)
    df["is_replicate"] = df["kin"] > DUPLICATE_MIN_KIN
    df = df.drop_duplicates(subset=["person_i", "person_j"])
    return df.sort_values(["kin"], ascending=False).reset_index(drop=True)


def strict_ancestry_map(cohort, min_prop=0.9):
    """{person_id: ANCESTRY} keeping ONLY people whose own admixture proportion for their predicted
    label is >= min_prop (Table 4 columns `p_afr`..`p_sas`). Everyone else is deliberately absent:
    an HWE-style expectation computed from an admixed person's "assigned" population is not a
    like-for-like expectation, and the obs/exp homozygosity ratio in step 5 is exactly the place
    where that mistake would masquerade as allele dropout."""
    out = {}
    for row in cohort.itertuples(index=False):
        anc = getattr(row, "ancestry_pred", None)
        if not isinstance(anc, str) or anc not in vc.ANCESTRY_ORDER:
            continue
        col = f"p_{anc.lower()}"
        prop = getattr(row, col, None)
        if prop is None or pd.isna(prop):
            continue
        if float(prop) >= min_prop:
            out[row.person_id] = anc
    return out


def load_year_of_birth(path):
    """Optional {person_id: year_of_birth}. Absent file -> empty dict, and every age-gap number is
    then skipped and SAID to be skipped in the report (never silently omitted)."""
    if not path or not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t", dtype={"person_id": str})
    ycol = "year_of_birth" if "year_of_birth" in df.columns else df.columns[-1]
    df[ycol] = pd.to_numeric(df[ycol], errors="coerce")
    return {r.person_id: float(getattr(r, ycol)) for r in df.itertuples(index=False)
            if not pd.isna(getattr(r, ycol))}


def age_gap_class(person_i, person_j, yob, gap=AGE_GAP_PARENT_CHILD):
    """'likely_parent_child' (|dYOB| >= gap), 'likely_sibling' (< gap), or 'unknown' when either
    year of birth is missing. A proxy for the pedigree the AoU relatedness table does not carry
    (11: no IBD0 column), not a claim about the real relationship."""
    a, b = yob.get(person_i), yob.get(person_j)
    if a is None or b is None:
        return "unknown"
    return "likely_parent_child" if abs(a - b) >= gap else "likely_sibling"


# ===========================================================================
# Per-person sequence + contig assembly (thin I/O + a pure reshaper)
# ===========================================================================
TABLE1_COLUMNS = ["person_id", "hap", "contig", "gene", "gene_class", "consensus",
                   "copy_index", "is_novel", "novelty_depth", "gene_start"]


def load_table1_lean(path):
    """Table 1, but only the ten columns this script actually uses, with the low-cardinality ones
    as `category`. At full scale Table 1 is ~1.5M rows x 25 columns and the VM has 31GB shared with
    everything else; `vc.load_table1` would pull in `cds_mut`, `alleles` and friends (long strings)
    for no reason. Same derivations as `vc.load_table1` otherwise, so downstream code is identical.

    `person_id`, `hap` and `contig` are deliberately NOT categorical: `derive_canonical_order`
    groups by all three, and a categorical groupby materializes the full cartesian product of
    categories (12k x 2 x n_contigs empty groups)."""
    if not os.path.exists(path):
        sys.exit(f"FATAL: Table 1 (hla_calls_rich.tsv) not found at {path!r}. "
                 f"Run 01_extract_rich.py first.")
    head = pd.read_csv(path, sep="\t", nrows=0)
    usecols = [c for c in TABLE1_COLUMNS if c in head.columns]
    missing = {"person_id", "hap", "contig", "gene", "consensus"} - set(usecols)
    if missing:
        sys.exit(f"FATAL: Table 1 at {path!r} is missing required columns {sorted(missing)} -- "
                 f"schema changed? See SCHEMA.md Table 1.")
    df = pd.read_csv(path, sep="\t", usecols=usecols,
                     dtype={"person_id": str, "hap": str, "contig": str,
                            "gene": "category", "gene_class": "category", "consensus": str})
    if "is_novel" in df.columns:
        df["is_novel"] = vc.parse_bool_column(df["is_novel"])
    for c in ("novelty_depth", "copy_index", "gene_start"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["gene_bare"] = df["gene"].astype(str).map(vc.gene_bare)
    return df


def seq_digest(seq):
    """Short stable digest of a CDS. Homozygosity and fragmentation only ever ask "is this the
    same sequence as that one", so the cohort-wide step keeps 16 hex chars per hap instead of a
    ~1-3kb string -- the difference between ~10MB and ~2GB at 12k people x 8 genes."""
    return hashlib.sha1(seq.encode()).hexdigest()[:16]


def load_person_cds(people_root, person_id, genes):
    """{(hap, gene_bare): {copy_index: seq}} for one person, using 17's `cds.fa.gz` parser (which
    already handles the real header format `>{contig}_{gene}_{i} {splice_sites}`). copy_index is
    kept here -- unlike 17's own loader, which flattens it -- because step 3 is defined only on
    copy_index==1 rows and a DRB paralog's second copy must not be silently mixed in."""
    want = {g: ("HLA-" + g if g not in mod17._MIC_TAP else g) for g in genes}
    out = {}
    for hap in ("hap1", "hap2"):
        path = os.path.join(people_root, person_id, "immuannot_output", hap, "cds.fa.gz")
        idx = parse_cds_fasta(path)
        for g, key in want.items():
            recs = idx.get(key)
            if recs:
                # REAL-DATA TRAP (found 2026-09-16 against the production cohort): the trailing
                # integer in a cds.fa.gz header `>{contig}_{gene}_{i}` is the record's ORDINAL
                # WITHIN THE FILE (HLA-E_1, HLA-L_2, HLA-K_3 ... across different genes), NOT the
                # gene's copy_index. Keying on it and taking "copy 1" silently dropped the
                # sequence of every gene except whichever happened to be the file's first record,
                # which made 78% of relative gene-comparisons look missing/discordant.
                # A gene with exactly one record on this hap is unambiguous; a gene with several
                # records is a real multi-copy case (DRB paralogs, segmental duplication) that
                # cannot be mapped onto Table 1's copy_index from this file alone, so it is kept
                # as a list and excluded by callers that are defined on copy_index==1.
                out[(hap, g)] = [seq for _, seq in recs]
    return out


def build_person_map(cds, contig_lookup, name_lookup, person_id, genes, digest=False):
    """Pure reshaper: -> {gene: {"seq": {hap: seq}, "contig": {hap: contig}, "name": {hap: str}}}
    restricted to copy_index==1 (SCHEMA Table 1 grain; >1 is real segmental-duplication copy
    number, a different question). Contigs and 2-field names come from Table 1, which is the only
    place the contig of a (person, hap, gene) is recorded."""
    pmap = {}
    for gene in genes:
        entry = {"seq": {}, "contig": {}, "name": {}}
        for hap in ("hap1", "hap2"):
            recs = cds.get((hap, gene)) or []
            seq = recs[0] if len(recs) == 1 else None  # >1 record: ambiguous copy, see load_person_cds
            if seq is not None:
                entry["seq"][hap] = seq_digest(seq) if digest else seq
            ctg = contig_lookup.get((person_id, gene, hap))
            if ctg is not None:
                entry["contig"][hap] = ctg
            nm = name_lookup.get((person_id, gene, hap))
            if nm is not None:
                entry["name"][hap] = nm
        if entry["seq"] or entry["name"]:
            pmap[gene] = entry
    return pmap


def build_table1_lookups(table1, resolution=2):
    """(contig_lookup, name_lookup, gene_class, novelty_lookup), all keyed
    (person_id, gene_bare, hap), copy_index==1 only."""
    sub = table1[table1["copy_index"].fillna(1) == 1]
    sub = sub[sub["hap"].isin(["hap1", "hap2"])]
    contig_lookup, name_lookup, novelty = {}, {}, {}
    for row in sub.itertuples(index=False):
        key = (row.person_id, row.gene_bare, row.hap)
        contig_lookup[key] = row.contig
        nm = vc.to_nfield(row.consensus, resolution)
        if nm is not None:
            name_lookup[key] = nm
        novelty[key] = (bool(getattr(row, "is_novel", False)),
                        getattr(row, "novelty_depth", None))
    gene_class = dict(zip(sub["gene_bare"], sub["gene_class"]))
    return contig_lookup, name_lookup, gene_class, novelty


def observed_seqs(pmap, gene):
    """The person's observed CDS multiset at `gene`, as a set of distinct sequences (1 if
    homozygous-by-sequence, 2 if heterozygous, 0 if the gene is absent from both haps)."""
    return set((pmap.get(gene, {}).get("seq") or {}).values())


def is_sequence_homozygous(pmap, gene):
    """True when the person's observed CDS at this gene is a single distinct sequence -- either
    both haps carry byte-identical CDS, or only one hap carries the gene at all. The second case
    is included on purpose: a collapsed heterozygous locus is exactly what allele dropout looks
    like, and it presents as "one hap has it, the other doesn't" just as often as it presents as
    two identical copies."""
    seqs = observed_seqs(pmap, gene)
    return len(seqs) == 1


def pair_gene_shared(a_map, b_map, gene):
    """Phase-agnostic sharing on the EXACT sequence (not the name): do the two people's observed
    CDS multisets intersect at all?"""
    a, b = observed_seqs(a_map, gene), observed_seqs(b_map, gene)
    if not a or not b:
        return None                      # not comparable; step 3 categorizes this, step 2 skips it
    return bool(a & b)


# ===========================================================================
# STEP 2 -- pair-level sharing fraction, reported not selected on
# ===========================================================================
def pair_sharing(a_map, b_map, sharing_genes):
    """(n_shared, n_compared, fraction) over `sharing_genes` where BOTH people have a sequence.
    Genes missing on either side are excluded HERE (this is an IBD-state estimate, and a missing
    gene carries no IBD information) but are counted, never dropped, in step 3's decomposition --
    that split is the whole point of audit finding 2."""
    n_shared = n_compared = 0
    for gene in sharing_genes:
        shared = pair_gene_shared(a_map, b_map, gene)
        if shared is None:
            continue
        n_compared += 1
        n_shared += int(shared)
    frac = (n_shared / n_compared) if n_compared else float("nan")
    return n_shared, n_compared, frac


def sharing_threshold_sweep(pair_rows, thresholds=SHARING_THRESHOLDS):
    """For each candidate "this pair is IBD>=1" threshold, what would the headline discordance rate
    have been? Shows the sensitivity of v1's hidden 0.8 cut instead of hiding it.

    `pair_rows`: list of dicts with sharing_frac, n_gene_comparisons, n_discordant (the step-3
    counts for that pair). Returns a DataFrame."""
    out = []
    for t in thresholds:
        kept = [r for r in pair_rows
                if not pd.isna(r["sharing_frac"]) and r["sharing_frac"] >= t]
        n_cmp = sum(r["n_gene_comparisons"] for r in kept)
        n_dis = sum(r["n_discordant"] for r in kept)
        p, lo, hi = _wilson(n_dis, n_cmp) if n_cmp else (float("nan"),) * 3
        out.append({"threshold": t, "n_pairs": len(kept), "n_gene_comparisons": n_cmp,
                    "n_discordant": n_dis, "discordance_rate": p, "ci_lo": lo, "ci_hi": hi})
    return pd.DataFrame(out)


# ===========================================================================
# STEP 3 -- mutually exclusive discordance decomposition
# ===========================================================================
def singleton_contig_genes(pmap, classical=tuple(CLASSICAL)):
    """Genes whose contig (on at least one hap) carries NO other classical gene -- i.e. the gene
    sits on its own assembly fragment, where nothing anchors it and a mis-assignment or a dropped
    haplotype is most likely. Used for the `fragmentation_candidate` category."""
    per_contig = defaultdict(set)
    for gene in classical:
        for hap, ctg in (pmap.get(gene, {}).get("contig") or {}).items():
            if ctg is not None:
                per_contig[(hap, ctg)].add(gene)
    out = set()
    for gene in classical:
        for hap, ctg in (pmap.get(gene, {}).get("contig") or {}).items():
            if ctg is not None and len(per_contig[(hap, ctg)]) <= 1:
                out.add(gene)
    return out


def classify_gene_comparison(a_map, b_map, gene, a_singleton, b_singleton,
                             point_diff_max=POINT_DIFF_MAX_BASES):
    """The one category this (pair, gene) comparison falls into. Exactly one, always, with the
    documented precedence (WS2 step 3):

        concordant
        > missing_both  > missing_one          (presence beats everything: no sequence, no verdict)
        > point_diff                           (<= `point_diff_max` bases in ONE block)
        > dropout_candidate                    (either person looks homozygous-by-sequence here)
        > fragmentation_candidate              (gene alone on a contig in either person)
        > other_multiblock                     (the residual -- the audit's "real concern" bucket)

    Returns (category, n_diff_bases, n_diff_blocks); the diff numbers are NaN where no pairwise
    comparison was possible.

    Precedence rationale: `point_diff` outranks `dropout_candidate` because a 1-3 base difference
    at a homozygous locus is still overwhelmingly consensus error (audit finding 5: implied
    QV ~46-47, inside the normal HiFi range), and calling it dropout would inflate the dropout
    estimate with plain sequencing noise."""
    a_seqs, b_seqs = observed_seqs(a_map, gene), observed_seqs(b_map, gene)
    if not a_seqs and not b_seqs:
        return "missing_both", float("nan"), float("nan")
    if not a_seqs or not b_seqs:
        return "missing_one", float("nan"), float("nan")
    cmp = compare_gene_seqs(sorted(a_seqs), sorted(b_seqs))
    if cmp["raw_shared"]:
        return "concordant", 0.0, 0.0
    n_bases, n_blocks = float(cmp["n_diff_bases"]), float(cmp["n_diff_blocks"])
    if n_bases <= point_diff_max and n_blocks <= 1:
        return "point_diff", n_bases, n_blocks
    if is_sequence_homozygous(a_map, gene) or is_sequence_homozygous(b_map, gene):
        return "dropout_candidate", n_bases, n_blocks
    if gene in a_singleton or gene in b_singleton:
        return "fragmentation_candidate", n_bases, n_blocks
    return "other_multiblock", n_bases, n_blocks


def decompose_pair(a_map, b_map, genes, classical=tuple(CLASSICAL)):
    """[{gene, category, n_diff_bases, n_diff_blocks}] for every gene in `genes`."""
    a_singleton = singleton_contig_genes(a_map, classical)
    b_singleton = singleton_contig_genes(b_map, classical)
    rows = []
    for gene in genes:
        cat, nb, nbl = classify_gene_comparison(a_map, b_map, gene, a_singleton, b_singleton)
        rows.append({"gene": gene, "category": cat, "n_diff_bases": nb, "n_diff_blocks": nbl})
    return rows


def aggregate_decomposition(records, gene_class):
    """records: iterable of dicts with gene + category (across all pairs). Returns a per-gene
    DataFrame of category counts with the discordance rate and its Wilson CI. `missing_both` is
    kept in the table but excluded from the rate's denominator -- a gene neither person called
    says nothing about agreement between them (it is a *coverage* result, reported separately)."""
    counts = defaultdict(Counter)
    for r in records:
        counts[r["gene"]][r["category"]] += 1
    rows = []
    for gene, c in counts.items():
        total = sum(c.values())
        comparable = total - c["missing_both"]
        discordant = comparable - c["concordant"]
        p, lo, hi = _wilson(discordant, comparable) if comparable else (float("nan"),) * 3
        rows.append({"gene": gene, "gene_class": gene_class.get(gene, "other"),
                     "n_total": total, "n_comparable": comparable, "n_discordant": discordant,
                     "discordance_rate": p, "ci_lo": lo, "ci_hi": hi,
                     **{cat: c[cat] for cat in DISCORDANCE_CATEGORIES}})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values("gene").reset_index(drop=True)


def aggregate_by_gene_class(decomp_df):
    if decomp_df.empty:
        return decomp_df
    agg = decomp_df.groupby("gene_class")[
        ["n_total", "n_comparable", "n_discordant"] + DISCORDANCE_CATEGORIES].sum().reset_index()
    p, lo, hi = vc.wilson_ci(agg["n_discordant"].values, agg["n_comparable"].values)
    agg["discordance_rate"], agg["ci_lo"], agg["ci_hi"] = p, lo, hi
    return agg


# ===========================================================================
# STEP 4 -- technical replicates (duplicate / MZ twin pairs) -> a direct QV
# ===========================================================================
def replicate_gene_diff(a_map, b_map, gene):
    """Zygosity-aware, hap-resolved comparison for a pair that should be the SAME genome. Tries
    both hap-to-hap assignments and keeps the one with fewer total differing bases (hap1/hap2
    labels are arbitrary between two independent assemblies). Returns
    {n_hap_compared, n_diff_bases, n_aligned_bases, both_equal} or None if not comparable on both
    haps."""
    a, b = (a_map.get(gene, {}).get("seq") or {}), (b_map.get(gene, {}).get("seq") or {})
    if len(a) < 1 or len(b) < 1:
        return None
    a_list = [a.get("hap1"), a.get("hap2")]
    b_list = [b.get("hap1"), b.get("hap2")]
    best = None
    for order in ((0, 1), (1, 0)):
        total_bases = total_len = n_cmp = 0
        for i, j in enumerate(order):
            x, y = a_list[i], b_list[j]
            if x is None or y is None:
                continue
            nb, _ = best_alignment_diff(x, y)
            total_bases += nb
            total_len += max(len(x), len(y))
            n_cmp += 1
        if n_cmp == 0:
            continue
        cand = {"n_hap_compared": n_cmp, "n_diff_bases": total_bases,
                "n_aligned_bases": total_len, "both_equal": total_bases == 0 and n_cmp == 2}
        if best is None or cand["n_diff_bases"] < best["n_diff_bases"]:
            best = cand
    return best


def replicate_qv(diff_records):
    """Phred-scale QV from replicate disagreements: err = sum(diff bases) / sum(aligned bases),
    QV = -10 log10(err). Two independent assemblies of the same genome differ by the sum of BOTH
    assemblies' errors, so this is a conservative (pessimistic) per-assembly QV, not an optimistic
    one. Returns (qv, err, n_diff, n_aligned); qv is inf when zero differences were seen."""
    n_diff = sum(r["n_diff_bases"] for r in diff_records)
    n_aln = sum(r["n_aligned_bases"] for r in diff_records)
    if n_aln == 0:
        return float("nan"), float("nan"), 0, 0
    err = n_diff / n_aln
    qv = float("inf") if n_diff == 0 else -10.0 * math.log10(err)
    return qv, err, n_diff, n_aln


# ===========================================================================
# STEP 5 -- homozygosity obs/exp, like-for-like across LR-CDS / LR-2field / SR-2field
# ===========================================================================
def expected_homozygosity(allele_counts):
    """sum_i x_i(x_i-1) / (n(n-1)) -- the finite-sample-unbiased estimator of sum_i p_i^2 (Nei
    1978's correction, same algebra). `allele_counts`: iterable of per-allele copy counts. NaN
    when fewer than 2 copies were observed."""
    x = np.asarray([c for c in allele_counts if c > 0], dtype=float)
    n = x.sum()
    if n < 2:
        return float("nan")
    return float(np.sum(x * (x - 1)) / (n * (n - 1)))


def _obs_exp_from_calls(calls):
    """calls: list of (a1, a2) allele identities (one per person). Returns (obs, exp, n_people)."""
    if not calls:
        return float("nan"), float("nan"), 0
    obs = float(np.mean([1.0 if a1 == a2 else 0.0 for a1, a2 in calls]))
    counter = Counter()
    for a1, a2 in calls:
        counter[a1] += 1
        counter[a2] += 1
    return obs, expected_homozygosity(counter.values()), len(calls)


def homozygosity_obs_exp(calls_by_group, n_boot=200, seed=0, min_people=MIN_ANCESTRY_PEOPLE):
    """calls_by_group: {(measure, gene, ancestry): [(a1, a2), ...]} -- one entry per PERSON, each
    an unordered pair of allele identities at the same resolution as that measure.

    For each group: observed homozygosity (both identities equal), expected homozygosity from that
    group's OWN allele-copy counts (finite-sample corrected), their ratio, and a percentile
    bootstrap CI of the ratio resampling PEOPLE (not alleles -- the unit of independence is a
    person, and resampling alleles would understate the uncertainty by treating a person's two
    copies as two independent draws). Deterministic given `seed`.

    Groups with fewer than `min_people` people are returned with thin_n=True, never dropped: a
    dropout signal in a small ancestry group is exactly the finding that must not be hidden."""
    rows = []
    for (measure, gene, anc), calls in sorted(calls_by_group.items()):
        obs, exp, n = _obs_exp_from_calls(calls)
        ratio = obs / exp if (exp and not pd.isna(exp) and exp > 0) else float("nan")
        lo = hi = float("nan")
        if n >= 2 and n_boot > 0:
            # zlib.crc32, not hash(): Python salts str hashes per process, which would make
            # the "deterministic" bootstrap differ between runs.
            key = zlib.crc32(f"{measure}|{gene}|{anc}".encode())
            rng = np.random.default_rng([int(seed), key])
            ratios = []
            idx_space = np.arange(n)
            for _ in range(n_boot):
                idx = rng.choice(idx_space, size=n, replace=True)
                b_obs, b_exp, _ = _obs_exp_from_calls([calls[i] for i in idx])
                if b_exp and not pd.isna(b_exp) and b_exp > 0:
                    ratios.append(b_obs / b_exp)
            if len(ratios) >= 10:
                lo, hi = (float(np.percentile(ratios, 2.5)),
                          float(np.percentile(ratios, 97.5)))
        rows.append({"measure": measure, "gene": gene, "ancestry": anc, "n_people": n,
                     "obs_hom": obs, "exp_hom": exp, "obs_exp_ratio": ratio,
                     "ci_lo": lo, "ci_hi": hi, "thin_n": n < min_people})
    return pd.DataFrame(rows)


def fragmentation_class(pmap, classical=tuple(CLASSICAL)):
    """How many distinct contigs carry this person's classical genes, per hap, expressed as the
    max over the two haps: '1', '2' or '>=3'. Returns None when no classical gene has a contig.
    Fragmentation is the competing explanation for excess homozygosity (a gene on its own fragment
    can simply be missing from the other hap), so homozygosity is reported split by this."""
    per_hap = defaultdict(set)
    for gene in classical:
        for hap, ctg in (pmap.get(gene, {}).get("contig") or {}).items():
            if ctg is not None:
                per_hap[hap].add(ctg)
    if not per_hap:
        return None
    n = max(len(v) for v in per_hap.values())
    return "1" if n <= 1 else ("2" if n == 2 else ">=3")


def homozygosity_by_fragmentation(frag_calls):
    """frag_calls: {(gene, frag_class): [(a1, a2), ...]}. Observed homozygosity only (an expected
    value per fragmentation class would be circular: fragmentation changes which alleles are
    observable at all)."""
    rows = []
    for (gene, frag), calls in sorted(frag_calls.items()):
        obs, _, n = _obs_exp_from_calls(calls)
        rows.append({"gene": gene, "fragmentation_class": frag, "n_people": n, "obs_hom": obs})
    return pd.DataFrame(rows)


# ===========================================================================
# STEP 6 -- the switch test WITH a denominator and an empirical power estimate
# ===========================================================================
def resolve_states(a_map, b_map, gene_order, genes=None):
    """16's linkage-state walk, re-expressed on exact CDS SEQUENCES instead of allele names.

    A gene gets a state (a_hap_index, b_hap_index) only when the two people share exactly ONE
    distinct sequence and that sequence sits on exactly one hap in each person. Genes with no
    state carry no phase information (they are not evidence of a switch by themselves). Each state
    also records which contig backed the chosen hap on each side, so `count_switches` can refuse
    to compare state across a contig boundary."""
    genes = genes if genes is not None else sorted(set(a_map) | set(b_map))
    ordered = sorted(genes, key=lambda g: gene_order.get(g, float("inf")))
    detail = []
    for gene in ordered:
        a_seq = (a_map.get(gene, {}).get("seq") or {})
        b_seq = (b_map.get(gene, {}).get("seq") or {})
        shared_seqs = set(a_seq.values()) & set(b_seq.values())
        state = a_contig = b_contig = None
        if len(shared_seqs) == 1:
            s = next(iter(shared_seqs))
            a_haps = [h for h, v in a_seq.items() if v == s]
            b_haps = [h for h, v in b_seq.items() if v == s]
            if len(a_haps) == 1 and len(b_haps) == 1:
                state = (a_haps[0], b_haps[0])
                a_contig = (a_map[gene].get("contig") or {}).get(a_haps[0])
                b_contig = (b_map[gene].get("contig") or {}).get(b_haps[0])
        detail.append({"gene": gene, "shared": bool(shared_seqs), "state": state,
                       "a_contig": a_contig, "b_contig": b_contig})
    return detail


def testable_transitions(detail):
    """THE MISSING DENOMINATOR (audit finding 3). Number of adjacent resolved-gene pairs where
    both people's resolved hap stayed on the same contig -- i.e. the number of places where a
    switch COULD have been observed at all. Also returns the resolved-gene count."""
    resolved = [d for d in detail if d["state"] is not None]
    n_testable = sum(
        1 for prev, cur in zip(resolved, resolved[1:])
        if cur["a_contig"] == prev["a_contig"] and cur["b_contig"] == prev["b_contig"])
    return n_testable, len(resolved)


def apply_synthetic_switch(pmap, detail, rng):
    """Insert ONE synthetic hap-label swap into `pmap` (a copy is returned; the input is not
    mutated), at a random breakpoint inside a contig -- the ground truth for the power estimate.

    Mechanics: find maximal runs of consecutive *resolved* genes (in gene order) that sit on one
    contig for both people, pick a run with >= 2 genes, pick a breakpoint strictly inside it, and
    exchange hap1/hap2 for every gene from the breakpoint onward that lies on that same contig in
    this person. That is precisely what an assembly switch error does to a phased contig.

    Returns (new_pmap, info) or (None, None) when the pair has no eligible run."""
    resolved = [d for d in detail if d["state"] is not None]
    runs, cur = [], []
    for d in resolved:
        if cur and (d["a_contig"] == cur[-1]["a_contig"] and d["b_contig"] == cur[-1]["b_contig"]):
            cur.append(d)
        else:
            if len(cur) >= 2:
                runs.append(cur)
            cur = [d]
    if len(cur) >= 2:
        runs.append(cur)
    if not runs:
        return None, None

    run = runs[int(rng.integers(len(runs)))]
    bp = int(rng.integers(1, len(run)))            # 1..len-1 -> strictly inside the run
    target_contig = run[bp]["b_contig"]
    order = {d["gene"]: i for i, d in enumerate(detail)}
    bp_rank = order[run[bp]["gene"]]

    new = {g: {"seq": dict(v["seq"]), "contig": dict(v["contig"]), "name": dict(v["name"])}
           for g, v in pmap.items()}
    n_swapped = 0
    for gene, entry in new.items():
        if order.get(gene, -1) < bp_rank:
            continue
        if target_contig is not None and target_contig not in entry["contig"].values():
            continue
        for key in ("seq", "contig", "name"):
            d = entry[key]
            h1, h2 = d.get("hap1"), d.get("hap2")
            d.pop("hap1", None)
            d.pop("hap2", None)
            if h2 is not None:
                d["hap1"] = h2
            if h1 is not None:
                d["hap2"] = h1
        n_swapped += 1
    return new, {"breakpoint_gene": run[bp]["gene"], "contig": target_contig,
                 "n_genes_swapped": n_swapped, "run_length": len(run)}


def switch_power(pairs, gene_order, genes, seed=0):
    """Empirical sensitivity of the switch detector. `pairs`: list of (a_map, b_map). For each
    pair with >= 1 testable transition, insert one synthetic swap into person B and ask whether
    the SAME detector (`count_switches`, reused from 16) now reports >= 1 switch.

    Returns (DataFrame per pair, summary dict). Deterministic given `seed`."""
    rows = []
    for i, (a_map, b_map) in enumerate(pairs):
        detail = resolve_states(a_map, b_map, gene_order, genes)
        n_sw, _ = count_switches(detail)
        n_testable, n_resolved = testable_transitions(detail)
        row = {"pair_index": i, "n_resolved_genes": n_resolved,
               "n_testable_transitions": n_testable, "n_switches_observed": n_sw,
               "eligible": n_testable >= 1, "detected": None, "n_genes_swapped": None}
        if n_testable >= 1:
            rng = np.random.default_rng(int(seed) + i)
            sim_map, info = apply_synthetic_switch(b_map, detail, rng)
            if sim_map is not None:
                sim_detail = resolve_states(a_map, sim_map, gene_order, genes)
                sim_sw, _ = count_switches(sim_detail)
                row["detected"] = bool(sim_sw >= 1)
                row["n_genes_swapped"] = info["n_genes_swapped"]
                row["n_switches_simulated"] = sim_sw
        rows.append(row)
    df = pd.DataFrame(rows)
    elig = df[df["detected"].notna()] if not df.empty else df
    summary = {
        "n_pairs": int(len(df)),
        "n_eligible_pairs": int(len(elig)),
        "detection_rate": float(elig["detected"].mean()) if len(elig) else float("nan"),
        "total_testable_transitions": int(df["n_testable_transitions"].sum()) if len(df) else 0,
        "median_testable_transitions": (float(df["n_testable_transitions"].median())
                                        if len(df) else float("nan")),
        "median_resolved_genes": (float(df["n_resolved_genes"].median())
                                  if len(df) else float("nan")),
        "n_pairs_zero_testable": int((df["n_testable_transitions"] == 0).sum()) if len(df) else 0,
        "n_pairs_with_observed_switch": int((df["n_switches_observed"] > 0).sum()) if len(df) else 0,
    }
    return df, summary


# ===========================================================================
# STEP 7 -- orthogonal SR-vs-LR 2-field concordance (no relatives involved)
# ===========================================================================
def lr_field_class(novelty_by_hap):
    """LR field class for one person-gene, from Table 1's `novelty_depth` encoding (SCHEMA Table
    1): not novel -> `known`; depth 4 -> `f4_noncoding`; 3 -> `f3_synonymous`; 2 -> `f2_protein`;
    depth 1 or undetermined -> None, meaning EXCLUDE and count the exclusion (the first field
    itself is unresolved, so a 2-field comparison is not defined).

    A person-gene whose 2-field name is missing entirely (the depth-1 case: even field 1 is
    unresolved) is likewise excluded-and-counted by the caller.

    `novelty_by_hap`: {hap: (is_novel, novelty_depth)}. When the two haps disagree the more severe
    class wins -- the concordance of a person-gene is limited by its worse allele."""
    severity = {"f2_protein": 3, "f3_synonymous": 2, "f4_noncoding": 1, "known": 0}
    classes = []
    for is_novel, depth in novelty_by_hap.values():
        if not is_novel:
            classes.append("known")
            continue
        try:
            d = int(depth)
        except (TypeError, ValueError):
            return None
        if d == 4:
            classes.append("f4_noncoding")
        elif d == 3:
            classes.append("f3_synonymous")
        elif d == 2:
            classes.append("f2_protein")
        else:
            return None                                  # depth 1 == undetermined -> excluded
    if not classes:
        return None
    return max(classes, key=lambda c: severity[c])


def genotype_match_count(lr_names, sr_names):
    """Alleles matching at 2 fields between two unordered genotypes, as a MULTISET intersection
    size: 0, 1 or 2. Both sides must already be normalized through `vc.to_nfield(..., 2)` (07's
    normalization, mirrored) so that '02:01:01:01', 'A*02:01' and '02:01' all compare equal."""
    a, b = Counter(lr_names), Counter(sr_names)
    return int(sum((a & b).values()))


def sr_lr_concordance_rows(lr_by_person_gene, sr_by_person_gene, ancestry, novelty_by_person_gene):
    """One row per (person, gene) present on BOTH sides, with n_match in {0,1,2}, the person's
    strict ancestry, and the LR field class. Person ids are used as keys only and never returned.

    Returns (DataFrame, n_excluded_undetermined)."""
    rows, n_excluded = [], 0
    for (pid, gene), lr_names in sorted(lr_by_person_gene.items()):
        sr_names = sr_by_person_gene.get((pid, gene))
        anc = ancestry.get(pid)
        if not sr_names or anc is None:
            continue                     # nothing to compare against / no strict ancestry label
        fc = lr_field_class(novelty_by_person_gene.get((pid, gene), {}))
        if fc is None or not lr_names:
            # depth-1 / undetermined: the first field itself is unresolved, so a 2-field
            # comparison is undefined. Counted, never silently dropped.
            n_excluded += 1
            continue
        rows.append({"gene": gene, "ancestry": anc, "field_class": fc,
                     "n_lr_alleles": len(lr_names), "n_sr_alleles": len(sr_names),
                     "n_match": genotype_match_count(lr_names, sr_names)})
    return pd.DataFrame(rows), n_excluded


def aggregate_sr_lr(conc_df, min_people=MIN_ANCESTRY_PEOPLE):
    """Mean matching alleles (out of 2) and the fraction fully concordant, by gene x ancestry x
    field class. `n_person_genes` is a person-level count, so it is disclosed via safe_count in
    every written artifact."""
    if conc_df.empty:
        return conc_df
    agg = (conc_df.groupby(["gene", "ancestry", "field_class"])
           .agg(n_person_genes=("n_match", "size"),
                mean_alleles_matching=("n_match", "mean"),
                frac_both_match=("n_match", lambda s: float((s == 2).mean())),
                frac_zero_match=("n_match", lambda s: float((s == 0).mean())))
           .reset_index())
    agg["thin_n"] = agg["n_person_genes"] < min_people
    return agg


# ===========================================================================
# Figures (5 per WS2 step 8; dpi 200, Okabe-Ito colorblind-safe)
# ===========================================================================
DPI = 200


def plot_sharing_histogram(pairs_df, out_path, thresholds=SHARING_THRESHOLDS):
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    sub = pairs_df.dropna(subset=["sharing_frac"])
    groups = [("first_degree", "#0072B2"), ("replicate", "#D55E00")]
    data = {"first_degree": sub.loc[~sub["is_replicate"], "sharing_frac"].values,
            "replicate": sub.loc[sub["is_replicate"], "sharing_frac"].values}
    bins = np.linspace(0, 1, 21)
    bottom = np.zeros(len(bins) - 1)
    for label, color in groups:
        vals = data[label]
        if len(vals) == 0:
            continue
        counts, _ = np.histogram(vals, bins=bins)
        ax.bar(bins[:-1], counts, width=np.diff(bins), align="edge", bottom=bottom,
               color=color, edgecolor="white", linewidth=0.4, label=f"{label} (n={len(vals)})")
        bottom += counts
    for t in thresholds:
        ax.axvline(t, color="#555555", linestyle=":", linewidth=0.9)
        ax.text(t, ax.get_ylim()[1], f" {t}", fontsize=7, va="top", color="#555555")
    ax.set_xlabel("Fraction of comparable genes where the pair shares an exact CDS sequence")
    ax.set_ylabel("Number of pairs")
    ax.set_title("Step 2: HLA sharing across ALL kinship-selected relative pairs\n"
                 "(no pair was selected on this quantity;\ndotted lines = the reported threshold "
                 "sweep)", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    vc.savefig(fig, out_path, dpi=DPI)


def plot_decomposition_by_gene(decomp_df, out_path, top_n=30):
    if decomp_df.empty:
        return
    df = decomp_df.copy()
    df["n_plot"] = df["n_comparable"]
    df = df.sort_values("n_plot", ascending=False).head(top_n)
    df = df.sort_values("discordance_rate", ascending=False)
    cats = [c for c in DISCORDANCE_CATEGORIES if c not in ("concordant", "missing_both")]
    frac = df[cats].div(df["n_comparable"].replace(0, np.nan), axis=0).fillna(0.0) * 100
    x = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(max(8, 0.42 * len(df)), 5.2))
    bottom = np.zeros(len(df))
    for cat in cats:
        vals = frac[cat].values
        ax.bar(x, vals, bottom=bottom, color=CATEGORY_COLORS[cat], width=0.78,
               label=cat, edgecolor="white", linewidth=0.3)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([vc.gene_display(g) for g in df["gene"]], rotation=90, fontsize=7)
    ax.set_ylabel("% of comparable gene-comparisons")
    ax.set_title("Step 3: what the discordance actually IS, per gene\n(mutually exclusive "
                 "categories;\ndenominator INCLUDES genes missing in one person)", fontsize=10)
    ax.legend(fontsize=7.5, ncol=3, frameon=False, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    vc.savefig(fig, out_path, dpi=DPI)


def plot_homozygosity_obs_exp(hom_df, out_path):
    if hom_df.empty:
        return
    measures = ["LR_CDS_exact", "LR_2field", "SR_2field"]
    markers = {"LR_CDS_exact": "o", "LR_2field": "s", "SR_2field": "^"}
    genes = sorted(hom_df["gene"].unique())
    ancestries = [a for a in vc.ANCESTRY_ORDER if a in set(hom_df["ancestry"])]
    fig, axes = plt.subplots(1, max(1, len(ancestries)), figsize=(
        max(7.5, 3.2 * max(1, len(ancestries))), 4.8), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, anc in zip(axes, ancestries or [None]):
        sub = hom_df[hom_df["ancestry"] == anc]
        for mi, measure in enumerate(measures):
            s = sub[sub["measure"] == measure]
            if s.empty:
                continue
            xs = [genes.index(g) + (mi - 1) * 0.22 for g in s["gene"]]
            ax.errorbar(xs, s["obs_exp_ratio"],
                        yerr=[np.clip(s["obs_exp_ratio"] - s["ci_lo"], 0, None),
                              np.clip(s["ci_hi"] - s["obs_exp_ratio"], 0, None)],
                        fmt=markers[measure], markersize=4.5, capsize=2, elinewidth=0.8,
                        color=vc.ANCESTRY_COLORS.get(anc, "#333333"),
                        markerfacecolor="white" if measure != "LR_CDS_exact" else None,
                        label=measure if ax is axes[0] else None)
            thin = s[s["thin_n"]]
            if not thin.empty:
                ax.scatter([genes.index(g) + (mi - 1) * 0.22 for g in thin["gene"]],
                           thin["obs_exp_ratio"], marker="x", s=28, color="#000000",
                           label="thin n (<50 people)" if ax is axes[0] and mi == 0 else None)
        ax.axhline(1.0, color="#555555", linewidth=0.9, linestyle="--")
        ax.set_xticks(range(len(genes)))
        ax.set_xticklabels([vc.gene_display(g) for g in genes], rotation=90, fontsize=7)
        ax.set_title(anc if anc else "(no strict-ancestry group)", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("observed / expected homozygosity")
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(handles, labels, fontsize=7, frameon=False)
    # Line breaks are explicit: savefig(bbox_inches='tight') grows the canvas to fit an overlong
    # single-line suptitle, which leaves the panels marooned in whitespace.
    fig.suptitle("Step 5: excess homozygosity is the allele-dropout signature\n"
                 "1.0 = Hardy-Weinberg expectation from that group's own frequencies;\n"
                 "LR sitting above SR on the SAME people is what collapse looks like",
                 fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    vc.savefig(fig, out_path, dpi=DPI)


def plot_switch_power(power_df, power_summary, out_path):
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0))
    ax = axes[0]
    if not power_df.empty:
        vals = power_df["n_resolved_genes"].values
        ax.hist(vals, bins=np.arange(-0.5, max(vals.max(), 1) + 1.5, 1), color="#0072B2",
                edgecolor="white", linewidth=0.4)
    ax.set_xlabel("Resolved genes per pair")
    ax.set_ylabel("Number of pairs")
    ax.set_title("How many genes carry phase info", fontsize=9)

    ax = axes[1]
    if not power_df.empty:
        vals = power_df["n_testable_transitions"].values
        ax.hist(vals, bins=np.arange(-0.5, max(vals.max(), 1) + 1.5, 1), color="#CC79A7",
                edgecolor="white", linewidth=0.4)
    ax.set_xlabel("TESTABLE same-contig transitions per pair")
    ax.set_ylabel("Number of pairs")
    ax.set_title("The denominator v1 never reported", fontsize=9)

    ax = axes[2]
    rate = power_summary.get("detection_rate", float("nan"))
    n_elig = power_summary.get("n_eligible_pairs", 0)
    n_zero = power_summary.get("n_pairs_zero_testable", 0)
    ax.bar(["detected", "missed", "untestable\n(0 transitions)"],
           [rate * n_elig if not pd.isna(rate) else 0,
            (1 - rate) * n_elig if not pd.isna(rate) else 0, n_zero],
           color=["#009E73", "#D55E00", "#999999"])
    ax.set_ylabel("Number of pairs")
    ax.set_title(f"Empirical power: {0 if pd.isna(rate) else rate*100:.1f}% of synthetic\n"
                 f"switches detected ({n_elig} eligible pairs)", fontsize=9)
    for a in axes:
        a.spines[["top", "right"]].set_visible(False)
        a.yaxis.set_major_locator(MaxNLocator(integer=True))   # these are pair COUNTS
    fig.suptitle("Step 6: the switch test's own sensitivity, measured rather than assumed",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    vc.savefig(fig, out_path, dpi=DPI)


def plot_sr_lr_concordance(sr_lr_agg, out_path):
    if sr_lr_agg.empty:
        return
    # Pool ancestries back together for the figure by re-weighting on person-gene counts.
    # (Plain arithmetic rather than groupby.apply: `include_groups` semantics differ across the
    # pandas versions present in this project's pixi envs.)
    tmp = sr_lr_agg.copy()
    tmp["_num"] = tmp["mean_alleles_matching"] * tmp["n_person_genes"]
    grp = tmp.groupby(["gene", "field_class"])[["_num", "n_person_genes"]].sum().reset_index()
    grp["mean_alleles_matching"] = grp["_num"] / grp["n_person_genes"].replace(0, np.nan)
    df = grp[["gene", "field_class", "mean_alleles_matching", "n_person_genes"]]
    genes = sorted(df["gene"].unique())
    classes = [c for c in FIELD_CLASS_ORDER if c in set(df["field_class"])]
    width = 0.8 / max(len(classes), 1)
    fig, ax = plt.subplots(figsize=(max(7.5, 1.1 * len(genes)), 4.4))
    for i, fc in enumerate(classes):
        sub = df[df["field_class"] == fc].set_index("gene").reindex(genes)
        xs = np.arange(len(genes)) + (i - (len(classes) - 1) / 2) * width
        ax.bar(xs, sub["mean_alleles_matching"].fillna(0), width=width,
               color=FIELD_CLASS_COLORS[fc], label=fc, edgecolor="white", linewidth=0.3)
        for x, v, n in zip(xs, sub["mean_alleles_matching"], sub["n_person_genes"]):
            if pd.isna(v):
                continue
            ax.text(x, v + 0.03, safe_count(n), ha="center", fontsize=5.5, rotation=90)
    ax.axhline(2.0, color="#555555", linestyle="--", linewidth=0.9)
    ax.set_xticks(np.arange(len(genes)))
    ax.set_xticklabels([vc.gene_display(g) for g in genes], rotation=0, fontsize=8)
    ax.set_ylim(0, 2.35)
    ax.set_ylabel("Mean alleles matching short-read calls (of 2)")
    ax.set_title("Step 7: orthogonal check -- long-read vs AoU short-read calls,\n"
                 "same people, 2-field. Prediction: noncoding-only novelty must NOT lower this;\n"
                 "protein-altering novelty should. Bar labels = person-gene counts.", fontsize=9)
    ax.legend(fontsize=8, frameon=False, ncol=len(classes))
    ax.spines[["top", "right"]].set_visible(False)
    vc.savefig(fig, out_path, dpi=DPI)


# ===========================================================================
# Report
# ===========================================================================
def build_report(ctx):
    """Plain-language report: every number is stated together with what it would mean if it were
    different, so a non-specialist can read the evidence rather than the jargon."""
    L = ["# Are the long-read HLA calls actually right? (QC v2)\n",
         "This is the load-bearing QC for the paper's central claim. Version 1 of this analysis "
         "(scripts 16/17) answered a question it had partly assumed the answer to; this version "
         "fixes each of those problems. Everything below is computed on the **exact assembled DNA "
         "sequence** of each gene, not on allele names, except where it says otherwise.\n"]

    L.append("## 1. Which pairs of people we looked at, and why that choice is now fair\n")
    L.append(f"We took **every** pair of people in the long-read cohort whom All of Us's own "
             f"genome-wide relatedness table calls first-degree relatives or closer "
             f"(kinship >= {FIRST_DEGREE_MIN_KIN}). That is "
             f"**{safe_count(ctx['n_pairs'])}** pairs, of which "
             f"**{safe_count(ctx['n_replicate_pairs'])}** are duplicates or identical twins "
             f"(kinship > {DUPLICATE_MIN_KIN}).\n")
    L.append("Version 1 instead kept only pairs that already agreed at >=80% of their genes, and "
             "then reported how often those same pairs disagreed. That is circular: the pairs "
             "were chosen using the very number being reported. Nothing here is chosen using HLA "
             "data at all.\n")

    L.append("## 2. How much HLA the relatives share (a distribution, not a filter)\n")
    L.append("Two first-degree relatives do not have to share HLA. A parent and child always "
             "share one copy, but two siblings can easily inherit different copies from each "
             "parent and share none. So the honest picture is the whole spread, not a cutoff:\n")
    L.append(f"- Median fraction of genes where the pair shares an exact sequence: "
             f"**{ctx['sharing_median']:.3f}**\n")
    L.append(f"- Pairs sharing at almost every gene (>=0.95): "
             f"**{safe_count(ctx['n_sharing_high'])}**; pairs sharing at almost none (<=0.1): "
             f"**{safe_count(ctx['n_sharing_low'])}**. Two modes here is the expected biology, "
             f"not a problem.\n")
    L.append("\nBecause any cutoff is a judgement call, the headline disagreement rate is "
             "reported at several cutoffs. If the number moves a lot across this table, the "
             "cutoff was doing the work rather than the data:\n")
    L.append(_md_table(ctx["sweep_display"]) + "\n")
    if ctx["yob_available"]:
        L.append("\nUsing year of birth, pairs at least "
                 f"{AGE_GAP_PARENT_CHILD} years apart (so probably parent-child, who must share "
                 "one copy at every gene) are reported separately:\n")
        L.append(_md_table(ctx["age_gap_display"]) + "\n")
    else:
        L.append("\n**Age-gap stratification was SKIPPED**: no `--yob-tsv` was supplied, so we "
                 "cannot separate probable parent-child pairs (which must share a copy at every "
                 "gene) from siblings (which need not). Supply `--yob-tsv` with columns "
                 "`person_id`, `year_of_birth` to enable it.\n")

    L.append("\n## 3. What the disagreements actually are\n")
    L.append("Every single gene comparison is placed in exactly one box. Crucially, a gene that "
             "one person has and the other does not is now COUNTED (as `missing_one`) instead of "
             "being quietly deleted from the denominator, which is what hid allele dropout and "
             "assembly gaps in version 1.\n")
    L.append("- `concordant` - both people have at least one byte-identical copy of this gene.\n"
             "- `missing_one` / `missing_both` - the gene was not assembled in one / either "
             "person. A coverage result, not an accuracy result.\n"
             f"- `point_diff` - up to {POINT_DIFF_MAX_BASES} differing bases in one place. This is "
             "what ordinary long-read consensus error looks like.\n"
             "- `dropout_candidate` - a real difference, and at least one of the two people looks "
             "homozygous here, which is what a collapsed (lost) second copy looks like.\n"
             "- `fragmentation_candidate` - a real difference at a gene that sits alone on its own "
             "assembly fragment, where there is nothing to anchor it.\n"
             "- `other_multiblock` - a real, scattered difference with none of the above excuses. "
             "This is the bucket that would genuinely worry us.\n")
    L.append("\n" + _md_table(ctx["class_display"]) + "\n")
    L.append(f"\nAcross all genes, **{ctx['overall_rate']*100:.2f}%** of comparable comparisons "
             f"disagree (95% CI {ctx['overall_lo']*100:.2f}-{ctx['overall_hi']*100:.2f}%), and "
             f"**{ctx['pct_point']:.1f}%** of those disagreements are just 1-"
             f"{POINT_DIFF_MAX_BASES} bases.\n")

    L.append("\n## 4. The cleanest possible check: the same person sequenced twice\n")
    if ctx["replicate_qv"] is None:
        L.append("No duplicate/identical-twin pairs were available in this run, so no direct "
                 "per-base quality estimate could be made.\n")
    else:
        qv = ctx["replicate_qv"]
        L.append(f"Duplicate/identical-twin pairs should be byte-identical. Comparing both "
                 f"haplotypes gene by gene gives **{ctx['replicate_n_diff']} differing bases in "
                 f"{ctx['replicate_n_aligned']} compared bases**, i.e. a per-base quality of "
                 f"**QV {'inf' if qv == float('inf') else f'{qv:.1f}'}** "
                 f"(higher is better; 40-50 is the normal range for this kind of assembly). "
                 f"This number involves no relatives, no inheritance assumptions and no allele "
                 f"names - it is the most direct accuracy statement in this report.\n")

    L.append("\n## 5. Are we silently losing one copy of a gene? (homozygosity check)\n")
    L.append("If an assembly ever collapses two different copies of a gene into one, people will "
             "look homozygous (two identical copies) more often than genetics allows. We compare "
             "how often people ARE homozygous with how often they SHOULD be, given the allele "
             "frequencies of their own ancestry group. A ratio near 1.0 is healthy; clearly above "
             "1.0 means copies are going missing.\n")
    L.append("We do this three ways on the SAME people - long-read exact sequence, long-read "
             "allele names, and All of Us's independent short-read allele names - so the "
             "comparison is like for like. Only people whose ancestry is unambiguous "
             f"(own admixture proportion >= {ctx['strict_min']}) are used, because the expectation "
             "depends on which population's frequencies apply.\n")
    L.append("Expected homozygosity uses the finite-sample-corrected formula "
             "`sum x(x-1)/(n(n-1))`; the naive version is biased upward in small groups, which "
             "would have hidden exactly this signal where it matters most.\n")
    if ctx["hom_display"] is not None:
        L.append("\n" + _md_table(ctx["hom_display"]) + "\n")
    L.append(f"\nGroups with fewer than {MIN_ANCESTRY_PEOPLE} people are shown but flagged "
             f"`thin_n`; treat them as suggestive only.\n")
    if ctx["frag_display"] is not None:
        L.append("\nSplit by how fragmented the person's assembly is (how many separate contigs "
                 "carry their 8 classical genes):\n")
        L.append(_md_table(ctx["frag_display"]) + "\n")

    L.append("\n## 6. Does the phase-switch test have any power? (the missing denominator)\n")
    ps = ctx["power_summary"]
    L.append("Version 1's headline was '0 of 545 pairs show a phase switch'. That is only "
             "meaningful if a switch COULD have been seen. A switch is only visible between two "
             "neighbouring genes that (a) both carry unambiguous shared-copy information and "
             "(b) sit on the same assembly fragment in both people. We now count those "
             "opportunities:\n")
    L.append(f"- Median genes carrying phase information per pair: "
             f"**{ps['median_resolved_genes']}**\n")
    L.append(f"- Median *testable* neighbouring-gene transitions per pair: "
             f"**{ps['median_testable_transitions']}**; total across all pairs: "
             f"**{ps['total_testable_transitions']}**\n")
    L.append(f"- Pairs where NO switch could ever have been detected (zero testable transitions): "
             f"**{safe_count(ps['n_pairs_zero_testable'])}** of {safe_count(ps['n_pairs'])}\n")
    L.append(f"- Pairs where a switch WAS observed: "
             f"**{safe_count(ps['n_pairs_with_observed_switch'])}**\n")
    rate = ps["detection_rate"]
    L.append(f"\nThen we measured the test instead of trusting it: in each eligible pair we "
             f"deliberately broke the phasing of one person at a random point inside one contig "
             f"(a synthetic switch error, with a fixed random seed so this is reproducible) and "
             f"asked the same detector to find it. It found "
             f"**{'NA' if pd.isna(rate) else f'{rate*100:.1f}%'}** of them across "
             f"{safe_count(ps['n_eligible_pairs'])} eligible pairs.\n")
    L.append("Read this as the ceiling on what '0 switches' can mean: if the detection rate is "
             "low, a clean result is mostly a statement about power, not about the assemblies.\n")

    L.append("\n## 7. An independent opinion: All of Us's own short-read calls\n")
    L.append("Steps 1-6 all use relatives. This step does not. For every person who has both "
             "kinds of data we ask how many of their two alleles per gene agree with All of Us's "
             "independent short-read HLA calls, at two-field resolution (both sides normalized the "
             "same way).\n")
    L.append("The sharp prediction: a long-read call that is novel only OUTSIDE the coding "
             "sequence (`f4_noncoding`) is invisible at two fields and must NOT agree less often. "
             "A novel call that changes the protein (`f2_protein`) is a two-field difference by "
             "construction and SHOULD agree less often. If `f4` novelty also lowered agreement, "
             "our novelty calls would look like noise.\n")
    if ctx["sr_lr_display"] is not None:
        L.append("\n" + _md_table(ctx["sr_lr_display"]) + "\n")
    L.append(f"\nPerson-genes excluded because the long-read call's very first field was "
             f"unresolved (`novelty_depth` 1 / undetermined, where a two-field comparison is not "
             f"defined): **{safe_count(ctx['sr_lr_excluded'])}**.\n")

    L.append("\n## What must still be checked before quoting any of this\n")
    for item in ctx["caveats"]:
        L.append(f"- {item}\n")
    L.append("\n## Disclosure discipline\n")
    L.append(f"No person identifiers appear in this directory. Any person or pair count between 1 "
             f"and {SMALL_COUNT_CEILING - 1} is written as `<{SMALL_COUNT_CEILING}`. Person-level "
             f"rows, if written at all, went only to the separate local directory.\n")
    return L


CAVEATS = [
    "Kinship alone cannot tell a parent-child pair from a full sibling pair (the All of Us "
    "relatedness table has no IBD0 column - see script 11). Without `--yob-tsv` the "
    "parent-child-specific 'must share a copy at every gene' logic cannot be applied to any pair.",
    "`dropout_candidate` is a *shape*, not a proof: a genuinely homozygous person who differs "
    "from their relative lands in it too. Its value is as an upper bound, read together with the "
    "homozygosity obs/exp ratio in step 5.",
    "`fragmentation_candidate` uses the contig labels in Table 1; a gene alone on a contig is "
    "suspicious but not wrong by itself.",
    "The synthetic-switch power estimate measures sensitivity to ONE clean hap-label swap inside "
    "one contig. It does not estimate sensitivity to short switch segments, to errors that also "
    "corrupt the sequence, or to switches that coincide with a contig break (which are "
    "undetectable by construction).",
    "Different-length sequences fall back to difflib alignment (inherited from script 17), which "
    "can misalign in repetitive stretches; equal-length comparisons are position-anchored and "
    "exact.",
    "copy_index > 1 rows (segmental duplication / DRB copy number) are excluded everywhere here.",
]


# ===========================================================================
# Orchestration (thin I/O wrappers only)
# ===========================================================================
def build_parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=os.path.join(vc.DEFAULT_OUTROOT, "hla_calls_rich.tsv"))
    ap.add_argument("--cohort-membership",
                    default=os.path.join(vc.DEFAULT_OUTROOT, "cohort_membership.tsv"))
    ap.add_argument("--relatedness-table", default=DEFAULT_RELATEDNESS)
    ap.add_argument("--sr-genotypes", default=DEFAULT_SR_GENOTYPES,
                    help="AoU-native short-read HLA calls (research_id + <gene>_1/<gene>_2), "
                         "read from the auto-mounted controlled-tier dataset.")
    ap.add_argument("--people-root", default=DEFAULT_PEOPLE_ROOT,
                    help="Root holding <person_id>/immuannot_output/hap{1,2}/cds.fa.gz.")
    ap.add_argument("--yob-tsv", default=None,
                    help="Optional TSV with person_id, year_of_birth. Absent -> the age-gap "
                         "stratification is skipped and the report says so.")
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR,
                    help="Aggregate-only outputs. No person ids, counts 1-19 written as '<20'.")
    ap.add_argument("--local-dir", default=DEFAULT_LOCAL_DIR,
                    help="The ONLY place anything person-level is written (pair rows). Set to "
                         "'' to write nothing person-level at all.")
    ap.add_argument("--limit-pairs", type=int, default=None,
                    help="Pilot mode: only analyze the first N pairs (highest kinship first).")
    ap.add_argument("--threads", type=int, default=8,
                    help="Threads used to read per-person cds.fa.gz files (I/O bound). The "
                         "Workbench VM has 4 vCPUs / 31GB, so 8 keeps the gzip reads overlapped "
                         "without thrashing; raise it only on a bigger machine.")
    ap.add_argument("--resolution", type=int, default=2, choices=[2, 3, 4],
                    help="Field depth for the NAME-based comparisons only (SR-vs-LR concordance "
                         "and LR-2field homozygosity). Sequence comparisons ignore it.")
    ap.add_argument("--strict-ancestry-min", type=float, default=0.9,
                    help="Minimum own-admixture proportion for a person to count under their "
                         "ancestry label.")
    ap.add_argument("--bootstrap", type=int, default=200,
                    help="Bootstrap replicates (over people) for homozygosity obs/exp CIs.")
    ap.add_argument("--seed", type=int, default=20260916,
                    help="Seed for the switch simulation and the bootstrap (deterministic).")
    ap.add_argument("--min-ancestry-n", type=int, default=MIN_ANCESTRY_PEOPLE)
    ap.add_argument("--limit-hom-people", type=int, default=None,
                    help="Pilot mode for steps 5/7: cap how many people (of the whole LR-and-SR "
                         "overlap) are streamed for the homozygosity and SR-concordance steps.")
    return ap


def _load_maps(people_root, person_ids, genes, contig_lookup, name_lookup, threads,
               digest=False):
    """Parallel (I/O-bound) read of every needed person's cds.fa.gz, then the pure reshape.
    `digest=True` keeps only sequence hashes -- used by the cohort-wide homozygosity step, which
    never needs to align anything, only to ask whether two copies are the same."""
    maps = {}

    def one(pid):
        cds = load_person_cds(people_root, pid, genes)
        return pid, build_person_map(cds, contig_lookup, name_lookup, pid, genes, digest=digest)

    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        for pid, pmap in ex.map(one, person_ids):
            maps[pid] = pmap
    return maps


def main_from_args(args):
    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading Table 1...", file=sys.stderr)
    table1 = load_table1_lean(args.table1)
    contig_lookup, name_lookup, gene_class, novelty = build_table1_lookups(
        table1, args.resolution)
    gene_order = derive_canonical_order(table1)
    all_genes = sorted(table1["gene_bare"].dropna().unique())
    sharing_genes = [g for g in all_genes
                     if gene_class.get(g) in SHARING_GENE_CLASSES]
    classical = [g for g in CLASSICAL if g in set(all_genes)] or CLASSICAL

    print("Loading cohort membership + relatedness...", file=sys.stderr)
    cohort = vc.load_cohort_membership(args.cohort_membership)
    rel = mod11.load_relatedness(args.relatedness_table)
    pairs = build_pair_universe(rel, cohort)
    if args.limit_pairs:
        pairs = pairs.head(args.limit_pairs).copy()
    print(f"  {len(pairs)} kinship-selected pairs, both in LR cohort.", file=sys.stderr)
    strict_anc = strict_ancestry_map(cohort, args.strict_ancestry_min)
    yob = load_year_of_birth(args.yob_tsv)

    people = sorted(set(pairs["person_i"]) | set(pairs["person_j"]))
    print(f"Reading cds.fa.gz for {len(people)} people ({args.threads} threads)...",
          file=sys.stderr)
    maps = _load_maps(args.people_root, people, all_genes, contig_lookup, name_lookup,
                      args.threads)

    # ---- steps 2/3/4/6 per pair -------------------------------------------------
    pair_rows, decomp_records, replicate_diffs, power_pairs = [], [], [], []
    for row in pairs.itertuples(index=False):
        a_map, b_map = maps.get(row.person_i, {}), maps.get(row.person_j, {})
        n_shared, n_cmp, frac = pair_sharing(a_map, b_map, sharing_genes)
        recs = decompose_pair(a_map, b_map, all_genes, classical)
        n_comparable = sum(1 for r in recs if r["category"] != "missing_both")
        n_discordant = sum(1 for r in recs if r["category"] not in ("concordant", "missing_both"))
        decomp_records.extend(recs)
        pair_rows.append({
            "person_i": row.person_i, "person_j": row.person_j, "kin": float(row.kin),
            "kin_degree": row.kin_degree, "is_replicate": bool(row.is_replicate),
            "n_shared_genes": n_shared, "n_sharing_genes_compared": n_cmp,
            "sharing_frac": frac, "n_gene_comparisons": n_comparable,
            "n_discordant": n_discordant,
            "age_gap_class": age_gap_class(row.person_i, row.person_j, yob),
        })
        if row.is_replicate:
            for gene in all_genes:
                d = replicate_gene_diff(a_map, b_map, gene)
                if d:
                    replicate_diffs.append({"gene": gene, **d})
        else:
            power_pairs.append((a_map, b_map))

    pairs_out = pd.DataFrame(pair_rows)
    decomp_df = aggregate_decomposition(decomp_records, gene_class)
    class_df = aggregate_by_gene_class(decomp_df)
    sweep_df = sharing_threshold_sweep(pair_rows)
    power_df, power_summary = switch_power(power_pairs, gene_order, all_genes, seed=args.seed)

    # ---- step 5: homozygosity on people present in BOTH cohorts ---------------
    print("Loading SR genotypes...", file=sys.stderr)
    sr_long = vc.load_sr_genotypes(args.sr_genotypes)
    # People for steps 5 and 7: everyone with BOTH kinds of data and an unambiguous ancestry --
    # not just the relative pairs. These two steps are cohort-wide statements, and restricting
    # them to relatives would re-introduce a selection the audit told us to remove.
    both_cohort = set(cohort.loc[cohort["in_lr"].astype(bool) & cohort["in_sr"].astype(bool),
                                  "person_id"])
    hom_people = sorted(both_cohort & set(strict_anc))
    if args.limit_hom_people:
        hom_people = hom_people[:args.limit_hom_people]
    if not hom_people:
        print("  WARNING: no people are in both cohorts with a strict ancestry label; steps 5 and "
              "7 will be empty.", file=sys.stderr)
    # Filter the (potentially ~500k x 16) SR table down to those people BEFORE the per-allele
    # normalization, so the expensive string work runs on the rows we keep, not on all of them.
    sr_long = sr_long[sr_long["person_id"].isin(set(hom_people))].copy()
    sr_long["allele_n"] = sr_long["allele"].map(lambda a: vc.to_nfield(a, args.resolution))
    sr_by_pg = defaultdict(list)
    for r in sr_long.dropna(subset=["allele_n"]).itertuples(index=False):
        sr_by_pg[(r.person_id, r.gene_bare)].append(r.allele_n)
    del sr_long

    calls_by_group = defaultdict(list)
    frag_calls = defaultdict(list)
    lr_by_pg, nov_by_pg = {}, {}
    # Streamed in chunks, keeping only sequence DIGESTS: at 12k people this step would otherwise
    # hold every classical-gene CDS in memory at once, for no benefit (it only tests equality).
    print(f"Homozygosity + SR concordance over {len(hom_people)} people (streamed, hashed)...",
          file=sys.stderr)
    CHUNK = 500
    for start in range(0, len(hom_people), CHUNK):
        chunk = hom_people[start:start + CHUNK]
        hom_maps = _load_maps(args.people_root, chunk, classical, contig_lookup, name_lookup,
                              args.threads, digest=True)
        for pid in chunk:
            anc = strict_anc[pid]
            pmap = hom_maps.get(pid, {})
            frag = fragmentation_class(pmap, classical)
            for gene in classical:
                seqs = (pmap.get(gene, {}).get("seq") or {})
                names = (pmap.get(gene, {}).get("name") or {})
                sr_names = sr_by_pg.get((pid, gene), [])
                if len(seqs) == 2:
                    pair_seq = tuple(sorted(seqs.values()))
                    calls_by_group[("LR_CDS_exact", gene, anc)].append(pair_seq)
                    if frag is not None:
                        frag_calls[(gene, frag)].append(pair_seq)
                if len(names) == 2:
                    calls_by_group[("LR_2field", gene, anc)].append(tuple(sorted(names.values())))
                if len(sr_names) == 2:
                    calls_by_group[("SR_2field", gene, anc)].append(tuple(sorted(sr_names)))
                if gene in pmap:
                    # Registered even when the 2-field name is MISSING (a depth-1 "undetermined"
                    # call has a sequence but no resolvable first field): step 7 must count that
                    # exclusion rather than never seeing the person-gene at all.
                    haps = sorted(set(seqs) | set(names))
                    lr_by_pg[(pid, gene)] = list(names.values())
                    nov_by_pg[(pid, gene)] = {
                        hap: novelty.get((pid, gene, hap), (False, None)) for hap in haps}
        del hom_maps
    hom_df = homozygosity_obs_exp(calls_by_group, n_boot=args.bootstrap, seed=args.seed,
                                  min_people=args.min_ancestry_n)
    frag_df = homozygosity_by_fragmentation(frag_calls)

    # ---- step 7: SR-vs-LR concordance ----------------------------------------
    sr_pairs_by_pg = {k: v for k, v in sr_by_pg.items() if len(v) >= 1}
    conc_df, n_excluded = sr_lr_concordance_rows(lr_by_pg, sr_pairs_by_pg, strict_anc, nov_by_pg)
    sr_lr_agg = aggregate_sr_lr(conc_df, args.min_ancestry_n)

    # ---- figures -------------------------------------------------------------
    print("Writing figures...", file=sys.stderr)
    plot_sharing_histogram(pairs_out, os.path.join(args.out_dir, "sharing_histogram.png"))
    plot_decomposition_by_gene(decomp_df,
                               os.path.join(args.out_dir, "decomposition_by_gene.png"))
    plot_homozygosity_obs_exp(hom_df, os.path.join(args.out_dir, "homozygosity_obs_exp.png"))
    plot_switch_power(power_df, power_summary, os.path.join(args.out_dir, "switch_power.png"))
    plot_sr_lr_concordance(sr_lr_agg, os.path.join(args.out_dir, "sr_lr_concordance.png"))

    # ---- aggregate tables ----------------------------------------------------
    decomp_df.to_csv(os.path.join(args.out_dir, "discordance_decomposition.tsv"),
                     sep="\t", index=False)
    sweep_df.to_csv(os.path.join(args.out_dir, "sharing_threshold_sweep.tsv"),
                    sep="\t", index=False)
    hom_out = hom_df.copy()
    if not hom_out.empty:
        hom_out["n_people"] = hom_out["n_people"].map(safe_count)
    hom_out.to_csv(os.path.join(args.out_dir, "homozygosity_obs_exp.tsv"), sep="\t", index=False)
    frag_out = frag_df.copy()
    if not frag_out.empty:
        frag_out["n_people"] = frag_out["n_people"].map(safe_count)
    frag_out.to_csv(os.path.join(args.out_dir, "homozygosity_by_fragmentation.tsv"),
                    sep="\t", index=False)
    sr_out = sr_lr_agg.copy()
    if not sr_out.empty:
        sr_out["n_person_genes"] = sr_out["n_person_genes"].map(safe_count)
    sr_out.to_csv(os.path.join(args.out_dir, "sr_lr_concordance.tsv"), sep="\t", index=False)
    power_out = power_df.drop(columns=["pair_index"], errors="ignore")
    power_out.to_csv(os.path.join(args.out_dir, "switch_power_per_pair.tsv"),
                     sep="\t", index=False)

    # ---- report + summary.json ----------------------------------------------
    comparable = int(class_df["n_comparable"].sum()) if not class_df.empty else 0
    discordant = int(class_df["n_discordant"].sum()) if not class_df.empty else 0
    o_p, o_lo, o_hi = _wilson(discordant, comparable) if comparable else (float("nan"),) * 3
    n_point = sum(1 for r in decomp_records if r["category"] == "point_diff")
    qv, err, n_diff, n_aln = (replicate_qv(replicate_diffs) if replicate_diffs
                              else (None, None, 0, 0))

    sweep_display = sweep_df.copy()
    sweep_display["n_pairs"] = sweep_display["n_pairs"].map(safe_count)
    sweep_display["discordance_rate_%"] = (sweep_display["discordance_rate"] * 100).round(3)
    sweep_display = sweep_display[["threshold", "n_pairs", "n_gene_comparisons",
                                   "n_discordant", "discordance_rate_%"]]

    class_display = pd.DataFrame()
    if not class_df.empty:
        class_display = class_df.copy()
        class_display["discordance_rate_%"] = (class_display["discordance_rate"] * 100).round(3)
        class_display = class_display[["gene_class", "n_total", "n_comparable", "n_discordant",
                                       "discordance_rate_%"] + DISCORDANCE_CATEGORIES]

    age_display = pd.DataFrame()
    if yob and not pairs_out.empty:
        age_display = (pairs_out.groupby("age_gap_class")
                       .agg(n_pairs=("kin", "size"),
                            median_sharing=("sharing_frac", "median"),
                            n_gene_comparisons=("n_gene_comparisons", "sum"),
                            n_discordant=("n_discordant", "sum")).reset_index())
        age_display["discordance_rate_%"] = (
            100 * age_display["n_discordant"] / age_display["n_gene_comparisons"]).round(3)
        age_display["n_pairs"] = age_display["n_pairs"].map(safe_count)
        age_display["median_sharing"] = age_display["median_sharing"].round(3)

    hom_display = None
    if not hom_df.empty:
        hom_display = hom_df.copy()
        for c in ("obs_hom", "exp_hom", "obs_exp_ratio", "ci_lo", "ci_hi"):
            hom_display[c] = hom_display[c].round(3)
        hom_display["n_people"] = hom_display["n_people"].map(safe_count)
        hom_display = hom_display[["measure", "gene", "ancestry", "n_people", "obs_hom",
                                   "exp_hom", "obs_exp_ratio", "ci_lo", "ci_hi", "thin_n"]]
    frag_display = None
    if not frag_df.empty:
        frag_display = frag_df.copy()
        frag_display["obs_hom"] = frag_display["obs_hom"].round(3)
        frag_display["n_people"] = frag_display["n_people"].map(safe_count)
    sr_lr_display = None
    if not sr_lr_agg.empty:
        sr_lr_display = sr_lr_agg.copy()
        for c in ("mean_alleles_matching", "frac_both_match", "frac_zero_match"):
            sr_lr_display[c] = sr_lr_display[c].round(3)
        sr_lr_display["n_person_genes"] = sr_lr_display["n_person_genes"].map(safe_count)

    sharing_vals = pairs_out["sharing_frac"].dropna() if not pairs_out.empty else pd.Series(
        dtype=float)
    ctx = {
        "n_pairs": len(pairs_out),
        "n_replicate_pairs": int(pairs_out["is_replicate"].sum()) if not pairs_out.empty else 0,
        "sharing_median": float(sharing_vals.median()) if len(sharing_vals) else float("nan"),
        "n_sharing_high": int((sharing_vals >= 0.95).sum()),
        "n_sharing_low": int((sharing_vals <= 0.10).sum()),
        "sweep_display": sweep_display,
        "yob_available": bool(yob),
        "age_gap_display": age_display,
        "class_display": class_display,
        "overall_rate": o_p, "overall_lo": o_lo, "overall_hi": o_hi,
        "pct_point": (100.0 * n_point / discordant) if discordant else float("nan"),
        "replicate_qv": qv, "replicate_n_diff": n_diff, "replicate_n_aligned": n_aln,
        "strict_min": args.strict_ancestry_min,
        "hom_display": hom_display, "frag_display": frag_display,
        "power_summary": power_summary,
        "sr_lr_display": sr_lr_display, "sr_lr_excluded": n_excluded,
        "caveats": CAVEATS,
    }
    vc.write_report(os.path.join(args.out_dir, "qc_relatives_v2_report.md"), build_report(ctx))

    summary = {
        "n_pairs": safe_count(len(pairs_out)),
        "n_replicate_pairs": safe_count(ctx["n_replicate_pairs"]),
        "sharing_frac_median": ctx["sharing_median"],
        "overall_discordance_rate": o_p,
        "overall_discordance_ci": [o_lo, o_hi],
        "n_gene_comparisons_comparable": comparable,
        "n_discordant": discordant,
        "category_counts": dict(Counter(r["category"] for r in decomp_records)),
        "threshold_sweep": sweep_df.assign(
            n_pairs=sweep_df["n_pairs"].map(safe_count)).to_dict(orient="records"),
        "replicate_qv": (None if qv is None else (None if qv == float("inf") else qv)),
        "replicate_zero_differences": bool(qv == float("inf")) if qv is not None else None,
        "replicate_n_diff_bases": n_diff, "replicate_n_aligned_bases": n_aln,
        "switch_power": power_summary,
        "homozygosity_n_groups": int(len(hom_df)),
        "homozygosity_median_ratio_by_measure": (
            hom_df.groupby("measure")["obs_exp_ratio"].median().to_dict()
            if not hom_df.empty else {}),
        "sr_lr_mean_alleles_matching_by_field_class": (
            conc_df.groupby("field_class")["n_match"].mean().to_dict()
            if not conc_df.empty else {}),
        "sr_lr_excluded_undetermined": safe_count(n_excluded),
        "yob_stratification": "applied" if yob else "skipped (no --yob-tsv)",
        "params": {"resolution": args.resolution, "seed": args.seed,
                   "bootstrap": args.bootstrap,
                   "strict_ancestry_min": args.strict_ancestry_min,
                   "min_ancestry_n": args.min_ancestry_n,
                   "limit_pairs": args.limit_pairs, "threads": args.threads,
                   "sharing_thresholds": list(SHARING_THRESHOLDS)},
    }
    with open(os.path.join(args.out_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # Person-level rows: the ONLY place they may go (never --out-dir).
    if args.local_dir:
        os.makedirs(args.local_dir, exist_ok=True)
        pairs_out.to_csv(os.path.join(args.local_dir, "26_pair_level_rows.tsv"),
                         sep="\t", index=False)
    print(f"Done. Aggregate outputs in {args.out_dir}", file=sys.stderr)
    return {"pairs": pairs_out, "decomposition": decomp_df, "class": class_df,
            "sweep": sweep_df, "homozygosity": hom_df, "fragmentation": frag_df,
            "power": power_df, "power_summary": power_summary, "sr_lr": sr_lr_agg,
            "summary": summary}


def main():
    main_from_args(build_parser().parse_args())


if __name__ == "__main__":
    main()
