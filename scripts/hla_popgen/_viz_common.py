#!/usr/bin/env python3
"""Shared constants + helpers for the three figure scripts in this directory
(`05_figures_frequency.py`, `06_figures_structure.py`, `07_figures_crosscohort.py`).

Not a SCHEMA.md table -- this is presentation-layer glue, kept in exactly ONE place so no two
figures can silently disagree on an ancestry's color, a gene list, or a confidence-interval
formula. The three figure scripts import this module; nothing under 00-04 imports or is imported
by it (those are owned by a different concurrent agent per the task boundary).

Read scripts/hla_popgen/SCHEMA.md and scripts/hla_popgen/research/VIZ_LIT.md before touching this
file -- every constant/helper here exists to satisfy a specific requirement documented there.
"""
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scipy.special import gammaln
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

try:
    import umap  # noqa: F401
    HAVE_UMAP = True
except ImportError:
    HAVE_UMAP = False

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_THIS_DIR))
DEFAULT_REPORT_ROOT = os.path.join(_REPO_ROOT, "reports", "hla_popgen")

# ---------------------------------------------------------------------------
# Ancestry palette -- VIZ_LIT.md Part 4.1: no verified official gnomAD/AoU hex-code ancestry
# palette could be found/confirmed. This is the Okabe-Ito colorblind-safe qualitative palette
# (Okabe & Ito 2008), used as a documented substitute, NOT a claim of matching gnomAD/AoU exactly.
# *** TRIVIAL TO SWAP: every figure in 05/06/07 imports this one dict. Replace the six hex values
# below the moment an official AoU/gnomAD ancestry palette is confirmed -- nothing else changes. ***
# ---------------------------------------------------------------------------
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
ANCESTRY_COLORS = {
    "AFR": "#D55E00",  # Okabe-Ito vermillion
    "AMR": "#E69F00",  # Okabe-Ito orange
    "EAS": "#009E73",  # Okabe-Ito bluish green
    "EUR": "#0072B2",  # Okabe-Ito blue
    "MID": "#CC79A7",  # Okabe-Ito reddish purple
    "SAS": "#56B4E9",  # Okabe-Ito sky blue
}
MISSING_COLOR = "#999999"

CLASSICAL_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
                    "HLA-DRB1"]
CLASSICAL_GENES_BARE = [g.replace("HLA-", "") for g in CLASSICAL_GENES]

# SCHEMA.md gene_class table -- the non-classical/framework genes VIZ_LIT.md 3.5 flags as unique
# to the long-read cohort (never analysed at population scale before this project).
# BARE names (no 'HLA-' prefix) -- bug found against the fixtures: every filter in this project
# matches against `gene_bare` (gene_bare() always strips 'HLA-'), so a mixed list with 'HLA-E' /
# 'HLA-DRB3' etc. silently matched ZERO rows for 8 of these 12 genes (only MICA/MICB/TAP1/TAP2,
# which never had a prefix to begin with, ever matched) -- the non-classical diversity map
# rendered those 8 genes as fully blank with no data, no warning.
NONCLASSICAL_GENES = ["MICA", "MICB", "TAP1", "TAP2", "DRB3", "DRB4", "DRB5",
                       "DQA2", "DQB2", "E", "F", "G"]

TD_STRATA = ["exact", "near", "loose", "distant"]

# Ancestry x gene (x allele) cells below this raw observation count are suppressed or visibly
# flagged, never plotted as if solid -- this repo has been bitten by noisy thin-N tails before
# (VIZ_LIT.md 4.2, task instructions "hard rule").
MIN_CELL_N = 5

# A SEPARATE, stricter threshold for statistics that are an AVERAGE ACROSS MANY ALLELES within one
# ancestry (e.g. mean SR-vs-LR frequency disagreement, resolution-cascade conditional entropy,
# template_distance summaries) rather than a single raw allele-copy count. Bug found in review: a
# per-ancestry bar chart used MIN_CELL_N (5) as its thin-N cutoff and MID/SAS (n_people=6-8 in the
# fixtures, genuinely the smallest AoU ancestry groups in real data too) cleared it, so their
# noisy, few-people-driven averages rendered as ordinary full-weight bars and visually dominated
# the well-powered EUR/AFR signal (n=95-110) the figure existed to show. 5 raw allele copies is a
# defensible floor for a single binomial proportion; it is not enough independent people to trust
# an average taken ACROSS an entire allele set for that ancestry. Never reuse MIN_CELL_N for this
# purpose -- pick this one instead.
MIN_CELL_N_PEOPLE = 10

NULL_TOKENS = {"", "NA", "NAN", "NONE", ".", "-", "UNDETERMINED"}


# ---------------------------------------------------------------------------
# Ancestry normalization -- ENVIRONMENT.md quirk #30
# ---------------------------------------------------------------------------
def normalize_ancestry(series):
    """Uppercase/strip a pandas Series of ancestry labels. ENVIRONMENT.md quirk #30:
    `ancestry_pred` is lowercase in every real source file -- normalize on load, always, from
    every column that carries it, independently."""
    out = series.astype(str).str.strip().str.upper()
    out = out.where(~out.isin({"NAN", "NONE", ""}), other=pd.NA)
    return out


def parse_bool_column(series):
    """TSV round-trips a bool column as the literal strings 'True'/'False' -- coerce robustly
    rather than relying on pandas' type inference (which is usually right but is not a contract)."""
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "t", "yes"})


# ---------------------------------------------------------------------------
# Allele string parsing (consensus / SR call strings)
# ---------------------------------------------------------------------------
def parse_allele_fields(consensus):
    """`consensus` like 'HLA-A*02:01:01:01', possibly with 'new' spliced in at some depth
    (SCHEMA.md Table 1 novelty encoding), or an SR call like 'A*02:01'. Returns the list of
    colon-delimited fields AFTER the gene*/allele-group prefix, stopping before any 'new' token.
    Returns [] for null/undetermined/missing calls."""
    if consensus is None or (isinstance(consensus, float) and pd.isna(consensus)):
        return []
    s = str(consensus).strip()
    if not s or s.upper() in NULL_TOKENS:
        return []
    if "*" in s:
        s = s.split("*", 1)[1]
    fields = [f for f in s.split(":") if f != ""]
    clean = []
    for f in fields:
        if f.strip().lower() == "new":
            break
        clean.append(f)
    return clean


def to_nfield(consensus, n):
    """First n colon-delimited fields, joined, or None if fewer than n resolved fields exist."""
    fields = parse_allele_fields(consensus)
    if len(fields) < n:
        return None
    return ":".join(fields[:n])


_MIC_TAP = {"MICA", "MICB", "TAP1", "TAP2"}


def gene_display(gene_bare_name):
    """Cosmetic inverse of gene_bare() for axis labels: MIC/TAP genes never had an 'HLA-' prefix
    to begin with; every other bare name (classical or non-classical) is shown with it restored,
    matching SCHEMA.md's gene_class table naming."""
    if gene_bare_name in _MIC_TAP:
        return gene_bare_name
    return f"HLA-{gene_bare_name}"


def gene_bare(gene):
    """Strip the 'HLA-' prefix Table 1 uses; a no-op on already-bare gene names (SR calls, MICA,
    TAP1, ...)."""
    if not isinstance(gene, str):
        return gene
    return gene[4:] if gene.startswith("HLA-") else gene


# ---------------------------------------------------------------------------
# Loaders -- Table 1 / Table 2 / Table 4 per SCHEMA.md, plus the AoU-native SR genotype file
# (NOT a SCHEMA.md table: 02_build_cohorts.py deliberately keeps only *membership* from
# hla_genotypes.tsv for Table 4, per its own docstring, so the actual SR allele calls needed for
# cross-cohort frequency comparison (VIZ_LIT.md 3.1, the headline novel figure) have to be read
# directly from the raw file here).
# ---------------------------------------------------------------------------
def load_table1(path):
    if not os.path.exists(path):
        sys.exit(f"FATAL: Table 1 (hla_calls_rich.tsv) not found at {path!r}. "
                 f"Run 01_extract_rich.py first.")
    df = pd.read_csv(path, sep="\t", dtype={"person_id": str, "consensus": str, "gene": str,
                                             "contig": str, "hap": str})
    for c in ("is_novel", "has_warning"):
        if c in df.columns:
            df[c] = parse_bool_column(df[c])
    for c in ("template_distance", "template_distance_per_kb", "cds_distance", "n_aa_changes",
              "n_fields", "n_tied", "gene_start", "gene_end", "copy_index"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["gene_bare"] = df["gene"].map(gene_bare)
    return df


def load_table2(path):
    if not os.path.exists(path):
        sys.exit(f"FATAL: Table 2 (hla_cis_pairs.tsv) not found at {path!r}. "
                 f"Run 01_extract_rich.py first.")
    df = pd.read_csv(path, sep="\t", dtype={"person_id": str, "hap": str, "contig": str,
                                             "pair": str, "allele_a": str, "allele_b": str,
                                             "haplotype_label": str})
    for c in ("both_exact", "either_novel"):
        if c in df.columns:
            df[c] = parse_bool_column(df[c])
    return df


def load_cohort_membership(path):
    if not os.path.exists(path):
        sys.exit(f"FATAL: Table 4 (cohort_membership.tsv) not found at {path!r}. "
                 f"Run 02_build_cohorts.py first.")
    df = pd.read_csv(path, sep="\t", dtype={"person_id": str})
    if "ancestry_pred" in df.columns:
        df["ancestry_pred"] = normalize_ancestry(df["ancestry_pred"])
    for c in ("in_sr", "in_lr"):
        if c in df.columns:
            df[c] = parse_bool_column(df[c])
    for c in ("max_template_distance", "mean_template_distance", "n_genes_called",
              "n_genes_exact") + tuple(f"p_{a.lower()}" for a in ANCESTRY_ORDER):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load_sr_genotypes(path):
    """Wide -> long: research_id + <gene>_1/<gene>_2 columns -> one row per (person_id, gene,
    copy, allele). Vectorized (pd.melt) rather than iterrows -- this file is ~500K rows x 16
    columns at full AoU scale."""
    if not os.path.exists(path):
        sys.exit(f"FATAL: AoU-native SR genotype file not found at {path!r} (expected "
                 f"research_id + <gene>_1/<gene>_2 columns for the 8 classical genes, e.g. "
                 f"v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv).")
    df = pd.read_csv(path, sep="\t", dtype=str)
    id_col = "research_id" if "research_id" in df.columns else "person_id"
    if id_col not in df.columns:
        sys.exit(f"FATAL: {path} has neither 'research_id' nor 'person_id'. "
                 f"Actual columns: {list(df.columns)}")
    allele_cols = [c for c in df.columns if c != id_col and c.rsplit("_", 1)[-1] in ("1", "2")]
    if not allele_cols:
        sys.exit(f"FATAL: {path} has no <gene>_1/<gene>_2 columns. Actual: {list(df.columns)}")
    long_df = df.melt(id_vars=[id_col], value_vars=allele_cols, var_name="col", value_name="allele")
    long_df = long_df.rename(columns={id_col: "person_id"})
    split = long_df["col"].str.rsplit("_", n=1, expand=True)
    long_df["gene_bare"] = split[0]
    long_df["copy"] = split[1].astype(int)
    long_df = long_df.drop(columns=["col"])
    return long_df


# ---------------------------------------------------------------------------
# The three-cohort selector -- SCHEMA.md: "Cohort 3 is a sweep, not a fixed threshold."
# ---------------------------------------------------------------------------
def select_cohort_people(cohort_df, cohort, td_max=None):
    """Returns (subset_of_cohort_membership, human_label) for one of the three cohorts.

    cohort: 'sr' (AoU-native short-read, ~500K) | 'lr' (Immuannot long-read, full, ~12K) |
    'lr_td' (long-read filtered on template_distance -- a SWEEP parameter, not a fixed cutoff;
    callers building a sweep figure should call this once per --td-max value in the sweep list,
    never hardcode a single threshold)."""
    if cohort == "sr":
        if "in_sr" not in cohort_df.columns:
            sys.exit("FATAL: cohort_membership.tsv has no 'in_sr' column.")
        return cohort_df[cohort_df["in_sr"]].copy(), "sr (AoU-native short-read)"
    if cohort == "lr":
        if "in_lr" not in cohort_df.columns:
            sys.exit("FATAL: cohort_membership.tsv has no 'in_lr' column.")
        return cohort_df[cohort_df["in_lr"]].copy(), "lr (Immuannot long-read, full)"
    if cohort == "lr_td":
        if td_max is None:
            sys.exit("FATAL: --cohort lr_td requires --td-max (e.g. --td-max 0, 1, 5). "
                     "This is a sweep parameter (SCHEMA.md) -- never hardcode one value across "
                     "the project; pass the specific cut you want for this run.")
        sub = cohort_df[cohort_df["in_lr"] & (cohort_df["max_template_distance"] <= td_max)].copy()
        return sub, f"lr_td{td_max} (Immuannot, max_template_distance<={td_max})"
    sys.exit(f"FATAL: unknown --cohort {cohort!r}; expected one of sr/lr/lr_td.")


def cohort_label_slug(cohort, td_max=None):
    return f"lr_td{td_max}" if cohort == "lr_td" else cohort


def add_common_args(ap):
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT,
                     help="Directory with hla_calls_rich.tsv, hla_cis_pairs.tsv, "
                          "cohort_membership.tsv, hla_genotypes.tsv (real run: ~/pipeline_outputs).")
    ap.add_argument("--table1", default=None, help="Override path to hla_calls_rich.tsv.")
    ap.add_argument("--table2", default=None, help="Override path to hla_cis_pairs.tsv.")
    ap.add_argument("--cohort-membership", default=None,
                     help="Override path to cohort_membership.tsv (Table 4).")
    ap.add_argument("--sr-genotypes", default=None,
                     help="Override path to hla_genotypes.tsv (AoU-native SR calls; not a "
                          "SCHEMA.md table -- read directly since Table 4 only carries SR "
                          "*membership*, not SR allele calls).")
    ap.add_argument("--cohort", choices=["sr", "lr", "lr_td"], required=True,
                     help="Which of the three cohorts to analyze: sr=AoU-native short-read "
                          "(~500K, 2-field, unphased, 8 classical genes); lr=Immuannot long-read, "
                          "full (~12K, phased, up to 4-field, 65 genes); lr_td=long-read filtered "
                          "on template_distance (requires --td-max; a sweep parameter).")
    ap.add_argument("--td-max", type=int, default=None,
                     help="Max per-person max_template_distance to keep, for --cohort lr_td.")
    ap.add_argument("--td-max-sweep", default="0,1,2,5,10",
                     help="Comma-separated max_template_distance cuts used by figures that show "
                          "how a statistic MOVES as the filter tightens (SCHEMA.md: 'the movement "
                          "... is itself the scientific result'). Never a single hardcoded value.")
    ap.add_argument("--out-dir", default=None,
                     help="Where to write PNGs + the markdown report. Default: "
                          "reports/hla_popgen/<script-name>/<cohort-label>/")
    ap.add_argument("--min-cell-n", type=int, default=MIN_CELL_N,
                     help="Ancestry x gene (x allele) cells with fewer than this many raw "
                          "observations are suppressed/flagged, never plotted as if solid.")
    ap.add_argument("--min-cell-n-people", type=int, default=MIN_CELL_N_PEOPLE,
                     help="Stricter threshold (people, not raw allele copies) for statistics "
                          "averaged ACROSS an ancestry's whole allele set -- e.g. SR-vs-LR "
                          "disagreement, resolution-cascade entropy, template_distance summaries. "
                          "See MIN_CELL_N_PEOPLE's docstring for why this must not reuse "
                          "--min-cell-n.")


def resolve_paths(args):
    table1 = args.table1 or os.path.join(args.outroot, "hla_calls_rich.tsv")
    table2 = args.table2 or os.path.join(args.outroot, "hla_cis_pairs.tsv")
    cohort_path = args.cohort_membership or os.path.join(args.outroot, "cohort_membership.tsv")
    sr_path = args.sr_genotypes or os.path.join(args.outroot, "hla_genotypes.tsv")
    return table1, table2, cohort_path, sr_path


def default_out_dir(script_name, cohort, td_max=None):
    return os.path.join(DEFAULT_REPORT_ROOT, script_name, cohort_label_slug(cohort, td_max))


# ---------------------------------------------------------------------------
# Statistics -- Wilson score CI (VIZ_LIT.md 4.3: preferred over Wald for rare alleles/small N),
# Shannon entropy, expected heterozygosity, Hurlbert rarefied allelic richness.
# ---------------------------------------------------------------------------
def wilson_ci(count, nobs, z=1.96):
    """Wilson score interval for a binomial proportion. Vectorized over numpy arrays.
    Returns (point_estimate, lo, hi), all NaN where nobs==0."""
    count = np.asarray(count, dtype=float)
    nobs = np.asarray(nobs, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(nobs > 0, count / np.where(nobs > 0, nobs, 1), np.nan)
        denom = 1 + z**2 / np.where(nobs > 0, nobs, 1)
        center = (p + z**2 / (2 * np.where(nobs > 0, nobs, 1))) / denom
        half = (z * np.sqrt(p * (1 - p) / np.where(nobs > 0, nobs, 1)
                             + z**2 / (4 * np.where(nobs > 0, nobs, 1)**2))) / denom
        lo = np.clip(center - half, 0, 1)
        hi = np.clip(center + half, 0, 1)
    lo = np.where(nobs > 0, lo, np.nan)
    hi = np.where(nobs > 0, hi, np.nan)
    return p, lo, hi


def shannon_entropy(freqs):
    freqs = np.asarray(freqs, dtype=float)
    freqs = freqs[freqs > 0]
    if len(freqs) == 0:
        return np.nan
    return float(-np.sum(freqs * np.log(freqs)))


def expected_heterozygosity(freqs):
    freqs = np.asarray(freqs, dtype=float)
    if len(freqs) == 0:
        return np.nan
    return float(1 - np.sum(freqs ** 2))


def _log_comb(n, k):
    if k < 0 or k > n or n < 0:
        return -np.inf
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


def rarefied_richness(counts, m):
    """Hurlbert (1971) rarefaction: expected number of distinct alleles in a sample of size m
    drawn WITHOUT replacement from a population whose per-allele copy counts are `counts`
    (summing to N total copies observed). VIZ_LIT.md 1.10: heterozygosity/richness are both
    sample-size sensitive, so cross-ancestry comparison needs a common rarefaction target rather
    than raw counts (AFR/AMR/EAS/MID/SAS will all be much smaller N than EUR in AoU)."""
    if not HAVE_SCIPY:
        return np.nan
    counts = np.asarray(counts, dtype=float)
    counts = counts[counts > 0]
    N = counts.sum()
    if N < m or m <= 0 or len(counts) == 0:
        return np.nan
    log_comb_N_m = _log_comb(N, m)
    total = 0.0
    for ni in counts:
        if N - ni < m:
            total += 1.0
        else:
            total += 1.0 - np.exp(_log_comb(N - ni, m) - log_comb_N_m)
    return float(total)


# ---------------------------------------------------------------------------
# Manual PCA (no sklearn dependency -- not in this project's stated dependency list).
# ---------------------------------------------------------------------------
def standardize(X):
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd = np.where(sd == 0, 1.0, sd)
    return (X - mu) / sd


def pca_fit_transform(X, n_components=10):
    """Simple mean/unit-variance-standardized PCA via SVD. Returns (scores, explained_var_ratio).
    X: (n_samples, n_features), already standardized by the caller if desired."""
    n_components = min(n_components, X.shape[0], X.shape[1])
    Xc = X - X.mean(axis=0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    scores = U[:, :n_components] * S[:n_components]
    var = (S ** 2) / max(X.shape[0] - 1, 1)
    var_ratio = var / var.sum()
    return scores, var_ratio[:n_components]


def silhouette_like(coords, labels):
    """Cheap silhouette-coefficient substitute (no sklearn): mean over points of
    (b - a) / max(a, b) where a = mean intra-label distance, b = mean nearest-other-label
    distance. Same [-1, 1] interpretation as sklearn's silhouette_score; not identical numerics
    (sklearn uses per-point nearest-cluster b, this does too, so it should track closely), kept
    dependency-free on purpose."""
    labels = np.asarray(labels)
    uniq = np.unique(labels)
    if len(uniq) < 2 or len(coords) < 3:
        return np.nan
    scores = []
    for i in range(len(coords)):
        own = labels[i]
        same = coords[labels == own]
        if len(same) > 1:
            a = np.mean(np.linalg.norm(same - coords[i], axis=1))
        else:
            a = 0.0
        b = np.inf
        for other in uniq:
            if other == own:
                continue
            pts = coords[labels == other]
            if len(pts) == 0:
                continue
            b = min(b, np.mean(np.linalg.norm(pts - coords[i], axis=1)))
        if not np.isfinite(b):
            continue
        scores.append((b - a) / max(a, b, 1e-12))
    return float(np.mean(scores)) if scores else np.nan


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------
def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def savefig(fig, path, dpi=150):
    ensure_dir(os.path.dirname(path))
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}", file=sys.stderr)


def write_report(path, lines):
    ensure_dir(os.path.dirname(path))
    text = "\n".join(lines)
    with open(path, "w") as f:
        f.write(text)
    print(f"(written report to {path})", file=sys.stderr)


def n_flag(n, min_cell_n):
    """Returns True if a cell's raw N is thin enough to suppress/flag."""
    return pd.isna(n) or n < min_cell_n
