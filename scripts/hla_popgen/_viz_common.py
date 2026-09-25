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


# ---------------------------------------------------------------------------
# Nature-style publication figures (S03 Task A). Extends this module, backward compatible --
# nothing above is touched, `savefig()` above still works exactly as before for the 05/06/07
# figure scripts. This section is for later, publication-grade figures (36_figure1_native.py's
# type scale was the model: one font size scale, no top/right spines, embedded fonts).
#
# Usage:
#     from _viz_common import nature_style, mm, save_fig, panel_letter, SUPPRESSED_COLOR, \
#         diverging_cmap, diverging_norm, hatch_suppressed, ANCESTRY_COLORS
#
#     with nature_style():
#         fig, ax = plt.subplots(figsize=(mm(NATURE_SINGLE_COL_MM), mm(60)))
#         ax.plot(...)
#         panel_letter(ax, "a")
#     save_fig(fig, "reports/hla_popgen/NN_x/my_panel")   # writes .pdf (embedded fonts) + .png (600 dpi)
#
# ANCESTRY_COLORS (defined near the top of this file) is reused here on purpose -- 36's ANC_COLORS
# dict has the identical six hex values, so this is the one place that dict should ever live;
# nothing here redefines it.
# ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
## cnsplots port (reference/FIGURE_STYLE.md, WS-C, S04). "Port, don't depend" --
## github.com/faridrashidi/cnsplots commit 634482c (v0.7.0), BSD-3-Clause, values only (rcParam
## defaults + palette hex), no runtime dependency added. What was ADOPTED vs where we deviate,
## with reasons, is spelled out here so a reader never has to diff against the upstream repo:
##
## Adopted verbatim (cnsplots `_settings.py` defaults):
##   - legend.frameon=False, legend's small unobtrusive glyphs (handled via LEGEND_KW below,
##     since matplotlib's rcParams has no direct markerscale/handlelength/handletextpad keys
##     that apply globally the way cnsplots' own legend() wrapper does -- see LEGEND_KW).
##   - axes.grid=False made EXPLICIT (was already true by matplotlib default + our own
##     nature_style unsetting nothing that turns it on -- cnsplots states it as a hard rule, so
##     we now pin it rather than rely on inherited default).
##   - axes.spines.top/right=False (already ours, unchanged).
##   - xtick.major.size / ytick.major.size = 2 (cnsplots: xtick_major_size=2) -- we previously
##     left tick length at matplotlib's default (3.5); shortened to match.
##   - mathtext.fontset="custom" (cnsplots) -- avoids matplotlib's default "dejavusans" mathtext
##     clashing visually with a Helvetica body font when a panel has a math expression.
##   - font.sans-serif extended with "Helvetica Neue", "Nimbus Sans", "Liberation Sans" (cnsplots'
##     fallback list) for cross-platform parity -- Linux CI boxes without Helvetica/Arial
##     installed now fall through to Nimbus Sans / Liberation Sans (metric-compatible clones)
##     before DejaVu Sans, rather than jumping straight to DejaVu.
##
## Deliberate DEVIATIONS from cnsplots (and why):
##   - savefig.dpi: we keep 600, not cnsplots' 288 (72x4). Nature's OWN artwork guidelines ask
##     for >= 300 dpi combination art and we've already shipped 600-dpi PNGs project-wide
##     (S03 CRITIC_2 checked figures at "native resolution" assuming 600 dpi); lowering to 288
##     now would be a visible regression for reviewers already looking at 600-dpi PNGs, for no
##     benefit (PDF vector export is unaffected either way -- dpi only touches the raster PNG).
##   - savefig.transparent: cnsplots defaults to True; we deliberately keep opaque white
##     backgrounds (save_fig() below does not set transparent=True) because these figures are
##     shared as flat PNGs in reports/Slack, where a transparent PNG on a dark viewer background
##     becomes illegible -- that's a distribution-context call, not a style disagreement.
##   - legend.fontsize stays 6 (already close to cnsplots' 7; unchanged rather than bumped, to
##     avoid a second silent change riding along on this port -- 6 vs 7 pt is not the thing
##     CRITIC_2 flagged, panel text sizes below the 5-7pt floor was).
##   - We do not adopt cnsplots' own plotting functions (barplot/heatmapplot/lineplot/...) or its
##     pvalue-star/forest/venn helpers -- FIGURE_STYLE.md's "port, don't depend" recommendation
##     covers rcParams + palette values only; our figures are hand-built matplotlib, not built on
##     cnsplots' API.
## ---------------------------------------------------------------------------
NATURE_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "Helvetica Neue", "Nimbus Sans",
                         "Liberation Sans", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 5,
    "ytick.labelsize": 5,
    "legend.fontsize": 6,
    "legend.frameon": False,     # cnsplots: legend_frameon=False -- no legend box, ever
    "axes.linewidth": 0.5,
    "axes.grid": False,          # cnsplots: axes_grid=False, made explicit (see note above)
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2,       # cnsplots: xtick_major_size=2 (was matplotlib default 3.5)
    "ytick.major.size": 2,
    "xtick.minor.width": 0.35,
    "ytick.minor.width": 0.35,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "mathtext.fontset": "custom",  # cnsplots: mathtext_fontset="custom"
    "pdf.fonttype": 42,   # embed as TrueType (Type 42), not Type 3 -- journals require this
    "ps.fonttype": 42,
    "savefig.dpi": 600,   # deviation from cnsplots' 288 -- see note above (we keep our 600 dpi)
}

# cnsplots legend_markerscale=0.5 / legend_handlelength=0.7 / legend_handletextpad=0.3 -- not
# global rcParams keys, so exposed as a kwargs dict callers can splat into ax.legend(**LEGEND_KW)
# (or fig.legend(**LEGEND_KW)) for the "small unobtrusive legend glyph" look, when a legend is
# unavoidable (prefer direct labelling first -- FIGURE_STYLE.md de-AI checklist item 8).
LEGEND_KW = {"frameon": False, "markerscale": 0.5, "handlelength": 0.7, "handletextpad": 0.3}

# Journal qualitative palettes (hex), via cnsplots `_palettes.py` (BSD-3-Clause, values only --
# not copyrightable expression, cited per FIGURE_STYLE.md). For any figure needing an alternate
# qualitative palette OTHER than ancestry -- e.g. distinguishing gene groups, cohorts, methods.
# Do NOT reuse these for ancestry: ANCESTRY_COLORS (Okabe-Ito, above) is the one true ancestry
# palette project-wide and must not be replaced by "Nature"/"Cell"/etc per FIGURE_STYLE.md.
JOURNAL_PALETTES = {
    "nature": ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2",
               "#DC0000", "#7E6148", "#B09C85"],
    "science": ["#3B4992", "#EE0000", "#008B45", "#631879", "#008280", "#BB0021", "#5F559B",
                "#A20056", "#808180", "#1B1919"],
    "lancet": ["#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A",
               "#ADB6B6", "#1B1919"],
    "nejm": ["#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91"],
    "cell": ["#C84C3A", "#2F7E8F", "#E1A22E", "#4E5A8A", "#5F9862", "#D07A6A", "#8B6FA8",
             "#7B8C9E", "#B85F7A", "#6B6B6B"],
}
# A restrained default accent for "one comparison, one hue" figures (de-AI checklist item 9) --
# picked from JOURNAL_PALETTES["nature"], not from ANCESTRY_COLORS, so it never gets confused
# with an ancestry encoding.
ACCENT_COLOR = JOURNAL_PALETTES["nature"][3]   # "#3C5488", a restrained slate blue

# Nature spec: panel letters are bold, lowercase, 8pt -- distinct from body text (5-7pt above).
PANEL_LETTER_FONTSIZE = 8

# Nature column widths in mm -- use with mm() for matplotlib figsize (always inches).
NATURE_SINGLE_COL_MM = 89
NATURE_DOUBLE_COL_MM = 183


def nature_style(apply=False):
    """Nature Publishing Group house style for matplotlib.

    Returns a context manager by default (`with nature_style(): ...`) that restores whatever
    rcParams were active on exit -- the normal, safe way to use this, since it can never leak
    style into a caller's other figures. Pass apply=True to instead push NATURE_RC onto the
    *global* rcParams with no restore, for the rare case where the style needs to outlive a
    `with` block (e.g. handed off to a helper in another module that draws later).
    """
    if apply:
        matplotlib.rcParams.update(NATURE_RC)
        return None
    return matplotlib.rc_context(rc=NATURE_RC)


def mm(x):
    """Millimetres -> inches, for matplotlib `figsize` (always inches). Nature widths: 89 mm
    single-column, 183 mm double-column/full-page (NATURE_SINGLE_COL_MM / NATURE_DOUBLE_COL_MM)."""
    return x / 25.4


def panel_letter(ax, letter, dx=-0.12, dy=1.05, fontsize=PANEL_LETTER_FONTSIZE):
    """Bold lowercase panel letter (Nature spec) just outside the axes, in axes-fraction
    coordinates. `letter` is lowercased automatically regardless of what's passed in."""
    ax.text(dx, dy, str(letter).lower(), transform=ax.transAxes, fontsize=fontsize,
            fontweight="bold", va="bottom", ha="left")


def save_fig(fig, path_stem, dpi=600, strict=True, reason=None):
    """Writes `<path_stem>.pdf` (vector, embedded TrueType/Type-42 fonts) and `<path_stem>.png`
    (raster, `dpi`, default 600 -- Nature's minimum for combination art), creating the parent
    directory if needed, then closes `fig`. Returns (pdf_path, png_path).

    Distinct from the plain `savefig()` above (which is PNG-only, 150 dpi, for the 05/06/07
    exploratory figures) -- this is the publication pair.

    Runs `check_layout(fig)` first (S04 WS-C phase 2 layout linter). By default (`strict=True`)
    any error-severity violation (text/text overlap, text or axes clipped beyond the figure bbox,
    a `mark_marginal()`-registered axes misaligned with its main axes) raises `RuntimeError`
    instead of writing a figure known to be broken. Warning-severity violations (excess
    whitespace) are printed but never block the save.

    Pass `strict=False` to save anyway despite error-severity violations -- this REQUIRES a
    `reason=` string (raises `ValueError` otherwise) so a deliberate override is never silent;
    the reason is printed alongside the suppressed violations.
    """
    violations = check_layout(fig)
    errors = [v for v in violations if v["severity"] == "error"]
    warns = [v for v in violations if v["severity"] == "warning"]
    for w in warns:
        print("  [layout WARNING] %s: %s" % (w["type"], w["detail"]), file=sys.stderr)
    if errors:
        msg = "\n".join("  - %s: %s" % (e["type"], e["detail"]) for e in errors)
        if strict:
            raise RuntimeError(
                "check_layout() found %d layout violation(s) for %r; refusing to save "
                "(strict=True, the default). Fix the layout, or call save_fig(..., "
                "strict=False, reason='...') to override deliberately:\n%s"
                % (len(errors), path_stem, msg))
        if not reason:
            raise ValueError(
                "save_fig(strict=False) requires a mandatory reason= string explaining why "
                "these %d layout violation(s) are being saved anyway." % len(errors))
        print("  [layout OVERRIDE] strict=False, reason=%r. %d violation(s) suppressed:\n%s"
              % (reason, len(errors), msg), file=sys.stderr)

    path_stem = str(path_stem)
    parent = os.path.dirname(path_stem)
    if parent:
        ensure_dir(parent)
    pdf_path = path_stem + ".pdf"
    png_path = path_stem + ".png"
    with matplotlib.rc_context(rc={"pdf.fonttype": 42, "ps.fonttype": 42}):
        fig.savefig(pdf_path, dpi=dpi, bbox_inches="tight")
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {pdf_path}", file=sys.stderr)
    print(f"  wrote {png_path}", file=sys.stderr)
    return pdf_path, png_path


# ---------------------------------------------------------------------------
# Diverging colormap for signed statistics (e.g. signed LD, a difference-from-zero), always used
# with a norm centered at 0 (diverging_norm), plus a fixed "not observed / suppressed" color +
# hatch helper for AoU small-cell (<20) cells -- these must read as visibly distinct from real
# data, never folded into a color scale as if they were an ordinary low value (feedback:
# "Suppressed counts are not zero").
# ---------------------------------------------------------------------------
SUPPRESSED_COLOR = "#D9D9D9"   # light grey -- reserved; not a valid data color on any scale here


def diverging_cmap():
    """Blue-white-red diverging colormap for signed statistics. Returns the colormap only, not
    a norm -- pair it with diverging_norm(vmin, vmax) so the actual data range is centered at 0
    correctly for each figure (the range differs per figure; the colormap doesn't need to)."""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(
        "nature_diverging_bwr", ["#2166AC", "#F7F7F7", "#B2182B"], N=256)


def diverging_norm(vmin=-1.0, vmax=1.0):
    """TwoSlopeNorm centered at 0, for use with diverging_cmap(). vmin/vmax should bracket the
    actual data range; defaults to [-1, 1] for correlation-like statistics (e.g. signed D')."""
    from matplotlib.colors import TwoSlopeNorm
    return TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)


# ---------------------------------------------------------------------------
# Layout linter (S04 WS-C phase 2, reference/FIGURE_STYLE.md). CRITIC_WSC.md: the rcParam port
# above never touched panel/subplot *layout* (spacing, alignment, sizing) -- the thing Marc
# actually flagged ("overlapping text is clearly inadmissible", marginals not aligned to their
# heatmap, colorbar floating with huge gaps). This section catches those mechanically, at
# save_fig() time, instead of relying on a human eyeballing a thumbnail.
# ---------------------------------------------------------------------------
def mark_marginal(ax_margin, ax_main, axis="x"):
    """Mark `ax_margin` as a marginal plot of `ax_main` along `axis` ('x' or 'y'), so
    check_layout()'s alignment check (c) can verify the two axes' data-to-display mapping
    actually coincides (i.e. the marginal's bars really do sit over/beside the main axes' cells)
    rather than trusting that `sharex`/`sharey` (or matching gridspec widths) were set correctly.
    Call once, after both axes' data limits are finalized (post-plotting, pre-save_fig)."""
    ax_margin._layout_marginal_of = (ax_main, axis)


def mark_decoration(artist):
    """Mark a Line2D/Patch (or any artist with `get_window_extent`) as a structural "decoration"
    -- an annotation bracket, a group divider, a manually-drawn connector -- so check_layout()'s
    text-vs-decoration check (a2) treats it as something text must never sit on top of, the same
    way it already treats two Text artists.

    Why this exists (S04 WS-C phase 2 linter false negative, `reports/hla_popgen/
    37_dq_g1g2_signed_ld/fig_dq_g1g2_committed_MAIN_POOLED.png`): `check_layout()`'s original
    overlap check (a) only ever compared Text against Text. The G1/G2 brackets in
    `37c_dq_g1g2_from_committed.py` (`_bracket_v`/`_bracket_h`) are drawn as plain `ax.plot(...)`
    Line2D segments in their own thin axes -- never Text -- so the rotated "DQA1" ylabel drawn
    straight through the vertical bracket line, and the "DQB1" xlabel colliding with the bracket's
    horizontal cap, were both structurally invisible to the linter: it had nothing in its `texts`
    list to compare either collision against. `check_layout()` returned zero violations for a
    figure Marc rejected on sight for exactly this overlap.

    Callers opt in explicitly (rather than the linter guessing which Line2D/Patch objects are
    "decorative") because most Line2D/Patch objects in a figure ARE the data -- a bar chart's
    Rectangle patches, a scatter's markers, a line plot's series -- and text legitimately sits
    near or on top of those (a bar's value label, a heatmap annotation) with no layout fault.
    Blanket-flagging every Line2D/Patch against every Text would drown real faults in false
    positives. A bracket/divider/connector is drawn purely to guide the eye, never to carry data
    a text label would legitimately overlap, so it opts in.

    Usage: `ln, = ax.plot(...); vc.mark_decoration(ln)`, or `mark_decoration(some_patch)`."""
    artist._layout_decoration = True


def mark_label(text_artist):
    """Mark a Text artist as a direct data-series label (e.g. an end-of-line ancestry label on a
    saturation curve, placed via `ax.annotate`/`ax.text`) so `check_layout()`'s (a4) check verifies
    it does not sit on top of any Line2D *data* series in its OWN axes.

    Why this exists (Figure 1 v5 panel e, S04 WS-C, 2026-09-25): direct end-of-line labels were
    placed at each curve's own endpoint with a small vertical-repulsion pass against each OTHER
    LABEL, but never checked against the actual drawn LINES -- when one ancestry's curve (e.g.
    AFR) is plotted well past another ancestry's endpoint (e.g. AMR's line stops earlier in x),
    the still-rising AFR line can pass directly through where AMR's label was placed. Two check
    categories already existed (`mark_decoration()` for "this Line2D/Patch must never be covered
    by text" and the plain text-vs-text check) but neither covers "this Text must never be covered
    by an ORDINARY plotted data line" -- opting the label in explicitly (rather than checking every
    Text in every axes against every line, which would flag ordinary in-line annotations, e.g. a
    value printed inside its own bar) keeps this scoped to the direct-labeling use case it exists
    for.

    Usage: `t = ax.annotate("AFR", (x, y), ...); vc.mark_label(t)`."""
    text_artist._layout_direct_label = True


def check_layout(fig, tol_overlap_px=2.0, tol_clip_frac=0.08, tol_margin_px=0.5,
                  whitespace_warn_frac=0.55, grid_n=48):
    """Mechanical layout linter, run against a fully-drawn figure. Returns a list of violation
    dicts `{"type": str, "severity": "error"|"warning", "detail": str}`. Checks:

    (a) pairwise overlaps among all visible Text artists (titles, axis labels, tick labels,
        legend text, annotations) -- small tolerance (`tol_overlap_px`) so touching-but-not-
        overlapping text doesn't false-positive.
    (a2) any visible Text overlapping an artist explicitly registered via `mark_decoration()`
        (a Line2D/Patch that is a structural bracket/divider/connector, never plotted data) --
        added in S04 WS-C phase 2 after check (a) alone missed a rotated axis label drawn
        straight through a G1/G2 annotation bracket (a Line2D in a neighboring axes, so no Text
        object existed for it to collide with under check (a)). See `mark_decoration()`'s
        docstring for the full false-negative writeup.
    (a3) any visible Text overlapping a visible spine of an Axes it does NOT itself belong to --
        e.g. a panel's title or label creeping into the neighboring panel's frame. A spine is
        exempted against Text that belongs to its OWN Axes (tick/axis labels are expected to sit
        right at their own axes' edge; that's normal, not a fault).
    (a4) any Text registered via `mark_label()` (a direct end-of-line data-series label) that sits
        on top of a plotted Line2D data series in its OWN axes -- catches a direct label placed at
        its own curve's endpoint that lands on top of a DIFFERENT series' line still passing
        through that region (e.g. two converging discovery curves where one series' line extends
        well past another's endpoint). Segment-vs-bbox intersection, not just endpoint/vertex
        containment, so a label sitting mid-segment between two vertices is still caught. See
        `mark_label()`'s docstring.
    (b) any text or Axes extending far beyond the figure's own (nominal, pre-`bbox_inches=
        "tight"`) bbox -- tolerance is `tol_clip_frac` of the figure's width/height (default 5%),
        not a fixed pixel count, since an ordinary axis label dipping a few px past the nominal
        canvas edge is normal matplotlib behavior that `save_fig()`'s `bbox_inches="tight"`
        handles cleanly. What this catches is a label/axes placed with a large, layout-breaking
        offset (e.g. a hardcoded negative-axes-fraction position well outside the panel; default
        8% tolerance covers ordinary tick/axis-label overflow, confirmed empirically against a
        plain `plt.subplots()` figure with default margins) that
        forces `bbox_inches="tight"` to balloon the saved canvas with mostly-empty space --
        exactly the "sits far outside the plot" / "floats with huge empty gaps" fault this linter
        exists to catch, which (d)'s whitespace check corroborates independently.
    (c) marginal alignment: for any axes registered with `mark_marginal()`, the data-to-display
        mapping along the shared axis must coincide with the main axes' mapping within
        `tol_margin_px` (default 0.5px) at the marginal's own xlim/ylim endpoints and midpoint --
        catches a marginal whose bars don't actually line up with the main heatmap's columns/rows
        even when nothing text-related overlaps.
    (d) whitespace: fraction of the figure's area not covered by any Axes or Text bbox, on a
        coarse `grid_n`x`grid_n` raster (cheap, no exact polygon union needed). Above
        `whitespace_warn_frac` this is a WARNING only (large gaps aren't automatically wrong --
        e.g. a deliberately spacious legend column -- but worth a human glance), never an error.

    All of (a), (a2), (a3), (a4), (b), (c) are "error" severity; (d) is "warning" only. Call via `save_fig()`, which
    raises on any error-severity violation by default (see its `strict`/`reason` kwargs).
    """
    import matplotlib.text as mtext
    from matplotlib.transforms import Bbox

    # Draw at (at least) save_fig()'s actual raster dpi, not whatever `fig.dpi` happens to be
    # (matplotlib's figure default is 100 unless set at construction). Empirically, on at least
    # one dev machine, rasterizing certain glyphs (e.g. a superscript/mathtext character) at
    # dpi<=100 raises `RuntimeError: failed to load glyph` from freetype, while the exact same
    # figure draws fine at dpi>=150 -- i.e. this is a font-rasterization floor, not a matplotlib
    # API restriction. save_fig() always rasterizes the PNG at 600 dpi, so checking at a lower
    # preview dpi would both risk that spurious crash AND measure pixel positions that don't match
    # what's actually saved. Bump (never lower) fig.dpi for the duration of this check, restore
    # after -- callers of check_layout() must never observe a mutated fig.dpi.
    _orig_dpi = fig.dpi
    try:
        if fig.dpi < 150:
            fig.dpi = 150
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        fig_bbox = fig.bbox
        violations = []

        # ---- tick-label Text objects sitting outside their own axis's view limits ----
        # matplotlib's default tick Locator deliberately places one extra major (and sometimes
        # minor) tick just past each end of the data range -- that tick's Text is already
        # invisible in the rendered PNG/PDF (clip_on=True, clipped to the axes bbox) but
        # `get_window_extent()` returns its UNCLIPPED geometry, so treating every Text with
        # `get_visible()==True` as "on-canvas" made these phantom off-view labels collide with a
        # panel letter or a neighboring subplot's own off-view label sharing the same dead space
        # -- a linter false positive, not a real visual fault (found by a concurrent S04 WS-C
        # session against 39/43b/38b, worked around per-script there with a local
        # `_prune_offview_ticklabels()` helper before this fix centralized it here). `sorted()`,
        # not unpacked-in-order: an inverted axis returns get_xlim()/get_ylim() as (high, low).
        offview_tick_ids = set()
        for ax in fig.axes:
            for axis_obj, get_lim in ((ax.xaxis, ax.get_xlim), (ax.yaxis, ax.get_ylim)):
                lo, hi = sorted(get_lim())
                for minor in (False, True):
                    try:
                        locs = axis_obj.get_ticklocs(minor=minor)
                        labels = axis_obj.get_ticklabels(minor=minor)
                    except Exception:
                        continue
                    for loc, t in zip(locs, labels):
                        if loc < lo - 1e-9 or loc > hi + 1e-9:
                            offview_tick_ids.add(id(t))

        # ---- collect visible, non-empty Text artists with a real on-canvas extent ----
        texts = []
        for t in fig.findobj(mtext.Text):
            if not t.get_visible():
                continue
            if id(t) in offview_tick_ids:
                continue
            s = t.get_text()
            if s is None or s.strip() == "":
                continue
            try:
                bb = t.get_window_extent(renderer=renderer)
            except Exception:
                continue
            if bb.width <= 0 or bb.height <= 0:
                continue
            texts.append((t, bb))

        # ---- (a) pairwise text/text overlap ----
        pad = -tol_overlap_px / 2.0
        for i in range(len(texts)):
            t1, b1 = texts[i]
            b1p = b1.padded(pad)
            for j in range(i + 1, len(texts)):
                t2, b2 = texts[j]
                inter = Bbox.intersection(b1p, b2.padded(pad))
                if inter is not None and inter.width > 0 and inter.height > 0:
                    violations.append({
                        "type": "text_overlap", "severity": "error",
                        "detail": "%r overlaps %r (overlap area=%.1fpx^2)"
                                  % (t1.get_text(), t2.get_text(), inter.width * inter.height),
                    })

        # ---- (a2) text vs registered decoration (mark_decoration()) ----
        decorations = []
        for a in fig.findobj(lambda art: getattr(art, "_layout_decoration", False)):
            if not a.get_visible():
                continue
            try:
                dbb = a.get_window_extent(renderer=renderer)
            except Exception:
                continue
            if dbb.width < 0 or dbb.height < 0:
                continue
            decorations.append((a, dbb))

        for t, tb in texts:
            tbp = tb.padded(pad)
            for d, dbb in decorations:
                # A zero-area bbox (e.g. a perfectly vertical/horizontal line segment) still
                # needs a real overlap test -- pad it out to a hairline width/height first so
                # Bbox.intersection can register a genuine crossing, not just touch it.
                dbbp = Bbox.from_extents(
                    dbb.x0 - max(tol_overlap_px, 0.5), dbb.y0 - max(tol_overlap_px, 0.5),
                    dbb.x1 + max(tol_overlap_px, 0.5), dbb.y1 + max(tol_overlap_px, 0.5))
                inter = Bbox.intersection(tbp, dbbp)
                if inter is not None and inter.width > 0 and inter.height > 0:
                    violations.append({
                        "type": "text_decoration_overlap", "severity": "error",
                        "detail": "text %r overlaps decoration artist %r (overlap area=%.1fpx^2)"
                                  % (t.get_text(), d, inter.width * inter.height),
                    })

        # ---- (a3) text vs a visible spine of an Axes it does not belong to ----
        for ax in fig.axes:
            for side, spine in ax.spines.items():
                if not spine.get_visible():
                    continue
                try:
                    sbb = spine.get_window_extent(renderer=renderer)
                except Exception:
                    continue
                if sbb.width < 0 or sbb.height < 0:
                    continue
                sbbp = Bbox.from_extents(
                    sbb.x0 - max(tol_overlap_px, 0.5), sbb.y0 - max(tol_overlap_px, 0.5),
                    sbb.x1 + max(tol_overlap_px, 0.5), sbb.y1 + max(tol_overlap_px, 0.5))
                for t, tb in texts:
                    if getattr(t, "axes", None) is ax:
                        continue  # a spine touching its own tick/axis labels is normal
                    inter = Bbox.intersection(tb.padded(pad), sbbp)
                    if inter is not None and inter.width > 0 and inter.height > 0:
                        violations.append({
                            "type": "text_foreign_spine_overlap", "severity": "error",
                            "detail": "text %r overlaps spine %r of a different axes "
                                      "(overlap area=%.1fpx^2)"
                                      % (t.get_text(), side, inter.width * inter.height),
                        })

        # ---- (a4) direct-label text (mark_label()) vs plotted data Line2D in the SAME axes ----
        def _seg_intersects_bbox(x0, y0, x1, y1, bb):
            # Liang-Barsky-lite: clip the segment parametrically against the box's 4 half-planes;
            # if any parametric range survives, the segment crosses (or starts/ends inside) bb.
            dx, dy = x1 - x0, y1 - y0
            tmin, tmax = 0.0, 1.0
            for p, q in ((-dx, x0 - bb.x0), (dx, bb.x1 - x0),
                        (-dy, y0 - bb.y0), (dy, bb.y1 - y0)):
                if p == 0:
                    if q < 0:
                        return False
                    continue
                r = q / p
                if p < 0:
                    tmin = max(tmin, r)
                else:
                    tmax = min(tmax, r)
                if tmin > tmax:
                    return False
            return True

        for t, tb in texts:
            if not getattr(t, "_layout_direct_label", False):
                continue
            ax = getattr(t, "axes", None)
            if ax is None:
                continue
            tbp = tb.padded(pad)
            for ln in ax.get_lines():
                if getattr(ln, "_layout_decoration", False):
                    continue  # a bracket/divider line -- (a2) already governs this one
                try:
                    verts = ln.get_transform().transform(ln.get_path().vertices)
                except Exception:
                    continue
                if len(verts) < 2:
                    continue
                hit = False
                for (x0, y0), (x1, y1) in zip(verts[:-1], verts[1:]):
                    if _seg_intersects_bbox(x0, y0, x1, y1, tbp):
                        hit = True
                        break
                if hit:
                    violations.append({
                        "type": "label_over_line_data", "severity": "error",
                        "detail": "direct label %r sits on top of a plotted data line %r"
                                  % (t.get_text(), ln.get_label()),
                    })

        # ---- (b) clipped beyond figure bbox ----
        pad_x = max(2.0, tol_clip_frac * fig_bbox.width)
        pad_y = max(2.0, tol_clip_frac * fig_bbox.height)
        fb = Bbox.from_extents(fig_bbox.x0 - pad_x, fig_bbox.y0 - pad_y,
                                fig_bbox.x1 + pad_x, fig_bbox.y1 + pad_y)

        def _within(bb):
            return fb.x0 <= bb.x0 and bb.x1 <= fb.x1 and fb.y0 <= bb.y0 and bb.y1 <= fb.y1

        for t, bb in texts:
            if not _within(bb):
                violations.append({
                    "type": "text_clipped", "severity": "error",
                    "detail": "%r extends beyond figure bbox (text=%s, fig=%s)"
                              % (t.get_text(), tuple(round(v, 1) for v in bb.bounds),
                                 tuple(round(v, 1) for v in fig_bbox.bounds)),
                })
        for ax in fig.axes:
            try:
                abb = ax.get_window_extent(renderer=renderer)
            except Exception:
                continue
            if not _within(abb):
                violations.append({
                    "type": "axes_clipped", "severity": "error",
                    "detail": "axes %r extends beyond figure bbox" % (ax.get_label() or ax,),
                })

        # ---- (c) marginal alignment ----
        for ax in fig.axes:
            info = getattr(ax, "_layout_marginal_of", None)
            if info is None:
                continue
            ax_main, axis = info
            if axis == "x":
                lo, hi = ax.get_xlim()
                for d in (lo, (lo + hi) / 2.0, hi):
                    p_margin = ax.transData.transform((d, 0))[0]
                    p_main = ax_main.transData.transform((d, 0))[0]
                    delta = abs(p_margin - p_main)
                    if delta > tol_margin_px:
                        violations.append({
                            "type": "marginal_misaligned", "severity": "error",
                            "detail": "x=%.4g: marginal axes px=%.2f vs main axes px=%.2f "
                                      "(delta=%.2fpx > tol=%.2fpx)"
                                      % (d, p_margin, p_main, delta, tol_margin_px),
                        })
            elif axis == "y":
                lo, hi = ax.get_ylim()
                for d in (lo, (lo + hi) / 2.0, hi):
                    p_margin = ax.transData.transform((0, d))[1]
                    p_main = ax_main.transData.transform((0, d))[1]
                    delta = abs(p_margin - p_main)
                    if delta > tol_margin_px:
                        violations.append({
                            "type": "marginal_misaligned", "severity": "error",
                            "detail": "y=%.4g: marginal axes px=%.2f vs main axes px=%.2f "
                                      "(delta=%.2fpx > tol=%.2fpx)"
                                      % (d, p_margin, p_main, delta, tol_margin_px),
                        })
            else:
                raise ValueError("mark_marginal axis must be 'x' or 'y', got %r" % (axis,))

        # ---- (d) whitespace (coarse raster coverage, warning only) ----
        fx0, fy0, fx1, fy1 = fig_bbox.x0, fig_bbox.y0, fig_bbox.x1, fig_bbox.y1
        fw, fh = max(fx1 - fx0, 1e-9), max(fy1 - fy0, 1e-9)
        covered = np.zeros((grid_n, grid_n), dtype=bool)

        def _mark(bb):
            x0 = int(np.clip((bb.x0 - fx0) / fw * grid_n, 0, grid_n))
            x1 = int(np.clip(np.ceil((bb.x1 - fx0) / fw * grid_n), 0, grid_n))
            y0 = int(np.clip((bb.y0 - fy0) / fh * grid_n, 0, grid_n))
            y1 = int(np.clip(np.ceil((bb.y1 - fy0) / fh * grid_n), 0, grid_n))
            if x1 > x0 and y1 > y0:
                covered[y0:y1, x0:x1] = True

        for ax in fig.axes:
            try:
                _mark(ax.get_window_extent(renderer=renderer))
            except Exception:
                pass
        for _t, bb in texts:
            _mark(bb)
        frac_white = 1.0 - covered.mean()
        if frac_white > whitespace_warn_frac:
            violations.append({
                "type": "excess_whitespace", "severity": "warning",
                "detail": "%.1f%% of figure area not covered by any axes/text (warn threshold "
                          "%.0f%%)" % (frac_white * 100, whitespace_warn_frac * 100),
            })

        return violations
    finally:
        # Never let callers observe a mutated fig.dpi -- restore it even though (unusually for
        # this codebase) `save_fig()`'s own subsequent `fig.savefig(..., dpi=600)` calls pass
        # dpi explicitly and wouldn't be affected either way; the invariant is still worth holding
        # since `check_layout()` is also called standalone (report-only mode, this session).
        fig.dpi = _orig_dpi


def hatch_suppressed(ax, x, y, width, height, **kwargs):
    """Draws a light-grey hatched rectangle (SUPPRESSED_COLOR fill) marking one 'not observed /
    below n=20' cell on a heatmap or grid, so it reads as visibly different from anywhere on a
    color scale rather than as if it were a real low value. x, y, width, height are in data
    coordinates by default (pass transform=ax.transAxes via kwargs for axes-fraction placement),
    matching matplotlib.patches.Rectangle's constructor. Returns the patch (already added to ax)."""
    from matplotlib.patches import Rectangle
    zorder = kwargs.pop("zorder", 4)
    rect = Rectangle((x, y), width, height, facecolor=SUPPRESSED_COLOR, edgecolor="#999999",
                     hatch="////", linewidth=0.3, zorder=zorder, **kwargs)
    ax.add_patch(rect)
    return rect
