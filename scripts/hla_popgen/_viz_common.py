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


def save_fig(fig, path_stem, dpi=600):
    """Writes `<path_stem>.pdf` (vector, embedded TrueType/Type-42 fonts) and `<path_stem>.png`
    (raster, `dpi`, default 600 -- Nature's minimum for combination art), creating the parent
    directory if needed, then closes `fig`. Returns (pdf_path, png_path).

    Distinct from the plain `savefig()` above (which is PNG-only, 150 dpi, for the 05/06/07
    exploratory figures) -- this is the publication pair.
    """
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
