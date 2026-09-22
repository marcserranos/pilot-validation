#!/usr/bin/env python3
"""Deletion supplement figure -- re-rendered from the COMMITTED 30_hla_sv tables only.

Cole (call #8, sprint S03): "simplify the colour scheme -- drop the orange/green distinction,
keep everything else." Script 30's `fig_deletion_rates()` (scripts/hla_popgen/30_hla_structural_
variation.py, ~line 314) colour-coded three classes: green = positive control (known CNV genes),
grey = negative control (should never be deleted), orange = everything else. That maps biological
category onto colour, which is the exact thing Cole asked to remove. This script keeps the
statistic (bridged-deletion %, Wilson CIs) and the control framing, but marks the controls with
annotation instead of colour: one neutral colour for every gene, and the two control groups called
out with a bracket/label above the bars.

This script is deliberately offline-only. It reads the already-committed, already small-cell-
suppressed TSVs under reports/hla_popgen/30_hla_sv/ -- NOT the VM, NOT hla_calls_rich.tsv, NOT
anything requiring a pipeline_outputs mount. That means:

  * `pct_deleted_bridged` is trusted as-is: script 30 computed it from the EXACT (unsuppressed)
    counts before writing the table, so the percentage itself is not censored even where the
    underlying `n_deleted_disp` is written `<20`.
  * `n_deleted` for the Wilson CI is *reconstructed* from `pct_deleted_bridged` and
    `n_bridged_disp` (round(pct/100 * n_bridged)) wherever `n_deleted_disp` is a real number, and
    from the committed percentage directly wherever it is `<20` -- this discloses nothing beyond
    what's already computable by inverting the committed percentage, since both numbers are
    already public in the TSV. `<20` is parsed as censored (NaN), NEVER as 0 -- see
    scripts/hla_popgen/32_novel_allele_callouts.py's `parse_count()` for the same convention.
  * If `n_bridged_disp` itself is `<20` (a thin denominator), the gene/ancestry cell is dropped
    from the plot and counted in the README's caveats rather than rendered with a fabricated CI.

Panels:
  (i)  per-gene bridged-deletion % across every gene in the committed table, one neutral colour,
       Wilson 95% CIs, genes grouped by SCHEMA.md's gene_class (a documented proxy for genomic
       region absent a genomic-position column in the committed table), expected-CNV genes
       (DRB3/4/5, C4A/B) and negative-control genes (A/B/C/DRA/DQA1/DQB1/DPA1/DPB1/DRB1) marked
       by a bracket + label above their bars rather than by colour.
  (ii) deletion frequency by ancestry for the CNV genes (DRB3, DRB4, DRB5, C4A, C4B) as a
       dot-plot with Wilson 95% CIs, ancestry palette (ANCESTRY_COLORS from _viz_common.py).

Usage (offline, no VM):
    python3 scripts/hla_popgen/38b_deletion_supplement_fig.py

Run: scripts/hla_popgen/tests/test_deletion_supplement_fig.py
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import _viz_common as vc

_REPO_ROOT = os.path.dirname(os.path.dirname(_THIS_DIR))
DEFAULT_SV_DIR = os.path.join(_REPO_ROOT, "reports", "hla_popgen", "30_hla_sv")
DEFAULT_OUT_DIR = os.path.join(_REPO_ROOT, "reports", "hla_popgen", "38_hla_a_deletion_validation")

# SCHEMA.md "Gene classification" table order -- used as the genomic-position proxy since the
# committed 30_hla_sv tables carry gene_class but not gene_start/gene_end.
GENE_CLASS_ORDER = [
    "classical_I", "classical_II", "class_II_accessory", "class_II_paralog",
    "nonclassical_I", "pseudogene_I", "mic_tap", "complement", "other",
]

NEGATIVE_CONTROL_GENES = ["A", "B", "C", "DRA", "DQA1", "DQB1", "DPA1", "DPB1", "DRB1"]
POSITIVE_CONTROL_GENES = ["DRB3", "DRB4", "DRB5"]
CNV_GENES = ["DRB3", "DRB4", "DRB5", "C4A", "C4B"]

NEUTRAL_COLOR = "#4C72B0"   # single neutral colour for panel (i) -- replaces the green/orange split

SUPPRESS_BELOW = 20


def parse_count(v):
    """A suppressed cell ('<20') -> NaN, never 0 and never 20. Same convention as
    32_novel_allele_callouts.py's parse_count() -- NaN forces every downstream use to be
    explicit about censoring rather than silently treating a censored cell as absent."""
    if v is None:
        return float("nan")
    s = str(v).strip()
    if s == "" or s.upper() in {"NA", "NAN", "NONE"}:
        return float("nan")
    if s.startswith("<"):
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def is_censored(v):
    return str(v).strip().startswith("<")


def load_gene_rates(sv_dir):
    path = os.path.join(sv_dir, "deletion_rates_by_gene.tsv")
    if not os.path.exists(path):
        sys.exit(f"FATAL: {path} not found. This script reads the COMMITTED 30_hla_sv tables "
                 f"only -- run 30_hla_structural_variation.py on the VM first if it's missing.")
    df = pd.read_csv(path, sep="\t", dtype=str)
    df["pct_deleted_bridged"] = pd.to_numeric(df["pct_deleted_bridged"], errors="coerce")
    df["n_bridged"] = df["n_bridged_disp"].map(parse_count)
    df["n_bridged_censored"] = df["n_bridged_disp"].map(is_censored)
    df["n_deleted_disp_censored"] = df["n_deleted_disp"].map(is_censored)
    return df


def load_ancestry_rates(sv_dir):
    path = os.path.join(sv_dir, "deletion_rates_by_gene_ancestry.tsv")
    if not os.path.exists(path):
        sys.exit(f"FATAL: {path} not found. This script reads the COMMITTED 30_hla_sv tables "
                 f"only -- run 30_hla_structural_variation.py on the VM first if it's missing.")
    df = pd.read_csv(path, sep="\t", dtype=str)
    df["pct_deleted_bridged"] = pd.to_numeric(df["pct_deleted_bridged"], errors="coerce")
    df["n_bridged"] = df["n_bridged_disp"].map(parse_count)
    df["n_bridged_censored"] = df["n_bridged_disp"].map(is_censored)
    return df


def add_wilson_ci(df):
    """Reconstructs n_deleted from the committed (uncensored) percentage and n_bridged, then
    attaches a Wilson 95% CI. Rows with a censored (thin) denominator get NaN CIs and are flagged
    `plottable=False` -- the caller drops them from the figure and reports them in the README."""
    df = df.copy()
    df["plottable"] = ~df["n_bridged_censored"] & df["n_bridged"].notna() & (df["n_bridged"] > 0)
    n_deleted_est = np.where(
        df["plottable"], np.round(df["pct_deleted_bridged"] / 100.0 * df["n_bridged"]), np.nan)
    df["n_deleted_est"] = n_deleted_est
    p, lo, hi = vc.wilson_ci(np.where(df["plottable"], n_deleted_est, np.nan),
                              np.where(df["plottable"], df["n_bridged"], np.nan))
    # Use the committed exact percentage as the point estimate (it's exact; the Wilson center
    # from the *rounded* reconstructed count would needlessly disagree with it in the last digit).
    df["ci_lo"] = lo * 100.0
    df["ci_hi"] = hi * 100.0
    return df


# ---------------------------------------------------------------------------
# panel (i): per-gene deletion rate, one colour, controls annotated not coloured
# ---------------------------------------------------------------------------
def draw_panel_gene(ax, gene_rates):
    df = gene_rates.copy()
    df["class_rank"] = df["gene_class"].map(
        {c: i for i, c in enumerate(GENE_CLASS_ORDER)}).fillna(len(GENE_CLASS_ORDER))
    df = df.sort_values(["class_rank", "gene"]).reset_index(drop=True)
    plottable = df[df["plottable"]].reset_index(drop=True)
    dropped = df[~df["plottable"]]

    x = np.arange(len(plottable))
    yerr_lo = plottable["pct_deleted_bridged"] - plottable["ci_lo"]
    yerr_hi = plottable["ci_hi"] - plottable["pct_deleted_bridged"]
    ax.bar(x, plottable["pct_deleted_bridged"], color=NEUTRAL_COLOR, width=0.72, zorder=3)
    ax.errorbar(x, plottable["pct_deleted_bridged"], yerr=[yerr_lo.clip(lower=0), yerr_hi.clip(lower=0)],
                fmt="none", ecolor="#222222", elinewidth=0.5, capsize=1.5, capthick=0.5, zorder=4)

    ax.set_xticks(x)
    ax.set_xticklabels(plottable["gene"], rotation=90, fontsize=5)
    ax.set_ylabel("bridged haplotypes lacking the gene (%)")
    ax.grid(axis="y", lw=0.3, alpha=0.35, zorder=0)

    # class boundaries -- thin vertical separators, not colour
    classes = plottable["gene_class"].tolist()
    for i in range(1, len(classes)):
        if classes[i] != classes[i - 1]:
            ax.axvline(i - 0.5, color="#BBBBBB", lw=0.4, zorder=1)

    # controls: bracket + label above the bars, per instructions ("subtle annotation... rather
    # than colour")
    def _bracket(gene_list, label, y_frac):
        idx = [i for i, g in enumerate(plottable["gene"]) if g in gene_list]
        if not idx:
            return
        ymax = ax.get_ylim()[1]
        y = ymax * y_frac
        for lo_i, hi_i in _contiguous_runs(idx):
            ax.plot([lo_i - 0.3, lo_i - 0.3, hi_i + 0.3, hi_i + 0.3],
                    [y * 0.94, y, y, y * 0.94], color="#333333", lw=0.4, zorder=5)
            ax.text((lo_i + hi_i) / 2.0, y * 1.02, label, ha="center", va="bottom",
                    fontsize=4.6, color="#333333")

    ax.set_ylim(0, max(plottable["ci_hi"].max(), plottable["pct_deleted_bridged"].max()) * 1.30)
    _bracket(POSITIVE_CONTROL_GENES, "expected CNV", 0.86)
    _bracket(NEGATIVE_CONTROL_GENES, "negative control", 0.72)

    ax.set_title("Bridged-deletion rate per gene (all genes, one colour)", fontsize=7, loc="left")
    if len(dropped):
        ax.text(0.99, 0.97, f"{len(dropped)} gene(s) omitted: denominator <20",
                transform=ax.transAxes, ha="right", va="top", fontsize=4.6, color="#888888",
                style="italic")
    return dropped


def _contiguous_runs(idx):
    idx = sorted(idx)
    runs = []
    start = prev = idx[0]
    for i in idx[1:]:
        if i == prev + 1:
            prev = i
            continue
        runs.append((start, prev))
        start = prev = i
    runs.append((start, prev))
    return runs


# ---------------------------------------------------------------------------
# panel (ii): CNV genes by ancestry, dot plot, ancestry palette
# ---------------------------------------------------------------------------
def draw_panel_ancestry(ax, anc_rates):
    df = anc_rates[anc_rates["gene"].isin(CNV_GENES)].copy()
    genes = [g for g in CNV_GENES if g in set(df["gene"])]
    ancs = [a for a in vc.ANCESTRY_ORDER if a in set(df["ancestry"])]

    y_positions = {}
    y = 0
    yticks, yticklabels = [], []
    n_suppressed = 0
    for g in genes:
        for a in ancs:
            y_positions[(g, a)] = y
            y += 1
        yticks.append(np.mean([y_positions[(g, a)] for a in ancs]))
        # gene_display() assumes an HLA gene (restores 'HLA-'); C4A/C4B are not HLA genes, so
        # the raw gene name is used here rather than vc.gene_display(g).
        yticklabels.append(g)
        y += 1  # gap between gene groups

    for a in ancs:
        sub = df[df["ancestry"] == a]
        ys, xs, xerr_lo, xerr_hi = [], [], [], []
        for _, r in sub.iterrows():
            key = (r["gene"], a)
            if key not in y_positions:
                continue
            if not r["plottable"]:
                n_suppressed += 1
                ax.text(0.5, y_positions[key], "suppressed (n<20)", fontsize=4.4,
                        color=vc.SUPPRESSED_COLOR if hasattr(vc, "SUPPRESSED_COLOR") else "#999999",
                        va="center", ha="left",
                        bbox=dict(boxstyle="round,pad=0.15", fc="#EEEEEE", ec="none"))
                continue
            ys.append(y_positions[key])
            xs.append(r["pct_deleted_bridged"])
            xerr_lo.append(max(0.0, r["pct_deleted_bridged"] - r["ci_lo"]))
            xerr_hi.append(max(0.0, r["ci_hi"] - r["pct_deleted_bridged"]))
        if ys:
            ax.errorbar(xs, ys, xerr=[xerr_lo, xerr_hi], fmt="o", ms=2.6,
                        color=vc.ANCESTRY_COLORS.get(a, "#777777"),
                        ecolor=vc.ANCESTRY_COLORS.get(a, "#777777"), elinewidth=0.6, capsize=1.5,
                        capthick=0.6, label=a, zorder=3)

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("bridged haplotypes lacking the gene (%)")
    ax.grid(axis="x", lw=0.3, alpha=0.35, zorder=0)
    ax.legend(frameon=False, fontsize=5, ncol=3, loc="lower right", handletextpad=0.3,
              columnspacing=0.8)
    ax.set_title("CNV-gene deletion frequency by ancestry", fontsize=7, loc="left")
    return n_suppressed


def make_figure(gene_rates, anc_rates, out_stem):
    with vc.nature_style():
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(190)),
            gridspec_kw={"height_ratios": [1.15, 1.0], "hspace": 0.55})
        dropped = draw_panel_gene(ax1, gene_rates)
        vc.panel_letter(ax1, "a")
        n_suppressed = draw_panel_ancestry(ax2, anc_rates)
        vc.panel_letter(ax2, "b")
        fig.subplots_adjust(left=0.09, right=0.98, top=0.96, bottom=0.14)
    pdf_path, png_path = vc.save_fig(fig, out_stem)
    return pdf_path, png_path, dropped, n_suppressed


README_TEMPLATE = """# 38 -- deletion supplement figure (HLA-A / MHC gene deletions, style pass)

Source: script `38b_deletion_supplement_fig.py`, re-rendered from the already-committed tables
in `reports/hla_popgen/30_hla_sv/` (`deletion_rates_by_gene.tsv`,
`deletion_rates_by_gene_ancestry.tsv`). No VM access, no re-computation from raw calls -- every
number here was already public.

Cole (call #8): "simplify the colour scheme -- drop the orange/green distinction, keep
everything else." Script 30's original `fig_deletion_rates()` colour-coded genes green
(positive control) / grey (negative control) / orange (everything else). This version uses one
neutral colour for every gene; the two control groups are called out with a bracket + label
above their bars instead.

## Panel (a) -- per-gene bridged-deletion rate

Every gene in the committed table, bridged-deletion % (Wilson 95% CI), grouped by SCHEMA.md's
`gene_class` (used here as a documented stand-in for genomic position, since the committed table
carries gene class but not gene_start/gene_end coordinates). Thin grey vertical lines mark class
boundaries. The two control groups from script 30 -- `{pos}` (expected copy-number genes) and
`{neg}` (negative controls that should essentially never be deleted) -- are marked with a bracket
and label above their bars, not by colour.

**How to read it:** a real, high bridged-deletion rate at DRB3/4/5 is the expected biology (DR51/
52/53 haplotype groups routinely lack the paralogous DRB locus). A near-zero rate at the negative
controls is the method's false-positive floor; anything meaningfully above zero there is measuring
assembly/detection error, not biology.

**Caveat:** {n_dropped} gene(s) were omitted from this panel because their bridged-haplotype
denominator itself is small-cell suppressed (`<20`) in the committed table -- see
`deletion_rates_by_gene.tsv`.

## Panel (b) -- CNV-gene deletion frequency by ancestry

The five expected-CNV genes ({cnv}) as a dot-plot: bridged-deletion % with Wilson 95% CI, one dot
per (gene, ancestry), coloured by the project's standard ancestry palette
(`ANCESTRY_COLORS` in `_viz_common.py`).

**How to read it:** DRB3/4/5 frequencies track DRB1 haplotype-group frequencies, which differ by
ancestry for well-established population-genetic reasons (not an artifact of this pipeline); C4A/
C4B copy number is independently variable. A cell drawn as "suppressed (n<20)" text instead of a
dot means the (gene, ancestry) bridged-haplotype count itself was `<20` in the committed table --
{n_suppressed} such cell(s) here -- and no rate is shown for it, per the AoU small-cell rule.
`<20` was parsed as censored throughout (never as 0).

## Caveats

- This figure cannot detect anything script 30 did not already detect: it is a restyle of
  committed, already-suppressed tables, not a new analysis.
- `pct_deleted_bridged` is the EXACT percentage script 30 computed before applying small-cell
  suppression to the underlying counts, so it is trusted as-is; the Wilson CI's count term is
  reconstructed from that exact percentage and the (uncensored) bridged-haplotype denominator,
  which discloses nothing beyond what the committed TSV already implies algebraically.
- "Deleted" here means *bridged absence* (script 30's stricter definition: a single contig spans
  the gene's position and the gene is not there) -- not raw absence, which conflates real
  deletions with assembly fragmentation. See `30_hla_structural_variation.py`'s module docstring.
"""


def write_readme(out_dir, gene_rates, dropped, n_suppressed):
    text = README_TEMPLATE.format(
        pos=", ".join(POSITIVE_CONTROL_GENES), neg=", ".join(NEGATIVE_CONTROL_GENES),
        cnv=", ".join(CNV_GENES), n_dropped=len(dropped), n_suppressed=n_suppressed)
    path = os.path.join(out_dir, "README.md")
    with open(path, "w") as f:
        f.write(text)
    print(f"  wrote {path}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sv-dir", default=DEFAULT_SV_DIR,
                     help="Directory with the committed 30_hla_sv tables.")
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    args = ap.parse_args()

    vc.ensure_dir(args.out_dir)

    gene_rates = add_wilson_ci(load_gene_rates(args.sv_dir))
    anc_rates = add_wilson_ci(load_ancestry_rates(args.sv_dir))

    out_stem = os.path.join(args.out_dir, "supp_deletions")
    pdf_path, png_path, dropped, n_suppressed = make_figure(gene_rates, anc_rates, out_stem)
    write_readme(args.out_dir, gene_rates, dropped, n_suppressed)

    print(f"[38b] wrote {pdf_path}, {png_path}, and README.md to {args.out_dir}")


if __name__ == "__main__":
    main()
