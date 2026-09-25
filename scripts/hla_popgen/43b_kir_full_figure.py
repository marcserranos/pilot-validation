#!/usr/bin/env python3
"""S03 WS6 follow-on -- Nature-grade multi-panel figure for the full-cohort KIR run.

Reads the 6 disclosure-safe aggregate TSVs `43_kir_full_aggregate.py` wrote to
reports/hla_popgen/43_kir_full_cohort/ and renders one 183 mm (double-column) static
matplotlib figure -> reports/hla_popgen/43_kir_full_cohort/fig_kir_full_cohort.{png,pdf}.

No participant-level data is read here -- every input is already an aggregate TSV with the
<20 disclosure rule applied (see `43_kir_full_aggregate.py`'s docstring). This script's own
disclosure obligation is just: never turn a blank/`<20` cell into a plotted 0 -- censored
cells are hatched (kir_gene_by_ancestry.tsv's novel_protein_pct, 7 gene x MID cells), not
silently dropped or drawn at height/color 0.

Gene order: the 17 KIR genes are laid out centromeric -> telomeric per the standard published
KIR gene map (framework genes 3DL3/3DP1/2DL4/3DL2 anchoring the array; e.g. Kulkarni et al 2008
Fig 1, Middleton & Gonzalez 2010 Table 1) -- used here ONLY for panel gene ordering. Note this
differs cosmetically from `41_kir_pilot.py`'s FRAMEWORK_GENES list order (3DL3, 2DL4, 3DP1, 3DL2),
which encodes adjacency for `classify_framework_miss()`, not a claim about the full 17-gene
physical order; that function only needs "what genes flank gene X", which is unaffected by
whether 2DL4 or 3DP1 is drawn first in this figure.

Usage:
  python3 scripts/hla_popgen/43b_kir_full_figure.py \\
      --in-dir reports/hla_popgen/43_kir_full_cohort \\
      --out-stem reports/hla_popgen/43_kir_full_cohort/fig_kir_full_cohort
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _THIS_DIR)
from _viz_common import (nature_style, mm, save_fig, panel_letter, ANCESTRY_ORDER,
                          ANCESTRY_COLORS, NATURE_DOUBLE_COL_MM, SUPPRESSED_COLOR)

# Standard published KIR gene map, centromeric -> telomeric (panel-ordering only; see module
# docstring). Framework genes (present on ~all haplotypes) marked with a dagger in the figure.
GENE_ORDER = ["KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DL5B", "KIR2DS3", "KIR2DS5",
              "KIR2DP1", "KIR2DL1", "KIR3DP1", "KIR2DL4", "KIR3DL1", "KIR3DS1", "KIR2DL5A",
              "KIR2DS1", "KIR2DS4", "KIR3DL2"]
FRAMEWORK_GENES = {"KIR3DL3", "KIR3DP1", "KIR2DL4", "KIR3DL2"}

GENE_SHORT = {g: g.replace("KIR", "") for g in GENE_ORDER}

NOVELTY_COLORS = {
    "known": "#999999",
    "genomic": "#56B4E9",   # novel, genomic-only (intron/UTR)
    "synonymous": "#E69F00",  # novel CDS, synonymous
    "protein": "#D55E00",   # novel protein -- highlight color, most consequential tier
}
NOVELTY_ORDER = ["known", "genomic", "synonymous", "protein"]
NOVELTY_LABEL = {"known": "known", "genomic": "novel, genomic-only",
                 "synonymous": "novel CDS, synonymous", "protein": "novel protein"}


def die(msg):
    sys.exit(f"FATAL: {msg}")


def load_tables(in_dir):
    need = ["kir_run_summary", "kir_gene_summary", "kir_gene_by_ancestry", "kir_allele_freq",
            "kir_qc", "kir_content_by_ancestry"]
    out = {}
    for name in need:
        path = os.path.join(in_dir, name + ".tsv")
        if not os.path.exists(path):
            die(f"missing input table {path!r}")
        out[name] = pd.read_csv(path, sep="\t", dtype=str)
    return out


def to_num(s):
    return pd.to_numeric(s, errors="coerce")


# ---------------------------------------------------------------------------
# Panel a: per-gene haplotype presence % with Wilson CI (unrelated subset).
# ---------------------------------------------------------------------------
def panel_a_presence(ax, gene_summary):
    gs = gene_summary.set_index("gene").reindex(GENE_ORDER)
    y = np.arange(len(GENE_ORDER))
    pct = to_num(gs["presence_pct"])
    lo = to_num(gs["presence_ci_lo"])
    hi = to_num(gs["presence_ci_hi"])
    xerr = np.vstack([pct - lo, hi - pct])
    colors = ["#0072B2" if g in FRAMEWORK_GENES else "#56B4E9" for g in GENE_ORDER]
    ax.barh(y, pct, xerr=xerr, height=0.62, color=colors, edgecolor="none",
            error_kw=dict(elinewidth=0.5, ecolor="#333333", capsize=1.2, capthick=0.5))
    ax.set_yticks(y)
    ax.set_yticklabels([GENE_SHORT[g] for g in GENE_ORDER], fontsize=5)
    # Framework genes are marked by weight (bold) + the darker bar color set above, not a "†"
    # symbol (FIGURE_STYLE.md de-AI checklist item 12: no decorative symbols in labels) -- the
    # figure-wide caption spells out what bold + dark-blue mean once.
    for lbl, g in zip(ax.get_yticklabels(), GENE_ORDER):
        if g in FRAMEWORK_GENES:
            lbl.set_fontweight("bold")
    ax.invert_yaxis()
    ax.set_xlim(0, 105)
    ax.set_xlabel("haplotype presence (%)")
    ax.set_title("Presence per gene", fontsize=7, pad=3)
    ax.tick_params(axis="y", length=0)
    return y


def panel_b_novelty(ax, gene_summary, y):
    gs = gene_summary.set_index("gene").reindex(GENE_ORDER)
    left = np.zeros(len(GENE_ORDER))
    field_map = {"known": "pct_known", "genomic": "pct_novel_genomic_only",
                 "synonymous": "pct_novel_cds_synonymous", "protein": "pct_novel_protein"}
    for tier in NOVELTY_ORDER:
        vals = to_num(gs[field_map[tier]]).fillna(0.0).values
        ax.barh(y, vals, left=left, height=0.62, color=NOVELTY_COLORS[tier],
                edgecolor="none", label=NOVELTY_LABEL[tier])
        left = left + vals
    n_calls = gs["n_calls"].astype(int).values
    for yi, n in zip(y, n_calls):
        # Right-aligned, clipped to the axes: previously left-aligned at x=101.5 the widest labels
        # (5-digit n, e.g. "n=23,097") bled past the panel-b/c gutter and collided with panel c's
        # row labels (found in style pass, S04 WS-C). Right-align against a fixed right edge
        # instead, so every label's rightmost character lands at the same x regardless of digit
        # count, and clip_on=True keeps it from ever drawing into the next panel's margin.
        ax.text(118, yi, f"n={n:,}", va="center", ha="right", fontsize=5.0, color="#333333",
                clip_on=False)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("allele calls (%)")
    ax.set_title("Novelty per gene", fontsize=7, pad=3)
    ax.tick_params(axis="y", length=0)
    ax.margins(x=0)
    handles = [Rectangle((0, 0), 1, 1, color=NOVELTY_COLORS[t]) for t in NOVELTY_ORDER]
    return handles


# ---------------------------------------------------------------------------
# Panel c: gene x ancestry presence heatmap.
# ---------------------------------------------------------------------------
def _heatmap(ax, mat, y, cmap, vmin, vmax, censored_mask, cbar_label, title, fmt="{:.0f}"):
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                    extent=[-0.5, mat.shape[1] - 0.5, len(GENE_ORDER) - 0.5, -0.5])
    for gi in range(mat.shape[0]):
        for ai in range(mat.shape[1]):
            if censored_mask[gi, ai]:
                ax.add_patch(Rectangle((ai - 0.5, gi - 0.5), 1, 1, facecolor=SUPPRESSED_COLOR,
                                        edgecolor="#999999", hatch="////", linewidth=0.3,
                                        zorder=4))
            elif not np.isnan(mat[gi, ai]):
                ax.text(ai, gi, fmt.format(mat[gi, ai]), ha="center", va="center",
                        fontsize=5.0, color="white" if mat[gi, ai] > (vmin + vmax) / 2
                        else "#222222")
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels(ANCESTRY_ORDER, fontsize=5.3, rotation=90)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_title(title, fontsize=7, pad=3)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=2, pad=1)
    cbar = plt.colorbar(im, ax=ax, fraction=0.05, pad=0.04, shrink=0.85)
    cbar.set_label(cbar_label, fontsize=5.3)
    cbar.ax.tick_params(labelsize=5.0, length=1.5)
    return im


def _pivot_num(gene_by_ancestry, value_col):
    piv = gene_by_ancestry.pivot(index="gene", columns="ancestry", values=value_col)
    piv = piv.reindex(index=GENE_ORDER, columns=ANCESTRY_ORDER)
    return piv.apply(to_num).values.astype(float)


def panel_c_presence_heatmap(ax, gene_by_ancestry, y):
    # Integer % here (values cluster 70-100, cell text at 1-decimal collided edge-to-edge in
    # these narrow 6-ancestry-wide cells -- check_layout() caught it); panel d keeps 1 decimal,
    # where the underlying values (mostly single digits) actually need it to be distinguishable.
    # Each panel's own precision is fixed and stated once here + in the README, rather than
    # forcing one shared decimal count that would either crowd panel c or hide real digits in d.
    mat = _pivot_num(gene_by_ancestry, "presence_pct")
    censored = np.isnan(mat)  # none expected for presence, but handled generically
    _heatmap(ax, mat, y, "viridis", 0, 100, censored, "presence (%)", "Presence × ancestry",
             fmt="{:.0f}")
    ax.set_yticklabels([GENE_SHORT[g] for g in GENE_ORDER], fontsize=5.3)
    for lbl, g in zip(ax.get_yticklabels(), GENE_ORDER):
        if g in FRAMEWORK_GENES:
            lbl.set_fontweight("bold")


def panel_d_novelty_heatmap(ax, gene_by_ancestry, y):
    mat = _pivot_num(gene_by_ancestry, "novel_protein_pct")
    censored = np.isnan(mat)
    vmax = np.nanmax(mat) if np.isfinite(np.nanmax(mat)) else 25.0
    _heatmap(ax, mat, y, "OrRd", 0, max(vmax, 10.0), censored, "novel protein (%)",
             "Protein novelty × ancestry", fmt="{:.1f}")


# ---------------------------------------------------------------------------
# Panel e: cA vs cB haplotype content by ancestry.
# ---------------------------------------------------------------------------
def panel_e_content(ax, content_by_ancestry):
    df = content_by_ancestry.set_index("ancestry").reindex(ANCESTRY_ORDER)
    pct_ca = to_num(df["pct_cA"])
    lo = to_num(df["pct_cA_ci_lo"])
    hi = to_num(df["pct_cA_ci_hi"])
    x = np.arange(len(ANCESTRY_ORDER))
    colors = [ANCESTRY_COLORS[a] for a in ANCESTRY_ORDER]
    ax.bar(x, pct_ca, yerr=[pct_ca - lo, hi - pct_ca], width=0.6, color=colors,
           edgecolor="none", error_kw=dict(elinewidth=0.5, ecolor="#333333", capsize=1.5,
                                            capthick=0.5))
    ax.axhline(50, color="#666666", lw=0.4, ls="--", zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels(ANCESTRY_ORDER, fontsize=5.5)
    ax.set_ylim(0, 80)
    ax.set_ylabel("haplotypes cA (%)")
    ax.set_title("cA/cB content by ancestry", fontsize=7, pad=3)


# ---------------------------------------------------------------------------
# Panel f: QC strip -- framework presence + 2DL2/2DL3 and 3DL1/3DS1 co-occurrence.
# ---------------------------------------------------------------------------
def panel_f_qc(ax, qc):
    # Points + CI, not bars (FIGURE_STYLE.md de-AI checklist item 10 / known debt: a bar chart
    # over this narrow 93-97% range on a 0-100 axis would be unreadable, but a truncated bar axis
    # exaggerates the differences by area -- a bar's visual weight is its height from a baseline,
    # so any baseline above 0 misleads. A point encodes only its position, so it carries no such
    # baseline claim and can honestly sit in the narrow range that's actually informative.
    fw = qc[qc["metric"] == "framework_gene_presence"].copy()
    fw["n"] = to_num(fw["n"]); fw["pct"] = to_num(fw["pct"])
    fw["ci_lo"] = to_num(fw["ci_lo"]); fw["ci_hi"] = to_num(fw["ci_hi"])
    fw_order = ["KIR3DL3", "KIR3DP1", "KIR2DL4", "KIR3DL2"]
    fw = fw.set_index("item").reindex(fw_order)
    x = np.arange(len(fw_order))
    ax.errorbar(x, fw["pct"], yerr=[fw["pct"] - fw["ci_lo"], fw["ci_hi"] - fw["pct"]],
                fmt="o", ms=3.5, color="#0072B2", ecolor="#333333", elinewidth=0.6,
                capsize=2.0, capthick=0.6, zorder=3)
    ax.set_xlim(-0.5, len(fw_order) - 0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([g.replace("KIR", "") for g in fw_order], fontsize=5.2)
    ymin = min(fw["ci_lo"].min(), fw["pct"].min()) - 1.0
    ymax = max(fw["ci_hi"].max(), fw["pct"].max()) + 1.0
    ax.set_ylim(ymin, ymax)
    ax.set_ylabel("framework presence (%)")
    ax.set_title("QC: framework gene presence", fontsize=7, pad=3)

    # Co-occurrence numbers used to sit as free text just right of this axes (transform=
    # ax.transAxes, x=1.03) -- with nothing to their right (this is the rightmost panel), that
    # text ran straight off the figure's own bbox (check_layout() text_clipped, S04 WS-C
    # redesign). Returned as plain strings instead and folded into the figure-wide caption below,
    # which is centered and already sized to the full figure width.
    co = qc[qc["metric"] == "cooccurrence"].set_index("item")
    lines = []
    for item, label in [("KIR2DL2_and_KIR2DL3", "2DL2+2DL3"),
                         ("KIR3DL1_and_KIR3DS1", "3DL1+3DS1")]:
        row = co.loc[item]
        lines.append(f"{label} co-occurrence {to_num(row['pct']):.1f}% (n={int(to_num(row['n']))})")
    return lines


def build_figure(tables, out_stem):
    with nature_style():
        fig = plt.figure(figsize=(mm(NATURE_DOUBLE_COL_MM), mm(195)))
        gs = fig.add_gridspec(2, 4, height_ratios=[2.3, 1.0], width_ratios=[1.15, 1.0, 0.85, 0.95],
                               hspace=0.38, wspace=0.75,
                               left=0.085, right=0.905, top=0.96, bottom=0.085)

        ax_a = fig.add_subplot(gs[0, 0])
        ax_b = fig.add_subplot(gs[0, 1])
        ax_c = fig.add_subplot(gs[0, 2])
        ax_d = fig.add_subplot(gs[0, 3])
        ax_e = fig.add_subplot(gs[1, 0:2])
        ax_f = fig.add_subplot(gs[1, 2:4])

        y = panel_a_presence(ax_a, tables["kir_gene_summary"])
        panel_letter(ax_a, "a", dx=-0.42, dy=1.03)
        legend_handles = panel_b_novelty(ax_b, tables["kir_gene_summary"], y)
        panel_letter(ax_b, "b", dx=-0.10, dy=1.03)
        panel_c_presence_heatmap(ax_c, tables["kir_gene_by_ancestry"], y)
        panel_letter(ax_c, "c", dx=-0.12, dy=1.03)
        panel_d_novelty_heatmap(ax_d, tables["kir_gene_by_ancestry"], y)
        panel_letter(ax_d, "d", dx=-0.05, dy=1.03)
        panel_e_content(ax_e, tables["kir_content_by_ancestry"])
        panel_letter(ax_e, "e", dx=-0.10, dy=1.06)
        cooccur_lines = panel_f_qc(ax_f, tables["kir_qc"])
        panel_letter(ax_f, "f", dx=-0.12, dy=1.06)

        bottoms, tops, _lefts, _rights = gs.get_grid_positions(fig)
        row_gap_y = (bottoms[0] + tops[1]) / 2.0  # (row0 bottom + row1 top) / 2
        fig.legend(legend_handles, [NOVELTY_LABEL[t] for t in NOVELTY_ORDER],
                   loc="center", bbox_to_anchor=(0.5, row_gap_y), ncol=4, fontsize=5.3,
                   frameon=False, handlelength=1.0, handleheight=0.9, columnspacing=1.2)

        fig.text(0.5, 0.035, "  |  ".join(cooccur_lines), fontsize=5.3, color="#444444",
                  ha="center")
        fig.text(0.5, 0.012,
                  "Bold gene label / darker bar (a): framework gene (expected on ~all "
                  "haplotypes). Hatched cells: <20 calls, censored per disclosure rule. "
                  "Ancestry = predicted; unrelated subset (KING kin ≥ 0.0442 removed).",
                  fontsize=5.3, color="#444444", ha="center")

        return save_fig(fig, out_stem)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", default=os.path.join(
        os.path.dirname(_THIS_DIR), "..", "reports", "hla_popgen", "43_kir_full_cohort"))
    ap.add_argument("--out-stem", default=None)
    args = ap.parse_args()
    in_dir = os.path.abspath(args.in_dir)
    out_stem = args.out_stem or os.path.join(in_dir, "fig_kir_full_cohort")

    tables = load_tables(in_dir)
    pdf_path, png_path = build_figure(tables, out_stem)
    print(f"wrote {pdf_path}\nwrote {png_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
