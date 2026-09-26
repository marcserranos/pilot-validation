#!/usr/bin/env python3
"""WS-D (Sprint S04) -- figures for `48_kir_hla_ligand_cooccurrence.py`'s aggregate TSVs.

Local-only: reads `reports/hla_popgen/48_kir_hla_ligand_cooccurrence/*.tsv` (pulled-back
aggregate exports, never raw per-person data) and renders one figure with two panels:

  (a) Forest plot of KIR-HLA receptor-ligand odds ratios (log scale, reference line at OR=1),
      one row per (pair, ancestry), rows grouped and ordered by biology -- inhibitory pairs first
      (2DL1xC2, 2DL2xC1, 2DL3xC1, 3DL1xBw4), then the one activating pair (3DS1xBw4, tests the
      epidemiologic Bw4-80I association, not a confirmed receptor-ligand binding pair -- see 48's
      own README caveats). A dotted reference line marks the Bonferroni threshold context (30
      pairs tested, alpha=0.05/30 ~= 0.0017) as a text note; the two nominal (uncorrected)
      signals are starred directly on their rows.
  (b) A compact C1/C2/Bw4 carrier-frequency-by-ancestry dot plot (`epitope_freq_by_ancestry.tsv`),
      for context alongside the forest plot's receptor-side numbers.

Disclosure (hard rule, `AGENT_PREAMBLE.md`): any 2x2 cell in [1, 19] is already masked upstream in
`kir_hla_ligand_cooccurrence.tsv` (blank `odds_ratio`/CI/p columns). This script never re-derives
or guesses an OR for a masked row -- it draws a hatched grey marker labelled "censored (<20)"
in that row instead, never a point at OR=1 or any other value that could read as real data.

Usage:
    python3 scripts/hla_popgen/48b_ligand_figure.py \\
        --in-dir reports/hla_popgen/48_kir_hla_ligand_cooccurrence \\
        --out-dir reports/hla_popgen/48_kir_hla_ligand_cooccurrence
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

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)
import _viz_common as vc

# Biological ordering: inhibitory pairs first, then the one activating/epidemiologic pair.
PAIR_ORDER = ["KIR2DL1xC2", "KIR2DL2xC1", "KIR2DL3xC1", "KIR3DL1xBw4", "KIR3DS1xBw4"]
PAIR_LABEL = {
    "KIR2DL1xC2": "2DL1 x C2 (inhibitory)",
    "KIR2DL2xC1": "2DL2 x C1 (inhibitory)",
    "KIR2DL3xC1": "2DL3 x C1 (inhibitory)",
    "KIR3DL1xBw4": "3DL1 x Bw4 (inhibitory)",
    "KIR3DS1xBw4": "3DS1 x Bw4-80I (activating, epidemiologic)",
}
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
N_PAIRS_TESTED = 30
BONFERRONI_ALPHA = 0.05 / N_PAIRS_TESTED
# The two nominal (uncorrected) single-ancestry signals called out in 48's README -- neither
# survives Bonferroni; starred on the plot, never presented as a finding.
NOMINAL_SIGNALS = {("EAS", "KIR2DL2xC1"), ("EUR", "KIR2DL3xC1")}


def load_tables(in_dir):
    cooc = pd.read_csv(os.path.join(in_dir, "kir_hla_ligand_cooccurrence.tsv"), sep="\t")
    cooc["pair"] = cooc["kir_gene"] + "x" + cooc["ligand"]
    epi = pd.read_csv(os.path.join(in_dir, "epitope_freq_by_ancestry.tsv"), sep="\t")
    return cooc, epi


def _row_positions(pairs, ancestries, group_gap=1.4):
    """Returns {(pair, ancestry): y} with a visual gap between pair-groups, y increasing top to
    bottom in plot order (first pair drawn at the top)."""
    y = 0.0
    pos = {}
    order = []
    for pair in pairs:
        for anc in ancestries:
            pos[(pair, anc)] = y
            order.append((pair, anc))
            y += 1.0
        y += group_gap
    return pos, order


def fig_forest_and_epitopes(cooc, epi, out_stem):
    pairs = [p for p in PAIR_ORDER if p in set(cooc["pair"])]
    pos, order = _row_positions(pairs, ANCESTRY_ORDER)
    ymax = max(pos.values())

    with vc.nature_style():
        fig = plt.figure(figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(185)),
                          constrained_layout=True)
        gs = fig.add_gridspec(1, 4)
        ax_a = fig.add_subplot(gs[0, :3])
        ax_b = fig.add_subplot(gs[0, 3])

        cooc_idx = cooc.set_index(["pair", "ancestry"])
        for pair, anc in order:
            y = ymax - pos[(pair, anc)]  # flip so first pair/ancestry is at the TOP
            if (pair, anc) not in cooc_idx.index:
                continue
            row = cooc_idx.loc[(pair, anc)]
            color = vc.ANCESTRY_COLORS.get(anc, vc.ACCENT_COLOR)
            censored = pd.isna(row["odds_ratio"])
            if censored:
                rect = Rectangle((0.55, y - 0.32), 1.2, 0.64, facecolor=vc.SUPPRESSED_COLOR,
                                  edgecolor="#999999", hatch="////", linewidth=0.4)
                ax_a.add_patch(rect)
                vc.mark_decoration(rect)
                t = ax_a.annotate("censored (<20 in a cell)", (1.85, y), va="center",
                                   ha="left", fontsize=4.6, color="#666666", style="italic")
                vc.mark_label(t)
                continue
            orv = float(row["odds_ratio"])
            lo, hi = float(row["or_ci_lo"]), float(row["or_ci_hi"])
            ax_a.plot([lo, hi], [y, y], color=color, lw=0.8, zorder=2)
            ax_a.scatter([orv], [y], s=10, color=color, zorder=3, edgecolor="white",
                         linewidth=0.3)
            if (anc, pair) in NOMINAL_SIGNALS:
                t = ax_a.annotate("*", (hi, y), xytext=(2, 0), textcoords="offset points",
                                   va="center", ha="left", fontsize=7, color=color,
                                   fontweight="bold")
                vc.mark_label(t)

        ax_a.axvline(1.0, color="#444444", lw=0.6, ls=(0, (2, 1.5)), zorder=1)
        ax_a.set_xscale("log")
        ax_a.set_xlim(0.35, 3.2)
        ax_a.set_xlabel("odds ratio (KIR present x ligand present), log scale")
        # Log-scale default tick formatting renders exponents via mathtext, which crashes under
        # this environment's font stack + mathtext.fontset="custom" (FT_Load_Glyph division-by-
        # zero) -- use plain decimal tick labels instead (no mathtext involved).
        xticks = [0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0]
        ax_a.set_xticks(xticks)
        ax_a.set_xticks([], minor=True)
        ax_a.xaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))

        # y ticks: ancestry labels colored per point; pair group labels to the left of those.
        yticks, yticklabels, ycolors = [], [], []
        for pair, anc in order:
            y = ymax - pos[(pair, anc)]
            yticks.append(y)
            yticklabels.append(anc)
            ycolors.append(vc.ANCESTRY_COLORS.get(anc, vc.ACCENT_COLOR))
        ax_a.set_yticks(yticks)
        ax_a.set_yticklabels(yticklabels, fontsize=4.6)
        for tick_lbl, c in zip(ax_a.get_yticklabels(), ycolors):
            tick_lbl.set_color(c)
        ax_a.set_ylim(-0.8, ymax + 0.8)

        for pair in pairs:
            group_ys = [ymax - pos[(pair, anc)] for anc in ANCESTRY_ORDER
                        if (pair, anc) in pos]
            y_top, y_bot = max(group_ys), min(group_ys)
            t = ax_a.annotate(PAIR_LABEL[pair], (-0.30, (y_top + y_bot) / 2),
                               xycoords=("axes fraction", "data"), ha="right", va="center",
                               fontsize=5.4, fontweight="bold", annotation_clip=False)
            vc.mark_label(t)

        note = ax_a.annotate(
            f"{N_PAIRS_TESTED} tests; Bonferroni p<{BONFERRONI_ALPHA:.4f}. * = nominal p<0.05\n"
            f"(uncorrected), not significant after correction. Dashed = OR 1. Hatched = a "
            f"cell <20, masked.",
            (0.0, 1.05), xycoords="axes fraction", ha="left", va="bottom", fontsize=4.6,
            color="#666666")
        vc.mark_label(note)
        vc.panel_letter(ax_a, "a")

        # --- Panel b: C1/C2/Bw4 carrier frequency by ancestry, compact dot plot ---
        epitopes = ["C1", "C2", "Bw4"]
        y_epi = {e: i for i, e in enumerate(epitopes)}
        for anc in ANCESTRY_ORDER:
            sub = epi[epi["ancestry"] == anc]
            if sub.empty:
                continue
            color = vc.ANCESTRY_COLORS.get(anc, vc.ACCENT_COLOR)
            for e in epitopes:
                r = sub[sub["epitope"] == e]
                if r.empty:
                    continue
                ax_b.scatter([float(r["pct"].iloc[0])], [y_epi[e]], s=10, color=color,
                             zorder=3, edgecolor="white", linewidth=0.3)
        ax_b.set_yticks(list(y_epi.values()))
        ax_b.set_yticklabels(epitopes, fontsize=5.5)
        ax_b.set_ylim(-0.6, len(epitopes) - 0.4)
        ax_b.set_xlim(0, 100)
        ax_b.set_xlabel("% carriers")
        ax_b.set_title("epitope carrier freq.\nby ancestry (dot color\nas in panel a)",
                        fontsize=5.4)
        vc.panel_letter(ax_b, "b")

        vc.save_fig(fig, out_stem)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", default="reports/hla_popgen/48_kir_hla_ligand_cooccurrence")
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args()
    out_dir = args.out_dir or args.in_dir
    os.makedirs(out_dir, exist_ok=True)

    cooc, epi = load_tables(args.in_dir)
    fig_forest_and_epitopes(cooc, epi, os.path.join(out_dir, "fig_kir_hla_ligand_forest"))
    print(f"[48b] 1 figure written to {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
