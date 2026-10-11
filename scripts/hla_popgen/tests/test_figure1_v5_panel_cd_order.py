#!/usr/bin/env python3
"""Regression test for a real correctness bug (orchestrator review, 2026-09-26): after giving
panel d's bar-chart axes `sharey=ax_c` (Figure 1 v5, `40_figure1_v5.py`) for structural row
alignment, the OLD code's `ax_bar.invert_yaxis(); ax_hm.invert_yaxis()` pair -- previously two
INDEPENDENT single inversions, one per (unshared) axes -- became two inversions of the SAME shared
view limits, which cancel out. Net effect verified against the pre-existing (separate-axes) code at
commit 64716c0: with `CLASSICAL = [A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1]` (index 0 = A), the OLD
render had `ax_hm`'s ylim ascending (A at bottom, DRB1 at top -- correct) but `ax_bar`'s ylim
DESCENDING (A at top, DRB1 at bottom) -- i.e. panel d's bars were already silently swapped
end-to-end relative to their own row labels in the pre-existing figure, not introduced by the
`sharey=` change. The current (shared-axes) code's net-zero double inversion happens to leave both
axes ascending -- gene labels and bar values correctly paired -- which this test locks in
mechanically rather than relying on eyeballing rendered pixel values against a TSV again.

Run: python3 -m pytest scripts/hla_popgen/tests/test_figure1_v5_panel_cd_order.py -q
"""
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)
if HLA_POPGEN_DIR not in sys.path:
    sys.path.insert(0, HLA_POPGEN_DIR)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

spec = importlib.util.spec_from_file_location(
    "figure1_v5", os.path.join(HLA_POPGEN_DIR, "40_figure1_v5.py"))
m40 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m40)


def _make_axes_like_compose_layout(fig):
    """Replicates compose_layout()'s ax_c/ax_d creation (sharey + private Ticker) exactly, so this
    test exercises the real wiring rather than a simplified stand-in."""
    ax_c = fig.add_subplot(1, 2, 1)
    ax_d = fig.add_subplot(1, 2, 2, sharey=ax_c)
    from matplotlib.axis import Ticker
    from matplotlib.ticker import NullFormatter, NullLocator
    ax_d.yaxis.major = Ticker()
    ax_d.yaxis.set_major_locator(NullLocator())
    ax_d.yaxis.set_major_formatter(NullFormatter())
    return ax_c, ax_d


def _distinct_grid_and_dtab():
    """Gene x ancestry values and per-gene bar totals, each a DISTINCT number so a row/column
    mixup is unambiguous (unlike using real data, where two genes' values could coincidentally
    collide)."""
    genes, ancs = m40.CLASSICAL, m40.ANC
    value = pd.DataFrame(
        {a: [100 * i + j for i in range(len(genes))] for j, a in enumerate(ancs)},
        index=genes, dtype=float)
    hatched = pd.DataFrame(False, index=genes, columns=ancs)
    text = value.round(0).astype(int).astype(str)
    grid = {"value": value, "hatched": hatched, "text": text}
    # Distinct, strictly increasing totals by CLASSICAL list position (A=10, B=20, ..., DRB1=80).
    dtab = pd.DataFrame({
        "seen once": [10 * (i + 1) for i in range(len(genes))],
        "2–19 unrelated people": [0] * len(genes),
    }, index=genes, dtype=float)
    return grid, dtab, genes, ancs


def test_panel_c_heatmap_value_matches_gene_row():
    """Each heatmap cell's rendered text must appear at the y-position whose tick label is that
    cell's own gene -- not a mirrored/reversed row."""
    grid, dtab, genes, ancs = _distinct_grid_and_dtab()
    fig = plt.figure(figsize=(6, 4))
    ax_c, ax_d = _make_axes_like_compose_layout(fig)
    m40.draw_panel_cd(ax_c, ax_d, grid, "any_field", dtab, cax=None)
    fig.canvas.draw()

    tick_gene_by_y = {}
    for loc, lab in zip(ax_c.get_yticks(), ax_c.get_yticklabels()):
        tick_gene_by_y[round(float(loc))] = lab.get_text()

    import matplotlib.text as mtext
    cell_texts = [t for t in ax_c.findobj(mtext.Text)
                 if t.get_text().lstrip("-").isdigit()]
    assert cell_texts, "no per-cell value text found on the heatmap"
    for t in cell_texts:
        x_data, y_data = t.get_position()
        y_idx = round(y_data)
        gene_bare = tick_gene_by_y[y_idx]
        gene = "HLA-" + gene_bare if not gene_bare.startswith("HLA-") else gene_bare
        expected = int(grid["value"].loc[gene].iloc[round(x_data)])
        assert int(t.get_text()) == expected, (
            f"row y={y_idx} labelled {gene_bare!r} shows {t.get_text()!r}, "
            f"expected {expected} (the TSV's own value for {gene!r})")
    plt.close(fig)


def test_panel_d_bar_total_matches_gene_row():
    """Each bar's total width (stacked recurrence-class segments) must equal the TSV's own total
    for the gene at that SAME y-position on panel c -- this is the exact fault: panel d's bars
    were reversed end-to-end relative to their row labels in the pre-`sharey` figure."""
    grid, dtab, genes, ancs = _distinct_grid_and_dtab()
    fig = plt.figure(figsize=(6, 4))
    ax_c, ax_d = _make_axes_like_compose_layout(fig)
    m40.draw_panel_cd(ax_c, ax_d, grid, "any_field", dtab, cax=None)
    fig.canvas.draw()

    tick_gene_by_y = {}
    for loc, lab in zip(ax_c.get_yticks(), ax_c.get_yticklabels()):
        tick_gene_by_y[round(float(loc))] = lab.get_text()

    from collections import defaultdict
    bar_total_by_y = defaultdict(float)
    for p in ax_d.patches:
        y_center = round(p.get_y() + p.get_height() / 2)
        bar_total_by_y[y_center] += p.get_width()

    assert set(bar_total_by_y) == set(tick_gene_by_y), (
        "panel d bar y-positions don't match panel c's tick y-positions at all -- "
        "sharey link is broken")
    for y_idx, gene_bare in tick_gene_by_y.items():
        gene = "HLA-" + gene_bare if not gene_bare.startswith("HLA-") else gene_bare
        expected = float(dtab.loc[gene, "seen once"] + dtab.loc[gene, "2–19 unrelated people"])
        got = bar_total_by_y[y_idx]
        assert got == pytest.approx(expected), (
            f"row y={y_idx} labelled {gene_bare!r}: bar total={got}, "
            f"expected {expected} (dtab's own total for {gene!r}) -- "
            f"panel d's bars are mismatched to panel c's gene labels")
    plt.close(fig)


def test_panel_c_and_d_agree_on_gene_order_top_to_bottom():
    """Sanity check on the intended visual order: DRB1 at the top, HLA-A at the bottom (matching
    the committed figure), for BOTH panels simultaneously."""
    grid, dtab, genes, ancs = _distinct_grid_and_dtab()
    fig = plt.figure(figsize=(6, 4))
    ax_c, ax_d = _make_axes_like_compose_layout(fig)
    m40.draw_panel_cd(ax_c, ax_d, grid, "any_field", dtab, cax=None)
    fig.canvas.draw()

    ylim_c = ax_c.get_ylim()
    ylim_d = ax_d.get_ylim()
    assert ylim_c == ylim_d, "panel c/d ylim diverged despite sharey="
    lo, hi = sorted(ylim_c)
    top_idx = round(hi) if ylim_c[0] < ylim_c[1] else round(lo)
    # Whichever tick is closest to the TOP of the drawn axes must be DRB1 (last in CLASSICAL).
    y_top_data = ylim_c[1]  # transform: matplotlib places ylim[1] at the top of the axes
    ticks = sorted(ax_c.get_yticks(), key=lambda v: abs(v - y_top_data))
    top_tick = ticks[0]
    top_label = dict(zip(ax_c.get_yticks(), [t.get_text() for t in ax_c.get_yticklabels()]))[top_tick]
    assert top_label == "DRB1", f"expected DRB1 at the top row, found {top_label!r}"
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
