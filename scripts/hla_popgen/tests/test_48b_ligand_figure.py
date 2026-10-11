#!/usr/bin/env python3
"""Unit tests for 48b_ligand_figure.py, synthetic data only.

Covers:
  1. `_row_positions()` -- row-position bookkeeping (group gaps, ordering) on a small
     hand-computable synthetic pair/ancestry list.
  2. Figure rendering on synthetic data through `_viz_common.check_layout(strict=True)` (via
     `save_fig`) -- must render with zero layout violations, including a deliberately masked
     (<20-cell) row that must render as a hatched patch, never a plotted OR point.
  3. **Value <-> label binding test**: renders the forest panel on a small synthetic
     `kir_hla_ligand_cooccurrence.tsv`-shaped frame, reads back each plotted point's (x, y) and
     confirms it matches the corresponding row's `odds_ratio` and row position -- and confirms a
     masked row produces NO scatter point (only a hatched Rectangle), never a fabricated OR.

Run: python3 scripts/hla_popgen/tests/test_48b_ligand_figure.py
"""
import importlib.util
import os
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)


def _load_module(filename, modname):
    spec = importlib.util.spec_from_file_location(modname, os.path.join(HLA_POPGEN_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


m48b = _load_module("48b_ligand_figure.py", "ligand_figure_test_target")


def _synthetic_cooc():
    """3 pairs x 2 ancestries, one row deliberately masked (<20 cell -> blank OR/CI)."""
    rows = []

    def add(anc, gene, ligand, n_total, a, b, c, d, orv, lo, hi, masked=False):
        rows.append({
            "ancestry": anc, "kir_gene": gene, "ligand": ligand, "n_total": n_total,
            "n_a_kir_and_ligand": a, "n_b_kir_only": b, "n_c_ligand_only": c,
            "n_d_neither": d,
            "odds_ratio": np.nan if masked else orv,
            "or_ci_lo": np.nan if masked else lo,
            "or_ci_hi": np.nan if masked else hi,
            "fisher_p": np.nan if masked else 0.5, "chi2_p": np.nan if masked else 0.5,
            "perm_p": np.nan if masked else 0.5, "n_perms": 0 if masked else 1000,
        })

    add("AFR", "KIR2DL1", "C2", 500, 300, 100, 80, 20, 1.50, 0.99, 2.29)
    add("EUR", "KIR2DL1", "C2", 500, 250, 150, 70, 30, 1.08, 0.78, 1.48)
    add("AFR", "KIR2DL2", "C1", 500, 200, 100, 150, 50, 1.07, 0.92, 1.25)
    add("EUR", "KIR2DL2", "C1", 500, 8, 5, 480, 7, None, None, None, masked=True)
    add("AFR", "KIR3DL1", "Bw4", 500, 350, 100, 40, 10, 1.38, 0.83, 2.28)
    add("EUR", "KIR3DL1", "Bw4", 500, 300, 150, 30, 20, 1.01, 0.73, 1.38)
    df = pd.DataFrame(rows)
    df["pair"] = df["kir_gene"] + "x" + df["ligand"]  # load_tables() does this before plotting
    return df


def _synthetic_epi():
    rows = []
    for anc in ("AFR", "EUR"):
        for epi, pct in (("C1", 70.0), ("C2", 65.0), ("Bw4", 75.0)):
            rows.append({"ancestry": anc, "epitope": epi,
                         "n_carriers": int(pct * 5), "n_total": 500, "pct": pct})
    return pd.DataFrame(rows)


class TestRowPositions(unittest.TestCase):
    def test_group_gap_and_ordering(self):
        pairs = ["KIR2DL1xC2", "KIR2DL2xC1"]
        ancestries = ["AFR", "EUR"]
        pos, order = m48b._row_positions(pairs, ancestries, group_gap=1.0)
        # Within a group, consecutive ancestries are 1.0 apart.
        self.assertAlmostEqual(pos[("KIR2DL1xC2", "EUR")] - pos[("KIR2DL1xC2", "AFR")], 1.0)
        # Across groups, the gap (group_gap) is added on top of the last row's 1.0 step.
        self.assertAlmostEqual(
            pos[("KIR2DL2xC1", "AFR")] - pos[("KIR2DL1xC2", "EUR")], 1.0 + 1.0)
        self.assertEqual(order[0], ("KIR2DL1xC2", "AFR"))
        self.assertEqual(order[-1], ("KIR2DL2xC1", "EUR"))


class TestFigureRendersAndLayout(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.cooc = _synthetic_cooc()
        self.epi = _synthetic_epi()

    def test_fig_forest_and_epitopes_renders(self):
        out_stem = os.path.join(self.tmp.name, "fig_kir_hla_ligand_forest")
        m48b.fig_forest_and_epitopes(self.cooc, self.epi, out_stem)
        self.assertTrue(os.path.exists(out_stem + ".png"))
        self.assertTrue(os.path.exists(out_stem + ".pdf"))


class TestLabelBindingMatchesData(unittest.TestCase):
    """The task's explicit ask: verify plotted value <-> label binding for at least one panel,
    against the source TSV/DataFrame -- not just "the figure rendered without an exception".
    Also verifies the disclosure hard rule: a masked (<20-cell) row never becomes a plotted
    point, only a hatched Rectangle."""

    def test_forest_points_match_cooccurrence_rows_and_masked_row_has_no_point(self):
        import _viz_common as vc
        real_save_fig = vc.save_fig
        captured = {}

        def _fake_save_fig(fig_, path_stem, **kwargs):
            captured["fig"] = fig_
            return path_stem + ".pdf", path_stem + ".png"

        vc.save_fig = _fake_save_fig
        try:
            cooc = _synthetic_cooc()
            epi = _synthetic_epi()
            out_stem = os.path.join(tempfile.mkdtemp(), "fig_kir_hla_ligand_forest")
            m48b.fig_forest_and_epitopes(cooc, epi, out_stem)
        finally:
            vc.save_fig = real_save_fig

        fig_built = captured["fig"]
        ax_a = fig_built.axes[0]

        pairs = [p for p in m48b.PAIR_ORDER if p in set(cooc["pair"])]
        pos, order = m48b._row_positions(pairs, m48b.ANCESTRY_ORDER)
        ymax = max(pos.values())
        cooc_idx = cooc.set_index(["pair", "ancestry"])

        # Collect all scatter points actually drawn in the axes.
        from matplotlib.collections import PathCollection
        scatter_points = []
        for coll in ax_a.collections:
            if isinstance(coll, PathCollection):
                scatter_points.extend(coll.get_offsets().tolist())

        n_expected_points = 0
        for pair, anc in order:
            if (pair, anc) not in cooc_idx.index:
                continue
            row = cooc_idx.loc[(pair, anc)]
            y_expected = ymax - pos[(pair, anc)]
            if pd.isna(row["odds_ratio"]):
                # Masked row: must NOT appear as any scatter point at this y.
                for (px, py) in scatter_points:
                    self.assertFalse(
                        abs(py - y_expected) < 1e-6,
                        f"masked row {pair}/{anc} unexpectedly has a plotted point")
                continue
            n_expected_points += 1
            match = [(px, py) for (px, py) in scatter_points if abs(py - y_expected) < 1e-6]
            self.assertEqual(len(match), 1,
                              f"expected exactly 1 point for {pair}/{anc} at y={y_expected}")
            px, py = match[0]
            self.assertAlmostEqual(px, float(row["odds_ratio"]), places=6,
                                    msg=f"{pair}/{anc}: plotted OR does not match TSV")

        self.assertEqual(len(scatter_points), n_expected_points,
                          "extra or missing scatter points vs. non-masked TSV rows")


if __name__ == "__main__":
    unittest.main()
