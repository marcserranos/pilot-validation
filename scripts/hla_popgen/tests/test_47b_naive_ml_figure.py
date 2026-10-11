#!/usr/bin/env python3
"""Unit tests for 47b_naive_ml_figure.py, synthetic data only.

Covers:
  1. `null_band()` -- the symmetric permutation-null band reconstruction from
     perm_auroc_mean/perm_auroc_p95, on hand-computable synthetic values (including the floor at
     0.5, since a mean very close to 0.5 with a slightly-above-mean p95 must never produce a
     lower edge below 0.5).
  2. Figure rendering on synthetic data through `_viz_common.check_layout(strict=True)` (via
     `save_fig`) -- must render with zero layout violations.
  3. **Value <-> label binding test**: renders panel (a)'s AUROC scatter+labels on a small
     synthetic `naive_ml_metrics.tsv`-shaped frame, reads back every `mark_label`-registered
     Annotation's anchor point (`xy`) and rounds its text back to a float, and asserts each one
     matches the corresponding task's `lr_auroc` value -- a regression guard against a
     task/AUROC mismatch (e.g. a sort-order bug between the scatter and the text loop).

Run: python3 scripts/hla_popgen/tests/test_47b_naive_ml_figure.py
"""
import importlib.util
import os
import sys
import tempfile
import unittest

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)


def _load_module(filename, modname):
    spec = importlib.util.spec_from_file_location(modname, os.path.join(HLA_POPGEN_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


m47b = _load_module("47b_naive_ml_figure.py", "naive_ml_figure_test_target")


def _synthetic_metrics():
    rows = [
        {"task": "ancestry_AFR", "n_total": 1000, "n_pos": 300, "lr_auroc": 0.900,
         "tree_auroc": 0.70, "n_perms": 20, "perm_auroc_mean": 0.500, "perm_auroc_p95": 0.520,
         "perm_pvalue_lr": 0.048},
        {"task": "ancestry_EUR", "n_total": 1000, "n_pos": 400, "lr_auroc": 0.850,
         "tree_auroc": 0.65, "n_perms": 20, "perm_auroc_mean": 0.501, "perm_auroc_p95": 0.515,
         "perm_pvalue_lr": 0.048},
        {"task": "platform_revio", "n_total": 1000, "n_pos": 900, "lr_auroc": 0.560,
         "tree_auroc": 0.51, "n_perms": 20, "perm_auroc_mean": 0.500, "perm_auroc_p95": 0.518,
         "perm_pvalue_lr": 0.048},
        {"task": "cB_from_HLA", "n_total": 1000, "n_pos": 600, "lr_auroc": 0.540,
         "tree_auroc": 0.52, "n_perms": 20, "perm_auroc_mean": 0.4995, "perm_auroc_p95": 0.513,
         "perm_pvalue_lr": 0.048},
    ]
    return pd.DataFrame(rows)


def _synthetic_features():
    rows = [
        {"task": "ancestry_AFR", "model": "l1_logreg", "feature": "HLA-B*07:05", "weight": -2.0},
        {"task": "ancestry_AFR", "model": "l1_logreg", "feature": "HLA-C*16:02", "weight": -1.5},
        {"task": "ancestry_EUR", "model": "l1_logreg", "feature": "HLA-DRB1*14:02", "weight": -3.0},
        {"task": "ancestry_EUR", "model": "l1_logreg", "feature": "HLA-B*39:10", "weight": -2.2},
    ]
    return pd.DataFrame(rows)


class TestNullBand(unittest.TestCase):
    def test_typical_band_fallback(self):
        # No perm_auroc_p05 column -> fallback path. mean 0.55, p95 0.60 -> mirrored lower edge
        # 2*0.55-0.60=0.50 (not floored, exactly at it).
        row = pd.Series({"perm_auroc_mean": 0.55, "perm_auroc_p95": 0.60})
        lo, hi, exact = m47b.null_band(row)
        self.assertFalse(exact)
        self.assertAlmostEqual(hi, 0.60)
        self.assertAlmostEqual(lo, 0.50)

        # mean 0.56, p95 0.60 -> mirrored lower edge 2*0.56-0.60=0.52, above the 0.5 floor
        row2 = pd.Series({"perm_auroc_mean": 0.56, "perm_auroc_p95": 0.60})
        lo2, hi2, exact2 = m47b.null_band(row2)
        self.assertFalse(exact2)
        self.assertAlmostEqual(hi2, 0.60)
        self.assertAlmostEqual(lo2, 0.52)

    def test_floors_at_point_five_fallback(self):
        # mean 0.50, p95 0.502 -> mirrored lower edge 0.498, must floor to 0.5.
        row = pd.Series({"perm_auroc_mean": 0.50, "perm_auroc_p95": 0.502})
        lo, hi, exact = m47b.null_band(row)
        self.assertFalse(exact)
        self.assertEqual(lo, 0.5)
        self.assertAlmostEqual(hi, 0.502)

    def test_p95_below_mean_still_gives_hi_ge_mean_fallback(self):
        # Defensive case: noisy small-sample p95 could in principle sit below the mean.
        row = pd.Series({"perm_auroc_mean": 0.51, "perm_auroc_p95": 0.505})
        lo, hi, exact = m47b.null_band(row)
        self.assertGreaterEqual(hi, row["perm_auroc_mean"])
        self.assertGreaterEqual(hi, lo)

    def test_exact_band_used_when_percentile_columns_present(self):
        # perm_auroc_p05 present and non-empty -> exact 5-95% band, ignoring the mean/mirroring.
        row = pd.Series({"perm_auroc_mean": 0.501, "perm_auroc_p05": 0.47, "perm_auroc_p95": 0.53})
        lo, hi, exact = m47b.null_band(row)
        self.assertTrue(exact)
        self.assertAlmostEqual(lo, 0.47)
        self.assertAlmostEqual(hi, 0.53)

    def test_exact_band_ignores_blank_percentile_column(self):
        # A blank string (as pandas.read_csv would give for an empty TSV cell) must not be
        # mistaken for a real 0.0 percentile -- falls back to the mean/p95 approximation.
        row = pd.Series({"perm_auroc_mean": 0.55, "perm_auroc_p05": "", "perm_auroc_p95": 0.60})
        lo, hi, exact = m47b.null_band(row)
        self.assertFalse(exact)


class TestWhiskers(unittest.TestCase):
    def test_returns_none_without_min_max_columns(self):
        row = pd.Series({"perm_auroc_mean": 0.5, "perm_auroc_p95": 0.52})
        self.assertIsNone(m47b.whiskers(row))

    def test_returns_min_max_when_present(self):
        row = pd.Series({"perm_auroc_min": 0.40, "perm_auroc_max": 0.62})
        self.assertEqual(m47b.whiskers(row), (0.40, 0.62))

    def test_blank_min_treated_as_absent(self):
        row = pd.Series({"perm_auroc_min": "", "perm_auroc_max": 0.62})
        self.assertIsNone(m47b.whiskers(row))


class TestFigureRendersAndLayout(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.metrics = _synthetic_metrics()
        self.feats = _synthetic_features()

    def test_fig_auroc_vs_null_renders(self):
        out_stem = os.path.join(self.tmp.name, "fig_naive_ml_auroc")
        m47b.fig_auroc_vs_null(self.metrics, self.feats, out_stem)
        self.assertTrue(os.path.exists(out_stem + ".png"))
        self.assertTrue(os.path.exists(out_stem + ".pdf"))


class TestLabelBindingMatchesData(unittest.TestCase):
    """The task's explicit ask: verify plotted value <-> label binding for at least one panel,
    against the source TSV/DataFrame -- not just "the figure rendered without an exception"."""

    def test_panel_a_auroc_labels_match_metrics_rows(self):
        metrics = _synthetic_metrics()
        m = metrics.set_index("task")
        tasks = [t for t in m47b.TASK_ORDER if t in set(metrics["task"])]

        # Monkeypatch save_fig to avoid writing/closing so we can inspect the axes afterwards --
        # _viz_common is already on sys.path/sys.modules from loading m47b above (it does its own
        # sys.path.insert(0, _THIS_DIR) at import time), so this is the same module instance.
        import _viz_common as vc
        real_save_fig = vc.save_fig
        captured = {}

        def _fake_save_fig(fig_, path_stem, **kwargs):
            captured["fig"] = fig_
            return path_stem + ".pdf", path_stem + ".png"

        vc.save_fig = _fake_save_fig
        try:
            out_stem = os.path.join(tempfile.mkdtemp(), "fig_naive_ml_auroc")
            m47b.fig_auroc_vs_null(metrics, _synthetic_features(), out_stem)
        finally:
            vc.save_fig = real_save_fig

        fig_built = captured["fig"]
        ax_a = fig_built.axes[0]
        found = {}
        for child in ax_a.get_children():
            if hasattr(child, "xy") and hasattr(child, "get_text"):
                txt = child.get_text()
                try:
                    val = float(txt)
                except ValueError:
                    continue
                x, y = child.xy
                found[round(x)] = (txt, y)

        for xi, t in enumerate(tasks):
            self.assertIn(xi, found, f"no numeric label found at x={xi} ({t})")
            txt, y = found[xi]
            expected_auroc = float(m.loc[t, "lr_auroc"])
            self.assertAlmostEqual(y, expected_auroc, places=6,
                                    msg=f"{t}: label y does not match metrics table")
            self.assertAlmostEqual(float(txt), expected_auroc, places=3,
                                    msg=f"{t}: label text does not match metrics table")


if __name__ == "__main__":
    unittest.main()
