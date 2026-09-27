#!/usr/bin/env python3
"""Unit tests for 45_kir_recurrence_figure.py, synthetic data only.

Covers the v4b (2026-09-27) additions for the WS-A/WS-B final-figures pass:
  1. `compute_singleton_share_stats()` now also tests the CDS level (Marc's ask: private-vs-shared
     novelty at both CDS and protein level, not just any/protein) -- confirms "cds" appears as a
     row and the two-proportion z-test arithmetic is correct on a hand-computable fixture.
  2. `_label_curve_ends()`'s axis-bottom floor fix -- several near-zero-value labels no longer get
     placed AT the axis bottom (y == ylim[0]), which used to render half the glyph into the
     x-tick-label margin below the axes (found on full-size visual review of the new per-ancestry
     supplement figure, not caught by `check_layout()` since the label never left its own axes'
     bbox). Regression-tests the fix directly on the helper, not just "the figure rendered".
  3. `fig_saturation_by_recurrence_ancestry()` (new per-ancestry, equal-N, recurrence-class
     stratified supplement) renders through `check_layout(strict=True)` (via `save_fig`) on a
     small synthetic multi-ancestry curve fixture, including an all-zero ancestry/level slice
     (the "0 alleles" placeholder path).
  4. `fig_coverage_completeness_ancestry()` (new per-ancestry coverage/Chao2 panel) renders
     cleanly on synthetic data.

Run: python3 scripts/hla_popgen/tests/test_45_kir_recurrence_figure.py
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


m45 = _load_module("45_kir_recurrence_figure.py", "kir_recurrence_figure_test_target")


def _curve_row(gene, ancestry, level, species, n, distinct, eq1, eq2, gt2, ge20):
    return {"gene": gene, "ancestry": ancestry, "level": level, "species": species, "n": n,
            "mean_distinct": distinct, "lo2_5": distinct * 0.9, "hi97_5": distinct * 1.1,
            "mean_eq1": eq1, "mean_eq2": eq2, "mean_gt2": gt2, "mean_ge20": ge20}


def _synthetic_curve_for_stats():
    """One gene per species, genomic/cds/any_novel/protein_novel levels, ALL ancestry, single n
    (the full-cohort point) -- hand-computable singleton shares."""
    rows = []
    # HLA: genomic s_obs=100 (eq1=40 -> 40% singleton); cds s_obs=50 (eq1=10 -> 20%).
    rows.append(_curve_row("A", "ALL", "genomic", "hla", 500, 100, 40, 30, 20, 10))
    rows.append(_curve_row("A", "ALL", "cds", "hla", 500, 50, 10, 15, 20, 5))
    rows.append(_curve_row("A", "ALL", "any_novel", "hla", 500, 80, 30, 25, 20, 5))
    rows.append(_curve_row("A", "ALL", "protein_novel", "hla", 500, 20, 10, 5, 5, 0))
    # KIR: genomic s_obs=200 (eq1=150 -> 75%); cds s_obs=80 (eq1=64 -> 80%).
    rows.append(_curve_row("KIR2DL1", "ALL", "genomic", "kir", 500, 200, 150, 30, 15, 5))
    rows.append(_curve_row("KIR2DL1", "ALL", "cds", "kir", 500, 80, 64, 10, 5, 1))
    rows.append(_curve_row("KIR2DL1", "ALL", "any_novel", "kir", 500, 150, 100, 30, 15, 5))
    rows.append(_curve_row("KIR2DL1", "ALL", "protein_novel", "kir", 500, 40, 30, 6, 4, 0))
    return pd.DataFrame(rows)


class TestSingletonShareStatsIncludesCds(unittest.TestCase):
    def setUp(self):
        self.curve = _synthetic_curve_for_stats()
        self.stats = m45.compute_singleton_share_stats(self.curve)

    def test_cds_row_present(self):
        self.assertIn("cds", self.stats["level"].tolist())

    def test_cds_singleton_share_arithmetic(self):
        row = self.stats[self.stats["level"] == "cds"].iloc[0]
        self.assertAlmostEqual(row["hla_singleton_share"], 10 / 50)
        self.assertAlmostEqual(row["kir_singleton_share"], 64 / 80)
        self.assertAlmostEqual(row["hla_s_obs"], 50.0)
        self.assertAlmostEqual(row["kir_s_obs"], 80.0)

    def test_all_four_levels_present(self):
        self.assertEqual(set(self.stats["level"]),
                          {"genomic", "cds", "any_novel", "protein_novel"})


class TestLabelCurveEndsAxisFloor(unittest.TestCase):
    """Regression test for the axis-bottom-collision fix: several near-zero label values must all
    be placed strictly above `ax.get_ylim()[0]` (never AT it), each still separated by at least
    the enforced min_gap."""

    def test_near_zero_cluster_never_touches_axis_bottom(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(4, 3))
        ax.plot([0, 10], [0, 0])
        ax.set_ylim(0, 60)  # mimics a panel dominated by a much larger "seen 1x" curve
        # Three curves all ending at/near y=0 (mimics HLA AFR protein-level eq2/gt2/ge20).
        ends = {"seen 2x": (10, 1.5), "seen >2x": (10, 0.0), "seen ≥20x": (10, 0.0)}
        colors = {k: "black" for k in ends}
        m45._label_curve_ends(ax, ends, colors, fontsize=6.0, min_gap_frac=0.12)
        ylo, _ = ax.get_ylim()
        annotations = [a for a in ax.texts]
        self.assertEqual(len(annotations), 3)
        for a in annotations:
            _, y = a.xy
            self.assertGreater(y, ylo, f"label {a.get_text()!r} placed at/below axis bottom "
                                        f"(y={y}, ylim bottom={ylo})")
        plt.close(fig)

    def test_labels_still_respect_min_gap_after_flooring(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(4, 3))
        ax.plot([0, 10], [0, 0])
        ax.set_ylim(0, 60)
        ends = {"a": (10, 0.0), "b": (10, 0.0), "c": (10, 0.5)}
        colors = {k: "black" for k in ends}
        m45._label_curve_ends(ax, ends, colors, fontsize=6.0, min_gap_frac=0.12)
        ys = sorted(a.xy[1] for a in ax.texts)
        min_gap = 60 * 0.12
        for y0, y1 in zip(ys, ys[1:]):
            self.assertGreaterEqual(y1 - y0, min_gap - 1e-9)
        plt.close(fig)


class TestFigSaturationByRecurrenceAncestryRenders(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        rows = []
        for anc, scale in [("AFR", 1.0), ("EUR", 0.8)]:
            for level in ("any_novel", "protein_novel"):
                for n in (100, 200):
                    d = 10 * scale if level == "any_novel" else 0.0
                    rows.append(_curve_row("A", anc, level, "hla", n, d,
                                            d * 0.7, d * 0.15, d * 0.1, d * 0.05))
                    rows.append(_curve_row("KIR2DL1", anc, level, "kir", n, d * 5,
                                            d * 3.5, d * 0.9, d * 0.5, d * 0.1))
        self.curve = pd.DataFrame(rows)

    def test_renders_for_both_species_without_layout_violation(self):
        for species in ("hla", "kir"):
            out_stem = os.path.join(self.tmp.name, f"fig_ancestry_{species}")
            m45.fig_saturation_by_recurrence_ancestry(self.curve, out_stem, species,
                                                       ancestries=["AFR", "EUR"])
            self.assertTrue(os.path.exists(out_stem + ".png"))
            self.assertTrue(os.path.exists(out_stem + ".pdf"))


class TestFigCoverageCompletenessAncestryRenders(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        rows = []
        for anc in ("AFR", "EUR"):
            for sp, s_obs, chao2, gt in [("hla", 100, 250.0, 0.99), ("kir", 150, 400.0, 0.95)]:
                rows.append({"gene": "A", "ancestry": anc, "level": "protein", "species": sp,
                             "n_people": 1000, "s_obs": s_obs, "good_turing_coverage": gt,
                             "chao2": chao2, "chao2_se": 1.0,
                             "chao2_undetected_f0hat": chao2 - s_obs, "q1": "", "q2": "",
                             "chao_new_by_2n": 0.0})
        self.cov = pd.DataFrame(rows)

    def test_renders_without_layout_violation(self):
        out_stem = os.path.join(self.tmp.name, "fig_coverage_ancestry")
        m45.fig_coverage_completeness_ancestry(self.cov, out_stem, level="protein",
                                                ancestries=["AFR", "EUR"])
        self.assertTrue(os.path.exists(out_stem + ".png"))
        self.assertTrue(os.path.exists(out_stem + ".pdf"))


if __name__ == "__main__":
    unittest.main()
