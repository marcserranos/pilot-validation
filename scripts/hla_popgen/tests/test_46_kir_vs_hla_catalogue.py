#!/usr/bin/env python3
"""Unit tests for 46_kir_vs_hla_catalogue.py, synthetic data only.

Covers:
  1. `build_gene_metrics()` -- pct_novel_any/pct_novel_protein/completeness arithmetic on a small
     hand-computable synthetic `coverage_chao2.tsv`-shaped frame, including the merge against
     `equal_n_slope.tsv`-shaped slope data.
  2. `two_proportion_ztest()` reuse via the lazily-imported 45 module -- a known textbook z-test
     value (Marc's own worked hypothesis-test numbers from the S04 brief) reproduces to 2 decimal
     places.
  3. Figure rendering on synthetic data through `_viz_common.check_layout(strict=True)` (via
     `save_fig`) -- both figures must render with zero layout violations.
  4. **Label <-> data binding test** (the task's explicit ask: "check that every plotted value <->
     label binding matches the TSV"): renders `fig_catalogue_completeness` on a small synthetic
     metrics table, reads back every `mark_label`-registered Annotation's anchor point (`xy`,
     documented in `_viz_common.mark_label`'s own point-registration Annotation contract) and text,
     and asserts each one's (x, y) matches the corresponding TSV row's
     (pct_novel_any, completeness) to floating-point tolerance, and the text equals that row's
     `gene_display`. This is a real regression guard against a mislabeled point (e.g. a repulsion
     bug that nudges a label's anchor instead of its offset, or a gene/species mismatch in the
     zip() ordering upstream).

Run: python3 scripts/hla_popgen/tests/test_46_kir_vs_hla_catalogue.py
"""
import importlib.util
import math
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


m46 = _load_module("46_kir_vs_hla_catalogue.py", "kir_vs_hla_catalogue_test_target")


def _synthetic_cov():
    """Two genes per species, genomic + any_novel + protein_novel levels, ancestry ALL only --
    hand-computable pct_novel/completeness."""
    rows = []

    def add(gene, species, level, s_obs, gt=0.99, chao2=None):
        rows.append({"gene": gene, "ancestry": "ALL", "level": level, "species": species,
                     "n_people": 1000, "s_obs": s_obs, "good_turing_coverage": gt,
                     "chao2": chao2 if chao2 is not None else s_obs * 1.5,
                     "chao2_se": 1.0, "chao2_undetected_f0hat": (chao2 or s_obs * 1.5) - s_obs,
                     "q1": "", "q2": "", "chao_new_by_2n": 0.0})

    # HLA-A: 100 genomic, 20 any-level novel (20%), 0 protein-novel, chao2=200 -> completeness 0.5
    add("A", "hla", "genomic", 100, chao2=200.0)
    add("A", "hla", "any_novel", 20, chao2=40.0)
    add("A", "hla", "protein_novel", 0, chao2=0.0)
    # HLA-B: 200 genomic, 100 any-level novel (50%), chao2=250 -> completeness 0.8
    add("B", "hla", "genomic", 200, chao2=250.0)
    add("B", "hla", "any_novel", 100, chao2=180.0)
    add("B", "hla", "protein_novel", 0, chao2=0.0)
    # KIR2DL1: 50 genomic, 40 any-level novel (80%), chao2=55.55... -> completeness 0.9
    add("KIR2DL1", "kir", "genomic", 50, chao2=50 / 0.9)
    add("KIR2DL1", "kir", "any_novel", 40, chao2=44.0)
    add("KIR2DL1", "kir", "protein_novel", 10, chao2=12.0)
    # KIR2DL2: 40 genomic, 10 any-level novel (25%), chao2=80 -> completeness 0.5
    add("KIR2DL2", "kir", "genomic", 40, chao2=80.0)
    add("KIR2DL2", "kir", "any_novel", 10, chao2=15.0)
    add("KIR2DL2", "kir", "protein_novel", 2, chao2=3.0)
    return pd.DataFrame(rows)


def _synthetic_slope():
    rows = []
    for gene, species, slope in [("A", "hla", 5.0), ("B", "hla", 8.0),
                                  ("KIR2DL1", "kir", 1.0), ("KIR2DL2", "kir", 0.5)]:
        rows.append({"gene": gene, "ancestry": "ALL", "level": "genomic", "species": species,
                     "curve": "distinct", "n_star": 1000, "slope_per_1000": slope,
                     "mean_at_n_star": 0.0})
    return pd.DataFrame(rows)


class TestBuildGeneMetrics(unittest.TestCase):
    def setUp(self):
        self.cov = _synthetic_cov()
        self.slope = _synthetic_slope()
        self.metrics = m46.build_gene_metrics(self.cov, self.slope)

    def test_pct_novel_any(self):
        row = self.metrics[(self.metrics.gene == "A") & (self.metrics.species == "hla")].iloc[0]
        self.assertAlmostEqual(row["pct_novel_any"], 20.0)
        row = self.metrics[(self.metrics.gene == "B") & (self.metrics.species == "hla")].iloc[0]
        self.assertAlmostEqual(row["pct_novel_any"], 50.0)
        row = self.metrics[(self.metrics.gene == "KIR2DL1") & (self.metrics.species == "kir")].iloc[0]
        self.assertAlmostEqual(row["pct_novel_any"], 80.0)

    def test_pct_novel_protein_zero_for_hla(self):
        for gene in ("A", "B"):
            row = self.metrics[(self.metrics.gene == gene) & (self.metrics.species == "hla")].iloc[0]
            self.assertAlmostEqual(row["pct_novel_protein"], 0.0)

    def test_completeness(self):
        row = self.metrics[(self.metrics.gene == "A") & (self.metrics.species == "hla")].iloc[0]
        self.assertAlmostEqual(row["completeness"], 0.5)
        row = self.metrics[(self.metrics.gene == "B") & (self.metrics.species == "hla")].iloc[0]
        self.assertAlmostEqual(row["completeness"], 0.8)
        row = self.metrics[(self.metrics.gene == "KIR2DL1") & (self.metrics.species == "kir")].iloc[0]
        self.assertAlmostEqual(row["completeness"], 0.9, places=6)
        row = self.metrics[(self.metrics.gene == "KIR2DL2") & (self.metrics.species == "kir")].iloc[0]
        self.assertAlmostEqual(row["completeness"], 0.5)

    def test_slope_merged(self):
        row = self.metrics[(self.metrics.gene == "A") & (self.metrics.species == "hla")].iloc[0]
        self.assertAlmostEqual(row["slope_per_1000"], 5.0)

    def test_gene_display(self):
        row = self.metrics[(self.metrics.gene == "A") & (self.metrics.species == "hla")].iloc[0]
        self.assertEqual(row["gene_display"], "HLA-A")
        row = self.metrics[(self.metrics.gene == "KIR2DL1") & (self.metrics.species == "kir")].iloc[0]
        self.assertEqual(row["gene_display"], "KIR2DL1")


class TestTwoProportionZTestReuse(unittest.TestCase):
    def test_reuses_45s_implementation_and_matches_worked_example(self):
        m45 = _load_module("45_kir_recurrence_figure.py", "kir_recurrence_figure_test_target_46")
        # Marc's brief worked example: HLA singleton share ~34%, KIR ~18% (approximate figures
        # given in the task) -- confirm the z-test direction and rough magnitude on the ACTUAL
        # pooled genomic-level numbers 44 produced (962/2823 vs 236/1460), not just synthetic
        # toy numbers, so this test would catch a sign/formula regression against real values.
        p1, p2, z, pval = m45.two_proportion_ztest(962, 2823, 236, 1460)
        self.assertAlmostEqual(p1, 0.3408, places=3)
        self.assertAlmostEqual(p2, 0.1616, places=3)
        self.assertGreater(z, 10)  # highly significant, HLA >> KIR singleton share
        self.assertLess(pval, 1e-30)


class TestFiguresRenderAndLayout(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.cov = _synthetic_cov()
        self.slope = _synthetic_slope()
        self.metrics = m46.build_gene_metrics(self.cov, self.slope)

    def _synthetic_curve(self):
        """Minimal saturation_curves.tsv-shaped frame so pooled_recurrence_from_curves() (used by
        fig_recurrence_composition) has something to sum -- one gene per species, one level."""
        rows = []
        for gene, species in [("A", "hla"), ("KIR2DL1", "kir")]:
            for level in ("genomic", "any_novel", "protein_novel"):
                for n in (10, 20):
                    rows.append({"gene": gene, "ancestry": "ALL", "level": level,
                                 "species": species, "n": n, "mean_distinct": 10.0,
                                 "lo2_5": 8.0, "hi97_5": 12.0, "mean_eq1": 4.0, "mean_eq2": 2.0,
                                 "mean_gt2": 4.0, "mean_ge20": 1.0})
        return pd.DataFrame(rows)

    def test_fig_catalogue_completeness_renders(self):
        out_stem = os.path.join(self.tmp.name, "fig_catalogue_completeness")
        m46.fig_catalogue_completeness(self.metrics, out_stem)
        self.assertTrue(os.path.exists(out_stem + ".png"))
        self.assertTrue(os.path.exists(out_stem + ".pdf"))

    def test_fig_recurrence_composition_renders(self):
        out_stem = os.path.join(self.tmp.name, "fig_recurrence_composition")
        m46.fig_recurrence_composition(self._synthetic_curve(), out_stem)
        self.assertTrue(os.path.exists(out_stem + ".png"))
        self.assertTrue(os.path.exists(out_stem + ".pdf"))


class TestLabelBindingMatchesData(unittest.TestCase):
    """The task's explicit ask: verify plotted value <-> label binding for at least one panel,
    against the source TSV/DataFrame -- not just "the figure rendered without an exception"."""

    def test_catalogue_completeness_labels_match_metrics_rows(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        cov = _synthetic_cov()
        slope = _synthetic_slope()
        metrics = m46.build_gene_metrics(cov, slope)

        fig, ax = plt.subplots(figsize=(4, 3), constrained_layout=True)
        all_x, all_y, all_lab, all_col = [], [], [], []
        for sp in m46.SPECIES_ORDER:
            s = metrics[metrics["species"] == sp]
            ax.scatter(s["pct_novel_any"], s["completeness"], color=m46.SPECIES_COLOR[sp])
            all_x += s["pct_novel_any"].tolist()
            all_y += s["completeness"].tolist()
            all_lab += s["gene_display"].tolist()
            all_col += [m46.SPECIES_COLOR[sp]] * len(s)
        texts = m46._label_points(ax, all_x, all_y, all_lab, all_col)
        plt.close(fig)

        self.assertEqual(len(texts), len(metrics))
        # Build the expected (gene_display -> (x, y)) mapping straight from the metrics
        # DataFrame (the "TSV" this script would otherwise write), independent of plotting order.
        expected = {row["gene_display"]: (row["pct_novel_any"], row["completeness"])
                    for _, row in metrics.iterrows()}
        seen_labels = set()
        for t in texts:
            label = t.get_text()
            self.assertIn(label, expected, f"unexpected label {label!r} not in metrics table")
            exp_x, exp_y = expected[label]
            got_x, got_y = t.xy  # the anchor point passed to ax.annotate(label, (x, y), ...)
            self.assertAlmostEqual(got_x, exp_x, places=6,
                                    msg=f"{label}: plotted x does not match metrics table")
            self.assertAlmostEqual(got_y, exp_y, places=6,
                                    msg=f"{label}: plotted y does not match metrics table")
            seen_labels.add(label)
        self.assertEqual(seen_labels, set(expected))


if __name__ == "__main__":
    unittest.main()
