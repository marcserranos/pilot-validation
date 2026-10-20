#!/usr/bin/env python3
"""Smoke test: 72 draws every figure from synthetic aggregate tables and the layout linter accepts them.
Run: python3 -m pytest scripts/hla_popgen/tests/test_72_call_quality_figures.py -q"""
import importlib.util
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
sys.path.insert(0, HERE)
import make_s06_fixture as fx  # noqa: E402

spec = importlib.util.spec_from_file_location("m72", os.path.join(SRC, "72_call_quality_figures.py"))
m = importlib.util.module_from_spec(spec)
sys.modules["m72"] = m
spec.loader.exec_module(m)


class T(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = fx.build(tempfile.mkdtemp())
        cls.out = tempfile.mkdtemp()
        m.main(["--results", cls.res, "--out", cls.out])

    def test_all_figures_written(self):
        names = os.listdir(self.out)
        for stem in ("fig_depth_distributions", "fig_gene_depth_profile", "fig_exon_depth", "fig_error_curves_pair",
                     "fig_error_curves_trio_child", "fig_auc", "fig_metric_correlation", "fig_pruning",
                     "fig_grid_rel_cds__template_distance", "fig_contig_truth"):
            self.assertIn(stem + ".png", names, stem)
            self.assertIn(stem + ".pdf", names, stem)


if __name__ == "__main__":
    unittest.main()
