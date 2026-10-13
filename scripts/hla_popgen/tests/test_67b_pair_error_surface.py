#!/usr/bin/env python3
"""Synthetic tests for 67b_pair_error_surface.py (no VM data).
Run: python3 -m pytest scripts/hla_popgen/tests/test_67b_pair_error_surface.py -q"""
import importlib.util
import os
import sys
import unittest

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location("m67b", os.path.join(SRC, "67b_pair_error_surface.py"))
m = importlib.util.module_from_spec(spec)
sys.modules["m67b"] = m
spec.loader.exec_module(m)


def synth(n=8000, seed=0):
    rng = np.random.default_rng(seed)
    a, b = rng.lognormal(0, 0.25, n), rng.lognormal(0, 0.25, n)
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    bad = (rng.random(n) < np.where(lo < 0.75, 0.12, 0.03)).astype(int)
    return pd.DataFrame({"dep_lo": lo, "dep_hi": hi, "bad": bad})


class T(unittest.TestCase):
    def test_surface_triangle_and_signal(self):
        s = m.surface(synth())
        self.assertTrue((s["dep_hi"] >= s["dep_lo"]).all())
        v = s[s["error_pct"] != ""].copy()
        v["e"] = v["error_pct"].astype(float)
        self.assertGreater(v[v["dep_lo"] < 0.7]["e"].mean(), v[v["dep_lo"] > 0.9]["e"].mean() + 3)

    def test_sparse_points_blank_and_eff_rounded(self):
        s = m.surface(synth())
        far = s[(s["dep_lo"] > 1.7)]
        self.assertTrue((far["error_pct"] == "").all())
        shown = s[s["eff_n"] != ""]
        self.assertTrue((shown["eff_n"].astype(int) % 10 == 0).all())
        self.assertTrue((shown["eff_n"].astype(int) >= m.MIN_EFF).all())


if __name__ == "__main__":
    unittest.main()
