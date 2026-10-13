#!/usr/bin/env python3
"""Synthetic tests for 56d_score_vs_template_distance.py (no VM data).
Run: python3 -m pytest scripts/hla_popgen/tests/test_56d_score_vs_template_distance.py -q"""
import importlib.util
import os
import sys
import unittest

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
import _s05_common as C  # noqa: E402

spec = importlib.util.spec_from_file_location("m56d", os.path.join(SRC, "56d_score_vs_template_distance.py"))
m = importlib.util.module_from_spec(spec)
sys.modules["m56d"] = m
spec.loader.exec_module(m)


class T(unittest.TestCase):
    def test_call_distance_max_and_missing(self):
        cs = pd.DataFrame({"person_id": ["p", "p", "q", "q"], "gene": ["A", "A", "A", "A"],
                           "template_distance": ["1", "5", None, "2"]})
        d = m.call_distance(cs)
        self.assertEqual(d[("p", "A")], 5)
        self.assertEqual(d[("q", "A")], np.inf)

    def test_records_skip_ibd0_pairs_and_take_max(self):
        G8 = ["A", "B"]
        geno = {("i", "A"): (("A*01:01", "A*02:01"), (), 5), ("j", "A"): (("A*01:01", "A*03:01"), (), 3),
                ("i", "B"): (("B*07:02",), (), 6), ("j", "B"): (("B*08:01",), (), 6),
                ("x", "A"): (("A*01:01",), (), 6), ("y", "A"): (("A*02:01",), (), 6),
                ("x", "B"): (("B*07:02",), (), 6), ("y", "B"): (("B*08:01",), (), 6)}
        dist = {k: float(n) for n, k in enumerate(geno)}
        fd = pd.DataFrame({"i": ["i", "x"], "j": ["j", "y"]})
        r = m.records(fd, geno, dist, G8, G8)
        self.assertEqual(len(r), 2)  # pair x/y fails both genes -> dropped
        a = r[r["gene"] == "A"].iloc[0]
        self.assertEqual((a["bad"], a["score"], a["dist"]), (0, 3, 1.0))
        self.assertEqual(r[r["gene"] == "B"].iloc[0]["bad"], 1)

    def test_sweep_shapes_and_masking(self):
        rng = np.random.default_rng(0)
        n = 4000
        dist = rng.integers(0, 30, n).astype(float)
        bad = (rng.random(n) < 0.01 + 0.004 * dist).astype(int)
        r = pd.DataFrame({"gene": "A", "bad": bad, "score": rng.integers(0, 8, n), "dist": dist})
        s = m.sweep(r, C.rate_cells)
        self.assertEqual(set(s["filter"]), {"score_min", "template_distance_max"})
        allrow = s[(s["filter"] == "template_distance_max") & (s["threshold"] == "all")].iloc[0]
        self.assertEqual(allrow["retained_pct"], "100.0")
        d = s[s["filter"] == "template_distance_max"]
        ok = d[d["pct"] != ""]
        self.assertLess(float(ok.iloc[0]["pct"]), float(ok.iloc[-1]["pct"]))  # tighter distance -> fewer errors
        self.assertFalse(any(k.isdigit() and 0 < int(k) < 20 for k in s["k"]))

    def test_drop_differencing(self):
        rows = [{"_n": n, "_k": k, "id": i} for i, (n, k) in enumerate([(6455, 264), (6449, 261), (6389, 242),
                                                                        (5273, 104), (5023, 87), (40, 0)])]
        kept = {r["id"] for r in m.drop_differencing(rows)}
        self.assertNotIn(1, kept)  # 6 records / 3 errors from row 0
        self.assertNotIn(3, kept)  # 17 errors above row 4
        self.assertTrue({0, 2, 4, 5} <= kept)


if __name__ == "__main__":
    unittest.main()
