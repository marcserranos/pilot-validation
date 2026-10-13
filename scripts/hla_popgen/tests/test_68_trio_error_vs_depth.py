"""68: synthetic checks of the prune curve, the 2D grid and the logistic table (no VM inputs)."""
import importlib.util
import os
import sys
import unittest

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
_spec = importlib.util.spec_from_file_location("_s68", os.path.join(os.path.dirname(HERE), "68_trio_error_vs_depth.py"))
s = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s)


def synth(n=4000, dep_effect=True, seed=1):
    r = np.random.default_rng(seed)
    dc, dp = r.normal(1, 0.15, n), r.normal(1, 0.15, n)
    lp = -3 + (8 * (1 - dc) if dep_effect else 0)
    y = (r.random(n) < 1 / (1 + np.exp(-lp))).astype(int)
    return pd.DataFrame(dict(child=[str(i // 20) for i in range(n)], d_child=dc, d_tpar=dp, mhc_child=dc,
                             mhc_tpar=dp, bad2=y, badfull=y, score=r.integers(0, 7, n)))


class Test68(unittest.TestCase):
    def test_prune_beats_random_only_with_effect(self):
        rng = np.random.default_rng(0)
        for eff, expect in ((True, True), (False, False)):
            d = synth(dep_effect=eff)
            rows = s.prune_curve(d, "bad2", "d_child", rng, "x", n_rand=200)
            r30 = [r for r in rows if r["prune_fraction"] == 0.3][0]
            self.assertEqual(r30["beats_random"], expect)

    def test_grid_and_logit(self):
        d = synth()
        g, edges = s.mendel_2d(d, "bad2")
        self.assertEqual(len(g), 9)
        self.assertFalse((g["error_pct"] == "").sum() == 1)
        M59 = s.load("m59", "59_depth_outcome_join.py")
        t = s.logit_table(d, M59)
        child = t[(t["metric"] == "bad2") & (t["depth"] == "gene_depth") & (t["term"] == "child")
                  & (~t["adjusted_for_score"])].iloc[0]
        par = t[(t["metric"] == "bad2") & (t["depth"] == "gene_depth") & (t["term"] == "transmitting_parent")
                & (~t["adjusted_for_score"])].iloc[0]
        self.assertGreater(child["or_per_0_1_lower"], 1.5)
        self.assertLess(abs(par["or_per_0_1_lower"] - 1), 0.2)

    def test_find_trios(self):
        pairs = [("c", "p", 0.25), ("c", "m", 0.25), ("c", "s", 0.25), ("p", "s", 0.25)]
        self.assertEqual(s.find_trios(["c"], pairs)["c"], ("m", "p"))


if __name__ == "__main__":
    unittest.main()
