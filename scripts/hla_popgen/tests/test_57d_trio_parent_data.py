#!/usr/bin/env python3
"""Synthetic tests for 57d_trio_parent_data.py (no VM data).
Run: python3 -m pytest scripts/hla_popgen/tests/test_57d_trio_parent_data.py -q"""
import importlib.util
import os
import sys
import unittest

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
spec = importlib.util.spec_from_file_location("m57d", os.path.join(SRC, "57d_trio_parent_data.py"))
m = importlib.util.module_from_spec(spec)
sys.modules["m57d"] = m
spec.loader.exec_module(m)


def synth(n=60, seed=1):
    """n trios; parent depth drives the error (thin parent -> more error), child platform does not."""
    rng = np.random.default_rng(seed)
    trio, dep, tr, o = [], {}, [], []
    for i in range(n):
        c, a, b = f"c{i}", f"a{i}", f"b{i}"
        cp = "revio" if i % 2 else "sequel2e"
        pa, pb = ("revio", "revio") if i % 3 else ("revio", "sequel2e")
        da, db, dc = rng.uniform(8, 40), rng.uniform(8, 40), rng.uniform(15, 35)
        for p, d in ((a, da), (b, db), (c, dc)):
            dep[(p, "window_MHC")] = d
            dep[(p, "window_KIR")] = d
            dep[(p, "control")] = d
        trio.append({"set": "57", "child": c, "pa": a, "pb": b, "child_platform": cp, "pa_plat": pa, "pb_plat": pb,
                     "parents_mix": m.parents_mix(pa, pb), "n_nonrevio_parents": int(pa != "revio") + int(pb != "revio")})
        err = 0.02 * (10 / min(da, db)) ** 2
        for t in m.WINDOW:
            for h in ("1", "2"):
                tr.append({"child": c, "child_platform": cp, "target": t, "hap": h,
                           "sw_err": int(1000 * err), "sw_tot": 1000, "hm_err": int(5000 * err), "hm_tot": 5000})
        for g in range(10):  # deeper parent gets the call; thinner parent -> ambiguous more often
            amb = rng.random() < 0.5 * 10 / min(da, db)
            og = "amb" if amb else ("pat" if da > db else "mat")
            o.append({"child": c, "hap": "hap1", "gene": f"G{g}", "pat_k": int(da), "mat_k": int(db), "origin": og})
    return pd.DataFrame(trio), dep, pd.DataFrame(tr), pd.DataFrame(o)


class T(unittest.TestCase):
    def test_parents_mix(self):
        self.assertEqual(m.parents_mix("revio", "revio"), "both_revio")
        self.assertEqual(m.parents_mix("sequel2e", "revio"), "mixed_revio_sequel2e")
        self.assertEqual(m.parents_mix("sequel2e", "sequel2"), "any_sequel2")

    def test_depth_index(self):
        d = pd.DataFrame({"person_id": ["x", "x", "x"], "region": ["window_MHC", "window_KIR", "HLA-A"],
                          "mapq": [0, 0, 20], "mean_depth": [12.0, 30.0, 5.0], "control_depth": [20.0] * 3})
        dep = m.depth_index(d)
        self.assertEqual(dep[("x", "window_MHC")], 12.0)
        self.assertEqual(dep[("x", "control")], 20.0)
        self.assertNotIn(("x", "HLA-A"), dep)

    def test_region_table_finds_parent_effect(self):
        trio, dep, tr, _ = synth()
        r = m.region_table(m.per_trio_rates(tr), trio, dep)
        row = r[(r["target"] == "MHC") & (r["metric"] == "hamming_error")].iloc[0]
        self.assertLess(float(row["rho_vs_thinner_parent_depth"]), -0.8)
        self.assertLess(float(row["hi_per_doubling_thinner_parent_depth"]), 0)
        self.assertEqual(row["n_trios"], "60")

    def test_origin_table_lean_and_tertiles(self):
        trio, dep, _, o = synth(90)
        r = m.origin_table(o, trio, dep)
        lean = r[r["test"] == "resolved_calls_assigned_to_deeper_parent"].iloc[0]
        self.assertGreater(float(lean["value"]), 95)
        t = r[r["test"] == "ambiguous_by_thinner_parent_depth_tertile"].set_index("group")
        self.assertGreater(float(t.loc["tertile_1", "value"]), float(t.loc["tertile_3", "value"]))

    def test_mix_table_masks_small(self):
        trio, _, _, _ = synth(60)
        mt = m.mix_table(trio)
        self.assertTrue(set(mt["n_trios"]) <= {"0", "<20", "20"} | {str(x) for x in range(20, 61)})
        self.assertFalse(any(v.isdigit() and 0 < int(v) < 20 for v in mt["n_trios"]))


if __name__ == "__main__":
    unittest.main()
