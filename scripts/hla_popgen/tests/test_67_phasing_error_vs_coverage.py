#!/usr/bin/env python3
"""Synthetic tests for 67_phasing_error_vs_coverage.py (no VM data).
Run: python3 -m pytest scripts/hla_popgen/tests/test_67_phasing_error_vs_coverage.py -q"""
import importlib.util
import os
import sys
import unittest

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
spec = importlib.util.spec_from_file_location("m67", os.path.join(SRC, "67_phasing_error_vs_coverage.py"))
m = importlib.util.module_from_spec(spec)
sys.modules["m67"] = m
spec.loader.exec_module(m)


class T(unittest.TestCase):
    def test_rnd_and_rate_mask_small(self):
        self.assertEqual(m.rnd(0), "0")
        self.assertEqual(m.rnd(7), "<20")
        self.assertEqual(m.rnd(1234), "1230")
        self.assertEqual(m.rate(5, 100), ("", "", ""))
        self.assertEqual(m.rate(0, 100)[0], "0.0")
        self.assertEqual(m.rate(30, 100)[0], "30.0")

    def test_gene_depth_normalised_per_platform(self):
        d = pd.DataFrame({"person_id": ["a", "b", "c", "d"], "region": ["HLA-A"] * 4, "mapq": [0] * 4,
                          "rel_depth": [0.8, 1.2, 0.5, 1.5]})
        dep = m.gene_depth_index(d, {"a": "revio", "b": "revio", "c": "sequel2e", "d": "sequel2e"})
        self.assertAlmostEqual(dep[("a", "A")], 0.8)
        self.assertAlmostEqual(dep[("c", "A")], 0.5)

    def test_pair_records_filter_and_depth(self):
        geno = {}
        for p in ("i", "j"):
            for g in ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]:
                geno[(p, g)] = (("x",), ("x",), 6)
        geno[("j", "A")] = (("y",), ("y",), 3)          # one failing gene: pair kept, A is bad
        dep = {(p, g): 1.0 for p in ("i", "j") for g in ["A", "B"]}
        dep[("j", "A")] = 0.6
        fd = pd.DataFrame({"i": ["i"], "j": ["j"], "pp": ["revio-revio"]})
        r = m.pair_records(fd, geno, dep, {"G8": ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"],
                                           "ALL": ["A", "B"]})
        self.assertEqual(len(r), 2)
        a = r[r["gene"] == "A"].iloc[0]
        self.assertEqual((a["bad"], a["score"], a["dep_lo"], a["dep_hi"]), (1, 3, 0.6, 1.0))
        geno[("j", "B")] = (("y",), ("y",), 3)          # two failing genes: pair not MHC-sharing
        self.assertEqual(len(m.pair_records(fd, geno, dep, {"G8": ["A", "B"], "ALL": ["A", "B"]})), 0)

    def test_grid_masks_small_cells(self):
        rng = np.random.default_rng(1)
        lo = rng.uniform(0.5, 1.0, 600)
        r = pd.DataFrame({"dep_lo": lo, "dep_hi": lo + rng.uniform(0, 0.5, 600), "bad": (lo < 0.6).astype(int)})
        e = m.quantile_edges(pd.concat([r["dep_lo"], r["dep_hi"]]), 3)
        rows = pd.DataFrame(m.depth_grid(r, "dep_lo", "dep_hi", e, "pair"))
        self.assertTrue((rows["y_bin"] >= rows["x_bin"]).all())
        for _, x in rows.iterrows():
            if x["n_records"] == "<20":
                self.assertEqual(x["error_pct"], "")

    def test_switch_intervals(self):
        o = pd.DataFrame({"child": ["c"] * 3, "hap": ["hap1"] * 3, "gene": ["A", "B", "DRB1"],
                          "contig": ["ctg1"] * 3, "origin": ["pat", "mat", "mat"]})
        calls = pd.DataFrame({"person_id": ["c"] * 3, "hap": ["hap1"] * 3, "gene": ["A", "B", "DRB1"],
                              "contig": ["ctg1:1-500000"] * 3, "gene_start": [100, 200000, 400000],
                              "copy_index": [1, 1, 1]})
        pos = {"A": 29_900_000, "B": 31_350_000, "DRB1": 32_580_000}
        iv = m.switch_intervals(o, calls, pos)
        self.assertEqual(list(iv["switch"]), [1, 0])
        self.assertEqual(list(iv["interval"]), ["classI", "classI->classII_DRDQ"])

    def test_prune_curve_removes_worst_and_null(self):
        n = 400
        pp = pd.DataFrame({"person_id": [f"p{k}" for k in range(n)],
                           "norm_window_MHC_q0": np.linspace(0.5, 1.5, n),
                           "hla_absent": [1 if k < 40 else 0 for k in range(n)],
                           "ancestry": ["AFR"] * n, "platform": ["revio"] * n})
        out = m.prune_curve(pp, [("MHC", "norm_window_MHC_q0", False)], everyone=pp, draws=3)
        row10 = out[out["prune_pct"] == 10].iloc[0]
        self.assertEqual(row10["hla_absent_pct"], "0.00")     # the 40 events were the 40 lowest-depth people
        self.assertGreater(float(row10["hla_absent_pct_random"]), 5)


if __name__ == "__main__":
    unittest.main()
