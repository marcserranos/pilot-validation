#!/usr/bin/env python3
"""Run: python3 -m pytest scripts/hla_popgen/tests/test_71d_support_tables.py -q"""
import importlib.util
import os
import sys
import unittest

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
spec = importlib.util.spec_from_file_location("m71d", os.path.join(SRC, "71d_support_tables.py"))
m = importlib.util.module_from_spec(spec)
sys.modules["m71d"] = m
spec.loader.exec_module(m)


class T(unittest.TestCase):
    def test_overlap(self):
        pair = pd.DataFrame({"person": ["c", "c", "a"], "other": ["p1", "p2", "b"], "gene": ["A", "A", "A"], "bad": ["0", "0", "1"]})
        trio = pd.DataFrame({"person": ["c", "c"], "other": ["p1", "p2"], "gene": ["A", "A"], "hap": ["1", "2"], "bad": ["0", "0"]})
        o = dict(zip(*[m.overlap_table(pair, trio)[c] for c in ("quantity", "value")]))
        self.assertEqual(o["pairs_distinct"], 3)
        self.assertEqual(o["edges_in_both"], 2)
        self.assertEqual(o["pair_records_on_trio_edges"], 2)
        self.assertEqual(o["people_in_both"], 3)

    def test_named_effect(self):
        ct = pd.DataFrame({"person_id": ["c", "c", "d"], "hap": ["1", "2", "1"], "gene": ["A", "A", "A"], "named": ["False", "True", "True"]})
        trio = pd.DataFrame({"person": ["c", "c"], "other": ["d", "d"], "gene": ["A", "A"], "hap": ["1", "2"], "bad": ["1", "0"]})
        pair = pd.DataFrame({"person": ["c"], "other": ["d"], "gene": ["A"], "bad": ["1"]})
        e = dict(m.named_effect(ct, trio, pair, trio.iloc[1:], pair.iloc[0:0]).values)
        self.assertEqual(e["calls_unnamed"], 1)
        self.assertEqual(e["trio_any_unnamed_records"], 1)
        self.assertEqual(e["trio_any_unnamed_errors"], 1)
        self.assertEqual(e["trio_named_records"], 1)
        self.assertEqual(e["pair_any_unnamed_errors"], 1)


if __name__ == "__main__":
    unittest.main()
