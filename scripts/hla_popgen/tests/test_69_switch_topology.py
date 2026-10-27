"""69: switch calling on synthetic window origins, and PAF window mapping."""
import importlib.util
import os
import sys
import tempfile
import unittest

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
_spec = importlib.util.spec_from_file_location("_s69", os.path.join(os.path.dirname(HERE), "69_switch_topology.py"))
s = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s)


def win(seq, ctg="c1"):
    return pd.DataFrame([dict(child="k", plat="revio", hap="hap1", ctg=ctg, qa=i * 10_000 + 1, qb=(i + 1) * 10_000,
                              ref_pos=30_000_000 + i * 10_000, pat_k=0, mat_k=0, origin=o)
                         for i, o in enumerate(seq)])


class Test69(unittest.TestCase):
    def test_block_and_isolated(self):
        ev, adj = s.call_switches(win(["pat", "pat", "amb", "mat", "pat", "pat", "mat", "mat", "mat"]))
        self.assertEqual(list(ev.sort_values("pos")["kind"]), ["isolated", "block"])
        self.assertEqual(len(adj), 7)
        iso = ev[ev["kind"] == "isolated"].iloc[0]
        self.assertEqual(iso["pos"], 30_030_000)

    def test_no_switch(self):
        ev, _ = s.call_switches(win(["pat"] * 5 + ["amb"] * 3 + ["pat"] * 2))
        self.assertEqual(len(ev), 0)

    def test_paf_windows_minus_strand(self):
        tmp = tempfile.mkdtemp()
        paf = os.path.join(tmp, "x.paf")
        open(paf, "w").write("\t".join(["ctg", "200000", "0", "100000", "-", "chr6", "171000000", "30000000",
                                        "30100000", "99000", "100000", "60", "tp:A:P"]) + "\n")
        ws = s.windows_from_paf(paf, 10_000)
        self.assertEqual(len(ws), 10)
        first = min(ws, key=lambda w: w[1])
        self.assertGreater(first[3], 30_090_000)  # contig start maps to the right end on '-'


if __name__ == "__main__":
    unittest.main()
