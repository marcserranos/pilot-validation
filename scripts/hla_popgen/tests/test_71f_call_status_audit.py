#!/usr/bin/env python3
"""Run: python3 -m pytest scripts/hla_popgen/tests/test_71f_call_status_audit.py -q"""
import importlib.util
import os
import sys
import unittest

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
spec = importlib.util.spec_from_file_location("m71f", os.path.join(SRC, "71f_call_status_audit.py"))
m = importlib.util.module_from_spec(spec)
sys.modules["m71f"] = m
spec.loader.exec_module(m)


class T(unittest.TestCase):
    def test_hap_status_and_taxonomy(self):
        rows = []
        # p1: A valid on both haps. p2: A valid on hap1, hap2 has other genes but no A. p3: A valid on hap1, A undetermined on hap2.
        # p4: only hap1 has any rows at all.
        for p, h, g, c in [("p1", "hap1", "HLA-A", "HLA-A*01:01:01:01"), ("p1", "hap2", "HLA-A", "HLA-A*02:01:01:01"),
                           ("p2", "hap1", "HLA-A", "HLA-A*01:01:01:01"), ("p2", "hap2", "HLA-B", "HLA-B*07:02:01:01"),
                           ("p3", "hap1", "HLA-A", "HLA-A*01:01:01:01"), ("p3", "hap2", "HLA-A", "undetermined"),
                           ("p4", "hap1", "HLA-A", "HLA-A*03:01:01:01")]:
            rows.append(dict(person_id=p, hap=h, gene=m.C.norm_gene(g), consensus=c))
        w = m.hap_status_table(pd.DataFrame(rows), ["A"]).set_index("person_id")
        self.assertEqual((w.loc["p1", "hap1"], w.loc["p1", "hap2"]), ("valid", "valid"))
        self.assertEqual(w.loc["p2", "hap2"], "no_gene_row")
        self.assertEqual(w.loc["p3", "hap2"], "undetermined")
        self.assertEqual(w.loc["p4", "hap2"], "no_haplotype")
        t = m.taxonomy(w.reset_index())
        r = t[(t["scope"] == "CLASSICAL8") & (t["status_pair"] == "valid + no_haplotype")]
        self.assertEqual(int(r["n_person_genes"].iloc[0]), 1)

    def test_trio_parent_status_counts_zero_calls(self):
        ct = pd.DataFrame([dict(person_id="c", hap="1", gene="A", tf="A*01:01", named=True),
                           dict(person_id="c", hap="2", gene="A", tf="A*02:01", named=True),
                           dict(person_id="mom", hap="1", gene="A", tf="A*02:01", named=True)])         # mom has ONE call; dad has none
        origin = pd.DataFrame([dict(child="c", hap="1", gene="A", origin="pat"), dict(child="c", hap="2", gene="A", origin="mat")])
        out = m.trio_parent_status(origin, ct, {"c": ("dad", "mom")})
        r0 = out[(out["scope"] == "ALL") & (out["parent_calls_at_gene"] == 0)]
        r1 = out[(out["scope"] == "ALL") & (out["parent_calls_at_gene"] == 1)]
        self.assertEqual(int(r0["child_calls"].iloc[0]), 1)
        self.assertEqual(int(r0["child_allele_missing_from_parent"].iloc[0]), 1)     # dad has no call: the child's allele is "missing" by construction
        self.assertEqual(int(r1["child_allele_missing_from_parent"].iloc[0]), 0)


if __name__ == "__main__":
    unittest.main()
