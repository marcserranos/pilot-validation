#!/usr/bin/env python3
"""Synthetic tests for 75_trio_error_anatomy.py. Run: python3 -m pytest scripts/hla_popgen/tests/test_75_trio_error_anatomy.py -q"""
import importlib.util
import os
import sys
import unittest

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
spec = importlib.util.spec_from_file_location("m75", os.path.join(SRC, "75_trio_error_anatomy.py"))
m = importlib.util.module_from_spec(spec)
sys.modules["m75"] = m
spec.loader.exec_module(m)

HOMO = "HLA-A*01:01:01:01|:300+aa:100|K(AAG)<E(GAG)"          # one inserted a: homopolymer indel
OTHER = "HLA-A*01:01:01:01|:300+acg:100|K(AAG)<E(GAG)"          # mixed bases inserted
SUBS = "HLA-A*01:01:01:01|:300*ag:100|K(AAG)<E(GAG)"


def call(p, h, g, cons):
    return dict(person_id=p, hap=h, gene=g, consensus=cons, tf=cons)


class T(unittest.TestCase):
    def test_mut_class(self):
        self.assertEqual(m.mut_class(HOMO), "homopolymer_indel")
        self.assertEqual(m.mut_class(OTHER), "other_indel")
        self.assertEqual(m.mut_class(SUBS), "substitution")
        self.assertEqual(m.mut_class(None), "none")

    def test_closest_and_protein_novel(self):
        self.assertEqual(m.closest_two_field(HOMO), "A*01:01")
        self.assertTrue(m.protein_novel("HLA-A*01:new"))
        self.assertFalse(m.protein_novel("HLA-A*01:01:new"))      # novel only at depth 3: the 2-field name is intact
        self.assertFalse(m.protein_novel("HLA-A*01:01:01:02"))

    def _run(self, child_cons, child_mut, parent_cons):
        ct = pd.DataFrame([call("c", "hap1", "A", child_cons), call("c", "hap2", "A", "HLA-A*02:01:01:01"),
                           call("f", "hap1", "A", parent_cons[0]), call("f", "hap2", "A", parent_cons[1])])
        t1 = pd.DataFrame([dict(person_id="c", hap="hap1", gene="A", copy_index=1, cds_mut=child_mut)])
        calls = m.call_frame(ct, t1)
        rec = pd.DataFrame([dict(person="c", other="f", gene="A", hap="hap1", bad=0)])
        return m.evaluate(rec, calls)

    def test_depth3_novelty_is_not_an_error(self):
        defs, an, dn, bg, dis = self._run("HLA-A*01:01:new", None, ["HLA-A*01:01:01:01", "HLA-A*03:01:01:01"])
        self.assertEqual(int(defs.iloc[0]["n_errors"]), 0)

    def test_frameshift_novel_is_error_under_p0_and_collapsed_by_p1(self):
        defs, an, dn, bg, dis = self._run("HLA-A*01:new", HOMO, ["HLA-A*01:01:01:01", "HLA-A*03:01:01:01"])
        d = defs.set_index(defs["definition"].str[:2])["n_errors"]
        self.assertEqual(int(d["P0"]), 1)
        self.assertEqual(int(d["P1"]), 0)        # homopolymer-only: collapsed to A*01:01, which the parent carries
        self.assertEqual(int(d["P2"]), 0)
        self.assertEqual(an.iloc[0]["category"], "child protein-novel: homopolymer_indel")

    def test_substitution_novel_is_collapsed_only_by_p2(self):
        defs, an, dn, bg, dis = self._run("HLA-A*01:new", SUBS, ["HLA-A*01:01:01:01", "HLA-A*03:01:01:01"])
        d = defs.set_index(defs["definition"].str[:2])["n_errors"]
        self.assertEqual((int(d["P0"]), int(d["P1"]), int(d["P2"])), (1, 1, 0))

    def test_true_allele_difference_stays_an_error(self):
        defs, an, dn, bg, dis = self._run("HLA-A*01:01:01:01", None, ["HLA-A*03:01:01:01", "HLA-A*11:01:01:01"])
        d = defs.set_index(defs["definition"].str[:2])["n_errors"]
        self.assertEqual((int(d["P0"]), int(d["P1"]), int(d["P2"])), (1, 1, 1))
        self.assertEqual(an.iloc[0]["category"], "both named, different 2-field allele")


if __name__ == "__main__":
    unittest.main()
