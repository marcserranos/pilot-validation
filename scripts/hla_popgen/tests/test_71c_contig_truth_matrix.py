#!/usr/bin/env python3
"""Synthetic tests for 71c. Run: python3 -m pytest scripts/hla_popgen/tests/test_71c_contig_truth_matrix.py -q"""
import importlib.util
import os
import sys
import unittest

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
spec = importlib.util.spec_from_file_location("m71c", os.path.join(SRC, "71c_contig_truth_matrix.py"))
m = importlib.util.module_from_spec(spec)
sys.modules["m71c"] = m
spec.loader.exec_module(m)


def row(child, hap, gene, contig, origin):
    return dict(child=child, hap=hap, gene=gene, contig=contig, origin=origin)


class T(unittest.TestCase):
    def setUp(self):
        # child c1 hap1: A, B on ctg1 both from dad; DRB1 on ctg2 from dad (phase split by the assembly but truly dad)
        # child c1 hap2: A on ctg3 from mom, B on ctg3 from dad (a switch inside one contig), DRB1 ambiguous
        self.o = pd.DataFrame([
            row("c1", "hap1", "A", "ctg1", "pat"), row("c1", "hap1", "B", "ctg1", "pat"),
            row("c1", "hap1", "DRB1", "ctg2", "pat"),
            row("c1", "hap2", "A", "ctg3", "mat"), row("c1", "hap2", "B", "ctg3", "pat"),
            row("c1", "hap2", "DRB1", "ctg4", "amb"),
        ])

    def test_truth_counts_for_a_b(self):
        t = m.truth_table(self.o, genes=["A", "B", "DRB1"])
        r = t[(t["stratum"] == "all") & (t["gene_a"] == "A") & (t["gene_b"] == "B")].iloc[0]
        self.assertEqual(r["n_haplotypes"], 2)
        self.assertAlmostEqual(r["colocation_pct"], 100.0)           # both haplotypes have A and B on one contig
        self.assertEqual(r["n_resolved"], 2)
        self.assertAlmostEqual(r["same_parent_all_pct"], 50.0)       # hap1 same parent, hap2 switched
        self.assertAlmostEqual(r["same_parent_samecontig_pct"], 50.0)
        self.assertTrue(pd.isna(r["same_parent_diffcontig_pct"]))    # no different-contig pairs for A~B

    def test_split_contig_still_true_phase(self):
        t = m.truth_table(self.o, genes=["A", "B", "DRB1"])
        r = t[(t["stratum"] == "all") & (t["gene_a"] == "A") & (t["gene_b"] == "DRB1")].iloc[0]
        self.assertEqual(r["n_haplotypes"], 2)
        self.assertAlmostEqual(r["colocation_pct"], 0.0)
        self.assertEqual(r["n_resolved"], 1)                          # hap2's DRB1 is ambiguous
        self.assertAlmostEqual(r["same_parent_all_pct"], 100.0)
        self.assertAlmostEqual(r["same_parent_diffcontig_pct"], 100.0)
        self.assertTrue(pd.isna(r["same_parent_samecontig_pct"]))

    def test_matrix_is_symmetric(self):
        t = m.truth_table(self.o, genes=["A", "B", "DRB1"]).set_index(["stratum", "gene_a", "gene_b"])
        for (s, a, b) in t.index:
            x, y = t.loc[(s, a, b)], t.loc[(s, b, a)]
            self.assertEqual(x["n_resolved"], y["n_resolved"])

    def test_strata(self):
        t = m.truth_table(self.o, genes=["A", "B"], strata={"c1": "revio"})
        self.assertEqual(set(t["stratum"]), {"all", "revio"})

    def test_colocation_cohort(self):
        calls = pd.DataFrame([dict(person_id="p", hap="1", gene=g, contig=c, copy_index=1, consensus="x*01:01")
                              for g, c in (("A", "c1"), ("B", "c1"), ("DRB1", "c2"))]
                             + [dict(person_id="p", hap="2", gene="A", contig="c9", copy_index=1, consensus="undetermined"),
                                dict(person_id="p", hap="2", gene="B", contig="c9", copy_index=1, consensus="x*02:01")])
        hg = m.hap_gene_contigs(calls)
        t = m.colocation_table(hg, genes=["A", "B", "DRB1"])
        r = t[(t["gene_a"] == "A") & (t["gene_b"] == "B")].iloc[0]
        self.assertEqual(r["n_haplotypes"], 1)                        # hap2's A is undetermined, so not called
        self.assertAlmostEqual(r["colocation_pct"], 100.0)

    def test_gene_resolution(self):
        r = m.gene_resolution(self.o).set_index("gene")
        self.assertAlmostEqual(r.loc["DRB1", "resolved_pct"], 50.0)


if __name__ == "__main__":
    unittest.main()
