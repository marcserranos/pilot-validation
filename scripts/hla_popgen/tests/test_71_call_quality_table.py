#!/usr/bin/env python3
"""Synthetic tests for 71_call_quality_table.py. Run: python3 -m pytest scripts/hla_popgen/tests/test_71_call_quality_table.py -q"""
import importlib.util
import os
import sys
import unittest

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
spec = importlib.util.spec_from_file_location("m71", os.path.join(SRC, "71_call_quality_table.py"))
m = importlib.util.module_from_spec(spec)
sys.modules["m71"] = m
spec.loader.exec_module(m)

FEATURES = pd.DataFrame([
    dict(name="A|gene", gene="A", kind="gene", length=1000, ars=False, exon_number=0, cds_ord=0),
    dict(name="A|exons1", gene="A", kind="exons", length=100, ars=False, exon_number=1, cds_ord=0),
    dict(name="A|exons2", gene="A", kind="exons", length=300, ars=False, exon_number=2, cds_ord=0),
    dict(name="A|cds1", gene="A", kind="cds", length=100, ars=False, exon_number=1, cds_ord=0),
    dict(name="A|cds2", gene="A", kind="cds", length=300, ars=False, exon_number=2, cds_ord=0),
    dict(name="A|ex1", gene="A", kind="ex", length=100, ars=False, exon_number=1, cds_ord=1),
    dict(name="A|ex2", gene="A", kind="ex", length=300, ars=True, exon_number=2, cds_ord=2),
])


def feature_depth(person, depths, mapq=0):
    return [dict(person_id=person, feature=n, mapq=mapq, sum_depth=d * FEATURES.set_index("name").loc[n, "length"],
                 n_reads=5, mean_depth=d) for n, d in depths.items()]


class T(unittest.TestCase):
    def setUp(self):
        d = {"A|gene": 20, "A|exons1": 10, "A|exons2": 30, "A|cds1": 10, "A|cds2": 30, "A|ex1": 10, "A|ex2": 30}
        rows = feature_depth("p1", d, 0) + feature_depth("p1", {k: v * 0.9 for k, v in d.items()}, 20)
        self.fd = pd.DataFrame(rows)
        self.ctrl = {("p1", 0): 20.0, ("p1", 20): 18.0}
        self.pb = pd.DataFrame([dict(person_id="p1", gene="A", mapq=0, cds_p10=5.0, cds_min=1.0, ars_min=25.0,
                                     cds_share_lt10=0.1, ars_share_lt10=0.0, cds_median=28.0, cds_mean=25.0)])

    def test_gene_depth_is_length_weighted_and_relative(self):
        g = m.gene_depth_table(self.fd, FEATURES, self.ctrl, self.pb).iloc[0]
        self.assertAlmostEqual(g["abs_exons"], (10 * 100 + 30 * 300) / 400)
        self.assertAlmostEqual(g["rel_cds"], 25.0 / 20.0)
        self.assertAlmostEqual(g["rel_ars"], 30.0 / 20.0)           # only the ARS exon
        self.assertAlmostEqual(g["rel_exon_min"], 10.0 / 20.0)      # shallowest coding exon
        self.assertAlmostEqual(g["mapq_frac"], 0.9)
        self.assertAlmostEqual(g["cds_p10_rel"], 5.0 / 20.0)

    def _calls(self):
        return pd.DataFrame([
            dict(person_id="p1", hap="1", contig="c1", gene="A", copy_index=1, consensus="HLA-A*01:01:01:01", tf="A*01:01",
                 plat="revio", template_distance=4, gene_start=1, gene_end=1000, cds_mut=None, template_warning="NA",
                 novelty_depth=np.nan, nctg=1, score=7, n_tied=1, cds_distance=np.nan),
            dict(person_id="p1", hap="2", contig="c2", gene="A", copy_index=1, consensus="HLA-A*02:01:01:01", tf="A*02:01",
                 plat="revio", template_distance=40, gene_start=1, gene_end=1000, cds_mut="X|:3*ag|K(AAG)<E(GAG)",
                 template_warning="partial_CDS", novelty_depth=2.0, nctg=3, score=2, n_tied=3, cds_distance=1),
        ])

    def test_call_table_orientation_and_na_semantics(self):
        gd = m.gene_depth_table(self.fd, FEATURES, self.ctrl, self.pb)
        t = m.build_call_table(self._calls(), gd, {"p1": 0.8}).set_index("hap")
        self.assertEqual(t.loc["1", "x_has_warning"], 0)             # the literal "NA" means no warning (quirk)
        self.assertEqual(t.loc["2", "x_has_warning"], 1)
        self.assertEqual(t.loc["1", "x_cds_distance"], 0)            # absent CDS distance = perfect template = 0
        self.assertAlmostEqual(t.loc["1", "x_template_distance_per_kb"], 4.0)
        self.assertEqual(t.loc["2", "x_is_novel"], 1)
        # oriented: depth and score flip sign so that higher is always worse
        self.assertAlmostEqual(t.loc["1", "o_rel_cds"], -1.25)
        self.assertEqual(t.loc["1", "o_score7"], -7)
        self.assertGreater(t.loc["2", "o_score7"], t.loc["1", "o_score7"])
        self.assertGreater(t.loc["2", "o_template_distance"], t.loc["1", "o_template_distance"])

    def test_trio_records(self):
        rows = []
        for p, tf in (("child", "A*01:01"), ("dad", "A*01:01"), ("mom", "A*03:01")):
            rows.append(dict(person_id=p, hap="1", contig="c", gene="A", consensus="HLA-" + tf + ":01", tf=tf, plat="revio",
                             **{"o_" + k: 1.0 for k in m.KEYS}))
        ct = pd.DataFrame(rows)
        origin = pd.DataFrame([dict(child="child", hap="1", gene="A", origin="pat"),
                               dict(child="child", hap="1", gene="A", origin="amb")])
        r = m.trio_records(origin, ct, {"child": ("dad", "mom")})
        self.assertEqual(len(r), 1)                                  # the ambiguous call is not a record
        self.assertEqual(r["bad"].iloc[0], 0)
        origin2 = pd.DataFrame([dict(child="child", hap="1", gene="A", origin="mat")])
        r2 = m.trio_records(origin2, ct, {"child": ("dad", "mom")})
        self.assertEqual(r2["bad"].iloc[0], 1)                       # mom lacks the child's allele
        self.assertIn("w_rel_cds", r2.columns)

    def test_pair_records_use_worse_person(self):
        base = {"o_" + k: 0.0 for k in m.KEYS}
        ct = pd.DataFrame([dict(person_id="i", hap="1", gene="A", **{**base, "o_rel_cds": -1.0}),
                           dict(person_id="j", hap="1", gene="A", **{**base, "o_rel_cds": -0.3})])
        geno = {("i", "A"): (("A*01:01",), ("c",), 5), ("j", "A"): (("A*02:01",), ("d",), 5)}
        pairs = pd.DataFrame({"i": ["i"], "j": ["j"]})
        r = m.pair_records(pairs, geno, ct, {"G8": ["A"], "ALL": ["A"]}, "twin")
        self.assertEqual(r["bad"].iloc[0], 1)
        self.assertAlmostEqual(r["w_rel_cds"].iloc[0], -0.3)          # the lower depth (higher oriented value) wins

    def test_hist_counts_sum_to_people(self):
        gd = pd.DataFrame({"person_id": [f"p{i}" for i in range(50)], "gene": "A",
                           **{f: np.linspace(0.2, 1.4, 50) for f in m.DEPTH_FEATURES}})
        h = m.hist_table(gd, {f"p{i}": "revio" for i in range(50)})
        s = h[(h["platform"] == "all") & (h["feature"] == "rel_cds")]["n"].sum()
        self.assertEqual(s, 50)
        self.assertEqual(len(h[(h["platform"] == "all") & (h["feature"] == "rel_cds")]), 200)   # 0.01-wide bins

    def test_metric_corr_symmetric_with_unit_diagonal(self):
        rng = np.random.default_rng(0)
        x = rng.normal(size=500)
        ct = pd.DataFrame({"o_" + k: x + (0 if i == 0 else rng.normal(size=500)) for i, k in enumerate(m.KEYS)})
        r = m.metric_corr(ct)
        d = r[r["metric_a"] == r["metric_b"]]["spearman"]
        self.assertTrue(np.allclose(d, 1.0))

    def test_priority_ids(self):
        pairs = pd.DataFrame({"i": ["a", "b"], "j": ["c", "d"]})
        self.assertEqual(m.priority_ids(pairs, [("k", "m1", "m2")]), ["a", "b", "c", "d", "k", "m1", "m2"])


if __name__ == "__main__":
    unittest.main()
