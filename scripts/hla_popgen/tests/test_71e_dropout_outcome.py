#!/usr/bin/env python3
"""Run: python3 -m pytest scripts/hla_popgen/tests/test_71e_dropout_outcome.py -q"""
import importlib.util
import os
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
spec = importlib.util.spec_from_file_location("m71e", os.path.join(SRC, "71e_dropout_outcome.py"))
m = importlib.util.module_from_spec(spec)
sys.modules["m71e"] = m
spec.loader.exec_module(m)
Q = m.Q


def calls(n_people=1500, seed=3):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_people):
        for g in ("A", "B", "DRB1"):
            lat = rng.normal()
            two = rng.random() > 1 / (1 + np.exp(-(-2.4 + 1.3 * lat)))   # latent thinness makes a second call less likely
            for hap in (("1", "2") if two else ("1",)):
                r = dict(person_id=f"p{i}", hap=hap, gene=g, plat="revio")
                for k in Q.KEYS:
                    r["o_" + k] = rng.normal()
                r["o_rel_cds"] = lat + rng.normal(scale=0.4)
                r["o_cds_share_lt10"] = lat + rng.normal(scale=0.5)
                rows.append(r)
    return pd.DataFrame(rows)


class T(unittest.TestCase):
    def test_person_gene_table(self):
        ct = calls(100)
        pg = m.person_gene_table(ct)
        self.assertEqual(len(pg), 300)
        self.assertEqual(int(pg["incomplete"].sum()), int((pg["n_called"] < 2).sum()))
        self.assertTrue(set(pg["n_called"]) <= {1, 2})

    def test_depth_predicts_incompleteness_and_outputs_written(self):
        pg = m.person_gene_table(calls())
        d = tempfile.mkdtemp()
        m.analyse_dropout(pg, d)
        a = pd.read_csv(os.path.join(d, "dropout_auc.tsv"), sep="\t").set_index("metric")
        self.assertGreater(a.loc["rel_cds", "auc"], 0.7)
        self.assertLess(abs(a.loc["rel_exons", "auc"] - 0.5), 0.1)
        self.assertGreater(a.loc["rule:worst_rank_depth", "auc"], 0.65)
        for f in ("dropout_prune.tsv", "dropout_curves.tsv", "dropout_rate_by_gene.tsv"):
            self.assertTrue(os.path.exists(os.path.join(d, f)), f)
        pr = pd.read_csv(os.path.join(d, "dropout_prune.tsv"), sep="\t")
        r = pr[(pr["metric"] == "rel_cds") & (pr["prune_pct"] == 10.0)].iloc[0]
        base = pr[(pr["metric"] == "rel_cds") & (pr["prune_pct"] == 0.0)].iloc[0]["error_pct_kept"]
        self.assertLess(r["error_pct_kept"], base)


if __name__ == "__main__":
    unittest.main()
