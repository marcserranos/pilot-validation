#!/usr/bin/env python3
"""Synthetic tests for 71b_call_quality_stats.py. Run: python3 -m pytest scripts/hla_popgen/tests/test_71b_call_quality_stats.py -q"""
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
spec = importlib.util.spec_from_file_location("m71b", os.path.join(SRC, "71b_call_quality_stats.py"))
m = importlib.util.module_from_spec(spec)
sys.modules["m71b"] = m
spec.loader.exec_module(m)


def planted(n=4000, seed=1, beta=1.5):
    """Error probability rises with a latent 'badness' that the metric `o_rel_cds` measures with noise."""
    rng = np.random.default_rng(seed)
    bad_lat = rng.normal(size=n)
    p = 1 / (1 + np.exp(-(-3.5 + beta * bad_lat)))
    y = (rng.random(n) < p).astype(int)
    v = pd.DataFrame({"person": [f"p{i // 2}" for i in range(n)], "gene": rng.choice(["A", "B", "DRB1"], n),
                      "plat": rng.choice(["revio", "sequel2e"], n), "bad": y})
    for k in m.KEYS:
        v["o_" + k] = rng.normal(size=n)                        # pure noise metrics
    v["o_rel_cds"] = bad_lat + rng.normal(scale=0.5, size=n)    # informative
    v["o_has_warning"] = (bad_lat + rng.normal(size=n) > 1).astype(float)   # informative, binary
    return v


class T(unittest.TestCase):
    def test_auc_basics(self):
        self.assertAlmostEqual(m.auc_score([1, 2, 3, 4], [0, 0, 1, 1]), 1.0)
        self.assertAlmostEqual(m.auc_score([4, 3, 2, 1], [0, 0, 1, 1]), 0.0)
        self.assertAlmostEqual(m.auc_score([1, 1, 1, 1], [0, 0, 1, 1]), 0.5)      # ties
        self.assertTrue(np.isnan(m.auc_score([1, 2], [0, 0])))

    def test_informative_metric_beats_noise(self):
        v = planted()
        a_inf = m.auc_score(v["o_rel_cds"], v["bad"])
        a_noise = m.auc_score(v["o_n_tied"], v["bad"])
        self.assertGreater(a_inf, 0.75)
        self.assertLess(abs(a_noise - 0.5), 0.06)
        lo, hi = m.cluster_boot_auc(v["o_rel_cds"], v["bad"], v["person"], B=100)
        self.assertLess(lo, a_inf)
        self.assertGreater(hi, a_inf)

    def test_within_gene_auc_ignores_gene_confounding(self):
        rng = np.random.default_rng(2)
        n = 3000
        gene = rng.choice(["A", "B"], n)
        y = (rng.random(n) < np.where(gene == "A", 0.3, 0.02)).astype(int)
        o = np.where(gene == "A", 5.0, 0.0) + rng.normal(size=n)       # metric differs by gene, not by error
        self.assertGreater(m.auc_score(o, y), 0.7)                      # pooled AUC is fooled
        self.assertLess(abs(m.auc_within_gene(o, y, gene) - 0.5), 0.08)  # within-gene is not

    def test_logit_cluster_recovers_slope(self):
        rng = np.random.default_rng(3)
        n = 6000
        x = rng.normal(size=n)
        y = (rng.random(n) < 1 / (1 + np.exp(-(-1 + 0.8 * x)))).astype(int)
        b, se = m.logit_cluster(np.column_stack([np.ones(n), x]), y, np.arange(n))
        self.assertAlmostEqual(b[1], 0.8, delta=0.12)
        self.assertGreater(se[1], 0)

    def test_adjusted_or_direction(self):
        v = planted()
        orr, lo, hi, p = m.adjusted_or(v["o_rel_cds"], v["bad"], v["gene"], v["plat"], v["person"])
        self.assertGreater(orr, 1.5)
        self.assertLess(p, 1e-6)
        orn = m.adjusted_or(v["o_n_tied"], v["bad"], v["gene"], v["plat"], v["person"])
        self.assertTrue(orn[1] < 1.0 < orn[2] or orn[3] > 0.01)

    def test_prune_curve_improves_for_informative_and_is_flat_for_noise(self):
        v = planted()
        pc = m.prune_curve(v["o_rel_cds"], v["bad"])
        self.assertAlmostEqual(pc.loc[pc["prune_pct"] == 0, "error_pct_kept"].iloc[0], 100 * v["bad"].mean(), places=6)
        e10 = pc.loc[pc["prune_pct"] == 10, "error_pct_kept"].iloc[0]
        self.assertLess(e10, 100 * v["bad"].mean() * 0.8)
        self.assertGreater(pc.loc[pc["prune_pct"] == 10, "captured_pct"].iloc[0], 20)
        pn = m.prune_curve(v["o_n_tied"], v["bad"])
        r = pn.loc[pn["prune_pct"] == 10].iloc[0]
        self.assertTrue(r["random_lo"] <= r["error_pct_kept"] <= r["random_hi"] + 0.5)   # inside the random band
        # captured error share is monotone non-decreasing
        self.assertTrue((np.diff(pc["captured_pct"]) >= -1e-9).all())

    def test_prune_ties_removed_proportionally(self):
        o = np.r_[np.ones(100), np.zeros(900)]                  # a binary metric: 10% flagged
        y = np.r_[np.ones(50), np.zeros(50), np.zeros(900)]     # all errors are among the flagged
        pc = m.prune_curve(o, y, grid=np.array([0.0, 5.0, 10.0, 20.0]))
        r = pc.set_index("prune_pct")
        self.assertAlmostEqual(r.loc[5.0, "errors_kept"], 25.0)   # half of the flagged group removed, half its errors
        self.assertAlmostEqual(r.loc[10.0, "errors_kept"], 0.0)
        self.assertAlmostEqual(r.loc[20.0, "errors_kept"], 0.0)

    def test_cohort_pct_removed(self):
        o = np.arange(1000, dtype=float)
        y = (o > 900).astype(int)
        cohort = np.arange(10000, dtype=float)
        pc = m.prune_curve(o, y, grid=np.array([0.0, 10.0]), cohort_o=cohort).set_index("prune_pct")
        self.assertEqual(pc.loc[0.0, "cohort_pct_removed"], 0.0)
        self.assertAlmostEqual(pc.loc[10.0, "cohort_pct_removed"], 100 * (10000 - 900) / 10000, delta=0.2)

    def test_within_gene_percentile(self):
        ct = pd.DataFrame({"gene": ["A"] * 100 + ["B"] * 100, "o_rel_cds": np.r_[np.arange(100), np.arange(100) + 1000.0]})
        p = m.within_gene_percentile([99.0, 1099.0, 50.0], ["A", "B", "A"], ct, "rel_cds")
        self.assertGreater(p[0], 98)
        self.assertGreater(p[1], 98)
        self.assertAlmostEqual(p[2], 50.5, delta=1.0)

    def test_rolling_curve_discrete_and_continuous(self):
        v = planted()
        d = m.rolling_curve(v["o_has_warning"], v["bad"], "has_warning")
        self.assertEqual(len(d), 2)
        self.assertGreater(d["error_pct"].iloc[1], d["error_pct"].iloc[0])
        c = m.rolling_curve(v["o_rel_cds"], v["bad"], "rel_cds")
        self.assertGreater(len(c), 100)
        self.assertTrue((c["lo"] <= c["error_pct"] + 1e-9).all())
        # x is reported in the metric's own units (depth: sign flipped back), so error falls as x rises
        self.assertLess(np.corrcoef(c["x_median"], c["error_pct"])[0, 1], -0.5)

    def test_models_oof_auc(self):
        v = planted()
        X, names = m.design_matrix(v, m.MODEL_SETS["model:all"])
        r = m.oof_risk(X, v["bad"].to_numpy(), v["person"])
        self.assertTrue(np.isfinite(r).all())
        self.assertGreater(m.auc_score(r, v["bad"]), 0.75)
        cf = m.model_coefficients(X, names, v["bad"].to_numpy())
        self.assertGreater(cf.set_index("term").loc["rel_cds", "coef"], 0.5)

    def test_grid2d_counts_sum(self):
        v = planted(n=2000)
        ct = planted(n=5000, seed=9)
        ct["gene"] = "A"
        g = m.grid2d(v, "rel_cds", "n_tied", ct)
        self.assertEqual(g["n"].sum(), 2000)
        self.assertEqual(g["k"].sum(), int(v["bad"].sum()))

    def _cohort(self, n=6000, seed=5):
        c = planted(n=n, seed=seed)
        c["gene"] = np.random.default_rng(seed).choice(["A", "B", "DRB1"], n)
        return c

    def test_tail_composites_find_outliers_that_a_bulk_model_misses(self):
        # errors only happen in the far tail of one metric (worst 3%): a bulk-weighted linear model sees almost nothing, a tail feature sees it
        rng = np.random.default_rng(11)
        n = 8000
        v = pd.DataFrame({"person": [f"p{i}" for i in range(n)], "gene": "A", "plat": "revio"})
        for k in m.KEYS:
            v["o_" + k] = rng.normal(size=n)
        x = rng.normal(size=n)
        v["o_rel_cds"] = x
        v["bad"] = ((x > np.quantile(x, 0.97)) & (rng.random(n) < 0.5)).astype(int)
        ct = v.copy()
        pcts = m.metric_percentiles(v, ct, m.ALL_KEYS)
        self.assertGreater(m.auc_score(pcts["rel_cds"], v["bad"]), 0.9)
        Xt, names = m.tail_design(pcts, m.ALL_KEYS, v)
        risk = m.oof_risk(Xt, v["bad"].to_numpy(), v["person"])
        self.assertGreater(m.auc_score(risk, v["bad"]), 0.85)
        self.assertIn("rel_cds>=p95", names)
        pc = m.prune_curve(pcts["rel_cds"], v["bad"].to_numpy(), grid=np.array([0.0, 3.0]))
        self.assertGreater(pc.set_index("prune_pct").loc[3.0, "captured_pct"], 80)

    def test_gbm_oof_runs(self):
        v = planted(n=3000)
        X = np.column_stack([v["o_" + k] for k in m.ALL_KEYS])
        r = m.gbm_oof(X, v["bad"].to_numpy(), v["person"])
        self.assertGreater(m.auc_score(r, v["bad"]), 0.7)

    def test_value_hist_counts(self):
        v = planted(n=3000)
        h = m.value_hist(v, "n_tied")
        self.assertEqual(h["n"].sum(), 3000)
        self.assertEqual(h["k"].sum(), int(v["bad"].sum()))
        hc = m.value_hist(v, "rel_cds")
        self.assertEqual(hc["n"].sum(), 3000)
        self.assertTrue((hc["x_lo"] <= hc["x_hi"]).all())

    def test_complete_flag_and_strata(self):
        v = planted(n=2000)
        v["o_n_called"] = np.where(np.arange(2000) % 5 == 0, -1.0, -2.0)
        rec = pd.DataFrame({"person": v["person"], "gene": v["gene"], "bad": v["bad"], "plat": v["plat"], "w_n_called": v["o_n_called"]})
        for k in m.KEYS:
            rec["c_" + k] = v["o_" + k]
        mv = m.make_view(rec, "c")
        self.assertEqual(int(mv["complete"].sum()), 1600)

    def test_analyse_end_to_end_writes_tables(self):
        v = planted(n=1500)
        ct = planted(n=3000, seed=5)
        ct["gene"] = "A"
        d = tempfile.mkdtemp()
        m.analyse({"pair": v}, ct, d, do_boot=False)
        for f in ("stat_error_baseline.tsv", "stat_curves.tsv", "stat_auc.tsv", "stat_prune.tsv", "stat_model_coef.tsv", "stat_metric_hist.tsv",
                  "stat_complete_strata.tsv"):
            self.assertTrue(os.path.exists(os.path.join(d, f)), f)
        auc = pd.read_csv(os.path.join(d, "stat_auc.tsv"), sep="\t").set_index("metric")
        self.assertGreater(auc.loc["rel_cds", "auc"], 0.75)
        self.assertIn("model:all", auc.index)
        pr = pd.read_csv(os.path.join(d, "stat_prune.tsv"), sep="\t")
        self.assertEqual(set(pr["scale"]), {"pooled", "within_gene"})


if __name__ == "__main__":
    unittest.main()
