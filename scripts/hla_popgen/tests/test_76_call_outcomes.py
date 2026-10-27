#!/usr/bin/env python3
"""Run: python3 -m pytest scripts/hla_popgen/tests/test_76_call_outcomes.py -q   (76 and 76b, synthetic)"""
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


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(SRC, fn))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


m = load("m76", "76_call_outcomes.py")
b = load("m76b", "76b_pruning_evidence.py")
Q = m.Q


def call_rows(person, hap, gene, cons, ctg="c1", depth=1.0, plat="revio", **kw):
    r = dict(person_id=person, hap=hap, contig=ctg, gene=gene, consensus=cons, plat=plat)
    r["tf"] = m.C.two_field(cons)
    r["named"] = r["tf"] is not None
    for k in Q.KEYS:
        r["x_" + k] = 0.0
        r["o_" + k] = 0.0
    r["x_rel_cds"] = depth
    r["o_rel_cds"] = -depth
    r["x_rel_mhc"] = 0.9
    r["x_cds_share_lt10"] = 0.0
    r["x_cds_min_rel"] = depth
    r.update(kw)
    return r


class StreamA(unittest.TestCase):
    def test_slot_classes(self):
        rows = []
        # p1 hap1: A,B on ctg1, C missing but flanked by A and B on one contig -> inside_contig; hap2 absent entirely
        for g, c in (("A", "x1"), ("B", "x1"), ("DRB1", "x1")):
            rows.append(dict(person_id="p1", hap="hap1", contig=c + ":1-100", gene=g, consensus="HLA-%s*01:01:01:01" % g))
        # p2 hap1: A on ctg x1, B on ctg x2 (C between): between_contigs; hap2: only DPB1 (A has no flank on the left, anchors only on the right)
        rows += [dict(person_id="p2", hap="hap1", contig="x1:1-9", gene="A", consensus="HLA-A*01:01:01:01"),
                 dict(person_id="p2", hap="hap1", contig="x2:1-9", gene="B", consensus="HLA-B*07:02:01:01"),
                 dict(person_id="p2", hap="hap2", contig="y:1-9", gene="DPB1", consensus="HLA-DPB1*04:01:01:01")]
        s = m.slot_table(pd.DataFrame(rows), ["A", "C", "B"])
        g = s.set_index(["person", "gene", "hap"])
        self.assertEqual(g.loc[("p1", "C", "1"), "hap_class"], "inside_contig")
        self.assertEqual(g.loc[("p1", "C", "1"), "status"], "missing")
        self.assertEqual(g.loc[("p1", "A", "1"), "status"], "called")
        self.assertEqual(g.loc[("p1", "A", "2"), "hap_class"], "no_haplotype")
        self.assertEqual(g.loc[("p2", "C", "1"), "hap_class"], "between_contigs")
        self.assertEqual(g.loc[("p2", "A", "2"), "hap_class"], "region_end")

    def test_analyse_slots_runs(self):
        rng = np.random.default_rng(1)
        rows = []
        for i in range(300):
            for g in ("A", "B", "DRB3"):
                for h in ("1", "2"):
                    miss = rng.random() < (0.3 if i % 5 == 0 else 0.05)
                    rows.append(dict(person="p%d" % i, gene=g, hap=h, status="missing" if miss else "called",
                                     hap_class=rng.choice(["inside_contig", "between_contigs", "region_end", "no_haplotype"]),
                                     hap_n_genes=int(rng.integers(0, 17)), hap_n_contigs=int(rng.integers(1, 4)),
                                     x_rel_cds=float(rng.normal(0.9, 0.2)), x_rel_mhc=0.9, x_cds_min_rel=0.5, x_cds_share_lt10=0.0,
                                     other_called=int(rng.random() < 0.9), missing=int(miss)))
        with tempfile.TemporaryDirectory() as d:
            m.analyse_slots(pd.DataFrame(rows), d, boot=30)
            for f in ("A_missing_by_gene.tsv", "A_missing_by_hap_class.tsv", "A_missing_by_feature_bins.tsv", "A_missing_feature_auc.tsv", "A_missing_explained.tsv"):
                self.assertTrue(os.path.exists(os.path.join(d, f)), f)
            ex = pd.read_csv(os.path.join(d, "A_missing_explained.tsv"), sep="\t")
            self.assertAlmostEqual(ex[ex["scope"] == "classical8"]["pct_of_missing"].sum(), 100.0, places=3)


def trio_fixture():
    """c1 (child) with parents pa (pat), ma (mat); gene A: child hap1 from pat, hap2 from mat; more genes to trigger statuses."""
    calls = []

    def add(p, g, a1, a2=None, **kw):
        calls.append(call_rows(p, "hap1", g, "HLA-%s*%s:01:01:01" % (g, a1), **kw))
        if a2:
            calls.append(call_rows(p, "hap2", g, "HLA-%s*%s:01:01:01" % (g, a2), **kw))

    add("c1", "A", "01", "02")
    add("pa", "A", "01", "03")
    add("ma", "A", "02", "04")                # A: both inherited correctly
    add("c1", "B", "07", "08")
    add("pa", "B", "07", "09")
    add("ma", "B", "35", "44")                # B hap2 (mat): 08 absent from ma, absent from pa -> neither
    add("c1", "C", "01", "02")
    add("pa", "C", "01", "02")
    add("ma", "C", "05", "06")                # C hap2 (mat): 02 absent from ma but present in pa -> other_parent_only
    add("c1", "DRB1", "01")                   # child one call only -> child_incomplete
    add("pa", "DRB1", "01", "03")
    add("ma", "DRB1", "04", "07")
    add("c1", "DQB1", "02", "03")
    add("pa", "DQB1", "02")                   # parent one call -> parent_incomplete
    add("ma", "DQB1", "03", "05")
    ct = pd.DataFrame(calls)
    # a novel-in-field-2 call
    ct = pd.concat([ct, pd.DataFrame([call_rows("c1", "hap1", "DPB1", "HLA-DPB1*04:new"), call_rows("c1", "hap2", "DPB1", "HLA-DPB1*05:01:01:01"),
                                      call_rows("pa", "hap1", "DPB1", "HLA-DPB1*04:01:01:01"), call_rows("pa", "hap2", "DPB1", "HLA-DPB1*02:01:01:01"),
                                      call_rows("ma", "hap1", "DPB1", "HLA-DPB1*05:01:01:01"), call_rows("ma", "hap2", "DPB1", "HLA-DPB1*03:01:01:01")])],
                   ignore_index=True)
    og = []
    for g in ("A", "B", "C", "DRB1", "DQB1", "DPB1"):
        og.append(dict(child="c1", hap="hap1", gene=g, contig="k1", origin="pat", pat_k=40, mat_k=2))
        og.append(dict(child="c1", hap="hap2", gene=g, contig="k2", origin="mat", pat_k=3, mat_k=30))
    origin = pd.DataFrame(og)
    # DQB1 on hap1 is a switch suspect: neighbours DQA1 (next anchor) has the opposite origin
    origin = pd.concat([origin, pd.DataFrame([dict(child="c1", hap="hap1", gene="DQA1", contig="k1", origin="mat", pat_k=1, mat_k=9)])], ignore_index=True)
    return ct, origin


class StreamB(unittest.TestCase):
    def setUp(self):
        self.ct, self.origin = trio_fixture()
        rel = pd.DataFrame(dict(i=["c1", "c1", "pa", "x"], j=["pa", "ma", "ma", "y"], kin=[.25, .25, 0, .25], ibs0=np.nan, ibd0=[0.0, 0.01, 0.9, 0.0]))
        self.trio_of, self.meta = m.trio_meta([("c1", "pa", "ma"), ("c1", "pa", "zz")], rel)

    def test_trio_meta(self):
        r = self.meta.iloc[0]
        self.assertEqual(r["n_candidates"], 2)
        self.assertEqual(r["verified"], "yes")
        self.assertFalse(r["primary"])             # two candidate trios
        self.assertEqual(self.trio_of["c1"], ("pa", "ma"))
        bad = m.trio_meta([("c1", "pa", "ma")], pd.DataFrame(dict(i=["c1", "c1"], j=["pa", "ma"], kin=[.25, .25], ibs0=np.nan, ibd0=[0.0, 0.3])))[1]
        self.assertEqual(bad.iloc[0]["verified"], "no")
        none = m.trio_meta([("c1", "pa", "ma")], None)[1]
        self.assertEqual(none.iloc[0]["verified"], "unverifiable")

    def test_classification(self):
        rec = m.classify_trio_calls(self.ct, self.origin, self.trio_of, self.meta)
        r = rec.set_index(["gene", "hap"])
        self.assertEqual(r.loc[("A", "hap1"), "cls"], "match")
        self.assertEqual(r.loc[("A", "hap2"), "cls"], "match")
        self.assertEqual(r.loc[("B", "hap2"), "cls"], "neither")
        self.assertEqual(r.loc[("B", "hap2"), "bad"], 1)
        self.assertEqual(r.loc[("C", "hap2"), "cls"], "other_parent_only")
        self.assertEqual(r.loc[("C", "hap2"), "bad"], 0)
        self.assertEqual(r.loc[("DRB1", "hap1"), "status"], "child_incomplete")
        self.assertEqual(r.loc[("DQB1", "hap1"), "status"], "parent_incomplete")
        self.assertEqual(r.loc[("DPB1", "hap1"), "status"], "novel_unlabelled")

    def test_origin_robustness(self):
        rec = m.classify_trio_calls(self.ct, self.origin, self.trio_of, self.meta)
        r = rec.set_index(["gene", "hap"])
        self.assertEqual(r.loc[("A", "hap1"), "strong"], 1)       # 40 vs 2 k-mers
        self.assertEqual(r.loc[("A", "hap2"), "strong"], 1)       # 30 vs 3 k-mers
        self.assertEqual(r.loc[("DQB1", "hap1"), "switch_suspect"], 1)    # next gene on the same contig says mat
        self.assertEqual(r.loc[("A", "hap1"), "switch_suspect"], 0)

    def test_positions_switch(self):
        o = pd.DataFrame([dict(child="c", hap="hap1", gene="A", contig="k", origin="pat"), dict(child="c", hap="hap1", gene="B", contig="k", origin="mat"),
                          dict(child="c", hap="hap1", gene="DPB1", contig="k", origin="mat")])
        pos = {("c", "hap1", "A"): ("k", 1000), ("c", "hap1", "B"): ("k", 50_000), ("c", "hap1", "DPB1"): ("k", 5_000_000)}
        sus = m.switch_suspects(o, pos)
        self.assertEqual(sus, {("c", "hap1", "A"), ("c", "hap1", "B")})     # DPB1 is 5 Mb away from A: not a neighbour of it

    def test_summary_and_units(self):
        rec = m.classify_trio_calls(self.ct, self.origin, self.trio_of, self.meta)
        rec.loc[:, "primary"] = True
        with tempfile.TemporaryDirectory() as d:
            fam = m.summarise_trio(rec, d)
            cc = pd.read_csv(os.path.join(d, "B_class_counts.tsv"), sep="\t")
            row = cc[(cc["trio_set"] == "all_trios") & (cc["scope"] == "all_fixed_genes") & (cc["origin_subset"] == "default") & (cc["gene"] == "ALL")].iloc[0]
            self.assertEqual(row["n_neither"], 1)
            self.assertEqual(row["n_other_parent_only"], 1)
            self.assertEqual(fam["n_neither"].sum(), 1)
            for f in ("B_status_counts.tsv", "B_origin_margin.tsv", "B_family_clustering.tsv", "B_name_vs_cds.tsv", "B_by_child_homozygosity.tsv"):
                self.assertTrue(os.path.exists(os.path.join(d, f)), f)
        pgu = m.person_gene_units(self.ct, ["A"])
        comb = m.combined_units(rec, pgu)
        # c1: DRB1 child_incomplete -> missing=1; B -> wrong=1 ; A -> clean (bad 0)
        c = comb.set_index("gene")
        self.assertEqual(c.loc["DRB1", "missing"], 1)
        self.assertEqual(c.loc["B", "wrong"], 1)
        self.assertEqual(c.loc["A", "bad"], 0)
        self.assertNotIn("DQB1", c.index)          # parent incomplete and child complete: unlabelled, excluded

    def test_pair_units(self):
        d = pd.DataFrame(dict(person=["a", "b", "c"], other=["b", "a", "d"], gene=["A", "DRB3", "B"], bad=[0, 1, 1], w_n_called=[-2, -2, -1]))
        pa, pc = m.pair_units(d)
        self.assertEqual(len(pa), 2)                 # DRB3 dropped
        self.assertEqual(len(pc), 1)
        self.assertNotEqual(pa["family"].iloc[0], pa["family"].iloc[1])

    def test_three_field(self):
        self.assertEqual(m.three_field("HLA-A*01:01:01:02N"), "A*01:01:01")
        self.assertIsNone(m.three_field("HLA-A*01:01:new"))


class Pruning(unittest.TestCase):
    def test_removed_curve_and_ties(self):
        o = np.array([5.0, 4, 3, 2, 1, 0, 0, 0, 0, 0])
        y = np.array([1, 0, 1, 0, 0, 0, 0, 0, 0, 0])
        grid = np.array([0, 10, 20, 100.0])
        nr, er, n, K = b.removed_curve(o, y, grid)
        self.assertEqual(list(nr), [0, 1, 2, 10])
        self.assertEqual(list(er[:, 0]), [0, 1, 1, 2])
        # all ties: proportional removal
        nr, er, n, K = b.removed_curve(np.zeros(10), y, np.array([50.0]))
        self.assertAlmostEqual(er[0, 0], 1.0)

    def test_x_opt(self):
        grid = np.arange(0, 10.01, 0.25)
        nr = grid * 10                                     # n = 1000
        er = np.minimum(grid, 4) * 3.0                      # errors saturate at 4%
        self.assertLessEqual(b.x_opt_from(nr, er, grid, cost=25, window=2.0), 6.0)
        self.assertGreaterEqual(b.x_opt_from(nr, er, grid, cost=25, window=2.0), 3.5)

    def test_threshold_and_apply(self):
        o = np.arange(100.0)
        thr = b.threshold_at(o, 10)
        self.assertEqual(thr, 90.0)
        r = b.apply_threshold(o, (o >= 95).astype(int), thr)
        self.assertEqual(r["n_pruned"], 10)
        self.assertEqual(r["errors_removed"], 5)
        self.assertAlmostEqual(r["captured_pct"], 100.0)
        self.assertAlmostEqual(r["enrichment"], 10.0)

    def test_power(self):
        p = b.power_row("v", 3000, 60)
        self.assertGreater(p["auc_min_detectable"], 0.5)
        self.assertLess(p["auc_min_detectable"], 0.7)

    def test_end_to_end(self):
        rng = np.random.default_rng(5)
        n = 1500
        fam = rng.integers(0, 300, n)
        lat = rng.normal(size=n)
        bad = (rng.random(n) < 1 / (1 + np.exp(-(-4 + 1.5 * lat)))).astype(int)

        def unit(extra=None):
            d = pd.DataFrame(dict(person=["p%d" % i for i in range(n)], gene="A", bad=bad, family=fam))
            for k in Q.KEYS:
                d["o_" + k] = rng.normal(size=n)
            d["o_cds_share_lt10"] = lat + rng.normal(scale=0.5, size=n)
            return d
        with tempfile.TemporaryDirectory() as d:
            i, o = os.path.join(d, "i"), os.path.join(d, "o")
            os.makedirs(i)
            a = unit()
            a["scope_classical8"] = True
            a.to_csv(os.path.join(i, "units_A.INTERNAL.tsv"), sep="\t", index=False)
            bt = unit()
            for pre in ("c", "p", "w"):
                for k in Q.KEYS:
                    bt[pre + "_" + k] = bt["o_" + k]
            bt.to_csv(os.path.join(i, "units_B.INTERNAL.tsv"), sep="\t", index=False)
            pr = unit().rename(columns={"o_" + k: "w_" + k for k in Q.KEYS})
            pr["other"] = "q"
            pr.to_csv(os.path.join(i, "units_pair_complete.INTERNAL.tsv"), sep="\t", index=False)
            comb = unit()
            comb["missing"] = bad
            comb["wrong"] = 0
            comb.to_csv(os.path.join(i, "units_combined.INTERNAL.tsv"), sep="\t", index=False)
            b.run(i, o, grid_max=10.0, boot=25, folds=3, reps=3)
            for f in ("prune_curves.tsv", "prune_summary.tsv", "auc.tsv", "power.tsv", "cv_by_metric.tsv", "cv_selected.tsv", "threshold_stability.tsv", "transfer.tsv"):
                self.assertTrue(os.path.exists(os.path.join(o, f)), f)
            s = pd.read_csv(os.path.join(o, "prune_summary.tsv"), sep="\t")
            self.assertEqual(set(s["view"]), {"A_all", "A_classical8", "B_child", "B_worst", "combined_child", "pair_complete"})
            top = s[s["view"] == "A_all"].sort_values("auc", ascending=False).iloc[0]
            self.assertEqual(top["metric"], "cds_share_lt10")
            self.assertGreater(top["auc"], 0.7)
            cv = pd.read_csv(os.path.join(o, "cv_by_metric.tsv"), sep="\t")
            r = cv[(cv["view"] == "A_all") & (cv["metric"] == "cds_share_lt10")]
            self.assertGreater(r["enrichment"].iloc[0], 1.0)
            cur = pd.read_csv(os.path.join(o, "prune_curves.tsv"), sep="\t")
            cc = cur[(cur["view"] == "combined_child") & (cur["metric"] == "cds_share_lt10")]
            self.assertIn("captured_pct_missing", cc.columns)

    def test_small_event_views_do_not_crash(self):
        """Regression: with <5 training events per split cv_run produced empty frames and `groupby("metric")` raised KeyError."""
        rng = np.random.default_rng(2)
        n = 400
        def view(k_events):
            d = pd.DataFrame(dict(person=["p%d" % i for i in range(n)], gene="A", family=np.arange(n) % 50,
                                  bad=np.r_[np.ones(k_events, int), np.zeros(n - k_events, int)][rng.permutation(n)]))
            for k in Q.KEYS:
                d["o_" + k] = rng.normal(size=n)
            return d
        with tempfile.TemporaryDirectory() as d:
            i, o = os.path.join(d, "i"), os.path.join(d, "o")
            os.makedirs(i)
            for name, k in (("units_B", 6), ("units_B_all", 2)):
                x = view(k)
                for pre in ("c", "w"):
                    x = pd.concat([x, pd.DataFrame({pre + "_" + q: x["o_" + q] for q in Q.KEYS})], axis=1)
                x.to_csv(os.path.join(i, name + ".INTERNAL.tsv"), sep="\t", index=False)
            a = view(30)
            a["scope_classical8"] = True
            a.to_csv(os.path.join(i, "units_A.INTERNAL.tsv"), sep="\t", index=False)
            b.run(i, o, grid_max=5.0, boot=20, folds=5, reps=2)
            s = pd.read_csv(os.path.join(o, "prune_summary.tsv"), sep="\t")
            self.assertIn("B_child", set(s["view"]))
            self.assertNotIn("B_all_child", set(s["view"]))      # 2 events: skipped, logged, no crash
            pw = pd.read_csv(os.path.join(o, "power.tsv"), sep="\t")
            self.assertIn("B_all_child", set(pw["view"]))


if __name__ == "__main__":
    unittest.main()
