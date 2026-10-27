#!/usr/bin/env python3
"""Synthetic tests for 79 (final S06 report): full build, pending placeholders, masking, think-tank placeholder. python3 tests/test_79_s06_final_report.py"""
import importlib.util
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
spec = importlib.util.spec_from_file_location("m79", os.path.join(SRC, "79_s06_final_report.py"))
M = importlib.util.module_from_spec(spec)
sys.modules["m79"] = M
spec.loader.exec_module(M)

G8 = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]


def w(root, sub, name, df):
    os.makedirs(os.path.join(root, sub), exist_ok=True)
    df.to_csv(os.path.join(root, sub, name), sep="\t", index=False)


def make_tables(root, with_prune=True):
    rng = np.random.default_rng(1)
    rows = []
    for sc in ("classical8", "all_fixed_genes"):
        for g in G8:
            n, k = 2000, int(rng.integers(40, 120))
            rows.append(dict(scope=sc, gene=g, n=n, k=k, pct=100 * k / n, lo=100 * k / n - 0.4, hi=100 * k / n + 0.4))
    w(root, "76", "A_missing_by_gene.tsv", pd.DataFrame(rows))
    hc = []
    for sc in ("classical8", "all_fixed_genes"):
        for c, n, p in (("inside_contig", 100000, 0.03), ("between_contigs", 40000, 10.5), ("region_end", 50000, 5.0), ("no_haplotype", 100, 100.0)):
            hc.append(dict(scope=sc, hap_class=c, n=n, k=int(n * p / 100), pct=p, lo=p, hi=p))
    w(root, "76", "A_missing_by_hap_class.tsv", pd.DataFrame(hc))
    ex = []
    for sc in ("classical8", "all_fixed_genes"):
        for c, p in (("assembly_break_at_gene", 56.6), ("region_end", 41.6), ("haplotype_absent", 1.4), ("inside_contig_depth_ok", 0.3), ("inside_contig_low_depth", 0.1)):
            ex.append(dict(scope=sc, explanation=c, n_missing_slots=int(p * 100), pct_of_missing=p))
    w(root, "76", "A_missing_explained.tsv", pd.DataFrame(ex))
    fb = []
    for f in ("hap_n_genes", "hap_n_contigs", "x_rel_cds", "x_rel_mhc", "x_cds_min_rel", "x_cds_share_lt10"):
        for i in range(10):
            fb.append(dict(scope="classical8", feature=f, bin_lo=i * 0.1, bin_hi=i * 0.1 + 0.1, n=100, k=10 + i, pct=10 + i, lo=8 + i, hi=12 + i))
    w(root, "76", "A_missing_by_feature_bins.tsv", pd.DataFrame(fb))
    fa = [dict(scope=sc, feature=f, n=1000, n_missing=100, auc=a, lo=a - 0.02, hi=a + 0.02) for sc in ("classical8", "all_fixed_genes")
          for f, a in (("hap_n_genes", 0.87), ("hap_n_contigs", 0.69), ("x_rel_cds", 0.6))]
    w(root, "76", "A_missing_feature_auc.tsv", pd.DataFrame(fa))
    cc = []
    for ts in ("primary", "all_trios"):
        for os_ in ("default", "strong_origin", "no_switch_suspect", "strong_and_no_switch"):
            cc.append(dict(trio_set=ts, scope="classical8", origin_subset=os_, gene="ALL", n_eligible=1800, n_match=1790, n_other_parent_only=3, n_neither=7, pct=0.39, lo=0.2, hi=0.8))
        for g in G8:
            cc.append(dict(trio_set=ts, scope="classical8", origin_subset="default", gene=g, n_eligible=230, n_match=228, n_other_parent_only=1, n_neither=1, pct=0.4, lo=0.07, hi=2.4))
    w(root, "76", "B_class_counts.tsv", pd.DataFrame(cc))
    w(root, "76", "B_family_event_hist.tsv", pd.DataFrame(dict(neither_events_in_family=[0, 1], families=[130, 6])))
    om = [dict(cls=c, feature=f, bin_lo=i, bin_hi=i + 1, n=int(rng.integers(0, 50))) for c in ("match", "other_parent_only", "neither") for f in ("log2_kmer_ratio", "winner_kmers") for i in range(10)]
    w(root, "76", "B_origin_margin.tsv", pd.DataFrame(om))
    w(root, "76", "B_family_clustering.tsv", pd.DataFrame([dict(families=136, n_neither=7)]))
    w(root, "76", "C_novel_rates_by_gene.tsv", pd.DataFrame([dict(gene=g, n_calls=1000, novel_any=60, novel_unnamed_field12=20, novel_any_pct=6, novel_any_lo=5, novel_any_hi=7, unnamed_pct=2,
                                                                  unnamed_lo=1.5, unnamed_hi=2.5) for g in G8]))
    w(root, "76", "C_novel_by_depth.tsv", pd.DataFrame([dict(outcome=o, feature="x_rel_cds", bin_lo=i * 0.1, bin_hi=i * 0.1 + 0.1, n=50, k=5, pct=10, lo=5, hi=15) for o in ("novel_any", "novel_unnamed") for i in range(10)]))
    w(root, "76", "definition_card.tsv", pd.DataFrame([dict(stream="A", term="slot", definition="person x gene x haplotype", n=100)]))
    w(root, "71f", "trio_parent_status.tsv", pd.DataFrame([
        dict(scope="CLASSICAL8", parent_calls_at_gene=1, child_call_named=True, child_calls=76, child_allele_missing_from_parent=32, missing_but_in_other_parent=2, missing_pct=42.1),
        dict(scope="CLASSICAL8", parent_calls_at_gene=2, child_call_named=True, child_calls=2161, child_allele_missing_from_parent=36, missing_but_in_other_parent=4, missing_pct=1.67),
        dict(scope="ALL", parent_calls_at_gene=1, child_call_named=True, child_calls=261, child_allele_missing_from_parent=46, missing_but_in_other_parent=5, missing_pct=17.6),
        dict(scope="ALL", parent_calls_at_gene=2, child_call_named=True, child_calls=3247, child_allele_missing_from_parent=54, missing_but_in_other_parent=9, missing_pct=1.66)]))
    w(root, "71f", "trio_parent_status_v4.tsv", pd.DataFrame([
        dict(scope="CLASSICAL8", parent_calls_at_gene=1, child_calls_at_gene=1, n_records=3, n_child_allele_absent_from_parent=1, pct=33.3),
        dict(scope="CLASSICAL8", parent_calls_at_gene=1, child_calls_at_gene=2, n_records=71, n_child_allele_absent_from_parent=29, pct=40.8),
        dict(scope="CLASSICAL8", parent_calls_at_gene=2, child_calls_at_gene=1, n_records=59, n_child_allele_absent_from_parent=0, pct=0.0),
        dict(scope="CLASSICAL8", parent_calls_at_gene=2, child_calls_at_gene=2, n_records=2025, n_child_allele_absent_from_parent=8, pct=0.4)]))
    w(root, "77c", "kir_label_instability.tsv", pd.DataFrame([dict(stratum="tx_unnamed", value="1", n_records=314, n_neither=85, rate=0.27), dict(stratum="tx_unnamed", value="0", n_records=4572, n_neither=43, rate=0.0094),
                                                              dict(stratum="tx_n_calls", value="1", n_records=866, n_neither=56, rate=0.065), dict(stratum="tx_n_calls", value="2", n_records=3827, n_neither=70, rate=0.018)]))
    w(root, "71f", "call_status_taxonomy.tsv", pd.DataFrame([dict(scope="gene", gene=g, status_pair=s, n_person_genes=100, pct=p) for g in G8
                                                              for s, p in (("valid + valid", 92.0), ("valid + no_gene_row", 7.0), ("no_gene_row + no_gene_row", 1.0))]))
    genes = ["KIR2DL1", "KIR3DL3", "KIR2DS4"]
    w(root, "77", "hap_class_counts.tsv", pd.DataFrame(dict(hap_class=["complete", "empty", "multi_contig", "one_contig_open_end"], n_haps=[18000, 100, 4000, 1700])))
    w(root, "77", "gene_n_called.tsv", pd.DataFrame([dict(gene=g, n_called=c, n_people=n) for g in genes for c, n in ((1, 300), (2, 700))]))
    w(root, "77", "gene_by_novelty_tier.tsv", pd.DataFrame([dict(gene=g, novelty_tier=t, n=n) for g in genes for t, n in (("known", 100), ("novel_protein", 20))]))
    w(root, "77b", "gene_origin_summary_kir.tsv", pd.DataFrame([dict(gene=g, n_windows=300, n_children=170, ambiguous_pct=20, ambiguous_n=60, nowin_pct=0, strong_of_resolved_pct=85, kmers_median=30,
                                                                    kmers_p10=2, kmers_p90=200) for g in ["ALL"] + genes]))
    w(root, "77b", "kmers_per_window_kir.tsv", pd.DataFrame([dict(gene=g, bin=b, n=n) for g in genes for b, n in (("0", "<20"), ("1-4", 40), ("5-9", 50), ("10+", 200))]))
    tg = []
    for v, n, k in (("all", 2433, 67), ("tx_named", 2268, 23), ("strong_origin", 2131, 61), ("complete_only", 1020, "<20")):
        tg.append(dict(variant=v, gene="ALL", n_records=n, n_informative=n, neither_n=k, neither_pct=(100 * k / n if k != "<20" else ""), lo=1, hi=3))
    w(root, "77c", "trio_classes_group.tsv", pd.DataFrame(tg))
    w(root, "77c", "pair_classes_group.tsv", pd.DataFrame([dict(kind="pair", variant=v, locus=l, n_records=500, n_informative=480, no_share_n=40, no_share_pct=8, lo=6, hi=10)
                                                           for v in ("all", "complete_only") for l in ("ALL", "2DL1", "3DL3")]))
    if with_prune:
        for sub, views in (("76b", ["A_classical8", "B_child", "pair_complete"]), ("78b", ["B_child", "H_hap"])):
            cur, auc, summ, cv, pw = [], [], [], [], []
            for v in views:
                for m, a in (("n_called", 0.85), ("rel_cds", 0.7), ("template_distance", 0.5)):
                    K = 40
                    n = 4000
                    for x in np.arange(0, 20.01, 0.25):
                        nr = n * x / 100
                        e = min(K, K * (x / 100) * (1 + (a - 0.5) * 6) if x > 0 else 0)
                        cur.append(dict(view=v, metric=m, prune_pct=x, n_units=n, n_errors=K, n_pruned=nr, errors_removed=e, captured_pct=100 * e / K,
                                        residual_error_pct=100 * (K - e) / (n - nr), base_error_pct=100 * K / n, enrichment=(e / K) / (nr / n) if nr else np.nan,
                                        nnp=nr / e if e else np.inf, marginal_nnp=(5 + x * 2) if x >= 2 else np.nan, rand_removed_lo=K * x / 100 * 0.5, rand_removed_hi=K * x / 100 * 1.6,
                                        p_tail=0.001, captured_lo=100 * e / K - 5, captured_hi=100 * e / K + 5, residual_lo=0, residual_hi=2))
                    auc.append(dict(view=v, metric=m, n_units=n, n_errors=K, auc=a, auc_lo=a - 0.05, auc_hi=a + 0.05))
                    summ.append(dict(view=v, metric=m, x_opt=5.0, p_tail_at_5=0.0004 if a > 0.6 else 0.4))
                    cv.append(dict(view=v, metric=m, reps=20, pruned_pct=5, captured_pct=20, enrichment=1 + (a - 0.5) * 6, enrichment_lo=0.8, enrichment_hi=3, share_reps_enrichment_gt1=0.95 if a > 0.6 else 0.4))
                pw.append(dict(view=v, n_units=4000, n_errors=40, auc_min_detectable=0.6))
            w(root, sub, "prune_curves.tsv", pd.DataFrame(cur))
            w(root, sub, "auc.tsv", pd.DataFrame(auc))
            w(root, sub, "prune_summary.tsv", pd.DataFrame(summ))
            w(root, sub, "cv_by_metric.tsv", pd.DataFrame(cv))
            w(root, sub, "power.tsv", pd.DataFrame(pw))


class T79(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def run_build(self, with_prune=True):
        tabs = os.path.join(self.d, "agg")
        make_tables(tabs, with_prune)
        out = os.path.join(self.d, "out")
        self.assertEqual(M.main(["--tables", tabs, "--out", out]), 0)
        return open(os.path.join(out, "S06_FINAL_REPORT.html")).read(), out

    def test_full_build(self):
        h, out = self.run_build(True)
        self.assertIn("id='think-tank'", h)
        self.assertIn("data:image/png;base64", h)
        self.assertNotIn("nan%", h)
        for nm in ("fig_missing_overview", "fig_stream_b", "fig_prune_hla", "fig_prune_kir"):
            self.assertTrue(os.path.exists(os.path.join(out, "figs", nm + ".png")), nm)
        self.assertIn("evidence", h)
        self.assertIn("strong", h)           # the synthetic n_called curve is strongly enriched
        self.assertNotIn("figure not produced", h.split("5. Raw evidence")[0])

    def test_pending_without_prune_tables(self):
        h, _ = self.run_build(False)
        self.assertIn("Pending:", h)
        self.assertIn("76b/prune_curves.tsv", h)
        self.assertIn("id='think-tank'", h)

    def test_numbers_are_computed(self):
        tabs = os.path.join(self.d, "agg")
        make_tables(tabs, False)
        n = M.N(M.Tables(tabs))
        self.assertAlmostEqual(n.p1[2], 100 * 30 / 74)
        self.assertIn("30/74", n.sentence_parent())
        self.assertIn("7 carry an allele found in neither parent", n.sentence_b())

    def test_masked_cells_stay_nan(self):
        tabs = os.path.join(self.d, "agg")
        make_tables(tabs, False)
        t = M.Tables(tabs).rd("77b", "kmers_per_window_kir.tsv")
        self.assertTrue(t["n"].isna().any())

    def test_xopt_rule(self):
        d = pd.DataFrame(dict(prune_pct=[0, 1, 2, 3, 4], marginal_nnp=[np.nan, np.nan, 5, 12, 40]))
        self.assertEqual(M.xopt_at_cost(d, 25), 3.0)
        self.assertEqual(M.xopt_at_cost(d, 100), 4.0)

    def test_evidence_labels(self):
        self.assertTrue(M.evidence_label(0.0001, 0.95, 40, 0.8, 0.6).startswith("strong"))
        self.assertTrue(M.evidence_label(0.4, 0.4, 40, 0.5, 0.6).startswith("weak"))
        self.assertTrue(M.evidence_label(0.01, 0.8, 3, 0.7, 0.6).startswith("none"))


if __name__ == "__main__":
    unittest.main()
