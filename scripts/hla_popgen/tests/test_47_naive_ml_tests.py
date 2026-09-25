#!/usr/bin/env python3
"""Unit tests for 47_naive_ml_tests.py (S04 WS-D), synthetic data only.

Covers:
  1. hla_two_field() -- colon-field truncation, 'new'-token handling, null tokens.
  2. build_hla_carriage() -- rare-allele (<20 carriers) drop, matrix shape.
  3. run_binary_task() -- underpowered guard (n_pos or n_neg < 20) returns None; a planted
     signal produces AUROC well above 0.5 and above its own permutation null; top features are
     drawn only from the already-filtered (>=20-carrier) feature set.
  4. suppressed() -- the 1-19 -> '<20' / true-0 -> '0' rule.
  5. End-to-end: --synthetic run via make_synthetic()+run_pipeline() writes the three aggregate
     TSVs and a STATUS.txt, with the planted ancestry signal detectable and the (unplanted)
     cB-from-HLA signal near chance.

Run: python3 scripts/hla_popgen/tests/test_47_naive_ml_tests.py
"""
import importlib.util
import os
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)


def _load_module(filename, modname):
    spec = importlib.util.spec_from_file_location(modname, os.path.join(HLA_POPGEN_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


m = _load_module("47_naive_ml_tests.py", "naive_ml_test_target")


class TestHlaTwoField(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(m.hla_two_field("HLA-A*02:01:01:01"), "02:01")

    def test_new_token_truncates(self):
        self.assertEqual(m.hla_two_field("HLA-A*02:new"), None)  # <2 resolved fields
        self.assertEqual(m.hla_two_field("HLA-A*02:01:new"), "02:01")

    def test_null_tokens(self):
        self.assertIsNone(m.hla_two_field("undetermined"))
        self.assertIsNone(m.hla_two_field(None))
        self.assertIsNone(m.hla_two_field(float("nan")))

    def test_single_field_is_none(self):
        self.assertIsNone(m.hla_two_field("HLA-A*02"))


class TestBuildHlaCarriage(unittest.TestCase):
    def test_rare_allele_dropped(self):
        rows = []
        # 25 people carry HLA-A*01:01 (common), 3 people carry HLA-A*99:99 (rare).
        for i in range(25):
            rows.append({"person_id": f"p{i}", "hap": "hap1", "gene": "HLA-A",
                        "consensus": "HLA-A*01:01"})
        for i in range(25, 28):
            rows.append({"person_id": f"p{i}", "hap": "hap1", "gene": "HLA-A",
                        "consensus": "HLA-A*99:99"})
        df = pd.DataFrame(rows)
        mat, n_dropped = m.build_hla_carriage(df, min_carriers=20)
        self.assertIn("HLA-A*01:01", mat.columns)
        self.assertNotIn("HLA-A*99:99", mat.columns)
        self.assertEqual(n_dropped, 1)
        self.assertTrue(mat.loc["p0", "HLA-A*01:01"])
        self.assertFalse(mat.loc["p25", "HLA-A*01:01"])


class TestSuppressed(unittest.TestCase):
    def test_true_zero_stays_zero(self):
        self.assertEqual(m.suppressed(0), "0")

    def test_small_count_masked(self):
        for n in (1, 5, 19):
            self.assertEqual(m.suppressed(n), "<20")

    def test_large_count_passthrough(self):
        self.assertEqual(m.suppressed(20), "20")
        self.assertEqual(m.suppressed(137), "137")


class TestRunBinaryTask(unittest.TestCase):
    def test_underpowered_returns_none(self):
        rng = np.random.default_rng(0)
        X = rng.random((50, 5))
        y = np.zeros(50, dtype=int)
        y[:5] = 1  # only 5 positives, < MIN_CARRIERS_FOR_FEATURE
        self.assertIsNone(m.run_binary_task(X, y, [f"f{i}" for i in range(5)], n_perms=5))

    def test_planted_signal_beats_chance_and_permutation_null(self):
        rng = np.random.default_rng(1)
        n = 300
        y = (rng.random(n) < 0.4).astype(int)
        # Feature 0 is strongly informative; features 1-4 are pure noise.
        X = rng.random((n, 5))
        X[:, 0] = y * 3 + rng.normal(0, 0.5, size=n)
        feat_names = ["informative", "noise1", "noise2", "noise3", "noise4"]
        r = m.run_binary_task(X, y, feat_names, n_perms=20, seed=0)
        self.assertIsNotNone(r)
        self.assertGreater(r["lr_auroc"], 0.85)
        self.assertGreater(r["lr_auroc"], r["perm_auroc_p95"])
        self.assertEqual(r["top_lr_features"][0][0], "informative")

    def test_no_signal_hovers_near_chance(self):
        rng = np.random.default_rng(2)
        n = 300
        y = (rng.random(n) < 0.4).astype(int)
        X = rng.random((n, 5))  # pure noise, unrelated to y
        r = m.run_binary_task(X, y, [f"f{i}" for i in range(5)], n_perms=20, seed=0)
        self.assertIsNotNone(r)
        self.assertLess(abs(r["lr_auroc"] - 0.5), 0.15)


class TestEndToEndSynthetic(unittest.TestCase):
    def test_pipeline_writes_expected_files_and_detects_planted_ancestry_signal(self):
        table1_df, cohort_df, kir_content_true = m.make_synthetic(n_people=300, seed=42)
        pids = sorted(cohort_df["person_id"])
        kir_presence, kir_content = m.build_kir_features_synthetic(pids, kir_content_true, seed=42)
        with tempfile.TemporaryDirectory() as tmp:
            status_path = os.path.join(tmp, "STATUS.txt")
            m.run_pipeline(table1_df, cohort_df, kir_presence, kir_content, tmp,
                          "platform", n_perms=10, pca_bins=10, status_path=status_path)
            self.assertTrue(os.path.exists(os.path.join(tmp, "naive_ml_metrics.tsv")))
            self.assertTrue(os.path.exists(os.path.join(tmp, "naive_ml_top_features.tsv")))
            self.assertTrue(os.path.exists(status_path))
            metrics = pd.read_csv(os.path.join(tmp, "naive_ml_metrics.tsv"), sep="\t")
            self.assertGreater(len(metrics), 0)
            # Ancestry tasks should show real signal (planted association).
            anc_rows = metrics[metrics["task"].str.startswith("ancestry_")]
            self.assertGreater(len(anc_rows), 0)
            self.assertTrue((anc_rows["lr_auroc"] > 0.55).any())
            # cB_from_HLA has no planted association -- should not blow past chance.
            cb_rows = metrics[metrics["task"] == "cB_from_HLA"]
            if len(cb_rows):
                self.assertLess(abs(cb_rows.iloc[0]["lr_auroc"] - 0.5), 0.2)

    def test_missing_platform_column_exits_loudly(self):
        table1_df, cohort_df, kir_content_true = m.make_synthetic(n_people=60, seed=1)
        pids = sorted(cohort_df["person_id"])
        kir_presence, kir_content = m.build_kir_features_synthetic(pids, kir_content_true, seed=1)
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(SystemExit):
                m.run_pipeline(table1_df, cohort_df, kir_presence, kir_content, tmp,
                              "totally_not_a_real_column", n_perms=5, pca_bins=5,
                              status_path=os.path.join(tmp, "STATUS.txt"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
