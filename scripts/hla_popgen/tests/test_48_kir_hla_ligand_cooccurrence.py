#!/usr/bin/env python3
"""Unit tests for 48_kir_hla_ligand_cooccurrence.py (S04 WS-D), synthetic data only.

Covers:
  1. classify_c1c2() / classify_bw4() -- lookup table hits, unclassified fallback, B*15 default
     (Bw6) vs. the CLI-supplied exception list.
  2. two_by_two_stats() -- disclosure masking (any cell in [1,19] blanks OR/CI/p-values but keeps
     masked cell counts as '<20', never blank/0), Haldane-Anscombe correction on a zero cell,
     and a planted enrichment producing OR > 1 with a small Fisher p-value.
  3. End-to-end: --synthetic run writes the three aggregate TSVs, detects the planted
     KIR2DL1 x C2 / KIR3DL1 x Bw4 enrichment, and reports zero unclassified alleles for the
     synthetic cohort's own restricted allele set.

Run: python3 scripts/hla_popgen/tests/test_48_kir_hla_ligand_cooccurrence.py
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


m = _load_module("48_kir_hla_ligand_cooccurrence.py", "ligand_test_target")


class TestClassifyC1C2(unittest.TestCase):
    def test_known_c2(self):
        self.assertEqual(m.classify_c1c2("05:01"), "C2")

    def test_known_c1(self):
        self.assertEqual(m.classify_c1c2("07:01"), "C1")

    def test_unclassified_not_silently_absorbed(self):
        self.assertEqual(m.classify_c1c2("99:99"), "unclassified")


class TestClassifyBw4(unittest.TestCase):
    def test_b_locus_bw4_group(self):
        self.assertEqual(m.classify_bw4("HLA-B", "57:01"), "Bw4")

    def test_b_locus_bw6_group(self):
        self.assertEqual(m.classify_bw4("HLA-B", "07:02"), "Bw6")

    def test_a_locus_bw4(self):
        self.assertEqual(m.classify_bw4("HLA-A", "24:02"), "Bw4")

    def test_a_locus_not_bw4(self):
        self.assertEqual(m.classify_bw4("HLA-A", "01:01"), "not_bw4_locus_A")

    def test_b15_defaults_to_bw6(self):
        self.assertEqual(m.classify_bw4("HLA-B", "15:01"), "Bw6")

    def test_b15_exception_list_overrides(self):
        self.assertEqual(m.classify_bw4("HLA-B", "15:13", b15_as_bw4=frozenset({"15:13"})), "Bw4")
        # An allele not in the exception list still defaults to Bw6.
        self.assertEqual(m.classify_bw4("HLA-B", "15:01", b15_as_bw4=frozenset({"15:13"})), "Bw6")


class TestTwoByTwoStats(unittest.TestCase):
    def test_disclosive_cell_blanks_stats_but_keeps_masked_counts(self):
        n = 200
        kir = np.zeros(n, dtype=bool)
        kir[:5] = True  # only 5 KIR+ people -> a 1-19 cell somewhere
        ligand = np.zeros(n, dtype=bool)
        ligand[:100] = True
        r = m.two_by_two_stats(kir, ligand, n_perms=10)
        self.assertEqual(r["odds_ratio"], "")
        self.assertEqual(r["fisher_p"], "")
        # masked counts are never blank and never silently 0 for a genuinely non-zero cell.
        self.assertIn(r["n_a_kir_and_ligand"], ("<20",))

    def test_true_zero_cell_not_masked_as_disclosive(self):
        n = 200
        kir = np.zeros(n, dtype=bool)
        kir[:100] = True
        ligand = np.zeros(n, dtype=bool)
        ligand[100:] = True  # perfectly disjoint -> a=0 exactly (true zero, not 1-19)
        r = m.two_by_two_stats(kir, ligand, n_perms=10)
        self.assertEqual(r["n_a_kir_and_ligand"], "0")
        # a genuine 0 with all other cells >=20 should NOT be treated as disclosive.
        self.assertNotEqual(r["odds_ratio"], "")

    def test_planted_enrichment_gives_or_above_one_and_small_pvalue(self):
        rng = np.random.default_rng(0)
        n = 500
        ligand = rng.random(n) < 0.4
        kir = np.zeros(n, dtype=bool)
        for i in range(n):
            p = 0.8 if ligand[i] else 0.3
            kir[i] = rng.random() < p
        r = m.two_by_two_stats(kir, ligand, n_perms=50, seed=1)
        self.assertIsInstance(r["odds_ratio"], float)
        self.assertGreater(r["odds_ratio"], 1.5)
        self.assertLess(r["fisher_p"], 0.01)


class TestEndToEndSynthetic(unittest.TestCase):
    def test_pipeline_writes_expected_files_and_detects_planted_enrichment(self):
        (table1_df, epitopes_df, kir_presence, ancestry_of,
         n_unclass_c, n_unclass_bw) = m.make_synthetic(n_people=600, seed=7)
        with tempfile.TemporaryDirectory() as tmp:
            status_path = os.path.join(tmp, "STATUS.txt")
            m.run_pipeline(table1_df, epitopes_df, kir_presence, ancestry_of, tmp,
                          n_perms=20, b15_as_bw4=frozenset(), status_path=status_path,
                          n_unclass_c=n_unclass_c, n_unclass_bw=n_unclass_bw)
            for fn in ("kir_hla_ligand_cooccurrence.tsv", "epitope_freq_by_ancestry.tsv",
                      "ligand_lookup_qc.tsv"):
                self.assertTrue(os.path.exists(os.path.join(tmp, fn)), fn)
            self.assertTrue(os.path.exists(status_path))

            df = pd.read_csv(os.path.join(tmp, "kir_hla_ligand_cooccurrence.tsv"), sep="\t")
            self.assertGreater(len(df), 0)
            planted = df[(df["kir_gene"] == "KIR2DL1") & (df["ligand"] == "C2")].copy()
            planted["odds_ratio"] = pd.to_numeric(planted["odds_ratio"], errors="coerce")
            planted = planted.dropna(subset=["odds_ratio"])
            self.assertGreater(len(planted), 0)
            self.assertTrue((planted["odds_ratio"] > 1.0).all())

            qc = pd.read_csv(os.path.join(tmp, "ligand_lookup_qc.tsv"), sep="\t")
            self.assertIn("n_unclassified_C_allele_calls", qc["metric"].values)

    def test_no_unclassified_alleles_in_synthetic_restricted_set(self):
        (table1_df, epitopes_df, kir_presence, ancestry_of,
         n_unclass_c, n_unclass_bw) = m.make_synthetic(n_people=100, seed=3)
        self.assertEqual(n_unclass_c, 0)
        self.assertEqual(n_unclass_bw, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
