#!/usr/bin/env python3
"""Unit tests for 44_kir_recurrence_saturation.py, synthetic data only.

Covers:
  1. Recurrence classing (classify_recurrence) -- exactness of eq1/eq2/gt2/ge20, including the
     documented gt2/ge20 overlap.
  2. Good-Turing incidence coverage (Chao & Jost 2012) -- Q2>0 and Q2==0 branches, and the
     "no unseen mass" edge case (coverage == 1).
  3. Chao2 SE -- non-negative, and the Q2==0 bias-corrected branch doesn't blow up / go negative.
  4. Equal-N rarefaction curves with FIXED, precomputed orders -- reproducibility (same orders =>
     identical curves) and correctness against a hand-countable tiny example.
  5. Equal-N slope sign/magnitude sanity on a synthetic saturating vs. still-climbing curve.
  6. Disclosure: build_recurrence_row/build_coverage_chao2_row mask every count/rate < 20,
     never emit an allele name next to a small count, and blank q1/q2 together (not one alone).
  7. A synthetic end-to-end run: fake gzip'd KIR hap GTFs (41/43's own fixture-writing convention)
     through parse_all -> kir_person_gene_sets -> run_gene_level -> the four TSVs, paired with a
     synthetic HLA `calls` DataFrame (the shape 39_saturation_by_ancestry.build_labeled_calls
     would hand back) through hla_person_gene_sets -> run_gene_level -- exercising the exact same
     code paths main() drives, without needing the VM-only refdata/cohort/relatedness inputs
     build_labeled_calls() itself requires. Then feeds the written TSVs through 45's figure
     functions and confirms they render without error.

Run: python3 scripts/hla_popgen/tests/test_44_kir_recurrence_saturation.py
"""
import gzip
import importlib.util
import os
import sys
import tempfile
import unittest
from collections import Counter

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


kir41 = _load_module("41_kir_pilot.py", "kir_pilot_test_target_44")
agg43 = _load_module("43_kir_full_aggregate.py", "kir_full_aggregate_test_target_44")
agg43._kir = kir41
m44 = _load_module("44_kir_recurrence_saturation.py", "recurrence_saturation_44_test_target")
# Point 44's lazy module loaders at the already-loaded 41/43 so tests don't re-import from disk.
m44._m41 = kir41
m44._m43 = agg43


def _gtf_transcript_line(gene, consensus):
    attrs = (f'gene_id "IG{gene}"; transcript_id "IAT{gene}.1"; gene_name "{gene}"; '
             f'consensus "{consensus}"; alleles "{consensus}";')
    return f"chr19_synth\tImmuannot\ttranscript\t100\t200\t.\t+\t.\t{attrs}\n"


def _write_gtf_gz(path, gene_consensus_pairs):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as f:
        f.write("##synthetic test gtf\n")
        for gene, consensus in gene_consensus_pairs:
            f.write(_gtf_transcript_line(gene, consensus))


class TestRecurrenceClassing(unittest.TestCase):
    def test_classify_recurrence_exact_and_overlap(self):
        # alleles: a (1 carrier), b (2 carriers), c (5 carriers), d (25 carriers)
        counts = Counter({"a": 1, "b": 2, "c": 5, "d": 25})
        cls = m44.classify_recurrence(counts)
        self.assertEqual(cls["eq1"], 1)
        self.assertEqual(cls["eq2"], 1)
        self.assertEqual(cls["gt2"], 2)   # c and d (>=3)
        self.assertEqual(cls["ge20"], 1)  # d only -- subset of gt2, by design

    def test_classify_recurrence_empty(self):
        self.assertEqual(m44.classify_recurrence(Counter()),
                         {"eq1": 0, "eq2": 0, "gt2": 0, "ge20": 0})

    def test_carrier_counts_counts_distinct_people_not_occurrences(self):
        unit_sets = {"p1": {"x", "y"}, "p2": {"x"}, "p3": {"x", "y"}}
        c = m44.carrier_counts(unit_sets)
        self.assertEqual(c["x"], 3)
        self.assertEqual(c["y"], 2)


class TestGoodTuringCoverage(unittest.TestCase):
    def test_no_unseen_mass_gives_full_coverage(self):
        # No singletons at all -> nothing suggests undetected species -> coverage == 1.
        Q = Counter({2: 5, 3: 2})
        cov = m44.good_turing_coverage(Q, m=50, n_incidences=200)
        self.assertAlmostEqual(cov, 1.0)

    def test_many_singletons_low_coverage(self):
        # Mostly singletons (rare alleles) -> low coverage.
        Q = Counter({1: 40, 2: 2})
        n_incidences = 40 * 1 + 2 * 2
        cov = m44.good_turing_coverage(Q, m=50, n_incidences=n_incidences)
        self.assertTrue(0.0 <= cov < 0.6)

    def test_q2_zero_branch_does_not_crash(self):
        Q = Counter({1: 3})
        cov = m44.good_turing_coverage(Q, m=10, n_incidences=3)
        self.assertTrue(0.0 <= cov <= 1.0)

    def test_zero_incidences_is_nan(self):
        self.assertTrue(np.isnan(m44.good_turing_coverage(Counter(), m=10, n_incidences=0)))


class TestChao2SE(unittest.TestCase):
    def test_se_nonnegative_q2_positive(self):
        se = m44.chao2_se(S_obs=100, Q=Counter({1: 10, 2: 4}), m=50)
        self.assertGreaterEqual(se, 0.0)

    def test_se_nonnegative_q2_zero(self):
        se = m44.chao2_se(S_obs=80, Q=Counter({1: 6}), m=30)
        self.assertGreaterEqual(se, 0.0)

    def test_se_zero_when_no_rarity(self):
        se = m44.chao2_se(S_obs=50, Q=Counter(), m=30)
        self.assertEqual(se, 0.0)

    def test_se_single_unit_is_zero(self):
        self.assertEqual(m44.chao2_se(S_obs=5, Q=Counter({1: 5}), m=1), 0.0)


class TestEqualNCurves(unittest.TestCase):
    def _small_unit_sets_list(self):
        # 4 "people", hand-countable alleles.
        return [{"a", "b"}, {"a"}, {"c"}, {"a", "c"}]

    def test_orders_are_reproducible_with_same_seed(self):
        o1 = m44.make_orders(4, 5, seed=7)
        o2 = m44.make_orders(4, 5, seed=7)
        for a, b in zip(o1, o2):
            np.testing.assert_array_equal(a, b)

    def test_rarefaction_matches_hand_count_full_order(self):
        unit_sets_list = self._small_unit_sets_list()
        order = np.array([0, 1, 2, 3])
        curves = m44.rarefaction_curves_fixed_orders(unit_sets_list, [order], thresholds=[1, 2, 3])
        # After all 4 people (order 0,1,2,3): a seen in units 0,1,3 (3x); b in unit0 (1x); c in
        # units 2,3 (2x). Distinct alleles overall = 3 (a,b,c).
        self.assertEqual(curves["distinct"][0, -1], 3)
        # threshold-1 ("reached >=1 carrier") at the end = 3 (all three alleles).
        self.assertEqual(curves[1][0, -1], 3)
        # threshold-2 ("reached >=2 carriers") at the end: a (3) and c (2) qualify -> 2.
        self.assertEqual(curves[2][0, -1], 2)
        # threshold-3: only a (3 carriers) -> 1.
        self.assertEqual(curves[3][0, -1], 1)

    def test_distinct_curve_is_monotone_nondecreasing(self):
        unit_sets_list = self._small_unit_sets_list()
        orders = m44.make_orders(4, 10, seed=1)
        curves = m44.rarefaction_curves_fixed_orders(unit_sets_list, orders, thresholds=[1, 2])
        diffs = np.diff(curves["distinct"], axis=1)
        self.assertTrue((diffs >= 0).all())

    def test_class_curves_from_thresholds_sums_consistently(self):
        unit_sets_list = self._small_unit_sets_list()
        order = np.array([0, 1, 2, 3])
        curves = m44.rarefaction_curves_fixed_orders(unit_sets_list, [order], m44.KIR_THRESHOLDS)
        classes = m44.class_curves_from_thresholds(curves)
        # eq1 + eq2 + gt2 must equal the >=1 "reached" curve at every step (partition of carriers).
        total = classes["eq1"] + classes["eq2"] + classes["gt2"]
        np.testing.assert_array_equal(total, curves[1])
        # ge20 must never exceed gt2 (subset).
        self.assertTrue((classes["ge20"] <= classes["gt2"]).all())

    def test_equal_n_slope_positive_for_still_rising_curve(self):
        steps = np.arange(1, 101)
        mean = steps * 2.0  # perfectly linear, still rising
        slope = m44.equal_n_slope(steps, mean, n_star=100)
        self.assertAlmostEqual(slope, 2000.0, delta=1.0)

    def test_equal_n_slope_near_zero_for_flat_curve(self):
        steps = np.arange(1, 101)
        mean = np.full(100, 50.0)  # fully saturated, flat
        slope = m44.equal_n_slope(steps, mean, n_star=100)
        self.assertAlmostEqual(slope, 0.0, delta=1e-6)

    def test_equal_n_slope_nan_when_n_star_out_of_range(self):
        steps = np.arange(1, 11)
        mean = steps.astype(float)
        self.assertTrue(np.isnan(m44.equal_n_slope(steps, mean, n_star=999)))

    def test_pick_n_star_excludes_underpowered_and_mid(self):
        n_by_anc = {"AFR": 500, "AMR": 400, "MID": 50, "SAS": 90}
        # SAS below MIN_PEOPLE_PER_ANCESTRY (100) -> excluded from the well-powered set too.
        n_star = m44.pick_n_star(n_by_anc, min_people=100, exclude=("MID",))
        self.assertEqual(n_star, 400)  # min(AFR=500, AMR=400)

    def test_pick_n_star_none_when_nobody_qualifies(self):
        self.assertIsNone(m44.pick_n_star({"MID": 10}, min_people=100))


class TestDisclosure(unittest.TestCase):
    def test_recurrence_row_masks_small_counts(self):
        unit_sets = {f"p{i}": {"only_allele"} for i in range(5)}  # 5 carriers, one allele
        row = m44.build_recurrence_row("KIR3DL1", "MID", "any_novel", "kir", unit_sets)
        self.assertEqual(row["n_distinct_alleles"], "<20")
        self.assertEqual(row["eq1"], "0")   # genuine zero stays "0"
        self.assertEqual(row["gt2"], "<20")  # the one allele has 5 carriers -> gt2 bucket, masked

    def test_recurrence_row_true_zero_not_masked(self):
        row = m44.build_recurrence_row("KIR2DS3", "EAS", "protein_novel", "kir", {})
        self.assertEqual(row["n_distinct_alleles"], "0")
        self.assertEqual(row["eq1"], "0")

    def test_coverage_row_blanks_q1_q2_together_when_either_is_small(self):
        unit_sets = {f"p{i}": {f"a{i}"} for i in range(5)}  # 5 singletons, q1=5<20
        row = m44.build_coverage_chao2_row("A", "SAS", "all", "hla", unit_sets)
        self.assertEqual(row["q1"], "")
        self.assertEqual(row["q2"], "")

    def test_coverage_row_keeps_large_q1_q2(self):
        # 30 people, each with a unique allele (q1=30) plus one allele shared by exactly 2 people
        # in every pair grouping -> construct q2 >= 20 too.
        unit_sets = {}
        for i in range(30):
            unit_sets[f"s{i}"] = {f"singleton{i}"}
        for i in range(20):
            unit_sets[f"d{i}a"] = unit_sets.get(f"d{i}a", set()) | {f"dup{i}"}
            unit_sets[f"d{i}b"] = unit_sets.get(f"d{i}b", set()) | {f"dup{i}"}
        row = m44.build_coverage_chao2_row("A", "AFR", "all", "hla", unit_sets)
        self.assertNotEqual(row["q1"], "")
        self.assertNotEqual(row["q2"], "")


class TestKirLevelMask(unittest.TestCase):
    def test_levels(self):
        self.assertTrue(m44.kir_level_mask("known", "all"))
        self.assertFalse(m44.kir_level_mask("undetermined", "all"))
        self.assertFalse(m44.kir_level_mask("known", "any_novel"))
        self.assertTrue(m44.kir_level_mask("novel_genomic_known_cds", "any_novel"))
        self.assertTrue(m44.kir_level_mask("novel_protein", "protein_novel"))
        self.assertFalse(m44.kir_level_mask("novel_cds_synonymous", "protein_novel"))


class TestHlaLevelMask(unittest.TestCase):
    def _frame(self):
        return pd.DataFrame({
            "keep_clean": [True, True, True, False],
            "field_class": ["known", "f4_noncoding", "f2_protein", "f2_protein"],
            "seq_class": ["cds_known", "cds_known", "novel_protein", "novel_protein"],
        })

    def test_all_level_is_keep_clean_only(self):
        df = self._frame()
        mask = m44.hla_level_mask(df, "all")
        self.assertEqual(list(mask), [True, True, True, False])

    def test_any_novel_excludes_known(self):
        df = self._frame()
        mask = m44.hla_level_mask(df, "any_novel")
        self.assertEqual(list(mask), [False, True, True, False])

    def test_protein_novel_is_strict(self):
        df = self._frame()
        mask = m44.hla_level_mask(df, "protein_novel")
        self.assertEqual(list(mask), [False, False, True, False])  # row 3 dropped: not keep_clean

    def test_id_col_choice(self):
        self.assertEqual(m44.hla_id_col("protein_novel"), "prot_id")
        self.assertEqual(m44.hla_id_col("all"), "cds_id")
        self.assertEqual(m44.hla_id_col("any_novel"), "cds_id")


class TestEndToEndSynthetic(unittest.TestCase):
    """Fake KIR GTFs -> parse_all -> kir_person_gene_sets -> run_gene_level, paired with a fake
    HLA `calls` frame (shape build_labeled_calls would produce) -> hla_person_gene_sets ->
    run_gene_level, feeding the same four TSVs 44's main() writes, then rendering 45's figures on
    them. Deliberately bypasses build_labeled_calls()/cohort_membership/relatedness-table I/O
    (VM-only inputs) -- see module docstring."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.outroot = os.path.join(self.tmp.name, "pipeline_outputs_kir")
        # 6 synthetic KIR people: 2 per ancestry (AFR, EUR), enough for a tiny curve.
        self.people = {
            "111111": ("AFR", [("KIR3DL1", "KIR3DL1*001"), ("KIR2DL1", "KIR2DL1*003new")]),
            "222222": ("AFR", [("KIR3DL1", "KIR3DL1*001"), ("KIR2DL1", "KIR2DL1*003new")]),
            "333333": ("AFR", [("KIR3DL1", "KIR3DL1*002new")]),
            "444444": ("EUR", [("KIR3DL1", "KIR3DL1*001")]),
            "555555": ("EUR", [("KIR3DL1", "KIR3DL1*001")]),
            "666666": ("EUR", [("KIR2DL1", "KIR2DL1*003new")]),
        }
        for pid, (anc, calls) in self.people.items():
            hap_dir = os.path.join(self.outroot, pid, "immuannot_output")
            _write_gtf_gz(os.path.join(hap_dir, "hap1.gtf.gz"), calls)
            _write_gtf_gz(os.path.join(hap_dir, "hap2.gtf.gz"), [])
        self.ancestry_of = {pid: anc for pid, (anc, _) in self.people.items()}

    def tearDown(self):
        self.tmp.cleanup()

    def _fake_hla_calls(self):
        # Mirrors build_labeled_calls()'s output columns for the "A" gene, same 6 people, with a
        # couple of any_novel / protein_novel calls so the HLA side has non-trivial curves too.
        rows = []
        data = {
            "111111": ("A*01:01", "known", "cds_known", "A*01:01", "A_cds_x1"),
            "222222": ("A*01:01", "known", "cds_known", "A*01:01", "A_cds_x1"),
            "333333": ("A_prot_ab12cd34", "f2_protein", "novel_protein", "A_prot_ab12cd34", "A_cds_novel1"),
            "444444": ("A*02:01", "known", "cds_known", "A*02:01", "A_cds_x2"),
            "555555": ("A*02:01", "known", "cds_known", "A*02:01", "A_cds_x2"),
            "666666": ("A_noncoding_new", "f4_noncoding", "cds_known", "A*02:01", "A_cds_x2b"),
        }
        for pid, (_cons, fc, sc, prot_id, cds_id) in data.items():
            rows.append({"person_id": pid, "gene_b": "A", "field_class": fc, "seq_class": sc,
                        "prot_id": prot_id, "cds_id": cds_id, "keep_clean": True})
        return pd.DataFrame(rows)

    def test_end_to_end_writes_four_tsvs_and_45_renders(self):
        kir = kir41
        agg = agg43
        pids = agg.discover_people(self.outroot)
        self.assertEqual(len(pids), 6)
        results = agg.parse_all(pids, self.outroot, workers=1,
                                kir_script_path=os.path.join(HLA_POPGEN_DIR, "41_kir_pilot.py"))
        results_by_pid = {r["person_id"]: r for r in results}

        calls = self._fake_hla_calls()

        pids_by_ancestry = {
            "ALL": pids,
            "AFR": [p for p in pids if self.ancestry_of[p] == "AFR"],
            "EUR": [p for p in pids if self.ancestry_of[p] == "EUR"],
        }
        n_permutations = 5
        orders_by_ancestry = {a: m44.make_orders(len(ps), n_permutations, seed=44 + hash(a) % 1000)
                              for a, ps in pids_by_ancestry.items()}
        n_star_by_ancestry = {"ALL": 6, "AFR": 3, "EUR": 3}

        all_rec, all_curve, all_slope, all_cov = [], [], [], []
        for level in m44.LEVELS:
            kir_sets_all = m44.kir_person_gene_sets(results_by_pid, pids, level, kir)
            hla_sets_all = m44.hla_person_gene_sets(calls, level, ["A"])
            for gene in ["KIR3DL1", "KIR2DL1"]:
                species_sets = {"kir": kir_sets_all.get(gene, {})}
                rec, curve, slope, cov = m44.run_gene_level(
                    gene, level, species_sets, pids_by_ancestry, orders_by_ancestry,
                    n_star_by_ancestry, stride=1, extrapolate_2n=True)
                all_rec += rec; all_curve += curve; all_slope += slope; all_cov += cov
            for gene in ["A"]:
                species_sets = {"hla": hla_sets_all.get(gene, {})}
                rec, curve, slope, cov = m44.run_gene_level(
                    gene, level, species_sets, pids_by_ancestry, orders_by_ancestry,
                    n_star_by_ancestry, stride=1, extrapolate_2n=True)
                all_rec += rec; all_curve += curve; all_slope += slope; all_cov += cov

        out_dir = os.path.join(self.tmp.name, "results_44")
        os.makedirs(out_dir, exist_ok=True)
        pd.DataFrame(all_rec).to_csv(os.path.join(out_dir, "recurrence_classes.tsv"),
                                     sep="\t", index=False)
        pd.DataFrame(all_curve).to_csv(os.path.join(out_dir, "saturation_curves.tsv"),
                                       sep="\t", index=False)
        pd.DataFrame(all_slope).to_csv(os.path.join(out_dir, "equal_n_slope.tsv"),
                                       sep="\t", index=False)
        pd.DataFrame(all_cov).to_csv(os.path.join(out_dir, "coverage_chao2.tsv"),
                                     sep="\t", index=False)
        for fname in ("recurrence_classes.tsv", "saturation_curves.tsv", "equal_n_slope.tsv",
                     "coverage_chao2.tsv"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, fname)))

        # Disclosure sanity: no small count leaked as a bare integer < 20 in recurrence table.
        rec_df = pd.read_csv(os.path.join(out_dir, "recurrence_classes.tsv"), sep="\t", dtype=str)
        for col in ("n_distinct_alleles", "eq1", "eq2", "gt2", "ge20"):
            for v in rec_df[col]:
                if v not in ("0", "<20"):
                    self.assertGreaterEqual(int(v), 20, f"{col}={v} should have been masked")

        # 45's figure functions must render on this output without raising.
        m45 = _load_module("45_kir_recurrence_figure.py", "kir_recurrence_figure_45_test_target")
        rec, curve, cov = m45.load_tables(out_dir)
        m45.fig_saturation_paired(curve, os.path.join(out_dir, "fig_saturation_paired"))
        m45.fig_recurrence_classes(rec, os.path.join(out_dir, "fig_recurrence_classes"))
        m45.fig_coverage_chao2(cov, os.path.join(out_dir, "fig_coverage_chao2"))
        for stem in ("fig_saturation_paired", "fig_recurrence_classes", "fig_coverage_chao2"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, stem + ".png")))
            self.assertTrue(os.path.exists(os.path.join(out_dir, stem + ".pdf")))


class TestJointUnrelatedSet(unittest.TestCase):
    def test_intersects_before_removing_relatives(self):
        kir_pids = ["1", "2", "3", "4"]
        hla_pids = ["2", "3", "4", "5"]
        pairs = [("2", "3", 0.10)]  # related pair, both in the intersection {2,3,4}
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
            f.write("i.s\tj.s\tkin\n2\t3\t0.10\n")
            path = f.name
        try:
            kept, n_removed, n_common = m44.joint_unrelated_set(kir_pids, hla_pids, path, 0.0442)
            self.assertEqual(n_common, 3)   # {2,3,4}
            self.assertEqual(n_removed, 1)  # one of 2/3 dropped
            self.assertNotIn("1", kept)     # not in HLA pool
            self.assertNotIn("5", kept)     # not in KIR pool
            self.assertIn("4", kept)
        finally:
            os.remove(path)


if __name__ == "__main__":
    unittest.main()
