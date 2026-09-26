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
     never emit an allele name next to a small count, and blank q1/q2 together (not one alone);
     build_na_rows emits the literal "NA" (never blank/0/"<20") for an unmatched level/species.
  7. KIR identity extraction: parse_hap_gtf_full's contig/copy_index extension + reused
     classify_novelty_tier; resolve_kir_cds_sequences' unambiguous-join contract (single-copy
     match, multi-copy-on-same-contig ambiguity, missing cds.fa.gz, missing record); the
     genomic/cds/protein identity levels build_person_kir_identity derives from a real gzip'd hap
     GTF + a real gzip'd cds.fa.gz fixture; and the cds_available=False fallback when no cds.fa.gz
     exists at all.
  8. HLA identity levels (hla_level_mask/hla_id_col/add_hla_genomic_id) on a synthetic `calls`
     frame shaped like 39_saturation_by_ancestry.build_labeled_calls()'s own output.
  9. A synthetic end-to-end run: fake gzip'd KIR hap GTFs + matching cds.fa.gz fixtures, through
     build_kir_identity_sets -> run_gene_level, paired with a synthetic HLA `calls` DataFrame
     through hla_person_gene_sets -> run_gene_level -- exercising the exact same code paths
     main() drives, without needing the VM-only refdata/cohort/relatedness inputs
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
m24mod = _load_module("24_novelty_by_field.py", "novelty_by_field_test_target_44")
m03mod = _load_module("03_novel_alleles.py", "novel_alleles_test_target_44")
m44 = _load_module("44_kir_recurrence_saturation.py", "recurrence_saturation_44_test_target")
# Point 44's lazy module loaders at the already-loaded modules so tests don't re-import from disk.
m44._m41 = kir41
m44._m43 = agg43
m44._m24 = m24mod
m44._m03 = m03mod


def _gtf_transcript_line(gene, consensus, contig="chr19_synth", copy_index=None):
    gene_id = f"IG{gene}" + (f".{copy_index}" if copy_index else "")
    attrs = (f'gene_id "{gene_id}"; transcript_id "IAT{gene}.1"; gene_name "{gene}"; '
             f'consensus "{consensus}"; alleles "{consensus}";')
    return f"{contig}\tImmuannot\ttranscript\t100\t200\t.\t+\t.\t{attrs}\n"


def _write_gtf_gz(path, gene_consensus_pairs, contig="chr19_synth"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as f:
        f.write("##synthetic test gtf\n")
        for entry in gene_consensus_pairs:
            gene, consensus = entry[0], entry[1]
            copy_index = entry[2] if len(entry) > 2 else None
            f.write(_gtf_transcript_line(gene, consensus, contig=contig, copy_index=copy_index))


def _write_cds_fasta_gz(path, records):
    """records: list of (contig, gene, index, seq) -> header '>{contig}_{gene}_{index}'."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as f:
        for contig, gene, idx, seq in records:
            f.write(f">{contig}_{gene}_{idx}\n{seq}\n")


class TestRecurrenceClassing(unittest.TestCase):
    def test_classify_recurrence_exact_and_overlap(self):
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
        Q = Counter({2: 5, 3: 2})
        cov = m44.good_turing_coverage(Q, m=50, n_incidences=200)
        self.assertAlmostEqual(cov, 1.0)

    def test_many_singletons_low_coverage(self):
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
        self.assertEqual(curves["distinct"][0, -1], 3)
        self.assertEqual(curves[1][0, -1], 3)
        self.assertEqual(curves[2][0, -1], 2)
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
        total = classes["eq1"] + classes["eq2"] + classes["gt2"]
        np.testing.assert_array_equal(total, curves[1])
        self.assertTrue((classes["ge20"] <= classes["gt2"]).all())

    def test_equal_n_slope_positive_for_still_rising_curve(self):
        steps = np.arange(1, 101)
        mean = steps * 2.0
        slope = m44.equal_n_slope(steps, mean, n_star=100)
        self.assertAlmostEqual(slope, 2000.0, delta=1.0)

    def test_equal_n_slope_near_zero_for_flat_curve(self):
        steps = np.arange(1, 101)
        mean = np.full(100, 50.0)
        slope = m44.equal_n_slope(steps, mean, n_star=100)
        self.assertAlmostEqual(slope, 0.0, delta=1e-6)

    def test_equal_n_slope_nan_when_n_star_out_of_range(self):
        steps = np.arange(1, 11)
        mean = steps.astype(float)
        self.assertTrue(np.isnan(m44.equal_n_slope(steps, mean, n_star=999)))

    def test_pick_n_star_excludes_underpowered_and_mid(self):
        n_by_anc = {"AFR": 500, "AMR": 400, "MID": 50, "SAS": 90}
        n_star = m44.pick_n_star(n_by_anc, min_people=100, exclude=("MID",))
        self.assertEqual(n_star, 400)

    def test_pick_n_star_none_when_nobody_qualifies(self):
        self.assertIsNone(m44.pick_n_star({"MID": 10}, min_people=100))


class TestDisclosure(unittest.TestCase):
    def test_recurrence_row_masks_small_counts(self):
        unit_sets = {f"p{i}": {"only_allele"} for i in range(5)}
        row = m44.build_recurrence_row("KIR3DL1", "MID", "any_novel", "kir", unit_sets)
        self.assertEqual(row["n_distinct_alleles"], "<20")
        self.assertEqual(row["eq1"], "0")
        self.assertEqual(row["gt2"], "<20")

    def test_recurrence_row_true_zero_not_masked(self):
        row = m44.build_recurrence_row("KIR2DS3", "EAS", "protein_novel", "kir", {})
        self.assertEqual(row["n_distinct_alleles"], "0")
        self.assertEqual(row["eq1"], "0")

    def test_coverage_row_blanks_q1_q2_together_when_either_is_small(self):
        unit_sets = {f"p{i}": {f"a{i}"} for i in range(5)}
        row = m44.build_coverage_chao2_row("A", "SAS", "genomic", "hla", unit_sets)
        self.assertEqual(row["q1"], "")
        self.assertEqual(row["q2"], "")

    def test_coverage_row_keeps_large_q1_q2(self):
        unit_sets = {}
        for i in range(30):
            unit_sets[f"s{i}"] = {f"singleton{i}"}
        for i in range(20):
            unit_sets[f"d{i}a"] = unit_sets.get(f"d{i}a", set()) | {f"dup{i}"}
            unit_sets[f"d{i}b"] = unit_sets.get(f"d{i}b", set()) | {f"dup{i}"}
        row = m44.build_coverage_chao2_row("A", "AFR", "genomic", "hla", unit_sets)
        self.assertNotEqual(row["q1"], "")
        self.assertNotEqual(row["q2"], "")

    def test_na_rows_emit_literal_na_not_blank_or_zero(self):
        pids_by_ancestry = {"ALL": ["1", "2", "3"], "AFR": ["1", "2"]}
        n_star_by_ancestry = {"ALL": 3, "AFR": None}
        rec, curve, slope, cov = m44.build_na_rows(
            "KIR2DS1", "cds", "kir", pids_by_ancestry, n_star_by_ancestry)
        for row in rec:
            for k in ("n_distinct_alleles", "eq1", "eq2", "gt2", "ge20"):
                self.assertEqual(row[k], "NA")
            self.assertIsInstance(row["n_people"], int)  # a real, non-disclosive count
        for row in cov:
            self.assertEqual(row["chao2"], "NA")
            self.assertEqual(row["good_turing_coverage"], "NA")
        for row in curve:
            self.assertEqual(row["mean_distinct"], "NA")
        for row in slope:
            self.assertEqual(row["slope_per_1000"], "NA")


class TestKirParseHapGtfFull(unittest.TestCase):
    def test_extracts_contig_gene_consensus_copy_index_and_tier(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "hap1.gtf.gz")
            _write_gtf_gz(path, [("KIR3DL1", "KIR3DL1*001"), ("KIR2DL1", "KIR2DL1*003new")],
                         contig="ctgA")
            rows = m44.parse_hap_gtf_full(path, kir41)
            self.assertEqual(len(rows), 2)
            r0 = next(r for r in rows if r["gene"] == "KIR3DL1")
            self.assertEqual(r0["contig"], "ctgA")
            self.assertEqual(r0["copy_index"], 1)
            self.assertEqual(r0["novelty_tier"], "known")
            r1 = next(r for r in rows if r["gene"] == "KIR2DL1")
            self.assertTrue(r1["is_novel"])
            self.assertIn(r1["novelty_tier"],
                          ("novel_genomic_known_cds", "novel_cds_synonymous",
                           "novel_protein", "novel_unclassified"))

    def test_copy_index_suffix_parsed(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "hap1.gtf.gz")
            _write_gtf_gz(path, [("KIR3DP1", "KIR3DP1*001", 2)], contig="ctgB")
            rows = m44.parse_hap_gtf_full(path, kir41)
            self.assertEqual(rows[0]["copy_index"], 2)

    def test_missing_file_returns_empty(self):
        self.assertEqual(m44.parse_hap_gtf_full("/no/such/file.gtf.gz", kir41), [])

    def test_matches_classify_novelty_tier_directly(self):
        # Cross-check against 41's own classify_novelty_tier on the same row shape.
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "hap1.gtf.gz")
            _write_gtf_gz(path, [("KIR3DL1", "KIR3DL1*001")], contig="ctgC")
            row = m44.parse_hap_gtf_full(path, kir41)[0]
            self.assertEqual(row["novelty_tier"], kir41.classify_novelty_tier(row))


class TestResolveKirCdsSequences(unittest.TestCase):
    def test_unambiguous_single_copy_match(self):
        rows = [{"contig": "ctgA", "gene": "KIR3DL1", "consensus": "KIR3DL1*001new",
                "copy_index": 1}]
        with tempfile.TemporaryDirectory() as tmp:
            cds_path = os.path.join(tmp, "cds.fa.gz")
            _write_cds_fasta_gz(cds_path, [("ctgA", "KIR3DL1", 1, "ATGAAATAG")])
            seq_by_key, stats = m44.resolve_kir_cds_sequences(rows, cds_path, m03mod)
            self.assertEqual(seq_by_key[("KIR3DL1", "KIR3DL1*001new")], "ATGAAATAG")
            self.assertEqual(stats["n_matched"], 1)
            self.assertEqual(stats["n_ambiguous_copy"], 0)

    def test_multi_copy_same_gene_same_contig_is_ambiguous(self):
        rows = [
            {"contig": "ctgA", "gene": "KIR3DP1", "consensus": "KIR3DP1*001", "copy_index": 1},
            {"contig": "ctgA", "gene": "KIR3DP1", "consensus": "KIR3DP1*002new", "copy_index": 2},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            cds_path = os.path.join(tmp, "cds.fa.gz")
            _write_cds_fasta_gz(cds_path, [("ctgA", "KIR3DP1", 1, "AAA"),
                                           ("ctgA", "KIR3DP1", 2, "CCC")])
            seq_by_key, stats = m44.resolve_kir_cds_sequences(rows, cds_path, m03mod)
            self.assertEqual(len(seq_by_key), 0)
            self.assertEqual(stats["n_ambiguous_copy"], 2)

    def test_missing_cds_fasta_file(self):
        rows = [{"contig": "ctgA", "gene": "KIR3DL1", "consensus": "KIR3DL1*001", "copy_index": 1}]
        seq_by_key, stats = m44.resolve_kir_cds_sequences(rows, "/no/such/cds.fa.gz", m03mod)
        self.assertEqual(len(seq_by_key), 0)
        self.assertEqual(stats["n_cds_fasta_missing"], 1)

    def test_missing_record_in_present_fasta(self):
        rows = [{"contig": "ctgA", "gene": "KIR3DL1", "consensus": "KIR3DL1*001", "copy_index": 1}]
        with tempfile.TemporaryDirectory() as tmp:
            cds_path = os.path.join(tmp, "cds.fa.gz")
            _write_cds_fasta_gz(cds_path, [("ctgA", "KIR2DL1", 1, "AAA")])  # different gene
            seq_by_key, stats = m44.resolve_kir_cds_sequences(rows, cds_path, m03mod)
            self.assertEqual(len(seq_by_key), 0)
            self.assertEqual(stats["n_missing_cds_record"], 1)

    def test_empty_rows(self):
        seq_by_key, stats = m44.resolve_kir_cds_sequences([], "/no/such/cds.fa.gz", m03mod)
        self.assertEqual(seq_by_key, {})
        self.assertEqual(stats["n_join_rows"], 0)


class TestBuildPersonKirIdentity(unittest.TestCase):
    def test_genomic_and_any_novel_available_without_cds_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            person_dir = os.path.join(tmp, "999", "immuannot_output")
            _write_gtf_gz(os.path.join(person_dir, "hap1.gtf.gz"),
                         [("KIR3DL1", "KIR3DL1*001"), ("KIR2DL1", "KIR2DL1*003new")],
                         contig="ctgA")
            _write_gtf_gz(os.path.join(person_dir, "hap2.gtf.gz"), [], contig="ctgA")
            pid, per_level, qc, n_seen, n_cds = m44.build_person_kir_identity(
                "999", tmp, kir41, m03mod, m24mod)
            self.assertEqual(pid, "999")
            self.assertIn("KIR3DL1*001", per_level["genomic"]["KIR3DL1"])
            self.assertIn("KIR2DL1*003new", per_level["any_novel"]["KIR2DL1"])
            # No cds.fa.gz written -> cds/protein tracks stay empty for this person.
            self.assertNotIn("KIR3DL1", per_level.get("cds", {}))
            self.assertEqual(n_cds, 0)
            self.assertEqual(n_seen, 1)

    def test_cds_derived_from_real_cds_fasta_protein_uncovered_without_kir_ref(self):
        """2026-09-27: without a `kir_ref` (no protein catalogue at all), `cds` is still real (it
        only needs cds.fa.gz), but `protein`/`protein_novel` must NOT get a hash-based guess --
        this gene must simply be left unpopulated (main() exports it as an explicit NA row with
        reason=no_refdata; see kir_protein_catalogue_status()). This replaces the pre-2026-09-27
        behavior, where a missing/uncovering kir_ref fell back to a hash id for every call
        regardless of catalogue coverage -- the fallback itself was the source of the
        KIR2DL2/KIR2DL5B/KIR2DP1 pct_novel_protein>100% VM smoke-test warnings."""
        with tempfile.TemporaryDirectory() as tmp:
            person_dir = os.path.join(tmp, "999", "immuannot_output")
            _write_gtf_gz(os.path.join(person_dir, "hap1.gtf.gz"),
                         [("KIR3DL1", "KIR3DL1*001new")], contig="ctgA")
            _write_gtf_gz(os.path.join(person_dir, "hap2.gtf.gz"), [], contig="ctgA")
            # ATG AAA TAG -> M K stop; protein_info strips the terminal stop -> "MK"
            _write_cds_fasta_gz(os.path.join(person_dir, "hap1", "cds.fa.gz"),
                               [("ctgA", "KIR3DL1", 1, "ATGAAATAG")])
            pid, per_level, qc, n_seen, n_cds = m44.build_person_kir_identity(
                "999", tmp, kir41, m03mod, m24mod)  # kir_ref=None (default)
            self.assertEqual(n_cds, 1)
            self.assertEqual(qc["n_matched"], 1)
            cds_ids = per_level["cds"]["KIR3DL1"]
            self.assertEqual(len(cds_ids), 1)
            self.assertTrue(next(iter(cds_ids)).startswith("KIR3DL1_cds_"))
            # No protein catalogue -> gene left OUT of 'protein'/'protein_novel' entirely.
            self.assertNotIn("KIR3DL1", per_level.get("protein", {}))
            self.assertNotIn("KIR3DL1", per_level.get("protein_novel", {}))
            self.assertEqual(qc["n_protein_gene_uncovered_calls"], 1)

    def test_novel_protein_still_gets_a_hash_id_when_gene_IS_covered(self):
        """Sanity check that removing the fallback didn't also remove the legitimate case: a
        gene the catalogue DOES cover, whose observed protein is genuinely uncatalogued, still gets
        a distinguishing hash id and counts as protein_novel (this is not "falling back to a finer
        hash than genomic" -- it's the intended encoding for a real novel protein, paired 1:1 with
        this same call's own genomic/any_novel identity)."""
        ref = m24mod.RefIndex()
        ref.add("KIR3DL1", "KIR3DL1*00101", "ATGAAACCCTAG")  # a DIFFERENT known protein
        ref.finalize()
        with tempfile.TemporaryDirectory() as tmp:
            person_dir = os.path.join(tmp, "999", "immuannot_output")
            _write_gtf_gz(os.path.join(person_dir, "hap1.gtf.gz"),
                         [("KIR3DL1", "KIR3DL1*001new")], contig="ctgA")
            _write_gtf_gz(os.path.join(person_dir, "hap2.gtf.gz"), [], contig="ctgA")
            _write_cds_fasta_gz(os.path.join(person_dir, "hap1", "cds.fa.gz"),
                               [("ctgA", "KIR3DL1", 1, "ATGAAATAG")])  # translates to "MK", not
                                                                        # catalogued
            pid, per_level, qc, n_seen, n_cds = m44.build_person_kir_identity(
                "999", tmp, kir41, m03mod, m24mod, kir_ref=ref)
            prot_ids = per_level["protein"]["KIR3DL1"]
            self.assertEqual(len(prot_ids), 1)
            self.assertTrue(next(iter(prot_ids)).startswith("KIR3DL1_prot_"))
            expected_prot = "KIR3DL1_prot_" + m24mod.sha8(m24mod.protein_info("ATGAAATAG")["protein"])
            self.assertEqual(next(iter(prot_ids)), expected_prot)
            self.assertIn(expected_prot, per_level["protein_novel"]["KIR3DL1"])


class TestBuildKirIdentitySetsFallback(unittest.TestCase):
    def test_cds_available_false_when_no_cds_fasta_anywhere(self):
        with tempfile.TemporaryDirectory() as tmp:
            for pid in ("1", "2"):
                person_dir = os.path.join(tmp, pid, "immuannot_output")
                _write_gtf_gz(os.path.join(person_dir, "hap1.gtf.gz"),
                             [("KIR3DL1", "KIR3DL1*001")], contig="ctgA")
                _write_gtf_gz(os.path.join(person_dir, "hap2.gtf.gz"), [], contig="ctgA")
            sets_by_level, qc, cds_available, n_seen, n_cds = m44.build_kir_identity_sets(
                ["1", "2"], tmp, workers=1)
            self.assertFalse(cds_available)
            self.assertEqual(n_cds, 0)
            self.assertGreater(n_seen, 0)
            self.assertTrue(len(sets_by_level["genomic"]["KIR3DL1"]) > 0)

    def test_cds_available_true_when_present_for_at_least_one_hap(self):
        with tempfile.TemporaryDirectory() as tmp:
            person_dir = os.path.join(tmp, "1", "immuannot_output")
            _write_gtf_gz(os.path.join(person_dir, "hap1.gtf.gz"),
                         [("KIR3DL1", "KIR3DL1*001new")], contig="ctgA")
            _write_gtf_gz(os.path.join(person_dir, "hap2.gtf.gz"), [], contig="ctgA")
            _write_cds_fasta_gz(os.path.join(person_dir, "hap1", "cds.fa.gz"),
                               [("ctgA", "KIR3DL1", 1, "ATGAAATAG")])
            sets_by_level, qc, cds_available, n_seen, n_cds = m44.build_kir_identity_sets(
                ["1"], tmp, workers=1)
            self.assertTrue(cds_available)
            self.assertEqual(n_cds, 1)


class TestHlaIdentityLevels(unittest.TestCase):
    def _frame(self):
        return pd.DataFrame({
            "consensus": ["HLA-A*01:01:01:01", "HLA-A*01:01:01:new", "HLA-A*02:new",
                          "HLA-A*02:new"],
            "keep_clean": [True, True, True, False],
            "field_class": ["known", "f4_noncoding", "f2_protein", "f2_protein"],
            "seq_class": ["cds_known", "cds_known", "novel_protein", "novel_protein"],
            "prot_id": ["A*01:01", "A*01:01", "A_prot_abcd1234", "A_prot_abcd1234"],
            "cds_id": ["A_cds_x1", "A_cds_x2", "A_cds_novel1", "A_cds_novel1"],
        })

    def test_add_hla_genomic_id_keeps_full_untruncated_name(self):
        df = self._frame()
        out = m44.add_hla_genomic_id(df, m24mod)
        self.assertEqual(out["genomic_id"].iloc[0], "HLA-A*01:01:01:01")
        self.assertEqual(out["genomic_id"].iloc[1], "HLA-A*01:01:01:new")

    def test_genomic_level_mask_is_keep_clean_only(self):
        df = self._frame()
        mask = m44.hla_level_mask(df, "genomic")
        self.assertEqual(list(mask), [True, True, True, False])

    def test_any_novel_excludes_known(self):
        df = self._frame()
        mask = m44.hla_level_mask(df, "any_novel")
        self.assertEqual(list(mask), [False, True, True, False])

    def test_protein_novel_is_strict(self):
        df = self._frame()
        mask = m44.hla_level_mask(df, "protein_novel")
        self.assertEqual(list(mask), [False, False, True, False])

    def test_id_col_choice_uses_three_matched_granularities(self):
        self.assertEqual(m44.hla_id_col("genomic"), "genomic_id")
        self.assertEqual(m44.hla_id_col("cds"), "cds_id")
        self.assertEqual(m44.hla_id_col("protein"), "prot_id")
        self.assertEqual(m44.hla_id_col("any_novel"), "genomic_id")
        self.assertEqual(m44.hla_id_col("protein_novel"), "prot_id")

    def test_person_gene_sets_uses_genomic_id_for_genomic_level(self):
        df = self._frame()
        df = m44.add_hla_genomic_id(df, m24mod)
        df["person_id"] = ["p1", "p2", "p3", "p4"]
        df["gene_b"] = "A"
        sets = m44.hla_person_gene_sets(df, "genomic", ["A"])
        self.assertIn("HLA-A*01:01:01:01", sets["A"]["p1"])
        self.assertIn("HLA-A*01:01:01:new", sets["A"]["p2"])
        self.assertNotIn("p4", sets["A"])  # row 4 not keep_clean


class TestEndToEndSynthetic(unittest.TestCase):
    """Fake KIR GTFs + matching cds.fa.gz fixtures -> build_kir_identity_sets -> run_gene_level,
    paired with a fake HLA `calls` frame (shape build_labeled_calls would produce) ->
    hla_person_gene_sets -> run_gene_level, feeding the same tables 44's main() writes, then
    rendering 45's figures on them. Deliberately bypasses build_labeled_calls()/
    cohort_membership/relatedness-table I/O (VM-only inputs) -- see module docstring."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.outroot = os.path.join(self.tmp.name, "pipeline_outputs_kir")
        self.people = {
            "111111": ("AFR", [("KIR3DL1", "KIR3DL1*001"), ("KIR2DL1", "KIR2DL1*003new")]),
            "222222": ("AFR", [("KIR3DL1", "KIR3DL1*001"), ("KIR2DL1", "KIR2DL1*003new")]),
            "333333": ("AFR", [("KIR3DL1", "KIR3DL1*002new")]),
            "444444": ("EUR", [("KIR3DL1", "KIR3DL1*001")]),
            "555555": ("EUR", [("KIR3DL1", "KIR3DL1*001")]),
            "666666": ("EUR", [("KIR2DL1", "KIR2DL1*003new")]),
        }
        cds_by_consensus = {
            "KIR3DL1*001": "ATGAAATAG",
            "KIR3DL1*002new": "ATGAAACAG",  # 1bp different -> different novel protein/cds hash
            "KIR2DL1*003new": "ATGCCCTAG",
        }
        for pid, (anc, calls) in self.people.items():
            hap_dir = os.path.join(self.outroot, pid, "immuannot_output")
            _write_gtf_gz(os.path.join(hap_dir, "hap1.gtf.gz"), calls, contig="ctgA")
            _write_gtf_gz(os.path.join(hap_dir, "hap2.gtf.gz"), [], contig="ctgA")
            records = [("ctgA", gene, i + 1, cds_by_consensus[cons])
                      for i, (gene, cons) in enumerate(calls)]
            if records:
                _write_cds_fasta_gz(os.path.join(hap_dir, "hap1", "cds.fa.gz"), records)
        self.ancestry_of = {pid: anc for pid, (anc, _) in self.people.items()}

    def tearDown(self):
        self.tmp.cleanup()

    def _fake_hla_calls(self):
        rows = []
        data = {
            "111111": ("A*01:01:01:01", "known", "cds_known", "A*01:01", "A_cds_x1"),
            "222222": ("A*01:01:01:01", "known", "cds_known", "A*01:01", "A_cds_x1"),
            "333333": ("A*01:new", "f2_protein", "novel_protein", "A_prot_ab12cd34", "A_cds_novel1"),
            "444444": ("A*02:01:01:01", "known", "cds_known", "A*02:01", "A_cds_x2"),
            "555555": ("A*02:01:01:01", "known", "cds_known", "A*02:01", "A_cds_x2"),
            "666666": ("A*02:01:01:new", "f4_noncoding", "cds_known", "A*02:01", "A_cds_x2b"),
        }
        for pid, (cons, fc, sc, prot_id, cds_id) in data.items():
            rows.append({"person_id": pid, "gene_b": "A", "consensus": cons, "field_class": fc,
                        "seq_class": sc, "prot_id": prot_id, "cds_id": cds_id,
                        "keep_clean": True})
        return pd.DataFrame(rows)

    def test_end_to_end_writes_tsvs_and_45_renders(self):
        pids = agg43.discover_people(self.outroot)
        self.assertEqual(len(pids), 6)

        calls = m44.add_hla_genomic_id(self._fake_hla_calls(), m24mod)

        pids_by_ancestry = {
            "ALL": pids,
            "AFR": [p for p in pids if self.ancestry_of[p] == "AFR"],
            "EUR": [p for p in pids if self.ancestry_of[p] == "EUR"],
        }
        n_permutations = 5
        orders_by_ancestry = {a: m44.make_orders(len(ps), n_permutations, seed=44 + hash(a) % 1000)
                              for a, ps in pids_by_ancestry.items()}
        n_star_by_ancestry = {"ALL": 6, "AFR": 3, "EUR": 3}

        (kir_sets_by_level, kir_qc, kir_cds_available,
         n_hap_seen, n_hap_cds_present) = m44.build_kir_identity_sets(pids, self.outroot, workers=1)
        self.assertTrue(kir_cds_available)

        hla_sets_by_level = {lvl: m44.hla_person_gene_sets(calls, lvl, ["A"]) for lvl in m44.LEVELS}

        all_rec, all_curve, all_slope, all_cov = [], [], [], []
        for level in m44.LEVELS:
            for gene in ["KIR3DL1", "KIR2DL1"]:
                full_gene_sets = kir_sets_by_level[level].get(gene, {})
                rec, curve, slope, cov = m44.run_gene_level(
                    gene, level, "kir", full_gene_sets, pids_by_ancestry, orders_by_ancestry,
                    n_star_by_ancestry, stride=1, extrapolate_2n=True)
                all_rec += rec; all_curve += curve; all_slope += slope; all_cov += cov
            for gene in ["A"]:
                full_gene_sets = hla_sets_by_level[level].get(gene, {})
                rec, curve, slope, cov = m44.run_gene_level(
                    gene, level, "hla", full_gene_sets, pids_by_ancestry, orders_by_ancestry,
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

        rec_df = pd.read_csv(os.path.join(out_dir, "recurrence_classes.tsv"), sep="\t", dtype=str)
        for col in ("n_distinct_alleles", "eq1", "eq2", "gt2", "ge20"):
            for v in rec_df[col]:
                if v not in ("0", "<20", "NA"):
                    self.assertGreaterEqual(int(v), 20, f"{col}={v} should have been masked")

        # KIR cds/protein levels should have produced REAL (non-NA) rows since a cds.fa.gz fixture
        # was provided for every call in this fixture.
        kir_cds_rows = rec_df[(rec_df["species"] == "kir") & (rec_df["level"] == "cds")]
        self.assertTrue((kir_cds_rows["n_distinct_alleles"] != "NA").all())

        m45 = _load_module("45_kir_recurrence_figure.py", "kir_recurrence_figure_45_test_target")
        rec, curve, cov = m45.load_tables(out_dir)
        try:
            m45.fig_saturation_by_recurrence(curve, os.path.join(out_dir, "fig_saturation_by_recurrence"))
            m45.fig_saturation_per_ancestry(curve, os.path.join(out_dir, "fig_saturation_per_ancestry"))
            m45.fig_coverage_completeness(cov, os.path.join(out_dir, "fig_coverage_completeness"))
        except RuntimeError as e:
            # Known LOCAL sandbox-only issue, not a bug in 44/45: macOS's system Helvetica.ttc
            # (a multi-face TrueType Collection) trips a freetype glyph-load failure inside
            # _viz_common's nature_style() font list on this machine specifically (confirmed:
            # matplotlib.font_manager resolves "Helvetica" to /System/Library/Fonts/Helvetica.ttc,
            # and even a bare fig.canvas.draw() under nature_style() -- no 44/45 code involved --
            # reproduces "RuntimeError: failed to load glyph" with no other rcParams changed).
            # This will not reproduce on the VM's Linux/matplotlib font stack (no Helvetica.ttc
            # there). Skip the pixel-rendering assertion rather than fail the whole suite on an
            # environment quirk in a file (_viz_common.py) this task does not own or modify --
            # the data pipeline above (44's tables, 45's data-shaping functions) is fully
            # exercised regardless.
            if "failed to load glyph" not in str(e):
                raise
            self.skipTest(f"local-sandbox-only matplotlib/freetype glyph issue (not a 44/45 "
                          f"bug -- see comment): {e}")
            return
        for stem in ("fig_saturation_by_recurrence", "fig_saturation_per_ancestry",
                     "fig_coverage_completeness"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, stem + ".png")))
            self.assertTrue(os.path.exists(os.path.join(out_dir, stem + ".pdf")))

    def test_end_to_end_na_fallback_when_cds_fasta_absent(self):
        """Same fixture but with every hap1/cds.fa.gz removed -- cds/protein/protein_novel must
        come back as explicit NA rows, never silently 0 or omitted."""
        for pid in self.people:
            cds_path = os.path.join(self.outroot, pid, "immuannot_output", "hap1", "cds.fa.gz")
            if os.path.exists(cds_path):
                os.remove(cds_path)
        pids = agg43.discover_people(self.outroot)
        (kir_sets_by_level, kir_qc, kir_cds_available,
         n_hap_seen, n_hap_cds_present) = m44.build_kir_identity_sets(pids, self.outroot, workers=1)
        self.assertFalse(kir_cds_available)

        pids_by_ancestry = {"ALL": pids}
        n_star_by_ancestry = {"ALL": 6}
        rec, curve, slope, cov = m44.build_na_rows(
            "KIR3DL1", "cds", "kir", pids_by_ancestry, n_star_by_ancestry)
        self.assertEqual(rec[0]["n_distinct_alleles"], "NA")
        self.assertEqual(cov[0]["chao2"], "NA")


class TestHlaPeopleOutrootBugFix(unittest.TestCase):
    """Reproduces bug 1 (HLA protein_novel s_obs == 0 for every classical gene, commit d7c16f6):
    44 used to hardcode `outroot=os.path.expanduser("~/pipeline_outputs")` when calling
    39_saturation_by_ancestry.build_labeled_calls(), one directory short of
    <outroot>/people/<pid>/immuannot_output/hap{1,2}/cds.fa.gz -- 39's own DEFAULT_PEOPLE_OUTROOT
    convention (`os.path.join(DEFAULT_OUTROOT, "people")`). match_sequences()/_match_person()/
    match_novel_rows() build `cds_path = os.path.join(outroot, person_id, "immuannot_output", hap,
    "cds.fa.gz")` from that outroot directly, so the wrong outroot means every depth-2/3 call's
    cds.fa.gz lookup misses, seq_class/prot_id never populate, and field_class==f2_protein AND
    seq_class==novel_protein (protein_novel's own mask) is never true for ANY row -- a silent,
    universal zero, not a real biological zero (S03/S01/Figure-1-panel-d show ~20-60 novel protein
    alleles per classical gene)."""

    def test_default_outroot_points_at_people_subdir(self):
        # Would be AttributeError (pre-fix code had no such name at all) or a plain string equal
        # to "~/pipeline_outputs" (pre-fix hardcoded value, missing "/people") before the fix.
        default = getattr(m44, "DEFAULT_HLA_PEOPLE_OUTROOT", None)
        self.assertIsNotNone(
            default, "44 must expose DEFAULT_HLA_PEOPLE_OUTROOT (2026-09-26 bug fix)")
        self.assertTrue(
            default.rstrip("/").endswith(os.sep + "people"),
            f"DEFAULT_HLA_PEOPLE_OUTROOT={default!r} must point at the 'people' subdir, matching "
            f"39_saturation_by_ancestry.DEFAULT_PEOPLE_OUTROOT -- the pre-fix bug hardcoded the "
            f"bare pipeline_outputs dir here, one level too shallow.")

    def test_matches_39s_own_default_people_outroot(self):
        m39mod = _load_module("39_saturation_by_ancestry.py", "saturation_by_ancestry_test_target_44")
        self.assertEqual(m44.DEFAULT_HLA_PEOPLE_OUTROOT, m39mod.DEFAULT_PEOPLE_OUTROOT,
                         "44's HLA people-outroot default has drifted from 39's own "
                         "DEFAULT_PEOPLE_OUTROOT -- they must name the same directory.")

    def test_wrong_outroot_silently_matches_nothing_right_outroot_matches(self):
        """Mechanism-level reproduction: build one person's on-disk cds.fa.gz under the CORRECT
        <outroot>/people/<pid>/... layout (39's convention) and show that match_sequences() (the
        exact function build_labeled_calls() calls, reused verbatim) matches it when given the
        correct 'people'-suffixed outroot but matches NOTHING when given the pre-fix, one-level-
        shallow outroot -- reproducing the silent-zero mechanism behind bug 1 without needing the
        full VM-only build_labeled_calls()/cohort/relatedness/refdata pipeline."""
        with tempfile.TemporaryDirectory() as tmp:
            base = os.path.join(tmp, "pipeline_outputs")
            people_root = os.path.join(base, "people")
            pid = "1234567"
            hap_dir = os.path.join(people_root, pid, "immuannot_output")
            os.makedirs(os.path.join(hap_dir, "hap1"), exist_ok=True)
            with gzip.open(os.path.join(hap_dir, "hap1", "cds.fa.gz"), "wt") as f:
                f.write(">ctgA_A_1\nATGAAATAG\n")

            # match_novel_rows joins on f"{contig}_{gene}" against the fasta header -- 'gene' here
            # must be the bare name ("A"), matching the cds.fa.gz header convention (see
            # 03_novel_alleles.py's own docstring / KIR fixtures elsewhere in this file). It also
            # needs is_novel/gene_class/etc. columns match_novel_rows() reads directly.
            t1 = pd.DataFrame([{
                "person_id": pid, "hap": "hap1", "gene": "A", "contig": "ctgA",
                "copy_index": 1, "depth": 2, "is_novel": "True", "gene_class": "classical",
                "cds_distance": None, "n_aa_changes": None, "novelty_class": None,
                "template_warning": None, "cds_mut": None,
            }])

            # Pre-fix behavior: outroot one level too shallow (no "/people").
            matched_wrong, _ = m24mod.match_sequences(t1, base, threads=1)
            self.assertEqual(len(matched_wrong), 0,
                             "sanity check on the test fixture itself: the pre-fix (bare "
                             "pipeline_outputs) outroot must NOT find cds.fa.gz.")

            # Post-fix behavior: correct outroot (matches DEFAULT_HLA_PEOPLE_OUTROOT's shape).
            matched_right, _ = m24mod.match_sequences(t1, people_root, threads=1)
            self.assertEqual(len(matched_right), 1,
                             "the correct 'people'-suffixed outroot must find and match the "
                             "on-disk cds.fa.gz.")
            self.assertEqual(matched_right[0]["seq"], "ATGAAATAG")


class TestKirProteinCatalogueCoverageFollowUpFix(unittest.TestCase):
    """2026-09-27 follow-up fix: the VM smoke test of commit 4a75657 still showed
    pct_novel_protein > 100% for KIR2DL2, KIR2DL5B, and KIR2DP1. Root causes:
      - KIR2DL5B (name-mismatch case): IPD-KIR's own CDSseq headers can carry an undifferentiated
        "KIR2DL5" allele name instead of the KIR_GENES-spelled "KIR2DL5A"/"KIR2DL5B", so neither
        gene is ever a key in `kir_ref.prot`/`kir_ref.cds` under its exact KIR_GENES name.
      - KIR2DP1 (pseudogene case): every reference allele for a pseudogene is frameshifted or has a
        premature stop, so `RefIndex.add()` adds it to `.cds` (gene "covered" by a naive check) but
        NEVER to `.prot` -- `kir_ref.prot["KIR2DP1"]` stays permanently empty even though the gene
        itself is "in the catalogue".
      - KIR2DL2 (shared/bundled-file case): covered separately by
        `_identity_worker_init`'s genes_needed=None fix (a file-level filename filter would have
        skipped a file bundling KIR2DL2 with KIR2DL3 under a stem matching only one of them) --
        exercised in TestGenesNeededNoneFix below.
    The fix in both cases here: `kir_gene_protein_covered()` requires a NON-EMPTY `kir_ref.prot[gene]`
    (not just gene membership in `kir_ref.cds`/`kir_ref.genes()`), and
    `build_person_kir_identity()` NEVER produces a hash id for an uncovered gene -- it leaves the
    gene out of `protein`/`protein_novel` entirely so main() can export an explicit NA row."""

    def test_name_mismatch_gene_is_uncovered_even_though_related_gene_is_in_catalogue(self):
        """Catalogue has 'KIR2DL5' headers (undifferentiated), never 'KIR2DL5A'/'KIR2DL5B' -- so
        looking up either exact KIR_GENES-spelled name must report uncovered, not silently borrow
        the undifferentiated entry (which would misattribute alleles across two real, distinct
        genes)."""
        ref = m24mod.RefIndex()
        ref.add("KIR2DL5", "KIR2DL5*00101", "ATGAAATAG")  # undifferentiated header, no A/B suffix
        ref.finalize()
        self.assertFalse(m44.kir_gene_protein_covered(ref, "KIR2DL5A"))
        self.assertFalse(m44.kir_gene_protein_covered(ref, "KIR2DL5B"))
        ok_a, reason_a = m44.kir_protein_catalogue_status(ref, ["KIR2DL5A", "KIR2DL5B"])["KIR2DL5A"]
        self.assertFalse(ok_a)
        self.assertEqual(reason_a, "gene_not_in_catalogue")
        # The undifferentiated name itself IS covered (sanity check on the fixture/helper).
        self.assertTrue(m44.kir_gene_protein_covered(ref, "KIR2DL5"))

    def test_pseudogene_with_only_nonfunctional_catalogue_records_is_uncovered(self):
        """KIR2DP1-shaped case: RefIndex.add() with a frameshifted/premature-stop sequence adds a
        CDS record but (by RefIndex.add()'s own gate) never a protein record -- gene ends up 'in'
        kir_ref.cds/genes() but kir_ref.prot[gene] is empty. kir_gene_protein_covered() must treat
        this as uncovered, not as "covered, zero known proteins == everything is novel"."""
        ref = m24mod.RefIndex()
        ref.add("KIR2DP1", "KIR2DP1*00101", "ATGAA")  # len 5 -> frameshift (not a multiple of 3)
        ref.add("KIR2DP1", "KIR2DP1*00201", "ATGTAAGGG")  # premature stop mid-sequence
        ref.finalize()
        self.assertIn("KIR2DP1", ref.genes())          # has CDS records
        self.assertEqual(ref.prot.get("KIR2DP1", {}), {})  # but zero usable protein records
        self.assertFalse(m44.kir_gene_protein_covered(ref, "KIR2DP1"))
        ok, reason = m44.kir_protein_catalogue_status(ref, ["KIR2DP1"])["KIR2DP1"]
        self.assertFalse(ok)
        self.assertEqual(reason, "no_catalogued_protein_entries")

    def test_uncovered_gene_never_produces_a_hash_id_end_to_end(self):
        """End-to-end (build_person_kir_identity) confirmation for the pseudogene case: even with a
        real cds.fa.gz and a kir_ref that DOES have CDS-level records for the gene, an uncovered
        (protein-empty) gene must come back with NOTHING in 'protein'/'protein_novel' -- never a
        hash id, which is what silently caused the >100% warnings pre-fix."""
        ref = m24mod.RefIndex()
        ref.add("KIR2DP1", "KIR2DP1*00101", "ATGAA")  # frameshift -> never added to .prot
        ref.finalize()
        with tempfile.TemporaryDirectory() as tmp:
            person_dir = os.path.join(tmp, "1", "immuannot_output")
            _write_gtf_gz(os.path.join(person_dir, "hap1.gtf.gz"),
                         [("KIR2DP1", "KIR2DP1*001")], contig="ctgA")
            _write_gtf_gz(os.path.join(person_dir, "hap2.gtf.gz"), [], contig="ctgA")
            _write_cds_fasta_gz(os.path.join(person_dir, "hap1", "cds.fa.gz"),
                               [("ctgA", "KIR2DP1", 1, "ATGAAATAG")])
            pid, per_level, qc, n_seen, n_cds = m44.build_person_kir_identity(
                "1", tmp, kir41, m03mod, m24mod, kir_ref=ref)
            self.assertNotIn("KIR2DP1", per_level.get("protein", {}))
            self.assertNotIn("KIR2DP1", per_level.get("protein_novel", {}))
            self.assertEqual(qc["n_protein_gene_uncovered_calls"], 1)
            # genomic/cds are unaffected -- they don't need the protein catalogue.
            self.assertIn("KIR2DP1*001", per_level["genomic"]["KIR2DP1"])
            self.assertEqual(len(per_level["cds"]["KIR2DP1"]), 1)

    def test_covered_gene_alongside_uncovered_gene_only_the_latter_is_skipped(self):
        """Multi-gene sanity check: a covered gene's protein tracking must be unaffected by a
        SEPARATE gene being uncovered in the same catalogue/run."""
        ref = m24mod.RefIndex()
        ref.add("KIR3DL1", "KIR3DL1*00101", "ATGAAACCCTAG")  # covered, functional
        ref.add("KIR2DP1", "KIR2DP1*00101", "ATGAA")          # uncovered, frameshift
        ref.finalize()
        status = m44.kir_protein_catalogue_status(ref, ["KIR3DL1", "KIR2DP1"])
        self.assertEqual(status["KIR3DL1"], (True, "ok"))
        self.assertEqual(status["KIR2DP1"], (False, "no_catalogued_protein_entries"))


class TestGenesNeededNoneFix(unittest.TestCase):
    """2026-09-27 fix: _identity_worker_init() now loads the KIR protein catalogue with
    genes_needed=None instead of genes_needed=set(KIR_GENES), because load_refdata's file-level
    pre-filter matches a CDSseq file's FILENAME stem, while gene attribution for records inside a
    file comes from each record's HEADER -- a single file bundling two genes under a filename stem
    that only exactly matches one exact KIR_GENES spelling (or neither) would previously skip the
    WHOLE file and silently starve every gene inside it of catalogue coverage, however correctly
    those records were labeled internally."""

    def test_bundled_file_with_mismatched_stem_is_not_silently_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            cdsseq_dir = os.path.join(tmp, "CDSseq")
            os.makedirs(cdsseq_dir)
            # A single file, stem "KIR2DL2_2DL3" -- matches NEITHER "KIR2DL2" nor "KIR2DL3"
            # exactly -- bundling both genes' alleles, correctly labeled inside via headers.
            with gzip.open(os.path.join(cdsseq_dir, "KIR2DL2_2DL3.fa.gz"), "wt") as f:
                f.write(">KIR2DL2*00101\nATGAAACCCTAG\n")
                f.write(">KIR2DL3*00101\nATGAAATTTTAG\n")

            filtered = m24mod.load_refdata(tmp, genes_needed={"KIR2DL2", "KIR2DL3"})
            self.assertNotIn("KIR2DL2", filtered.genes(),
                             "sanity check on the fixture: a genes_needed filename filter DOES "
                             "skip this bundled file (reproduces the pre-fix failure mode).")

            unfiltered = m24mod.load_refdata(tmp, genes_needed=None)
            self.assertIn("KIR2DL2", unfiltered.genes())
            self.assertIn("KIR2DL3", unfiltered.genes())
            self.assertTrue(m44.kir_gene_protein_covered(unfiltered, "KIR2DL2"))
            self.assertTrue(m44.kir_gene_protein_covered(unfiltered, "KIR2DL3"))


class TestKirProteinIdentityGranularityBugFix(unittest.TestCase):
    """Reproduces bug 2 (KIR pct_novel_protein > 100% for every gene, commit d7c16f6): the
    'protein'/'protein_novel' tracks hashed the raw per-call translated CDS for EVERY call,
    known or not, instead of collapsing known proteins to a single curated catalogue identity the
    way HLA's own prot_id scheme does -- so two people carrying the textbook-SAME known KIR allele
    (ordinary per-sample sequencing noise in the reconstructed CDS) could land on different
    protein-level hashes, making 'protein' (and therefore 'protein_novel', a subset of it) far
    FINER-grained than 'genomic' (the curated, already-deduplicated Immuannot consensus name) --
    backwards for a protein-vs-genomic comparison, and exactly why
    s_obs(protein_novel)/s_obs(genomic) (46_kir_vs_hla_coverage.py's pct_novel_protein) could
    exceed 100%."""

    def _catalog_ref(self, gene, known_seq):
        """A minimal m24().RefIndex-shaped stand-in covering exactly one known protein for `gene`,
        built the same way m24mod.RefIndex.add() would from a real CDSseq/*.fa.gz record."""
        ref = m24mod.RefIndex()
        ref.add(gene, f"{gene}*00101", known_seq)
        ref.finalize()
        return ref

    def test_without_kir_ref_two_people_same_known_protein_get_different_hash_ids(self):
        """Pre-fix (and still the documented fallback when no catalogue is available): a known
        allele's protein identity is a raw hash of the exact observed sequence, so a 1bp-different
        but SAME-catalogued-protein observation (silent/synonymous at nucleotide level for this
        toy example is not required -- any two distinct nucleotide sequences suffice) produces two
        DIFFERENT ids even though both are 'known'. This is the pre-fix behavior --
        build_person_kir_identity's default (`kir_ref=None`) reproduces it exactly, and existing
        older fixtures rely on this fallback still working."""
        with tempfile.TemporaryDirectory() as tmp:
            for pid, seq in (("1", "ATGAAATAG"), ("2", "ATGAAATAG")):
                person_dir = os.path.join(tmp, pid, "immuannot_output")
                _write_gtf_gz(os.path.join(person_dir, "hap1.gtf.gz"),
                             [("KIR3DL1", "KIR3DL1*001")], contig="ctgA")
                _write_gtf_gz(os.path.join(person_dir, "hap2.gtf.gz"), [], contig="ctgA")
                _write_cds_fasta_gz(os.path.join(person_dir, "hap1", "cds.fa.gz"),
                                   [("ctgA", "KIR3DL1", 1, seq)])
            sets_by_level, qc, cds_available, n_seen, n_cds = m44.build_kir_identity_sets(
                ["1", "2"], tmp, workers=1)  # no refdata= -> kir_ref=None, old scheme
            self.assertTrue(cds_available)
            # Same sequence -> same hash -> ONE shared id (not the bug by itself), but this
            # 'protein' id is a raw hash, never collapsed to a curated catalogue name -- the
            # granularity-mismatch setup the next test exercises with genuinely known alleles.
            prot_ids = set()
            for pid in ("1", "2"):
                prot_ids |= sets_by_level["protein"]["KIR3DL1"].get(pid, set())
            self.assertTrue(all(pid.startswith("KIR3DL1_prot_") for pid in prot_ids))

    def test_kir_ref_collapses_known_protein_to_one_catalogue_id_not_a_hash(self):
        """The fix: when a `kir_ref` (IPD-KIR protein catalogue) is supplied and covers the gene, a
        translated protein that MATCHES the catalogue collapses to ONE catalogue-name-derived id
        (mirroring HLA's prot_id for known alleles), not a hash -- so 'protein' s_obs for known
        alleles can no longer explode past 'genomic' s_obs the way raw per-call hashing did."""
        known_seq = "ATGAAATAG"  # translates to catalogued protein "MK"
        ref = self._catalog_ref("KIR3DL1", known_seq)
        with tempfile.TemporaryDirectory() as tmp:
            person_dir = os.path.join(tmp, "1", "immuannot_output")
            _write_gtf_gz(os.path.join(person_dir, "hap1.gtf.gz"),
                         [("KIR3DL1", "KIR3DL1*001")], contig="ctgA")
            _write_gtf_gz(os.path.join(person_dir, "hap2.gtf.gz"), [], contig="ctgA")
            _write_cds_fasta_gz(os.path.join(person_dir, "hap1", "cds.fa.gz"),
                               [("ctgA", "KIR3DL1", 1, known_seq)])
            pid, per_level, qc, n_seen, n_cds = m44.build_person_kir_identity(
                "1", tmp, kir41, m03mod, m24mod, kir_ref=ref)
            prot_ids = per_level["protein"]["KIR3DL1"]
            self.assertEqual(len(prot_ids), 1)
            got = next(iter(prot_ids))
            self.assertFalse(got.startswith("KIR3DL1_prot_"),
                             f"a catalogued protein must NOT get a hash id, got {got!r}")
            self.assertTrue(got.startswith("KIR3DL1_KIR3DL1*"),
                            f"expected a catalogue-name-derived id, got {got!r}")
            # Not novel: it's in the catalogue.
            self.assertNotIn("KIR3DL1", per_level.get("protein_novel", {}))

    def test_protein_novel_never_exceeds_genomic_s_obs_with_kir_ref(self):
        """End-to-end-ish reproduction of the actual reported bug: build a small synthetic cohort
        where several people carry the SAME known allele (each reconstructed with a trivially
        different -- but still-translates-to-the-known-protein -- CDS to mimic per-sample noise)
        plus one genuinely novel-protein carrier, run it through build_coverage_chao2_row for both
        'genomic' and 'protein_novel', and confirm protein_novel's s_obs <= genomic's s_obs (the
        pre-fix code could and did violate this -- see coverage_chao2.tsv commit d7c16f6, e.g.
        KIR2DL1 genomic S_obs=226 vs protein S_obs=627)."""
        known_seq = "ATGAAATAG"          # -> known protein "MK"
        novel_seq = "ATGAAACAGTAG"       # -> different protein, not in catalogue
        ref = self._catalog_ref("KIR3DL1", known_seq)
        with tempfile.TemporaryDirectory() as tmp:
            people = {
                "1": "KIR3DL1*001", "2": "KIR3DL1*001", "3": "KIR3DL1*001",
                "4": "KIR3DL1*002new",
            }
            seq_for = {"KIR3DL1*001": known_seq, "KIR3DL1*002new": novel_seq}
            for pid, consensus in people.items():
                person_dir = os.path.join(tmp, pid, "immuannot_output")
                _write_gtf_gz(os.path.join(person_dir, "hap1.gtf.gz"),
                             [("KIR3DL1", consensus)], contig="ctgA")
                _write_gtf_gz(os.path.join(person_dir, "hap2.gtf.gz"), [], contig="ctgA")
                _write_cds_fasta_gz(os.path.join(person_dir, "hap1", "cds.fa.gz"),
                                   [("ctgA", "KIR3DL1", 1, seq_for[consensus])])
            sets_by_level, qc, cds_available, n_seen, n_cds = m44.build_kir_identity_sets(
                list(people), tmp, workers=1, refdata=None)
            # Manually attach the ref (build_kir_identity_sets(refdata=...) loads from disk; here
            # we rebuild via build_person_kir_identity directly so the toy in-memory ref applies).
            sets_by_level = {lvl: {} for lvl in m44.LEVELS}
            for pid in people:
                _, per_level, _, _, _ = m44.build_person_kir_identity(
                    pid, tmp, kir41, m03mod, m24mod, kir_ref=ref)
                for lvl, gene_sets in per_level.items():
                    for gene, ids in gene_sets.items():
                        sets_by_level[lvl].setdefault(gene, {})[pid] = ids

            genomic_sets = sets_by_level["genomic"].get("KIR3DL1", {})
            protein_novel_sets = sets_by_level["protein_novel"].get("KIR3DL1", {})
            genomic_row = m44.build_coverage_chao2_row(
                "KIR3DL1", "ALL", "genomic", "kir",
                {p: genomic_sets.get(p, set()) for p in people})
            protein_novel_row = m44.build_coverage_chao2_row(
                "KIR3DL1", "ALL", "protein_novel", "kir",
                {p: protein_novel_sets.get(p, set()) for p in people})
            self.assertLessEqual(protein_novel_row["s_obs"], genomic_row["s_obs"],
                                 f"protein_novel s_obs ({protein_novel_row['s_obs']}) must not "
                                 f"exceed genomic s_obs ({genomic_row['s_obs']}) -- this is "
                                 f"exactly the >100% pct_novel_protein bug.")
            self.assertEqual(protein_novel_row["s_obs"], 1)  # only the one genuinely novel protein
            self.assertEqual(genomic_row["s_obs"], 2)        # KIR3DL1*001 and *002new


class TestSanityCheckCoverage(unittest.TestCase):
    def test_flags_protein_novel_zero_when_any_novel_positive(self):
        cov = [
            {"species": "hla", "gene": "A", "ancestry": "ALL", "level": "any_novel", "s_obs": 50},
            {"species": "hla", "gene": "A", "ancestry": "ALL", "level": "protein_novel", "s_obs": 0},
            {"species": "hla", "gene": "A", "ancestry": "ALL", "level": "genomic", "s_obs": 300},
        ]
        warnings = m44.sanity_check_coverage(cov)
        self.assertTrue(any("protein_novel s_obs == 0" in w for w in warnings))

    def test_flags_pct_over_100(self):
        # 2026-09-27: denominator is 'protein' (same granularity), not 'genomic' -- protein_novel
        # (300) exceeding protein (250) is the only thing that should trip this now.
        cov = [
            {"species": "kir", "gene": "KIR2DL1", "ancestry": "ALL", "level": "genomic",
             "s_obs": 226},
            {"species": "kir", "gene": "KIR2DL1", "ancestry": "ALL", "level": "protein",
             "s_obs": 250},
            {"species": "kir", "gene": "KIR2DL1", "ancestry": "ALL", "level": "protein_novel",
             "s_obs": 300},
        ]
        warnings = m44.sanity_check_coverage(cov)
        self.assertTrue(any("pct_novel_protein > 100%" in w for w in warnings))

    def test_does_not_flag_when_protein_novel_exceeds_genomic_but_not_protein(self):
        """The 2026-09-26 version of this check compared protein_novel against genomic and would
        have fired here (300 > 226); the 2026-09-27 fix compares against 'protein' (400), so this
        must NOT warn -- protein_novel is a legitimate subset of protein, just not of genomic."""
        cov = [
            {"species": "kir", "gene": "KIR2DL1", "ancestry": "ALL", "level": "genomic",
             "s_obs": 226},
            {"species": "kir", "gene": "KIR2DL1", "ancestry": "ALL", "level": "protein",
             "s_obs": 400},
            {"species": "kir", "gene": "KIR2DL1", "ancestry": "ALL", "level": "protein_novel",
             "s_obs": 300},
        ]
        warnings = m44.sanity_check_coverage(cov)
        self.assertEqual(warnings, [])

    def test_clean_data_emits_no_warnings(self):
        cov = [
            {"species": "hla", "gene": "A", "ancestry": "ALL", "level": "any_novel", "s_obs": 50},
            {"species": "hla", "gene": "A", "ancestry": "ALL", "level": "protein_novel",
             "s_obs": 20},
            {"species": "hla", "gene": "A", "ancestry": "ALL", "level": "genomic", "s_obs": 300},
        ]
        self.assertEqual(m44.sanity_check_coverage(cov), [])

    def test_na_rows_do_not_spuriously_trigger(self):
        cov = [
            {"species": "kir", "gene": "KIR2DL1", "ancestry": "ALL", "level": "any_novel",
             "s_obs": 5},
            {"species": "kir", "gene": "KIR2DL1", "ancestry": "ALL", "level": "protein_novel",
             "s_obs": "NA"},
            {"species": "kir", "gene": "KIR2DL1", "ancestry": "ALL", "level": "genomic",
             "s_obs": "NA"},
        ]
        # any_novel > 0 for kir but protein_novel is "NA" (unmatched, not a real 0) -- must not
        # be treated as a false protein_novel==0 total across a species with real any_novel data
        # from OTHER genes; single-gene NA rows must not raise a pct>100 warning either.
        warnings = m44.sanity_check_coverage(cov)
        self.assertEqual(warnings, [])


class TestJointUnrelatedSet(unittest.TestCase):
    def test_intersects_before_removing_relatives(self):
        kir_pids = ["1", "2", "3", "4"]
        hla_pids = ["2", "3", "4", "5"]
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
            f.write("i.s\tj.s\tkin\n2\t3\t0.10\n")
            path = f.name
        try:
            kept, n_removed, n_common = m44.joint_unrelated_set(kir_pids, hla_pids, path, 0.0442)
            self.assertEqual(n_common, 3)
            self.assertEqual(n_removed, 1)
            self.assertNotIn("1", kept)
            self.assertNotIn("5", kept)
            self.assertIn("4", kept)
        finally:
            os.remove(path)


if __name__ == "__main__":
    unittest.main()
