#!/usr/bin/env python3
"""Unit tests for 48_kir_hla_ligand_cooccurrence.py (S04 WS-D), synthetic data only.

Covers:
  1. Sequence-derived assignment primitives: translate_cds() (frame offset, stop-codon
     truncation), mature_protein_from_cds() (leader stripping), classify_c1c2_seq() /
     classify_bw4_seq() (residue rules, pinned to the known C1/C2 and Bw4/Bw6 reference alleles
     named in the module docstring), using short synthetic sequences that mimic the CDSseq
     layout rather than real IPD-IMGT/HLA data.
  2. CDS reference loading + allele mapping: parse_cds_header(), load_cds_fasta_gz() (synthetic
     gzip fixtures), map_called_allele_to_cds() (exact match, two-field fallback, novel-call and
     missing-reference unresolved paths).
  3. classify_c1c2()/classify_bw4() -- the lookup-table FALLBACK: table hits, unclassified
     fallback, B*15 default (Bw6) vs. the CLI-supplied exception list.
  4. two_by_two_stats() -- disclosure masking (any cell in [1,19] blanks OR/CI/p-values but keeps
     masked cell counts as '<20', never blank/0), Haldane-Anscombe correction on a zero cell,
     and a planted enrichment producing OR > 1 with a small Fisher p-value.
  5. build_person_epitopes() end-to-end with a synthetic CDS reference: sequence path used when a
     reference resolves, lookup fallback otherwise, and the seq-vs-lookup cross-check rows.
  6. End-to-end: --synthetic run writes the four aggregate TSVs, detects the planted
     KIR2DL1 x C2 / KIR3DL1 x Bw4 enrichment, and reports zero unclassified alleles for the
     synthetic cohort's own restricted allele set.

Run: python3 scripts/hla_popgen/tests/test_48_kir_hla_ligand_cooccurrence.py
"""
import gzip
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


# ---------------------------------------------------------------------------
# Helpers for building short synthetic protein / CDS sequences that mimic the real layout
# (24-aa leader + mature protein) without using any real IPD-IMGT/HLA sequence.
# ---------------------------------------------------------------------------
_AA_TO_CODON = {  # one representative codon per amino acid (enough to build a synthetic CDS)
    "A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT", "Q": "CAA", "E": "GAA",
    "G": "GGT", "H": "CAT", "I": "ATT", "L": "CTT", "K": "AAA", "M": "ATG", "F": "TTT",
    "P": "CCT", "S": "TCT", "T": "ACT", "W": "TGG", "Y": "TAT", "V": "GTT",
}


def _make_mature(overrides, length=90, filler="A"):
    """1-based-position overrides on a filler-only mature protein of the given length."""
    arr = [filler] * length
    for pos, aa in overrides.items():
        arr[pos - 1] = aa
    return "".join(arr)


def _protein_to_nt(protein, stop=True):
    nt = "".join(_AA_TO_CODON[aa] for aa in protein)
    if stop:
        nt += "TAA"
    return nt


def _make_synthetic_cds(mature_overrides, leader_len=24, mature_length=90, frame=1,
                        junk_prefix_len=0):
    """Builds a full synthetic nucleotide CDS: [junk prefix] + [leader, 'A' * leader_len] +
    [mature protein with overrides] + [stop codon]. junk_prefix_len (in nt, < 3) tests the
    frame= offset."""
    leader = "A" * leader_len
    mature = _make_mature(mature_overrides, length=mature_length)
    nt = ("N" * junk_prefix_len) + _protein_to_nt(leader, stop=False) + _protein_to_nt(mature)
    return nt, frame


class TestTranslateCds(unittest.TestCase):
    def test_basic_translation_and_stop_codon_truncation(self):
        nt = _protein_to_nt("MAK", stop=True) + "GGG"  # trailing codon after stop must be dropped
        self.assertEqual(m.translate_cds(nt, frame=1), "MAK")

    def test_frame_offset(self):
        protein = "MAK"
        nt = "NN" + _protein_to_nt(protein, stop=True)  # 2 junk bases -> frame=3 (0-based skip 2)
        self.assertEqual(m.translate_cds(nt, frame=3), protein)

    def test_unknown_codon_is_x_not_dropped(self):
        nt = "ATG" + "NNN" + "AAA" + "TAA"
        self.assertEqual(m.translate_cds(nt, frame=1), "MXK")


class TestMatureProteinFromCds(unittest.TestCase):
    def test_leader_stripped(self):
        nt, frame = _make_synthetic_cds({80: "K"})
        mature = m.mature_protein_from_cds(nt, frame)
        self.assertEqual(len(mature), 90)
        self.assertEqual(mature[79], "K")

    def test_too_short_returns_none(self):
        nt = _protein_to_nt("MAK", stop=True)  # shorter than the 24-aa leader
        self.assertIsNone(m.mature_protein_from_cds(nt, 1))


class TestClassifyC1C2Seq(unittest.TestCase):
    """Pinned to the module docstring's worked example: mature residue 80 is N (with S77) for
    C*01:02/C*03:04/C*07:01 (C1), and K for C*02:02/C*04:01/C*05:01/C*06:02 (C2)."""

    def test_c1_pattern_asn80_ser77(self):
        mature = _make_mature({77: "S", 80: "N"})
        self.assertEqual(m.classify_c1c2_seq(mature), "C1")

    def test_c2_pattern_lys80(self):
        mature = _make_mature({80: "K"})
        self.assertEqual(m.classify_c1c2_seq(mature), "C2")

    def test_asn80_without_ser77_is_other_not_silently_c1(self):
        mature = _make_mature({77: "D", 80: "N"})
        self.assertEqual(m.classify_c1c2_seq(mature), "other")

    def test_neither_pattern_is_other(self):
        mature = _make_mature({77: "D", 80: "R"})
        self.assertEqual(m.classify_c1c2_seq(mature), "other")

    def test_too_short_is_none(self):
        self.assertIsNone(m.classify_c1c2_seq("A" * 50))

    def test_none_input_is_none(self):
        self.assertIsNone(m.classify_c1c2_seq(None))


class TestClassifyBw4Seq(unittest.TestCase):
    def test_bw4_80i(self):
        mature = _make_mature({77: "N", 80: "I", 81: "A", 82: "V", 83: "R"})
        self.assertEqual(m.classify_bw4_seq(mature), ("Bw4", "80I"))

    def test_bw4_80t(self):
        mature = _make_mature({77: "N", 80: "T", 81: "A", 82: "V", 83: "R"})
        self.assertEqual(m.classify_bw4_seq(mature), ("Bw4", "80T"))

    def test_bw6_reference_pattern(self):
        mature = _make_mature({77: "S", 80: "N", 81: "L", 82: "R", 83: "G"})
        self.assertEqual(m.classify_bw4_seq(mature), ("Bw6", None))

    def test_neither_pattern_is_other(self):
        mature = _make_mature({77: "D", 80: "Q", 81: "L", 82: "R", 83: "G"})
        self.assertEqual(m.classify_bw4_seq(mature), ("other", None))

    def test_too_short_is_none(self):
        self.assertEqual(m.classify_bw4_seq("A" * 50), (None, None))


class TestCdsHeaderAndLoading(unittest.TestCase):
    def test_parse_cds_header(self):
        allele, frame = m.parse_cds_header(">HLA-A*01:01:01:01 HLA00001 frame=1 1098bp")
        self.assertEqual(allele, "HLA-A*01:01:01:01")
        self.assertEqual(frame, 1)

    def test_parse_cds_header_frame2(self):
        allele, frame = m.parse_cds_header(">HLA-B*07:02:01:01 HLA00002 frame=2 1101bp")
        self.assertEqual(frame, 2)

    def test_load_cds_fasta_gz_roundtrip(self):
        nt, frame = _make_synthetic_cds({77: "S", 80: "N"})  # C*01:02-like, C1
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "C.fa.gz")
            with gzip.open(path, "wt") as f:
                f.write(f">HLA-C*01:02:01:01 HLA00099 frame={frame} {len(nt)}bp\n")
                f.write(nt + "\n")
            idx = m.load_cds_fasta_gz(path)
            self.assertIn("HLA-C*01:02:01:01", idx)
            seq, fr = idx["HLA-C*01:02:01:01"]
            self.assertEqual(seq, nt)
            self.assertEqual(fr, frame)

    def test_load_cds_fasta_gz_missing_file_returns_empty(self):
        self.assertEqual(m.load_cds_fasta_gz("/nonexistent/path/C.fa.gz"), {})


class TestMapCalledAlleleToCds(unittest.TestCase):
    def setUp(self):
        nt_c1, frame = _make_synthetic_cds({77: "S", 80: "N"})
        self.cds_index = {"HLA-C*01:02:01:01": (nt_c1, frame)}
        self.cds_2f = m.build_cds_two_field_index(self.cds_index)

    def test_exact_match(self):
        mature, level = m.map_called_allele_to_cds("HLA-C*01:02:01:01", self.cds_index,
                                                    self.cds_2f)
        self.assertEqual(level, "exact")
        self.assertEqual(m.classify_c1c2_seq(mature), "C1")

    def test_two_field_fallback(self):
        mature, level = m.map_called_allele_to_cds("HLA-C*01:02", self.cds_index, self.cds_2f)
        self.assertEqual(level, "two_field_fallback")
        self.assertEqual(m.classify_c1c2_seq(mature), "C1")

    def test_novel_call_is_unresolved(self):
        mature, level = m.map_called_allele_to_cds("HLA-C*01:02:01:new", self.cds_index,
                                                    self.cds_2f)
        self.assertEqual(level, "unresolved")
        self.assertIsNone(mature)

    def test_no_reference_allele_is_unresolved(self):
        mature, level = m.map_called_allele_to_cds("HLA-C*99:99", self.cds_index, self.cds_2f)
        self.assertEqual(level, "unresolved")
        self.assertIsNone(mature)

    def test_empty_cds_index_is_unresolved(self):
        mature, level = m.map_called_allele_to_cds("HLA-C*01:02:01:01", {}, {})
        self.assertEqual(level, "unresolved")
        self.assertIsNone(mature)


class TestClassifyC1C2LookupFallback(unittest.TestCase):
    def test_known_c2(self):
        self.assertEqual(m.classify_c1c2("05:01"), "C2")

    def test_known_c1(self):
        self.assertEqual(m.classify_c1c2("07:01"), "C1")

    def test_unclassified_not_silently_absorbed(self):
        self.assertEqual(m.classify_c1c2("99:99"), "unclassified")


class TestClassifyBw4LookupFallback(unittest.TestCase):
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


class TestBuildPersonEpitopesWithCds(unittest.TestCase):
    def setUp(self):
        # Two synthetic HLA-C reference alleles: one C1-pattern, one C2-pattern, both agreeing
        # with C1C2_TABLE's real-world group assignment for these two-field names.
        nt_c1, f1 = _make_synthetic_cds({77: "S", 80: "N"})   # C*01:02 is C1 in C1C2_TABLE
        nt_c2, f2 = _make_synthetic_cds({80: "K"})            # C*05:01 is C2 in C1C2_TABLE
        # A deliberately mismatched fixture: C1C2_TABLE says "07:04" is C2, but this synthetic
        # sequence is built with the C1 residue pattern -- exercises the disagreement path
        # (a synthetic-fixture test, not a real biological claim about C*07:04).
        nt_mismatch, fmis = _make_synthetic_cds({77: "S", 80: "N"})
        self.cds_by_gene = {
            "HLA-A": {},
            "HLA-B": {},
            "HLA-C": {
                "HLA-C*01:02:01:01": (nt_c1, f1),
                "HLA-C*05:01:01:01": (nt_c2, f2),
                "HLA-C*07:04:01:01": (nt_mismatch, fmis),
            },
        }

    def test_sequence_path_used_when_reference_resolves(self):
        table1_df = pd.DataFrame([
            {"person_id": "p1", "hap": "hap1", "gene": "HLA-C", "consensus": "HLA-C*01:02:01:01"},
            {"person_id": "p1", "hap": "hap2", "gene": "HLA-C", "consensus": "HLA-C*05:01:01:01"},
        ])
        epitopes_df, qc, crosscheck = m.build_person_epitopes(table1_df,
                                                              cds_by_gene=self.cds_by_gene)
        self.assertTrue(epitopes_df.loc["p1", "has_C1"])
        self.assertTrue(epitopes_df.loc["p1", "has_C2"])
        self.assertEqual(qc.get("n_C_source_sequence", 0), 2)
        self.assertEqual(qc.get("n_C_source_lookup_fallback", 0), 0)

    def test_falls_back_to_lookup_when_no_cds_dir(self):
        table1_df = pd.DataFrame([
            {"person_id": "p1", "hap": "hap1", "gene": "HLA-C", "consensus": "HLA-C*01:02:01:01"},
        ])
        epitopes_df, qc, crosscheck = m.build_person_epitopes(table1_df, cds_by_gene=None)
        self.assertTrue(epitopes_df.loc["p1", "has_C1"])
        self.assertEqual(qc.get("n_C_source_sequence", 0), 0)
        self.assertEqual(qc.get("n_C_source_lookup_fallback", 0), 1)

    def test_crosscheck_agreement_recorded(self):
        table1_df = pd.DataFrame([
            {"person_id": "p1", "hap": "hap1", "gene": "HLA-C", "consensus": "HLA-C*01:02:01:01"},
        ])
        _, _, crosscheck = m.build_person_epitopes(table1_df, cds_by_gene=self.cds_by_gene)
        rows = [r for r in crosscheck if r["two_field"] == "01:02"]
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["seq_label"], "C1")
        self.assertEqual(rows[0]["lookup_label"], "C1")
        self.assertEqual(rows[0]["agree"], True)

    def test_crosscheck_disagreement_recorded(self):
        table1_df = pd.DataFrame([
            {"person_id": "p1", "hap": "hap1", "gene": "HLA-C", "consensus": "HLA-C*07:04:01:01"},
        ])
        _, _, crosscheck = m.build_person_epitopes(table1_df, cds_by_gene=self.cds_by_gene)
        rows = [r for r in crosscheck if r["two_field"] == "07:04"]
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["seq_label"], "C1")
        self.assertEqual(rows[0]["lookup_label"], "C2")  # C1C2_TABLE's real-world assignment
        self.assertEqual(rows[0]["agree"], False)

    def test_novel_call_unresolved_and_falls_back(self):
        table1_df = pd.DataFrame([
            {"person_id": "p1", "hap": "hap1", "gene": "HLA-C", "consensus": "HLA-C*01:02:01:new"},
        ])
        epitopes_df, qc, crosscheck = m.build_person_epitopes(table1_df,
                                                              cds_by_gene=self.cds_by_gene)
        # Novel call -> sequence unresolved -> falls back to the lookup table on its two-field
        # prefix (01:02 -> C1), and is never silently dropped.
        self.assertTrue(epitopes_df.loc["p1", "has_C1"])
        self.assertEqual(qc.get("n_C_source_lookup_fallback", 0), 1)


class TestEndToEndSynthetic(unittest.TestCase):
    def test_pipeline_writes_expected_files_and_detects_planted_enrichment(self):
        (table1_df, epitopes_df, kir_presence, ancestry_of,
         qc, crosscheck_rows) = m.make_synthetic(n_people=600, seed=7)
        with tempfile.TemporaryDirectory() as tmp:
            status_path = os.path.join(tmp, "STATUS.txt")
            m.run_pipeline(table1_df, epitopes_df, kir_presence, ancestry_of, tmp,
                          n_perms=20, b15_as_bw4=frozenset(), status_path=status_path,
                          qc=qc, crosscheck_rows=crosscheck_rows)
            for fn in ("kir_hla_ligand_cooccurrence.tsv", "epitope_freq_by_ancestry.tsv",
                      "ligand_lookup_qc.tsv", "ligand_seq_vs_lookup_crosscheck.tsv"):
                self.assertTrue(os.path.exists(os.path.join(tmp, fn)), fn)
            self.assertTrue(os.path.exists(status_path))

            df = pd.read_csv(os.path.join(tmp, "kir_hla_ligand_cooccurrence.tsv"), sep="\t")
            self.assertGreater(len(df), 0)
            planted = df[(df["kir_gene"] == "KIR2DL1") & (df["ligand"] == "C2")].copy()
            planted["odds_ratio"] = pd.to_numeric(planted["odds_ratio"], errors="coerce")
            planted = planted.dropna(subset=["odds_ratio"])
            self.assertGreater(len(planted), 0)
            self.assertTrue((planted["odds_ratio"] > 1.0).all())

            qc_df = pd.read_csv(os.path.join(tmp, "ligand_lookup_qc.tsv"), sep="\t")
            self.assertIn("n_unclassified_C_allele_calls", qc_df["metric"].values)
            self.assertIn("n_two_field_groups_seq_vs_lookup_compared", qc_df["metric"].values)

    def test_no_unclassified_alleles_in_synthetic_restricted_set(self):
        (table1_df, epitopes_df, kir_presence, ancestry_of,
         qc, crosscheck_rows) = m.make_synthetic(n_people=100, seed=3)
        self.assertEqual(qc.get("n_unclassified_C_allele_calls", 0), 0)
        self.assertEqual(qc.get("n_unclassified_Bw_allele_calls", 0), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
