#!/usr/bin/env python3
"""Unit tests for 41_kir_pilot.py (WS5, S03 call #8 KIR scoping pilot), synthetic data only.

Covers the statistics that would be silently wrong and hard to notice in a real run:
  1. greedy_unrelated() -- same algorithm as 24_novelty_by_field.py's, re-tested here so a future
     edit to this script's own copy doesn't quietly diverge from the shared "unrelated" definition.
  2. parse_hap_gtf() -- novelty detection (":new" suffix on consensus) and "undetermined" handling,
     against a hand-built synthetic gtf.gz (not a real Immuannot output, but the exact attribute
     grammar documented in reference/IMMUANNOT_GTF_SPEC.md).
  3. aggregate_quality() -- framework-gene presence %, A/B content classification, and the
     KIR2DL2/KIR2DL3 + KIR3DL1/KIR3DS1 mutual-exclusivity rates, on a small synthetic 2-person
     (4-haplotype) cohort with known ground truth, run against real gzip'd gtf files on disk (the
     function's actual I/O path) so a path-construction bug would also be caught.

Run: python3 scripts/hla_popgen/tests/test_kir_pilot.py
"""
import gzip
import importlib.util
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)


def _load_module(filename, modname):
    spec = importlib.util.spec_from_file_location(modname, os.path.join(HLA_POPGEN_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


kir = _load_module("41_kir_pilot.py", "kir_pilot_test_target")


def _gtf_transcript_line(gene, consensus):
    attrs = (f'gene_id "IG{gene}"; transcript_id "IAT{gene}.1"; gene_name "{gene}"; '
             f'consensus "{consensus}"; alleles "{consensus}";')
    return f"chr19_synth\tImmuannot\ttranscript\t100\t200\t.\t+\t.\t{attrs}\n"


def _write_gtf_gz(path, gene_consensus_pairs):
    with gzip.open(path, "wt") as f:
        f.write("##synthetic test gtf\n")
        for gene, consensus in gene_consensus_pairs:
            f.write(_gtf_transcript_line(gene, consensus))


class TestGreedyUnrelated(unittest.TestCase):
    def test_removes_minimum_to_break_all_pairs_above_threshold(self):
        # A-B related (kin above threshold), C isolated. Must drop exactly one of A/B.
        pairs = [("A", "B", 0.10), ("B", "C", 0.01)]  # B-C below KIN_MIN, not a real edge
        kept, removed = kir.greedy_unrelated(["A", "B", "C"], pairs, kin_min=kir.KIN_MIN)
        self.assertEqual(len(removed), 1)
        self.assertIn(removed[0], ("A", "B"))
        self.assertEqual(len(kept), 2)
        self.assertIn("C", kept)

    def test_below_threshold_pair_is_not_removed(self):
        pairs = [("A", "B", 0.01)]  # below KIN_MIN (0.0442)
        kept, removed = kir.greedy_unrelated(["A", "B"], pairs, kin_min=kir.KIN_MIN)
        self.assertEqual(removed, [])
        self.assertEqual(kept, {"A", "B"})

    def test_hub_person_removed_first(self):
        # D related to A, B, and C; nobody else related -- must remove D alone.
        pairs = [("D", "A", 0.2), ("D", "B", 0.2), ("D", "C", 0.2)]
        kept, removed = kir.greedy_unrelated(["A", "B", "C", "D"], pairs, kin_min=kir.KIN_MIN)
        self.assertEqual(removed, ["D"])
        self.assertEqual(kept, {"A", "B", "C"})


class TestParseHapGtf(unittest.TestCase):
    def test_novel_and_known_and_undetermined(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "hap1.gtf.gz")
            _write_gtf_gz(path, [
                ("KIR3DL3", "KIR3DL3*00101"),
                ("KIR2DS4", "KIR2DS4*00301:new"),  # HLA-style colon convention (also must work)
                ("KIR3DL2", "KIR3DL2*00201new"),   # REAL KIR convention -- no colon before "new"
                ("KIR2DL5A", "undetermined"),
            ])
            rows = kir.parse_hap_gtf(path)
            self.assertEqual(len(rows), 4)
            by_gene = {r["gene"]: r for r in rows}
            self.assertFalse(by_gene["KIR3DL3"]["is_novel"])
            self.assertFalse(by_gene["KIR3DL3"]["is_undetermined"])
            self.assertTrue(by_gene["KIR2DS4"]["is_novel"])
            self.assertTrue(by_gene["KIR3DL2"]["is_novel"])
            self.assertTrue(by_gene["KIR2DL5A"]["is_undetermined"])
            self.assertFalse(by_gene["KIR2DL5A"]["is_novel"])

    def test_missing_file_returns_empty(self):
        self.assertEqual(kir.parse_hap_gtf("/no/such/path.gtf.gz"), [])


class TestClassifyNovoltyTier(unittest.TestCase):
    def _row(self, consensus, cds_distance=None, cds_mut=None):
        return {"consensus": consensus, "cds_distance": cds_distance, "cds_mut": cds_mut,
                "is_novel": consensus.endswith("new") and consensus != "undetermined",
                "is_undetermined": consensus == "undetermined"}

    def test_known(self):
        row = self._row("KIR3DL3*00101")
        self.assertEqual(kir.classify_novelty_tier(row), "known")

    def test_undetermined(self):
        row = self._row("undetermined")
        self.assertEqual(kir.classify_novelty_tier(row), "undetermined")

    def test_novel_genomic_known_cds(self):
        row = self._row("KIR3DL2*00201new", cds_distance=0)
        self.assertEqual(kir.classify_novelty_tier(row), "novel_genomic_known_cds")

    def test_novel_cds_synonymous(self):
        # Real example from pilot output: Gly(GGG)<Gly(GGC) -- same amino acid both sides.
        row = self._row("KIR2DL2*004new", cds_distance=1,
                         cds_mut="KIR2DL2*004|:35*cg:1011|Gly(GGG)<Gly(GGC);")
        self.assertEqual(kir.classify_novelty_tier(row), "novel_cds_synonymous")

    def test_novel_protein(self):
        # Real example: Gln(CAA)<Rrg(CGA) -- different amino acid.
        row = self._row("KIR3DL2*00801new", cds_distance=1,
                         cds_mut="KIR3DL2*0080101|:121*ga:1246|Gln(CAA)<Rrg(CGA);")
        self.assertEqual(kir.classify_novelty_tier(row), "novel_protein")

    def test_novel_protein_wins_if_any_candidate_nonsynonymous(self):
        # Multiple tied candidates, comma-joined -- one synonymous, one not: must be novel_protein.
        row = self._row("KIR2DP1*001new", cds_distance=1,
                         cds_mut="A|:1*ac:1|Gly(GGG)<Gly(GGC),B|:2*ac:2|Gln(CAA)<Rrg(CGA);")
        self.assertEqual(kir.classify_novelty_tier(row), "novel_protein")

    def test_missing_cds_distance_is_unclassified(self):
        row = self._row("KIR2DL1*001new")
        self.assertEqual(kir.classify_novelty_tier(row), "novel_unclassified")


class TestClassifyFrameworkMiss(unittest.TestCase):
    def test_centromeric_end_missing_nothing_after(self):
        # KIR3DL3 missing, and none of 2DL4/3DP1/3DL2 called either -> edge fragmented.
        self.assertEqual(kir.classify_framework_miss("KIR3DL3", []), "edge_fragmented")

    def test_centromeric_end_missing_with_telomeric_genes_present(self):
        self.assertEqual(kir.classify_framework_miss("KIR3DL3", ["KIR3DL2"]), "flanked")

    def test_telomeric_end_missing_nothing_before(self):
        self.assertEqual(kir.classify_framework_miss("KIR3DL2", []), "edge_fragmented")

    def test_middle_gene_flanked_both_sides(self):
        self.assertEqual(kir.classify_framework_miss("KIR2DL4", ["KIR3DL3", "KIR3DL2"]), "flanked")

    def test_middle_gene_flanked_one_side_only(self):
        self.assertEqual(kir.classify_framework_miss("KIR3DP1", ["KIR3DL3"]),
                          "ambiguous_partial_flank")


class TestAggregateQuality(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.outroot = self.tmpdir.name

        class Args:
            outroot = self.outroot
        self.args = Args()

        # Person 1: cA content on both haps (framework genes + KIR2DL3 + KIR3DL1, no B genes).
        p1 = os.path.join(self.outroot, "1", "immuannot_output")
        os.makedirs(p1, exist_ok=True)
        cA_genes = [("KIR3DL3", "KIR3DL3*001"), ("KIR2DL4", "KIR2DL4*001"),
                    ("KIR3DP1", "KIR3DP1*001"), ("KIR3DL2", "KIR3DL2*001"),
                    ("KIR2DL3", "KIR2DL3*001"), ("KIR3DL1", "KIR3DL1*001")]
        _write_gtf_gz(os.path.join(p1, "hap1.gtf.gz"), cA_genes)
        _write_gtf_gz(os.path.join(p1, "hap2.gtf.gz"), cA_genes)

        # Person 2: cB content on hap1 (adds KIR2DS1, KIR3DS1 -- so 3DL1+3DS1 co-occur here),
        # hap2 has no KIR calls at all (simulates a failed/empty haplotype).
        p2 = os.path.join(self.outroot, "2", "immuannot_output")
        os.makedirs(p2, exist_ok=True)
        cB_genes = cA_genes + [("KIR2DS1", "KIR2DS1*001"), ("KIR3DS1", "KIR3DS1*001:new")]
        _write_gtf_gz(os.path.join(p2, "hap1.gtf.gz"), cB_genes)
        _write_gtf_gz(os.path.join(p2, "hap2.gtf.gz"), [])  # empty -> no KIR transcript rows

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_hap_success_rate_and_framework_presence(self):
        q = kir.aggregate_quality(["1", "2"], self.args)
        self.assertEqual(q["n_hap_attempted"], 4)
        self.assertEqual(q["n_hap_with_kir_calls"], 3)  # person 2's hap2 is empty
        self.assertAlmostEqual(q["hap_success_rate_pct"], 75.0)
        # Framework genes present on all 3 successful haps.
        for g in kir.FRAMEWORK_GENES:
            self.assertAlmostEqual(q["framework_gene_presence_pct"][g], 100.0)

    def test_ab_content_classification(self):
        q = kir.aggregate_quality(["1", "2"], self.args)
        # 2 cA haps (person 1) + 1 cB hap (person 2's hap1).
        self.assertEqual(q["hap_content_A_vs_B"].get("cA"), 2)
        self.assertEqual(q["hap_content_A_vs_B"].get("cB"), 1)

    def test_mutual_exclusivity_rates(self):
        q = kir.aggregate_quality(["1", "2"], self.args)
        # KIR2DL2/KIR2DL3: only KIR2DL3 ever appears (no KIR2DL2 anywhere) -> 0% co-occurrence.
        self.assertEqual(q["pct_haps_2dl2_and_2dl3_cooccur"], 0.0)
        # KIR3DL1/KIR3DS1: all 3 successful haps carry KIR3DL1; only person 2 hap1 also carries
        # KIR3DS1 -> 1/3 haps with either gene co-occur.
        self.assertAlmostEqual(q["pct_haps_3dl1_and_3ds1_cooccur"], round(100.0 / 3, 1), places=1)

    def test_novelty_fraction(self):
        q = kir.aggregate_quality(["1", "2"], self.args)
        # 6+6+8 = 20 total calls, exactly 1 tagged ":new" (KIR3DS1 on person 2 hap1).
        self.assertEqual(q["total_kir_calls"], 20)
        self.assertAlmostEqual(q["pct_novel"], 5.0)


if __name__ == "__main__":
    unittest.main()
