#!/usr/bin/env python3
"""Unit tests for 43_kir_full_aggregate.py, synthetic data only.

Covers what would be silently wrong and hard to notice on a real ~12k-person run:
  1. Disclosure primitives (suppressed()/rate_row()) -- the hard <20 masking rule, both for raw
     counts and for rates whose numerator OR denominator is small.
  2. kir_allele_freq.tsv's pooling rule -- alleles with count >= 20 kept individually, alleles
     with count < 20 pooled into one "other (<20 each)" row, and that pooled row itself follows
     the same disclosure rule if its summed count is still < 20.
  3. End-to-end small run (discover_people -> parse_all -> compute_bundle, on real gzip'd
     synthetic hap*.gtf.gz files, same fixture style as test_kir_pilot.py) feeding
     kir_gene_by_ancestry.tsv's construction, confirming ancestry grouping + presence/novelty
     rates + the small-cell blanking behave correctly together, not just in isolation.

Run: python3 scripts/hla_popgen/tests/test_kir_full_aggregate.py
"""
import gzip
import importlib.util
import os
import sys
import tempfile
import unittest
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)


def _load_module(filename, modname):
    spec = importlib.util.spec_from_file_location(modname, os.path.join(HLA_POPGEN_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


kir = _load_module("41_kir_pilot.py", "kir_pilot_test_target_43")
agg = _load_module("43_kir_full_aggregate.py", "kir_full_aggregate_test_target")
agg._kir = kir  # module-level global parse_person()/compute_bundle() callers expect


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


class TestDisclosurePrimitives(unittest.TestCase):
    def test_suppressed_masks_1_to_19(self):
        self.assertEqual(agg.suppressed(0), "0")  # genuine zero, not disclosive -- must not mask
        for n in (1, 5, 19):
            self.assertEqual(agg.suppressed(n), "<20")
        self.assertEqual(agg.suppressed(20), "20")
        self.assertEqual(agg.suppressed(1000), "1000")

    def test_rate_row_zero_numerator_is_not_blanked(self):
        # 0/50 is a genuine, non-disclosive zero rate -- must NOT be blanked like a 1-19 count.
        n, d, pct, lo, hi = agg.rate_row(0, 50)
        self.assertEqual(n, "0")
        self.assertEqual(d, "50")
        self.assertEqual(pct, "0.0")

    def test_rate_row_zero_denominator_is_blank(self):
        n, d, pct, lo, hi = agg.rate_row(0, 0)
        self.assertEqual(n, "0")
        self.assertEqual(d, "0")
        self.assertEqual(pct, "")

    def test_rate_row_blank_when_numerator_small(self):
        n, d, pct, lo, hi = agg.rate_row(5, 100)
        self.assertEqual(n, "<20")
        self.assertEqual(d, "100")
        self.assertEqual(pct, "")
        self.assertEqual(lo, "")
        self.assertEqual(hi, "")

    def test_rate_row_blank_when_denominator_small(self):
        n, d, pct, lo, hi = agg.rate_row(15, 18)
        self.assertEqual(n, "<20")
        self.assertEqual(d, "<20")
        self.assertEqual(pct, "")

    def test_rate_row_populated_when_both_clear_threshold(self):
        n, d, pct, lo, hi = agg.rate_row(25, 100)
        self.assertEqual(n, "25")
        self.assertEqual(d, "100")
        self.assertEqual(pct, "25.0")
        self.assertLess(float(lo), 25.0)
        self.assertGreater(float(hi), 25.0)

    def test_rate_row_exactly_20_is_not_censored(self):
        n, d, pct, lo, hi = agg.rate_row(20, 40)
        self.assertEqual(n, "20")
        self.assertEqual(pct, "50.0")

    def test_wilson_ci_matches_hand_computed_bounds(self):
        # p=0.5, n=100 -- Wilson interval should be roughly symmetric around 0.5, within (0.4,0.6).
        p, lo, hi = agg.wilson_ci(50, 100)
        self.assertAlmostEqual(p, 0.5)
        self.assertGreater(lo, 0.39)
        self.assertLess(hi, 0.61)

    def test_wilson_ci_zero_nobs_is_nan(self):
        p, lo, hi = agg.wilson_ci(0, 0)
        self.assertTrue(p != p)  # NaN != NaN


class TestAlleleFreqPooling(unittest.TestCase):
    def _bundle_with_known_alleles(self, allele_counts):
        """Builds a minimal bundle dict shaped like compute_bundle()'s output, with only the
        gene_known_allele_counts field populated (the only field build_allele_freq() reads)."""
        b = defaultdict(lambda: defaultdict(dict))
        b["gene_known_allele_counts"] = {"KIR3DL3": dict(allele_counts)}
        # build_allele_freq iterates kir.KIR_GENES -- give every other gene an empty dict so the
        # loop is well-defined (matches what compute_bundle() would actually produce).
        for g in kir.KIR_GENES:
            b["gene_known_allele_counts"].setdefault(g, {})
        return b

    def test_alleles_ge20_kept_individually(self):
        b = self._bundle_with_known_alleles({"A*001": 25, "A*002": 30})
        df = agg.build_allele_freq(b, kir)
        rows = df[df["gene"] == "KIR3DL3"]
        alleles = set(rows["allele"])
        self.assertIn("A*001", alleles)
        self.assertIn("A*002", alleles)
        self.assertNotIn("other (<20 each)", alleles)
        row_a1 = rows[rows["allele"] == "A*001"].iloc[0]
        self.assertEqual(row_a1["n_haplotypes"], "25")
        self.assertEqual(row_a1["n_known_calls_total_gene"], "55")
        self.assertEqual(row_a1["freq_pct"], str(round(100 * 25 / 55, 1)))

    def test_small_alleles_pooled_and_disclosed_if_sum_clears_20(self):
        # A: 25 (kept). B: 12, C: 9 -> pooled sum = 21 (>=20, so freq computed for the pooled row).
        b = self._bundle_with_known_alleles({"A*001": 25, "B*001": 12, "C*001": 9})
        df = agg.build_allele_freq(b, kir)
        rows = df[df["gene"] == "KIR3DL3"]
        alleles = set(rows["allele"])
        self.assertEqual(alleles, {"A*001", "other (<20 each)"})
        pooled = rows[rows["allele"] == "other (<20 each)"].iloc[0]
        self.assertEqual(pooled["n_haplotypes"], "21")
        self.assertNotEqual(pooled["freq_pct"], "")

    def test_pooled_row_itself_censored_if_sum_still_below_20(self):
        # A: 25 (kept). B: 8, C: 6 -> pooled sum = 14 (<20) -> pooled row must read "<20", blank pct.
        b = self._bundle_with_known_alleles({"A*001": 25, "B*001": 8, "C*001": 6})
        df = agg.build_allele_freq(b, kir)
        rows = df[df["gene"] == "KIR3DL3"]
        pooled = rows[rows["allele"] == "other (<20 each)"].iloc[0]
        self.assertEqual(pooled["n_haplotypes"], "<20")
        self.assertEqual(pooled["freq_pct"], "")

    def test_gene_with_no_known_alleles_produces_no_rows(self):
        b = self._bundle_with_known_alleles({})
        df = agg.build_allele_freq(b, kir)
        self.assertEqual(len(df), 0)


class TestEndToEndAncestryTable(unittest.TestCase):
    """Builds a small synthetic --outroot (real gzip'd gtf.gz files) for people spread across two
    ancestry groups of very different size, runs the real discover_people/parse_all/compute_bundle
    path, and checks kir_gene_by_ancestry.tsv's construction: the large group shows real numbers,
    the small group is blanked, and per-gene presence is counted correctly."""

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.outroot = os.path.join(self.tmpdir.name, "outroot")
        os.makedirs(self.outroot, exist_ok=True)

        # AFR: 12 people (24 haps), all carrying KIR3DL3 (known allele) + KIR2DL4 (framework).
        self.afr_ids = [f"afr{i}" for i in range(12)]
        for pid in self.afr_ids:
            genes = [("KIR3DL3", "KIR3DL3*001"), ("KIR2DL4", "KIR2DL4*001")]
            for hap in ("hap1", "hap2"):
                _write_gtf_gz(os.path.join(self.outroot, pid, "immuannot_output", f"{hap}.gtf.gz"),
                             genes)

        # EAS: 3 people (6 haps) -- deliberately below the 20-haplotype disclosure floor.
        self.eas_ids = [f"eas{i}" for i in range(3)]
        for pid in self.eas_ids:
            genes = [("KIR3DL3", "KIR3DL3*001")]
            for hap in ("hap1", "hap2"):
                _write_gtf_gz(os.path.join(self.outroot, pid, "immuannot_output", f"{hap}.gtf.gz"),
                             genes)

        self.ancestry_of = {pid: "AFR" for pid in self.afr_ids}
        self.ancestry_of.update({pid: "EAS" for pid in self.eas_ids})

    def tearDown(self):
        self.tmpdir.cleanup()

    def _bundles_by_ancestry(self):
        pids = agg.discover_people(self.outroot)
        self.assertEqual(len(pids), 15)
        results = agg.parse_all(pids, self.outroot, workers=1, kir_script_path=None)
        results_by_pid = {r["person_id"]: r for r in results}
        by_anc = defaultdict(list)
        for p in pids:
            by_anc[self.ancestry_of[p]].append(p)
        return {a: agg.compute_bundle(ps, results_by_pid, kir) for a, ps in by_anc.items()}

    def test_afr_group_has_ge20_haps_and_full_presence(self):
        bundles = self._bundles_by_ancestry()
        self.assertEqual(bundles["AFR"]["n_hap_with_calls"], 24)
        self.assertEqual(bundles["AFR"]["gene_hap_presence"]["KIR3DL3"], 24)

    def test_gene_by_ancestry_table_blanks_the_small_group(self):
        bundles = self._bundles_by_ancestry()
        df = agg.build_gene_by_ancestry(bundles, kir)
        row_afr = df[(df["gene"] == "KIR3DL3") & (df["ancestry"] == "AFR")].iloc[0]
        self.assertEqual(row_afr["n_haplotypes_carrying"], "24")
        self.assertEqual(row_afr["presence_pct"], "100.0")

        row_eas = df[(df["gene"] == "KIR3DL3") & (df["ancestry"] == "EAS")].iloc[0]
        self.assertEqual(row_eas["n_haplotypes_carrying"], "<20")
        self.assertEqual(row_eas["n_haplotypes_total"], "<20")
        self.assertEqual(row_eas["presence_pct"], "")

    def test_gene_absent_in_all_haplotypes_still_has_a_row_with_zero_presence(self):
        bundles = self._bundles_by_ancestry()
        df = agg.build_gene_by_ancestry(bundles, kir)
        row = df[(df["gene"] == "KIR3DL1") & (df["ancestry"] == "AFR")].iloc[0]
        self.assertEqual(row["n_haplotypes_carrying"], "0")


class TestDiscoverPeople(unittest.TestCase):
    def test_limit_is_deterministic_prefix(self):
        with tempfile.TemporaryDirectory() as d:
            for pid in ("c", "a", "b"):
                os.makedirs(os.path.join(d, pid))
            all_pids = agg.discover_people(d)
            self.assertEqual(all_pids, ["a", "b", "c"])
            self.assertEqual(agg.discover_people(d, limit=2), ["a", "b"])

    def test_missing_outroot_returns_empty(self):
        self.assertEqual(agg.discover_people("/no/such/outroot"), [])


class TestParsePersonMissingOrCorrupt(unittest.TestCase):
    def test_missing_gtf_files_count_as_invalid_not_crash(self):
        with tempfile.TemporaryDirectory() as d:
            r = agg.parse_person("ghost", d)
            self.assertEqual(r["haps"], [])
            self.assertEqual(r["n_invalid_haps"], 2)

    def test_corrupt_gzip_counts_as_invalid(self):
        with tempfile.TemporaryDirectory() as d:
            bad_path = os.path.join(d, "p1", "immuannot_output", "hap1.gtf.gz")
            os.makedirs(os.path.dirname(bad_path), exist_ok=True)
            with open(bad_path, "wb") as f:
                f.write(b"not actually gzip data")
            good_path = os.path.join(d, "p1", "immuannot_output", "hap2.gtf.gz")
            _write_gtf_gz(good_path, [("KIR3DL3", "KIR3DL3*001")])
            r = agg.parse_person("p1", d)
            self.assertEqual(len(r["haps"]), 1)
            self.assertEqual(r["n_invalid_haps"], 1)


if __name__ == "__main__":
    unittest.main()
