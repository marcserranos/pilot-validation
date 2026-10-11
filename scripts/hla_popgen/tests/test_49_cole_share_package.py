#!/usr/bin/env python3
"""Unit tests for 49_cole_share_package.py (S04 WS-E), synthetic data only -- no VM, no AoU data.

Covers:
  1. hla_calls: schema, novelty_tier mapping, already_shared flag, two-field truncation.
  2. kir_calls: synthetic hap*.gtf.gz -> known/novel rows, copy_index, novelty_tier vocabulary.
  3. persons.tsv: hla_called/kir_called flags and exclusion_reason (incl. the sequel2 KIR gap).
  4. MANIFEST.tsv: md5/sha256 correctness, and --verify catching a corrupted file.
  5. --dry-run: prints only aggregate counts (no person_id/allele strings) and writes nothing.
  6. Guard rail: refuses to write under ~/pipeline_outputs* / ~/pipeline_outputs_kir.

Run:
    cd scripts/hla_popgen && python3 tests/test_49_cole_share_package.py
    (or) python3 -m pytest scripts/hla_popgen/tests/test_49_cole_share_package.py -q
"""
import contextlib
import gzip
import importlib.util
import io
import os
import sys
import tempfile
import unittest

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)
if HLA_POPGEN_DIR not in sys.path:
    sys.path.insert(0, HLA_POPGEN_DIR)


def _load(filename, modname):
    spec = importlib.util.spec_from_file_location(modname, os.path.join(HLA_POPGEN_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


pkg = _load("49_cole_share_package.py", "cole_share_package_test_target")

TABLE1_COLS = ["person_id", "hap", "contig", "gene", "copy_index", "gene_class", "consensus",
               "is_novel", "novelty_depth", "novelty_class", "template_warning", "cds_mut",
               "alleles", "template_allele", "cds_distance", "n_aa_changes"]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
def make_table1_tsv(path):
    """Two people, one known HLA-A call, one novel protein-altering DQB1 call (person p2), and
    one DQA1 call for p2 on the same hap/contig (to exercise the DQ/DP phasing-confidence join)."""
    rows = [
        # person_id, hap, contig, gene, copy_index, gene_class, consensus, is_novel,
        # novelty_depth, novelty_class, template_warning, cds_mut, alleles, template_allele,
        # cds_distance, n_aa_changes
        ["p1", "hap1", "ctgA", "HLA-A", 1, "classical_I", "HLA-A*01:01:01:01", False, "", "",
         "NA", "", "", "", "", ""],
        ["p2", "hap1", "ctgB", "HLA-DQB1", 1, "classical_II", "HLA-DQB1*02:01:01:new", True, 2,
         "protein_altering", "NA", "Ref|cs|Gly(GGG)<Ser(AGC)", "", "", "0", "1"],
        ["p2", "hap1", "ctgB", "HLA-DQA1", 1, "classical_II", "HLA-DQA1*01:01:01:01", False, "",
         "", "NA", "", "", "", "", ""],
    ]
    df = pd.DataFrame(rows, columns=TABLE1_COLS)
    df.to_csv(path, sep="\t", index=False)


def _gtf_transcript_line(gene, consensus, gene_id_suffix="", contig="chr19_synth",
                          cds_distance=None, cds_mut=None):
    attrs = (f'gene_id "IG{gene}{gene_id_suffix}"; transcript_id "IAT{gene}.1"; '
             f'gene_name "{gene}"; consensus "{consensus}"; alleles "{consensus}";')
    if cds_distance is not None:
        attrs += f' cds_distance {cds_distance};'
    if cds_mut is not None:
        attrs += f' cds_mut "{cds_mut}";'
    return f"{contig}\tImmuannot\ttranscript\t100\t200\t.\t+\t.\t{attrs}\n"


def write_gtf_gz(path, lines):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as f:
        f.write("##synthetic\n")
        for line in lines:
            f.write(line)


def make_kir_outroot(root):
    """p1: hap1 known KIR2DL1, hap2 nothing. p2 (sequel2-style gap simulation handled at the
    persons level, not here): hap1 known + novel_protein KIR3DL1, hap2 known KIR2DL1."""
    write_gtf_gz(os.path.join(root, "p1", "immuannot_output", "hap1.gtf.gz"),
                 [_gtf_transcript_line("KIR2DL1", "KIR2DL1*00101")])
    write_gtf_gz(os.path.join(root, "p1", "immuannot_output", "hap2.gtf.gz"), [])
    write_gtf_gz(os.path.join(root, "p2", "immuannot_output", "hap1.gtf.gz"),
                 [_gtf_transcript_line("KIR3DL1", "KIR3DL1*00101"),
                  _gtf_transcript_line("KIR2DL1", "KIR2DL1*00201new", cds_distance=1,
                                        cds_mut="Ref|cs|Gly(GGG)<Ser(AGC)")])
    write_gtf_gz(os.path.join(root, "p2", "immuannot_output", "hap2.gtf.gz"),
                 [_gtf_transcript_line("KIR2DL1", "KIR2DL1*00101")])
    # p3 has a directory but both haps are empty/unparseable -> kir_called False, not a
    # "dir_missing" case.
    write_gtf_gz(os.path.join(root, "p3", "immuannot_output", "hap1.gtf.gz"), [])
    write_gtf_gz(os.path.join(root, "p3", "immuannot_output", "hap2.gtf.gz"), [])


def make_cohort_membership_tsv(path):
    rows = [
        {"person_id": "p1", "ancestry_pred": "eur", "platform": "revio",
         "p_afr": 0.01, "p_amr": 0.01, "p_eas": 0.01, "p_eur": 0.95, "p_mid": 0.01, "p_sas": 0.01},
        {"person_id": "p2", "ancestry_pred": "afr", "platform": "revio",
         "p_afr": 0.9, "p_amr": 0.02, "p_eas": 0.02, "p_eur": 0.02, "p_mid": 0.02, "p_sas": 0.02},
        {"person_id": "p3", "ancestry_pred": "eas", "platform": "sequel2",
         "p_afr": 0.02, "p_amr": 0.02, "p_eas": 0.9, "p_eur": 0.02, "p_mid": 0.02, "p_sas": 0.02},
        {"person_id": "p4", "ancestry_pred": "sas", "platform": "revio",
         "p_afr": 0.02, "p_amr": 0.02, "p_eas": 0.02, "p_eur": 0.02, "p_mid": 0.02, "p_sas": 0.9},
    ]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


# ---------------------------------------------------------------------------
# 1. hla_calls
# ---------------------------------------------------------------------------
class TestHlaCalls(unittest.TestCase):
    def test_schema_and_novelty_and_already_shared(self):
        with tempfile.TemporaryDirectory() as d:
            t1_path = os.path.join(d, "hla_calls_rich.tsv")
            make_table1_tsv(t1_path)
            df, stats = pkg.build_hla_calls(t1_path, phasing_confidence_path=None)

            expected_cols = {"person_id", "hap", "contig", "gene", "gene_class", "copy_index",
                              "allele", "two_field", "novelty_tier", "artifact_label",
                              "dq_dp_phasing_confidence", "already_shared_two_field_common"}
            self.assertEqual(expected_cols, set(df.columns))
            self.assertEqual(len(df), 3)

            known_row = df[(df.person_id == "p1") & (df.gene == "HLA-A")].iloc[0]
            self.assertEqual(known_row["novelty_tier"], "known")
            self.assertEqual(known_row["two_field"], "HLA-A*01:01")
            self.assertTrue(bool(known_row["already_shared_two_field_common"]))

            novel_row = df[(df.person_id == "p2") & (df.gene == "HLA-DQB1")].iloc[0]
            self.assertEqual(novel_row["novelty_tier"], "protein")
            self.assertFalse(bool(novel_row["already_shared_two_field_common"]))

            self.assertFalse(stats["phasing_confidence_available"])

    def test_phasing_confidence_join_dq_dp_only(self):
        with tempfile.TemporaryDirectory() as d:
            t1_path = os.path.join(d, "hla_calls_rich.tsv")
            make_table1_tsv(t1_path)
            phase_path = os.path.join(d, "phasing_confidence_per_person.tsv")
            pd.DataFrame([
                {"person_id": "p2", "hap": "hap1", "gene_a": "DQA1", "gene_b": "DQB1",
                 "status": "physical", "contig": "ctgB"},
            ]).to_csv(phase_path, sep="\t", index=False)

            df, stats = pkg.build_hla_calls(t1_path, phasing_confidence_path=phase_path)
            self.assertTrue(stats["phasing_confidence_available"])
            dqb1 = df[(df.person_id == "p2") & (df.gene == "HLA-DQB1")].iloc[0]
            self.assertIn("DQA1~DQB1:physical", dqb1["dq_dp_phasing_confidence"])
            # HLA-A (not a DQ/DP gene) must not get a phasing flag.
            a_row = df[(df.person_id == "p1") & (df.gene == "HLA-A")].iloc[0]
            self.assertTrue(pd.isna(a_row["dq_dp_phasing_confidence"]))

    def test_missing_phasing_file_does_not_crash(self):
        with tempfile.TemporaryDirectory() as d:
            t1_path = os.path.join(d, "hla_calls_rich.tsv")
            make_table1_tsv(t1_path)
            df, stats = pkg.build_hla_calls(
                t1_path, phasing_confidence_path=os.path.join(d, "does_not_exist.tsv"))
            self.assertFalse(stats["phasing_confidence_available"])
            self.assertTrue(df["dq_dp_phasing_confidence"].isna().all())


# ---------------------------------------------------------------------------
# 2. kir_calls
# ---------------------------------------------------------------------------
class TestKirCalls(unittest.TestCase):
    def test_row_counts_and_novelty_tiers(self):
        with tempfile.TemporaryDirectory() as d:
            kir_root = os.path.join(d, "pipeline_outputs_kir")
            make_kir_outroot(kir_root)
            df, stats, kir_call_status = pkg.build_kir_calls(kir_root)

            self.assertEqual(set(df.columns), {
                "person_id", "hap", "contig", "gene", "copy_index", "allele", "novelty_tier",
                "cds_sha8", "protein_sha8", "cds_available"})
            # p1: 1 call (hap1 only). p2: 3 calls (hap1 x2 + hap2 x1). p3: 0 calls.
            self.assertEqual(len(df), 4)
            self.assertEqual(stats["n_persons_attempted"], 3)
            self.assertEqual(stats["n_persons_with_calls"], 2)
            self.assertFalse(kir_call_status["p3"])
            self.assertTrue(kir_call_status["p1"])

            novel = df[(df.person_id == "p2") & (df.gene == "KIR2DL1") & (df.hap == "hap1")]
            self.assertEqual(len(novel), 1)
            self.assertEqual(novel.iloc[0]["novelty_tier"], "novel_protein")

            known = df[(df.person_id == "p1") & (df.gene == "KIR2DL1")].iloc[0]
            self.assertEqual(known["novelty_tier"], "known")
            # No cds.fa.gz was written in this fixture -> cds_available False, hashes None.
            self.assertFalse(bool(known["cds_available"]))
            self.assertIsNone(known["cds_sha8"])
            self.assertFalse(stats["cds_join_available"])


# ---------------------------------------------------------------------------
# 3. persons.tsv
# ---------------------------------------------------------------------------
class TestPersons(unittest.TestCase):
    def test_hla_kir_flags_and_exclusion_reason(self):
        with tempfile.TemporaryDirectory() as d:
            cohort_path = os.path.join(d, "cohort_membership.tsv")
            make_cohort_membership_tsv(cohort_path)
            kir_call_status = {"p1": True, "p2": True, "p3": False}
            hla_person_ids = {"p1", "p2"}  # p3, p4 have no HLA rows in this fixture

            persons = pkg.build_persons(cohort_path, relatedness_path=None, kir_call_status=kir_call_status,
                                        hla_person_ids=hla_person_ids, skip_relatedness=True)

            self.assertEqual(set(persons.columns), {
                "person_id", "ancestry_pred", "platform", "p_afr", "p_amr", "p_eas", "p_eur",
                "p_mid", "p_sas", "unrelated", "hla_called", "kir_called", "exclusion_reason"})

            p1 = persons.set_index("person_id").loc["p1"]
            self.assertTrue(p1["hla_called"])
            self.assertTrue(p1["kir_called"])
            self.assertEqual(p1["exclusion_reason"], "")

            p3 = persons.set_index("person_id").loc["p3"]  # platform sequel2, no KIR calls
            self.assertFalse(p3["kir_called"])
            self.assertIn("kir_excluded_sequel2_tier3_self_align", p3["exclusion_reason"])
            self.assertIn("hla_not_called", p3["exclusion_reason"])

            p4 = persons.set_index("person_id").loc["p4"]  # not in kir_call_status at all
            self.assertIn("kir_dir_missing", p4["exclusion_reason"])


# ---------------------------------------------------------------------------
# 4. Manifest
# ---------------------------------------------------------------------------
class TestManifest(unittest.TestCase):
    def test_manifest_md5_sha256_correct_and_verify_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            f1 = os.path.join(d, "a.tsv")
            f2 = os.path.join(d, "b.tsv")
            with open(f1, "w") as f:
                f.write("col1\tcol2\nx\ty\n")
            with open(f2, "w") as f:
                f.write("colA\ncontent\n")
            manifest_path, manifest = pkg.write_manifest(d, [("a.tsv", 1), ("b.tsv", 1)])
            self.assertTrue(os.path.exists(manifest_path))

            import hashlib
            with open(f1, "rb") as fh:
                expected_md5 = hashlib.md5(fh.read()).hexdigest()
            row = manifest[manifest.file == "a.tsv"].iloc[0]
            self.assertEqual(row["md5"], expected_md5)
            self.assertEqual(row["bytes"], os.path.getsize(f1))

            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                ok = pkg.verify_manifest(d)
            self.assertTrue(ok)
            self.assertNotIn("content", buf.getvalue())  # file contents never printed

    def test_verify_detects_corruption(self):
        with tempfile.TemporaryDirectory() as d:
            f1 = os.path.join(d, "a.tsv")
            with open(f1, "w") as f:
                f.write("col1\ncontent\n")
            pkg.write_manifest(d, [("a.tsv", 1)])
            with open(f1, "a") as f:
                f.write("tampered\n")
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                ok = pkg.verify_manifest(d)
            self.assertFalse(ok)
            self.assertIn("MISMATCH", buf.getvalue())


# ---------------------------------------------------------------------------
# 5. --dry-run and guard rail, end to end via main()
# ---------------------------------------------------------------------------
class TestMainDryRunAndGuardRail(unittest.TestCase):
    def _build_inputs(self, d):
        t1_path = os.path.join(d, "hla_calls_rich.tsv")
        make_table1_tsv(t1_path)
        kir_root = os.path.join(d, "pipeline_outputs_kir")
        make_kir_outroot(kir_root)
        cohort_path = os.path.join(d, "cohort_membership.tsv")
        make_cohort_membership_tsv(cohort_path)
        return t1_path, kir_root, cohort_path

    def test_dry_run_writes_nothing_and_prints_no_ids(self):
        with tempfile.TemporaryDirectory() as d:
            t1_path, kir_root, cohort_path = self._build_inputs(d)
            out_dir = os.path.join(d, "s04_out")
            argv = ["--hla-table1", t1_path, "--kir-outroot", kir_root,
                    "--cohort-membership", cohort_path, "--out-dir", out_dir,
                    "--skip-relatedness", "--dry-run"]
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                pkg.main(argv)
            output = buf.getvalue()
            self.assertFalse(os.path.exists(out_dir))
            for pid in ("p1", "p2", "p3", "p4"):
                self.assertNotIn(pid, output)
            for allele in ("HLA-A*01:01:01:01", "KIR2DL1*00101", "KIR3DL1*00101"):
                self.assertNotIn(allele, output)
            self.assertIn("DRY RUN", output)
            self.assertIn("hla_calls", output)
            self.assertIn("kir_calls", output)

    def test_build_then_verify_end_to_end(self):
        with tempfile.TemporaryDirectory() as d:
            t1_path, kir_root, cohort_path = self._build_inputs(d)
            out_dir = os.path.join(d, "s04", "share_release_test")
            argv = ["--hla-table1", t1_path, "--kir-outroot", kir_root,
                    "--cohort-membership", cohort_path, "--out-dir", out_dir,
                    "--skip-relatedness"]
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                pkg.main(argv)
            for fname in ("hla_calls.tsv.gz", "kir_calls.tsv.gz", "persons.tsv", "README.md",
                          "SCHEMA.md", "MANIFEST.tsv"):
                self.assertTrue(os.path.exists(os.path.join(out_dir, fname)), fname)

            argv_verify = ["--out-dir", out_dir, "--verify"]
            buf2 = io.StringIO()
            with contextlib.redirect_stdout(buf2), self.assertRaises(SystemExit) as cm:
                pkg.main(argv_verify)
            self.assertEqual(cm.exception.code, 0)
            self.assertIn("all files match manifest", buf2.getvalue())

    def test_refuses_to_write_under_pipeline_outputs(self):
        with tempfile.TemporaryDirectory() as d:
            t1_path, kir_root, cohort_path = self._build_inputs(d)
            bad_out = os.path.expanduser("~/pipeline_outputs/should_not_write_here")
            argv = ["--hla-table1", t1_path, "--kir-outroot", kir_root,
                    "--cohort-membership", cohort_path, "--out-dir", bad_out,
                    "--skip-relatedness"]
            with self.assertRaises(SystemExit):
                pkg.main(argv)
            self.assertFalse(os.path.exists(bad_out))


if __name__ == "__main__":
    unittest.main()
