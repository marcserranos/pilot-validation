#!/usr/bin/env python3
"""Tests for 62_cole_haplotype_units.py -- synthetic Table 1 only, no AoU data.

Covers: the pair-rule evidence levels (physical / flagged / not allowed / split / incomplete),
configurable rules, DQ G1/G2 classes, the DQ spacing and multi-copy anomaly flags, the coverage
flag (pass / low / not measured / file absent), the unrelated flag, PUBLIC masking and
complementary suppression, no person ids in PUBLIC outputs, INTERNAL refusing repo and reports
paths, and the committed dictionary matching the script.

Run: python3 -m pytest scripts/hla_popgen/tests/test_62_cole_haplotype_units.py -q
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
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
spec = importlib.util.spec_from_file_location("hu62", os.path.join(SRC, "62_cole_haplotype_units.py"))
hu = importlib.util.module_from_spec(spec)
sys.modules["hu62"] = hu
spec.loader.exec_module(hu)

RULES = hu.parse_pair_rules(hu.DEFAULT_PAIR_RULES)


def t1_rows(pid, hap, calls):
    """calls: list of (gene, consensus, contig, start, copy)."""
    return [{"person_id": pid, "hap": hap, "contig": c, "gene": f"HLA-{g}", "copy_index": cp,
             "consensus": f"HLA-{cons}", "novelty_depth": np.nan, "template_distance": 0,
             "gene_start": s, "gene_end": s + 5000, "cds_mut": np.nan, "template_warning": np.nan}
            for g, cons, c, s, cp in calls]


def make_inputs(tmp, n=400, seed=3):
    """n people x 2 haps. 70% G1 haplotype DRB1*03:01~DQA1*05:01~DQB1*02:01, 30% G2
    DRB1*15:01~DQA1*01:02~DQB1*06:02; DQ spacing 20 kb. Planted: person 0 hap1 split contigs,
    person 1 hap1 DQB1 missing, person 2 hap1 spacing 150 kb, person 3 hap1 second DQB1 copy,
    person 4 hap1 cross-group pair on one contig. Ancestry EUR/AFR + 5 MID."""
    rng = np.random.default_rng(seed)
    rows, coh = [], []
    for i in range(n):
        pid = str(7_000_000 + i)
        anc = "MID" if i >= n - 5 else ("EUR" if i % 2 else "AFR")
        coh.append({"person_id": pid, "ancestry_pred": anc.lower(),
                    "platform": "revio" if i % 3 else "sequel2e", "in_lr": True})
        for hap in ("hap1", "hap2"):
            g1 = rng.random() < 0.7
            drb, dqa, dqb = (("DRB1*03:01:01:01", "DQA1*05:01:01:02", "DQB1*02:01:01")
                             if g1 else ("DRB1*15:01:01:01", "DQA1*01:02:01:01", "DQB1*06:02:01"))
            ctg = f"{hap}_c1:1-900000"
            calls = [("DRB1", drb, ctg, 100_000, 1), ("DQA1", dqa, ctg, 200_000, 1),
                     ("DQB1", dqb, ctg, 220_000, 1), ("A", "A*02:01:01:01", f"{hap}_c0", 10, 1)]
            if hap == "hap1":
                if i == 0:
                    calls[2] = ("DQB1", dqb, f"{hap}_c2:1-50000", 1_000, 1)
                elif i == 1:
                    calls = [c for c in calls if c[0] != "DQB1"]
                elif i == 2:
                    calls[2] = ("DQB1", dqb, ctg, 350_000, 1)
                elif i == 3:
                    calls.append(("DQB1", "DQB1*06:02:01", ctg, 260_000, 2))
                elif i == 4:
                    calls[1] = ("DQA1", "DQA1*01:02:01:01", ctg, 200_000, 1)
                    calls[2] = ("DQB1", "DQB1*02:01:01", ctg, 220_000, 1)
            rows += t1_rows(pid, hap, calls)
    t1 = os.path.join(tmp, "t1.tsv")
    pd.DataFrame(rows).to_csv(t1, sep="\t", index=False)
    cohp = os.path.join(tmp, "cohort.tsv")
    pd.DataFrame(coh).to_csv(cohp, sep="\t", index=False)
    rel = os.path.join(tmp, "rel.tsv")
    pd.DataFrame({"i.s": ["7000010"], "j.s": ["7000011"], "kin": [0.25]}).to_csv(
        rel, sep="\t", index=False)
    depth = os.path.join(tmp, "per_person_depth.tsv")
    with open(depth, "w") as fh:
        fh.write("# INTERNAL — contains n<20 cells, do not export\n")
        fh.write("person_id\tregion\tmapq\tmean_depth\tn_reads\tcontrol_depth\trel_depth\n")
        for i in range(n - 50):  # last 50 people not measured
            rel_d = 0.3 if i in (5, 6) else 0.9
            fh.write(f"{7_000_000 + i}\twindow_MHC\t0\t10\t100\t12\t{rel_d}\n")
            fh.write(f"{7_000_000 + i}\twindow_MHC\t20\t9\t90\t11\t0.1\n")
    return t1, cohp, rel, depth


class TestRules(unittest.TestCase):
    def test_qc_summary_cross_field_suppression(self):
        """AMR-like stratum: G1 + G2 shown, predicted_incompatible 1-19, not_called equal to the
        unit's incomplete count (shown in evidence_level). Without cross-field suppression the
        incompatible count is of_haps - G1 - G2 - incomplete."""
        rows = []
        spec = [("physical", "G1", 600), ("physical", "G2", 300), ("physical", "predicted_incompatible", 7),
                ("incomplete", "not_called", 40), ("split_contigs", "G1", 25)]
        for ev, dq, k in spec:
            for u, _ in hu.UNITS:
                rows += [dict(unit=u, evidence_level=ev, dq_g1g2=dq, dq_structure="canonical",
                              coverage_flag="pass", platform="revio", unrelated="True",
                              ancestry_pred="AMR")] * k
        q = hu.public_qc_summary(pd.DataFrame(rows))
        q = q[(q["unit"] == "DQA1~DQB1") & (q["ancestry"] == "AMR")]
        cell = {(r.field, r.value): r.n_haps for r in q.itertuples()}
        self.assertEqual(cell[("dq_g1g2", "predicted_incompatible")], "<20")
        hidden = [v for (f, v), d in cell.items() if f == "dq_g1g2" and not str(d).isdigit()]
        self.assertGreaterEqual(len(hidden), 2)
        # not_called and incomplete are the same haplotypes: both shown or both hidden
        self.assertEqual(str(cell[("dq_g1g2", "not_called")]).isdigit(),
                         str(cell[("evidence_level", "incomplete")]).isdigit())

    def test_evidence_levels(self):
        g = ("DRB1", "DQA1", "DQB1")
        self.assertEqual(hu.evidence_for_unit(g, {"DRB1": "c", "DQA1": "c", "DQB1": "c"}, RULES)[0],
                         "physical")
        self.assertEqual(hu.evidence_for_unit(g, {"DRB1": "c", "DQA1": "c", "DQB1": "d"}, RULES)[0],
                         "split_contigs")
        self.assertEqual(hu.evidence_for_unit(g, {"DRB1": "c", "DQA1": "c"}, RULES)[0], "incomplete")
        flagged = hu.parse_pair_rules("DQA1~DQB1=flagged,DRB1~DQA1=allowed,DRB1~DQB1=allowed")
        self.assertEqual(hu.evidence_for_unit(("DQA1", "DQB1"), {"DQA1": "c", "DQB1": "c"},
                                              flagged)[0], "physical_flagged")
        self.assertEqual(hu.evidence_for_unit(("A", "B"), {"A": "c", "B": "c"}, RULES)[0],
                         "not_allowed_by_rule")
        self.assertEqual(hu.evidence_for_unit(("B", "C"), {"B": "c", "C": "c"}, RULES)[0], "physical")
        self.assertEqual(hu.pair_status("DPB1", "DQB1", RULES), "flagged")  # order-free
        self.assertEqual(hu.pair_status("X", "Y", RULES), "not_allowed")
        with self.assertRaises(ValueError):
            hu.parse_pair_rules("A~B=maybe")

    def test_dq(self):
        self.assertEqual(hu.classify_dq("DQA1*05:01", "DQB1*02:01"), "G1")
        self.assertEqual(hu.classify_dq("DQA1*01:02", "DQB1*06:02"), "G2")
        self.assertEqual(hu.classify_dq("DQA1*01:02", "DQB1*02:01"), "predicted_incompatible")
        self.assertEqual(hu.classify_dq(None, "DQB1*02:01"), "not_called")
        a = {"contig": "c", "start": 0}
        self.assertEqual(hu.dq_structure(a, {"contig": "c", "start": 20_000}, False, 100)[0],
                         "canonical")
        self.assertEqual(hu.dq_structure(a, {"contig": "c", "start": 150_000}, False, 100)[0],
                         "anomalous_spacing")
        self.assertEqual(hu.dq_structure(a, {"contig": "c", "start": 150_000}, True, 100)[0],
                         "anomalous_spacing;anomalous_multicopy")
        self.assertEqual(hu.dq_structure(a, {"contig": "d", "start": 1}, False, 100)[0],
                         "not_assessable")

    def test_coverage_flag(self):
        self.assertEqual(hu.coverage_flag(0.3, 0.5), "low_coverage")
        self.assertEqual(hu.coverage_flag(0.5, 0.5), "pass")
        self.assertEqual(hu.coverage_flag(None, 0.5), "not_measured")

    def test_normalise_depth_per_platform(self):
        depth = {f"r{i}": 0.8 for i in range(30)}
        depth.update({f"s{i}": 0.6 for i in range(30)})
        depth.update({"x1": 0.7, "x2": 0.3})  # tiny platform -> all-platform median
        plat = {**{f"r{i}": "revio" for i in range(30)}, **{f"s{i}": "sequel2e" for i in range(30)},
                "x1": "odd", "x2": "odd"}
        n = hu.normalise_depth(depth, plat)
        self.assertAlmostEqual(n["r0"], 1.0)
        self.assertAlmostEqual(n["s0"], 1.0)
        allmed = float(np.median(list(depth.values())))
        self.assertAlmostEqual(n["x2"], 0.3 / allmed)


class TestEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.t1, cls.coh, cls.rel, cls.depth = make_inputs(cls.tmp)
        cls.internal = os.path.join(cls.tmp, "internal")
        cls.out = os.path.join(cls.tmp, "public")
        hu.main(["--table1", cls.t1, "--cohort-membership", cls.coh, "--relatedness-table",
                 cls.rel, "--depth-table", cls.depth, "--internal-dir", cls.internal,
                 "--out-dir", cls.out])
        path = os.path.join(cls.internal, "haplotype_units.tsv.gz")
        with gzip.open(path, "rt") as fh:
            cls.first_line = fh.readline()
        cls.u = pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)

    def row(self, i, hap="hap1", unit="DRB1~DQA1~DQB1"):
        u = self.u
        return u[(u["person_id"] == str(7_000_000 + i)) & (u["hap"] == hap) & (u["unit"] == unit)].iloc[0]

    def test_internal_header_and_grain(self):
        self.assertTrue(self.first_line.startswith("# INTERNAL"))
        self.assertEqual(len(self.u), 400 * 2 * 2)
        self.assertFalse(self.u.duplicated(["person_id", "hap", "unit"]).any())
        self.assertTrue(os.path.exists(os.path.join(self.internal, "MANIFEST.tsv")))

    def test_planted_cases(self):
        self.assertEqual(self.row(0)["evidence_level"], "split_contigs")
        self.assertEqual(self.row(0, unit="DQA1~DQB1")["dq_structure"], "not_assessable")
        self.assertEqual(self.row(1)["evidence_level"], "incomplete")
        self.assertEqual(self.row(1)["dq_g1g2"], "not_called")
        self.assertEqual(self.row(2)["dq_structure"], "anomalous_spacing")
        self.assertAlmostEqual(float(self.row(2)["dq_spacing_kb"]), 150.0)
        self.assertEqual(self.row(3)["dq_structure"], "anomalous_multicopy")
        self.assertEqual(self.row(3)["allele_DQB1"].count("/"), 1)
        self.assertEqual(self.row(4)["dq_g1g2"], "predicted_incompatible")
        r = self.row(7)
        self.assertEqual(r["evidence_level"], "physical")
        self.assertEqual(r["dq_structure"], "canonical")
        self.assertIn(r["unit_two_field"], ("DRB1*03:01~DQA1*05:01~DQB1*02:01",
                                            "DRB1*15:01~DQA1*01:02~DQB1*06:02"))
        self.assertEqual(r["allele_DRB1"].count(":"), 3)
        dq = self.row(7, unit="DQA1~DQB1")
        self.assertEqual(dq["allele_DRB1"], "")
        self.assertEqual(dq["unit_two_field"].count("~"), 1)

    def test_coverage_and_flags(self):
        self.assertEqual(self.row(5)["coverage_flag"], "low_coverage")
        self.assertEqual(self.row(7)["coverage_flag"], "pass")
        self.assertEqual(self.row(399)["coverage_flag"], "not_measured")
        self.assertEqual(self.row(7)["coverage_cutoff_status"], "provisional_until_59")
        # fallback is normalised: typical people sit at 1.0, cutoff 0.6 on that scale
        self.assertEqual(self.row(7)["coverage_scale"], "norm")
        self.assertAlmostEqual(float(self.row(7)["coverage_value"]), 1.0, places=3)
        self.assertAlmostEqual(float(self.row(7)["mhc_rel_depth"]), 0.9)
        self.assertEqual(float(self.row(7)["coverage_cutoff"]), 0.6)
        self.assertEqual(self.row(7)["unrelated"], "True")
        rel = {self.row(10)["unrelated"], self.row(11)["unrelated"]}
        self.assertEqual(rel, {"True", "False"})
        self.assertIn(self.row(9)["platform"], ("revio", "sequel2e"))
        self.assertEqual(self.row(9)["ancestry_pred"], "EUR")

    def test_public_outputs(self):
        ids = set(self.u["person_id"])
        for f in os.listdir(self.out):
            txt = open(os.path.join(self.out, f)).read()
            for pid in ids:
                self.assertNotIn(pid, txt, f)
        fr = pd.read_csv(os.path.join(self.out, "unit_freq_by_ancestry.tsv"), sep="\t", dtype=str,
                         keep_default_na=False)
        named = fr[fr["rank"] != ""]
        self.assertTrue((named["n_haps"].astype(int) >= 20).all())
        eur = fr[(fr["unit"] == "DQA1~DQB1") & (fr["ancestry"] == "EUR") & (fr["rank"] == "1")].iloc[0]
        self.assertEqual(eur["haplotype"], "DQA1*05:01~DQB1*02:01")
        self.assertAlmostEqual(float(eur["freq_pct"]), 70, delta=6)
        mid = fr[fr["ancestry"] == "MID"]
        self.assertTrue(mid["n_resolved_haps"].isin(["<20", ""]).all())
        qc = pd.read_csv(os.path.join(self.out, "unit_qc_summary.tsv"), sep="\t", dtype=str,
                         keep_default_na=False)
        for v in qc["n_haps"]:
            if v.isdigit():
                self.assertFalse(0 < int(v) < 20)
        # complementary suppression: a group with one primary-masked cell has a second masked cell
        for _, g in qc.groupby(["unit", "field", "ancestry"]):
            m = g["n_haps"].isin(["<20", "<20*"])
            if (g["n_haps"] == "<20").sum() == 1 and len(g) > 1:
                self.assertGreaterEqual(m.sum(), 2)

    def test_no_depth_file(self):
        internal = os.path.join(self.tmp, "internal_nodepth")
        hu.main(["--table1", self.t1, "--cohort-membership", self.coh, "--relatedness-table",
                 os.path.join(self.tmp, "missing_rel.tsv"), "--depth-table",
                 os.path.join(self.tmp, "nope.tsv"), "--internal-dir", internal,
                 "--out-dir", os.path.join(self.tmp, "pub2")])
        u = pd.read_csv(os.path.join(internal, "haplotype_units.tsv.gz"), sep="\t", comment="#",
                        dtype=str)
        self.assertTrue((u["coverage_flag"] == "not_measured").all())
        self.assertTrue((u["unrelated"] == "not_assessed").all())

    def test_missing_table_uses_cached_pair_list(self):
        from unittest import mock
        internal = os.path.join(self.tmp, "internal_fb")
        with mock.patch.object(hu.C, "FALLBACK_RELATEDNESS", self.rel), \
                mock.patch.object(hu.C.load_relatedness_pairs, "__defaults__", (self.rel,)):
            hu.main(["--table1", self.t1, "--cohort-membership", self.coh, "--relatedness-table",
                     os.path.join(self.tmp, "missing_rel.tsv"), "--internal-dir", internal,
                     "--out-dir", os.path.join(self.tmp, "pub3")])
        read = lambda d: pd.read_csv(os.path.join(d, "haplotype_units.tsv.gz"), sep="\t", comment="#",
                                     dtype=str).drop_duplicates("person_id").set_index("person_id")["unrelated"]
        fb, ref = read(internal), read(self.internal)
        self.assertTrue(fb.isin(["True", "False"]).all())
        self.assertTrue((fb == ref.loc[fb.index]).all())   # same result as passing the table directly


class TestCoverageFrom59(unittest.TestCase):
    def test_exclude_people_file_sets_flag(self):
        with tempfile.TemporaryDirectory() as d:
            t1, coh, rel, depth = make_inputs(d, n=60)
            ex = os.path.join(d, "exclude_people.tsv")
            rows = []
            for i in range(50):  # people 50-59 absent -> not_measured
                pid = str(7_000_000 + i)
                mhc = 0.5 if i == 7 else 1.0  # person 7 thin on 59's scale; person 5 fine here
                kir = 0.4 if i == 8 else 1.0  # person 8 excluded for KIR only
                rows += [dict(person_id=pid, region="window_KIR", mapq=0, scale="norm", cutoff=0.65,
                              value=kir, below_cutoff=kir < 0.65, exclude=(mhc < 0.6 or kir < 0.65)),
                         dict(person_id=pid, region="window_MHC", mapq=0, scale="norm", cutoff=0.6,
                              value=mhc, below_cutoff=mhc < 0.6, exclude=(mhc < 0.6 or kir < 0.65))]
            with open(ex, "w") as fh:
                fh.write("# INTERNAL — contains n<20 cells, do not export\n")
                pd.DataFrame(rows).to_csv(fh, sep="\t", index=False)
            internal = os.path.join(d, "int")
            hu.main(["--table1", t1, "--cohort-membership", coh, "--relatedness-table", rel,
                     "--depth-table", depth, "--exclude-people", ex, "--internal-dir", internal,
                     "--out-dir", os.path.join(d, "pub")])
            u = pd.read_csv(os.path.join(internal, "haplotype_units.tsv.gz"), sep="\t", comment="#",
                            dtype=str, keep_default_na=False).drop_duplicates("person_id").set_index("person_id")
            self.assertEqual(u.loc["7000007", "coverage_flag"], "low_coverage")
            self.assertEqual(u.loc["7000005", "coverage_flag"], "pass")   # 59 overrides the 58 fallback
            self.assertEqual(u.loc["7000008", "coverage_flag"], "pass")   # KIR-only exclusion ignored
            self.assertEqual(u.loc["7000055", "coverage_flag"], "not_measured")
            self.assertEqual(set(u["coverage_cutoff_status"]), {"set_by_59"})
            self.assertEqual(u.loc["7000007", "coverage_cutoff"], "0.6")
            self.assertEqual(u.loc["7000007", "coverage_scale"], "norm")
            # 59 proposed no MHC cutoff -> everyone measured passes
            pd.DataFrame(rows).assign(cutoff=lambda x: np.where(x["region"] == "window_MHC", np.nan,
                                                                 x["cutoff"])).to_csv(ex, sep="\t", index=False)
            hu.main(["--table1", t1, "--cohort-membership", coh, "--relatedness-table", rel,
                     "--depth-table", depth, "--exclude-people", ex, "--internal-dir", internal,
                     "--out-dir", os.path.join(d, "pub")])
            u = pd.read_csv(os.path.join(internal, "haplotype_units.tsv.gz"), sep="\t", comment="#",
                            dtype=str, keep_default_na=False).drop_duplicates("person_id").set_index("person_id")
            self.assertEqual(u.loc["7000007", "coverage_flag"], "pass")
            self.assertEqual(u.loc["7000007", "coverage_cutoff"], "")


class TestInternalRefusal(unittest.TestCase):
    def test_refuses_repo_and_reports(self):
        with tempfile.TemporaryDirectory() as d:
            t1, coh, rel, depth = make_inputs(d, n=30)
            for bad in (os.path.join(SRC, "x_internal"), os.path.join(d, "reports", "x")):
                with self.assertRaises(ValueError):
                    hu.main(["--table1", t1, "--cohort-membership", coh, "--relatedness-table",
                             rel, "--depth-table", depth, "--internal-dir", bad,
                             "--out-dir", os.path.join(d, "pub")])
            self.assertFalse(os.path.exists(os.path.join(SRC, "x_internal")))
            with self.assertRaises(ValueError):
                hu.write_internal_gz(pd.DataFrame({"a": [1]}), os.path.join(SRC, "u.tsv.gz"))


class TestDictionary(unittest.TestCase):
    def test_repo_readme_matches_script(self):
        if not os.path.exists(hu.REPO_DICTIONARY):
            self.skipTest("dictionary not yet written")
        self.assertEqual(open(hu.REPO_DICTIONARY).read(), hu.dictionary_markdown())

    def test_dictionary_covers_columns(self):
        md = hu.dictionary_markdown()
        for c in hu.COLUMNS:
            key = c if not c.split("_")[-1] in ("DRB1", "DQA1", "DQB1") else c.rsplit("_", 1)[0] + "_<GENE>"
            self.assertIn(f"`{key}`", md, c)


if __name__ == "__main__":
    unittest.main()
