#!/usr/bin/env python3
"""Tests for 61_kir_by_ancestry_panel.py -- synthetic KIR GTFs only, no AoU data.

Covers: protein-field truncation and per-tier allele labels, complementary suppression, the
complete-haplotype restriction (a gene seen only on fragmented haplotypes must not count),
relatives and --exclude-people removal, PUBLIC masking (<20, no named allele under 20), no person
ids in any output, the local draft from the committed 51/43 tables, and figure rendering.

Run: python3 -m pytest scripts/hla_popgen/tests/test_61_kir_by_ancestry_panel.py -q
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
spec = importlib.util.spec_from_file_location("kp61", os.path.join(SRC, "61_kir_by_ancestry_panel.py"))
kp = importlib.util.module_from_spec(spec)
sys.modules["kp61"] = kp
spec.loader.exec_module(kp)

FRAME = ["KIR3DL3", "KIR2DL4", "KIR3DP1", "KIR3DL2"]


def gtf_line(contig, gene, consensus, extra=""):
    attrs = (f'gene_id "IG{gene}"; transcript_id "IAT{gene}.1"; gene_name "{gene}"; '
             f'consensus "{consensus}";{extra}')
    return f"{contig}\tImmuannot\ttranscript\t100\t200\t.\t+\t.\t{attrs}\n"


def write_gtf(path, lines):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write("##synthetic\n")
        fh.writelines(lines)


def make_cohort(tmp, seed=7):
    """300 EUR + 300 AFR + 10 MID people. KIR2DS1 carried by 40% of EUR / 10% of AFR haps.
    20% of haplotypes are fragmented (3DL2 on a second contig) and carry KIR2DS3 -- a gene that
    therefore never appears on a complete haplotype. Five related pairs; 5 people excluded."""
    rng = np.random.default_rng(seed)
    root = os.path.join(tmp, "kir")
    people = []
    anc = ["EUR"] * 300 + ["AFR"] * 300 + ["MID"] * 10
    for i, a in enumerate(anc):
        pid = str(9_000_000 + i)
        people.append((pid, a))
        for hap in ("hap1", "hap2"):
            frag = rng.random() < 0.2
            lines = []
            for g in FRAME:
                ctg = "ctgB" if (frag and g == "KIR3DL2") else "ctgA"
                lines.append(gtf_line(ctg, g, f"{g}*0010101"))
            p = 0.4 if a == "EUR" else 0.1
            if rng.random() < p:
                lines.append(gtf_line("ctgA", "KIR2DS1", "KIR2DS1*0020103"))
            if frag:
                lines.append(gtf_line("ctgA", "KIR2DS3", "KIR2DS3*0010301"))
            allele = "KIR2DL1*0030204" if rng.random() < 0.7 else "KIR2DL1*0040101"
            if rng.random() < 0.1:
                lines.append(gtf_line("ctgA", "KIR2DL1", "KIR2DL1*00302new",
                                      " cds_distance 0;"))
            else:
                lines.append(gtf_line("ctgA", "KIR2DL1", allele))
            write_gtf(os.path.join(root, pid, "immuannot_output", f"{hap}.gtf.gz"), lines)
    coh = pd.DataFrame({"person_id": [p for p, _ in people],
                        "ancestry_pred": [a.lower() for _, a in people],
                        "platform": "revio", "in_lr": True})
    coh_path = os.path.join(tmp, "cohort.tsv")
    coh.to_csv(coh_path, sep="\t", index=False)
    rel = pd.DataFrame({"i.s": [people[2 * k][0] for k in range(5)],
                        "j.s": [people[2 * k + 1][0] for k in range(5)], "kin": 0.25})
    rel_path = os.path.join(tmp, "rel.tsv")
    rel.to_csv(rel_path, sep="\t", index=False)
    excl_path = os.path.join(tmp, "exclude.tsv")
    pd.DataFrame({"person_id": [people[i][0] for i in range(100, 105)],
                  "exclude": True}).to_csv(excl_path, sep="\t", index=False)
    return root, coh_path, rel_path, excl_path, people


class TestHelpers(unittest.TestCase):
    def test_protein_field(self):
        self.assertEqual(kp.protein_field("KIR2DL1*0030204"), "KIR2DL1*003")
        self.assertEqual(kp.protein_field("KIR3DL2*00201new"), "KIR3DL2*002")
        self.assertEqual(kp.protein_field("KIR2DL1*01102"), "KIR2DL1*011")
        self.assertIsNone(kp.protein_field("undetermined"))

    def test_allele_label_tiers(self):
        row = {"consensus": "KIR2DL1*0030204", "novelty_tier": "known"}
        self.assertEqual(kp.allele_label(row, "protein"), "KIR2DL1*003")
        self.assertEqual(kp.allele_label(row, "full_known"), "KIR2DL1*0030204")
        syn = {"consensus": "KIR2DL1*00302new", "novelty_tier": "novel_cds_synonymous"}
        self.assertEqual(kp.allele_label(syn, "protein"), "KIR2DL1*003")
        self.assertIsNone(kp.allele_label(syn, "full_known"))
        prot = {"consensus": "KIR2DL1*00302new", "novelty_tier": "novel_protein"}
        self.assertEqual(kp.allele_label(prot, "protein"), kp.NOVEL_PROTEIN_LABEL)
        und = {"consensus": "undetermined", "novelty_tier": "undetermined"}
        self.assertIsNone(kp.allele_label(und, "protein"))

    def test_suppress_partition(self):
        self.assertEqual(kp.suppress_partition({"a": 30, "b": 100}), set())
        self.assertEqual(kp.suppress_partition({"a": 5, "b": 100}), {"a", "b"})
        self.assertEqual(kp.suppress_partition({"a": 5, "b": 10, "c": 30, "d": 100}),
                         {"a", "b", "c"})
        self.assertEqual(kp.suppress_partition({"a": 0, "b": 100}), set())

    def test_allele_rows_never_name_small(self):
        rows = kp.allele_rows("KIR2DL1", "EUR", {"X*001": 100, "X*002": 25, "X*003": 10},
                              "protein", top_n=5)
        named = {r["allele"] for r in rows}
        self.assertNotIn("X*003", named)
        other = [r for r in rows if r["allele"] == "other"]
        # remainder 10 would be recoverable -> X*002 folded into other (35)
        self.assertEqual(other[0]["n_haps"], "35")
        self.assertNotIn("X*002", named)

    def test_exclude_people_formats(self):
        with tempfile.TemporaryDirectory() as d:
            p1 = os.path.join(d, "a.txt")
            open(p1, "w").write("123\n456\n")
            self.assertEqual(kp.read_exclude_people(p1), {"123", "456"})
            p2 = os.path.join(d, "b.tsv")
            pd.DataFrame({"person_id": ["1", "2"], "exclude": [True, False]}).to_csv(
                p2, sep="\t", index=False)
            self.assertEqual(kp.read_exclude_people(p2), {"1"})
            # 59's INTERNAL exclude_people.tsv: header comment, two rows per person
            p3 = os.path.join(d, "exclude_people.tsv")
            with open(p3, "w") as fh:
                fh.write("# INTERNAL — contains n<20 cells, do not export\n")
                pd.DataFrame({"person_id": ["7", "7", "8", "8"],
                              "region": ["window_KIR", "window_MHC"] * 2, "mapq": 0,
                              "scale": "norm", "cutoff": 0.6, "value": [0.5, 0.9, 0.9, 0.9],
                              "below_cutoff": [True, False, False, False],
                              "exclude": [True, True, False, False]}).to_csv(fh, sep="\t", index=False)
            self.assertEqual(kp.read_exclude_people(p3), {"7"})


class TestVmMode(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        root, coh, rel, excl, cls.people = make_cohort(cls.tmp)
        cls.out = os.path.join(cls.tmp, "out")
        kp.main(["--mode", "vm", "--kir-outroot", root, "--cohort-membership", coh,
                 "--relatedness-table", rel, "--exclude-people", excl, "--out-dir", cls.out,
                 "--workers", "1", "--no-figure"])
        cls.t = {f[:-4]: pd.read_csv(os.path.join(cls.out, f), sep="\t", dtype=str)
                 for f in os.listdir(cls.out) if f.endswith(".tsv")}

    def test_complete_only(self):
        g = self.t["gene_freq_by_ancestry"]
        ds3 = g[(g["gene"] == "KIR2DS3") & (g["ancestry"] == "ALL")].iloc[0]
        self.assertEqual(ds3["n_haps_carrying"], "0")  # only on fragmented haplotypes
        l3 = g[(g["gene"] == "KIR3DL3") & (g["ancestry"] == "EUR")].iloc[0]
        self.assertEqual(float(l3["hap_freq_pct"]), 100.0)
        comp = self.t["completeness_by_ancestry"]
        allrow = comp[comp["ancestry"] == "ALL"].iloc[0]
        self.assertAlmostEqual(float(allrow["complete_pct"]), 80, delta=4)

    def test_ancestry_gradient(self):
        g = self.t["gene_freq_by_ancestry"]
        eur = float(g[(g["gene"] == "KIR2DS1") & (g["ancestry"] == "EUR")]["hap_freq_pct"].iloc[0])
        afr = float(g[(g["gene"] == "KIR2DS1") & (g["ancestry"] == "AFR")]["hap_freq_pct"].iloc[0])
        self.assertAlmostEqual(eur, 40, delta=6)
        self.assertAlmostEqual(afr, 10, delta=4)
        cb = self.t["content_by_ancestry"].set_index("ancestry")
        self.assertGreater(float(cb.loc["EUR", "cB_pct"]), float(cb.loc["AFR", "cB_pct"]))

    def test_relatives_and_exclusions_dropped(self):
        rs = dict(self.t["run_summary"][["metric", "value"]].values)
        self.assertEqual(rs["people_kept_unrelated_not_excluded"], str(610 - 5 - 5))
        self.assertEqual(rs["people_excluded_by_list"], "<20")

    def test_masking(self):
        g = self.t["gene_freq_by_ancestry"]
        mid = g[(g["ancestry"] == "MID")]
        self.assertTrue((mid["n_complete_haps"] == "<20").all())
        self.assertTrue(mid["hap_freq_pct"].isna().all())
        for name, df in self.t.items():
            for col in df.columns:
                if col.startswith("n_") and col != "n_genes":
                    for v in df[col].dropna():
                        if str(v).isdigit():
                            self.assertFalse(0 < int(v) < 20, f"{name}.{col}={v}")
        al = self.t["top_alleles_by_ancestry"]
        named = al[~al["allele"].isin(["other", kp.NOVEL_PROTEIN_LABEL])]
        self.assertTrue(named["n_haps"].str.isdigit().all())
        self.assertTrue((named["n_haps"].astype(int) >= 20).all())
        prot = al[(al["resolution"] == "protein") & (al["gene"] == "KIR2DL1") & (al["ancestry"] == "ALL")]
        top = prot[prot["rank"] == "1"].iloc[0]
        self.assertEqual(top["allele"], "KIR2DL1*003")
        self.assertGreater(float(top["freq_pct"]), 70)  # genomic-novel 003 calls keep protein 003

    def test_no_person_ids(self):
        ids = {p for p, _ in self.people}
        for f in os.listdir(self.out):
            txt = open(os.path.join(self.out, f), errors="ignore").read()
            for pid in list(ids)[:50] + [self.people[0][0], self.people[-1][0]]:
                self.assertNotIn(pid, txt, f)

    def test_figure_mode(self):
        kp.main(["--mode", "figure", "--tables-dir", self.out])
        self.assertTrue(os.path.exists(os.path.join(self.out, "fig_kir_by_ancestry_panel.png")))


class TestLocalMode(unittest.TestCase):
    def test_local_draft_from_committed(self):
        root = os.path.join(os.path.dirname(os.path.dirname(SRC)), "reports", "hla_popgen")
        if not os.path.exists(os.path.join(root, "51_kir_structure_qc",
                                           "gene_content_haplotypes.tsv")):
            self.skipTest("committed 51 tables not present")
        t = kp.local_tables(root)
        g = t["gene_freq_by_ancestry"]
        self.assertTrue((g[g["gene"] == "KIR3DL3"]["hap_freq_pct"] == "100.0").all())
        lo = pd.to_numeric(g["hap_freq_lower_bound"])
        hi = pd.to_numeric(g["hap_freq_upper_bound"])
        est = pd.to_numeric(g["hap_freq_pct"])
        self.assertTrue(((lo <= est + 1e-9) & (est <= hi + 1e-9)).all())
        cb = t["content_by_ancestry"].set_index("ancestry")
        self.assertEqual(cb.loc["SAS", "cB_pct"], "57.1")
        self.assertTrue(t["gene_freq_by_ancestry"]["basis"].str.startswith("draft").all())
        with tempfile.TemporaryDirectory() as d:
            kp.main(["--mode", "local", "--report-root", root, "--out-dir", d])
            self.assertTrue(os.path.exists(os.path.join(d, "fig_kir_by_ancestry_panel_DRAFT.pdf")))


if __name__ == "__main__":
    unittest.main()
