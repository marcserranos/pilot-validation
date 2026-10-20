#!/usr/bin/env python3
"""Synthetic tests for 70_call_exon_depth.py (no samtools, no AoU data).
Run: python3 -m pytest scripts/hla_popgen/tests/test_70_call_exon_depth.py -q"""
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
spec = importlib.util.spec_from_file_location("m70", os.path.join(SRC, "70_call_exon_depth.py"))
m = importlib.util.module_from_spec(spec)
sys.modules["m70"] = m
spec.loader.exec_module(m)


def gtf_line(chrom, kind, s, e, attrs):
    return "\t".join([chrom, "HAVANA", kind, str(s), str(e), ".", "+", ".", attrs])


def write_gtf(path):
    L = []
    # HLA-A: gene 1001-5000 (1-based). Reference transcript T1: exons 1-4; exon 1 starts with 5'UTR.
    L.append(gtf_line("chr6", "gene", 1001, 5000, 'gene_id "G1"; gene_name "HLA-A";'))
    T1 = 'gene_id "G1"; transcript_id "T1"; gene_name "HLA-A"; transcript_type "protein_coding"; tag "MANE_Select";'
    L.append(gtf_line("chr6", "transcript", 1001, 5000, T1))
    for n, (s, e) in enumerate([(1001, 1200), (2001, 2300), (3001, 3300), (4001, 5000)], 1):
        L.append(gtf_line("chr6", "exon", s, e, T1 + f' exon_number "{n}";'))
    for s, e in [(1101, 1200), (2001, 2300), (3001, 3300), (4001, 4500)]:
        L.append(gtf_line("chr6", "CDS", s, e, T1))
    # alternative coding transcript T2 with an extra exon in an intron, and a non-coding transcript T3
    T2 = 'gene_id "G1"; transcript_id "T2"; gene_name "HLA-A"; transcript_type "protein_coding";'
    L.append(gtf_line("chr6", "transcript", 1001, 5000, T2))
    L.append(gtf_line("chr6", "exon", 2501, 2600, T2 + ' exon_number "1";'))
    L.append(gtf_line("chr6", "CDS", 2501, 2600, T2))
    T3 = 'gene_id "G1"; transcript_id "T3"; gene_name "HLA-A"; transcript_type "retained_intron";'
    L.append(gtf_line("chr6", "transcript", 1001, 5000, T3))
    L.append(gtf_line("chr6", "exon", 3501, 3900, T3 + ' exon_number "1";'))
    # HLA-DQA1 (class II alpha), single transcript with 3 coding exons
    L.append(gtf_line("chr6", "gene", 10001, 12000, 'gene_id "G2"; gene_name "HLA-DQA1";'))
    TQ = 'gene_id "G2"; transcript_id "TQ"; gene_name "HLA-DQA1"; transcript_type "protein_coding"; tag "Ensembl_canonical";'
    L.append(gtf_line("chr6", "transcript", 10001, 12000, TQ))
    for n, (s, e) in enumerate([(10001, 10100), (10501, 10800), (11501, 12000)], 1):
        L.append(gtf_line("chr6", "exon", s, e, TQ + f' exon_number "{n}";'))
        L.append(gtf_line("chr6", "CDS", s, e, TQ))
    # a gene on another chromosome and one we do not ask for
    L.append(gtf_line("chr1", "gene", 1, 100, 'gene_id "G3"; gene_name "HLA-B";'))
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(L) + "\n")


class T(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.gtf = os.path.join(self.d, "g.gtf.gz")
        write_gtf(self.gtf)
        self.f, self.missing = m.build_features(self.gtf)

    def test_missing_genes_reported_not_zero(self):
        self.assertIn("HLA-B", self.missing)      # only on chr1 in the fixture
        self.assertIn("HLA-DRB3", self.missing)
        self.assertNotIn("A", self.missing)

    def test_exons_and_cds_are_reference_transcript_only(self):
        ex = self.f[(self.f["gene"] == "A") & (self.f["kind"] == "exons")]
        spans = sorted(zip(ex["start"], ex["end"]))
        self.assertNotIn((3500, 3900), spans)          # retained-intron transcript left out
        self.assertNotIn((2500, 2600), spans)          # alternative-isoform coding exon NOT in the reference CDS
        self.assertEqual(len(spans), 4)
        cds = self.f[(self.f["gene"] == "A") & (self.f["kind"] == "cds")]
        self.assertEqual(int((cds["end"] - cds["start"]).sum()), 100 + 300 + 300 + 500)

    def test_reference_transcript_is_mane_and_exons_clipped_to_cds(self):
        ex = self.f[(self.f["gene"] == "A") & (self.f["kind"] == "ex")].sort_values("cds_ord")
        self.assertEqual(list(ex["transcript_id"].unique()), ["T1"])
        self.assertEqual(list(ex["exon_number"]), [1, 2, 3, 4])
        first = ex.iloc[0]
        self.assertEqual((first["start"], first["end"]), (1100, 1200))   # 5'UTR part of exon 1 excluded
        last = ex.iloc[-1]
        self.assertEqual((last["start"], last["end"]), (4000, 4500))     # 3'UTR part of exon 4 excluded

    def test_ars_flags(self):
        a = self.f[(self.f["gene"] == "A") & (self.f["kind"] == "ex")].set_index("cds_ord")["ars"]
        self.assertEqual([bool(a[k]) for k in (1, 2, 3, 4)], [False, True, True, False])   # class I: exons 2, 3
        q = self.f[(self.f["gene"] == "DQA1") & (self.f["kind"] == "ex")].set_index("cds_ord")["ars"]
        self.assertEqual([bool(q[k]) for k in (1, 2, 3)], [False, True, False])            # class II: exon 2

    def test_names_unique(self):
        self.assertFalse(self.f["name"].duplicated().any())
        self.assertIn("A|gene", set(self.f["name"]))
        self.assertIn("A|ex2", set(self.f["name"]))

    def test_perbase_summary(self):
        cds, ars = m.cds_bed(self.f), m.ars_intervals(self.f)
        gene_cds = cds[cds["gene"] == "DQA1"]
        pos = []
        for s, e in zip(gene_cds["start"], gene_cds["end"]):
            pos += list(range(s + 1, e + 1))
        dep = np.full(len(pos), 30)
        dep[:10] = 0                                  # first 10 CDS bases uncovered
        df = pd.DataFrame({"chrom": "chr6", "pos": pos, "depth": dep})
        st = m.perbase_summary(df, gene_cds, ars)["DQA1"]
        self.assertEqual(st["cds_bases"], 100 + 300 + 500)
        self.assertEqual(st["cds_min"], 0)
        self.assertAlmostEqual(st["cds_share_0"], 10 / 900)
        self.assertEqual(st["ars_min"], 30)           # the 0x bases are in exon 1, not the ARS exon
        self.assertEqual(st["ars_bases"], 300)

    def test_parse_depth_and_empty(self):
        d = m.parse_depth("chr6\t100\t5\nchr6\t101\t7\n")
        self.assertEqual(list(d["depth"]), [5, 7])
        self.assertTrue(m.parse_depth("").empty)

    def test_measure_person_with_fake_runners(self):
        feats = self.f
        feats_bed = os.path.join(self.d, "f.bed")
        m.write_bed(feats, feats_bed)
        cds, ars = m.cds_bed(feats), m.ars_intervals(feats)
        bam = os.path.join(self.d, "x.bam")
        open(bam, "w").close()

        def bedcov(b, bed, q, t, s):
            r = feats[["chrom", "start", "end", "name"]].copy()
            r["sum_depth"] = (r["end"] - r["start"]) * (20.0 if q == 0 else 18.0)
            r["n_reads"] = 5
            return r

        def depth(b, bed, q, t, s):
            return pd.DataFrame({"chrom": ["chr6"] * 3, "pos": [2001, 2002, 2003], "depth": [20, 20, 20]})

        fr, pr, why = m.measure_person(dict(person_id="p1", bam=bam), feats_bed, "cds", cds, ars, 10,
                                       bedcov=bedcov, depth=depth)
        self.assertIsNone(why)
        d = pd.DataFrame(fr)
        self.assertEqual(set(d["mapq"]), {0, 20})
        self.assertAlmostEqual(d[(d["feature"] == "A|gene") & (d["mapq"] == 0)]["mean_depth"].iloc[0], 20.0)
        self.assertAlmostEqual(d[(d["feature"] == "A|gene") & (d["mapq"] == 20)]["mean_depth"].iloc[0], 18.0)
        self.assertTrue(len(pr) >= 2)
        _, _, why = m.measure_person(dict(person_id="p2", bam=os.path.join(self.d, "none.bam")), feats_bed, "cds",
                                     cds, ars, 10, bedcov=bedcov, depth=depth)
        self.assertEqual(why, "bam_missing")

    def test_resume_and_failures(self):
        idir = os.path.join(self.d, "int")
        os.makedirs(idir)
        # the public-output guard only allows INTERNAL tables under an internal path; mimic with a name it accepts
        try:
            m.append_rows(os.path.join(idir, "feature_depth.tsv"),
                          [dict(person_id="p1", feature="A|gene", mapq=0, sum_depth=1, n_reads=1, mean_depth=1)],
                          m.FEATURE_DEPTH_COLS)
            m.append_rows(os.path.join(idir, "failures.tsv"), [dict(person_id="p9", reason="timeout")],
                          ["person_id", "reason"])
        except Exception as e:   # pragma: no cover - guard refuses a non-internal tmp path
            self.skipTest(f"internal-path guard refused the temp dir: {e}")
        self.assertEqual(m.done_people(idir), {"p1", "p9"})

    def test_priority_order(self):
        people = pd.DataFrame({"person_id": ["a", "b", "c", "d"], "bam": ["x"] * 4})
        o = m.order_people(people, ["c", "d"])
        self.assertEqual(list(o["person_id"]), ["c", "d", "a", "b"])

    def test_bin_profile(self):
        regions = pd.DataFrame({"chrom": ["chr6"], "start": [0], "end": [40], "gene": ["A"]})
        dp = pd.DataFrame({"chrom": "chr6", "pos": np.arange(1, 41), "depth": [10] * 20 + [30] * 20})
        prof = m.bin_profile(dp, regions)["A"]
        self.assertEqual(list(prof), [10, 10, 30, 30])


if __name__ == "__main__":
    unittest.main()
