#!/usr/bin/env python3
"""Synthetic tests for 78_kir_gene_depth.py and 78b_kir_pruning_evidence.py (no samtools, no AoU data).
Run: python3 -m pytest scripts/hla_popgen/tests/test_78_kir_depth_pruning.py -q"""
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


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(SRC, fn))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


d78 = load("m78", "78_kir_gene_depth.py")
p78 = load("m78b", "78b_kir_pruning_evidence.py")


def gtf_line(chrom, kind, s, e, attrs, strand="+"):
    return "\t".join([chrom, "HAVANA", kind, str(s), str(e), ".", strand, ".", attrs])


def write_gtf(path):
    L = []
    T = 'gene_id "G1"; transcript_id "T1"; gene_name "KIR3DL1"; transcript_type "protein_coding"; tag "MANE_Select";'
    L.append(gtf_line("chr19", "gene", 54801001, 54805000, 'gene_id "G1"; gene_name "KIR3DL1";'))
    L.append(gtf_line("chr19", "transcript", 54801001, 54805000, T))
    for n, (s, e) in enumerate([(54801001, 54801200), (54802001, 54802300), (54804001, 54805000)], 1):
        L.append(gtf_line("chr19", "exon", s, e, T + ' exon_number "%d";' % n))
    for s, e in [(54801101, 54801200), (54802001, 54802300), (54804001, 54804500)]:
        L.append(gtf_line("chr19", "CDS", s, e, T))
    # same gene name on chr6 must be ignored; an HLA gene on chr19 must be ignored
    L.append(gtf_line("chr6", "gene", 100, 500, 'gene_id "G9"; gene_name "KIR3DL1";'))
    L.append(gtf_line("chr19", "gene", 54900001, 54901000, 'gene_id "G2"; gene_name "LILRB1";'))
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


class DepthTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.gtf = os.path.join(self.tmp, "g.gtf")
        write_gtf(self.gtf)

    def test_features_chr19_reference_transcript(self):
        feats, missing = d78.build_kir_features(self.gtf)
        self.assertEqual(set(feats["gene"]), {"KIR3DL1"})
        self.assertIn("KIR3DX1", missing)
        cds = feats[feats["kind"] == "cds"]
        self.assertEqual(int((cds["end"] - cds["start"]).sum()), 100 + 300 + 500)
        self.assertFalse(feats["ars"].astype(bool).any())          # no peptide-binding-exon convention for KIR
        self.assertTrue((feats["chrom"] == "chr19").all())
        self.assertTrue(d78.M70.parse_gtf.__name__ != "on_chr19")   # patch restored

    def test_regions_and_controls(self):
        self.assertEqual(d78.parse_regions("chr19:51,000,001-51,100,000, chr19:5-9"), [("chr19", 51000001, 51100000), ("chr19", 5, 9)])
        paths = d78.write_control_beds(self.tmp, d78.parse_regions(d78.DEFAULT_CONTROL))
        self.assertEqual(len(paths), 2)
        self.assertEqual(open(paths[0][0]).read().split("\t")[1], "51000000")

    def fake_depth(self, level0=40.0, level20=30.0, bad_gene=False):
        def f(bam, bed, mapq, timeout, samtools, region):
            rows = []
            for ln in open(bed):
                c, s, e = ln.split("\t")[:3]
                for p in range(int(s) + 1, int(e) + 1):
                    rows.append((c, p, level20 if mapq == 20 else level0))
            return pd.DataFrame(rows, columns=["chrom", "pos", "depth"])
        return f

    def test_measure_and_metrics(self):
        feats, _ = d78.build_kir_features(self.gtf)
        regions = d78.M70.gene_regions(feats)
        bed = os.path.join(self.tmp, "gene_regions.bed")
        regions.to_csv(bed, sep="\t", header=False, index=False)
        ctrl_beds = d78.write_control_beds(self.tmp, [("chr19", 51000001, 51000100)])
        bam = os.path.join(self.tmp, "x.bam")
        open(bam, "w").write("x")
        person = dict(person_id="p1", bam=bam)
        fr, pr, cr, why = d78.measure_person(person, bed, regions, feats, ctrl_beds, 10, depth=self.fake_depth())
        self.assertIsNone(why)
        self.assertEqual({r["mapq"] for r in cr}, {0, 20})
        self.assertAlmostEqual([r["control_depth"] for r in cr if r["mapq"] == 20][0], 30.0)
        gm = d78.gene_metrics(pd.DataFrame(fr), pd.DataFrame(pr), pd.DataFrame(cr), feats)
        r = gm.iloc[0]
        self.assertAlmostEqual(r["rel_cds"], 1.0)
        self.assertAlmostEqual(r["mapq_frac"], 30 / 40.0)
        self.assertAlmostEqual(r["rel_cds_q0"], 40 / 40.0)      # control also 40 at MAPQ0
        self.assertEqual(r["cds_share_lt10"], 0.0)
        missing = d78.measure_person(dict(person_id="p2", bam=os.path.join(self.tmp, "none.bam")), bed, regions, feats, ctrl_beds, 10)
        self.assertEqual(missing[3], "bam_missing")
        dist = d78.distribution_table(pd.concat([gm.assign(person_id=str(i)) for i in range(25)]))
        self.assertTrue((dist["gene"] == "KIR3DL1").all())

    def test_absent_genes_listed(self):
        self.assertTrue(set(d78.NOT_ON_GRCH38).isdisjoint(d78.KIR_ON_CHR19))
        self.assertEqual(len(d78.KIR_ON_CHR19) + len(d78.NOT_ON_GRCH38), 18)


def synth(tmp, n_trio=420, seed=3):
    """Write 77 / 77c / 78 INTERNAL files for n_trio synthetic trios; errors are enriched at high template distance."""
    rng = np.random.default_rng(seed)
    genes = ["KIR3DL3", "KIR2DL4", "KIR3DL1", "KIR2DL2", "KIR2DS4"]
    calls, haps, trio, pair = [], [], [], []
    for t in range(n_trio):
        c, p, m = "c%d" % t, "p%d" % t, "m%d" % t
        for who in (c, p, m):
            for hp in ("hap1", "hap2"):
                cls = "complete" if rng.random() > 0.2 else "multi_contig"
                haps.append(dict(person_id=who, hap=hp, hap_class=cls, n_calls=5, n_genes=5, n_contigs=1 if cls == "complete" else 2))
                for g in genes:
                    td = rng.gamma(2.0, 20.0)
                    calls.append(dict(person_id=who, hap=hp, contig="x", gene=g, copy_index=1, consensus="%s*0010101" % g,
                                      novelty_tier="known", label3="%s*001" % g, label5="%s*00101" % g, gene_start=100, gene_end=900,
                                      template_distance=td, cds_distance=td / 10, cds_mut="", template_warning="", n_tied=1,
                                      hap_class=cls, n_contigs_hap=1 if cls == "complete" else 2, contig_len=100000,
                                      dist_start=99, dist_end=99000, n_called=2))
        for g in genes:
            for hp in ("hap1", "hap2"):
                td = [x["template_distance"] for x in calls if x["person_id"] == c and x["hap"] == hp and x["gene"] == g][0]
                bad = int(rng.random() < (0.30 if td > 110 else 0.01))
                trio.append(dict(child=c, tx=p, other=m, gene=g, locus=g, mode="gene", hap=hp, copy_index=1, origin="pat",
                                 strong=int(rng.random() < 0.8), stream="named", cls3="neither" if bad else "match",
                                 tx_complete=True, tx_n_calls=2))
        pair.append(dict(person=c, other=p, locus=genes[t % 5], mode="gene", kind="pair",
                         cls="no_share" if rng.random() < 0.05 else "share", bad=0, both_complete=True))
    pair = pd.DataFrame(pair)
    pair["bad"] = (pair["cls"] == "no_share").astype(int)
    i77, i77c, i78 = (os.path.join(tmp, x) for x in ("i77", "i77c", "i78"))
    for d in (i77, i77c, i78):
        os.makedirs(d)
    pd.DataFrame(calls).to_csv(os.path.join(i77, "kir_calls.INTERNAL.tsv"), sep="\t", index=False)
    pd.DataFrame(haps).to_csv(os.path.join(i77, "kir_haps.INTERNAL.tsv"), sep="\t", index=False)
    pd.DataFrame(trio).to_csv(os.path.join(i77c, "records_trio_kir.INTERNAL.tsv"), sep="\t", index=False)
    pair.to_csv(os.path.join(i77c, "records_pair_kir.INTERNAL.tsv"), sep="\t", index=False)
    return i77, i77c, i78


class PruningTests(unittest.TestCase):
    def test_helpers(self):
        self.assertEqual(p78.count_aa("A|cs|aa,B|cs|aa"), 2)
        self.assertEqual(p78.count_aa(""), 0)
        self.assertEqual(p78.has_warning("."), 0)
        self.assertEqual(p78.has_warning("low_identity"), 1)
        comp = p78.components([("a", "b"), ("b", "c"), ("x", "y")])
        self.assertEqual(comp["a"], comp["c"])
        self.assertNotEqual(comp["a"], comp["x"])
        m = p78.mask_frame(pd.DataFrame({"n_units": [5, 25, 0], "captured_pct": [1.0, 2.0, 3.0], "errors_removed": [2.5, 40.0, 0.0]}))
        self.assertEqual(list(m["n_units"]), ["<20", "25", "0"])
        self.assertEqual(m["errors_removed"][0], "<20")
        self.assertEqual(list(m["captured_pct"]), [1.0, 2.0, 3.0])

    def test_views_and_end_to_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            i77, i77c, i78 = synth(tmp)
            views, ct = p78.build_views(i77, i77c, i78)
            self.assertEqual(set(views), {"B_child", "B_child_strong", "B_child_anytx", "B_worst", "B_worst_strong", "P_pair_all",
                                          "P_pair_complete", "H_hap"})
            b = views["B_child"]
            self.assertGreater(b["bad"].sum(), 20)
            # depth was not provided: depth metrics are NaN, not 0, and unusable
            self.assertTrue(b["o_rel_cds"].isna().all())
            self.assertNotIn("rel_cds", p78.usable_metrics(b))
            self.assertIn("template_distance", p78.usable_metrics(b))
            # the planted signal: template distance predicts `neither`
            self.assertGreater(p78.S.auc_score(b["o_template_distance"], b["bad"]), 0.65)
            # worst >= child on every row where both exist
            w = views["B_worst"]
            self.assertTrue((w["o_template_distance"].fillna(0) >= views["B_child"]["o_template_distance"].fillna(0)).all())
            # the hap view must not use the outcome-defining metrics
            h = views["H_hap"]
            self.assertFalse(any(c in h.columns for c in ("o_n_contigs_hap", "o_end_proximity")))
            self.assertTrue(set(h["bad"]) <= {0, 1})
            out, inte = os.path.join(tmp, "out"), os.path.join(tmp, "int")
            p78.run(i77, i77c, i78, inte, out, boot=20, folds=3, reps=2, grid_max=10.0, step=0.5)
            for fn in ("prune_curves.tsv", "prune_summary.tsv", "auc.tsv", "power.tsv", "per_stratum_auc.tsv", "definition_card_kir.tsv",
                       "POWER_STATEMENT.txt", "depth_metric_coverage.tsv"):
                self.assertTrue(os.path.exists(os.path.join(out, fn)), fn)
            pc = pd.read_csv(os.path.join(out, "prune_curves.tsv"), sep="\t", low_memory=False)
            self.assertIn("B_child", set(pc["view"]))
            self.assertTrue(os.path.exists(os.path.join(inte, "unmasked", "prune_curves.tsv")))
            au = pd.read_csv(os.path.join(inte, "unmasked", "auc.tsv"), sep="\t")
            a = au[(au["view"] == "B_child") & (au["metric"] == "template_distance")]["auc"].iloc[0]
            self.assertGreater(a, 0.65)
            self.assertIn("min detectable AUC", open(os.path.join(out, "POWER_STATEMENT.txt")).read())


if __name__ == "__main__":
    unittest.main()
