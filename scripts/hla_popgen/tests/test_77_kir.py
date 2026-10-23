#!/usr/bin/env python3
"""Synthetic tests for 77 / 77b / 77c (KIR call table, trio origin helpers, Mendelian records). python3 tests/test_77_kir.py"""
import gzip
import importlib.util
import os
import shutil
import sys
import tempfile
import unittest

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


m77, m77b, m77c = load("m77", "77_kir_call_table.py"), load("m77b", "77b_kir_trio_origin.py"), load("m77c", "77c_kir_trio_records.py")


def gtf(path, calls):
    """calls: (contig, gene, start, end, consensus, copy, cds_distance, cds_mut)"""
    with gzip.open(path, "wt") as f:
        f.write("# synthetic\n")
        for ctg, g, s, e, cons, cp, cd, mut in calls:
            gid = g if cp == 1 else "%s.%d" % (g, cp)
            f.write('%s\tImmuannot\tgene\t%d\t%d\t.\t+\t.\tgene_id "%s"; gene_name "%s"; template_distance 2;\n' % (ctg, s, e, gid, g))
            a = 'gene_id "%s"; transcript_id "%s.t"; gene_name "%s"; consensus "%s"; alleles "%s,%s"; template_distance 2;' % (
                gid, gid, g, cons, cons, cons)
            if cd is not None:
                a += ' cds_distance %d; cds_mut "%s";' % (cd, mut or "")
            f.write("%s\tImmuannot\ttranscript\t%d\t%d\t.\t+\t.\t%s\n" % (ctg, s, e, a))


class T(unittest.TestCase):
    def test_labels(self):
        self.assertEqual(m77.protein_field("KIR2DL1*0030204"), "KIR2DL1*003")
        self.assertEqual(m77.protein_field("KIR3DL2*00201new"), "KIR3DL2*002")
        self.assertEqual(m77.coding_field("KIR2DL1*0030204"), "KIR2DL1*00302")
        self.assertIsNone(m77.coding_field("KIR2DL1*003"))
        self.assertEqual(m77.labels("KIR2DL1*0030204", "known"), ("KIR2DL1*003", "KIR2DL1*00302"))
        self.assertEqual(m77.labels("KIR2DL1*0030204", "novel_cds_synonymous"), ("KIR2DL1*003", None))
        self.assertEqual(m77.labels("KIR2DL1*0030204", "novel_protein"), (None, None))
        self.assertEqual(m77.labels("undetermined", "undetermined"), (None, None))

    def test_call_table(self):
        d = tempfile.mkdtemp()
        try:
            for pid, h1, h2 in (("1", [("c1:1-300000", "KIR3DL3", 1000, 3000, "KIR3DL3*0030101", 1, None, None),
                                       ("c1:1-300000", "KIR3DL2", 200000, 202000, "KIR3DL2*00201new", 1, 0, None)],
                                      [("c2:5-1000", "KIR3DL3", 50, 900, "KIR3DL3*0010101", 1, None, None)]),
                                ("2", [], [])):
                os.makedirs(os.path.join(d, pid, "immuannot_output"))
                gtf(os.path.join(d, pid, "immuannot_output", "hap1.gtf.gz"), h1)
                gtf(os.path.join(d, pid, "immuannot_output", "hap2.gtf.gz"), h2)
            calls, haps = m77.run_chunk((["1", "2"], d))
            self.assertEqual(len(calls), 3)
            hc = dict(zip(zip(haps.person_id, haps.hap), haps.hap_class))
            self.assertEqual(hc[("1", "hap1")], "complete")
            self.assertEqual(hc[("1", "hap2")], "one_contig_open_end")
            self.assertEqual(hc[("2", "hap1")], "empty")
            calls = m77.add_person_level(calls, haps)
            r = calls[(calls.gene == "KIR3DL3")]
            self.assertEqual(set(r.n_called), {2})
            self.assertEqual(set(calls.person_n_haps), {2})
            x = calls[calls.gene == "KIR3DL2"].iloc[0]
            self.assertEqual(x.label3, "KIR3DL2*002")
            self.assertEqual(x.novelty_tier, "novel_genomic_known_cds")
            self.assertEqual(int(x.dist_start), 199999)
            self.assertEqual(int(x.n_tied), 2)
        finally:
            shutil.rmtree(d)

    def test_origin_and_windows(self):
        self.assertEqual(m77b.origin(10, 2), "pat")
        self.assertEqual(m77b.origin(2, 10), "mat")
        self.assertEqual(m77b.origin(4, 0), "amb")
        self.assertEqual(m77b.origin(9, 4), "amb")
        self.assertEqual(m77b.strong(12, 2), 1)
        self.assertEqual(m77b.strong(9, 0), 0)
        seqs = {"c1:1-100": "ACGT" * 5000}
        txt, meta, nowin = m77b.make_windows(seqs, [("KIR2DL1", 1, "c1:1-100", 100, 300), ("KIR2DL3", 1, "zz", 1, 2)])
        self.assertEqual(meta, [("KIR2DL1", 1, "c1:1-100", 1, 10300)])
        self.assertEqual(nowin, [("KIR2DL3", 1, "zz")])
        self.assertTrue(txt.startswith(">w0\n"))

    def test_classes(self):
        def e(*l):
            return {"labels": set(l), "n": len(l) or 1, "unnamed": False}
        self.assertEqual(m77c.classify("a", e("a", "b"), e("c")), "match")
        self.assertEqual(m77c.classify("a", e("b"), e("a")), "only_other_parent")
        self.assertEqual(m77c.classify("a", e("b"), e("c")), "neither")
        self.assertEqual(m77c.classify("a", None, None), "absent_in_parents")
        self.assertEqual(m77c.classify("a", None, e("a")), "tx_absent_other_has")

    def _calls(self):
        rows = []
        def add(p, hap, gene, lab, cls="complete"):
            rows.append(dict(person_id=p, hap=hap, gene=gene, copy_index=1, label3=lab, label5=None, novelty_tier="known",
                             hap_class=cls, template_distance=0, cds_distance=0, n_called=1))
        add("c", "hap1", "KIR2DL1", "KIR2DL1*003"); add("c", "hap2", "KIR2DL1", "KIR2DL1*004")
        add("c", "hap1", "KIR3DL1", "KIR3DL1*001")          # parents lack 3DL1 but have 3DS1 (group mode)
        add("c", "hap2", "KIR2DL3", None)                  # novel
        add("p", "hap1", "KIR2DL1", "KIR2DL1*003"); add("p", "hap2", "KIR2DL1", "KIR2DL1*009")
        add("p", "hap1", "KIR3DS1", "KIR3DS1*013")
        add("m", "hap1", "KIR2DL1", "KIR2DL1*004"); add("m", "hap2", "KIR2DL1", "KIR2DL1*010")
        return pd.DataFrame(rows)

    def test_trio_records(self):
        calls = self._calls()
        origin = pd.DataFrame(dict(child=["c"] * 4, hap=["hap1", "hap2", "hap1", "hap2"], gene=["KIR2DL1", "KIR2DL1", "KIR3DL1", "KIR2DL3"],
                                   copy_index=1, pat_k=[20, 20, 20, 20], mat_k=[0, 0, 0, 0], origin=["pat"] * 4, strong=[1] * 4))
        comp = {"c": True, "p": True, "m": False}
        r = m77c.trio_records(calls, origin, {"c": ("p", "m")}, comp, "gene")
        by = {(x.gene, x.hap): x for x in r.itertuples()}
        self.assertEqual(by[("KIR2DL1", "hap1")].cls3, "match")
        self.assertEqual(by[("KIR2DL1", "hap2")].cls3, "only_other_parent")   # 004 only in the other parent
        self.assertEqual(by[("KIR3DL1", "hap1")].cls3, "absent_in_parents")
        self.assertEqual(by[("KIR2DL3", "hap2")].stream, "novel")
        self.assertFalse(by[("KIR2DL1", "hap1")].other_complete)
        g = m77c.trio_records(calls, origin, {"c": ("p", "m")}, comp, "group")
        self.assertEqual({(x.gene): x.cls3 for x in g.itertuples() if x.hap == "hap1"}["KIR3DL1"], "neither")
        s = m77c.summarize_trio(r)
        self.assertIn("complete_only", set(s.variant))

    def test_pairs(self):
        calls = self._calls()
        pr = m77c.pair_records(calls, [("c", "p"), ("c", "m"), ("p", "m")], {}, "gene", "pair")
        cl = {(x.person, x.other, x.locus): x.cls for x in pr.itertuples()}
        self.assertEqual(cl[("c", "p", "KIR2DL1")], "share")
        self.assertEqual(cl[("c", "p", "KIR3DL1")], "presence_discordant")
        self.assertEqual(cl[("p", "m", "KIR2DL1")], "no_share")
        a, b = m77c.select_pairs({"1", "2", "3"}, [("1", "2", 0.25), ("1", "3", 0.5), ("2", "3", 0.01), ("1", "9", .25)])
        self.assertEqual((a, b), ([("1", "2")], [("1", "3")]))


if __name__ == "__main__":
    unittest.main()
