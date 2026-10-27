import gzip
import importlib.util
import os
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("e105", os.path.join(HERE, "..", "105_novel_allele_evidence.py"))
e = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e)


class T(unittest.TestCase):
    def test_plus_strand(self):
        m = {"strand": "+", "exons": [(100, 109), (200, 209)]}
        self.assertEqual(e.cds_to_genomic(m, 1), 100)
        self.assertEqual(e.cds_to_genomic(m, 10), 109)
        self.assertEqual(e.cds_to_genomic(m, 11), 200)
        self.assertEqual(e.cds_to_genomic(m, 20), 209)
        self.assertIsNone(e.cds_to_genomic(m, 21))

    def test_minus_strand(self):
        m = {"strand": "-", "exons": [(200, 209), (100, 109)]}    # transcript order = descending
        self.assertEqual(e.cds_to_genomic(m, 1), 209)
        self.assertEqual(e.cds_to_genomic(m, 10), 200)
        self.assertEqual(e.cds_to_genomic(m, 11), 109)
        self.assertEqual(e.genomic_bases(m, "A", "G"), ("T", "C"))
        self.assertEqual(e.genomic_bases({"strand": "+"}, "A", "G"), ("A", "G"))

    def test_pileup_parse(self):
        c = e.parse_pileup_bases(".,..^~,GgG$.+2AC-1T*")
        self.assertEqual((c["ref"], c["G"], c["del"]), (6, 3, 1))
        self.assertEqual(c["A"], 0)

    def test_models(self):
        with tempfile.NamedTemporaryFile("wb", suffix=".gtf.gz", delete=False) as f:
            pass
        with gzip.open(f.name, "wt") as g:
            for s, en, tr in ((100, 109, "T1"), (200, 209, "T1"), (100, 105, "T2")):
                g.write('chr6\tHAVANA\tCDS\t%d\t%d\t.\t+\t0\tgene_name "HLA-F"; transcript_id "%s";\n' % (s, en, tr))
        m = e.load_cds_models(f.name, ["HLA-F"])
        self.assertEqual(m["HLA-F"]["transcript"], "T1")
        self.assertEqual(m["HLA-F"]["cds_len"], 20)


if __name__ == "__main__":
    unittest.main()
