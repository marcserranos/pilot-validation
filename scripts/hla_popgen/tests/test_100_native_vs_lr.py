import importlib.util
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("n100", os.path.join(HERE, "..", "100_native_vs_lr.py"))
n100 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(n100)


class T(unittest.TestCase):
    def test_two(self):
        self.assertEqual(n100.two("A*01:01:01"), "01:01")
        self.assertEqual(n100.two("HLA-A*01:01:01:02"), "01:01")
        self.assertEqual(n100.two("01:01"), "01:01")
        self.assertIsNone(n100.two("A*01"))
        self.assertIsNone(n100.two("HLA-F*01:new"))
        self.assertIsNone(n100.two(None))

    def test_compare_kinds(self):
        c = n100.compare(("01:01", "02:01"), ("02:01", "01:01"))
        self.assertEqual(c["kind"], "concordant")
        self.assertEqual(n100.compare(("01:01", "02:01"), ("01:01", "01:01"))["kind"], "lr_hom_native_het")
        self.assertEqual(n100.compare(("01:01", "01:01"), ("01:01", "02:01"))["kind"], "lr_het_native_hom")
        s = n100.compare(("01:01", "03:01"), ("01:01", "02:01"))
        self.assertEqual((s["kind"], s["conf"], s["alleles_ok"]), ("allele_swap", [("03:01", "02:01")], 1))
        self.assertEqual(n100.compare((None, "01:01"), ("01:01", "01:01"))["kind"], "uncallable")

    def test_td(self):
        self.assertEqual(n100.td_stratum([0, 0], False), "td0")
        self.assertEqual(n100.td_stratum([0, 2], False), "td1-2")
        self.assertEqual(n100.td_stratum([0, 5], False), "td3+")
        self.assertEqual(n100.td_stratum([0, 0], True), "warned")


if __name__ == "__main__":
    unittest.main()
