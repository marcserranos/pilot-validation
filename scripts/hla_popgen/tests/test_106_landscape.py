import importlib.util
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("l106", os.path.join(HERE, "..", "106_nonclassical_landscape.py"))
l = importlib.util.module_from_spec(spec)
spec.loader.exec_module(l)


class T(unittest.TestCase):
    def test_two_field(self):
        self.assertEqual(l.two_field("HLA-F*01:01:01:01"), "F*01:01")
        self.assertEqual(l.two_field("G*01:05N"), "G*01:05")
        self.assertIsNone(l.two_field("HLA-F*01:new"))
        self.assertIsNone(l.two_field(None))

    def test_null_like(self):
        self.assertTrue(l.null_like("HLA-G*01:05N"))
        self.assertTrue(l.null_like("HLA-B*44:02:01:02S"))
        self.assertFalse(l.null_like("HLA-A*02:01:01:01"))

    def test_het(self):
        self.assertAlmostEqual(l.heterozygosity(["a", "a", "b", "b"]), 0.5)
        self.assertEqual(l.heterozygosity(["a", "a"]), 0.0)


if __name__ == "__main__":
    unittest.main()
