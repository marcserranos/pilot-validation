import importlib.util
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("c95", os.path.join(HERE, "..", "95_signature_context.py"))
c95 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c95)


class T(unittest.TestCase):
    def test_two(self):
        self.assertEqual(c95.two("HLA-A*02:01:01:01"), "A*02:01")
        self.assertEqual(c95.two("A*24:02"), "A*24:02")
        self.assertIsNone(c95.two("HLA-F*01:new"))
        self.assertIsNone(c95.two(None))

    def test_eff_number(self):
        self.assertAlmostEqual(c95.eff_number(["a", "a", "a", "a"]), 1.0)
        self.assertAlmostEqual(c95.eff_number(["a", "b", "c", "d"]), 4.0)
        self.assertEqual(c95.eff_number([]), 0.0)

    def test_tersile(self):
        lo, hi = c95.tersile_bounds(list(range(9)))
        self.assertEqual((lo, hi), (3, 6))


if __name__ == "__main__":
    unittest.main()
