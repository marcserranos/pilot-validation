import importlib.util
import os
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("e102", os.path.join(HERE, "..", "102_expression_dosage_tests.py"))
e102 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e102)


class T(unittest.TestCase):
    def test_ols_recovers_slope(self):
        rng = np.random.default_rng(1)
        x = rng.integers(0, 3, 500).astype(float)
        c = rng.normal(size=500)
        y = 1.0 + 0.5 * c + 2.0 * x + rng.normal(scale=0.3, size=500)
        X = np.column_stack([np.ones(500), c, x])
        b, se, p = e102.ols(y, X)
        self.assertAlmostEqual(b, 2.0, delta=0.1)
        self.assertLess(p, 1e-10)

    def test_ols_null(self):
        rng = np.random.default_rng(2)
        x = rng.integers(0, 3, 500).astype(float)
        y = rng.normal(size=500)
        b, se, p = e102.ols(y, np.column_stack([np.ones(500), x]))
        self.assertGreater(p, 0.001)

    def test_supp(self):
        self.assertEqual(e102.supp(5), "<20")
        self.assertEqual(e102.supp(0), "0")


if __name__ == "__main__":
    unittest.main()
