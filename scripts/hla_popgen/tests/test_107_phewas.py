import importlib.util
import os
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("p107", os.path.join(HERE, "..", "107_phewas_engine.py"))
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)


class T(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(4)
        n = 4000
        self.Z = np.column_stack([np.ones(n), rng.integers(0, 2, n).astype(float), rng.normal(size=n)])
        self.X = rng.binomial(2, 0.2, size=(n, 5)).astype(float)
        lin = -2.0 + 0.5 * self.Z[:, 1] + 0.9 * self.X[:, 0]           # planted effect of exposure 0 on outcome 0
        y0 = rng.binomial(1, 1 / (1 + np.exp(-lin))).astype(float)
        y1 = rng.binomial(1, 0.1, n).astype(float)                     # null outcome
        self.Y = np.column_stack([y0, y1])

    def test_planted_effect_found(self):
        U, I = p.score_tests(self.X, self.Y, self.Z)
        z = p.zmat(U, I)
        self.assertGreater(abs(z[0, 0]), 8)
        self.assertLess(abs(z[1, 0]), 4)
        self.assertAlmostEqual(U[0, 0] / I[0, 0], 0.9, delta=0.2)
        self.assertLess(np.abs(z[:, 1]).max(), 4.5)

    def test_near_constant_exposure_stays_calibrated(self):
        rng = np.random.default_rng(9)
        n = 6000
        Z = np.column_stack([np.ones(n), rng.normal(size=n)])
        X = np.column_stack([2.0 - rng.binomial(1, 0.01, n), rng.binomial(2, 0.3, n).astype(float)])   # first: allele at 99% frequency
        Y = rng.binomial(1, 0.01, size=(n, 40)).astype(float)
        z = p.zmat(*p.score_tests(X, Y, Z))
        self.assertLess(np.abs(z).max(), 5.0)
        self.assertLess(abs(np.median(z ** 2) / 0.4549 - 1), 0.5)

    def test_bh_monotone(self):
        q = p.bh(np.array([0.001, 0.01, 0.04, 0.5]))
        self.assertTrue(all(q[i] <= q[i + 1] for i in range(3)))
        self.assertLessEqual(q.max(), 1)


if __name__ == "__main__":
    unittest.main()
