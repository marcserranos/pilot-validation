import importlib.util
import os
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("w103", os.path.join(HERE, "..", "103_expression_within_g_group.py"))
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)


class T(unittest.TestCase):
    def test_three_field(self):
        self.assertEqual(w.three_field("C*07:01:01:02"), "C*07:01:01")

    def test_within_stat_detects_class_effect(self):
        res = np.array([0.0] * 30 + [1.0] * 30)
        groups = {"G1": {"a": list(range(30)), "b": list(range(30, 60))}}
        self.assertGreater(w.within_stat(res, groups), 10)
        flat = np.zeros(60)
        self.assertEqual(w.within_stat(flat, groups), 0.0)

    def test_residualise(self):
        X = np.column_stack([np.ones(10), np.arange(10.0)])
        y = 2 + 3 * np.arange(10.0)
        self.assertLess(abs(w.residualise(y, X)).max(), 1e-9)


if __name__ == "__main__":
    unittest.main()
