import importlib.util
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("s93", os.path.join(HERE, "..", "93_shortlist_ipd_release.py"))
s93 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s93)


class T(unittest.TestCase):
    def test_status(self):
        g = {"01:01": "AAAACCCC", "02:01": "CCCC"}
        self.assertEqual(s93.status(g, "AAAACCCC"), ("exact", ["01:01"]))
        self.assertEqual(s93.status(g, "TTCCCCTT")[0], "partial")
        self.assertEqual(s93.status(g, "GGGGGGGG")[0], "still_novel")


if __name__ == "__main__":
    unittest.main()
