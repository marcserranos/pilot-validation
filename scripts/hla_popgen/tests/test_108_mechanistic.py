import importlib.util
import os
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("m108", os.path.join(HERE, "..", "108_mechanistic_rna_tests.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class T(unittest.TestCase):
    def test_leader_map(self):
        d = tempfile.mkdtemp()
        lines = [" Prot -24 1", " B*07:02:01:01 MLVMAPRTVLLLL", " B*08:01:01:01 -R-T---------"]
        open(os.path.join(d, "B_prot.txt"), "w").write("\n".join(lines) + "\n")
        for g in ("A", "C", "G"):
            open(os.path.join(d, g + "_prot.txt"), "w").write(" Prot -24 1\n %s*01:01:01:01 MAVMAPRTLLLLL\n" % g)
        mp = m.leader_m21_map(d)
        self.assertTrue(mp["B*07:02"])
        self.assertFalse(mp["B*08:01"])
        self.assertTrue(mp["A*01:01"])


if __name__ == "__main__":
    unittest.main()
