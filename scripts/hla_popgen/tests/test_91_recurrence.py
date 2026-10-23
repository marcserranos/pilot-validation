"""Unit tests for 91 novel-call recurrence helpers (synthetic data only)."""
import importlib.util
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("r91", os.path.join(HERE, "..", "91_novel_recurrence.py"))
r91 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r91)


class T(unittest.TestCase):
    def test_supp(self):
        self.assertEqual(r91.supp(0), "0")
        self.assertEqual(r91.supp(5), "<20")
        self.assertEqual(r91.supp(20), "20")

    def test_greedy_unrelated(self):
        kin = {"a": {"b"}, "b": {"a", "c"}, "c": {"b"}}
        keep = r91.greedy_unrelated(["a", "b", "c", "d"], kin)
        self.assertEqual(sorted(keep), ["a", "c", "d"])

    def test_parents_of(self):
        first = {"c": {"p1", "p2", "s"}, "p1": {"c"}, "p2": {"c"}, "s": {"c"}}
        anyk = {"c": {"p1", "p2", "s"}, "s": {"c", "p1", "p2"}, "p1": {"c", "s"}, "p2": {"c", "s"}}
        self.assertEqual(r91.parents_of("c", first, anyk), ("p1", "p2"))  # p1/p2 have no kin edge between them
        anyk["p1"].add("p2")
        self.assertIsNone(r91.parents_of("c", first, anyk), "related first-degree pair is not a parent pair")

    def test_signature_snv_only(self):
        core = r91.load_core()
        ref = "ATGGCCGTCATGGCGCCCCGAACC"
        cs = ":20*ga:3"  # core.parse_cs reads *<obs><ref>
        s = r91.signature(core, ref, "REF*01|%s|" % cs)
        self.assertIsNotNone(s)
        self.assertEqual(s[0], "REF*01")
        self.assertEqual(s[1], ("21:A>G",))
        self.assertIsNone(r91.signature(core, ref, "REF*01|:20+a:4|"))  # indel -> None


if __name__ == "__main__":
    unittest.main()
