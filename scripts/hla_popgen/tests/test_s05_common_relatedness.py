#!/usr/bin/env python3
"""_s05_common.load_relatedness_pairs: AoU table, cached LR pair-list fallback, and no-file case."""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import _s05_common as C  # noqa: E402


class TestRelatedness(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.aou = os.path.join(self.tmp, "samples_relatedness.tsv")
        with open(self.aou, "w") as fh:
            fh.write("i.s\tj.s\tkin\n1\t2\t0.25\n2\t3\t0.12\n")
        self.cache = os.path.join(self.tmp, "relatedness_lr_overlap_pairs.tsv")
        with open(self.cache, "w") as fh:
            fh.write("person_i\tperson_j\tkin\tkin_degree\tbucket\n1\t2\t0.25\tfirst\tboth_lr\n")

    def test_aou_table_preferred(self):
        pairs = C.load_relatedness_pairs(self.aou, fallback=self.cache)
        self.assertEqual(len(pairs), 2)

    def test_fallback_used_when_table_missing(self):
        pairs = C.load_relatedness_pairs(os.path.join(self.tmp, "gone.tsv"), fallback=self.cache)
        self.assertEqual(pairs, [("1", "2", 0.25)])
        self.assertEqual(len(C.greedy_unrelated(["1", "2", "4"], pairs)), 2)

    def test_nothing_available(self):
        self.assertEqual(C.load_relatedness_pairs(os.path.join(self.tmp, "gone.tsv"),
                                                  fallback=os.path.join(self.tmp, "gone2.tsv")), [])


if __name__ == "__main__":
    unittest.main()
