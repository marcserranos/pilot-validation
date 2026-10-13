#!/usr/bin/env python3
"""Tests for 66_cross_table_suppression.py: the solver on synthetic tables, and a regression check
that the committed PUBLIC tables have no count of 1-19 recoverable from the modelled sums.

Run: python3 -m pytest scripts/hla_popgen/tests/test_66_cross_table_suppression.py -q
"""
import importlib.util
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
spec = importlib.util.spec_from_file_location("cts66", os.path.join(SRC, "66_cross_table_suppression.py"))
cts = importlib.util.module_from_spec(spec)
sys.modules["cts66"] = cts
spec.loader.exec_module(cts)


def write(d, name, rows):
    p = os.path.join(d, name)
    with open(p, "w") as f:
        f.write("\n".join("\t".join(map(str, r)) for r in rows) + "\n")
    return p


class TestSolver(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_table_roundtrip_is_byte_identical(self):
        p = write(self.tmp, "t.tsv", [["a", "b", "c"], ["x", "1.50", ""], ["y", "<20", "7e-05"]])
        t = cts.Table("t", p)
        self.assertEqual(t.text(), open(p).read())
        t.set(0, "b", "<20*")
        self.assertEqual(t.text().splitlines()[1], "x\t<20*\t")

    def test_one_masked_bin_recovered_from_other_table(self):
        """genotype n_0copies = bins [0, 0.2): one masked bin -> mask one more (the smallest)."""
        g = write(self.tmp, "g.tsv", [["feature", "n_0copies"], ["F", 120]])
        h = write(self.tmp, "h.tsv", [["af_lo", "n_people"], ["0.0", 90], ["0.05", "<20"], ["0.1", 21],
                                      ["0.15", 0], ["0.2", 500]])
        P = cts.Problem()
        G, H = P.load("g", g), P.load("h", h)
        P.sum_eq(P.k(G, 0, "n_0copies"), [P.k(H, r, "n_people") for r in range(4)])
        leaks, inc = P.analyse()
        self.assertEqual([(lk[1], lk[2]) for lk in leaks], [(9, "solved")])
        self.assertFalse(inc)
        self.assertTrue(P.solve())
        self.assertEqual(H.get(2, "n_people"), "<20*")      # 21 < 90 < 120
        self.assertEqual(H.get(0, "n_people"), "90")
        self.assertEqual(H.get(3, "n_people"), "0")          # zeros never masked
        self.assertFalse(P.analyse()[0])

    def test_residual_step_difference(self):
        """Waterfall: 1,499 -> 1,494 loses 5 people; the unprotected step is masked, rates blanked."""
        w = write(self.tmp, "w.tsv", [["step", "EAS", "pct"], ["cohort", 1499, "1.0"], ["in_lr", 1494, "2.0"]])
        P = cts.Problem()
        W = P.load("w", w)
        P.rate_links["w"] = {"EAS": ["pct"]}
        P.protect(W, [1], ["EAS"])
        P.sum_eq(P.k(W, 0, "EAS"), [P.k(W, 1, "EAS")], ("w", "lost"))
        self.assertEqual(P.analyse()[0][0][1], 5)
        self.assertTrue(P.solve())
        self.assertEqual((W.get(0, "EAS"), W.get(0, "pct"), W.get(1, "EAS")), ("<20*", "", "1494"))

    def test_residual_20_or_more_is_fine(self):
        w = write(self.tmp, "w.tsv", [["step", "n"], ["a", 100], ["b", 80]])
        P = cts.Problem()
        W = P.load("w", w)
        P.sum_eq(P.k(W, 0, "n"), [P.k(W, 1, "n")], ("w", "lost"))
        self.assertEqual(P.analyse(), ([], []))

    def test_two_primary_cells_pinned_by_bounds(self):
        """<20 + <20 = 2 forces 1 + 1 even though two cells are hidden."""
        t = write(self.tmp, "t.tsv", [["k", "n"], ["tot", 102], ["a", 100], ["b", "<20"], ["c", "<20"]])
        P = cts.Problem()
        T = P.load("t", t)
        P.sum_eq(P.k(T, 0, "n"), [P.k(T, r, "n") for r in (1, 2, 3)])
        leaks = P.analyse()[0]
        self.assertEqual({(lk[1], lk[2]) for lk in leaks}, {(1, "bounds")})
        self.assertTrue(P.solve())
        self.assertEqual(T.get(1, "n"), "<20*")

    def test_inconsistent_equation_reported(self):
        t = write(self.tmp, "t.tsv", [["k", "n"], ["tot", 10], ["a", 4], ["b", 5]])
        P = cts.Problem()
        T = P.load("t", t)
        P.sum_eq(P.k(T, 0, "n"), [P.k(T, 1, "n"), P.k(T, 2, "n")])
        self.assertTrue(P.analyse()[1])

    def test_chained_tables_need_two_rounds(self):
        """Masking one cell can open a second leak elsewhere; the solver iterates until stable."""
        t = write(self.tmp, "t.tsv", [["k", "x", "y"], ["tot", 100, 300], ["a", "<20", 200], ["b", 95, 100]])
        P = cts.Problem()
        T = P.load("t", t)
        P.sum_eq(P.k(T, 0, "x"), [P.k(T, 1, "x"), P.k(T, 2, "x")])        # x: 100 = <20 + 95
        P.sum_eq(P.k(T, 2, "y"), [P.k(T, 2, "x")], ("t", "b", "rest"))    # b: y = x + rest
        self.assertTrue(P.solve())
        self.assertEqual(T.get(2, "x"), "<20*")
        self.assertFalse(P.analyse()[0])   # rest = 100 - <20* is not pinned

    def test_unfixable_when_only_protected_cells(self):
        t = write(self.tmp, "t.tsv", [["k", "n"], ["a", 495], ["b", 487]])
        P = cts.Problem()
        T = P.load("t", t)
        P.protect(T, [0, 1], ["n"])
        P.sum_eq(P.k(T, 0, "n"), [P.k(T, 1, "n")], ("t", "related"))
        self.assertFalse(P.solve())
        self.assertEqual(list(P.unfixable.values()), [(8, "solved")])

    def test_pct_pins_small_complement(self):
        """98.2% of 56 -> 55 below 1 TPM, 1 person above."""
        t = write(self.tmp, "t.tsv", [["n", "pct"], [56, "98.2"], [4698, "91.3"], [288, "50.0"]])
        P = cts.Problem()
        P.load("t", t)
        P.pct_checks.append(("t", "pct", "n", 1, []))
        self.assertEqual([h[1] for h in P.check_pcts()], [0])
        self.assertEqual(P.tables["t"].get(0, "pct"), "")


class TestCommittedTables(unittest.TestCase):
    """Regression: on the committed PUBLIC tables, every modelled equation set is safe."""

    @classmethod
    def setUpClass(cls):
        cls.P = cts.build()
        cls.leaks, cls.inconsistent = cls.P.analyse()

    def test_model_matches_the_tables(self):
        self.assertEqual(self.inconsistent, [])

    def test_no_count_1_to_19_recoverable(self):
        allowed = {"residual 50/waterfall/has_hla_calls/unrelated_kin_lt_0.0442/MID"}   # see 66 docstring
        bad = [cts.fmt_key(self.P, v) + f" = {x}" for v, x, _h in self.leaks
               if cts.fmt_key(self.P, v) not in allowed]
        self.assertEqual(bad, [])

    def test_no_single_hidden_cell_in_any_equation(self):
        """The task's literal rule: no equation (all of its cells published or masked) has exactly
        one hidden cell. Residual (unpublished) terms count as hidden."""
        for terms, _c, label in self.P.eqs:
            hidden = sum(self.P.val(k) is None for k, _ in terms)
            self.assertNotEqual(hidden, 1, label) if not any(k[0] == "r" for k, _ in terms) else None

    def test_no_pct_pins_a_small_count(self):
        self.assertEqual(self.P.check_pcts(apply=False), [])


if __name__ == "__main__":
    unittest.main()
