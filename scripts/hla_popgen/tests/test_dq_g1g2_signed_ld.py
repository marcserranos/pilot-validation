#!/usr/bin/env python3
"""Unit tests for scripts/hla_popgen/37_dq_g1g2_signed_ld.py's statistics, on synthetic data.

Covers: signed D' at the three textbook extremes (perfect repulsion -> -1, independence -> 0,
perfect coupling -> +1), the Petersdorf G1/G2 classifier, and a sanity check on the
observed/expected purge statistic (independence -> O/E ~= 1).

Run: python3 scripts/hla_popgen/tests/test_dq_g1g2_signed_ld.py
"""
import importlib.util
import os
import random
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)
if HLA_POPGEN_DIR not in sys.path:
    sys.path.insert(0, HLA_POPGEN_DIR)


def _load_module(filename, modname):
    path = os.path.join(HLA_POPGEN_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


mod = _load_module("37_dq_g1g2_signed_ld.py", "dq_g1g2_signed_ld")
m29 = mod.m29()

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print("  [%s] %s%s" % (status, name, (" -- " + detail) if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


def _table(alleles_a, alleles_b):
    return m29.contingency(alleles_a, alleles_b)


def test_perfect_repulsion():
    # Two alleles per locus, exactly 50 haplotypes each of A1~B2 and A2~B1, ZERO of A1~B1/A2~B2:
    # A and B are perfectly anti-correlated -> D' = -1 for A1~B1 (and A2~B2).
    a = ["A1"] * 50 + ["A2"] * 50
    b = ["B2"] * 50 + ["B1"] * 50
    M, a_lab, b_lab = _table(a, b)
    rows = mod.signed_dprime_table(M, a_lab, b_lab, min_haps=20)
    look = {(r["allele_a"], r["allele_b"]): r["signed_Dprime"] for r in rows}
    check("perfect repulsion A1~B1 D'=-1", abs(look[("A1", "B1")] - (-1.0)) < 1e-9,
          detail=str(look.get(("A1", "B1"))))
    check("perfect repulsion A2~B2 D'=-1", abs(look[("A2", "B2")] - (-1.0)) < 1e-9,
          detail=str(look.get(("A2", "B2"))))
    check("perfect repulsion A1~B2 D'=+1 (coupled)", abs(look[("A1", "B2")] - 1.0) < 1e-9)


def test_perfect_coupling():
    # A1 always with B1, A2 always with B2 -> D'=+1 on the diagonal.
    a = ["A1"] * 60 + ["A2"] * 40
    b = ["B1"] * 60 + ["B2"] * 40
    M, a_lab, b_lab = _table(a, b)
    rows = mod.signed_dprime_table(M, a_lab, b_lab, min_haps=20)
    look = {(r["allele_a"], r["allele_b"]): r["signed_Dprime"] for r in rows}
    check("perfect coupling A1~B1 D'=+1", abs(look[("A1", "B1")] - 1.0) < 1e-9,
          detail=str(look.get(("A1", "B1"))))
    check("perfect coupling A2~B2 D'=+1", abs(look[("A2", "B2")] - 1.0) < 1e-9)


def test_independence():
    # Alleles assigned independently at random with fixed marginal frequencies -> D' ~= 0.
    rng = random.Random(42)
    n = 20000
    a = rng.choices(["A1", "A2"], weights=[0.6, 0.4], k=n)
    b = rng.choices(["B1", "B2"], weights=[0.3, 0.7], k=n)
    M, a_lab, b_lab = _table(a, b)
    rows = mod.signed_dprime_table(M, a_lab, b_lab, min_haps=20)
    look = {(r["allele_a"], r["allele_b"]): r["signed_Dprime"] for r in rows}
    d = look[("A1", "B1")]
    check("independence D' near 0 (|D'|<0.05 at n=20000)", abs(d) < 0.05, detail=str(d))


def test_min_haps_floor_excludes_rare_alleles():
    # A3 carried by only 5 haplotypes must not appear as a row even though A1/A2 pass.
    a = ["A1"] * 50 + ["A2"] * 50 + ["A3"] * 5
    b = ["B1"] * 50 + ["B2"] * 50 + ["B1"] * 5
    M, a_lab, b_lab = _table(a, b)
    rows = mod.signed_dprime_table(M, a_lab, b_lab, min_haps=20)
    alleles_seen = {r["allele_a"] for r in rows}
    check("rare allele A3 (n=5) excluded by the 20-haplotype floor", "A3" not in alleles_seen,
          detail=str(alleles_seen))


def test_g1g2_classification():
    check("DQA1*01:01 x DQB1*05:01 is G2", mod.classify_pair("DQA1*01:01", "DQB1*05:01") == "G2")
    check("DQA1*03:01 x DQB1*03:02 is G1", mod.classify_pair("DQA1*03:01", "DQB1*03:02") == "G1")
    check("DQA1*01:01 x DQB1*03:01 is predicted_incompatible (G2 A x G1 B)",
          mod.classify_pair("DQA1*01:01", "DQB1*03:01") == "predicted_incompatible")
    check("DQA1*05:01 x DQB1*05:01 is predicted_incompatible (G1 A x G2 B)",
          mod.classify_pair("DQA1*05:01", "DQB1*05:01") == "predicted_incompatible")
    check("malformed allele -> unclassified", mod.classify_pair("garbage", "DQB1*05:01") == "unclassified")


def test_observed_expected_independence():
    # Draw DQA1/DQB1 alleles independently from group-only marginals -> O/E should hover near 1.
    rng = random.Random(7)
    n = 8000
    a_alleles = rng.choices(["DQA1*01:01", "DQA1*03:01"], weights=[0.5, 0.5], k=n)  # G2, G1
    b_alleles = rng.choices(["DQB1*05:01", "DQB1*03:01"], weights=[0.5, 0.5], k=n)  # G2, G1
    oe = mod.observed_expected_incompatible(a_alleles, b_alleles, n_bootstrap=200,
                                            rng=random.Random(1))
    check("independent groups: O/E within 0.85-1.15", 0.85 <= oe["oe_ratio"] <= 1.15,
          detail=str(oe["oe_ratio"]))
    check("O/E CI brackets 1.0 under independence",
          oe["oe_ci_lo"] <= 1.0 <= oe["oe_ci_hi"],
          detail="[%s, %s]" % (oe["oe_ci_lo"], oe["oe_ci_hi"]))


def test_observed_expected_purge():
    # Construct a population with the cross-group combination systematically absent (the
    # "purged" biology G1/G2 predicts) -> O/E should be near 0, well below the CI's lower bound
    # under independence.
    n = 4000
    a_alleles = ["DQA1*01:01"] * (n // 2) + ["DQA1*03:01"] * (n // 2)
    b_alleles = ["DQB1*05:01"] * (n // 2) + ["DQB1*03:01"] * (n // 2)  # G2-G2 and G1-G1 only
    oe = mod.observed_expected_incompatible(a_alleles, b_alleles, n_bootstrap=200,
                                            rng=random.Random(1))
    check("purged population: O/E far below 1", oe["oe_ratio"] < 0.2, detail=str(oe["oe_ratio"]))
    check("purged population: observed count is 0", oe["observed"] == 0)


def main():
    print("test_dq_g1g2_signed_ld.py")
    test_perfect_repulsion()
    test_perfect_coupling()
    test_independence()
    test_min_haps_floor_excludes_rare_alleles()
    test_g1g2_classification()
    test_observed_expected_independence()
    test_observed_expected_purge()
    if FAILURES:
        print("\n%d FAILURE(S): %s" % (len(FAILURES), ", ".join(FAILURES)))
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
