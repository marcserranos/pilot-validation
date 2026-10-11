#!/usr/bin/env python3
"""Unit tests for 38_hla_a_deletion_validation.py, on synthetic data.

Focus: the new statistics this script adds on top of 30's bridged_absences (already tested in
test_hla_structural_variation.py) -- the HWE biallelic test, within-person concordance, and the
short-read contradiction test.

Run: python3 scripts/hla_popgen/tests/test_hla_a_deletion_validation.py
"""
import importlib.util
import os
import sys

import pandas as pd

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


m38 = _load_module("38_hla_a_deletion_validation.py", "m38")

FAILS = []


def check(name, cond):
    status = "OK" if cond else "FAIL"
    print("[%s] %s" % (status, name))
    if not cond:
        FAILS.append(name)


def _calls(rows):
    """rows: list of (person_id, hap, gene, status)."""
    return pd.DataFrame(rows, columns=["person_id", "hap", "gene", "status"])


# ---------------------------------------------------------------------------
# per_person_hap_status / HWE / concordance
# ---------------------------------------------------------------------------
def test_per_person_hap_status_drops_unbridged():
    calls = _calls([
        ("p1", "hap1", "A", "present"), ("p1", "hap2", "A", "deleted_bridged"),
        ("p2", "hap1", "A", "deleted_bridged"), ("p2", "hap2", "A", "absent_unbridged"),
    ])
    status = m38.per_person_hap_status(calls, "A")
    # p2 has one unbridged haplotype -> excluded entirely, not counted as a half-observation
    check("per_person_hap_status excludes people with any unbridged hap",
          set(status["person_id"]) == {"p1"})
    check("per_person_hap_status marks p1 as one haplotype deleted",
          int(status.loc[status["person_id"] == "p1", "n_deleted_haps"].iloc[0]) == 1)


def test_hwe_matches_hand_calc():
    # 10 people fully bridged, 4 haplotype-deletions total among the 20 haplotypes -> p = 0.2
    # 1 person has both haplotypes deleted.
    rows = []
    dels = [("p1", "hap1"), ("p1", "hap2"), ("p2", "hap1"), ("p3", "hap1")]
    for i in range(1, 11):
        pid = "p%d" % i
        for hap in ("hap1", "hap2"):
            status = "deleted_bridged" if (pid, hap) in dels else "present"
            rows.append((pid, hap, "A", status))
    calls = _calls(rows)
    status_df = m38.per_person_hap_status(calls, "A")
    hwe = m38.hwe_biallelic_test(status_df)
    check("hwe p_hap_deleted == 0.2", abs(hwe["p_hap_deleted"] - 0.2) < 1e-9)
    check("hwe n_both_deleted == 1", hwe["n_both_deleted"] == 1)
    check("hwe expected == N*p^2 == 0.4", abs(hwe["expected_both_deleted"] - 0.4) < 1e-9)
    check("hwe reports a p-value", "poisson_exact_p" in hwe)


def test_within_person_concordance_independent_case():
    # Deletions on hap1 and hap2 chosen independently (no person has both) -> OR should not blow
    # up; table sums correctly.
    rows = []
    for i in range(1, 21):
        pid = "p%d" % i
        s1 = "deleted_bridged" if i <= 5 else "present"
        s2 = "deleted_bridged" if 6 <= i <= 10 else "present"
        rows.append((pid, "hap1", "A", s1))
        rows.append((pid, "hap2", "A", s2))
    calls = _calls(rows)
    status_df = m38.per_person_hap_status(calls, "A")
    conc = m38.within_person_concordance(status_df)
    total = sum(sum(row) for row in conc["table"])
    check("concordance table covers all 20 people", total == 20)
    check("concordance table has zero (1,1) cell (no overlap by construction)",
          conc["table"][1][1] == 0)


# ---------------------------------------------------------------------------
# SR contradiction test
# ---------------------------------------------------------------------------
def _sr_long(pairs):
    """pairs: dict person_id -> (a1, a2) for gene 'A'."""
    rows = []
    for pid, (a1, a2) in pairs.items():
        rows.append((pid, "A", 1, a1))
        rows.append((pid, "A", 2, a2))
    return pd.DataFrame(rows, columns=["person_id", "gene_bare", "copy", "allele"])


def test_sr_contradiction_flags_heterozygous_carrier():
    # p1: LR hemizygous-deleted at A (hap1 deleted, hap2 present). SR reports TWO distinct
    # alleles -> contradiction. p2: LR hemizygous-deleted, SR reports ONE effective allele
    # (homozygous-looking) -> no contradiction. p3/p4: LR non-carriers (both present).
    calls = _calls([
        ("p1", "hap1", "A", "deleted_bridged"), ("p1", "hap2", "A", "present"),
        ("p2", "hap1", "A", "deleted_bridged"), ("p2", "hap2", "A", "present"),
        ("p3", "hap1", "A", "present"), ("p3", "hap2", "A", "present"),
        ("p4", "hap1", "A", "present"), ("p4", "hap2", "A", "present"),
    ])
    status_df = m38.per_person_hap_status(calls, "A")
    sr_long = _sr_long({
        "p1": ("A*02:01", "A*03:01"),   # heterozygous SR -> contradicts the LR deletion
        "p2": ("A*02:01", "A*02:01"),   # homozygous-looking SR -> consistent with hemizygosity
        "p3": ("A*01:01", "A*24:02"),
        "p4": ("A*01:01", "A*01:01"),
    })
    sr_map = m38.sr_alleles_per_person(sr_long, "A")
    anc_of = {p: "EUR" for p in ("p1", "p2", "p3", "p4")}
    res = m38.sr_contradiction_test(status_df, sr_map, anc_of, "A", t1=pd.DataFrame(), calls=calls)
    check("carrier het rate is 1/2 = 50%", abs(res["carrier"]["pct_het"] - 50.0) < 1e-9)
    check("noncarrier het rate is 1/2 = 50%", abs(res["noncarrier"]["pct_het"] - 50.0) < 1e-9)
    check("sr_gene_typed is True", res["sr_gene_typed"] is True)


def test_sr_gene_not_typed_reported():
    status_df = m38.per_person_hap_status(_calls([
        ("p1", "hap1", "DRB4", "deleted_bridged"), ("p1", "hap2", "DRB4", "present"),
    ]), "DRB4")
    res = m38.sr_contradiction_test(status_df, {}, {"p1": "EUR"}, "DRB4", t1=pd.DataFrame(),
                                    calls=pd.DataFrame(), positive_control=True)
    check("untyped gene reports sr_gene_typed == False", res["sr_gene_typed"] is False)


def test_first2_raw_handles_formats():
    check("first2_raw strips gene prefix", m38.first2_raw("A*02:01:01:01") == "02:01")
    check("first2_raw handles bare fields", m38.first2_raw("02:01") == "02:01")
    check("first2_raw None on non-str", m38.first2_raw(None) is None)


if __name__ == "__main__":
    test_per_person_hap_status_drops_unbridged()
    test_hwe_matches_hand_calc()
    test_within_person_concordance_independent_case()
    test_sr_contradiction_flags_heterozygous_carrier()
    test_sr_gene_not_typed_reported()
    test_first2_raw_handles_formats()
    if FAILS:
        print("\n%d test(s) FAILED: %s" % (len(FAILS), FAILS))
        sys.exit(1)
    print("\nAll tests passed.")
