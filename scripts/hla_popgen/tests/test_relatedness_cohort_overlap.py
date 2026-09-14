#!/usr/bin/env python3
"""Regression test for 11_relatedness_cohort_overlap.py -- specifically the coverage-counting bug
caught during fixture testing 2026-09-08: naively unioning both sides of an `lr_sr` pair counted the
short-read-only relative as if they were part of the long-read cohort, inflating "coverage" past
100%. `lr_people_covered` must only ever contain ids that are actually in `lr_ids`."""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import importlib
mod = importlib.import_module("11_relatedness_cohort_overlap")


def _cohort():
    return pd.DataFrame({
        "person_id": ["1", "2", "3", "4", "5", "6", "7", "8", "9"],
        "in_lr":     [True, True, True, False, True, True, False, False, False],
        "in_sr":     [True, True, True, True, True, True, True, True, False],
    })


def test_bucket_classification():
    rel = pd.DataFrame({
        "person_i": ["1", "2", "3", "5", "6"],
        "person_j": ["2", "4", "7", "8", "9"],
        "kin": [0.25, 0.25, 0.12, 0.02, 0.20],
    })
    out = mod.classify_pairs(rel, _cohort())
    got = dict(zip(zip(out["person_i"], out["person_j"]), out["bucket"]))
    assert got[("1", "2")] == "both_lr"       # both lr
    assert got[("2", "4")] == "lr_sr"         # 2 is lr, 4 is sr-only
    assert got[("3", "7")] == "lr_sr"         # 3 is lr, 7 is sr-only
    assert got[("5", "8")] == "lr_sr"         # 5 is lr, 8 is sr-only
    assert got[("6", "9")] == "lr_but_relative_uncalled"  # 6 is lr, 9 is neither


def test_coverage_never_counts_sr_only_relative_as_lr_covered():
    rel = pd.DataFrame({
        "person_i": ["1", "2"],
        "person_j": ["2", "4"],   # 2~4: 2 is lr, 4 is sr-only
        "kin": [0.25, 0.25],
    })
    cohort = _cohort()
    rel = mod.classify_pairs(rel, cohort)
    lr_ids = set(cohort.loc[cohort["in_lr"], "person_id"])
    _, _, lr_people_covered = mod.summarize(rel, lr_ids)
    # person "4" (sr-only) must never appear in the lr_sr coverage set
    assert "4" not in lr_people_covered["lr_sr"]
    assert lr_people_covered["lr_sr"] == {"2"}
    assert lr_people_covered["both_lr"] == {"1", "2"}


def test_kin_degree_bins():
    assert mod.kin_degree(0.5) == "duplicate_or_MZ_twin"
    assert mod.kin_degree(0.25) == "first_degree_parentchild_or_sibling"
    assert mod.kin_degree(0.12) == "second_degree"
    assert mod.kin_degree(0.05) == "third_degree"
    assert mod.kin_degree(0.01) == "below_third_degree"


if __name__ == "__main__":
    test_bucket_classification()
    test_coverage_never_counts_sr_only_relative_as_lr_covered()
    test_kin_degree_bins()
    print("OK")
