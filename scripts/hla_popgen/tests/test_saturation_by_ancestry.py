#!/usr/bin/env python3
"""Unit tests for 39_saturation_by_ancestry.py's core statistics, on synthetic data only.

Covers: rarefaction_threshold_curves (distinct-count + carrier-threshold accounting),
summarize_curve (mean/percentile bands), equal_n_comparison (bootstrap difference direction),
and fit_and_extrapolate (Clench/Chao2 wiring against 04_allele_saturation.py).

Run: python3 scripts/hla_popgen/tests/test_saturation_by_ancestry.py
"""
import importlib.util
import os
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


sat = _load_module("39_saturation_by_ancestry.py", "saturation_by_ancestry")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


def test_rarefaction_final_step_matches_true_counts():
    # 5 "people", each contributing a small allele set. At the FINAL step (all 5 included, in any
    # order), distinct-count and >=k-carrier counts must equal the ground truth regardless of the
    # random ordering -- this is the one thing that must hold no matter how permutation shuffles.
    unit_sets = [
        {"A*01:01", "A*02:01"},
        {"A*01:01"},
        {"A*01:01", "A*03:01"},
        {"A*02:01"},
        {"A*01:01", "A*02:01", "A*03:01"},
    ]
    # ground truth carrier counts: A*01:01 -> 4, A*02:01 -> 3, A*03:01 -> 2
    curves = sat.rarefaction_threshold_curves(unit_sets, [1, 2, 10], n_permutations=30, seed=1)
    final_distinct = curves["distinct"][:, -1]
    check("final-step distinct count is always 3 (ground truth), across all permutations",
          bool(np.all(final_distinct == 3)), str(final_distinct))
    final_ge1 = curves[1][:, -1]
    final_ge2 = curves[2][:, -1]
    final_ge10 = curves[10][:, -1]
    check(">=1 carrier count at final step == 3 (all 3 alleles seen at least once)",
          bool(np.all(final_ge1 == 3)), str(final_ge1))
    check(">=2 carrier count at final step == 3 (all 3 alleles reach >=2 carriers)",
          bool(np.all(final_ge2 == 3)), str(final_ge2))
    check(">=10 carrier count at final step == 0 (no allele reaches 10 carriers)",
          bool(np.all(final_ge10 == 0)), str(final_ge10))


def test_distinct_curve_is_monotone_nondecreasing():
    rng = np.random.default_rng(0)
    alphabet = [f"X*{i:02d}" for i in range(20)]
    unit_sets = [set(rng.choice(alphabet, size=rng.integers(1, 4), replace=False))
                 for _ in range(60)]
    curves = sat.rarefaction_threshold_curves(unit_sets, [1, 2], n_permutations=10, seed=2)
    ok = True
    for p in range(curves["distinct"].shape[0]):
        diffs = np.diff(curves["distinct"][p])
        if np.any(diffs < 0):
            ok = False
    check("cumulative distinct-allele curve never decreases within a permutation", ok)


def test_summarize_curve_bounds():
    arr = np.array([[1, 2, 3], [1, 2, 5], [1, 2, 1]])
    steps, mean, lo, hi = sat.summarize_curve(arr)
    check("summarize_curve returns steps 1..m", list(steps) == [1, 2, 3], str(steps))
    check("mean at step 3 is between min and max of that column",
          arr[:, 2].min() <= mean[2] <= arr[:, 2].max(), str(mean))
    check("lo <= mean <= hi at every step", bool(np.all(lo <= mean + 1e-9))
          and bool(np.all(mean <= hi + 1e-9)))


def test_equal_n_comparison_detects_a_real_difference():
    # AFR carries systematically more distinct alleles per person than EUR -> AFR should show up
    # with a positive, statistically distinguishable difference at equal N.
    rng = np.random.default_rng(3)
    afr_pool = [f"A*{i:02d}" for i in range(40)]   # bigger allele pool
    eur_pool = [f"A*{i:02d}" for i in range(10)]   # smaller allele pool -> saturates faster
    afr_units = [set(rng.choice(afr_pool, size=2, replace=False)) for _ in range(150)]
    eur_units = [set(rng.choice(eur_pool, size=2, replace=False)) for _ in range(150)]
    eq_df, diff_df = sat.equal_n_comparison(
        {"AFR": afr_units, "EUR": eur_units}, n_min=150, n_bootstrap=100, seed=5)
    check("equal_n_comparison returns one row per ancestry", len(eq_df) == 2, str(eq_df))
    afr_mean = eq_df.set_index("ancestry").loc["AFR", "mean_distinct_at_n_min"]
    eur_mean = eq_df.set_index("ancestry").loc["EUR", "mean_distinct_at_n_min"]
    check("AFR (larger pool) shows more distinct alleles at equal N than EUR",
          afr_mean > eur_mean, f"AFR={afr_mean}, EUR={eur_mean}")
    row = diff_df[diff_df["comparison"] == "AFR - EUR"].iloc[0]
    check("AFR - EUR difference is positive", row["mean_diff"] > 0, str(row["mean_diff"]))
    check("bootstrap p-value is well-formed (in [0, 1])",
          0.0 <= row["p_two_sided_bootstrap"] <= 1.0, str(row["p_two_sided_bootstrap"]))


def test_fit_and_extrapolate_runs_on_synthetic_saturating_data():
    # A population with only 5 alleles total, sampled with replacement many times: the discovery
    # curve should visibly saturate, and the Clench fit + Chao2 should both run without error.
    rng = np.random.default_rng(7)
    pool = [f"A*0{i}" for i in range(5)]
    unit_sets = [set(rng.choice(pool, size=1)) for _ in range(200)]
    curves = sat.rarefaction_threshold_curves(unit_sets, [1], n_permutations=25, seed=9)
    steps, mean, lo, hi = sat.summarize_curve(curves["distinct"])
    out = sat.fit_and_extrapolate(steps, mean, lo, hi, unit_sets)
    check("s_obs_now caps near the true pool size (5) for a saturating population",
          out["s_obs_now"] <= 5.5, str(out["s_obs_now"]))
    check("chao2_richness is finite and >= observed richness",
          out["chao2_richness"] >= out["s_obs_now"] - 1e-6, str(out))
    if out["clench_fit_ok"]:
        check("expected_new_per_1000 is small for an already-saturated population",
              out["expected_new_per_1000"] < 1.0, str(out["expected_new_per_1000"]))


def test_nonparametric_extrapolation_is_nonnegative_and_agrees_with_curve_direction():
    # Orchestrator review 2026-09-23: the Clench fit produced NEGATIVE "expected new alleles"
    # on real data, which is impossible -- the replacement non-parametric estimators must never
    # do that on a still-rising curve, and the empirical rate must reflect the curve's own slope.
    import pandas as pd
    rows = []
    for n in range(1, 101):
        rows.append({"scheme": "pred", "ancestry": "AFR", "gene_group": "classical_pooled",
                     "category": "all", "n": n, "mean_distinct": 5 * np.sqrt(n),
                     "lo2_5": 5 * np.sqrt(n) - 1, "hi97_5": 5 * np.sqrt(n) + 1})
        rows.append({"scheme": "pred", "ancestry": "AFR", "gene_group": "classical_pooled",
                     "category": "carriers_ge_1", "n": n, "mean_distinct": 5 * np.sqrt(n),
                     "lo2_5": 0, "hi97_5": 0})
        rows.append({"scheme": "pred", "ancestry": "AFR", "gene_group": "classical_pooled",
                     "category": "carriers_ge_2", "n": n, "mean_distinct": 4 * np.sqrt(n),
                     "lo2_5": 0, "hi97_5": 0})
    curve_df = pd.DataFrame(rows)
    fit_df = pd.DataFrame([{"scheme": "pred", "ancestry": "AFR", "gene_group": "classical_pooled",
                            "category": "all", "n_now": 100, "s_obs_now": 50.0,
                            "chao2_richness": 80.0, "expected_new_per_1000": -5.0,
                            "clench_s_max": 60.0}])
    out = sat.nonparametric_extrapolation(curve_df, fit_df, scheme="pred")
    row = out[out["ancestry"] == "AFR"].iloc[0]
    check("empirical_rate_per_1000 is positive for a still-rising curve",
          row["empirical_rate_per_1000"] > 0, str(row["empirical_rate_per_1000"]))
    check("chao_new_alleles_by_2n is non-negative (never a negative 'discovery')",
          row["chao_new_alleles_by_2n"] >= 0, str(row["chao_new_alleles_by_2n"]))
    check("Clench value is preserved verbatim as the supplementary column, not silently dropped",
          row["clench_expected_new_per_1000_supplementary"] == -5.0,
          str(row["clench_expected_new_per_1000_supplementary"]))


def test_equal_n_descriptive_excludes_requested_ancestries():
    import pandas as pd
    curve_df = pd.DataFrame([
        {"scheme": "pred", "ancestry": "AFR", "gene_group": "classical_pooled", "category": "all",
         "n": 500, "mean_distinct": 300.0, "lo2_5": 290.0, "hi97_5": 310.0},
        {"scheme": "pred", "ancestry": "MID", "gene_group": "classical_pooled", "category": "all",
         "n": 500, "mean_distinct": 277.0, "lo2_5": 277.0, "hi97_5": 277.0},
    ])
    out = sat.equal_n_descriptive(curve_df, 500, exclude=("MID",), scheme="pred")
    check("MID is excluded when requested", "MID" not in set(out["ancestry"]), str(out))
    check("AFR is still present", "AFR" in set(out["ancestry"]), str(out))


def test_format_p_never_prints_a_literal_zero():
    check("a zero bootstrap p-value is reported as a resolution bound, not '0.0000'",
          sat.format_p(0.0, 200) == "<0.005", sat.format_p(0.0, 200))
    check("a nonzero p-value above the resolution floor prints normally",
          sat.format_p(0.24, 200) == "0.2400", sat.format_p(0.24, 200))


def main():
    test_rarefaction_final_step_matches_true_counts()
    test_distinct_curve_is_monotone_nondecreasing()
    test_summarize_curve_bounds()
    test_equal_n_comparison_detects_a_real_difference()
    test_fit_and_extrapolate_runs_on_synthetic_saturating_data()
    test_nonparametric_extrapolation_is_nonnegative_and_agrees_with_curve_direction()
    test_equal_n_descriptive_excludes_requested_ancestries()
    test_format_p_never_prints_a_literal_zero()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("All tests passed.")


if __name__ == "__main__":
    main()
