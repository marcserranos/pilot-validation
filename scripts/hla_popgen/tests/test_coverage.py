#!/usr/bin/env python3
"""Tests for scripts/hla_popgen/_coverage.py -- sample-coverage / Chao & Jost estimators.

Runs standalone (no pytest, per repo convention -- see e.g. tests/test_novel.py):

    cd scripts/hla_popgen && python3 -c "import sys; sys.path.insert(0,'tests'); \\
        import test_coverage as t; [getattr(t,n)() for n in dir(t) if n.startswith('test_')]; \\
        print('ok')"

Also collectible by pytest if it's ever installed (plain `assert`-based test_* functions).
"""
import itertools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)
if HLA_POPGEN_DIR not in sys.path:
    sys.path.insert(0, HLA_POPGEN_DIR)

import _coverage as cov


def close(a, b, tol=1e-9):
    return abs(a - b) <= tol


# ---------------------------------------------------------------------------
# 1. sample_coverage
# ---------------------------------------------------------------------------

def test_sample_coverage_hand_computed():
    X = [1] * 10 + [5] * 10
    n = 10 * 1 + 10 * 5  # 60
    f1, f2 = 10, 0
    # f2 == 0 branch: A = (n-1)*f1 / ((n-1)*f1 + 2)
    A = (n - 1) * f1 / ((n - 1) * f1 + 2)
    expected = 1 - (f1 / n) * A
    got = cov.sample_coverage(X)
    assert close(got, expected), (got, expected)
    assert close(got, 1 - (10 / 60) * (590 / 592))


def test_sample_coverage_no_singletons():
    X = [2] * 5 + [3] * 5
    assert cov.sample_coverage(X) == 1.0


def test_sample_coverage_f2_zero_branch():
    X = [1, 1, 1, 4, 4]
    n = sum(X)
    f1, f2 = 3, 0
    A = (n - 1) * f1 / ((n - 1) * f1 + 2)
    expected = 1 - (f1 / n) * A
    assert close(cov.sample_coverage(X), expected)


def test_sample_coverage_f2_positive():
    X = [1, 1, 2, 2, 10, 10]
    n = sum(X)
    f1, f2 = 2, 2
    A = (n - 1) * f1 / ((n - 1) * f1 + 2 * f2)
    expected = 1 - (f1 / n) * A
    assert close(cov.sample_coverage(X), expected)


# ---------------------------------------------------------------------------
# 2. coverage_at_size: equality at m==n, monotonicity, m=1 closed form, brute force
# ---------------------------------------------------------------------------

def test_coverage_at_size_equals_sample_coverage_at_n():
    X = [3, 2, 1, 1, 5, 8]
    n = sum(X)
    assert close(cov.coverage_at_size(X, n), cov.sample_coverage(X))


def test_coverage_at_size_monotone():
    X = [3, 2, 1, 1, 5, 8, 1, 1, 2]
    n = sum(X)
    ms = np.arange(1, n + 1)
    vals = cov.coverage_at_size(X, ms)
    diffs = np.diff(vals)
    assert np.all(diffs >= -1e-9), diffs.min()


def test_coverage_at_size_m1_closed_form():
    # Hand-derived: coverage of a size-1 sample = P(next single draw, from the remaining
    # n-1 haplotypes, repeats the allele of the one just drawn) = sum_i p_i*(X_i-1)/(n-1).
    X = np.array([3.0, 2.0, 1.0, 1.0, 5.0, 8.0])
    n = X.sum()
    expected = float(np.sum((X / n) * (X - 1) / (n - 1)))
    got = cov.coverage_at_size(X, 1)
    assert close(got, expected), (got, expected)


def _brute_force_interp_coverage(X, m, n_reps=None):
    """Exact enumeration: average, over ALL C(n,m) subsets of the n haplotype "slots", of
    the fraction of the (n-m) HELD-OUT slots whose allele already appears in the m-subset.
    This is the quantity coverage_at_size's interpolation formula computes (see module
    docstring's item 2) -- NOT "sum of full-sample allele frequencies present in the
    subset", which is a different quantity (also checked below, and shown to differ)."""
    X = np.asarray(X, dtype=int)
    n = int(X.sum())
    labels = np.repeat(np.arange(len(X)), X)  # one label per haplotype "slot"
    idx = np.arange(n)
    total = 0.0
    count = 0
    for combo in itertools.combinations(idx, m):
        combo = set(combo)
        held_out = [i for i in idx if i not in combo]
        sub_labels = set(labels[list(combo)])
        if not held_out:
            continue
        n_repeat = sum(1 for i in held_out if labels[i] in sub_labels)
        total += n_repeat / len(held_out)
        count += 1
    return total / count


def _brute_force_presence_weighted(X, m):
    """The DIFFERENT quantity: average over subsets of size m of
    sum_i (X_i/n) * 1(allele i present in the subset). Documented in _coverage.py as the
    quantity coverage_at_size does NOT compute."""
    X = np.asarray(X, dtype=int)
    n = int(X.sum())
    labels = np.repeat(np.arange(len(X)), X)
    idx = np.arange(n)
    total = 0.0
    count = 0
    for combo in itertools.combinations(idx, m):
        sub_labels = set(labels[list(combo)])
        s = sum((X[i] / n) for i in range(len(X)) if i in sub_labels)
        total += s
        count += 1
    return total / count


def test_interpolation_brute_force():
    X = [3, 2, 1, 1]
    n = sum(X)
    for m in (1, 2, 3, 5):
        expected = _brute_force_interp_coverage(X, m)
        got = cov.coverage_at_size(X, m)
        assert close(got, expected, tol=1e-9), (m, got, expected)


def test_interpolation_differs_from_presence_weighted_definition():
    # Confirms the two candidate "rarefied coverage" definitions are genuinely different,
    # and that coverage_at_size matches the one documented (interp_coverage), not the other.
    X = [3, 2, 1, 1]
    m = 1
    interp = _brute_force_interp_coverage(X, m)
    presence = _brute_force_presence_weighted(X, m)
    assert not close(interp, presence, tol=1e-6), "expected these to differ, they didn't"
    assert close(cov.coverage_at_size(X, m), interp)


def test_coverage_at_size_extrapolation_matches_formula():
    X = [1, 1, 1, 2, 2, 10]
    n = sum(X)
    f1, f2 = 3, 2
    A = cov._A_hat(n, f1, f2)
    for m_star in (1, 5, 20):
        expected = 1 - (f1 / n) * A ** (m_star + 1)
        got = cov.coverage_at_size(X, n + m_star)
        assert close(got, expected), (m_star, got, expected)


# ---------------------------------------------------------------------------
# 3. richness_at_size
# ---------------------------------------------------------------------------

def _brute_force_richness(X, m):
    X = np.asarray(X, dtype=int)
    n = int(X.sum())
    labels = np.repeat(np.arange(len(X)), X)
    idx = np.arange(n)
    total = 0.0
    count = 0
    for combo in itertools.combinations(idx, m):
        sub_labels = set(labels[list(combo)])
        total += len(sub_labels)
        count += 1
    return total / count


def test_richness_at_size_brute_force():
    X = [3, 2, 1, 1]
    for m in (1, 2, 3, 6):
        expected = _brute_force_richness(X, m)
        got = cov.richness_at_size(X, m)
        assert close(got, expected, tol=1e-9), (m, got, expected)


def test_richness_at_size_extrapolation_mstar_zero():
    X = [1, 1, 1, 2, 2, 10]
    n = sum(X)
    S_obs = len(X)
    got = cov.richness_at_size(X, n)  # m == n is the interpolation branch, equals S_obs
    assert close(got, S_obs)
    got_ext = cov.richness_at_size(X, n + 1e-9)
    assert close(got_ext, S_obs, tol=1e-6)


# ---------------------------------------------------------------------------
# 4. size_for_coverage
# ---------------------------------------------------------------------------

def test_size_for_coverage_inverse_consistency():
    X = [1, 1, 1, 2, 2, 10, 20, 30]
    C_hat = cov.sample_coverage(X)
    for target in (0.3, 0.6, min(0.95, C_hat) if C_hat > 0.05 else 0.3):
        m, _ = cov.size_for_coverage(X, target)
        assert np.isfinite(m)
        m = int(m)
        assert cov.coverage_at_size(X, m) >= target - 1e-9
        if m > 1:
            assert cov.coverage_at_size(X, m - 1) < target + 1e-9


def test_size_for_coverage_extrapolation_branch():
    X = [1, 1, 1, 2, 2, 10, 20, 30]
    C_hat = cov.sample_coverage(X)
    target = min(0.999, (C_hat + 1.0) / 2)
    if target <= C_hat:
        return  # C_hat already very high in this fixture; nothing to test here
    m, reliable = cov.size_for_coverage(X, target)
    n = sum(X)
    if np.isfinite(m):
        assert m > n
        assert cov.coverage_at_size(X, m) >= target - 1e-6
        assert reliable == (m <= 2 * n)


# ---------------------------------------------------------------------------
# 5. cross_coverage
# ---------------------------------------------------------------------------

def test_cross_coverage_identical_populations_full():
    counts = {"a": 5, "b": 3, "c": 2}
    cc, unseen = cov.cross_coverage(counts, counts, m_ref=sum(counts.values()))
    assert close(cc, 1.0)
    expected_unseen = 1 - cov.sample_coverage(list(counts.values()))
    assert close(unseen, expected_unseen)


def test_cross_coverage_disjoint():
    ref = {"a": 5, "b": 3}
    target = {"x": 4, "y": 6}
    cc, _ = cov.cross_coverage(ref, target)
    assert close(cc, 0.0)


def _brute_force_cross_coverage(counts_ref, counts_target, m_ref):
    allele_ids = sorted(set(counts_ref) | set(counts_target))
    n_A = sum(counts_ref.values())
    labels = []
    for i, a in enumerate(allele_ids):
        labels += [a] * counts_ref.get(a, 0)
    idx = np.arange(n_A)
    n_B = sum(counts_target.values())
    total = 0.0
    count = 0
    for combo in itertools.combinations(idx, m_ref):
        observed = set(labels[i] for i in combo)
        cc = sum(counts_target.get(a, 0) / n_B for a in observed)
        total += cc
        count += 1
    return total / count


def test_cross_coverage_brute_force_partial_subsample():
    counts_ref = {"a": 2, "b": 1, "c": 1}
    counts_target = {"a": 3, "b": 1, "d": 2}
    m_ref = 2
    expected = _brute_force_cross_coverage(counts_ref, counts_target, m_ref)
    got, _ = cov.cross_coverage(counts_ref, counts_target, m_ref=m_ref)
    assert close(got, expected, tol=1e-9), (got, expected)


# ---------------------------------------------------------------------------
# 6. pooled_design_coverage
# ---------------------------------------------------------------------------

def test_pooled_design_single_population_matches_cross_coverage():
    counts = {"a": 5, "b": 3, "c": 2}
    s = 6
    cov_by_pop, global_cov = cov.pooled_design_coverage(
        {"A": counts}, {"A": s}, {"A": 1.0})
    expected, _ = cov.cross_coverage(counts, counts, m_ref=s)
    assert close(cov_by_pop["A"], expected)
    assert close(global_cov, expected)


def test_pooled_design_adding_population_never_decreases_coverage():
    countsA = {"a": 5, "b": 3, "c": 2}
    countsB = {"a": 1, "d": 4, "e": 3}
    cov_A_only, _ = cov.pooled_design_coverage(
        {"A": countsA}, {"A": 5}, {"A": 1.0})
    cov_A_and_B, _ = cov.pooled_design_coverage(
        {"A": countsA, "B": countsB}, {"A": 5, "B": 4}, {"A": 1.0})
    assert cov_A_and_B["A"] >= cov_A_only["A"] - 1e-9


def test_pooled_design_rejects_oversized_subsample():
    counts = {"a": 5, "b": 3}
    try:
        cov.pooled_design_coverage({"A": counts}, {"A": 100}, {"A": 1.0})
        assert False, "expected ValueError"
    except ValueError:
        pass


# ---------------------------------------------------------------------------
# 7. bootstrap_ci
# ---------------------------------------------------------------------------

def test_bootstrap_ci_contains_point_estimate_roughly():
    X = [1, 1, 1, 2, 2, 5, 5, 5, 10, 20]
    point = cov.sample_coverage(X)
    lo, hi = cov.bootstrap_ci(X, cov.sample_coverage, n_boot=300, seed=1)
    assert lo <= hi
    # Not a strict guarantee for every seed/statistic, but sample_coverage's bootstrap
    # distribution should be tight around the point estimate for this fixture.
    assert lo - 0.2 <= point <= hi + 0.2


def test_bootstrap_ci_reproducible_with_seed():
    X = [1, 1, 1, 2, 2, 5, 5, 5, 10, 20]
    r1 = cov.bootstrap_ci(X, cov.sample_coverage, n_boot=50, seed=42)
    r2 = cov.bootstrap_ci(X, cov.sample_coverage, n_boot=50, seed=42)
    assert r1 == r2


if __name__ == "__main__":
    names = [n for n in dir(sys.modules[__name__]) if n.startswith("test_")]
    failures = []
    for name in names:
        fn = globals()[name]
        try:
            fn()
            print(f"  PASS  {name}")
        except AssertionError as e:
            print(f"  FAIL  {name} -- {e}")
            failures.append(name)
    if failures:
        sys.exit(f"{len(failures)} test(s) failed: {failures}")
    print(f"ok ({len(names)} tests)")
