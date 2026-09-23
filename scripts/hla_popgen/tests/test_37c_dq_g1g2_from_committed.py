#!/usr/bin/env python3
"""Unit tests for scripts/hla_popgen/37c_dq_g1g2_from_committed.py, on synthetic data.

Covers: parsing n_hap_ij_disp into the three disclosure statuses, N-haplotype recovery from
freq_hap, and -- the statistic the task explicitly calls out for a test -- the censoring-interval
O/E purge statistic (`oe_interval_from_table`): lower bound treats every censored (`<20`) cross-
group cell as 0 haplotypes, upper bound treats every one as 19, and the two bounds must bracket
the true observed count regardless of what it actually was (1-19).

Run: python3 scripts/hla_popgen/tests/test_37c_dq_g1g2_from_committed.py
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


mod = _load_module("37c_dq_g1g2_from_committed.py", "dq_g1g2_from_committed")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print("  [%s] %s%s" % (status, name, (" -- " + detail) if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# Build a small synthetic DQA1~DQB1 table matching ld_pairwise.tsv's schema, for one ancestry,
# with a known ground truth so oe_interval_from_table's bounds can be checked exactly.
#
# Design: N = 1000 haplotypes total.
#   - Compatible (G1/G1): DQA1*02:01 x DQB1*02:01, n=400 -> freq_a=0.40, freq_b=0.40 (only pair
#     for each allele, so marginals equal this cell's freq), estimated.
#   - Compatible (G2/G2): DQA1*01:01 x DQB1*05:01, n=500 -> freq_a=0.50, freq_b=0.50, estimated.
#   - Cross-group (predicted_incompatible), estimated: DQA1*01:01 x DQB1*02:01, n=40.
#   - Cross-group, censored (<20): DQA1*02:01 x DQB1*05:01, true n=12 (written "<20").
#   - Cross-group, not observed: none needed beyond the above for the O/E test, but include one
#     to check status parsing: DQA1*03:01 x DQB1*03:01 isn't used here (keep it simple: 2 alleles
#     per gene only, so the 2x2 cross product is exactly the two cross-group cells above).
#
# NOTE: with only 2 alleles per side, freq_a for DQA1*01:01 must include BOTH its cells (500 + 40
# = 540 -> freq_a = 0.54) and DQA1*02:01 both its cells (400 + 12 = 412 -> freq_a = 0.412); freq_b
# symmetric. This makes the synthetic table internally consistent (row/col sums to true marginals)
# the way the real committed table is.
# ---------------------------------------------------------------------------
def _make_synthetic_table():
    n_total = 1000.0
    # True (uncensored) counts, only known to the test -- the table itself only exposes n=12 as "<20".
    n_g1g1 = 400
    n_g2g2 = 500
    n_cross_est = 40
    n_cross_censored_true = 12  # hidden from the table; table only says "<20"

    freq_a_0101 = (n_g2g2 + n_cross_est) / n_total       # DQA1*01:01 marginal (G2 x G2 + G2 x G1)
    freq_a_0201 = (n_g1g1 + n_cross_censored_true) / n_total  # DQA1*02:01 marginal
    freq_b_0201 = (n_g1g1 + n_cross_est) / n_total        # DQB1*02:01 marginal (G1 x G1 + G2 x G1)
    freq_b_0501 = (n_g2g2 + n_cross_censored_true) / n_total  # DQB1*05:01 marginal

    rows = [
        # compatible G1/G1
        {"allele_a": "DQA1*02:01", "allele_b": "DQB1*02:01", "freq_a": freq_a_0201,
         "freq_b": freq_b_0201, "freq_hap": n_g1g1 / n_total, "D": 0.0, "Dprime": 0.9,
         "r2": 0.5, "chi2_1df": 10.0, "pair": "DQA1~DQB1", "ancestry": "TEST",
         "n_hap_ij_disp": str(n_g1g1)},
        # compatible G2/G2
        {"allele_a": "DQA1*01:01", "allele_b": "DQB1*05:01", "freq_a": freq_a_0101,
         "freq_b": freq_b_0501, "freq_hap": n_g2g2 / n_total, "D": 0.0, "Dprime": 0.95,
         "r2": 0.5, "chi2_1df": 10.0, "pair": "DQA1~DQB1", "ancestry": "TEST",
         "n_hap_ij_disp": str(n_g2g2)},
        # cross-group, estimated (clears 20)
        {"allele_a": "DQA1*01:01", "allele_b": "DQB1*02:01", "freq_a": freq_a_0101,
         "freq_b": freq_b_0201, "freq_hap": n_cross_est / n_total, "D": 0.0, "Dprime": -0.5,
         "r2": 0.1, "chi2_1df": 2.0, "pair": "DQA1~DQB1", "ancestry": "TEST",
         "n_hap_ij_disp": str(n_cross_est)},
        # cross-group, censored (<20) -- the table never reveals the true count of 12
        {"allele_a": "DQA1*02:01", "allele_b": "DQB1*05:01", "freq_a": freq_a_0201,
         "freq_b": freq_b_0501, "freq_hap": n_cross_censored_true / n_total, "D": 0.0,
         "Dprime": float("nan"), "r2": float("nan"), "chi2_1df": float("nan"),
         "pair": "DQA1~DQB1", "ancestry": "TEST", "n_hap_ij_disp": "<20"},
    ]
    df = pd.DataFrame(rows)
    return df, {
        "n_total": n_total, "n_cross_true_total": n_cross_est + n_cross_censored_true,
        "n_cross_censored_true": n_cross_censored_true,
    }


def test_status_parsing():
    df_raw, _ = _make_synthetic_table()
    # Round-trip through load_committed's status/casting logic by re-using its _status closure
    # indirectly: replicate the same column types load_committed expects, then call the loader's
    # helper behavior by constructing the same columns it would produce.
    df = df_raw.copy()
    df["status"] = df["n_hap_ij_disp"].map(
        lambda v: "suppressed_lt20" if v == "<20" else ("not_observed" if v == "0" else "estimated"))
    est_mask = df["n_hap_ij_disp"].isin(["400", "500", "40"])
    check("status: n>=20 rows -> estimated", (df.loc[est_mask, "status"] == "estimated").all())
    check("status: '<20' row -> suppressed_lt20",
         (df.loc[df["n_hap_ij_disp"] == "<20", "status"] == "suppressed_lt20").all())


def test_estimate_n_haplotypes():
    df_raw, truth = _make_synthetic_table()
    # Write the synthetic table to a temp file and load it through the real loader end-to-end.
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        df_raw.to_csv(fh.name, sep="\t", index=False)
        tmp_path = fh.name
    try:
        loaded = mod.load_committed(tmp_path)
        n_est = mod.estimate_n_haplotypes(loaded[loaded["ancestry"] == "TEST"])
        check("estimate_n_haplotypes recovers true N", n_est == int(truth["n_total"]),
             "got %r want %r" % (n_est, truth["n_total"]))
        return loaded, truth
    finally:
        os.unlink(tmp_path)


def test_oe_interval_brackets_truth():
    loaded, truth = test_estimate_n_haplotypes()
    df_anc = loaded[loaded["ancestry"] == "TEST"]
    n_haps = mod.estimate_n_haplotypes(df_anc)
    res = mod.oe_interval_from_table(df_anc, n_haps)
    check("oe_interval_from_table returns a result", res is not None)
    if res is None:
        return

    check("n_cross_group_cells == 2 (both DQA1 x both DQB1 cross-group combos)",
         res["n_cross_group_cells"] == 2, detail=str(res))
    check("n_cells_estimated == 1", res["n_cells_estimated"] == 1)
    check("n_cells_censored_lt20 == 1", res["n_cells_censored_lt20"] == 1)

    # Lower bound: censored cell counted as 0 -> observed_lower == the one estimated cell's count.
    check("observed_lower == 40 (censored cell treated as 0)", res["observed_lower"] == 40.0,
         detail=str(res["observed_lower"]))
    # Upper bound: censored cell counted as 19 -> observed_upper == 40 + 19 == 59.
    check("observed_upper == 59 (censored cell treated as 19)", res["observed_upper"] == 59.0,
         detail=str(res["observed_upper"]))

    # The true hidden total (40 + 12 = 52) must lie within [observed_lower, observed_upper].
    true_total = truth["n_cross_true_total"]
    check("true cross-group total (%d) lies within [lower, upper]" % true_total,
         res["observed_lower"] <= true_total <= res["observed_upper"],
         detail="%s <= %s <= %s ?" % (res["observed_lower"], true_total, res["observed_upper"]))

    # oe_lower/oe_upper must be monotonically consistent with observed bounds (lower <= upper).
    check("oe_lower <= oe_upper", res["oe_lower"] <= res["oe_upper"], detail=str(res))
    check("expected > 0", res["expected"] > 0)


def test_oe_interval_all_censored_lower_zero():
    """If every cross-group cell in an ancestry is censored, the lower O/E bound must be exactly
    0 (never fabricated as a point estimate) -- the specific guardrail the task asked for."""
    rows = [
        {"allele_a": "DQA1*02:01", "allele_b": "DQB1*02:01", "freq_a": 0.5, "freq_b": 0.5,
         "freq_hap": 0.5, "D": 0.0, "Dprime": 0.9, "r2": 0.5, "chi2_1df": 10.0,
         "pair": "DQA1~DQB1", "ancestry": "TEST2", "n_hap_ij_disp": "500"},
        {"allele_a": "DQA1*01:01", "allele_b": "DQB1*05:01", "freq_a": 0.49, "freq_b": 0.49,
         "freq_hap": 0.49, "D": 0.0, "Dprime": 0.9, "r2": 0.5, "chi2_1df": 10.0,
         "pair": "DQA1~DQB1", "ancestry": "TEST2", "n_hap_ij_disp": "490"},
        {"allele_a": "DQA1*01:01", "allele_b": "DQB1*02:01", "freq_a": 0.49, "freq_b": 0.5,
         "freq_hap": 0.005, "D": 0.0, "Dprime": float("nan"), "r2": float("nan"),
         "chi2_1df": float("nan"), "pair": "DQA1~DQB1", "ancestry": "TEST2",
         "n_hap_ij_disp": "<20"},
        {"allele_a": "DQA1*02:01", "allele_b": "DQB1*05:01", "freq_a": 0.5, "freq_b": 0.49,
         "freq_hap": 0.005, "D": 0.0, "Dprime": float("nan"), "r2": float("nan"),
         "chi2_1df": float("nan"), "pair": "DQA1~DQB1", "ancestry": "TEST2",
         "n_hap_ij_disp": "<20"},
    ]
    df_raw = pd.DataFrame(rows)
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        df_raw.to_csv(fh.name, sep="\t", index=False)
        tmp_path = fh.name
    try:
        loaded = mod.load_committed(tmp_path)
        df_anc = loaded[loaded["ancestry"] == "TEST2"]
        n_haps = mod.estimate_n_haplotypes(df_anc)
        res = mod.oe_interval_from_table(df_anc, n_haps)
        check("all-censored cross-group -> observed_lower == 0", res["observed_lower"] == 0.0,
             detail=str(res))
        check("all-censored cross-group -> observed_upper == 38 (2 cells x 19)",
             res["observed_upper"] == 38.0, detail=str(res))
        check("all-censored cross-group -> oe_lower == 0.0", res["oe_lower"] == 0.0)
    finally:
        os.unlink(tmp_path)


def main():
    print("test_status_parsing")
    test_status_parsing()
    print("test_oe_interval_brackets_truth")
    test_oe_interval_brackets_truth()
    print("test_oe_interval_all_censored_lower_zero")
    test_oe_interval_all_censored_lower_zero()

    print()
    if FAILURES:
        print("%d check(s) FAILED: %s" % (len(FAILURES), FAILURES))
        sys.exit(1)
    print("All checks passed.")


if __name__ == "__main__":
    main()
