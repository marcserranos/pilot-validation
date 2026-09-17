#!/usr/bin/env python3
"""Unit tests for 33_figure1_v3_compose.py.

These exist because of a bug this script actually shipped with for one run:
`pd.to_numeric(...).fillna(0)` turned every AoU-suppressed `<20` cell into a zero. 3,675 of the
8,784 cells in the source table are suppressed, so the result was a Figure 1 panel reporting a
0.00% artifact rate for Middle Eastern ancestry across all eight genes -- a flat, clean-looking,
entirely fabricated row. Nothing about the plot looked wrong.

So the tests pin the censoring contract:
  1. `<20` contributes [0, 19], never 0 and never 20.
  2. A cell that is entirely censored yields a lower bound of 0 AND an upper bound that is
     visibly non-zero, so the interval shows the ignorance instead of hiding it.
  3. The rate's lower bound pairs the smallest numerator with the largest denominator (and vice
     versa), so the interval really brackets the truth rather than being two independent guesses.
  4. The SFS bin order is the numeric one, not the lexicographic one that would put 129-256
     before 2.

Run: python3 scripts/hla_popgen/tests/test_figure1_v3_compose.py
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


m = _load_module("33_figure1_v3_compose.py", "fig1v3")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


def _counts(rows):
    return pd.DataFrame(rows, columns=["gene", "gene_class", "ancestry_scheme", "ancestry",
                                       "unrelated_only", "field_class", "artifact_label",
                                       "n_haplotypes"])


# ---------------------------------------------------------------------------
# censoring
# ---------------------------------------------------------------------------
def test_bounds_of_a_censored_cell():
    lo, hi, n = m._bounds(["<20"])
    check("a censored cell's lower bound is 0", lo == 0.0, str(lo))
    check("...its upper bound is 19, not 20", hi == 19.0, str(hi))
    check("...and it is counted as censored", n == 1, str(n))


def test_bounds_of_plain_counts():
    lo, hi, n = m._bounds(["10", "5"])
    check("plain counts give an exact interval", (lo, hi) == (15.0, 15.0), f"{lo},{hi}")
    check("no censored cells are reported", n == 0, str(n))


def test_bounds_mixed():
    lo, hi, n = m._bounds(["100", "<20", "<20"])
    check("mixed cells give lower = disclosable sum", lo == 100.0, str(lo))
    check("...and upper = disclosable + 19 per censored cell", hi == 138.0, str(hi))


def test_bounds_ignores_junk_without_crashing():
    lo, hi, n = m._bounds(["NA", "", "7"])
    check("unparseable cells are skipped, not counted as zero-with-certainty",
          (lo, hi, n) == (7.0, 7.0, 0), f"{lo},{hi},{n}")


def test_fully_censored_numerator_produces_a_visible_interval():
    """THE regression. Every artifact cell suppressed -> the reported rate must be an interval
    whose upper bound is clearly non-zero, not a confident 0.00%."""
    counts = _counts([
        ("HLA-B", "classical_I", "strict", "MID", "True", "known", "clean", "520"),
        ("HLA-B", "classical_I", "strict", "MID", "True", "f2_protein", "clean", "<20"),
        ("HLA-B", "classical_I", "strict", "MID", "True", "f2_protein", "partial_cds", "<20"),
        ("HLA-B", "classical_I", "strict", "MID", "True", "f2_protein", "homopolymer_indel", "<20"),
    ])
    r = m.rate_by_gene_ancestry(counts, "strict", True, ["HLA-B"])
    art = r[r["field_class"] == "artifact_control"].iloc[0]
    check("the artifact lower bound is 0 (all its cells are suppressed)",
          art["pct_lower"] == 0.0, str(art["pct_lower"]))
    check("...but the upper bound is visibly non-zero, so the ignorance is shown",
          art["pct_upper"] > 1.0, str(art["pct_upper"]))
    check("...and the number of censored cells is recorded",
          art["n_censored_cells"] == 2, str(art["n_censored_cells"]))


def test_interval_brackets_correctly():
    """Lower bound = smallest numerator over largest denominator; upper = the reverse. Pairing
    them the other way round would produce an interval that does not contain the truth."""
    counts = _counts([
        ("HLA-A", "classical_I", "strict", "AFR", "True", "known", "clean", "900"),
        ("HLA-A", "classical_I", "strict", "AFR", "True", "f2_protein", "clean", "<20"),
    ])
    r = m.rate_by_gene_ancestry(counts, "strict", True, ["HLA-A"])
    row = r[r["field_class"] == "f2_protein"].iloc[0]
    # numerator in [0,19]; denominator in [900, 919]
    check("lower bound uses min numerator / max denominator", row["pct_lower"] == 0.0,
          str(row["pct_lower"]))
    expected_hi = 100.0 * 19.0 / 900.0
    check("upper bound uses max numerator / min denominator",
          abs(row["pct_upper"] - expected_hi) < 1e-9,
          f"{row['pct_upper']} vs {expected_hi}")


def test_uncalled_rows_are_excluded_from_the_denominator():
    """A gene that simply fails to be called more often in one ancestry must not look like it has
    lower novelty there."""
    counts = _counts([
        ("HLA-A", "classical_I", "strict", "EUR", "True", "known", "clean", "100"),
        ("HLA-A", "classical_I", "strict", "EUR", "True", "f2_protein", "clean", "100"),
        ("HLA-A", "classical_I", "strict", "EUR", "True", "uncalled", "clean", "800"),
    ])
    r = m.rate_by_gene_ancestry(counts, "strict", True, ["HLA-A"])
    row = r[r["field_class"] == "f2_protein"].iloc[0]
    check("the rate is 50% (100 of 200 called), not 10% (100 of 1000 rows)",
          abs(row["pct_lower"] - 50.0) < 1e-9, str(row["pct_lower"]))


def test_ancestry_scheme_filter_is_respected():
    counts = _counts([
        ("HLA-A", "classical_I", "pred", "EUR", "True", "known", "clean", "100"),
        ("HLA-A", "classical_I", "pred", "EUR", "True", "f2_protein", "clean", "100"),
    ])
    r = m.rate_by_gene_ancestry(counts, "strict", True, ["HLA-A"])
    check("asking for strict ancestry does not silently return the pred rows", r.empty,
          str(len(r)))


def test_unrelated_filter_is_respected():
    counts = _counts([
        ("HLA-A", "classical_I", "strict", "EUR", "False", "known", "clean", "100"),
    ])
    r = m.rate_by_gene_ancestry(counts, "strict", True, ["HLA-A"])
    check("related-inclusive rows are excluded when unrelated_only is requested", r.empty,
          str(len(r)))


# ---------------------------------------------------------------------------
# SFS
# ---------------------------------------------------------------------------
def test_sfs_bins_are_in_numeric_order():
    idx = m.SFS_BIN_ORDER
    check("bin '2' comes before bin '129-256'", idx.index("2") < idx.index("129-256"))
    check("the singleton bin is first", idx[0] == "1")
    check("the overflow bin is last", idx[-1] == ">1024")


def test_sfs_density_divides_by_bin_width():
    sfs = pd.DataFrame([
        {"gene": "ALL_CLASSICAL", "series": "known", "occurrence_bin": "1", "n_alleles": 10},
        {"gene": "ALL_CLASSICAL", "series": "known", "occurrence_bin": "3-4", "n_alleles": 10},
    ])
    counts, dens = m.sfs_known_vs_novel(sfs, "ALL_CLASSICAL")
    check("a width-1 bin keeps its count", dens.loc["1", "known"] == 10.0,
          str(dens.loc["1", "known"]))
    check("a width-2 bin is halved, removing the sawtooth",
          dens.loc["3-4", "known"] == 5.0, str(dens.loc["3-4", "known"]))


def test_sfs_flagged_artifacts_are_not_counted_as_novel():
    sfs = pd.DataFrame([
        {"gene": "ALL_CLASSICAL", "series": "flagged_artifact", "occurrence_bin": "1",
         "n_alleles": 99},
        {"gene": "ALL_CLASSICAL", "series": "cds_changing", "occurrence_bin": "1", "n_alleles": 1},
    ])
    counts, _ = m.sfs_known_vs_novel(sfs, "ALL_CLASSICAL")
    check("flagged artifacts form their own series", counts.loc["1", "flagged"] == 99.0,
          str(counts.loc["1"].to_dict()))
    check("...and are not folded into 'novel'", counts.loc["1", "novel"] == 1.0,
          str(counts.loc["1"].to_dict()))


def main():
    print("test_figure1_v3_compose.py")
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
