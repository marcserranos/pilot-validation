#!/usr/bin/env python3
"""Unit tests for 31_aa_diversity_selection.py.

What can go silently wrong here, and is therefore what is tested:

  1. **Alignment frame.** If a protein with a deletion is compared position-by-position without
     alignment, every residue after the deletion is compared to the wrong one and the whole
     downstream half of the gene looks hyper-diverse. The result would look like a dramatic
     finding, not like a bug.
  2. **Gaps counted as a 21st amino acid.** That inflates pi exactly where alignment is least
     certain. Gaps must be excluded from the denominator instead.
  3. **The unbiased correction.** pi must be n/(n-1) * (1 - sum p^2), not the raw 1 - sum p^2.
  4. **Split codons at exon boundaries.** HLA exon junctions fall mid-codon. A codon must go to
     the exon contributing most of its bases, not be dropped and not default to the later exon.
  5. **Fst vs G'st.** The whole reason G'st is reported is that raw Fst is bounded near zero when
     within-population heterozygosity is high. A test asserts that bound is real, so nobody
     "simplifies" the code by deleting G'st later.

Run: python3 scripts/hla_popgen/tests/test_aa_diversity_selection.py
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


m = _load_module("31_aa_diversity_selection.py", "aa_div")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# alignment
# ---------------------------------------------------------------------------
def test_equal_length_is_identity():
    check("equal-length proteins map position-for-position",
          m.align_to_reference("ACDEF", "AGDEF") == "AGDEF")


def test_deletion_is_placed_as_a_gap_not_a_frame_shift():
    """THE test. The observed protein is missing the reference's 3rd residue. Everything after it
    must still line up with the reference; a naive comparison would shift by one and report four
    differences instead of one."""
    ref = "ACDEFGH"
    obs = "ACEFGH"          # 'D' deleted
    al = m.align_to_reference(ref, obs)
    check("the aligned string has one column per reference residue", len(al) == len(ref), al)
    check("the deletion appears as a gap at the right position", al[2] == "-", al)
    check("residues after the deletion stay in frame", al[3:] == "EFGH", al)
    diffs = sum(1 for a, b in zip(ref, al) if a != b and b != "-")
    check("no spurious downstream differences are created", diffs == 0, f"{al} vs {ref}")


def test_insertion_has_no_reference_column():
    ref = "ACDEF"
    obs = "ACWDEF"          # 'W' inserted after C
    al = m.align_to_reference(ref, obs)
    check("an inserted residue does not create a reference column", len(al) == len(ref), al)
    check("the reference residues still match", al == ref, f"{al} vs {ref}")


# ---------------------------------------------------------------------------
# pi
# ---------------------------------------------------------------------------
def test_pi_is_zero_for_a_monomorphic_position():
    pi, cov = m.pi_per_codon([("AA", 10), ("AA", 30)], 2)
    check("a monomorphic position has pi = 0", list(pi) == [0.0, 0.0], str(pi))
    check("coverage is the summed weight", list(cov) == [40.0, 40.0], str(cov))


def test_pi_matches_hand_computed_unbiased_value():
    # 50/50 split over 2 amino acids, n=100: 1 - (0.25+0.25) = 0.5; x 100/99 = 0.505050...
    pi, _ = m.pi_per_codon([("A", 50), ("C", 50)], 1)
    check("pi applies the n/(n-1) unbiased correction",
          abs(pi[0] - (100 / 99) * 0.5) < 1e-12, str(pi[0]))


def test_gaps_are_excluded_from_the_denominator_not_counted_as_a_state():
    """Two real states plus a gap. If the gap were a third state, pi would rise; it must instead
    shrink the denominator, leaving the two real states at a 50/50 split."""
    with_gap, cov = m.pi_per_codon([("A", 50), ("C", 50), ("-", 50)], 1)
    without, _ = m.pi_per_codon([("A", 50), ("C", 50)], 1)
    check("a gap does not add diversity", abs(with_gap[0] - without[0]) < 1e-12,
          f"{with_gap[0]} vs {without[0]}")
    check("a gap is excluded from coverage", cov[0] == 100.0, str(cov[0]))


def test_unknown_residue_is_also_excluded():
    pi, cov = m.pi_per_codon([("A", 50), ("C", 50), ("X", 50)], 1)
    check("'X' is treated as missing, not as a state", cov[0] == 100.0, str(cov[0]))


def test_single_observation_position_has_undefined_pi_reported_as_zero():
    pi, _ = m.pi_per_codon([("A", 1)], 1)
    check("n < 2 gives pi = 0 rather than a negative or a crash", pi[0] == 0.0, str(pi[0]))


# ---------------------------------------------------------------------------
# exon mapping
# ---------------------------------------------------------------------------
def test_codon_to_exon_simple():
    # exon 1 = 6bp (2 codons), exon 2 = 9bp (3 codons)
    check("codons map to their exon when boundaries are codon-aligned",
          m.codon_to_exon([(1, 6), (2, 9)]) == [1, 1, 2, 2, 2],
          str(m.codon_to_exon([(1, 6), (2, 9)])))


def test_split_codon_goes_to_the_majority_exon():
    """Exon 1 contributes 7 bases, so codon 3 is 1 base of exon 1 and 2 of exon 2 -> exon 2.
    Dropping split codons instead would shorten the groove by a residue at each junction."""
    out = m.codon_to_exon([(1, 7), (2, 8)])
    check("a codon split 1:2 across a boundary goes to the exon with 2 bases",
          out == [1, 1, 2, 2, 2], str(out))
    out2 = m.codon_to_exon([(1, 8), (2, 7)])
    check("a codon split 2:1 goes to the exon with 2 bases", out2[2] == 1, str(out2))


def test_no_codon_is_lost_at_a_boundary():
    out = m.codon_to_exon([(1, 7), (2, 8)])
    check("total codons equal floor(total_bases / 3)", len(out) == (7 + 8) // 3, str(out))


# ---------------------------------------------------------------------------
# permutation test
# ---------------------------------------------------------------------------
def test_permutation_detects_real_enrichment():
    rng = np.random.default_rng(0)
    pi = np.concatenate([np.full(20, 0.8), np.zeros(80)])
    mask = np.zeros(100, dtype=bool)
    mask[:20] = True
    diff, p, ratio = m.permutation_test(pi, mask, 500, rng)
    check("a genuinely enriched region gives a small p", p < 0.01, str(p))
    check("the ratio is reported", ratio > 10, str(ratio))


def test_permutation_finds_nothing_when_there_is_nothing():
    rng = np.random.default_rng(0)
    pi = rng.random(200)
    mask = np.zeros(200, dtype=bool)
    mask[:50] = True
    _, p, _ = m.permutation_test(pi, mask, 500, rng)
    check("a randomly placed region does not give a significant p", p > 0.05, str(p))


def test_permutation_handles_an_all_or_nothing_mask():
    rng = np.random.default_rng(0)
    pi = np.ones(10)
    _, p, _ = m.permutation_test(pi, np.ones(10, dtype=bool), 100, rng)
    check("a mask covering everything returns NaN rather than crashing", np.isnan(p), str(p))


# ---------------------------------------------------------------------------
# differentiation
# ---------------------------------------------------------------------------
def test_fst_zero_when_populations_are_identical():
    f = m.allele_fst({"AFR": {"a": 0.5, "b": 0.5}, "EUR": {"a": 0.5, "b": 0.5}},
                     {"AFR": 100, "EUR": 100})
    check("identical frequencies give Fst = 0", abs(f["fst"]) < 1e-12, str(f))


def test_fst_is_bounded_near_zero_at_high_heterozygosity():
    """The Brandt et al. 2018 point, as a test. Two populations sharing NO alleles at all -- total
    differentiation -- but each highly polymorphic. Raw Fst stays small; G'st goes to 1. If someone
    later deletes G'st as redundant, this fails."""
    afr = {"a%d" % i: 0.05 for i in range(20)}
    eur = {"e%d" % i: 0.05 for i in range(20)}
    f = m.allele_fst({"AFR": afr, "EUR": eur}, {"AFR": 1000, "EUR": 1000})
    check("completely disjoint allele pools still give a small raw Fst",
          f["fst"] < 0.1, str(f["fst"]))
    check("...while Hedrick's G'st correctly reports near-complete differentiation",
          f["gst_prime"] > 0.9, str(f["gst_prime"]))
    check("within-population heterozygosity is what causes the bound",
          f["hs"] > 0.9, str(f["hs"]))


def test_fst_needs_two_populations():
    f = m.allele_fst({"AFR": {"a": 1.0}}, {"AFR": 100})
    check("a single population yields NaN, not a number", np.isnan(f["fst"]), str(f))


# ---------------------------------------------------------------------------
# misc
# ---------------------------------------------------------------------------
def test_two_field_and_suppression():
    check("two_field truncates to the protein", m.two_field("HLA-DRB1*15:01:01:02") == "DRB1*15:01")
    check("a field-2-novel call is unresolved", m.two_field("HLA-DRB1*15:new") is None)
    check("1-19 is suppressed", m.suppress(5) == "<20")
    check("0 is disclosable", m.suppress(0) == "0")


def main():
    print("test_aa_diversity_selection.py")
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
