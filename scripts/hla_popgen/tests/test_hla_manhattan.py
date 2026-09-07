#!/usr/bin/env python3
"""Unit tests for 15_hla_manhattan.py.

Focus is the logic that 11_gene_diversity_track.py got wrong or never had:
  1. STRAND handling in the canonical-coordinate cs walk. `11` ignored strand entirely; real PAF
     rows for one gene come back on both strands, so mishandling silently mirrors ~half the
     haplotypes. The reverse-strand test below is the one that would have caught it.
  2. cs operation orientation (query = canonical allele, target = contig), i.e. that a `*xy`
     substitution reports the PERSON's base (x, the contig base), not the canonical base.
  3. Unbiased per-site pi over allele labels, including the "everyone matches" and
     "coverage < 2" edge cases.

Run: python3 scripts/hla_popgen/tests/test_hla_manhattan.py
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)


def _load_module(filename, modname):
    path = os.path.join(HLA_POPGEN_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


m = _load_module("15_hla_manhattan.py", "hla_manhattan")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# cs walk, forward strand
# ---------------------------------------------------------------------------
def test_forward_substitution_position_and_base():
    # ':10*ag:5' -- 10 matches, then contig 'a' vs canonical 'g', then 5 matches.
    obs, ins = m.walk_cs_canonical(":10*ag:5", qstart=0, qend=16, strand="+")
    check("forward: substitution lands at canonical position 10", set(obs) == {10}, str(obs))
    check("forward: records the PERSON's base (contig base 'a'), not the canonical 'g'",
          obs.get(10) == "a", str(obs))
    check("forward: no insertions recorded", ins == [], str(ins))


def test_forward_qstart_offset_applied():
    # Alignment starts at canonical position 100, so the substitution is at 100+10.
    obs, _ = m.walk_cs_canonical(":10*ag:5", qstart=100, qend=116, strand="+")
    check("forward: qstart offsets canonical positions", set(obs) == {110}, str(obs))


def test_deletion_of_canonical_bases_spans_positions():
    # '+cc' = bases present in canonical, absent from contig -> person DELETED them.
    obs, ins = m.walk_cs_canonical(":10+cc:5", qstart=0, qend=17, strand="+")
    check("'+' marks each canonical position the person lacks as DEL",
          obs == {10: "DEL", 11: "DEL"}, str(obs))
    check("'+' is not an insertion event", ins == [], str(ins))


def test_insertion_has_no_canonical_coordinate():
    # '-tt' = bases present in contig, absent from canonical -> no canonical column exists.
    obs, ins = m.walk_cs_canonical(":10-tt:5", qstart=0, qend=15, strand="+")
    check("'-' contributes no per-site allele label", obs == {}, str(obs))
    check("'-' is recorded as an insertion anchored at the current canonical position",
          ins == [10], str(ins))


# ---------------------------------------------------------------------------
# cs walk, reverse strand -- the case 11_gene_diversity_track.py never handled
# ---------------------------------------------------------------------------
def test_reverse_strand_mirrors_position():
    # Same cs string, reverse strand, full-length span [0, 16). The k-th consumed query base maps
    # to qend-1-k, so the substitution at k=10 lands at 16-1-10 = 5, NOT at 10.
    obs, _ = m.walk_cs_canonical(":10*ag:5", qstart=0, qend=16, strand="-")
    check("reverse: substitution maps to qend-1-k (5), not qstart+k (10)",
          set(obs) == {5}, str(obs))


def test_reverse_strand_differs_from_forward():
    fwd, _ = m.walk_cs_canonical(":10*ag:5", qstart=0, qend=16, strand="+")
    rev, _ = m.walk_cs_canonical(":10*ag:5", qstart=0, qend=16, strand="-")
    check("reverse and forward give different canonical positions (regression guard for the "
          "strand bug in 11_gene_diversity_track.py)", set(fwd) != set(rev),
          f"fwd={fwd} rev={rev}")


def test_reverse_strand_deletion_span_runs_backwards():
    obs, _ = m.walk_cs_canonical(":10+cc:5", qstart=0, qend=17, strand="-")
    # k = 10 and 11 -> 17-1-10 = 6 and 17-1-11 = 5
    check("reverse: DEL span walks backwards down the canonical axis",
          obs == {6: "DEL", 5: "DEL"}, str(obs))


# ---------------------------------------------------------------------------
# per-site pi
# ---------------------------------------------------------------------------
def test_pi_all_reference_is_zero():
    pi, n_diff = m.per_site_pi([10, 10], {})
    check("no observed differences -> pi = 0 at every site", pi == [0.0, 0.0], str(pi))
    check("no observed differences -> n_diff = 0", n_diff == [0, 0], str(n_diff))


def test_pi_two_equal_alleles_is_one_half_corrected():
    # n=2 covering, 1 alt: labels {alt:1, REF:1} -> 1 - (0.25+0.25) = 0.5; x n/(n-1)=2 -> 1.0
    pi, n_diff = m.per_site_pi([2], {0: {"a": 1}})
    check("n=2 with a 1:1 split gives unbiased pi = 1.0", abs(pi[0] - 1.0) < 1e-12, str(pi))
    check("n_diff counts the alt-carrying haplotypes", n_diff == [1], str(n_diff))


def test_pi_matches_hand_computed_value():
    # n=4: {a:1, g:1, REF:2} -> 1 - (0.0625+0.0625+0.25) = 0.625; x 4/3 = 0.8333...
    pi, _ = m.per_site_pi([4], {0: {"a": 1, "g": 1}})
    check("n=4 three-allele site matches hand-computed unbiased pi (0.8333)",
          abs(pi[0] - (4 / 3) * 0.625) < 1e-12, str(pi))


def test_pi_requires_at_least_two_covering_haplotypes():
    pi, _ = m.per_site_pi([1], {0: {"a": 1}})
    check("a site covered by a single haplotype yields pi = 0 (undefined, not negative)",
          pi == [0.0], str(pi))


def test_pi_never_negative_even_if_counts_exceed_coverage():
    # Defensive: alt count > coverage must not produce a negative REF count / pi.
    pi, _ = m.per_site_pi([2], {0: {"a": 5}})
    check("alt count exceeding coverage is clamped, pi stays in [0, 1]",
          0.0 <= pi[0] <= 1.0, str(pi))


# ---------------------------------------------------------------------------
# annotation parsing + CDS mask
# ---------------------------------------------------------------------------
def test_parse_ranges():
    check("parse_ranges reads comma-separated 1-based inclusive ranges",
          m.parse_ranges("301..373,504..773") == [(301, 373), (504, 773)],
          str(m.parse_ranges("301..373,504..773")))
    check("parse_ranges tolerates junk", m.parse_ranges("") == [])


def test_cds_mask_is_zero_based_inclusive_of_both_ends():
    # 1-based inclusive (2,3) -> 0-based indices 1 and 2.
    mask = m.in_ranges_mask(5, [(2, 3)])
    check("in_ranges_mask converts 1-based inclusive to 0-based correctly",
          mask == [False, True, True, False, False], str(mask))


def test_sliding_mean_smooths_without_changing_length():
    vals = [0.0, 0.0, 1.0, 0.0, 0.0]
    sm = m.sliding_mean(vals, 3)
    check("sliding_mean preserves length", len(sm) == len(vals), str(sm))
    check("sliding_mean spreads a spike into its neighbours",
          sm[1] > 0 and sm[3] > 0 and sm[2] > 0, str(sm))
    check("sliding_mean with window=1 is a no-op", m.sliding_mean(vals, 1) == vals)


def main():
    print("test_hla_manhattan.py")
    test_forward_substitution_position_and_base()
    test_forward_qstart_offset_applied()
    test_deletion_of_canonical_bases_spans_positions()
    test_insertion_has_no_canonical_coordinate()
    test_reverse_strand_mirrors_position()
    test_reverse_strand_differs_from_forward()
    test_reverse_strand_deletion_span_runs_backwards()
    test_pi_all_reference_is_zero()
    test_pi_two_equal_alleles_is_one_half_corrected()
    test_pi_matches_hand_computed_value()
    test_pi_requires_at_least_two_covering_haplotypes()
    test_pi_never_negative_even_if_counts_exceed_coverage()
    test_parse_ranges()
    test_cds_mask_is_zero_based_inclusive_of_both_ends()
    test_sliding_mean_smooths_without_changing_length()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
