#!/usr/bin/env python3
"""Unit tests for 11_gene_diversity_track.py's cs-string walker and target<->query boundary
translation -- the two pieces of new logic this script introduces (PAF/GTF join was validated live
against real VM data instead; see NEEDLE_VIEW_BRIEF.md's "Live VM validation" section).

Run: python3 scripts/hla_popgen/tests/test_gene_diversity_track.py
"""
import importlib.util
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)


def _load_module(filename, modname):
    path = os.path.join(HLA_POPGEN_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


m = _load_module("11_gene_diversity_track.py", "gene_diversity_track")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


def test_walk_cs_pure_match():
    variants, bp = m.walk_cs(":50")
    check("pure match: no variants", variants == [], detail=str(variants))
    check("pure match: breakpoints span full length",
          bp[0] == (0, 0) and bp[-1] == (50, 50), detail=str(bp))


def test_walk_cs_substitution():
    variants, _ = m.walk_cs(":10*ag:10")
    check("substitution at query offset 10", variants == [10], detail=str(variants))


def test_walk_cs_insertion_query_extra():
    # '+cc' = query (reference) has 2 extra bases the contig doesn't -- query advances, target
    # doesn't. Every one of those reference positions is disturbed.
    variants, bp = m.walk_cs(":10+cc:10")
    check("insertion touches both query positions it spans",
          variants == [10, 11], detail=str(variants))
    check("insertion advances query but not target",
          bp[2] == (10, 12), detail=str(bp))


def test_walk_cs_deletion_target_extra():
    # '-tt' = target (contig) has 2 extra bases the reference doesn't -- target advances, query
    # doesn't. One tick recorded at the anchor query offset.
    variants, bp = m.walk_cs(":10-tt:10")
    check("deletion anchors one variant tick at query offset 10",
          variants == [10], detail=str(variants))
    check("deletion advances target but not query",
          bp[2] == (12, 10), detail=str(bp))


def test_walk_cs_combined_matches_manual_trace():
    # ':10*ag:5+cc:5-tt:10' traced by hand in NEEDLE_VIEW_BRIEF.md-adjacent scratch work:
    # sub at q=10; insertion (2bp) touches q=16,17; deletion (2bp) anchors at q=23.
    variants, bp = m.walk_cs(":10*ag:5+cc:5-tt:10")
    check("combined cs string variant positions", variants == [10, 16, 17, 23],
          detail=str(variants))
    check("combined cs string final breakpoint (target=33, query=33)", bp[-1] == (33, 33),
          detail=str(bp))


def test_target_to_query_within_match_run():
    _, bp = m.walk_cs(":10*ag:10")  # match(10) sub match(10): t=0..21, q=0..21 in lockstep after
    # the substitution too (1-for-1), so target and query offsets stay equal throughout.
    check("midpoint of a pure-match run maps 1:1", m.target_to_query(bp, 5) == 5)
    check("offset just after the substitution maps 1:1", m.target_to_query(bp, 15) == 15)


def test_target_to_query_across_deletion():
    # Target has 2 extra bases (a deletion from the person's view) between t=10 and t=12; query
    # stays at 10 throughout that span. A target offset landing inside the deleted span should
    # anchor to the query offset at the start of the span (best a breakpoint map can do).
    _, bp = m.walk_cs(":10-tt:10")
    check("target offset at the deletion's start maps to query 10",
          m.target_to_query(bp, 10) == 10, detail=str(bp))


def test_clip_ranges_and_range_bp():
    ranges = [(0, 100), (200, 300)]
    clipped = m.clip_ranges(ranges, 50, 250)
    check("clip drops/truncates ranges outside [lo, hi]",
          clipped == [(50, 100), (200, 250)], detail=str(clipped))
    check("range_bp sums correctly", m.range_bp(clipped) == 100, detail=str(m.range_bp(clipped)))
    check("clip that excludes everything returns empty",
          m.clip_ranges(ranges, 500, 600) == [])


def test_build_manhattan_track_stratifies_correctly():
    # Two haplotypes on the SAME template (directly comparable positions), one on a DIFFERENT
    # template (must be excluded from this call -- caller is responsible for pre-filtering to the
    # dominant template only; this test locks in that build_manhattan_track itself just aggregates
    # whatever records it's given, treating them as already-comparable).
    records = [
        ({10, 20}, [(5, 15)], 100),
        ({20, 30}, [(5, 15)], 100),
    ]
    counts, cds_ranges, qlen = m.build_manhattan_track(records)
    check("position 20 seen in both haplotypes", counts[20] == 2, detail=str(counts))
    check("position 10 seen in one haplotype", counts[10] == 1, detail=str(counts))
    check("cds_ranges taken from first record", cds_ranges == [(5, 15)])
    check("qlen taken from first record", qlen == 100)


def test_build_manhattan_track_empty():
    counts, cds_ranges, qlen = m.build_manhattan_track([])
    check("empty input -> empty counts, no crash", counts == Counter() and cds_ranges == []
          and qlen is None)


def test_consensus_cds_ranges_basic():
    ranges = [[(100, 200), (300, 400)], [(101, 199), (299, 401)], [(99, 201), (301, 399)]]
    consensus, spreads = m.consensus_cds_ranges(ranges)
    check("consensus has 2 segments (matches all 3 inputs)", len(consensus) == 2,
          detail=str(consensus))
    check("consensus segment 1 is the median (100, 200)", consensus[0] == (100, 200),
          detail=str(consensus))
    check("spread is small for tightly-clustered synthetic boundaries",
          all(s[0] <= 2 and s[1] <= 2 for s in spreads), detail=str(spreads))


def test_consensus_cds_ranges_empty():
    consensus, spreads = m.consensus_cds_ranges([])
    check("empty input -> empty output, no crash", consensus == [] and spreads == [])


def main():
    print("test_gene_diversity_track.py")
    test_walk_cs_pure_match()
    test_walk_cs_substitution()
    test_walk_cs_insertion_query_extra()
    test_walk_cs_deletion_target_extra()
    test_walk_cs_combined_matches_manual_trace()
    test_target_to_query_within_match_run()
    test_target_to_query_across_deletion()
    test_clip_ranges_and_range_bp()
    test_build_manhattan_track_stratifies_correctly()
    test_build_manhattan_track_empty()
    test_consensus_cds_ranges_basic()
    test_consensus_cds_ranges_empty()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
