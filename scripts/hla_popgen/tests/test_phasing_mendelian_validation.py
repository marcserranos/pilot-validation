#!/usr/bin/env python3
"""Fixture tests for 12_phasing_mendelian_validation.py -- constructs a tiny synthetic cohort with
KNOWN ground truth (a clean parent-child transmission, an injected single-gene mismatch, a
simulated crossover, and an IBD0-style unrelated pair) and asserts the pipeline recovers exactly
the expected sharing fractions, mismatch locations, and switch counts. This is the correctness gate
for a script whose whole point is measuring a real error rate -- a silent bug here would misreport
the actual accuracy of the long-read pipeline, so every core function is checked against hand-
computed expectations before this is trusted against real data."""
import importlib
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
mod = importlib.import_module("12_phasing_mendelian_validation")

GENES = [f"G{i}" for i in range(1, 11)]
GENE_CLASS = {**{f"G{i}": "classical_I" for i in range(1, 6)},
              **{f"G{i}": "pseudogene_I" for i in range(6, 11)}}


def _row(person_id, hap, gene, allele, gene_start, contig="ctg1", copy_index=1):
    return {"person_id": person_id, "hap": hap, "contig": contig, "gene": f"HLA-{gene}",
            "gene_bare": gene, "gene_class": GENE_CLASS[gene], "consensus": allele,
            "copy_index": copy_index, "gene_start": gene_start, "gene_end": gene_start + 50}


def _make_table1():
    rows = []
    # P1: the shared "parent" -- hap1 = A-series, hap2 = B-series
    for i, g in enumerate(GENES, start=1):
        rows.append(_row("P1", "hap1", g, f"A{i}:01", i * 100))
        rows.append(_row("P1", "hap2", g, f"B{i}:01", i * 100))
    # P2: clean child -- hap1 is an exact copy of P1's hap1 all the way through; hap2 unrelated.
    for i, g in enumerate(GENES, start=1):
        rows.append(_row("P2", "hap1", g, f"A{i}:01", i * 100))
        rows.append(_row("P2", "hap2", g, f"Y{i}:01", i * 100))
    # P3: child with ONE injected mismatch at G6 (a pseudogene) -- otherwise a clean copy of hap1.
    for i, g in enumerate(GENES, start=1):
        allele = "WRONG:99" if g == "G6" else f"A{i}:01"
        rows.append(_row("P3", "hap1", g, allele, i * 100))
        rows.append(_row("P3", "hap2", g, f"Z{i}:01", i * 100))
    # P4: simulated crossover -- G1-G5 copy P1's hap1, G6-G10 copy P1's hap2 (one switch expected).
    for i, g in enumerate(GENES, start=1):
        allele = f"A{i}:01" if i <= 5 else f"B{i}:01"
        rows.append(_row("P4", "hap1", g, allele, i * 100))
        rows.append(_row("P4", "hap2", g, f"X{i}:01", i * 100))
    # P5: unrelated (IBD0-style) -- shares nothing with P1.
    for i, g in enumerate(GENES, start=1):
        rows.append(_row("P5", "hap1", g, f"Q{i}:01", i * 100))
        rows.append(_row("P5", "hap2", g, f"R{i}:01", i * 100))
    # P6: the SAME real transmitted chromosome as P2 (no true recombination) for the whole region,
    # but assembled as two separate contigs with the hap1/hap2 LABEL flipped on the second contig --
    # simulates the exact artifact a real run of this script caught 2026-09-08 (arbitrary hap-label
    # reshuffling between independently-assembled fragments). G1-G5 on ctgA with the transmitted
    # allele on hap1 (matches P1's hap1, like P2); G6-G10 on ctgB with the SAME transmitted allele
    # now on hap2 instead. Zero real switches should be reported once contig is respected.
    for i, g in enumerate(GENES, start=1):
        if i <= 5:
            rows.append(_row("P6", "hap1", g, f"A{i}:01", i * 100, contig="ctgA"))
            rows.append(_row("P6", "hap2", g, f"Y{i}:01", i * 100, contig="ctgA"))
        else:
            rows.append(_row("P6", "hap2", g, f"A{i}:01", i * 100, contig="ctgB"))
            rows.append(_row("P6", "hap1", g, f"Y{i}:01", i * 100, contig="ctgB"))
    return pd.DataFrame(rows)


def test_build_allele_lookup_and_canonical_order():
    table1 = _make_table1()
    lookup, gene_class, contig_lookup = mod.build_allele_lookup(table1, resolution=2)
    assert lookup[("P1", "G1", "hap1")] == "A1:01"
    assert lookup[("P1", "G1", "hap2")] == "B1:01"
    assert gene_class["G6"] == "pseudogene_I"

    order = mod.derive_canonical_order(table1)
    ordered_genes = sorted(order, key=order.get)
    assert ordered_genes == GENES  # must recover the true gene_start order exactly


def test_clean_parent_child_pair_zero_mismatch_zero_switches():
    table1 = _make_table1()
    lookup, _, contig_lookup = mod.build_allele_lookup(table1, resolution=2)
    order = mod.derive_canonical_order(table1)
    result = mod.analyze_pair("P1", "P2", GENES, lookup, order, contig_lookup)
    assert result["n_compared"] == 10
    assert result["n_shared"] == 10
    assert result["sharing_frac"] == 1.0
    n_switches, n_states = mod.count_switches(result["detail"])
    assert n_switches == 0
    assert n_states == 10


def test_injected_mismatch_is_localized_to_the_right_gene():
    table1 = _make_table1()
    lookup, _, contig_lookup = mod.build_allele_lookup(table1, resolution=2)
    order = mod.derive_canonical_order(table1)
    result = mod.analyze_pair("P1", "P3", GENES, lookup, order, contig_lookup)
    assert result["n_compared"] == 10
    assert result["n_shared"] == 9          # every gene except G6
    assert result["sharing_frac"] == 0.9
    mismatched_genes = [d["gene"] for d in result["detail"] if not d["shared"]]
    assert mismatched_genes == ["G6"]


def test_simulated_crossover_gives_exactly_one_switch():
    table1 = _make_table1()
    lookup, _, contig_lookup = mod.build_allele_lookup(table1, resolution=2)
    order = mod.derive_canonical_order(table1)
    result = mod.analyze_pair("P1", "P4", GENES, lookup, order, contig_lookup)
    assert result["sharing_frac"] == 1.0    # every gene shares SOME allele with P1
    n_switches, n_states = mod.count_switches(result["detail"])
    assert n_switches == 1
    assert n_states == 10


def test_unrelated_pair_is_low_sharing():
    table1 = _make_table1()
    lookup, _, contig_lookup = mod.build_allele_lookup(table1, resolution=2)
    order = mod.derive_canonical_order(table1)
    result = mod.analyze_pair("P1", "P5", GENES, lookup, order, contig_lookup)
    assert result["sharing_frac"] == 0.0
    cls = mod.classify_pair(result["sharing_frac"], result["n_compared"],
                             min_genes=5, high_thresh=0.8, low_thresh=0.3)
    assert cls == "low_sharing"


def test_contig_boundary_suppresses_spurious_switch():
    """The exact artifact caught in the first real VM run (2026-09-08): P6 carries the SAME real
    transmitted chromosome as P1's hap1 for the whole region, but the hap1/hap2 LABEL flips between
    two independently-assembled contigs (ctgA for G1-G5, ctgB for G6-G10). Without the contig guard,
    this reads as a switch at the G5->G6 boundary; with it, that transition must be suppressed
    (real switches only count within one contig for both people)."""
    table1 = _make_table1()
    lookup, _, contig_lookup = mod.build_allele_lookup(table1, resolution=2)
    order = mod.derive_canonical_order(table1)
    result = mod.analyze_pair("P1", "P6", GENES, lookup, order, contig_lookup)
    assert result["sharing_frac"] == 1.0
    n_switches, n_states = mod.count_switches(result["detail"])
    assert n_switches == 0
    assert n_states == 10


def test_classify_pair_thresholds():
    assert mod.classify_pair(0.9, 10, 5, 0.8, 0.3) == "high_sharing"
    assert mod.classify_pair(0.2, 10, 5, 0.8, 0.3) == "low_sharing"
    assert mod.classify_pair(0.5, 10, 5, 0.8, 0.3) == "mixed"
    assert mod.classify_pair(0.9, 2, 5, 0.8, 0.3) == "insufficient_data"


def test_run_all_pairs_aggregates_mismatch_to_the_correct_gene_and_class():
    table1 = _make_table1()
    lookup, gene_class, contig_lookup = mod.build_allele_lookup(table1, resolution=2)
    order = mod.derive_canonical_order(table1)
    pairs_df = pd.DataFrame([
        {"person_i": "P1", "person_j": "P2", "kin": 0.25,
         "kin_degree": "first_degree_parentchild_or_sibling"},
        {"person_i": "P1", "person_j": "P3", "kin": 0.25,
         "kin_degree": "first_degree_parentchild_or_sibling"},
        {"person_i": "P1", "person_j": "P4", "kin": 0.25,
         "kin_degree": "first_degree_parentchild_or_sibling"},
        {"person_i": "P1", "person_j": "P5", "kin": 0.25,
         "kin_degree": "first_degree_parentchild_or_sibling"},
    ])
    pairs_out, gene_out, switches_out, mixed_examples = mod.run_all_pairs(
        pairs_df, lookup, order, GENES, min_genes=5, high_thresh=0.8, low_thresh=0.3,
        contig_lookup=contig_lookup)

    classes = dict(zip(zip(pairs_out["person_i"], pairs_out["person_j"]), pairs_out["pair_class"]))
    assert classes[("P1", "P2")] == "high_sharing"
    assert classes[("P1", "P3")] == "high_sharing"
    assert classes[("P1", "P4")] == "high_sharing"
    assert classes[("P1", "P5")] == "low_sharing"

    g6 = gene_out[gene_out["gene"] == "G6"].iloc[0]
    # P2 and P4 both match cleanly at G6; only P3 mismatches -- 3 high_sharing pairs compared G6.
    assert g6["n_compared"] == 3
    assert g6["n_mismatch"] == 1
    for g in ["G1", "G2", "G3", "G4", "G5", "G7", "G8", "G9", "G10"]:
        row = gene_out[gene_out["gene"] == g].iloc[0]
        assert row["n_mismatch"] == 0

    # P1~P4 (the simulated crossover) is the only switch-carrying pair.
    switch_map = dict(zip(zip(switches_out["person_i"], switches_out["person_j"]),
                           switches_out["n_switches"]))
    assert switch_map[("P1", "P4")] == 1
    assert switch_map[("P1", "P2")] == 0
    assert switch_map[("P1", "P3")] == 0


if __name__ == "__main__":
    test_build_allele_lookup_and_canonical_order()
    test_clean_parent_child_pair_zero_mismatch_zero_switches()
    test_injected_mismatch_is_localized_to_the_right_gene()
    test_simulated_crossover_gives_exactly_one_switch()
    test_unrelated_pair_is_low_sharing()
    test_contig_boundary_suppresses_spurious_switch()
    test_classify_pair_thresholds()
    test_run_all_pairs_aggregates_mismatch_to_the_correct_gene_and_class()
    print("OK")
