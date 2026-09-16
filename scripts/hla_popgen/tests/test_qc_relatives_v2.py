#!/usr/bin/env python3
"""Fixture tests for 26_qc_relatives_v2.py -- the load-bearing QC script, so every step of the
WS2 spec is checked against a hand-computed expectation before it is trusted on real AoU data.

Covered here:
  * step 1: the pair universe is selected on kinship ONLY (the fix for the audit's circularity
    finding) -- a pair that shares no HLA at all must still be in the universe;
  * step 2: sharing is decided on the exact CDS sequence, so two people with the SAME allele name
    but different sequences are NOT shared, and different names with the same sequence ARE;
  * step 3: one test per discordance category, plus the documented precedence;
  * step 4: replicate QV arithmetic;
  * step 5: the finite-sample-corrected expected homozygosity, against a hand-computed value, and
    the thin-n flag;
  * step 6: testable-transition counting, that assembly fragmentation REDUCES it, and that a
    known synthetic hap-label swap is actually detected;
  * step 7: 2-field multiset matching, field-class severity, and the exclusion of
    undetermined (novelty_depth 1) calls;
  * step 8: an end-to-end run on a complete tiny synthetic dataset (Table 1 + cds.fa.gz +
    relatedness + cohort_membership + SR genotypes), asserting all promised outputs exist, that
    no person id leaks into --out-dir, and that the report says so when the optional year-of-birth
    file is absent.

Runs without pytest (pytest is not installed in the target env): every test is a plain function
taking an optional `tmp_path`.
"""
import gzip
import importlib
import json
import os
import sys
import tempfile
import pathlib

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
mod = importlib.import_module("26_qc_relatives_v2")

CLASSICAL = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
GENE_CLASS = {"A": "classical_I", "B": "classical_I", "C": "classical_I",
              "DPA1": "classical_II", "DPB1": "classical_II", "DQA1": "classical_II",
              "DQB1": "classical_II", "DRB1": "classical_II", "H": "pseudogene_I"}
BASE = ("ACGTACGTAC" * 6)          # 60bp toy CDS


def _seq(tag, n_sub=0, positions=()):
    """A deterministic 60bp toy CDS: `tag` picks a distinct base sequence, `positions` substitutes
    single bases (to build known point / scattered differences)."""
    # sha1 of the tag, not sum(ord(...)): the first version of this helper collided (tags whose
    # character sums differed by 4 produced byte-identical "different" alleles), which silently
    # turned a sharing test into a tautology.
    import hashlib
    h = hashlib.sha1(tag.encode()).digest()
    s = list(BASE)
    for i in range(0, 60, 3):
        s[i] = "ACGT"[h[(i // 3) % len(h)] % 4]
    for p in positions:
        s[p] = "T" if s[p] != "T" else "G"
    return "".join(s)


def _entry(seq1=None, seq2=None, ctg1="c1", ctg2="c1", name1=None, name2=None):
    e = {"seq": {}, "contig": {}, "name": {}}
    if seq1 is not None:
        e["seq"]["hap1"] = seq1
        e["contig"]["hap1"] = ctg1
        if name1:
            e["name"]["hap1"] = name1
    if seq2 is not None:
        e["seq"]["hap2"] = seq2
        e["contig"]["hap2"] = ctg2
        if name2:
            e["name"]["hap2"] = name2
    return e


# ---------------------------------------------------------------------------
# step 1 -- pair universe / strict ancestry
# ---------------------------------------------------------------------------
def test_pair_universe_is_selected_on_kinship_only():
    rel = pd.DataFrame({"person_i": ["P1", "P2", "P3", "P4"],
                        "person_j": ["P2", "P3", "P4", "P5"],
                        "kin": [0.50, 0.25, 0.12, 0.30]})
    cohort = pd.DataFrame({"person_id": ["P1", "P2", "P3", "P4", "P5"],
                           "in_lr": [True, True, True, True, False]})
    out = mod.build_pair_universe(rel, cohort)
    got = set(zip(out["person_i"], out["person_j"]))
    # P3-P4 drops out on kinship (second degree); P4-P5 drops out because P5 has no long-read
    # assembly. Nothing is filtered on HLA sharing -- that is the whole point of the fix.
    assert got == {("P1", "P2"), ("P2", "P3")}, got
    assert list(out["is_replicate"]) == [True, False]
    assert out.iloc[0]["kin_degree"] == "duplicate_or_MZ_twin"


def test_strict_ancestry_requires_own_proportion():
    cohort = pd.DataFrame({
        "person_id": ["A", "B", "C", "D"],
        "ancestry_pred": ["EUR", "EUR", "AFR", "XXX"],
        "p_eur": [0.95, 0.60, 0.10, 0.99], "p_afr": [0.02, 0.35, 0.93, 0.0],
        "p_amr": [0.0] * 4, "p_eas": [0.0] * 4, "p_mid": [0.0] * 4, "p_sas": [0.0] * 4})
    got = mod.strict_ancestry_map(cohort, 0.9)
    assert got == {"A": "EUR", "C": "AFR"}, got   # B is admixed, D's label is not a real ancestry


def test_age_gap_class_needs_both_years():
    yob = {"A": 1950.0, "B": 1980.0, "C": 1978.0}
    assert mod.age_gap_class("A", "B", yob) == "likely_parent_child"
    assert mod.age_gap_class("B", "C", yob) == "likely_sibling"
    assert mod.age_gap_class("A", "Z", yob) == "unknown"


# ---------------------------------------------------------------------------
# step 2 -- sharing is by SEQUENCE, never by name
# ---------------------------------------------------------------------------
def test_sharing_uses_sequence_not_allele_name():
    s1, s2 = _seq("x1"), _seq("x2")
    # Same allele NAME on both sides, but different sequences -> not shared.
    a = {"A": _entry(s1, _seq("x3"), name1="02:01", name2="11:01")}
    b = {"A": _entry(s2, _seq("x4"), name1="02:01", name2="24:02")}
    assert mod.pair_gene_shared(a, b, "A") is False
    # Different names, identical sequence -> shared (Immuannot naming-collapse case from 17).
    b2 = {"A": _entry(s1, _seq("x4"), name1="02:09", name2="24:02")}
    assert mod.pair_gene_shared(a, b2, "A") is True
    assert mod.pair_gene_shared(a, {}, "A") is None          # not comparable


def test_pair_sharing_fraction_and_threshold_sweep():
    shared, other = _seq("s"), _seq("o")
    a = {g: _entry(shared, _seq("a" + g)) for g in CLASSICAL}
    b = {g: _entry(shared if g != "DRB1" else other, _seq("b" + g)) for g in CLASSICAL}
    n_shared, n_cmp, frac = mod.pair_sharing(a, b, CLASSICAL)
    assert (n_shared, n_cmp) == (7, 8) and abs(frac - 7 / 8) < 1e-12
    sweep = mod.sharing_threshold_sweep(
        [{"sharing_frac": 0.9, "n_gene_comparisons": 10, "n_discordant": 1},
         {"sharing_frac": 0.4, "n_gene_comparisons": 10, "n_discordant": 6}],
        thresholds=(0.5, 0.8))
    assert list(sweep["n_pairs"]) == [1, 1]
    assert abs(sweep.iloc[0]["discordance_rate"] - 0.1) < 1e-12


# ---------------------------------------------------------------------------
# step 3 -- one test per mutually exclusive category + precedence
# ---------------------------------------------------------------------------
def _classify(a, b, gene="A", a_singleton=(), b_singleton=()):
    return mod.classify_gene_comparison(a, b, gene, set(a_singleton), set(b_singleton))[0]


def test_category_concordant():
    s = _seq("c")
    assert _classify({"A": _entry(s, _seq("z1"))}, {"A": _entry(s, _seq("z2"))}) == "concordant"


def test_category_missing_one():
    assert _classify({"A": _entry(_seq("m"), None)}, {}) == "missing_one"


def test_category_missing_both():
    assert _classify({}, {}) == "missing_both"


def test_category_point_diff():
    s = _seq("p")
    near = _seq("p", positions=(30,))                       # exactly 1 base, 1 block
    a = {"A": _entry(s, _seq("q1"))}
    b = {"A": _entry(near, _seq("q2"))}
    assert _classify(a, b) == "point_diff"


def test_category_dropout_candidate():
    s = _seq("d")
    scattered = _seq("d", positions=(0, 10, 20, 30, 40))    # 5 bases, 5 blocks
    a = {"A": _entry(s, s)}                                 # homozygous by sequence
    b = {"A": _entry(scattered, _seq("d2"))}
    assert _classify(a, b) == "dropout_candidate"
    # Only ONE hap called is also "homozygous by observation" -- the dropout shape.
    a_one = {"A": _entry(s, None)}
    assert _classify(a_one, b) == "dropout_candidate"


def test_category_fragmentation_candidate():
    s = _seq("f")
    scattered = _seq("f", positions=(0, 10, 20, 30, 40))
    a = {"A": _entry(s, _seq("f2"))}                        # heterozygous -> not dropout
    b = {"A": _entry(scattered, _seq("f3"))}
    assert _classify(a, b, a_singleton={"A"}) == "fragmentation_candidate"


def test_category_other_multiblock():
    s = _seq("o1")
    scattered = _seq("o1", positions=(0, 10, 20, 30, 40))
    a = {"A": _entry(s, _seq("o2"))}
    b = {"A": _entry(scattered, _seq("o3"))}
    assert _classify(a, b) == "other_multiblock"


def test_category_precedence_point_diff_beats_dropout():
    """A 1-base difference at a homozygous locus is consensus error (audit finding 5), not
    dropout -- the documented precedence must keep it out of the dropout bucket."""
    s = _seq("pd")
    a = {"A": _entry(s, s)}                                 # homozygous AND singleton contig
    b = {"A": _entry(_seq("pd", positions=(15,)), _seq("pd2"))}
    assert _classify(a, b, a_singleton={"A"}) == "point_diff"


def test_categories_are_mutually_exclusive_and_exhaustive():
    s = _seq("e")
    a = {g: _entry(s, _seq("ea" + g)) for g in CLASSICAL}
    b = {g: _entry(s if g != "B" else _seq("eb"), _seq("eb" + g)) for g in CLASSICAL}
    recs = mod.decompose_pair(a, b, CLASSICAL + ["H"], tuple(CLASSICAL))
    assert len(recs) == 9                                    # exactly one row per gene, no drops
    assert all(r["category"] in mod.DISCORDANCE_CATEGORIES for r in recs)
    agg = mod.aggregate_decomposition(recs, GENE_CLASS)
    # n_total == sum of the category columns, for every gene: partition, not overlapping tags.
    cats = agg[mod.DISCORDANCE_CATEGORIES].sum(axis=1)
    assert (cats == agg["n_total"]).all()
    # 'H' was never called -> missing_both -> excluded from the rate's denominator.
    h = agg[agg["gene"] == "H"].iloc[0]
    assert h["missing_both"] == 1 and h["n_comparable"] == 0


def test_singleton_contig_genes():
    pmap = {"A": _entry(_seq("1"), _seq("2"), ctg1="cA", ctg2="cA"),
            "B": _entry(_seq("3"), _seq("4"), ctg1="cA", ctg2="cA"),
            "C": _entry(_seq("5"), _seq("6"), ctg1="cLONE", ctg2="cLONE")}
    assert mod.singleton_contig_genes(pmap, ("A", "B", "C")) == {"C"}


# ---------------------------------------------------------------------------
# step 4 -- replicates / QV
# ---------------------------------------------------------------------------
def test_replicate_diff_tries_both_hap_assignments_and_qv():
    s1, s2 = _seq("r1"), _seq("r2")
    a = {"A": _entry(s1, s2)}
    b = {"A": _entry(s2, s1)}                               # same genome, hap labels swapped
    d = mod.replicate_gene_diff(a, b, "A")
    assert d["n_diff_bases"] == 0 and d["both_equal"] is True
    b_err = {"A": _entry(s2, _seq("r1", positions=(5,)))}    # one base wrong
    d2 = mod.replicate_gene_diff(a, b_err, "A")
    assert d2["n_diff_bases"] == 1
    qv, err, n_diff, n_aln = mod.replicate_qv([{"n_diff_bases": 1, "n_aligned_bases": 1000}])
    assert (n_diff, n_aln) == (1, 1000) and abs(qv - 30.0) < 1e-9
    qv0, _, _, _ = mod.replicate_qv([{"n_diff_bases": 0, "n_aligned_bases": 1000}])
    assert qv0 == float("inf")


# ---------------------------------------------------------------------------
# step 5 -- homozygosity obs/exp
# ---------------------------------------------------------------------------
def test_expected_homozygosity_finite_sample_correction():
    # counts 6 and 4, n=10: sum x(x-1)/(n(n-1)) = (30 + 12)/90 = 0.4667, strictly BELOW the
    # biased naive sum p^2 = 0.36+0.16 = 0.52? -- no: the correction is smaller than naive here,
    # which is the documented direction (naive is biased upward).
    got = mod.expected_homozygosity([6, 4])
    assert abs(got - (30 + 12) / 90) < 1e-12
    assert got < 0.36 + 0.16
    assert np.isnan(mod.expected_homozygosity([1]))


def test_homozygosity_obs_exp_flags_thin_groups_and_is_deterministic():
    calls = {("LR_CDS_exact", "A", "EUR"): [("s1", "s1")] * 3 + [("s1", "s2")] * 3}
    d1 = mod.homozygosity_obs_exp(calls, n_boot=50, seed=7, min_people=50)
    d2 = mod.homozygosity_obs_exp(calls, n_boot=50, seed=7, min_people=50)
    row = d1.iloc[0]
    assert row["n_people"] == 6 and abs(row["obs_hom"] - 0.5) < 1e-12
    assert bool(row["thin_n"]) is True
    assert row["ci_lo"] <= row["obs_exp_ratio"] <= row["ci_hi"]
    pd.testing.assert_frame_equal(d1, d2)          # seeded bootstrap, identical across calls


def test_fragmentation_class():
    one = {g: _entry(_seq(g), _seq(g + "b"), ctg1="c1", ctg2="c1") for g in CLASSICAL}
    assert mod.fragmentation_class(one, tuple(CLASSICAL)) == "1"
    three = dict(one)
    three["A"] = _entry(_seq("A"), _seq("Ab"), ctg1="c2", ctg2="c1")
    three["B"] = _entry(_seq("B"), _seq("Bb"), ctg1="c3", ctg2="c1")
    assert mod.fragmentation_class(three, tuple(CLASSICAL)) == ">=3"
    assert mod.fragmentation_class({}, tuple(CLASSICAL)) is None


# ---------------------------------------------------------------------------
# step 6 -- switch test denominator + empirical power
# ---------------------------------------------------------------------------
def _switch_pair(contigs_b=None):
    """A clean 'parent/child' pair over the 8 classical genes: person B's hap1 is an exact copy of
    person A's hap1 at every gene, everything else distinct. `contigs_b` optionally fragments B."""
    order = {g: i for i, g in enumerate(CLASSICAL)}
    a = {g: _entry(_seq("sh" + g), _seq("aa" + g)) for g in CLASSICAL}
    b = {}
    for g in CLASSICAL:
        ctg = (contigs_b or {}).get(g, "c1")
        b[g] = _entry(_seq("sh" + g), _seq("bb" + g), ctg1=ctg, ctg2=ctg)
    return a, b, order


def test_testable_transitions_and_fragmentation_reduces_them():
    a, b, order = _switch_pair()
    detail = mod.resolve_states(a, b, order, CLASSICAL)
    n_testable, n_resolved = mod.testable_transitions(detail)
    assert n_resolved == 8 and n_testable == 7          # 8 resolved genes -> 7 transitions
    assert mod.count_switches(detail)[0] == 0           # clean pair, no switch

    # Same genotypes, but B's assembly is broken into 4 fragments: the SAME data now supports
    # fewer testable transitions. This is the audit's "power is unknown" point, made concrete.
    frag = {"A": "c1", "B": "c1", "C": "c2", "DPA1": "c2", "DPB1": "c3", "DQA1": "c3",
            "DQB1": "c4", "DRB1": "c4"}
    a2, b2, _ = _switch_pair(frag)
    detail2 = mod.resolve_states(a2, b2, order, CLASSICAL)
    n_testable2, n_resolved2 = mod.testable_transitions(detail2)
    assert n_resolved2 == 8
    assert n_testable2 == 4 and n_testable2 < n_testable


def test_synthetic_switch_with_known_breakpoint_is_detected():
    a, b, order = _switch_pair()
    detail = mod.resolve_states(a, b, order, CLASSICAL)
    sim, info = mod.apply_synthetic_switch(b, detail, np.random.default_rng(3))
    assert info["n_genes_swapped"] >= 1
    sim_detail = mod.resolve_states(a, sim, order, CLASSICAL)
    assert mod.count_switches(sim_detail)[0] >= 1
    assert b["A"]["seq"]["hap1"] == _seq("shA")        # input was not mutated

    # Hand-built swap, independent of the simulator: flip hap labels from DPB1 onward.
    manual = {g: {"seq": dict(v["seq"]), "contig": dict(v["contig"]), "name": dict(v["name"])}
              for g, v in b.items()}
    for g in CLASSICAL[4:]:
        s = manual[g]["seq"]
        s["hap1"], s["hap2"] = s["hap2"], s["hap1"]
    assert mod.count_switches(mod.resolve_states(a, manual, order, CLASSICAL))[0] == 1


def test_switch_power_is_deterministic_and_reports_the_denominator():
    a, b, order = _switch_pair()
    frag = {g: ("c1" if i < 2 else f"c{i}") for i, g in enumerate(CLASSICAL)}
    a2, b2, _ = _switch_pair(frag)
    df1, s1 = mod.switch_power([(a, b), (a2, b2)], order, CLASSICAL, seed=11)
    df2, s2 = mod.switch_power([(a, b), (a2, b2)], order, CLASSICAL, seed=11)
    pd.testing.assert_frame_equal(df1, df2)
    assert s1 == s2
    assert s1["n_pairs"] == 2 and s1["n_eligible_pairs"] >= 1
    assert s1["total_testable_transitions"] == int(df1["n_testable_transitions"].sum())
    assert s1["detection_rate"] == 1.0
    # A pair with no same-contig neighbours is ineligible and must be reported, not dropped.
    lone = {g: _entry(_seq("l" + g), _seq("m" + g), ctg1=f"x{g}", ctg2=f"x{g}")
            for g in CLASSICAL}
    df3, s3 = mod.switch_power([(a, lone)], order, CLASSICAL, seed=5)
    assert s3["n_eligible_pairs"] == 0 and s3["n_pairs_zero_testable"] == 1


# ---------------------------------------------------------------------------
# step 7 -- SR vs LR concordance
# ---------------------------------------------------------------------------
def test_genotype_match_count_is_a_multiset_intersection():
    assert mod.genotype_match_count(["02:01", "11:01"], ["02:01", "11:01"]) == 2
    assert mod.genotype_match_count(["02:01", "11:01"], ["02:01", "24:02"]) == 1
    assert mod.genotype_match_count(["02:01", "02:01"], ["02:01", "11:01"]) == 1
    assert mod.genotype_match_count(["02:01", "11:01"], ["24:02", "31:01"]) == 0


def test_lr_field_class_severity_and_undetermined_exclusion():
    assert mod.lr_field_class({"hap1": (False, None), "hap2": (False, None)}) == "known"
    assert mod.lr_field_class({"hap1": (True, 4), "hap2": (False, None)}) == "f4_noncoding"
    assert mod.lr_field_class({"hap1": (True, 3), "hap2": (True, 4)}) == "f3_synonymous"
    assert mod.lr_field_class({"hap1": (True, 2), "hap2": (True, 4)}) == "f2_protein"
    assert mod.lr_field_class({"hap1": (True, 1), "hap2": (False, None)}) is None
    assert mod.lr_field_class({"hap1": (True, None)}) is None


def test_sr_lr_concordance_rows_counts_and_excludes():
    lr = {("P1", "A"): ["02:01", "11:01"], ("P2", "A"): ["02:01", "11:01"],
          ("P3", "A"): ["02:01", "11:01"]}
    sr = {("P1", "A"): ["02:01", "11:01"], ("P2", "A"): ["02:01", "24:02"],
          ("P3", "A"): ["02:01", "11:01"]}
    anc = {"P1": "EUR", "P2": "EUR", "P3": "EUR"}
    nov = {("P1", "A"): {"hap1": (False, None), "hap2": (True, 4)},
           ("P2", "A"): {"hap1": (True, 2), "hap2": (False, None)},
           ("P3", "A"): {"hap1": (True, 1), "hap2": (False, None)}}   # undetermined -> excluded
    df, n_excluded = mod.sr_lr_concordance_rows(lr, sr, anc, nov)
    assert n_excluded == 1 and len(df) == 2
    by_class = dict(zip(df["field_class"], df["n_match"]))
    assert by_class["f4_noncoding"] == 2 and by_class["f2_protein"] == 1
    agg = mod.aggregate_sr_lr(df, min_people=50)
    assert set(agg["field_class"]) == {"f4_noncoding", "f2_protein"}
    assert bool(agg["thin_n"].all()) is True


def test_safe_count_masks_small_person_counts():
    assert mod.safe_count(0) == "0"
    assert mod.safe_count(1) == "<20"
    assert mod.safe_count(19) == "<20"
    assert mod.safe_count(20) == "20"
    assert mod.safe_count(None) == "NA"


# ---------------------------------------------------------------------------
# step 8 -- end-to-end on a complete tiny synthetic dataset
# ---------------------------------------------------------------------------
PEOPLE = [f"P{i:02d}" for i in range(1, 13)]
# A tiny allele pool per gene so that homozygotes occur and frequencies are estimable.
ALLELE_POOL = {g: [f"{g}*01:01", f"{g}*02:01", f"{g}*03:01"] for g in CLASSICAL}


def _person_genotype(pid, gene):
    """Deterministic pseudo-genotype: person index mod pool size, with every 4th person made
    homozygous so the homozygosity step has something to measure."""
    i = PEOPLE.index(pid)
    pool = ALLELE_POOL[gene]
    a1 = pool[(i + CLASSICAL.index(gene)) % len(pool)]
    a2 = a1 if i % 4 == 0 else pool[(i + CLASSICAL.index(gene) + 1) % len(pool)]
    return a1, a2


def _allele_seq(allele):
    return _seq(allele)


def _write_cds(root, pid, hap, records):
    path = os.path.join(root, "people", pid, "immuannot_output", hap, "cds.fa.gz")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as f:
        for contig, gene, idx, seq in records:
            # Real header format: '>{contig}_{gene}_{i} {splice_sites}' (see 17's parser).
            f.write(f">{contig}_{gene}_{idx} GT_AG:GT_AG\n{seq}\n")


def _build_dataset(root):
    os.makedirs(root, exist_ok=True)
    t1_rows, sr_rows, cohort_rows = [], [], []
    for pid in PEOPLE:
        i = PEOPLE.index(pid)
        hap_records = {"hap1": [], "hap2": []}
        for gi, gene in enumerate(CLASSICAL):
            a1, a2 = _person_genotype(pid, gene)
            # P05 is deliberately fragmented: its DRB1 sits alone on a second contig.
            contig = "ctg_b" if (pid == "P05" and gene == "DRB1") else "ctg_a"
            for hap, allele in (("hap1", a1), ("hap2", a2)):
                hap_records[hap].append((contig, f"HLA-{gene}", 1, _allele_seq(allele)))
                # Novel (depth 4 / noncoding-only) for one gene in one person, and a
                # depth-1 undetermined call for another, to exercise step 7's field classes.
                consensus = allele.split("*")[1]
                is_novel, depth = False, ""
                if pid == "P02" and gene == "A":
                    consensus = consensus + ":01:new"
                    is_novel, depth = True, 4
                if pid == "P03" and gene == "B":
                    consensus = "new"
                    is_novel, depth = True, 1
                t1_rows.append({
                    "person_id": pid, "hap": hap, "contig": contig, "gene": f"HLA-{gene}",
                    "copy_index": 1, "gene_class": GENE_CLASS[gene],
                    "consensus": f"{gene}*{consensus}", "n_fields": 2,
                    "is_novel": is_novel, "novelty_depth": depth,
                    "gene_start": 1000 * (gi + 1), "gene_end": 1000 * (gi + 1) + 500,
                    "template_distance": 0, "has_warning": False})
            sr_rows.append((pid, gene, a1, a2))
        for hap, recs in hap_records.items():
            _write_cds(root, pid, hap, recs)
        cohort_rows.append({"person_id": pid, "in_sr": True, "in_lr": True,
                            "ancestry_pred": "eur",                  # lowercase, as in real files
                            "p_afr": 0.01, "p_amr": 0.01, "p_eas": 0.01,
                            "p_eur": 0.95 if i % 5 else 0.5,         # one admixed person
                            "p_mid": 0.01, "p_sas": 0.01,
                            "max_template_distance": 0, "mean_template_distance": 0.0,
                            "n_genes_called": len(CLASSICAL), "n_genes_exact": len(CLASSICAL),
                            "td_stratum": "exact"})

    pd.DataFrame(t1_rows).to_csv(os.path.join(root, "hla_calls_rich.tsv"), sep="\t", index=False)
    pd.DataFrame(cohort_rows).to_csv(os.path.join(root, "cohort_membership.tsv"),
                                     sep="\t", index=False)
    # SR genotype file: wide, research_id + <gene>_1/<gene>_2 (07 / _viz_common.load_sr_genotypes).
    wide = {"research_id": PEOPLE}
    for gene in CLASSICAL:
        wide[f"{gene}_1"] = [_person_genotype(p, gene)[0] for p in PEOPLE]
        wide[f"{gene}_2"] = [_person_genotype(p, gene)[1] for p in PEOPLE]
    pd.DataFrame(wide).to_csv(os.path.join(root, "hla_genotypes.tsv"), sep="\t", index=False)
    # Relatedness: one duplicate/MZ pair, three first-degree pairs, one second-degree (excluded).
    rel = pd.DataFrame({"i.s": ["P01", "P03", "P05", "P07", "P09"],
                        "j.s": ["P02", "P04", "P06", "P08", "P10"],
                        "kin": [0.49, 0.25, 0.22, 0.20, 0.12]})
    rel.to_csv(os.path.join(root, "samples_relatedness.tsv"), sep="\t", index=False)
    return root


def _args(root, out_dir, **over):
    ns = mod.build_parser().parse_args([
        "--table1", os.path.join(root, "hla_calls_rich.tsv"),
        "--cohort-membership", os.path.join(root, "cohort_membership.tsv"),
        "--relatedness-table", os.path.join(root, "samples_relatedness.tsv"),
        "--sr-genotypes", os.path.join(root, "hla_genotypes.tsv"),
        "--people-root", os.path.join(root, "people"),
        "--out-dir", out_dir,
        "--local-dir", os.path.join(root, "local"),
        "--bootstrap", "20", "--threads", "2"])
    for k, v in over.items():
        setattr(ns, k, v)
    return ns


def test_end_to_end_run_produces_every_promised_output(tmp_path):
    root = _build_dataset(os.path.join(str(tmp_path), "data"))
    out_dir = os.path.join(str(tmp_path), "out")
    res = mod.main_from_args(_args(root, out_dir))

    expected = ["qc_relatives_v2_report.md", "discordance_decomposition.tsv",
                "sharing_histogram.png", "decomposition_by_gene.png",
                "homozygosity_obs_exp.png", "switch_power.png",
                "sr_lr_concordance.tsv", "sr_lr_concordance.png", "summary.json"]
    for name in expected:
        p = os.path.join(out_dir, name)
        assert os.path.exists(p) and os.path.getsize(p) > 0, f"missing/empty output: {name}"

    # 4 pairs pass the first-degree kinship floor (the 0.12 pair is correctly excluded).
    assert len(res["pairs"]) == 4
    assert int(res["pairs"]["is_replicate"].sum()) == 1
    # Every gene comparison landed in exactly one category.
    decomp = res["decomposition"]
    assert (decomp[mod.DISCORDANCE_CATEGORIES].sum(axis=1) == decomp["n_total"]).all()
    # Homozygosity ran for all three like-for-like measures on the same people.
    assert set(res["homozygosity"]["measure"]) == {"LR_CDS_exact", "LR_2field", "SR_2field"}
    # The admixed person (p_eur == 0.5) must be absent from the strict-ancestry step.
    assert res["homozygosity"]["n_people"].max() <= len(PEOPLE) - 2
    # SR and LR are the same calls here, so 2-field concordance must be perfect for known alleles.
    sr_lr = res["sr_lr"]
    known = sr_lr[sr_lr["field_class"] == "known"]
    assert not known.empty and abs(known["mean_alleles_matching"].min() - 2.0) < 1e-9
    # The depth-1 undetermined call (P03 / HLA-B) was excluded and counted, not silently dropped.
    summary = json.load(open(os.path.join(out_dir, "summary.json")))
    assert summary["sr_lr_excluded_undetermined"] == "<20"
    assert summary["yob_stratification"] == "skipped (no --yob-tsv)"
    assert summary["switch_power"]["n_pairs"] == 3           # replicate pair excluded from step 6
    # Person-level rows went ONLY to the local dir.
    assert os.path.exists(os.path.join(root, "local", "26_pair_level_rows.tsv"))


def test_end_to_end_writes_no_person_ids_and_masks_small_counts(tmp_path):
    root = _build_dataset(os.path.join(str(tmp_path), "data"))
    out_dir = os.path.join(str(tmp_path), "out")
    mod.main_from_args(_args(root, out_dir))
    for name in sorted(os.listdir(out_dir)):
        if name.endswith(".png"):
            continue
        text = open(os.path.join(out_dir, name)).read()
        for pid in PEOPLE:
            assert pid not in text, f"person id {pid} leaked into {name}"
    report = open(os.path.join(out_dir, "qc_relatives_v2_report.md")).read()
    assert "SKIPPED" in report and "--yob-tsv" in report      # says so, never silently omits
    assert "<20" in report                                    # small counts are masked


def test_end_to_end_with_yob_stratifies_by_age_gap(tmp_path):
    root = _build_dataset(os.path.join(str(tmp_path), "data"))
    out_dir = os.path.join(str(tmp_path), "out")
    yob_path = os.path.join(root, "yob.tsv")
    pd.DataFrame({"person_id": PEOPLE,
                  "year_of_birth": [1950 + (30 if i % 2 else 0) for i in range(len(PEOPLE))]
                  }).to_csv(yob_path, sep="\t", index=False)
    mod.main_from_args(_args(root, out_dir, yob_tsv=yob_path))
    report = open(os.path.join(out_dir, "qc_relatives_v2_report.md")).read()
    assert "likely_parent_child" in report and "SKIPPED" not in report
    summary = json.load(open(os.path.join(out_dir, "summary.json")))
    assert summary["yob_stratification"] == "applied"


def test_end_to_end_limit_pairs_pilot_mode(tmp_path):
    root = _build_dataset(os.path.join(str(tmp_path), "data"))
    out_dir = os.path.join(str(tmp_path), "out")
    res = mod.main_from_args(_args(root, out_dir, limit_pairs=2))
    assert len(res["pairs"]) == 2


if __name__ == "__main__":
    import inspect
    for name, fn in inspect.getmembers(sys.modules[__name__], inspect.isfunction):
        if not name.startswith("test_"):
            continue
        kw = ({"tmp_path": pathlib.Path(tempfile.mkdtemp())}
              if "tmp_path" in inspect.signature(fn).parameters else {})
        fn(**kw)
        print(f"ok {name}")
    print("ok")


def test_cds_header_ordinal_is_not_copy_index(tmp_path):
    """REGRESSION (real-data bug, 2026-09-16): the trailing integer in a cds.fa.gz header runs
    across the whole file (HLA-E_1, HLA-L_2, HLA-A_3 ...), it is NOT the gene's copy_index.
    Treating it as copy_index and keeping only '1' dropped every gene but the file's first record,
    which made 78% of relative gene-comparisons look missing. Each gene here appears once, with
    ordinals 1..3, and all three sequences must come back."""
    import gzip
    root = tmp_path / "people" / "P1" / "immuannot_output" / "hap1"
    root.mkdir(parents=True)
    with gzip.open(str(root / "cds.fa.gz"), "wt") as f:
        for ordinal, (gene, seq) in enumerate(
                [("HLA-E", "AAA"), ("HLA-A", "CCC"), ("HLA-DRB1", "GGG")], start=1):
            f.write(">ctg1_%s_%d GT_AG:GT_AG\n%s\n" % (gene, ordinal, seq))
    cds = mod.load_person_cds(str(tmp_path / "people"), "P1", ["A", "E", "DRB1"])
    assert cds[("hap1", "A")] == ["CCC"]
    assert cds[("hap1", "E")] == ["AAA"]
    assert cds[("hap1", "DRB1")] == ["GGG"]
    pmap = mod.build_person_map(cds, {}, {}, "P1", ["A", "E", "DRB1"])
    assert pmap["A"]["seq"]["hap1"] == "CCC"
    assert pmap["DRB1"]["seq"]["hap1"] == "GGG"


def test_multi_record_gene_is_treated_as_ambiguous(tmp_path):
    """Two records for the same gene on one hap (DRB paralog / segmental duplication) cannot be
    mapped to Table 1's copy_index from this file, so the gene carries no copy_index==1 sequence."""
    import gzip
    root = tmp_path / "people" / "P2" / "immuannot_output" / "hap1"
    root.mkdir(parents=True)
    with gzip.open(str(root / "cds.fa.gz"), "wt") as f:
        f.write(">ctg1_HLA-DRB1_1 GT_AG\nAAA\n>ctg2_HLA-DRB1_2 GT_AG\nTTT\n")
    cds = mod.load_person_cds(str(tmp_path / "people"), "P2", ["DRB1"])
    assert sorted(cds[("hap1", "DRB1")]) == ["AAA", "TTT"]
    pmap = mod.build_person_map(cds, {}, {}, "P2", ["DRB1"])
    assert "hap1" not in pmap.get("DRB1", {}).get("seq", {})
