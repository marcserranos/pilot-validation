#!/usr/bin/env python3
"""Unit tests for 29_hla_ld_by_ancestry.py.

The tests that matter here are the ones that would catch a *silently wrong* LD number:

  1. r^2 against hand-computed values on a 2x2 with known D, including perfect LD and perfect
     linkage equilibrium. If this is wrong, every figure is wrong and nothing looks odd.
  2. D' vs r^2 divergence -- the classic case where complete LD (D'=1) coexists with a small r^2
     because the marginal frequencies differ. This is exactly the trap Cole's "just compute r^2"
     instruction walks into, so the script must get it right and report both.
  3. The cis key. A haplotype may only be formed from two genes on the SAME contig. A test where
     the two genes sit on different contigs of the same assembly haplotype must yield ZERO
     haplotypes -- not a pairing.
  4. Multi-copy genes are dropped, not silently collapsed to copy 1 (the S01 trap that made 78% of
     relative comparisons look discordant was exactly this class of mistake).
  5. Bergsma's bias correction actually reduces the table-size inflation it exists to remove.

Run: python3 scripts/hla_popgen/tests/test_hla_ld_by_ancestry.py
"""
import importlib.util
import math
import os
import random
import sys

import numpy as np
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


m = _load_module("29_hla_ld_by_ancestry.py", "hla_ld")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


def _rec(recs, a, b):
    for r in recs:
        if r["allele_a"] == a and r["allele_b"] == b:
            return r
    return None


# ---------------------------------------------------------------------------
# two_field
# ---------------------------------------------------------------------------
def test_two_field_truncation():
    check("two_field keeps exactly two fields",
          m.two_field("HLA-DQB1*02:01:01:02") == "DQB1*02:01",
          m.two_field("HLA-DQB1*02:01:01:02"))
    check("two_field strips the HLA- prefix",
          m.two_field("HLA-A*01:01") == "A*01:01")
    check("a call novel at field 2 has no protein identity and must be dropped",
          m.two_field("HLA-A*01:new") is None, str(m.two_field("HLA-A*01:new")))
    check("a call novel at field 1 must be dropped",
          m.two_field("HLA-A*new") is None)
    check("field-3 novelty keeps its protein identity and is kept",
          m.two_field("HLA-A*01:01:new") == "A*01:01",
          str(m.two_field("HLA-A*01:01:new")))
    check("undetermined is dropped", m.two_field("undetermined") is None)
    check("a one-field call is dropped", m.two_field("HLA-A*01") is None)


# ---------------------------------------------------------------------------
# pairwise r^2 / D'
# ---------------------------------------------------------------------------
def test_perfect_ld_gives_r2_one():
    # 50 haplotypes a1~b1, 50 haplotypes a2~b2. Perfect correlation.
    M = np.array([[50.0, 0.0], [0.0, 50.0]])
    recs = m.pairwise_ld(M, ["a1", "a2"], ["b1", "b2"], min_haps=1)
    r = _rec(recs, "a1", "b1")
    check("perfect LD on equal frequencies gives r^2 = 1", abs(r["r2"] - 1.0) < 1e-12, str(r["r2"]))
    check("perfect LD gives D' = 1", abs(r["Dprime"] - 1.0) < 1e-12, str(r["Dprime"]))


def test_linkage_equilibrium_gives_r2_zero():
    # Product of marginals exactly: D = 0 everywhere.
    M = np.array([[25.0, 25.0], [25.0, 25.0]])
    recs = m.pairwise_ld(M, ["a1", "a2"], ["b1", "b2"], min_haps=1)
    r = _rec(recs, "a1", "b1")
    check("linkage equilibrium gives r^2 = 0", abs(r["r2"]) < 1e-12, str(r["r2"]))
    check("linkage equilibrium gives D = 0", abs(r["D"]) < 1e-12, str(r["D"]))


def test_r2_matches_hand_computed_value():
    # n=100. p(a1)=0.5, p(b1)=0.5, p(a1b1)=0.4 -> D = 0.4 - 0.25 = 0.15
    # r^2 = D^2 / (0.5*0.5*0.5*0.5) = 0.0225 / 0.0625 = 0.36
    M = np.array([[40.0, 10.0], [10.0, 40.0]])
    recs = m.pairwise_ld(M, ["a1", "a2"], ["b1", "b2"], min_haps=1)
    r = _rec(recs, "a1", "b1")
    check("r^2 matches the hand-computed 0.36", abs(r["r2"] - 0.36) < 1e-12, str(r["r2"]))
    check("chi2 on 1 df equals n * r^2", abs(r["chi2_1df"] - 36.0) < 1e-9, str(r["chi2_1df"]))


def test_complete_ld_with_small_r2():
    """The trap. a1 occurs on 10% of haplotypes and ALWAYS with b1; b1 is on 50%. There is no
    a1~b2 haplotype at all, so D' = 1 (complete LD), yet r^2 is only ~0.11. Reporting r^2 alone
    would call this weak linkage; it is the strongest linkage the frequencies permit."""
    #          b1    b2
    #  a1      10     0
    #  a2      40    50
    M = np.array([[10.0, 0.0], [40.0, 50.0]])
    recs = m.pairwise_ld(M, ["a1", "a2"], ["b1", "b2"], min_haps=1)
    r = _rec(recs, "a1", "b1")
    check("complete LD at unequal frequencies still gives D' = 1",
          abs(r["Dprime"] - 1.0) < 1e-12, str(r["Dprime"]))
    check("...while r^2 stays small, which is why both are reported",
          r["r2"] < 0.2, str(r["r2"]))


def test_min_allele_floor_excludes_rare_alleles():
    M = np.array([[50.0, 50.0], [1.0, 1.0]])
    recs = m.pairwise_ld(M, ["common", "rare"], ["b1", "b2"], min_haps=20)
    check("an allele below the haplotype floor produces no rows",
          all(r["allele_a"] != "rare" for r in recs), str([r["allele_a"] for r in recs]))


# ---------------------------------------------------------------------------
# multi-allelic statistics
# ---------------------------------------------------------------------------
def test_multiallelic_dprime_one_on_perfect_association():
    M = np.eye(4) * 25.0
    st = m.multiallelic_stats(M)
    check("a perfectly associated 4x4 table gives multi-allelic D' = 1",
          abs(st["Dprime_multi"] - 1.0) < 1e-9, str(st["Dprime_multi"]))
    check("...and Cramer's V = 1", abs(st["cramers_v"] - 1.0) < 1e-9, str(st["cramers_v"]))


def test_multiallelic_dprime_zero_at_equilibrium():
    M = np.full((4, 4), 25.0)
    st = m.multiallelic_stats(M)
    check("an independent table gives multi-allelic D' = 0",
          abs(st["Dprime_multi"]) < 1e-9, str(st["Dprime_multi"]))
    check("...and NMI = 0", abs(st["nmi"]) < 1e-9, str(st["nmi"]))


def test_bergsma_correction_removes_table_size_inflation():
    """The reason the bias-corrected V exists. Draw two INDEPENDENT loci, one with few alleles and
    one with many, at the same n. Raw Cramer's V is visibly inflated for the many-allele table;
    the corrected V is near zero for both. Without this, a cross-ancestry comparison would report
    'more LD in AFR' purely because AFR has more alleles."""
    rng = np.random.default_rng(7)
    n = 400

    def v_pair(k):
        a = rng.integers(0, k, n)
        b = rng.integers(0, k, n)
        M, al, bl = m.contingency([f"a{x}" for x in a], [f"b{x}" for x in b])
        return m.multiallelic_stats(M)

    small = v_pair(3)
    large = v_pair(20)
    check("raw Cramer's V is inflated by table size under independence",
          large["cramers_v"] > small["cramers_v"] + 0.02,
          f"small={small['cramers_v']:.3f} large={large['cramers_v']:.3f}")
    check("bias-corrected V stays near zero for both",
          large["cramers_v_bc"] < 0.1 and small["cramers_v_bc"] < 0.1,
          f"small={small['cramers_v_bc']:.3f} large={large['cramers_v_bc']:.3f}")


def test_rarefaction_is_without_replacement_and_respects_size():
    pairs = [("a%d" % (i % 5), "b%d" % (i % 5)) for i in range(200)]
    rng = random.Random(1)
    out = m.rarefied_stats(pairs, 50, 20, rng)
    check("rarefied_stats reports the size it used", out["rarefied_to"] == 50, str(out))
    check("perfectly linked pairs stay at D' = 1 after rarefaction",
          abs(out["Dprime_multi_mean"] - 1.0) < 1e-6, str(out["Dprime_multi_mean"]))
    check("rarefaction refuses to extrapolate beyond the sample",
          m.rarefied_stats(pairs, 500, 5, rng) is None)


def test_rarefaction_lowers_apparent_allele_count():
    """Sanity check that rarefaction is doing the thing it is for: a subsample of 30 haplotypes
    from a 50-allele locus cannot show 50 alleles."""
    pairs = [("a%d" % i, "b%d" % i) for i in range(200)]
    out = m.rarefied_stats(pairs, 30, 20, random.Random(3))
    check("a 30-haplotype subsample shows at most 30 alleles",
          out["n_alleles_a_mean"] <= 30.0, str(out["n_alleles_a_mean"]))


# ---------------------------------------------------------------------------
# haplotype extraction -- the cis key
# ---------------------------------------------------------------------------
def _t1(rows):
    return pd.DataFrame(rows, columns=["person_id", "hap", "contig", "gene", "copy_index",
                                       "consensus"])


def test_same_contig_yields_a_haplotype():
    t1 = _t1([
        ("p1", "hap1", "ctg1", "HLA-DQA1", 1, "HLA-DQA1*01:01:01"),
        ("p1", "hap1", "ctg1", "HLA-DQB1", 1, "HLA-DQB1*05:01:01"),
    ])
    haps, diag = m.extract_haplotypes(t1, "DQA1", "DQB1")
    check("two genes on one contig give exactly one phased haplotype", len(haps) == 1, str(haps))
    check("the alleles are truncated to two fields",
          list(haps.iloc[0][["allele_a", "allele_b"]]) == ["DQA1*01:01", "DQB1*05:01"],
          str(haps.iloc[0].to_dict()))


def test_different_contigs_yield_nothing():
    """THE test. Same person, same assembly haplotype, but the two genes landed on different
    contigs -- so their phase relative to each other is unknown. Pairing them would invent a
    haplotype that was never observed, and it would look completely normal in the output."""
    t1 = _t1([
        ("p1", "hap1", "ctg1", "HLA-DQA1", 1, "HLA-DQA1*01:01"),
        ("p1", "hap1", "ctg2", "HLA-DQB1", 1, "HLA-DQB1*05:01"),
    ])
    haps, diag = m.extract_haplotypes(t1, "DQA1", "DQB1")
    check("genes on different contigs of the same hap yield NO haplotype", len(haps) == 0,
          str(haps))
    check("...and the fragmentation is counted, not silently dropped",
          diag.get("contig_carries_only_one_gene", 0) == 2, str(diag))
    check("...while the assembly is still counted as carrying both genes somewhere",
          diag.get("assemblies_with_both_genes_somewhere", 0) == 1, str(diag))


def test_multi_copy_gene_is_dropped_not_collapsed():
    t1 = _t1([
        ("p1", "hap1", "ctg1", "HLA-DQA1", 1, "HLA-DQA1*01:01"),
        ("p1", "hap1", "ctg1", "HLA-DQA1", 2, "HLA-DQA1*05:01"),
        ("p1", "hap1", "ctg1", "HLA-DQB1", 1, "HLA-DQB1*05:01"),
    ])
    haps, diag = m.extract_haplotypes(t1, "DQA1", "DQB1")
    check("an ambiguous two-copy gene yields no haplotype", len(haps) == 0, str(haps))
    check("...and is counted as multi_copy_ambiguous",
          diag.get("multi_copy_ambiguous", 0) == 1, str(diag))


def test_unresolved_call_is_dropped_and_counted():
    t1 = _t1([
        ("p1", "hap1", "ctg1", "HLA-DQA1", 1, "HLA-DQA1*new"),
        ("p1", "hap1", "ctg1", "HLA-DQB1", 1, "HLA-DQB1*05:01"),
    ])
    haps, diag = m.extract_haplotypes(t1, "DQA1", "DQB1")
    check("a field-1-novel call yields no haplotype", len(haps) == 0, str(haps))
    check("...and is counted separately from fragmentation",
          diag.get("unresolved_2field_call", 0) == 1, str(diag))


def test_both_assembly_haplotypes_contribute():
    t1 = _t1([
        ("p1", "hap1", "ctg1", "HLA-DPA1", 1, "HLA-DPA1*01:03"),
        ("p1", "hap1", "ctg1", "HLA-DPB1", 1, "HLA-DPB1*04:01"),
        ("p1", "hap2", "ctg9", "HLA-DPA1", 1, "HLA-DPA1*02:01"),
        ("p1", "hap2", "ctg9", "HLA-DPB1", 1, "HLA-DPB1*02:01"),
    ])
    haps, _ = m.extract_haplotypes(t1, "DPA1", "DPB1")
    check("one person contributes two independent phased haplotypes", len(haps) == 2, str(haps))


# ---------------------------------------------------------------------------
# disclosure
# ---------------------------------------------------------------------------
def test_suppression_rule():
    check("zero is disclosable", m.suppress(0) == "0")
    check("1-19 is suppressed", m.suppress(19) == "<20", m.suppress(19))
    check("20 and above is disclosed", m.suppress(20) == "20")


def main():
    print("test_hla_ld_by_ancestry.py")
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
