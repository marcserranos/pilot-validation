#!/usr/bin/env python3
"""Unit tests for 32_novel_allele_callouts.py.

This script turns S01's cluster table into a list someone will paste into a manuscript, so the
failure modes are all "a wrong number looks completely normal":

  1. `<20` read as 0 would erase real carriers; read as 20 would invent them and let a censored
     allele outrank a disclosable one in the main-text ranking. It must be NaN.
  2. A synonymous cluster (same protein, new CDS) must never end up in a "novel proteins" count.
  3. Ancestry concentration must be computed over disclosable cells only, and the censored
     ancestries must be carried through, or the list will overstate how population-specific an
     allele is.
  4. The catalogue-gap control must be attached per gene -- without it, "TAP1 tops the list" reads
     as biology when it may be reference depth.

Run: python3 scripts/hla_popgen/tests/test_novel_allele_callouts.py
"""
import importlib.util
import math
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


m = _load_module("32_novel_allele_callouts.py", "novel_callouts")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# censoring
# ---------------------------------------------------------------------------
def test_censored_cell_is_nan_not_zero_and_not_twenty():
    v = m.parse_count("<20")
    check("'<20' parses to NaN", math.isnan(v), str(v))
    check("...which is not 0", v != 0)
    check("...and not 20", v != 20)
    check("is_censored recognises it", m.is_censored("<20"))
    check("a plain number is not censored", not m.is_censored("42"))


def test_plain_counts_parse():
    check("a plain integer parses", m.parse_count("94") == 94.0)
    check("zero parses to zero, not NaN", m.parse_count("0") == 0.0)
    check("NA parses to NaN", math.isnan(m.parse_count("NA")))
    check("an empty cell parses to NaN", math.isnan(m.parse_count("")))


def test_censored_allele_cannot_outrank_a_disclosable_one():
    """The ranking guard. A cluster with a censored count must never be sorted above one with a
    real count of 20, which is what would happen if '<20' were read as 20 (ties) or worse."""
    clusters = pd.DataFrame([
        {"cluster_id": "c_big", "cluster_type": "novel_protein", "gene": "TAP1",
         "n_persons_unrelated_clean": "20", "clean_recurrent_unrelated": "True"},
        {"cluster_id": "c_censored", "cluster_type": "novel_protein", "gene": "TAP1",
         "n_persons_unrelated_clean": "<20", "clean_recurrent_unrelated": "True"},
    ])
    _, tier1, tier2 = m.build_tiers(clusters, {}, {})
    check("only the disclosable cluster reaches tier 1", list(tier1["cluster_id"]) == ["c_big"],
          str(list(tier1["cluster_id"])))
    check("the censored one falls to tier 2", list(tier2["cluster_id"]) == ["c_censored"],
          str(list(tier2["cluster_id"])))


def test_non_recurrent_is_excluded_from_both_tiers():
    clusters = pd.DataFrame([
        {"cluster_id": "singleton", "cluster_type": "novel_protein", "gene": "TAP1",
         "n_persons_unrelated_clean": "50", "clean_recurrent_unrelated": "False"},
    ])
    _, tier1, tier2 = m.build_tiers(clusters, {}, {})
    check("a non-recurrent cluster is in neither tier", tier1.empty and tier2.empty,
          f"t1={len(tier1)} t2={len(tier2)}")


# ---------------------------------------------------------------------------
# ancestry
# ---------------------------------------------------------------------------
def test_ancestry_profile_separates_disclosable_from_censored():
    row = {"n_persons_strict_AFR": "80", "n_persons_strict_EUR": "<20",
           "n_persons_strict_EAS": "0", "n_persons_strict_AMR": "NA",
           "n_persons_strict_MID": "0", "n_persons_strict_SAS": "0"}
    counts, censored = m.ancestry_profile(row)
    check("disclosable non-zero counts are kept", counts == {"AFR": 80.0}, str(counts))
    check("a censored ancestry is recorded, not dropped silently", censored == ["EUR"],
          str(censored))


def test_concentration_uses_disclosable_counts_only():
    share, top = m.concentration({"AFR": 90.0, "EUR": 10.0})
    check("the top ancestry is identified", top == "AFR")
    check("the share is over disclosable carriers", abs(share - 0.9) < 1e-12, str(share))
    share2, top2 = m.concentration({})
    check("no disclosable carriers gives NaN, not 1.0", math.isnan(share2), str(share2))
    check("...and no top ancestry", top2 is None)


# ---------------------------------------------------------------------------
# the reference-depth control
# ---------------------------------------------------------------------------
def test_catalogue_gap_is_haplotype_weighted(tmpdir=None):
    import tempfile
    d = tempfile.mkdtemp()
    p = os.path.join(d, "imgt.tsv")
    # 1000 haplotypes at 90% catalogued + 1000 at 100% -> gap = 5%, not the unweighted 5% by luck:
    # use unequal n to make weighting observable.
    pd.DataFrame([
        {"gene": "TAP1", "resolution": "protein", "ancestry": "AFR", "n_haplotypes": 3000,
         "frac_in_imgt_today": 0.90},
        {"gene": "TAP1", "resolution": "protein", "ancestry": "EUR", "n_haplotypes": 1000,
         "frac_in_imgt_today": 1.00},
    ]).to_csv(p, sep="\t", index=False)
    gap = m.catalogue_gap_by_gene(p)
    # weighted: 1 - (0.9*3000 + 1.0*1000)/4000 = 1 - 0.925 = 0.075 -> 7.5%
    check("the catalogue gap is weighted by haplotype count, not averaged over ancestries",
          abs(gap["TAP1"] - 7.5) < 1e-6, str(gap))


def test_catalogue_gap_is_attached_to_each_cluster():
    clusters = pd.DataFrame([
        {"cluster_id": "c1", "cluster_type": "novel_protein", "gene": "TAP1",
         "n_persons_unrelated_clean": "426", "clean_recurrent_unrelated": "True"},
    ])
    _, tier1, _ = m.build_tiers(clusters, {}, {"TAP1": 8.0})
    check("the gene's catalogue gap travels with the callout",
          tier1.iloc[0]["pct_gene_haplotypes_not_in_imgt"] == 8.0,
          str(tier1.iloc[0].to_dict()))


def test_missing_catalogue_gap_is_none_not_zero():
    """A missing control must not read as 'this gene has no catalogue gap'."""
    clusters = pd.DataFrame([
        {"cluster_id": "c1", "cluster_type": "novel_protein", "gene": "HLA-A",
         "n_persons_unrelated_clean": "30", "clean_recurrent_unrelated": "True"},
    ])
    _, tier1, _ = m.build_tiers(clusters, {}, {})
    check("an absent gap is None, never 0", tier1.iloc[0]["pct_gene_haplotypes_not_in_imgt"] is None,
          str(tier1.iloc[0]["pct_gene_haplotypes_not_in_imgt"]))


# ---------------------------------------------------------------------------
# rates
# ---------------------------------------------------------------------------
def test_carrier_rate_uses_people_not_haplotypes():
    clusters = pd.DataFrame([
        {"cluster_id": "c1", "cluster_type": "novel_protein", "gene": "TAP1",
         "n_persons_unrelated_clean": "100", "clean_recurrent_unrelated": "True",
         "n_persons_strict_AFR": "100"},
    ])
    # 4000 haplotypes -> 2000 people -> 100/2000 = 50 per 1000
    _, tier1, _ = m.build_tiers(clusters, {"TAP1": {"AFR": 2000.0}}, {})
    check("the rate is carriers per 1,000 people, with haplotypes halved upstream",
          abs(tier1.iloc[0]["carriers_per_1000_in_top_ancestry"] - 50.0) < 1e-6,
          str(tier1.iloc[0]["carriers_per_1000_in_top_ancestry"]))


def test_load_denominators_halves_haplotypes(tmpdir=None):
    import tempfile
    d = tempfile.mkdtemp()
    p = os.path.join(d, "cov.tsv")
    pd.DataFrame([{"gene": "HLA-A", "resolution": "protein", "ancestry": "AFR",
                   "n_haplotypes": 4294}]).to_csv(p, sep="\t", index=False)
    den = m.load_denominators(p)
    check("haplotypes are converted to people", den["A"]["AFR"] == 2147.0, str(den))
    check("the HLA- prefix is stripped so gene keys match", "A" in den, str(den.keys()))


def main():
    print("test_novel_allele_callouts.py")
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
