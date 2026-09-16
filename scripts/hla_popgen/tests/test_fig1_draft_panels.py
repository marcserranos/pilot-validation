#!/usr/bin/env python3
"""Fixture tests for 23_fig1_draft_panels.py -- the three-bucket split, the per-ancestry rate
denominators, the SFS binning, and the disclosure floor on participant counts."""
import importlib
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
mod = importlib.import_module("23_fig1_draft_panels")


def _novel():
    rows = [
        # gene, novelty_class, tier, n_haplotypes, n_persons, ancestry_counts
        ("HLA-A", "protein_altering", "recurrent", 2, 2, {"AFR": 2}),
        ("HLA-A", "beyond_cds", "high", 30, 25, {"AFR": 10, "EUR": 15}),
        ("HLA-A", "protein_altering", "flagged_artifact", 1, 1, {"EUR": 1}),
        ("HLA-A", "synonymous", "singleton_clean", 1, 1, {"EAS": 1}),
        ("HLA-A", "undetermined", "singleton_clean", 1, 1, {"EAS": 1}),
    ]
    df = pd.DataFrame([{
        "novel_id": f"n{i}", "gene": g, "nearest_allele": "HLA-A*01:01", "novelty_class": nc,
        "confidence_tier": t, "cds_distance": 1.0, "n_aa_changes": 1.0, "n_haplotypes": h,
        "n_persons": p, "n_ancestries": len(ac), "ancestry_counts": json.dumps(ac),
    } for i, (g, nc, t, h, p, ac) in enumerate(rows)])
    return df


def _richness():
    rows = []
    for anc, n in [("AFR", 100), ("AMR", 50), ("EAS", 50), ("EUR", 200), ("MID", 10),
                   ("SAS", 40), ("POOLED", 450)]:
        rows.append({"gene": "HLA-A", "category": "all", "ancestry": anc, "n_haplotypes": n})
    rows.append({"gene": "HLA-A", "category": "novel", "ancestry": "POOLED", "n_haplotypes": 35})
    return pd.DataFrame(rows)


def _loaded(tmp_path):
    p = tmp_path / "novel.tsv"
    _novel().to_csv(p, sep="\t", index=False)
    return mod.load_novel(str(p))


def test_buckets_keep_artifacts_and_noncoding_apart(tmp_path):
    novel = _loaded(tmp_path)
    assert list(novel["bucket"][:4]) == ["cds_changing", "noncoding_only", "flagged_artifact",
                                         "cds_changing"]
    assert pd.isna(novel["bucket"].iloc[4])  # undetermined stays out of every bucket


def test_rates_use_per_ancestry_haplotype_denominators(tmp_path):
    rates = mod.novel_rate_table(_loaded(tmp_path), _richness(), ["HLA-A"]).set_index(
        ["ancestry", "bucket"])
    assert rates.loc[("AFR", "cds_changing"), "per_100_haplotypes"] == 2.0
    assert rates.loc[("EUR", "noncoding_only"), "per_100_haplotypes"] == 7.5
    assert rates.loc[("EUR", "flagged_artifact"), "per_100_haplotypes"] == 0.5
    # a bucket with no carriers in an ancestry is a real zero, not a missing row
    assert rates.loc[("MID", "cds_changing"), "per_100_haplotypes"] == 0.0
    assert len(rates) == 6 * 3


def test_reconciliation_counts_undetermined_separately(tmp_path):
    recon = mod.reconciliation_table(_loaded(tmp_path), _richness(), ["HLA-A"])
    r = recon.loc["HLA-A"]
    assert (r["cds_changing"], r["noncoding_only"], r["flagged_artifact"], r["undetermined"]) == \
        (3, 30, 1, 1)
    assert r["sum_all_buckets"] == 35 == r["richness_novel_haplotypes"]


def test_sfs_bins_are_left_closed():
    s = pd.Series({1: 5, 2: 3, 4: 1, 5: 2, 8: 1, 9: 4, 1024: 1, 1025: 7})
    b = mod.bin_sfs(s)
    assert b["1"] == 5 and b["2"] == 3 and b["3-4"] == 1 and b["5-8"] == 3
    assert b["9-16"] == 4 and b["513-1024"] == 1 and b[">1024"] == 7
    assert b.sum() == s.sum()
    assert len(mod.SFS_BIN_WIDTHS) == len(mod.SFS_BIN_LABELS)


def test_participant_counts_below_floor_are_suppressed(tmp_path):
    lst = mod.top_novel_list(_loaded(tmp_path), ["HLA-A"], "noncoding_only", min_persons=11)
    row = lst.iloc[0]
    assert row["n_persons"] == "25"
    assert row["n_AFR"] == "<20" and row["n_EUR"] == "<20" and row["n_MID"] == "0"
    assert mod.suppress(0) == "0" and mod.suppress(19) == "<20" and mod.suppress(20) == "20"


def test_top_list_excludes_singletons_and_artifacts(tmp_path):
    lst = mod.top_novel_list(_loaded(tmp_path), ["HLA-A"], "cds_changing", min_persons=1)
    assert list(lst["novel_id"]) == ["n0"]
