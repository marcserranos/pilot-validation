#!/usr/bin/env python3
"""Unit tests for scripts/hla_popgen/_disclosure.py -- the shared PUBLIC/INTERNAL disclosure
policy module for S04 (reference/AOU_SMALL_CELL_POLICY.md, context/DECISIONS.md "two disclosure
versions"). Synthetic values only, no real cohort data.

Run: pytest scripts/hla_popgen/tests/test_disclosure.py -q
"""
import importlib.util
import os
import sys

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)
if HLA_POPGEN_DIR not in sys.path:
    sys.path.insert(0, HLA_POPGEN_DIR)

import _disclosure as d


# ---------------------------------------------------------------------------
# mask(): each count_type x mode
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("count_type", [d.PARTICIPANT, d.CARRIER_NAMED_ALLELE, d.RATE_NUM_DENOM])
def test_masked_types_public_censors_1_to_19(count_type):
    for n in range(1, 20):
        r = d.mask(n, count_type, d.PUBLIC)
        assert r.display == "<20"
        assert r.lt20 is True
        assert r.review_flag is False


@pytest.mark.parametrize("count_type", [d.PARTICIPANT, d.CARRIER_NAMED_ALLELE, d.RATE_NUM_DENOM])
def test_masked_types_public_keeps_true_zero(count_type):
    r = d.mask(0, count_type, d.PUBLIC)
    assert r.display == "0"
    assert r.lt20 is False


@pytest.mark.parametrize("count_type", [d.PARTICIPANT, d.CARRIER_NAMED_ALLELE, d.RATE_NUM_DENOM])
def test_masked_types_public_exact_at_or_above_20(count_type):
    for n in (20, 21, 500):
        r = d.mask(n, count_type, d.PUBLIC)
        assert r.display == str(n)
        assert r.lt20 is False


@pytest.mark.parametrize("count_type", [d.PARTICIPANT, d.CARRIER_NAMED_ALLELE, d.RATE_NUM_DENOM])
def test_masked_types_internal_always_exact_with_lt20_flag(count_type):
    for n in (0, 1, 5, 19, 20, 100):
        r = d.mask(n, count_type, d.INTERNAL)
        assert r.display == str(n)
        assert r.lt20 == (0 < n < 20)
        assert r.review_flag is False


@pytest.mark.parametrize("count_type", [d.ALLELE_DISTINCT, d.RICHNESS, d.QC_TALLY])
@pytest.mark.parametrize("mode", [d.PUBLIC, d.INTERNAL])
def test_exact_types_always_exact_in_both_modes(count_type, mode):
    for n in (0, 1, 5, 19, 20, 500):
        r = d.mask(n, count_type, mode)
        assert r.display == str(n)
        assert r.lt20 == (0 < n < 20)
        assert r.review_flag is False


def test_recurrence_class_exact_in_both_modes():
    for mode in (d.PUBLIC, d.INTERNAL):
        for n in (0, 1, 19, 20):
            r = d.mask(n, d.RECURRENCE_CLASS, mode, total_carriers=1000)
            assert r.display == str(n)


def test_recurrence_class_review_flag_public_only():
    # thin stratum (< floor) -> flagged, but ONLY in PUBLIC mode
    r_public_thin = d.mask(3, d.RECURRENCE_CLASS, d.PUBLIC, total_carriers=50)
    assert r_public_thin.review_flag is True
    r_public_wide = d.mask(3, d.RECURRENCE_CLASS, d.PUBLIC, total_carriers=1000)
    assert r_public_wide.review_flag is False
    r_internal_thin = d.mask(3, d.RECURRENCE_CLASS, d.INTERNAL, total_carriers=50)
    assert r_internal_thin.review_flag is False  # review_flag is a PUBLIC-only concept


def test_recurrence_class_review_flag_conservative_default_when_total_missing():
    r = d.mask(3, d.RECURRENCE_CLASS, d.PUBLIC)  # no total_carriers given
    assert r.review_flag is True


def test_mask_missing_value_is_blank_in_every_type_and_mode():
    for count_type in d.COUNT_TYPES:
        for mode in d.MODES:
            assert d.mask(None, count_type, mode).display == ""
            assert d.mask(float("nan"), count_type, mode).display == ""


def test_mask_rejects_unknown_count_type_or_mode():
    with pytest.raises(ValueError):
        d.mask(5, "not_a_type", d.PUBLIC)
    with pytest.raises(ValueError):
        d.mask(5, d.PARTICIPANT, "not_a_mode")


def test_mask_rejects_negative_value():
    with pytest.raises(ValueError):
        d.mask(-1, d.PARTICIPANT, d.PUBLIC)


# ---------------------------------------------------------------------------
# mask_rate()
# ---------------------------------------------------------------------------
def test_mask_rate_public_blanks_on_small_numerator_or_denominator():
    r = d.mask_rate(5, 30, d.PUBLIC)
    assert r["n"] == "<20" and r["d"] == "30"
    assert r["pct"] == "" and r["lo"] == "" and r["hi"] == ""
    assert r["lt20_n"] is True and r["lt20_d"] is False

    r2 = d.mask_rate(25, 10, d.PUBLIC)
    assert r2["n"] == "25" and r2["d"] == "<20"
    assert r2["pct"] == ""


def test_mask_rate_public_blanks_on_zero_denominator():
    r = d.mask_rate(0, 0, d.PUBLIC)
    assert r["pct"] == ""


def test_mask_rate_public_reports_real_zero_numerator():
    r = d.mask_rate(0, 100, d.PUBLIC)
    assert r["n"] == "0" and r["d"] == "100"
    assert r["pct"] == "0.0"


def test_mask_rate_public_exact_when_both_ge20():
    r = d.mask_rate(40, 100, d.PUBLIC)
    assert r["n"] == "40" and r["d"] == "100"
    assert r["pct"] == "40.0"
    assert 0.0 < float(r["lo"]) < float(r["hi"]) < 100.0


def test_mask_rate_internal_never_blanks():
    r = d.mask_rate(5, 30, d.INTERNAL)
    assert r["n"] == "5" and r["d"] == "30"
    assert r["pct"] != ""
    assert r["lt20_n"] is True and r["lt20_d"] is False

    r2 = d.mask_rate(0, 0, d.INTERNAL)
    assert r2["n"] == "0" and r2["d"] == "0"
    assert r2["pct"] == ""  # still undefined for d==0, but not "blanked for disclosure"


# ---------------------------------------------------------------------------
# INTERNAL path refusal (task item 2)
# ---------------------------------------------------------------------------
def test_internal_refuses_repo_working_tree_path():
    repo_path = os.path.join(HLA_POPGEN_DIR, "some_internal_output.tsv")
    with pytest.raises(ValueError):
        d.assert_internal_path_allowed(repo_path)


def test_internal_refuses_reports_path_even_outside_repo(tmp_path):
    reports_path = tmp_path / "reports" / "hla_popgen" / "foo.tsv"
    with pytest.raises(ValueError):
        d.assert_internal_path_allowed(str(reports_path))


def test_internal_allows_vm_style_path():
    d.assert_internal_path_allowed(os.path.expanduser("~/s04/internal/44/recurrence_classes.tsv"))


def test_write_tsv_internal_refuses_repo_path(tmp_path, monkeypatch):
    df = pd.DataFrame({"a": [1, 2]})
    bad_path = os.path.join(HLA_POPGEN_DIR, "should_not_write.tsv")
    with pytest.raises(ValueError):
        d.write_tsv(df, bad_path, d.INTERNAL)
    assert not os.path.exists(bad_path)


def test_write_tsv_internal_prepends_header_comment(tmp_path):
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    out_path = tmp_path / "internal" / "x.tsv"
    d.write_tsv(df, str(out_path), d.INTERNAL)
    text = out_path.read_text()
    lines = text.splitlines()
    assert lines[0] == d.INTERNAL_HEADER_COMMENT
    # body still parses as a normal TSV once the header comment line is skipped
    reread = pd.read_csv(out_path, sep="\t", skiprows=1)
    assert list(reread.columns) == ["a", "b"]
    assert reread["a"].tolist() == [1, 2]


def test_write_tsv_public_has_no_header_comment(tmp_path):
    df = pd.DataFrame({"a": [1, 2]})
    out_path = tmp_path / "public" / "x.tsv"
    d.write_tsv(df, str(out_path), d.PUBLIC)
    text = out_path.read_text()
    assert not text.startswith("#")
    reread = pd.read_csv(out_path, sep="\t")
    assert reread["a"].tolist() == [1, 2]


def test_default_out_dir_routes_internal_away_from_public_default():
    pub = d.default_out_dir(d.PUBLIC, "44", "/some/public/dir")
    assert pub == "/some/public/dir"
    internal = d.default_out_dir(d.INTERNAL, "44", "/some/public/dir")
    assert "internal" in internal
    assert "44" in internal
    assert internal != "/some/public/dir"


# ---------------------------------------------------------------------------
# End-to-end: a synthetic PUBLIC table never contains a bare 1-19 value in a
# participant/carrier/rate column.
# ---------------------------------------------------------------------------
def test_public_output_table_has_no_1_to_19_values_in_disclosive_columns():
    rows = []
    rng = np.random.default_rng(0)
    for _ in range(200):
        participant_n = int(rng.integers(0, 60))
        allele_n = int(rng.integers(0, 60))
        rows.append({
            "n_carriers": d.mask(participant_n, d.PARTICIPANT, d.PUBLIC).display,
            "n_distinct_alleles": d.mask(allele_n, d.ALLELE_DISTINCT, d.PUBLIC).display,
        })
    df = pd.DataFrame(rows)

    def _is_1_to_19(cell):
        try:
            v = int(cell)
        except (TypeError, ValueError):
            return False
        return 0 < v < 20

    assert not df["n_carriers"].map(_is_1_to_19).any()
    # allele_distinct is intentionally exact -- 1-19 values ARE allowed here
    assert df["n_distinct_alleles"].map(_is_1_to_19).any()


def test_internal_output_table_lt20_flag_matches_masked_public_cells():
    rng = np.random.default_rng(1)
    for _ in range(100):
        n = int(rng.integers(0, 60))
        pub = d.mask(n, d.PARTICIPANT, d.PUBLIC)
        internal = d.mask(n, d.PARTICIPANT, d.INTERNAL)
        assert internal.display == str(n)
        assert internal.lt20 == (pub.display == "<20")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
