#!/usr/bin/env python3
"""Tests for scripts/hla_popgen/27_allele_space_coverage.py.

Runs standalone (no pytest, per repo convention -- see e.g. tests/test_coverage.py):

    cd scripts/hla_popgen && python3 -c "import sys,tempfile,pathlib,inspect; \\
        sys.path.insert(0,'tests'); import test_allele_space_coverage as t; \\
        [f(**({'tmp_path': pathlib.Path(tempfile.mkdtemp())} \\
              if 'tmp_path' in inspect.signature(f).parameters else {})) \\
         for n,f in inspect.getmembers(t, inspect.isfunction) if n.startswith('test_')]; \\
        print('ok')"

Also collectible by pytest if it's ever installed (plain `assert`-based test_* functions).
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)
if HLA_POPGEN_DIR not in sys.path:
    sys.path.insert(0, HLA_POPGEN_DIR)

import _coverage as cov  # noqa: F401 -- imported for parity with test_coverage.py's style

# The module file is "27_allele_space_coverage.py" -- not a valid bare import name
# (leading digit), so load it explicitly by path, matching this repo's convention for
# numerically-prefixed script modules (see e.g. test_novelty_by_field.py's importlib use).
import importlib.util

_SPEC = importlib.util.spec_from_file_location(
    "s27", os.path.join(HLA_POPGEN_DIR, "27_allele_space_coverage.py"))
m27 = importlib.util.module_from_spec(_SPEC)
sys.modules["s27"] = m27
_SPEC.loader.exec_module(m27)

REQUIRED_COLUMNS = m27.REQUIRED_COLUMNS


# ---------------------------------------------------------------------------
# Synthetic fixture builder
# ---------------------------------------------------------------------------
def _rows_for_gene_resolution(gene, resolution, ancestry_counts, exclude_flagged="clean_only",
                               ancestry_scheme="strict", unrelated_only=True):
    """ancestry_counts: {ancestry: {allele_id: n_haplotypes}}."""
    rows = []
    for anc, counts in ancestry_counts.items():
        for allele_id, n in counts.items():
            rows.append({"gene": gene, "resolution": resolution,
                         "exclude_flagged": exclude_flagged,
                         "ancestry_scheme": ancestry_scheme, "ancestry": anc,
                         "unrelated_only": unrelated_only, "allele_id": allele_id,
                         "n_haplotypes": n})
    return rows


def make_synthetic_counts():
    rows = []

    # Gene HLA-A: EUR all singletons (low coverage), AFR all doubletons (high coverage).
    eur_counts = {f"HLA-A_prot_singleton{i:03d}": 1 for i in range(30)}
    afr_counts = {f"HLA-A*{i:02d}:01": 2 for i in range(15)}
    for resolution in ("protein", "cds"):
        rows += _rows_for_gene_resolution(
            "HLA-A", resolution, {"EUR": eur_counts, "AFR": afr_counts})

    # Gene HLA-B: EUR and AFR share a common allele pool plus each has private alleles,
    # for the cross-ancestry test (private alleles => diagonal >= off-diagonal).
    shared = {f"HLA-B*{i:02d}:01": 20 for i in range(5)}
    eur_private = {f"HLA-B_prot_eur_priv{i:03d}": 10 for i in range(10)}
    afr_private = {f"HLA-B_prot_afr_priv{i:03d}": 10 for i in range(10)}
    eur_b = dict(shared)
    eur_b.update(eur_private)
    afr_b = dict(shared)
    afr_b.update(afr_private)
    for resolution in ("protein", "cds"):
        rows += _rows_for_gene_resolution("HLA-B", resolution, {"EUR": eur_b, "AFR": afr_b})

    # Gene HLA-C: three ancestries, each with mostly-private alleles, for design-scenario
    # tests (equal allocation should beat a 94%-EUR design on the cross-ancestry mean).
    eur_c = {f"HLA-C*{i:02d}:01": 50 for i in range(3)}
    eur_c.update({f"HLA-C_prot_eur_priv{i:03d}": 5 for i in range(30)})
    afr_c = {f"HLA-C_prot_afr_priv{i:03d}": 5 for i in range(60)}
    eas_c = {f"HLA-C_prot_eas_priv{i:03d}": 5 for i in range(40)}
    for resolution in ("protein", "cds"):
        rows += _rows_for_gene_resolution(
            "HLA-C", resolution, {"EUR": eur_c, "AFR": afr_c, "EAS": eas_c})

    return pd.DataFrame(rows, columns=REQUIRED_COLUMNS)


def write_synthetic_counts(path):
    df = make_synthetic_counts()
    df.to_csv(path, sep="\t", index=False)
    return df


# ---------------------------------------------------------------------------
# 1. loading / fatal-missing-column behavior
# ---------------------------------------------------------------------------
def test_load_allele_counts_ok(tmp_path):
    path = str(tmp_path / "counts.tsv")
    write_synthetic_counts(path)
    df = m27.load_allele_counts(path)
    for c in REQUIRED_COLUMNS:
        assert c in df.columns
    assert df["n_haplotypes"].dtype.kind in "if"
    assert set(df["ancestry"].unique()) >= {"EUR", "AFR", "EAS"}


def test_load_allele_counts_missing_column_exits(tmp_path):
    path = str(tmp_path / "bad_counts.tsv")
    df = make_synthetic_counts().drop(columns=["allele_id"])
    df.to_csv(path, sep="\t", index=False)
    try:
        m27.load_allele_counts(path)
        raised = False
    except SystemExit as e:
        raised = True
        assert "allele_id" in str(e)
    assert raised, "expected SystemExit listing the missing column"


def test_missing_file_exits():
    try:
        m27.load_allele_counts("/nonexistent/path/does_not_exist.tsv")
        raised = False
    except SystemExit:
        raised = True
    assert raised


# ---------------------------------------------------------------------------
# 2. singletons -> low coverage, doubletons -> high coverage
# ---------------------------------------------------------------------------
def test_singletons_low_doubletons_high(tmp_path):
    path = str(tmp_path / "counts.tsv")
    df = write_synthetic_counts(path)
    scoped = m27.filter_scope(df, "protein", "clean_only", "strict", True)
    table = m27.build_coverage_table(scoped, "protein", n_boot=20)
    eur_row = table[(table["gene"] == "HLA-A") & (table["ancestry"] == "EUR")].iloc[0]
    afr_row = table[(table["gene"] == "HLA-A") & (table["ancestry"] == "AFR")].iloc[0]
    assert eur_row["f1"] == 30 and eur_row["f2"] == 0
    assert afr_row["f1"] == 0 and afr_row["f2"] == 15  # doubletons only -> f1=0, f2=all alleles
    assert eur_row["coverage"] < 0.5, "all-singleton spectrum should give low coverage"
    assert afr_row["coverage"] == 1.0, "no singletons (f1=0) -> sample_coverage is exactly 1.0"
    assert afr_row["coverage"] > eur_row["coverage"]


# ---------------------------------------------------------------------------
# 3. cross-ancestry diagonal >= off-diagonal with private alleles
# ---------------------------------------------------------------------------
def test_cross_ancestry_diagonal_ge_offdiagonal(tmp_path):
    path = str(tmp_path / "counts.tsv")
    df = write_synthetic_counts(path)
    scoped = m27.filter_scope(df, "protein", "clean_only", "strict", True)
    matrix = m27.build_cross_ancestry_matrix(scoped, "protein")
    b = matrix[(matrix["gene"] == "HLA-B")]
    diag = b[b["is_diagonal"]]["cross_coverage"]
    offdiag = b[~b["is_diagonal"]]["cross_coverage"]
    assert len(diag) and len(offdiag)
    assert diag.min() >= offdiag.max() - 1e-9, (
        f"diagonal {diag.tolist()} should be >= off-diagonal {offdiag.tolist()} "
        f"when populations carry private alleles")


# ---------------------------------------------------------------------------
# 4. design putting everything into one ancestry maximizes that ancestry's coverage
# ---------------------------------------------------------------------------
def test_design_all_in_one_ancestry_maximizes_its_coverage():
    counts_by_pop = {
        "EUR": {f"e{i}": 5 for i in range(20)},
        "AFR": {f"a{i}": 5 for i in range(20)},
    }
    weights = {"EUR": 0.5, "AFR": 0.5}
    full_eur_sizes = {"EUR": 100, "AFR": 0}
    spread_sizes = {"EUR": 20, "AFR": 80}
    cov_full, _ = m27.pooled_design_coverage_safe(counts_by_pop, full_eur_sizes, weights)
    cov_spread, _ = m27.pooled_design_coverage_safe(counts_by_pop, spread_sizes, weights)
    assert cov_full["EUR"] >= cov_spread["EUR"], (
        "putting the whole design budget into EUR's own subsample should not decrease "
        "EUR's own coverage relative to a design that spends less of the budget on EUR")


# ---------------------------------------------------------------------------
# 5. equal allocation beats 94%-EUR design on the mean, with private alleles
# ---------------------------------------------------------------------------
def test_equal_allocation_beats_eur_heavy_design(tmp_path):
    path = str(tmp_path / "counts.tsv")
    df = write_synthetic_counts(path)
    scoped = m27.filter_scope(df, "protein", "clean_only", "strict", True)
    design_df = m27.build_design_scenarios(scoped, "protein")
    c = design_df[(design_df["gene"] == "HLA-C") & (design_df["target_ancestry"] == "GLOBAL_MEAN")]
    equal_cov = c[c["design"] == "equal_allocation"]["coverage"].iloc[0]
    eur_heavy_cov = c[c["design"] == "94pct_eur_ukb_like"]["coverage"].iloc[0]
    assert equal_cov > eur_heavy_cov, (
        f"equal allocation ({equal_cov}) should beat the 94%-EUR design ({eur_heavy_cov}) "
        f"on mean coverage when AFR/EAS carry mostly-private alleles starved by the "
        f"EUR-heavy design")


# ---------------------------------------------------------------------------
# 6. suppression
# ---------------------------------------------------------------------------
def test_suppression_applied():
    assert m27.suppress(5) == "<20"
    assert m27.suppress(19) == "<20"
    assert m27.suppress(20) == "20"
    assert m27.suppress(0) == "0"
    df = pd.DataFrame({"n_haplotypes": [3, 25, 19, 20]})
    out = m27.suppress_df(df, ["n_haplotypes"])
    assert list(out["n_haplotypes"]) == ["<20", "25", "<20", "20"]


# ---------------------------------------------------------------------------
# 7. IMGT-known-allele convention
# ---------------------------------------------------------------------------
def test_is_known_allele_convention():
    assert m27.is_known_allele("HLA-A*02:01")
    assert m27.is_known_allele("HLA-B*07:02:01")
    assert not m27.is_known_allele("A_prot_deadbeef")
    assert not m27.is_known_allele("A_cds_deadbeef")
    assert not m27.is_known_allele("A_nonfunc_deadbeef")
    assert not m27.is_known_allele(None)


# ---------------------------------------------------------------------------
# 8. end-to-end run writes every file
# ---------------------------------------------------------------------------
EXPECTED_OUTPUTS = [
    "coverage_by_gene_ancestry.tsv",
    "coverage_curves.tsv",
    "coverage_curves_protein.png",
    "coverage_curves_cds.png",
    "cross_ancestry_matrix.tsv",
    "cross_ancestry_matrix.png",
    "imgt_coverage_today.tsv",
    "imgt_coverage_today.png",
    "design_scenarios.tsv",
    "design_scenarios.png",
    "summary.json",
    "allele_space_coverage_report.md",
]


def test_end_to_end(tmp_path):
    counts_path = str(tmp_path / "allele_counts_by_resolution.tsv")
    write_synthetic_counts(counts_path)
    out_dir = str(tmp_path / "out")
    argv = ["--allele-counts", counts_path, "--out-dir", out_dir, "--n-boot", "10"]
    m27.main(argv)
    for fname in EXPECTED_OUTPUTS:
        fpath = os.path.join(out_dir, fname)
        assert os.path.exists(fpath), f"missing expected output {fname}"
        assert os.path.getsize(fpath) > 0, f"{fname} is empty"

    # No allele-level rows anywhere: the allele_id column must never appear in any
    # committed TSV (only aggregate gene/ancestry/design columns).
    for fname in EXPECTED_OUTPUTS:
        if fname.endswith(".tsv"):
            with open(os.path.join(out_dir, fname)) as f:
                header = f.readline()
            assert "allele_id" not in header, f"{fname} leaks allele-level column"

    # Suppression made it into the written coverage table for thin cells.
    cov_tsv = pd.read_csv(os.path.join(out_dir, "coverage_by_gene_ancestry.tsv"), sep="\t",
                           dtype=str)
    assert cov_tsv["n_haplotypes"].isin(["<20"]).any() or (
        pd.to_numeric(cov_tsv["n_haplotypes"], errors="coerce").min() >= 20), (
        "expected at least one suppressed thin cell, or all cells to be >= 20")

    # --limit-genes actually limits scope.
    out_dir2 = str(tmp_path / "out_limited")
    m27.main(["--allele-counts", counts_path, "--out-dir", out_dir2, "--n-boot", "10",
              "--limit-genes", "1"])
    cov2 = pd.read_csv(os.path.join(out_dir2, "coverage_by_gene_ancestry.tsv"), sep="\t")
    assert cov2["gene"].nunique() <= 1
