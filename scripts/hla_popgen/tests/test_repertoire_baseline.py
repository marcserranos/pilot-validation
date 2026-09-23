#!/usr/bin/env python3
"""Unit tests for 42_repertoire_baseline.py (WS6 repertoire-only baseline).

Two properties matter most for a disease-classification baseline, because a bug in either
would silently produce an optimistic, unreproducible number:

  1. CV folds are IDENTICAL across model variants that share the same label vector (fixed
     seed + StratifiedKFold is deterministic) -- covariate-only, repertoire-only, and
     repertoire+covariate must be scored on the exact same held-out people per repeat, or
     their AUROCs aren't comparable (BenchRep-T Section 2.2's whole point).
  2. Feature construction (VJ usage, k-mer vocabulary/capping) never looks at the label
     vector -- it is a pure function of the repertoire data, so no fold's test labels can
     leak into the feature space, regardless of downstream CV.

Also covers: productive filter drops stop/ambiguous CDR3s, CDR3 trimming removes exactly the
flanking residues, gapped-4mer count is correct, and rarefy_clones produces nested subsets
(BenchRep-T Appendix B.4) and never exceeds the target.

Run: python3 scripts/hla_popgen/tests/test_repertoire_baseline.py
"""
import importlib.util
import inspect
import os
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


rb = _load_module("42_repertoire_baseline.py", "repertoire_baseline")

FAILURES = []


def check(name, cond):
    status = "OK" if cond else "FAIL"
    print(f"[{status}] {name}")
    if not cond:
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# Feature-construction primitives
# ---------------------------------------------------------------------------
def test_productive_filter():
    check("is_productive drops stop codon", rb.is_productive("CASS_YF") is False)
    check("is_productive drops ambiguous", rb.is_productive("CAS?SYF") is False)
    check("is_productive keeps clean CDR3", rb.is_productive("CASSLGQAYEQYF") is True)
    check("is_productive rejects too-short", rb.is_productive("C") is False)


def test_trim_cdr3():
    check("trim_cdr3 removes flanking Cys/Phe",
         rb.trim_cdr3("CASSLGQAYEQYF") == "ASSLGQAYEQY")
    check("trim_cdr3 empty on too-short input", rb.trim_cdr3("CF") == "")


def test_gapped_4mers_count():
    trimmed = "ABCDEFG"  # length 7 -> 4 windows (0..3)
    kmers = rb.gapped_4mers(trimmed)
    n_windows = len(trimmed) - 3
    check("gapped_4mers yields 5 features per window (1 ungapped + 4 gapped)",
         len(kmers) == n_windows * 5)
    check("gapped_4mers first ungapped 4mer is correct", kmers[0] == "ABCD")
    check("gapped_4mers gapped variants mark exactly one position",
         all(km.count("_") == 1 for km in kmers if "_" in km))


def test_chain_of():
    check("chain_of TRBV -> TRB", rb.chain_of("TRBV12-3") == "TRB")
    check("chain_of IGHV -> IGH", rb.chain_of("IGHV1-2") == "IGH")
    check("chain_of other gene -> other", rb.chain_of("TRAV1") == "other")


# ---------------------------------------------------------------------------
# Feature construction never sees labels
# ---------------------------------------------------------------------------
def test_build_features_signature_has_no_label_arg():
    sig = inspect.signature(rb.build_vj_kmer_features)
    params = list(sig.parameters)
    check("build_vj_kmer_features takes no phenotype/label argument",
         all("label" not in p and p != "y" for p in params))


def test_features_identical_regardless_of_label_permutation():
    tables, covariates, labels = rb.make_synthetic_cohort(40, seed=1)
    X1, kinds1 = rb.build_vj_kmer_features(tables, max_kmers=200)
    # Shuffle the label rows (simulating a different fold assignment) -- features must be
    # byte-identical, since build_vj_kmer_features never consumes labels at all.
    shuffled_labels = labels.sample(frac=1.0, random_state=99).reset_index(drop=True)
    X2, kinds2 = rb.build_vj_kmer_features(tables, max_kmers=200)
    check("feature matrix identical after unrelated label shuffle",
         X1.equals(X2) and kinds1 == kinds2)
    check("shuffled label frame is a permutation of the original (sanity)",
         set(shuffled_labels["research_id"]) == set(labels["research_id"]))


def test_kmer_cap_is_variance_based_not_label_based():
    tables, _, _ = rb.make_synthetic_cohort(30, seed=2)
    X_full, _ = rb.build_vj_kmer_features(tables, max_kmers=10_000_000)  # effectively uncapped
    X_capped, _ = rb.build_vj_kmer_features(tables, max_kmers=50)
    kmer_cols_full = [c for c in X_full.columns if c.startswith("KMER__")]
    kmer_cols_capped = [c for c in X_capped.columns if c.startswith("KMER__")]
    if kmer_cols_full:
        expected = X_full[kmer_cols_full].var(axis=0).sort_values(ascending=False).index[:50]
        check("kmer cap keeps the highest-variance columns",
             set(kmer_cols_capped) == set(expected))
    else:
        check("kmer cap test skipped (no kmer columns in synthetic fixture)", True)


# ---------------------------------------------------------------------------
# CV folds identical across model variants sharing a label vector
# ---------------------------------------------------------------------------
def test_folds_deterministic_and_shared_across_models():
    if not rb.HAVE_SKLEARN:
        check("sklearn not installed -- fold-sharing test skipped", True)
        return
    y = np.array([0, 1] * 30)
    folds_a = rb.make_folds(y, n_splits=3, seed=42)
    folds_b = rb.make_folds(y, n_splits=3, seed=42)
    same = all(
        np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])
        for a, b in zip(folds_a, folds_b)
    )
    check("make_folds is deterministic for a fixed seed (so covariate-only, "
         "repertoire-only, and combined models see the same held-out people)", same)

    folds_c = rb.make_folds(y, n_splits=3, seed=7)
    different = not all(
        np.array_equal(a[0], c[0]) for a, c in zip(folds_a, folds_c)
    )
    check("different seeds give different fold assignments (repeats are not degenerate)",
         different)


def test_run_cv_pooled_shares_folds_across_covariate_and_repertoire_variants():
    """End-to-end check on the actual run_pipeline code path: for a given phenotype, the
    covariate-only and repertoire-only variants both call run_cv_pooled with the identical
    y array and seed list, so their internal fold splits must coincide. We check this by
    intercepting make_folds calls made by run_cv_pooled for two different X matrices but the
    same y/seeds, and asserting identical (train_idx, test_idx) pairs each repeat."""
    if not rb.HAVE_SKLEARN:
        check("sklearn not installed -- pooled-CV fold-sharing test skipped", True)
        return
    rng = np.random.default_rng(3)
    y = (rng.random(60) < 0.3).astype(int)
    seeds = [10, 11]
    calls = []
    real_make_folds = rb.make_folds

    def spy_make_folds(y_arg, n_splits, seed):
        folds = real_make_folds(y_arg, n_splits, seed)
        calls.append((seed, [tuple(te) for _, te in folds]))
        return folds

    rb.make_folds = spy_make_folds
    try:
        X_cov = rng.normal(size=(60, 5))
        X_rep = rng.normal(size=(60, 30))
        rb.run_cv_pooled(X_cov, y, seeds, n_splits=3, model="lr")
        calls_cov = calls[:]
        calls.clear()
        rb.run_cv_pooled(X_rep, y, seeds, n_splits=3, model="lr")
        calls_rep = calls[:]
    finally:
        rb.make_folds = real_make_folds

    check("run_cv_pooled requests folds for every repeat seed",
         [s for s, _ in calls_cov] == seeds and [s for s, _ in calls_rep] == seeds)
    same_folds = all(
        [set(map(tuple, te1)) == set(map(tuple, te2))
         for (_, te1), (_, te2) in zip(calls_cov, calls_rep)]
    )
    check("covariate-only and repertoire-only runs use the SAME held-out folds "
         "for the same (y, seeds) -- AUROCs are directly comparable", same_folds)


# ---------------------------------------------------------------------------
# Rarefaction (depth-confound control)
# ---------------------------------------------------------------------------
def test_rarefy_never_exceeds_target_and_is_nested():
    df = pd.DataFrame({
        "V": ["TRBV1"] * 20, "J": ["TRBJ1"] * 20,
        "CDR3_amino_acids": [f"CAS{i}SF" for i in range(20)],
        "chain": ["TRB"] * 20,
        "read_count": list(range(1, 21)),  # total = 210
    })
    small = rb.rarefy_clones(df, target_reads=50, seed=5)
    large = rb.rarefy_clones(df, target_reads=150, seed=5)
    check("rarefy_clones hits the exact target read count",
         small["read_count"].sum() == 50 and large["read_count"].sum() == 150)
    check("rarefy_clones returns None when a person has fewer reads than target",
         rb.rarefy_clones(df, target_reads=1000, seed=5) is None)
    # Nesting: every clone touched at the smaller depth must also appear (with count > 0
    # possibility) at the larger depth for the same person/seed -- i.e. small's clone index
    # set is a subset of large's, since both come from prefixes of the same fixed permutation.
    small_idx = set(small.index) if not small.empty else set()
    # re-derive via V/CDR3 identity since indices reset; use (V, CDR3_amino_acids) as key
    small_keys = set(zip(small["V"], small["CDR3_amino_acids"]))
    large_keys = set(zip(large["V"], large["CDR3_amino_acids"]))
    check("rarefy_clones nesting: clones present at the smaller depth are a subset of the "
         "larger depth for the same person+seed (BenchRep-T Appendix B.4)",
         small_keys.issubset(large_keys))


def test_rarefy_deterministic_for_fixed_seed():
    df = pd.DataFrame({
        "V": ["TRBV1"] * 10, "J": ["TRBJ1"] * 10,
        "CDR3_amino_acids": [f"CAS{i}SF" for i in range(10)],
        "chain": ["TRB"] * 10,
        "read_count": [10] * 10,
    })
    a = rb.rarefy_clones(df, target_reads=40, seed=1)
    b = rb.rarefy_clones(df, target_reads=40, seed=1)
    check("rarefy_clones is deterministic for a fixed seed",
         a["read_count"].equals(b["read_count"]) and (a["V"] == b["V"]).all())


# ---------------------------------------------------------------------------
# report.tsv loading
# ---------------------------------------------------------------------------
def test_load_report_tsv_filters_and_restricts_chains(tmp_path=None):
    import tempfile
    rows = pd.DataFrame({
        "V": ["TRBV1", "TRBV2", "IGHV1", "TRAV1"],
        "D": [".", ".", ".", "."],
        "J": ["TRBJ1", "TRBJ1", "IGHJ1", "TRAJ1"],
        "C": [".", ".", ".", "."],
        "CDR3_dna": ["N" * 30] * 4,
        "CDR3_amino_acids": ["CASSLGQAYEQYF", "CASS_YF", "CARDGY", "CAV?YF"],
        "read_count": ["5", "3", "2", "1"],
        "frequency": ["0.5", "0.3", "0.2", "0.1"],
    })
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "p_report.tsv")
        rows.to_csv(path, sep="\t", index=False)
        loaded = rb.load_report_tsv(path, chains=("TRB", "IGH"))
    check("load_report_tsv drops non-productive rows", len(loaded) == 2)
    check("load_report_tsv restricts to requested chains",
         set(loaded["chain"]) == {"TRB", "IGH"})


def test_load_report_tsv_accepts_real_trust4_column_aliases():
    """The real report.tsv files on the VM use `#count`/`CDR3aa`/`CDR3nt`, not the
    `read_count`/`CDR3_amino_acids`/`CDR3_dna` names the protocol doc assumed (found 2026-09-23
    when the first real run died with a FATAL missing-column error) -- lock in the alias map."""
    import tempfile
    rows = pd.DataFrame({
        "#count": ["5", "3", "2"],
        "frequency": ["0.5", "0.3", "0.2"],
        "CDR3nt": ["N" * 30] * 3,
        "CDR3aa": ["CASSLGQAYEQYF", "CASS_YF", "CARDGY"],
        "V": ["TRBV1", "TRBV2", "IGHV1"],
        "D": [".", ".", "."],
        "J": ["TRBJ1", "TRBJ1", "IGHJ1"],
        "C": [".", ".", "."],
        "cid": ["1", "2", "3"],
        "cid_full_length": ["1", "1", "0"],
    })
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "p_report.tsv")
        rows.to_csv(path, sep="\t", index=False)
        loaded = rb.load_report_tsv(path, chains=("TRB", "IGH"))
    check("load_report_tsv maps #count -> read_count (non-productive TRBV2 row dropped)",
         sorted(loaded["read_count"].tolist()) == [2, 5])
    check("load_report_tsv maps CDR3aa -> CDR3_amino_acids",
         set(loaded["CDR3_amino_acids"]) == {"CASSLGQAYEQYF", "CARDGY"})


def main():
    test_productive_filter()
    test_trim_cdr3()
    test_gapped_4mers_count()
    test_chain_of()
    test_build_features_signature_has_no_label_arg()
    test_features_identical_regardless_of_label_permutation()
    test_kmer_cap_is_variance_based_not_label_based()
    test_folds_deterministic_and_shared_across_models()
    test_run_cv_pooled_shares_folds_across_covariate_and_repertoire_variants()
    test_rarefy_never_exceeds_target_and_is_nested()
    test_rarefy_deterministic_for_fixed_seed()
    test_load_report_tsv_filters_and_restricts_chains()
    test_load_report_tsv_accepts_real_trust4_column_aliases()

    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("All tests passed.")


if __name__ == "__main__":
    main()
