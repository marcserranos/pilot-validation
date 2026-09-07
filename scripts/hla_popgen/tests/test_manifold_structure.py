#!/usr/bin/env python3
"""Real assertions against 14_manifold_structure.py's pure-math functions, per SCHEMA.md hard rule
#6 ("every script here must run against synthetic fixtures before it is handed to Marc").

These are unit tests against small hand-built numpy/pandas inputs (not the full fixture pipeline --
the functions here operate on an already-built dosage matrix + ancestry labels, the same shape
`08_embedding_compare.py`'s `embedding_cache.pkl` provides, so there is nothing pipeline-specific
left to exercise). All of them run with only numpy/pandas/matplotlib -- no scipy or umap-learn
needed locally, matching this repo's HAVE_UMAP-guarded-import convention.

Run:
    python3 scripts/hla_popgen/tests/test_manifold_structure.py
"""
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)


def _load_module(filename, modname):
    path = os.path.join(HLA_POPGEN_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ms = _load_module("14_manifold_structure.py", "hla_popgen_manifold_structure")

FAILED = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILED.append(name)


def test_clean_ancestry():
    print("test_clean_ancestry")
    # Real cache dicts hand back pd.NA (not plain None) for people missing an ancestry call --
    # this crashed np.unique/`==`/`or` in the real VM run (TypeError: boolean value of NA is
    # ambiguous) because none of the local synthetic tests below ever built an array containing
    # pd.NA until this regression test was added.
    check("pd.NA maps to None", ms.clean_ancestry(pd.NA) is None)
    check("None maps to None", ms.clean_ancestry(None) is None)
    check("float nan maps to None", ms.clean_ancestry(float("nan")) is None)
    check("a real label passes through as str", ms.clean_ancestry("AFR") == "AFR")
    # The actual failure mode: building a labels array with a pd.NA mixed in, then using it the
    # way ancestry_f_stat/residualize_ancestry/the KNN test all do. `==` against a cleaned array
    # is safe; np.unique on a None-containing object array is NOT (None vs str has no ordering),
    # which is exactly why ancestry_f_stat filters None out before calling np.unique rather than
    # after -- this checks both halves of that pattern.
    mixed = np.array([ms.clean_ancestry(v) for v in ["AFR", pd.NA, "EUR", None, "AFR"]],
                     dtype=object)
    try:
        _ = mixed == "AFR"
        eq_ok = True
    except TypeError:
        eq_ok = False
    check("== on a cleaned mixed-None array doesn't crash", eq_ok)
    try:
        uniq = np.unique(mixed[np.array([v is not None for v in mixed])])
        filtered_unique_ok = True
    except TypeError:
        filtered_unique_ok = False
    check("np.unique on the None-filtered array doesn't crash", filtered_unique_ok,
          f"uniq={list(uniq) if filtered_unique_ok else None}")


def test_ancestry_f_stat():
    print("test_ancestry_f_stat")
    rng = np.random.default_rng(0)
    n = 300
    ancestry = np.array((["AFR"] * (n // 3)) + (["EUR"] * (n // 3)) + (["EAS"] * (n // 3)))
    # Column with a real between-group mean shift.
    col_diff = np.concatenate([rng.normal(0, 0.5, n // 3), rng.normal(3, 0.5, n // 3),
                                rng.normal(6, 0.5, n // 3)])
    # Column with no relationship to ancestry at all.
    col_same = rng.normal(0, 0.5, n)
    f_diff = ms.ancestry_f_stat(col_diff, ancestry)
    f_same = ms.ancestry_f_stat(col_same, ancestry)
    check("differentiated column has much higher F-stat than flat column",
          f_diff > 20 * max(f_same, 1e-6), f"f_diff={f_diff:.2f} f_same={f_same:.2f}")
    # <2 usable groups (one singleton ancestry) should degrade gracefully to 0.0, not crash.
    ancestry_thin = np.array(["AFR"] * (n - 1) + ["EAS"])
    f_thin = ms.ancestry_f_stat(col_diff, ancestry_thin)
    check("single-member group is dropped, doesn't crash", f_thin == 0.0, f"f_thin={f_thin}")


def test_spearman_corr():
    print("test_spearman_corr")
    rng = np.random.default_rng(1)
    x = rng.normal(size=200)
    y_monotonic = x ** 3 + rng.normal(0, 0.01, 200)  # monotonic-ish, not linear
    y_random = rng.normal(size=200)
    r_mono = ms.spearman_corr(x, y_monotonic)
    r_rand = ms.spearman_corr(x, y_random)
    check("monotonically related series has high rank correlation", r_mono > 0.9,
          f"r_mono={r_mono:.3f}")
    check("unrelated series has low rank correlation", abs(r_rand) < 0.2,
          f"r_rand={r_rand:.3f}")
    const = np.zeros(50)
    check("constant input doesn't crash (returns nan)",
          np.isnan(ms.spearman_corr(const, rng.normal(size=50))))


def test_residualize_ancestry():
    print("test_residualize_ancestry")
    rng = np.random.default_rng(2)
    n, p = 150, 5
    ancestry_by_person = {}
    ids = [f"p{i}" for i in range(n)]
    groups = (["AFR"] * 50) + (["EUR"] * 50) + (["EAS"] * 50)
    for pid, g in zip(ids, groups):
        ancestry_by_person[pid] = g
    ancestry_by_person[ids[0]] = pd.NA  # real cache dicts carry pd.NA for missing ancestry calls
    # Column 0 has a strong ancestry-driven offset; column 1 has none.
    X = rng.normal(0, 1, (n, p))
    offset = np.array([0.0 if g == "AFR" else (5.0 if g == "EUR" else -5.0) for g in groups])
    X[:, 0] += offset
    mat = pd.DataFrame(X, index=ids, columns=[f"col{i}" for i in range(p)])

    resid = ms.residualize_ancestry(mat, ancestry_by_person)
    check("residual has same shape as input", resid.shape == mat.shape)

    group_means_before = mat.groupby([ancestry_by_person[i] for i in mat.index])["col0"].mean()
    group_means_after = resid.groupby([ancestry_by_person[i] for i in resid.index])["col0"].mean()
    check("per-ancestry-group mean of the ancestry-driven column collapses toward 0 after "
          "residualizing", group_means_after.abs().max() < 1e-9,
          f"before={group_means_before.to_dict()} after={group_means_after.to_dict()}")
    check("overall variance of the ancestry-driven column drops substantially",
          resid["col0"].var() < 0.3 * mat["col0"].var(),
          f"before_var={mat['col0'].var():.2f} after_var={resid['col0'].var():.2f}")
    check("a column unrelated to ancestry is only mean-centered per group, not distorted",
          resid["col1"].var() > 0.7 * mat["col1"].var(),
          f"before_var={mat['col1'].var():.2f} after_var={resid['col1'].var():.2f}")


def test_knn_label_enrichment():
    print("test_knn_label_enrichment")
    rng = np.random.default_rng(3)
    n_per = 100
    cluster_a = rng.normal(loc=[0, 0], scale=0.3, size=(n_per, 2))
    cluster_b = rng.normal(loc=[10, 10], scale=0.3, size=(n_per, 2))
    coords = np.vstack([cluster_a, cluster_b])
    labels_matching = np.array([0] * n_per + [1] * n_per)  # perfectly aligned with clusters
    rng.shuffle(coords)  # shuffling both together keeps alignment; do it via index instead
    # (undo: re-derive aligned arrays properly instead of the shuffle above)
    coords = np.vstack([cluster_a, cluster_b])
    idx_all = ms.knn_indices(coords, k=10)
    obs, null_mean, null_std, p = ms.knn_label_enrichment(idx_all, labels_matching, n_perm=200,
                                                          seed=0)
    check("labels aligned with spatial clusters show near-total same-label enrichment",
          obs > 0.95, f"obs={obs:.3f}")
    check("enrichment far exceeds the permutation null", obs > null_mean + 5 * max(null_std, 1e-6),
          f"obs={obs:.3f} null_mean={null_mean:.3f} null_std={null_std:.4f}")
    check("p-value is significant for real spatial structure", p < 0.01, f"p={p:.4g}")

    labels_random = rng.integers(0, 2, size=2 * n_per)  # unrelated to spatial position
    obs_r, null_mean_r, null_std_r, p_r = ms.knn_label_enrichment(idx_all, labels_random,
                                                                  n_perm=200, seed=0)
    check("labels unrelated to spatial position land within the null band",
          abs(obs_r - null_mean_r) < 5 * max(null_std_r, 1e-6),
          f"obs_r={obs_r:.3f} null_mean_r={null_mean_r:.3f} null_std_r={null_std_r:.4f}")
    check("p-value is non-significant for unrelated labels", p_r > 0.05, f"p_r={p_r:.4g}")


def test_variance_loadings_report():
    print("test_variance_loadings_report")
    rng = np.random.default_rng(4)
    n, p = 200, 20
    ids = [f"p{i}" for i in range(n)]
    groups = (["AFR"] * (n // 2)) + (["EUR"] * (n // 2))
    ancestry_by_person = dict(zip(ids, groups))
    ancestry_by_person[ids[-1]] = pd.NA  # missing-ancestry-call case, must not crash np.unique
    X = rng.integers(0, 2, size=(n, p)).astype(float)
    cols = [f"HLA-B*{i:02d}" for i in range(p)]
    # Make columns 0-2 strongly ancestry-differentiated (like real HLA background variation) but
    # keep a little within-group noise so ss_within isn't exactly 0 (degenerate F-stat guard).
    for j in range(3):
        X[:n // 2, j] = rng.normal(0, 0.05, n // 2)
        X[n // 2:, j] = 1 + rng.normal(0, 0.05, n - n // 2)
    mat = pd.DataFrame(X, index=ids, columns=cols)
    disease_cols_by_label = {"fake disease": [cols[5], cols[6]]}

    df, corr = ms.variance_loadings_report(mat, ancestry_by_person, disease_cols_by_label,
                                            pca_denoise_k=5)
    check("returns one row per matrix column", len(df) == p, f"len(df)={len(df)}")
    check("disease-proxy flag matches the columns passed in",
          set(df.loc[df["is_disease_proxy"], "column"]) == {cols[5], cols[6]})
    check("ancestry-differentiated columns get a much higher F-stat than the rest",
          df.loc[df["column"].isin(cols[:3]), "ancestry_f_stat"].mean() >
          10 * df.loc[~df["column"].isin(cols[:3]), "ancestry_f_stat"].mean())
    check("correlation is a finite float (or nan on degenerate input, never a crash)",
          isinstance(corr, float))


if __name__ == "__main__":
    test_clean_ancestry()
    test_ancestry_f_stat()
    test_spearman_corr()
    test_residualize_ancestry()
    test_knn_label_enrichment()
    test_variance_loadings_report()
    print()
    if FAILED:
        print(f"FAILED: {len(FAILED)} check(s): {FAILED}")
        sys.exit(1)
    print("All checks passed.")
