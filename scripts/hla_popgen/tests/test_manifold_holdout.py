#!/usr/bin/env python3
"""Real assertions against 15_manifold_holdout_validation.py's pure-math functions, per SCHEMA.md
hard rule #6 ("every script here must run against synthetic fixtures before it is handed to Marc").

Unit tests against small hand-built numpy/pandas inputs -- same rationale as
test_manifold_structure.py: these functions operate on an already-built embedding + labels, so
there's nothing pipeline-specific left to exercise, and everything here runs with only
numpy/pandas -- no umap-learn needed locally (the actual UMAP fit/transform calls are
HAVE_UMAP-guarded and only exercised on the VM).

Run:
    python3 scripts/hla_popgen/tests/test_manifold_holdout.py
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


mh = _load_module("15_manifold_holdout_validation.py", "hla_popgen_manifold_holdout")

FAILED = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILED.append(name)


def test_stratified_split():
    print("test_stratified_split")
    rng = np.random.default_rng(0)
    n = 1000
    labels = (rng.random(n) < 0.09).astype(int)  # ~9% prevalence, like the real any_diag rate
    train_idx, test_idx = mh.stratified_split(labels, test_frac=0.3, seed=1)
    check("train + test cover every index exactly once",
          set(train_idx) | set(test_idx) == set(range(n)) and
          len(set(train_idx) & set(test_idx)) == 0)
    check("test split is roughly the requested fraction",
          abs(len(test_idx) / n - 0.3) < 0.02, f"got {len(test_idx) / n:.3f}")
    train_prev = labels[train_idx].mean()
    test_prev = labels[test_idx].mean()
    overall_prev = labels.mean()
    check("prevalence is preserved in both splits (stratified, not just random)",
          abs(train_prev - overall_prev) < 0.01 and abs(test_prev - overall_prev) < 0.02,
          f"overall={overall_prev:.3f} train={train_prev:.3f} test={test_prev:.3f}")


def test_knn_predict_score():
    print("test_knn_predict_score")
    rng = np.random.default_rng(1)
    n_per = 150
    train_coords = np.vstack([rng.normal([0, 0], 0.3, (n_per, 2)),
                              rng.normal([10, 10], 0.3, (n_per, 2))])
    train_labels = np.array([0] * n_per + [1] * n_per)
    query_near_neg = rng.normal([0, 0], 0.3, (20, 2))
    query_near_pos = rng.normal([10, 10], 0.3, (20, 2))
    scores_neg = mh.knn_predict_score(train_coords, train_labels, query_near_neg, k=15)
    scores_pos = mh.knn_predict_score(train_coords, train_labels, query_near_pos, k=15)
    check("queries near the negative cluster score close to 0",
          scores_neg.mean() < 0.1, f"mean={scores_neg.mean():.3f}")
    check("queries near the positive cluster score close to 1",
          scores_pos.mean() > 0.9, f"mean={scores_pos.mean():.3f}")


def test_auroc():
    print("test_auroc")
    rng = np.random.default_rng(2)
    labels = np.array([0] * 100 + [1] * 100)
    perfect_scores = np.concatenate([rng.uniform(0, 0.4, 100), rng.uniform(0.6, 1.0, 100)])
    reversed_scores = np.concatenate([rng.uniform(0.6, 1.0, 100), rng.uniform(0, 0.4, 100)])
    random_scores = rng.uniform(0, 1, 200)
    check("perfectly separating scores give AUROC ~= 1.0",
          mh.auroc(perfect_scores, labels) > 0.99,
          f"auroc={mh.auroc(perfect_scores, labels):.4f}")
    check("perfectly reversed scores give AUROC ~= 0.0",
          mh.auroc(reversed_scores, labels) < 0.01,
          f"auroc={mh.auroc(reversed_scores, labels):.4f}")
    r = mh.auroc(random_scores, labels)
    check("unrelated scores give AUROC near 0.5", 0.35 < r < 0.65, f"auroc={r:.4f}")
    check("single-class labels return nan, not a crash", np.isnan(mh.auroc(random_scores,
                                                                            np.ones(200))))


def test_bootstrap_auc_ci():
    print("test_bootstrap_auc_ci")
    rng = np.random.default_rng(3)
    labels = np.array([0] * 100 + [1] * 100)
    scores = np.concatenate([rng.uniform(0, 0.4, 100), rng.uniform(0.6, 1.0, 100)])
    point, lo, hi = mh.bootstrap_auc_ci(scores, labels, n_boot=300, seed=4)
    check("point estimate matches plain auroc()", abs(point - mh.auroc(scores, labels)) < 1e-9)
    check("CI bounds bracket the point estimate", lo <= point <= hi, f"lo={lo:.3f} hi={hi:.3f}")
    check("CI lower bound is well above 0.5 for a strongly separating signal", lo > 0.8,
          f"lo={lo:.3f}")


def test_auroc_permutation_pvalue():
    print("test_auroc_permutation_pvalue")
    rng = np.random.default_rng(5)
    labels = np.array([0] * 100 + [1] * 100)
    strong_scores = np.concatenate([rng.uniform(0, 0.4, 100), rng.uniform(0.6, 1.0, 100)])
    obs, null_mean, null_std, p = mh.auroc_permutation_pvalue(strong_scores, labels, n_perm=500,
                                                               seed=6)
    check("null mean sits near 0.5 for a permutation null", abs(null_mean - 0.5) < 0.05,
          f"null_mean={null_mean:.3f}")
    check("strong signal gets a significant p-value", p < 0.01, f"p={p:.4g}")

    rng2 = np.random.default_rng(7)
    random_scores = rng2.uniform(0, 1, 200)
    _obs_r, _null_mean_r, _null_std_r, p_r = mh.auroc_permutation_pvalue(random_scores, labels,
                                                                         n_perm=500, seed=8)
    check("unrelated signal gets a non-significant p-value", p_r > 0.05, f"p_r={p_r:.4g}")


if __name__ == "__main__":
    test_stratified_split()
    test_knn_predict_score()
    test_auroc()
    test_bootstrap_auc_ci()
    test_auroc_permutation_pvalue()
    print()
    if FAILED:
        print(f"FAILED: {len(FAILED)} check(s): {FAILED}")
        sys.exit(1)
    print("All checks passed.")
