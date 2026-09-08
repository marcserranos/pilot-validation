#!/usr/bin/env python3
"""Hypothesis #1 from `research/ANCESTRY_VS_DISEASE_MANIFOLD.md`'s "New hypotheses" section:
`14_manifold_structure.py`'s supervised-UMAP result (obs=0.97 same-label KNN fraction vs a 0.83
null, p=0.002) is a real finding but an incomplete one -- `umap.UMAP(..., y=labels)` is explicitly
optimized to pull same-label points together, so a high score on the SAME data it was fit and
labeled on is expected almost by construction. It shows the allele-dosage space *contains*
disease-relevant information; it does not show that information *generalizes* to people the fit
never saw. This script is the direct fix: a proper train/held-out-test split.

## Design

1. Stratified train/test split (default 70/30) on `any_diag`, preserving prevalence in both halves.
2. Fit UMAP TWICE on the TRAIN split only, using train-only mean/std (no leakage from test people
   into standardization or the fit):
   - `supervised`: `y=y_train` passed to `umap.UMAP.fit`.
   - `unsupervised`: same train data, no `y` -- the baseline. If supervision isn't actually adding
     predictive information beyond whatever structure already exists unsupervised, this baseline
     should score just as well.
3. `.transform()` the TEST split into each fitted embedding (a real out-of-sample projection --
   umap-learn's `transform` never sees test labels, only test allele-dosage rows).
4. For each test person, score = fraction of their k=15 nearest TRAIN neighbors (in embedding
   space) who are diagnosed. This is a k-NN classifier operating entirely on the frozen train-fit
   embedding.
5. Evaluate with AUROC (rank-based, no sklearn) between predicted score and true test-set
   diagnosis, plus:
   - a bootstrap 95% CI (resampling test people) for estimation uncertainty, and
   - a label-permutation p-value (shuffling TEST labels only, scores fixed -- valid because the
     scores never depended on test labels in the first place, so this is a cheap, honest null:
     "would an AUROC this high arise if these held-out labels were unrelated to the scores?").

If `supervised`'s held-out AUROC clears its bootstrap CI above 0.5 AND beats `unsupervised`'s, that
is real evidence the signal generalizes, not an artifact of supervised UMAP's training-time
optimization. If it collapses to ~0.5 out of sample, the original supervised-UMAP finding was
overfit to the labels it was given, and the honest conclusion reverts to "no evidence of disease
structure" pending better features (see checklist items 2-4 in the research doc).

Usage (fixtures -- pure-math functions only, no umap-learn needed locally):
    python3 scripts/hla_popgen/tests/test_manifold_holdout.py

Real run (VM, after 08 and 12 have both run at least once):
    python3 scripts/hla_popgen/15_manifold_holdout_validation.py
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _viz_common as vc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


# ---------------------------------------------------------------------------
# Split
# ---------------------------------------------------------------------------
def stratified_split(labels, test_frac=0.3, seed=0):
    """Returns (train_idx, test_idx) -- positions into `labels`, each label value split
    independently so both halves preserve the input prevalence (important here: any_diag is only
    ~9% positive, and an unstratified split could easily leave the test set with very few or zero
    positives)."""
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels)
    train_idx, test_idx = [], []
    for lbl in np.unique(labels):
        grp = np.where(labels == lbl)[0]
        rng.shuffle(grp)
        n_test = int(round(len(grp) * test_frac))
        test_idx.extend(grp[:n_test])
        train_idx.extend(grp[n_test:])
    return np.array(sorted(train_idx)), np.array(sorted(test_idx))


# ---------------------------------------------------------------------------
# UMAP fit/transform (train-only stats, real out-of-sample projection)
# ---------------------------------------------------------------------------
def fit_umap_train(X_train_std, y, seed, umap_neighbors, umap_min_dist):
    """y=None fits unsupervised; y=array fits supervised. Returns (reducer, train_coords), or
    (None, None) if umap-learn isn't installed (same HAVE_UMAP guard used throughout this repo)."""
    if not vc.HAVE_UMAP:
        return None, None
    import umap
    reducer = umap.UMAP(n_neighbors=umap_neighbors, min_dist=umap_min_dist, random_state=seed)
    coords = reducer.fit_transform(X_train_std, y=y) if y is not None else \
        reducer.fit_transform(X_train_std)
    return reducer, coords


# ---------------------------------------------------------------------------
# k-NN-in-embedding classifier (train neighbors only -- true held-out scoring)
# ---------------------------------------------------------------------------
def knn_predict_score(train_coords, train_labels, query_coords, k=15, chunk=500):
    """For each query (test) point, the fraction of its k nearest TRAIN-embedding neighbors that
    are label==1. A continuous [0,1] score, evaluated against true query labels via AUROC below."""
    train_labels = np.asarray(train_labels).astype(float)
    n_q = len(query_coords)
    k = min(k, len(train_coords))
    scores = np.empty(n_q)
    for start in range(0, n_q, chunk):
        end = min(start + chunk, n_q)
        d = np.linalg.norm(query_coords[start:end, None, :] - train_coords[None, :, :], axis=2)
        idx = np.argpartition(d, k - 1, axis=1)[:, :k]
        scores[start:end] = train_labels[idx].mean(axis=1)
    return scores


# ---------------------------------------------------------------------------
# AUROC + uncertainty (no sklearn/scipy: rank-based Mann-Whitney U formula)
# ---------------------------------------------------------------------------
def auroc(scores, labels):
    scores = np.asarray(scores)
    labels = np.asarray(labels).astype(int)
    n_pos = int(labels.sum())
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = pd.Series(scores).rank().values
    sum_ranks_pos = ranks[labels == 1].sum()
    return float((sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def bootstrap_auc_ci(scores, labels, n_boot=1000, seed=0, alpha=0.05):
    """Resamples TEST people (not train) with replacement; skips a resample if it happens to be
    single-class (AUROC undefined). Returns (point_estimate, ci_lo, ci_hi)."""
    rng = np.random.default_rng(seed)
    scores, labels = np.asarray(scores), np.asarray(labels)
    n = len(scores)
    point = auroc(scores, labels)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        b_labels = labels[idx]
        if 0 < b_labels.sum() < len(b_labels):
            boots.append(auroc(scores[idx], b_labels))
    if not boots:
        return point, np.nan, np.nan
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return point, float(lo), float(hi)


def auroc_permutation_pvalue(scores, labels, n_perm=2000, seed=0):
    """Shuffles TEST labels only (scores never depended on them, so this is a cheap, honest null --
    no re-fitting/re-transforming needed, unlike a label-permutation test on the UMAP fit itself
    would require). Precomputes ranks once since scores are fixed across permutations. Returns
    (observed, null_mean, null_std, p) -- p is one-sided (null >= observed)."""
    scores, labels = np.asarray(scores), np.asarray(labels).astype(int)
    n_pos = int(labels.sum())
    n_neg = len(labels) - n_pos
    observed = auroc(scores, labels)
    if n_pos == 0 or n_neg == 0:
        return observed, np.nan, np.nan, np.nan
    ranks = pd.Series(scores).rank().values
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for p in range(n_perm):
        perm_labels = rng.permutation(labels)
        sum_ranks_pos = ranks[perm_labels == 1].sum()
        null[p] = (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    p_value = (float((null >= observed).sum()) + 1) / (n_perm + 1)
    return observed, float(null.mean()), float(null.std()), p_value


# ---------------------------------------------------------------------------
# Plotting -- where do held-out people actually land in a mapping they never trained on?
# ---------------------------------------------------------------------------
LABEL_COLORS = {0: "#9a9a9a", 1: "#1b9e77"}  # neutral gray (undiagnosed) / teal (diagnosed)


def plot_train_test_panel(ax, train_coords, y_train, test_coords, y_test, title):
    """Faint small dots = train people (the fit never optimizes anything about the held-out
    people). Larger black-edged triangles = test people, `.transform()`-projected into that same
    frozen mapping using ONLY their allele dosage -- their color (true diagnosis) was never seen
    by the fit or the transform. If diagnosed test triangles land inside the diagnosed train
    region rather than scattered randomly, that is the visual version of the held-out AUROC."""
    for lbl in (0, 1):
        m = y_train == lbl
        ax.scatter(train_coords[m, 0], train_coords[m, 1], s=6, alpha=0.25,
                   color=LABEL_COLORS[lbl], linewidths=0)
    for lbl in (0, 1):
        m = y_test == lbl
        ax.scatter(test_coords[m, 0], test_coords[m, 1], s=26, alpha=0.95,
                   color=LABEL_COLORS[lbl], edgecolors="black", linewidths=0.5, marker="^")
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])


def build_train_test_figure(panels, out_dir):
    """panels: {name: (train_coords, y_train, test_coords, y_test)}."""
    if not panels:
        return None
    fig, axes = plt.subplots(1, len(panels), figsize=(5.5 * len(panels), 5))
    if len(panels) == 1:
        axes = [axes]
    for ax, (name, (tr_c, tr_y, te_c, te_y)) in zip(axes, panels.items()):
        plot_train_test_panel(ax, tr_c, tr_y, te_c, te_y,
                              f"{name}\n(small = train, triangle = held-out test)")
    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=LABEL_COLORS[0],
                  markersize=6, alpha=0.5, label="train, undiagnosed"),
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=LABEL_COLORS[1],
                  markersize=6, alpha=0.5, label="train, diagnosed"),
        plt.Line2D([0], [0], marker="^", color="none", markerfacecolor=LABEL_COLORS[0],
                  markeredgecolor="black", markersize=8, label="held-out test, undiagnosed"),
        plt.Line2D([0], [0], marker="^", color="none", markerfacecolor=LABEL_COLORS[1],
                  markeredgecolor="black", markersize=8, label="held-out test, diagnosed"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8, frameon=False,
              bbox_to_anchor=(0.5, -0.02))
    path = os.path.join(out_dir, "holdout_train_test_embedding.png")
    vc.savefig(fig, path)
    return path


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--embedding-cache",
                     default=os.path.expanduser("~/repos/pilot-validation/reports/hla_popgen/"
                                                 "08_embedding_compare/lr/embedding_cache.pkl"))
    ap.add_argument("--wide-diagnosis-labels",
                     default=os.path.expanduser("~/pipeline_outputs/rnaseq/pheno/"
                                                 "person_disease_labels_wide.tsv"))
    ap.add_argument("--out-dir", default=None,
                     help="Default: reports/hla_popgen/15_manifold_holdout_validation/")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--test-frac", type=float, default=0.3)
    ap.add_argument("--umap-neighbors", type=int, default=15)
    ap.add_argument("--umap-min-dist", type=float, default=0.1)
    ap.add_argument("--knn-k", type=int, default=15)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--n-perm", type=int, default=2000)
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(os.path.dirname(here))
    out_dir = args.out_dir or os.path.join(repo_root, "reports", "hla_popgen",
                                            "15_manifold_holdout_validation")
    os.makedirs(out_dir, exist_ok=True)

    import pickle
    print(f"Loading {args.embedding_cache!r} ...", file=sys.stderr)
    with open(os.path.expanduser(args.embedding_cache), "rb") as f:
        cache = pickle.load(f)
    mat_full = pd.DataFrame(cache["mat_full_values"], index=cache["mat_full_index"],
                             columns=cache["mat_full_columns"])

    print(f"Loading {args.wide_diagnosis_labels!r} ...", file=sys.stderr)
    diag = pd.read_csv(os.path.expanduser(args.wide_diagnosis_labels), sep="\t",
                        dtype={"person_id": str}).set_index("person_id")
    diag = diag.reindex(mat_full.index).fillna(False)
    any_diag = (diag[diag.columns].sum(axis=1) > 0).astype(int).values

    print(f"Splitting {mat_full.shape[0]} people {1 - args.test_frac:.0%}/{args.test_frac:.0%} "
          f"train/test, stratified on any_diag ...", file=sys.stderr)
    train_idx, test_idx = stratified_split(any_diag, test_frac=args.test_frac, seed=args.seed)
    y_train, y_test = any_diag[train_idx], any_diag[test_idx]
    report = ["# Held-out validation of the supervised-UMAP disease signal "
              "(`15_manifold_holdout_validation.py`)\n",
              "Tackles hypothesis #1 from `research/ANCESTRY_VS_DISEASE_MANIFOLD.md`: is "
              "`14_manifold_structure.py`'s supervised-UMAP result (obs=0.97 vs null=0.83 same-"
              "label KNN fraction) real generalizing signal, or an artifact of `umap.UMAP(y=...)` "
              "being fit and scored on the same people? See this script's module docstring for "
              "the full design.\n",
              f"\n## Split\n{len(train_idx)} train / {len(test_idx)} test "
              f"({100 * y_train.mean():.1f}% / {100 * y_test.mean():.1f}% diagnosed -- stratified "
              f"split, prevalence preserved in both halves).\n"]

    X = mat_full.values.astype(float)
    X_train, X_test = X[train_idx], X[test_idx]
    mu, sd = X_train.mean(axis=0), X_train.std(axis=0)
    sd_safe = np.where(sd == 0, 1.0, sd)
    X_train_std = (X_train - mu) / sd_safe
    X_test_std = (X_test - mu) / sd_safe  # train-only stats applied to test -- no leakage

    results = {}
    panels = {}
    for name, y_fit in [("supervised", y_train), ("unsupervised", None)]:
        print(f"Fitting {name} UMAP on train split only ...", file=sys.stderr)
        reducer, train_coords = fit_umap_train(X_train_std, y_fit, args.seed, args.umap_neighbors,
                                                args.umap_min_dist)
        if reducer is None:
            report.append(f"\n## {name}\numap-learn not installed -- skipped.\n")
            continue
        test_coords = reducer.transform(X_test_std)
        panels[name] = (train_coords, y_train, test_coords, y_test)
        scores = knn_predict_score(train_coords, y_train, test_coords, k=args.knn_k)
        point, ci_lo, ci_hi = bootstrap_auc_ci(scores, y_test, n_boot=args.n_boot, seed=args.seed)
        obs, null_mean, null_std, p = auroc_permutation_pvalue(scores, y_test,
                                                                n_perm=args.n_perm, seed=args.seed)
        results[name] = {"auroc": point, "ci_lo": ci_lo, "ci_hi": ci_hi, "p": p}
        report.append(
            f"\n## {name} (train-fit, held-out test scored)\n"
            f"Held-out AUROC: **{point:.3f}** (bootstrap 95% CI {ci_lo:.3f}-{ci_hi:.3f}, "
            f"{args.n_boot} resamples).\n"
            f"Label-permutation null (test labels shuffled, scores fixed, {args.n_perm} perms): "
            f"null mean {null_mean:.3f} +/- {null_std:.3f}, **p={p:.4g}**.\n")

    fig_path = build_train_test_figure(panels, out_dir)
    if fig_path:
        report.append(f"\n## Where held-out people land\nSmall dots are train people (colored by "
                      f"true diagnosis, used to fit/supervise the mapping); triangles are held-out "
                      f"test people `.transform()`-projected into that frozen mapping from their "
                      f"allele dosage alone -- their color was never seen by the fit or the "
                      f"transform. Figure: `{fig_path}`\n")

    if "supervised" in results and "unsupervised" in results:
        s, u = results["supervised"], results["unsupervised"]
        verdict = (
            "Supervised beats 0.5 and beats the unsupervised baseline on held-out data -- the "
            "original supervised-UMAP finding generalizes; this is real disease-relevant "
            "structure, not a training-time artifact."
            if s["auroc"] > 0.5 and s["ci_lo"] > 0.5 and s["auroc"] > u["auroc"] else
            "Supervised does NOT clear its own bootstrap CI above 0.5, or does not beat the "
            "unsupervised baseline, on held-out data -- the original supervised-UMAP finding "
            "does not generalize under this test; treat it as overfit to the labels it was "
            "given, not evidence of real disease structure."
        )
        report.append(f"\n## Verdict\nSupervised held-out AUROC {s['auroc']:.3f} "
                      f"(CI {s['ci_lo']:.3f}-{s['ci_hi']:.3f}) vs. unsupervised baseline "
                      f"{u['auroc']:.3f} (CI {u['ci_lo']:.3f}-{u['ci_hi']:.3f}).\n\n{verdict}\n")

    report_path = os.path.join(out_dir, "holdout_validation_report.md")
    vc.write_report(report_path, report)
    print(f"Wrote {report_path!r}.", file=sys.stderr)


if __name__ == "__main__":
    main()
