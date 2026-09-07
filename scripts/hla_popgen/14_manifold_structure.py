#!/usr/bin/env python3
"""Why does the allele-space UMAP/PCA show ancestry structure but not disease structure, even
though `13_disease_allele_association.py` finds a real, strong, statistically robust association
(B*27/ankylosing spondylitis)? See `research/ANCESTRY_VS_DISEASE_MANIFOLD.md` for the full
mechanistic writeup this script tests. Four experiments, all reusing `08_embedding_compare.py`'s
already-cached `mat_full` (no re-extraction needed):

1. Variance/loadings audit: how much of PC1-PC3's loading weight sits on disease-proxy allele
   columns vs. ancestry-differentiated columns (measured by a between/within-ancestry variance
   ratio, i.e. a one-way-ANOVA-style F-statistic, per column).
2. Ancestry-residualized re-embedding: subtract each column's per-ancestry-group mean, re-fit
   PCA + UMAP on the residual, color by disease burden -- does a second-order manifold appear once
   the dominant ancestry axis is removed?
3. Supervised UMAP: fit with disease-diagnosis/burden as the `y=` target (umap-learn's
   `target_metric`) -- the sharpest test of whether ANY combination of allele-dosage space
   separates cases from controls.
4. KNN-label-enrichment permutation test: turns "visible clustering or not" into a p-value, run on
   all three embeddings (locked full_enriched, residualized, supervised) x two label sets (ancestry
   as a positive control, disease burden as the real test). The KNN graph per embedding is computed
   ONCE; permutations only reshuffle labels and re-look-up neighbor labels, so this stays cheap even
   at n~9,000.

Disease burden score: count of `12_disease_phenotypes.py`'s HLA_LINKED labels matched per person
(0-11), from the wide per-person file -- read VM-side only, never committed (same rule as 12/13).

Usage (fixtures -- pure-math functions only, no umap-learn/scipy needed locally):
    python3 scripts/hla_popgen/tests/test_manifold_structure.py

Real run (VM, after 08 and 12 have both run at least once):
    python3 scripts/hla_popgen/14_manifold_structure.py
"""
import argparse
import importlib.util
import os
import pickle
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _viz_common as vc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def _load_module(filename, modname):
    path = os.path.join(HERE, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def clean_ancestry(v):
    """`ancestry_by_person.get(pid)` can come back `pd.NA` (not plain `None`) for people missing
    an ancestry call -- `pd.NA`'s `__bool__` raises, which breaks `np.unique`, `==` comparisons,
    and `or` fallbacks alike (same failure mode as 08_embedding_compare.py's ancestry-mask crash).
    Coerce to a plain str-or-None up front so every downstream comparison is safe."""
    return None if pd.isna(v) else str(v)


# ---------------------------------------------------------------------------
# 1. Variance / loadings audit
# ---------------------------------------------------------------------------
def ancestry_f_stat(col_values, ancestry_labels):
    """One-way-ANOVA-style between/within variance ratio for a single column, grouped by ancestry.
    Higher = more ancestry-differentiated. Groups with <2 members are dropped (can't estimate
    within-group variance from 1 point); returns 0.0 if fewer than 2 usable groups remain."""
    labels = np.asarray(ancestry_labels)
    # Filter None out BEFORE np.unique -- np.unique sorts its input, and comparing None to a str
    # raises TypeError, so filtering after the fact (only checking `g is not None` on the unique
    # result) is already too late.
    non_null = np.unique(labels[np.array([lbl is not None for lbl in labels])])
    uniq = [g for g in non_null if (labels == g).sum() >= 2]
    if len(uniq) < 2:
        return 0.0
    grand_mean = col_values.mean()
    ss_between, ss_within, df_between, df_within = 0.0, 0.0, len(uniq) - 1, 0
    for g in uniq:
        vals = col_values[labels == g]
        ss_between += len(vals) * (vals.mean() - grand_mean) ** 2
        ss_within += ((vals - vals.mean()) ** 2).sum()
        df_within += len(vals) - 1
    if df_within <= 0 or ss_within == 0:
        return 0.0
    return (ss_between / df_between) / (ss_within / df_within)


def variance_loadings_report(mat_full, ancestry_by_person, disease_columns_by_label,
                              pca_denoise_k=15):
    """Returns a per-column DataFrame (raw_var_frac, std_var_frac, pc1..pc3 abs loading,
    ancestry_f_stat, is_disease_proxy) plus the Spearman-free rank correlation (via numpy, no
    scipy) between |PC1 loading| and ancestry_f_stat across all columns."""
    X = mat_full.values.astype(float)
    ancestry = np.array([clean_ancestry(ancestry_by_person.get(pid)) for pid in mat_full.index],
                        dtype=object)

    raw_var = X.var(axis=0)
    raw_var_frac = raw_var / max(raw_var.sum(), 1e-12)

    mu, sd = X.mean(axis=0), X.std(axis=0)
    sd_safe = np.where(sd == 0, 1.0, sd)
    Xs = (X - mu) / sd_safe
    std_var = Xs.var(axis=0)
    std_var_frac = std_var / max(std_var.sum(), 1e-12)

    Xc = Xs - Xs.mean(axis=0)
    _U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    k = min(pca_denoise_k, Vt.shape[0])
    loadings = np.abs(Vt[:k])  # (k, n_features)

    f_stats = np.array([ancestry_f_stat(X[:, j], ancestry) for j in range(X.shape[1])])

    disease_cols = set()
    for cols in disease_columns_by_label.values():
        disease_cols.update(cols)

    df = pd.DataFrame({
        "column": mat_full.columns,
        "raw_var_frac": raw_var_frac,
        "std_var_frac": std_var_frac,
        "pc1_abs_loading": loadings[0] if k >= 1 else np.nan,
        "pc2_abs_loading": loadings[1] if k >= 2 else np.nan,
        "pc3_abs_loading": loadings[2] if k >= 3 else np.nan,
        "ancestry_f_stat": f_stats,
        "is_disease_proxy": [c in disease_cols for c in mat_full.columns],
    })

    corr = spearman_corr(df["pc1_abs_loading"].values, df["ancestry_f_stat"].values)
    return df, corr


def spearman_corr(a, b):
    """Rank correlation without scipy: Pearson correlation of the ranks."""
    ra = pd.Series(a).rank().values
    rb = pd.Series(b).rank().values
    if ra.std() == 0 or rb.std() == 0:
        return np.nan
    return float(np.corrcoef(ra, rb)[0, 1])


# ---------------------------------------------------------------------------
# 2. Ancestry residualization
# ---------------------------------------------------------------------------
def residualize_ancestry(mat_full, ancestry_by_person):
    """Subtracts each column's per-ancestry-group mean (removes between-ancestry variance only,
    leaves within-ancestry variance untouched). People with unknown/missing ancestry are left
    un-residualized (subtract the grand mean instead) rather than dropped, so the cohort stays the
    same size as `mat_full` for every downstream comparison."""
    X = mat_full.values.astype(float)
    ancestry = np.array([clean_ancestry(ancestry_by_person.get(pid)) for pid in mat_full.index],
                        dtype=object)
    grand_mean = X.mean(axis=0)
    resid = X - grand_mean
    for g in pd.unique(ancestry):
        if g is None:
            continue
        mask = ancestry == g
        if mask.sum() < 2:
            continue
        resid[mask] = X[mask] - X[mask].mean(axis=0)
    return pd.DataFrame(resid, index=mat_full.index, columns=mat_full.columns)


# ---------------------------------------------------------------------------
# 3. Supervised UMAP
# ---------------------------------------------------------------------------
def fit_supervised_umap(mat_full, y, seed, umap_neighbors, umap_min_dist):
    """y: array-like of int/categorical labels, same order as mat_full.index. Returns None if
    umap-learn isn't installed (same HAVE_UMAP guard as 08_embedding_compare.py)."""
    if not vc.HAVE_UMAP:
        return None
    import umap
    X = mat_full.values.astype(float)
    mu, sd = X.mean(axis=0), X.std(axis=0)
    Xs = (X - mu) / np.where(sd == 0, 1.0, sd)
    reducer = umap.UMAP(n_neighbors=umap_neighbors, min_dist=umap_min_dist, random_state=seed)
    return reducer.fit_transform(Xs, y=np.asarray(y))


# ---------------------------------------------------------------------------
# 4. KNN-label-enrichment permutation test (no scipy dependency; chunked to stay memory-safe)
# ---------------------------------------------------------------------------
def knn_indices(coords, k=15, chunk=500):
    n = len(coords)
    k = min(k, n - 1)
    idx_all = np.empty((n, k), dtype=int)
    for start in range(0, n, chunk):
        end = min(start + chunk, n)
        d = np.linalg.norm(coords[start:end, None, :] - coords[None, :, :], axis=2)
        for local_i, global_i in enumerate(range(start, end)):
            d[local_i, global_i] = np.inf
        idx_all[start:end] = np.argpartition(d, k - 1, axis=1)[:, :k]
    return idx_all


def knn_label_enrichment(idx_all, labels, n_perm=500, seed=0):
    """Mean fraction of a point's k nearest neighbors (in embedding space) sharing its label,
    against a null of label-permutation (neighbor graph fixed -- only labels reshuffled, so this
    stays cheap even for n~9,000). Returns (observed, null_mean, null_std, p_value) where p_value
    is one-sided (fraction of permutations with enrichment >= observed)."""
    labels = np.asarray(labels)

    def mean_same_frac(lab):
        neighbor_labels = lab[idx_all]
        return float((neighbor_labels == lab[:, None]).mean())

    observed = mean_same_frac(labels)
    rng = np.random.default_rng(seed)
    perm_scores = np.empty(n_perm)
    for p in range(n_perm):
        perm_scores[p] = mean_same_frac(rng.permutation(labels))
    p_value = (float((perm_scores >= observed).sum()) + 1) / (n_perm + 1)
    return observed, float(perm_scores.mean()), float(perm_scores.std()), p_value


# ---------------------------------------------------------------------------
# Disease burden
# ---------------------------------------------------------------------------
def load_burden(wide_diagnosis_path, person_ids):
    diag = pd.read_csv(wide_diagnosis_path, sep="\t", dtype={"person_id": str}) \
        .set_index("person_id")
    diag = diag.reindex(person_ids).fillna(False)
    label_cols = [c for c in diag.columns]
    burden = diag[label_cols].sum(axis=1).astype(int)
    any_diag = (burden > 0).astype(int)
    return burden.values, any_diag.values, label_cols


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_embedding(ax, coords, color_vals, title, cmap="viridis", discrete=False):
    if coords is None:
        ax.text(0.5, 0.5, "umap-learn not installed", ha="center", va="center", fontsize=9,
                color="#999999", transform=ax.transAxes)
        ax.set_title(title, fontsize=9)
        return None
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=color_vals, cmap=cmap, s=4, alpha=0.6,
                     linewidths=0)
    ax.set_title(title, fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
    return sc


def build_figure(embeddings, burden, out_dir, suffix):
    fig, axes = plt.subplots(1, len(embeddings), figsize=(5 * len(embeddings), 4.5))
    if len(embeddings) == 1:
        axes = [axes]
    sc = None
    for ax, (name, coords) in zip(axes, embeddings.items()):
        sc = plot_embedding(ax, coords, burden, f"{name}\n(colored by disease burden)")
    if sc is not None:
        fig.colorbar(sc, ax=axes, shrink=0.7, label="# HLA-linked diagnoses (0-11)")
    path = os.path.join(out_dir, f"manifold_structure_embeddings{suffix}.png")
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
                     help="Default: reports/hla_popgen/14_manifold_structure/")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--umap-neighbors", type=int, default=15)
    ap.add_argument("--umap-min-dist", type=float, default=0.1)
    ap.add_argument("--pca-denoise-k", type=int, default=15)
    ap.add_argument("--knn-k", type=int, default=15)
    ap.add_argument("--n-perm", type=int, default=500)
    ap.add_argument("--sample", action="store_true")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(os.path.dirname(here))
    out_dir = args.out_dir or os.path.join(repo_root, "reports", "hla_popgen",
                                            "14_manifold_structure")
    os.makedirs(out_dir, exist_ok=True)
    suffix = ".sample" if args.sample else ""

    print(f"Loading {args.embedding_cache!r} ...", file=sys.stderr)
    with open(os.path.expanduser(args.embedding_cache), "rb") as f:
        cache = pickle.load(f)
    mat_full = pd.DataFrame(cache["mat_full_values"], index=cache["mat_full_index"],
                             columns=cache["mat_full_columns"])
    ancestry_by_person = cache["ancestry_by_person"]
    locked_raw_coords = cache["coords_by_cell"].get(("full_enriched", "raw_ancestry"))

    print(f"Loading {args.wide_diagnosis_labels!r} ...", file=sys.stderr)
    burden, any_diag, label_cols = load_burden(os.path.expanduser(args.wide_diagnosis_labels),
                                                mat_full.index)

    emb08 = _load_module("08_embedding_compare.py", "hla_popgen_embedding_compare")

    report = ["# Ancestry vs. disease manifold structure (`14_manifold_structure.py`)\n",
              "See `research/ANCESTRY_VS_DISEASE_MANIFOLD.md` for the full hypothesis and "
              "experiment checklist this implements.\n",
              f"Cohort: {mat_full.shape[0]} people x {mat_full.shape[1]} allele-identity columns "
              f"(same `full_enriched` matrix as `08_embedding_compare.py`). "
              f"{int((burden > 0).sum())} ({100 * (burden > 0).mean():.1f}%) matched >=1 of "
              f"{len(label_cols)} HLA-linked diagnoses.\n"]

    # --- 1. Variance / loadings audit ---
    print("1/4: variance/loadings audit ...", file=sys.stderr)
    disease_groups = emb08.DISEASE_ALLELE_GROUPS
    disease_cols_by_label = {}
    for label, gene, allele_group in disease_groups:
        prefix = f"{vc.gene_display(gene)}*{allele_group.split('*', 1)[1]}" if "*" in \
            allele_group else None
        cols = [c for c in mat_full.columns if prefix and c.split(":", 1)[-1].startswith(prefix)]
        disease_cols_by_label[label] = cols
    loadings_df, corr = variance_loadings_report(mat_full, ancestry_by_person,
                                                  disease_cols_by_label, args.pca_denoise_k)
    loadings_path = os.path.join(out_dir, f"variance_loadings{suffix}.tsv")
    loadings_df.to_csv(loadings_path, sep="\t", index=False)
    disease_rows = loadings_df[loadings_df["is_disease_proxy"]]
    other_rows = loadings_df[~loadings_df["is_disease_proxy"]]
    report.append(
        f"\n## 1. Variance / loadings audit\n"
        f"Rank correlation (|PC1 loading| vs. ancestry F-stat, all "
        f"{len(loadings_df)} columns): **{corr:.3f}**.\n\n"
        f"Disease-proxy columns ({len(disease_rows)}): mean |PC1 loading| "
        f"{disease_rows['pc1_abs_loading'].mean():.4f}, mean ancestry F-stat "
        f"{disease_rows['ancestry_f_stat'].mean():.2f}.\n"
        f"All other columns ({len(other_rows)}): mean |PC1 loading| "
        f"{other_rows['pc1_abs_loading'].mean():.4f}, mean ancestry F-stat "
        f"{other_rows['ancestry_f_stat'].mean():.2f}.\n\n"
        f"Full per-column table: `{loadings_path}`\n")

    # --- 2. Ancestry residualization ---
    print("2/4: residualizing ancestry + re-fitting PCA/UMAP ...", file=sys.stderr)
    mat_resid = residualize_ancestry(mat_full, ancestry_by_person)
    X_resid_std = (mat_resid.values - mat_resid.values.mean(axis=0)) / \
        np.where(mat_resid.values.std(axis=0) == 0, 1.0, mat_resid.values.std(axis=0))
    pca_resid, _var, _mu, _vt = emb08.fit_pca_with_loadings(X_resid_std,
                                                             n_components=args.pca_denoise_k)
    resid_raw_coords, _resid_pca_coords = emb08.fit_umap_pair(
        X_resid_std, pca_resid, args.pca_denoise_k, args.seed, args.umap_neighbors,
        args.umap_min_dist)
    report.append(
        "\n## 2. Ancestry-residualized re-embedding\n"
        "Each column's per-ancestry-group mean subtracted before re-fitting PCA + UMAP from "
        "scratch (removes between-ancestry variance, keeps within-ancestry variance). Compare "
        "`resid_raw` panel below against `full_enriched`'s `raw_ancestry` panel in "
        "`08_embedding_compare`'s matrix figure.\n")

    # --- 3. Supervised UMAP ---
    print("3/4: supervised UMAP (target=disease burden) ...", file=sys.stderr)
    supervised_coords = fit_supervised_umap(mat_full, any_diag, args.seed, args.umap_neighbors,
                                             args.umap_min_dist)
    report.append(
        "\n## 3. Supervised UMAP (target = any HLA-linked diagnosis)\n"
        "`umap-learn` fit with `y=any_diag` -- the sharpest test of whether any combination of "
        "the allele-dosage space separates diagnosed from undiagnosed people.\n")

    fig_path = build_figure(
        {"locked_raw (unsupervised)": locked_raw_coords,
         "ancestry_residualized": resid_raw_coords,
         "supervised (y=any_diag)": supervised_coords},
        burden, out_dir, suffix)
    report.append(f"\nFigure: `{fig_path}`\n")

    # --- 4. KNN-label-enrichment permutation test ---
    print("4/4: KNN-label-enrichment permutation tests ...", file=sys.stderr)
    ancestry_labels = np.array([clean_ancestry(ancestry_by_person.get(pid)) or "unknown" for pid
                                 in mat_full.index])
    embeddings_for_test = {"locked_raw": locked_raw_coords, "ancestry_residualized":
                            resid_raw_coords, "supervised": supervised_coords}
    label_sets = {"ancestry (positive control)": ancestry_labels,
                  "any_HLA_linked_diagnosis": any_diag}
    rows = []
    for emb_name, coords in embeddings_for_test.items():
        if coords is None:
            continue
        idx_all = knn_indices(coords, k=args.knn_k)
        for label_name, labels in label_sets.items():
            obs, null_mean, null_std, p = knn_label_enrichment(idx_all, labels,
                                                                n_perm=args.n_perm, seed=args.seed)
            rows.append({"embedding": emb_name, "label_set": label_name,
                         "observed_same_label_frac": round(obs, 4),
                         "null_mean": round(null_mean, 4), "null_std": round(null_std, 4),
                         "p_value": p})
    enrich_df = pd.DataFrame(rows)
    enrich_path = os.path.join(out_dir, f"knn_label_enrichment{suffix}.tsv")
    enrich_df.to_csv(enrich_path, sep="\t", index=False)
    report.append(
        f"\n## 4. KNN-label-enrichment permutation test (k={args.knn_k}, "
        f"{args.n_perm} permutations)\n"
        f"Mean fraction of a point's k nearest embedding-space neighbors sharing its label, vs. "
        f"a label-permutation null (neighbor graph fixed). p-value is one-sided "
        f"(observed >= null).\n")
    report.append("| embedding | label set | observed | null mean +/- sd | p |")
    report.append("|---|---|---|---|---|")
    for r in rows:
        report.append(f"| {r['embedding']} | {r['label_set']} | "
                      f"{r['observed_same_label_frac']:.4f} | "
                      f"{r['null_mean']:.4f} +/- {r['null_std']:.4f} | {r['p_value']:.4g} |")
    report.append(f"\nFull table: `{enrich_path}`\n")

    report_path = os.path.join(out_dir, f"manifold_structure_report{suffix}.md")
    vc.write_report(report_path, report)
    print(f"Wrote {report_path!r}.", file=sys.stderr)


if __name__ == "__main__":
    main()
