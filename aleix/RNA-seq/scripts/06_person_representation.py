#!/usr/bin/env python3
"""06 -- What is the best way to turn a person's TRB repertoire into a vector?
(report: reports/06_person_representation/)

Report 03 found the representation we had been using -- the mean of a person's SCEPTR
clonotype vectors -- is not stable: recomputing it from that person's top 250 clonotypes
instead of their top 500 gives pairwise distances correlating only 0.35 with the original,
and read-weighting moves it as much again. A representation that changes that much under an
arbitrary choice cannot support anything downstream. This benchmarks the alternatives.

Judged on four axes, none of which needs a disease label:

  1  IDENTIFIABILITY (the primary axis). Split each person's clonotypes at random into two
     halves, build the representation independently on each, and ask whether a person's own
     second half is their nearest neighbour among all people. This is the fingerprinting
     test used for the same purpose in connectomics and microbiome work (Finn et al. 2015
     Nat Neurosci). It measures whether a representation captures stable individual identity
     or sampling noise, and it is exactly what the top-250-vs-500 instability violated.
  2  BIOLOGICAL SIGNAL. Out-of-fold prediction of age (ridge, R^2 and MAE) and of sex
     (logistic, AUROC). Age is the right probe because 02 established a strong, real
     age effect in this cohort (diversity -3.3%/decade, expansion +1.9 pp/decade).
  3  NUISANCE. Out-of-fold prediction of log sequencing depth. A good representation should
     carry LESS of this. Reported against signal, because the two trade off.
  4  ANCESTRY. Macro one-vs-rest AUROC -- reported as a property to be aware of, not a
     target: 03 showed plain TRBV usage already reaches it, so it is largely germline
     V-gene composition.

Representations compared, built from 02's full clonotype cache: TRBV usage, CDR3-length profile, CDR3 3-mer composition, mean and
read-weighted mean SCEPTR vectors, SCEPTR metacluster abundance profiles at several
cluster counts, their centred-log-ratio transforms, and concatenations. A label-permuted
control is run through the identical pipeline.

Which clonotypes: every full-cache clonotype that has a SCEPTR vector. With
embed_full_cache.py's output present (the default since 2026-10-04) that is every clonotype
with a usable TRBV; without it, vectors are looked up from the top-500 embedding pool, which
covers only clonotypes whose sequence is in somebody's top 500 (biased toward expanded and
public sequences). The coverage actually achieved is written to summary.csv.

Identifiability is reported three ways (2026-10-04): on the raw vectors (cosine), on
z-scored features (each column standardised with full-data statistics, so cosine becomes a
correlation and a shared mean direction cannot dominate), and on z-scored features with the
linear effect of that half's log sequencing depth removed (how much identity survives once
depth is taken out). Concatenated representations are given both as a raw np.hstack and
block-balanced (each block standardised and scaled by 1/sqrt(width), repfig.balanced_concat),
because a raw hstack lets the block with the larger magnitude dominate the cosine.

Compositional profiles are also given as CLR (centred log-ratio) because abundance profiles
are compositional: raw proportions are constrained to sum to one, which induces spurious
negative correlations between components (Aitchison 1982).

Usage (from aleix/RNA-seq/, after 02):
  pixi run python3 -u scripts/06_person_representation.py [--n-people 4000] [--seed 0]
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402
from sklearn.cluster import MiniBatchKMeans  # noqa: E402
from sklearn.decomposition import TruncatedSVD  # noqa: E402
from sklearn.linear_model import LogisticRegression, Ridge  # noqa: E402
from sklearn.metrics import mean_absolute_error, r2_score, roc_auc_score  # noqa: E402
from sklearn.model_selection import KFold, StratifiedKFold  # noqa: E402
from sklearn.pipeline import make_pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

import repfig as R  # noqa: E402

MIN_CLONOTYPES = 100      # need enough to split in half
KS = (30, 100, 300)       # metacluster counts to test
SVD_DIM = 100
COMBOS = {                # concatenations: raw hstack here, block-balanced added after build
    "SCEPTR mean + TRBV": ("SCEPTR mean", "TRBV usage"),
    "TRBV + clusters k=100 (CLR)": ("_TRBV usage, CLR", "SCEPTR clusters k=100, CLR"),
}


def clr(P, eps=1e-6):
    """Centred log-ratio: proportions -> unconstrained coordinates (Aitchison 1982)."""
    L = np.log(P + eps)
    return L - L.mean(axis=1, keepdims=True)


def profile(codes, n_cat, person_idx, n_people, weights=None):
    """Sparse person x category count matrix, row-normalised to proportions."""
    w = np.ones(len(codes)) if weights is None else weights
    M = csr_matrix((w, (person_idx, codes)), shape=(n_people, n_cat)).toarray()
    return M / np.maximum(M.sum(axis=1, keepdims=True), 1e-9)


def mean_pool(embs, person_idx, n_people, weights=None):
    w = np.ones(len(person_idx)) if weights is None else weights
    out = np.zeros((n_people, embs.shape[1]))
    np.add.at(out, person_idx, embs * w[:, None])
    denom = np.bincount(person_idx, weights=w, minlength=n_people)
    return out / np.maximum(denom, 1e-9)[:, None]


def identifiability(A, B):
    """Split-half fingerprinting. A[i], B[i] are two independent halves of person i.
    Returns (top-1 accuracy, mean percentile rank of the correct match)."""
    A = A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-9)
    B = B / np.maximum(np.linalg.norm(B, axis=1, keepdims=True), 1e-9)
    S = A @ B.T
    self_s = np.diag(S).copy()
    better = (S > self_s[:, None]).sum(axis=1)
    n = len(A)
    return float(np.mean(better == 0)), float(np.mean(1 - better / (n - 1)))


def cv_regress(X, y, seed, folds=5):
    pred = np.zeros(len(y))
    for tr, te in KFold(folds, shuffle=True, random_state=seed).split(X):
        m = make_pipeline(StandardScaler(), Ridge(alpha=10.0)).fit(X[tr], y[tr])
        pred[te] = m.predict(X[te])
    return r2_score(y, pred), mean_absolute_error(y, pred)


def cv_auroc_binary(X, y, seed, folds=5):
    pred = np.zeros(len(y))
    for tr, te in StratifiedKFold(folds, shuffle=True, random_state=seed).split(X, y):
        m = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=3000)).fit(X[tr], y[tr])
        pred[te] = m.predict_proba(X[te])[:, 1]
    return roc_auc_score(y, pred)


def cv_auroc_macro(X, y, classes, seed, folds=5):
    proba = np.zeros((len(y), len(classes)))
    for tr, te in StratifiedKFold(folds, shuffle=True, random_state=seed).split(X, y):
        m = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=3000)).fit(X[tr], y[tr])
        cols = [list(m.classes_).index(c) for c in classes]
        proba[te] = m.predict_proba(X[te])[:, cols]
    return float(np.mean([roc_auc_score(y == c, proba[:, i])
                          for i, c in enumerate(classes)]))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", default="~/pipeline_outputs/rnaseq/cohort_full.tsv")
    ap.add_argument("--pheno-dir", default="~/pipeline_outputs/rnaseq/pheno")
    ap.add_argument("--cache-dir", default="~/pipeline_outputs/rnaseq/atlas")
    ap.add_argument("--embeddings-dir", default="~/pipeline_outputs/rnaseq/embeddings")
    ap.add_argument("--tag", default="cohort_full_vcdr3")
    ap.add_argument("--n-people", type=int, default=4000,
                    help="people subsampled for the benchmark (0 = all)")
    ap.add_argument("--outdir", default="~/pipeline_outputs/rnaseq/reports/06_person_representation")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    outdir = os.path.expanduser(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # ---------------- clonotypes, with a SCEPTR vector for each ----------------
    clono = pd.read_pickle(os.path.join(os.path.expanduser(args.cache_dir),
                                        "trb_clonotypes_full.pkl"))
    clono["research_id"] = clono["research_id"].astype(str)
    clono["trbv"] = clono["trbv"].astype(str)
    edir = os.path.expanduser(args.embeddings_dir)
    full_npy = os.path.join(edir, "embeddings_sceptr_fullcache.npy")
    full_keys = os.path.join(edir, "keys_sceptr_fullcache.tsv")
    if os.path.exists(full_npy) and os.path.exists(full_keys):
        embs_pool = np.load(full_npy, mmap_mode="r")
        pool = pd.read_csv(full_keys, sep="\t", keep_default_na=False)
        vector_source = "full cache (embed_full_cache.py)"
    else:
        embs_pool = np.load(os.path.join(edir, f"embeddings_sceptr_{args.tag}.npy"))
        pool = pd.read_csv(os.path.join(edir, f"pool_sceptr_{args.tag}.tsv"), sep="\t",
                           dtype={"research_id": str}, keep_default_na=False)
        pool["trbv"] = pool["v_gene"].astype(str).str.split("*").str[0]
        vector_source = "top-500 embedding pool (lookup)"

    # One vector per distinct (TRBV, CDR3): SCEPTR is deterministic, so the embedding of a
    # sequence is the same whoever carries it.
    key = pool["trbv"] + "|" + pool["cdr3aa"]
    first = ~key.duplicated()
    lut = pd.Series(np.flatnonzero(first), index=key[first])
    ck = clono["trbv"] + "|" + clono["cdr3aa"]
    row = lut.reindex(ck).to_numpy()
    have = ~pd.isna(row)
    coverage = float(have.mean())
    print(f"{len(clono):,} clonotypes; {have.sum():,} ({100 * coverage:.1f}%) have a "
          f"SCEPTR vector from the {vector_source}", file=sys.stderr)
    clono = clono[have].copy()
    emb_idx = row[have].astype(int)

    # ---------------- people ----------------
    people = R.load_people(args.cohort, args.pheno_dir)
    n_cl = clono.groupby("research_id").size()
    ok = n_cl[n_cl >= MIN_CLONOTYPES].index
    people = people[people["research_id"].isin(ok)].dropna(subset=["age", "sex"])
    people = people[people["age"] >= 18]
    if args.n_people and len(people) > args.n_people:
        people = people.sample(args.n_people, random_state=args.seed)
    people = people.sort_values("research_id").reset_index(drop=True)
    pid = pd.Series(np.arange(len(people)), index=people["research_id"])
    keep = clono["research_id"].isin(pid.index)
    clono, emb_idx = clono[keep], emb_idx[keep.to_numpy()]
    person_idx = pid.reindex(clono["research_id"]).to_numpy()
    n_people = len(people)
    print(f"Benchmark cohort: {n_people:,} people, {len(clono):,} clonotypes",
          file=sys.stderr)

    embs = np.asarray(embs_pool[emb_idx], dtype=np.float32)
    reads = clono["reads"].to_numpy(float)
    half = rng.random(len(clono)) < 0.5      # split-half assignment, fixed across methods

    # ---------------- feature spaces ----------------
    vcode, vcats = pd.factorize(clono["trbv"])
    lcode, lcats = pd.factorize(clono["cdr3aa"].str.len())
    t0 = time.time()
    clusters = {}
    for k in KS:
        clusters[k] = MiniBatchKMeans(n_clusters=k, random_state=args.seed,
                                      batch_size=10_000, n_init=3).fit_predict(embs)
    print(f"k-means at k={KS}: {time.time() - t0:.0f}s", file=sys.stderr)

    uniq_seq, inv = np.unique(clono["cdr3aa"].to_numpy(), return_inverse=True)
    vocab, rows, cols = {}, [], []
    for i, s in enumerate(uniq_seq):
        for j in range(len(s) - 2):
            cols.append(vocab.setdefault(s[j:j + 3], len(vocab)))
            rows.append(i)
    Kmer = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(uniq_seq), len(vocab)))
    Kmer = Kmer.multiply(1 / np.maximum(np.sqrt(Kmer.multiply(Kmer).sum(1)).A1, 1e-9)[:, None])
    kmer_svd = TruncatedSVD(SVD_DIM, random_state=args.seed).fit_transform(Kmer.tocsr())
    kmer_per_clono = kmer_svd[inv]
    print(f"3-mer vocabulary {len(vocab)}, SVD -> {SVD_DIM}", file=sys.stderr)

    def build(mask):
        """All representations, computed on the clonotypes selected by `mask`."""
        m = np.flatnonzero(mask)
        pi, w = person_idx[m], reads[m]
        out = {
            "TRBV usage": profile(vcode[m], len(vcats), pi, n_people),
            "CDR3 length profile": profile(lcode[m], len(lcats), pi, n_people),
            "CDR3 3-mer (SVD)": mean_pool(kmer_per_clono[m], pi, n_people),
            "SCEPTR mean": mean_pool(embs[m], pi, n_people),
            "SCEPTR mean, read-wt": mean_pool(embs[m], pi, n_people, w),
        }
        for k in KS:
            P = profile(clusters[k][m], k, pi, n_people)
            out[f"SCEPTR clusters k={k}"] = P
            out[f"SCEPTR clusters k={k}, CLR"] = clr(P)
            if k == KS[1]:
                out[f"SCEPTR clusters k={k}, read-wt"] = profile(clusters[k][m], k, pi,
                                                                 n_people, w)
        out["_TRBV usage, CLR"] = clr(profile(vcode[m], len(vcats), pi, n_people))
        out["_log depth"] = np.log(np.maximum(
            np.bincount(pi, weights=w, minlength=n_people), 1.0))
        for name, parts in COMBOS.items():
            out[name] = np.hstack([out[p_] for p_ in parts])
        return out

    t0 = time.time()
    full = build(np.ones(len(clono), bool))
    A, B = build(half), build(~half)
    # block-balanced versions of the concatenations (statistics from the full build)
    for name, parts in COMBOS.items():
        for D_ in (full, A, B):
            D_[f"{name} (balanced)"] = R.balanced_concat([D_[p_] for p_ in parts],
                                                         [full[p_] for p_ in parts])
    depth_half = {"A": A.pop("_log depth"), "B": B.pop("_log depth")}
    full.pop("_log depth")
    for D_ in (full, A, B):
        D_.pop("_TRBV usage, CLR")
    print(f"Built {len(full)} representations x3: {time.time() - t0:.0f}s", file=sys.stderr)

    # ---------------- evaluate ----------------
    age = people["age"].to_numpy(float)
    sex = (people["sex"] == "Female").to_numpy().astype(int)
    anc = people["ancestry"].to_numpy()
    classes = [a for a in R.ANCESTRY_ORDER if np.sum(anc == a) >= R.MIN_PEOPLE]
    anc_ok = np.isin(anc, classes)
    depth = np.log(clono.groupby("research_id")["reads"].sum()
                        .reindex(people["research_id"]).to_numpy())

    rows_out = []
    for name in full:
        t1, pr = identifiability(A[name], B[name])
        Az, Bz = R.zscore_like(A[name], full[name]), R.zscore_like(B[name], full[name])
        t1z, prz = identifiability(Az, Bz)
        t1zd, _ = identifiability(R.residualize(Az, depth_half["A"]),
                                  R.residualize(Bz, depth_half["B"]))
        r2_age, mae_age = cv_regress(full[name], age, args.seed)
        r2_depth, _ = cv_regress(full[name], depth, args.seed)
        auc_sex = cv_auroc_binary(full[name], sex, args.seed)
        auc_anc = cv_auroc_macro(full[name][anc_ok], anc[anc_ok], classes, args.seed)
        rows_out.append((name, full[name].shape[1], t1, pr, t1z, prz, t1zd, r2_age, mae_age,
                         auc_sex, auc_anc, r2_depth))
        print(f"  {name:34s} dim {full[name].shape[1]:4d}  ID {100 * t1:5.1f}%  "
              f"age R2 {r2_age:+.3f}  sex {auc_sex:.3f}  depth R2 {r2_depth:+.3f}",
              file=sys.stderr)

    # permuted-label control on the best representation by identifiability
    res = pd.DataFrame(rows_out, columns=["representation", "dim", "identifiability_top1",
                                          "identifiability_pct_rank", "identifiability_top1_z",
                                          "identifiability_pct_rank_z",
                                          "identifiability_top1_z_depthres", "age_r2",
                                          "age_mae", "sex_auroc", "ancestry_auroc",
                                          "depth_r2"])
    for col in ("identifiability_top1", "identifiability_top1_z",
                "identifiability_top1_z_depthres"):
        _, lo, hi = R.wilson(np.round(res[col].to_numpy() * n_people), n_people)
        res[f"{col}_lo"], res[f"{col}_hi"] = lo, hi
    best = res.sort_values("identifiability_top1", ascending=False).iloc[0]["representation"]
    best_z = res.sort_values("identifiability_top1_z", ascending=False).iloc[0]["representation"]
    perm = rng.permutation(n_people)
    ctrl = {
        "age_r2": cv_regress(full[best], age[perm], args.seed)[0],
        "sex_auroc": cv_auroc_binary(full[best], sex[perm], args.seed),
        "identifiability_top1": identifiability(A[best], B[best][perm])[0],
    }
    res.to_csv(os.path.join(outdir, "representation_benchmark.csv"), index=False)
    pd.Series(ctrl).to_csv(os.path.join(outdir, "permuted_control.csv"))

    # ---------------- figure ----------------
    import cnsplots as cns
    mp = cns.multipanel(max_width=540)
    pal = R.nature()
    W, H, GAP = 105, 105, 34

    def fam(n):
        if " + " in n:
            return "combined"
        if n.startswith("SCEPTR clusters"):
            return "SCEPTR clusters"
        if n.startswith("SCEPTR mean"):
            return "SCEPTR mean"
        return "simple baseline"
    fcol = {"SCEPTR clusters": pal[0], "SCEPTR mean": pal[3], "combined": pal[2],
            "simple baseline": "0.65"}
    res["family"] = res["representation"].map(fam)
    short = {n: n.replace("SCEPTR ", "").replace(" profile", "").replace("(SVD)", "")
             for n in res["representation"]}

    r = res.sort_values("identifiability_top1")
    ax = mp.panel("a", width=180, height=H + 30, margin_right=12, margin_bottom=GAP)
    y = np.arange(len(r))
    ax.barh(y, 100 * r["identifiability_top1"], height=0.68,
            color=[fcol[f] for f in r["family"]], linewidth=0, label="raw cosine")
    ax.scatter(100 * r["identifiability_top1_z"], y, s=9, color="0.1", zorder=3,
               linewidths=0, label="z-scored")
    ax.scatter(100 * r["identifiability_top1_z_depthres"], y, s=9, facecolors="none",
               edgecolors="0.1", linewidths=0.6, zorder=3, label="z-scored, depth removed")
    ax.legend(frameon=False, fontsize=5, loc="lower right")
    ax.axvline(100 / n_people, color="0.5", lw=0.5, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([short[n] for n in r["representation"]], fontsize=5)
    ax.set_title("Split-half identifiability", pad=11)
    ax.text(0.5, 1.02, f"own second half is nearest of {n_people:,}",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=5.5)
    R.style(ax, None, "Correct (%)", None)

    ax = mp.panel("b", width=W, height=H + 30, margin_right=12, margin_bottom=GAP)
    r2 = res.set_index("representation").loc[r["representation"]].reset_index()  # a's order
    y = np.arange(len(r2))
    ax.barh(y, r2["age_r2"], height=0.68, color=[fcol[f] for f in r2["family"]],
            linewidth=0)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.axvline(0, color="0.5", lw=0.5)
    ax.set_title("Age prediction", pad=11)
    ax.text(0.5, 1.02, "5-fold out-of-fold", transform=ax.transAxes, ha="center",
            va="bottom", fontsize=5.5)
    R.style(ax, None, "R²", None)

    ax = mp.panel("c", width=W, height=H, margin_right=60, margin_bottom=GAP)
    for f, g_ in res.groupby("family"):
        ax.scatter(g_["depth_r2"], g_["age_r2"], s=16, linewidths=0, color=fcol[f], label=f)
    ax.legend(frameon=False)
    cns.take_legend_out(title="Family", ax=ax)
    ax.set_title("Signal vs nuisance", pad=11)
    ax.text(0.5, 1.02, "up and to the left is better", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=5.5)
    R.style(ax, None, "Sequencing depth R² (nuisance)", "Age R² (signal)")

    ax = mp.panel("d", width=W, height=H, margin_right=12)
    kk = [k for k in KS]
    raw = [res.loc[res["representation"] == f"SCEPTR clusters k={k}",
                   "identifiability_top1"].iat[0] for k in kk]
    cl_ = [res.loc[res["representation"] == f"SCEPTR clusters k={k}, CLR",
                   "identifiability_top1"].iat[0] for k in kk]
    ax.plot(kk, 100 * np.array(raw), "-o", ms=3, lw=0.9, color=pal[0], label="proportions")
    ax.plot(kk, 100 * np.array(cl_), "-o", ms=3, lw=0.9, color=pal[1], label="CLR")
    ax.set_xscale("log")
    ax.set_xticks(kk)
    ax.set_xticklabels([str(k) for k in kk])
    ax.minorticks_off()
    ax.legend(frameon=False, fontsize=5.5)
    ax.set_title("Cluster count", pad=11)
    R.style(ax, None, "Metaclusters (k)", "Identifiability (%)")

    ax = mp.panel("e", width=W, height=H, margin_right=12)
    rs = res.sort_values("sex_auroc")
    y = np.arange(len(rs))
    ax.barh(y - 0.2, rs["sex_auroc"] - 0.5, left=0.5, height=0.38, color=pal[3],
            linewidth=0, label="Sex")
    ax.barh(y + 0.2, rs["ancestry_auroc"] - 0.5, left=0.5, height=0.38, color=pal[1],
            linewidth=0, label="Ancestry (macro)")
    ax.axvline(0.5, color="0.5", lw=0.5, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([short[n] for n in rs["representation"]], fontsize=4.5)
    ax.legend(frameon=False, fontsize=5.5, loc="lower center", ncol=2,
              bbox_to_anchor=(0.5, 1.0))
    ax.set_title("Demographics", pad=16)
    R.style(ax, None, "AUROC", None)

    ax = mp.panel("f", width=W, height=H, margin_right=12)
    cols = ["identifiability_top1", "age_r2", "sex_auroc", "depth_r2"]
    labs = ["Identity", "Age", "Sex", "Depth"]
    Z = res.set_index("representation")[cols]
    Z = (Z - Z.mean()) / Z.std(ddof=0)
    Z = Z.loc[res.sort_values("identifiability_top1", ascending=False)["representation"]]
    im = ax.imshow(Z.to_numpy(), cmap="bwr", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(len(labs)))
    ax.set_xticklabels(labs, rotation=45, ha="right")
    ax.set_yticks(range(len(Z)))
    ax.set_yticklabels([short[n] for n in Z.index], fontsize=4.5)
    R.colorbar(ax, im, "z-score across representations")
    ax.set_title("All metrics", pad=11)
    R.save("fig_person_representation", outdir)

    summary = pd.concat([
        pd.DataFrame({"metric": ["people", "clonotypes", "vector_source",
                                 "frac_cache_clonotypes_with_vector",
                                 "best_by_identifiability", "best_by_identifiability_z"],
                      "value": [n_people, len(clono), vector_source, coverage, best, best_z]}),
        pd.DataFrame({"metric": [f"permuted_{k}" for k in ctrl],
                      "value": list(ctrl.values())})])
    summary.to_csv(os.path.join(outdir, "summary.csv"), index=False)
    print("\n" + res.to_string(index=False))
    print(f"\nBest by identifiability: {best}")
    print(f"Permuted-label control: {ctrl}")
    print(f"\nAll outputs in {outdir}")


if __name__ == "__main__":
    main()
