#!/usr/bin/env python3
"""03 -- Why TCRs are public, and how robust the person-level picture is
(report: reports/03_publicness_and_robustness/).

  a  fraction of clonotypes shared by >= 2 people, by CDR3 length
  b  generation probability (OLGA, human TRB default model) by sharing level
  c  CDR3 length by sharing level
  d  person-vector robustness to the per-person cap (top-k clonotypes by reads) and to
     read-weighting, vs the reference representation (top-500, unweighted mean): Spearman
     correlation of all pairwise person distances, and preservation of each person's 10
     nearest neighbours
  e  how much genetic ancestry the person representation carries: cross-validated
     one-vs-rest AUROC, SCEPTR person vector vs plain TRBV-usage profile vs permuted labels

Known biology tested in a-c: public TCRs are public largely because V(D)J recombination
makes them likely -- short, near-germline junctions with high generation probability
(Venturi et al. 2006; Elhanati et al. 2018 Immunol Rev; Sethna et al. 2019 Bioinformatics,
OLGA). If sharing in this cohort is real convergent recombination rather than artifact,
Pgen must rise monotonically with the number of people sharing a clonotype.

d answers "is Fig 2 of report 01 an artifact of top-500 / unweighted mean?". e turns report
01's "ancestries overlap heavily" into a number, and asks whether any ancestry signal lives
in V-gene usage (germline TRBV variation) or beyond it -- the confounder the HLA phase has
to handle.

Inputs (VM-local): 02's trb_sharing.pkl; embed_cdr3s.py's pool + embeddings (--tag).
Outputs: figure + de-identified CSVs to --outdir.

Usage (from aleix/RNA-seq/, after 02; needs `pixi run pip install olga`):
  pixi run python3 -u scripts/03_publicness_and_robustness.py [--per-bin 5000] [--workers 16]
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.spatial.distance import pdist, squareform  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import repfig as R  # noqa: E402

SHARE_BINS = [(1, 1, "1"), (2, 4, "2–4"), (5, 19, "5–19"), (20, 99, "20–99"),
              (100, 10 ** 9, "≥100")]
CAPS = [25, 50, 100, 250, 500]
def person_vectors(pool, embs, k, weighted=False):
    """Mean of each person's first k pool rows (pool is ranked by reads within person)."""
    rank = pool.groupby("research_id", sort=False).cumcount().to_numpy()
    keep = rank < k
    codes, people = pd.factorize(pool["research_id"].to_numpy()[keep])
    w = pool["reads"].to_numpy(float)[keep] if weighted else np.ones(keep.sum())
    vec = np.zeros((len(people), embs.shape[1]))
    np.add.at(vec, codes, embs[keep] * w[:, None])
    vec /= np.bincount(codes, weights=w)[:, None]
    return pd.DataFrame(vec, index=people)


def knn_sets(D, k=10):
    D = D.copy()
    np.fill_diagonal(D, np.inf)
    return np.argsort(D, axis=1)[:, :k]


def cv_auroc(X, y, classes, seed, folds=5):
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    proba = np.zeros((len(y), len(classes)))
    for tr, te in StratifiedKFold(folds, shuffle=True, random_state=seed).split(X, y):
        clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0))
        clf.fit(X[tr], y[tr])
        proba[te] = clf.predict_proba(X[te])[:, [list(clf.classes_).index(c) for c in classes]]
    return {c: roc_auc_score(y == c, proba[:, i]) for i, c in enumerate(classes)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache-dir", default="~/pipeline_outputs/rnaseq/atlas")
    ap.add_argument("--embeddings-dir", default="~/pipeline_outputs/rnaseq/embeddings")
    ap.add_argument("--tag", default="cohort_full_vcdr3")
    ap.add_argument("--outdir",
                    default="~/pipeline_outputs/rnaseq/reports/03_publicness_and_robustness")
    ap.add_argument("--per-bin", type=int, default=5000,
                    help="unique clonotypes sampled per sharing bin for Pgen")
    ap.add_argument("--n-people-dist", type=int, default=1500,
                    help="people subsampled for the pairwise-distance robustness check")
    ap.add_argument("--permutations", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    args = ap.parse_args()
    outdir = os.path.expanduser(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    share_path = os.path.join(os.path.expanduser(args.cache_dir), "trb_sharing.pkl")
    if not os.path.exists(share_path):
        R.die(f"missing {share_path} -- run 02_repertoire_atlas.py first")
    share = pd.read_pickle(share_path)
    share["len"] = share["cdr3aa"].str.len()
    share["bin"] = pd.cut(share["n_people"], [lo - 0.5 for lo, _, _ in SHARE_BINS] + [1e12],
                          labels=[lab for *_, lab in SHARE_BINS])
    labels = [lab for *_, lab in SHARE_BINS]

    # ---------------- a: public fraction by length ----------------
    by_len = (share.groupby("len")
                   .agg(n=("n_people", "size"), public=("n_people", lambda s: np.mean(s >= 2)))
                   .loc[8:22])
    by_len.to_csv(os.path.join(outdir, "public_fraction_by_length.csv"))

    # ---------------- b/c: Pgen by sharing level ----------------
    samp = (share.sample(frac=1, random_state=args.seed)
                 .groupby("bin", observed=True).head(args.per_bin).copy())
    t0 = time.time()
    samp["pgen"] = R.olga_pgen(samp["cdr3aa"], args.workers)
    print(f"OLGA Pgen for {len(samp):,} clonotypes: {time.time() - t0:.0f}s", file=sys.stderr)
    samp["log10_pgen"] = np.log10(samp["pgen"].where(samp["pgen"] > 0))
    zero_frac = float((samp["pgen"] <= 0).mean())
    ok = samp.dropna(subset=["log10_pgen"])
    rho, rho_p = spearmanr(ok["n_people"], ok["log10_pgen"])
    pg_sum = R.quantile_summary(ok, "bin", "log10_pgen", min_n=1, order=labels)
    len_sum = R.quantile_summary(samp, "bin", "len", min_n=1, order=labels)
    pd.concat({"log10_pgen": pg_sum, "cdr3_length": len_sum}).to_csv(
        os.path.join(outdir, "pgen_and_length_by_sharing.csv"))

    # ---------------- d/e: robustness of the person representation ----------------
    edir = os.path.expanduser(args.embeddings_dir)
    embs = np.load(os.path.join(edir, f"embeddings_sceptr_{args.tag}.npy"))
    pool = pd.read_csv(os.path.join(edir, f"pool_sceptr_{args.tag}.tsv"), sep="\t",
                       dtype={"research_id": str}, keep_default_na=False)
    pool["reads"] = pd.to_numeric(pool["reads"], errors="coerce").fillna(1)
    ref = person_vectors(pool, embs, 500)
    sub = rng.choice(ref.index.to_numpy(), size=min(args.n_people_dist, len(ref)), replace=False)
    d_ref = squareform(pdist(ref.loc[sub].to_numpy()))
    nn_ref = knn_sets(d_ref)
    rob = []
    for label, k, w in [(f"top-{k}", k, False) for k in CAPS] + [("top-500, read-weighted", 500, True)]:
        v = person_vectors(pool, embs, k, weighted=w).reindex(sub)
        d = squareform(pdist(v.to_numpy()))
        iu = np.triu_indices(len(sub), 1)
        r = spearmanr(d[iu], d_ref[iu])[0]
        nn = knn_sets(d)
        overlap = np.mean([len(set(a) & set(b)) / 10 for a, b in zip(nn, nn_ref)])
        rob.append((label, k, w, r, overlap))
    rob = pd.DataFrame(rob, columns=["representation", "k", "read_weighted",
                                     "spearman_vs_ref", "knn10_overlap_vs_ref"])
    rob.to_csv(os.path.join(outdir, "person_vector_robustness.csv"), index=False)

    # ---------------- f: ancestry signal ----------------
    anc_of = pool.drop_duplicates("research_id").set_index("research_id")["ancestry"]
    y = anc_of.reindex(ref.index).to_numpy()
    classes = [a for a in R.ANCESTRY_ORDER if np.sum(y == a) >= R.MIN_PEOPLE]
    keep = np.isin(y, classes)
    pool["trbv"] = pool["v_gene"].astype(str).str.split("*").str[0]
    vuse = (pool.groupby(["research_id", "trbv"]).size().unstack(fill_value=0))
    vuse = vuse.div(vuse.sum(axis=1), axis=0).reindex(ref.index).fillna(0)
    feats = {"SCEPTR person vector": ref.to_numpy(), "TRBV usage": vuse.to_numpy()}
    auc = {name: cv_auroc(X[keep], y[keep], classes, args.seed) for name, X in feats.items()}
    null = []
    for i in range(args.permutations):
        yp = rng.permutation(y[keep])
        null.append(cv_auroc(feats["SCEPTR person vector"][keep], yp, classes, args.seed + i))
    null = pd.DataFrame(null)
    auc_df = pd.DataFrame(auc)
    auc_df["null_mean"] = null.mean()
    auc_df["null_sd"] = null.std()
    auc_df.index.name = "ancestry"
    auc_df.to_csv(os.path.join(outdir, "ancestry_auroc.csv"))

    # ---------------- figure ----------------
    import cnsplots as cns
    mp = cns.multipanel(max_width=540)
    pal = R.nature()
    W, H, GAP = 100, 100, 30
    share_cols = [pal[i] for i in (6, 1, 3, 0, 7)]

    ax = mp.panel("a", width=W, height=H, margin_right=15, margin_bottom=GAP)
    ax.plot(by_len.index, 100 * by_len["public"], "-o", color=pal[3], ms=2.5, lw=0.9)
    R.style(ax, "Publicness by length", "CDR3 length (aa)", "Shared by ≥2 people (%)")

    ax = mp.panel("b", width=W, height=H, margin_right=15, margin_bottom=GAP)
    R.pointrange(ax, pg_sum, colors=share_cols[:len(pg_sum)])
    ax.set_title("Generation probability", pad=11)
    ax.text(0.5, 1.02, f"Spearman ρ = {rho:.2f}", transform=ax.transAxes, ha="center",
            va="bottom", fontsize=5.5)
    R.style(ax, None, "People sharing", "log$_{10}$ P$_{gen}$")

    ax = mp.panel("c", width=W, height=H, margin_right=15, margin_bottom=GAP)
    R.pointrange(ax, len_sum, colors=share_cols[:len(len_sum)])
    R.style(ax, "CDR3 length", "People sharing", "Length (aa)")

    ax = mp.panel("d", width=W, height=H, margin_right=95)
    unw = rob[~rob["read_weighted"]]
    wt = rob[rob["read_weighted"]]
    from matplotlib.lines import Line2D
    handles = []
    for col, lab, colr in (("spearman_vs_ref", "Pairwise distances (ρ)", pal[3]),
                           ("knn10_overlap_vs_ref", "10 nearest neighbours", pal[0])):
        ax.plot(unw["k"], unw[col], "-o", color=colr, ms=2.5, lw=0.9)
        ax.scatter(wt["k"], wt[col], marker="D", s=14, facecolors="white", edgecolors=colr,
                   linewidths=0.9, zorder=3)
        handles.append(Line2D([], [], color=colr, marker="o", ms=2.5, lw=0.9, label=lab))
    handles.append(Line2D([], [], color="0.3", marker="D", ms=3, mfc="white", lw=0,
                          label="Read-weighted, k=500"))
    ax.set_xscale("log")
    ax.set_xticks(CAPS)
    ax.set_xticklabels([str(k) for k in CAPS])
    ax.minorticks_off()
    ax.set_ylim(0, 1.02)
    ax.legend(handles=handles, frameon=False)
    cns.take_legend_out(title="Agreement with k=500", ax=ax)
    R.style(ax, "Robustness of person vector", "Clonotypes per person (k)", "Agreement")

    ax = mp.panel("e", width=150, height=H, margin_right=70)
    x = np.arange(len(classes))
    bw = 0.36
    for j, (name, colr) in enumerate((("SCEPTR person vector", pal[3]),
                                      ("TRBV usage", pal[1]))):
        ax.bar(x + (j - 0.5) * bw, auc_df.loc[classes, name], width=bw, color=colr,
               linewidth=0, label=name)
    ax.fill_between([-0.6, len(classes) - 0.4],
                    [auc_df["null_mean"].mean() - 2 * auc_df["null_sd"].max()] * 2,
                    [auc_df["null_mean"].mean() + 2 * auc_df["null_sd"].max()] * 2,
                    color="0.85", linewidth=0, zorder=0, label="Permuted labels (±2 s.d.)")
    ax.axhline(0.5, color="0.5", lw=0.5, ls="--", zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels(classes)
    ax.set_xlim(-0.6, len(classes) - 0.4)
    ax.set_ylim(0.4, 1.0)
    ax.legend(frameon=False)
    cns.take_legend_out(title="Features", ax=ax)
    R.style(ax, "Ancestry predictable from repertoire?", "Genetic ancestry",
            "One-vs-rest AUROC (5-fold CV)")
    R.save("fig_publicness_and_robustness", outdir)

    rows = [("clonotypes_unique", len(share)), ("pgen_sampled", len(samp)),
            ("pgen_zero_frac", zero_frac), ("spearman_sharing_vs_log10_pgen", rho),
            ("spearman_p", rho_p), ("people_distance_subsample", len(sub))]
    rows += [(f"robust_spearman_{r.representation}", r.spearman_vs_ref) for r in rob.itertuples()]
    rows += [(f"robust_knn10_{r.representation}", r.knn10_overlap_vs_ref) for r in rob.itertuples()]
    rows += [(f"auroc_macro_{n}", float(np.mean(list(a.values())))) for n, a in auc.items()]
    rows += [("auroc_macro_null_mean", float(null.mean(axis=1).mean()))]
    summary = pd.DataFrame(rows, columns=["metric", "value"])
    summary.to_csv(os.path.join(outdir, "summary.csv"), index=False)
    print("\n" + summary.to_string(index=False))
    print(f"\nAll outputs in {outdir}")


if __name__ == "__main__":
    main()
