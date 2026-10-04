#!/usr/bin/env python3
"""05 -- Metacluster structure of the SCEPTR TRB space, and what defines each cluster
(report: reports/05_embedding_clusters/).

Report 01 showed the clonotype embedding breaks into V-gene islands. This asks the next
question quantitatively: if the space is partitioned into k metaclusters, what distinguishes
them, and do TCRs of KNOWN antigen specificity concentrate in particular clusters?

  a  cluster map (PC1-PC2 of the 64-d centroids), point size = clonotypes, colour = dominant
     TRBV family
  b  the properties that separate clusters: CDR3 length, generation probability, publicness,
     clonal expansion (each cluster one point)
  c  VDJdb specificity enrichment: for every (cluster, pathogen) pair, odds ratio of exact
     VDJdb matches falling in that cluster vs the rest of the space, Fisher exact,
     Benjamini-Hochberg across all tested pairs
  d  the strongest enriched cluster-pathogen pairs, as odds ratio with 95% CI

Also writes `clusters_3d.csv` -- cluster centroids in 3-d PCA space with every per-cluster
statistic -- which is what the interactive 3-d cluster map is built from. Only cluster-level
aggregates leave the VM; every reported cluster pools >= 20 distinct participants.

Clustering: MiniBatchKMeans on all embedded clonotypes (k-means on a learned metric space;
no claim that clusters are discrete biological types -- they are a partition used to make
the space describable). Stability is reported as the mean Jaccard of cluster membership
against a second seed.

Usage (from aleix/RNA-seq/, after 02 and 04's VDJdb download):
  pixi run python3 -u scripts/05_embedding_clusters.py [--k 60] [--workers 16]
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402
from sklearn.cluster import MiniBatchKMeans  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402

import repfig as R  # noqa: E402

MIN_MATCHES_PER_TEST = 20   # VDJdb-matched clonotypes needed before testing a pathogen
TOP_PAIRS = 12


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--embeddings-dir", default="~/pipeline_outputs/rnaseq/embeddings")
    ap.add_argument("--tag", default="cohort_full_vcdr3")
    ap.add_argument("--cache-dir", default="~/pipeline_outputs/rnaseq/atlas")
    ap.add_argument("--vdjdb", default=None, help="vdjdb.slim.txt or .zip (default: the copy "
                    "04 downloaded into --cache-dir)")
    ap.add_argument("--k", type=int, default=60)
    ap.add_argument("--pgen-sample", type=int, default=60_000,
                   help="clonotypes sampled for per-cluster generation probability")
    ap.add_argument("--outdir", default="~/pipeline_outputs/rnaseq/reports/05_embedding_clusters")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    args = ap.parse_args()
    outdir = os.path.expanduser(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    edir = os.path.expanduser(args.embeddings_dir)
    embs = np.load(os.path.join(edir, f"embeddings_sceptr_{args.tag}.npy"))
    pool = pd.read_csv(os.path.join(edir, f"pool_sceptr_{args.tag}.tsv"), sep="\t",
                       dtype={"research_id": str}, keep_default_na=False)
    if len(embs) != len(pool):
        R.die("embeddings and pool are not row-aligned")
    pool["trbv"] = pool["v_gene"].astype(str).str.split("*").str[0]
    pool["reads"] = pd.to_numeric(pool["reads"], errors="coerce").fillna(1)
    pool["len"] = pool["cdr3aa"].str.len()
    print(f"{len(pool):,} clonotypes x {embs.shape[1]} dims, "
          f"{pool['research_id'].nunique():,} people", file=sys.stderr)

    # ---------------- cluster ----------------
    t0 = time.time()
    km = MiniBatchKMeans(n_clusters=args.k, random_state=args.seed, batch_size=10_000,
                         n_init=3, max_iter=300)
    pool["cluster"] = km.fit_predict(embs)
    print(f"MiniBatchKMeans k={args.k}: {time.time() - t0:.0f}s", file=sys.stderr)
    km2 = MiniBatchKMeans(n_clusters=args.k, random_state=args.seed + 1, batch_size=10_000,
                          n_init=3, max_iter=300).fit_predict(embs)
    ct = pd.crosstab(pool["cluster"], km2).to_numpy()
    jac = ct / (ct.sum(0, keepdims=True) + ct.sum(1, keepdims=True) - ct)
    stability = float(np.mean(jac.max(axis=1)))
    print(f"Cluster stability (mean best Jaccard vs second seed): {stability:.3f}",
          file=sys.stderr)

    # ---------------- per-cluster statistics ----------------
    sample = pool.sample(min(args.pgen_sample, len(pool)), random_state=args.seed).copy()
    sample["pgen"] = R.olga_pgen(sample["cdr3aa"], args.workers)
    lp = (sample[sample["pgen"] > 0].groupby("cluster")["pgen"]
                .apply(lambda s: float(np.median(np.log10(s)))).rename("median_log10_pgen"))

    share = (pool.groupby(["trbv", "cdr3aa"], observed=True)["research_id"].transform("nunique"))
    pool["n_people_sharing"] = share
    g = pool.groupby("cluster")
    stats = pd.DataFrame({
        "n_clonotypes": g.size(),
        "n_people": g["research_id"].nunique(),
        "mean_cdr3_len": g["len"].mean(),
        "mean_log10_reads": g["reads"].apply(lambda s: float(np.mean(np.log10(np.maximum(s, 1))))),
        "public_frac": g["n_people_sharing"].apply(lambda s: float(np.mean(s >= 2))),
        "top_trbv": g["trbv"].agg(lambda s: s.value_counts().idxmax()),
        "top_trbv_frac": g["trbv"].agg(lambda s: s.value_counts(normalize=True).iat[0]),
    }).join(lp)
    stats["trbv_family"] = stats["top_trbv"].str.split("-").str[0]
    kept = stats["n_people"] >= R.MIN_PEOPLE
    print(f"{int((~kept).sum())} of {args.k} clusters pool < {R.MIN_PEOPLE} people and are "
          f"dropped", file=sys.stderr)

    # ---------------- 3-d / 2-d centroid layout ----------------
    pca = PCA(n_components=3, random_state=args.seed).fit(km.cluster_centers_)
    xyz = pca.transform(km.cluster_centers_)
    for i, c in enumerate("xyz"):
        stats[f"pc{i + 1}"] = xyz[:, i]
    var = pca.explained_variance_ratio_

    # ---------------- VDJdb specificity enrichment ----------------
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from prep_vdjdb_pool import DEFAULT_URL

    from importlib import import_module
    load_vdjdb = import_module("04_antigen_specificity").load_vdjdb
    vdj = load_vdjdb(os.path.expanduser(args.vdjdb) if args.vdjdb else None, DEFAULT_URL,
                     os.path.expanduser(args.cache_dir))
    hit = (pool[["cluster", "trbv", "cdr3aa"]]
           .merge(vdj[["trbv", "cdr3aa", "pathogen", "epitope"]].drop_duplicates(
               ["trbv", "cdr3aa"]), on=["trbv", "cdr3aa"]))
    print(f"{len(hit):,} cohort clonotypes exactly match a VDJdb TCR", file=sys.stderr)

    # Crude test: Fisher exact, cluster vs rest. V-adjusted test: Mantel-Haenszel across TRBV
    # strata, i.e. "within the same V gene, are matches over-represented in this cluster?".
    # Clusters are ~74% one V family and VDJdb's V usage is skewed, so a crude enrichment can
    # be V-gene bias; the stratified one is the CDR3-motif part.
    ct_all = pd.crosstab(pool["cluster"], pool["trbv"])
    rows = []
    tot = len(pool)
    for path, h in hit.groupby("pathogen"):
        if len(h) < MIN_MATCHES_PER_TEST:
            continue
        per = h["cluster"].value_counts()
        ct_hit = pd.crosstab(h["cluster"], h["trbv"]).reindex(
            index=ct_all.index, columns=ct_all.columns, fill_value=0)
        H_v, N_v = ct_hit.sum(0).to_numpy(), ct_all.sum(0).to_numpy()
        for cl in stats.index[kept]:
            a = int(per.get(cl, 0))
            if a == 0:
                continue
            n_cl = int(stats.loc[cl, "n_clonotypes"])
            table = [[a, len(h) - a], [n_cl - a, tot - n_cl - (len(h) - a)]]
            orr, p = fisher_exact(table, alternative="two-sided")
            se = np.sqrt(sum(1 / max(v, 0.5) for row in table for v in row))
            av, cv = ct_hit.loc[cl].to_numpy(), ct_all.loc[cl].to_numpy()
            or_mh, p_mh = R.mantel_haenszel(av, H_v - av, cv - av, N_v - cv - (H_v - av))
            rows.append((path, cl, a, len(h), n_cl, orr, np.exp(np.log(max(orr, 1e-9)) - 1.96 * se),
                         np.exp(np.log(max(orr, 1e-9)) + 1.96 * se), p, or_mh, p_mh))
    enr = pd.DataFrame(rows, columns=["pathogen", "cluster", "matches_in_cluster",
                                      "matches_total", "cluster_clonotypes", "odds_ratio",
                                      "or_lo", "or_hi", "p", "or_mh_trbv", "p_cmh_trbv"])
    if len(enr):
        enr["q"] = bh(enr["p"])
        enr["q_cmh_trbv"] = bh(enr["p_cmh_trbv"].fillna(1.0))
        enr = enr.sort_values("p")
    enr.to_csv(os.path.join(outdir, "vdjdb_cluster_enrichment.csv"), index=False)

    out = stats[kept].copy()
    out.index.name = "cluster"
    out.to_csv(os.path.join(outdir, "clusters_3d.csv"))

    # ---------------- figure ----------------
    import cnsplots as cns
    mp = cns.multipanel(max_width=540)
    pal = R.nature()
    W, H, GAP = 110, 110, 34
    s = out
    fams = s["trbv_family"].value_counts().head(8).index.tolist()
    fcol = dict(zip(fams, pal[:len(fams)]))
    size = 6 + 34 * (s["n_clonotypes"] / s["n_clonotypes"].max())

    ax = mp.panel("a", width=W, height=H, margin_right=70, margin_bottom=GAP)
    for f in fams + ["other"]:
        m = s["trbv_family"] == f if f != "other" else ~s["trbv_family"].isin(fams)
        if m.any():
            ax.scatter(s.loc[m, "pc1"], s.loc[m, "pc2"], s=size[m], linewidths=0,
                       color=fcol.get(f, "0.7"), label=f.replace("TRBV", "V"), alpha=0.9)
    ax.legend(frameon=False, markerscale=0.6)
    cns.take_legend_out(title="Dominant TRBV", ax=ax)
    ax.set_title(f"{len(s)} metaclusters", pad=11)
    ax.text(0.5, 1.02, f"dot area ∝ clonotypes; PC1–2 = {100 * var[:2].sum():.0f}% of "
            f"centroid variance", transform=ax.transAxes, ha="center", va="bottom",
            fontsize=5.5)
    R.style(ax, None, "Centroid PC1", "Centroid PC2")

    ax = mp.panel("b", width=W, height=H, margin_right=55, margin_bottom=GAP)
    sc = ax.scatter(s["mean_cdr3_len"], s["median_log10_pgen"], s=size,
                    c=100 * s["public_frac"], cmap="parula", linewidths=0)
    R.colorbar(ax, sc, "Shared by ≥2 people (%)")
    R.style(ax, "What separates clusters", "Mean CDR3 length (aa)",
            "Median log$_{10}$ P$_{gen}$")

    ax = mp.panel("c", width=W, height=H, margin_right=15)
    if len(enr):
        sig = enr[enr["q"] < 0.05]
        ax.scatter(np.log2(enr["odds_ratio"].clip(0.05, 20)), -np.log10(enr["p"].clip(1e-300)),
                   s=8, linewidths=0, color="0.75", label="n.s.")
        if len(sig):
            ax.scatter(np.log2(sig["odds_ratio"].clip(0.05, 20)),
                       -np.log10(sig["p"].clip(1e-300)), s=10, linewidths=0, color=pal[0],
                       label="BH q < 0.05")
        ax.axvline(0, color="0.6", lw=0.5, ls="--")
        ax.legend(frameon=False, fontsize=5.5, loc="upper left")
        ax.set_title("Specificity enrichment", pad=11)
        ax.text(0.5, 1.02, f"{len(enr):,} cluster × pathogen tests, "
                f"{len(sig):,} significant", transform=ax.transAxes, ha="center",
                va="bottom", fontsize=5.5)
    R.style(ax, None, "log$_2$ odds ratio", "−log$_{10}$ P")

    ax = mp.panel("d", width=190, height=H, margin_right=15)
    if len(enr):
        top = enr[enr["odds_ratio"] > 1].head(TOP_PAIRS).iloc[::-1]
        y = np.arange(len(top))
        pc = {p: c for p, c in zip(sorted(top["pathogen"].unique()), pal)}
        for i, r in enumerate(top.itertuples()):
            c = pc[r.pathogen]
            ax.plot([r.or_lo, r.or_hi], [i, i], color=c, lw=1.0)
            ax.scatter([r.odds_ratio], [i], s=16, color=c, zorder=3)
            if np.isfinite(r.or_mh_trbv) and r.or_mh_trbv > 0:
                ax.scatter([r.or_mh_trbv], [i], s=14, facecolors="none", edgecolors=c,
                           marker="D", linewidths=0.7, zorder=4)
        ax.scatter([], [], s=16, color="0.3", label="crude")
        ax.scatter([], [], s=14, facecolors="none", edgecolors="0.3", marker="D",
                   linewidths=0.7, label="within TRBV gene")
        ax.legend(frameon=False, fontsize=5, loc="lower right")
        ax.axvline(1, color="0.6", lw=0.5, ls="--")
        ax.set_xscale("log")
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r.pathogen} · cluster {r.cluster} "
                            f"({r.matches_in_cluster}/{r.matches_total})"
                            for r in top.itertuples()], fontsize=5.5)
        ax.minorticks_off()
    R.style(ax, "Strongest enrichments", "Odds ratio (95% CI)", None)
    R.save("fig_embedding_clusters", outdir)

    rows = [("clonotypes", len(pool)), ("k", args.k), ("clusters_kept", int(kept.sum())),
            ("cluster_stability_jaccard", stability),
            ("centroid_pc1_var", var[0]), ("centroid_pc2_var", var[1]),
            ("centroid_pc3_var", var[2]),
            ("median_cluster_size", float(s["n_clonotypes"].median())),
            ("mean_top_trbv_frac", float(s["top_trbv_frac"].mean())),
            ("vdjdb_exact_matched_clonotypes", len(hit)),
            ("enrichment_tests", len(enr)),
            ("enrichment_significant_q05", int((enr["q"] < 0.05).sum()) if len(enr) else 0),
            ("enrichment_significant_q05_within_trbv",
             int((enr["q_cmh_trbv"] < 0.05).sum()) if len(enr) else 0),
            ("enrichment_significant_both",
             int(((enr["q"] < 0.05) & (enr["q_cmh_trbv"] < 0.05)).sum()) if len(enr) else 0)]
    if len(enr):
        for path in ("HomoSapiens", "CMV", "EBV", "InfluenzaA"):
            e_ = enr[(enr["pathogen"] == path) & (enr["odds_ratio"] > 1)]
            if len(e_):
                top_ = e_.sort_values("p").iloc[0]
                rows += [(f"top_or_{path}", top_["odds_ratio"]),
                         (f"top_or_{path}_within_trbv", top_["or_mh_trbv"]),
                         (f"top_q_{path}_within_trbv", top_["q_cmh_trbv"])]
    summary = pd.DataFrame(rows, columns=["metric", "value"])
    summary.to_csv(os.path.join(outdir, "summary.csv"), index=False)
    print("\n" + summary.to_string(index=False))
    print(f"\nAll outputs in {outdir}")


if __name__ == "__main__":
    main()
