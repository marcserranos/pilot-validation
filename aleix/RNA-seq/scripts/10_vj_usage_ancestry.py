#!/usr/bin/env python3
"""10 -- W1: TRB V and J gene usage by genetic ancestry (report: reports/10_vj_usage_ancestry/).

Cole's request (call 2026-10-06): a model-free view of whether repertoire composition
differs by ancestry. Per person, the frequency of each TRBV and TRBJ gene among their unique
TRB clonotypes, then PCA, coloured by ancestry. Plus three numbers that make it more than a
picture: how much ancestry each principal component carries, which genes differ most, and
how well V/J usage predicts ancestry against a permuted control and a covariate-only model.

  a  PC1 vs PC2 of V+J usage, ancestry drawn as 50% / 80% density contours (people with
     ancestry probability >= 0.90; no individual drawn)
  b  per component: variance explained, and the share of it explained by ancestry (eta^2),
     next to its correlation with depth and age -- a PC that tracks depth is not ancestry
  c  the genes that differ most by ancestry: median usage per group, scaled per gene
  d  ancestry prediction AUROC from V usage, J usage, V+J, a covariate-only model (age, sex,
     log depth) and permuted labels; grouped CV by family, so relatives never straddle folds

CHOICES (DECISIONS.md, 2026-10-08):
  - Unit = unique clonotypes, not reads. Read counts measure expansion and transcription;
    clonotype counts measure which genes recombine and survive selection, which is what a
    germline or HLA effect on usage would change.
  - Common genes only: mean usage >= 0.1% and present in >= 50% of people. Rare and
    unusable symbols are reported, not modelled.
  - Compositional data: centred log-ratio per composition (V separately from J), counts
    + 0.5 pseudocount, then PCA on the centred CLR matrix without scaling (Aitchison PCA).
  - Subcohorts: features and PCA on everyone passing the read floor (7,527); the figure on
    the >= 0.90-purity subset; prediction on read floor + unrelated, grouped by family.
  - Compositional echo: usage shares sum to one, so one gene that is strongly up in a group
    pushes every other gene of the same composition down a little in that group. A gene
    list ranked by effect size therefore mixes primary effects with their echoes; read the
    top of the list, not the tail (seen on planted test data: +25% TRBV9 in one group made
    every other TRBV look depleted there).
  - Germline caveat: TRBV deletion/loss-of-function alleles differ in frequency across
    populations (Mantena et al. 2026), so V-usage differences by ancestry can be genetic
    without involving HLA. This report describes; it does not attribute.

Inputs (VM): analysis/person_table.pkl (report 09), <research_id>/<rid>_cdr3.out.
Outputs: aggregate CSVs + figure to --outdir; VM-local analysis/vj_usage.pkl (person x gene
counts, reused by W3 feature sets).

Usage (from aleix/RNA-seq/, after 09):
  pixi run python3 -u scripts/10_vj_usage_ancestry.py [--workers 16]
"""
import argparse
import importlib
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import repfig as R  # noqa: E402
from embed_cdr3s import load_person_cdr3s  # noqa: E402

ANC = R.ANCESTRY_ORDER
MIN_MEAN_USAGE = 0.001
MIN_PRESENT = 0.5
N_PCS = 10


def gene(g):
    return str(g).split("*")[0] if g and str(g) not in ("nan", "") else ""


def usage_one(job):
    trust4_dir, rid = job
    df, _ = load_person_cdr3s(trust4_dir, rid, 0.02, False, False, 30, keep_j=True)
    if df is None:
        return rid, None, None
    d = df[df["chain"] == "TRB"]
    v = d["v_gene"].map(gene)
    j = d["j_gene"].map(gene)
    return rid, v[v.str.startswith("TRBV")].value_counts(), j[j.str.startswith("TRBJ")].value_counts()


def clr_counts(C, pseudo=0.5):
    L = np.log(C + pseudo)
    return L - L.mean(axis=1, keepdims=True)


def eta2(y, groups):
    """Share of variance in y explained by group membership (one-way ANOVA R^2)."""
    y = np.asarray(y, float)
    g = pd.Series(y).groupby(np.asarray(groups))
    return float(sum(len(v) * (v.mean() - y.mean()) ** 2 for _, v in g) / np.sum((y - y.mean()) ** 2))


def cv_ancestry_auroc(X, y, groups, seed, folds=5):
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedGroupKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    classes = [a for a in ANC if np.sum(y == a) >= R.MIN_PEOPLE]
    ok = np.isin(y, classes)
    X, y, groups = X[ok], y[ok], groups[ok]
    proba = np.zeros((len(y), len(classes)))
    cv = StratifiedGroupKFold(folds, shuffle=True, random_state=seed)
    for tr, te in cv.split(X, y, groups):
        m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=3000)).fit(X[tr], y[tr])
        proba[te] = m.predict_proba(X[te])[:, [list(m.classes_).index(c) for c in classes]]
    per = {c: roc_auc_score(y == c, proba[:, i]) for i, c in enumerate(classes)}
    return float(np.mean(list(per.values()))), per


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    base = "~/pipeline_outputs/rnaseq"
    ap.add_argument("--person-table", default=f"{base}/analysis/person_table.pkl")
    ap.add_argument("--trust4-dir", default=base)
    ap.add_argument("--local-out", default=f"{base}/analysis")
    ap.add_argument("--outdir", default=f"{base}/reports/10_vj_usage_ancestry")
    ap.add_argument("--purity", type=float, default=0.90)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--reread", action="store_true",
                    help="rebuild usage counts from cdr3.out even if analysis/vj_usage.pkl exists")
    ap.add_argument("--indel-genes", default="TRBV4-3,TRBV6-2,TRBV3-2",
                    help="genes in the common TRB-locus insertion/deletion polymorphism; the "
                         "sensitivity PCA (panel e) is recomputed without them")
    args = ap.parse_args()
    outdir, local = os.path.expanduser(args.outdir), os.path.expanduser(args.local_out)
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    pt = pd.read_pickle(os.path.expanduser(args.person_table))
    pt = pt[pt["pass_read_floor"]].reset_index(drop=True)
    print(f"Read floor subcohort: {len(pt):,} people", file=sys.stderr)

    # ---------------- usage counts ----------------
    cache = os.path.join(local, "vj_usage.pkl")
    if os.path.exists(cache) and not args.reread:
        uc = pd.read_pickle(cache)
        Vc = uc["V"].reindex(pt["research_id"]).fillna(0)
        Jc = uc["J"].reindex(pt["research_id"]).fillna(0)
        print(f"Loaded usage counts from {cache} (--reread to rebuild)", file=sys.stderr)
    else:
        t0 = time.time()
        jobs = [(os.path.expanduser(args.trust4_dir), r) for r in pt["research_id"]]
        V, J = {}, {}
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            for i, (rid, v, j) in enumerate(ex.map(usage_one, jobs, chunksize=16), 1):
                if v is not None:
                    V[rid], J[rid] = v, j
                if i % 1000 == 0:
                    print(f"  read {i:,}/{len(jobs):,} ({time.time() - t0:.0f}s)", file=sys.stderr)
        Vc = pd.DataFrame(V).T.fillna(0).reindex(pt["research_id"]).fillna(0)
        Jc = pd.DataFrame(J).T.fillna(0).reindex(pt["research_id"]).fillna(0)
        pd.to_pickle({"V": Vc, "J": Jc}, cache)

    def common(C):
        P = C.div(C.sum(1).replace(0, np.nan), axis=0)
        keep = (P.mean() >= MIN_MEAN_USAGE) & ((C > 0).mean() >= MIN_PRESENT)
        return C.loc[:, keep].reindex(sorted(C.columns[keep], key=natkey), axis=1), int((~keep).sum())

    def natkey(g):
        import re
        return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", g)]

    Vk, v_drop = common(Vc)
    Jk, j_drop = common(Jc)
    print(f"Common genes: {Vk.shape[1]} TRBV ({v_drop} rare/unused dropped), "
          f"{Jk.shape[1]} TRBJ ({j_drop} dropped)", file=sys.stderr)
    CV, CJ = clr_counts(Vk.to_numpy(float)), clr_counts(Jk.to_numpy(float))
    X = np.hstack([CV, CJ])
    genes = list(Vk.columns) + list(Jk.columns)
    anc = pt["ancestry"].to_numpy()

    # ---------------- PCA ----------------
    from sklearn.decomposition import PCA
    pca = PCA(n_components=N_PCS, random_state=args.seed).fit(X - X.mean(0))
    S = pca.transform(X - X.mean(0))
    from scipy.stats import spearmanr
    lr = np.log10(pt["trb_reads"].to_numpy(float))
    rows = []
    for k in range(N_PCS):
        rows.append({"pc": f"PC{k + 1}", "var_explained": pca.explained_variance_ratio_[k],
                     "ancestry_eta2": eta2(S[:, k], anc),
                     "spearman_log_reads": spearmanr(S[:, k], lr)[0],
                     "spearman_age": spearmanr(S[:, k], pt["age"], nan_policy="omit")[0],
                     "spearman_aou_PC1": spearmanr(S[:, k], pt["PC1"], nan_policy="omit")[0]
                     if "PC1" in pt else np.nan})
    pcs = pd.DataFrame(rows)
    pcs.to_csv(os.path.join(outdir, "pca_components.csv"), index=False)
    load = pd.DataFrame(pca.components_[:4].T, index=genes, columns=["PC1", "PC2", "PC3", "PC4"])
    load.to_csv(os.path.join(outdir, "pca_loadings.csv"))
    print(pcs.round(3).to_string(index=False), file=sys.stderr)

    # ---------------- genes by ancestry ----------------
    from scipy.stats import kruskal
    P_all = pd.concat([Vk.div(Vk.sum(1), axis=0), Jk.div(Jk.sum(1), axis=0)], axis=1)
    g_rows = []
    for gi, g in enumerate(genes):
        x = X[:, gi]
        grp = [x[anc == a] for a in ANC if np.sum(anc == a) >= R.MIN_PEOPLE]
        row = {"gene": g, "eta2": eta2(x, anc), "kruskal_p": kruskal(*grp).pvalue,
               "median_usage_pct_all": 100 * P_all[g].median()}
        for a in ANC:
            m = anc == a
            if m.sum() >= R.MIN_PEOPLE:
                row[f"median_usage_pct_{a}"] = 100 * P_all.loc[m, g].median()
        g_rows.append(row)
    gt = pd.DataFrame(g_rows)
    from statsmodels.stats.multitest import multipletests
    gt["kruskal_q"] = multipletests(gt["kruskal_p"], method="fdr_bh")[1]
    gt = gt.sort_values("eta2", ascending=False)
    gt.to_csv(os.path.join(outdir, "genes_by_ancestry.csv"), index=False)

    # ---------------- genotype signature: people who never use a gene ----------------
    # Zero usage happens by chance when a person has few clonotypes and the gene is rare. The
    # expected zero fraction if everyone carried the gene: mean over people of (1 - p)^n_i,
    # p = the gene's median usage among people who use it, n_i = the person's clonotype count.
    # Zeros well above that expectation point to people lacking the gene (germline).
    z_rows = []
    for comp, C in (("V", Vc), ("J", Jc)):
        tot = C.sum(1).to_numpy(float)
        for g in C.columns:
            c = C[g].to_numpy(float)
            present = c > 0
            if present.sum() < R.MIN_PEOPLE:
                continue
            pg = float(np.median(c[present] / tot[present]))
            exp0 = (1 - pg) ** tot
            row = {"gene": g, "median_usage_pct_if_used": 100 * pg,
                   "obs_zero_frac": float(np.mean(~present)), "exp_zero_frac": float(np.mean(exp0))}
            for a in ANC:
                m = anc == a
                if m.sum() >= R.MIN_PEOPLE:
                    row[f"obs_zero_{a}"] = float(np.mean(~present[m]))
                    row[f"exp_zero_{a}"] = float(np.mean(exp0[m]))
            z_rows.append(row)
    zt = pd.DataFrame(z_rows)
    zt["excess_zero_frac"] = zt["obs_zero_frac"] - zt["exp_zero_frac"]
    # if zeros are homozygous absence and alleles are in Hardy-Weinberg equilibrium,
    # absence-allele frequency q ~ sqrt(excess zero fraction)
    zt["implied_absence_allele_freq"] = np.sqrt(zt["excess_zero_frac"].clip(lower=0))
    zt = zt.sort_values("excess_zero_frac", ascending=False)
    zt.to_csv(os.path.join(outdir, "zero_usage_by_ancestry.csv"), index=False)
    print("Genes with excess zero usage (observed - expected if everyone carried them):\n"
          + zt.head(8)[["gene", "obs_zero_frac", "exp_zero_frac", "excess_zero_frac",
                        "implied_absence_allele_freq"]].round(3).to_string(index=False),
          file=sys.stderr)

    # ---------------- sensitivity: PCA without the indel genes ----------------
    indel = [g for g in args.indel_genes.split(",") if g in Vk.columns]
    Vk2 = Vk.drop(columns=indel)
    X2 = np.hstack([clr_counts(Vk2.to_numpy(float)), CJ])
    pca2 = PCA(n_components=N_PCS, random_state=args.seed).fit(X2 - X2.mean(0))
    S2 = pca2.transform(X2 - X2.mean(0))
    pcs2 = pd.DataFrame({"pc": [f"PC{k + 1}" for k in range(N_PCS)],
                         "var_explained": pca2.explained_variance_ratio_,
                         "ancestry_eta2": [eta2(S2[:, k], anc) for k in range(N_PCS)]})
    pcs2.to_csv(os.path.join(outdir, "pca_components_without_indel_genes.csv"), index=False)

    # ---------------- ancestry prediction ----------------
    sub = pt["unrelated"].to_numpy()
    groups = pt["family_id"].to_numpy()
    cov = np.column_stack([pt["age"].fillna(pt["age"].median()), lr,
                           (pt["sex_model"] == "Female").astype(float),
                           (pt["sex_model"] == "Unknown").astype(float)])
    sets = {"TRBV usage": CV, "TRBJ usage": CJ, "TRBV + TRBJ usage": X,
            "covariates only (age, sex, depth)": cov,
            "TRBV + TRBJ + covariates": np.hstack([X, cov]),
            "TRBV + TRBJ without indel genes": X2}
    a_rows = []
    for name, F in sets.items():
        macro, per = cv_ancestry_auroc(F[sub], anc[sub], groups[sub], args.seed)
        a_rows.append({"features": name, "macro_auroc": macro, **{f"auroc_{k}": v for k, v in per.items()}})
        print(f"  {name:36s} macro AUROC {macro:.3f}", file=sys.stderr)
    perm = []
    for r in range(5):
        macro, _ = cv_ancestry_auroc(X[sub], rng.permutation(anc[sub]), groups[sub], args.seed + r)
        perm.append(macro)
    a_rows.append({"features": "TRBV + TRBJ, permuted labels (mean of 5)", "macro_auroc": float(np.mean(perm))})
    au = pd.DataFrame(a_rows)
    au.to_csv(os.path.join(outdir, "ancestry_auroc.csv"), index=False)

    # ---------------- figure ----------------
    import cnsplots as cns
    m01 = importlib.import_module("01_sceptr_embedding_viz")
    col = R.ancestry_colors()
    mp = cns.multipanel(max_width=540)
    W, H, GAP = 120, 120, 36

    pure = (pt["anc_maxprob"] >= args.purity).to_numpy() if "anc_maxprob" in pt else np.ones(len(pt), bool)
    ax = mp.panel("a", width=W, height=H, margin_right=75, margin_bottom=GAP)
    m01.hdr_contours(ax, S[pure][:, :2], anc[pure], ANC, [col[a] for a in ANC],
                     f"V+J usage, ancestry prob ≥ {args.purity:.2f}", "Ancestry")
    ax.set_xlabel(f"PC1 ({100 * pca.explained_variance_ratio_[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({100 * pca.explained_variance_ratio_[1]:.1f}%)")

    ax = mp.panel("b", width=W, height=H, margin_right=20, margin_bottom=GAP)
    x = np.arange(N_PCS)
    ax.bar(x - 0.2, 100 * pcs["var_explained"], width=0.4, color="0.7", lw=0, label="variance explained")
    ax.bar(x + 0.2, 100 * pcs["var_explained"] * pcs["ancestry_eta2"], width=0.4,
           color=R.nature()[0], lw=0, label="of which ancestry")
    for k in range(N_PCS):
        if abs(pcs.loc[k, "spearman_log_reads"]) >= 0.3 or abs(pcs.loc[k, "spearman_age"]) >= 0.3:
            ax.text(k, 100 * pcs.loc[k, "var_explained"] + 0.3, "d" if abs(pcs.loc[k, "spearman_log_reads"]) >= 0.3 else "a",
                    ha="center", fontsize=5.5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{k + 1}" for k in x])
    ax.legend(frameon=False, fontsize=5.5)
    ax.set_title("What each component carries", pad=11)
    ax.text(0.5, 1.02, "letters: |ρ| ≥ 0.3 with log depth (d) or age (a)", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=5.5)
    R.style(ax, None, "Principal component", "% of total variance")

    ax = mp.panel("c", width=W + 10, height=H + 20, margin_right=60, margin_bottom=GAP)
    top = gt.head(20)
    M = np.array([[top.iloc[i].get(f"median_usage_pct_{a}", np.nan) for a in ANC] for i in range(len(top))])
    allmed = top["median_usage_pct_all"].to_numpy()[:, None]
    Z = np.log2(np.maximum(M, 1e-6) / np.maximum(allmed, 1e-6))   # fold change vs whole cohort
    lim = float(min(1.0, max(0.25, np.nanquantile(np.abs(Z), 0.98))))
    im = ax.imshow(Z, cmap="bwr", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(len(ANC)))
    ax.set_xticklabels(ANC)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([f"{g}  (η²={e:.2f})" for g, e in zip(top["gene"], top["eta2"])], fontsize=4.5)
    R.colorbar(ax, im, "log$_2$ median usage vs cohort")
    ax.set_title("Genes that differ most", pad=4)

    ax = mp.panel("d", width=W, height=H, margin_right=12, margin_bottom=GAP)
    order = list(au["features"])
    y = np.arange(len(order))[::-1]
    pal = R.nature()
    fcol = {"TRBV usage": pal[0], "TRBJ usage": pal[1], "TRBV + TRBJ usage": pal[3],
            "covariates only (age, sex, depth)": "0.6", "TRBV + TRBJ + covariates": pal[2],
            "TRBV + TRBJ without indel genes": pal[4]}
    ax.barh(y, au["macro_auroc"] - 0.5, left=0.5, height=0.6,
            color=[fcol.get(f, "0.85") for f in order], lw=0)
    ax.axvline(0.5, color="0.5", lw=0.5, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=5)
    ax.set_xlim(0.45, 1.0)
    R.style(ax, "Predicting ancestry", "Macro one-vs-rest AUROC", None)

    ax = mp.panel("e", width=W, height=H, margin_right=75, margin_bottom=GAP)
    m01.hdr_contours(ax, S2[pure][:, :2], anc[pure], ANC, [col[a] for a in ANC],
                     "Without " + ", ".join(indel), "Ancestry")
    ax.set_xlabel(f"PC1 ({100 * pca2.explained_variance_ratio_[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({100 * pca2.explained_variance_ratio_[1]:.1f}%)")

    ax = mp.panel("f", width=W, height=H, margin_right=12, margin_bottom=GAP)
    show = [g for g in ("TRBV4-3", "TRBV6-2", "TRBV3-2") if g in set(zt["gene"])][:2]
    if not show:
        show = list(zt["gene"].head(2))
    w = 0.38
    for k, g in enumerate(show):
        r = zt.set_index("gene").loc[g]
        xs = np.arange(len(ANC)) + (k - 0.5) * w
        obs = [100 * r.get(f"obs_zero_{a}", np.nan) for a in ANC]
        exp_ = [100 * r.get(f"exp_zero_{a}", np.nan) for a in ANC]
        ax.bar(xs, obs, width=w, color=pal[k], lw=0, label=f"{g} observed")
        ax.scatter(xs, exp_, s=10, color="black", zorder=3, label="expected by chance" if k == 0 else None)
    ax.set_xticks(range(len(ANC)))
    ax.set_xticklabels(ANC)
    ax.legend(frameon=False, fontsize=5.5)
    R.style(ax, "People who never use the gene", None, "% of people")
    R.save("fig_vj_usage_ancestry", outdir)

    summary = [("people_read_floor", len(pt)), ("people_unrelated_for_auroc", int(sub.sum())),
               ("people_pure_figure", int(pure.sum())), ("trbv_common", Vk.shape[1]),
               ("trbj_common", Jk.shape[1]), ("trbv_dropped", v_drop), ("trbj_dropped", j_drop)]
    summary += [(f"auroc_{r.features}", r.macro_auroc) for r in au.itertuples()]
    summary += [(f"{r.pc}_var", r.var_explained) for r in pcs.itertuples()]
    summary += [(f"{r.pc}_ancestry_eta2", r.ancestry_eta2) for r in pcs.itertuples()]
    pd.DataFrame(summary, columns=["metric", "value"]).to_csv(os.path.join(outdir, "summary.csv"), index=False)
    print(f"\n{au.round(3).to_string(index=False)}\n\nTop genes by eta2:\n"
          f"{gt.head(12)[['gene', 'eta2', 'kruskal_q']].round(4).to_string(index=False)}\n"
          f"All outputs in {outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
