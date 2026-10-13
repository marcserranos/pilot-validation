#!/usr/bin/env python3
"""14 -- W5: disease classification v0 from the TRB repertoire
(report: reports/14_disease_classification/).

Cole's request (call 2026-10-06): TCR-beta only, SCEPTR embeddings + k-mer features,
XGBoost and LASSO, "get something working". The question is not "can we predict disease"
(age alone predicts most chronic disease) but "what does the repertoire ADD to age, sex,
ancestry and depth". So every model is reported as AUROC of base + repertoire features
against AUROC of base alone, on the same people and the same folds.

Diseases: >= 150 cases in report 11. Feature sets (report 12): base; base + diversity;
base + V/J usage; base + SCEPTR (mean + TRBV); base + 3-mers; base + all.
Models: L1-penalised logistic regression (penalty chosen by inner 3-fold CV) and XGBoost
(fixed, conservative settings: 300 trees, depth 3, learning rate 0.05, row and column
subsampling 0.8; no tuning in v0). Evaluation: 5-fold stratified CV, out-of-fold
predictions, AUROC; the gain over base gets a 95% CI from a paired bootstrap over people
(500 resamples of the out-of-fold predictions). analysis_main is unrelated, so plain
stratified folds cannot leak relatives.

Usage (from aleix/RNA-seq/, after 11 and 12; needs `pip install xgboost`):
  pixi run python3 -u scripts/14_disease_classification.py [--workers 16]
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import repfig as R  # noqa: E402

MIN_CASES = 150
N_BOOT = 500
SETS = [("base", ["base"]), ("+ diversity", ["base", "diversity"]), ("+ V/J usage", ["base", "vj"]),
        ("+ SCEPTR", ["base", "sceptr"]), ("+ 3-mers", ["base", "kmer"]),
        ("+ all", ["base", "diversity", "vj", "sceptr", "kmer"])]


def models(seed, workers):
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegressionCV
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    out = {"LASSO": lambda: make_pipeline(
        SimpleImputer(strategy="median"), StandardScaler(),
        LogisticRegressionCV(Cs=10, cv=3, penalty="l1", solver="liblinear", scoring="roc_auc",
                             max_iter=2000, random_state=seed, n_jobs=workers))}
    try:
        from xgboost import XGBClassifier
        out["XGBoost"] = lambda: XGBClassifier(n_estimators=300, max_depth=3, learning_rate=0.05,
                                               subsample=0.8, colsample_bytree=0.8, n_jobs=workers,
                                               random_state=seed, eval_metric="logloss")
    except ImportError:
        from sklearn.ensemble import HistGradientBoostingClassifier
        print("!! xgboost not installed: using sklearn HistGradientBoosting as the tree model",
              file=sys.stderr)
        out["GradientBoosting"] = lambda: HistGradientBoostingClassifier(
            max_iter=300, max_depth=3, learning_rate=0.05, random_state=seed)
    return out


def plot(res, outdir):
    """Figure from the results table alone, so it can be redrawn without refitting:
    a  out-of-fold AUROC, base vs base + all features, per disease (LASSO)
    b  gain over base in AUROC points, disease x feature set, one column block per model;
       ring = 95% CI of the gain excludes 0"""
    import cnsplots as cns
    pal = R.nature()
    mods = list(dict.fromkeys(res["model"]))
    sets = [s_ for s_, _ in SETS if s_ != "base"]
    base_ = res[(res["features"] == "base") & (res["model"] == mods[0])].sort_values("auroc")
    dis = list(base_["disease"])
    H = max(110, 10 * len(dis))
    mp = cns.multipanel(max_width=540)
    ax = mp.panel("a", width=110, height=H, margin_right=12, margin_bottom=48)
    y = np.arange(len(dis))
    for sname, c in (("base", "0.6"), ("+ all", pal[0])):
        s_ = res[(res["features"] == sname) & (res["model"] == mods[0])].set_index("disease").loc[dis]
        ax.scatter(s_["auroc"], y, s=12, color=c, zorder=3, label=f"{sname} ({mods[0]})")
    ax.axvline(0.5, color="0.5", lw=0.5, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{d} ({c})" for d, c in zip(dis, base_["cases"])], fontsize=4.8)
    ax.legend(frameon=False, fontsize=5, loc="lower right")
    R.style(ax, "Out-of-fold AUROC", "AUROC", None)

    cols = [(m, s_) for m in mods for s_ in sets]
    G = np.array([[res[(res["disease"] == d) & (res["model"] == m) & (res["features"] == s_)]["gain_over_base"].iat[0]
                   for m, s_ in cols] for d in dis]) * 100
    sig = np.array([[bool(res[(res["disease"] == d) & (res["model"] == m) & (res["features"] == s_)]["gain_ci_excludes_0"].iat[0])
                     for m, s_ in cols] for d in dis])
    ax = mp.panel("b", width=170, height=H, margin_left=95, margin_right=60, margin_bottom=48)
    lim = max(1.0, float(np.nanquantile(np.abs(G), 0.95)))
    im = ax.imshow(G, cmap="bwr", vmin=-lim, vmax=lim, aspect="auto", origin="lower")  # row 0 at the bottom, as in panel a
    yy, xx = np.nonzero(sig)
    ax.scatter(xx, yy, s=16, facecolors="none", edgecolors="black", linewidths=0.6)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([s_ for _, s_ in cols], rotation=60, ha="right", fontsize=5)
    for k, m in enumerate(mods):
        ax.text(k * len(sets) + (len(sets) - 1) / 2, len(dis) - 0.3, m, ha="center", va="bottom", fontsize=6)
        if k:
            ax.axvline(k * len(sets) - 0.5, color="black", lw=0.6)
    ax.set_yticks(range(len(dis)))
    ax.set_yticklabels(dis, fontsize=4.8)
    R.colorbar(ax, im, "gain over base (AUROC points)")
    ax.set_title("What the repertoire adds", pad=14)
    R.save("fig_disease_classification", outdir)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    base = "~/pipeline_outputs/rnaseq"
    ap.add_argument("--labels", default=f"{base}/analysis/disease_labels.pkl")
    ap.add_argument("--features", default=f"{base}/analysis/feature_sets.pkl")
    ap.add_argument("--counts", default=f"{base}/reports/11_disease_labels/case_counts.csv")
    ap.add_argument("--outdir", default=f"{base}/reports/14_disease_classification")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--plot-only", action="store_true",
                    help="redraw the figure from <outdir>/classification.csv without refitting")
    args = ap.parse_args()
    if args.plot_only:
        plot(pd.read_csv(os.path.join(os.path.expanduser(args.outdir), "classification.csv")),
             os.path.expanduser(args.outdir))
        return
    outdir = os.path.expanduser(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold

    L = pd.read_pickle(os.path.expanduser(args.labels))
    F = pd.read_pickle(os.path.expanduser(args.features))
    names = pd.read_csv(os.path.expanduser(args.counts)).set_index("key")["disease"].to_dict()
    M = models(args.seed, args.workers)
    rng = np.random.default_rng(args.seed)
    rows = []
    t_all = time.time()
    for key in L.columns:
        lab = L[key]
        ok = lab.notna()
        y = lab[ok].to_numpy(int)
        if y.sum() < MIN_CASES:
            continue
        ids = lab.index[ok]
        folds = list(StratifiedKFold(5, shuffle=True, random_state=args.seed).split(np.zeros(len(y)), y))
        boot = [rng.integers(len(y), size=len(y)) for _ in range(N_BOOT)]
        t0 = time.time()
        for mname, make in M.items():
            oof = {}
            for sname, parts in SETS:
                X = pd.concat([F[p].reindex(ids) for p in parts], axis=1).to_numpy(float)
                pred = np.zeros(len(y))
                for tr, te in folds:
                    m = make().fit(X[tr], y[tr])
                    pred[te] = m.predict_proba(X[te])[:, 1]
                oof[sname] = pred
            a_base = roc_auc_score(y, oof["base"])
            for sname, _ in SETS:
                a = roc_auc_score(y, oof[sname])
                bs, bd = [], []
                for b in boot:
                    yb = y[b]
                    if 0 < yb.sum() < len(yb):
                        ab = roc_auc_score(yb, oof[sname][b])
                        bs.append(ab)
                        bd.append(ab - roc_auc_score(yb, oof["base"][b]))
                rows.append({"key": key, "disease": names.get(key, key), "model": mname,
                             "features": sname, "cases": int(y.sum()), "controls": int((y == 0).sum()),
                             "auroc": a, "auroc_lo": np.quantile(bs, .025), "auroc_hi": np.quantile(bs, .975),
                             "gain_over_base": a - a_base, "gain_lo": np.quantile(bd, .025),
                             "gain_hi": np.quantile(bd, .975)})
        print(f"  {names.get(key, key):40s} cases {int(y.sum()):5d}  ({time.time() - t0:.0f}s)", file=sys.stderr)
    res = pd.DataFrame(rows)
    res["gain_ci_excludes_0"] = res["gain_lo"] > 0
    res.to_csv(os.path.join(outdir, "classification.csv"), index=False)

    plot(res, outdir)

    best = (res[res["features"] != "base"].sort_values("gain_over_base", ascending=False)
               .groupby("disease").head(1).sort_values("gain_over_base", ascending=False))
    print("\nBase AUROC and best gain per disease:\n" + best.merge(
        res[(res["features"] == "base")][["disease", "model", "auroc"]].rename(columns={"auroc": "base_auroc"}),
        on=["disease", "model"])[["disease", "cases", "model", "base_auroc", "features", "auroc",
                                  "gain_over_base", "gain_lo", "gain_hi"]].round(3).to_string(index=False),
          file=sys.stderr)
    pd.DataFrame([("diseases", res["key"].nunique()), ("models", ",".join(dict.fromkeys(res["model"]))),
                  ("gains_ci_excludes_0", int(res["gain_ci_excludes_0"].sum())),
                  ("runtime_s", round(time.time() - t_all))],
                 columns=["metric", "value"]).to_csv(os.path.join(outdir, "summary.csv"), index=False)


if __name__ == "__main__":
    main()
