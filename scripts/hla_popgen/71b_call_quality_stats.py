#!/usr/bin/env python3
"""S06 71b: how well does each per-call quality metric predict typing error, and what does pruning the worst calls do? (VM)

Reads 71's INTERNAL tables (call_table, records_trio / pair / twin) and writes AGGREGATE tables with fine resolution
(rolling windows, 0.25% pruning steps, 2.5-percentile grids), from which 72 draws the figures. Nothing person-level
leaves the Workbench: every row here is a bin, a window or a summary.

Error label for the trios = the Mendelian mismatch only: the parent that yak says transmitted the child's haplotype does not
carry the child's 2-field allele. Whose metrics to attribute it to is a modelling choice, so three views are kept: the
child's own call (`trio_child`, primary), the transmitting parent's call (`trio_parent`, if parent metrics predict the mismatch
as well as the child's, attribution is ambiguous) and the worse of the two (`trio_worst`).

Views of a truth set. Each record has the oriented metrics of the child (`c_`), the other person (`p_`: transmitting
parent, or the relative) and the worse of the two (`w_`; higher = worse). Trio records are analysed on the child's own call
(`trio_child`, the call being judged) and on the worse of child and parent (`trio_worst`, because a Mendelian mismatch can
come from the parent's call); pair and twin records on the worse of the two people (`w_`).

1. stat_error_baseline.tsv  records, errors and error rate per truth view x gene
2. stat_curves.tsv          error rate in rolling windows along each metric (window = n/25 records, step = n/300; discrete
                            metrics: one point per value), with Wilson 95% intervals, x in the metric's own units
3. stat_auc.tsv             per truth view x metric: AUC of the metric for error (cluster bootstrap over people, 300
                            draws, 95% CI), within-gene AUC (genes weighted by their positive x negative pairs), Spearman rho
                            with the 0/1 error, logistic odds ratio per SD adjusted for gene and platform with person-
                            clustered robust SE. Also the multi-metric models (below) and the old 0-7 score.
4. stat_prune.tsv           pruning: remove the worst x% of records by the metric (x = 0 to 50 in 0.25 steps; ties removed
                            proportionally), error among the rest with Wilson 95% interval, the exact 95% band of removing
                            the same number AT RANDOM (hypergeometric), the share of all errors removed, the metric threshold
                            at the cut and the share of ALL COHORT calls beyond that threshold. `scale` is `pooled` (the
                            metric as is) or `within_gene` (percentile within the gene, over all cohort calls).
5. stat_model_coef.tsv      standardised coefficients of the multi-metric logistic models (all records, for reading)
6. stat_grid_<a>__<b>.tsv   error rate on a 40 x 40 grid of cohort percentile ranks of two metrics

Multi-metric models (L2 logistic, 5-fold CV grouped by person, out-of-fold risk is the "metric"):
  model:immuannot   Immuannot numbers only          model:depth     read depth only
  model:all         depth + Immuannot + assembly    (+ gene indicators in every model)

VM:  pixi run -e spechla -- python3 scripts/hla_popgen/71b_call_quality_stats.py
Tests: scripts/hla_popgen/tests/test_71b_call_quality_stats.py
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _s05_common as C  # noqa: E402


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


Q = sys.modules.get("m71") or _load("m71", "71_call_quality_table.py")
KEYS, SIGN, FAMILY = Q.KEYS, Q.SIGN, Q.FAMILY
INT = os.path.expanduser("~/s06/internal/71/")
OUT = os.path.expanduser("~/s06/results/71/")
PRUNE_GRID = np.round(np.arange(0, 50.0001, 0.25), 2)
VIEWS = [("trio_child", "trio", "c"), ("trio_parent", "trio", "p"), ("trio_worst", "trio", "w"), ("pair", "pair", "w"), ("twin", "twin", "w")]
MODEL_SETS = {
    "model:immuannot": [k for k in KEYS if Q.FAMILY[k] == "immuannot"],
    "model:depth": [k for k in KEYS if Q.FAMILY[k] == "depth"],
    "model:all": [k for k in KEYS if Q.FAMILY[k] in ("immuannot", "depth", "assembly")],
}
LOGY = {"template_distance", "template_distance_per_kb", "cds_distance", "n_aa_changes", "n_tied"}


def log(msg):
    C.log("[71b] " + msg)


# ---------------------------------------------------------------------------
# views
# ---------------------------------------------------------------------------
def make_view(rec, prefix):
    """Records -> DataFrame with person, gene, plat, bad and o_<metric> (higher = worse) for the chosen prefix."""
    v = rec[["person", "gene", "bad"]].copy()
    v["plat"] = rec["plat"] if "plat" in rec else "NA"
    for k in KEYS:
        col = f"{prefix}_{k}"
        v["o_" + k] = rec[col] if col in rec else np.nan
    return v.reset_index(drop=True)


def raw_value(k, o):
    """Oriented value back to the metric's own units."""
    return SIGN[k] * np.asarray(o, float)


# ---------------------------------------------------------------------------
# primitives
# ---------------------------------------------------------------------------
def auc_score(v, y):
    """P(error call has the higher value than a correct call) + 0.5 ties. v oriented (higher = worse)."""
    v, y = np.asarray(v, float), np.asarray(y, int)
    ok = np.isfinite(v)
    v, y = v[ok], y[ok]
    n1, n0 = int(y.sum()), int((1 - y).sum())
    if n1 == 0 or n0 == 0:
        return np.nan
    r = stats.rankdata(v)
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def cluster_boot_auc(v, y, groups, B=300, seed=71):
    rng = np.random.default_rng(seed)
    v, y = np.asarray(v, float), np.asarray(y, int)
    g = pd.factorize(np.asarray(groups))[0]
    order = np.argsort(g, kind="stable")
    v, y, g = v[order], y[order], g[order]
    starts = np.r_[0, np.flatnonzero(np.diff(g)) + 1]
    ends = np.r_[starts[1:], len(g)]
    ng = len(starts)
    out = []
    for _ in range(B):
        pick = rng.integers(0, ng, ng)
        idx = np.concatenate([np.arange(starts[i], ends[i]) for i in pick])
        a = auc_score(v[idx], y[idx])
        if np.isfinite(a):
            out.append(a)
    if len(out) < 20:
        return (np.nan, np.nan)
    return (float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5)))


def auc_within_gene(v, y, genes):
    v, y, genes = np.asarray(v, float), np.asarray(y, int), np.asarray(genes)
    num = den = 0.0
    for g in np.unique(genes):
        m = genes == g
        a = auc_score(v[m], y[m])
        if np.isfinite(a):
            w = float(y[m & np.isfinite(v)].sum() * (1 - y[m & np.isfinite(v)]).sum())
            num += a * w
            den += w
    return num / den if den else np.nan


def logit_cluster(X, y, groups, ridge=1e-4, iters=50):
    """Logistic regression by IRLS with a person-clustered sandwich covariance. X includes the intercept."""
    X, y = np.asarray(X, float), np.asarray(y, float)
    b = np.zeros(X.shape[1])
    pen = ridge * np.eye(X.shape[1])
    pen[0, 0] = 0
    for _ in range(iters):
        eta = X @ b
        p = 1 / (1 + np.exp(-np.clip(eta, -30, 30)))
        W = np.clip(p * (1 - p), 1e-9, None)
        H = X.T @ (X * W[:, None]) + pen
        step = np.linalg.solve(H, X.T @ (y - p) - pen @ b)
        b = b + step
        if np.max(np.abs(step)) < 1e-8:
            break
    p = 1 / (1 + np.exp(-np.clip(X @ b, -30, 30)))
    W = np.clip(p * (1 - p), 1e-9, None)
    bread = np.linalg.inv(X.T @ (X * W[:, None]) + pen)
    g = pd.factorize(np.asarray(groups))[0]
    scores = (X * (y - p)[:, None])
    meat = np.zeros((X.shape[1],) * 2)
    for gi in np.unique(g):
        s = scores[g == gi].sum(axis=0)
        meat += np.outer(s, s)
    cov = bread @ meat @ bread
    return b, np.sqrt(np.clip(np.diag(cov), 0, None))


def adjusted_or(v, y, genes, plat, groups, min_events=5):
    """OR per SD of the oriented metric, adjusted for gene and platform dummies (levels with >= min_events errors)."""
    v = np.asarray(v, float)
    ok = np.isfinite(v)
    if ok.sum() < 50 or np.asarray(y)[ok].sum() < 5:
        return (np.nan,) * 4
    z = (v[ok] - v[ok].mean()) / (v[ok].std() or 1.0)
    cols = [np.ones(ok.sum()), z]
    ys = np.asarray(y)[ok]
    for name, arr in (("gene", np.asarray(genes)[ok]), ("plat", np.asarray(plat)[ok])):
        lev = pd.Series(ys).groupby(arr).sum()
        for l in sorted(lev[lev >= min_events].index)[1:]:
            cols.append((arr == l).astype(float))
    try:
        b, se = logit_cluster(np.column_stack(cols), ys, np.asarray(groups)[ok])
    except np.linalg.LinAlgError:
        return (np.nan,) * 4
    zv = b[1] / se[1] if se[1] > 0 else np.nan
    p = float(2 * stats.norm.sf(abs(zv))) if np.isfinite(zv) else np.nan
    return float(np.exp(b[1])), float(np.exp(b[1] - 1.96 * se[1])), float(np.exp(b[1] + 1.96 * se[1])), p


def wilson_arr(k, n):
    k, n = np.asarray(k, float), np.asarray(n, float)
    z = 1.96
    with np.errstate(divide="ignore", invalid="ignore"):
        p = k / n
        den = 1 + z * z / n
        c = (p + z * z / (2 * n)) / den
        h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, np.clip(c - h, 0, 1), np.clip(c + h, 0, 1)


# ---------------------------------------------------------------------------
# rolling curve
# ---------------------------------------------------------------------------
def rolling_curve(o, y, k, max_discrete=25):
    """Error rate along the metric. Discrete metrics: one row per value. Continuous: rolling windows of n/25 records,
    stepping n/300, ordered from best to worst."""
    o, y = np.asarray(o, float), np.asarray(y, int)
    ok = np.isfinite(o)
    o, y = o[ok], y[ok]
    n = len(o)
    if n < 50:
        return pd.DataFrame()
    u = np.unique(o)
    rows = []
    if len(u) <= max_discrete:
        for val in u:
            m = o == val
            rows.append((raw_value(k, val), raw_value(k, val), raw_value(k, val), int(m.sum()), int(y[m].sum())))
    else:
        order = np.argsort(o, kind="stable")
        so, sy = o[order], y[order]
        w = max(60, n // 25)
        step = max(1, n // 300)
        cs = np.r_[0, np.cumsum(sy)]
        for s in range(0, n - w + 1, step):
            e = s + w
            rows.append((raw_value(k, np.median(so[s:e])), raw_value(k, so[s]), raw_value(k, so[e - 1]), w,
                         int(cs[e] - cs[s])))
    d = pd.DataFrame(rows, columns=["x_median", "x_a", "x_b", "n", "k"])
    p, lo, hi = wilson_arr(d["k"], d["n"])
    d["error_pct"], d["lo"], d["hi"] = 100 * p, 100 * lo, 100 * hi
    return d


# ---------------------------------------------------------------------------
# pruning
# ---------------------------------------------------------------------------
def prune_curve(o, y, grid=PRUNE_GRID, cohort_o=None, k_name=None, band_draws=True):
    """Remove the worst x% of records (highest oriented value). Tied values are removed proportionally, so a binary
    metric is not cut at an arbitrary row. -> DataFrame."""
    o, y = np.asarray(o, float), np.asarray(y, int)
    ok = np.isfinite(o)
    o, y = o[ok], y[ok]
    n, K = len(o), int(y.sum())
    if n < 50 or K == 0:
        return pd.DataFrame()
    vals, inv = np.unique(-o, return_inverse=True)          # ascending in -o = worst first
    cnt = np.bincount(inv)
    err = np.bincount(inv, weights=y)
    ccnt, cerr = np.cumsum(cnt), np.cumsum(err)
    rows = []
    cs = None
    if cohort_o is not None:
        cs = np.sort(np.asarray(cohort_o, float)[np.isfinite(cohort_o)])
    for x in grid:
        kk = n * x / 100.0
        g = int(np.searchsorted(ccnt, kk, side="right"))     # whole tie groups removed
        removed_n = ccnt[g - 1] if g > 0 else 0.0
        removed_e = cerr[g - 1] if g > 0 else 0.0
        if g < len(cnt) and kk > removed_n:
            f = (kk - removed_n) / cnt[g]
            removed_n, removed_e = removed_n + f * cnt[g], removed_e + f * err[g]
        keep_n, keep_e = n - removed_n, K - removed_e
        rate = keep_e / keep_n if keep_n > 0 else np.nan
        p, lo, hi = wilson_arr([keep_e], [keep_n])
        row = dict(prune_pct=float(x), n_removed=float(removed_n), n_kept=float(keep_n), errors_kept=float(keep_e),
                   error_pct_kept=100 * rate, lo=100 * float(lo[0]), hi=100 * float(hi[0]),
                   captured_pct=100 * removed_e / K, thr_oriented=float(-vals[min(g, len(vals) - 1)]))
        if band_draws and 0 < kk < n:
            m = int(round(kk))
            rl, rh = stats.hypergeom.ppf([0.025, 0.975], n, K, m)      # errors among the removed at random
            row["random_lo"] = 100 * (K - rh) / (n - m)
            row["random_hi"] = 100 * (K - rl) / (n - m)
        else:
            row["random_lo"] = row["random_hi"] = 100 * K / n
        row["random_mean"] = 100 * K / n
        if cs is not None:
            row["cohort_pct_removed"] = 100 * (len(cs) - np.searchsorted(cs, row["thr_oriented"], side="left")) / len(cs) \
                if x > 0 else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def within_gene_percentile(o, genes, cohort_ct, k):
    """Percentile (0-100, higher = worse) of each record's oriented value among ALL cohort calls of its gene."""
    out = np.full(len(o), np.nan)
    o = np.asarray(o, float)
    genes = np.asarray(genes)
    col = "o_" + k
    for g in np.unique(genes):
        ref = np.sort(cohort_ct.loc[cohort_ct["gene"] == g, col].dropna().to_numpy(float))
        m = genes == g
        if len(ref) < 20:
            continue
        lo = np.searchsorted(ref, o[m], side="left")
        hi = np.searchsorted(ref, o[m], side="right")
        out[m] = 100.0 * (lo + hi) / 2.0 / len(ref)
    return out


# ---------------------------------------------------------------------------
# multi-metric models
# ---------------------------------------------------------------------------
def design_matrix(view, keys, cohort_median=None):
    """Oriented metrics (log1p for distance counts), median-imputed, standardised, + gene dummies. -> X, names."""
    cols, names = [], []
    for k in keys:
        v = view["o_" + k].to_numpy(float)
        if k in LOGY:
            v = np.log1p(np.clip(v, 0, None))
        med = np.nanmedian(v) if np.isfinite(v).any() else 0.0
        v = np.where(np.isfinite(v), v, med)
        sd = v.std()
        if sd == 0:
            continue
        cols.append((v - v.mean()) / sd)
        names.append(k)
    for g in sorted(view["gene"].unique()):
        d = (view["gene"] == g).to_numpy(float)
        if 0 < d.sum() < len(d):
            cols.append(d)
            names.append("gene:" + g)
    return np.column_stack(cols), names


def oof_risk(X, y, groups, folds=5, C_=0.3, seed=71):
    """Out-of-fold predicted risk with folds grouped by person. sklearn when available, else the IRLS fit above."""
    y = np.asarray(y, int)
    gid = pd.factorize(np.asarray(groups))[0]
    rng = np.random.default_rng(seed)
    fold_of_group = rng.integers(0, folds, gid.max() + 1)
    fold = fold_of_group[gid]
    risk = np.full(len(y), np.nan)
    try:
        from sklearn.linear_model import LogisticRegression
        for f in range(folds):
            tr, te = fold != f, fold == f
            if y[tr].sum() < 3 or te.sum() == 0:
                continue
            clf = LogisticRegression(C=C_, max_iter=2000).fit(X[tr], y[tr])
            risk[te] = clf.predict_proba(X[te])[:, 1]
    except ImportError:
        for f in range(folds):
            tr, te = fold != f, fold == f
            if y[tr].sum() < 3 or te.sum() == 0:
                continue
            Xt = np.column_stack([np.ones(tr.sum()), X[tr]])
            b, _ = logit_cluster(Xt, y[tr], np.arange(tr.sum()), ridge=1.0 / C_)
            risk[te] = 1 / (1 + np.exp(-(np.column_stack([np.ones(te.sum()), X[te]]) @ b)))
    return risk


def model_coefficients(X, names, y):
    """Standardised coefficients from an L2 logistic fit on all records, with person-free robust SE omitted (reading only)."""
    Xt = np.column_stack([np.ones(len(X)), X])
    b, se = logit_cluster(Xt, y, np.arange(len(y)), ridge=1.0 / 0.3)
    return pd.DataFrame({"term": ["intercept"] + names, "coef": b, "se_naive": se, "odds_ratio": np.exp(b)})


# ---------------------------------------------------------------------------
# 2-D grids on cohort percentile ranks
# ---------------------------------------------------------------------------
def percentile_rank(values, ref):
    ref = np.sort(np.asarray(ref, float)[np.isfinite(ref)])
    v = np.asarray(values, float)
    lo, hi = np.searchsorted(ref, v, side="left"), np.searchsorted(ref, v, side="right")
    return np.where(np.isfinite(v), 100.0 * (lo + hi) / 2.0 / max(len(ref), 1), np.nan)


def grid2d(view, ka, kb, cohort_ct, bins=40):
    a = percentile_rank(view["o_" + ka], cohort_ct["o_" + ka])
    b = percentile_rank(view["o_" + kb], cohort_ct["o_" + kb])
    ok = np.isfinite(a) & np.isfinite(b)
    edges = np.linspace(0, 100, bins + 1)
    ia = np.clip(np.digitize(a[ok], edges) - 1, 0, bins - 1)
    ib = np.clip(np.digitize(b[ok], edges) - 1, 0, bins - 1)
    y = view["bad"].to_numpy(int)[ok]
    n = np.zeros((bins, bins))
    k = np.zeros((bins, bins))
    np.add.at(n, (ia, ib), 1)
    np.add.at(k, (ia, ib), y)
    rows = [dict(a_lo=edges[i], b_lo=edges[j], n=int(n[i, j]), k=int(k[i, j]))
            for i in range(bins) for j in range(bins) if n[i, j] > 0]
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
def analyse(views, cohort_ct, out_dir, do_boot=True):
    w = lambda df, name: C.write_tsv(df, os.path.join(out_dir, name))  # noqa: E731
    base, curves, aucs, prunes, coefs = [], [], [], [], []
    for vname, v in views.items():
        if not len(v):
            continue
        y = v["bad"].to_numpy(int)
        for g, d in [("ALL", v)] + list(v.groupby("gene")):
            base.append(dict(view=vname, gene=g, n_records=len(d), n_errors=int(d["bad"].sum()),
                             error_pct=100 * d["bad"].mean()))
        scores = {k: v["o_" + k].to_numpy(float) for k in KEYS}
        # multi-metric out-of-fold risks
        for mname, ks in MODEL_SETS.items():
            if y.sum() < 10:
                continue
            X, names = design_matrix(v, ks)
            scores[mname] = oof_risk(X, y, v["person"])
            cf = model_coefficients(X, names, y)
            cf.insert(0, "model", mname)
            cf.insert(0, "view", vname)
            coefs.append(cf)
        for k, o in scores.items():
            is_model = k.startswith("model:")
            fam = "model" if is_model else FAMILY[k]
            lo, hi = cluster_boot_auc(o, y, v["person"]) if do_boot and np.isfinite(o).sum() > 50 else (np.nan, np.nan)
            a = auc_score(o, y)
            okm = np.isfinite(o)
            rho, pr = stats.spearmanr(o[okm], y[okm]) if okm.sum() > 30 and y[okm].std() > 0 else (np.nan, np.nan)
            orr = adjusted_or(o, y, v["gene"], v["plat"], v["person"]) if not is_model else (np.nan,) * 4
            aucs.append(dict(view=vname, metric=k, family=fam, n=int(okm.sum()), n_errors=int(y[okm].sum()), auc=a,
                             auc_lo=lo, auc_hi=hi, auc_within_gene=auc_within_gene(o, y, v["gene"]), spearman=rho,
                             spearman_p=pr, or_per_sd=orr[0], or_lo=orr[1], or_hi=orr[2], or_p=orr[3]))
            if not is_model:
                c = rolling_curve(o, y, k)
                if len(c):
                    c.insert(0, "metric", k)
                    c.insert(0, "view", vname)
                    curves.append(c)
            cohort = cohort_ct["o_" + k].to_numpy(float) if (not is_model and ("o_" + k) in cohort_ct) else None
            for scale in ("pooled", "within_gene"):
                if scale == "within_gene":
                    if is_model or cohort is None:
                        continue
                    oo = within_gene_percentile(o, v["gene"], cohort_ct, k)
                    pc = prune_curve(oo, y, cohort_o=None)
                else:
                    pc = prune_curve(o, y, cohort_o=cohort)
                if len(pc):
                    pc.insert(0, "scale", scale)
                    pc.insert(0, "metric", k)
                    pc.insert(0, "view", vname)
                    prunes.append(pc)
    w(pd.DataFrame(base), "stat_error_baseline.tsv")
    w(pd.concat(curves, ignore_index=True) if curves else pd.DataFrame(), "stat_curves.tsv")
    w(pd.DataFrame(aucs), "stat_auc.tsv")
    w(pd.concat(prunes, ignore_index=True) if prunes else pd.DataFrame(), "stat_prune.tsv")
    w(pd.concat(coefs, ignore_index=True) if coefs else pd.DataFrame(), "stat_model_coef.tsv")
    for vname, ka, kb in (("pair", "rel_cds", "template_distance"), ("pair", "rel_ars", "cds_distance"),
                          ("pair", "rel_cds", "rel_mhc")):
        if vname in views and len(views[vname]):
            w(grid2d(views[vname], ka, kb, cohort_ct), f"stat_grid_{ka}__{kb}.tsv")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--internal-dir", default=INT)
    ap.add_argument("--out-dir", default=OUT)
    ap.add_argument("--no-boot", action="store_true", help="skip the cluster bootstrap CIs (quick look)")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    ct = Q.read_int(a.internal_dir + "call_table.INTERNAL.tsv", dtype={"person_id": str})
    recs = {}
    for name in ("trio", "pair", "twin"):
        p = a.internal_dir + f"records_{name}.INTERNAL.tsv"
        recs[name] = pd.read_csv(p, sep="\t", dtype={"person": str, "other": str}) if os.path.exists(p) else pd.DataFrame()
    views = {vn: make_view(recs[src], pre) for vn, src, pre in VIEWS if len(recs[src])}
    analyse(views, ct, a.out_dir, do_boot=not a.no_boot)
    log("DONE 71b")
    return 0


if __name__ == "__main__":
    sys.exit(main())
