#!/usr/bin/env python3
"""S06 76b: pruning evidence per outcome, with transfer and robustness (VM; reads 76's INTERNAL unit tables, writes aggregate curves).

Question: if the worst x% of calls by a metric are removed, how many problems go with them, at what cost, and does a rule chosen on one set of
families / one outcome hold on another? Nothing here picks "the" answer; every number needed to choose is tabulated.

Views (each a table of units with `bad`, `family`/group and oriented o_<metric> columns, higher = worse):
  A_classical8, A_all     person x gene units, bad = incomplete genotype (stream A)       [metrics: one random call, 71e]
  B_child, B_worst        trio calls eligible for labelling, bad = `neither` class         [metrics: child's call / worse of child and parent]
  B_child_all             same with all trios (not only the primary set)
  combined_child          child person x gene, bad = missing OR wrong; missing / wrong curves kept apart
  pair_complete, pair_all pair records, bad = no shared 2-field allele                    [metrics: worse of the two people]

Per view x metric (prune_curves.tsv, 0.25% steps up to --max-prune):
  n_pruned, errors_removed, captured_pct (errors removed / all errors), residual_error_pct, base_error_pct, enrichment (captured fraction / pruned fraction),
  nnp = calls pruned per error removed (cumulative), marginal_nnp = calls pruned per error removed over the last --nnp-window percent,
  random band for errors removed (exact hypergeometric 95%), p_tail (P[random pruning removes >= as many errors]), person/family-cluster bootstrap
  95% intervals for captured_pct and residual_error_pct.
"Optimal" is NOT chosen by the script; for convenience prune_summary.tsv gives x_opt = the largest x before the marginal NNP exceeds --cost-ratio
(default 25 calls per error removed), with its bootstrap spread, next to the full curves.
Robustness: cv_by_metric.tsv / cv_selected.tsv (family-grouped K-fold repeated: threshold learned on training families, applied to held-out families),
transfer.tsv (threshold learned on one view applied to another), power.tsv (events, minimal detectable AUC), threshold_stability.tsv.

VM:  pixi run -e spechla -- python3 scripts/hla_popgen/76b_pruning_evidence.py --internal-dir ~/s06/internal/76/ --out-dir ~/s06/results/76b/
Tests: scripts/hla_popgen/tests/test_76b_pruning_evidence.py
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


S = sys.modules.get("m71b") or _load("m71b", "71b_call_quality_stats.py")
Q = S.Q
INT = os.path.expanduser("~/s06/internal/76/")
OUT = os.path.expanduser("~/s06/results/76b/")
FIXED_X = (1.0, 2.5, 5.0, 10.0)           # pre-specified tail tests / reporting points
TRANSFERS = [("A_classical8", "B_child"), ("B_child", "A_classical8"), ("A_all", "B_child"), ("B_child", "pair_complete"),
             ("pair_complete", "B_child"), ("B_worst", "pair_complete"), ("pair_complete", "B_worst"), ("A_all", "pair_all"),
             ("pair_all", "A_all")]


def log(m):
    C.log("[76b] " + m)


# ---------------------------------------------------------------------------
# curve primitives
# ---------------------------------------------------------------------------
def removed_curve(o, Y, grid):
    """Remove the worst x% (highest o) for each x in grid. Tied values are removed proportionally (a binary metric is not cut at an arbitrary row).
    o: (n,), Y: (n, m) label columns. -> n_removed (len(grid),), e_removed (len(grid), m)."""
    o = np.asarray(o, float)
    Y = np.asarray(Y, float)
    if Y.ndim == 1:
        Y = Y[:, None]
    ok = np.isfinite(o)
    o, Y = o[ok], Y[ok]
    n = len(o)
    vals, inv = np.unique(-o, return_inverse=True)
    cnt = np.bincount(inv).astype(float)
    err = np.column_stack([np.bincount(inv, weights=Y[:, j]) for j in range(Y.shape[1])])
    ccnt, cerr = np.cumsum(cnt), np.cumsum(err, axis=0)
    nr = np.zeros(len(grid))
    er = np.zeros((len(grid), Y.shape[1]))
    for i, x in enumerate(grid):
        kk = n * x / 100.0
        g = int(np.searchsorted(ccnt, kk, side="right"))
        rn = ccnt[g - 1] if g > 0 else 0.0
        re = cerr[g - 1].copy() if g > 0 else np.zeros(Y.shape[1])
        if g < len(cnt) and kk > rn:
            f = (kk - rn) / cnt[g]
            rn, re = rn + f * cnt[g], re + f * err[g]
        nr[i], er[i] = rn, re
    return nr, er, n, Y.sum(axis=0)


def x_opt_from(nr, er_main, grid, cost, window):
    """Largest x reached before the marginal NNP (calls pruned per error removed over the last `window` percent) first exceeds `cost`."""
    step = grid[1] - grid[0]
    w = max(int(round(window / step)), 1)
    best = 0.0
    for i in range(w, len(grid)):
        de = er_main[i] - er_main[i - w]
        dn = nr[i] - nr[i - w]
        mn = dn / de if de > 0 else np.inf
        if mn > cost:
            return float(grid[i - w])
        best = float(grid[i])
    return best


def cluster_index(groups):
    g = pd.factorize(np.asarray(groups))[0]
    order = np.argsort(g, kind="stable")
    gs = g[order]
    starts = np.r_[0, np.flatnonzero(np.diff(gs)) + 1]
    ends = np.r_[starts[1:], len(gs)]
    return order, starts, ends


def curves_for(o, Y, groups, grid, cost, window, boot, seed):
    """Point curve plus cluster-bootstrap draws of captured_pct / residual_error_pct / x_opt."""
    nr, er, n, K = removed_curve(o, Y, grid)
    out = dict(nr=nr, er=er, n=n, K=K)
    if boot and K[0] > 0:
        rng = np.random.default_rng(seed)
        order, starts, ends = cluster_index(groups)
        o_, Y_ = np.asarray(o, float)[order], np.asarray(Y, float).reshape(len(o), -1)[order]
        ng = len(starts)
        cap, res, xo = [], [], []
        for _ in range(boot):
            pick = rng.integers(0, ng, ng)
            idx = np.concatenate([np.arange(starts[i], ends[i]) for i in pick])
            nr_b, er_b, n_b, K_b = removed_curve(o_[idx], Y_[idx], grid)
            if K_b[0] == 0:
                continue
            cap.append(100 * er_b[:, 0] / K_b[0])
            res.append(100 * (K_b[0] - er_b[:, 0]) / np.maximum(n_b - nr_b, 1))
            xo.append(x_opt_from(nr_b, er_b[:, 0], grid, cost, window))
        if len(cap) >= 20:
            out["cap_ci"] = (np.percentile(cap, 2.5, axis=0), np.percentile(cap, 97.5, axis=0))
            out["res_ci"] = (np.percentile(res, 2.5, axis=0), np.percentile(res, 97.5, axis=0))
            out["xopt_boot"] = np.array(xo)
    return out


def curve_rows(view, metric, cv, grid, window, label_names):
    nr, er, n, K = cv["nr"], cv["er"], cv["n"], cv["K"]
    rows = []
    step = grid[1] - grid[0]
    w = max(int(round(window / step)), 1)
    for i, x in enumerate(grid):
        e = er[i, 0]
        keep_n = n - nr[i]
        row = dict(view=view, metric=metric, prune_pct=float(x), n_units=int(n), n_errors=int(K[0]), n_pruned=float(nr[i]), errors_removed=float(e),
                   captured_pct=100 * e / K[0], residual_error_pct=100 * (K[0] - e) / keep_n if keep_n > 0 else np.nan,
                   base_error_pct=100 * K[0] / n,
                   enrichment=(e / K[0]) / (nr[i] / n) if nr[i] > 0 else np.nan, nnp=nr[i] / e if e > 0 else np.inf)
        if i >= w:
            de, dn = e - er[i - w, 0], nr[i] - nr[i - w]
            row["marginal_nnp"] = dn / de if de > 0 else np.inf
        else:
            row["marginal_nnp"] = np.nan
        m = int(round(nr[i]))
        if 0 < m < n and K[0] > 0:
            row["rand_removed_lo"] = float(stats.hypergeom.ppf(0.025, n, int(K[0]), m))
            row["rand_removed_hi"] = float(stats.hypergeom.ppf(0.975, n, int(K[0]), m))
            row["p_tail"] = float(stats.hypergeom.sf(int(round(e)) - 1, n, int(K[0]), m))
        else:
            row["rand_removed_lo"] = row["rand_removed_hi"] = row["p_tail"] = np.nan
        if "cap_ci" in cv:
            row["captured_lo"], row["captured_hi"] = float(cv["cap_ci"][0][i]), float(cv["cap_ci"][1][i])
            row["residual_lo"], row["residual_hi"] = float(cv["res_ci"][0][i]), float(cv["res_ci"][1][i])
        for j, nm in enumerate(label_names[1:], start=1):
            row["captured_pct_" + nm] = 100 * er[i, j] / K[j] if K[j] > 0 else np.nan
            row["errors_removed_" + nm] = float(er[i, j])
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# views
# ---------------------------------------------------------------------------
def read_view(path, rename=None):
    d = Q.read_int(path)
    return d


def make_views(internal_dir):
    """-> {name: DataFrame(group, bad, [missing, wrong], o_<metric>...)}"""
    v = {}

    def get(name):
        p = os.path.join(internal_dir, name)
        return S.Q.read_int(p, dtype={"person": str, "other": str}) if os.path.exists(p) else None

    a = get("units_A.INTERNAL.tsv")
    if a is not None:
        a["group"] = a["person"]
        v["A_all"] = a
        v["A_classical8"] = a[a["scope_classical8"].astype(str).str.lower() == "true"].reset_index(drop=True)
    b = get("units_B.INTERNAL.tsv")
    ball = get("units_B_all.INTERNAL.tsv")
    for nm, d in (("B", b), ("B_all", ball)):
        if d is None:
            continue
        for pre, tag in (("c", "child"), ("w", "worst")):
            x = d[["person", "gene", "bad", "family"]].copy()
            x["group"] = d["family"]
            for k in Q.KEYS:
                x["o_" + k] = d.get(pre + "_" + k, np.nan)
            v["%s_%s" % (nm, tag)] = x
    c = get("units_combined.INTERNAL.tsv")
    if c is not None:
        c["group"] = c["family"]
        v["combined_child"] = c
    for nm in ("all", "complete"):
        p = get("units_pair_%s.INTERNAL.tsv" % nm)
        if p is not None:
            x = p[["person", "gene", "bad", "family"]].copy()
            x["group"] = p["family"]
            for k in Q.KEYS:
                x["o_" + k] = p.get("w_" + k, np.nan)
            v["pair_" + nm] = x
    return v


def usable_metrics(d, min_finite=50):
    return [k for k in Q.KEYS if ("o_" + k) in d and np.isfinite(d["o_" + k].to_numpy(float)).sum() >= min_finite and k != "n_called"]


def power_row(view, n, K):
    n0 = n - K
    if K < 1 or n0 < 1:
        return dict(view=view, n_units=n, n_errors=K, auc_min_detectable=np.nan)
    se = np.sqrt((n + 1) / (12.0 * K * n0))
    return dict(view=view, n_units=n, n_errors=K, auc_se_null=se, auc_min_detectable=0.5 + (1.645 + 0.842) * se,
                note="one-sided alpha 0.05, power 0.8, independent units (clustering makes real power lower)")


def threshold_at(o, x):
    """Oriented value of the boundary call when the worst x% are pruned (prune rows with o >= thr)."""
    v = np.sort(np.asarray(o, float)[np.isfinite(o)])[::-1]
    if not len(v) or x <= 0:
        return np.inf
    k = int(np.ceil(len(v) * x / 100.0)) - 1
    return float(v[min(max(k, 0), len(v) - 1)])


def apply_threshold(o, y, thr):
    o = np.asarray(o, float)
    ok = np.isfinite(o)
    pr = ok & (o >= thr)
    n, K = int(ok.sum()), int(np.asarray(y)[ok].sum())
    e = int(np.asarray(y)[pr].sum())
    m = int(pr.sum())
    p = float(stats.hypergeom.sf(e - 1, n, K, m)) if 0 < m < n and K > 0 else np.nan
    return dict(n=n, n_errors=K, n_pruned=m, pruned_pct=100 * m / n if n else np.nan, errors_removed=e,
                captured_pct=100 * e / K if K else np.nan, enrichment=(e / K) / (m / n) if K and m else np.nan, p_tail=p)


def learn_x(o, y, grid, cost, window):
    nr, er, n, K = removed_curve(o, np.asarray(y, float), grid)
    if K[0] < 5:
        return np.nan
    return x_opt_from(nr, er[:, 0], grid, cost, window)


def cv_run(d, metrics, grid, cost, window, folds, reps, seed, select_by="auc"):
    """Family-grouped K-fold repeated `reps` times. In every training split: per metric learn x_opt and its value threshold; per split also pick the
    best metric by training AUC. On the held-out families: prune rows at/above the threshold and score."""
    groups = d["group"].to_numpy()
    ug = np.unique(groups)
    y = d["bad"].to_numpy(int)
    O = {k: d["o_" + k].to_numpy(float) for k in metrics}
    per_metric, sel_rows, thr_rows = [], [], []
    for r in range(reps):
        rng = np.random.default_rng(seed + r)
        perm = rng.permutation(ug)
        fold_of = {g: i % folds for i, g in enumerate(perm)}
        f = np.array([fold_of[g] for g in groups])
        acc = {k: np.zeros(4) for k in metrics}     # test n, test errors, pruned, errors removed
        sacc = np.zeros(4)
        chosen = []
        for fo in range(folds):
            tr, te = f != fo, f == fo
            if y[tr].sum() < 5 or y[te].sum() < 1:
                continue
            aucs = {}
            for k in metrics:
                aucs[k] = S.auc_score(O[k][tr], y[tr])
                x = learn_x(O[k][tr], y[tr], grid, cost, window)
                if not np.isfinite(x) or x <= 0:
                    thr = np.inf
                else:
                    thr = threshold_at(O[k][tr], x)
                res = apply_threshold(O[k][te], y[te], thr)
                acc[k] += [res["n"], res["n_errors"], res["n_pruned"], res["errors_removed"]]
                thr_rows.append(dict(rep=r, fold=fo, metric=k, x_opt=x, thr=thr))
            best = max((k for k in metrics if np.isfinite(aucs[k])), key=lambda k: aucs[k], default=None)
            if best is not None:
                thr_b = [t for t in thr_rows[-len(metrics):] if t["metric"] == best][0]["thr"]
                res = apply_threshold(O[best][te], y[te], thr_b)
                sacc += [res["n"], res["n_errors"], res["n_pruned"], res["errors_removed"]]
                chosen.append(best)
        for k in metrics:
            n, K, m, e = acc[k]
            per_metric.append(dict(rep=r, metric=k, test_units=n, test_errors=K, pruned=m, errors_removed=e))
        sel_rows.append(dict(rep=r, test_units=sacc[0], test_errors=sacc[1], pruned=sacc[2], errors_removed=sacc[3], chosen=";".join(chosen).replace(";", ",")))
    pm = pd.DataFrame(per_metric)
    by = []
    for k, g in pm.groupby("metric"):
        g = g[g["test_errors"] > 0]
        if g.empty:
            continue
        pr = g["pruned"] / g["test_units"]
        cp = g["errors_removed"] / g["test_errors"]
        enr = cp / pr.replace(0, np.nan)
        by.append(dict(metric=k, reps=len(g), pruned_pct=100 * pr.mean(), captured_pct=100 * cp.mean(), captured_pct_lo=100 * cp.quantile(0.025),
                       captured_pct_hi=100 * cp.quantile(0.975), enrichment=enr.mean(), enrichment_lo=enr.quantile(0.025),
                       enrichment_hi=enr.quantile(0.975), share_reps_enrichment_gt1=float((enr > 1).mean())))
    sel = pd.DataFrame(sel_rows)
    sel = sel[sel["test_errors"] > 0]
    sel_out = pd.DataFrame()
    if len(sel):
        pr = sel["pruned"] / sel["test_units"]
        cp = sel["errors_removed"] / sel["test_errors"]
        enr = cp / pr.replace(0, np.nan)
        freq = pd.Series(",".join(sel["chosen"]).split(",")).value_counts()
        sel_out = pd.DataFrame([dict(reps=len(sel), pruned_pct=100 * pr.mean(), captured_pct=100 * cp.mean(), captured_pct_lo=100 * cp.quantile(0.025),
                                     captured_pct_hi=100 * cp.quantile(0.975), enrichment=enr.mean(), enrichment_lo=enr.quantile(0.025),
                                     enrichment_hi=enr.quantile(0.975), top_chosen_metrics=",".join("%s:%d" % (a, b) for a, b in freq.head(5).items()))])
    th = pd.DataFrame(thr_rows)
    stab = []
    for k, g in th.groupby("metric"):
        t = g["thr"].replace(np.inf, np.nan).dropna()
        x = g["x_opt"].dropna()
        stab.append(dict(metric=k, folds=len(g), folds_with_prune=int(len(t)), x_opt_median=x.median() if len(x) else np.nan,
                         x_opt_iqr_lo=x.quantile(0.25) if len(x) else np.nan, x_opt_iqr_hi=x.quantile(0.75) if len(x) else np.nan,
                         thr_median=t.median() if len(t) else np.nan, thr_iqr_lo=t.quantile(0.25) if len(t) else np.nan,
                         thr_iqr_hi=t.quantile(0.75) if len(t) else np.nan))
    return pd.DataFrame(by), sel_out, pd.DataFrame(stab)


def run(internal_dir, out_dir, grid_max=20.0, step=0.25, cost=25.0, window=2.0, boot=200, folds=5, reps=20, seed=76, views=None):
    os.makedirs(out_dir, exist_ok=True)
    w = lambda df, name: C.write_tsv(df, os.path.join(out_dir, name))  # noqa: E731
    grid = np.round(np.arange(0, grid_max + 1e-9, step), 4)
    V = make_views(internal_dir)
    if views:
        V = {k: d for k, d in V.items() if k in views}
    curves, summ, aucs, power, cvm, cvs, stab = [], [], [], [], [], [], []
    for vn, d in V.items():
        d = d.reset_index(drop=True)
        y = d["bad"].to_numpy(int)
        power.append(power_row(vn, len(d), int(y.sum())))
        if y.sum() < 3:
            log("%s: %d events, skipped" % (vn, int(y.sum())))
            continue
        metrics = usable_metrics(d)
        labels = ["bad"] + [c for c in ("missing", "wrong") if c in d]
        Y = d[labels].to_numpy(float)
        for k in metrics:
            o = d["o_" + k].to_numpy(float)
            cv = curves_for(o, Y, d["group"].to_numpy(), grid, cost, window, boot, seed)
            curves += curve_rows(vn, k, cv, grid, window, labels)
            a = S.auc_score(o, y)
            lo, hi = S.cluster_boot_auc(o, y, d["group"].to_numpy(), B=boot) if np.isfinite(a) else (np.nan, np.nan)
            aucs.append(dict(view=vn, metric=k, n_units=int(np.isfinite(o).sum()), n_errors=int(y[np.isfinite(o)].sum()), auc=a, auc_lo=lo, auc_hi=hi))
            xo = x_opt_from(cv["nr"], cv["er"][:, 0], grid, cost, window)
            row = dict(view=vn, metric=k, n_units=cv["n"], n_errors=int(cv["K"][0]), auc=a, auc_lo=lo, auc_hi=hi, x_opt=xo)
            if "xopt_boot" in cv:
                xb = cv["xopt_boot"]
                row.update(x_opt_boot_median=float(np.median(xb)), x_opt_boot_q25=float(np.percentile(xb, 25)), x_opt_boot_q75=float(np.percentile(xb, 75)))
            for x in FIXED_X:
                i = int(np.argmin(np.abs(grid - x)))
                if x <= grid_max:
                    e = cv["er"][i, 0]
                    row["captured_pct_at_%g" % x] = 100 * e / cv["K"][0]
                    row["enrichment_at_%g" % x] = (e / cv["K"][0]) / (cv["nr"][i] / cv["n"]) if cv["nr"][i] > 0 else np.nan
                    m = int(round(cv["nr"][i]))
                    row["p_tail_at_%g" % x] = float(stats.hypergeom.sf(int(round(e)) - 1, cv["n"], int(cv["K"][0]), m)) if 0 < m < cv["n"] else np.nan
            i = int(np.argmin(np.abs(grid - xo)))
            row["captured_pct_at_xopt"] = 100 * cv["er"][i, 0] / cv["K"][0]
            row["residual_error_pct_at_xopt"] = 100 * (cv["K"][0] - cv["er"][i, 0]) / (cv["n"] - cv["nr"][i])
            summ.append(row)
        by, sel, st = cv_run(d, metrics, grid, cost, window, folds, reps, seed)
        for df, lst in ((by, cvm), (sel, cvs), (st, stab)):
            if len(df):
                df.insert(0, "view", vn)
                lst.append(df)
        log("%s: %d units, %d events, %d metrics" % (vn, len(d), int(y.sum()), len(metrics)))
    w(pd.DataFrame(curves), "prune_curves.tsv")
    w(pd.DataFrame(summ), "prune_summary.tsv")
    w(pd.DataFrame(aucs), "auc.tsv")
    w(pd.DataFrame(power), "power.tsv")
    if cvm:
        w(pd.concat(cvm, ignore_index=True), "cv_by_metric.tsv")
    if cvs:
        w(pd.concat(cvs, ignore_index=True), "cv_selected.tsv")
    if stab:
        w(pd.concat(stab, ignore_index=True), "threshold_stability.tsv")
    # transfer: threshold learned on the source view (x_opt rule and a fixed 5%), applied to the target view
    tr = []
    for src, dst in TRANSFERS:
        if src not in V or dst not in V:
            continue
        ds, dd = V[src], V[dst]
        for k in sorted(set(usable_metrics(ds)) & set(usable_metrics(dd))):
            os_, ys = ds["o_" + k].to_numpy(float), ds["bad"].to_numpy(int)
            if ys.sum() < 5 or dd["bad"].sum() < 3:
                continue
            for rule, x in (("x_opt", learn_x(os_, ys, grid, cost, window)), ("fixed_5pct", 5.0)):
                if not np.isfinite(x) or x <= 0:
                    continue
                thr = threshold_at(os_, x)
                r = apply_threshold(dd["o_" + k].to_numpy(float), dd["bad"].to_numpy(int), thr)
                tr.append(dict(source=src, target=dst, metric=k, rule=rule, source_x=x, thr_oriented=thr,
                               source_pruned_pct=apply_threshold(os_, ys, thr)["pruned_pct"], **{"target_" + a: b for a, b in r.items()}))
    if tr:
        w(pd.DataFrame(tr), "transfer.tsv")
    log("DONE 76b")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--internal-dir", default=INT)
    ap.add_argument("--out-dir", default=OUT)
    ap.add_argument("--max-prune", type=float, default=20.0)
    ap.add_argument("--step", type=float, default=0.25)
    ap.add_argument("--cost-ratio", type=float, default=25.0, help="calls one is willing to prune per error removed")
    ap.add_argument("--nnp-window", type=float, default=2.0, help="percent of units over which the marginal NNP is measured")
    ap.add_argument("--boot", type=int, default=200)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--views", nargs="*", default=None)
    a = ap.parse_args(argv)
    run(a.internal_dir, a.out_dir, a.max_prune, a.step, a.cost_ratio, a.nnp_window, a.boot, a.folds, a.reps, views=a.views)
    return 0


if __name__ == "__main__":
    sys.exit(main())
