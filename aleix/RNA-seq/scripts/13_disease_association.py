#!/usr/bin/env python3
"""13 -- W4: repertoire diversity vs chronic disease, at fixed depth
(report: reports/13_disease_association/).

The depth-controlled version of Cole's entropy-wide association (which used raw Shannon
entropy ~ phenotype + PCs + age + sex + total reads). For every disease with >= 50 cases
(report 11) and every repertoire measure:

    measure (standardised) ~ case + age + age^2 + sex + PC1-10 + log10 TRB reads

OLS with heteroscedasticity-robust (HC3) errors; effect = difference in SD units between
cases and controls. Measures: richness, Shannon, inverse Simpson and public fraction at a
fixed 823 reads, top-10 clonal share, and RAW Shannon (Cole's measure, for direct
comparison). Multiple testing: Benjamini-Hochberg and Bonferroni across all tests of the
main model.

One sensitivity model adds log(1 + conditions per EHR-year): people who are seen more often
accumulate more diagnoses of everything, and acute illness at a visit lowers circulating T
cells, so a disease effect that disappears with this term is about health-care contact, not
the disease.

  a  effect heatmap, disease x measure (dot = BH q < 0.05, ring = Bonferroni)
  b  raw vs fixed-depth Shannon effects, per disease: what depth control changes
  c  inverse Simpson effects with 95% CI, main vs utilisation-adjusted

Usage (from aleix/RNA-seq/, after 11 and 12):
  pixi run python3 -u scripts/13_disease_association.py
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import repfig as R  # noqa: E402

MIN_CASES = 50
MEASURES = [("shannon_at_depth", "Shannon (823 reads)"),
            ("inv_simpson_at_depth", "Inverse Simpson (823)"),
            ("richness_at_depth", "Richness (823)"),
            ("public_frac_at_depth", "Public fraction (823)"),
            ("top10_frac", "Top-10 clonal share"),
            ("shannon_raw", "Shannon, raw")]


def fit(y, X):
    import statsmodels.api as sm
    m = sm.OLS(y, sm.add_constant(X, has_constant="add")).fit(cov_type="HC3")
    return m.params["case"], m.bse["case"], m.pvalues["case"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    base = "~/pipeline_outputs/rnaseq"
    ap.add_argument("--person-table", default=f"{base}/analysis/person_table.pkl")
    ap.add_argument("--labels", default=f"{base}/analysis/disease_labels.pkl")
    ap.add_argument("--features", default=f"{base}/analysis/feature_sets.pkl")
    ap.add_argument("--counts", default=f"{base}/reports/11_disease_labels/case_counts.csv")
    ap.add_argument("--outdir", default=f"{base}/reports/13_disease_association")
    args = ap.parse_args()
    outdir = os.path.expanduser(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    from statsmodels.stats.multitest import multipletests

    pt = pd.read_pickle(os.path.expanduser(args.person_table)).set_index("research_id")
    L = pd.read_pickle(os.path.expanduser(args.labels))
    Fb = pd.read_pickle(os.path.expanduser(args.features))["base"]
    names = pd.read_csv(os.path.expanduser(args.counts)).set_index("key")["disease"].to_dict()
    util = np.log1p(pt["conditions_per_ehr_year"]).reindex(Fb.index)

    rows = []
    for key in L.columns:
        lab = L[key]
        n_case = int((lab == 1).sum())
        if n_case < MIN_CASES:
            continue
        for col, mname in MEASURES:
            y = pt[col].reindex(L.index)
            ok = lab.notna() & y.notna() & Fb.notna().all(1)
            yz = (y[ok] - y[ok].mean()) / y[ok].std()
            X = Fb[ok].copy()
            X.insert(0, "case", lab[ok])
            b, se, p = fit(yz, X)
            X2 = X.copy()
            X2["utilisation"] = util[ok].fillna(util.median())
            b2, se2, p2 = fit(yz, X2)
            rows.append({"key": key, "disease": names.get(key, key), "measure": col,
                         "measure_label": mname, "cases": int((lab[ok] == 1).sum()),
                         "controls": int((lab[ok] == 0).sum()), "effect_sd": b,
                         "lo": b - 1.96 * se, "hi": b + 1.96 * se, "p": p,
                         "effect_sd_util_adj": b2, "lo_util_adj": b2 - 1.96 * se2,
                         "hi_util_adj": b2 + 1.96 * se2, "p_util_adj": p2})
    res = pd.DataFrame(rows)
    res["q_bh"] = multipletests(res["p"], method="fdr_bh")[1]
    res["bonferroni"] = res["p"] < 0.05 / len(res)
    res["q_bh_util_adj"] = multipletests(res["p_util_adj"], method="fdr_bh")[1]
    res.to_csv(os.path.join(outdir, "association.csv"), index=False)

    # ---------------- figure ----------------
    import cnsplots as cns
    pal = R.nature()
    mp = cns.multipanel(max_width=540)
    piv = res.pivot(index="disease", columns="measure", values="effect_sd")
    q = res.pivot(index="disease", columns="measure", values="q_bh")
    bf = res.pivot(index="disease", columns="measure", values="bonferroni")
    order_m = [m for m, _ in MEASURES]
    piv, q, bf = piv[order_m], q[order_m], bf[order_m]
    piv = piv.loc[piv["inv_simpson_at_depth"].sort_values().index]
    q, bf = q.loc[piv.index], bf.loc[piv.index]
    H = max(120, 9 * len(piv))
    ax = mp.panel("a", width=150, height=H, margin_right=70, margin_bottom=40)
    lim = float(np.nanquantile(np.abs(piv.to_numpy()), 0.98)) or 0.3
    im = ax.imshow(piv.to_numpy(), cmap="bwr", vmin=-lim, vmax=lim, aspect="auto")
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            if q.iat[i, j] < 0.05:
                ax.scatter(j, i, s=4, color="black")
            if bool(bf.iat[i, j]):
                ax.scatter(j, i, s=16, facecolors="none", edgecolors="black", linewidths=0.5)
    ax.set_xticks(range(len(order_m)))
    ax.set_xticklabels([dict(MEASURES)[m] for m in order_m], rotation=45, ha="right", fontsize=5)
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels(piv.index, fontsize=5)
    R.colorbar(ax, im, "cases − controls (SD)")
    ax.set_title("Diversity vs disease, adjusted", pad=4)

    ax = mp.panel("b", width=110, height=110, margin_right=12, margin_bottom=40)
    w = res.pivot(index="disease", columns="measure", values="effect_sd")
    qq = res.pivot(index="disease", columns="measure", values="q_bh")
    sig_raw, sig_fix = qq["shannon_raw"] < 0.05, qq["shannon_at_depth"] < 0.05
    colr = [pal[0] if r_ and f_ else pal[3] if r_ else pal[2] if f_ else "0.75"
            for r_, f_ in zip(sig_raw, sig_fix)]
    ax.scatter(w["shannon_raw"], w["shannon_at_depth"], s=12, c=colr, linewidths=0)
    lo_, hi_ = np.nanmin(w[["shannon_raw", "shannon_at_depth"]].to_numpy()), np.nanmax(w[["shannon_raw", "shannon_at_depth"]].to_numpy())
    ax.plot([lo_, hi_], [lo_, hi_], color="0.6", lw=0.5, ls="--")
    ax.axhline(0, color="0.8", lw=0.4)
    ax.axvline(0, color="0.8", lw=0.4)
    for lab_, c_ in (("both q<0.05", pal[0]), ("raw only", pal[3]), ("fixed depth only", pal[2])):
        ax.scatter([], [], s=12, color=c_, label=lab_)
    ax.legend(frameon=False, fontsize=5)
    R.style(ax, "Raw vs fixed-depth Shannon", "Effect, raw Shannon (SD)", "Effect at 823 reads (SD)")

    ax = mp.panel("c", width=150, height=H, margin_right=12, margin_bottom=40)
    s = res[res["measure"] == "inv_simpson_at_depth"].sort_values("effect_sd").reset_index(drop=True)
    y = np.arange(len(s))
    ax.errorbar(s["effect_sd"], y + 0.15, xerr=[s["effect_sd"] - s["lo"], s["hi"] - s["effect_sd"]],
                fmt="o", ms=2.5, color=pal[0], elinewidth=0.6, label="main model")
    ax.errorbar(s["effect_sd_util_adj"], y - 0.15,
                xerr=[s["effect_sd_util_adj"] - s["lo_util_adj"], s["hi_util_adj"] - s["effect_sd_util_adj"]],
                fmt="o", ms=2.5, color=pal[3], elinewidth=0.6, label="+ health-care utilisation")
    ax.axvline(0, color="0.5", lw=0.5, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{d} ({c})" for d, c in zip(s["disease"], s["cases"])], fontsize=4.8)
    ax.legend(frameon=False, fontsize=5, loc="lower right")
    R.style(ax, "Inverse Simpson (dominant clones)", "Cases − controls (SD, 95% CI)", None)
    R.save("fig_disease_association", outdir)

    top = res.sort_values("p").head(25)
    print(top[["disease", "measure", "cases", "effect_sd", "q_bh", "bonferroni",
               "effect_sd_util_adj", "q_bh_util_adj"]].round(4).to_string(index=False), file=sys.stderr)
    pd.DataFrame([("tests", len(res)), ("diseases", res["key"].nunique()),
                  ("bh_q05", int((res["q_bh"] < 0.05).sum())), ("bonferroni", int(res["bonferroni"].sum())),
                  ("bh_q05_util_adj", int((res["q_bh_util_adj"] < 0.05).sum()))],
                 columns=["metric", "value"]).to_csv(os.path.join(outdir, "summary.csv"), index=False)


if __name__ == "__main__":
    main()
