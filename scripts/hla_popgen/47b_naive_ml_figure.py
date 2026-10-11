#!/usr/bin/env python3
"""WS-D (Sprint S04) -- figure for `47_naive_ml_tests.py`'s aggregate TSVs.

Local-only: reads `reports/hla_popgen/47_naive_ml_tests/*.tsv` (pulled-back aggregate exports,
never raw per-person data) and renders one figure with two panels:

  (a) Per-task observed AUROC (a point) against its permutation-null RANGE (grey band, min-max
      over the `--n-perms` shuffles), all nine tasks on one shared 0.5-1.0 axis, grouped into
      three blocks (ancestry / platform / cB-from-HLA) with direct labels -- no legend. The
      message this panel exists to carry: ancestry is strongly encoded in HLA+KIR carriage,
      platform weakly, and HLA-derived KIR gene content not at all.
  (b) The top ancestry-informative HLA alleles (by |L1-logistic-regression coefficient|),
      restricted to features with >=20 carriers (the same floor `47` already enforced before
      modeling -- re-stated here as a belt-and-suspenders filter against a stale/edited TSV).

Interpretation guard (do not violate when extending this script or writing prose from it):
  - `47` now stores the full permutation-null distribution per task (min/5th/25th/50th/75th/95th
    percentile/max, `--n-perms` default 200), not just a mean/p95 pair. This script's `null_band()`
    prefers the exact 5-95th-percentile columns (`perm_auroc_p05`/`perm_auroc_p95`) when present,
    with min-max drawn as whiskers -- both the band and the whiskers are then the real shuffle
    distribution, not an approximation.
  - Fallback: an older `naive_ml_metrics.tsv` (pre-percentile-columns run) only has
    `perm_auroc_mean`/`perm_auroc_p95`. In that case `null_band()` reconstructs a defensible band
    using `perm_auroc_p95` as the upper edge and `2*perm_auroc_mean - perm_auroc_p95` (mirrored
    around the mean) as a conservative lower edge, floored at 0.5 -- a symmetric approximation of
    the null's spread, not the true range. The panel's legend text says which case it is.
  - Even with 200 shuffles the empirical p-value floors at 1/201 ~= 0.005 -- still not a substitute
    for reading the AUROC directly against its null band, which is what this panel plots; it never
    plots or claims a p-value.

Usage:
    python3 scripts/hla_popgen/47b_naive_ml_figure.py \\
        --in-dir reports/hla_popgen/47_naive_ml_tests \\
        --out-dir reports/hla_popgen/47_naive_ml_tests
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)
import _viz_common as vc

MIN_CARRIERS_FOR_FEATURE = 20

# NOTE: platform_sequel2e[_EUR/_AFR] are omitted from the panel on purpose -- within a
# two-valued platform subset, a one-vs-rest AUROC for "is revio" and "is sequel2e" are the same
# ROC curve read from opposite ends and are numerically identical (verified: both 0.573 overall,
# both 0.558 within EUR, both 0.489 within AFR in the v2 run). Plotting both would visually
# double-count one effect; the README states the equivalence and both rows are still in
# naive_ml_metrics.tsv for anyone who wants sequel2e's own row.
TASK_ORDER = ["ancestry_AFR", "ancestry_AMR", "ancestry_EAS", "ancestry_EUR", "ancestry_MID",
              "ancestry_SAS", "platform_revio", "platform_revio_EUR", "platform_revio_AFR",
              "cB_from_HLA", "cB_from_HLA_EUR"]
TASK_LABEL = {
    "ancestry_AFR": "ancestry: AFR", "ancestry_AMR": "ancestry: AMR",
    "ancestry_EAS": "ancestry: EAS", "ancestry_EUR": "ancestry: EUR",
    "ancestry_MID": "ancestry: MID", "ancestry_SAS": "ancestry: SAS",
    "platform_revio": "platform\n(all)", "platform_revio_EUR": "platform\n(EUR only)",
    "platform_revio_AFR": "platform\n(AFR only)",
    "cB_from_HLA": "KIR cB\n(all)", "cB_from_HLA_EUR": "KIR cB\n(EUR only)",
}
TASK_GROUP = {
    "ancestry_AFR": "ancestry", "ancestry_AMR": "ancestry", "ancestry_EAS": "ancestry",
    "ancestry_EUR": "ancestry", "ancestry_MID": "ancestry", "ancestry_SAS": "ancestry",
    "platform_revio": "platform", "platform_revio_EUR": "platform", "platform_revio_AFR": "platform",
    "cB_from_HLA": "cB_from_HLA", "cB_from_HLA_EUR": "cB_from_HLA",
}
GROUP_COLOR = {"ancestry": vc.JOURNAL_PALETTES["nature"][0],
               "platform": vc.JOURNAL_PALETTES["nature"][1],
               "cB_from_HLA": vc.JOURNAL_PALETTES["nature"][3]}


def load_tables(in_dir):
    metrics = pd.read_csv(os.path.join(in_dir, "naive_ml_metrics.tsv"), sep="\t")
    feats = pd.read_csv(os.path.join(in_dir, "naive_ml_top_features.tsv"), sep="\t")
    return metrics, feats


def _has(row, col):
    return col in row.index if hasattr(row, "index") else col in row


def null_band(row):
    """Returns (lo, hi, exact) for the permutation-null band.

    Preferred path: `naive_ml_metrics.tsv` now stores the actual permutation-null
    5th/95th-percentile columns (`perm_auroc_p05`/`perm_auroc_p95`) alongside min/max -- when
    present, the band is the exact 5-95% range (`exact=True`) and whiskers are drawn out to the
    stored min/max.

    Fallback (older TSVs without those columns, e.g. from before this run's --n-perms increase):
    reconstruct a symmetric approximation from `perm_auroc_mean`/`perm_auroc_p95` alone --
    `exact=False` -- see module docstring's Interpretation guard for why this is only an
    approximation, not the true shuffle range."""
    if _has(row, "perm_auroc_p05") and row.get("perm_auroc_p05") == row.get("perm_auroc_p05") \
            and str(row.get("perm_auroc_p05")) != "":
        lo = float(row["perm_auroc_p05"])
        hi = float(row["perm_auroc_p95"])
        return lo, hi, True
    mean = float(row["perm_auroc_mean"])
    p95 = float(row["perm_auroc_p95"])
    hi = max(p95, mean)
    lo = max(0.5, 2 * mean - hi)
    return lo, hi, False


def whiskers(row):
    """Returns (min, max) of the true permutation-null draws when the min/max columns are
    present, else None -- used to draw min-max whiskers on top of the 5-95% band (exact case
    only; the mirrored-approximation fallback has no real min/max to draw)."""
    if _has(row, "perm_auroc_min") and row.get("perm_auroc_min") == row.get("perm_auroc_min") \
            and str(row.get("perm_auroc_min")) != "":
        return float(row["perm_auroc_min"]), float(row["perm_auroc_max"])
    return None


def fig_auroc_vs_null(metrics, feats, out_stem, min_carriers=MIN_CARRIERS_FOR_FEATURE):
    tasks = [t for t in TASK_ORDER if t in set(metrics["task"])]
    m = metrics.set_index("task")

    with vc.nature_style():
        fig, (ax_a, ax_b) = plt.subplots(
            2, 1, figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(150)),
            constrained_layout=True, gridspec_kw={"height_ratios": [1.1, 1.0]})

        # --- Panel a: AUROC vs. permutation-null band, all 9 tasks, one shared 0.5-1.0 axis ---
        x = np.arange(len(tasks))
        any_exact = False
        for xi, t in enumerate(tasks):
            row = m.loc[t]
            lo, hi, exact = null_band(row)
            any_exact = any_exact or exact
            grp = TASK_GROUP[t]
            color = GROUP_COLOR[grp]
            ax_a.add_patch(matplotlib.patches.Rectangle(
                (xi - 0.32, lo), 0.64, hi - lo, facecolor=vc.SUPPRESSED_COLOR,
                edgecolor="none", zorder=1))
            wk = whiskers(row)
            if wk is not None:
                wmin, wmax = wk
                ax_a.plot([xi, xi], [wmin, lo], color="#999999", lw=0.5, zorder=1.3)
                ax_a.plot([xi, xi], [hi, wmax], color="#999999", lw=0.5, zorder=1.3)
                ax_a.plot([xi - 0.14, xi + 0.14], [wmin, wmin], color="#999999", lw=0.5,
                         zorder=1.3)
                ax_a.plot([xi - 0.14, xi + 0.14], [wmax, wmax], color="#999999", lw=0.5,
                         zorder=1.3)
            ax_a.plot([xi - 0.32, xi + 0.32], [0.5, 0.5], color="#BBBBBB", lw=0.4,
                      zorder=1.2, ls=(0, (1, 1)))
            auroc = float(row["lr_auroc"])
            ax_a.scatter([xi], [auroc], s=16, color=color, zorder=3, edgecolor="white",
                         linewidth=0.4)
            # Offset sideways (not straight up) so the label text doesn't sit on top of the
            # vertical min/max whisker line, which is drawn at the same x position as the point.
            t_lab = ax_a.annotate(f"{auroc:.3f}", (xi, auroc), xytext=(9, 0),
                                   textcoords="offset points", ha="left", va="center",
                                   fontsize=5.2, color=color, fontweight="bold")
            vc.mark_label(t_lab)
        ax_a.set_xlim(-0.6, len(tasks) - 0.4)
        ax_a.set_ylim(0.45, 1.03)
        ax_a.set_xticks(x)
        ax_a.set_xticklabels([TASK_LABEL[t] for t in tasks], fontsize=5.0, rotation=0)
        for tick_lbl, t in zip(ax_a.get_xticklabels(), tasks):
            tick_lbl.set_color(GROUP_COLOR[TASK_GROUP[t]])
            tick_lbl.set_fontweight("bold")
        ax_a.set_ylabel("cross-validated AUROC")
        ax_a.axhline(0.5, color="#999999", lw=0.5, zorder=0.5)
        band_desc = ("grey band = permutation-null 5-95% range, whiskers = min-max" if any_exact
                    else "grey band = permutation-null AUROC range (mirrored mean/p95 "
                         "approximation)")
        legend_txt = ax_a.annotate(
            f"{band_desc}; point = observed. Label color: carriage->ancestry, "
            "carriage->platform, HLA->KIR cB (matches point color).",
            (0.0, 1.10), xycoords="axes fraction", ha="left", va="bottom", fontsize=4.6,
            color="#666666")
        vc.mark_label(legend_txt)
        vc.panel_letter(ax_a, "a")

        # --- Panel b: top ancestry-informative alleles, LR coefficients, >=20-carrier features ---
        anc_tasks = [t for t in tasks if TASK_GROUP[t] == "ancestry"]
        sub = feats[(feats["task"].isin(anc_tasks)) & (feats["model"] == "l1_logreg")].copy()
        sub["ancestry"] = sub["task"].str.replace("ancestry_", "", regex=False)
        sub["abs_weight"] = sub["weight"].abs()
        top = (sub.sort_values("abs_weight", ascending=False)
                  .groupby("ancestry", group_keys=False).head(3))
        top = top.sort_values(["ancestry", "abs_weight"], ascending=[True, True])
        top["row_label"] = top["ancestry"] + ": " + top["feature"]
        yy = np.arange(len(top))
        colors = [vc.ANCESTRY_COLORS.get(a, vc.ACCENT_COLOR) for a in top["ancestry"]]
        ax_b.barh(yy, top["weight"], color=colors, height=0.65)
        ax_b.axvline(0, color="black", lw=0.5)
        ax_b.set_yticks(yy)
        ax_b.set_yticklabels(top["row_label"], fontsize=5.0)
        ax_b.set_xlabel("L1 logistic-regression coefficient (sign = direction of association)")
        ax_b.set_title(
            f"top 3 ancestry-informative HLA alleles per ancestry (≥20 carriers)",
            fontsize=6.5)
        vc.panel_letter(ax_b, "b")

        vc.save_fig(fig, out_stem)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", default="reports/hla_popgen/47_naive_ml_tests")
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args()
    out_dir = args.out_dir or args.in_dir
    os.makedirs(out_dir, exist_ok=True)

    metrics, feats = load_tables(args.in_dir)
    fig_auroc_vs_null(metrics, feats, os.path.join(out_dir, "fig_naive_ml_auroc"))
    print(f"[47b] 1 figure written to {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
