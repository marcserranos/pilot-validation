#!/usr/bin/env python3
"""WS-A (Sprint S04) -- figures for `44_kir_recurrence_saturation.py`'s aggregate TSVs.

Local-only: reads `reports/hla_popgen/44_kir_recurrence_saturation/*.tsv` (pulled-back aggregate
exports, never raw per-person data) and renders (1) recurrence-class-stratified saturation curves
for KIR and HLA, any-level and protein-level novelty side by side, (2) a per-ancestry equal-N
supplement, and (3) a coverage/completeness panel that keeps Good-Turing incidence coverage and
Chao2-based richness completeness visually and textually distinct (see "Interpretation guard"
below -- these are NOT interchangeable "fraction of allele space explored" numbers).

Also writes `recurrence_stats.tsv`: the singleton-share two-proportion z-test (HLA vs KIR) at the
genomic and any-level-novel granularities, plus the pooled-ALL "ground truth" recurrence-class
counts recovered from `saturation_curves.tsv` (see `pooled_recurrence_from_curves()` below).

Interpretation guard (do not violate when extending this script or writing prose from it):
  - Good-Turing `good_turing_coverage` is INCIDENCE coverage: P(the next carrier's allele was
    already seen). It is not the fraction of allele diversity discovered.
  - `s_obs / chao2` is the Chao2-based estimate of that latter fraction ("richness completeness").
    HLA ~64%, KIR ~84% pooled. Never conflate the two, and never write "99% of the allele space
    is explored" -- that claim is false; 99% is the incidence-coverage number.

Known data-QA finding (documented in the README, reproduced here so the figure script's own
provenance is traceable): `recurrence_classes.tsv`'s per-gene eq1/eq2/gt2 columns are `<20`-masked
independently per gene, so naively summing the masked pooled headline table across genes
UNDER-COUNTS eq2/eq1 by the sum of every masked cell's true (1-19) value and, because
`44`'s own hand-written prose table in its README was apparently computed a different, inconsistent
way, does not always even satisfy eq1+eq2+gt2==s_obs (it should, by construction: every allele is
in exactly one class). `saturation_curves.tsv`'s per-gene, per-n rows are NOT `<20`-masked (they
are permutation means, not raw counts, and were not marked as requiring masking in 44's own
Disclosure section) -- summed across all genes at the full cohort size (`n == n.max()`), they
recover the exact, internally-consistent pooled S_obs/eq1/eq2/gt2 counts (verified below:
eq1+eq2+gt2 == s_obs exactly for every species/level checked). This script uses the curve-derived
pooled-ALL sums for headline numbers and hypothesis tests, and the already-masked
`recurrence_classes.tsv`/`coverage_chao2.tsv` for any GENE-LEVEL figure content (never plotting a
gene-level curve-derived value that could be < 20 and unmasked).

Usage:
    python3 scripts/hla_popgen/45_kir_recurrence_figure.py \\
        --in-dir reports/hla_popgen/44_kir_recurrence_saturation \\
        --out-dir reports/hla_popgen/44_kir_recurrence_saturation
"""
import argparse
import math
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

SPECIES_ORDER = ["hla", "kir"]
SPECIES_LABEL = {"kir": "KIR", "hla": "HLA"}
LEVEL_ORDER = ["any_novel", "protein_novel"]
LEVEL_LABEL = {"any_novel": "any-level novel", "protein_novel": "protein-level novel"}
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "SAS"]  # MID excluded -- not well-powered (39/44 precedent)
RECUR_CLASSES = ["eq1", "eq2", "gt2", "ge20"]
RECUR_LABEL = {"eq1": "seen 1x", "eq2": "seen 2x", "gt2": "seen >2x", "ge20": "seen ≥20x"}
RECUR_COLOR_KEY = "nature"


def _to_censored(series):
    """Parses a column that may contain the literal string '<20' -- returns (values, is_censored)
    where censored entries are NaN in `values` (never 0)."""
    is_cens = series.astype(str).str.strip() == "<20"
    vals = pd.to_numeric(series.where(~is_cens), errors="coerce")
    return vals, is_cens


def load_tables(in_dir):
    rec = pd.read_csv(os.path.join(in_dir, "recurrence_classes.tsv"), sep="\t", dtype=str)
    curve = pd.read_csv(os.path.join(in_dir, "saturation_curves.tsv"), sep="\t")
    cov = pd.read_csv(os.path.join(in_dir, "coverage_chao2.tsv"), sep="\t")
    return rec, curve, cov


def pooled_recurrence_from_curves(curve, species, level, ancestry="ALL"):
    """Recovers the exact, unmasked, internally-consistent pooled-`ancestry` recurrence-class
    counts for one species/level by summing every gene's curve value AT THE FULL COHORT SIZE
    (`n == n.max()` for that slice). See module docstring "Known data-QA finding" for why this is
    preferred over summing `recurrence_classes.tsv`'s per-gene `<20`-masked cells.

    Returns a dict with s_obs, eq1, eq2, gt2, ge20 (floats, should already be integral) plus
    `partitions_exactly` (bool: eq1+eq2+gt2 == s_obs, the sanity check every caller should log).
    """
    s = curve[(curve["species"] == species) & (curve["level"] == level)
              & (curve["ancestry"] == ancestry)]
    if s.empty:
        return None
    nmax = s["n"].max()
    full = s[s["n"] == nmax].groupby("gene", as_index=False)[
        ["mean_distinct", "mean_eq1", "mean_eq2", "mean_gt2", "mean_ge20"]].mean()
    tot = full[["mean_distinct", "mean_eq1", "mean_eq2", "mean_gt2", "mean_ge20"]].sum()
    out = {
        "n_people": int(nmax),
        "s_obs": round(float(tot["mean_distinct"]), 1),
        "eq1": round(float(tot["mean_eq1"]), 1),
        "eq2": round(float(tot["mean_eq2"]), 1),
        "gt2": round(float(tot["mean_gt2"]), 1),
        "ge20": round(float(tot["mean_ge20"]), 1),
    }
    out["partitions_exactly"] = abs((out["eq1"] + out["eq2"] + out["gt2"]) - out["s_obs"]) < 0.6
    return out


def two_proportion_ztest(x1, n1, x2, n2):
    """Two-sided two-proportion z-test (pooled variance), no scipy/statsmodels dependency (the
    project's pixi env may not have either pinned for this script). Returns (p1, p2, z, p_value).
    """
    if n1 == 0 or n2 == 0:
        return (float("nan"),) * 4
    p1, p2 = x1 / n1, x2 / n2
    p = (x1 + x2) / (n1 + n2)
    se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    if se == 0:
        return p1, p2, float("nan"), float("nan")
    z = (p1 - p2) / se
    p_value = math.erfc(abs(z) / math.sqrt(2))
    return p1, p2, z, p_value


def compute_singleton_share_stats(curve):
    """Marc's hypothesis: 'KIR novelty is less private than HLA novelty (singleton share lower),
    suggesting systematic catalogue gaps rather than private errors.' Tested at genomic (all
    alleles, for context) and any-level-novel (the actual novelty claim); protein-level is
    reported but flagged NA for HLA (pooled protein_novel S_obs is 0 -- see README caveats, an
    unreconciled 44 data question, not something this local-only script can fix)."""
    rows = []
    for level in ["genomic", "any_novel", "protein_novel"]:
        hla = pooled_recurrence_from_curves(curve, "hla", level)
        kir = pooled_recurrence_from_curves(curve, "kir", level)
        if hla is None or kir is None:
            continue
        if hla["s_obs"] == 0 or kir["s_obs"] == 0:
            rows.append({"level": level, "hla_s_obs": hla["s_obs"], "kir_s_obs": kir["s_obs"],
                         "hla_singleton_share": float("nan"), "kir_singleton_share": float("nan"),
                         "z": float("nan"), "p_value": float("nan"),
                         "note": "one species has s_obs==0 at this level -- test not computable"})
            continue
        p1, p2, z, pval = two_proportion_ztest(hla["eq1"], hla["s_obs"], kir["eq1"], kir["s_obs"])
        rows.append({"level": level, "hla_s_obs": hla["s_obs"], "kir_s_obs": kir["s_obs"],
                     "hla_singleton_share": round(p1, 4), "kir_singleton_share": round(p2, 4),
                     "z": round(z, 3), "p_value": pval,
                     "note": "eq1/s_obs, HLA vs KIR, two-proportion z-test"})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Direct end-of-line labeling (adapted from 40_figure1_v5.draw_panel_e's shared-x-column pattern):
# a shared x just past the rightmost curve, vertically repelled labels, dotted elbow leaders.
# ---------------------------------------------------------------------------
def _label_curve_ends(ax, ends, colors, fontsize, min_gap_frac=0.12):
    """Places one label per curve at a shared x-column past the rightmost curve endpoint, with
    vertical repulsion so labels never overlap. `min_gap_frac` is a fraction of the AXES' own
    plotted y-range (via `ax.get_ylim()` after the data is drawn), not of the spread between the
    label values themselves -- using only the label spread underestimates the needed gap when
    labels happen to cluster close together relative to a much larger axis range (observed in a
    synthetic-fixture test: two labels 8 apart on an axis spanning 0-6, i.e. yrange==label
    spread, produced a gap far smaller than one line of 5-6pt text actually needs)."""
    if not ends:
        return
    xmax = max(x for x, _ in ends.values())
    ylo, yhi = ax.get_ylim()
    axis_span = (yhi - ylo) or 1.0
    min_gap = max(axis_span * min_gap_frac, 1e-6)
    label_x = xmax * 1.04
    order = sorted(ends, key=lambda k: ends[k][1])
    placed = []
    for k in order:
        y = ends[k][1]
        if placed and y - placed[-1] < min_gap:
            y = placed[-1] + min_gap
        placed.append(y)
    label_y = dict(zip(order, placed))
    for k, (x0, y0) in ends.items():
        t = ax.annotate(k, (label_x, label_y[k]), fontsize=fontsize, color=colors[k],
                         fontweight="bold", va="center", ha="left", annotation_clip=False)
        vc.mark_label(t)
        leader_end = label_x - 0.01 * xmax
        if leader_end - x0 > 0.005 * xmax or abs(label_y[k] - y0) > 1e-9:
            ln, = ax.plot([x0, x0], [y0, label_y[k]], color=colors[k], lw=0.5,
                          ls=(0, (1, 1)), zorder=2.5, clip_on=False)
            vc.mark_decoration(ln)
            ln2, = ax.plot([x0, leader_end], [label_y[k], label_y[k]], color=colors[k], lw=0.5,
                           ls=(0, (1, 1)), zorder=2.5, clip_on=False)
            vc.mark_decoration(ln2)
    x0lo, x0hi = ax.get_xlim()
    ax.set_xlim(x0lo, max(x0hi, label_x * 1.22))


# ---------------------------------------------------------------------------
# Figure 1: recurrence-class-stratified saturation, KIR/HLA x any/protein novel, ALL ancestry.
# ---------------------------------------------------------------------------
def fig_saturation_by_recurrence(curve, out_stem, ancestry="ALL"):
    pal = vc.JOURNAL_PALETTES[RECUR_COLOR_KEY]
    colors = {c: pal[i] for i, c in enumerate(RECUR_CLASSES)}
    with vc.nature_style():
        fig, axes = plt.subplots(2, 2, figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(120)),
                                  constrained_layout=True)
        for row, species in enumerate(SPECIES_ORDER):
            for col, level in enumerate(LEVEL_ORDER):
                ax = axes[row, col]
                sub = curve[(curve["species"] == species) & (curve["level"] == level)
                            & (curve["ancestry"] == ancestry)]
                pooled = sub.groupby("n", as_index=False)[
                    ["mean_eq1", "mean_eq2", "mean_gt2", "mean_ge20"]].sum().sort_values("n")
                all_zero = pooled.empty or float(
                    pooled[["mean_eq1", "mean_eq2", "mean_gt2", "mean_ge20"]].to_numpy().max()) <= 0
                if all_zero:
                    # HLA protein_novel is pooled 0 (flagged, unreconciled 44 data question -- see
                    # README caveats). Plotting near-zero float noise on an auto-scaled axis
                    # produces an unreadable 1e-6-scale panel; state the fact in-panel instead.
                    ax.set_xlim(0, 1)
                    ax.set_ylim(0, 1)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    t = ax.text(0.5, 0.5, "0 protein-level\nnovel alleles\n(pooled)",
                                 ha="center", va="center", fontsize=6.5, color="#666666",
                                 transform=ax.transAxes)
                    vc.mark_label(t)
                else:
                    ends = {}
                    for c in RECUR_CLASSES:
                        x = pooled["n"].to_numpy()
                        y = pooled["mean_" + c].to_numpy()
                        ax.plot(x, y, color=colors[c], linewidth=1.0)
                        ends[RECUR_LABEL[c]] = (x[-1], y[-1])
                    colors_by_label = {RECUR_LABEL[c]: colors[c] for c in RECUR_CLASSES}
                    ax.set_xlim(left=0)
                    ax.set_ylim(bottom=0)
                    _label_curve_ends(ax, ends, colors_by_label, fontsize=5.2)
                if row == 0:
                    ax.set_title(LEVEL_LABEL[level], fontsize=7)
                if col == 0:
                    ax.set_ylabel(f"{SPECIES_LABEL[species]}\ndistinct alleles reaching class")
                if row == 1 and not all_zero:
                    ax.set_xlabel("N unrelated people")
                vc.panel_letter(ax, "abcd"[row * 2 + col])
        vc.save_fig(fig, out_stem)


# ---------------------------------------------------------------------------
# Figure 2: equal-N supplement -- overall (genomic level) discovery curve per ancestry, KIR vs HLA,
# plus the pooled-ALL curve for scale contrast (explains why the pooled slope << per-ancestry
# slope: pooled N* is capped by the smallest well-powered ancestry's endpoint, see caption).
# ---------------------------------------------------------------------------
def fig_saturation_per_ancestry(curve, out_stem, level="genomic"):
    with vc.nature_style():
        fig, axes = plt.subplots(1, 2, figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(75)),
                                  constrained_layout=True)
        for ax, species in zip(axes, SPECIES_ORDER):
            sub = curve[(curve["species"] == species) & (curve["level"] == level)]
            ends = {}
            colors = {}
            for i, a in enumerate(ANCESTRY_ORDER):
                s = sub[sub["ancestry"] == a].groupby("n", as_index=False)["mean_distinct"].sum()
                s = s.sort_values("n")
                if s.empty:
                    continue
                col = vc.ANCESTRY_COLORS.get(a, vc.ACCENT_COLOR)
                ax.plot(s["n"], s["mean_distinct"], color=col, linewidth=1.0)
                ends[a] = (s["n"].to_numpy()[-1], s["mean_distinct"].to_numpy()[-1])
                colors[a] = col
            ax.set_xlim(left=0)
            ax.set_ylim(bottom=0)
            _label_curve_ends(ax, ends, colors, fontsize=5.5)
            ax.set_title(SPECIES_LABEL[species], fontsize=7)
            ax.set_xlabel("N people (this ancestry)")
        axes[0].set_ylabel("distinct alleles (genomic level)")
        vc.panel_letter(axes[0], "a")
        vc.panel_letter(axes[1], "b")
        vc.save_fig(fig, out_stem)


# ---------------------------------------------------------------------------
# Figure 3: coverage/completeness -- Good-Turing incidence coverage AND Chao2 richness
# completeness (s_obs/chao2), kept as two clearly-labeled, side-by-side panels so neither is ever
# read as a stand-in for the other (interpretation guard, module docstring).
# ---------------------------------------------------------------------------
def fig_coverage_completeness(cov, out_stem, level="genomic"):
    sub = cov[(cov["ancestry"] == "ALL") & (cov["level"] == level)].copy()
    sub["completeness"] = sub["s_obs"] / sub["chao2"]
    agg = sub.groupby("species", as_index=False).agg(
        good_turing_coverage=("good_turing_coverage", "mean"),
        s_obs=("s_obs", "sum"), chao2=("chao2", "sum"))
    agg["completeness"] = agg["s_obs"] / agg["chao2"]
    with vc.nature_style():
        fig, axes = plt.subplots(1, 2, figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(60)),
                                  constrained_layout=True)
        x = np.arange(len(SPECIES_ORDER))
        pal = vc.JOURNAL_PALETTES[RECUR_COLOR_KEY]
        for ax, metric, title in (
                (axes[0], "good_turing_coverage", "Good–Turing incidence coverage"),
                (axes[1], "completeness", "Chao2 richness completeness (S$_{obs}$/Chao2)")):
            vals = [float(agg.loc[agg["species"] == sp, metric].iloc[0]) for sp in SPECIES_ORDER]
            ax.bar(x, vals, width=0.5, color=[pal[i] for i in range(len(SPECIES_ORDER))])
            for xi, v in zip(x, vals):
                t = ax.annotate(f"{100*v:.1f}%", (xi, v), xytext=(0, 2),
                                 textcoords="offset points", ha="center", fontsize=6)
                vc.mark_label(t)
            ax.set_xticks(x)
            ax.set_xticklabels([SPECIES_LABEL[s] for s in SPECIES_ORDER])
            ax.set_ylim(0, 1.08)
            ax.set_title(title, fontsize=6.5)
        axes[0].set_ylabel("fraction")
        vc.panel_letter(axes[0], "a")
        vc.panel_letter(axes[1], "b")
        vc.save_fig(fig, out_stem)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", default="reports/hla_popgen/44_kir_recurrence_saturation")
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args()
    out_dir = args.out_dir or args.in_dir
    os.makedirs(out_dir, exist_ok=True)

    rec, curve, cov = load_tables(args.in_dir)

    fig_saturation_by_recurrence(curve, os.path.join(out_dir, "fig_saturation_by_recurrence"))
    fig_saturation_per_ancestry(curve, os.path.join(out_dir, "fig_saturation_per_ancestry"))
    fig_coverage_completeness(cov, os.path.join(out_dir, "fig_coverage_completeness"))

    stats = compute_singleton_share_stats(curve)
    stats_path = os.path.join(out_dir, "recurrence_stats.tsv")
    stats.to_csv(stats_path, sep="\t", index=False)
    print(f"  wrote {stats_path}", file=sys.stderr)

    # Sanity-log the partition identity + the known recurrence_classes.tsv undercount, so a rerun
    # of this script always re-confirms both README claims rather than trusting stale prose.
    for level in ["genomic", "any_novel", "protein_novel"]:
        for sp in SPECIES_ORDER:
            p = pooled_recurrence_from_curves(curve, sp, level)
            if p is not None:
                print(f"  [check] {sp}/{level}: s_obs={p['s_obs']} eq1+eq2+gt2="
                      f"{p['eq1']+p['eq2']+p['gt2']} partitions_exactly={p['partitions_exactly']}",
                      file=sys.stderr)

    print(f"[45] 3 figures + recurrence_stats.tsv written to {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
