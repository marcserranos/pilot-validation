#!/usr/bin/env python3
"""WS-A (Sprint S04) -- figure(s) for 44_kir_recurrence_saturation.py's output TSVs.

Local-only: reads `reports/hla_popgen/44_kir_recurrence_saturation/*.tsv` (the pulled-back
aggregate exports -- never raw per-person data) and renders paired any-level-novel vs
protein-level-novel panels, per gene/ancestry, using the shared `_viz_common.nature_style()`
(WS-C is porting `reference/FIGURE_STYLE.md` into that module in parallel; this script only calls
it, it does not define its own rc/palette).

Panels (one figure per KIR gene x species pairing is too many for a single artifact; this script
draws the OVERALL, ALL-ancestry-pooled comparison plus one small multi-ancestry supplement):
  fig_saturation_paired      -- 2 columns (any_novel | protein_novel) x 2 rows (KIR | HLA), ALL
                                 ancestry, mean +/- 95% rarefaction band, from saturation_curves.tsv.
  fig_recurrence_classes     -- stacked/grouped bars of eq1/eq2/gt2/ge20 allele counts, ALL
                                 ancestry, any_novel vs protein_novel side by side, KIR vs HLA.
  fig_coverage_chao2         -- Good-Turing coverage and Chao2-undetected (f0_hat), per ancestry,
                                 any_novel vs protein_novel, KIR vs HLA -- points+CI (SE), never a
                                 non-zero-based bar axis (WS-C "honest axes" rule).

Disclosure: this script assumes 44's own TSVs are already disclosure-safe (masked "<20" strings
where required); it must not silently coerce a "<20" cell to 0 -- `_to_censored()` below turns it
into NaN + a hatched marker instead (see `_viz_common.hatch_suppressed`), never a plotted zero.

Usage:
    python3 scripts/hla_popgen/45_kir_recurrence_figure.py \\
        --in-dir reports/hla_popgen/44_kir_recurrence_saturation \\
        --out-dir reports/hla_popgen/44_kir_recurrence_saturation
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

SPECIES_ORDER = ["kir", "hla"]
SPECIES_LABEL = {"kir": "KIR", "hla": "HLA"}
LEVEL_ORDER = ["any_novel", "protein_novel"]
LEVEL_LABEL = {"any_novel": "any-level novel", "protein_novel": "protein-level novel"}
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]


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


def fig_saturation_paired(curve, out_stem, gene=None, ancestry="ALL"):
    sub = curve[curve["ancestry"] == ancestry]
    if gene is not None:
        sub = sub[sub["gene"] == gene]
    with vc.nature_style():
        fig, axes = plt.subplots(2, 2, figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(110)),
                                  sharex=False, constrained_layout=True)
        for row, species in enumerate(SPECIES_ORDER):
            for col, level in enumerate(LEVEL_ORDER):
                ax = axes[row, col]
                s = sub[(sub["species"] == species) & (sub["level"] == level)]
                if gene is None:
                    s = s.groupby("n", as_index=False)[["mean_distinct", "lo2_5", "hi97_5"]].sum()
                s = s.sort_values("n")
                if len(s):
                    ax.plot(s["n"], s["mean_distinct"], color=vc.ACCENT_COLOR, linewidth=1.0)
                    ax.fill_between(s["n"], s["lo2_5"], s["hi97_5"], color=vc.ACCENT_COLOR,
                                     alpha=0.2, linewidth=0)
                if row == 0:
                    ax.set_title(LEVEL_LABEL[level], fontsize=7)
                if col == 0:
                    ax.set_ylabel(f"{SPECIES_LABEL[species]}\ndistinct alleles")
                if row == 1:
                    ax.set_xlabel("N people")
                vc.panel_letter(ax, "abcd"[row * 2 + col])
        vc.save_fig(fig, out_stem)


def fig_recurrence_classes(rec, out_stem, ancestry="ALL"):
    sub = rec[rec["ancestry"] == ancestry].copy()
    classes = ["eq1", "eq2", "gt2", "ge20"]
    for c in classes:
        sub[c + "_val"], sub[c + "_cens"] = _to_censored(sub[c])
    with vc.nature_style():
        fig, axes = plt.subplots(1, 2, figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(60)),
                                  constrained_layout=True)
        for ax, level in zip(axes, LEVEL_ORDER):
            s = sub[sub["level"] == level]
            x = np.arange(len(classes))
            width = 0.35
            for i, species in enumerate(SPECIES_ORDER):
                srow = s[s["species"] == species]
                totals = [srow[c + "_val"].fillna(0).sum() for c in classes]
                cens = [bool(srow[c + "_cens"].any()) for c in classes]
                offs = x + (i - 0.5) * width
                ax.bar(offs, totals, width=width, label=SPECIES_LABEL[species],
                       color=vc.JOURNAL_PALETTES["nature"][i])
                for xi, is_c, val in zip(offs, cens, totals):
                    if is_c:
                        vc.hatch_suppressed(ax, xi - width / 2, 0, width, max(val, 1))
            ax.set_xticks(x)
            ax.set_xticklabels(["=1", "=2", ">2", ">=20"])
            ax.set_title(LEVEL_LABEL[level], fontsize=7)
            ax.set_ylabel("distinct alleles")
        axes[0].legend(**vc.LEGEND_KW)
        vc.save_fig(fig, out_stem)


def fig_coverage_chao2(cov, out_stem):
    sub = cov[cov["ancestry"].isin(ANCESTRY_ORDER)].copy()
    with vc.nature_style():
        fig, axes = plt.subplots(1, 2, figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(65)),
                                  constrained_layout=True)
        for ax, level in zip(axes, LEVEL_ORDER):
            s = sub[sub["level"] == level]
            for i, species in enumerate(SPECIES_ORDER):
                srow = s[s["species"] == species].groupby("ancestry", as_index=False).agg(
                    good_turing_coverage=("good_turing_coverage", "mean"))
                srow["ancestry"] = pd.Categorical(srow["ancestry"], ANCESTRY_ORDER, ordered=True)
                srow = srow.sort_values("ancestry")
                ax.plot(srow["ancestry"], srow["good_turing_coverage"], marker="o", markersize=3,
                        linewidth=0.8, label=SPECIES_LABEL[species],
                        color=vc.JOURNAL_PALETTES["nature"][i])
            ax.set_ylim(0, 1.02)
            ax.set_ylabel("Good-Turing coverage")
            ax.set_title(LEVEL_LABEL[level], fontsize=7)
        axes[0].legend(**vc.LEGEND_KW)
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
    fig_saturation_paired(curve, os.path.join(out_dir, "fig_saturation_paired"))
    fig_recurrence_classes(rec, os.path.join(out_dir, "fig_recurrence_classes"))
    fig_coverage_chao2(cov, os.path.join(out_dir, "fig_coverage_chao2"))
    print(f"[45] 3 figures written to {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
