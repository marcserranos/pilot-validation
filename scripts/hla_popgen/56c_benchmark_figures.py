#!/usr/bin/env python3
"""Figures for the 56/56b phasing & typing benchmark, drawn from the committed PUBLIC tables
(no VM, no participant data). Counts 1-19 are never drawn; masked cells are hatched.

    python3 scripts/hla_popgen/56c_benchmark_figures.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(REPO, "reports", "hla_popgen", "56_phasing_typing_benchmark")
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7", "#008300", "#e34948"]
BLUE, ORANGE, GREEN, RED = SERIES[0], SERIES[1], SERIES[2], SERIES[7]
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e6e5e0"

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "axes.axisbelow": True, "figure.dpi": 150, "savefig.bbox": "tight", "legend.frameon": False,
    "font.family": "DejaVu Sans", "hatch.linewidth": 0.8})


def num(s):
    return pd.to_numeric(s, errors="coerce")


def save(fig, name):
    fig.savefig(os.path.join(D, name))
    plt.close(fig)
    print("wrote", os.path.join("reports/hla_popgen/56_phasing_typing_benchmark", name))


def sweep():
    return pd.read_csv(f"{D}/refined_sweep.tsv", sep="\t", dtype=str)


def row(sw, check, item, ms):
    r = sw[(sw.check == check) & (sw["item"] == item) & (sw.min_score == str(ms))]
    return r.iloc[0]


def fig_sweep_typing():
    sw = sweep()
    xs = list(range(0, 8))  # score 7 = all seven points; very few pair-genes reach it
    rows = [row(sw, "gene_nonsharing", "ALL", t) for t in xs]
    pct = np.array([num(r.pct) for r in rows])
    lo = np.array([num(r.lo) for r in rows])
    hi = np.array([num(r.hi) for r in rows])
    n = np.array([num(r.n) for r in rows])
    fig, (a, b) = plt.subplots(1, 2, figsize=(8.2, 3.2), gridspec_kw={"width_ratios": [1.35, 1]})
    a.errorbar(xs[:7], pct[:7], yerr=[pct[:7] - lo[:7], hi[:7] - pct[:7]], fmt="o-", color=BLUE, capsize=3,
               lw=1.6, ms=5)
    # score 7: 0 errors in a tiny set; its upper CI runs off the axis, so draw it hollow with an arrow
    ytop = 5.6
    a.plot([7], [pct[7]], "o", mfc="none", mec=BLUE, ms=5)
    a.annotate("", xy=(7, ytop), xytext=(7, pct[7]), arrowprops=dict(arrowstyle="->", color=BLUE, lw=1))
    a.text(6.85, 4.1, f"0 of {int(n[7])}\nCI to {hi[7]:.1f}%", ha="right", va="top", fontsize=8, color=INK2)
    a.set_ylim(0, ytop)
    a.set_xticks(xs)
    a.set_xlabel("minimum per-call confidence score (0-7)")
    a.set_ylabel("pair-genes with no shared allele (%)")
    a.set_title("Typing discordance in MHC-sharing relatives")
    a.annotate("4.1%", (0, pct[0]), xytext=(0.3, pct[0] + 0.85), ha="center", fontsize=9, color=INK)
    a.annotate("1.15%", (6, pct[6]), xytext=(5.3, 2.6), ha="center", fontsize=9, color=INK,
               arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    a.text(3.6, 4.9, "4.1% → 1.15% (about 3.5× lower)", ha="center", fontsize=9, color=BLUE)
    ret = n / n[0] * 100
    b.bar(xs, ret, color=GREEN, width=0.65)
    for x, r in zip(xs, ret):
        b.text(x, r + 2, f"{r:.0f}" if r >= 1 else f"{r:.1f}", ha="center", fontsize=8, color=INK)
    b.set_xticks(xs)
    b.set_ylim(0, 118)
    b.set_xlabel("minimum per-call confidence score (0-7)")
    b.set_ylabel("pair-genes retained (%)")
    b.set_title("Share of pair-genes kept")
    b.grid(axis="x", visible=False)
    fig.tight_layout(w_pad=2)
    save(fig, "fig_sweep_typing.png")


def hatched(ax, y, h, label_txt):
    ax.barh(y, 100, height=h, color="none", edgecolor=MUTED, hatch="////", lw=0.8)


def fig_haplotype_transmission():
    sw = sweep()
    pairs = ["DQA1~DQB1", "DPA1~DPB1", "DRB1~DQB1", "B~C", "A~B"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 3.4), gridspec_kw={"width_ratios": [1.9, 1]})
    y = np.arange(len(pairs))[::-1]
    xmax = 80
    for yi, p in zip(y, pairs):
        r = row(sw, "hap_not_transmitted", p, 0)
        if r.pct in ("", "nan") or pd.isna(num(r.pct)):
            a.barh(yi, 6, height=0.6, color="none", edgecolor=MUTED, hatch="////", lw=0.8)
            a.text(7, yi, f"<20 failures (masked), of {int(num(r.n))}", va="center", fontsize=8, color=INK2)
        else:
            v, l, h = num(r.pct), num(r.lo), num(r.hi)
            a.barh(yi, v, height=0.6, color=RED if p == "A~B" else BLUE)
            a.errorbar(v, yi, xerr=[[v - l], [h - v]], fmt="none", ecolor=INK, capsize=3, lw=1)
            a.text(h + 1.5, yi, f"{v:.1f}%  (n={int(num(r.n))})", va="center", fontsize=8, color=INK)
    a.set_yticks(y, pairs)
    a.set_xlim(0, xmax)
    a.set_xlabel("haplotype not transmitted to relative (%), score ≥ 0")
    a.set_title("By gene pair (both relatives same-contig)")
    a.grid(axis="y", visible=False)
    cats = [("Revio–Revio", "platform:revio-revio", BLUE), ("Sequel IIe–Sequel IIe", "platform:sequel2e-sequel2e", ORANGE)]
    for i, (lab, it, col) in enumerate(cats):
        r = row(sw, "hap_not_transmitted", "ALL", it)
        v, l, h = num(r.pct), num(r.lo), num(r.hi)
        b.bar(i, v, color=col, width=0.55)
        b.errorbar(i, v, yerr=[[v - l], [h - v]], fmt="none", ecolor=INK, capsize=3, lw=1)
        b.text(i, h + 0.5, f"{v:.1f}%\nn={int(num(r.n)):,}", ha="center", fontsize=8, color=INK)
    b.set_xticks([0, 1], [c[0].replace("–", "–\n") for c in cats])
    b.set_xlim(-0.6, 1.6)
    b.set_ylim(0, 18)
    b.set_ylabel("not transmitted (%)")
    b.set_title("By platform pair")
    b.grid(axis="x", visible=False)
    fig.tight_layout(w_pad=2)
    save(fig, "fig_haplotype_transmission.png")


def fig_pair_sharing():
    d = pd.read_csv(f"{D}/pair_sharing_distribution.tsv", sep="\t", dtype=str)
    vals = {r.n_classical_genes_no_shared_allele: r.n_pairs for _, r in d.iterrows()}
    labs = ["0", "1", "2", "3+ (<20 each)"]
    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    for i, k in enumerate(["0", "1", "2"]):
        v = int(vals[k])
        ax.bar(i, v, color=BLUE, width=0.65)
        ax.text(i, v + 8, f"{v}", ha="center", fontsize=9, color=INK)
    ax.bar(3, 20, color="none", edgecolor=MUTED, hatch="////", width=0.65, lw=0.8)
    ax.text(3, 28, "<20 each\n(masked)", ha="center", fontsize=8, color=INK2)
    ax.set_xticks(range(4), labs)
    ax.set_ylim(0, 440)
    ax.set_xlabel("classical genes (of 8) with no shared allele in the pair")
    ax.set_ylabel("first-degree pairs")
    ax.set_title("Most first-degree pairs share an MHC haplotype")
    ax.grid(axis="x", visible=False)
    ax.annotate("IBD0 siblings would sit at 5+", xy=(3.0, 75), xytext=(2.4, 200), fontsize=9, color=RED,
                ha="center", arrowprops=dict(arrowstyle="->", color=RED, lw=1))
    fig.tight_layout()
    save(fig, "fig_pair_sharing.png")


def fig_platform_quality():
    sw = sweep()
    cats = [("Revio–Revio", "platform:revio-revio", BLUE), ("Sequel IIe–Sequel IIe", "platform:sequel2e-sequel2e", ORANGE)]
    fig, axs = plt.subplots(1, 2, figsize=(8.2, 3.3))
    for ax, check, title, yl in [(axs[0], "gene_nonsharing", "Typing: no shared allele", "pair-genes with no shared allele (%)"),
                                 (axs[1], "hap_not_transmitted", "Phasing: haplotype not transmitted", "haplotypes not transmitted (%)")]:
        top = 0
        for i, (lab, it, col) in enumerate(cats):
            r = row(sw, check, "ALL", it)
            v, l, h = num(r.pct), num(r.lo), num(r.hi)
            ax.bar(i, v, color=col, width=0.55)
            ax.errorbar(i, v, yerr=[[v - l], [h - v]], fmt="none", ecolor=INK, capsize=3, lw=1)
            ax.text(i, h + 0.25, f"{v:.2f}%\nn={int(num(r.n)):,}", ha="center", fontsize=8, color=INK)
            top = max(top, h)
        ax.set_xticks([0, 1], [c[0].replace("–", "–\n") for c in cats])
        ax.set_xlim(-0.6, 1.6)
        ax.set_ylim(0, top * 1.3)
        ax.set_ylabel(yl)
        ax.set_title(title)
        ax.grid(axis="x", visible=False)
    fig.tight_layout(w_pad=2.5)
    save(fig, "fig_platform_quality.png")


if __name__ == "__main__":
    fig_sweep_typing()
    fig_haplotype_transmission()
    fig_pair_sharing()
    fig_platform_quality()
