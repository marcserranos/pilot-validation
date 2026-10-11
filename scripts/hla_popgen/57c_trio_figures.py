#!/usr/bin/env python3
"""Figures for 57/57b (yak trio benchmark), drawn from the committed PUBLIC tables only
(no VM, no participant data). Masked cells (<20) are not drawn.

    python3 scripts/hla_popgen/57c_trio_figures.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(REPO, "reports", "hla_popgen", "57_yak_trioeval")
D56 = os.path.join(REPO, "reports", "hla_popgen", "56_phasing_typing_benchmark")
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7", "#008300", "#e34948"]
BLUE, ORANGE, GREEN = SERIES[0], SERIES[1], SERIES[2]
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e6e5e0"
PLAT = {"revio": ("Revio", BLUE), "sequel2e": ("Sequel IIe", ORANGE), "all": ("all", INK2)}
TARGET_LABEL = {"MHC_classII_DR_DQ": "class II\nDR-DQ\n0.3 Mb", "KIR": "KIR\n\n0.2 Mb",
                "MHC": "whole\nMHC\n4 Mb", "MHC_classI_A_to_B": "class I\nA..B\n1.5 Mb",
                "chr20_control": "chr20\ncontrol\n55 Mb"}
FLOOR = 0.02  # zero medians/quartiles are drawn at this floor on the log axis (open marker)

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "axes.axisbelow": True, "figure.dpi": 150, "savefig.bbox": "tight", "legend.frameon": False,
    "font.family": "DejaVu Sans"})


def num(s):
    return pd.to_numeric(s, errors="coerce")


def save(fig, name):
    fig.savefig(os.path.join(D, name))
    plt.close(fig)
    print("wrote", os.path.join("reports/hla_popgen/57_yak_trioeval", name))


def fig_regions():
    s = pd.read_csv(f"{D}/trioeval_summary.tsv", sep="\t", dtype=str)
    tg = [t for t in TARGET_LABEL if t in set(s["target"])]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.9))
    for ax, k, title in ((axes[0], "sw", "a  switch error: adjacent hap-specific k-mers flip parent"),
                         (axes[1], "hm", "b  hamming error: hap-specific k-mers from the contig's minority parent")):
        for j, pl in enumerate(("revio", "sequel2e")):
            lab, col = PLAT[pl]
            for i, t in enumerate(tg):
                r = s[(s["target"] == t) & (s["child_platform"] == pl)]
                if r.empty:
                    continue
                r = r.iloc[0]
                x = i + (j - 0.5) * 0.3
                med, lo, hi, po = (num(r[f"{k}_{q}_pct"]) for q in ("median", "p25", "p75", "pooled"))
                if np.isnan(med):
                    continue
                m2, l2, h2 = max(med, FLOOR), max(lo, FLOOR), max(hi, FLOOR)
                ax.plot([x, x], [l2, h2], color=col, lw=1.4, solid_capstyle="butt")
                ax.plot(x, m2, "o", color=col, ms=6, mfc=col if med > 0 else "white", mew=1.5)
                ax.plot(x, po, "_", color=col, ms=13, mew=2.2)
        ax.set_xticks(range(len(tg)))
        ax.set_xticklabels([TARGET_LABEL[t] for t in tg], fontsize=8)
        ax.set_yscale("log")
        ax.set_ylim(FLOOR * 0.8, 60 if k == "hm" else 3)
        ax.set_title(title, fontsize=9, loc="left", color=INK)
        ax.set_ylabel("% (log scale)")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=BLUE, marker="o", lw=1.4, label="Revio child: median and IQR (25 trios)"),
         Line2D([], [], color=ORANGE, marker="o", lw=1.4, label="Sequel IIe child: median and IQR (25 trios)"),
         Line2D([], [], color=INK2, marker="_", ms=13, mew=2.2, lw=0, label="pooled over all k-mers"),
         Line2D([], [], color=INK2, marker="o", mfc="white", lw=0, label=f"median 0 (drawn at {FLOOR}%)")]
    fig.legend(handles=h, loc="lower center", ncol=4, fontsize=7.5, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("yak trioeval, 50 AoU trios: child hifiasm haplotypes scored against parental k-mers",
                 x=0.01, ha="left", fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    save(fig, "fig_trioeval_regions.png")


def fig_pairs():
    s = pd.read_csv(f"{D}/pair_phase_truth.tsv", sep="\t", dtype=str)
    s = s[s["child_platform"] == "all"].set_index("pair")
    ref = pd.read_csv(f"{D56}/refined_sweep.tsv", sep="\t", dtype=str)
    ref = ref[(ref["check"] == "hap_not_transmitted") & (ref["min_score"] == "0")].set_index("item")
    order = [x for x in ["DQA1~DQB1", "DRB1~DQA1", "DRB1~DQB1", "DPA1~DPB1", "B~C", "DQB1~DPB1", "A~C", "A~B"]
             if x in s.index]
    fig, ax = plt.subplots(figsize=(8.4, 3.6))
    x = np.arange(len(order))
    for xi, k in zip(x, order):
        r = s.loc[k]
        n = num(r["n_same_contig_haps"])
        if r["n_parent_mismatch"] == "<20":  # masked: draw the upper bound 19/n, hatched
            ub = 100 * 19 / n
            ax.bar(xi - 0.19, ub, 0.36, color="white", edgecolor=BLUE, hatch="////", lw=1)
            ax.text(xi - 0.19, ub + 0.8, f"<20/{int(n)}", ha="center", va="bottom", fontsize=7, color=INK2)
        else:
            pct, lo, hi = num(r["mismatch_pct"]), num(r["lo"]), num(r["hi"])
            ax.bar(xi - 0.19, pct, 0.36, color=BLUE)
            ax.errorbar(xi - 0.19, pct, yerr=[[pct - lo], [hi - pct]], fmt="none", ecolor=INK2, capsize=3, lw=1)
            ax.text(xi - 0.19, hi + 0.8, f"{r['n_parent_mismatch']}/{int(n)}", ha="center", va="bottom", fontsize=7,
                    color=INK2)
        if k in ref.index and not np.isnan(num(ref.loc[k, "pct"])):
            v, lo2, hi2 = num(ref.loc[k, "pct"]), num(ref.loc[k, "lo"]), num(ref.loc[k, "hi"])
            ax.bar(xi + 0.19, v, 0.36, color=MUTED)
            ax.errorbar(xi + 0.19, v, yerr=[[v - lo2], [hi2 - v]], fmt="none", ecolor=INK2, capsize=3, lw=1)
        elif k in ref.index:
            ax.text(xi + 0.19, 1, "<20", ha="center", va="bottom", fontsize=7, color=MUTED)
    ax.set_xticks(x)
    ax.set_xticklabels(order, fontsize=8)
    ax.set_ylabel("% same-contig haplotypes")
    ax.set_ylim(0, 72)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=BLUE, label="57b trio truth: the two genes come from different parents"),
                       Patch(facecolor="white", edgecolor=BLUE, hatch="////", label="57b, <20 events: upper bound 19/n"),
                       Patch(color=MUTED, label="56b relatives: one-contig haplotype not shared (indirect)")],
              fontsize=7, loc="upper left")
    ax.set_title("Physical phasing truth, 179 trios: does one contig carry one parent's haplotype?", fontsize=9,
                 loc="left", color=INK)
    fig.tight_layout()
    save(fig, "fig_pair_phase_truth.png")


def fig_mendel():
    s = pd.read_csv(f"{D}/mendel_typing.tsv", sep="\t", dtype=str)
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 3.4), gridspec_kw={"width_ratios": [1, 1.6]})
    for g, col in (("ALL", BLUE), ("classical8", GREEN)):
        d = s[s["gene"] == g]
        xs = num(d["min_score"]).to_numpy()
        for c, ls, lab in (("mismatch_2field_pct", "-", "2-field"), ("mismatch_full_pct", "--", "full resolution")):
            a.plot(xs, num(d[c]), ls, marker="o", color=col, lw=1.6, ms=5, label=f"{'all 15 genes' if g == 'ALL' else 'classical 8'}, {lab}")
        d2 = d.dropna(subset=["lo"])
        a.fill_between(num(d2["min_score"]), num(d2["lo"]), num(d2["hi"]), color=col, alpha=0.12, lw=0)
    a.set_xticks([0, 4, 6])
    a.set_xlabel("minimum per-call confidence score")
    a.set_ylabel("% child calls absent from\nthe transmitting parent")
    a.set_ylim(0, 14)
    a.legend(fontsize=6.8, loc="center left", bbox_to_anchor=(0.02, 0.5))
    a.set_title("a  Mendelian typing error vs score", fontsize=9, loc="left", color=INK)
    d = s[(s["min_score"] == "0") & ~s["gene"].isin(["ALL", "classical8"])].copy()
    d["full"] = num(d["mismatch_full_pct"])
    d = d.dropna(subset=["full"]).sort_values("full")
    b.barh(range(len(d)), d["full"], color=BLUE, height=0.65)
    b.set_yticks(range(len(d)))
    b.set_yticklabels(d["gene"], fontsize=8)
    for i, (v, n) in enumerate(zip(d["full"], d["n_calls"])):
        b.text(v + 0.3, i, f"{v:.1f}%  (n={n})", va="center", fontsize=7, color=INK2)
    b.set_xlim(0, 26)
    b.set_xlabel("% full-resolution allele absent from the transmitting parent (score ≥0)")
    b.set_title("b  per gene, full resolution (2-field per-gene counts <20 → not drawn)", fontsize=9, loc="left",
                color=INK)
    fig.suptitle("Typing truth, 179 trios: is the child's allele carried by the parent it came from?", x=0.01,
                 ha="left", fontsize=10, color=INK)
    fig.tight_layout()
    save(fig, "fig_mendel_typing.png")


if __name__ == "__main__":
    fig_regions()
    fig_pairs()
    fig_mendel()
