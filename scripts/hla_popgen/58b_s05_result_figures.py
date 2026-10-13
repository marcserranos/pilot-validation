#!/usr/bin/env python3
"""Figures for the S05 VM results (58 regional depth, 64 novel-allele validation), drawn from the
committed PUBLIC tables only (no VM, no participant data). Masked cells (<20, <20*) are drawn as
hatched placeholders or left out, never as numbers.

    python3 scripts/hla_popgen/58b_s05_result_figures.py

Writes PNGs next to the tables:
  reports/hla_popgen/58_regional_depth/fig_window_depth.png
  reports/hla_popgen/58_regional_depth/fig_trb_locus_profile.png
  reports/hla_popgen/58_regional_depth/fig_trb_histograms.png
  reports/hla_popgen/64_novel_allele_validation/fig_funnel_read_evidence.png
  reports/hla_popgen/64_novel_allele_validation/fig_hp_runlength.png
  reports/hla_popgen/64_novel_allele_validation/fig_call_outcomes.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D58 = os.path.join(REPO, "reports", "hla_popgen", "58_regional_depth")
D64 = os.path.join(REPO, "reports", "hla_popgen", "64_novel_allele_validation")
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7", "#008300", "#e34948"]
BLUE, ORANGE, GREEN, YELLOW, PINK, PURPLE, _, RED = SERIES
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e6e5e0"
ANC = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "axes.axisbelow": True, "figure.dpi": 150, "savefig.bbox": "tight", "legend.frameon": False,
    "font.family": "DejaVu Sans"})


def num(s):
    return pd.to_numeric(s, errors="coerce")


def save(fig, d, name):
    fig.savefig(os.path.join(d, name))
    plt.close(fig)
    print("wrote", os.path.relpath(os.path.join(d, name), REPO))


# ---------------------------------------------------------------- 58

def fig_window_depth():
    s = pd.read_csv(f"{D58}/depth_summary_by_region.tsv", sep="\t", dtype=str)
    s = s[(s["mapq"] == "0") & s["region"].str.startswith("window_")]
    strata = [("ALL", "ALL")] + [(a, a) for a in ANC] + [("revio", "Revio"), ("sequel2e", "Sequel IIe")]
    wins = [("window_MHC", "MHC  chr6:29.5-33.5 Mb", BLUE), ("window_KIR", "KIR region  chr19", GREEN),
            ("window_TRB", "TRB locus  chr7", PURPLE)]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, (reg, title, col) in zip(axes, wins):
        for i, (key, lab) in enumerate(strata):
            r = s[(s["region"] == reg) & (s["stratum"] == key)]
            if r.empty:
                continue
            r = r.iloc[0]
            p5, p25, p50, p75, p95 = (float(r[f"rel_p{q}"]) for q in ("05", "25", "50", "75", "95"))
            c = col if key not in ("revio", "sequel2e") else (BLUE if key == "revio" else ORANGE)
            ax.plot([i, i], [p5, p95], color=c, lw=1)
            ax.add_patch(plt.Rectangle((i - 0.28, p25), 0.56, p75 - p25, facecolor=c, alpha=0.25,
                                       edgecolor=c, lw=0.8))
            ax.plot([i - 0.28, i + 0.28], [p50, p50], color=c, lw=2)
            if key == "ALL":
                ax.text(i, p95 + 0.012, f"{p50:.3f}", ha="center", va="bottom", fontsize=8, color=INK)
        ax.axhline(1.0, color=INK2, lw=0.8, ls="--")
        ax.axvline(0.5, color=GRID, lw=0.8)
        ax.axvline(len(ANC) + 0.5, color=GRID, lw=0.8)
        ax.set_xticks(range(len(strata)))
        ax.set_xticklabels([lab for _, lab in strata], rotation=45, ha="right")
        ax.set_title(title)
        ax.set_xlim(-0.6, len(strata) - 0.4)
    axes[0].set_ylabel("depth relative to the person's\ngenome-wide depth (MAPQ >= 0)")
    axes[0].set_ylim(0.6, 1.15)
    fig.text(0.01, -0.12, "box = 25th-75th percentile, thick line = median, whisker = 5th-95th percentile, "
             "12,233 people. Dashed line = covered like the rest of the genome.", fontsize=8, color=INK2)
    fig.suptitle("Long-read depth over the three immune regions, relative to each person's genome", y=1.02)
    save(fig, D58, "fig_window_depth.png")


def _gene_class(g):
    if g in ("TRBV4-2", "TRBV6-2"):
        return "ins", GREEN
    if g in ("TRBV5-7", "TRBV6-5", "TRBV6-8", "TRBV7-4", "TRBV7-5", "TRBV7-7", "TRBV6-6"):
        return "high", PURPLE
    if g.startswith(("TRBD", "TRBJ", "TRBC")):
        return "djc", ORANGE
    return "v", BLUE


def fig_trb_profile():
    s = pd.read_csv(f"{D58}/depth_summary_by_region.tsv", sep="\t", dtype=str)
    s = s[(s["mapq"] == "0") & (s["region_class"] == "gene_TRB") & (s["stratum"] == "ALL")]
    bed = pd.read_csv(f"{D58}/regions.bed", sep="\t", header=None,
                      names=["chrom", "start", "end", "name", "cls"], dtype=str)
    pos = {r["name"]: (int(r["start"]) + int(r["end"])) / 2e6 for _, r in bed.iterrows()}
    s = s[s["region"].isin(pos)].copy()
    s["x"] = s["region"].map(pos)
    for q in ("25", "50", "75"):
        s[q] = num(s[f"rel_p{q}"])
    s = s.sort_values("x")
    fig, ax = plt.subplots(figsize=(11, 4.2))
    for _, r in s.iterrows():
        _, c = _gene_class(r["region"])
        ax.plot([r["x"], r["x"]], [r["25"], r["75"]], color=c, lw=1, alpha=0.6)
        ax.plot(r["x"], r["50"], "o", color=c, ms=4)
    for g, dy, ha in (("TRBV4-2", 0.08, "right"), ("TRBV6-2", 0.08, "left"), ("TRBJ1-1", -0.09, "center"),
                      ("TRBV5-7", 0.07, "center"), ("TRBV30", 0.07, "center"), ("TRBV2", 0.07, "center"),
                      ("TRBC2", 0.06, "left")):
        r = s[s["region"] == g]
        if not r.empty:
            r = r.iloc[0]
            ax.annotate(g, (r["x"], r["50"]), xytext=(r["x"], r["50"] + dy), ha=ha, fontsize=7.5,
                        color=INK2, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5))
    ax.axhline(1.0, color=INK2, lw=0.8, ls="--")
    ax.set_xlabel("position on chr7 (Mb, GRCh38 primary)")
    ax.set_ylabel("relative depth (median; line = IQR)")
    ax.set_ylim(0.55, 1.75)
    ax.legend(handles=[
        Patch(color=BLUE, label="TRBV genes (typical)"),
        Patch(color=GREEN, label="TRBV4-2, TRBV6-2: paralogs of the 21.6 kb insertion (extra copies map here)"),
        Patch(color=PURPLE, label="other TRBV genes running high (not explained yet)"),
        Patch(color=ORANGE, label="TRBD, TRBJ, TRBC: depth falls toward D1-J1 (pattern of T-cell V(D)J deletions)")],
        loc="upper left", fontsize=7.5)
    ax.set_title("TRB locus: median relative depth per gene across 12,233 people")
    save(fig, D58, "fig_trb_locus_profile.png")


def fig_trb_histograms():
    h = pd.read_csv(f"{D58}/trb_rel_depth_histogram.tsv", sep="\t", dtype=str)
    panels = [("TRBV3-1", "a  TRBV3-1: a typical gene", BLUE),
              ("TRBJ1-1", "b  TRBJ1-1: shifted down,\none broad mode", ORANGE),
              ("TRBV4-2", "c  TRBV4-2: insertion paralog,\npile-up in the open last bin", GREEN)]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.3), sharey=True)
    for ax, (g, title, col) in zip(axes, panels):
        r = h[h["gene"] == g].copy()
        r["lo"] = num(r["bin_lo"])
        masked = r["n_people"].str.startswith("<")
        r["n"] = num(r["n_people"]).fillna(0)
        ax.bar(r["lo"][~masked], r["n"][~masked], width=0.05, align="edge", color=col, alpha=0.8)
        ax.bar(r["lo"][masked], 10, width=0.05, align="edge", facecolor="none", edgecolor=MUTED,
               hatch="////", lw=0.4)
        for cut in (0.25, 0.75):
            ax.axvline(cut, color=RED, lw=0.8, ls=":")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("relative depth (0.05 bins; last is >= 1.45)")
        ax.set_xlim(0, 1.5)
    axes[0].set_ylabel("people")
    axes[0].text(0.26, axes[0].get_ylim()[1] * 0.92, "copy-number cuts\n0.25 / 0.75", color=RED, fontsize=7)
    fig.text(0.01, -0.08, "Hatched = 1-19 people (masked). Expected modes for a germline deletion: 0, 0.5 "
             "(one copy) and 1.0 (two copies). None of the genes shows a separate 0.5 mode.",
             fontsize=8, color=INK2)
    save(fig, D58, "fig_trb_histograms.png")


# ---------------------------------------------------------------- 64

STEP_LABEL = {
    "distinct": "distinct protein-novel sequences",
    "ge2_unrelated_carriers": ">= 2 unrelated carriers",
    "not_homopolymer_only": "not a homopolymer-only indel",
    "no_template_warning": "no template warning",
    "ge2_platforms": "seen on >= 2 platforms",
    "ge2_ancestries": "seen in >= 2 ancestries",
    "ge20_carriers_clean": ">= 20 clean carriers"}


def fig_funnel():
    f = pd.read_csv(f"{D64}/tables/funnel_with_evidence.tsv", sep="\t", dtype=str)
    f = f[(f["novelty_class"] == "protein") & (f["gene"] == "ALL")]
    f = f.set_index("step").loc[list(STEP_LABEL)]
    cats = [("n_pooled_supported", "reads support the novel allele", GREEN),
            ("n_pooled_contradicted", "reads carry the known allele instead", RED),
            ("n_pooled_ambiguous", "reads undecided", YELLOW),
            ("n_pooled_untested", "no pooled read verdict", "#c9c7c0")]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 3.6), gridspec_kw=dict(width_ratios=[1.25, 1]),
                                 sharey=True)
    y = np.arange(len(f))[::-1]
    for yi, (step, r) in zip(y, f.iterrows()):
        tot = int(r["n_keys"])
        if tot:
            ax.barh(yi, tot, color=BLUE, height=0.6)
        nk = int(r["n_now_known"])
        ax.text(max(tot, 0.6) * 1.25, yi, f"{tot:,}" + (f"  ({nk} now in IPD-IMGT/HLA)" if nk else ""),
                va="center", fontsize=8, color=INK)
        left = 0.0
        for col, _, c in cats:
            v = int(r[col])
            if tot and v:
                bx.barh(yi, 100 * v / tot, left=left, color=c, height=0.6)
                if 100 * v / tot >= 9:
                    bx.text(left + 50 * v / tot, yi, str(v), ha="center", va="center", fontsize=7.5,
                            color="white" if c != "#c9c7c0" else INK)
                left += 100 * v / tot
    ax.set_yticks(y)
    ax.set_yticklabels([STEP_LABEL[s] for s in f.index])
    ax.set_xscale("log")
    ax.set_xlim(0.6, 3e5)
    ax.set_xticks([1, 10, 100, 1000])
    ax.set_xticklabels(["1", "10", "100", "1,000"])
    ax.set_xlabel("distinct sequences (keys), log scale")
    ax.set_title("a  keys left after each filter", loc="left", fontsize=9)
    bx.set_xlim(0, 100)
    bx.set_xlabel("% of the step's keys, by pooled long-read verdict")
    bx.set_title("b  what the long reads say about them", loc="left", fontsize=9)
    bx.legend(handles=[Patch(color=c, label=l) for _, l, c in cats], loc="upper center",
              bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=7.5)
    fig.suptitle("Classical protein-novel HLA alleles (8 genes, unrelated people)", y=1.02)
    save(fig, D64, "fig_funnel_read_evidence.png")


def fig_hp_runlength():
    t = pd.read_csv(f"{D64}/tables/hp_runlength_by_gene.tsv", sep="\t", dtype=str)
    t = t.iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 3.2))
    y = np.arange(len(t))
    left = np.zeros(len(t))
    for col, lab, c in (("pct_called_len", "run length the assembly called (novel)", GREEN),
                        ("pct_nearest_len", "run length of the nearest known allele", RED),
                        ("pct_other_len", "another length", "#c9c7c0")):
        v = num(t[col]).to_numpy()
        ax.barh(y, v, left=left, color=c, height=0.6, label=lab)
        for yi, (l0, vi) in enumerate(zip(left, v)):
            if vi >= 8:
                ax.text(l0 + vi / 2, yi, f"{vi:.0f}%", ha="center", va="center", fontsize=7.5, color="white")
        left += v
    ax.set_yticks(y)
    ax.set_yticklabels([f"{g}  ({int(n):,} reads)" for g, n in zip(t["gene"], t["n_read_observations"])])
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of HiFi reads spanning the homopolymer run")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=7.5)
    ax.set_title("At homopolymer-indel novel calls, most reads carry the known run length")
    save(fig, D64, "fig_hp_runlength.png")


def fig_call_outcomes():
    t = pd.read_csv(f"{D64}/tables/failures_summary.tsv", sep="\t", dtype=str)
    get = {(r["stage"], r["reason"]): r["n_calls"] for _, r in t.iterrows()}
    bars = [("tested calls\n(4,616)", "call_status", [
                ("ok", "CDS rebuilt = called", BLUE), ("cds_disagree", "Immuannot CDS record differs", ORANGE),
                ("nearest_too_distant", "> 40 differences", PURPLE),
                ("indel_spans_exon_junction", "indel across exon junction", PINK)]),
            ("long-read test of\nthe 4,066 ok calls", "long_read", [
                ("ran", "clear verdict", GREEN), ("mixed", "reads mixed", YELLOW),
                ("low_reads", "too few reads", "#c9c7c0"), ("no_window_spanned", "no read spans", MUTED)])]
    fig, ax = plt.subplots(figsize=(9, 2.6))
    handles = []
    for yi, (lab, stage, parts) in enumerate(bars[::-1]):
        left = 0
        for key, plab, c in parts:
            raw = get.get((stage, key), "0")
            masked = raw.startswith("<")
            v = 10 if masked else int(raw)
            ax.barh(yi, v, left=left, color=c, height=0.55, hatch="////" if masked else None)
            if v > 250:
                ax.text(left + v / 2, yi, f"{v:,}", ha="center", va="center", fontsize=8, color="white")
            left += v
            handles.append(Patch(color=c, label=plab + (f" ({raw}, masked)" if masked else "")))
    ax.set_yticks([0, 1])
    ax.set_yticklabels([b[0] for b in bars[::-1]])
    ax.set_xlabel("novel-allele calls")
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=4, fontsize=7.5)
    ax.set_title("What happened to each tested novel call")
    save(fig, D64, "fig_call_outcomes.png")


if __name__ == "__main__":
    fig_window_depth()
    fig_trb_profile()
    fig_trb_histograms()
    fig_funnel()
    fig_hp_runlength()
    fig_call_outcomes()
