#!/usr/bin/env python3
"""Panel-level overview of the per-gene diversity results from `15_hla_manhattan.py`.

The per-gene Manhattan plots show *where* diversity sits along one gene. This shows the comparison
ACROSS genes, which is where the actual argument lives:

  Top panel  -- absolute mean pi inside vs outside the CDS, log scale. This is the control result:
                coding diversity spans ~55x across the panel, ordered exactly as known biology
                (classical class I highest, HLA-DRA / E / F lowest). A method returning uniform
                diversity would be measuring noise; this one reproduces a gradient it was never
                tuned on.

  Bottom     -- the CDS / non-CDS ratio, diverging around 1.0. Ratio < 1 is the DEFAULT for an
                ordinary gene: purifying selection removes coding variation, pushing pi_CDS below
                the roughly-neutral intronic background. Ratio > 1 is the anomaly, and it is the
                signature of balancing selection maintaining coding variation. So the interesting
                claim is not "CDS is diverse" but "CDS is diverse ONLY in the genes where balancing
                selection is expected".

Genes flagged with an asterisk carry a data-quality caveat recorded in the report and are drawn
hatched so they are not read as clean support:
  HLA-E    -- ratio rests on tiny absolute values and an intronic pi ~30x below its paralogs F/G
  HLA-DPB1 -- NM max 1608 vs 51-412 elsewhere, i.e. some haplotypes align poorly to its canonical

Usage:
    python3 scripts/hla_popgen/16_panel_overview.py \\
        --summary reports/hla_popgen/15_hla_manhattan/panel_summary.tsv \\
        --out reports/hla_popgen/15_hla_manhattan/panel_overview.png
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Gene classes drive the grouping and colour, and are the whole point of the control arm.
GENE_CLASS = {
    "HLA-A": "classical I", "HLA-B": "classical I", "HLA-C": "classical I",
    "HLA-DRB1": "classical II", "HLA-DQB1": "classical II", "HLA-DPB1": "classical II",
    "HLA-DQA1": "classical II", "HLA-DPA1": "classical II",
    "HLA-DRA": "conserved", "HLA-E": "conserved", "HLA-F": "conserved", "HLA-G": "conserved",
}
CLASS_COLOR = {
    "classical I": "#B3341F",
    "classical II": "#E08A2E",
    "conserved": "#3E6E8E",
}
CAVEATED = {"HLA-E", "HLA-DPB1"}


def load(summary_path):
    rows = []
    with open(summary_path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append({
                "gene": r["gene"],
                "pi_cds": float(r["mean_pi_cds"]),
                "pi_non": float(r["mean_pi_noncds"]),
                "ratio": float(r["pi_ratio_cds_over_noncds"]),
                "n_haps": int(r["n_haps"]),
            })
    rows.sort(key=lambda d: -d["pi_cds"])
    return rows


def plot(rows, out_path):
    genes = [r["gene"] for r in rows]
    short = [g.replace("HLA-", "") for g in genes]
    classes = [GENE_CLASS.get(g, "conserved") for g in genes]
    colors = [CLASS_COLOR[c] for c in classes]
    x = np.arange(len(genes))

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(13, 8.2), dpi=180, sharex=True,
        gridspec_kw={"height_ratios": [1.35, 1.0], "hspace": 0.13})
    fig.patch.set_facecolor("white")

    # ---- Top: absolute pi, CDS vs non-CDS, log scale ----
    w = 0.38
    cds_vals = [r["pi_cds"] for r in rows]
    non_vals = [r["pi_non"] for r in rows]
    b1 = ax1.bar(x - w / 2, cds_vals, w, color=colors, edgecolor="black", linewidth=0.5,
                 label="inside CDS")
    b2 = ax1.bar(x + w / 2, non_vals, w, color="#D8D5DE", edgecolor="black", linewidth=0.5,
                 label="outside CDS (intron / UTR)")
    for bar_cds, bar_non, g in zip(b1, b2, genes):
        if g in CAVEATED:
            bar_cds.set_hatch("///")
            bar_non.set_hatch("///")

    ax1.set_yscale("log")
    ax1.set_ylabel("mean $\\pi$ per site  (log scale)", fontsize=11)
    ax1.grid(axis="y", color="#E8E8E8", lw=0.8)
    ax1.set_axisbelow(True)
    ax1.spines[["top", "right"]].set_visible(False)
    ax1.legend(frameon=False, fontsize=9.5, loc="upper right", ncol=2)

    fold = max(cds_vals) / min(cds_vals)
    ax1.set_title(
        f"Coding diversity spans {fold:.0f}× across the panel, ordered as known biology\n"
        f"(classical class I highest → HLA-DRA / E / F lowest; ~24,000 haplotypes per gene)",
        fontsize=12, fontweight="bold", pad=10)

    # ---- Bottom: the CDS/non-CDS ratio, diverging about 1.0 ----
    # Colour by gene CLASS, not by whether the ratio clears 1. Colouring by ratio would contradict
    # the class key below the figure (DPB1 is class II but sits above 1), and the axhline already
    # carries the above/below-1 signal.
    ratios = [r["ratio"] for r in rows]
    b3 = ax2.bar(x, ratios, 0.62, color=colors, edgecolor="black", linewidth=0.5)
    for bar, g in zip(b3, genes):
        if g in CAVEATED:
            bar.set_hatch("///")

    ax2.axhline(1.0, color="black", lw=1.4, zorder=5)
    # Park the annotation over the sub-1 genes: everything above the line is empty there, whereas
    # both ends of the axis are occupied (enriched bars on the left, HLA-E's 3.69 on the right).
    first_below = next((i for i, v in enumerate(ratios) if v < 1), len(genes) - 1)
    ax2.text(first_below - 0.3, 1.06, "ratio = 1  (no enrichment)", ha="left", va="bottom",
             fontsize=9, color="black")
    for xi, v in zip(x, ratios):
        ax2.text(xi, v + (0.09 if v > 1 else 0.06), f"{v:.2f}", ha="center", va="bottom",
                 fontsize=9)

    ax2.set_ylabel("$\\pi_{CDS}$ / $\\pi_{non-CDS}$", fontsize=11)
    ax2.set_ylim(0, max(ratios) * 1.22)
    ax2.grid(axis="y", color="#E8E8E8", lw=0.8)
    ax2.set_axisbelow(True)
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.set_xticks(x)
    ax2.set_xticklabels(short, fontsize=11)
    ax2.set_title(
        "Enrichment appears ONLY where balancing selection is expected — "
        "ratio < 1 is the normal purifying-selection case",
        fontsize=11.5, pad=8)

    # Two legends: the CDS/non-CDS pairing lives on ax1, the colour key is a FIGURE legend below
    # everything. Putting the colour key on ax1 with a negative bbox draws it into the gap where
    # ax2's title sits, where it is silently overdrawn.
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=CLASS_COLOR[c], edgecolor="black", lw=0.5)
               for c in ["classical I", "classical II", "conserved"]]
    labels = ["classical class I", "classical class II", "conserved / non-classical"]
    handles.append(plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="black", lw=0.5,
                                 hatch="///"))
    labels.append("data-quality caveat (see report)")
    fig.legend(handles=handles, labels=labels, frameon=False, fontsize=9.5,
               loc="lower center", bbox_to_anchor=(0.5, -0.015), ncol=4)

    fig.savefig(out_path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", out_path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument("--summary",
                     default=os.path.join(here, "reports/hla_popgen/15_hla_manhattan/panel_summary.tsv"))
    ap.add_argument("--out",
                     default=os.path.join(here, "reports/hla_popgen/15_hla_manhattan/panel_overview.png"))
    args = ap.parse_args()
    rows = load(args.summary)
    plot(rows, args.out)


if __name__ == "__main__":
    main()
