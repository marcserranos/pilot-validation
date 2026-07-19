#!/usr/bin/env python3
"""LOCAL-ONLY plotting script -- run on your Mac, not the VM. No pixi, no gcsfuse mount,
no AoU access needed: it only reads the two aggregate CSVs already downloaded from the VM.

Inputs (download via the JupyterLab file browser, see reports/aou_callset_validation.md):
  ~/pipeline_outputs/experiment_a2/allele_freq_by_ancestry.csv  (gene, ancestry, allele, freq, n_individuals)
  ~/pipeline_outputs/experiment_a3/allele_tree_data.csv         (gene, group_2field, catalogued_3field_variants,
                                                                   observed_3field_variants, n_individuals,
                                                                   dominant_ancestry, dominant_ancestry_share)

Setup (one time):
  pip3 install pandas matplotlib

Usage:
  python3 local_plot_allele_ancestry.py \
      --freq-csv ~/Downloads/allele_freq_by_ancestry.csv \
      --tree-csv ~/Downloads/allele_tree_data.csv \
      --gene DRB1 \
      --out-dir ~/Downloads/aou_plots

Produces 3 PNGs in --out-dir:
  1. <gene>_ancestry_stacked.png   -- top alleles for --gene, stacked bar = ancestry SHARE of
     carriers (not a single dominant-ancestry color -- the full composition).
  2. <gene>_top_per_ancestry.png   -- most common alleles WITHIN each of the 6 ancestry groups,
     small multiples, for --gene.
  3. allele_space_sunburst.png     -- two-level sunburst, ALL 8 genes: inner ring = gene (sized
     by total carriers), outer ring = each gene's top allele groups (sized by carriers, colored
     by dominant ancestry, remainder folded into "other").
"""
import argparse
import math
import os

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Wedge
import pandas as pd

ANCESTRY_COLORS = {
    "afr": "#2a78d6", "amr": "#008300", "eas": "#e87ba4",
    "eur": "#eda100", "mid": "#1baf7a", "sas": "#eb6834",
}
GENE_GRAY = "#888780"
ANCESTRY_ORDER = ["afr", "amr", "eas", "eur", "mid", "sas"]


def plot_ancestry_stacked(freq_df, gene, out_path, top_n=12):
    """Top alleles for `gene` by total carrier copies, stacked bar = full ancestry composition
    (share of that allele's copies coming from each ancestry group) -- not a single color."""
    g = freq_df[freq_df.gene == gene].copy()
    g["copies"] = g["freq"] * 2 * g["n_individuals"]
    totals = g.groupby("allele")["copies"].sum().sort_values(ascending=False)
    top_alleles = totals.head(top_n).index.tolist()

    pivot = g[g.allele.isin(top_alleles)].pivot_table(
        index="allele", columns="ancestry", values="copies", fill_value=0)
    pivot = pivot.reindex(top_alleles)
    for a in ANCESTRY_ORDER:
        if a not in pivot.columns:
            pivot[a] = 0.0
    pivot = pivot[ANCESTRY_ORDER]
    shares = pivot.div(pivot.sum(axis=1), axis=0) * 100

    fig, ax = plt.subplots(figsize=(9, 0.5 * len(top_alleles) + 2))
    left = pd.Series(0.0, index=shares.index)
    for anc in ANCESTRY_ORDER:
        ax.barh(shares.index, shares[anc], left=left, color=ANCESTRY_COLORS[anc], label=anc)
        left += shares[anc]
    ax.invert_yaxis()
    ax.set_xlabel("Share of carriers, by ancestry (%)")
    ax.set_title(f"{gene} -- top {len(top_alleles)} alleles, ancestry composition")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=6, frameon=False)
    ax.set_xlim(0, 100)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"wrote {out_path}")


def plot_top_per_ancestry(freq_df, gene, out_path, top_n=5):
    """Most common alleles WITHIN each ancestry group (not overall) -- small multiples."""
    g = freq_df[freq_df.gene == gene]
    fig, axes = plt.subplots(1, len(ANCESTRY_ORDER), figsize=(3 * len(ANCESTRY_ORDER), 4), sharex=True)
    for ax, anc in zip(axes, ANCESTRY_ORDER):
        sub = g[g.ancestry == anc].sort_values("freq", ascending=False).head(top_n)
        colors = [ANCESTRY_COLORS[anc]] * len(sub)
        ax.barh(sub["allele"], sub["freq"] * 100, color=colors)
        ax.invert_yaxis()
        ax.set_title(anc)
        ax.set_xlabel("Freq (%)")
    fig.suptitle(f"{gene} -- most common allele within each ancestry group")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"wrote {out_path}")


def _pt(deg, r):
    rad = math.radians(deg)
    return r * math.cos(rad), r * math.sin(rad)


def plot_sunburst(tree_df, out_path, top_n_per_gene=6):
    """Two-level sunburst: inner ring = genes (sized by total carriers), outer ring = each
    gene's top allele groups (sized by carriers, colored by dominant ancestry, 'other' folds
    the remainder). Pure matplotlib (Wedge patches) -- no extra dependencies."""
    genes = sorted(tree_df["gene"].unique())
    gene_totals = tree_df.groupby("gene")["n_individuals"].sum()
    total = gene_totals.sum()

    fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=dict(aspect="equal"))
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.axis("off")

    inner_r, outer_r_in, outer_r_out = 0.55, 0.6, 1.0
    angle = 0.0
    for gene in genes:
        gtotal = gene_totals[gene]
        gspan = 360.0 * gtotal / total

        ax.add_patch(Wedge((0, 0), inner_r, angle, angle + gspan,
                           facecolor=GENE_GRAY, edgecolor="white", linewidth=1))
        mid = angle + gspan / 2
        x, y = _pt(mid, inner_r * 0.6)
        ax.text(x, y, gene, ha="center", va="center", fontsize=10, fontweight="bold")

        sub = tree_df[tree_df.gene == gene].sort_values("n_individuals", ascending=False)
        top = sub.head(top_n_per_gene)
        other_n = sub["n_individuals"].iloc[top_n_per_gene:].sum()
        segs = list(zip(top["group_2field"], top["n_individuals"], top["dominant_ancestry"]))
        if other_n > 0:
            segs.append(("other", other_n, None))

        a2 = angle
        for label, n, dom_anc in segs:
            span2 = gspan * n / gtotal if gtotal else 0
            color = ANCESTRY_COLORS.get(dom_anc, "#c3c2b7")
            ax.add_patch(Wedge((0, 0), outer_r_out, a2, a2 + span2,
                               width=outer_r_out - outer_r_in,
                               facecolor=color, edgecolor="white", linewidth=0.5))
            if span2 > 6:  # only label wedges wide enough to read
                mx, my = _pt(a2 + span2 / 2, (outer_r_in + outer_r_out) / 2)
                ax.text(mx, my, label, ha="center", va="center", fontsize=7)
            a2 += span2

        angle += gspan

    handles = [mpatches.Patch(color=c, label=a) for a, c in ANCESTRY_COLORS.items()]
    handles.append(mpatches.Patch(color=GENE_GRAY, label="gene ring"))
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(-0.05, -0.05),
             fontsize=9, title="Outer ring = dominant ancestry", frameon=False)
    ax.set_title("AoU-native HLA calls -- allele space by gene and ancestry", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--freq-csv", required=True)
    ap.add_argument("--tree-csv", required=True)
    ap.add_argument("--gene", default="DRB1", help="gene for the two per-gene plots")
    ap.add_argument("--out-dir", default=os.path.expanduser("~/Downloads/aou_plots"))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    freq_df = pd.read_csv(os.path.expanduser(args.freq_csv))
    tree_df = pd.read_csv(os.path.expanduser(args.tree_csv))

    plot_ancestry_stacked(freq_df, args.gene,
                          os.path.join(args.out_dir, f"{args.gene}_ancestry_stacked.png"))
    plot_top_per_ancestry(freq_df, args.gene,
                          os.path.join(args.out_dir, f"{args.gene}_top_per_ancestry.png"))
    plot_sunburst(tree_df, os.path.join(args.out_dir, "allele_space_sunburst.png"))


if __name__ == "__main__":
    main()
