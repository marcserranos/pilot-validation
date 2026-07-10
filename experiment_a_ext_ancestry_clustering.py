#!/usr/bin/env python3
"""Experiment A extension: locus x ancestry homozygosity heatmap (all 8 genes,
not just DQA1) + a PCA scatter using AoU's own already-computed pca_features,
colored by genetic ancestry and by per-person cross-locus homozygosity count.
See EXPERIMENTS.md 2026-07-10 for design. Writes PNGs to disk (view via the
Jupyter file browser on the VM) -- prints only aggregate stats to stdout, no
per-person data, and nothing here should be downloaded off the Workbench."""
import argparse
import ast
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CLASSICAL_GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
ANCESTRY_ORDER = ["afr", "amr", "eas", "eur", "mid", "sas", "oth"]


def locus_ancestry_homozygosity(df, ancestry_col):
    rows = []
    for gene in CLASSICAL_GENES:
        c1, c2 = f"{gene}_1", f"{gene}_2"
        if c1 not in df.columns or c2 not in df.columns:
            continue
        for anc, sub in df.groupby(ancestry_col):
            resolved = sub.dropna(subset=[c1, c2])
            if len(resolved) == 0:
                continue
            hz = (resolved[c1] == resolved[c2]).sum() / len(resolved) * 100
            rows.append({"gene": gene, "ancestry": anc, "n": len(sub), "homozygosity_pct": round(hz, 2)})
    return pd.DataFrame(rows)


def print_heatmap_table(hz_df):
    pivot = hz_df.pivot(index="gene", columns="ancestry", values="homozygosity_pct")
    cols = [a for a in ANCESTRY_ORDER if a in pivot.columns]
    pivot = pivot.reindex(index=CLASSICAL_GENES, columns=cols)
    print("## Locus x ancestry homozygosity (%)\n")
    print(pivot.to_markdown())
    print()
    return pivot


def save_heatmap_png(pivot, outpath):
    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(pivot.values, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.values[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=8)
    ax.set_title("Homozygosity rate (%) by gene x genetic ancestry")
    fig.colorbar(im, ax=ax, label="% homozygous")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def build_pca_scatter(df, ancestry_col, pca_col, sample_n, seed, outpath):
    n_total = len(df)
    sample = df.sample(n=sample_n, random_state=seed).copy() if n_total > sample_n else df.copy()
    n_sampled = len(sample)

    def parse_pcs(v):
        if pd.isna(v):
            return None
        try:
            parsed = ast.literal_eval(str(v))
            return parsed[0], parsed[1]
        except Exception:
            return None

    parsed = sample[pca_col].map(parse_pcs)
    valid = parsed.notna()
    n_failed = int((~valid).sum())
    sample = sample.loc[valid].copy()
    sample["pc1"] = parsed.loc[valid].map(lambda t: t[0])
    sample["pc2"] = parsed.loc[valid].map(lambda t: t[1])

    hz_count = pd.Series(0, index=sample.index)
    for gene in CLASSICAL_GENES:
        c1, c2 = f"{gene}_1", f"{gene}_2"
        if c1 in sample.columns and c2 in sample.columns:
            hz_count += ((sample[c1].notna()) & (sample[c2].notna()) & (sample[c1] == sample[c2])).astype(int)
    sample["hz_count"] = hz_count

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    ancestries = [a for a in ANCESTRY_ORDER if a in sample[ancestry_col].unique()]
    cmap = plt.get_cmap("tab10")
    for i, anc in enumerate(ancestries):
        pts = sample[sample[ancestry_col] == anc]
        axes[0].scatter(pts["pc1"], pts["pc2"], s=4, alpha=0.35, color=cmap(i), label=f"{anc} (n={len(pts)})")
    axes[0].set_title("PC1 vs PC2, colored by genetic ancestry")
    axes[0].set_xlabel("PC1")
    axes[0].set_ylabel("PC2")
    axes[0].legend(markerscale=3, fontsize=8, loc="best")

    sc = axes[1].scatter(sample["pc1"], sample["pc2"], s=4, alpha=0.5, c=sample["hz_count"], cmap="viridis")
    axes[1].set_title("PC1 vs PC2, colored by cross-locus homozygosity count (of 8)")
    axes[1].set_xlabel("PC1")
    axes[1].set_ylabel("PC2")
    fig.colorbar(sc, ax=axes[1], label="homozygous loci (of 8)")

    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

    return {"n_total": n_total, "n_sampled": n_sampled, "n_parse_failed": n_failed, "n_plotted": len(sample)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tsv", required=True, help="Path to hla_genotypes.tsv (mounted path)")
    ap.add_argument("--ancestry-tsv", required=True, help="Path to ancestry_preds.tsv (mounted path)")
    ap.add_argument("--id-col", default="research_id")
    ap.add_argument("--outdir", default=os.path.expanduser("~/pipeline_outputs/experiment_a_ext"))
    ap.add_argument("--sample-n", type=int, default=50000,
                     help="Subsample size for the PCA scatter -- 535k points renders as an unreadable blob; "
                          "a random subsample of this size keeps real proportions across ancestry groups "
                          "while staying fast to render. Raise it if smaller groups (e.g. mid, n~2100) look "
                          "too sparse.")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.tsv, sep="\t")
    anc = pd.read_csv(args.ancestry_tsv, sep="\t")

    anc_id_col = next((c for c in anc.columns if c.lower() in ("research_id", "person_id", "s")), None)
    ancestry_col = next((c for c in anc.columns if "ancestry" in c.lower() and "pca" not in c.lower()), None)
    pca_col = next((c for c in anc.columns if "pca" in c.lower()), None)
    if anc_id_col is None or ancestry_col is None:
        print(f"Could not identify id/ancestry columns. Ancestry TSV columns: {list(anc.columns)}", file=sys.stderr)
        sys.exit(1)

    keep_cols = [c for c in [anc_id_col, ancestry_col, pca_col] if c is not None]
    merged = df.merge(anc[keep_cols], left_on=args.id_col, right_on=anc_id_col, how="inner")
    print(f"Joined {len(merged)} / {len(df)} rows to an ancestry label.\n")

    hz_df = locus_ancestry_homozygosity(merged, ancestry_col)
    pivot = print_heatmap_table(hz_df)
    heatmap_path = os.path.join(args.outdir, "locus_ancestry_homozygosity_heatmap.png")
    save_heatmap_png(pivot, heatmap_path)
    print(f"Heatmap saved to {heatmap_path}\n")

    if pca_col is None:
        print("No pca_features-like column found in the ancestry TSV -- skipping PCA scatter.")
        print(f"Ancestry TSV columns were: {list(anc.columns)}")
        return

    scatter_path = os.path.join(args.outdir, "pca_ancestry_homozygosity_scatter.png")
    stats = build_pca_scatter(merged, ancestry_col, pca_col, args.sample_n, args.seed, scatter_path)
    print(f"PCA scatter saved to {scatter_path}")
    print(f"  total joined rows: {stats['n_total']}, subsampled to: {stats['n_sampled']}, "
          f"failed to parse pca_features: {stats['n_parse_failed']}, plotted: {stats['n_plotted']}")


if __name__ == "__main__":
    main()
