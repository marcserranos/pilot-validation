#!/usr/bin/env python3
"""S06 72: paper-grade figures for the call-quality sprint, drawn from the AGGREGATE tables of 70, 71, 71b and 71c.

Nothing here reads a person-level table. Resolution is the point (Marc, 2026-10-20): histograms in 0.01 bins, rolling error
curves, 0.25% pruning steps, 40 x 40 percentile grids; no deciles, no sparse matrices. 600 dpi PNG plus vector PDF, Nature
widths, the repo's layout linter (`save_fig` refuses overlapping or clipped text).

  fig_depth_distributions   per-gene depth of exons / CDS / peptide-binding exons, per-person MHC depth, gene x gene correlation
  fig_gene_depth_profile    depth along each classical gene (population median and 10-90th percentile) under the gene model
  fig_exon_depth            depth of each coding exon of each gene (ARS exons marked)
  fig_error_curves_<view>   typing error along every quality metric (rolling windows, Wilson 95% bands)
  fig_auc                   AUC of each metric for typing error, with person-clustered bootstrap CIs; within-gene AUC
  fig_metric_correlation    metric x metric Spearman correlation, clustered
  fig_pruning               remove the worst x% of calls by a metric: error among the rest, and the share of errors removed
  fig_grid_<a>__<b>         error on a 40 x 40 grid of two metrics' cohort percentile ranks
  fig_contig_truth          MHC gene-pair matrices: contig co-location, same-parent truth over all haplotypes, same-parent
                            truth restricted to pairs on one contig

Usage (laptop, after pulling the VM aggregate tables):
    python3 72_call_quality_figures.py --results ~/WORK/.../s06_results --out reports/hla_popgen/72_call_quality
"""
import argparse
import os
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _viz_common as V  # noqa: E402

REPO = os.path.dirname(os.path.dirname(HERE))
CLASSICAL = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
FAM_COLOR = {"depth": "#3C5488", "immuannot": "#E64B35", "assembly": "#00A087", "legacy": "#7E6148", "model": "#1B1919"}
VIEW_LABEL = {"trio_child": "Trios: child's call", "trio_parent": "Trios: transmitting parent's call",
              "trio_worst": "Trios: worse of child and parent (primary)", "pair": "First-degree pairs (worse call)",
              "twin": "Twins / duplicates (worse call)"}
METRIC_LABEL = {
    "template_distance": "Template distance (gene span)", "template_distance_per_kb": "Template distance per kb",
    "cds_distance": "CDS distance", "n_aa_changes": "Amino-acid changes", "has_warning": "Template warning",
    "n_tied": "Alleles tied", "is_novel": "Novel allele", "rel_gene": "Gene depth", "rel_exons": "Exon depth",
    "rel_cds": "CDS depth", "rel_ars": "Peptide-binding exon depth", "rel_exon_min": "Shallowest exon depth",
    "cds_p10_rel": "CDS 10th-percentile depth", "cds_min_rel": "CDS minimum depth",
    "ars_min_rel": "Peptide-binding minimum depth", "cds_share_lt10": "CDS bases under 10x",
    "ars_share_lt10": "Peptide-binding bases under 10x", "mapq_frac": "Unique-mapping fraction",
    "rel_cds_n": "CDS depth (gene-normalised)", "rel_ars_n": "Peptide-binding depth (gene-normalised)",
    "rel_exon_min_n": "Shallowest exon (exon-normalised)", "rel_mhc": "MHC window depth", "mhc_contigs": "MHC contigs", "diff_contig": "Off partner's contig",
    "score7": "Old 0-7 score", "model:immuannot": "Model: Immuannot", "model:depth": "Model: depth",
    "model:all": "Model: all metrics",
}
CURVE_METRICS = ["rel_cds", "rel_ars", "rel_ars_n", "rel_exon_min_n", "cds_p10_rel", "cds_share_lt10", "mapq_frac", "rel_mhc", "template_distance",
                 "cds_distance", "n_aa_changes", "n_tied", "has_warning", "mhc_contigs", "diff_contig", "score7"]
LOG_X = {"template_distance", "template_distance_per_kb", "cds_distance"}


def log(m):
    print("[72] " + m, file=sys.stderr, flush=True)


SCOPE = ""      # "" (all genes) or "classical8": stat_* tables are read from 71/classical8 first


def load(res, *names):
    for n in names:
        subs = (["71/" + SCOPE] if SCOPE and n.startswith("stat_") else []) + ["", "71", "70", "71c"]
        for sub in subs:
            p = os.path.join(res, sub, n)
            if os.path.exists(p):
                return pd.read_csv(p, sep="\t")
    return None


def style():
    """The repo's Nature style with Arial first: macOS' Helvetica.ttc cannot be loaded by this matplotlib build (every text
    draw fails with 'failed to load glyph'); Arial.ttf is the same look and equally acceptable to journals."""
    return matplotlib.rc_context(rc=dict(V.NATURE_RC, **{"font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                                                         "axes.unicode_minus": False}))


# ---------------------------------------------------------------------------
def fig_depth_distributions(res, out):
    hist = load(res, "dist_gene_depth_hist.tsv")
    q = load(res, "dist_gene_depth_quantiles.tsv")
    mh = load(res, "dist_person_mhc_hist.tsv")
    grid = load(res, "dist_mhc_vs_genome_grid.tsv")
    corr = load(res, "corr_gene_depth.tsv")
    if hist is None or q is None:
        return
    with style():
        fig = plt.figure(figsize=(V.mm(183), V.mm(150)))
        gs = fig.add_gridspec(2, 3, hspace=0.5, wspace=0.42, left=0.07, right=0.98, top=0.95, bottom=0.08)
        ax = fig.add_subplot(gs[0, 0])
        h = hist[(hist["platform"] == "all") & (hist["feature"] == "rel_cds") & hist["gene"].isin(CLASSICAL)]
        for i, g in enumerate(CLASSICAL[::-1]):
            d = h[h["gene"] == g].sort_values("bin_lo")
            y = d["n"].to_numpy(float)
            y = y / max(y.max(), 1) * 0.95
            x = d["bin_lo"].to_numpy() + 0.005
            ax.fill_between(x, i, i + y, color=V.JOURNAL_PALETTES["nature"][3], alpha=0.35, lw=0)
            ax.plot(x, i + y, color=V.JOURNAL_PALETTES["nature"][3], lw=0.5)
            ax.text(-0.02, i + 0.1, g, ha="right", va="bottom", fontsize=5.5)
        ax.set_xlim(0, 2)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_xlabel("CDS depth / genome depth")
        ax.set_title("Per-gene CDS depth (0.01 bins)", loc="left", fontsize=7)
        V.panel_letter(ax, "a")
        ax = fig.add_subplot(gs[0, 1])
        feats = [("rel_gene", "gene span"), ("rel_exons", "exons"), ("rel_cds", "CDS"), ("rel_ars", "peptide-binding exons")]
        qa = q[(q["platform"] == "all") & q["gene"].isin(CLASSICAL)]
        cols = V.JOURNAL_PALETTES["nature"]
        for j, (f, lab) in enumerate(feats):
            d = qa[qa["feature"] == f].set_index("gene").reindex(CLASSICAL)
            y = np.arange(len(CLASSICAL)) + (j - 1.5) * 0.18
            ax.hlines(y, d["p05"], d["p95"], color=cols[j], lw=0.5)
            ax.hlines(y, d["p25"], d["p75"], color=cols[j], lw=1.6)
            ax.plot(d["p50"], y, "o", ms=2, color=cols[j], label=lab)
        ax.set_yticks(range(len(CLASSICAL)))
        ax.set_yticklabels(CLASSICAL)
        ax.invert_yaxis()
        ax.set_xlabel("depth / genome depth (median, IQR, 5-95%)")
        ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4, columnspacing=0.8, **V.LEGEND_KW)
        V.panel_letter(ax, "b", dy=1.12)
        ax = fig.add_subplot(gs[0, 2])
        if mh is not None:
            d = mh[(mh["scale"] == "relative") & (mh["stratum_type"] == "platform")]
            for k, (pl, g) in enumerate(d.groupby("stratum")):
                if g["n"].sum() < 50:
                    continue
                ax.plot((g["bin_lo"] + 0.005), g["n"] / g["n"].sum() / 0.01, lw=0.8, label=pl,
                        color=cols[k])
            ax.set_xlim(0.4, 1.4)
            ax.set_xlabel("MHC-window depth / genome depth")
            ax.set_ylabel("density")
            ax.legend(**V.LEGEND_KW)
        ax.set_title("Per-person MHC depth", loc="left", fontsize=7)
        V.panel_letter(ax, "c")
        ax = fig.add_subplot(gs[1, 0])
        if grid is not None and len(grid):
            G = np.zeros((60, 60))
            for r in grid.itertuples():
                G[int(r.mhc_lo), int(r.genome_lo)] = r.n
            ax.pcolormesh(np.arange(61), np.arange(61), np.where(G > 0, G, np.nan), cmap="viridis", norm=matplotlib.colors.LogNorm(),
                          rasterized=True)
            ax.plot([0, 60], [0, 60], color="0.5", lw=0.4, ls="--")
            ax.set_xlim(5, 45)
            ax.set_ylim(5, 45)
            ax.set_xlabel("genome control depth (x)")
            ax.set_ylabel("MHC-window depth (x)")
        V.panel_letter(ax, "d")
        ax = fig.add_subplot(gs[1, 1])
        if mh is not None:
            d = mh[(mh["scale"] == "absolute") & (mh["stratum_type"] == "ancestry")]
            for k, (a, g) in enumerate(d.groupby("stratum")):
                if g["n"].sum() < 50:
                    continue
                ax.plot(g["bin_lo"] + 0.25, g["n"] / g["n"].sum() / 0.5, lw=0.8, label=a,
                        color=V.ANCESTRY_COLORS.get(a, "0.5"))
            ax.set_xlim(0, 45)
            ax.set_xlabel("MHC-window depth (x)")
            ax.set_ylabel("density")
            ax.legend(ncol=2, **V.LEGEND_KW)
        V.panel_letter(ax, "e")
        ax = fig.add_subplot(gs[1, 2])
        if corr is not None:
            gl = [g for g in ["A", "B", "C", "MICA", "DRA", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1", "TAP1", "TAP2"]
                  if g in set(corr["gene_a"])]
            M = corr.pivot(index="gene_a", columns="gene_b", values="spearman").reindex(index=gl, columns=gl)
            im = ax.imshow(M.to_numpy(float), cmap="RdBu_r", vmin=-1, vmax=1)
            ax.set_xticks(range(len(gl)))
            ax.set_xticklabels(gl, rotation=90)
            ax.set_yticks(range(len(gl)))
            ax.set_yticklabels(gl)
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cb.set_label("Spearman, CDS depth", fontsize=6)
            cb.ax.tick_params(labelsize=5)
        V.panel_letter(ax, "f")
        V.save_fig(fig, os.path.join(out, "fig_depth_distributions"))


def gene_model_axes(ax, feats, gene):
    f = feats[feats["gene"] == gene]
    g = f[f["kind"] == "gene"].iloc[0]
    strand = g["strand"]
    y0 = 0.5
    ax.hlines(y0, g["start"], g["end"], color="0.3", lw=0.6)
    for r in f[f["kind"] == "exons"].itertuples():
        ax.add_patch(Rectangle((r.start, y0 - 0.14), r.end - r.start, 0.28, fc="#B7C4D8", ec="none"))
    for r in f[f["kind"] == "cds"].itertuples():
        ax.add_patch(Rectangle((r.start, y0 - 0.28), r.end - r.start, 0.56, fc="#3C5488", ec="none"))
    for r in f[(f["kind"] == "ex") & f["ars"].astype(bool)].itertuples():
        ax.add_patch(Rectangle((r.start, y0 - 0.28), r.end - r.start, 0.56, fc="#E64B35", ec="none"))
    ax.annotate("", xy=(g["end"] if strand == "+" else g["start"], 0.05), xytext=(g["start"] if strand == "+" else g["end"], 0.05),
                arrowprops=dict(arrowstyle="->", lw=0.4, color="0.4"))
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    for s in ("left", "top", "right"):
        ax.spines[s].set_visible(False)


def fig_gene_depth_profile(res, out):
    prof = load(res, "profile_by_gene.tsv")
    feats = load(res, "features.tsv")
    if prof is None or feats is None:
        return
    with style():
        fig = plt.figure(figsize=(V.mm(183), V.mm(170)))
        outer = fig.add_gridspec(4, 2, hspace=0.55, wspace=0.28, left=0.07, right=0.985, top=0.96, bottom=0.075)
        for i, gene in enumerate(CLASSICAL):
            inner = outer[i // 2, i % 2].subgridspec(2, 1, height_ratios=[4, 1], hspace=0.05)
            a1, a2 = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
            d = prof[prof["gene"] == gene].sort_values("bin_start")
            x = (d["bin_start"] + d["bin_end"]) / 2
            a1.fill_between(x, d["rel_p10"], d["rel_p90"], color="#8491B4", alpha=0.3, lw=0)
            a1.fill_between(x, d["rel_p25"], d["rel_p75"], color="#3C5488", alpha=0.35, lw=0)
            a1.plot(x, d["rel_p50"], color="#3C5488", lw=0.6)
            a1.set_xlim(x.min(), x.max())
            a1.set_ylim(0, max(1.6, d["rel_p90"].max() * 1.05))
            a1.set_title(f"HLA-{gene}" if gene in CLASSICAL else gene, loc="left", fontsize=7)
            a1.tick_params(labelbottom=False)
            if i % 2 == 0:
                a1.set_ylabel("depth / genome")
            gene_model_axes(a2, feats, gene)
            a2.set_xlim(x.min(), x.max())
            a2.set_xlabel("chr6 (GRCh38)" if i >= 6 else "")
            if i == 0:
                a1.text(0.99, 0.97, "median, IQR, 10-90%", transform=a1.transAxes, ha="right", va="top", fontsize=5)
        fig.text(0.07, 0.012, "gene model: light = exons incl. UTR, dark blue = CDS, red = peptide-binding exons (reference transcript)",
                 fontsize=5, ha="left")
        V.save_fig(fig, os.path.join(out, "fig_gene_depth_profile"))


def fig_exon_depth(res, out):
    ex = load(res, "dist_exon_depth.tsv")
    if ex is None:
        return
    genes = [g for g in CLASSICAL if g in set(ex["gene"])]
    with style():
        fig, axes = plt.subplots(2, 4, figsize=(V.mm(183), V.mm(95)), sharey=True)
        for ax, g in zip(axes.ravel(), genes):
            d = ex[ex["gene"] == g].sort_values("cds_ord")
            x = np.arange(len(d))
            col = np.where(d["ars"].astype(bool), "#E64B35", "#3C5488")
            ax.vlines(x, d["p05"], d["p95"], color=col, lw=0.5)
            ax.vlines(x, d["p25"], d["p75"], color=col, lw=2.2)
            ax.scatter(x, d["p50"], c=col, s=7, zorder=3)
            ax.set_xticks(x)
            ax.set_xticklabels([str(int(v)) for v in d["exon_number"]])
            ax.set_title(f"HLA-{g}", loc="left", fontsize=7)
            ax.set_xlabel("exon")
            ax.axhline(1.0, color="0.7", lw=0.4, ls="--")
        axes[0, 0].set_ylabel("depth / genome")
        axes[1, 0].set_ylabel("depth / genome")
        fig.text(0.99, 0.01, "red: peptide-binding exons; median, IQR, 5-95%", ha="right", fontsize=5)
        fig.tight_layout(rect=(0, 0.03, 1, 1))
        V.save_fig(fig, os.path.join(out, "fig_exon_depth"))


# ---------------------------------------------------------------------------
def fig_error_curves(res, out, view):
    cur = load(res, "stat_curves.tsv")
    base = load(res, "stat_error_baseline.tsv")
    if cur is None or not len(cur) or view not in set(cur["view"]):
        return
    c = cur[cur["view"] == view]
    b = base[(base["view"] == view) & (base["gene"] == "ALL")]["error_pct"].iloc[0]
    metrics = [m for m in CURVE_METRICS if m in set(c["metric"])]
    nc = 4
    nr = int(np.ceil(len(metrics) / nc))
    with style():
        fig, axes = plt.subplots(nr, nc, figsize=(V.mm(183), V.mm(40 * nr)))
        axes = np.atleast_2d(axes)
        for ax, m in zip(axes.ravel(), metrics):
            d = c[c["metric"] == m].sort_values("x_median")
            col = FAM_COLOR.get({"score7": "legacy", "mhc_contigs": "assembly", "diff_contig": "assembly"}.get(m, None) or (
                "depth" if m.startswith(("rel_", "cds_p", "cds_min", "ars_", "cds_share", "mapq")) else "immuannot"))
            if len(d) <= 25:
                ax.errorbar(d["x_median"], d["error_pct"], yerr=[np.clip(d["error_pct"] - d["lo"], 0, None), np.clip(d["hi"] - d["error_pct"], 0, None)],
                            fmt="o", ms=2.5, color=col, lw=0.6, capsize=1)
            else:
                ax.fill_between(d["x_median"], d["lo"], d["hi"], color=col, alpha=0.25, lw=0)
                ax.plot(d["x_median"], d["error_pct"], color=col, lw=0.8)
            ax.axhline(b, color="0.5", lw=0.4, ls="--")
            if m in LOG_X and (d["x_median"] > 0).any():
                ax.set_xscale("function", functions=(lambda x: np.log1p(np.clip(x, 0, None)), np.expm1))
                hi = d["x_median"].max()
                ticks = [t for t in (0, 1, 2, 5, 10, 20, 50, 100, 200, 500) if t <= hi]
                ax.set_xticks(ticks)
                ax.set_xticklabels([str(t) for t in ticks])
                ax.minorticks_off()
            ax.set_xlabel(METRIC_LABEL.get(m, m))
            ax.set_ylim(bottom=0)
        for ax in axes.ravel()[len(metrics):]:
            ax.axis("off")
        for ax in axes[:, 0]:
            ax.set_ylabel("typing error (%)")
        fig.suptitle(VIEW_LABEL.get(view, view) + f"; dashed: overall error {b:.2f}%", x=0.01, ha="left", fontsize=7)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        V.save_fig(fig, os.path.join(out, f"fig_error_curves_{view}"))


def fig_auc(res, out):
    a = load(res, "stat_auc.tsv")
    if a is None or not len(a):
        return
    views = [v for v in ("trio_worst", "trio_child", "pair", "twin") if v in set(a["view"])]
    prim = "pair" if "pair" in views else views[0]
    order = a[a["view"] == prim].sort_values("auc")["metric"].tolist()
    mk = {"trio_child": "o", "pair": "s", "twin": "^", "trio_parent": "v", "trio_worst": "D"}
    with style():
        fig, axes = plt.subplots(1, 2, figsize=(V.mm(183), V.mm(max(80, 4.2 * len(order)))), sharey=True)
        for ax, col, title in ((axes[0], "auc", "AUC for typing error"), (axes[1], "auc_within_gene", "AUC within gene")):
            for j, v in enumerate(views):
                d = a[a["view"] == v].set_index("metric").reindex(order)
                y = np.arange(len(order)) + (j - (len(views) - 1) / 2) * 0.22
                fam = [d["family"].get(m, "depth") if m in d.index else "depth" for m in order]
                colors = [FAM_COLOR.get(f, "0.4") for f in fam]
                if col == "auc":
                    ax.hlines(y, d["auc_lo"], d["auc_hi"], color=colors, lw=0.6)
                ax.scatter(d[col], y, c=colors, marker=mk[v], s=9, zorder=3, label=VIEW_LABEL[v], edgecolors="none")
            ax.axvline(0.5, color="0.5", lw=0.4, ls="--")
            ax.set_xlabel(title)
            ax.set_xlim(0.4, 1.0)
        axes[0].set_yticks(range(len(order)))
        axes[0].set_yticklabels([METRIC_LABEL.get(m, m) for m in order])
        h, l = axes[1].get_legend_handles_labels()
        h = [plt.Line2D([], [], marker=mk[v], color="0.3", ls="", ms=3) for v in views]
        fam = [plt.Line2D([], [], marker="o", color=c, ls="", ms=3) for c in FAM_COLOR.values()]
        axes[1].legend(h + fam, l + list(FAM_COLOR), loc="lower right", **V.LEGEND_KW)
        fig.tight_layout()
        V.save_fig(fig, os.path.join(out, "fig_auc"))


def fig_metric_correlation(res, out):
    c = load(res, "corr_metric_metric.tsv")
    if c is None:
        return
    from scipy.cluster.hierarchy import leaves_list, linkage
    M = c.pivot(index="metric_a", columns="metric_b", values="spearman")
    M = M.fillna(0)
    z = linkage(1 - np.abs(M.to_numpy()), method="average")
    order = [M.index[i] for i in leaves_list(z)]
    M = M.loc[order, order]
    with style():
        fig, ax = plt.subplots(figsize=(V.mm(120), V.mm(110)))
        im = ax.imshow(M.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([METRIC_LABEL.get(m, m) for m in order], rotation=90)
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels([METRIC_LABEL.get(m, m) for m in order])
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        cb.set_label("Spearman (higher = worse call for both)", fontsize=6)
        cb.ax.tick_params(labelsize=5)
        fig.tight_layout()
        V.save_fig(fig, os.path.join(out, "fig_metric_correlation"))


def fig_pruning(res, out):
    p = load(res, "stat_prune.tsv")
    if p is None or not len(p):
        return
    p = p[p["scale"] == "pooled"]
    views = [v for v in ("trio_worst", "trio_child", "pair", "twin") if v in set(p["view"])]
    show = [("rel_cds", "CDS depth"), ("rel_ars", "Peptide-binding depth"), ("cds_p10_rel", "CDS 10th-pct depth"),
            ("template_distance", "Template distance"), ("cds_distance", "CDS distance"), ("n_tied", "Alleles tied"),
            ("score7", "Old 0-7 score"), ("model:all", "Model: all metrics")]
    cols = {"rel_cds": "#3C5488", "rel_ars": "#4DBBD5", "cds_p10_rel": "#8491B4", "template_distance": "#E64B35",
            "cds_distance": "#F39B7F", "n_tied": "#DC0000", "score7": "#7E6148", "model:all": "#1B1919"}
    with style():
        fig, axes = plt.subplots(2, len(views), figsize=(V.mm(183), V.mm(105)), squeeze=False)
        for j, v in enumerate(views):
            ax, ax2 = axes[0, j], axes[1, j]
            d0 = p[(p["view"] == v)]
            for m, lab in show:
                d = d0[d0["metric"] == m].sort_values("prune_pct")
                if not len(d):
                    continue
                ax.plot(d["prune_pct"], d["error_pct_kept"], color=cols[m], lw=0.8, label=lab)
                ax2.plot(d["prune_pct"], d["captured_pct"], color=cols[m], lw=0.8)
            any_d = d0[d0["metric"] == d0["metric"].iloc[0]].sort_values("prune_pct")
            ax.fill_between(any_d["prune_pct"], any_d["random_lo"], any_d["random_hi"], color="0.8", lw=0, label="random removal (95%)")
            ax2.plot([0, 50], [0, 100], color="0.6", lw=0.5, ls="--")
            ax.set_title(textwrap.fill(VIEW_LABEL[v], 24), loc="left", fontsize=6.5)
            ax.set_ylabel("typing error among calls kept (%)" if j == 0 else "")
            ax2.set_ylabel("share of all errors removed (%)" if j == 0 else "")
            ax2.set_xlabel("worst calls removed (% of records)")
            ax.set_ylim(bottom=0)
            ax.set_xlim(0, 50)
            ax2.set_xlim(0, 50)
            ax2.set_ylim(0, 100)
        h, l = axes[0, 0].get_legend_handles_labels()
        fig.legend(h, l, loc="lower center", ncol=5, **V.LEGEND_KW)
        fig.tight_layout(rect=(0, 0.06, 1, 1))
        V.save_fig(fig, os.path.join(out, "fig_pruning"))


def fig_grid(res, out, name):
    g = load(res, name)
    if g is None or not len(g):
        return
    a, b = name.replace("stat_grid_", "").replace(".tsv", "").split("__")
    from scipy.ndimage import gaussian_filter
    N = np.zeros((40, 40))
    K = np.zeros((40, 40))
    for r in g.itertuples():
        i, j = int(round(r.a_lo / 2.5)), int(round(r.b_lo / 2.5))
        N[i, j], K[i, j] = r.n, r.k
    Ns, Ks = gaussian_filter(N, 1.2), gaussian_filter(K, 1.2)
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.where(Ns >= 8, 100 * Ks / Ns, np.nan)
    with style():
        fig, ax = plt.subplots(figsize=(V.mm(89), V.mm(80)))
        im = ax.pcolormesh(np.linspace(0, 100, 41), np.linspace(0, 100, 41), R.T, cmap="magma_r", rasterized=True)
        ax.set_xlabel(METRIC_LABEL.get(a, a) + " (cohort percentile, 100 = worst)")
        ax.set_ylabel(METRIC_LABEL.get(b, b) + " (cohort percentile, 100 = worst)")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label("typing error (%), kernel-smoothed", fontsize=6)
        cb.ax.tick_params(labelsize=5)
        fig.tight_layout()
        V.save_fig(fig, os.path.join(out, f"fig_grid_{a}__{b}"))


GENE_ORDER = ["A", "C", "B", "MICA", "C4A", "DRA", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1", "TAP1"]


def matrix(df, col, genes):
    M = pd.DataFrame(np.nan, index=genes, columns=genes)
    for r in df.itertuples():
        if r.gene_a in M.index and r.gene_b in M.columns:
            M.loc[r.gene_a, r.gene_b] = getattr(r, col)
    return M


def fig_contig_truth(res, out):
    t = load(res, "contig_matrix_trio.tsv")
    c = load(res, "contig_matrix_cohort.tsv")
    if t is None:
        return
    t = t[t["stratum"] == "all"]
    genes = [g for g in GENE_ORDER if g in set(t["gene_a"])]
    panels = [("colocation_pct", "Contig co-location\n(% haplotypes, genes on one contig)", t),
              ("same_parent_all_pct", "Same parent, trio truth\n(% haplotypes, any contig)", t),
              ("same_parent_samecontig_pct", "Same parent, trio truth\n(% of pairs on one contig only)", t)]
    cmap = LinearSegmentedColormap.from_list("w2b", ["#F7FBFF", "#08306B"])
    with style():
        fig, axes = plt.subplots(1, 3, figsize=(V.mm(183), V.mm(66)))
        for k, (ax, (col, title, df)) in enumerate(zip(axes, panels)):
            M = matrix(df, col, genes)
            arr = M.to_numpy(float).copy()
            np.fill_diagonal(arr, np.nan)
            ax.imshow(arr, cmap=cmap, vmin=0, vmax=100)
            for i in range(len(genes)):
                for j in range(len(genes)):
                    if i != j and np.isfinite(arr[i, j]):
                        ax.text(j, i, f"{arr[i, j]:.0f}", ha="center", va="center", fontsize=3.8,
                                color="white" if arr[i, j] > 55 else "black")
            ax.set_xticks(range(len(genes)))
            ax.set_xticklabels(genes, rotation=90)
            ax.set_yticks(range(len(genes)))
            ax.set_yticklabels(genes if k == 0 else [])
            ax.set_title(title, fontsize=6.5, loc="left")
            V.panel_letter(ax, "abc"[k], dy=1.18)
        sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.Normalize(0, 100))
        cb = fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.01)
        cb.set_label("%", fontsize=6)
        cb.ax.tick_params(labelsize=5)
        V.save_fig(fig, os.path.join(out, "fig_contig_truth"))


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, help="directory with the pulled aggregate TSVs (or its 70/71/71c subfolders)")
    ap.add_argument("--out", default=os.path.join(REPO, "reports", "hla_popgen", "72_call_quality"))
    ap.add_argument("--only", help="comma list of figure names")
    ap.add_argument("--scope", choices=["all", "classical8"], default="all")
    a = ap.parse_args(argv)
    global SCOPE
    SCOPE = "" if a.scope == "all" else "classical8"
    if SCOPE:
        a.out = os.path.join(a.out, SCOPE)
    os.makedirs(a.out, exist_ok=True)
    jobs = {
        "depth_distributions": lambda: fig_depth_distributions(a.results, a.out),
        "gene_depth_profile": lambda: fig_gene_depth_profile(a.results, a.out),
        "exon_depth": lambda: fig_exon_depth(a.results, a.out),
        "error_curves": lambda: [fig_error_curves(a.results, a.out, v) for v in VIEW_LABEL],
        "auc": lambda: fig_auc(a.results, a.out),
        "metric_correlation": lambda: fig_metric_correlation(a.results, a.out),
        "pruning": lambda: fig_pruning(a.results, a.out),
        "grids": lambda: [fig_grid(a.results, a.out, n) for n in sorted(os.listdir(os.path.join(a.results, "71"))
                                                                         if os.path.isdir(os.path.join(a.results, "71")) else os.listdir(a.results))
                          if n.startswith("stat_grid_")],
        "contig_truth": lambda: fig_contig_truth(a.results, a.out),
    }
    want = a.only.split(",") if a.only else list(jobs)
    for k in want:
        log(f"figure {k}")
        jobs[k]()
    return 0


if __name__ == "__main__":
    sys.exit(main())
