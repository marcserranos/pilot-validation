#!/usr/bin/env python3
"""Allele-frequency-by-ancestry figures, reproducible across all three cohorts (`--cohort
sr|lr|lr_td`). Implements VIZ_LIT.md's ranked-list items 1, 3, and pieces of 5 (research/VIZ_LIT.md
Part 5, "if we only build 8 figures"):

  1. Allele x ancestry frequency heatmap, hierarchically clustered on both axes (VIZ_LIT.md 1.2,
     "HLA map of the world" template) -- one panel per gene.
  2. Grouped allele-frequency bar chart per gene, with Wilson-score CIs (VIZ_LIT.md 1.1 + 4.3).
  3. Rank-frequency / allele-frequency-spectrum (Zipf) curves per gene, log-log, small multiples
     across genes, one line per ancestry (VIZ_LIT.md 1.3).
  4. Diversity indices (expected heterozygosity, Shannon entropy, rarefied allelic richness) per
     gene per ancestry, rarefied to a common copy-count so AFR/AMR/EAS/MID/SAS (small N in AoU)
     are fairly comparable to EUR (VIZ_LIT.md 1.10).

All frequencies computed at 2-field resolution (the SR cohort's ceiling; VIZ_LIT.md 3.1's
"coarsest common resolution" rule -- 2-field is also what the LR cohort's cross-cohort comparison
figure in 07_figures_crosscohort.py needs, so this script's 2-field choice stays consistent with
that one rather than picking a different resolution here).

Thin-N discipline (task hard rule): any (gene, ancestry[, allele]) cell with fewer than
--min-cell-n raw allele-copy observations is greyed out / flagged in every figure, never plotted
as a solid, trustworthy-looking estimate.

Usage (fixtures):
    python3 scripts/hla_popgen/tests/make_fixtures.py --outroot /tmp/hla_fixtures -n 300
    python3 scripts/hla_popgen/01_extract_rich.py --outroot /tmp/hla_fixtures --sample
    python3 scripts/hla_popgen/02_build_cohorts.py --outroot /tmp/hla_fixtures \\
        --table1 /tmp/hla_fixtures/hla_calls_rich.sample.tsv \\
        --cohort-full /tmp/hla_fixtures/immuannot_cohort_full.tsv \\
        --ancestry-preds /tmp/hla_fixtures/ancestry_preds.tsv \\
        --hla-genotypes /tmp/hla_fixtures/hla_genotypes.tsv --skip-mount-check --sample
    python3 scripts/hla_popgen/05_figures_frequency.py --outroot /tmp/hla_fixtures \\
        --table1 /tmp/hla_fixtures/hla_calls_rich.sample.tsv \\
        --cohort-membership /tmp/hla_fixtures/cohort_membership.sample.tsv \\
        --sr-genotypes /tmp/hla_fixtures/hla_genotypes.tsv --cohort lr

Real run (VM): python3 scripts/hla_popgen/05_figures_frequency.py --cohort lr
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _viz_common as vc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

TOP_N_HEATMAP = 15
TOP_N_BAR = 8


def build_long_frame(args, cohort_people, genes):
    """One row per (person_id, gene_bare, allele_2field) copy, joined to ancestry_pred. This is
    the single frame every figure in this script pivots off of."""
    person_ids = set(cohort_people["person_id"])
    anc_map = dict(zip(cohort_people["person_id"], cohort_people["ancestry_pred"]))
    if args.cohort == "sr":
        sr_path = args.sr_genotypes or os.path.join(args.outroot, "hla_genotypes.tsv")
        sr_long = vc.load_sr_genotypes(sr_path)
        sub = sr_long[sr_long["person_id"].isin(person_ids) & sr_long["gene_bare"].isin(genes)].copy()
        sub["allele_2field"] = sub["allele"].map(lambda a: vc.to_nfield(a, 2))
    else:
        table1_path = args.table1 or os.path.join(args.outroot, "hla_calls_rich.tsv")
        table1 = vc.load_table1(table1_path)
        sub = table1[table1["person_id"].isin(person_ids) & table1["gene_bare"].isin(genes)].copy()
        sub["allele_2field"] = sub["consensus"].map(lambda a: vc.to_nfield(a, 2))
    sub = sub.dropna(subset=["allele_2field"])
    sub["ancestry_pred"] = sub["person_id"].map(anc_map)
    sub = sub[sub["ancestry_pred"].isin(vc.ANCESTRY_ORDER)]
    sub["allele_label"] = sub["gene_bare"] + "*" + sub["allele_2field"]
    return sub[["person_id", "gene_bare", "allele_2field", "allele_label", "ancestry_pred"]]


def freq_matrix_for_gene(long_df, gene, min_cell_n):
    """Returns (freq_df, n_df) both indexed [allele_label] x columns ANCESTRY_ORDER; n_df holds
    raw copy counts (for thin-N flagging), freq_df holds allele frequency (copies/total copies in
    that ancestry x gene cell), NaN where the ancestry has zero total copies at this gene."""
    g = long_df[long_df["gene_bare"] == gene]
    if g.empty:
        return pd.DataFrame(), pd.DataFrame()
    counts = g.groupby(["allele_label", "ancestry_pred"]).size().unstack(fill_value=0)
    counts = counts.reindex(columns=vc.ANCESTRY_ORDER, fill_value=0)
    totals = counts.sum(axis=0)
    freq = counts.div(totals.replace(0, np.nan), axis=1)
    return freq, counts


# ---------------------------------------------------------------------------
# Figure 1: allele x ancestry heatmap, double-clustered per gene (VIZ_LIT.md 1.2)
# ---------------------------------------------------------------------------
def _average_linkage_order(dist):
    """Minimal average-linkage hierarchical clustering leaf order, via scipy if present, else the
    identity order (degrade gracefully rather than hard-crash if scipy is somehow unavailable)."""
    try:
        from scipy.cluster.hierarchy import linkage, dendrogram
        from scipy.spatial.distance import squareform
        if dist.shape[0] < 2:
            return list(range(dist.shape[0])), None
        condensed = squareform(dist, checks=False)
        Z = linkage(condensed, method="average")
        order = dendrogram(Z, no_plot=True)["leaves"]
        return order, Z
    except Exception as e:  # pragma: no cover - defensive only
        print(f"  WARNING: hierarchical clustering unavailable ({e}); using input order.",
              file=sys.stderr)
        return list(range(dist.shape[0])), None


def _draw_dendrogram_rescaled(ax, Z, n_leaves, orientation):
    """Draw dendrogram line segments rescaled into imshow's own 0..n_leaves-1 index coordinate
    frame, so it lines up with an adjacent heatmap WITHOUT matplotlib axis sharing (see the layout
    bug note at the call site). orientation='left' draws leaves top-to-bottom along the y-axis
    (rows); orientation='top' draws leaves left-to-right along the x-axis (columns)."""
    if Z is None or n_leaves < 2:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return
    from scipy.cluster.hierarchy import dendrogram
    dres = dendrogram(Z, no_plot=True)
    icoord = np.array(dres["icoord"], dtype=float)  # leaf spacing of 10, leaf i at 5+10*i
    dcoord = np.array(dres["dcoord"], dtype=float)
    leaf_pos = (icoord - 5.0) / 10.0  # rescaled to 0..n_leaves-1, matching imshow row/col index
    dmax = dcoord.max() if dcoord.size else 1.0
    for xs, ys in zip(leaf_pos, dcoord):
        if orientation == "top":
            ax.plot(xs, ys, color="#888888", linewidth=1)
        else:
            ax.plot(ys, xs, color="#888888", linewidth=1)
    if orientation == "top":
        ax.set_xlim(-0.5, n_leaves - 0.5)
        ax.set_ylim(0, dmax * 1.05 if dmax > 0 else 1)
    else:
        ax.set_xlim(dmax * 1.05 if dmax > 0 else 1, 0)  # leaves (x=0) adjacent to the heatmap
        ax.set_ylim(n_leaves - 0.5, -0.5)  # leaf 0 at top, matching imshow's row 0 at top


def plot_heatmap_gene(gene, freq, counts, min_cell_n, out_path):
    if freq.empty:
        return False
    # Restrict to top-N alleles by pooled frequency so the figure stays legible (VIZ_LIT.md 4.2).
    pooled = counts.sum(axis=1).sort_values(ascending=False)
    top_alleles = list(pooled.head(TOP_N_HEATMAP).index)
    freq = freq.loc[top_alleles]
    counts = counts.loc[top_alleles]
    thin = counts < min_cell_n

    # Cluster rows (alleles) and columns (ancestries) on Euclidean distance of frequency vectors,
    # sqrt-scaled per VIZ_LIT.md 1.2 ("linear color scale hides rare-allele signal").
    vals = freq.fillna(0.0).values
    if vals.shape[0] >= 2:
        row_dist = np.sqrt(((vals[:, None, :] - vals[None, :, :]) ** 2).sum(axis=2))
        row_order, row_Z = _average_linkage_order(row_dist)
    else:
        row_order, row_Z = list(range(vals.shape[0])), None
    if vals.shape[1] >= 2:
        col_dist = np.sqrt(((vals.T[:, None, :] - vals.T[None, :, :]) ** 2).sum(axis=2))
        col_order, col_Z = _average_linkage_order(col_dist)
    else:
        col_order, col_Z = list(range(vals.shape[1])), None

    alleles_ord = [top_alleles[i] for i in row_order]
    ancestries_ord = [vc.ANCESTRY_ORDER[i] for i in col_order]
    disp = freq.loc[alleles_ord, ancestries_ord]
    thin_disp = thin.loc[alleles_ord, ancestries_ord]

    fig = plt.figure(figsize=(9, 1.4 + 0.32 * len(alleles_ord)))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 5], height_ratios=[1, 6],
                           wspace=0.02, hspace=0.05)
    ax_heat = fig.add_subplot(gs[1, 1])
    ax_row = fig.add_subplot(gs[1, 0])
    ax_col = fig.add_subplot(gs[0, 1])
    ax_row.axis("off")
    ax_col.axis("off")

    # NOTE (layout bug found + fixed against the fixtures): matplotlib's `sharex`/`sharey` with a
    # raw scipy dendrogram() call breaks this figure -- scipy's dendrogram coordinate system uses
    # leaf spacing of 10 (leaf i at 5 + 10*i), and sharing that axis's scale with imshow's native
    # 0..N-1 index scale makes autoscale blow the shared axis out to the dendrogram's range,
    # squishing the heatmap into an invisible sliver (confirmed: rendered blank/white, labels
    # overlapping into unreadable text). Fix: draw the dendrogram lines manually, rescaled into
    # imshow's own index coordinate frame (leaf i -> position i), with NO axis sharing.
    _draw_dendrogram_rescaled(ax_row, row_Z, len(alleles_ord), orientation="left")
    _draw_dendrogram_rescaled(ax_col, col_Z, len(ancestries_ord), orientation="top")

    disp_sqrt = np.sqrt(disp.fillna(0.0).values)
    im = ax_heat.imshow(disp_sqrt, aspect="auto", cmap="viridis", vmin=0)
    ax_heat.set_xticks(range(len(ancestries_ord)))
    ax_heat.set_xticklabels(ancestries_ord, fontsize=9)
    # Layout bug found + fixed against the fixtures: row labels drawn on the LEFT (matplotlib's
    # default) render into the same horizontal band as the row dendrogram immediately to their
    # left, striking through both the labels and the dendrogram lines. Moving labels to the RIGHT
    # edge of the heatmap (opposite the dendrogram) removes the collision entirely.
    ax_heat.yaxis.tick_right()
    ax_heat.set_yticks(range(len(alleles_ord)))
    ax_heat.set_yticklabels(alleles_ord, fontsize=8)
    # Same bug, second half: a title set via ax_heat.set_title() lands directly above ax_heat,
    # which is exactly where the (small height-ratio) column dendrogram axis sits -- the title
    # text rendered on top of the dendrogram lines. fig.suptitle() draws outside the gridspec
    # entirely, well clear of ax_col's rendered content.
    fig.suptitle(f"{gene} allele frequency by ancestry (sqrt color scale, "
                 f"top {len(alleles_ord)} alleles, hierarchically clustered)", fontsize=10, y=1.04)

    # Thin-N flag: hatch/annotate cells below min_cell_n rather than plotting a bare solid color.
    for yi in range(len(alleles_ord)):
        for xi in range(len(ancestries_ord)):
            if thin_disp.values[yi, xi]:
                ax_heat.add_patch(plt.Rectangle((xi - 0.5, yi - 0.5), 1, 1, fill=False,
                                                 hatch="////", edgecolor="white", linewidth=0))
                ax_heat.text(xi, yi, "N<%d" % min_cell_n, ha="center", va="center",
                             fontsize=5.5, color="white")

    # Layout bug found + fixed against the fixtures: the row labels now render on ax_heat's RIGHT
    # edge (moved there above, away from the row dendrogram) -- a small colorbar `pad` (measured
    # from ax_heat's own edge, not from the label text's rendered extent) put the colorbar right
    # on top of those labels. A generous pad clears even the longest label ("HLA-DRB1*NN:NN").
    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.025, pad=0.22)
    cbar.set_label("sqrt(allele frequency)", fontsize=8)
    vc.savefig(fig, out_path)
    return True


# ---------------------------------------------------------------------------
# Figure 2: grouped frequency bar chart per gene with Wilson CIs (VIZ_LIT.md 1.1 + 4.3)
# ---------------------------------------------------------------------------
def plot_bar_gene(gene, freq, counts, totals, min_cell_n, out_path):
    if freq.empty:
        return False
    pooled = counts.sum(axis=1).sort_values(ascending=False)
    top_alleles = list(pooled.head(TOP_N_BAR).index)

    fig, ax = plt.subplots(figsize=(max(7, 1.1 * len(top_alleles)), 5))
    width = 0.8 / len(vc.ANCESTRY_ORDER)
    x = np.arange(len(top_alleles))
    for i, anc in enumerate(vc.ANCESTRY_ORDER):
        c = np.array([counts.loc[a, anc] if a in counts.index else 0 for a in top_alleles])
        n = totals.get(anc, 0)
        p, lo, hi = vc.wilson_ci(c, np.full_like(c, n, dtype=float))
        thin = c < min_cell_n
        offs = x + (i - (len(vc.ANCESTRY_ORDER) - 1) / 2) * width
        colors = [vc.ANCESTRY_COLORS[anc] if not t else "#DDDDDD" for t in thin]
        ax.bar(offs, np.nan_to_num(p), width=width * 0.95, color=colors,
               edgecolor=vc.ANCESTRY_COLORS[anc], linewidth=0.6)
        # Bug (crash) found against updated fixtures: zero-filling NaN p for the BAR (n=0 cell,
        # no data at all) but computing err_lo/err_hi from the un-zeroed p/lo/hi meant a plotted
        # y of 0 paired with a finite ci_lo produced err_lo = 0 - ci_lo < 0 -- matplotlib's
        # errorbar hard-rejects negative yerr. A cell with no data should have NO error bar, not
        # a zero-height one glued to a zero-height bar -- mask those cells out of errorbar
        # entirely rather than zero-filling them.
        has_data = ~np.isnan(p)
        err_lo = (p - lo)[has_data]
        err_hi = (hi - p)[has_data]
        assert np.all(err_lo >= -1e-9), f"{gene}/{anc}: negative err_lo {err_lo}"
        assert np.all(err_hi >= -1e-9), f"{gene}/{anc}: negative err_hi {err_hi}"
        if has_data.any():
            ax.errorbar(offs[has_data], p[has_data], yerr=[np.clip(err_lo, 0, None),
                                                             np.clip(err_hi, 0, None)],
                         fmt="none", ecolor="black", elinewidth=0.6, capsize=1.5, alpha=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([a.split("*", 1)[1] for a in top_alleles], rotation=45, ha="right")
    ax.set_ylabel("Allele frequency (Wilson 95% CI)")
    ax.set_title(f"{gene}: top {len(top_alleles)} alleles by pooled frequency, by ancestry\n"
                 f"(grey bars = fewer than {min_cell_n} observed copies in that ancestry)",
                 fontsize=10)
    # Bug found against the fixtures: a single ax.bar() call per ancestry mixes solid-colored and
    # thin-N-grey patches in one BarContainer. get_legend_handles_labels() only ever exposes ONE
    # representative color per label -- whichever bar happened to be first in that container -- so
    # an ancestry whose most-common-overall allele happened to be rare in ITS OWN group (small
    # ancestries like MID/SAS routinely hit this) showed a grey legend swatch, implying "this
    # ancestry's color is grey" instead of "this specific bar is thin". Explicit Patch proxies,
    # one per ancestry in its TRUE color, fix this regardless of which individual bars are thin.
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=vc.ANCESTRY_COLORS[a],
                              edgecolor=vc.ANCESTRY_COLORS[a]) for a in vc.ANCESTRY_ORDER]
    ax.legend(handles, vc.ANCESTRY_ORDER, fontsize=7, ncol=3, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    vc.savefig(fig, out_path)
    return True


# ---------------------------------------------------------------------------
# Figure 3: rank-frequency / Zipf curves (VIZ_LIT.md 1.3)
# ---------------------------------------------------------------------------
def plot_rank_frequency(long_df, genes, out_path):
    ncols = 4
    nrows = int(np.ceil(len(genes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.4 * nrows), squeeze=False)
    for idx, gene in enumerate(genes):
        ax = axes.flat[idx]
        g = long_df[long_df["gene_bare"] == gene]
        any_line = False
        for anc in vc.ANCESTRY_ORDER:
            sub = g[g["ancestry_pred"] == anc]
            if sub.empty:
                continue
            counts = sub["allele_label"].value_counts()
            total = counts.sum()
            freqs = (counts / total).sort_values(ascending=False).values
            if len(freqs) == 0:
                continue
            ax.step(range(1, len(freqs) + 1), freqs, where="mid",
                     color=vc.ANCESTRY_COLORS[anc], label=f"{anc} (n={total})", linewidth=1.3)
            any_line = True
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(gene, fontsize=10)
        ax.set_xlabel("allele rank", fontsize=8)
        ax.set_ylabel("frequency", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        if not any_line:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
    for idx in range(len(genes), nrows * ncols):
        axes.flat[idx].axis("off")
    handles, labels_ = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels_, loc="upper center", ncol=6, bbox_to_anchor=(0.5, 1.03),
               fontsize=8, frameon=False)
    fig.suptitle("Allele-frequency spectrum (rank-frequency / Zipf curve) by ancestry", y=1.06)
    fig.tight_layout()
    vc.savefig(fig, out_path)


# ---------------------------------------------------------------------------
# Figure 4: diversity indices, rarefied (VIZ_LIT.md 1.10)
# ---------------------------------------------------------------------------
def compute_diversity(long_df, genes, min_cell_n):
    """Returns a tidy DataFrame: gene, ancestry, n_copies, n_people, heterozygosity, shannon,
    richness_observed, richness_rarefied, rarefaction_target, thin_n."""
    rows = []
    for gene in genes:
        g = long_df[long_df["gene_bare"] == gene]
        if g.empty:
            continue
        anc_totals = g.groupby("ancestry_pred").size()
        valid_ancs = [a for a in vc.ANCESTRY_ORDER if anc_totals.get(a, 0) >= min_cell_n]
        if not valid_ancs:
            rare_target = None
        else:
            rare_target = int(anc_totals[valid_ancs].min())
        for anc in vc.ANCESTRY_ORDER:
            sub = g[g["ancestry_pred"] == anc]
            n_copies = len(sub)
            n_people = sub["person_id"].nunique()
            thin = n_copies < min_cell_n
            if n_copies == 0:
                rows.append(dict(gene=gene, ancestry=anc, n_copies=0, n_people=0,
                                  heterozygosity=np.nan, shannon=np.nan,
                                  richness_observed=0, richness_rarefied=np.nan,
                                  rarefaction_target=rare_target, thin_n=True))
                continue
            counts = sub["allele_label"].value_counts()
            freqs = (counts / n_copies).values
            het = vc.expected_heterozygosity(freqs)
            sh = vc.shannon_entropy(freqs)
            richness_obs = len(counts)
            richness_raref = (vc.rarefied_richness(counts.values, rare_target)
                               if rare_target and n_copies >= rare_target else np.nan)
            rows.append(dict(gene=gene, ancestry=anc, n_copies=n_copies, n_people=n_people,
                              heterozygosity=het, shannon=sh, richness_observed=richness_obs,
                              richness_rarefied=richness_raref, rarefaction_target=rare_target,
                              thin_n=thin))
    return pd.DataFrame(rows)


def plot_diversity(div_df, genes, out_path):
    metrics = [("heterozygosity", "Expected heterozygosity (1-Σp²)"),
               ("shannon", "Shannon entropy (nats)"),
               ("richness_rarefied", "Rarefied allelic richness")]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.2))
    for ax, (col, label) in zip(axes, metrics):
        for gi, gene in enumerate(genes):
            sub = div_df[div_df["gene"] == gene]
            for anc in vc.ANCESTRY_ORDER:
                row = sub[sub["ancestry"] == anc]
                if row.empty:
                    continue
                r = row.iloc[0]
                val = r[col]
                if pd.isna(val):
                    continue
                marker = "x" if r["thin_n"] else "o"
                ax.scatter(gi + (vc.ANCESTRY_ORDER.index(anc) - 2.5) * 0.08, val,
                           color=vc.ANCESTRY_COLORS[anc], marker=marker, s=28,
                           label=anc if gi == 0 else None)
        ax.set_xticks(range(len(genes)))
        ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
        ax.set_title(label, fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels_ = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels_, loc="upper center", ncol=6, bbox_to_anchor=(0.5, 1.08),
               fontsize=8, frameon=False)
    fig.suptitle("Diversity indices by gene x ancestry ('x' marker = thin-N cell, "
                 "richness rarefied to the smallest sufficiently-sampled ancestry's copy count)",
                 y=1.15, fontsize=10)
    fig.tight_layout()
    vc.savefig(fig, out_path)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    vc.add_common_args(ap)
    ap.add_argument("--genes", default=None,
                     help="Comma-separated bare gene names to restrict to (default: the 8 "
                          "classical genes -- non-classical genes are 07_figures_crosscohort.py's "
                          "job, LR-only anyway).")
    args = ap.parse_args()

    genes = (args.genes.split(",") if args.genes else vc.CLASSICAL_GENES_BARE)

    cohort_df = vc.load_cohort_membership(
        args.cohort_membership or os.path.join(args.outroot, "cohort_membership.tsv"))
    cohort_people, cohort_label = vc.select_cohort_people(cohort_df, args.cohort, args.td_max)
    print(f"Cohort: {cohort_label} -- {len(cohort_people)} people.", file=sys.stderr)
    if len(cohort_people) < 20:
        sys.exit(f"FATAL: only {len(cohort_people)} people in cohort {cohort_label!r} -- too few "
                 f"for a meaningful frequency figure.")

    out_dir = args.out_dir or vc.default_out_dir("05_figures_frequency", args.cohort, args.td_max)
    vc.ensure_dir(out_dir)

    long_df = build_long_frame(args, cohort_people, genes)
    if long_df.empty:
        sys.exit(f"FATAL: 0 usable (gene, ancestry-labeled) allele calls for cohort {cohort_label!r} "
                 f"restricted to genes {genes}. Check --table1/--sr-genotypes point at real data.")

    report = [f"# Allele-frequency figures -- cohort `{cohort_label}`\n",
              f"N people = {len(cohort_people)}. Genes analyzed: {', '.join(genes)}. "
              f"Resolution: 2-field. Thin-N threshold: {args.min_cell_n} raw allele copies.\n"]

    fig_paths = []
    for gene in genes:
        freq, counts = freq_matrix_for_gene(long_df, gene, args.min_cell_n)
        if freq.empty:
            print(f"  WARNING: no calls for gene {gene} in this cohort -- skipping its figures.",
                  file=sys.stderr)
            continue
        totals = counts.sum(axis=0)

        heat_path = os.path.join(out_dir, f"heatmap_{gene}.png")
        if plot_heatmap_gene(gene, freq, counts, args.min_cell_n, heat_path):
            fig_paths.append(heat_path)

        bar_path = os.path.join(out_dir, f"bar_{gene}.png")
        if plot_bar_gene(gene, freq, counts, totals, args.min_cell_n, bar_path):
            fig_paths.append(bar_path)

        report.append(f"\n## {gene}\n")
        report.append(f"Total copies observed by ancestry: "
                       f"{ {k: int(v) for k, v in totals.items()} }\n")
        report.append("| Allele | " + " | ".join(vc.ANCESTRY_ORDER) + " | pooled n |")
        report.append("|---|" + "---|" * (len(vc.ANCESTRY_ORDER) + 1))
        pooled = counts.sum(axis=1).sort_values(ascending=False)
        for allele in pooled.head(TOP_N_BAR).index:
            cells = []
            for anc in vc.ANCESTRY_ORDER:
                n = counts.loc[allele, anc] if allele in counts.index else 0
                p, lo, hi = vc.wilson_ci(np.array([n]), np.array([totals.get(anc, 0)],
                                                                   dtype=float))
                flag = " *thin*" if n < args.min_cell_n else ""
                if totals.get(anc, 0) == 0:
                    cells.append("--")
                else:
                    cells.append(f"{p[0]:.1%} [{lo[0]:.1%}-{hi[0]:.1%}]{flag}")
            report.append(f"| {allele} | " + " | ".join(cells) + f" | {int(pooled[allele])} |")

    rank_path = os.path.join(out_dir, "rank_frequency_spectrum.png")
    plot_rank_frequency(long_df, genes, rank_path)
    fig_paths.append(rank_path)

    div_df = compute_diversity(long_df, genes, args.min_cell_n)
    div_path = os.path.join(out_dir, "diversity_indices.png")
    plot_diversity(div_df, genes, div_path)
    fig_paths.append(div_path)

    report.append("\n## Diversity indices (rarefied)\n")
    report.append("| Gene | Ancestry | N copies | N people | Heterozygosity | Shannon | "
                   "Richness (observed) | Richness (rarefied) | Rarefaction target | Thin N |")
    report.append("|---|---|---|---|---|---|---|---|---|---|")
    for _, r in div_df.iterrows():
        raref_str = "" if pd.isna(r["richness_rarefied"]) else f"{r['richness_rarefied']:.2f}"
        het_str = "" if pd.isna(r["heterozygosity"]) else f"{r['heterozygosity']:.3f}"
        shannon_str = "" if pd.isna(r["shannon"]) else f"{r['shannon']:.3f}"
        report.append(f"| {r['gene']} | {r['ancestry']} | {int(r['n_copies'])} | "
                       f"{int(r['n_people'])} | {het_str} | {shannon_str} | "
                       f"{int(r['richness_observed'])} | {raref_str} | "
                       f"{r['rarefaction_target']} | {'YES' if r['thin_n'] else ''} |")

    report.append("\n## Figures\n")
    for p in fig_paths:
        report.append(f"- `{p}`")

    vc.write_report(os.path.join(out_dir, "frequency_report.md"), report)


if __name__ == "__main__":
    main()
