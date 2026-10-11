#!/usr/bin/env python3
"""S03 WS6 follow-on -- Nature-grade multi-panel figure for the full-cohort KIR run.

Reads the 6 disclosure-safe aggregate TSVs `43_kir_full_aggregate.py` wrote to
reports/hla_popgen/43_kir_full_cohort/ and renders one 183 mm (double-column) static
matplotlib figure -> reports/hla_popgen/43_kir_full_cohort/fig_kir_full_cohort.{png,pdf}.

No participant-level data is read here -- every input is already an aggregate TSV with the
<20 disclosure rule applied (see `43_kir_full_aggregate.py`'s docstring). This script's own
disclosure obligation is just: never turn a blank/`<20` cell into a plotted 0 -- censored
cells are hatched (kir_gene_by_ancestry.tsv's novel_protein_pct, 7 gene x MID cells), not
silently dropped or drawn at height/color 0.

Gene order: the 17 KIR genes are laid out centromeric -> telomeric per the standard published
KIR gene map (framework genes 3DL3/3DP1/2DL4/3DL2 anchoring the array; e.g. Kulkarni et al 2008
Fig 1, Middleton & Gonzalez 2010 Table 1) -- used here ONLY for panel gene ordering. Note this
differs cosmetically from `41_kir_pilot.py`'s FRAMEWORK_GENES list order (3DL3, 2DL4, 3DP1, 3DL2),
which encodes adjacency for `classify_framework_miss()`, not a claim about the full 17-gene
physical order; that function only needs "what genes flank gene X", which is unaffected by
whether 2DL4 or 3DP1 is drawn first in this figure.

Usage:
  python3 scripts/hla_popgen/43b_kir_full_figure.py \\
      --in-dir reports/hla_popgen/43_kir_full_cohort \\
      --out-stem reports/hla_popgen/43_kir_full_cohort/fig_kir_full_cohort
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _THIS_DIR)
from _viz_common import (nature_style, mm, save_fig, panel_letter, ANCESTRY_ORDER,
                          ANCESTRY_COLORS, NATURE_DOUBLE_COL_MM, SUPPRESSED_COLOR,
                          diverging_cmap, diverging_norm, add_internal_watermark)
import _disclosure as _disc

# Standard published KIR gene map, centromeric -> telomeric (panel-ordering only; see module
# docstring). Framework genes (present on ~all haplotypes) marked with a dagger in the figure.
GENE_ORDER = ["KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DL5B", "KIR2DS3", "KIR2DS5",
              "KIR2DP1", "KIR2DL1", "KIR3DP1", "KIR2DL4", "KIR3DL1", "KIR3DS1", "KIR2DL5A",
              "KIR2DS1", "KIR2DS4", "KIR3DL2"]
FRAMEWORK_GENES = {"KIR3DL3", "KIR3DP1", "KIR2DL4", "KIR3DL2"}

GENE_SHORT = {g: g.replace("KIR", "") for g in GENE_ORDER}

NOVELTY_COLORS = {
    "known": "#999999",
    "genomic": "#56B4E9",   # novel, genomic-only (intron/UTR)
    "synonymous": "#E69F00",  # novel CDS, synonymous
    "protein": "#D55E00",   # novel protein -- highlight color, most consequential tier
}
NOVELTY_ORDER = ["known", "genomic", "synonymous", "protein"]
NOVELTY_LABEL = {"known": "known", "genomic": "novel, genomic-only",
                 "synonymous": "novel CDS, synonymous", "protein": "novel protein"}


def die(msg):
    sys.exit(f"FATAL: {msg}")


def _text_color_for_rgba(rgba):
    """Perceived-luminance (ITU-R BT.601 luma) text colour choice: black on light cells, white on
    dark ones, computed from the ACTUAL rendered RGBA of a colormap sample rather than a fixed
    "above the midpoint of vmin/vmax is dark" assumption. The previous version used the latter,
    which breaks for any non-monotonic-lightness colormap and, worse, for a diverging map where
    BOTH ends are dark and the middle is light -- exactly panel c's new deviation-from-pooled map.
    Orchestrator review (2026-09-25, full-size pass): dark text on viridis's dark-purple cells was
    unreadable (e.g. AFR/3DS1 "8"); OrRd's darkest red cells had the same problem for "38.9"."""
    r, g, b = rgba[0], rgba[1], rgba[2]
    luma = 0.299 * r + 0.587 * g + 0.114 * b
    return "#111111" if luma > 0.55 else "#FFFFFF"


def load_tables(in_dir):
    need = ["kir_run_summary", "kir_gene_summary", "kir_gene_by_ancestry", "kir_allele_freq",
            "kir_qc", "kir_content_by_ancestry"]
    out = {}
    for name in need:
        path = os.path.join(in_dir, name + ".tsv")
        if not os.path.exists(path):
            die(f"missing input table {path!r}")
        out[name] = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    return out


def to_num(s):
    return pd.to_numeric(s, errors="coerce")


# ---------------------------------------------------------------------------
# Panel a: per-gene haplotype presence % with Wilson CI (unrelated subset).
# ---------------------------------------------------------------------------
def panel_a_presence(ax, gene_summary):
    gs = gene_summary.set_index("gene").reindex(GENE_ORDER)
    y = np.arange(len(GENE_ORDER))
    pct = to_num(gs["presence_pct"])
    lo = to_num(gs["presence_ci_lo"])
    hi = to_num(gs["presence_ci_hi"])
    xerr = np.vstack([pct - lo, hi - pct])
    colors = ["#0072B2" if g in FRAMEWORK_GENES else "#56B4E9" for g in GENE_ORDER]
    ax.barh(y, pct, xerr=xerr, height=0.62, color=colors, edgecolor="none",
            error_kw=dict(elinewidth=0.5, ecolor="#333333", capsize=1.2, capthick=0.5))
    ax.set_yticks(y)
    ax.set_yticklabels([GENE_SHORT[g] for g in GENE_ORDER], fontsize=5)
    # Framework genes are marked by weight (bold) + the darker bar color set above, not a "†"
    # symbol (FIGURE_STYLE.md de-AI checklist item 12: no decorative symbols in labels) -- the
    # figure-wide caption spells out what bold + dark-blue mean once.
    for lbl, g in zip(ax.get_yticklabels(), GENE_ORDER):
        if g in FRAMEWORK_GENES:
            lbl.set_fontweight("bold")
    ax.invert_yaxis()
    ax.set_xlim(0, 105)
    ax.set_xlabel("haplotype presence (%)")
    ax.set_title("Presence per gene", fontsize=7, pad=3)
    ax.tick_params(axis="y", length=0)
    return y


def panel_b_novelty(ax, gene_summary, y):
    gs = gene_summary.set_index("gene").reindex(GENE_ORDER)
    left = np.zeros(len(GENE_ORDER))
    field_map = {"known": "pct_known", "genomic": "pct_novel_genomic_only",
                 "synonymous": "pct_novel_cds_synonymous", "protein": "pct_novel_protein"}
    for tier in NOVELTY_ORDER:
        vals = to_num(gs[field_map[tier]]).fillna(0.0).values
        ax.barh(y, vals, left=left, height=0.62, color=NOVELTY_COLORS[tier],
                edgecolor="none", label=NOVELTY_LABEL[tier])
        left = left + vals
    n_calls = gs["n_calls"].astype(int).values
    # `n=` counts as a genuine separate right-hand COLUMN, not a right-aligned label floated past
    # the bar end: the previous right-aligned version anchored each string's right edge at a fixed
    # x, which means a long string (e.g. "n=23,142") extends LEFTWARD by its own rendered width --
    # in a panel this narrow (100 data-units across a ~1/4-page-wide axes), that reliably reached
    # back far enough to sit ON TOP of the bar's own end (orchestrator review, full-size pass:
    # "text over data ... inadmissible"). Left-aligning at a fixed x instead means every label
    # starts at the same, deterministic position clear of x=100 regardless of digit count -- the
    # string can only grow rightward, into blank axes margin that exists for exactly this purpose.
    # The axes' own xlim is widened to make room; a thin grey rule marks where the bars end and
    # the label column begins, so the column reads as a deliberate table, not a stray overflow.
    N_COL_X = 108
    ax.axvline(100, color="#CCCCCC", lw=0.4, zorder=1)
    for yi, n in zip(y, n_calls):
        ax.text(N_COL_X, yi, f"n={n:,}", va="center", ha="left", fontsize=4.6, color="#888888",
                clip_on=False)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("allele calls (%)")
    ax.set_title("Novelty per gene", fontsize=7, pad=3)
    ax.tick_params(axis="y", length=0)
    ax.margins(x=0)
    handles = [Rectangle((0, 0), 1, 1, color=NOVELTY_COLORS[t]) for t in NOVELTY_ORDER]
    return handles


# ---------------------------------------------------------------------------
# Panels c/d: gene x ancestry heatmaps.
# ---------------------------------------------------------------------------
def _signed_cell_label(value):
    """Round-then-format a signed statistic for a heatmap cell: zero (after rounding, from either
    side) prints as a plain "0" with no sign; anything else keeps an explicit +/- prefix. Fixes
    "+0"/"0"/"-0" all appearing in the same panel (orchestrator review, pass 3) -- the previous
    version formatted the RAW float with "{:+.0f}" and then string-patched the one "-0.0" spelling
    it happened to notice, which still left an unsigned "0" (from a value that happened to round
    to +0 in a way "-0.0" didn't match) sitting next to signed "+0"s elsewhere in the same row."""
    rounded = int(round(value))
    return "0" if rounded == 0 else f"{rounded:+d}"


def _heatmap(ax, mat, y, cmap, norm, censored_mask, cbar_label, title, fmt="{:.0f}", signed=False,
             cbar=True, cbar_fraction=0.055, cbar_pad=0.08):
    cmap = plt.get_cmap(cmap) if isinstance(cmap, str) else cmap
    im = ax.imshow(mat, aspect="auto", cmap=cmap, norm=norm,
                    extent=[-0.5, mat.shape[1] - 0.5, len(GENE_ORDER) - 0.5, -0.5])
    for gi in range(mat.shape[0]):
        for ai in range(mat.shape[1]):
            if censored_mask[gi, ai]:
                ax.add_patch(Rectangle((ai - 0.5, gi - 0.5), 1, 1, facecolor=SUPPRESSED_COLOR,
                                        edgecolor="#999999", hatch="////", linewidth=0.3,
                                        zorder=4))
            elif not np.isnan(mat[gi, ai]):
                # Luminance-aware text colour computed from the CELL'S OWN rendered RGBA, not a
                # linear vmin/vmax midpoint assumption -- the latter is wrong for a diverging map
                # (panel c) where both ends are dark and the middle is light, and orchestrator
                # review (full-size pass) also caught it failing on viridis's/OrRd's darkest cells
                # under the old assumption. See `_text_color_for_rgba()`.
                rgba = cmap(norm(mat[gi, ai]))
                label = _signed_cell_label(mat[gi, ai]) if signed else fmt.format(mat[gi, ai])
                ax.text(ai, gi, label, ha="center", va="center",
                        fontsize=5.0, color=_text_color_for_rgba(rgba))
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels(ANCESTRY_ORDER, fontsize=5.3, rotation=90)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_title(title, fontsize=7, pad=3)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=2, pad=1)
    if cbar:
        # Slimmer, closer-fitted colorbar (orchestrator review: "colourbar must not crowd the
        # heatmap" -- previously a large fraction/shrink pairing left an oversized colorbar
        # column relative to a now-wider heatmap; a slimmer fraction + modest pad keeps a clean,
        # even gap). Panel c instead passes cbar=False and gets its own thin HORIZONTAL colorbar
        # in a reserved row below it (see build_figure) -- a vertical one here was still visually
        # crowding panel d's own column (orchestrator review, pass 3).
        cb = plt.colorbar(im, ax=ax, fraction=cbar_fraction, pad=cbar_pad, shrink=0.82, aspect=22)
        cb.set_label(cbar_label, fontsize=5.3)
        cb.ax.tick_params(labelsize=5.0, length=1.5)
    return im


def _pivot_num(gene_by_ancestry, value_col):
    piv = gene_by_ancestry.pivot(index="gene", columns="ancestry", values=value_col)
    piv = piv.reindex(index=GENE_ORDER, columns=ANCESTRY_ORDER)
    return piv.apply(to_num).values.astype(float)


def panel_c_deviation_heatmap(ax, gene_by_ancestry, gene_summary, y):
    """Redesigned from a flat presence-per-ancestry heatmap (which was almost entirely redundant
    with panel a's pooled presence per gene -- orchestrator review) to each ancestry's DEVIATION
    from the pooled presence in panel a: mat[gene, ancestry] = presence_pct[gene, ancestry] -
    presence_pct_pooled[gene]. Diverging, centred at 0 (`_viz_common.diverging_cmap/_norm`, the
    project's one signed-statistic colormap, already used for signed LD elsewhere) -- a cell near
    white means "this ancestry looks like the pooled average for this gene"; blue/red means
    under-/over-represented relative to pooled, which is the actually informative comparison once
    panel a already establishes the pooled value itself.

    Returns the `im` (AxesImage) so the caller can draw its OWN thin horizontal colorbar in a
    reserved row below this panel, rather than a vertical one wedged between panels c and d
    (orchestrator review, pass 3: "nudge c's colourbar so it doesn't crowd d")."""
    mat_pct = _pivot_num(gene_by_ancestry, "presence_pct")
    pooled = to_num(gene_summary.set_index("gene").reindex(GENE_ORDER)["presence_pct"]).values
    mat = mat_pct - pooled[:, None]
    censored = np.isnan(mat_pct)  # a censored ancestry cell has no presence_pct to begin with
    vmax = np.nanmax(np.abs(mat)) if np.isfinite(np.nanmax(np.abs(mat))) else 20.0
    vmax = max(5.0, np.ceil(vmax / 5.0) * 5.0)  # round up to a clean bracket for the colorbar
    im = _heatmap(ax, mat, y, diverging_cmap(), diverging_norm(-vmax, vmax), censored,
                  "deviation from pooled (pct pts)", "Presence: deviation from pooled",
                  signed=True, cbar=False)
    ax.set_yticklabels([GENE_SHORT[g] for g in GENE_ORDER], fontsize=5.3)
    for lbl, g in zip(ax.get_yticklabels(), GENE_ORDER):
        if g in FRAMEWORK_GENES:
            lbl.set_fontweight("bold")
    return im


def panel_d_novelty_heatmap(ax, gene_by_ancestry, y):
    # Integer %, per orchestrator review ("use integer percentages"): the previous 1-decimal
    # version's digits were hard to read even before the luminance fix, and integers are enough
    # precision at this cell size -- the underlying TSV keeps the full-precision numbers.
    mat = _pivot_num(gene_by_ancestry, "novel_protein_pct")
    censored = np.isnan(mat)
    vmax = np.nanmax(mat) if np.isfinite(np.nanmax(mat)) else 25.0
    vmax = max(vmax, 10.0)
    from matplotlib.colors import Normalize
    _heatmap(ax, mat, y, "OrRd", Normalize(vmin=0, vmax=vmax), censored, "novel protein (%)",
             "Protein novelty × ancestry", fmt="{:.0f}")


# ---------------------------------------------------------------------------
# Panel e: cA vs cB haplotype content by ancestry.
# ---------------------------------------------------------------------------
def panel_e_content(ax, content_by_ancestry):
    df = content_by_ancestry.set_index("ancestry").reindex(ANCESTRY_ORDER)
    pct_ca = to_num(df["pct_cA"])
    lo = to_num(df["pct_cA_ci_lo"])
    hi = to_num(df["pct_cA_ci_hi"])
    x = np.arange(len(ANCESTRY_ORDER))
    colors = [ANCESTRY_COLORS[a] for a in ANCESTRY_ORDER]
    ax.bar(x, pct_ca, yerr=[pct_ca - lo, hi - pct_ca], width=0.6, color=colors,
           edgecolor="none", error_kw=dict(elinewidth=0.5, ecolor="#333333", capsize=1.5,
                                            capthick=0.5))
    # No dashed 50% reference line (orchestrator review: an unexplained reference line invites a
    # "why 50?" the caption doesn't answer). The bar chart's own zero baseline is the only
    # reference needed to read cA-content magnitude and cross-ancestry differences.
    ax.set_xticks(x)
    ax.set_xticklabels(ANCESTRY_ORDER, fontsize=5.5)
    ax.set_ylim(0, 80)
    ax.set_ylabel("haplotypes cA (%)")
    ax.set_title("cA/cB content by ancestry", fontsize=7, pad=3)


# ---------------------------------------------------------------------------
# Panel f: QC strip -- framework presence + 2DL2/2DL3 and 3DL1/3DS1 co-occurrence.
# ---------------------------------------------------------------------------
def panel_f_qc(ax, qc):
    # Points + CI, not bars (FIGURE_STYLE.md de-AI checklist item 10 / known debt: a bar chart
    # over this narrow 93-97% range on a 0-100 axis would be unreadable, but a truncated bar axis
    # exaggerates the differences by area -- a bar's visual weight is its height from a baseline,
    # so any baseline above 0 misleads. A point encodes only its position, so it carries no such
    # baseline claim and can honestly sit in the narrow range that's actually informative.
    fw = qc[qc["metric"] == "framework_gene_presence"].copy()
    fw["n"] = to_num(fw["n"]); fw["pct"] = to_num(fw["pct"])
    fw["ci_lo"] = to_num(fw["ci_lo"]); fw["ci_hi"] = to_num(fw["ci_hi"])
    fw_order = ["KIR3DL3", "KIR3DP1", "KIR2DL4", "KIR3DL2"]
    fw = fw.set_index("item").reindex(fw_order)
    x = np.arange(len(fw_order))
    ax.errorbar(x, fw["pct"], yerr=[fw["pct"] - fw["ci_lo"], fw["ci_hi"] - fw["pct"]],
                fmt="o", ms=3.5, color="#0072B2", ecolor="#333333", elinewidth=0.6,
                capsize=2.0, capthick=0.6, zorder=3)
    ax.set_xlim(-0.5, len(fw_order) - 0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([g.replace("KIR", "") for g in fw_order], fontsize=5.2)
    ymin = min(fw["ci_lo"].min(), fw["pct"].min()) - 1.0
    ymax = max(fw["ci_hi"].max(), fw["pct"].max()) + 1.0
    ax.set_ylim(ymin, ymax)
    ax.set_ylabel("framework presence (%)")
    ax.set_title("QC: framework gene presence", fontsize=7, pad=3)
    # Co-occurrence numbers (KIR2DL2/KIR2DL3, KIR3DL1/KIR3DS1) previously sat as free text right
    # of this axes (ran off the figure's own edge -- check_layout() text_clipped, S04 WS-C, then
    # moved into an in-figure caption line). Orchestrator review (full-size pass): in-figure
    # caption text belongs in the README/FIGURES_INDEX caption, not baked into the raster --
    # removed here entirely; the numbers are already reported in
    # `reports/hla_popgen/43_kir_full_cohort/README.md`'s "QC sanity checks" section.


def build_figure(tables, out_stem, disclosure=_disc.PUBLIC):
    with nature_style():
        fig = plt.figure(figsize=(mm(NATURE_DOUBLE_COL_MM), mm(183)))
        # Pass 3 (orchestrator review, 2026-09-26): passes 1-2 anchored the panel-b legend by a
        # computed figure-fraction y-offset below panel b's axes. At thumbnail scale that looked
        # clear of row 2; at full resolution it actually overlapped panel e's title ("cA/cB
        # content by ancestry") -- a real overlap `check_layout()` did NOT catch (see this
        # script's module docstring note / FIGURES_INDEX.md for why: legend Text objects ARE
        # scanned by check_layout()'s Text-vs-Text pass, so the likely explanation is that the
        # offset math put the legend just outside the `tol_overlap_px` slack at whatever dpi the
        # linter drew at, while the final 600dpi raster's own text metrics differ by enough to
        # cross the actual pixel gap -- a computed offset is fundamentally fragile this way.
        # Fixed for real by giving the legend and panel c's colorbar their OWN reserved GridSpec
        # row (row 1, thin) between the main row and the e/f row -- a real Axes with real height
        # can never be "computed to just barely miss" the row below it, because GridSpec itself
        # guarantees non-overlapping rows.
        gs = fig.add_gridspec(3, 4, height_ratios=[2.3, 0.26, 1.0],
                               width_ratios=[1.15, 1.0, 0.85, 0.95],
                               hspace=0.30, wspace=0.75,
                               left=0.085, right=0.905, top=0.965, bottom=0.075)

        ax_a = fig.add_subplot(gs[0, 0])
        ax_b = fig.add_subplot(gs[0, 1])
        ax_c = fig.add_subplot(gs[0, 2])
        ax_d = fig.add_subplot(gs[0, 3])
        ax_legend = fig.add_subplot(gs[1, 1])
        ax_cbar_c = fig.add_subplot(gs[1, 2])
        ax_e = fig.add_subplot(gs[2, 0:2])
        ax_f = fig.add_subplot(gs[2, 2:4])

        y = panel_a_presence(ax_a, tables["kir_gene_summary"])
        panel_letter(ax_a, "a", dx=-0.42, dy=1.03)
        legend_handles = panel_b_novelty(ax_b, tables["kir_gene_summary"], y)
        panel_letter(ax_b, "b", dx=-0.10, dy=1.03)
        im_c = panel_c_deviation_heatmap(ax_c, tables["kir_gene_by_ancestry"],
                                          tables["kir_gene_summary"], y)
        panel_letter(ax_c, "c", dx=-0.12, dy=1.07)
        panel_d_novelty_heatmap(ax_d, tables["kir_gene_by_ancestry"], y)
        panel_letter(ax_d, "d", dx=-0.05, dy=1.07)
        panel_e_content(ax_e, tables["kir_content_by_ancestry"])
        panel_letter(ax_e, "e", dx=-0.10, dy=1.06)
        panel_f_qc(ax_f, tables["kir_qc"])
        panel_letter(ax_f, "f", dx=-0.12, dy=1.06)

        # Panel b's own legend, in its reserved row directly under it -- a real Axes (blanked to
        # a bare rectangle, no ticks/spines) rather than a floated `fig.legend()`, so its vertical
        # extent is accounted for by GridSpec like any other panel and cannot drift into row 2.
        # `axis("off")` alone leaves this Axes' default-view (0-1) tick-label Text objects
        # in existence with get_visible()==True -- matplotlib skips drawing them at draw time
        # (via the axes' own `axison=False`, a separate mechanism from a Text's own visibility
        # flag), so nothing is ever actually rendered, but `check_layout()`'s Text-vs-Text check
        # only looks at `get_visible()` and doesn't know about `axison` -- exactly the same class
        # of "phantom Text that never renders" false positive as the off-view-tick issue elsewhere
        # in this sprint (see 39_saturation_by_ancestry.py's `_prune_offview_ticklabels()`), just
        # via a different matplotlib mechanism. Found here as a real, reproducible check_layout()
        # false positive ('0.50' -- this Axes' default-view tick at the midpoint -- overlapping
        # the legend's own "novel, genomic-only" text). Fixed properly by removing the tick
        # Artists outright (`set_xticks([])`/`set_yticks([])`) instead of only hiding them.
        ax_legend.axis("off")
        ax_legend.set_xticks([])
        ax_legend.set_yticks([])
        ax_legend.legend(legend_handles, [NOVELTY_LABEL[t] for t in NOVELTY_ORDER],
                          loc="center", ncol=2, fontsize=5.3, frameon=False, handlelength=1.0,
                          handleheight=0.9, columnspacing=1.0, labelspacing=0.5,
                          bbox_to_anchor=(0.5, 0.5))

        # Panel c's colorbar: thin and HORIZONTAL, in its own reserved-row cell directly under
        # panel c, rather than a vertical bar squeezed into the c/d gutter (orchestrator review,
        # pass 3: "nudge c's colourbar so it doesn't crowd d"). `ax_cbar_c` fills its whole
        # GridSpec cell by default -- shrink it to an actual thin band (top-aligned, ~28% of the
        # cell's own height) with its own position override, or "thin" ends up meaning "as tall
        # as a normal small panel" just because the reserved row itself has that much height.
        cell_pos = ax_cbar_c.get_position()
        ax_cbar_c.set_position([cell_pos.x0, cell_pos.y0 + cell_pos.height * 0.62,
                                 cell_pos.width, cell_pos.height * 0.28])
        cbar_c = fig.colorbar(im_c, cax=ax_cbar_c, orientation="horizontal")
        cbar_c.set_label("deviation from pooled (pct pts)", fontsize=5.3, labelpad=2)
        cbar_c.ax.tick_params(labelsize=5.0, length=1.5)

        if disclosure == _disc.INTERNAL:
            add_internal_watermark(fig)
        return save_fig(fig, out_stem)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", default=os.path.join(
        os.path.dirname(_THIS_DIR), "..", "reports", "hla_popgen", "43_kir_full_cohort"))
    ap.add_argument("--out-stem", default=None)
    _disc.add_disclosure_arg(ap)
    args = ap.parse_args()
    in_dir = os.path.abspath(args.in_dir)
    if args.out_stem is not None:
        out_stem = args.out_stem
    elif args.disclosure == _disc.INTERNAL:
        out_stem = os.path.join(_disc.default_out_dir(_disc.INTERNAL, "43b", in_dir),
                                 "fig_kir_full_cohort")
    else:
        out_stem = os.path.join(in_dir, "fig_kir_full_cohort")
    if args.disclosure == _disc.INTERNAL:
        _disc.assert_internal_path_allowed(out_stem)
        os.makedirs(os.path.dirname(out_stem), exist_ok=True)

    tables = load_tables(in_dir)
    pdf_path, png_path = build_figure(tables, out_stem, disclosure=args.disclosure)
    print(f"wrote {pdf_path}\nwrote {png_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
