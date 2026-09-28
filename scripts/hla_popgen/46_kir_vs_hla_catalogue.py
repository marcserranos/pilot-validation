#!/usr/bin/env python3
"""WS-B (Sprint S04) -- KIR vs HLA reference-catalogue coverage, script 46.

Local-only: reads `44_kir_recurrence_saturation.py`'s already-pulled aggregate TSVs (never raw
per-person data) and puts KIR and HLA catalogue coverage on one comparable footing, per gene where
possible, per ORCHESTRATOR_HANDOFF.md WS-B. Candidate paper Figure 2 panel set.

Question: which region -- KIR (IPD-KIR 2.13.0) or HLA (IPD-IMGT/HLA 3.55.0) -- is better
represented by its reference catalogue, and does that differ gene by gene?

Method: for each gene (pooled ALL ancestry, the genomic identity granularity unless noted),
compute:
  - pct_novel_any    = 100 * s_obs(any_novel) / s_obs(genomic)      -- % of distinct alleles that
                         carry ANY sequence difference from the reference (novel at the genomic
                         level).
  - pct_novel_protein = 100 * s_obs(protein_novel) / s_obs(genomic) -- % that are novel at the
                         translated-protein level (a real amino-acid change).
  - good_turing_coverage -- incidence coverage, Chao & Jost (2012) (see 44's README + the
    interpretation guard reproduced below -- this is NOT "fraction of allele space explored").
  - completeness = s_obs / chao2                                    -- Chao2 (1987) richness
                         completeness, the fraction of ESTIMATED total allele richness already
                         observed. This is the metric that answers "how well is this gene's
                         reference catalogue represented in this cohort."
  - chao2_undetected_f0hat -- Chao2's estimated count of still-unseen alleles.
  - slope_per_1000    -- equal-N discovery slope (genomic level, pooled ALL ancestry, from
                         44's own equal_n_slope.tsv; see that file's own N* convention).

Interpretation guard (repeated from 45; do not violate when writing prose from this figure):
  Good-Turing coverage != Chao2 completeness. 99% incidence coverage does NOT mean "99% of the
  allele space explored" -- Chao2 completeness (s_obs/chao2) is the number that answers that
  question, and it is much lower (HLA ~64%, KIR ~84% pooled) than incidence coverage (~99% both).

Data-QA flag (not resolved by this script -- see README Caveats): `coverage_chao2.tsv`'s `s_obs`
column is NOT `<20`-masked per-gene (unlike `recurrence_classes.tsv`'s allele-class counts), even
though several KIR genes' `any_novel`/`protein_novel` s_obs are below 20 (e.g. KIR2DL5A=11,
KIR2DS3=16, KIR2DL2=16). This script uses these numbers as committed by 44 (VM-run, already local)
without re-masking them (re-masking someone else's committed aggregate file is out of scope for a
local-only figure script), but flags the question explicitly per CLAUDE.md's disclosure guard
("when in doubt... treat as withheld until the open question in context/DECISIONS.md is
resolved") -- Marc/Aleix should confirm this is consistent with 44's stated masking policy before
this figure is shown outside the team.

Deliverables:
  - `46_catalogue_metrics.tsv`: one row per gene x species, all metrics above.
  - `fig_catalogue_completeness.{png,pdf}`: Cleveland dot plot (rebuilt 2026-09-28 -- see
    `build_row_order()`/`fig_catalogue_completeness()` docstrings for the full rationale; the
    orchestrator rejected the earlier per-gene scatter version for leader-line spaghetti, an axis
    stretched past 100% to fit labels, a headline panel with no gene labels at all, and a tiny
    far-away legend). Rows = the 8 HLA classical genes + 17 KIR genes, in two blocks, ordered by
    protein-level Chao2 completeness -- the row IS the gene's label, so there are no leader lines
    anywhere. Panel a: Chao2 completeness per identity level (protein = headline, CDS, genomic =
    small grey "upper bound" tick), explained by a compact inline key, not a legend box. Panel b,
    sharing panel a's rows: % of distinct alleles novel at the protein level. Panel c: a
    per-ancestry strip -- HLA-minus-KIR protein completeness per ancestry plus pooled (the
    Simpson's-paradox point: HLA ahead in every single ancestry, near-tied only when pooled).
  - `fig_recurrence_composition.{png,pdf}`: paired recurrence-class composition bars (share of
    S_obs in each of eq1/eq2/gt2/ge20), KIR vs HLA, any-level vs protein-level novelty, using the
    corrected (curve-derived, exactly-partitioning) pooled counts from
    `45_kir_recurrence_figure.pooled_recurrence_from_curves()` -- reused, not re-derived, so the
    two scripts never disagree about the pooled recurrence split.

Usage:
    python3 scripts/hla_popgen/46_kir_vs_hla_catalogue.py \\
        --in-dir reports/hla_popgen/44_kir_recurrence_saturation \\
        --out-dir reports/hla_popgen/46_kir_vs_hla_catalogue
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

_m45 = None


def _load_45():
    """Lazily imports 45_kir_recurrence_figure.py by path (its filename isn't a valid Python
    module name) so this script reuses `pooled_recurrence_from_curves()`/`two_proportion_ztest()`
    verbatim rather than re-deriving the pooled-recurrence-from-curves logic a second time."""
    global _m45
    if _m45 is None:
        import importlib.util
        path = os.path.join(_THIS_DIR, "45_kir_recurrence_figure.py")
        spec = importlib.util.spec_from_file_location("kir_recurrence_figure_45", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _m45 = mod
    return _m45


SPECIES_ORDER = ["hla", "kir"]
SPECIES_LABEL = {"kir": "KIR", "hla": "HLA"}
SPECIES_COLOR = {"hla": vc.JOURNAL_PALETTES["nature"][0], "kir": vc.JOURNAL_PALETTES["nature"][3]}
LEVEL_ORDER = ["any_novel", "protein_novel"]
LEVEL_LABEL = {"any_novel": "any-level novel", "protein_novel": "protein-level novel"}
RECUR_CLASSES = ["eq1", "eq2", "gt2"]  # gt2 already includes ge20 -- avoid double-drawing a subset
RECUR_LABEL = {"eq1": "seen 1x", "eq2": "seen 2x", "gt2": "seen >2x (incl. ≥20x)"}


def load_tables(in_dir):
    rec = pd.read_csv(os.path.join(in_dir, "recurrence_classes.tsv"), sep="\t", dtype=str)
    curve = pd.read_csv(os.path.join(in_dir, "saturation_curves.tsv"), sep="\t")
    cov = pd.read_csv(os.path.join(in_dir, "coverage_chao2.tsv"), sep="\t")
    slope = pd.read_csv(os.path.join(in_dir, "equal_n_slope.tsv"), sep="\t")
    qc_path = os.path.join(in_dir, "kir_protein_catalogue_qc.tsv")
    protein_qc = pd.read_csv(qc_path, sep="\t") if os.path.exists(qc_path) else pd.DataFrame(
        columns=["gene", "protein_catalogue_covered", "reason"])
    return rec, curve, cov, slope, protein_qc


def build_gene_metrics(cov, slope, protein_qc=None, ancestry="ALL"):
    """Returns a tidy DataFrame, one row per (gene, species): pct_novel_any, pct_novel_protein,
    good_turing_coverage, completeness (s_obs/chao2), chao2_undetected_f0hat, slope_per_1000, all
    at the genomic identity granularity pooled over `ancestry` (default ALL) -- EXCEPT
    pct_novel_protein, whose denominator is `protein` s_obs (the SAME granularity as its own
    numerator, `protein_novel`), not `genomic` s_obs. 2026-09-27 fix (S04 coordinator item, VM
    smoke-test follow-up on commit 4a75657): `protein_novel` and `genomic` are two DIFFERENT
    identity granularities (a translated-protein identity vs. the curated full allele name) --
    dividing one by the other is comparing apples to oranges regardless of how correctly either
    side is computed, and was the proximate reason `pct_novel_protein` could exceed 100% even after
    44's own per-call identity-collapsing fix. `protein_novel` is BY CONSTRUCTION a subset of
    `protein` (44's `build_person_kir_identity`/`hla_level_mask` only ever add a call to
    `protein_novel` after first adding the same id to `protein`), so protein_novel/protein can only
    exceed 100% from an actual remaining bug, never a granularity artifact. Rows where the gene's
    `protein`/`protein_novel` level is 44's explicit `"NA"` (catalogue-uncovered gene -- see
    `kir_protein_catalogue_status()`) come through `pd.to_numeric(..., errors="coerce")` as NaN, so
    `pct_novel_protein` is NaN for that gene rather than a fabricated 0% or masked ratio."""
    sub = cov[cov["ancestry"] == ancestry].copy()
    for col in ("s_obs", "good_turing_coverage", "chao2", "chao2_undetected_f0hat"):
        sub[col] = pd.to_numeric(sub[col], errors="coerce")

    def s_obs_at(level):
        d = sub[sub["level"] == level][["gene", "species", "s_obs"]]
        return d.rename(columns={"s_obs": f"s_obs_{level}"})

    def chao2_at(level):
        d = sub[sub["level"] == level][["gene", "species", "chao2"]]
        return d.rename(columns={"chao2": f"chao2_{level}"})

    def gt_at(level):
        d = sub[sub["level"] == level][["gene", "species", "good_turing_coverage"]]
        return d.rename(columns={"good_turing_coverage": f"good_turing_coverage_{level}"})

    base = sub[sub["level"] == "genomic"][
        ["gene", "species", "s_obs", "good_turing_coverage", "chao2",
         "chao2_undetected_f0hat"]].rename(columns={"s_obs": "s_obs_genomic"})
    base = base.merge(s_obs_at("any_novel"), on=["gene", "species"], how="left")
    base = base.merge(s_obs_at("protein"), on=["gene", "species"], how="left")
    base = base.merge(s_obs_at("protein_novel"), on=["gene", "species"], how="left")
    base = base.merge(chao2_at("protein"), on=["gene", "species"], how="left")
    # CDS-level metrics (v4b, 2026-09-27 -- Marc's ask: base the headline KIR-vs-HLA comparison on
    # CDS and protein, NOT genomic. Genomic (true gene-span) identity is an UPPER BOUND, likely
    # inflated by person-specific span/UTR boundaries and intronic assembly noise the CDS-level
    # artifact gate can't see -- see README caveat 1). CDS has no "novel" sub-level of its own in
    # 44's exports (only genomic/any_novel and protein/protein_novel are novelty pairs); CDS's own
    # s_obs/chao2/coverage is exported here as a headline completeness metric in its own right.
    base = base.merge(s_obs_at("cds"), on=["gene", "species"], how="left")
    base = base.merge(chao2_at("cds"), on=["gene", "species"], how="left")
    base = base.merge(gt_at("cds"), on=["gene", "species"], how="left")

    base["pct_novel_any"] = 100.0 * base["s_obs_any_novel"] / base["s_obs_genomic"]
    base["pct_novel_protein"] = 100.0 * base["s_obs_protein_novel"] / base["s_obs_protein"]
    base["completeness"] = base["s_obs_genomic"] / base["chao2"]
    # Protein-level richness completeness -- answers WS-B's "which catalogue is better
    # represented" question AT THE PROTEIN LEVEL specifically (2026-09-27 coordinator ask). NaN
    # for a gene 44 marked catalogue-uncovered (its `protein`/`chao2_protein` are already NaN via
    # pd.to_numeric(..., errors="coerce") on 44's literal "NA") -- never divides to a fabricated 0.
    base["completeness_protein"] = base["s_obs_protein"] / base["chao2_protein"]
    # CDS-level richness completeness -- the OTHER headline granularity (v4b). Renamed the
    # pre-existing genomic-level `completeness` column is deliberately avoided (would break every
    # caller/test that already reads `completeness` as genomic) -- CDS gets its own column.
    base["completeness_cds"] = base["s_obs_cds"] / base["chao2_cds"]

    sl = slope[(slope["ancestry"] == ancestry) & (slope["level"] == "genomic")
               & (slope["curve"] == "distinct")][["gene", "species", "slope_per_1000"]]
    base = base.merge(sl, on=["gene", "species"], how="left")

    base["gene_display"] = np.where(base["species"] == "hla", "HLA-" + base["gene"], base["gene"])

    # Protein-catalogue coverage status (2026-09-27, `kir_protein_catalogue_qc.tsv`): HLA genes
    # have no such QC file (every HLA classical gene has a curated IPD-IMGT/HLA protein catalogue
    # entry) so they default to covered=True/"ok"; KIR genes are joined from 44's own QC table --
    # a gene 44 flagged uncovered (KIR2DP1/KIR3DP1, both pseudogenes with no catalogued reference
    # protein) must be rendered as an explicit NA/"pseudogene" case downstream, never plotted as
    # if its NaN pct_novel_protein/completeness_protein were an ordinary missing value.
    if protein_qc is None or protein_qc.empty:
        base["protein_catalogue_covered"] = True
        base["protein_catalogue_reason"] = "ok"
    else:
        qc = protein_qc.rename(columns={"protein_catalogue_covered": "_covered",
                                         "reason": "_reason"})
        base = base.merge(qc[["gene", "_covered", "_reason"]], on="gene", how="left")
        base["protein_catalogue_covered"] = base["_covered"].fillna(True)
        base["protein_catalogue_reason"] = base["_reason"].fillna("ok")
        base = base.drop(columns=["_covered", "_reason"])

    return base.sort_values(["species", "gene"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Cleveland dot plot (2026-09-28 rebuild -- orchestrator rejected the v4b spaghetti-scatter
# version of `fig_catalogue_completeness`: leader-line spaghetti, labels past 100% forcing a
# 140%-wide axis, no gene labels on the headline CDS panel, labels still touching dots, and a
# tiny far-away legend). Rows = genes; the y-axis tick label IS the gene identity, so there are
# NO leader lines anywhere and no free-repulsion label placement to get wrong -- a gene's row
# means the same thing in every panel that shares it, by construction, not by a repulsion pass
# that could disagree with the data.
# ---------------------------------------------------------------------------
GENE_BLOCKS = [("hla", "HLA  ·  8 classical genes"), ("kir", "KIR  ·  17 genes")]

# Marker SHAPE encodes identity level (consistent across species); marker COLOR encodes species
# (redundant with the block grouping, but keeps a reader oriented at a glance without a legend
# box). Genomic is a small grey tick, not a species color -- it is explicitly a context/upper-
# bound reference, not one of the two headline levels (44's Caveat 1; see README).
LEVEL_MARKER = {"protein": "o", "cds": "s", "genomic": "|"}
LEVEL_MARKERSIZE = {"protein": 4.6, "cds": 3.4, "genomic": 7.5}
GENOMIC_COLOR = vc.MISSING_COLOR  # "#999999" -- grey, deliberately not a species color


def build_row_order(metrics):
    """Returns the ordered list of row dicts used by EVERY reader of this figure's shared y-axis
    (the plotting code below, and this script's own binding test) -- one block-header dict per
    species (`{"kind": "header", "species": ..., "label": ...}`), followed by that species' gene
    rows sorted by protein-level Chao2 completeness (descending -- most-complete gene at the top
    of its block), NaN last (KIR2DP1/KIR3DP1, pseudogenes with no catalogued reference protein --
    never coerced to a fabricated 0, `na_position="last"` just puts them at the bottom of their
    own block, which is where a NaN protein-completeness gene belongs on a "sorted by protein
    completeness" axis). Plot and test share this ONE function so they cannot silently disagree
    about which row a gene is drawn in.
    """
    rows = []
    for species, header_label in GENE_BLOCKS:
        rows.append({"kind": "header", "species": species, "label": header_label})
        sub = metrics[metrics["species"] == species].sort_values(
            "completeness_protein", ascending=False, na_position="last")
        for _, r in sub.iterrows():
            rows.append({
                "kind": "gene", "species": species, "gene": r["gene"],
                "gene_display": r["gene_display"],
                "completeness": r["completeness"],
                "completeness_cds": r["completeness_cds"],
                "completeness_protein": r["completeness_protein"],
                "pct_novel_protein": r["pct_novel_protein"],
                "protein_catalogue_covered": bool(r["protein_catalogue_covered"]),
            })
    return rows


NO_PROTEIN_REF_SUFFIX = "  (no protein ref.)"


def row_label(row):
    """The single y-tick-label string for one `build_row_order()` row -- a block header's own
    label, or a gene row's `gene_display` (e.g. "HLA-A", "KIR2DL1"). This IS the row's identity
    label; there is no separate leader-lined text anywhere else in the figure.

    A gene with no catalogued reference protein (KIR2DP1/KIR3DP1, `protein_catalogue_covered=
    False`) gets `NO_PROTEIN_REF_SUFFIX` appended -- 2026-09-28 coordinator fix: an earlier
    version instead drew a HOLLOW marker at x=0 for these two genes, which a reader could still
    read as "completeness 0" at a glance despite the hollow styling; putting the caveat in the
    row's own label removes any ambiguity and needs no marker or in-panel note at all."""
    if row["kind"] == "header":
        return row["label"]
    if not row.get("protein_catalogue_covered", True):
        return row["gene_display"] + NO_PROTEIN_REF_SUFFIX
    return row["gene_display"]


def _style_row_axis(ax, rows, ys):
    """Shared y-axis setup for panels a/b: tick at every row (header rows included), header rows
    rendered bold/grey as a thin block label, gene rows left at normal weight -- and a thin grey
    rule drawn through each header row (inside the axes, so it never touches the tick-label text
    which lives outside the axes) as the visual block divider. A gene with no catalogued protein
    reference gets its WHOLE tick label (name + `NO_PROTEIN_REF_SUFFIX`) rendered in grey italic,
    distinct from both the normal black gene rows and the bold grey block headers."""
    ax.set_yticks(ys)
    ax.set_yticklabels([row_label(r) for r in rows])
    for tick_label, r in zip(ax.get_yticklabels(), rows):
        if r["kind"] == "header":
            tick_label.set_fontweight("bold")
            tick_label.set_fontsize(6.0)
            tick_label.set_color("#555555")
        elif not r.get("protein_catalogue_covered", True):
            tick_label.set_fontsize(5.6)
            tick_label.set_color("#888888")
            tick_label.set_fontstyle("italic")
        else:
            tick_label.set_fontsize(5.6)
    for r, y in zip(rows, ys):
        if r["kind"] == "header":
            ax.axhline(y, color="#DDDDDD", lw=0.7, zorder=0)
    # Headroom above the top row / below the bottom row so a header's bold tick label (drawn just
    # left of the axes, not inside it) never has to fight the panel's own inline key / x-tick row
    # for vertical space -- found necessary on full-size review of the first draft.
    ax.set_ylim(min(ys) - 0.7, max(ys) + 0.9)


def _row_shading(ax, rows, ys):
    """Very faint alternating row shading behind the markers (never a full gridline elsewhere in
    this figure) -- drawn identically on panels a and b from the SAME `rows`/`ys`, so a shaded (or
    unshaded) band falls on the exact same gene in both panels and the eye can track one gene
    across the shared y-axis without needing a leader line or a repeated label."""
    gene_i = 0
    for r, y in zip(rows, ys):
        if r["kind"] != "gene":
            continue
        if gene_i % 2 == 1:
            ax.axhspan(y - 0.5, y + 0.5, color="#F2F2F2", zorder=0, linewidth=0)
        gene_i += 1


def _inline_key(ax, items, y_offset_pt=24.0, fontsize=5.6, color="#444444", gap_pt=14.0,
                 marker_text_gap_pt=7.0, row_gap_pt=11.0):
    """A compact in-panel key for marker SHAPE, drawn along the top of `ax` (never a legend box,
    per FIGURE_STYLE.md's de-AI checklist item 3 / the task's own 'not a legend box' instruction)
    -- one real marker of each shape, followed by its label, each successive item starting where
    the previous item's rendered text actually ended (measured with the real renderer, the same
    draw-then-measure approach `check_layout()` itself uses) so the key never collides with itself
    regardless of font metrics or DPI.

    Added directly to the FIGURE (`fig.add_artist`/`fig.text`), not to `ax` (`ax.plot`/`ax.text`),
    and positioned in DISPLAY (pixel) coordinates converted from one final snapshot of `ax`'s own
    bbox -- deliberately not `ax.transAxes`. `constrained_layout` computes each axes' required
    margin from the tight bbox of everything BELONGING to that axes, including artists drawn with
    `transform=ax.transAxes` that stick out past its spines (e.g. this key's own rightmost items).
    Adding the key straight to `ax` therefore fed back into `constrained_layout`'s own next-draw
    margin computation, which kept shrinking `ax`'s width as more key items were added -- each
    item's fraction-space offset (correct against the `ax` bbox it was computed from) no longer
    matched the axes' NEW, narrower bbox by the time the next item (or `check_layout`'s own later
    redraw) measured it, producing real, reproducible text-on-text overlap despite every individual
    offset having been computed correctly at the time (found on first two render attempts: fixed
    key text was actually non-monotonic, then monotonic-but-still-overlapping, once other axes'
    content was also finalized after this key). A figure-level artist is invisible to
    `constrained_layout`'s per-axes margin accounting, so it cannot destabilize `ax`, and pixel
    coordinates (unlike `ax.transAxes` fractions) do not silently rescale if `ax`'s box does move
    for an unrelated reason.

    `items` is a list of (marker_kwargs, label) pairs, in the order plotted, left to right,
    starting just above `ax`'s own top-left corner, WRAPPING onto a new row (row spacing
    `row_gap_pt`) whenever the next item would extend past `ax`'s own right edge -- otherwise a
    long key row can spill rightward into the NEXT panel's own column and collide with that
    panel's title (both panels' titles sit at roughly the same absolute height, so this is a real
    cross-axes collision, found on full-size review: "novelty (headline)" (panel b's title)
    overlapped this key's own rightmost item once `gap_pt`/`marker_text_gap_pt` were widened to
    fix the earlier too-tight spacing).
    """
    import matplotlib.lines as mlines

    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax_bb = ax.get_window_extent(renderer=renderer)
    pt_to_px = fig.dpi / 72.0
    x0_px, y_px = ax_bb.x0, ax_bb.y1 + y_offset_pt * pt_to_px
    x_px = x0_px
    inv = fig.transFigure.inverted()
    for marker_kwargs, label in items:
        # Peek at this item's own rendered width BEFORE committing to placing it on the current
        # row, using an off-screen probe text -- wrapping needs to know the width in advance,
        # not just react after the fact.
        probe = fig.text(-1.0, -1.0, label, fontsize=fontsize)
        fig.canvas.draw()
        probe_w = probe.get_window_extent(renderer=fig.canvas.get_renderer()).width
        probe.remove()
        item_end_px = x_px + marker_text_gap_pt * pt_to_px + probe_w
        if x_px > x0_px and item_end_px > ax_bb.x1:
            x_px = x0_px
            y_px -= row_gap_pt * pt_to_px
        mx, my = inv.transform((x_px, y_px))
        ln = mlines.Line2D([mx], [my], transform=fig.transFigure, linestyle="none",
                            markeredgewidth=0.8, **marker_kwargs)
        fig.add_artist(ln)
        tx, ty = inv.transform((x_px + marker_text_gap_pt * pt_to_px, y_px))
        t = fig.text(tx, ty, label, fontsize=fontsize, color=color, va="center", ha="left")
        fig.canvas.draw()
        bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
        x_px = bb.x1 + gap_pt * pt_to_px


def _plot_completeness_panel(ax, rows, ys):
    """Panel a: Chao2 completeness (S_obs/Chao2), 0-1 axis. One marker per identity level: a
    small grey genomic tick (context, "upper bound" -- 44's Caveat 1), a CDS square, and the
    headline protein circle. KIR2DP1/KIR3DP1 (`protein_catalogue_covered=False`) get NO protein
    marker at all (2026-09-28 coordinator fix: an earlier hollow-marker-at-x=0 design still read
    as "completeness 0" at a glance) -- the row's own y-tick label already carries
    `NO_PROTEIN_REF_SUFFIX` (see `row_label()`), so no marker or in-panel note is needed here."""
    _row_shading(ax, rows, ys)
    for r, y in zip(rows, ys):
        if r["kind"] == "header":
            continue
        color = SPECIES_COLOR[r["species"]]
        ax.plot(r["completeness"], y, marker=LEVEL_MARKER["genomic"],
                 markersize=LEVEL_MARKERSIZE["genomic"], color=GENOMIC_COLOR,
                 markeredgewidth=1.1, zorder=2)
        ax.plot(r["completeness_cds"], y, marker=LEVEL_MARKER["cds"],
                 markersize=LEVEL_MARKERSIZE["cds"], color=color, markeredgecolor="white",
                 markeredgewidth=0.3, zorder=3)
        if r["protein_catalogue_covered"]:
            ax.plot(r["completeness_protein"], y, marker=LEVEL_MARKER["protein"],
                     markersize=LEVEL_MARKERSIZE["protein"], color=color,
                     markeredgecolor="white", markeredgewidth=0.3, zorder=4)
    ax.set_xlim(0.0, 0.72)
    ax.set_xlabel(r"Chao2 completeness (S$_{obs}$/Chao2)")


def _completeness_panel_key_items():
    """The (marker_kwargs, label) pairs for panel a's inline key -- factored out of
    `_plot_completeness_panel` so `fig_catalogue_completeness` can draw the key LAST, once every
    axes (a, b, c) has its final content: `_inline_key`'s draw-then-measure positioning needs
    constrained_layout to have already converged on ax_a's true final width, which it has not yet
    done while ax_b/ax_c are still empty (found on first render -- drawing the key immediately
    inside this function, before the other two panels existed, produced overlapping key text that
    only appeared once the OTHER axes were populated and the layout engine reflowed ax_a again)."""
    return [
        (dict(marker="o", markersize=LEVEL_MARKERSIZE["protein"], color="#444444"),
         "protein (headline)"),
        (dict(marker="s", markersize=LEVEL_MARKERSIZE["cds"], color="#444444"), "CDS"),
        (dict(marker="|", markersize=LEVEL_MARKERSIZE["genomic"], color=GENOMIC_COLOR),
         "genomic (upper bound)"),
    ]


def _plot_novelty_panel(ax, rows, ys):
    """Panel b, sharing panel a's rows: % of distinct alleles novel at the protein level (the
    headline novelty granularity). KIR2DP1/KIR3DP1 have no protein-level alleles to be novel
    among (`pct_novel_protein` is NaN by construction) -- left blank in this panel (never 0);
    the row's own y-tick label (shared with panel a) already carries the "(no protein ref.)"
    caveat, so this panel repeats nothing."""
    _row_shading(ax, rows, ys)
    for r, y in zip(rows, ys):
        if r["kind"] == "header" or not r["protein_catalogue_covered"]:
            continue
        color = SPECIES_COLOR[r["species"]]
        ax.plot(r["pct_novel_protein"], y, marker="o", markersize=LEVEL_MARKERSIZE["protein"],
                 color=color, markeredgecolor="white", markeredgewidth=0.3, zorder=3)
    ax.set_xlim(0.0, 100.0)
    ax.set_xlabel("% of distinct alleles\nnovel (protein level)")
    ax.set_title("novelty (headline)", fontsize=6.5, pad=14)


ANCESTRY_ORDER_SIMPSON = ["AFR", "AMR", "EAS", "EUR", "SAS"]


def _species_chao2_completeness(cov, ancestry, species, level="protein"):
    """Species-level (not per-gene) Chao2 completeness: sum(s_obs) / sum(chao2) across every gene
    of `species` at `level` within `ancestry`, straight from `coverage_chao2.tsv` -- the SAME
    aggregation the README's per-ancestry headline table is built from (confirmed by
    reproducing its numbers, e.g. AFR HLA protein 50.7-50.8%, to within rounding). Deliberately
    NOT an average of `build_gene_metrics()`'s per-gene completeness ratios -- that is a
    different (and not what this repo's committed headline numbers are) aggregation choice."""
    sub = cov[(cov["ancestry"] == ancestry) & (cov["species"] == species) & (cov["level"] == level)]
    s_obs = pd.to_numeric(sub["s_obs"], errors="coerce").sum()
    chao2 = pd.to_numeric(sub["chao2"], errors="coerce").sum()
    if chao2 == 0:
        return float("nan")
    return s_obs / chao2


def build_ancestry_simpson_rows(cov):
    """One row per well-powered ancestry (5, MID excluded per 44/46 precedent) plus a pooled
    "ALL" row: HLA-minus-KIR protein-level Chao2 completeness -- the Simpson's-paradox point from
    the README (HLA ahead of KIR in every single ancestry; pooled looks tied). Pooled is listed
    last and flagged `pooled=True` so the plotting code can set it apart from the 5 per-ancestry
    rows with its own thin divider, mirroring panel a/b's HLA/KIR block convention."""
    rows = []
    for anc in ANCESTRY_ORDER_SIMPSON:
        hla = _species_chao2_completeness(cov, anc, "hla")
        kir = _species_chao2_completeness(cov, anc, "kir")
        rows.append({"label": anc, "hla": hla, "kir": kir, "diff": hla - kir, "pooled": False})
    hla_all = _species_chao2_completeness(cov, "ALL", "hla")
    kir_all = _species_chao2_completeness(cov, "ALL", "kir")
    rows.append({"label": "Pooled (ALL)", "hla": hla_all, "kir": kir_all,
                 "diff": hla_all - kir_all, "pooled": True})
    return rows


PER_ANCESTRY_COLOR = "#555555"  # neutral dark grey -- NOT SPECIES_COLOR["kir"] (2026-09-28 fix:
# the per-ancestry points are a HLA-minus-KIR DIFFERENCE, not a KIR value, so coloring them KIR
# blue implied "this is about KIR" to a reader; only the pooled row gets the accent color).


def _plot_ancestry_simpson_panel(ax, sim_rows, y_top, ylim):
    """Panel c (optional per-ancestry strip): HLA-minus-KIR protein-level Chao2 completeness per
    ancestry plus pooled, 0 as the reference line -- "HLA more complete in every ancestry, tied
    when pooled." Per-ancestry points are a neutral dark grey (this is a DIFFERENCE, not a KIR
    value -- coloring them KIR blue would misleadingly imply otherwise); the pooled row alone gets
    the accent color, and is set off from the 5 per-ancestry rows by a thin divider, echoing panels
    a/b's HLA/KIR block convention.

    `y_top`/`ylim`: this panel's 6 rows start at `y_top` (the SAME y-value as panel a/b's HLA block
    header) and this axes is given panel a/b's own `ylim` verbatim -- not just a visually similar
    range -- so panel c's top row aligns EXACTLY with the top of the HLA block despite having far
    fewer rows than panels a/b share (2026-09-28 coordinator fix: previously each panel used its
    own independently-centered y-range, so panel c's content sat vertically centered in the middle
    of the figure, disconnected from panel a/b's top-aligned HLA block)."""
    n = len(sim_rows)
    ys = [y_top - i for i in range(n)]
    ax.axvline(0.0, color="#999999", lw=0.8, linestyle="--", zorder=1)
    for r, y in zip(sim_rows, ys):
        color = vc.ACCENT_COLOR if r["pooled"] else PER_ANCESTRY_COLOR
        ax.plot(r["diff"], y, marker="D" if r["pooled"] else "o", markersize=4.4,
                 color=color, markeredgecolor="white", markeredgewidth=0.3, zorder=3)
    divider_y = ys[-1] + 0.5
    ax.axhline(divider_y, color="#DDDDDD", lw=0.7, zorder=0)
    ax.set_yticks(ys)
    labels = ax.set_yticklabels([r["label"] for r in sim_rows])
    labels[-1].set_fontweight("bold")
    labels[-1].set_fontstyle("italic")
    ax.set_ylim(*ylim)
    diffs = [r["diff"] for r in sim_rows]
    pad = max(0.03, 0.15 * (max(diffs) - min(0.0, min(diffs))))
    ax.set_xlim(min(0.0, min(diffs)) - pad, max(diffs) + pad)
    ax.set_xlabel("HLA − KIR completeness\n(protein, Chao2)")
    ax.set_title("per-ancestry gap", fontsize=6.5, pad=14)
    # Short direction cue to the right of the zero line, just above the top ("AFR") row -- ASCII
    # "->" (not a unicode arrow glyph): this codebase has hit real font-rasterization crashes from
    # unicode arrows under this environment's Helvetica/mathtext stack before (see 47's own
    # ASCII-arrow fix, FIGURES_INDEX.md), so plain ASCII is used here too rather than re-risking it.
    ax.text(0.02, y_top + 0.55, "HLA more complete ->", fontsize=5.2, color="#666666",
             style="italic", ha="left", va="bottom")


def fig_catalogue_completeness(metrics, cov, out_stem):
    """Cleveland dot plot, rebuilt 2026-09-28 (orchestrator rejection of the v4b scatter version
    -- see the module-level comment above `GENE_BLOCKS`). Rows = the 8 HLA classical genes + 17
    KIR genes, grouped into two blocks and ordered within each block by protein-level Chao2
    completeness; the row itself is the gene's label, so there are no leader lines anywhere.

    Panel a: Chao2 completeness (S_obs/Chao2), 0-1 axis, one marker per identity level (protein =
    headline, CDS, genomic = small grey "upper bound" tick), explained by a compact inline key at
    the top of the panel rather than a legend box. KIR2DP1/KIR3DP1 (pseudogenes, no catalogued
    reference protein) get NO protein marker; their row's own y-tick label carries the caveat
    instead (see `row_label()`), never a fabricated 0.
    Panel b, sharing panel a's rows: % of distinct alleles novel at the protein level (the other
    v4b headline granularity), 0-100% axis.
    Panel c: a per-ancestry strip -- HLA-minus-KIR protein-level Chao2 completeness for the 5
    well-powered ancestries plus pooled, 0 as the reference line (the Simpson's-paradox point from
    the README: HLA ahead in every single ancestry, near-tied only when pooled). Top-aligned with
    panel a/b's HLA block, not vertically centered in its own axes.

    Width ratios ~2.2 : 1 : 1.1 with a small `wspace` (2026-09-28 coordinator fix): panels b and c
    were previously thin columns separated by large empty gaps at the figure's full 183mm width;
    the narrower gap and closer-to-content ratios make the three panels read as one figure.
    """
    rows = build_row_order(metrics)
    n = len(rows)
    ys = [n - 1 - i for i in range(n)]
    sim_rows = build_ancestry_simpson_rows(cov)

    fig_height_mm = max(150.0, n * 4.6 + 34.0)
    with vc.nature_style():
        fig, (ax_a, ax_b, ax_c) = plt.subplots(
            1, 3, figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(fig_height_mm)), dpi=150,
            gridspec_kw={"width_ratios": [2.2, 1.0, 1.1], "wspace": 0.28},
            constrained_layout=True)
        ax_b.sharey(ax_a)

        _style_row_axis(ax_a, rows, ys)
        _plot_completeness_panel(ax_a, rows, ys)
        _plot_novelty_panel(ax_b, rows, ys)
        plt.setp(ax_b.get_yticklabels(), visible=False)
        ax_b.tick_params(axis="y", length=0)
        # Panel c's top row (AFR) aligns with panel a/b's top row (the "HLA" block header, y=n-1)
        # -- pass that same y AND panel a/b's exact ylim (set inside `_style_row_axis`) through, so
        # the alignment is exact, not merely visually close.
        _plot_ancestry_simpson_panel(ax_c, sim_rows, y_top=n - 1, ylim=ax_a.get_ylim())

        vc.panel_letter(ax_a, "a", dx=-0.62, dy=1.16)
        vc.panel_letter(ax_b, "b", dy=1.16)
        vc.panel_letter(ax_c, "c", dy=1.16)
        # Drawn LAST, after every axes has its final content -- see `_completeness_panel_key_items`
        # docstring for why the key must not be placed while ax_b/ax_c are still empty.
        _inline_key(ax_a, _completeness_panel_key_items())
        vc.save_fig(fig, out_stem)


def fig_recurrence_composition(curve, out_stem):
    m45 = _load_45()
    with vc.nature_style():
        fig, axes = plt.subplots(1, 2, figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(65)),
                                  constrained_layout=True)
        pal = vc.JOURNAL_PALETTES["nature"]
        class_color = {c: pal[i + 4] for i, c in enumerate(RECUR_CLASSES)}
        for ax, level in zip(axes, LEVEL_ORDER):
            x = np.arange(len(SPECIES_ORDER))
            bottoms = np.zeros(len(SPECIES_ORDER))
            any_data = False
            for c in RECUR_CLASSES:
                shares = []
                for sp in SPECIES_ORDER:
                    p = m45.pooled_recurrence_from_curves(curve, sp, level)
                    if p is None or p["s_obs"] == 0:
                        shares.append(0.0)
                    else:
                        shares.append(100.0 * p[c] / p["s_obs"])
                        any_data = True
                ax.bar(x, shares, bottom=bottoms, width=0.5, color=class_color[c],
                       label=RECUR_LABEL[c])
                bottoms += np.array(shares)
            ax.set_xticks(x)
            ax.set_xticklabels([SPECIES_LABEL[s] for s in SPECIES_ORDER])
            ax.set_title(LEVEL_LABEL[level], fontsize=7)
            ax.set_ylim(0, 100)
            if not any_data:
                t = ax.text(0.5, 0.5, "no alleles\nat this level", ha="center", va="center",
                             fontsize=6, color="#666666", transform=ax.transAxes)
                vc.mark_label(t)
            # HLA protein-level bar is 0-height by construction (flagged data question) --
            # annotate rather than leave an unexplained blank bar.
            hla_p = m45.pooled_recurrence_from_curves(curve, "hla", level)
            if level == "protein_novel" and hla_p is not None and hla_p["s_obs"] == 0:
                t = ax.annotate("0 alleles", (0, 2), ha="center", va="bottom", fontsize=5.5,
                                 color="#666666")
                vc.mark_label(t)
        axes[0].set_ylabel("% of distinct novel alleles")
        axes[0].legend(loc="upper left", bbox_to_anchor=(0.0, -0.18), ncol=3, **vc.LEGEND_KW)
        vc.panel_letter(axes[0], "a")
        vc.panel_letter(axes[1], "b")
        vc.save_fig(fig, out_stem)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", default="reports/hla_popgen/44_kir_recurrence_saturation")
    ap.add_argument("--out-dir", default="reports/hla_popgen/46_kir_vs_hla_catalogue")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    rec, curve, cov, slope, protein_qc = load_tables(args.in_dir)
    metrics = build_gene_metrics(cov, slope, protein_qc)
    metrics_path = os.path.join(args.out_dir, "46_catalogue_metrics.tsv")
    metrics.to_csv(metrics_path, sep="\t", index=False)
    print(f"  wrote {metrics_path}", file=sys.stderr)

    # Per-ancestry gene metrics (Marc's ask: "check it per gene and per ancestry") -- one combined
    # tidy TSV with an `ancestry` column, not 5 separate files, so a reader can filter/pivot in one
    # place. MID excluded (39/44 precedent: not well-powered).
    ANCESTRY_ORDER_46 = ["AFR", "AMR", "EAS", "EUR", "SAS"]
    per_anc = []
    for anc in ANCESTRY_ORDER_46:
        m = build_gene_metrics(cov, slope, protein_qc, ancestry=anc)
        m.insert(0, "ancestry", anc)
        per_anc.append(m)
    metrics_ancestry = pd.concat(per_anc, ignore_index=True)
    metrics_ancestry_path = os.path.join(args.out_dir, "46_catalogue_metrics_by_ancestry.tsv")
    metrics_ancestry.to_csv(metrics_ancestry_path, sep="\t", index=False)
    print(f"  wrote {metrics_ancestry_path}", file=sys.stderr)

    fig_catalogue_completeness(metrics, cov,
                                os.path.join(args.out_dir, "fig_catalogue_completeness"))
    fig_recurrence_composition(curve, os.path.join(args.out_dir, "fig_recurrence_composition"))

    print(f"[46] 2 figures + 46_catalogue_metrics.tsv written to {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
