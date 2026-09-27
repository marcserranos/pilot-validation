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
  - `fig_catalogue_completeness.{png,pdf}`: per-gene scatter, x = % novel (any level), y = Chao2
    richness completeness, color = species, direct gene labels (no legend for genes; a small
    2-entry species legend since color alone needs a key).
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
# Direct-labeling helper for scattered points: places each label at a small fixed offset from its
# point, then iteratively pushes apart any pair of labels whose rendered bounding boxes actually
# overlap (measured via the real renderer, not a heuristic distance threshold -- a fixed-distance
# heuristic under/over-corrects depending on font metrics and DPI, and `check_layout()`'s own
# overlap check (S04 WS-C) is exactly this kind of real-bbox measurement, so this helper mirrors
# it rather than guessing at a "close enough" threshold that could still fail the linter).
# ---------------------------------------------------------------------------
def _nudge_text_points(text_artist, dy_points):
    x_off, y_off = text_artist.xyann
    text_artist.xyann = (x_off, y_off + dy_points)


def _label_points(ax, xs, ys, labels, colors, fontsize=5.0, dx_pt=3.0, max_iter=200):
    fig = ax.figure
    texts = []
    for x, y, lab, col in zip(xs, ys, labels, colors):
        t = ax.annotate(lab, (x, y), xytext=(dx_pt, 2.0), textcoords="offset points",
                         fontsize=fontsize, color=col, ha="left", va="center",
                         annotation_clip=False)
        vc.mark_label(t)
        texts.append(t)
    if not texts:
        return texts
    # Match _viz_common.check_layout()'s own dpi floor: on at least one dev machine, freetype
    # raises "RuntimeError: failed to load glyph" rasterizing at fig.dpi<150 (matplotlib's default
    # figure dpi is 100) but never at >=150 -- a font-rasterization floor, not an API misuse.
    # save_fig() -> check_layout() will redraw at >=150 dpi anyway, so measuring at that same
    # floor here also keeps the pixel geometry consistent with what check_layout() will verify.
    _orig_dpi = fig.dpi
    if fig.dpi < 150:
        fig.dpi = 150
    try:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        px_to_pt = 72.0 / fig.dpi
        _label_points_repel(texts, renderer, px_to_pt, max_iter)
    finally:
        fig.dpi = _orig_dpi
    return texts


def _label_points_repel(texts, renderer, px_to_pt, max_iter):
    fig = texts[0].get_figure()
    for _ in range(max_iter):
        moved = False
        boxes = [t.get_window_extent(renderer) for t in texts]
        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                bi, bj = boxes[i], boxes[j]
                if not bi.overlaps(bj):
                    continue
                overlap_y = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
                if overlap_y <= 0:
                    continue
                push_px = overlap_y / 2.0 + 6.0
                push_pt = push_px * px_to_pt
                if bi.y0 <= bj.y0:
                    _nudge_text_points(texts[i], -push_pt)
                    _nudge_text_points(texts[j], push_pt)
                else:
                    _nudge_text_points(texts[i], push_pt)
                    _nudge_text_points(texts[j], -push_pt)
                moved = True
        if not moved:
            break
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
    return texts


def _scatter_panel(ax, metrics, x_col, y_col, na_note=None, column_species=(),
                    column_fontsize=5.0, column_min_gap_frac=0.05, label_genes=True):
    """One completeness-vs-novelty scatter panel, direct-labeled, both species. Rows where
    `x_col`/`y_col` is NaN (a catalogue-uncovered gene, e.g. KIR2DP1/KIR3DP1 at the protein level)
    are NEVER silently dropped or plotted as 0 -- they are listed by name in an in-panel note
    (`na_note`, e.g. "KIR2DP1, KIR3DP1: no catalogue protein (pseudogene)"), so a reader sees
    explicitly which genes are missing and why, per the 2026-09-27 coordinator instruction.

    Labeling is done PER SPECIES so a tight cluster of one species doesn't force the other
    species' already-well-spaced labels into a shared far-away column. Species named in
    `column_species` use 45's shared-x-column end-of-line labeling (`_label_curve_ends`,
    originally built for curve endpoints but equally valid for any (x, y) point), positioned
    just past THAT SPECIES' OWN rightmost point -- i.e. the KIR protein-level cluster (~13 genes
    packed into x=90-100/y=0.06-0.25) gets a short local column with leader lines, not a
    figure-spanning one. Every other species uses `_label_points()`'s free per-point repulsion,
    which reads fine when points are already spread out (e.g. panel a, or panel b's HLA points).
    The free style was tried for the KIR cluster first and passed `check_layout(strict=True)`
    with zero TEXT-vs-TEXT violations yet still visually sat labels on top of marker dots --
    `check_layout()`'s overlap check is text-vs-text/decoration only, it does not compare a label
    against a plain scatter marker, so that failure mode is a real "text on data" violation
    (FIGURE_STYLE.md) the mechanical linter cannot catch; caught only by rendering and viewing the
    PNG at full size (2026-09-27 visual review)."""
    handles = []
    per_species = {}
    na_genes = []
    for sp in SPECIES_ORDER:
        s_all = metrics[metrics["species"] == sp]
        na_genes += s_all.loc[s_all[x_col].isna() | s_all[y_col].isna(), "gene_display"].tolist()
        s = s_all.dropna(subset=[x_col, y_col])
        h = ax.scatter(s[x_col], s[y_col], s=16, color=SPECIES_COLOR[sp],
                        label=SPECIES_LABEL[sp], zorder=3, edgecolors="white", linewidths=0.3)
        handles.append(h)
        per_species[sp] = (s[x_col].tolist(), s[y_col].tolist(), s["gene_display"].tolist())
    # Bottom padding (-0.04, not 0): a free-style label for a near-zero-completeness gene is
    # offset only 2pt above its own point (`_label_points`' default), which can otherwise render
    # low enough to collide with the x-tick-label row just below the axes' own spine (found on
    # full-size review: 'KIR2DP1' over the '100' x-tick in the CDS panel -- the two are plain Text
    # objects, and check (a3)'s own-axes-spine exemption doesn't cover a sibling tick LABEL). No
    # real data point is ever negative, so this buffer band is always empty of data.
    ax.set_ylim(-0.04, 1.05)
    # Right-pad xlim past the data max (v4b, 3-panel layout): a free-style label offset a few
    # points to the right of its point (`_label_points`' dx_pt=3.0) can otherwise land on top of
    # the rightmost x-tick label itself (found on full-size review: 'KIR2DP1' over the '100'
    # tick in the CDS panel) -- pad by 15% of the observed x-range so labels near the right edge
    # have somewhere to go.
    all_x_vals = [v for sp in SPECIES_ORDER for v in per_species[sp][0]]
    xmax_data = max(all_x_vals) if all_x_vals else 100.0
    ax.set_xlim(left=-2, right=max(xmax_data * 1.15, xmax_data + 10))
    if not label_genes:
        # Both species' any-level novelty ranges overlap substantially at the CDS granularity
        # (unlike the protein panel, where KIR's cluster and HLA's cluster sit at very different
        # x) -- neither the free-repulsion nor the column-leader style can place ~30 gene labels
        # here without collisions between species (tried both, 2026-09-27: free-vs-free overlaps,
        # and free HLA labels landing on KIR's column leader lines). Gene identities for this
        # panel are in `46_catalogue_metrics.tsv`/`_by_ancestry.tsv`; the caption says so.
        if na_genes:
            note = na_note or (", ".join(na_genes) + ": no catalogue data at this level")
            t = ax.text(0.02, 0.98, note, transform=ax.transAxes, fontsize=5.2, color="#666666",
                         ha="left", va="top", style="italic", wrap=True)
            vc.mark_label(t)
        return handles
    m45 = _load_45() if any(sp in column_species for sp in SPECIES_ORDER) else None
    # Free-style species are repelled together in ONE pass (not one call per species) -- a
    # per-species-only pass would resolve overlaps within each species but miss a label from one
    # species landing on a label from the other (caught by re-running check_layout, not by eye:
    # 'HLA-DRB1 overlaps KIR3DS1' in panel a, two close points from different species).
    free_x, free_y, free_lab, free_col = [], [], [], []
    for sp in SPECIES_ORDER:
        xs, ys, labs = per_species[sp]
        if not xs:
            continue
        if sp in column_species:
            ends = {lab: (x, y) for x, y, lab in zip(xs, ys, labs)}
            colors_by_label = {lab: SPECIES_COLOR[sp] for lab in labs}
            m45._label_curve_ends(ax, ends, colors_by_label, fontsize=column_fontsize,
                                   min_gap_frac=column_min_gap_frac)
        else:
            free_x += xs
            free_y += ys
            free_lab += labs
            free_col += [SPECIES_COLOR[sp]] * len(xs)
    if free_x:
        _label_points(ax, free_x, free_y, free_lab, free_col)
    if na_genes:
        note = na_note or (", ".join(na_genes) + ": no catalogue data at this level")
        # Top-left, not bottom-left: this data's y-range clusters low-to-mid (completeness
        # 0.08-0.57), so the top of the panel is the reliably empty corner -- checked against the
        # actual metrics range, not assumed.
        t = ax.text(0.02, 0.98, note, transform=ax.transAxes, fontsize=5.2, color="#666666",
                     ha="left", va="top", style="italic", wrap=True)
        vc.mark_label(t)
    return handles


def fig_catalogue_completeness(metrics, out_stem):
    """Three panels, left to right: (a) genomic (span-level) identity -- the UPPER-BOUND context
    view (see README caveat 1: person-specific span/UTR boundaries and intronic assembly noise the
    CDS-level artifact gate can't see likely inflate this one), (b) CDS-level completeness (a v4b
    headline granularity), (c) protein-level completeness (the other v4b headline granularity).
    Panel (b)'s x-axis reuses `pct_novel_any` (a genomic-level novelty measure) as the only novelty
    percentage 44 exports that pairs with the CDS s_obs/chao2 identity granularity -- CDS itself
    has no separate cds_novel/non-novel split in 44's tables, only its own identity-level
    s_obs/chao2 (see `build_gene_metrics` docstring)."""
    with vc.nature_style():
        fig, axes = plt.subplots(1, 3, figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(140)),
                                  constrained_layout=True)
        # Fix labels/limits BEFORE labeling points -- `_label_points()` measures real rendered
        # pixel positions, which shift if axis decorations are added afterward.
        axes[0].set_xlabel("% of distinct alleles novel (any level)")
        axes[0].set_ylabel(r"Chao2 richness completeness (S$_{obs}$/Chao2)")
        axes[0].set_title("genomic identity (upper bound)", fontsize=6.5)
        axes[1].set_xlabel("% of distinct alleles novel (any level)")
        axes[1].set_title("CDS identity (headline)", fontsize=6.5)
        axes[2].set_xlabel("% of distinct alleles novel (protein level)")
        axes[2].set_title("protein identity (headline)", fontsize=6.5)
        handles = _scatter_panel(axes[0], metrics, "pct_novel_any", "completeness")
        # CDS panel: both species' any-level-novelty values span a similar, overlapping range
        # here (unlike the protein panel), so gene labels are dropped for legibility (see
        # `_scatter_panel(..., label_genes=False)` docstring) -- values are in
        # `46_catalogue_metrics.tsv`/`_by_ancestry.tsv`, and the caption says so explicitly.
        _scatter_panel(axes[1], metrics, "pct_novel_any", "completeness_cds", label_genes=False)
        # Protein-level panel: KIR2DP1/KIR3DP1 (pseudogenes, no catalogued reference protein --
        # `kir_protein_catalogue_qc.tsv`) are NaN on both axes here and are named explicitly
        # rather than silently vanishing from the plot.
        _scatter_panel(axes[2], metrics, "pct_novel_protein", "completeness_protein",
                       na_note="KIR2DP1, KIR3DP1: pseudogenes, no catalogued reference protein "
                                "(excluded, not 0)", column_species=("kir",))
        axes[0].legend(handles=handles, loc="lower left", **vc.LEGEND_KW)
        vc.panel_letter(axes[0], "a")
        vc.panel_letter(axes[1], "b")
        vc.panel_letter(axes[2], "c")
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

    fig_catalogue_completeness(metrics, os.path.join(args.out_dir, "fig_catalogue_completeness"))
    fig_recurrence_composition(curve, os.path.join(args.out_dir, "fig_recurrence_composition"))

    print(f"[46] 2 figures + 46_catalogue_metrics.tsv written to {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
