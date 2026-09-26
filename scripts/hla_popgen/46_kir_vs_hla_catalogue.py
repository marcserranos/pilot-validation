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
    return rec, curve, cov, slope


def build_gene_metrics(cov, slope, ancestry="ALL"):
    """Returns a tidy DataFrame, one row per (gene, species): pct_novel_any, pct_novel_protein,
    good_turing_coverage, completeness (s_obs/chao2), chao2_undetected_f0hat, slope_per_1000, all
    at the genomic identity granularity pooled over `ancestry` (default ALL)."""
    sub = cov[cov["ancestry"] == ancestry].copy()
    for col in ("s_obs", "good_turing_coverage", "chao2", "chao2_undetected_f0hat"):
        sub[col] = pd.to_numeric(sub[col], errors="coerce")

    def s_obs_at(level):
        d = sub[sub["level"] == level][["gene", "species", "s_obs"]]
        return d.rename(columns={"s_obs": f"s_obs_{level}"})

    base = sub[sub["level"] == "genomic"][
        ["gene", "species", "s_obs", "good_turing_coverage", "chao2",
         "chao2_undetected_f0hat"]].rename(columns={"s_obs": "s_obs_genomic"})
    base = base.merge(s_obs_at("any_novel"), on=["gene", "species"], how="left")
    base = base.merge(s_obs_at("protein_novel"), on=["gene", "species"], how="left")

    base["pct_novel_any"] = 100.0 * base["s_obs_any_novel"] / base["s_obs_genomic"]
    base["pct_novel_protein"] = 100.0 * base["s_obs_protein_novel"] / base["s_obs_genomic"]
    base["completeness"] = base["s_obs_genomic"] / base["chao2"]

    sl = slope[(slope["ancestry"] == ancestry) & (slope["level"] == "genomic")
               & (slope["curve"] == "distinct")][["gene", "species", "slope_per_1000"]]
    base = base.merge(sl, on=["gene", "species"], how="left")

    base["gene_display"] = np.where(base["species"] == "hla", "HLA-" + base["gene"], base["gene"])
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
                push_px = overlap_y / 2.0 + 4.0
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


def fig_catalogue_completeness(metrics, out_stem):
    with vc.nature_style():
        fig, ax = plt.subplots(figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(115)),
                                constrained_layout=True)
        handles = []
        all_x, all_y, all_lab, all_col = [], [], [], []
        for sp in SPECIES_ORDER:
            s = metrics[metrics["species"] == sp].dropna(subset=["pct_novel_any", "completeness"])
            h = ax.scatter(s["pct_novel_any"], s["completeness"], s=16,
                            color=SPECIES_COLOR[sp], label=SPECIES_LABEL[sp], zorder=3,
                            edgecolors="white", linewidths=0.3)
            handles.append(h)
            all_x += s["pct_novel_any"].tolist()
            all_y += s["completeness"].tolist()
            all_lab += s["gene_display"].tolist()
            all_col += [SPECIES_COLOR[sp]] * len(s)
        # Fix axis limits and the legend BEFORE labeling -- `_label_points()` measures real
        # rendered pixel positions, which shift if xlim/ylim/legend are added afterward (each
        # changes the data<->pixel mapping or adds a sibling artist that reflows the layout).
        ax.set_xlabel("% of distinct alleles novel (any level)")
        ax.set_ylabel(r"Chao2 richness completeness (S$_{obs}$/Chao2)")
        ax.set_ylim(0, 1.05)
        ax.set_xlim(left=-2)
        ax.legend(handles=handles, loc="lower left", **vc.LEGEND_KW)
        # ALL points (both species) go through one repulsion pass -- resolving overlaps per
        # species separately would miss a KIR label landing on an HLA label placed nearby.
        _label_points(ax, all_x, all_y, all_lab, all_col)
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

    rec, curve, cov, slope = load_tables(args.in_dir)
    metrics = build_gene_metrics(cov, slope)
    metrics_path = os.path.join(args.out_dir, "46_catalogue_metrics.tsv")
    metrics.to_csv(metrics_path, sep="\t", index=False)
    print(f"  wrote {metrics_path}", file=sys.stderr)

    fig_catalogue_completeness(metrics, os.path.join(args.out_dir, "fig_catalogue_completeness"))
    fig_recurrence_composition(curve, os.path.join(args.out_dir, "fig_recurrence_composition"))

    print(f"[46] 2 figures + 46_catalogue_metrics.tsv written to {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
