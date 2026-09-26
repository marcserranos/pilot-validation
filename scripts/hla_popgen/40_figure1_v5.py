#!/usr/bin/env python3
"""Figure 1 v5 -- the introductory figure of the paper, redrawn to Cole's 2026-09-22 panel-by-panel
feedback (sprints/CALL_SUMMARY_2026-09-22.md Sec 6) and then to the orchestrator's 2026-09-23
review of the first v5 draft (fixed 183x150 mm canvas, censored-cell display bug, panel a rebuilt
from disclosure-safe bins, new panel f). Nature spec: double column 183 mm, <=~150 mm tall, 6 pt
body text, 7 pt axis titles, bold lowercase 8 pt panel letters, no in-panel caption text (the
legend lives in the README only).

Two-stage pipeline -- panel (a) needs a VM-side aggregation step of its own now (see
40a_admixture_bins.py), panel (b) needs a VM render of per-allele centroids (never per-person
data leaves the VM for either):

  Stage 0 (VM, `40a_admixture_bins.py`, separate script): bins unrelated people into >=20-person
  bins (sorted by predicted ancestry, then dominant-component probability) and exports per-bin
  MEAN admixture proportions -- disclosure-safe by construction, no bin can be smaller than the
  AoU small-cell floor. Writes panel_a_admixture_bins.tsv.

  Stage 1 (`--mode vm-panels`, VM): reuses 36_figure1_native.py's own `panel_b_data` (never
  re-derived) to export the per-allele ternary centroid table for panel (b). Only the aggregate
  table leaves the VM.

  Stage 2 (`--mode compose`, local): composes the full 183x150 mm figure from
    - panel_a_admixture_bins.tsv (binned admixture strip; placeholder box until the VM run lands)
    - panel_b_ternary_alleles.tsv (ternary, 0/50/100 ticks, light 25-step gridlines)
    - panel (c): gene x ancestry heatmap of the PROTEIN-LEVEL novelty rate (field 2 only -- same
      definition panel (d) uses, and the metric least inflated by non-coding diversity/artifacts;
      see README for why this was chosen over the any-field genomic-level rate 33/36 used).
      Every cell whose numerator has ANY censored (`<20`) contribution is drawn hatched and
      labelled as an upper-bound interval ('<=x%'), never as a bare point value -- this is the
      exact bug flagged in the 2026-09-23 review (MID x HLA-A rendered as a bogus '0').
    - panel (d): novel protein alleles per gene (same rows as c), stacked by recurrence class,
      as a marginal bar sharing c's y-axis. Recurrence-class COUNTS are allele-level, not
      participant-level, so they are not subject to the small-cell rule (established convention,
      36_figure1_native.py's own panel-d docstring) -- but a legend entry with zero total across
      the displayed genes is dropped rather than shown as a phantom category.
    - panel (e): per-ancestry allele-discovery curves from 39_saturation_by_ancestry, re-read at
      render time (safe against 39's concurrent extrapolation-fit correction -- only the raw
      curve points are plotted, not the Clench fit).
    - panel (f): DQ G1/G2 observed-vs-expected purge by ancestry, from 37_dq_g1g2_signed_ld's
      already-committed, disclosure-cleared oe_purge_committed.tsv. Chosen over a KIR preview
      because 41_kir_scoping has no aggregate data yet (scoping-only, pending the WS3 rerun),
      while the DQ O/E table is well powered (563-4,146 haplotypes/ancestry) and already public.

Usage (VM, run 40a first, then this):
    python3 scripts/hla_popgen/40a_admixture_bins.py --out ~/s03/results/40/panel_a_admixture_bins.tsv
    python3 scripts/hla_popgen/40_figure1_v5.py vm-panels --out-dir ~/s03/results/40

Usage (local):
    python3 scripts/hla_popgen/40_figure1_v5.py compose --pick A
"""
import argparse
import importlib.util
import json
import os
import sys

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)
import _viz_common as vc

_REPO_ROOT = os.path.dirname(os.path.dirname(_THIS_DIR))
DEFAULT_REPORTS = os.path.join(_REPO_ROOT, "reports", "hla_popgen")
DEFAULT_OUT_DIR = os.path.join(DEFAULT_REPORTS, "40_figure1_v5")

SUPPRESS_BELOW = 20
ANC = vc.ANCESTRY_ORDER  # ["AFR","AMR","EAS","EUR","MID","SAS"]
ANC_COLORS = vc.ANCESTRY_COLORS
CLASSICAL = ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"]
TERNARY = ["AFR", "EUR", "AMR"]

# Orchestrator spec: 6 pt body text, 7 pt axis titles, 8 pt bold panel letters.
FS_PANEL = 8
FS_TITLE = 7.0
FS_LABEL = 6.0
FS_TICK = 6.0
FS_LEG = 6.0
FS_ANNOT = 6.0

CANVAS_W_MM = 183.0
CANVAS_H_MM = 150.0


def _load(fn, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_THIS_DIR, fn))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def suppress(n):
    n = int(n)
    return "0" if n == 0 else ("<%d" % SUPPRESS_BELOW if n < SUPPRESS_BELOW else str(n))


# ===========================================================================================
# STAGE 1 -- VM only: panel (b) aggregate table
# ===========================================================================================
def run_vm_panels(args):
    """Panel (a) is now handled by 40a_admixture_bins.py (run separately). This mode only
    exports panel (b)'s per-allele ternary centroid table."""
    os.makedirs(args.out_dir, exist_ok=True)
    m36 = _load("36_figure1_native.py", "fig1_native_v36")
    m24 = _load("24_novelty_by_field.py", "novelty_by_field")

    print("[40 vm-panels] loading ...", flush=True)
    t1 = pd.read_csv(args.table1, sep="\t", dtype=str, low_memory=False)
    cohort = pd.read_csv(args.cohort_membership, sep="\t", dtype=str)

    people, n_removed = m24.build_people(
        sorted(set(t1["person_id"].astype(str))), args.cohort_membership,
        args.relatedness_table, args.kin_min, args.strict_threshold, args.skip_relatedness)
    keep = people[people["unrelated"]]
    anc_of = {p: a for p, a in zip(keep["person_id"].astype(str), keep["anc_strict"]) if a in ANC}
    people_keep = set(anc_of)
    print("[40] strict-ancestry (>=%.2f) unrelated=%d" % (args.strict_threshold, len(people_keep)),
          flush=True)

    pb = m36.panel_b_data(t1, cohort, people_keep, args.gene, args.min_carriers,
                          args.strict_threshold)
    print("[40] panel b: %d alleles (gene=%s)" % (len(pb), args.gene), flush=True)
    pb.drop(columns=["n_carriers"]).to_csv(
        os.path.join(args.out_dir, "panel_b_ternary_alleles.tsv"), sep="\t", index=False)
    with open(os.path.join(args.out_dir, "vm_panels_summary.json"), "w") as fh:
        json.dump({"n_people_strict_ancestry": len(people_keep), "n_alleles_panel_b": len(pb),
                   "gene_b": args.gene, "strict_threshold": args.strict_threshold,
                   "min_carriers": args.min_carriers, "n_removed_relatedness": n_removed}, fh,
                  indent=2)
    print("[40 vm-panels] done -> %s" % args.out_dir)


# ===========================================================================================
# STAGE 2 -- local: compose
# ===========================================================================================
def mm(x):
    return vc.mm(x)


# ---- panel a: binned admixture strip ------------------------------------------------------
def draw_panel_a(ax, bins_path):
    if not (bins_path and os.path.exists(bins_path)):
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.add_patch(__import__("matplotlib.pyplot", fromlist=["Rectangle"]).Rectangle(
            (0, 0), 1, 1, fill=False, edgecolor="#BBBBBB", lw=0.6, linestyle="--"))
        ax.text(0.5, 0.5, "panel a: awaiting 40a_admixture_bins.py VM run", ha="center",
               va="center", fontsize=FS_LABEL, color="#999999", transform=ax.transAxes)
        ax.axis("off")
        return None
    d = pd.read_csv(bins_path, sep="\t")
    d = d.sort_values(["anc", "order_within_anc"]).reset_index(drop=True)
    x = np.arange(len(d))
    bottom = np.zeros(len(d))
    for a in ANC:
        v = d["p_" + a.lower()].to_numpy(dtype=float)
        ax.fill_between(x, bottom, bottom + v, color=ANC_COLORS[a], linewidth=0)
        bottom += v
    codes = d["anc"].to_numpy()
    for i in range(1, len(codes)):
        if codes[i] != codes[i - 1]:
            ax.axvline(i, color="white", lw=0.4)
    ax.set_xlim(0, len(d))
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([0, 1.0])
    ax.set_ylabel("admixture\nproportion", fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    start = 0
    for i in range(1, len(codes) + 1):
        if i == len(codes) or codes[i] != codes[start]:
            mid = (start + i) / 2.0
            ax.text(mid, -0.10, codes[start], ha="center", va="top", fontsize=FS_TICK,
                   fontweight="bold", color=ANC_COLORS[codes[start]],
                   transform=ax.get_xaxis_transform())
            start = i
    ax.spines[["top", "right"]].set_visible(False)
    return d


# ---- panel b: ternary --------------------------------------------------------------------
def tern_xy(p):
    l, r, t = p
    return r + 0.5 * t, (np.sqrt(3) / 2.0) * t


def draw_ternary_grid(ax):
    """Light gridlines every 20% (orchestrator spec), numeric labels only at 0/50/100 per side,
    placed OUTSIDE the triangle so they never sit on top of a vertex's own bold corner label or
    on top of data points inside the triangle."""
    V = {"AFR": (1, 0, 0), "EUR": (0, 1, 0), "AMR": (0, 0, 1)}
    xy = {k: tern_xy(v) for k, v in V.items()}
    tri = np.array([xy["AFR"], xy["EUR"], xy["AMR"], xy["AFR"]])
    ax.plot(tri[:, 0], tri[:, 1], color="#333333", lw=0.6, zorder=2)
    for k in [20, 40, 60, 80]:
        kk = k / 100.0
        p0 = tern_xy((kk, 1 - kk, 0)); p1 = tern_xy((kk, 0, 1 - kk))
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#E8E8E8", lw=0.3, zorder=1)
        p0 = tern_xy((0, kk, 1 - kk)); p1 = tern_xy((1 - kk, kk, 0))
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#E8E8E8", lw=0.3, zorder=1)
        p0 = tern_xy((1 - kk, 0, kk)); p1 = tern_xy((0, 1 - kk, kk))
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#E8E8E8", lw=0.3, zorder=1)
    # Each side gets its own outward normal (perpendicular to that edge, away from the centroid),
    # distinct from the direction any vertex's corner label sits in -- so the 0/50/100 numerals
    # along an edge never stack on a corner label even where a tick coincides with a vertex.
    centroid = np.mean(tri[:3], axis=0)
    edges = [
        ("AFR", xy["AFR"], xy["EUR"], lambda kk: tern_xy((kk, 1 - kk, 0))),
        ("EUR", xy["EUR"], xy["AMR"], lambda kk: tern_xy((0, kk, 1 - kk))),
        ("AMR", xy["AMR"], xy["AFR"], lambda kk: tern_xy((1 - kk, 0, kk))),
    ]
    for _, p_a, p_b, pt_fn in edges:
        p_a, p_b = np.array(p_a), np.array(p_b)
        edge_vec = p_b - p_a
        normal = np.array([edge_vec[1], -edge_vec[0]])
        normal = normal / (np.linalg.norm(normal) + 1e-9)
        mid = (p_a + p_b) / 2.0
        if np.dot(normal, mid - centroid) < 0:
            normal = -normal
        # Only the midpoint (50) is labelled -- S04 WS-C phase 2 fix: the 0/100 ticks sit right on
        # top of a vertex, which is also where that vertex's own bold ancestry label (AFR/EUR/AMR)
        # is anchored, so "100"/"0" collided with "AFR"/"AMR" (`check_layout()` text_overlap).
        # The vertex label already states the 0/100 endpoint unambiguously; the numeral added
        # nothing but a collision.
        tp = np.array(pt_fn(0.5))
        lp = tp + normal * 0.075
        ax.text(lp[0], lp[1], "50", fontsize=FS_ANNOT, ha="center", va="center", color="#888888")
    return xy


def draw_panel_b(ax, pb, gene_b, min_carriers, strict):
    xy = draw_ternary_grid(ax)
    if not pb.empty:
        n = pb["n_carriers_disp"].apply(lambda s: 19 if str(s).startswith("<") else float(s))
        n = np.maximum(n.to_numpy(dtype=float), 1.0)
        lo, hi = np.log10(n.min()), np.log10(n.max())
        rng = max(1e-9, hi - lo)
        alpha = 0.25 + 0.65 * (np.log10(n) - lo) / rng
        size = 5 + 26 * (np.log10(n) - lo) / rng
        for (_, r), al, sz in zip(pb.iterrows(), alpha, size):
            xyp = tern_xy((r[TERNARY[0]], r[TERNARY[1]], r[TERNARY[2]]))
            novel = r["frac_novel_calls"] > 0.5
            dom = max(TERNARY, key=lambda a: r[a])
            ax.scatter([xyp[0]], [xyp[1]], s=sz, marker="^" if novel else "o",
                      facecolor=ANC_COLORS[dom], edgecolor="none", alpha=float(al), zorder=3)
    for lab, pos, ha, va, off in [("AFR", xy["AFR"], "right", "top", (-4, -6)),
                                  ("EUR", xy["EUR"], "left", "top", (4, -6)),
                                  ("AMR", xy["AMR"], "center", "bottom", (0, 6))]:
        ax.annotate(lab, pos, fontsize=FS_LABEL, fontweight="bold", ha=ha, va=va,
                   xytext=off, textcoords="offset points", color=ANC_COLORS[lab])
    ax.set_xlim(-0.16, 1.16)
    ax.set_ylim(-0.16, np.sqrt(3) / 2 + 0.12)
    ax.axis("off")
    # `axis("off")` hides the CURRENT axis decorations but does not survive a later
    # `fig.canvas.draw()` in this matplotlib version -- a redraw (e.g. `check_layout()`'s or
    # `save_fig()`'s) can still invoke this axes' default Locator/Formatter and grow its x/y
    # majorTicks from 1 to several, rendering plain default-formatted numbers ("0.0", "1.0", ...)
    # from its raw xlim/ylim (-0.16 to ~1.16) directly on top of this panel's own "catalogued"/
    # "first observed here" legend text (found by the orchestrator's re-review; same bug class as
    # `37c_dq_g1g2_from_committed.py`'s Ticker-replacement fix and Figure 1 v5's own ax_d fix --
    # see either docstring for the full writeup). Null out both axes explicitly so there is
    # nothing for a later draw to (re)populate.
    from matplotlib.ticker import NullLocator, NullFormatter
    for _axis in (ax.xaxis, ax.yaxis):
        _axis.set_major_locator(NullLocator())
        _axis.set_major_formatter(NullFormatter())
    from matplotlib.lines import Line2D
    # Marker legend anchored INSIDE the panel's own axes-fraction box (S04 WS-C phase 2 fix: the
    # previous bbox_to_anchor=(1.14, -0.08) sat outside the axes entirely, in the gap toward panel
    # c, where it collided with panel c's rotated y-axis label -- `check_layout()`'s
    # text_overlap check caught it, and check (a3)/text_foreign_spine would have too. (1.0, 0.0)
    # keeps the legend fully inside this axes' own box, tucked in the empty lower-right corner of
    # the triangle (between the AFR-AMR edge and the plot's own right boundary).
    # Below the triangle's own base edge (data only ever falls inside the triangle, so the
    # padding band below y=0 is guaranteed empty) rather than inside the plotted region -- the
    # (1.0, 0.34) placement passed `check_layout()` but sat visually on top of scattered data
    # points, which a Nature editor would still flag even though the mechanical linter can't see
    # marker/text overlap. ncol=2 keeps it to one compact row.
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", color="#666666", ms=2.4,
                              label="catalogued"),
                       Line2D([], [], marker="^", ls="none", color="#666666", ms=2.4,
                              label="first observed here")],
             frameon=False, fontsize=FS_LEG - 0.5, loc="upper center", handletextpad=0.3,
             borderaxespad=0.0, labelspacing=0.25, columnspacing=1.0, ncol=2,
             bbox_to_anchor=(0.5, -0.02))


# ---- panel c/d: novelty heatmap + recurrence marginal bar --------------------------------
# Short, single-line direct labels (de-AI checklist item 7/13: no sentence-length axis labels;
# state the unit, let the README carry the full metric definition and caveats). S04 WS-C phase 2
# fix: the previous two-line labels, rotated 90 for a y-axis label, laid their second line out
# WIDTH-wise (perpendicular to the vertical reading direction) rather than adding height -- that
# extra horizontal thickness is exactly what reached left into panel b's space and collided with
# the "EUR" corner label and the ternary legend (`check_layout()` text_overlap). A single line has
# no second line to add that width.
# Colorbar label (orchestrator review, 2026-09-26: panel c's y-axis label was removed -- its ticks
# are gene names, not a "novelty" scale -- and the quantity the heatmap color encodes moved onto
# the colorbar itself, which used to just say "% (exact)", not naming what was 100% of). One clear
# line, unit stated.
_METRIC_CBAR_LABEL = {
    "any_field": "any-field novelty, % of calls",
    "cds": "novel CDS, % of calls",
    "protein": "novel protein allele, % of calls",
}


def load_panel_cd(field_counts_path, clusters_path, clusters_totals_path=None,
                  ancestry_scheme="strict", c_metric="auto"):
    """Prefers a direct VM export of the TOTALS (panel_c_novelty_totals.tsv, from
    40b_novelty_rate_export.py) when present -- exact counts, suppressed exactly once, so a cell
    is hatched only when it is genuinely <20, not because it is the sum of several independently-
    censored sub-cells (the 2026-09-23 fix that made panel c almost entirely hatched). Falls back
    to the older, more-heavily-censored rendering built from 24/33's committed sub-split table
    when the totals export hasn't been run yet.

    Returns (grid, metric, source) where grid is a dict of gene x ancestry DataFrames: value (%,
    exact where not hatched), hatched (bool), text (str already formatted for the cell).
    """
    m34 = _load("34_novel_protein_recurrence.py", "novel_recurrence")
    cl = pd.read_csv(clusters_path, sep="\t", dtype=str)
    prot_cl = cl[cl["cluster_type"] == "novel_protein"]
    prot_cl = m34.classify(prot_cl)
    dtab = m34.gene_table(prot_cl).reindex(CLASSICAL).fillna(0)

    if clusters_totals_path and os.path.exists(clusters_totals_path):
        t = pd.read_csv(clusters_totals_path, sep="\t", dtype=str)
        metric = c_metric
        if metric == "auto":
            frac_ge20 = {}
            for key in ("protein", "cds", "any_field"):
                col = "n_novel_%s" % key
                if col in t.columns:
                    frac_ge20[key] = float((~t[col].astype(str).str.startswith("<")).mean())
            # Prefer the most specific (protein) metric that is usably powered (>=50% of cells
            # clear 20); fall back toward the broader metrics otherwise.
            metric = next((k for k in ("protein", "cds", "any_field")
                          if frac_ge20.get(k, 0) >= 0.5), "any_field")
        value = pd.DataFrame(index=CLASSICAL, columns=ANC, dtype=float)
        hatched = pd.DataFrame(index=CLASSICAL, columns=ANC, dtype=bool).fillna(True)
        text = pd.DataFrame(index=CLASSICAL, columns=ANC, dtype=object)
        for _, r in t.iterrows():
            g, a = r["gene"], r["ancestry"]
            if g not in CLASSICAL or a not in ANC:
                continue
            n_str = str(r.get("n_novel_%s" % metric, "<20"))
            n_total_str = str(r.get("n_total", "<20"))
            is_cens = n_str.startswith("<") or n_total_str.startswith("<")
            hatched.loc[g, a] = is_cens
            if is_cens:
                text.loc[g, a] = n_str if n_str.startswith("<") else "n/a"
                value.loc[g, a] = np.nan
            else:
                rate = r.get("rate_%s_pct" % metric, "")
                v = float(rate) if rate not in ("", None) and not pd.isna(rate) else np.nan
                value.loc[g, a] = v
                text.loc[g, a] = "%.0f" % v if not np.isnan(v) else n_str
        return {"value": value, "hatched": hatched, "text": text}, metric, "totals_export", dtab

    # ---- fallback: old, sub-split-and-summed rendering (heavier hatching) -------------------
    m33 = _load("33_figure1_v3_compose.py", "fig1v3_compose")
    counts = pd.read_csv(field_counts_path, sep="\t")
    rates = m33.rate_by_gene_ancestry(counts, ancestry_scheme, True, CLASSICAL)
    fields = rates[rates["field_class"] != "artifact_control"].copy()
    lo = fields.pivot_table(index="gene", columns="ancestry", values="pct_lower",
                            aggfunc="sum").reindex(index=CLASSICAL, columns=ANC)
    hi = fields.pivot_table(index="gene", columns="ancestry", values="pct_upper",
                            aggfunc="sum").reindex(index=CLASSICAL, columns=ANC)
    ncens = fields.pivot_table(index="gene", columns="ancestry", values="n_censored_cells",
                               aggfunc="sum").reindex(index=CLASSICAL, columns=ANC).fillna(0)
    called_hi = fields.pivot_table(index="gene", columns="ancestry", values="n_called_upper",
                                   aggfunc="max").reindex(index=CLASSICAL, columns=ANC)
    value = pd.DataFrame(index=CLASSICAL, columns=ANC, dtype=float)
    hatched = pd.DataFrame(index=CLASSICAL, columns=ANC, dtype=bool)
    text = pd.DataFrame(index=CLASSICAL, columns=ANC, dtype=object)
    for g in CLASSICAL:
        for a in ANC:
            ch = called_hi.loc[g, a] if (g in called_hi.index and a in called_hi.columns) else 0
            nc = ncens.loc[g, a] if (g in ncens.index and a in ncens.columns) else 0
            l = lo.loc[g, a] if (g in lo.index and a in lo.columns) else np.nan
            h = hi.loc[g, a] if (g in hi.index and a in hi.columns) else np.nan
            if pd.isna(l) or (ch is not None and ch < SUPPRESS_BELOW):
                hatched.loc[g, a] = True; value.loc[g, a] = np.nan; text.loc[g, a] = "n/a"
            elif nc and nc > 0:
                hatched.loc[g, a] = True; value.loc[g, a] = l; text.loc[g, a] = "≤%.0f" % h
            else:
                hatched.loc[g, a] = False; value.loc[g, a] = l; text.loc[g, a] = "%.0f" % l
    return {"value": value, "hatched": hatched, "text": text}, "any_field", "sub_split_fallback", dtab


def draw_panel_cd(ax_hm, ax_bar, grid, metric, dtab, cax=None):
    genes, ancs = CLASSICAL, ANC
    value = grid["value"].reindex(index=genes, columns=ancs)
    hatched = grid["hatched"].reindex(index=genes, columns=ancs)
    text = grid["text"].reindex(index=genes, columns=ancs)
    data = value.to_numpy(dtype=float)
    vmax = np.nanmax(data)
    vmax = vmax if vmax and vmax > 0 else 1.0
    im = ax_hm.imshow(data, cmap="YlOrRd", aspect="auto", vmin=0, vmax=vmax)
    for i, g in enumerate(genes):
        for j, a in enumerate(ancs):
            if bool(hatched.loc[g, a]):
                vc.hatch_suppressed(ax_hm, j - 0.5, i - 0.5, 1, 1, alpha=0.55)
                ax_hm.text(j, i, str(text.loc[g, a]), ha="center", va="center",
                          fontsize=FS_ANNOT - 0.5, color="#555555")
            else:
                ax_hm.text(j, i, str(text.loc[g, a]), ha="center", va="center",
                          fontsize=FS_ANNOT, color="#222222")
    ax_hm.set_xticks(range(len(ancs)))
    ax_hm.set_xticklabels(ancs, fontsize=FS_TICK)
    for tick, a in zip(ax_hm.get_xticklabels(), ancs):
        tick.set_color(ANC_COLORS[a])
        tick.set_fontweight("bold")
    ax_hm.set_yticks(range(len(genes)))
    ax_hm.set_yticklabels([g.replace("HLA-", "") for g in genes], fontsize=FS_TICK)
    # No y-axis label: the previous "any-field novelty (%)" label described the COLOR, not the
    # axis (whose ticks are gene names) -- a reader reads a y-axis label as "what varies down this
    # axis", which here is just "gene". The quantity the color encodes belongs on the colorbar,
    # which now carries it directly (orchestrator review, 2026-09-26).
    for sp in ax_hm.spines.values():
        sp.set_visible(False)
    if cax is not None:
        cb = ax_hm.figure.colorbar(im, cax=cax, orientation="vertical")
        cb.ax.tick_params(labelsize=FS_ANNOT, length=2)
        cb.set_label(_METRIC_CBAR_LABEL.get(metric, metric), fontsize=FS_ANNOT, labelpad=2)

    cols_all = [("seen once", "#C9CFD6"), ("2–19 unrelated people", "#5B8FBF"),
               ("≥20 unrelated people", "#B4472E")]
    cols = [(k, c) for k, c in cols_all
           if k in dtab.columns and float(dtab.reindex(genes)[k].sum()) > 0]
    y = np.arange(len(genes))
    left = np.zeros(len(genes))
    for key, col in cols:
        v = dtab.reindex(genes)[key].to_numpy(dtype=float)
        ax_bar.barh(y, v, left=left, height=0.72, color=col, label=key, edgecolor="white",
                   linewidth=0.3, zorder=3)
        left += v
    # No set_yticks/set_yticklabels([]) here: ax_bar's y-axis Ticker is a PRIVATE NullLocator/
    # NullFormatter (set at creation in compose_layout(), since sharey= makes the Ticker object
    # shared with ax_hm otherwise) -- calling set_yticklabels([]) on it would be a no-op for
    # display but would also needlessly reinstall a FixedFormatter, undoing that protection.
    ax_bar.tick_params(axis="y", left=False, labelleft=False)
    ax_bar.set_ylim(-0.5, len(genes) - 0.5)
    ax_bar.invert_yaxis()
    ax_hm.invert_yaxis()
    # Structural row alignment (orchestrator review, 2026-09-26: "make it visibly aligned to c's
    # rows") -- `ax_bar` is created with `sharey=ax_hm` in `compose_layout()`, which ties the two
    # axes' y data-to-display mapping exactly; `mark_marginal()` makes that guarantee mechanically
    # checked by `check_layout()` rather than merely asserted in a comment.
    vc.mark_marginal(ax_bar, ax_hm, axis="y")
    ax_bar.set_xlabel("novel protein alleles", fontsize=FS_TITLE, labelpad=2)
    ax_bar.tick_params(axis="x", labelsize=FS_TICK)
    ax_bar.spines[["top", "right"]].set_visible(False)
    # No in-panel legend (orchestrator review, 2026-09-26, 4th attempt -- moved to the report
    # README/caption instead): below-axis placement didn't physically fit this row (measured: the
    # gap needed ~9% of figure height, had ~7.5%); the bottom-right in-panel corner (tried next)
    # sat on top of the HLA-A bar, which is only short RELATIVE to the longest bars (B, C), not
    # short enough in absolute terms to clear a 2-entry legend box. `check_layout()` never caught
    # this because it has no check comparing legend/text Patches against ordinary data Patches
    # (bars) in general -- only against artists explicitly registered as decorations
    # (`mark_decoration()`) or plotted lines a label is registered against (`mark_label()`); a
    # blanket text-vs-every-Patch check would flag legitimate in-bar value labels elsewhere in
    # this codebase, so it was deliberately never added. The colour legend (grey="seen once",
    # blue="2-19 unrelated people") is stated in this report's README instead.


# ---- panel e: discovery curves --------------------------------------------------------------
def draw_panel_e(ax, curves_path, scheme="pred"):
    """Redesigned (orchestrator review, 2026-09-26, two passes):

    Pass 1 fixed a naive other-label-only repulsion (labels never checked against the actual
    LINES): "AMR" landed on top of AFR's still-rising curve (AFR extends further in x, so its line
    is still there, near its own plateau, at AMR's shorter endpoint), and AFR's own label, pushed
    up to clear AMR's, cleared the autoscaled ylim top and was silently clipped (no AFR label at
    all -- autoscale only ever saw the DATA, never the text).

    Pass 2: even with per-curve, line-aware repulsion, MID's label still landed on EAS/EUR/SAS's
    lines. Root cause is structural, not a tuning problem: MID's cohort is the smallest (N=487),
    so its curve ends VERY early in x, in the region where all six discovery curves are still
    close together (early cohort growth looks similar across ancestries before they diverge at
    larger N) -- there is no y position near MID's own endpoint x that clears 3 other lines at
    once. Fix: label all six at a SHARED x just past the LONGEST curve's endpoint (here AFR, the
    largest cohort) rather than each at its own endpoint -- no line is drawn beyond its own last
    point, so a shared x past the rightmost one guarantees zero line collisions by construction,
    at the cost of a short implicit "leader" gap between a shorter curve's true end and its label
    (standard practice for this style of chart, e.g. an Economist-style end-of-line legend
    column). Labels are still ordered/spaced by each curve's own final value.
    """
    df = pd.read_csv(curves_path, sep="\t")
    d = df[(df["scheme"] == scheme) & (df["gene_group"] == "classical_pooled")
          & (df["category"] == "all")].copy()
    curves_xy, ends = {}, {}
    xmax = 0.0
    for a in ANC:
        sub = d[d["ancestry"] == a].sort_values("n")
        if sub.empty:
            continue
        x = sub["n"].to_numpy(dtype=float)
        y = sub["mean_distinct"].to_numpy(dtype=float)
        lo_ = sub["lo2_5"].to_numpy(dtype=float)
        hi_ = sub["hi97_5"].to_numpy(dtype=float)
        ax.plot(x, y, color=ANC_COLORS[a], lw=0.9, zorder=3, label=a)
        ax.fill_between(x, lo_, hi_, color=ANC_COLORS[a], alpha=0.15, linewidth=0, zorder=2)
        curves_xy[a] = (x, y)
        ends[a] = (x[-1], y[-1])
        xmax = max(xmax, x[-1])
    if ends:
        label_x_shared = xmax * 1.03  # past every curve's own last point -- no line reaches here
        label_y = {a: ends[a][1] for a in ends}
        yrange = max(v[1] for v in ends.values()) - min(v[1] for v in ends.values())
        # 0.12 -> 0.16 (orchestrator, 2026-09-26): AFR (511) and AMR (503) start only 8 apart,
        # and 0.12*yrange wasn't quite enough separation once real font metrics were accounted
        # for -- they still touched.
        min_gap = max(1.0, yrange * 0.16)

        order = sorted(ends, key=lambda a: label_y[a])
        placed = []
        for a in order:
            y = label_y[a]
            if placed and y - placed[-1] < min_gap:
                y = placed[-1] + min_gap
            placed.append(y)
        for a, y in zip(order, placed):
            label_y[a] = y

        for a in ends:
            t = ax.annotate(a, (label_x_shared, label_y[a]), fontsize=FS_ANNOT,
                           color=ANC_COLORS[a], fontweight="bold", va="center", ha="left")
            vc.mark_label(t)
            # Short leader dash from the curve's true endpoint to the shared label column, for any
            # ancestry whose own end sits noticeably left of it (else the label would otherwise
            # look unconnected to its curve) -- a plain thin line, not requiring its own label, so
            # it is registered as a decoration rather than data the (a4) check would compare labels
            # against.
            x0, y0 = ends[a]
            leader_x_end = label_x_shared - 0.01 * xmax  # stop short of the label's own text bbox
            if leader_x_end - x0 > 0.01 * xmax:
                # Horizontal-only leader (orchestrator, 2026-09-26: the previous diagonal from the
                # curve's true endpoint (x0, y0) straight to the label's (possibly repulsion-
                # shifted) y looked like an odd, ungrounded diagonal, especially for AFR). Drawn
                # as an "elbow": a short vertical tick right at the curve's own end (from y0 up/
                # down to the label's row), then a horizontal run at the label's own y out to the
                # label -- the dominant, eye-tracing segment is horizontal, and the vertical part
                # is short enough to read as "this curve's end connects to this row."
                ln, = ax.plot([x0, x0], [y0, label_y[a]], color=ANC_COLORS[a], lw=0.5,
                             ls=(0, (1, 1)), zorder=2.5, clip_on=False)
                vc.mark_decoration(ln)
                ln2, = ax.plot([x0, leader_x_end], [label_y[a], label_y[a]], color=ANC_COLORS[a],
                              lw=0.5, ls=(0, (1, 1)), zorder=2.5, clip_on=False)
                vc.mark_decoration(ln2)

        # Extend the view to fit both the longest curve's x and every placed label -- autoscale
        # only ever accounted for the DATA (curves/fill_between), never these text annotations, so
        # a label pushed up by repulsion could land above the autoscaled top and be clipped
        # invisibly (the actual root cause of the very first pass's missing AFR label).
        y_lo_data = 0.0
        y_hi_data = max(v[1] for v in ends.values())
        y_lo = min([y_lo_data] + list(label_y.values()))
        y_hi = max([y_hi_data] + list(label_y.values()))
        pad = max(1.0, (y_hi - y_lo) * 0.04)
        ax.set_ylim(y_lo - pad * 0.3, y_hi + pad)
        # 1.14 -> 1.02 (orchestrator, 2026-09-26): the wider multiplier left a large empty band of
        # x-range (out to ~3500) with nothing in it -- the labels/leaders all fit within a much
        # smaller margin past label_x_shared.
        ax.set_xlim(0, label_x_shared * 1.02)
    ax.set_xlabel("people sampled (both haplotypes)", fontsize=FS_TITLE)
    ax.set_ylabel("distinct HLA protein alleles\n(8 classical genes, pooled)", fontsize=FS_TITLE)
    ax.tick_params(labelsize=FS_TICK)
    ax.spines[["top", "right"]].set_visible(False)


# ---- panel f: DQ G1/G2 observed-vs-expected purge -------------------------------------------
def draw_panel_f(ax, oe_path):
    d = pd.read_csv(oe_path, sep="\t")
    d = d[d["ancestry"].isin(ANC)].set_index("ancestry").reindex(ANC)
    x = np.arange(len(ANC))
    w = 0.34
    exp_v = d["expected"].to_numpy(dtype=float)
    obs_hi = d["observed_upper"].to_numpy(dtype=float)
    # observed_lower == observed_upper == 0 for every ancestry (n_cells_censored_lt20 is 0
    # throughout -- see 37_dq_g1g2_signed_ld/oe_purge_committed.tsv) -- this is an EXACT zero,
    # not a suppressed/upper-bound one, so it is labelled and coloured accordingly.
    ax.bar(x - w / 2, exp_v, width=w, color="#B0B0B0", label="expected", zorder=3)
    # No legend entry for the "observed" bars (orchestrator review, 2026-09-26): every observed
    # bar is ~0.6% of its expected bar's height -- visually indistinguishable from "invisible" at
    # this scale -- so a legend swatch for it is a color a reader can never actually match against
    # the plot. The per-ancestry coloured "0" already states observed=0 directly, once per bar;
    # `label="_nolegend_"` keeps the bars themselves (needed for the "0" annotations' x-position
    # and as a visual placeholder next to "expected") without a matching, unreadable legend key.
    ax.bar(x + w / 2, np.maximum(obs_hi, exp_v * 0.006), width=w, color="#DDDDDD",
          edgecolor="#999999", linewidth=0.4, label="_nolegend_", zorder=3)
    for i, a in enumerate(ANC):
        ax.text(i + w / 2, max(exp_v) * 0.02, "0", ha="center", va="bottom",
               fontsize=FS_ANNOT, color=ANC_COLORS[a], fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(ANC, fontsize=FS_TICK)
    for tick, a in zip(ax.get_xticklabels(), ANC):
        tick.set_color(ANC_COLORS[a])
        tick.set_fontweight("bold")
    # Shortened from "DQA1~DQB1 cis haplotypes\nin G1x G2 cross-group cells" (orchestrator review,
    # 2026-09-26) -- full definition stays in the report README.
    ax.set_ylabel("cross-group DQA1×DQB1\nhaplotypes", fontsize=FS_TITLE)
    ax.tick_params(axis="y", labelsize=FS_TICK)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=FS_LEG, loc="upper right", handlelength=1.0,
             title="observed = 0 in all ancestries", title_fontsize=FS_LEG - 0.5,
             alignment="left")


def compose_layout(layout, bins_path, pb, cd_data, curves_path, oe_path, gene_b, min_carriers,
                   strict):
    import matplotlib
    matplotlib.use("Agg")
    with vc.nature_style():
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        W, H = mm(CANVAS_W_MM), mm(CANVAS_H_MM)
        fig = plt.figure(figsize=(W, H))
        # Reduced from 0.85 -- most of that gap was pure whitespace below panel a's thin strip
        # (orchestrator: reduce whitespace between rows 1 and 2). Round 2 (orchestrator,
        # 2026-09-26): "tighten whitespace between the rows" -- 0.42 -> 0.30.
        outer = GridSpec(3, 1, height_ratios=[22, 60, 55], hspace=0.30, left=0.050, right=0.985,
                         top=0.98, bottom=0.075, figure=fig)

        ax_a = fig.add_subplot(outer[0])
        draw_panel_a(ax_a, bins_path)
        vc.panel_letter(ax_a, "a", dx=-0.028, dy=1.12)

        # Column order b / c / colorbar / d (was b / c / d / colorbar): orchestrator review,
        # 2026-09-26 -- the colorbar describes panel c's heatmap, but sitting AFTER panel d put it
        # far from c with a wide, empty-looking gap (d's own bars don't fill their whole column,
        # so the visual gap between d's plotted bars and the colorbar's thin column read as
        # "floating"). Placing it directly against c's right edge removes that gap and reads as
        # what it is -- c's colorbar, not d's.
        if layout == "A":
            gs2 = outer[1].subgridspec(1, 4, width_ratios=[55, 78, 4, 32], wspace=0.55)
        else:
            # layout B: swap b/c+d emphasis -- give the heatmap a touch more width and the
            # ternary a touch less, to see whether it reads better with 8 gene rows.
            gs2 = outer[1].subgridspec(1, 4, width_ratios=[50, 84, 4, 32], wspace=0.55)
        ax_b = fig.add_subplot(gs2[0, 0])
        ax_c = fig.add_subplot(gs2[0, 1])
        cax = fig.add_subplot(gs2[0, 2])
        # sharey=ax_c: structural row alignment for panel d (orchestrator: "make it visibly
        # aligned to c's rows") -- ties the y data-to-display mapping exactly, rather than relying
        # on both axes independently being told the same ylim/invert calls to stay in sync.
        ax_d = fig.add_subplot(gs2[0, 3], sharey=ax_c)
        # `sharey=` makes ax_d.yaxis.major (locator+formatter) the SAME OBJECT as ax_c.yaxis.major
        # (matplotlib's actual sharex/sharey behavior -- see _viz_common.mark_decoration()'s
        # docstring and 37c_dq_g1g2_from_committed.py's Ticker-replacement fix for the full
        # writeup of the bug class this avoids): ax_d never needs its own y tick labels (the gene
        # names live on ax_c only), so give it a private NullLocator/NullFormatter on the shared
        # axis right away -- this keeps the alignment (which only depends on shared ylim via the
        # Grouper, not on the Locator) while preventing a later `fig.canvas.draw()` from growing
        # ax_d's own tick-label objects and rendering a stray duplicate copy of the gene names.
        from matplotlib.axis import Ticker as _Ticker
        from matplotlib.ticker import NullLocator as _NullLocator, NullFormatter as _NullFormatter
        ax_d.yaxis.major = _Ticker()
        ax_d.yaxis.set_major_locator(_NullLocator())
        ax_d.yaxis.set_major_formatter(_NullFormatter())

        draw_panel_b(ax_b, pb, gene_b, min_carriers, strict)
        vc.panel_letter(ax_b, "b", dx=-0.10, dy=1.06)

        grid, metric, dtab = cd_data
        draw_panel_cd(ax_c, ax_d, grid, metric, dtab, cax=cax)
        vc.panel_letter(ax_c, "c", dx=-0.34, dy=1.05)
        vc.panel_letter(ax_d, "d", dx=-0.12, dy=1.05)
        # shrink the colorbar to ~3 x 30 mm, vertically centered on the heatmap, instead of the
        # full 62 mm row height a plain gridspec cell would give it.
        fig.canvas.draw()
        p_hm = ax_c.get_position()
        p_cax = cax.get_position()
        target_h = mm(30) / H
        cy = (p_hm.y0 + p_hm.y1) / 2.0
        cax.set_position([p_cax.x0, cy - target_h / 2, min(p_cax.width, mm(3) / W), target_h])

        gs3 = outer[2].subgridspec(1, 2, width_ratios=[1, 1], wspace=0.55)
        ax_e = fig.add_subplot(gs3[0, 0])
        ax_f = fig.add_subplot(gs3[0, 1])
        draw_panel_e(ax_e, curves_path)
        vc.panel_letter(ax_e, "e", dx=-0.15, dy=1.05)
        draw_panel_f(ax_f, oe_path)
        vc.panel_letter(ax_f, "f", dx=-0.16, dy=1.05)

        # Explicit edge alignment (orchestrator: "align the left edges of b/e and right edges of
        # d/f") -- b/c/d/cax and e/f come from two INDEPENDENT subgridspecs (gs2, gs3), each
        # spanning the outer row's full width on its own, so their fractional column widths need
        # not put a shared boundary at the same physical x. Read the actual drawn positions and
        # nudge e's left edge to match b's, and f's right edge to match d's, preserving each axes'
        # own width (same pattern already used above to reposition the colorbar).
        fig.canvas.draw()
        pb_, pe_ = ax_b.get_position(), ax_e.get_position()
        if abs(pe_.x0 - pb_.x0) > 1e-6:
            ax_e.set_position([pb_.x0, pe_.y0, pe_.width, pe_.height])
        pd_, pf_ = ax_d.get_position(), ax_f.get_position()
        target_right = pd_.x1
        if abs(pf_.x1 - target_right) > 1e-6:
            ax_f.set_position([target_right - pf_.width, pf_.y0, pf_.width, pf_.height])

        return fig


def run_compose(args):
    os.makedirs(args.out_dir, exist_ok=True)
    pb = pd.read_csv(args.panel_b_table, sep="\t")
    grid, metric, source, dtab = load_panel_cd(args.field_counts, args.clusters,
                                               args.panel_c_totals, args.ancestry_scheme,
                                               args.c_metric)
    print("[40 compose] panel c source=%s metric=%s" % (source, metric))
    cd_data = (grid, metric, dtab)

    figs = {}
    for layout in ["A", "B"]:
        fig = compose_layout(layout, args.panel_a_bins, pb, cd_data, args.curves, args.oe_table,
                             args.gene, args.min_carriers, args.strict_threshold)
        stem = os.path.join(args.out_dir, "layout_%s" % layout)
        vc.save_fig(fig, stem)
        figs[layout] = stem + ".png"

    chosen = args.pick
    import shutil
    shutil.copy(figs[chosen], os.path.join(args.out_dir, "figure1_v5.png"))
    shutil.copy(figs[chosen].replace(".png", ".pdf"),
               os.path.join(args.out_dir, "figure1_v5.pdf"))

    with open(os.path.join(args.out_dir, "compose_summary.json"), "w") as fh:
        json.dump({"layouts_built": ["A", "B"], "chosen": chosen, "gene_b": args.gene,
                   "min_carriers_b": args.min_carriers, "strict_threshold": args.strict_threshold,
                   "ancestry_scheme_cd": args.ancestry_scheme,
                   "panel_c_source": source, "panel_c_metric": metric,
                   "panel_a_bins_present": bool(args.panel_a_bins and
                                                os.path.exists(args.panel_a_bins)),
                   "canvas_mm": [CANVAS_W_MM, CANVAS_H_MM]}, fh, indent=2)
    print("[40 compose] wrote layouts A/B, chose %s -> %s" % (chosen, args.out_dir))


# ===========================================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)

    vmp = sub.add_parser("vm-panels")
    vmp.add_argument("--table1", default=os.path.expanduser("~/pipeline_outputs/hla_calls_rich.tsv"))
    vmp.add_argument("--cohort-membership",
                     default=os.path.expanduser("~/pipeline_outputs/cohort_membership.tsv"))
    vmp.add_argument("--relatedness-table", default=os.path.expanduser(
        "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/"
        "relatedness/samples_relatedness.tsv"))
    vmp.add_argument("--skip-relatedness", action="store_true")
    vmp.add_argument("--kin-min", type=float, default=0.0442)
    vmp.add_argument("--strict-threshold", type=float, default=0.98)
    vmp.add_argument("--min-carriers", type=int, default=20)
    vmp.add_argument("--gene", default="HLA-B")
    vmp.add_argument("--out-dir", default=os.path.expanduser("~/s03/results/40"))

    cmp_ = sub.add_parser("compose")
    cmp_.add_argument("--panel-a-bins",
                      default=os.path.join(DEFAULT_OUT_DIR, "panel_a_admixture_bins.tsv"))
    cmp_.add_argument("--panel-b-table",
                      default=os.path.join(DEFAULT_OUT_DIR, "panel_b_ternary_alleles.tsv"))
    cmp_.add_argument("--field-counts", default=os.path.join(
        DEFAULT_REPORTS, "24_novelty_by_field", "novelty_field_counts.tsv"))
    cmp_.add_argument("--clusters", default=os.path.join(
        DEFAULT_REPORTS, "24_novelty_by_field", "protein_level_novel_clusters.tsv"))
    cmp_.add_argument("--panel-c-totals", default=os.path.join(
        DEFAULT_OUT_DIR, "panel_c_novelty_totals.tsv"),
        help="From 40b_novelty_rate_export.py (VM). Falls back to the heavier-censored "
             "sub-split rendering if this file doesn't exist yet.")
    cmp_.add_argument("--c-metric", choices=["auto", "any_field", "cds", "protein"],
                      default="auto",
                      help="auto picks the most specific metric (protein > cds > any_field) "
                           "that has >=50%% of cells at >=20 in the totals export.")
    cmp_.add_argument("--curves", default=os.path.join(
        DEFAULT_REPORTS, "39_saturation_by_ancestry", "curves.tsv"))
    cmp_.add_argument("--oe-table", default=os.path.join(
        DEFAULT_REPORTS, "37_dq_g1g2_signed_ld", "oe_purge_committed.tsv"))
    cmp_.add_argument("--ancestry-scheme", choices=["strict", "pred"], default="strict")
    cmp_.add_argument("--gene", default="HLA-B")
    cmp_.add_argument("--min-carriers", type=int, default=20)
    cmp_.add_argument("--strict-threshold", type=float, default=0.98)
    cmp_.add_argument("--pick", choices=["A", "B"], default="A")
    cmp_.add_argument("--out-dir", default=DEFAULT_OUT_DIR)

    args = ap.parse_args(argv)
    if args.mode == "vm-panels":
        run_vm_panels(args)
    else:
        run_compose(args)


if __name__ == "__main__":
    main()
