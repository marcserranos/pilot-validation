#!/usr/bin/env python3
"""Figure 1 v5 — the introductory figure of the paper, redrawn to Cole's 2026-09-22 panel-by-panel
feedback (sprints/CALL_SUMMARY_2026-09-22.md Sec 6) and Nature's own figure spec (double column,
183 mm, <=~170 mm tall, 5-7 pt text, bold lowercase 8 pt panel letters).

Two-stage pipeline, because panels (a) and (b) are built from per-person admixture/carriage that
must never leave the VM (per-person renormalised admixture IS participant-level information):

  Stage 1 (`--mode vm-panels`, run on the VM): loads Table 1 + cohort membership, reuses
  36_figure1_native.py's own `panel_a_data`/`panel_b_data` (never re-derived -- this project's one
  definition of "unrelated, strict-ancestry, per-allele centroid"), draws panel (a) ALONE as a
  vector figure (no '>=98%' annotation, per Cole's ask), and exports two aggregate,
  disclosure-safe tables: panel_a_group_sizes.tsv (group sizes + fraction clearing the strict
  threshold) and panel_b_ternary_alleles.tsv (per-allele simplex centroid + suppressed carrier
  count, >=20-carrier floor already enforced by 36's own panel_b_data). Only the rendered panel
  (a) image and these two aggregate tables leave the VM.

  Stage 2 (`--mode compose`, run locally): composes the full figure from
    - the VM-rendered panel (a) image (embedded as-is, already disclosure-safe)
    - the pulled panel_b_ternary_alleles.tsv, redrawn locally with 0-100 tick marks on all three
      axes (Cole's ask) at the stricter ancestry filter already baked into the VM export
    - panel (c)+(d) merged: a gene x ancestry novelty-rate heatmap (this IS 33_figure1_v3's "the
      800,000 one" panel, Cole's explicit preferred replacement for the old reference-incompleteness
      bar chart), with panel (d) -- novel protein alleles per gene, labelled 'protein alleles not
      in IPD-IMGT/HLA/HLA' (not 'not in classical genes') -- as a marginal bar aligned to the same
      gene rows, per Cole's "d as a marginal bar on c's heatmap" suggestion. Restricted to the 8
      classical genes shown in the heatmap; the extended, ancestry-split gene list is the existing
      34_novel_recurrence/ supplement, referenced in the legend rather than duplicated here.
    - panel (e): per-ancestry allele-discovery curves from 39_saturation_by_ancestry, re-read at
      render time (another agent was correcting 39's extrapolation fits concurrently -- the curve
      *points* themselves (mean_distinct/lo2_5/hi97_5 per N) are independent of the Clench
      extrapolation fit and unaffected, but this script always reads curves.tsv fresh rather than
      caching it, so a corrected run is picked up automatically), with direct end-of-curve labels
      instead of a legend box.
  Produces two layout variants (A, B) and a `--pick {A,B}` final export.

Usage (VM):
    python3 scripts/hla_popgen/40_figure1_v5.py vm-panels --out-dir ~/s03/results/40

Usage (local, after pulling panel_a.png/.pdf + the two aggregate tables back):
    python3 scripts/hla_popgen/40_figure1_v5.py compose \\
        --panel-a-image reports/hla_popgen/40_figure1_v5/panel_a.png \\
        --panel-a-table reports/hla_popgen/40_figure1_v5/panel_a_group_sizes.tsv \\
        --panel-b-table reports/hla_popgen/40_figure1_v5/panel_b_ternary_alleles.tsv \\
        --pick A
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
CLASSICAL_BARE = [g.replace("HLA-", "") for g in CLASSICAL]
TERNARY = ["AFR", "EUR", "AMR"]

FS_PANEL = 8       # Nature spec: bold lowercase panel letters, 8pt
FS_TITLE = 6.5      # no in-panel titles per design brief -- used only for axis/cbar labels
FS_LABEL = 6.0
FS_TICK = 5.0
FS_LEG = 5.0
FS_ANNOT = 5.0


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
# STAGE 1 -- VM only: panel (a) rendered image + panel (a)/(b) aggregate tables
# ===========================================================================================
def draw_panel_a_alone(pa, strict, n_people, out_stem):
    """Panel (a) admixture barcode, standalone, vector. NO '>=98%' annotation under each block
    (Cole's ask, 2026-09-22 call, panel a row) -- the earlier 36_figure1_native.py version printed
    it directly on the panel; here it is dropped from the figure and kept only in the pulled
    panel_a_group_sizes.tsv table, where it belongs (a number, not a plot annotation)."""
    import matplotlib
    matplotlib.use("Agg")
    with vc.nature_style():
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(30)))
        x = np.arange(len(pa))
        bottom = np.zeros(len(pa))
        for a in ANC:
            v = pa["p_" + a.lower()].to_numpy(dtype=float)
            ax.fill_between(x, bottom, bottom + v, color=ANC_COLORS[a], linewidth=0, label=a)
            bottom += v
        codes = pa["anc"].to_numpy()
        for i in range(1, len(codes)):
            if codes[i] != codes[i - 1]:
                ax.axvline(i, color="white", lw=0.5)
        ax.set_xlim(0, len(pa))
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_ylabel("admixture\nproportion", fontsize=FS_LABEL)
        ax.tick_params(labelsize=FS_TICK)
        start = 0
        for i in range(1, len(codes) + 1):
            if i == len(codes) or codes[i] != codes[start]:
                mid = (start + i) / 2.0
                if (i - start) / float(len(codes)) > 0.02:
                    ax.text(mid, -0.06, codes[start], ha="center", va="top",
                            fontsize=FS_TICK, fontweight="bold", color="#333333",
                            transform=ax.get_xaxis_transform())
                start = i
        ax.legend(ncol=6, frameon=False, fontsize=FS_LEG, loc="upper center",
                  bbox_to_anchor=(0.5, -0.28), handlelength=1.0, columnspacing=1.0)
        vc.save_fig(fig, out_stem)


def run_vm_panels(args):
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
    all_unrelated = set(keep["person_id"].astype(str))
    print("[40] unrelated=%d, strict-ancestry (>=%.2f)=%d" %
          (len(all_unrelated), args.strict_threshold, len(people_keep)), flush=True)

    pa, pa_cols = m36.panel_a_data(cohort, all_unrelated, args.strict_threshold)
    pb = m36.panel_b_data(t1, cohort, people_keep, args.gene, args.min_carriers,
                          args.strict_threshold)
    print("[40] panel a n=%d | panel b %d alleles (gene=%s)" % (len(pa), len(pb), args.gene),
          flush=True)

    draw_panel_a_alone(pa, args.strict_threshold, len(pa),
                       os.path.join(args.out_dir, "panel_a"))

    grp = (pa.groupby("anc").agg(n_people=("anc", "size"),
                                 frac_passing_strict=("passes_strict", "mean")).reset_index())
    grp["n_people_disp"] = grp["n_people"].apply(suppress)
    grp.to_csv(os.path.join(args.out_dir, "panel_a_group_sizes.tsv"), sep="\t", index=False)
    pb.drop(columns=["n_carriers"]).to_csv(
        os.path.join(args.out_dir, "panel_b_ternary_alleles.tsv"), sep="\t", index=False)

    with open(os.path.join(args.out_dir, "vm_panels_summary.json"), "w") as fh:
        json.dump({"n_people_panel_a": len(pa), "n_people_strict_ancestry": len(people_keep),
                   "n_alleles_panel_b": len(pb), "gene_b": args.gene,
                   "strict_threshold": args.strict_threshold,
                   "min_carriers": args.min_carriers, "n_removed_relatedness": n_removed}, fh,
                  indent=2)
    print("[40 vm-panels] done -> %s (pull panel_a.png/.pdf, panel_a_group_sizes.tsv, "
          "panel_b_ternary_alleles.tsv, vm_panels_summary.json back)" % args.out_dir)


# ===========================================================================================
# STAGE 2 -- local: compose the full figure
# ===========================================================================================
def tern_xy(p):
    l, r, t = p
    return r + 0.5 * t, (np.sqrt(3) / 2.0) * t


def draw_ternary_grid(ax):
    """0-100 tick marks + gridlines on all three ternary axes (Cole's ask, panel b row)."""
    V = {"AFR": (1, 0, 0), "EUR": (0, 1, 0), "AMR": (0, 0, 1)}
    xy = {k: tern_xy(v) for k, v in V.items()}
    tri = np.array([xy["AFR"], xy["EUR"], xy["AMR"], xy["AFR"]])
    ax.plot(tri[:, 0], tri[:, 1], color="#333333", lw=0.6, zorder=2)
    ticks = [0, 20, 40, 60, 80, 100]
    # AFR axis (l): gridlines of constant l=k, parallel to EUR-AMR edge; ticks along the AFR-AMR
    # edge (r=0 side), reading 0 at AMR vertex up to 100 at AFR vertex.
    for k in ticks:
        kk = k / 100.0
        if 0 < k < 100:
            p0 = tern_xy((kk, 1 - kk, 0)); p1 = tern_xy((kk, 0, 1 - kk))
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#DDDDDD", lw=0.3, zorder=1)
        # tick label on the AFR-EUR edge (t=0 side): point where l=k, r=1-k, t=0
        tp = tern_xy((kk, 1 - kk, 0))
        ax.text(tp[0], tp[1] - 0.035, str(k), fontsize=FS_ANNOT, ha="center", va="top",
                color="#666666")
    # EUR axis (r): ticks along AMR-AFR edge is already the l axis; put EUR ticks along the
    # EUR-AMR edge (l=0 side): point where r=k, t=1-k, l=0
    for k in ticks:
        kk = k / 100.0
        tp = tern_xy((0, kk, 1 - kk))
        ax.text(tp[0] - 0.02, tp[1], str(k), fontsize=FS_ANNOT, ha="right", va="center",
                color="#666666")
    # AMR axis (t): ticks along the AFR-EUR edge is shared with l; put AMR ticks along the
    # EUR-AFR... use the remaining free edge: AFR-AMR edge (r=0 side): point where t=k, l=1-k, r=0
    for k in ticks:
        kk = k / 100.0
        tp = tern_xy((1 - kk, 0, kk))
        ax.text(tp[0] + 0.02, tp[1], str(k), fontsize=FS_ANNOT, ha="left", va="center",
                color="#666666")
    return xy


def draw_panel_b(ax, pb, gene_b, min_carriers, strict, n_alleles_note=True):
    xy = draw_ternary_grid(ax)
    if not pb.empty:
        n = pb["n_carriers_disp"].apply(lambda s: 19 if str(s).startswith("<") else float(s))
        n = n.to_numpy(dtype=float)
        n = np.maximum(n, 1.0)
        lo, hi = np.log10(n.min()), np.log10(n.max())
        rng = max(1e-9, hi - lo)
        alpha = 0.25 + 0.65 * (np.log10(n) - lo) / rng
        size = 6 + 34 * (np.log10(n) - lo) / rng
        for (_, r), al, sz in zip(pb.iterrows(), alpha, size):
            xyp = tern_xy((r[TERNARY[0]], r[TERNARY[1]], r[TERNARY[2]]))
            novel = r["frac_novel_calls"] > 0.5
            dom = max(TERNARY, key=lambda a: r[a])
            ax.scatter([xyp[0]], [xyp[1]], s=sz, marker="^" if novel else "o",
                      facecolor=ANC_COLORS[dom], edgecolor="none", alpha=float(al), zorder=3)
    for lab, pos, ha, va, off in [("AFR", xy["AFR"], "right", "top", (-3, -8)),
                                  ("EUR", xy["EUR"], "left", "top", (3, -8)),
                                  ("AMR", xy["AMR"], "center", "bottom", (0, 9))]:
        ax.annotate(lab, pos, fontsize=FS_LABEL, fontweight="bold", ha=ha, va=va,
                    xytext=off, textcoords="offset points")
    ax.set_xlim(-0.16, 1.16)
    ax.set_ylim(-0.14, np.sqrt(3) / 2 + 0.13)
    ax.axis("off")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", color="#666666", ms=3,
                              label="catalogued"),
                       Line2D([], [], marker="^", ls="none", color="#666666", ms=3,
                              label="first observed here")],
              frameon=False, fontsize=FS_LEG, loc="upper right", handletextpad=0.3,
              bbox_to_anchor=(1.05, 1.08))


def load_panel_cd(field_counts_path, clusters_path, ancestry_scheme="strict"):
    """(c) novelty-rate heatmap table (gene x ancestry, lower-bound %) + (d) marginal bar table
    (gene -> once/mid/high novel-protein-allele counts, classical genes only, pooled ancestry)."""
    m33 = _load("33_figure1_v3_compose.py", "fig1v3_compose")
    counts = pd.read_csv(field_counts_path, sep="\t")
    rates = m33.rate_by_gene_ancestry(counts, ancestry_scheme, True, CLASSICAL)
    pivot = (rates[rates["field_class"] != "artifact_control"]
             .pivot_table(index="gene", columns="ancestry", values="pct_lower", aggfunc="sum")
             .reindex(index=CLASSICAL, columns=ANC))
    n_cens = (rates.pivot_table(index="gene", columns="ancestry", values="n_censored_cells",
                                aggfunc="sum").reindex(index=CLASSICAL, columns=ANC))
    n_called = (rates.pivot_table(index="gene", columns="ancestry", values="n_called_upper",
                                  aggfunc="max").reindex(index=CLASSICAL, columns=ANC))

    m34 = _load("34_novel_protein_recurrence.py", "novel_recurrence")
    cl = pd.read_csv(clusters_path, sep="\t", dtype=str)
    prot = cl[cl["cluster_type"] == "novel_protein"]
    prot = m34.classify(prot)
    dtab = m34.gene_table(prot)
    dtab = dtab.reindex(CLASSICAL).fillna(0)
    return pivot, n_cens, n_called, dtab


def draw_panel_cd(ax_hm, ax_bar, pivot, n_cens, n_called, dtab, cax=None):
    genes = CLASSICAL
    ancs = ANC
    data = pivot.reindex(index=genes, columns=ancs).to_numpy(dtype=float)
    im = ax_hm.imshow(data, cmap="YlOrRd", aspect="auto", vmin=0,
                      vmax=np.nanmax(data) if np.nanmax(data) > 0 else 1)
    for i, g in enumerate(genes):
        for j, a in enumerate(ancs):
            nc = n_cens.loc[g, a] if (g in n_cens.index and a in n_cens.columns) else 0
            called_hi = n_called.loc[g, a] if (g in n_called.index and a in n_called.columns) else 0
            if pd.isna(data[i, j]) or (nc and nc > 0 and called_hi < SUPPRESS_BELOW):
                vc.hatch_suppressed(ax_hm, j - 0.5, i - 0.5, 1, 1)
            else:
                ax_hm.text(j, i, "%.0f" % data[i, j], ha="center", va="center",
                          fontsize=FS_ANNOT, color="#222222")
    ax_hm.set_xticks(range(len(ancs)))
    ax_hm.set_xticklabels(ancs, fontsize=FS_TICK)
    ax_hm.set_yticks(range(len(genes)))
    ax_hm.set_yticklabels([g.replace("HLA-", "") for g in genes], fontsize=FS_TICK)
    ax_hm.set_ylabel("% called haplotypes carrying sequence\nabsent from IPD-IMGT/HLA",
                     fontsize=FS_LABEL)
    for sp in ax_hm.spines.values():
        sp.set_visible(False)
    if cax is not None:
        cb = ax_hm.figure.colorbar(im, cax=cax, orientation="vertical")
        cb.ax.tick_params(labelsize=FS_TICK)
        cb.set_label("% novel haplotypes", fontsize=FS_LABEL)

    # marginal bar, same gene rows, horizontal, stacked by recurrence class
    cols = [("seen once", "#C9CFD6"), ("2–19 unrelated people", "#5B8FBF"),
            ("≥20 unrelated people", "#B4472E")]
    y = np.arange(len(genes))
    left = np.zeros(len(genes))
    for key, col in cols:
        v = dtab.reindex(genes)[key].to_numpy(dtype=float) if key in dtab.columns else np.zeros(len(genes))
        ax_bar.barh(y, v, left=left, height=0.72, color=col, label=key, edgecolor="white",
                   linewidth=0.3, zorder=3)
        left += v
    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels([])
    ax_bar.set_ylim(-0.5, len(genes) - 0.5)
    ax_bar.invert_yaxis()
    ax_hm.invert_yaxis()
    ax_bar.set_xlabel("protein alleles not in\nIPD-IMGT/HLA", fontsize=FS_LABEL)
    ax_bar.tick_params(axis="x", labelsize=FS_TICK)
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.legend(frameon=False, fontsize=FS_LEG, loc="upper center",
                 bbox_to_anchor=(0.5, -0.24), ncol=1, handlelength=1.0, labelspacing=0.25)


def draw_panel_e(ax, curves_path, scheme="pred"):
    df = pd.read_csv(curves_path, sep="\t")
    d = df[(df["scheme"] == scheme) & (df["gene_group"] == "classical_pooled")
          & (df["category"] == "all")].copy()
    ends = {}
    xmax = 0.0
    for a in ANC:
        sub = d[d["ancestry"] == a].sort_values("n")
        if sub.empty:
            continue
        x = sub["n"].to_numpy(dtype=float)
        y = sub["mean_distinct"].to_numpy(dtype=float)
        lo = sub["lo2_5"].to_numpy(dtype=float)
        hi = sub["hi97_5"].to_numpy(dtype=float)
        ax.plot(x, y, color=ANC_COLORS[a], lw=1.0, zorder=3)
        ax.fill_between(x, lo, hi, color=ANC_COLORS[a], alpha=0.15, linewidth=0, zorder=2)
        ends[a] = (x[-1], y[-1])
        xmax = max(xmax, x[-1])
    # direct end labels, spaced apart vertically so close-together curve ends don't collide --
    # sort by y, then greedily push labels apart by a minimum fraction of the y-range.
    if ends:
        yr = ax.get_ylim()
        ymin_data = min(v[1] for v in ends.values())
        ymax_data = max(v[1] for v in ends.values())
        min_gap = max(1.0, (ymax_data - ymin_data) * 0.09) if ymax_data > ymin_data else 1.0
        order = sorted(ends, key=lambda a: ends[a][1])
        placed = []
        for a in order:
            y = ends[a][1]
            if placed and y - placed[-1] < min_gap:
                y = placed[-1] + min_gap
            placed.append(y)
        for a, y_lab in zip(order, placed):
            x0, y0 = ends[a]
            ax.annotate(a, (x0, y_lab), xytext=(4, 0), textcoords="offset points",
                       fontsize=FS_ANNOT, color=ANC_COLORS[a], fontweight="bold", va="center")
        ax.set_xlim(0, xmax * 1.16)
    ax.set_xlabel("people sampled (both haplotypes)", fontsize=FS_LABEL)
    ax.set_ylabel("distinct HLA protein alleles\n(8 classical genes, pooled)", fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.spines[["top", "right"]].set_visible(False)


def compose_layout(layout, panel_a_img, pb, pivot, n_cens, n_called, dtab, curves_path,
                   n_people_a, gene_b, min_carriers, strict):
    import matplotlib
    matplotlib.use("Agg")
    with vc.nature_style():
        import matplotlib.pyplot as plt
        import matplotlib.image as mpimg
        from matplotlib.gridspec import GridSpec

        W, H = vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(168)
        fig = plt.figure(figsize=(W, H))

        if layout == "A":
            # full-width strip (a) on top, then b | c+d, then e | spare
            gs = GridSpec(3, 2, height_ratios=[0.62, 1.35, 1.15], width_ratios=[1.0, 1.35],
                         hspace=0.62, wspace=0.32, left=0.075, right=0.965, top=0.975,
                         bottom=0.10, figure=fig)
            ax_a = fig.add_subplot(gs[0, :])
            ax_b = fig.add_subplot(gs[1, 0])
            gs_cd = gs[1, 1].subgridspec(1, 3, width_ratios=[1.0, 0.5, 0.06], wspace=0.12)
            ax_c = fig.add_subplot(gs_cd[0, 0])
            ax_d = fig.add_subplot(gs_cd[0, 1])
            cax = fig.add_subplot(gs_cd[0, 2])
            ax_e = fig.add_subplot(gs[2, 0])
            ax_spare = fig.add_subplot(gs[2, 1])
            ax_spare.axis("off")
            ax_spare.text(0.02, 0.9, "(reserved: KIR panel, pending WS3 rerun)",
                          fontsize=FS_ANNOT, color="#999999", transform=ax_spare.transAxes)
        else:
            # layout B: (a) top strip; (b) and (e) share the left column stacked; (c+d) full
            # height on the right -- keeps the heatmap tall enough to read all 8 gene rows without
            # a marginal bar dominating the page, and puts the two point-cloud/line panels (b, e)
            # together since they share the same ancestry-colour reading convention.
            gs = GridSpec(3, 2, height_ratios=[0.62, 1.0, 1.0], width_ratios=[1.0, 1.15],
                         hspace=0.55, wspace=0.30, left=0.075, right=0.965, top=0.975,
                         bottom=0.10, figure=fig)
            ax_a = fig.add_subplot(gs[0, :])
            ax_b = fig.add_subplot(gs[1, 0])
            ax_e = fig.add_subplot(gs[2, 0])
            gs_cd = gs[1:, 1].subgridspec(1, 3, width_ratios=[1.0, 0.45, 0.055], wspace=0.14)
            ax_c = fig.add_subplot(gs_cd[0, 0])
            ax_d = fig.add_subplot(gs_cd[0, 1])
            cax = fig.add_subplot(gs_cd[0, 2])

        if panel_a_img and os.path.exists(panel_a_img):
            img = mpimg.imread(panel_a_img)
            ax_a.imshow(img, aspect="auto")
        else:
            ax_a.set_xlim(0, 1); ax_a.set_ylim(0, 1)
            ax_a.add_patch(plt.Rectangle((0, 0), 1, 1, fill=False, edgecolor="#BBBBBB", lw=0.6,
                                         linestyle="--"))
            ax_a.text(0.5, 0.5, "panel a pending VM re-render (see README)", ha="center",
                      va="center", fontsize=FS_LABEL, color="#999999", transform=ax_a.transAxes)
        ax_a.axis("off")
        vc.panel_letter(ax_a, "a", dx=-0.01, dy=1.02)

        draw_panel_b(ax_b, pb, gene_b, min_carriers, strict)
        vc.panel_letter(ax_b, "b", dx=-0.06)

        draw_panel_cd(ax_c, ax_d, pivot, n_cens, n_called, dtab, cax=cax)
        vc.panel_letter(ax_c, "c", dx=-0.30)
        vc.panel_letter(ax_d, "d", dx=-0.10)

        draw_panel_e(ax_e, curves_path)
        vc.panel_letter(ax_e, "e", dx=-0.14)

        note = ("a  all unrelated participants (n=%s), predicted ancestry, ordered least->most "
                "admixed within block; per-group %% clearing the strict ancestry threshold "
                "(>=%.2f) is in panel_a_group_sizes.tsv, not annotated on the panel.\n"
                "b  HLA-%s alleles, strict ancestry (probability >=%.2f), >=%d carriers; point "
                "area/opacity ~ carrier count (log).   c  exact lower-bound rate; hatched cell = "
                "denominator <%d haplotypes.   d  pooled across ancestry (ancestry split: "
                "reports/hla_popgen/34_novel_recurrence/, all genes).   e  25 permutations/point, "
                "mean +/- 95%% band, `pred` ancestry scheme; see reports/hla_popgen/"
                "39_saturation_by_ancestry/ for the strict95 sensitivity curve."
                % (n_people_a, strict, gene_b.replace("HLA-", ""), strict, min_carriers,
                   SUPPRESS_BELOW))
        fig.text(0.01, 0.006, note, fontsize=FS_ANNOT, color="#666666", ha="left", va="bottom",
                 linespacing=1.4)
        return fig


def run_compose(args):
    os.makedirs(args.out_dir, exist_ok=True)
    pb = pd.read_csv(args.panel_b_table, sep="\t")
    pa_groups = pd.read_csv(args.panel_a_table, sep="\t")
    n_people_a = int(pd.to_numeric(pa_groups["n_people"], errors="coerce").sum())

    pivot, n_cens, n_called, dtab = load_panel_cd(args.field_counts, args.clusters,
                                                  args.ancestry_scheme)

    figs = {}
    for layout in ["A", "B"]:
        fig = compose_layout(layout, args.panel_a_image, pb, pivot, n_cens, n_called, dtab,
                             args.curves, n_people_a, args.gene, args.min_carriers,
                             args.strict_threshold)
        stem = os.path.join(args.out_dir, "layout_%s" % layout)
        vc.save_fig(fig, stem)
        figs[layout] = stem + ".png"

    chosen = args.pick
    import shutil
    shutil.copy(figs[chosen] , os.path.join(args.out_dir, "figure1_v5.png"))
    shutil.copy(figs[chosen].replace(".png", ".pdf"),
               os.path.join(args.out_dir, "figure1_v5.pdf"))

    with open(os.path.join(args.out_dir, "compose_summary.json"), "w") as fh:
        json.dump({"layouts_built": ["A", "B"], "chosen": chosen,
                   "n_people_panel_a": n_people_a, "gene_b": args.gene,
                   "min_carriers_b": args.min_carriers,
                   "strict_threshold": args.strict_threshold,
                   "ancestry_scheme_cd": args.ancestry_scheme}, fh, indent=2)
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
    cmp_.add_argument("--panel-a-image", default=os.path.join(DEFAULT_OUT_DIR, "panel_a.png"))
    cmp_.add_argument("--panel-a-table",
                      default=os.path.join(DEFAULT_OUT_DIR, "panel_a_group_sizes.tsv"))
    cmp_.add_argument("--panel-b-table",
                      default=os.path.join(DEFAULT_OUT_DIR, "panel_b_ternary_alleles.tsv"))
    cmp_.add_argument("--field-counts", default=os.path.join(
        DEFAULT_REPORTS, "24_novelty_by_field", "novelty_field_counts.tsv"))
    cmp_.add_argument("--clusters", default=os.path.join(
        DEFAULT_REPORTS, "24_novelty_by_field", "protein_level_novel_clusters.tsv"))
    cmp_.add_argument("--curves", default=os.path.join(
        DEFAULT_REPORTS, "39_saturation_by_ancestry", "curves.tsv"))
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
