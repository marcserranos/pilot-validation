#!/usr/bin/env python3
"""Figure 1, drawn natively from data as a single figure.

Replaces script 35, which pasted four rendered PNGs into one sheet. That collage had four
different type scales, three internal working titles ("alpha ~ carrier count (no
minimum-carrier-count floor)"), duplicated panel letters, and an aspect ratio no journal would
take. This script draws all four panels in one matplotlib figure with one type scale, one legend
convention and one palette.

It must run **on the VM**, because panel (a) needs per-person admixture proportions and panel (b)
needs per-person allele carriage. Neither may leave the controlled environment. What leaves is the
rendered figure plus **aggregate, disclosure-safe tables** for panels b/c/d, so those three become
restyleable offline later — closing the gap WS4 flagged, where scripts 06 and 10 committed figures
but never their underlying tables.

Panels (Cole's spec, 2026-09-17 ask A1):
  a  cohort ancestry (admixture barcode)
  b  HLA-B per-allele ancestry centroids on the AFR/EUR/AMR simplex
  c  novelty by gene, ancestry and nomenclature field, with the artifact rate as negative control
  d  allele-frequency spectrum, novel vs already catalogued

Two things this version fixes beyond the layout:

* **Ask A9, the stricter admixture threshold.** Panels (a) and (b) now take `--strict-threshold`
  (default 0.98, which is what Cole and David asked for). Script 10 had no such parameter.
* **A disclosure hole in the old panel (b).** Script 10's ternary plotted a point per allele with
  *no minimum carrier count* — a point representing a single carrier is that person's own
  renormalised admixture, which is participant-level information. Here an allele must reach
  `--min-carriers` (default 20) to be drawn, which is both compliant and visually cleaner.

Exact counts are used for every rate, since this runs where the data is uncensored; no count is
ever written to an output file below the suppression threshold, and no rate is drawn from a
denominator under it.

Usage (VM):

    python3 scripts/hla_popgen/36_figure1_native.py --out-dir ~/results/36_figure1
"""
import argparse
import importlib.util
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)


def _load(fn, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_THIS_DIR, fn))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_TABLE1 = os.path.join(DEFAULT_OUTROOT, "hla_calls_rich.tsv")
DEFAULT_COHORT = os.path.join(DEFAULT_OUTROOT, "cohort_membership.tsv")
DEFAULT_CLUSTERS = os.path.expanduser(
    "~/results/24_novelty_by_field/protein_level_novel_clusters.tsv")
DEFAULT_RELATEDNESS = os.path.expanduser(
    "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/"
    "relatedness/samples_relatedness.tsv")
DEFAULT_OUT_DIR = os.path.expanduser("~/results/36_figure1")

SUPPRESS_BELOW = 20
KIN_MIN = 0.0442
ANC = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
ANC_COLORS = {"AFR": "#D55E00", "AMR": "#E69F00", "EAS": "#009E73",
              "EUR": "#0072B2", "MID": "#CC79A7", "SAS": "#56B4E9"}
CLASSICAL = ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"]
TERNARY = ["AFR", "EUR", "AMR"]

FIELD_ORDER = ["f4_noncoding", "f3_synonymous", "f2_protein"]
FIELD_LABEL = {"f4_noncoding": "new non-coding only",
               "f3_synonymous": "new synonymous CDS",
               "f2_protein": "new protein"}
# Hatching does not resolve at panel size in print; three alpha steps plus an outline on the
# rarest class read cleanly instead.
FIELD_ALPHA = {"f4_noncoding": 1.0, "f3_synonymous": 0.55, "f2_protein": 0.22}
FIELD_EDGE = {"f4_noncoding": "none", "f3_synonymous": "none", "f2_protein": "#333333"}

# One type scale for the whole figure.
FS_PANEL = 11     # panel letter
FS_TITLE = 9.0    # panel title
FS_LABEL = 8.0    # axis labels
FS_TICK = 7.0
FS_LEG = 7.0
FS_ANNOT = 6.2


def suppress(n):
    n = int(n)
    return "0" if n == 0 else ("<%d" % SUPPRESS_BELOW if n < SUPPRESS_BELOW else str(n))


def two_field(name):
    if not isinstance(name, str) or "*" not in name:
        return None
    g, rest = name.split("*", 1)
    f = rest.split(":")
    if len(f) < 2 or f[0].strip().lower() == "new" or f[1].strip().lower() == "new":
        return None
    return "%s*%s:%s" % (g.replace("HLA-", ""), f[0], f[1])


# ---------------------------------------------------------------------------
# panel data
# ---------------------------------------------------------------------------
def panel_a_data(cohort, people_all, strict):
    """Per-person admixture, ordered for the barcode.

    **Uses every unrelated person on their PREDICTED ancestry, not the strict subset.** The first
    version filtered this panel at the same 0.98 threshold as panels b and c, which made every bar
    a solid block of one colour — a plot of admixture with the admixture filtered out. The
    threshold belongs to the analyses that need clean ancestry labels, not to the panel whose job
    is to show how admixed the cohort actually is. The fraction of each group that survives the
    strict filter is annotated instead, which is the more useful statement.

    Sorted by predicted ancestry, then by that ancestry's own proportion descending, so each block
    runs from the least admixed person to the most.
    """
    d = cohort[cohort["person_id"].astype(str).isin(people_all)].copy()
    cols = ["p_" + a.lower() for a in ANC]
    for c in cols:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=cols)
    d["anc"] = d["ancestry_pred"].astype(str).str.upper()
    d = d[d["anc"].isin(ANC)]
    d["dom"] = [row["p_" + row["anc"].lower()] for _, row in d.iterrows()]
    d["passes_strict"] = d["dom"] >= strict
    d = d.sort_values(["anc", "dom"], ascending=[True, False])
    return d, cols


def panel_b_data(t1, cohort, people_keep, gene, min_carriers, strict):
    """Per-allele ancestry centroid on the AFR/EUR/AMR simplex.

    For each two-field allele, take everyone carrying it, renormalise each carrier's admixture
    proportions to the three displayed ancestries, and average. Renormalising **per person before
    averaging** (script 10's `per_person` mode) is the right order: it means a person who is 100%
    East Asian contributes nothing rather than a misleading value.
    """
    pcols = ["p_" + a.lower() for a in TERNARY]
    co = cohort.set_index(cohort["person_id"].astype(str))
    for c in pcols:
        co[c] = pd.to_numeric(co[c], errors="coerce")

    sub = t1[t1["gene"].astype(str) == gene].copy()
    sub = sub[sub["person_id"].astype(str).isin(people_keep)]
    sub["allele"] = sub["consensus"].map(two_field)
    sub = sub.dropna(subset=["allele"])
    sub["is_novel"] = sub["consensus"].astype(str).str.contains("new")

    rows = []
    for allele, grp in sub.groupby("allele"):
        pids = sorted(set(grp["person_id"].astype(str)))
        if len(pids) < min_carriers:
            continue
        P = co.reindex(pids)[pcols].to_numpy(dtype=float)
        s = P.sum(axis=1)
        ok = s > 0
        if ok.sum() < min_carriers:
            continue
        cent = (P[ok] / s[ok, None]).mean(axis=0)
        rows.append({"allele": allele, "n_carriers": len(pids),
                     "n_carriers_disp": suppress(len(pids)),
                     "frac_novel_calls": round(float(grp["is_novel"].mean()), 4),
                     TERNARY[0]: round(float(cent[0]), 5),
                     TERNARY[1]: round(float(cent[1]), 5),
                     TERNARY[2]: round(float(cent[2]), 5)})
    return pd.DataFrame(rows)


def panel_c_data(t1, anc_of, m24):
    """Novelty rate per gene x ancestry x nomenclature field, plus the artifact rate.

    Exact counts, because this runs where the table is uncensored. A gene x ancestry cell is only
    emitted when its denominator clears the suppression threshold, so no rate is ever drawn from a
    handful of haplotypes.
    """
    d = t1[t1["gene"].astype(str).isin(CLASSICAL)].copy()
    d["ancestry"] = d["person_id"].astype(str).map(anc_of)
    d = d[d["ancestry"].isin(ANC)]
    d = m24.add_call_labels(d)

    rows = []
    for (g, a), sub in d.groupby(["gene", "ancestry"]):
        called = sub[sub["field_class"] != "uncalled"]
        n = len(called)
        if n < SUPPRESS_BELOW:
            continue
        clean = called[called["artifact_label"] == "clean"]
        rec = {"gene": g, "ancestry": a, "n_called_disp": suppress(n)}
        for fc in FIELD_ORDER:
            rec[fc] = round(100.0 * int((clean["field_class"] == fc).sum()) / n, 4)
        rec["artifact"] = round(
            100.0 * int((called["artifact_label"] != "clean").sum()) / n, 4)
        rows.append(rec)
    return pd.DataFrame(rows)


def panel_d_data(t1, clusters, anc_of):
    """Novel protein alleles per gene, split by how many unrelated people carry them.

    **This replaced an allele-frequency spectrum, and the reason matters.** Script 24 applies the
    All of Us small-cell rule *when it writes its output*, so every carrier count from 1 to 19 is
    stored as the string `<20` even in the copy that lives on the VM. The uncensored counts are
    not recoverable from that file at all — a per-count spectrum or a 1,2,3…19 histogram would
    need script 24's clustering re-run with an uncensored internal dump. Rather than draw a novel
    series that is silently all zeros (which is exactly what the first version of this figure
    did), the panel shows the three-class split the data can actually support:

        seen once  |  2-19 unrelated people  |  >= 20 unrelated people

    which is the same information at the resolution the disclosure rule permits, and is per gene,
    so it carries the finding: the novelty is not in the classical genes.
    """
    cl = clusters[clusters["cluster_type"] == "novel_protein"].copy()

    def klass(r):
        rec = str(r.get("clean_recurrent_unrelated", "")).strip().lower() in {"true", "1", "yes"}
        raw = str(r.get("n_persons_unrelated_clean", "")).strip()
        if not rec:
            return "once"
        return "high" if (not raw.startswith("<") and raw.replace(".", "").isdigit()
                          and float(raw) >= SUPPRESS_BELOW) else "mid"

    cl["klass"] = cl.apply(klass, axis=1)
    t = (cl.groupby(["gene", "klass"]).size().unstack(fill_value=0)
         .reindex(columns=["once", "mid", "high"], fill_value=0))
    t["total"] = t.sum(axis=1)
    t = t.sort_values("total", ascending=False).reset_index()
    return t


def recurrence_histogram(clusters):
    """The per-count histogram that cannot be built off-VM: novel proteins seen in exactly
    1, 2, 3 ... 19, then >= 20 unrelated people. Only the binned counts leave; they are counts of
    ALLELES, not of participants, so they are not subject to the small-cell rule."""
    cl = clusters[clusters["cluster_type"] == "novel_protein"].copy()
    n = pd.to_numeric(cl.get("n_persons_unrelated_clean"), errors="coerce").fillna(0).astype(int)
    rows = []
    for k in range(1, SUPPRESS_BELOW):
        rows.append({"n_unrelated_carriers": str(k), "n_novel_proteins": int((n == k).sum())})
    rows.append({"n_unrelated_carriers": "≥%d" % SUPPRESS_BELOW,
                 "n_novel_proteins": int((n >= SUPPRESS_BELOW).sum())})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
def tern_xy(p):
    """Barycentric (left, right, top) -> cartesian, for an equilateral triangle of side 1."""
    l, r, t = p
    return r + 0.5 * t, (np.sqrt(3) / 2.0) * t


def draw(pa, pa_cols, pb, pc, pd_, out_png, out_pdf, gene_b, strict, min_carriers, n_people):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": 0.7,
                         "xtick.major.width": 0.7, "ytick.major.width": 0.7})

    fig = plt.figure(figsize=(13.2, 9.1))
    gs = GridSpec(2, 2, width_ratios=[1.42, 1.0], height_ratios=[1.0, 1.18],
                  hspace=0.40, wspace=0.20, left=0.06, right=0.985, top=0.95, bottom=0.11,
                  figure=fig)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    def panel_letter(ax, letter, dx=-0.085):
        ax.text(dx, 1.06, letter, transform=ax.transAxes, fontsize=FS_PANEL,
                fontweight="bold", va="bottom", ha="left")

    # ---- a: admixture barcode -----------------------------------------------------------
    x = np.arange(len(pa))
    bottom = np.zeros(len(pa))
    for a in ANC:
        v = pa["p_" + a.lower()].to_numpy(dtype=float)
        ax_a.fill_between(x, bottom, bottom + v, color=ANC_COLORS[a], linewidth=0, label=a)
        bottom += v
    # boundaries between predicted-ancestry blocks
    codes = pa["anc"].to_numpy()
    for i in range(1, len(codes)):
        if codes[i] != codes[i - 1]:
            ax_a.axvline(i, color="white", lw=0.8)
    ax_a.set_xlim(0, len(pa))
    ax_a.set_ylim(0, 1)
    ax_a.set_xticks([])
    ax_a.set_yticks([0, 0.5, 1.0])
    ax_a.tick_params(labelsize=FS_TICK)
    ax_a.set_ylabel("admixture proportion", fontsize=FS_LABEL)
    ax_a.set_title("Genetic ancestry of the long-read cohort (n = %s unrelated)"
                   % f"{n_people:,}", fontsize=FS_TITLE, loc="left", pad=8)
    # Block labels go BELOW the axis. Above it they collided with the title, which is exactly the
    # kind of overlap that makes a figure look unfinished.
    start = 0
    for i in range(1, len(codes) + 1):
        if i == len(codes) or codes[i] != codes[start]:
            mid = (start + i) / 2.0
            frac = float(pa["passes_strict"].to_numpy()[start:i].mean())
            ax_a.text(mid, -0.03, "%s" % codes[start], ha="center", va="top",
                      fontsize=FS_TICK, fontweight="bold", color="#333333")
            # A narrow block (MID is ~1% of the cohort) cannot hold this label without
            # colliding with its neighbour, so it is left to the exported table instead.
            if (i - start) / float(len(codes)) > 0.045:
                ax_a.text(mid, -0.115, "%.0f%% ≥%.2f" % (100 * frac, strict), ha="center",
                          va="top", fontsize=FS_ANNOT, color="#666666")
            start = i
    ax_a.legend(ncol=6, frameon=False, fontsize=FS_LEG, loc="lower center",
                bbox_to_anchor=(0.5, -0.27), handlelength=1.1, columnspacing=1.2)
    panel_letter(ax_a, "a")

    # ---- b: ternary ---------------------------------------------------------------------
    V = {k: tern_xy(v) for k, v in
         {"AFR": (1, 0, 0), "EUR": (0, 1, 0), "AMR": (0, 0, 1)}.items()}
    tri = np.array([V["AFR"], V["EUR"], V["AMR"], V["AFR"]])
    ax_b.plot(tri[:, 0], tri[:, 1], color="#333333", lw=0.9, zorder=2)
    if not pb.empty:
        n = pb["n_carriers"].to_numpy(dtype=float)
        alpha = 0.25 + 0.65 * (np.log10(n) - np.log10(n.min())) / max(
            1e-9, (np.log10(n.max()) - np.log10(n.min())))
        size = 8 + 42 * (np.log10(n) - np.log10(n.min())) / max(
            1e-9, (np.log10(n.max()) - np.log10(n.min())))
        for (_, r), al, sz in zip(pb.iterrows(), alpha, size):
            xy = tern_xy((r[TERNARY[0]], r[TERNARY[1]], r[TERNARY[2]]))
            novel = r["frac_novel_calls"] > 0.5
            # colour by which vertex dominates, so the eye reads structure not category
            dom = max(TERNARY, key=lambda a: r[a])
            ax_b.scatter([xy[0]], [xy[1]], s=sz, marker="^" if novel else "o",
                         facecolor=ANC_COLORS[dom], edgecolor="none", alpha=float(al), zorder=3)
        if not (pb["frac_novel_calls"] > 0.5).any():
            ax_b.text(0.5, 0.015, "no novel %s allele reaches %d carriers"
                      % (gene_b.replace("HLA-", ""), min_carriers),
                      transform=ax_b.transAxes, ha="center", va="bottom", fontsize=FS_ANNOT,
                      color="#B4472E", style="italic")
        # Only the six most common, with offsets alternating around the point, because eight
        # labels on a crowded simplex overlapped each other in the previous version.
        top = pb.nlargest(6, "n_carriers").copy()
        top["_xy"] = [tern_xy((r[TERNARY[0]], r[TERNARY[1]], r[TERNARY[2]]))
                      for _, r in top.iterrows()]
        top = top.sort_values("_xy", key=lambda c: [p[0] for p in c])
        offsets = [(5, 5), (5, -8), (-5, 6), (-5, -9), (6, 0), (-6, 0)]
        for (_, r), off in zip(top.iterrows(), offsets):
            ha = "left" if off[0] > 0 else "right"
            ax_b.annotate(r["allele"].split("*")[1], r["_xy"], fontsize=FS_ANNOT, ha=ha,
                          xytext=off, textcoords="offset points", color="#222222", zorder=5)
    for lab, pos, ha, va, off in [("AFR", V["AFR"], "right", "top", (-4, -4)),
                                  ("EUR", V["EUR"], "left", "top", (4, -4)),
                                  ("AMR", V["AMR"], "center", "bottom", (0, 5))]:
        ax_b.annotate(lab, pos, fontsize=FS_LABEL, fontweight="bold", ha=ha, va=va,
                      xytext=off, textcoords="offset points")
    ax_b.set_xlim(-0.09, 1.09)
    ax_b.set_ylim(-0.09, np.sqrt(3) / 2 + 0.11)
    ax_b.axis("off")
    ax_b.set_title("%s alleles by ancestry composition of their carriers"
                   % gene_b.replace("HLA-", "HLA-"), fontsize=FS_TITLE, loc="left", pad=6)
    ax_b.legend(handles=[Line2D([], [], marker="o", ls="none", color="#666666", ms=4,
                                label="catalogued in IPD-IMGT/HLA"),
                         Line2D([], [], marker="^", ls="none", color="#666666", ms=4,
                                label="first observed here")],
                frameon=False, fontsize=FS_LEG, loc="upper right", handletextpad=0.4)
    panel_letter(ax_b, "b", dx=-0.02)

    # ---- c: novelty by gene x ancestry --------------------------------------------------
    genes = [g for g in CLASSICAL if g in set(pc["gene"])]
    ancs = [a for a in ANC if a in set(pc["ancestry"])]
    w = 0.80 / max(1, len(ancs))
    for ai, a in enumerate(ancs):
        for gi, g in enumerate(genes):
            r = pc[(pc["gene"] == g) & (pc["ancestry"] == a)]
            if r.empty:
                continue
            r = r.iloc[0]
            xx = gi - 0.40 + w * (ai + 0.5)
            bot = 0.0
            for fc in FIELD_ORDER:
                ax_c.bar(xx, r[fc], bottom=bot, width=w * 0.88, color=ANC_COLORS[a],
                         alpha=FIELD_ALPHA[fc], edgecolor=FIELD_EDGE[fc],
                         linewidth=0.35 if FIELD_EDGE[fc] != "none" else 0, zorder=3)
                bot += r[fc]
            ax_c.plot([xx - w * 0.44, xx + w * 0.44], [r["artifact"], r["artifact"]],
                      color="#111111", lw=1.0, zorder=5, solid_capstyle="butt")
    ax_c.set_xticks(range(len(genes)))
    ax_c.set_xticklabels([g.replace("HLA-", "") for g in genes], fontsize=FS_TICK)
    ax_c.tick_params(axis="y", labelsize=FS_TICK)
    ax_c.set_ylabel("% of called haplotypes carrying\nsequence absent from IPD-IMGT/HLA",
                    fontsize=FS_LABEL)
    ax_c.set_title("Where the reference catalogue is incomplete", fontsize=FS_TITLE,
                   loc="left", pad=6)
    ax_c.grid(axis="y", lw=0.4, alpha=0.3, zorder=0)
    ax_c.spines[["top", "right"]].set_visible(False)
    # No ancestry legend here: panel (a) already establishes the colours, and repeating it was
    # what collided with this panel's title.
    h2 = [Patch(facecolor="#777777", alpha=FIELD_ALPHA[fc], edgecolor=FIELD_EDGE[fc],
                linewidth=0.35 if FIELD_EDGE[fc] != "none" else 0, label=FIELD_LABEL[fc])
          for fc in FIELD_ORDER]
    h2.append(Line2D([], [], color="#111111", lw=1.0, label="artifact rate (negative control)"))
    ax_c.legend(handles=h2, ncol=2, frameon=False, fontsize=FS_LEG, loc="upper left",
                handlelength=1.3, columnspacing=1.0, labelspacing=0.3)
    ax_c.set_ylim(0, max(1.0, pc[FIELD_ORDER].sum(axis=1).max()) * 1.34)
    panel_letter(ax_c, "c")

    # ---- d: novel proteins per gene, by recurrence ---------------------------------------
    top = pd_.head(16)
    xs = np.arange(len(top))
    cols = [("once", "#C9CFD6", "seen once"),
            ("mid", "#5B8FBF", "2–19 unrelated"),
            ("high", "#B4472E", "≥20 unrelated")]
    bottom = np.zeros(len(top))
    for key, col, lab in cols:
        v = top[key].to_numpy(dtype=float)
        ax_d.bar(xs, v, bottom=bottom, width=0.76, color=col, label=lab,
                 edgecolor="white", linewidth=0.4, zorder=3)
        bottom += v
    ax_d.set_xticks(xs)
    ax_d.set_xticklabels([g.replace("HLA-", "") for g in top["gene"]], rotation=90,
                         fontsize=FS_TICK)
    ax_d.tick_params(axis="y", labelsize=FS_TICK)
    ax_d.set_ylabel("distinct novel protein alleles", fontsize=FS_LABEL)
    ax_d.set_title("Novel proteins are rare, and not in the classical genes",
                   fontsize=FS_TITLE, loc="left", pad=6)
    ax_d.legend(frameon=False, fontsize=FS_LEG, loc="upper right")
    ax_d.grid(axis="y", lw=0.4, alpha=0.3, zorder=0)
    ax_d.spines[["top", "right"]].set_visible(False)
    panel_letter(ax_d, "d")

    note = ("a  all unrelated participants, predicted ancestry; the figure under each block is "
            "the share clearing the strict threshold.\n"
            "b, c  strict ancestry (probability ≥ %.2f); b shows only alleles with ≥ %d carriers.\n"
            "c  exact rates; a gene × ancestry cell with < %d called haplotypes is omitted.   "
            "d  carrier counts of 1–19 are suppressed at source, so 2–19 is one class."
            % (strict, min_carriers, SUPPRESS_BELOW))
    fig.text(0.02, 0.005, note, fontsize=FS_ANNOT, color="#666666", ha="left", va="bottom",
             linespacing=1.5)

    # No bbox_inches="tight": the GridSpec margins are explicit, and tight cropping was
    # re-flowing the canvas around the footnote and stretching every panel.
    fig.savefig(out_png, dpi=300, facecolor="white")
    fig.savefig(out_pdf, facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------------------
def run(args):
    os.makedirs(args.out_dir, exist_ok=True)
    m24 = _load("24_novelty_by_field.py", "novelty_by_field")

    print("[36] loading ...", flush=True)
    t1 = pd.read_csv(args.table1, sep="\t", dtype=str, low_memory=False)
    cohort = pd.read_csv(args.cohort_membership, sep="\t", dtype=str)
    clusters = pd.read_csv(args.clusters, sep="\t", dtype=str)

    people, n_removed = m24.build_people(
        sorted(set(t1["person_id"].astype(str))), args.cohort_membership,
        args.relatedness_table, args.kin_min, args.strict_threshold, args.skip_relatedness)
    keep = people[people["unrelated"]]
    anc_of = {p: a for p, a in zip(keep["person_id"].astype(str), keep["anc_strict"]) if a in ANC}
    people_keep = set(anc_of)
    print("[36] unrelated with strict ancestry (>=%.2f): %d" % (args.strict_threshold,
                                                                len(people_keep)), flush=True)

    all_unrelated = set(keep["person_id"].astype(str))
    pa, pa_cols = panel_a_data(cohort, all_unrelated, args.strict_threshold)
    pb = panel_b_data(t1, cohort, people_keep, args.gene, args.min_carriers,
                      args.strict_threshold)
    pc = panel_c_data(t1, anc_of, m24)
    pd_ = panel_d_data(t1, clusters, anc_of)
    print("[36] panels: a n=%d | b %d alleles | c %d cells | d %d genes"
          % (len(pa), len(pb), len(pc), len(pd_)), flush=True)

    # Aggregate, disclosure-safe exports so b/c/d can be restyled offline later.
    pb.drop(columns=["n_carriers"]).to_csv(
        os.path.join(args.out_dir, "panel_b_ternary_alleles.tsv"), sep="\t", index=False)
    pc.to_csv(os.path.join(args.out_dir, "panel_c_novelty_rates.tsv"), sep="\t", index=False)
    pd_.to_csv(os.path.join(args.out_dir, "panel_d_recurrence_by_gene.tsv"), sep="\t",
               index=False)
    (pa.groupby("anc").agg(n_people=("anc", "size"),
                           frac_passing_strict=("passes_strict", "mean")).reset_index()
       .to_csv(os.path.join(args.out_dir, "panel_a_group_sizes.tsv"), sep="\t", index=False))

    draw(pa, pa_cols, pb, pc, pd_,
         os.path.join(args.out_dir, "figure1.png"),
         os.path.join(args.out_dir, "figure1.pdf"),
         args.gene, args.strict_threshold, args.min_carriers, len(pa))

    with open(os.path.join(args.out_dir, "summary.json"), "w") as fh:
        json.dump({"n_people": len(pa), "strict_threshold": args.strict_threshold,
                   "gene_ternary": args.gene, "min_carriers_ternary": args.min_carriers,
                   "n_alleles_ternary": int(len(pb)),
                   "native_replot": True}, fh, indent=2)
    print("[36] done -> %s" % args.out_dir, flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=DEFAULT_TABLE1)
    ap.add_argument("--cohort-membership", default=DEFAULT_COHORT)
    ap.add_argument("--clusters", default=DEFAULT_CLUSTERS)
    ap.add_argument("--relatedness-table", default=DEFAULT_RELATEDNESS)
    ap.add_argument("--skip-relatedness", action="store_true")
    ap.add_argument("--kin-min", type=float, default=KIN_MIN)
    ap.add_argument("--strict-threshold", type=float, default=0.98,
                    help="ancestry-probability floor for panels a and b (ask A9)")
    ap.add_argument("--gene", default="HLA-B")
    ap.add_argument("--min-carriers", type=int, default=SUPPRESS_BELOW)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    run(ap.parse_args(argv))


if __name__ == "__main__":
    main()
