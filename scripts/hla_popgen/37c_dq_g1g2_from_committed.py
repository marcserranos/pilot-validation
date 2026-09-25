#!/usr/bin/env python3
"""37c -- DQ G1/G2 signed-D' figure, purge statistic (O/E), and bimodality check, computed
ENTIRELY from the already-committed `reports/hla_popgen/29_hla_ld/ld_pairwise.tsv`.

This is the `--from-committed` mode the task asked for: no VM access, no raw haplotype
recomputation. The committed table already carries, at 2-field (protein) resolution, per
ancestry, for DQA1~DQB1: `Dprime` (SIGNED -- verified empirically: D and Dprime always share
sign in this file, so 29's `pairwise_ld` was never sign-collapsed the way `37_dq_g1g2_signed_ld.py`
docstrings assumed when written without VM access), `freq_a`, `freq_b`, `freq_hap`, and
`n_hap_ij_disp` (`<20` = censored, `0` = genuinely not observed -- two distinct disclosure
statuses, never conflated). Critically, 29's `pairwise_ld` already required BOTH alleles in a row
to individually clear a 20-haplotype-carrier floor on their own marginal (`MIN_ALLELE_HAPS = 20`
in `29_hla_ld_by_ancestry.py`) before emitting any row for that pair -- so the "every row/col
allele needs >=20 carriers" floor this task asks for is already satisfied by construction; no
extra filtering is needed or applied here.

Reuses `dq_group`/`classify_pair` (Petersdorf 2022 G1/G2 rule) from
`37_dq_g1g2_signed_ld.py` verbatim via import -- does not redefine the rule.

Outputs (into reports/hla_popgen/37_dq_g1g2_signed_ld/):
  - fig_dq_g1g2_committed_MAIN_<ancestry>.png/.pdf -- Cole's target layout, main panel
    (ancestry with the most DQA1~DQB1 haplotypes among the estimable set).
  - fig_dq_g1g2_committed_supp_<ancestry>.png/.pdf -- one per ancestry in {AFR, AMR, EAS, EUR,
    SAS} (small multiple).
  - oe_purge_committed.tsv -- O/E interval (lower/upper bound from the censored-cell interval)
    per ancestry, cross-group ("predicted incompatible") haplotypes vs independence expectation.
  - recurrent_cross_group_pairs.tsv -- cross-group (predicted_incompatible) allele pairs that
    clear the 20-haplotype floor (status == estimated, i.e. Cole's "light blue" cells) in at
    least one ancestry, with which ancestries and their signed D'.
  - fig_bimodality_committed.png/.pdf -- histogram of signed D' within compatible (G1/G1, G2/G2)
    vs incompatible (cross-group) quadrants, pooled over all estimable ancestries.

Usage (local, no VM):
    python3 scripts/hla_popgen/37c_dq_g1g2_from_committed.py
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import _viz_common as vc  # noqa: E402


def _load_module(filename, modname):
    path = os.path.join(_THIS_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


_m37 = None


def m37():
    """Reuse dq_group/classify_pair (Petersdorf 2022 G1/G2 rule) verbatim, per task instructions."""
    global _m37
    if _m37 is None:
        _m37 = _load_module("37_dq_g1g2_signed_ld.py", "dq_g1g2_signed_ld")
    return _m37


DEFAULT_LD_PAIRWISE = os.path.join(
    _THIS_DIR, "..", "..", "reports", "hla_popgen", "29_hla_ld", "ld_pairwise.tsv")
DEFAULT_OUT_DIR = os.path.join(
    _THIS_DIR, "..", "..", "reports", "hla_popgen", "37_dq_g1g2_signed_ld")

SUPPRESS_BELOW = 20
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
SUPPLEMENT_ANCESTRIES = ["AFR", "AMR", "EAS", "EUR", "SAS"]  # MID excluded -- too few haplotypes
# to be usefully "estimable" per the task's own list; still reported in the O/E table with a
# caveat rather than silently dropped from every output.


# ---------------------------------------------------------------------------
# Load + annotate
# ---------------------------------------------------------------------------
def load_committed(path=DEFAULT_LD_PAIRWISE, pair="DQA1~DQB1"):
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df[df["pair"] == pair].copy()
    for col in ("freq_a", "freq_b", "freq_hap", "D", "Dprime", "r2", "chi2_1df"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    def _status(v):
        if v == "<20":
            return "suppressed_lt20"
        if v == "0":
            return "not_observed"
        return "estimated"

    df["status"] = df["n_hap_ij_disp"].map(_status)
    df["n_hap_ij"] = df["n_hap_ij_disp"].map(
        lambda v: np.nan if v == "<20" else int(v))
    df["signed_Dprime"] = df["Dprime"]
    return df


def estimate_n_haplotypes(df_anc):
    """N (total phased haplotypes) for this ancestry, recovered from n_hap_ij / freq_hap on any
    uncensored, nonzero row -- exact and internally consistent (verified: single value per
    ancestry across all its rows, to floating-point precision)."""
    sub = df_anc[(df_anc["status"] == "estimated") & (df_anc["freq_hap"] > 0)]
    if sub.empty:
        return None
    ests = sub["n_hap_ij"] / sub["freq_hap"]
    return int(round(ests.median()))


# ---------------------------------------------------------------------------
# Figure: Cole's target layout, built directly from the committed table's rows for one ancestry
# (or the pooled reconstruction -- see `pool_across_ancestries`).
#
# IMPORTANT correction (orchestrator review, 2026-09-23): D' is NOT undefined just because a
# pair's disclosed count is `0` ("not_observed") or `<20` ("suppressed_lt20"). `freq_hap`,
# `freq_a`, `freq_b` are disclosed as exact floats in the committed table regardless of the
# haplotype-count disclosure status (verified: every DQA1~DQB1 row, all 3 statuses, has a
# non-NaN `Dprime` -- 0 NaNs out of 927 rows), so `signed_Dprime` is a real, already-public value
# for every row that exists in the table at all. A cell should be hatched ONLY if the allele pair
# is entirely ABSENT from the table for this ancestry (both alleles individually clear the
# 20-haplotype floor but never appear together as a row -- this happens for a handful of the
# union-of-ancestries pairs used in the pooled reconstruction, never within a single ancestry's
# own dense grid). A `suppressed_lt20` cell IS coloured by its real D' but gets a small dot marker
# so the reader knows its count is disclosed only as "<20", not as an exact integer.
# ---------------------------------------------------------------------------
# Fill for the "predicted incompatible" quadrants (S04 WS-C phase 2 redesign). Distinct from
# vc.SUPPRESSED_COLOR ("#D9D9D9", always hatched, means "disclosure-censored data") -- this is a
# plain, unhatched, slightly lighter grey meaning "structurally uniform, nothing to see here": the
# data (recurrent_cross_group_pairs.tsv / oe_purge_committed.tsv) show D'=-1 (0 haplotypes
# observed) for every cross-group cell in every ancestry, so painting this quadrant on the same
# red-white-blue scale as the informative G1/G1, G2/G2 blocks wastes ink and a reader's attention
# on a foregone, uniform result. A pair that's absent from the table entirely still gets the usual
# hatched SUPPRESSED_COLOR treatment, drawn on top of this fill, so the two "nothing here" reasons
# (confirmed-zero vs never-cleared-the-disclosure-floor) stay visually distinct.
NEUTRAL_INCOMPAT_COLOR = "#EDEDED"


def _bracket_v(ax, y0, y1, label, x=0.30, cap=0.14, fontsize=5.5):
    """Vertical G1/G2 group bracket (a plain square-bracket '[' shape) spanning data-y range
    [y0, y1] on `ax`, with `label` set just to its left, rotated 90. Used for the DQA1 (row) axis,
    placed in a narrow dedicated axes to the LEFT of the heatmap so it sits directly adjacent to
    the row tick labels instead of floating off in a large negative-axes-fraction offset (the
    fault CRITIC_WSC.md flagged: 'sit far outside the plot'). `x` is deliberately on the FAR side
    of this axes (away from the heatmap, small x) -- the heatmap's own y tick labels hug the
    heatmap's left edge (large x in this axes' 0-1 coordinate space) and need the x>0.3 region of
    this column clear to render into without colliding with the bracket."""
    ax.plot([x, x], [y0, y1], color="black", lw=0.6, clip_on=False, solid_capstyle="butt")
    ax.plot([x, x - cap], [y0, y0], color="black", lw=0.6, clip_on=False)
    ax.plot([x, x - cap], [y1, y1], color="black", lw=0.6, clip_on=False)
    ax.text(x - cap - 0.10, (y0 + y1) / 2.0, label, ha="right", va="center", fontsize=fontsize,
            rotation=90)


def _bracket_h(ax, x0, x1, label, y=0.55, cap=0.16, fontsize=5.5):
    """Horizontal G1/G2 group bracket for the DQB1 (column) axis, in a narrow dedicated axes
    BELOW the heatmap (below the rotated column tick labels). `y` is on the FAR side of this
    axes (away from the heatmap, small-ish y after leaving headroom above for the tick labels
    that hug the heatmap's bottom edge, i.e. large y in this axes) -- same rationale as
    `_bracket_v`."""
    ax.plot([x0, x1], [y, y], color="black", lw=0.6, clip_on=False, solid_capstyle="butt")
    ax.plot([x0, x0], [y, y - cap], color="black", lw=0.6, clip_on=False)
    ax.plot([x1, x1], [y, y - cap], color="black", lw=0.6, clip_on=False)
    ax.text((x0 + x1) / 2.0, y - cap - 0.12, label, ha="center", va="top", fontsize=fontsize)


def fig_g1g2_from_table(df_anc, path_stem, corner_label="", incompatible_note=None):
    """Redesigned S04 WS-C phase 2 layout (was: independently-sized marginal axes that didn't
    line up with the heatmap columns/rows, a colorbar floating far to the right with large empty
    gaps, G1/G2 labels pushed to a large negative-axes-fraction offset outside the panel, a
    sentence-length corner label, and two solid dark-blue 'predicted incompatible' blocks that
    dominate the figure's ink despite carrying zero information -- see CRITIC_WSC.md).

    Fix: one GridSpec, sharex/sharey between the heatmap and both marginals (so alignment is
    structural, not a matched-width/height coincidence), a compact same-row colorbar axes, G1/G2
    brackets in dedicated thin axes adjacent to the tick labels, and the incompatible quadrants
    flattened to a plain neutral fill (see NEUTRAL_INCOMPAT_COLOR) so the eye goes to the
    informative G1/G1, G2/G2 blocks. `check_layout()` (via `save_fig`, strict=True) is the
    acceptance test for the layout fixes; visual encoding choices were iterated by re-reading the
    rendered PNG.
    """
    import matplotlib.pyplot as plt

    m37mod = m37()
    dq_group = m37mod.dq_group

    have = df_anc[df_anc["signed_Dprime"].notna()]
    if have.empty:
        return None

    a_alleles = sorted(
        have["allele_a"].unique(),
        key=lambda a: (dq_group("DQA1", a) != "G2",
                        -have.loc[have["allele_a"] == a, "freq_a"].max()))
    b_alleles = sorted(
        have["allele_b"].unique(),
        key=lambda b: (dq_group("DQB1", b) != "G2",
                        -have.loc[have["allele_b"] == b, "freq_b"].max()))
    n_a, n_b = len(a_alleles), len(b_alleles)
    n_a2 = sum(1 for a in a_alleles if dq_group("DQA1", a) == "G2")
    n_b2 = sum(1 for b in b_alleles if dq_group("DQB1", b) == "G2")

    M = np.full((n_a, n_b), np.nan)
    censored_dot = np.zeros_like(M, dtype=bool)
    look = {(r["allele_a"], r["allele_b"]): r["signed_Dprime"] for _, r in have.iterrows()}
    status_look = {(r["allele_a"], r["allele_b"]): r.get("status") for _, r in df_anc.iterrows()}
    for i, a in enumerate(a_alleles):
        for j, b in enumerate(b_alleles):
            v = look.get((a, b))
            if v is not None and not np.isnan(v):
                M[i, j] = v
                if status_look.get((a, b)) == "suppressed_lt20":
                    censored_dot[i, j] = True
            # else: genuinely absent from the table for this ancestry -- stays NaN -> hatched.

    # Incompatible-quadrant mask: (G2-alpha row, G1-beta col) or (G1-alpha row, G2-beta col).
    incompat = np.zeros_like(M, dtype=bool)
    incompat[:n_a2, n_b2:] = True
    incompat[n_a2:, :n_b2] = True
    Mdisplay = np.where(incompat, np.nan, M)  # incompatible cells never get D'-scale color

    marg_a = {a: have.loc[have["allele_a"] == a, "freq_a"].max() for a in a_alleles}
    marg_b = {b: have.loc[have["allele_b"] == b, "freq_b"].max() for b in b_alleles}

    width_mm = vc.NATURE_DOUBLE_COL_MM  # 183mm -- many alleles need the double-column width
    height_mm = float(np.clip(58 + 6.2 * n_a, 95, 165))

    with vc.nature_style():
        fig = plt.figure(figsize=(vc.mm(width_mm), vc.mm(height_mm)))
        gs = fig.add_gridspec(
            3, 4,
            width_ratios=[2.8, n_b, 1.9, 0.42],
            height_ratios=[1.9, n_a, 2.3],
            wspace=0.08, hspace=0.10,
            left=0.10, right=0.90, top=0.95, bottom=0.05)
        ax_heat = fig.add_subplot(gs[1, 1])
        ax_top = fig.add_subplot(gs[0, 1], sharex=ax_heat)
        ax_right = fig.add_subplot(gs[1, 2], sharey=ax_heat)
        ax_cbar = fig.add_subplot(gs[1, 3])
        ax_brk_y = fig.add_subplot(gs[1, 0], sharey=ax_heat)
        ax_brk_x = fig.add_subplot(gs[2, 1], sharex=ax_heat)
        vc.mark_marginal(ax_top, ax_heat, axis="x")
        vc.mark_marginal(ax_right, ax_heat, axis="y")

        cmap = vc.diverging_cmap()
        norm = vc.diverging_norm(-1.0, 1.0)

        # Neutral fill for the two incompatible quadrants FIRST (plain rectangles, not imshow --
        # a flat single patch per quadrant is less ink than per-cell color for a uniform result).
        if n_b2 < n_b and n_a2 > 0:
            ax_heat.add_patch(plt.Rectangle((n_b2 - 0.5, -0.5), n_b - n_b2, n_a2,
                                             facecolor=NEUTRAL_INCOMPAT_COLOR, edgecolor="none",
                                             zorder=1))
        if n_b2 > 0 and n_a2 < n_a:
            ax_heat.add_patch(plt.Rectangle((-0.5, n_a2 - 0.5), n_b2, n_a - n_a2,
                                             facecolor=NEUTRAL_INCOMPAT_COLOR, edgecolor="none",
                                             zorder=1))
        ax_heat.imshow(Mdisplay, cmap=cmap, norm=norm, aspect="auto", interpolation="none",
                       zorder=2)
        for i in range(n_a):
            for j in range(n_b):
                if np.isnan(M[i, j]):
                    vc.hatch_suppressed(ax_heat, j - 0.5, i - 0.5, 1, 1, zorder=3)
                elif censored_dot[i, j]:
                    ax_heat.plot(j, i, marker="o", markersize=1.3, color="black", zorder=6)

        ax_heat.set_xticks(range(n_b))
        ax_heat.set_xticklabels(b_alleles, rotation=90, fontsize=5)
        ax_heat.set_yticks(range(n_a))
        ax_heat.set_yticklabels(a_alleles, fontsize=5)
        ax_heat.tick_params(length=2)
        ax_heat.set_xlabel("DQB1", fontsize=6)
        ax_heat.set_ylabel("DQA1", fontsize=6)
        ax_heat.set_xlim(-0.5, n_b - 0.5)
        ax_heat.set_ylim(n_a - 0.5, -0.5)
        ax_heat.axhline(n_a2 - 0.5, color="black", lw=0.8, zorder=5)
        ax_heat.axvline(n_b2 - 0.5, color="black", lw=0.8, zorder=5)

        # One small annotation carries the headline (de-AI checklist item 8: direct label, no
        # legend) instead of duplicating "predicted incompatible" text across both quadrants.
        note = incompatible_note or "predicted incompatible\n(D′≈−1)"
        if n_b2 < n_b and n_a2 > 0:
            ax_heat.text((n_b2 + n_b) / 2 - 0.5, n_a2 / 2 - 0.5, note, ha="center", va="center",
                         fontsize=5, color="#555555", zorder=4)
        elif n_b2 > 0 and n_a2 < n_a:
            ax_heat.text(n_b2 / 2 - 0.5, (n_a2 + n_a) / 2 - 0.5, note, ha="center", va="center",
                         fontsize=5, color="#555555", zorder=4)

        # G1/G2 brackets, adjacent to the tick labels (dedicated thin axes, not a large offset).
        for a in (ax_brk_y, ax_brk_x):
            a.set_axis_off()
        ax_brk_y.set_xlim(0, 1)
        if n_a2 > 0:
            _bracket_v(ax_brk_y, -0.5, n_a2 - 0.5, "G2")
        if n_a2 < n_a:
            _bracket_v(ax_brk_y, n_a2 - 0.5, n_a - 0.5, "G1")
        ax_brk_x.set_ylim(0, 1)
        if n_b2 > 0:
            _bracket_h(ax_brk_x, -0.5, n_b2 - 0.5, "G2")
        if n_b2 < n_b:
            _bracket_h(ax_brk_x, n_b2 - 0.5, n_b - 0.5, "G1")

        ax_top.bar(range(n_b), [marg_b[b] for b in b_alleles], color="#AAAAAA",
                  width=0.7, linewidth=0)
        ax_top.set_ylabel("carrier\nfreq.", fontsize=4.5)
        top_max = max(0.31, max(marg_b.values()) * 1.05)
        ax_top.set_yticks([0, round(top_max / 2, 2), round(top_max, 2)])
        ax_top.set_ylim(0, top_max)
        ax_top.tick_params(labelbottom=False, bottom=False, labelsize=4)
        ax_top.spines[["top", "right"]].set_visible(False)

        ax_right.barh(range(n_a), [marg_a[a] for a in a_alleles], color="#AAAAAA",
                     height=0.7, linewidth=0)
        ax_right.set_xlabel("carrier\nfreq.", fontsize=4.5)
        right_max = max(0.31, max(marg_a.values()) * 1.05)
        ax_right.set_xticks([0, round(right_max / 2, 2), round(right_max, 2)])
        ax_right.set_xlim(0, right_max)
        ax_right.tick_params(labelleft=False, left=False, labelsize=4)
        ax_right.spines[["top", "right"]].set_visible(False)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cbar = fig.colorbar(sm, cax=ax_cbar)
        cbar.set_label("Signed phased D′", fontsize=5.5)
        cbar.ax.tick_params(labelsize=4.5, length=2)

        # Short panel label (not a sentence) -- N and full context live in the report README.
        if corner_label:
            fig.text(0.005, 0.995, corner_label, ha="left", va="top", fontsize=6,
                     fontweight="bold")
    return vc.save_fig(fig, path_stem)


# ---------------------------------------------------------------------------
# Purge statistic: O/E for cross-group ("predicted incompatible") haplotypes, censored-interval.
# ---------------------------------------------------------------------------
def oe_interval_from_table(df_anc, n_haps):
    m37mod = m37()
    classify_pair = m37mod.classify_pair

    df = df_anc.copy()
    df["group"] = [classify_pair(a, b) for a, b in zip(df["allele_a"], df["allele_b"])]
    cross = df[df["group"] == "predicted_incompatible"]
    if cross.empty or n_haps is None:
        return None

    n_estimated = int((cross["status"] == "estimated").sum())
    n_not_observed = int((cross["status"] == "not_observed").sum())
    n_censored = int((cross["status"] == "suppressed_lt20").sum())

    obs_estimated = cross.loc[cross["status"] == "estimated", "n_hap_ij"].sum()
    observed_lower = float(obs_estimated) + 0 * n_not_observed + 0 * n_censored
    observed_upper = float(obs_estimated) + 0 * n_not_observed + 19 * n_censored

    expected = float((cross["freq_a"] * cross["freq_b"]).sum() * n_haps)

    return {
        "n_haplotypes": n_haps,
        "n_cross_group_cells": len(cross),
        "n_cells_estimated": n_estimated,
        "n_cells_not_observed": n_not_observed,
        "n_cells_censored_lt20": n_censored,
        "observed_lower": observed_lower,
        "observed_upper": observed_upper,
        "expected": round(expected, 2),
        "oe_lower": round(observed_lower / expected, 4) if expected > 0 else float("nan"),
        "oe_upper": round(observed_upper / expected, 4) if expected > 0 else float("nan"),
    }


def recurrent_cross_group_pairs(df_all):
    """Cross-group (predicted_incompatible) allele pairs with status == estimated (i.e. an exact
    disclosed count of >=20 -- Cole's 'light blue' cells, defined here by the same disclosure
    threshold as everywhere else in this pipeline, NOT merely 'D' prime is not exactly -1', since
    at 2-field every not_observed cross-group cell has D'==-1.0 exactly -- see README) in >=1
    ancestry, tidy one-row-per-(pair, ancestry)."""
    m37mod = m37()
    classify_pair = m37mod.classify_pair
    df = df_all.copy()
    df["group"] = [classify_pair(a, b) for a, b in zip(df["allele_a"], df["allele_b"])]
    hits = df[(df["group"] == "predicted_incompatible") & (df["status"] == "estimated")]
    cols = ["allele_a", "allele_b", "ancestry", "n_hap_ij", "freq_a", "freq_b", "signed_Dprime"]
    return hits[cols].sort_values(["allele_a", "allele_b", "ancestry"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Pooled-across-ancestries reconstruction. D' depends only on allele/haplotype frequencies, and
# pooled frequencies are exact n-weighted means of the per-ancestry frequencies (both disclosed
# as floats regardless of the haplotype-count disclosure status) -- so pooled D' can be
# reconstructed exactly, ancestry-marginal-floor edge cases aside (see caveat below).
# ---------------------------------------------------------------------------
def pool_across_ancestries(df_all, ancestries=ANCESTRY_ORDER, n_by_anc=None):
    """Returns (pooled_df, pooled_N). pooled_df has one row per (allele_a, allele_b) pair that
    appears in >=1 ancestry's table, with pooled freq_a/freq_b/freq_hap (N_anc-weighted sums,
    ancestries where the pair/allele is absent contribute 0 -- a slight undercount only possible
    when an allele's own marginal in that ancestry is <20, i.e. the omitted contribution itself
    is <20/pooled_N, negligible at these cohort sizes) and pooled signed_Dprime recomputed from
    those pooled frequencies via the same Lewontin formula as 29/37. status is 'pooled' for every
    row (no per-cell count-based censoring applies to a quantity built from disclosed floats);
    n_ancestries_observed records how many ancestries actually contributed to each cell, for
    transparency."""
    if n_by_anc is None:
        n_by_anc = {a: estimate_n_haplotypes(df_all[df_all["ancestry"] == a]) for a in ancestries}
    sub = df_all[df_all["ancestry"].isin(ancestries)].copy()
    sub["N_anc"] = sub["ancestry"].map(n_by_anc)
    sub = sub[sub["N_anc"].notna()]
    pooled_N = float(sub.drop_duplicates("ancestry")["N_anc"].sum())
    if pooled_N <= 0:
        return pd.DataFrame(), 0

    a_marg = sub.drop_duplicates(["allele_a", "ancestry"])
    a_pool = (a_marg.assign(contrib=a_marg["N_anc"] * a_marg["freq_a"])
              .groupby("allele_a")["contrib"].sum() / pooled_N)
    b_marg = sub.drop_duplicates(["allele_b", "ancestry"])
    b_pool = (b_marg.assign(contrib=b_marg["N_anc"] * b_marg["freq_b"])
              .groupby("allele_b")["contrib"].sum() / pooled_N)

    grp = sub.assign(contrib=sub["N_anc"] * sub["freq_hap"]).groupby(["allele_a", "allele_b"])
    pooled_pij = grp["contrib"].sum() / pooled_N
    n_anc_obs = grp["ancestry"].nunique()

    rows = []
    for (a, b), pij in pooled_pij.items():
        pi, qj = a_pool.get(a), b_pool.get(b)
        if pi is None or qj is None or pi <= 0 or qj <= 0:
            continue
        D = pij - pi * qj
        if D >= 0:
            dmax = min(pi * (1 - qj), (1 - pi) * qj)
        else:
            dmax = min(pi * qj, (1 - pi) * (1 - qj))
        dprime = (D / dmax) if dmax > 0 else float("nan")
        rows.append({"allele_a": a, "allele_b": b, "freq_a": pi, "freq_b": qj, "freq_hap": pij,
                    "D": D, "signed_Dprime": dprime, "status": "pooled",
                    "n_ancestries_observed": int(n_anc_obs.loc[(a, b)])})
    return pd.DataFrame(rows), pooled_N


# ---------------------------------------------------------------------------
# Bimodality: signed D' within compatible vs incompatible quadrants, pooled over ancestries.
# ---------------------------------------------------------------------------
def fig_bimodality(df_all, path_stem, ancestries=SUPPLEMENT_ANCESTRIES):
    import matplotlib.pyplot as plt

    m37mod = m37()
    classify_pair = m37mod.classify_pair

    df = df_all[df_all["ancestry"].isin(ancestries) & df_all["signed_Dprime"].notna()].copy()
    df["group"] = [classify_pair(a, b) for a, b in zip(df["allele_a"], df["allele_b"])]
    compatible = df.loc[df["group"].isin(["G1", "G2"]), "signed_Dprime"].dropna().to_numpy()
    incompatible = df.loc[df["group"] == "predicted_incompatible", "signed_Dprime"].dropna().to_numpy()

    with vc.nature_style():
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(vc.mm(vc.NATURE_SINGLE_COL_MM), vc.mm(60)),
                                       sharey=True)
        bins = np.linspace(-1, 1, 41)
        ax1.hist(compatible, bins=bins, density=True, color="#B2182B")
        ax1.set_title("compatible\n(G1/G1, G2/G2)\nn=%d" % len(compatible), fontsize=5.5)
        ax1.set_ylabel("density", fontsize=5.5)
        if len(incompatible):
            ax2.hist(incompatible, bins=bins, density=True, color="#2166AC")
            ax2.set_title("predicted\nincompatible\nn=%d" % len(incompatible), fontsize=5.5)
        else:
            ax2.text(0.5, 0.5, "n=0 cells clear\nthe 20-haplotype\nfloor in any ancestry\n"
                    "(complete purge)", transform=ax2.transAxes, ha="center", va="center",
                    fontsize=5, color="#2166AC")
            ax2.set_title("predicted\nincompatible\nn=0", fontsize=5.5)
        for ax in (ax1, ax2):
            ax.set_xlabel("Signed D'", fontsize=5.5)
            # Fixed 3-tick scale (was matplotlib's auto ~7 ticks at 0.5 spacing, which crowded
            # into overlapping labels at this panel's narrow single-column width -- CRITIC_WSC.md
            # / S04 WS-C phase 2 fix, caught by check_layout()'s text_overlap check).
            ax.set_xticks([-1, 0, 1])
            ax.tick_params(labelsize=4.5)
            ax.spines[["top", "right"]].set_visible(False)
        # Short panel label, not a sentence (de-AI checklist item 7) -- the full description
        # ("pooled across AFR/AMR/EAS/EUR/SAS") lives in the report README/caption.
        fig.suptitle("DQ G1/G2 bimodality", fontsize=6.5, fontweight="bold", y=0.99)
        fig.subplots_adjust(top=0.74, wspace=0.15)
    return vc.save_fig(fig, path_stem), {
        "n_compatible": int(len(compatible)),
        "n_incompatible": int(len(incompatible)),
        "compatible_median": float(np.median(compatible)) if len(compatible) else float("nan"),
        "compatible_frac_near_plus1": float(np.mean(compatible >= 0.9)) if len(compatible) else float("nan"),
        "incompatible_median": float(np.median(incompatible)) if len(incompatible) else float("nan"),
        "incompatible_frac_near_minus1": float(np.mean(incompatible <= -0.9)) if len(incompatible) else float("nan"),
    }


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def run(args):
    vc.ensure_dir(args.out_dir)
    df_all = load_committed(args.ld_pairwise)

    n_by_anc = {anc: estimate_n_haplotypes(df_all[df_all["ancestry"] == anc])
               for anc in ANCESTRY_ORDER}
    print("[37c] N haplotypes per ancestry (DQA1~DQB1, 2-field):", n_by_anc, flush=True)

    main_anc = max((a for a in n_by_anc if n_by_anc[a]), key=lambda a: n_by_anc[a])
    print("[37c] most-haplotypes single ancestry: %s (N=%d)" % (main_anc, n_by_anc[main_anc]),
         flush=True)

    # Compute the O/E purge stats BEFORE the figures, so the pooled/per-ancestry headline
    # ("0 of N predicted-incompatible pairs observed") can be printed directly on the panel
    # instead of relying on a reader to cross-reference oe_purge_committed.tsv separately.
    oe_rows = []
    for anc in ANCESTRY_ORDER:
        res = oe_interval_from_table(df_all[df_all["ancestry"] == anc], n_by_anc.get(anc))
        if res:
            res["ancestry"] = anc
            oe_rows.append(res)
    oe_df = pd.DataFrame(oe_rows)[
        ["ancestry", "n_haplotypes", "n_cross_group_cells", "n_cells_estimated",
        "n_cells_not_observed", "n_cells_censored_lt20", "observed_lower", "observed_upper",
        "expected", "oe_lower", "oe_upper"]]
    total_estimated = int(oe_df["n_cells_estimated"].sum())
    total_cross = int(oe_df["n_cross_group_cells"].sum())

    pooled_df, pooled_N = pool_across_ancestries(df_all, ANCESTRY_ORDER, n_by_anc)
    if not pooled_df.empty:
        pooled_df["ancestry"] = "POOLED"
        main_label = "DQ G1/G2, pooled"
        main_stem = "fig_dq_g1g2_committed_MAIN_POOLED"
        print("[37c] main panel: pooled reconstruction, N=%d haplotypes, %d pairs (%d absent "
             "from every ancestry's own table -> hatched)"
             % (pooled_N, len(pooled_df),
                pooled_df["signed_Dprime"].isna().sum()), flush=True)
        pooled_note = ("predicted incompatible\n0/%d observed\n(all 6 ancestries)" % total_cross
                       if total_estimated == 0 else
                       "predicted incompatible\n%d/%d observed" % (total_estimated, total_cross))
        fig_g1g2_from_table(pooled_df, os.path.join(args.out_dir, main_stem),
                            corner_label=main_label, incompatible_note=pooled_note)
    else:
        fig_g1g2_from_table(df_all[df_all["ancestry"] == main_anc],
                            os.path.join(args.out_dir, "fig_dq_g1g2_committed_MAIN_%s" % main_anc),
                            corner_label="DQ G1/G2, %s" % main_anc)

    for anc in SUPPLEMENT_ANCESTRIES:
        anc_row = oe_df[oe_df["ancestry"] == anc]
        if not anc_row.empty:
            est = int(anc_row["n_cells_estimated"].iloc[0])
            cross = int(anc_row["n_cross_group_cells"].iloc[0])
            anc_note = "predicted incompatible\n%d/%d observed" % (est, cross)
        else:
            anc_note = None
        fig_g1g2_from_table(df_all[df_all["ancestry"] == anc],
                            os.path.join(args.out_dir, "fig_dq_g1g2_committed_supp_%s" % anc),
                            corner_label="DQ G1/G2, %s" % anc, incompatible_note=anc_note)

    oe_df.to_csv(os.path.join(args.out_dir, "oe_purge_committed.tsv"), sep="\t", index=False)
    print(oe_df.to_string(index=False), flush=True)

    recur = recurrent_cross_group_pairs(df_all)
    recur.to_csv(os.path.join(args.out_dir, "recurrent_cross_group_pairs.tsv"), sep="\t",
                index=False)
    print("[37c] recurrent cross-group (>=20 hap) pairs: %d rows across %d unique allele pairs"
         % (len(recur), recur[["allele_a", "allele_b"]].drop_duplicates().shape[0]), flush=True)

    _, bimod_stats = fig_bimodality(df_all, os.path.join(args.out_dir, "fig_bimodality_committed"))
    print("[37c] bimodality stats:", bimod_stats, flush=True)
    pd.Series(bimod_stats).to_json(os.path.join(args.out_dir, "bimodality_stats.json"), indent=2)

    print("[37c] done -> %s" % args.out_dir, flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ld-pairwise", default=DEFAULT_LD_PAIRWISE)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    args = ap.parse_args(argv)
    run(args)


if __name__ == "__main__":
    main()
