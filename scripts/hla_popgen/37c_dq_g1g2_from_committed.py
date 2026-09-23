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
def fig_g1g2_from_table(df_anc, path_stem, corner_label=""):
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
    n_a2 = sum(1 for a in a_alleles if dq_group("DQA1", a) == "G2")
    n_b2 = sum(1 for b in b_alleles if dq_group("DQB1", b) == "G2")

    M = np.full((len(a_alleles), len(b_alleles)), np.nan)
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

    marg_a = {a: have.loc[have["allele_a"] == a, "freq_a"].max() for a in a_alleles}
    marg_b = {b: have.loc[have["allele_b"] == b, "freq_b"].max() for b in b_alleles}

    with vc.nature_style():
        fig = plt.figure(figsize=(vc.mm(120), vc.mm(112)))
        gs = fig.add_gridspec(2, 2, width_ratios=[len(b_alleles), 5],
                               height_ratios=[3.5, len(a_alleles)], wspace=0.02, hspace=0.02)
        ax_heat = fig.add_subplot(gs[1, 0])
        ax_top = fig.add_subplot(gs[0, 0], sharex=ax_heat)
        ax_right = fig.add_subplot(gs[1, 1], sharey=ax_heat)

        cmap = vc.diverging_cmap()
        norm = vc.diverging_norm(-1.0, 1.0)
        ax_heat.imshow(M, cmap=cmap, norm=norm, aspect="auto", interpolation="none")
        for i in range(len(a_alleles)):
            for j in range(len(b_alleles)):
                if np.isnan(M[i, j]):
                    vc.hatch_suppressed(ax_heat, j - 0.5, i - 0.5, 1, 1)
                elif censored_dot[i, j]:
                    ax_heat.plot(j, i, marker="o", markersize=1.3, color="black", zorder=6)

        ax_heat.set_xticks(range(len(b_alleles)))
        ax_heat.set_xticklabels(b_alleles, rotation=90, fontsize=5)
        ax_heat.set_yticks(range(len(a_alleles)))
        ax_heat.set_yticklabels(a_alleles, fontsize=5)
        ax_heat.tick_params(length=2)
        ax_heat.set_xlabel("DQB1", fontsize=6)
        ax_heat.set_ylabel("DQA1", fontsize=6)
        ax_heat.set_xlim(-0.5, len(b_alleles) - 0.5)
        ax_heat.set_ylim(len(a_alleles) - 0.5, -0.5)

        ax_heat.axhline(n_a2 - 0.5, color="black", lw=0.8, zorder=5)
        ax_heat.axvline(n_b2 - 0.5, color="black", lw=0.8, zorder=5)

        if n_b2 < len(b_alleles) and n_a2 > 0:
            ax_heat.text((n_b2 + len(b_alleles)) / 2 - 0.5, n_a2 / 2 - 0.5,
                         "Predicted\nincompatible", ha="center", va="center", fontsize=5.5,
                         zorder=10, bbox=dict(boxstyle="round", fc="white", ec="none", alpha=0.85))
        if n_b2 > 0 and n_a2 < len(a_alleles):
            ax_heat.text(n_b2 / 2 - 0.5, (n_a2 + len(a_alleles)) / 2 - 0.5,
                         "Predicted\nincompatible", ha="center", va="center", fontsize=5.5,
                         zorder=10, bbox=dict(boxstyle="round", fc="white", ec="none", alpha=0.85))

        # Group-block axis labels (G1/G2), per orchestrator review. Pushed well clear of the
        # (variable-width) allele tick labels -- bbox_inches='tight' in save_fig expands the
        # canvas to fit, so a large negative offset costs nothing.
        label_x = -0.62
        ax_heat.text(label_x, (n_a2 - 1) / 2 if n_a2 else 0, "G2 α (DQA1*01)",
                    transform=ax_heat.get_yaxis_transform(), ha="center", va="center", fontsize=5,
                    rotation=90)
        if n_a2 < len(a_alleles):
            ax_heat.text(label_x, n_a2 + (len(a_alleles) - n_a2 - 1) / 2, "G1 α (DQA1*02–06)",
                        transform=ax_heat.get_yaxis_transform(), ha="center", va="center",
                        fontsize=5, rotation=90)
        ax_top.text((n_b2 - 1) / 2 if n_b2 else 0, 1.35, "G2 β (DQB1*05/06)",
                   transform=ax_top.get_xaxis_transform(), ha="center", va="bottom", fontsize=4.5)
        if n_b2 < len(b_alleles):
            ax_top.text(n_b2 + (len(b_alleles) - n_b2 - 1) / 2, 1.35, "G1 β (DQB1*02/03/04)",
                       transform=ax_top.get_xaxis_transform(), ha="center", va="bottom", fontsize=4.5)

        ax_top.bar(range(len(b_alleles)), [marg_b[b] for b in b_alleles], color="#AAAAAA",
                  width=0.6, linewidth=0)
        ax_top.set_ylabel("carrier\nfreq.", fontsize=4.5)
        ax_top.set_yticks([0, 0.1, 0.2, 0.3])
        ax_top.set_ylim(0, max(0.31, max(marg_b.values()) * 1.05))
        ax_top.tick_params(labelbottom=False, bottom=False, labelsize=4)
        ax_top.spines[["top", "right"]].set_visible(False)

        ax_right.barh(range(len(a_alleles)), [marg_a[a] for a in a_alleles], color="#AAAAAA",
                     height=0.6, linewidth=0)
        ax_right.set_xlabel("carrier\nfreq.", fontsize=4.5)
        ax_right.set_xticks([0, 0.1, 0.2, 0.3])
        ax_right.set_xlim(0, max(0.31, max(marg_a.values()) * 1.05))
        ax_right.tick_params(labelleft=False, left=False, labelsize=4)
        ax_right.spines[["top", "right"]].set_visible(False)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cbar = fig.colorbar(sm, ax=[ax_heat, ax_right], location="right", fraction=0.05,
                            pad=0.16, shrink=0.4, aspect=12)
        cbar.set_label("Signed phased D'", fontsize=5.5)
        cbar.ax.tick_params(labelsize=4.5)

        if corner_label:
            fig.text(0.01, 0.99, corner_label, ha="left", va="top", fontsize=5.5)
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
            ax.tick_params(labelsize=4.5)
            ax.spines[["top", "right"]].set_visible(False)
        fig.suptitle("Bimodality within compatible vs incompatible quadrants (pooled %s)"
                     % ", ".join(ancestries), fontsize=6.5)
        fig.subplots_adjust(top=0.78, wspace=0.15)
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

    pooled_df, pooled_N = pool_across_ancestries(df_all, ANCESTRY_ORDER, n_by_anc)
    if not pooled_df.empty:
        pooled_df["ancestry"] = "POOLED"
        main_label = "POOLED (all 6 ancestries, N=%d haplotypes)" % pooled_N
        main_stem = "fig_dq_g1g2_committed_MAIN_POOLED"
        print("[37c] main panel: pooled reconstruction, N=%d haplotypes, %d pairs (%d absent "
             "from every ancestry's own table -> hatched)"
             % (pooled_N, len(pooled_df),
                pooled_df["signed_Dprime"].isna().sum()), flush=True)
        fig_g1g2_from_table(pooled_df, os.path.join(args.out_dir, main_stem), corner_label=main_label)
    else:
        fig_g1g2_from_table(df_all[df_all["ancestry"] == main_anc],
                            os.path.join(args.out_dir, "fig_dq_g1g2_committed_MAIN_%s" % main_anc),
                            corner_label="%s, N=%d haplotypes" % (main_anc, n_by_anc[main_anc]))

    for anc in SUPPLEMENT_ANCESTRIES:
        fig_g1g2_from_table(df_all[df_all["ancestry"] == anc],
                            os.path.join(args.out_dir, "fig_dq_g1g2_committed_supp_%s" % anc),
                            corner_label="%s, N=%d haplotypes" % (anc, n_by_anc.get(anc) or 0))

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
