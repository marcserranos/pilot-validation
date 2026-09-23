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
# Figure: Cole's target layout, built directly from the committed table's rows for one ancestry.
# A local (not 37's) renderer, because 37's `fig_g1g2_signed_ld` has a latent gap for this input:
# a pair present in the table with status suppressed_lt20/not_observed is in `all_pairs` but not
# in `look`, so it is neither plotted with a real D' value NOR hatched -- it silently renders as
# an empty (NaN) cell. Task requires suppressed/not-estimable cells to be visibly hatched grey, so
# this renderer treats "not estimated" (whether censored, not_observed, or simply absent from the
# table -- e.g. an allele pair below the per-allele floor) as one hatched category.
# ---------------------------------------------------------------------------
def fig_g1g2_from_table(df_anc, path_stem, title_suffix=""):
    import matplotlib.pyplot as plt

    m37mod = m37()
    dq_group = m37mod.dq_group

    estimated = df_anc[df_anc["status"] == "estimated"]
    if estimated.empty:
        return None

    a_alleles = sorted(
        estimated["allele_a"].unique(),
        key=lambda a: (dq_group("DQA1", a) != "G2",
                        -estimated.loc[estimated["allele_a"] == a, "freq_a"].max()))
    b_alleles = sorted(
        estimated["allele_b"].unique(),
        key=lambda b: (dq_group("DQB1", b) != "G2",
                        -estimated.loc[estimated["allele_b"] == b, "freq_b"].max()))
    n_a2 = sum(1 for a in a_alleles if dq_group("DQA1", a) == "G2")
    n_b2 = sum(1 for b in b_alleles if dq_group("DQB1", b) == "G2")

    M = np.full((len(a_alleles), len(b_alleles)), np.nan)
    hatched = np.zeros_like(M, dtype=bool)
    look = {(r["allele_a"], r["allele_b"]): r["signed_Dprime"] for _, r in estimated.iterrows()}
    for i, a in enumerate(a_alleles):
        for j, b in enumerate(b_alleles):
            v = look.get((a, b))
            if v is not None and not np.isnan(v):
                M[i, j] = v
            else:
                hatched[i, j] = True  # covers suppressed_lt20, not_observed, and simply-absent

    marg_a = {a: estimated.loc[estimated["allele_a"] == a, "freq_a"].max() for a in a_alleles}
    marg_b = {b: estimated.loc[estimated["allele_b"] == b, "freq_b"].max() for b in b_alleles}

    with vc.nature_style():
        fig = plt.figure(figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(170)))
        gs = fig.add_gridspec(2, 2, width_ratios=[len(b_alleles), 6],
                               height_ratios=[4, len(a_alleles)], wspace=0.02, hspace=0.02)
        ax_heat = fig.add_subplot(gs[1, 0])
        ax_top = fig.add_subplot(gs[0, 0], sharex=ax_heat)
        ax_right = fig.add_subplot(gs[1, 1], sharey=ax_heat)

        cmap = vc.diverging_cmap()
        norm = vc.diverging_norm(-1.0, 1.0)
        ax_heat.imshow(M, cmap=cmap, norm=norm, aspect="auto", interpolation="none")
        for i in range(len(a_alleles)):
            for j in range(len(b_alleles)):
                if hatched[i, j]:
                    vc.hatch_suppressed(ax_heat, j - 0.5, i - 0.5, 1, 1)

        ax_heat.set_xticks(range(len(b_alleles)))
        ax_heat.set_xticklabels(b_alleles, rotation=90, fontsize=5)
        ax_heat.set_yticks(range(len(a_alleles)))
        ax_heat.set_yticklabels(a_alleles, fontsize=5)
        ax_heat.set_xlabel("DQB1")
        ax_heat.set_ylabel("DQA1")
        ax_heat.set_xlim(-0.5, len(b_alleles) - 0.5)
        ax_heat.set_ylim(len(a_alleles) - 0.5, -0.5)

        ax_heat.axhline(n_a2 - 0.5, color="black", lw=1.0, zorder=5)
        ax_heat.axvline(n_b2 - 0.5, color="black", lw=1.0, zorder=5)

        if n_b2 < len(b_alleles) and n_a2 > 0:
            ax_heat.text((n_b2 + len(b_alleles)) / 2 - 0.5, n_a2 / 2 - 0.5,
                         "Predicted\nincompatible", ha="center", va="center", fontsize=6,
                         zorder=10, bbox=dict(boxstyle="round", fc="white", ec="none", alpha=0.85))
        if n_b2 > 0 and n_a2 < len(a_alleles):
            ax_heat.text(n_b2 / 2 - 0.5, (n_a2 + len(a_alleles)) / 2 - 0.5,
                         "Predicted\nincompatible", ha="center", va="center", fontsize=6,
                         zorder=10, bbox=dict(boxstyle="round", fc="white", ec="none", alpha=0.85))

        ax_top.bar(range(len(b_alleles)), [marg_b[b] for b in b_alleles], color="#555555",
                  width=0.7)
        ax_top.set_ylabel("carrier\nfreq.", fontsize=5)
        ax_top.tick_params(labelbottom=False, bottom=False, labelsize=4)
        ax_top.spines[["top", "right"]].set_visible(False)

        ax_right.barh(range(len(a_alleles)), [marg_a[a] for a in a_alleles], color="#555555",
                     height=0.7)
        ax_right.set_xlabel("carrier\nfreq.", fontsize=5)
        ax_right.tick_params(labelleft=False, left=False, labelsize=4)
        ax_right.spines[["top", "right"]].set_visible(False)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cbar = fig.colorbar(sm, ax=[ax_heat, ax_right], location="right", fraction=0.05,
                            pad=0.14, shrink=0.6)
        cbar.set_label("Signed phased D' (2-field)", fontsize=6)

        fig.suptitle("DQA1-DQB1 signed phased D' -- from committed 29_hla_ld aggregates%s"
                     % (" (%s)" % title_suffix if title_suffix else ""), fontsize=8)
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
    """Cross-group (predicted_incompatible) allele pairs with status == estimated (i.e. clearing
    the 20-haplotype floor -- Cole's 'light blue' cells) in >=1 ancestry, tidy one-row-per-
    (pair, ancestry)."""
    m37mod = m37()
    classify_pair = m37mod.classify_pair
    df = df_all.copy()
    df["group"] = [classify_pair(a, b) for a, b in zip(df["allele_a"], df["allele_b"])]
    hits = df[(df["group"] == "predicted_incompatible") & (df["status"] == "estimated")]
    cols = ["allele_a", "allele_b", "ancestry", "n_hap_ij", "freq_a", "freq_b", "signed_Dprime"]
    return hits[cols].sort_values(["allele_a", "allele_b", "ancestry"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Bimodality: signed D' within compatible vs incompatible quadrants, pooled over ancestries.
# ---------------------------------------------------------------------------
def fig_bimodality(df_all, path_stem, ancestries=SUPPLEMENT_ANCESTRIES):
    import matplotlib.pyplot as plt

    m37mod = m37()
    classify_pair = m37mod.classify_pair

    df = df_all[df_all["ancestry"].isin(ancestries) & (df_all["status"] == "estimated")].copy()
    df["group"] = [classify_pair(a, b) for a, b in zip(df["allele_a"], df["allele_b"])]
    compatible = df.loc[df["group"].isin(["G1", "G2"]), "signed_Dprime"].dropna().to_numpy()
    incompatible = df.loc[df["group"] == "predicted_incompatible", "signed_Dprime"].dropna().to_numpy()

    with vc.nature_style():
        fig, ax = plt.subplots(figsize=(vc.mm(vc.NATURE_SINGLE_COL_MM), vc.mm(70)))
        bins = np.linspace(-1, 1, 41)
        ax.hist(compatible, bins=bins, alpha=0.7, density=True,
               label="compatible (G1/G1, G2/G2)\nn=%d" % len(compatible), color="#B2182B")
        if len(incompatible):
            ax.hist(incompatible, bins=bins, alpha=0.7, density=True,
                   label="predicted incompatible\nn=%d" % len(incompatible), color="#2166AC")
        else:
            ax.text(0.98, 0.92, "predicted incompatible: n=0 cells clear\nthe 20-haplotype floor "
                    "in any ancestry\n(complete purge at 2-field)", transform=ax.transAxes,
                    ha="right", va="top", fontsize=5, color="#2166AC")
        ax.set_xlabel("Signed phased D' (2-field, estimated cells only)")
        ax.set_ylabel("density")
        ax.legend(fontsize=5, frameon=False, loc="upper left")
        ax.set_title("Bimodality within the compatible quadrants\n(pooled %s)"
                     % ", ".join(ancestries), fontsize=7)
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
    print("[37c] main panel ancestry (most haplotypes): %s (N=%d)"
         % (main_anc, n_by_anc[main_anc]), flush=True)

    fig_g1g2_from_table(df_all[df_all["ancestry"] == main_anc],
                        os.path.join(args.out_dir, "fig_dq_g1g2_committed_MAIN_%s" % main_anc),
                        title_suffix="%s, main panel, N=%d haplotypes" % (main_anc, n_by_anc[main_anc]))

    for anc in SUPPLEMENT_ANCESTRIES:
        fig_g1g2_from_table(df_all[df_all["ancestry"] == anc],
                            os.path.join(args.out_dir, "fig_dq_g1g2_committed_supp_%s" % anc),
                            title_suffix="%s, N=%d haplotypes" % (anc, n_by_anc.get(anc) or 0))

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
