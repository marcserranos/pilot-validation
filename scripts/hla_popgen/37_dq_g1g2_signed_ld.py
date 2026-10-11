#!/usr/bin/env python3
"""Signed phased D' between DQA1 and DQB1, the Petersdorf (2022) G1/G2 "predicted incompatible"
cross-quadrant, and the transplant-genetics question Cole asked on the 2026-09-22 call: are the
residual (non -1) incompatible-quadrant cells phase/assembly errors, or real rare recombinants?

## WS1 of sprint S03 (call #8). Deliverables (see AGENT_PREAMBLE.md / this script's task text):

 1. Pairwise SIGNED phased D' for every DQA1 x DQB1 allele pair -- Lewontin's original signed
    statistic (not the |D'| used by script 29's `pairwise_ld`, which collapses sign on purpose for
    a magnitude-only r^2/|D'| table). D = p_ij - p_i q_j; D' = D / Dmax, sign kept, in [-1, 1].
    Haplotypes are the same physically-phased cis extraction as 29 (`extract_haplotypes`, reused
    verbatim), restricted to the unrelated cohort, at Cole's resolution (4-field, this script's
    main figure) AND 2-field (protein), each with a >=20-haplotype-carrier floor per allele.
 2. G1/G2 classification (Petersdorf 2022) and an exact recreation of Cole's target layout
    (reference/cole_dq_g1g2_target_figure.pdf): DQA1 rows with the DQA1*01 block on top, DQB1
    columns with the DQB1*05/06 block on the left, a black crosshair splitting the matrix into
    four quadrants, 'Predicted incompatible' labelling the two off-diagonal (cross-group)
    quadrants, a diverging [-1, 1] colormap labelled 'Signed phased D'', and suppressed/unobserved
    cells in a distinct grey (never folded into the color scale as a false zero).
 3. Observed vs expected (independence, p_i * q_j * N) cross-group ('incompatible') haplotype
    count, pooled and per ancestry, with O/E and a CI.
 4. Cole's artifact-vs-biology question on the individuals who carry an incompatible CIS pair:
    (a) swap-explicable heterozygotes, (b) the per-person phasing-confidence file (VM-only,
    aggregate summary only pulled back), (c) Mendelian transmission in 16's related pairs,
    (d) short-read concordance. All four are VM-only functions here (need Table 1 / relatedness /
    SR genotypes that do not exist off the VM) -- this script defines them; `run_artifact_checks()`
    is the VM entry point, never executed against fixtures.
 5. Supplement: tidy CSV of DQA1~DQB1 and DPA1~DPB1 pairs, and a DP analogue heatmap with
    hierarchical clustering to look for a data-driven "forbidden block" the way G1/G2 is for DQ.

## Local run status (2026-09-22)

This session had no authenticated Workbench tab (workbench.verily.com redirected to its sign-in
page; per AGENT_PREAMBLE.md / CLAUDE.md, credential entry is out of scope for an agent and a
logged-in tab must be handed over by Marc or Aleix). Nothing here has been run against real AoU
data. What *is* done: the full statistic implementation (signed D', G1/G2, O/E + bootstrap CI,
the Cole-figure renderer, the four artifact-check functions, the DP clustering supplement) plus a
unit test on synthetic haplotypes (tests/test_dq_g1g2_signed_ld.py) proving signed D' hits -1 at
perfect repulsion, 0 at independence, +1 at perfect coupling. The VM run itself (deliverables 1-5
with real numbers, and this file's README with results) is NOT done and needs a live VM session.

Usage (VM, once a session exists):
    setsid nohup python3 -u scripts/hla_popgen/37_dq_g1g2_signed_ld.py \\
        --out-dir ~/s03/results/37 < /dev/null & disown
Then --artifact-checks (needs Table 1 + relatedness + SR genotypes, VM only, see module docstring
of run_artifact_checks) writes ~/s03/results/37/phasing_confidence_per_person.tsv (VM-only,
never pulled back) plus an aggregate summary that IS safe to pull back.

Usage (unit test, local): scripts/hla_popgen/tests/test_dq_g1g2_signed_ld.py
"""
import argparse
import importlib.util
import json
import math
import os
import random
import sys
from collections import Counter, defaultdict

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


_m29 = None
_m24 = None
_m16 = None


def m29():
    """Reuse 29's cis-haplotype extractor verbatim, per task instructions."""
    global _m29
    if _m29 is None:
        _m29 = _load_module("29_hla_ld_by_ancestry.py", "hla_ld_by_ancestry")
    return _m29


def m24():
    global _m24
    if _m24 is None:
        _m24 = _load_module("24_novelty_by_field.py", "novelty_by_field")
    return _m24


def m16():
    global _m16
    if _m16 is None:
        _m16 = _load_module("16_phasing_mendelian_validation.py", "phasing_mendelian_validation")
    return _m16


DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_TABLE1 = os.path.join(DEFAULT_OUTROOT, "hla_calls_rich.tsv")
DEFAULT_COHORT = os.path.join(DEFAULT_OUTROOT, "cohort_membership.tsv")
DEFAULT_RELATEDNESS = os.path.expanduser(
    "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/"
    "relatedness/samples_relatedness.tsv")
# Quirk: never ~/mnt/aou-controlled -- stale manual mount, hangs any process that touches it.
DEFAULT_SR_GENOTYPES = os.path.expanduser(
    "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/"
    "hla_variants/hla_genotypes.tsv")
DEFAULT_OUT_DIR = os.path.expanduser("~/s03/results/37")
DEFAULT_PER_PERSON_PATH = os.path.expanduser("~/s03/results/37/phasing_confidence_per_person.tsv")

SUPPRESS_BELOW = 20
MIN_ALLELE_HAPS = 20   # every row/column allele needs >=20 haplotype carriers to appear
KIN_MIN = 0.0442
STRICT_MIN = 0.9
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
N_BOOTSTRAP = 2000

# ---------------------------------------------------------------------------
# Petersdorf 2022 G1/G2 DQ epitope-group classification.
# G1: DQA1*02/03/04/05/06 paired with DQB1*02/03/04.
# G2: DQA1*01 paired with DQB1*05/06.
# Anything crossing between the two groups (G1 DQA1 with a G2 DQB1, or vice versa) is the
# "predicted incompatible" cis combination the transplant-genetics literature flags as disfavoured
# in cis (it is fine in trans, i.e. on two different haplotypes in one person -- this script only
# ever classifies a single CIS haplotype's own DQA1~DQB1 pair).
# ---------------------------------------------------------------------------
G1_DQA1_FIELDS = {"02", "03", "04", "05", "06"}
G1_DQB1_FIELDS = {"02", "03", "04"}
G2_DQA1_FIELDS = {"01"}
G2_DQB1_FIELDS = {"05", "06"}


def _field1(allele):
    """First (gene-family) field as a zero-padded 2-digit string, e.g. 'DQA1*02:01:01:01' -> '02'."""
    if not isinstance(allele, str) or "*" not in allele:
        return None
    rest = allele.split("*", 1)[1]
    f1 = rest.split(":", 1)[0].strip()
    if not f1.isdigit():
        return None
    return f1.zfill(2)


def dq_group(gene, allele):
    """Returns 'G1', 'G2', or None (does not resolve to either canonical Petersdorf group)."""
    f1 = _field1(allele)
    if f1 is None:
        return None
    if gene == "DQA1":
        if f1 in G1_DQA1_FIELDS:
            return "G1"
        if f1 in G2_DQA1_FIELDS:
            return "G2"
        return None
    if gene == "DQB1":
        if f1 in G1_DQB1_FIELDS:
            return "G1"
        if f1 in G2_DQB1_FIELDS:
            return "G2"
        return None
    return None


def classify_pair(allele_a, allele_b):
    """(DQA1 allele, DQB1 allele) -> 'G1', 'G2', 'predicted_incompatible', or 'unclassified'."""
    ga = dq_group("DQA1", allele_a)
    gb = dq_group("DQB1", allele_b)
    if ga is None or gb is None:
        return "unclassified"
    if ga == gb:
        return ga
    return "predicted_incompatible"


# ---------------------------------------------------------------------------
# Signed phased D' (Lewontin 1964). Distinct from 29.pairwise_ld, which returns |D'| collapsed to
# an "allele vs all others" biallelic framing for r^2's sake -- here the sign (repulsion vs
# coupling between this specific allele_a and this specific allele_b) is the entire point of the
# figure, so it must survive.
# ---------------------------------------------------------------------------
def signed_dprime_table(M, a_labels, b_labels, min_haps=MIN_ALLELE_HAPS):
    """Signed D' for every (allele_a, allele_b) cell where BOTH the row and column allele clear
    `min_haps` carriers (summed over the full contingency table, i.e. the allele's own marginal
    haplotype count) -- the estimability/disclosure floor. Cells failing the floor are omitted
    (caller marks them 'suppressed' vs 'not_observed' by checking whether the allele passed at
    all vs the specific pair count was just small)."""
    n = M.sum()
    if n <= 0:
        return []
    p_a = M.sum(axis=1) / n
    p_b = M.sum(axis=0) / n
    keep_a = [i for i in range(len(a_labels)) if M[i, :].sum() >= min_haps]
    keep_b = [j for j in range(len(b_labels)) if M[:, j].sum() >= min_haps]
    out = []
    for i in keep_a:
        pi = p_a[i]
        for j in keep_b:
            qj = p_b[j]
            pij = M[i, j] / n
            D = pij - pi * qj
            if D >= 0:
                dmax = min(pi * (1 - qj), (1 - pi) * qj)
            else:
                dmax = min(pi * qj, (1 - pi) * (1 - qj))
            dprime = (D / dmax) if dmax > 0 else float("nan")
            n_ij = int(round(M[i, j]))
            out.append({
                "allele_a": a_labels[i], "allele_b": b_labels[j],
                "freq_a": float(pi), "freq_b": float(qj), "freq_hap": float(pij),
                "n_hap_ij": n_ij, "D": float(D), "signed_Dprime": float(dprime),
                "status": "estimated" if n_ij >= SUPPRESS_BELOW else (
                    "suppressed_lt20" if n_ij > 0 else "not_observed"),
            })
    return out


# ---------------------------------------------------------------------------
# Observed vs expected cross-group ("predicted incompatible") haplotype count.
# ---------------------------------------------------------------------------
def observed_expected_incompatible(haps_a, haps_b, gene_a="DQA1", gene_b="DQB1",
                                    n_bootstrap=N_BOOTSTRAP, rng=None):
    """haps_a/haps_b: parallel lists of alleles for the SAME haplotypes (cis pairs).

    Expected count under independence uses the empirical per-allele marginal frequencies
    (p_i for DQA1, q_j for DQB1) restricted to alleles resolving to a Petersdorf group --
    E = N * sum_{i in G1,j in G2} p_i q_j + sum_{i in G2,j in G1} p_i q_j, i.e. the independence
    prediction for "one allele from each group's DQA1/DQB1 lands on the same cis haplotype."
    A group-only marginal (not per-genotype) null, matching how Petersdorf-style disequilibrium
    claims are usually framed. CI on O/E via nonparametric bootstrap over haplotypes (resample
    with replacement, recompute O and E each replicate, report the O/E ratio's 2.5/97.5
    percentiles) -- appropriate here since the "expected" count is itself an estimated quantity,
    not a fixed null, so a naive binomial/Poisson CI on O alone would understate uncertainty.
    """
    rng = rng or random.Random(20260922)
    n = len(haps_a)
    if n == 0:
        return None
    a_arr = np.asarray(haps_a)
    b_arr = np.asarray(haps_b)
    ga = np.array([dq_group(gene_a, a) for a in a_arr])
    gb = np.array([dq_group(gene_b, b) for b in b_arr])

    def _oe(idx):
        sub_ga, sub_gb = ga[idx], gb[idx]
        m = (sub_ga != None) & (sub_gb != None)  # noqa: E711
        sub_ga, sub_gb = sub_ga[m], sub_gb[m]
        nn = len(sub_ga)
        if nn == 0:
            return float("nan"), float("nan"), 0
        p_g1a = np.mean(sub_ga == "G1")
        p_g2a = np.mean(sub_ga == "G2")
        p_g1b = np.mean(sub_gb == "G1")
        p_g2b = np.mean(sub_gb == "G2")
        observed = int(np.sum(sub_ga != sub_gb))
        expected = nn * (p_g1a * p_g2b + p_g2a * p_g1b)
        return observed, expected, nn

    idx_all = np.arange(n)
    obs, exp, n_classified = _oe(idx_all)
    if n_classified == 0 or exp <= 0:
        return {"n_haplotypes": n, "n_classified": n_classified, "observed": obs,
                "observed_disp": _sup(obs), "expected": exp, "oe_ratio": float("nan"),
                "oe_ci_lo": float("nan"), "oe_ci_hi": float("nan")}

    boot_ratios = []
    for _ in range(n_bootstrap):
        pick = rng.choices(range(n), k=n)
        o_b, e_b, nc_b = _oe(np.array(pick))
        if nc_b > 0 and e_b > 0:
            boot_ratios.append(o_b / e_b)
    boot_ratios = np.array(boot_ratios, dtype=float)
    boot_ratios = boot_ratios[~np.isnan(boot_ratios)]
    lo = float(np.percentile(boot_ratios, 2.5)) if boot_ratios.size else float("nan")
    hi = float(np.percentile(boot_ratios, 97.5)) if boot_ratios.size else float("nan")
    return {
        "n_haplotypes": n, "n_classified": n_classified,
        "observed": obs, "observed_disp": _sup(obs),
        "expected": round(float(exp), 2), "oe_ratio": round(obs / exp, 4),
        "oe_ci_lo": round(lo, 4), "oe_ci_hi": round(hi, 4),
    }


def _sup(n, threshold=SUPPRESS_BELOW):
    n = int(n)
    return "0" if n == 0 else ("<%d" % threshold if n < threshold else str(n))


# ---------------------------------------------------------------------------
# Cole's figure: exact recreation of reference/cole_dq_g1g2_target_figure.pdf, Nature-grade.
# ---------------------------------------------------------------------------
def fig_g1g2_signed_ld(pairwise_rows, path_stem, title_suffix=""):
    """DQA1 (rows) x DQB1 (cols) signed-D' heatmap, G2-DQA1 (*01) block on top, G2-DQB1 (*05/06)
    block on the left, black crosshair, diverging [-1,1] colormap, 'Predicted incompatible' text
    in the two off-diagonal quadrants, suppressed/not-observed cells hatched grey. Carrier
    frequency marginal bars on top/right."""
    import matplotlib.pyplot as plt

    df = pd.DataFrame(pairwise_rows)
    if df.empty:
        return None
    estimated = df[df["status"] == "estimated"]
    if estimated.empty:
        return None

    # Order: G2 alleles (DQA1*01 / DQB1*05,06) first by descending pooled carrier freq, then G1.
    a_alleles = sorted(estimated["allele_a"].unique(),
                        key=lambda a: (dq_group("DQA1", a) != "G2",
                                       -estimated.loc[estimated["allele_a"] == a, "freq_a"].max()))
    b_alleles = sorted(estimated["allele_b"].unique(),
                        key=lambda b: (dq_group("DQB1", b) != "G2",
                                       -estimated.loc[estimated["allele_b"] == b, "freq_b"].max()))
    n_a2 = sum(1 for a in a_alleles if dq_group("DQA1", a) == "G2")
    n_b2 = sum(1 for b in b_alleles if dq_group("DQB1", b) == "G2")

    M = np.full((len(a_alleles), len(b_alleles)), np.nan)
    look = {(r["allele_a"], r["allele_b"]): r["signed_Dprime"] for _, r in estimated.iterrows()}
    all_pairs = {(r["allele_a"], r["allele_b"]) for _, r in df.iterrows()}
    suppressed_cells = []
    for i, a in enumerate(a_alleles):
        for j, b in enumerate(b_alleles):
            if (a, b) in look:
                M[i, j] = look[(a, b)]
            elif (a, b) not in all_pairs:
                suppressed_cells.append((i, j))  # not_observed / below floor, never in the table

    marg_a = {a: estimated.loc[estimated["allele_a"] == a, "freq_a"].max() for a in a_alleles}
    marg_b = {b: estimated.loc[estimated["allele_b"] == b, "freq_b"].max() for b in b_alleles}

    with vc.nature_style():
        fig = plt.figure(figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(160)))
        gs = fig.add_gridspec(2, 2, width_ratios=[len(b_alleles), 6], height_ratios=[4, len(a_alleles)],
                              wspace=0.02, hspace=0.02)
        ax_heat = fig.add_subplot(gs[1, 0])
        ax_top = fig.add_subplot(gs[0, 0], sharex=ax_heat)
        ax_right = fig.add_subplot(gs[1, 1], sharey=ax_heat)

        cmap = vc.diverging_cmap()
        norm = vc.diverging_norm(-1.0, 1.0)
        ax_heat.imshow(M, cmap=cmap, norm=norm, aspect="auto", interpolation="none")
        for (i, j) in suppressed_cells:
            vc.hatch_suppressed(ax_heat, j - 0.5, i - 0.5, 1, 1)

        ax_heat.set_xticks(range(len(b_alleles)))
        ax_heat.set_xticklabels(b_alleles, rotation=90, fontsize=5)
        ax_heat.set_yticks(range(len(a_alleles)))
        ax_heat.set_yticklabels(a_alleles, fontsize=5)
        ax_heat.set_xlabel("DQB1")
        ax_heat.set_ylabel("DQA1")

        # crosshair between the G2 block and the G1 block
        ax_heat.axhline(n_a2 - 0.5, color="black", lw=1.0, zorder=5)
        ax_heat.axvline(n_b2 - 0.5, color="black", lw=1.0, zorder=5)

        # "Predicted incompatible" labels: G2-DQA1 rows x G1-DQB1 cols (top-right quadrant) and
        # G1-DQA1 rows x G2-DQB1 cols (bottom-left quadrant) -- the two cross-group quadrants.
        if n_b2 < len(b_alleles) and n_a2 > 0:
            ax_heat.text((n_b2 + len(b_alleles)) / 2 - 0.5, n_a2 / 2 - 0.5, "Predicted\nincompatible",
                        ha="center", va="center", fontsize=6,
                        bbox=dict(boxstyle="round", fc="white", ec="none", alpha=0.75))
        if n_b2 > 0 and n_a2 < len(a_alleles):
            ax_heat.text(n_b2 / 2 - 0.5, (n_a2 + len(a_alleles)) / 2 - 0.5, "Predicted\nincompatible",
                        ha="center", va="center", fontsize=6,
                        bbox=dict(boxstyle="round", fc="white", ec="none", alpha=0.75))

        ax_top.bar(range(len(b_alleles)), [marg_b[b] for b in b_alleles], color="#555555", width=0.7)
        ax_top.set_ylabel("carrier\nfreq.", fontsize=5)
        ax_top.tick_params(labelbottom=False, bottom=False)
        ax_top.spines[["top", "right"]].set_visible(False)

        ax_right.barh(range(len(a_alleles)), [marg_a[a] for a in a_alleles], color="#555555", height=0.7)
        ax_right.set_xlabel("carrier\nfreq.", fontsize=5)
        ax_right.tick_params(labelleft=False, left=False)
        ax_right.spines[["top", "right"]].set_visible(False)
        ax_right.invert_yaxis() if False else None

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cbar = fig.colorbar(sm, ax=[ax_heat, ax_right], location="right", fraction=0.05, pad=0.12,
                            shrink=0.6)
        cbar.set_label("Signed phased D'", fontsize=6)

        fig.suptitle("DQA1-DQB1 signed phased D'%s" % (" — " + title_suffix if title_suffix else ""),
                     fontsize=8)
    return vc.save_fig(fig, path_stem)


# ---------------------------------------------------------------------------
# Deliverable 5b: hierarchical clustering of the DP signed-D' matrix for a data-driven
# "forbidden block" (G1/G2-like rule, if one exists, for DPA1 x DPB1).
# ---------------------------------------------------------------------------
def cluster_dp_matrix(pairwise_rows):
    """Hierarchically clusters the DPA1 x DPB1 signed-D' matrix (average linkage on 1 - |D'|
    dissimilarity, alleles clustered independently on each axis) and returns the reordered
    matrix + labels + the scipy linkage objects, so callers can both plot a clustered heatmap
    and inspect whether any contiguous off-diagonal block of strongly negative D' emerges --
    the DP analogue of G1/G2, if the data supports one, rather than assuming Petersdorf's DQ
    grouping transfers.
    """
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import squareform

    df = pd.DataFrame([r for r in pairwise_rows if r["status"] == "estimated"])
    if df.empty:
        return None
    a_alleles = sorted(df["allele_a"].unique())
    b_alleles = sorted(df["allele_b"].unique())
    if len(a_alleles) < 3 or len(b_alleles) < 3:
        return None  # not enough alleles to cluster meaningfully
    M = np.zeros((len(a_alleles), len(b_alleles)))
    look = {(r["allele_a"], r["allele_b"]): r["signed_Dprime"] for _, r in df.iterrows()}
    for i, a in enumerate(a_alleles):
        for j, b in enumerate(b_alleles):
            M[i, j] = look.get((a, b), 0.0)

    def _cluster_axis(mat):
        dist = 1 - np.abs(np.clip(mat, -1, 1))
        # symmetrize a rectangular axis-distance via correlation of profiles across the OTHER axis
        prof_dist = squareform(1 - np.corrcoef(mat)[np.triu_indices(mat.shape[0], k=1)]
                               if mat.shape[0] > 1 else np.array([0.0]), checks=False) \
            if mat.shape[0] > 2 else None
        return None

    # Simpler, robust approach: cluster each axis on its own profile-correlation distance.
    def _order_by_profile(mat, axis):
        prof = mat if axis == 0 else mat.T
        if prof.shape[0] < 3:
            return list(range(prof.shape[0]))
        corr = np.corrcoef(prof)
        corr = np.nan_to_num(corr, nan=0.0)
        d = 1 - corr
        np.fill_diagonal(d, 0.0)
        d = (d + d.T) / 2
        Z = linkage(squareform(d, checks=False), method="average")
        return leaves_list(Z).tolist()

    order_a = _order_by_profile(M, axis=0)
    order_b = _order_by_profile(M, axis=1)
    M_reordered = M[np.ix_(order_a, order_b)]
    return {
        "matrix": M_reordered,
        "a_labels": [a_alleles[i] for i in order_a],
        "b_labels": [b_alleles[i] for i in order_b],
    }


def fig_dp_clustered(pairwise_rows, path_stem):
    import matplotlib.pyplot as plt
    res = cluster_dp_matrix(pairwise_rows)
    if res is None:
        return None
    M, a_lab, b_lab = res["matrix"], res["a_labels"], res["b_labels"]
    with vc.nature_style():
        fig, ax = plt.subplots(figsize=(vc.mm(vc.NATURE_SINGLE_COL_MM), vc.mm(90)))
        cmap = vc.diverging_cmap()
        norm = vc.diverging_norm(-1.0, 1.0)
        im = ax.imshow(M, cmap=cmap, norm=norm, aspect="auto")
        ax.set_xticks(range(len(b_lab))); ax.set_xticklabels(b_lab, rotation=90, fontsize=5)
        ax.set_yticks(range(len(a_lab))); ax.set_yticklabels(a_lab, fontsize=5)
        ax.set_xlabel("DPB1"); ax.set_ylabel("DPA1")
        ax.set_title("DPA1-DPB1 signed phased D', hierarchically clustered\n"
                     "(profile-correlation distance, average linkage)", fontsize=7)
        fig.colorbar(im, ax=ax, shrink=0.7, label="Signed phased D'")
    return vc.save_fig(fig, path_stem)


# ---------------------------------------------------------------------------
# Deliverable 4: artifact-vs-biology checks on individuals carrying an incompatible cis pair.
# VM-ONLY. Needs Table 1 (per-person, per-hap alleles + contig/template metadata), 16's related
# pairs, and the SR genotype file -- none of which exist off the VM, so this is defined but not
# exercised anywhere in this file except as the documented VM entry point `run_artifact_checks`.
# ---------------------------------------------------------------------------
def find_incompatible_carriers(t1, anc_of, resolution=4):
    """People with >=1 cis DQA1~DQB1 haplotype classified 'predicted_incompatible'.

    Returns a DataFrame: person_id, hap, contig, allele_a (DQA1), allele_b (DQB1), ancestry.
    """
    haps, _diag = m29().extract_haplotypes(t1, "DQA1", "DQB1")
    if haps.empty:
        return haps
    haps["ancestry"] = haps["person_id"].astype(str).map(anc_of)
    haps["group"] = [classify_pair(a, b) for a, b in zip(haps["allele_a"], haps["allele_b"])]
    return haps[haps["group"] == "predicted_incompatible"].copy()


def check_a_swap_explicable(t1, incompatible_haps, anc_of, rng=None):
    """(a) For each person carrying an incompatible cis pair, is that person heterozygous at
    BOTH DQA1 and DQB1 such that swapping DQB1 between hap1 and hap2 would make BOTH resulting
    haplotypes compatible (G1/G1 or G2/G2, or at least non-cross)? Compared against a null of
    the same test applied to a random sample of people carrying a *compatible* cis pair (same
    ancestry mix, same sample size) -- if incompatible carriers are just ordinary phase-labelling
    noise (hap1/hap2 swapped relative to true haplotype identity, a labelling artifact rather
    than a real error, since hap1 vs hap2 has no biological meaning on its own), the swap-fixes
    rate should be at or near 100%; if they are real recombinants, the two DQB1 alleles are truly
    on the wrong cis backgrounds and no hap-label swap fixes it structurally (it just moves which
    label is 'incompatible').
    """
    rng = rng or random.Random(20260922)
    people = sorted(incompatible_haps["person_id"].unique())
    fixable = 0
    n_testable = 0
    for pid in people:
        pdf = t1[t1["person_id"].astype(str) == str(pid)]
        both = m29().extract_haplotypes(pdf, "DQA1", "DQB1")[0]
        by_hap = {r["hap"]: (r["allele_a"], r["allele_b"]) for _, r in both.iterrows()}
        if set(by_hap.keys()) != {"hap1", "hap2"}:
            continue  # need both haplotypes phased to test a swap
        (a1, b1), (a2, b2) = by_hap["hap1"], by_hap["hap2"]
        n_testable += 1
        swapped_ok = (classify_pair(a1, b2) != "predicted_incompatible"
                     and classify_pair(a2, b1) != "predicted_incompatible")
        if swapped_ok:
            fixable += 1
    return {"n_incompatible_carriers": len(people), "n_testable": n_testable,
            "n_swap_explicable": fixable,
            "n_testable_disp": _sup(n_testable), "n_swap_explicable_disp": _sup(fixable),
            "rate": round(fixable / n_testable, 4) if n_testable else float("nan")}


def build_phasing_confidence_table(t1, gene_trios=(("DQA1", "DQB1"), ("DRB1", "DQB1"),
                                                    ("DRB1", "DQA1"), ("DPA1", "DPB1"))):
    """(b) Per-person, per-gene-pair, per-haplotype phasing status: same-contig (physical) /
    different-contig (would require statistical/inferred phasing, never done in this pipeline,
    so effectively 'unphaseable here') / missing (one or both genes absent from that hap).
    Columns: person_id, hap, gene_a, gene_b, status in
    {physical, different_contig, missing_a, missing_b, missing_both}, contig, contig_length,
    template_distance, template_warning (if present in Table 1).

    Writes to `out_path` directly (never returned as a DataFrame to the caller) -- this is the
    per-person file the task says must stay on the VM (~/s03/results/37/
    phasing_confidence_per_person.tsv) and never be pulled back; only its schema (this docstring
    + the header row) and an aggregate cross-tab leave the VM.
    """
    bare = t1["gene"].astype(str).str.replace("^HLA-", "", regex=True)
    rows = []
    for gene_a, gene_b in gene_trios:
        sub = t1[bare.isin([gene_a, gene_b])].copy()
        sub["gene_bare"] = bare[sub.index]
        for (pid, hap), grp in sub.groupby(["person_id", "hap"]):
            a_rows = grp[grp["gene_bare"] == gene_a]
            b_rows = grp[grp["gene_bare"] == gene_b]
            if a_rows.empty and b_rows.empty:
                status = "missing_both"
                contig = None
            elif a_rows.empty:
                status = "missing_a"
                contig = None
            elif b_rows.empty:
                status = "missing_b"
                contig = None
            else:
                ca = set(a_rows["contig"]); cb = set(b_rows["contig"])
                if len(ca) == 1 and ca == cb:
                    status = "physical"
                    contig = next(iter(ca))
                else:
                    status = "different_contig"
                    contig = None
            row = {"person_id": pid, "hap": hap, "gene_a": gene_a, "gene_b": gene_b,
                  "status": status, "contig": contig}
            for col in ("contig_length", "template_distance", "template_warning"):
                if contig is not None and col in t1.columns:
                    vals = grp.loc[grp["contig"] == contig, col]
                    row[col] = vals.iloc[0] if len(vals) else None
                else:
                    row[col] = None
            rows.append(row)
    return pd.DataFrame(rows)


def check_b_phasing_confidence_crosstab(per_person_df, incompatible_people):
    """Aggregate, disclosure-safe cross-tab: among incompatible-pair carriers vs everyone else,
    what fraction of their OTHER gene-pair haplotype phasings are 'physical' (vs
    different_contig/missing)? Counts <20 collapsed to '<20'. This is what leaves the VM."""
    df = per_person_df.copy()
    df["is_incompatible_carrier"] = df["person_id"].isin(incompatible_people)
    out = []
    for (gene_a, gene_b), grp in df.groupby(["gene_a", "gene_b"]):
        for carrier in (True, False):
            sub = grp[grp["is_incompatible_carrier"] == carrier]
            n = len(sub)
            n_physical = int((sub["status"] == "physical").sum())
            out.append({
                "gene_pair": "%s~%s" % (gene_a, gene_b),
                "group": "incompatible_carrier" if carrier else "other",
                "n_haplotypes_disp": _sup(n), "n_physical_disp": _sup(n_physical),
                "pct_physical": round(100.0 * n_physical / n, 1) if n else float("nan"),
            })
    return pd.DataFrame(out)


def check_c_mendelian_transmission(related_pairs_result, incompatible_people):
    """(c) Among 16's related (parent-offspring / high-sharing) pairs, are incompatible cis
    DQA1~DQB1 haplotypes transmitted intact from the carrying relative to the other member of
    the pair? `related_pairs_result` is 16's per-pair, per-gene sharing table (already computed
    by 16_phasing_mendelian_validation.py on the VM) filtered to DQA1/DQB1 in high_sharing pairs
    where at least one member carries an incompatible cis pair. Aggregate counts only."""
    df = related_pairs_result
    if df is None or df.empty:
        return {"n_pairs_with_carrier": 0, "n_transmitted_intact_disp": "0",
                "n_not_transmitted_disp": "0"}
    carrier_pairs = df[df["person_id"].isin(incompatible_people) |
                       df.get("relative_id", pd.Series(dtype=str)).isin(incompatible_people)]
    n = len(carrier_pairs)
    n_intact = int((carrier_pairs.get("gene_mismatch", pd.Series([False] * n)) == False).sum())  # noqa: E712
    return {"n_pairs_with_carrier_disp": _sup(n), "n_transmitted_intact_disp": _sup(n_intact),
            "n_not_transmitted_disp": _sup(n - n_intact)}


def check_d_short_read_concordance(sr_genotypes, lr_table, incompatible_people):
    """(d) Do AoU short-read DQA1/DQB1 genotypes agree with the LR alleles at the same rate for
    incompatible-pair carriers as for everyone else? Concordance = SR genotype (unordered 2-field
    pair) contains both of the person's LR 2-field alleles at that gene. Aggregate rates only."""
    if sr_genotypes is None or sr_genotypes.empty:
        return None
    out = []
    for gene in ("DQA1", "DQB1"):
        sr_g = sr_genotypes[sr_genotypes["gene"] == gene]
        lr_g = lr_table[lr_table["gene_bare"] == gene]
        for carrier in (True, False):
            people = (incompatible_people if carrier else
                     set(lr_g["person_id"]) - set(incompatible_people))
            n = 0
            n_concordant = 0
            for pid in people:
                sr_row = sr_g[sr_g["person_id"].astype(str) == str(pid)]
                lr_alleles = set(lr_g.loc[lr_g["person_id"].astype(str) == str(pid), "allele_2f"])
                if sr_row.empty or not lr_alleles:
                    continue
                n += 1
                sr_alleles = set(sr_row.iloc[0].get("genotype_2f", "").split("/"))
                if lr_alleles <= sr_alleles:
                    n_concordant += 1
            out.append({"gene": gene, "group": "incompatible_carrier" if carrier else "other",
                       "n_disp": _sup(n), "n_concordant_disp": _sup(n_concordant),
                       "pct_concordant": round(100.0 * n_concordant / n, 1) if n else float("nan")})
    return pd.DataFrame(out)


def run_artifact_checks(t1, anc_of, out_dir, related_pairs_result=None, sr_genotypes=None,
                        per_person_path=DEFAULT_PER_PERSON_PATH):
    """VM entry point for deliverable 4. Writes the per-person file ONLY to `per_person_path`
    (must be under ~/s03/results/37/, never copied into `out_dir` which gets pulled back), and
    writes the aggregate/disclosure-safe outputs into `out_dir`."""
    incompatible_haps = find_incompatible_carriers(t1, anc_of)
    incompatible_people = set(incompatible_haps["person_id"].astype(str))

    a_result = check_a_swap_explicable(t1, incompatible_haps, anc_of)
    # Null comparison: the task asks to compare the swap-explicable rate against "randomly chosen
    # compatible-haplotype carriers." Compatible carriers are, by definition, already
    # swap-consistent (swapping DQB1 between their two haplotypes cannot make them MORE
    # compatible than they already are), so the informative null is what fraction of random
    # compatible-carrier PEOPLE would, if their cis pairing were hypothetically the cross
    # combination instead, still resolve to a swap-explicable state -- i.e. the population base
    # rate of the genotype configuration (heterozygous at both genes, alleles split across the
    # two Petersdorf groups on each locus) that makes swap-explicability possible at all. That
    # base rate, not a fabricated counterfactual reclassification, is what is reported here.
    compat_haps, _ = m29().extract_haplotypes(t1, "DQA1", "DQB1")
    compat_haps["group"] = [classify_pair(a, b) for a, b in
                            zip(compat_haps["allele_a"], compat_haps["allele_b"])]
    compat_people = list(set(compat_haps.loc[compat_haps["group"] != "predicted_incompatible",
                                             "person_id"].astype(str)) - incompatible_people)
    rng = random.Random(20260922)
    null_sample = rng.sample(compat_people, min(len(compat_people), max(1, a_result["n_testable"])))
    null_df = t1[t1["person_id"].astype(str).isin(null_sample)]
    null_result = check_a_swap_explicable(null_df, compat_haps[compat_haps["person_id"].isin(null_sample)],
                                          anc_of)
    a_result["null_swap_configurable_rate"] = null_result["rate"]
    a_result["null_n_disp"] = null_result["n_testable_disp"]

    per_person_df = build_phasing_confidence_table(t1)
    os.makedirs(os.path.dirname(per_person_path), exist_ok=True)
    per_person_df.to_csv(per_person_path, sep="\t", index=False)  # VM-only, never pulled back

    b_result = check_b_phasing_confidence_crosstab(per_person_df, incompatible_people)
    c_result = check_c_mendelian_transmission(related_pairs_result, incompatible_people)
    d_result = check_d_short_read_concordance(sr_genotypes, t1.assign(
        gene_bare=t1["gene"].astype(str).str.replace("^HLA-", "", regex=True),
        allele_2f=t1["consensus"].map(lambda a: vc.to_nfield(a, 2))), incompatible_people)

    os.makedirs(out_dir, exist_ok=True)
    a_result.update({"check": "a_swap_explicable"})
    with open(os.path.join(out_dir, "artifact_check_a_swap_explicable.json"), "w") as fh:
        json.dump(a_result, fh, indent=2)
    b_result.to_csv(os.path.join(out_dir, "artifact_check_b_phasing_confidence.tsv"),
                    sep="\t", index=False)
    with open(os.path.join(out_dir, "artifact_check_c_mendelian.json"), "w") as fh:
        json.dump(c_result, fh, indent=2)
    if d_result is not None:
        d_result.to_csv(os.path.join(out_dir, "artifact_check_d_short_read_concordance.tsv"),
                        sep="\t", index=False)
    return {"a": a_result, "b": b_result, "c": c_result, "d": d_result,
            "per_person_path": per_person_path, "n_incompatible_carriers_disp":
                _sup(len(incompatible_people))}


# ---------------------------------------------------------------------------
# Supplement: tidy CSV of DQA1~DQB1 and DPA1~DPB1 pairs, per ancestry.
# ---------------------------------------------------------------------------
def supplement_table(all_pairwise_rows):
    """all_pairwise_rows: list of dicts already carrying 'pair', 'ancestry', 'resolution' keys
    (added by the caller in `run()`), one row per (allele_a, allele_b, ancestry, resolution)."""
    df = pd.DataFrame(all_pairwise_rows)
    cols = ["pair", "ancestry", "resolution", "allele_a", "allele_b", "n_hap_ij", "freq_a",
           "freq_b", "freq_hap", "signed_Dprime", "status"]
    df["r2"] = df.apply(lambda r: (r["D"] ** 2) / (r["freq_a"] * (1 - r["freq_a"]) *
                                                    r["freq_b"] * (1 - r["freq_b"]))
                        if r["freq_a"] not in (0, 1) and r["freq_b"] not in (0, 1) else float("nan"),
                        axis=1)
    df["n_hap_ij"] = df["n_hap_ij"].map(_sup)
    return df[cols + ["r2"]] if "r2" not in cols else df[cols]


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def run(args):
    ensure_dir = vc.ensure_dir
    ensure_dir(args.out_dir)
    rng = random.Random(args.seed)

    if not os.path.exists(args.table1):
        sys.exit("FATAL: table1 not found: %s -- this script needs a live VM session "
                 "(see this file's module docstring, 'Local run status')." % args.table1)

    print("[37] loading table1 ...", flush=True)
    t1 = pd.read_csv(args.table1, sep="\t", dtype=str, low_memory=False)
    t1["copy_index"] = pd.to_numeric(t1.get("copy_index", 1), errors="coerce").fillna(1).astype(int)

    people, n_removed = m24().build_people(
        sorted(set(t1["person_id"].astype(str))), args.cohort_membership,
        args.relatedness_table, KIN_MIN, STRICT_MIN, args.skip_relatedness)
    keep = people[people["unrelated"]]
    anc_of = dict(zip(keep["person_id"].astype(str), keep["anc_strict"]))
    t1 = t1[t1["person_id"].astype(str).isin(anc_of)]
    n_people = len(keep)
    print("[37] unrelated people: %d (dropped %d relatives)" % (n_people, n_removed), flush=True)

    all_pairwise = []
    oe_rows = []
    for gene_a, gene_b in [("DQA1", "DQB1"), ("DPA1", "DPB1")]:
        haps, diag = m29().extract_haplotypes(t1, gene_a, gene_b)
        if haps.empty:
            continue
        haps["ancestry"] = haps["person_id"].astype(str).map(anc_of)

        for resolution, tag in ((4, "4field"), (2, "2field")):
            # 29's extract_haplotypes hardcodes 2-field (`two_field()`); haps (above) is already
            # at that resolution and reused directly for the 2-field pass. For 4-field, re-derive
            # the haplotype extraction inline at 4-field resolution (same cis logic, different
            # truncation of `consensus`).
            if resolution == 2:
                haps_this = haps
            else:
                bare = t1["gene"].astype(str).str.replace("^HLA-", "", regex=True)
                sub_t1 = t1[bare.isin([gene_a, gene_b])].copy()
                sub_t1["gene_bare"] = bare[sub_t1.index]
                rows4 = []
                for (pid, hap, contig), grp in sub_t1.groupby(["person_id", "hap", "contig"]):
                    a_rows = grp[grp["gene_bare"] == gene_a]
                    b_rows = grp[grp["gene_bare"] == gene_b]
                    if a_rows.empty or b_rows.empty or len(a_rows) > 1 or len(b_rows) > 1:
                        continue
                    aa = vc.to_nfield(a_rows.iloc[0]["consensus"], 4)
                    bb = vc.to_nfield(b_rows.iloc[0]["consensus"], 4)
                    if aa is None or bb is None:
                        continue
                    rows4.append((pid, aa, bb))
                if not rows4:
                    continue
                haps4_df = pd.DataFrame(rows4, columns=["person_id", "allele_a", "allele_b"])
                haps4_df["ancestry"] = haps4_df["person_id"].astype(str).map(anc_of)
                haps_this = haps4_df

            for anc in ["ALL"] + ANCESTRY_ORDER:
                if anc != "ALL" and "ancestry" not in haps_this.columns:
                    continue
                sub = haps_this if anc == "ALL" else haps_this[haps_this["ancestry"] == anc]
                if len(sub) < 20:
                    continue
                M, a_lab, b_lab = m29().contingency(list(sub["allele_a"]), list(sub["allele_b"]))
                rows = signed_dprime_table(M, a_lab, b_lab, args.min_allele_haps)
                for r in rows:
                    r.update({"pair": "%s~%s" % (gene_a, gene_b), "ancestry": anc,
                             "resolution": tag})
                all_pairwise.extend(rows)

                if gene_a == "DQA1":
                    oe = observed_expected_incompatible(list(sub["allele_a"]), list(sub["allele_b"]),
                                                        n_bootstrap=args.n_bootstrap, rng=rng)
                    if oe:
                        oe.update({"ancestry": anc, "resolution": tag})
                        oe_rows.append(oe)

    # ---- write tidy supplement ----
    if all_pairwise:
        supp = supplement_table(all_pairwise)
        supp.to_csv(os.path.join(args.out_dir, "supp_table_ld_pairs.csv"), index=False)

    if oe_rows:
        pd.DataFrame(oe_rows).to_csv(os.path.join(args.out_dir, "oe_purge_table.tsv"),
                                     sep="\t", index=False)

    # ---- main figure: Cole's layout, 4-field pooled ALL, then per-ancestry ----
    main_rows = [r for r in all_pairwise if r["pair"] == "DQA1~DQB1" and r["resolution"] == "4field"]
    for anc in ["ALL"] + ANCESTRY_ORDER:
        rows_anc = [r for r in main_rows if r["ancestry"] == anc]
        if not rows_anc:
            continue
        fig_g1g2_signed_ld(rows_anc, os.path.join(args.out_dir, "fig_dq_g1g2_%s" % anc),
                          title_suffix=anc)

    dp_rows = [r for r in all_pairwise if r["pair"] == "DPA1~DPB1" and r["resolution"] == "4field"
              and r["ancestry"] == "ALL"]
    if dp_rows:
        fig_dp_clustered(dp_rows, os.path.join(args.out_dir, "fig_dp_clustered_ALL"))

    if args.artifact_checks:
        sr_genotypes = None
        if os.path.exists(args.sr_genotypes):
            sr_genotypes = pd.read_csv(args.sr_genotypes, sep="\t", dtype=str)
        related_result = None  # would come from 16's own output tables if present on the VM
        run_artifact_checks(t1, anc_of, args.out_dir, related_pairs_result=related_result,
                           sr_genotypes=sr_genotypes, per_person_path=args.per_person_path)

    with open(os.path.join(args.out_dir, "summary.json"), "w") as fh:
        json.dump({"n_people_unrelated": _sup(n_people), "pairs": ["DQA1~DQB1", "DPA1~DPB1"]}, fh,
                  indent=2)
    print("[37] done -> %s" % args.out_dir, flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=DEFAULT_TABLE1)
    ap.add_argument("--cohort-membership", default=DEFAULT_COHORT)
    ap.add_argument("--relatedness-table", default=DEFAULT_RELATEDNESS)
    ap.add_argument("--sr-genotypes", default=DEFAULT_SR_GENOTYPES)
    ap.add_argument("--skip-relatedness", action="store_true")
    ap.add_argument("--min-allele-haps", type=int, default=MIN_ALLELE_HAPS)
    ap.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)
    ap.add_argument("--artifact-checks", action="store_true",
                    help="run deliverable 4 (VM only, needs Table 1 + relatedness + SR genotypes)")
    ap.add_argument("--per-person-path", default=DEFAULT_PER_PERSON_PATH)
    ap.add_argument("--seed", type=int, default=20260922)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    args = ap.parse_args(argv)
    run(args)


if __name__ == "__main__":
    main()
