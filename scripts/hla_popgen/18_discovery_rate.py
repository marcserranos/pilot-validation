#!/usr/bin/env python3
"""Rate-of-discovery curves -- the derivative of the rarefaction/accumulation curve, per gene.

04_allele_saturation.py plots the CUMULATIVE curve S(N) (`discovery_curves_*.png`) and its Clench
asymptote (`discovery_rate_convergence.png`). Both are dominated by the fact that every classical
gene is still in the near-linear regime, which makes the interesting quantity -- how fast the rate
of discovery is actually decaying -- almost unreadable off a cumulative axis. This script plots
dS/dN directly: **expected new distinct alleles per additional haplotype sequenced.**

WHY THIS NEEDS NO VM ROUND-TRIP
-------------------------------
It reads only aggregate tables 04 already wrote to `reports/hla_popgen/`:
`frequency_spectrum.tsv` (the full incidence frequency-of-frequency spectrum Q_k), and
`allele_richness.tsv` (m = n_haplotypes, S_obs, Q1, Q2). No person-level data, no cds.fa.gz, no
`~/pipeline_outputs/`. It runs on a laptop.

THE IDENTITY THAT MAKES THAT POSSIBLE
-------------------------------------
04's `rarefaction_curve()` builds S(N) by permutation-averaging (200 random orderings of the
sampling units). It never writes those curve points out -- but it doesn't need to, because the
permutation average is a Monte-Carlo estimate of a quantity with a closed form in Q_k alone.

For incidence data with m sampling units (haplotypes), an allele seen in exactly k of them is absent
from a random subsample of size n with probability

    alpha_k(n) = C(m-k, n) / C(m, n)

so the expected rarefaction curve is  E[S(n)] = S_obs - sum_k Q_k * alpha_k(n)  (Colwell et al.
2012, eq. 17). The marginal discovery rate is its forward difference, and since
alpha_k(n+1) = alpha_k(n) * (m-k-n)/(m-n), that difference collapses to

    dS/dn = E[S(n+1)] - E[S(n)] = sum_k Q_k * alpha_k(n) * k / (m - n)          <-- exact, no MC

Two facts make this trustworthy rather than clever:
  * at n = m-1 it reduces to Q1/m -- exactly the Good-Turing `p(next haplotype = new)` that
    04_allele_saturation_report.md already reports. The curve's right-hand endpoint IS that
    published number; the curve is its full-trajectory generalization. `--self-check` asserts this.
  * it is an expectation over ALL orderings, not 200 sampled ones, so it is smoother than (and the
    limit of) the curve 04 plots.

BEYOND TODAY'S SAMPLE (the dashed segment)
------------------------------------------
Past N = m the exact form is unavailable (alpha_k is only defined for n <= m-k), so the projection
uses the Chao/Colwell (2012) incidence extrapolation, differentiated:

    S(m + n*) = S_obs + Q0_hat * (1 - (1-g)^n*),   g = Q1 / (m*Q0_hat + Q1)
    dS/dn*    = -Q0_hat * ln(1-g) * (1-g)^n*

which is continuous with the exact curve at n* = 0 (it tends to Q1/m). Q0_hat is the Chao2 unseen
count -- so the dashed segment inherits Chao2's known instability in the Q2 = 0 regime; genes
flagged `low_confidence_chao2` in `allele_richness.tsv` are drawn but marked. The Clench fit's
closed-form rate, dS/dN = S_max*b/(b+N)^2 (from `discovery_rate_fits.tsv`), is overlaid as the
independent second opinion -- the same disagreement 04's report asks the reader to weigh, shown on
the axis where the two estimators actually differ.

Outputs (reports/hla_popgen/):
  discovery_rate_per_gene.png      2x4 panel, one per classical gene, pooled + by ancestry
  discovery_rate_cross_gene.png    all 8 genes on one axis (pooled) -- the direct comparison
  discovery_rate_curve.tsv         the plotted curves, tidy, for reuse
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.special import gammaln

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_THIS_DIR))
DEFAULT_REPORT_DIR = os.path.join(_REPO_ROOT, "reports", "hla_popgen")

CLASSICAL_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
                   "HLA-DRB1"]
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
# Same palette as 04_allele_saturation.py's plot_discovery_curves, so the derivative figure is
# visually joinable to the cumulative one it differentiates.
ANCESTRY_COLOR = {
    "AFR": "#4C72B0", "AMR": "#DD8452", "EAS": "#55A868",
    "EUR": "#C44E52", "MID": "#8172B2", "SAS": "#937860", "POOLED": "#333333",
}
GENE_COLOR = {
    "HLA-A": "#4C72B0", "HLA-B": "#DD8452", "HLA-C": "#55A868", "HLA-DPA1": "#C44E52",
    "HLA-DPB1": "#8172B2", "HLA-DQA1": "#937860", "HLA-DQB1": "#DA8BC3", "HLA-DRB1": "#64B5CD",
}
CLENCH_COLOR = "#C44E52"
# Rate thresholds annotated on the cross-gene panel: "1 new allele per 10 / 100 / 1000 haplotypes".
RATE_MARKS = [(0.1, "1 per 10"), (0.01, "1 per 100"), (0.001, "1 per 1,000")]
# Below ~10 sampling units the exact curve is dominated by the trivial "almost everything is new"
# regime and just compresses the informative part of a log axis.
MIN_N_PLOTTED = 10


# ---------------------------------------------------------------------------
# The two rate curves
# ---------------------------------------------------------------------------
def exact_rate(Q, m, n_grid):
    """Exact expected marginal discovery rate dS/dn at each n in `n_grid`, for incidence data with
    m sampling units and frequency-of-frequency spectrum Q (dict k -> Q_k).

    dS/dn = sum_k Q_k * alpha_k(n) * k/(m-n), with alpha_k(n) = C(m-k,n)/C(m,n) evaluated in
    log-space via gammaln (m runs to ~15,000 haplotypes for DRB1; the direct binomial ratio
    overflows long before that). alpha_k(n) is identically 0 for n > m-k -- an allele present in k
    units cannot be missing from a subsample larger than m-k -- and that branch must be masked
    rather than computed, since gammaln of a negative argument is not the analytic continuation we
    want here."""
    n = np.asarray(n_grid, dtype=float)
    out = np.zeros_like(n)
    for k, q_k in Q.items():
        if q_k <= 0:
            continue
        ok = n <= (m - k)
        if not ok.any():
            continue
        nn = n[ok]
        log_alpha = ((gammaln(m - k + 1) - gammaln(m - k - nn + 1))
                     - (gammaln(m + 1) - gammaln(m - nn + 1)))
        out[ok] += q_k * np.exp(log_alpha) * k / (m - nn)
    return out


def chao_q0(Q1, Q2, m):
    """Chao2 estimate of the number of alleles present in the population but seen in zero sampling
    units. Bias-corrected form when Q2 > 0; the Q2 = 0 fallback is the standard substitution, and is
    exactly the regime 04's report flags as LOW-CONFIDENCE."""
    if Q1 <= 0:
        return 0.0
    if Q2 > 0:
        return ((m - 1) / m) * Q1 ** 2 / (2 * Q2)
    return ((m - 1) / m) * Q1 * (Q1 - 1) / 2


def chao_extrap_rate(Q1, Q2, m, n_star):
    """dS/dn* for the Chao/Colwell (2012) incidence extrapolation, n* haplotypes beyond the m
    already sampled. Continuous with exact_rate() at n* = 0 (both -> Q1/m)."""
    n_star = np.asarray(n_star, dtype=float)
    q0 = chao_q0(Q1, Q2, m)
    if q0 <= 0 or Q1 <= 0:
        return np.zeros_like(n_star)
    g = Q1 / (m * q0 + Q1)
    if not (0 < g < 1):
        return np.zeros_like(n_star)
    return -q0 * np.log1p(-g) * (1 - g) ** n_star


def clench_rate(n, s_max, b):
    """Closed-form rate of the Clench/Michaelis-Menten fit 04 already produced. Duplicated here
    (three lines) rather than imported, so this script stays runnable off the report TSVs alone."""
    return s_max * b / (b + np.asarray(n, dtype=float)) ** 2


def haplotypes_to_reach(Q1, Q2, m, target_rate):
    """Under the Chao extrapolation, total haplotypes N at which the discovery rate first falls to
    `target_rate`. Returns nan if the rate is already below target, or never reaches it."""
    q0 = chao_q0(Q1, Q2, m)
    if q0 <= 0 or Q1 <= 0:
        return float("nan")
    g = Q1 / (m * q0 + Q1)
    if not (0 < g < 1):
        return float("nan")
    r0 = -q0 * np.log1p(-g)          # rate at n* = 0
    if r0 <= target_rate:
        return float(m)
    return float(m + np.log(target_rate / r0) / np.log1p(-g))


# ---------------------------------------------------------------------------
# Assembling curves from the two report tables
# ---------------------------------------------------------------------------
def spectrum_dict(spec, gene, category, ancestry):
    sub = spec[(spec["gene"] == gene) & (spec["category"] == category)
               & (spec["ancestry"] == ancestry)]
    return {int(k): float(v) for k, v in zip(sub["occurrence_k"], sub["n_alleles"])}


def build_curve(rich, spec, gene, category, ancestry, n_points, extrap_factor):
    """One (gene, category, ancestry) rate curve: exact segment over [MIN_N_PLOTTED, m-1], then the
    Chao-extrapolated segment out to extrap_factor * m. Returns None when the stratum is missing or
    too small to say anything."""
    r = rich[(rich["gene"] == gene) & (rich["category"] == category)
             & (rich["ancestry"] == ancestry)]
    if r.empty:
        return None
    r = r.iloc[0]
    m, S_obs = int(r["n_haplotypes"]), int(r["S_obs"])
    Q1, Q2 = float(r["Q1_uniques"]), float(r["Q2_duplicates"])
    if m <= MIN_N_PLOTTED + 1:
        return None
    Q = spectrum_dict(spec, gene, category, ancestry)
    if not Q:
        return None
    n_obs = np.unique(np.round(np.geomspace(MIN_N_PLOTTED, m - 1, n_points)).astype(int)).astype(float)
    y_obs = exact_rate(Q, m, n_obs)
    n_ext = np.unique(np.round(np.geomspace(1, m * (extrap_factor - 1), n_points)).astype(int)).astype(float)
    y_ext = chao_extrap_rate(Q1, Q2, m, n_ext)
    return {
        "gene": gene, "category": category, "ancestry": ancestry, "m": m, "S_obs": S_obs,
        "Q1": Q1, "Q2": Q2, "good_turing": Q1 / m,
        "low_confidence": bool(r.get("low_confidence_chao2", False)) or Q2 == 0,
        "n_obs": n_obs, "y_obs": y_obs,
        "n_ext": m + n_ext, "y_ext": y_ext,
    }


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _shared_ylim(curves, genes):
    """A single y-range across panels, driven by the data actually drawn. The rate spans ~1.5
    decades here; letting matplotlib pick per-panel decade bounds wasted most of the axis on empty
    space below the lowest curve and made the cross-gene comparison harder to read, which is the
    whole point of the figure."""
    vals = []
    for gene in genes:
        for anc in ANCESTRY_ORDER + ["POOLED"]:
            c = curves.get((gene, anc))
            if c is None:
                continue
            vals.append(c["y_obs"][np.isfinite(c["y_obs"]) & (c["y_obs"] > 0)])
            if anc == "POOLED":
                vals.append(c["y_ext"][np.isfinite(c["y_ext"]) & (c["y_ext"] > 0)])
    allv = np.concatenate(vals)
    return float(allv.min()) / 1.6, float(allv.max()) * 1.4


def plot_per_gene(curves, fits, genes, out_path, category):
    fig, axes = plt.subplots(2, 4, figsize=(20, 9.5), sharey=True)
    ylo, yhi = _shared_ylim(curves, genes)
    for ax, gene in zip(axes.flat, genes):
        pooled = curves.get((gene, "POOLED"))
        if pooled is None:
            ax.axis("off")
            continue
        for anc in ANCESTRY_ORDER:
            c = curves.get((gene, anc))
            if c is None:
                continue
            ax.plot(c["n_obs"], c["y_obs"], color=ANCESTRY_COLOR[anc], lw=1.1, alpha=0.85,
                    label=anc if gene == genes[0] else None)
        ax.plot(pooled["n_obs"], pooled["y_obs"], color=ANCESTRY_COLOR["POOLED"], lw=2.4,
                label="POOLED (observed)" if gene == genes[0] else None)
        ax.plot(pooled["n_ext"], pooled["y_ext"], color=ANCESTRY_COLOR["POOLED"], lw=1.6, ls="--",
                label="Chao2 projection" if gene == genes[0] else None)
        fit = fits.get(gene)
        if fit is not None:
            n_all = np.geomspace(MIN_N_PLOTTED, pooled["n_ext"][-1], 300)
            ax.plot(n_all, clench_rate(n_all, fit["clench_s_max"], fit["clench_b"]),
                    color=CLENCH_COLOR, lw=1.4, ls=":",
                    label="Clench fit rate" if gene == genes[0] else None)
        # Today's sample: the point where "observed" stops and "projected" begins.
        ax.plot([pooled["m"]], [pooled["good_turing"]], "o", ms=6, color="#333333", zorder=5)
        ax.axvline(pooled["m"], color="#999999", lw=0.7, ls="-", alpha=0.6)
        ax.annotate(f"today: N={pooled['m']:,}\n{pooled['good_turing']:.3f} new/hap",
                    xy=(pooled["m"], pooled["good_turing"]), xytext=(-8, -10),
                    textcoords="offset points", fontsize=7.5, color="#333333",
                    ha="right", va="top")
        title = gene + ("  *" if pooled["low_confidence"] else "")
        ax.set_title(title, fontsize=11)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Haplotypes sampled (log)", fontsize=8)
        ax.set_ylabel("New alleles per haplotype (log)", fontsize=8)
        ax.set_ylim(ylo, yhi)
        for y, lab in RATE_MARKS:
            if ylo <= y <= yhi:
                ax.axhline(y, color="#CCCCCC", lw=0.6, zorder=0)
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=10, bbox_to_anchor=(0.5, 1.02),
               frameon=False, fontsize=9)
    fig.suptitle(
        f"Rate of discovery dS/dN -- the derivative of the accumulation curve -- {category} alleles, "
        "classical genes\n"
        "solid = exact expectation from the observed incidence spectrum; dashed = Chao2 projection "
        "beyond today's sample; dotted = Clench-fit rate;  * = Chao2 low-confidence (Q2=0)",
        y=1.08, fontsize=12.5)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_cross_gene(curves, genes, out_path, category):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))
    ax = axes[0]
    for gene in genes:
        c = curves.get((gene, "POOLED"))
        if c is None:
            continue
        ax.plot(c["n_obs"], c["y_obs"], color=GENE_COLOR[gene], lw=2.0, label=gene)
        ax.plot(c["n_ext"], c["y_ext"], color=GENE_COLOR[gene], lw=1.2, ls="--", alpha=0.8)
        ax.plot([c["m"]], [c["good_turing"]], "o", ms=5, color=GENE_COLOR[gene], zorder=5)
    ylo, yhi = _shared_ylim({k: v for k, v in curves.items() if k[1] == "POOLED"}, genes)
    ax.set_ylim(ylo, yhi)
    for y, lab in RATE_MARKS:
        if ylo <= y <= yhi:
            ax.axhline(y, color="#CCCCCC", lw=0.7, zorder=0)
            ax.annotate(lab, xy=(1.0, y), xycoords=("axes fraction", "data"), fontsize=7,
                        color="#888888", ha="right", va="bottom")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Haplotypes sampled (log)")
    ax.set_ylabel("New alleles discovered per additional haplotype (log)")
    ax.set_title("All classical genes, pooled cohort\n(dot = today's sample; dashed = projection)",
                 fontsize=11)
    ax.legend(fontsize=8, frameon=False, ncol=2)
    ax.spines[["top", "right"]].set_visible(False)

    # Right panel: the same information as a single decision-relevant number per gene -- how many
    # haplotypes until the gene stops paying out at a given rate.
    ax = axes[1]
    targets = [0.1, 0.01, 0.001]
    width = 0.26
    xs = np.arange(len(genes))
    for i, t in enumerate(targets):
        vals = []
        for gene in genes:
            c = curves.get((gene, "POOLED"))
            vals.append(np.nan if c is None
                        else haplotypes_to_reach(c["Q1"], c["Q2"], c["m"], t))
        ax.bar(xs + (i - 1) * width, vals, width, label=f"rate falls to {t:g}/hap",
               color=["#4C72B0", "#DD8452", "#55A868"][i])
    for gene, x in zip(genes, xs):
        c = curves.get((gene, "POOLED"))
        if c is not None:
            ax.plot([x - 1.6 * width, x + 1.6 * width], [c["m"], c["m"]], color="#333333", lw=1.4)
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels([g.replace("HLA-", "") for g in genes], fontsize=9)
    ax.set_ylabel("Haplotypes required (log)")
    ax.set_title("Sampling effort to reach a target discovery rate\n"
                 "(Chao2 projection; black bar = haplotypes sequenced today)", fontsize=11)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle(f"Rate of discovery across genes -- {category} alleles", fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def curves_to_frame(curves):
    rows = []
    for (gene, anc), c in curves.items():
        for seg, ns, ys in (("observed", c["n_obs"], c["y_obs"]),
                            ("chao2_projected", c["n_ext"], c["y_ext"])):
            for n, y in zip(ns, ys):
                rows.append({"gene": gene, "category": c["category"], "ancestry": anc,
                             "segment": seg, "n_haplotypes": int(n), "discovery_rate": float(y),
                             "n_haplotypes_observed": c["m"], "S_obs": c["S_obs"],
                             "good_turing_at_m": c["good_turing"],
                             "low_confidence_chao2": c["low_confidence"]})
    return pd.DataFrame(rows)


def self_check(curves):
    """The exact curve's right-hand endpoint must equal Good-Turing Q1/m -- the number
    04_allele_saturation_report.md already publishes. If this drifts, the closed form is wrong."""
    worst, worst_gene = 0.0, None
    for (gene, anc), c in curves.items():
        end = float(exact_rate(c["_Q"], c["m"], [c["m"] - 1.0])[0])
        d = abs(end - c["good_turing"])
        if d > worst:
            worst, worst_gene = d, (gene, anc)
    print(f"[self-check] max |curve endpoint - Good-Turing Q1/m| = {worst:.3e} "
          f"(worst: {worst_gene})", file=sys.stderr)
    assert worst < 1e-9, "closed-form discovery rate disagrees with Good-Turing at n = m-1"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report-dir", default=DEFAULT_REPORT_DIR,
                    help="directory holding 04's allele_richness.tsv / frequency_spectrum.tsv")
    ap.add_argument("--category", default="novel", choices=["novel", "known", "all"],
                    help="allele category to plot (default: novel, matching 04's headline)")
    ap.add_argument("--n-points", type=int, default=400, help="points per curve segment")
    ap.add_argument("--extrap-factor", type=float, default=5.0,
                    help="project the rate out to this multiple of today's haplotype count")
    ap.add_argument("--genes", nargs="*", default=CLASSICAL_GENES)
    ap.add_argument("--self-check", action="store_true",
                    help="assert the closed form reproduces the published Good-Turing number")
    args = ap.parse_args()

    rich = pd.read_csv(os.path.join(args.report_dir, "allele_richness.tsv"), sep="\t")
    spec = pd.read_csv(os.path.join(args.report_dir, "frequency_spectrum.tsv"), sep="\t")
    fits_path = os.path.join(args.report_dir, "discovery_rate_fits.tsv")
    fits = {}
    if os.path.exists(fits_path):
        fdf = pd.read_csv(fits_path, sep="\t")
        for _, r in fdf[(fdf["ancestry"] == "POOLED") & fdf["fit_ok"]].iterrows():
            fits[r["gene"]] = r

    curves = {}
    for gene in args.genes:
        for anc in ANCESTRY_ORDER + ["POOLED"]:
            c = build_curve(rich, spec, gene, args.category, anc, args.n_points,
                            args.extrap_factor)
            if c is not None:
                c["_Q"] = spectrum_dict(spec, gene, args.category, anc)
                curves[(gene, anc)] = c
    if not curves:
        sys.exit(f"no curves built for category={args.category} -- check --report-dir")

    if args.self_check:
        self_check(curves)

    suffix = "" if args.category == "novel" else f"_{args.category}"
    p1 = os.path.join(args.report_dir, f"discovery_rate_per_gene{suffix}.png")
    p2 = os.path.join(args.report_dir, f"discovery_rate_cross_gene{suffix}.png")
    p3 = os.path.join(args.report_dir, f"discovery_rate_curve{suffix}.tsv")
    plot_per_gene(curves, fits, args.genes, p1, args.category)
    plot_cross_gene(curves, args.genes, p2, args.category)
    curves_to_frame(curves).to_csv(p3, sep="\t", index=False)
    for p in (p1, p2, p3):
        print(f"wrote {p}", file=sys.stderr)


if __name__ == "__main__":
    main()
