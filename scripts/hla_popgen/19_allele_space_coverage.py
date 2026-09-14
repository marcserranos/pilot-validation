#!/usr/bin/env python3
"""How much of the allele space have we actually explored? -- per classical gene, with the evidence.

04_allele_saturation_report.md answers this twice and gets two different answers: Chao2 says 5.4%
of HLA-A's novel allele space, the Clench fit says 18.8%, and both rows carry a LOW-CONFIDENCE
flag. Neither number should be quoted on its own. This script plots the whole evidence base behind
the question instead of picking a winner.

WHY THIS NEEDS NO VM ROUND-TRIP
-------------------------------
Reads only aggregate tables 04 already wrote to `reports/hla_popgen/`: `allele_richness.tsv`
(S_obs, Q1, Q2, chao2/ace/jackknife), `frequency_spectrum.tsv` (full incidence spectrum Q_k), and
`discovery_rate_fits.tsv` (Clench S_max). Runs on a laptop.

THE TWO QUESTIONS THAT BOTH SOUND LIKE "WHAT % HAVE WE EXPLORED"
----------------------------------------------------------------
They have very different answers, and conflating them is the main way this number gets misreported:

  (1) RICHNESS coverage -- what fraction of the DISTINCT ALLELES that exist have we seen?
      S_obs / S_hat. This is the 5-55% number. It depends entirely on estimating how many alleles
      were missed, which is exactly what is hard here.

  (2) COPY coverage (Good-Turing / Chao-Jost sample coverage) -- what fraction of the ALLELE COPIES
      in the population belong to an allele we have already catalogued?
          C_hat = 1 - (Q1/U) * [ (m-1)Q1 / ((m-1)Q1 + 2*Q2) ]      (Chao & Jost 2012)
      with U = sum_k k*Q_k the total incidence count. This is the 55-96% number, and it is far more
      stable because it never has to guess the size of the unseen tail -- only its total weight.

Both are true at once: for HLA-A we have catalogued a small minority of the distinct novel alleles
that exist, but those already account for ~55% of the novel-allele copies people actually carry.
HLA is a long tail of rare alleles under balancing selection; that is precisely the regime where
these two numbers separate, so the figure reports both rather than collapsing them.

WHY THE RICHNESS ESTIMATORS DISAGREE BY 10x (and why the figure shows a band, not a point)
-------------------------------------------------------------------------------------------
Chao2 scales as Q1^2 / 2*Q2. The novel-allele spectrum here is nearly all singletons -- HLA-A has
Q1 = 652 against Q2 = 17, a ratio of ~38 where a well-behaved sample sits near 2-3. In that regime
Chao2's bias-corrected form grows quadratically in Q1 with almost nothing anchoring it, which is
the same degeneracy 04's report already flags for its Q2 = 0 rows, just short of the hard zero. The
jackknife estimators add at most ~Q1 to S_obs and so are bounded and conservative by construction.
The truth is bracketed, not pinned:

  * Chao2 is a LOWER bound on richness under its assumptions, so S_obs/Chao2 is formally an UPPER
    bound on the fraction discovered -- but a lower bound computed in a degenerate regime is not a
    usable point estimate, and its own bootstrap CI in `allele_richness.tsv` had to be widened to
    bracket its point estimate at all.
  * The jackknives are known to UNDERSTATE richness when a sample is far from coverage-complete,
    so their ~52% is best read as an optimistic ceiling.

So the defensible claim per gene is a range, and the figure draws that range. Panel B is the
control that says the spread is real signal rather than broken machinery: run the identical
estimators on KNOWN (already-in-IPD) alleles from the same haplotypes, where near-saturation is the
expected answer, and they land at 53-93% with a much tighter spread.

Outputs (reports/hla_popgen/):
  allele_space_explored.png     4-panel evidence figure
  allele_space_explored.tsv     every number plotted, tidy
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_THIS_DIR))
DEFAULT_REPORT_DIR = os.path.join(_REPO_ROOT, "reports", "hla_popgen")

CLASSICAL_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
                   "HLA-DRB1"]
# Estimator -> (column in allele_richness.tsv, marker, colour, short gloss for the legend).
# Ordered low-% to high-% so the legend reads in the same direction as the band.
ESTIMATORS = [
    ("Chao2",        "chao2",      "o", "#C44E52", "Q1^2/2Q2 -- degenerate when Q2 is tiny"),
    ("ACE",          "ace",        "s", "#DD8452", "abundance-coverage estimator"),
    ("Clench S_max", None,         "D", "#8172B2", "asymptote of the fitted accumulation curve"),
    ("Jackknife2",   "jackknife2", "^", "#55A868", "conservative; ~S_obs + 2Q1 - Q2"),
    ("Jackknife1",   "jackknife1", "v", "#4C72B0", "most conservative; ~S_obs + Q1"),
]
BAND_COLOR = "#B0B0B0"
CAT_COLOR = {"novel": "#C44E52", "known": "#4C72B0", "all": "#8172B2"}
# A frequency-of-frequency spectrum from a sample near coverage-completeness sits around Q1/Q2 ~ 2-3.
WELL_BEHAVED_Q1Q2 = (1.0, 4.0)


def sample_coverage(Q1, Q2, m, U):
    """Chao & Jost (2012) incidence sample-coverage estimator: the fraction of the population's
    allele copies that belong to an allele already seen at least once. Unlike a richness estimator
    it never has to guess how MANY alleles were missed, only their combined weight, which is why it
    is stable in exactly the singleton-dominated regime that breaks Chao2."""
    if U <= 0 or Q1 <= 0:
        return 1.0 if U > 0 else float("nan")
    denom = (m - 1) * Q1 + 2 * Q2
    if denom <= 0:
        return float("nan")
    return 1.0 - (Q1 / U) * ((m - 1) * Q1 / denom)


def collect(rich, spec, fits, genes, category):
    """One row per gene: S_obs, the spectrum diagnostics, every richness estimator's implied
    % discovered, and copy coverage."""
    rows = []
    for gene in genes:
        r = rich[(rich["gene"] == gene) & (rich["category"] == category)
                 & (rich["ancestry"] == "POOLED")]
        if r.empty:
            continue
        r = r.iloc[0]
        m, S_obs = int(r["n_haplotypes"]), int(r["S_obs"])
        Q1, Q2 = float(r["Q1_uniques"]), float(r["Q2_duplicates"])
        sp = spec[(spec["gene"] == gene) & (spec["category"] == category)
                  & (spec["ancestry"] == "POOLED")]
        U = float((sp["occurrence_k"] * sp["n_alleles"]).sum())
        row = {"gene": gene, "category": category, "n_haplotypes": m, "S_obs": S_obs,
               "Q1": Q1, "Q2": Q2, "U": U,
               "q1_q2_ratio": Q1 / Q2 if Q2 > 0 else np.inf,
               "copy_coverage": sample_coverage(Q1, Q2, m, U),
               "low_confidence_chao2": bool(r.get("low_confidence_chao2", False)),
               "chao2_ci_widened": bool(r.get("chao2_ci_widened", False))}
        for name, col, _, _, _ in ESTIMATORS:
            if col is not None:
                est = float(r[col])
            else:
                # Clench lives in the separate discovery_rate_fits.tsv, which 04 builds ONLY for
                # the novel category (`key = (gene, "novel", anc)`) and which carries no `category`
                # column to say so. Guard on category here: joining it to the known-category S_obs
                # silently produces a meaningless ratio (it read ~8% "discovered" for known
                # alleles, against ~55-93% from every other estimator).
                est = np.nan
                if category == "novel":
                    f = fits[(fits["gene"] == gene) & (fits["ancestry"] == "POOLED")]
                    if len(f) and bool(f["fit_ok"].iloc[0]):
                        est = float(f["clench_s_max"].iloc[0])
            row[f"S_hat_{name}"] = est
            row[f"pct_{name}"] = 100.0 * S_obs / est if est and est == est and est > 0 else np.nan
        pcts = [row[f"pct_{n}"] for n, *_ in ESTIMATORS if row[f"pct_{n}"] == row[f"pct_{n}"]]
        row["pct_min"], row["pct_max"] = (min(pcts), max(pcts)) if pcts else (np.nan, np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


def _forest(ax, df, genes, title, subtitle, show_legend):
    """Estimator-spread dot plot: one row per gene, a grey band spanning the estimators, a marker
    per estimator. Deliberately NOT a bar chart -- a bar implies a point estimate we do not have."""
    ys = np.arange(len(genes))[::-1]
    y_of = {g: y for g, y in zip(genes, ys)}
    for _, r in df.iterrows():
        y = y_of[r["gene"]]
        if r["pct_min"] == r["pct_min"]:
            ax.plot([r["pct_min"], r["pct_max"]], [y, y], color=BAND_COLOR, lw=7, alpha=0.55,
                    solid_capstyle="round", zorder=1)
            ax.annotate(f"{r['pct_min']:.0f}-{r['pct_max']:.0f}%",
                        xy=(r["pct_max"], y), xytext=(8, 0), textcoords="offset points",
                        va="center", fontsize=8.5, color="#555555")
        for name, _, marker, color, _ in ESTIMATORS:
            v = r.get(f"pct_{name}", np.nan)
            if v == v:
                ax.plot([v], [y], marker=marker, ms=6.5, color=color, mec="white", mew=0.6,
                        zorder=3, ls="none")
    ax.set_yticks(ys)
    ax.set_yticklabels([g.replace("HLA-", "") for g in genes], fontsize=10)
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of distinct alleles discovered  (S_obs / estimated total)", fontsize=9)
    ax.set_title(title, fontsize=11.5, loc="left", pad=20)
    ax.annotate(subtitle, xy=(0, 1.008), xycoords="axes fraction", ha="left", va="bottom",
                fontsize=8.5, color="#666666")
    ax.grid(axis="x", color="#EEEEEE", lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    if show_legend:
        handles = [Line2D([], [], marker=m, color=c, ls="none", ms=6.5, label=n)
                   for n, _, m, c, _ in ESTIMATORS]
        ax.legend(handles=handles, fontsize=7.5, frameon=False, loc="lower right", ncol=1)


def plot(novel, known, genes, out_path):
    fig, axes = plt.subplots(2, 2, figsize=(16.5, 11.5))

    # --- A: the headline. Novel alleles, estimator spread. ---
    _forest(axes[0, 0], novel, genes,
            "A.  Novel allele space: what fraction have we found?",
            "each marker = one estimator; grey band = their full range",
            show_legend=True)

    # --- B: the control. Same estimators, known alleles, where near-saturation is expected. ---
    _forest(axes[0, 1], known, genes,
            "B.  Control: known (already-in-IPD) alleles, same haplotypes",
            "estimators agree far better where the sample IS near-saturated",
            show_legend=False)

    # --- C: the other definition -- copy coverage. ---
    ax = axes[1, 0]
    xs = np.arange(len(genes))
    width = 0.38
    for i, (df, cat, lab) in enumerate([
            (novel, "novel", "novel alleles\n(among novel-carrying haplotypes)"),
            (known, "known", "known alleles\n(among all typed haplotypes)")]):
        vals = [100 * df.loc[df["gene"] == g, "copy_coverage"].iloc[0]
                if (df["gene"] == g).any() else np.nan for g in genes]
        bars = ax.bar(xs + (i - 0.5) * width, vals, width, label=lab, color=CAT_COLOR[cat],
                      alpha=0.9)
        for b, v in zip(bars, vals):
            if v == v:
                ax.annotate(f"{v:.0f}", xy=(b.get_x() + b.get_width() / 2, v), xytext=(0, 2),
                            textcoords="offset points", ha="center", fontsize=7.5,
                            color="#444444")
    ax.set_xticks(xs)
    ax.set_xticklabels([g.replace("HLA-", "") for g in genes], fontsize=9)
    ax.set_ylim(0, 112)
    ax.set_ylabel("% of allele copies already catalogued", fontsize=9)
    ax.set_title("C.  The other definition: copy coverage, not allele count",
                 fontsize=11.5, loc="left", pad=20)
    ax.annotate("Chao-Jost sample coverage -- stable because it weighs the unseen tail "
                "rather than counting it", xy=(0, 1.008), xycoords="axes fraction",
                ha="left", va="bottom", fontsize=8.5, color="#666666")
    # Every bar here is >=55% tall, so an in-axes legend sits on top of data wherever it goes.
    ax.legend(fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=2)
    ax.grid(axis="y", color="#EEEEEE", lw=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    # --- D: the diagnostic that explains panel A's width. ---
    ax = axes[1, 1]
    for i, (df, cat, lab) in enumerate([(novel, "novel", "novel alleles"),
                                        (known, "known", "known alleles")]):
        vals = [df.loc[df["gene"] == g, "q1_q2_ratio"].iloc[0]
                if (df["gene"] == g).any() else np.nan for g in genes]
        ax.bar(xs + (i - 0.5) * width, vals, width, label=lab, color=CAT_COLOR[cat], alpha=0.9)
    ax.axhspan(*WELL_BEHAVED_Q1Q2, color="#55A868", alpha=0.16, zorder=0)
    ax.annotate("Q1/Q2 ~ 1-4: regime where Chao2 is well behaved",
                xy=(len(genes) - 0.4, WELL_BEHAVED_Q1Q2[1]), xytext=(0, 4),
                textcoords="offset points", ha="right", fontsize=7.5, color="#3C7A50")
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels([g.replace("HLA-", "") for g in genes], fontsize=9)
    ax.set_ylabel("Q1 / Q2   (singletons per doubleton, log)", fontsize=9)
    ax.set_title("D.  Why panel A is a band and not a number", fontsize=11.5, loc="left", pad=20)
    ax.annotate("Chao2 scales as Q1^2/2Q2; the novel spectrum is almost all singletons",
                xy=(0, 1.008), xycoords="axes fraction", ha="left", va="bottom",
                fontsize=8.5, color="#666666")
    ax.legend(fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=2)
    ax.grid(axis="y", color="#EEEEEE", lw=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("How much of the HLA allele space have we explored?  "
                 "Classical genes, pooled long-read cohort",
                 fontsize=14, y=0.985)
    fig.text(0.5, 0.005,
             "Panels A/B: % discovered = S_obs / estimated total richness, one marker per "
             "estimator. Chao2 is a lower bound on richness under its assumptions (so an upper "
             "bound on % discovered) but is degenerate at the Q1/Q2 ratios in panel D; the "
             "jackknives understate richness far from coverage-completeness, so their end of the "
             "band is an optimistic ceiling. The honest per-gene claim is the band.",
             ha="center", fontsize=8, color="#666666", wrap=True)
    fig.tight_layout(rect=[0, 0.025, 1, 0.975])
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report-dir", default=DEFAULT_REPORT_DIR)
    ap.add_argument("--genes", nargs="*", default=CLASSICAL_GENES)
    args = ap.parse_args()

    rich = pd.read_csv(os.path.join(args.report_dir, "allele_richness.tsv"), sep="\t")
    spec = pd.read_csv(os.path.join(args.report_dir, "frequency_spectrum.tsv"), sep="\t")
    fits_path = os.path.join(args.report_dir, "discovery_rate_fits.tsv")
    fits = pd.read_csv(fits_path, sep="\t") if os.path.exists(fits_path) else pd.DataFrame(
        columns=["gene", "ancestry", "clench_s_max", "fit_ok"])

    novel = collect(rich, spec, fits, args.genes, "novel")
    known = collect(rich, spec, fits, args.genes, "known")
    if novel.empty:
        sys.exit("no novel-category rows found -- check --report-dir")

    out_png = os.path.join(args.report_dir, "allele_space_explored.png")
    out_tsv = os.path.join(args.report_dir, "allele_space_explored.tsv")
    plot(novel, known, args.genes, out_png)
    pd.concat([novel, known], ignore_index=True).to_csv(out_tsv, sep="\t", index=False)

    print(f"{'gene':10s}{'novel %discovered':>20s}{'novel copy cov':>16s}{'known %disc':>14s}",
          file=sys.stderr)
    for g in args.genes:
        n = novel[novel["gene"] == g]
        k = known[known["gene"] == g]
        if n.empty:
            continue
        n = n.iloc[0]
        rng = f"{n['pct_min']:.0f}-{n['pct_max']:.0f}%"
        kr = f"{k.iloc[0]['pct_min']:.0f}-{k.iloc[0]['pct_max']:.0f}%" if not k.empty else "-"
        print(f"{g:10s}{rng:>20s}{100 * n['copy_coverage']:15.1f}%{kr:>14s}", file=sys.stderr)
    for p in (out_png, out_tsv):
        print(f"wrote {p}", file=sys.stderr)


if __name__ == "__main__":
    main()
