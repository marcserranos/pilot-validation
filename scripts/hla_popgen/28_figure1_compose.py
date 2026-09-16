#!/usr/bin/env python3
"""Figure 1 for the long-read HLA paper, composed from committed aggregate results.

Runs on a laptop: every input is a committed file under `reports/hla_popgen/`, produced by
24_novelty_by_field.py (field classes, funnel counts), 25_noncoding_novelty_paf.py (non-coding
pooling) and 27_allele_space_coverage.py (coverage, cross-ancestry transfer). No VM, no
participant data.

THE STORY THE FIGURE TELLS (and why these four panels)
------------------------------------------------------
a. **What "novel" survives scrutiny.** 280,695 novel-flagged calls fall to 347 clean novel proteins
   recurrent in unrelated people, and to 8 in the classical genes. The panel is the audit trail:
   most novelty is non-coding, artifact-flagged, or already catalogued under another name.
b. **Novelty is a property of the gene, not of HLA in general.** Per-gene composition shows DRB1's
   60% non-coding novelty (an IPD genomic-reference gap) next to the classical genes' near-absent
   protein novelty, and the non-classical genes where protein novelty is real.
c. **Catalogues do not transfer across ancestries.** Equal-size catalogues built from one ancestry,
   used to type another: strongly asymmetric (AFR->EUR 97.3%, EUR->AFR 83.1%).
d. **How much is left.** Coverage vs catalogue size per ancestry: the classical protein space is
   ~99% covered at this cohort size, and the curves show what a smaller panel would have missed.

Panels a/b answer "what did we find", c/d answer "what would it take to find the rest" - the two
questions the supervisors asked, in the order a reader needs them.

Usage:
    python3 scripts/hla_popgen/28_figure1_compose.py [--report-root reports/hla_popgen]
                                                      [--out-dir reports/hla_popgen/28_figure1]
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _viz_common import (ANCESTRY_COLORS, ANCESTRY_ORDER, CLASSICAL_GENES,  # noqa: E402
                         DEFAULT_REPORT_ROOT, ensure_dir, plt)

FIELD_ORDER = ["known", "f4_noncoding", "f3_synonymous", "f2_protein"]
FIELD_LABEL = {"known": "catalogued allele", "f4_noncoding": "new non-coding sequence",
               "f3_synonymous": "new synonymous CDS", "f2_protein": "new protein"}
FIELD_COLOR = {"known": "#D9D9D9", "f4_noncoding": "#7FB3D5", "f3_synonymous": "#F5B041",
               "f2_protein": "#C0392B"}
SHOWCASE_NONCLASSICAL = ["HLA-E", "HLA-G", "MICA", "TAP1", "TAP2", "HLA-DRB3"]


def read_suppressed_int(series):
    """Counts written as '<20' by the disclosure filter become NaN, not 0 -- a suppressed cell is
    unknown-but-small, and silently zeroing it would understate rare classes."""
    return pd.to_numeric(series, errors="coerce")


def funnel_steps(summary24, summary25):
    """The audit trail of panel a, as (label, count, tone) from the committed summary of 24."""
    h = summary24["headline"]

    def n(key):
        return int(str(h[key]).replace("<20", "19"))

    return [
        ("flagged novel by the caller", n("n_novel_calls_any_depth"), "all"),
        ("new protein (field 2)", n("n_f2_protein_calls"), "all"),
        ("  not artifact-flagged", n("n_f2_protein_calls_clean"), "clean"),
        ("  CDS new, no frameshift", n("n_novel_protein_calls"), "clean"),
        ("distinct novel proteins", n("n_novel_protein_clusters"), "final"),
        ("  in \u22652 unrelated people", n("n_novel_protein_clusters_clean_recurrent_unrelated"),
         "final"),
    ]


def field_composition(counts, genes, scheme="pred", unrelated_only=True):
    """Fraction of haplotype-gene calls in each field class, plus the flagged share of each."""
    d = counts[(counts["ancestry"] == "POOLED") & (counts["ancestry_scheme"] == scheme)
               & (counts["unrelated_only"] == unrelated_only) & counts["gene"].isin(genes)].copy()
    d["n"] = read_suppressed_int(d["n_haplotypes"])
    tot = d.groupby("gene")["n"].sum()
    rows = []
    for gene in genes:
        if gene not in tot.index or not tot[gene]:
            continue
        sub = d[d["gene"] == gene]
        row = {"gene": gene, "total": float(tot[gene])}
        for fc in FIELD_ORDER:
            s = sub[sub["field_class"] == fc]
            row[fc] = float(s["n"].sum()) / float(tot[gene])
            flagged = float(s[s["artifact_label"] != "clean"]["n"].sum())
            row[fc + "_flagged_frac"] = (flagged / float(s["n"].sum())) if float(s["n"].sum()) else 0.0
        rows.append(row)
    return pd.DataFrame(rows).set_index("gene").reindex([g for g in genes if g in tot.index])


def transfer_matrix(cross, resolution="protein", genes=None):
    d = cross[(cross["resolution"] == resolution)
              & (cross["ancestry_ref"] != "UNASSIGNED") & (cross["ancestry_target"] != "UNASSIGNED")]
    if genes is not None:
        d = d[d["gene"].isin(genes)]
    m = d.pivot_table(index="ancestry_ref", columns="ancestry_target", values="cross_coverage")
    ancs = [a for a in ANCESTRY_ORDER if a in m.index]
    return m.reindex(index=ancs, columns=ancs) * 100.0, float(d["m_ref"].median())


def panel_a(ax, steps):
    labels = [s[0] for s in steps]
    vals = [max(s[1], 1) for s in steps]
    tone = {"all": "#95A5A6", "clean": "#5499C7", "final": "#C0392B"}
    y = np.arange(len(vals))[::-1]
    ax.barh(y, vals, color=[tone[s[2]] for s in steps], height=0.68)
    for yi, v, lab in zip(y, vals, labels):
        ax.text(v * 1.25, yi, f"{v:,}", va="center", ha="left", fontsize=6.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.set_xscale("log")
    ax.set_xlim(1, vals[0] * 12)
    ax.set_xlabel("haplotype-gene calls (log scale)", fontsize=7)
    ax.tick_params(axis="x", labelsize=6.5)
    ax.set_title("a   What survives scrutiny", loc="left", fontsize=8, fontweight="bold", pad=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_b(ax, comp):
    genes = list(comp.index)
    x = np.arange(len(genes))
    bottom = np.zeros(len(genes))
    for fc in FIELD_ORDER:
        vals = comp[fc].to_numpy(dtype=float)
        ax.bar(x, vals, bottom=bottom, color=FIELD_COLOR[fc], width=0.78,
               edgecolor="white", linewidth=0.3, label=FIELD_LABEL[fc])
        flag = comp[fc + "_flagged_frac"].to_numpy(dtype=float) * vals
        ax.bar(x, flag, bottom=bottom, color="none", width=0.78, hatch="////",
               edgecolor="#4D4D4D", linewidth=0.0)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([g.replace("HLA-", "") for g in genes], fontsize=6.5, rotation=60,
                       ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("share of calls", fontsize=7)
    ax.tick_params(axis="y", labelsize=6.5)
    # Divider only: the classical/non-classical split is stated in the caption rather than as
    # in-axes text, which collided with the title at this figure width.
    ax.axvline(len(CLASSICAL_GENES) - 0.5, color="#333333", lw=0.8, ls=":")
    ax.legend(fontsize=5.8, frameon=False, ncol=1, loc="center left",
              bbox_to_anchor=(1.01, 0.5), handlelength=1.2)
    ax.set_title("b   Novelty is gene-specific\n     (left of the line: classical genes;"
                 " hatch: artifact-flagged)", loc="left", fontsize=8, fontweight="bold", pad=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_c(ax, mat, m_ref):
    data = mat.to_numpy(dtype=float)
    im = ax.imshow(data, cmap="viridis", vmin=np.nanmin(data), vmax=100)
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, fontsize=6.5)
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels(mat.index, fontsize=6.5)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            if np.isnan(v):
                continue
            ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6,
                    color="white" if v < 92 else "black")
    ax.set_xlabel("population typed", fontsize=7)
    ax.set_ylabel("catalogue built from", fontsize=7)
    ax.set_title(f"c   Catalogues transfer asymmetrically\n     (equal size: {int(m_ref)} haplotypes)",
                 loc="left", fontsize=8, fontweight="bold", pad=6)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=6)
    cb.set_label("% of haplotypes typed", fontsize=6.5)


def panel_d(ax, curves, gene):
    d = curves[(curves["gene"] == gene) & (curves["resolution"] == "protein")]
    for anc in ANCESTRY_ORDER:
        s = d[d["ancestry"] == anc].sort_values("m")
        if s.empty:
            continue
        obs = s[~s["is_extrapolated"]]
        ext = s[s["is_extrapolated"]]
        ax.plot(obs["m"], obs["coverage"] * 100, color=ANCESTRY_COLORS[anc], lw=1.3, label=anc)
        if not ext.empty:
            ax.plot(ext["m"], ext["coverage"] * 100, color=ANCESTRY_COLORS[anc], lw=1.0, ls=":")
    ax.axhline(99, color="#666666", lw=0.7, ls="--")
    ax.text(1.2, 99.15, "99%", fontsize=6, color="#666666")
    ax.set_xscale("log")
    ax.set_ylim(60, 100.6)
    ax.set_xlabel("haplotypes in the catalogue (log scale)", fontsize=7)
    ax.set_ylabel("% of haplotypes typed", fontsize=7)
    ax.tick_params(labelsize=6.5)
    ax.legend(fontsize=5.8, frameon=False, ncol=1, loc="center left",
              bbox_to_anchor=(1.01, 0.5), handlelength=1.2)
    ax.set_title(f"d   How much is left\n     ({gene}, protein level; dotted = extrapolated)",
                 loc="left", fontsize=8, fontweight="bold", pad=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--report-root", default=DEFAULT_REPORT_ROOT)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--curve-gene", default="HLA-B")
    args = ap.parse_args()
    root = args.report_root
    out = args.out_dir or os.path.join(root, "28_figure1")
    ensure_dir(out)

    with open(os.path.join(root, "24_novelty_by_field", "summary.json")) as f:
        s24 = json.load(f)
    s25_path = os.path.join(root, "25_noncoding_novelty_paf", "summary.json")
    s25 = json.load(open(s25_path)) if os.path.exists(s25_path) else {}
    counts = pd.read_csv(os.path.join(root, "24_novelty_by_field", "novelty_field_counts.tsv"),
                         sep="\t", dtype=str)
    counts["unrelated_only"] = counts["unrelated_only"].astype(str).str.lower() == "true"
    cross = pd.read_csv(os.path.join(root, "27_allele_space_coverage", "cross_ancestry_matrix.tsv"),
                        sep="\t")
    curves = pd.read_csv(os.path.join(root, "27_allele_space_coverage", "coverage_curves.tsv"),
                         sep="\t")

    genes_b = CLASSICAL_GENES + [g for g in SHOWCASE_NONCLASSICAL]
    comp = field_composition(counts, genes_b)
    mat, m_ref = transfer_matrix(cross, "protein", CLASSICAL_GENES)

    fig = plt.figure(figsize=(7.4, 5.6))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.75,
                          left=0.135, right=0.88, top=0.90, bottom=0.10)
    panel_a(fig.add_subplot(gs[0, 0]), funnel_steps(s24, s25))
    panel_b(fig.add_subplot(gs[0, 1]), comp)
    panel_c(fig.add_subplot(gs[1, 0]), mat, m_ref)
    panel_d(fig.add_subplot(gs[1, 1]), curves, args.curve_gene)

    for ext in ("png", "pdf"):
        path = os.path.join(out, f"figure1.{ext}")
        fig.savefig(path, dpi=400 if ext == "png" else None, bbox_inches="tight")
        print(f"  wrote {path}", file=sys.stderr)
    plt.close(fig)


if __name__ == "__main__":
    main()
