#!/usr/bin/env python3
"""S06 73: methodology schematics for the call-quality report (public reference coordinates only, simulated reads).

  fig_schematic_depth    how the depth of a call is measured: one gene drawn to scale (GENCODE v46 reference transcript: exons,
                         CDS, peptide-binding exons), simulated long reads over it, the per-base depth track, and the measures
                         taken from that track (gene span, exons, CDS, peptide-binding exons, per-base minimum / 10th
                         percentile / share under 10x) with the normalisation by the genome control depth
  fig_schematic_truth    how a typing error is defined: child call -> yak k-mers say which parent transmitted the haplotype
                         -> Mendelian mismatch; and where each quality metric comes from

The reads are simulated for illustration (read lengths and starts are random); they are NOT data. The gene model is the real
one, read from 70's features.tsv (GENCODE v46, GRCh38).

    python3 73_depth_schematic.py --features ~/.../70/features.tsv --gene DQA1 --out reports/hla_popgen/72_call_quality
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _viz_common as V  # noqa: E402

REPO = os.path.dirname(os.path.dirname(HERE))
BLUE, RED, LIGHT, GREY = "#3C5488", "#E64B35", "#B7C4D8", "#7E6148"


def style():
    return matplotlib.rc_context(rc=dict(V.NATURE_RC, **{"font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                                                         "axes.unicode_minus": False}))


def simulate_reads(start, end, depth=17, seed=73, mean_len=14000, sd=4000):
    rng = np.random.default_rng(seed)
    L = end - start
    n = int(depth * (L + 2 * mean_len) / mean_len)
    lens = np.clip(rng.normal(mean_len, sd, n), 4000, 30000).astype(int)
    st = rng.integers(start - mean_len, end, n)
    keep = (st + lens > start) & (st < end)
    return st[keep], (st + lens)[keep]


def stack_rows(starts, ends, gap=300):
    rows, last = [], []
    for s, e in sorted(zip(starts, ends)):
        for i, le in enumerate(last):
            if s > le + gap:
                last[i] = e
                rows.append(i)
                break
        else:
            last.append(e)
            rows.append(len(last) - 1)
    return rows, len(last)


def draw_model(ax, f, y0, h, label=True):
    g = f[f["kind"] == "gene"].iloc[0]
    ax.hlines(y0, g["start"], g["end"], color="0.3", lw=0.7)
    for r in f[f["kind"] == "exons"].itertuples():
        ax.add_patch(Rectangle((r.start, y0 - h / 3), r.end - r.start, 2 * h / 3, fc=LIGHT, ec="none"))
    for r in f[f["kind"] == "cds"].itertuples():
        ax.add_patch(Rectangle((r.start, y0 - h / 2), r.end - r.start, h, fc=BLUE, ec="none"))
    for r in f[(f["kind"] == "ex") & f["ars"].astype(bool)].itertuples():
        ax.add_patch(Rectangle((r.start, y0 - h / 2), r.end - r.start, h, fc=RED, ec="none"))


def fig_schematic_depth(feats, gene, out):
    f = feats[feats["gene"] == gene]
    g = f[f["kind"] == "gene"].iloc[0]
    s, e = int(g["start"]) - 1500, int(g["end"]) + 1500
    rs, re_ = simulate_reads(s, e)
    rows, nrows = stack_rows(rs, re_)
    pos = np.arange(s, e)
    depth = np.zeros(len(pos))
    for a, b in zip(rs, re_):
        depth[max(a, s) - s: min(b, e) - s] += 1
    with style():
        fig = plt.figure(figsize=(V.mm(183), V.mm(150)))
        gs = fig.add_gridspec(3, 1, height_ratios=[1.5, 1.2, 1.0], hspace=0.12, left=0.08, right=0.64, top=0.94, bottom=0.07)
        a1, a2, a3 = (fig.add_subplot(gs[i]) for i in range(3))
        for a, b, r in zip(rs, re_, rows):
            a1.hlines(r, max(a, s), min(b, e), color="0.55", lw=0.9)
        a1.set_xlim(s, e)
        a1.set_ylim(-1, nrows + 1)
        a1.set_yticks([])
        a1.tick_params(labelbottom=False)
        a1.spines["left"].set_visible(False)
        a1.set_title(f"HLA-{gene}: simulated long reads (illustration, not data)", loc="left", fontsize=6.5)
        a2.fill_between(pos, depth, step="mid", color="#8491B4", alpha=0.6, lw=0)
        a2.plot(pos, depth, color=BLUE, lw=0.5)
        a2.axhline(10, color="0.5", lw=0.4, ls="--")
        a2.text(e, 10.4, "10x", fontsize=5, ha="right", color="0.4")
        a2.set_xlim(s, e)
        a2.set_ylabel("depth (x)")
        a2.tick_params(labelbottom=False)
        draw_model(a3, f, 0.5, 0.5)
        a3.set_xlim(s, e)
        a3.set_ylim(0, 1)
        a3.set_yticks([])
        for sp in ("left", "top", "right"):
            a3.spines[sp].set_visible(False)
        a3.set_xlabel(f"chr6 position (GRCh38), HLA-{gene}")
        a3.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: f"{x / 1e6:.3f}"))
        a3.set_xlabel(f"chr6 position (Mb, GRCh38), HLA-{gene}")
        # right-hand text: the measures
        t = fig.add_axes([0.67, 0.07, 0.32, 0.87])
        t.axis("off")
        t.set_xlim(0, 1)
        t.set_ylim(0, 1)
        items = [
            ("Where each base is counted", None, None),
            ("gene span", "mean depth over the whole gene, introns included (what 58 measured)", "0.3"),
            ("exons", "mean over the exons of the reference transcript, UTR included", LIGHT),
            ("CDS", "mean over the coding sequence only (GENCODE v46 CDS rows of the MANE / canonical transcript)", BLUE),
            ("peptide-binding exons", "class I exons 2-3, class II exon 2: the exons that decide 2-field typing", RED),
            ("per-base summaries", "minimum, 10th percentile and share of CDS bases under 10x, from the same track", "0.3"),
            ("relative depth", "every mean divided by the person's genome depth (median of five 2 Mb control windows), so people with"
                               " different overall coverage compare; MAPQ >= 0 and >= 20 both kept", "0.3"),
        ]
        y = 0.98
        for k, (head, body, col) in enumerate(items):
            if body is None:
                t.text(0, y, head, fontsize=7, fontweight="bold", va="top")
                y -= 0.07
                continue
            t.add_patch(Rectangle((0, y - 0.034), 0.05, 0.034, fc=col, ec="none"))
            t.text(0.09, y, head, fontsize=6.5, fontweight="bold", va="top")
            wrapped = _wrap(body, 44)
            t.text(0.09, y - 0.038, wrapped, fontsize=5.8, va="top", linespacing=1.25)
            y -= 0.045 + 0.027 * (wrapped.count("\n") + 1) + 0.03
        V.save_fig(fig, os.path.join(out, "fig_schematic_depth"))


def _wrap(s, n):
    words, lines, cur = s.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > n:
            lines.append(cur)
            cur = w
        else:
            cur = (cur + " " + w).strip()
    lines.append(cur)
    return "\n".join(lines)


def box(ax, x, y, w, h, text, fc="#F3F5F9", ec="0.4", fs=6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.005,rounding_size=0.01", fc=fc, ec=ec, lw=0.5))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, linespacing=1.25)


def arrow(ax, p, q):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=6, lw=0.6, color="0.3"))


def fig_schematic_truth(out):
    with style():
        fig, ax = plt.subplots(figsize=(V.mm(183), V.mm(98)))
        ax.axis("off")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        W, H = 0.27, 0.25
        box(ax, 0.01, 0.68, W, H, "Child's long-read assembly\n(two haplotype contigs, hifiasm)\nImmuannot types the gene\non each haplotype", fs=5.8)
        box(ax, 0.365, 0.68, W, H, "Parents' unique k-mers (yak)\ncounted on the child's gene\n+-10 kb: which parent does\nthis haplotype come from?", fs=5.8)
        box(ax, 0.72, 0.68, W - 0.01, H, "Resolved when one parent has\n>= 5 k-mers and >= 3x the other.\nDoes that parent carry the\nchild's 2-field allele?", fs=5.8)
        box(ax, 0.72, 0.32, W - 0.01, H, "no: typing error (Mendelian mismatch)\nyes: consistent\n(the label used throughout)", fc="#FBEDEA", fs=5.8)
        box(ax, 0.01, 0.32, W, H, "Child's HiFi reads on GRCh38:\nread depth over the gene, exons,\nCDS, peptide-binding exons,\nper base", fc="#EEF3FA", fs=5.8)
        box(ax, 0.365, 0.32, W, H, "Metrics of that call\nImmuannot: template / CDS distance,\nwarnings, ties.  Depth.  Assembly:\ncontigs, partner on the same contig", fc="#EEF3FA", fs=5.8)
        arrow(ax, (0.28, 0.805), (0.365, 0.805))
        arrow(ax, (0.635, 0.805), (0.72, 0.805))
        arrow(ax, (0.85, 0.68), (0.85, 0.57))
        arrow(ax, (0.22, 0.68), (0.42, 0.57))
        arrow(ax, (0.28, 0.445), (0.365, 0.445))
        arrow(ax, (0.635, 0.445), (0.72, 0.445))
        ax.text(0.68, 0.455, "predicts?", fontsize=5.5, ha="center", va="bottom")
        ax.text(0.01, 0.22, "Metrics are attributed to the child's call (primary), to the transmitting parent's call, or to the worse of the two. "
                            "Two further truth sets:\nfirst-degree pairs that share their MHC (no shared 2-field allele) and twin / duplicate pairs "
                            "(unordered genotypes differ).", fontsize=5.8, va="top", linespacing=1.3)
        ax.text(0.01, 0.07, "Pruning: rank calls from worst to best by a metric, remove the worst x%, recompute the error among the rest, "
                            "and compare with removing the same number at random.", fontsize=5.8, va="top")
        V.save_fig(fig, os.path.join(out, "fig_schematic_truth"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--features", required=True)
    ap.add_argument("--gene", default="DQA1")
    ap.add_argument("--out", default=os.path.join(REPO, "reports", "hla_popgen", "72_call_quality"))
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    feats = pd.read_csv(a.features, sep="\t")
    fig_schematic_depth(feats, a.gene, a.out)
    fig_schematic_truth(a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
