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



# CDS lengths (bp, with stop codon) of the modal full-length allele in IPD-IMGT/HLA 3.55.0, read from the `frame=1 NNNbp` headers of
# ~/tools/Immuannot_refdata/CDSseq/<gene>.fa.gz on the VM, 2026-10-20. Public reference data. MICA has no single modal length
# (1,158 / 1,156 / 1,155 bp are the three commonest), 1,155 is shown.
IPD_CDS_BP = {"A": 1098, "B": 1089, "C": 1101, "DPA1": 783, "DPB1": 777, "DQA1": 768, "DQB1": 786, "DRB1": 801, "DRA": 765,
              "MICA": 1155, "TAP1": 2247, "TAP2": 2061}
CLASS_I = {"A", "B", "C"}
CLASS_II_ALPHA = {"DQA1", "DPA1", "DRA"}
DOMAIN = {  # coding-exon ordinal -> protein domain, by class
    "I": {1: "leader", 2: "alpha1", 3: "alpha2", 4: "alpha3"},
    "IIb": {1: "leader", 2: "beta1", 3: "beta2"},
    "IIa": {1: "leader", 2: "alpha1", 3: "alpha2"},
}


def klass(g):
    return "I" if g in CLASS_I else ("IIa" if g in CLASS_II_ALPHA else "IIb")


def fig_schematic_cds(feats, out):
    genes = ["A", "B", "C", "DRB1", "DQB1", "DPB1", "DQA1", "DPA1"]
    with style():
        fig = plt.figure(figsize=(V.mm(183), V.mm(165)))
        gs = fig.add_gridspec(2, 2, height_ratios=[1.55, 1.0], width_ratios=[1.7, 1], hspace=0.22, wspace=0.3, left=0.07, right=0.985,
                              top=0.95, bottom=0.07)
        a = fig.add_subplot(gs[0, :])
        maxlen = 0
        for k, g in enumerate(genes):
            f = feats[(feats["gene"] == g) & (feats["kind"] == "ex")].sort_values("cds_ord")
            x = 0
            y = len(genes) - 1 - k
            for r in f.itertuples():
                col = RED if bool(r.ars) else BLUE
                a.add_patch(Rectangle((x, y - 0.32), r.length, 0.64, fc=col, ec="white", lw=0.4))
                dom = DOMAIN[klass(g)].get(int(r.cds_ord), "TM / cytoplasmic")
                if r.length >= 50:
                    lab = f"{int(r.exon_number)}\n{dom}" if (r.length > 180 or int(r.cds_ord) == 1) else f"{int(r.exon_number)}"
                    a.text(x + r.length / 2, y, lab, ha="center", va="center", fontsize=4.6, color="white", linespacing=1.3)
                x += r.length
            a.text(-12, y, f"HLA-{g}", ha="right", va="center", fontsize=6)
            a.text(x + 12, y, f"{x} bp", ha="left", va="center", fontsize=5, color="0.35")
            maxlen = max(maxlen, x)
        a.set_xlim(-170, maxlen + 130)
        a.set_ylim(-0.7, len(genes) - 0.3)
        a.axis("off")
        a.set_title("Coding exons of the reference transcript, concatenated (width = coding bases; number = GENCODE exon number)",
                    loc="left", fontsize=6.5)
        a.add_patch(Rectangle((maxlen * 0.50, -0.62), 40, 0.22, fc=RED, ec="none"))
        a.text(maxlen * 0.50 + 55, -0.51, "peptide-binding (groove) exon", fontsize=5.5, va="center")
        a.add_patch(Rectangle((maxlen * 0.78, -0.62), 40, 0.22, fc=BLUE, ec="none"))
        a.text(maxlen * 0.78 + 55, -0.51, "other coding exon", fontsize=5.5, va="center")
        V.panel_letter(a, "a", dx=-0.02, dy=1.0)
        b = fig.add_subplot(gs[1, 0])
        rows = []
        for g, ipd in IPD_CDS_BP.items():
            c = feats[(feats["gene"] == g) & (feats["kind"] == "cds")]["length"].sum()
            if c:
                rows.append((g, int(c) + 3, ipd))
        rows.sort(key=lambda r: r[2])
        for i, (g, gc, ipd) in enumerate(rows):
            ok = abs(gc - ipd) <= 3
            b.plot([ipd, gc], [i, i], color="0.7", lw=0.6)
            b.scatter([ipd], [i], s=14, marker="o", fc="white", ec="0.2", lw=0.6, zorder=3, label="IPD-IMGT/HLA 3.55.0, modal allele" if i == 0 else None)
            b.scatter([gc], [i], s=10, marker="s", fc=BLUE if ok else RED, ec="none", zorder=4,
                      label="GENCODE v46 reference CDS + stop codon" if i == 0 else None)
        b.set_yticks(range(len(rows)))
        b.set_yticklabels([f"HLA-{r[0]}" if r[0] not in ("MICA", "TAP1", "TAP2") else r[0] for r in rows])
        b.set_xscale("log")
        b.set_xticks([700, 1000, 1500, 2000])
        b.set_xticklabels(["700", "1000", "1500", "2000"])
        b.minorticks_off()
        b.set_xlabel("coding sequence length (bp)")
        b.legend(loc="lower right", **V.LEGEND_KW)
        b.set_title("CDS check against the IPD reference", loc="left", fontsize=6.5)
        V.panel_letter(b, "b")
        c = fig.add_subplot(gs[1, 1])
        c.axis("off")
        c.set_xlim(0, 1)
        c.set_ylim(0, 1)
        steps = [
            "1  Reference transcript of the gene: MANE Select, else Ensembl canonical, else longest CDS (GENCODE v46, GRCh38)",
            "2  CDS = its CDS rows. UTR, introns and alternative-isoform exons are not CDS",
            "3  Coding exons numbered in order; class I exons 2-3 and class II exon 2 are the groove (peptide-binding) exons",
            "4  Depth is counted per base over each set of positions",
        ]
        yy = 0.98
        for st in steps:
            w = _wrap(st[3:], 40)
            c.text(0.0, yy, st[0], fontsize=7, fontweight="bold", va="top", color=BLUE)
            c.text(0.08, yy, w, fontsize=5.6, va="top", linespacing=1.25)
            yy -= 0.045 + 0.058 * (w.count("\n") + 1)
        V.panel_letter(c, "c", dx=-0.02, dy=1.02)
        V.save_fig(fig, os.path.join(out, "fig_schematic_cds"))



def person(ax, x, y, label, fc="white", r=0.045):
    ax.add_patch(plt.Circle((x, y), r, fc=fc, ec="0.25", lw=0.7, zorder=3))
    ax.text(x, y, label, ha="center", va="center", fontsize=5.5, zorder=4)


def fig_schematic_sets(overlap, out):
    o = dict(zip(overlap["quantity"], overlap["value"]))
    with style():
        fig = plt.figure(figsize=(V.mm(183), V.mm(150)))
        gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.28, wspace=0.12, left=0.03, right=0.985, top=0.94, bottom=0.04)
        a = fig.add_subplot(gs[0, 0])
        a.set_xlim(0, 1)
        a.set_ylim(0, 1)
        a.axis("off")
        a.set_title("Who is in each set", loc="left", fontsize=6.5)
        # family: two parents, child, sibling
        pos = {"P1": (0.18, 0.78), "P2": (0.50, 0.78), "C": (0.34, 0.40), "S": (0.80, 0.40)}
        for k_, (x, y) in pos.items():
            pass
        a.plot([pos["P1"][0], pos["C"][0]], [pos["P1"][1], pos["C"][1]], color=BLUE, lw=2.0, zorder=1)
        a.plot([pos["P2"][0], pos["C"][0]], [pos["P2"][1], pos["C"][1]], color=BLUE, lw=2.0, zorder=1)
        a.plot([pos["C"][0], pos["S"][0]], [pos["C"][1], pos["S"][1]], color=RED, lw=1.4, ls="--", zorder=1)
        a.plot([pos["P1"][0], pos["P2"][0]], [pos["P1"][1], pos["P2"][1]], color="0.75", lw=0.8, zorder=1)
        a.plot([pos["P2"][0], pos["S"][0]], [pos["P2"][1], pos["S"][1]], color=BLUE, lw=0.0)
        person(a, *pos["P1"], "parent")
        person(a, *pos["P2"], "parent")
        person(a, *pos["C"], "child", fc="#EEF3FA")
        person(a, *pos["S"], "sibling", r=0.05)
        a.text(0.02, 0.12, "blue: a trio is a child with both parents (2 child-parent links)\n"
                           "pair set: every first-degree pair (child-parent links, and siblings, dashed red);\n"
                           "parent-child and sibling cannot be told apart here (no IBS0 in All of Us)",
               fontsize=5.4, va="top", linespacing=1.35)
        b = fig.add_subplot(gs[0, 1])
        b.set_xlim(-0.03, 1.03)
        b.set_ylim(0, 1)
        b.axis("off")
        b.set_title("What each test asks, for one gene", loc="left", fontsize=6.5)
        box(b, 0.0, 0.70, 0.46, 0.22, "child's two haplotypes carry X and Y", fs=5.6)
        box(b, 0.54, 0.70, 0.46, 0.22, "mother: Y, Z     father: A, B", fs=5.6)
        box(b, 0.0, 0.36, 0.46, 0.27, "TRIO: yak says haplotype 1 came from\nthe father, haplotype 2 from the mother.\nIs X in the father's alleles?  No: ERROR\n(each haplotype is checked on its own)",
            fc="#EEF3FA", fs=5.3)
        box(b, 0.54, 0.36, 0.46, 0.27, "PAIR (child-mother): do they share any\nallele? Y = Y: pass. A wrong X is\nhidden by the matching Y. Error only if\nNO allele is shared",
            fc="#FBEDEA", fs=5.3)
        b.text(0.0, 0.26, "Pairs also count recombination as an error, and siblings who share only part of their\nhaplotypes; pairs that share nothing across the whole MHC (haplotypes not shared) are\nleft out, so a pair record is one gene of a pair that should share an allele.",
               fontsize=5.3, va="top", linespacing=1.35)
        c = fig.add_subplot(gs[1, 0])
        c.set_xlim(0, 1.45)
        c.set_ylim(0, 1)
        c.set_aspect("equal", adjustable="box")
        c.axis("off")
        c.set_title("How much the two sets overlap (this analysis)", loc="left", fontsize=6.5)
        pairs, edges, both = o["pairs_distinct"], o["trio_edges"], o["edges_in_both"]
        rp = 0.34
        rt = rp * np.sqrt(edges / pairs)
        cx, cy = 0.40, 0.46
        c.add_patch(plt.Circle((cx, cy), rp, fc=RED, alpha=0.22, ec=RED, lw=0.7))
        c.add_patch(plt.Circle((cx - rp + rt + 0.03, cy), rt, fc=BLUE, alpha=0.30, ec=BLUE, lw=0.7))
        c.text(cx - rp + rt + 0.03, cy, f"{int(both)}\nin both", ha="center", va="center", fontsize=6, fontweight="bold")
        c.text(cx, cy + rp + 0.03, f"pair set: {int(pairs)} pairs", ha="center", fontsize=5.6, color=RED)
        c.text(cx, cy - rp - 0.08, f"trio child-parent links: {int(edges)}", ha="center", fontsize=5.6, color=BLUE)
        c.text(0.86, 0.80, _wrap(f"{int(pairs - both)} pairs are not in a trio: siblings, or a parent-child pair whose other parent is not in the data", 28),
               fontsize=5.3, va="top", linespacing=1.3)
        c.text(0.86, 0.42, _wrap(f"{int(edges - both)} trio links are not in the pair set (their MHC does not pass the pair filter, or a call is missing)", 28),
               fontsize=5.3, va="top", linespacing=1.3)
        d = fig.add_subplot(gs[1, 1])
        d.axis("off")
        d.set_xlim(0, 1)
        d.set_ylim(0, 1)
        d.set_title("Size of each set", loc="left", fontsize=6.5)
        rows = [("People", f"{int(o['pair_people'])}", f"{int(o['trio_people'])} ({int(o['trio_children'])} children, rest parents)"),
                ("Relationships", f"{int(o['pairs_distinct'])} pairs", f"{int(o['trio_edges'])} child-parent links"),
                ("Records", f"{int(o['pair_records'])} pair-genes", f"{int(o['trio_records'])} child haplotype calls"),
                ("People in both sets", f"{int(o['people_in_both'])}", ""),
                ("Pair-genes on a trio link", f"{int(o['pair_records_on_trio_edges'])} of {int(o['pair_records'])} ({100 * o['pair_records_on_trio_edges'] / o['pair_records']:.0f}%)", "")]
        d.text(0.30, 0.95, "pair set", fontsize=6, fontweight="bold", color=RED, va="top")
        d.text(0.66, 0.95, "trio set", fontsize=6, fontweight="bold", color=BLUE, va="top")
        y = 0.82
        for lab, x1, x2 in rows:
            d.text(0.0, y, lab, fontsize=5.6, va="top")
            d.text(0.30, y, x1, fontsize=5.6, va="top")
            d.text(0.66, y, _wrap(x2, 22) if x2 else "", fontsize=5.6, va="top", linespacing=1.2)
            y -= 0.15
        V.save_fig(fig, os.path.join(out, "fig_schematic_sets"))



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
    ap.add_argument("--overlap", default=None, help="pair_trio_overlap.tsv (default: next to --features, under 71/)")
    ap.add_argument("--out", default=os.path.join(REPO, "reports", "hla_popgen", "72_call_quality"))
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    feats = pd.read_csv(a.features, sep="\t")
    fig_schematic_depth(feats, a.gene, a.out)
    fig_schematic_cds(feats, a.out)
    fig_schematic_truth(a.out)
    ov = a.overlap or os.path.join(os.path.dirname(os.path.dirname(a.features)), "71", "pair_trio_overlap.tsv")
    if os.path.exists(ov):
        fig_schematic_sets(pd.read_csv(ov, sep="\t"), a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
