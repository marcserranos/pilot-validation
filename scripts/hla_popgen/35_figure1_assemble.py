#!/usr/bin/env python3
"""Figure 1 — assemble the four panels Cole specified into one sheet.

Cole's spec (2026-09-17 call, ask A1):

    (a) the admixture bar plot
    (b) a class I ternary, "the best one you can make" — he liked the HLA-B one as it is
    (c) novel-allele rate by ancestry, broken out by gene — "where is the reference data bias worst?"
    (d) the allele-frequency spectrum, novel vs known, two colours, one plot

Each panel already exists as its own figure. This script does no analysis: it reads the four
rendered PNGs and lays them out as a single figure sheet with panel letters, so there is one
artefact to put in a slide or a manuscript draft.

Why assemble from rendered PNGs rather than re-plot from data
-------------------------------------------------------------
Panels (c) and (d) come from script 33, which reads committed aggregates — those could be
re-plotted. Panels (a) and (b) come from scripts 06 and 10, which committed **figures but not the
underlying per-allele frequency tables**, so they genuinely cannot be re-plotted offline. Rather
than half-compose (two panels live, two pasted) this script treats all four the same way and is
honest that the output is a layout, not a re-analysis. The resolution ceiling is whatever the
source PNGs carry; for a final manuscript figure the two upstream scripts must be rerun to emit
vector output and their tables.

Known gap, deliberately not papered over: panels (a) and (b) are at the ORIGINAL admixture
threshold, not the stricter ~98% Cole and David asked for separately (ask A9). Scripts 06 and 10
have no such parameter today, so that is a code change, not a rerun. The caption says so.

Usage:

    python3 scripts/hla_popgen/35_figure1_assemble.py
"""
import argparse
import json
import os
import sys

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_THIS_DIR))
R = os.path.join(_REPO_ROOT, "reports", "hla_popgen")

PANELS = [
    ("a", os.path.join(R, "06_figures_structure", "lr_full", "admixture_barcode.png"),
     "Genetic ancestry of the long-read cohort"),
    ("b", os.path.join(R, "10_allele_ancestry_geometry", "gene_B", "ternary_B.png"),
     "HLA-B alleles are strongly structured by ancestry"),
    ("c", os.path.join(R, "33_figure1_v3", "figure1_v3_panels_cd.png"),
     "Where IPD-IMGT/HLA is incomplete, and what the novel alleles look like"),
]
DEFAULT_OUT_DIR = os.path.join(R, "35_figure1")


def assemble(panels, out_png, out_pdf, dpi=200):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    from matplotlib.gridspec import GridSpec

    imgs = []
    for letter, path, cap in panels:
        if not os.path.exists(path):
            sys.exit("FATAL: panel %s missing: %s" % (letter, path))
        imgs.append((letter, mpimg.imread(path), cap))

    # The three source images have very different aspect ratios: the admixture barcode is wide
    # and short, the ternary is square, the c+d sheet is wide and tall. Forcing (a) and (b) into
    # one row squashes (a) to a sliver. So each gets its own row, sized from its real aspect
    # ratio, and (b) is centred at the width that keeps it square without dominating the page.
    a_img, b_img, cd_img = imgs[0][1], imgs[1][1], imgs[2][1]
    fig_w = 13.0
    b_frac = 0.52                                   # fraction of page width given to the ternary
    a_h = fig_w * (a_img.shape[0] / a_img.shape[1])
    b_h = fig_w * b_frac * (b_img.shape[0] / b_img.shape[1])
    cd_h = fig_w * (cd_img.shape[0] / cd_img.shape[1])
    fig_h = a_h + b_h + cd_h + 0.5

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(3, 100, height_ratios=[a_h, b_h, cd_h], hspace=0.05, figure=fig)
    lo = int((100 - b_frac * 100) / 2)
    axes = [fig.add_subplot(gs[0, :]),
            fig.add_subplot(gs[1, lo:lo + int(b_frac * 100)]),
            fig.add_subplot(gs[2, :])]
    for ax, (letter, img, cap) in zip(axes, imgs):
        ax.imshow(img)
        ax.axis("off")
        ax.text(0.0 if letter != "b" else -0.06, 1.02, letter, transform=ax.transAxes,
                fontsize=15, fontweight="bold", va="bottom", ha="left", family="DejaVu Sans")

    fig.savefig(out_png, dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_readme(path):
    with open(path, "w") as fh:
        fh.write("""# 35 — Figure 1, assembled

The four panels Cole specified on 2026-09-17, laid out as one sheet.

| Panel | Content | Source |
|---|---|---|
| a | Genetic ancestry of the long-read cohort (admixture barcode) | script 06 |
| b | HLA-B allele ancestry simplex — the ternary Cole said he liked "as it is" | script 10 |
| c | Novelty by gene, ancestry and nomenclature field, with the artifact rate as a flat negative control | script 33 |
| d | Allele-frequency spectrum, novel vs known, density-normalised | script 33 |

*(c and d arrive as one wide sub-sheet from script 33 and occupy the bottom row together.)*

## Draft caption

**Figure 1 | Long-read HLA typing across 11,856 unrelated All of Us participants reveals where the
reference catalogue is incomplete.**
**(a)** Genetic-ancestry composition of the long-read cohort.
**(b)** Per-allele ancestry centroids for HLA-B on the AFR/EUR/AMR simplex; circles are alleles
catalogued in IPD-IMGT/HLA, triangles are alleles first observed here. Alleles sit close to the
vertices, i.e. HLA-B allele frequencies are strongly structured by ancestry.
**(c)** Percentage of called haplotypes carrying sequence absent from IPD-IMGT/HLA, per classical
gene and per ancestry, split by the nomenclature field at which the novelty appears (new protein,
new synonymous CDS, new non-coding only). Black bars show the assembly/annotation artifact rate,
which is flat across ancestries (0.90–3.22% in the five well-powered groups) and therefore acts as
a negative control: the ancestry gradient in the coloured stack is reference incompleteness, not
uneven assembly quality. Grey whiskers show how high the stack could reach if every
disclosure-suppressed count sat at its ceiling.
**(d)** Allele-frequency spectrum over all classical-gene alleles, novel versus already catalogued,
with counts divided by bin width. Novel alleles collapse after a few observations; catalogued
alleles reach hundreds of haplotypes.

## Honest limitations of this version

1. **This is a layout, not a re-analysis.** Panels are composited from rendered PNGs. Scripts 06
   and 10 committed figures but not their underlying per-allele frequency tables, so those two
   panels cannot be re-plotted offline. For a manuscript figure both must be rerun to emit vector
   output plus their tables.
2. **Panels (a) and (b) are at the original admixture threshold**, not the stricter ~98% Cole and
   David asked for (ask A9). Scripts 06 and 10 have no such parameter today — that is a code
   change, not just a rerun.
3. **Panel (d) predates S01's corrected novelty definition.** Its shape is right; the absolute
   novel-allele counts are over-estimates.
""")


def run(args):
    os.makedirs(args.out_dir, exist_ok=True)
    png = os.path.join(args.out_dir, "figure1.png")
    assemble(PANELS, png, os.path.join(args.out_dir, "figure1.pdf"), dpi=args.dpi)
    write_readme(os.path.join(args.out_dir, "README.md"))
    with open(os.path.join(args.out_dir, "summary.json"), "w") as fh:
        json.dump({"panels": [p[0] for p in PANELS],
                   "sources": [os.path.relpath(p[1], _REPO_ROOT) for p in PANELS],
                   "is_layout_not_reanalysis": True,
                   "strict_admixture_applied": False}, fh, indent=2)
    print("[35] Figure 1 assembled -> %s" % png)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    ap.add_argument("--dpi", type=int, default=200)
    run(ap.parse_args(argv))


if __name__ == "__main__":
    main()
