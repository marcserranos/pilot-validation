#!/usr/bin/env python3
"""S05 67b: a high-resolution view of 67's two-axis pair-error figure (Marc, 2026-10-13: "could we do this with
points, coloured by error or not?").

A scatter of the pair records (one point per relative pair x gene, coloured by disagreement) is individual-level
data, so it is drawn ON THE VM ONLY (--internal-scatter) for QC viewing in JupyterLab. The public version is a
smooth error surface: at every point of a fine grid over (lower depth of the pair, higher depth of the pair), the
Gaussian-kernel-weighted share of records that disagree. Each grid value averages many records, so no count is
shown; points whose effective number of records ((sum w)^2 / sum w^2) is under MIN_EFF are left blank.

Input (VM): ~/s05/internal/67/pair_records.INTERNAL.tsv (67 part A): dep_lo, dep_hi (gene depth / platform
median, MAPQ 0), bad (pair shares no 2-field allele).

VM:     python3 67b_pair_error_surface.py                     -> ~/s05/results/67b/67b_pair_error_surface.tsv
        python3 67b_pair_error_surface.py --internal-scatter  -> also ~/s05/internal/67b/pair_scatter.INTERNAL.png
Here:   python3 67b_pair_error_surface.py --figure            -> reports/.../fig_pair_error_surface.png
"""
import argparse
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
PUB = os.path.join(REPO, "reports", "hla_popgen", "67_phasing_error_vs_coverage")
REC = os.path.expanduser("~/s05/internal/67/pair_records.INTERNAL.tsv")
OUT = os.path.expanduser("~/s05/results/67b/")
INT = os.path.expanduser("~/s05/internal/67b/")
LO, HI, STEP = 0.30, 1.80, 0.05     # grid in normalised depth (1.0 = platform median at that gene)
BW = 0.08                           # kernel sd, same units
MIN_EFF = 100                       # blank a grid point resting on fewer effective records


def read_int(path):
    with open(path) as f:
        first = f.readline()
    return pd.read_csv(path, sep="\t", skiprows=1 if first.startswith("#") else 0)


def surface(r, lo=LO, hi=HI, step=STEP, bw=BW, min_eff=MIN_EFF):
    """Kernel-smoothed error over the (dep_lo <= dep_hi) triangle. Returns x, y, error_pct, eff_n (rounded to 10)."""
    x = pd.to_numeric(r["dep_lo"], errors="coerce").to_numpy(float)
    y = pd.to_numeric(r["dep_hi"], errors="coerce").to_numpy(float)
    b = pd.to_numeric(r["bad"], errors="coerce").to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(b)
    x, y, b = x[ok], y[ok], b[ok]
    g = np.round(np.arange(lo, hi + step / 2, step), 4)
    rows = []
    for gx in g:
        wx = np.exp(-0.5 * ((x - gx) / bw) ** 2)
        for gy in g:
            if gy < gx:
                continue
            w = wx * np.exp(-0.5 * ((y - gy) / bw) ** 2)
            sw, sw2 = w.sum(), (w * w).sum()
            eff = sw * sw / sw2 if sw2 > 0 else 0.0
            pct = 100 * (w * b).sum() / sw if eff >= min_eff else np.nan
            rows.append({"dep_lo": gx, "dep_hi": gy, "error_pct": "" if np.isnan(pct) else f"{pct:.2f}",
                         "eff_n": "" if eff < min_eff else str(int(round(eff / 10.0)) * 10)})
    return pd.DataFrame(rows)


def internal_scatter(r, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 6.4))
    ok, bad = r[r["bad"] == 0], r[r["bad"] == 1]
    ax.scatter(ok["dep_lo"], ok["dep_hi"], s=3, c="#9fb7d6", alpha=0.35, lw=0, label="pair agrees")
    ax.scatter(bad["dep_lo"], bad["dep_hi"], s=9, c="#e34948", alpha=0.9, lw=0, label="pair disagrees")
    ax.set(xlim=(0, 2.2), ylim=(0, 2.2), xlabel="lower depth of the pair (÷ platform median at the gene)",
           ylabel="higher depth of the pair", title="INTERNAL - individual-level, do not export")
    ax.legend(loc="lower right")
    fig.savefig(path, dpi=150, bbox_inches="tight")


def figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    s = pd.read_csv(os.path.join(PUB, "67b_pair_error_surface.tsv"), sep="\t")
    g = np.sort(s["dep_lo"].unique())
    Z = np.full((len(g), len(g)), np.nan)
    E = np.full((len(g), len(g)), np.nan)
    ix = {v: k for k, v in enumerate(g)}
    for a, b_, p, e in zip(s["dep_lo"], s["dep_hi"], s["error_pct"], s["eff_n"]):
        Z[ix[b_], ix[a]] = p
        E[ix[b_], ix[a]] = e
    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    vmax = np.nanpercentile(Z, 99)
    im = ax.pcolormesh(g, g, Z, cmap="magma_r", vmin=0, vmax=vmax, shading="nearest")
    cs = ax.contour(g, g, E, levels=[200, 500, 1000], colors="#52514e", linewidths=0.6, linestyles=":")
    ax.clabel(cs, fmt=lambda v: f"~{int(v)} records", fontsize=6)
    ax.plot([g[0], g[-1]], [g[0], g[-1]], color="#898781", lw=0.6)
    for v in (0.75, 1.0):
        ax.axvline(v, color="#898781", lw=0.5, ls="--")
    ax.text(0.755, 1.70, "0.75", fontsize=7, color="#52514e", va="top")
    cb = fig.colorbar(im, ax=ax, shrink=0.85)
    cb.set_label("pair-genes with no shared allele (%), smoothed")
    ax.set_xlabel("lower depth of the pair (÷ platform median at the gene)")
    ax.set_ylabel("higher depth of the pair")
    ax.set_title("Relative-pair disagreement over both members' depth (kernel sd 0.08)", fontsize=10)
    ax.set_xlim(0.38, 1.72)
    ax.set_ylim(0.58, 1.72)
    ax.set_aspect("equal")
    fig.savefig(os.path.join(PUB, "fig_pair_error_surface.png"), dpi=170, bbox_inches="tight")
    print("wrote", os.path.join(PUB, "fig_pair_error_surface.png"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figure", action="store_true")
    ap.add_argument("--internal-scatter", action="store_true")
    a = ap.parse_args(argv)
    if a.figure:
        return figure()
    r = read_int(REC)
    os.makedirs(OUT, exist_ok=True)
    s = surface(r)
    s.to_csv(OUT + "67b_pair_error_surface.tsv", sep="\t", index=False)
    print(f"records {len(r)}; grid points {len(s)}, shown {int((s['error_pct'] != '').sum())}")
    if a.internal_scatter:
        os.makedirs(INT, exist_ok=True)
        internal_scatter(r, INT + "pair_scatter.INTERNAL.png")
        print("wrote", INT + "pair_scatter.INTERNAL.png")


if __name__ == "__main__":
    main()
