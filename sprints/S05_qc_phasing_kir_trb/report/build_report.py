#!/usr/bin/env python3
"""Build the S05 sprint report: src/*.html (joined in name order) -> ../S05_REPORT.html (one self-contained file).

What the builder does:
  - <img data-src="reports/...png"> becomes an inline data: URI (PNG figures from the repo; wide
    images are downscaled to 1600 px);
  - <figure class="plate" id="f-..."> elements are numbered in order, and <a class="fref" href="#f-...">
    links get "Figure N";
  - <!--SVG:name--> placeholders become inline SVG charts drawn here from committed PUBLIC tables
    (theme-aware: colours come from the page's CSS tokens).

Only PUBLIC repo files are read. Run from anywhere:
    python3 sprints/S05_qc_phasing_kir_trb/report/build_report.py [--extra-root DIR ...]
--extra-root: extra directories searched for repo-relative paths that are not in the checkout yet
(e.g. result tables still on a branch).
"""
import argparse
import base64
import html
import io
import os
import re

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
SRC_DIR = os.path.join(HERE, "src")   # report source, one file per section, joined in name order
OUT = os.path.join(os.path.dirname(HERE), "S05_REPORT.html")
ANC = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
ROOTS = [REPO]


def find(rel):
    for r in ROOTS:
        p = os.path.join(r, rel)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(rel)


def img_uri(rel, max_w=1600):
    p = find(rel)
    data = open(p, "rb").read()
    try:
        from PIL import Image
        im = Image.open(io.BytesIO(data))
        if im.width > max_w:
            im = im.convert("RGB").resize((max_w, round(im.height * max_w / im.width)), Image.LANCZOS)
            buf = io.BytesIO()
            im.save(buf, "PNG", optimize=True)
            data = buf.getvalue()
    except ImportError:
        pass
    return "data:image/png;base64," + base64.b64encode(data).decode()


def esc(s):
    return html.escape(str(s), quote=True)


# ----------------------------------------------------------------- inline charts

def svg_hbars(rows, xmax, unit="%", width=640, row_h=26, label_w=150, ticks=None, ref=None,
              title=None):
    """rows: (label, value, css_class, text) ; value None = no bar (text only)."""
    h = row_h * len(rows) + 34
    plot_w = width - label_w - 70
    sx = lambda v: label_w + plot_w * v / xmax  # noqa: E731
    out = [f'<svg viewBox="0 0 {width} {h}" class="chart" role="img"'
           f'{f" aria-label=\"{esc(title)}\"" if title else ""}>']
    for t in ticks or []:
        x = sx(t)
        out.append(f'<line x1="{x:.1f}" y1="4" x2="{x:.1f}" y2="{h - 26}" class="grid"/>')
        out.append(f'<text x="{x:.1f}" y="{h - 10}" class="tick" text-anchor="middle">{t:g}{unit}</text>')
    if ref is not None:
        x = sx(ref[0])
        out.append(f'<line x1="{x:.1f}" y1="0" x2="{x:.1f}" y2="{h - 26}" class="refline"/>')
        out.append(f'<text x="{x + 4:.1f}" y="10" class="reflabel">{esc(ref[1])}</text>')
    for i, (lab, v, cls, txt) in enumerate(rows):
        y = 6 + i * row_h
        out.append(f'<text x="{label_w - 8}" y="{y + row_h / 2 + 4:.1f}" class="lab" text-anchor="end">{esc(lab)}</text>')
        if v is not None:
            out.append(f'<rect x="{label_w}" y="{y + 4}" width="{max(sx(v) - label_w, 1):.1f}" height="{row_h - 10}" class="{cls}"/>')
        tx = sx(v) + 6 if v is not None else label_w + 4
        out.append(f'<text x="{tx:.1f}" y="{y + row_h / 2 + 4:.1f}" class="val">{esc(txt)}</text>')
    out.append("</svg>")
    return "".join(out)


def chart_65_zero_copy():
    """65 reads mode, unrelated people: % with no TRBV4-3 copy, next to Aleix's 'never used' rates
    (aleix/hla-resolve-phase1 @ a9d404b, aleix/RNA-seq/context/EXPERIMENTS.md, run 887db29).
    Aggregates from the VM summary (logged in LOG.md 10-12); the PUBLIC tables follow."""
    zero = {"AFR": 38.5, "AMR": 30.1, "EAS": 37.1, "EUR": 19.5, "MID": 32.2, "SAS": 55.1}
    freq = {"AFR": 0.424, "AMR": 0.477, "EAS": 0.421, "EUR": 0.584, "MID": 0.466, "SAS": 0.287}
    never = {"SAS": 43.5, "AFR": 33.3, "EAS": 30.1, "MID": 25.8, "AMR": 25.4, "EUR": 16.7}
    order = sorted(ANC, key=lambda a: zero[a])
    rows = []
    for a in order:
        # MID counts are masked in the PUBLIC table (ALL minus the six ancestries gives the 23 unlisted people),
        # so MID is shown rounded to whole percent only
        txt = (f"{zero[a]:.0f}%  (rounded; counts masked)" if a == "MID"
               else f"{zero[a]:.1f}%  (insertion freq {freq[a]:.2f})")
        rows.append((f"{a}  no copy", zero[a], "bar-a", txt))
        rows.append((f"{a}  never used", never[a], "bar-b", f"{never[a]:.1f}%"))
    rows.append(("ALL  no copy", 33.0, "bar-a", "33.0%  (11,075 people called)"))
    rows.append(("ALL  never used", 27.3, "bar-b", "27.3%  (repertoire cohort)"))
    return svg_hbars(rows, 72, ticks=[0, 20, 40, 60], width=700, label_w=150,
                     title="People with no TRBV4-3 copy vs people who never use TRBV4-3, by ancestry")


def chart_59(rel):
    """59: outcome rate by normalised window depth (MAPQ 0, unrelated, ALL), with Wilson CIs.
    Bins with < 20 people are masked in the table and not drawn."""
    d = pd.read_csv(find(rel), sep="\t", dtype=str, keep_default_na=False)
    d = d[(d.scale == "norm") & (d.mapq == "0") & (d.stratum == "ALL")]
    panels = [("KIR window", [("kir_incomplete", "KIR haplotype incomplete", "s1"),
                              ("kir_fw_missing", "KIR framework gene missing", "s2")]),
              ("MHC window", [("hla_absent", "an HLA gene absent", "s1"),
                              ("mhc_fragmented", "MHC on ≥ 7 contigs", "s2")])]
    bins = ["0.7-0.8", "0.8-0.9", "0.9-1", "1-1.1", ">=1.1"]
    blab = ["0.7–0.8", "0.8–0.9", "0.9–1.0", "1.0–1.1", "≥ 1.1"]
    pw, ph, l, t, gap = 300, 190, 44, 34, 40
    w, h = 2 * (l + pw) + gap, t + ph + 70
    out = [f'<svg viewBox="0 0 {w} {h}" class="chart" role="img" aria-label="Outcome rate by normalised depth">']
    for k, (ptitle, series) in enumerate(panels):
        x0 = k * (l + pw + gap) + l
        sy = lambda v: t + ph * (1 - v / 70)  # noqa: E731
        sx = lambda i: x0 + pw * (i + 0.5) / len(bins)  # noqa: E731
        out.append(f'<text x="{x0}" y="16" class="lab strong">{ptitle}</text>')
        for v in (0, 20, 40, 60):
            out.append(f'<line x1="{x0}" y1="{sy(v):.1f}" x2="{x0 + pw}" y2="{sy(v):.1f}" class="grid"/>')
            out.append(f'<text x="{x0 - 6}" y="{sy(v) + 4:.1f}" class="tick" text-anchor="end">{v}%</text>')
        for i, bl in enumerate(blab):
            out.append(f'<text x="{sx(i):.1f}" y="{t + ph + 16}" class="tick" text-anchor="middle">{bl}</text>')
        out.append(f'<text x="{x0 + pw / 2}" y="{t + ph + 34}" class="tiny" text-anchor="middle">depth ÷ cohort median (1.0 = typical)</text>')
        for j, (oc, lab, cls) in enumerate(series):
            r = d[d.outcome == oc].set_index("depth_bin")
            pts = []
            for i, b in enumerate(bins):
                if b not in r.index or not r.at[b, "pct_events"]:
                    continue
                v, lo, hi = (float(r.at[b, c]) for c in ("pct_events", "pct_events_lo", "pct_events_hi"))
                pts.append((sx(i), sy(v)))
                out.append(f'<line x1="{sx(i):.1f}" y1="{sy(lo):.1f}" x2="{sx(i):.1f}" y2="{sy(hi):.1f}" class="ci {cls}"/>')
            out.append('<polyline points="' + " ".join(f"{x:.1f},{y:.1f}" for x, y in pts) + f'" class="ser {cls}"/>')
            for x, y in pts:
                out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.5" class="pt {cls}"/>')
            ly = t + ph + 50 + 0 * j
            lx = x0 + j * 160
            out.append(f'<circle cx="{lx + 4}" cy="{ly - 4}" r="4" class="pt {cls}"/>')
            out.append(f'<text x="{lx + 12}" y="{ly}" class="tiny">{lab}</text>')
    out.append("</svg>")
    return "".join(out)


def heat_class(v):
    if v is None:
        return "h0"
    for cut, c in ((90, "h6"), (70, "h5"), (50, "h4"), (30, "h3"), (15, "h2"), (5, "h1")):
        if v >= cut:
            return c
    return "h0"


def chart_61_heatmap(rel):
    t = pd.read_csv(find(rel), sep="\t", dtype=str, keep_default_na=False)
    order = ["KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DL5B", "KIR2DS3", "KIR2DP1", "KIR2DL1",
             "KIR3DP1", "KIR2DL4", "KIR3DL1", "KIR3DS1", "KIR2DL5A", "KIR2DS5", "KIR2DS1", "KIR2DS4", "KIR3DL2"]
    order = [g for g in order if g in set(t["gene"])]
    cw, rh, lw, top = 46, 22, 86, 26
    w, h = lw + cw * len(ANC) + 10, top + rh * len(order) + 6
    out = [f'<svg viewBox="0 0 {w} {h}" class="chart heat" role="img" aria-label="KIR gene frequency on complete haplotypes by ancestry">']
    for j, a in enumerate(ANC):
        out.append(f'<text x="{lw + cw * j + cw / 2}" y="16" class="tick" text-anchor="middle">{a}</text>')
    for i, g in enumerate(order):
        y = top + i * rh
        out.append(f'<text x="{lw - 6}" y="{y + rh / 2 + 4}" class="lab" text-anchor="end">{g.replace("KIR", "")}</text>')
        for j, a in enumerate(ANC):
            r = t[(t["gene"] == g) & (t["ancestry"] == a)]
            v = pd.to_numeric(r["hap_freq_pct"], errors="coerce").iloc[0] if len(r) else float("nan")
            v = None if pd.isna(v) else float(v)
            x = lw + cw * j
            out.append(f'<rect x="{x + 1}" y="{y + 1}" width="{cw - 2}" height="{rh - 2}" class="{heat_class(v)}"/>')
            txt = "—" if v is None else f"{v:.0f}"
            out.append(f'<text x="{x + cw / 2}" y="{y + rh / 2 + 4}" class="cell {heat_class(v)}t" text-anchor="middle">{txt}</text>')
    out.append("</svg>")
    return "".join(out)


def chart_trb_map():
    """Schematic of GRCh38 primary vs the ALT KI270803v1 over the TRB V cluster (from the probe-mapping table)."""
    w, h = 720, 150
    out = [f'<svg viewBox="0 0 {w} {h}" class="chart" role="img" aria-label="TRB locus: primary vs ALT">']
    out.append('<text x="10" y="34" class="lab">GRCh38 primary chr7</text>')
    out.append('<text x="10" y="104" class="lab">ALT KI270803v1</text>')
    out.append('<line x1="150" y1="30" x2="710" y2="30" class="axis"/>')
    out.append('<line x1="150" y1="100" x2="710" y2="100" class="axis"/>')
    genes_p = [("V2", 170), ("V3-1", 205), ("V4-1", 240), ("V5-1", 275), ("V6-1", 310), ("V4-2", 345),
               ("V6-2/3", 385), ("V7-2", 420), ("V7-3", 470), ("…", 600)]
    genes_a = [("V2", 170), ("V3-1", 205), ("V4-1", 240), ("V5-1", 275), ("V6-1", 310), ("V4-2", 345),
               ("V3-2", 380), ("V4-3", 412), ("V6-2/3", 444), ("V6-2/3", 478), ("V7-2", 512), ("V7-3", 552),
               ("V7-8", 640)]
    for lab, x in genes_p:
        out.append(f'<rect x="{x - 3}" y="24" width="6" height="12" class="gene"/>')
        out.append(f'<text x="{x}" y="52" class="tiny" text-anchor="middle">{lab}</text>')
    out.append('<rect x="364" y="82" width="132" height="36" class="insA"/>')
    out.append('<text x="430" y="78" class="tiny strong" text-anchor="middle">insertion A, 21.6 kb</text>')
    out.append('<rect x="612" y="82" width="56" height="36" class="insB"/>')
    out.append('<text x="640" y="78" class="tiny strong" text-anchor="middle">segment B, ~40 kb</text>')
    for lab, x in genes_a:
        out.append(f'<rect x="{x - 3}" y="94" width="6" height="12" class="gene"/>')
        out.append(f'<text x="{x}" y="124" class="tiny" text-anchor="middle">{lab}</text>')
    out.append('<path d="M365 36 L364 82" class="arr"/>')
    out.append('<text x="368" y="66" class="tiny">breakpoint chr7:142,346,200–142,352,104</text>')
    out.append('<text x="150" y="146" class="tiny">Not to scale. Genes after insertion A shift by +21.7 kb; TRBV7-8 has no ≥92% copy on primary.</text>')
    out.append("</svg>")
    return "".join(out)


def chart_65_af(rel):
    """65 reads mode, unrelated, ALL: people per 0.05 bin of the TRBV4-3 insertion allele fraction,
    with the 0.2/0.8 genotype cuts and the 0.5 a balanced heterozygote should sit at."""
    d = pd.read_csv(find(rel), sep="\t", dtype=str)
    d = d[(d.feature == "TRBV4-3") & (d.stratum_type == "ancestry") & (d.stratum == "ALL")]
    d = d.assign(lo=d.af_lo.astype(float)).sort_values("lo")
    vals = [(r.lo, r.n_people) for r in d.itertuples()]
    w, h, l, b, t = 720, 300, 56, 44, 16
    pw, ph = w - l - 14, h - b - t
    ymax = 4000
    sx = lambda x: l + pw * x  # noqa: E731
    sy = lambda y: t + ph * (1 - y / ymax)  # noqa: E731
    out = [f'<svg viewBox="0 0 {w} {h}" class="chart" role="img" aria-label="TRBV4-3 insertion allele '
           f'fraction histogram, unrelated people">']
    for y in range(0, ymax + 1, 1000):
        out.append(f'<line x1="{l}" y1="{sy(y):.1f}" x2="{w - 14}" y2="{sy(y):.1f}" class="grid"/>')
        out.append(f'<text x="{l - 6}" y="{sy(y) + 4:.1f}" class="tick" text-anchor="end">{y:,}</text>')
    for x in (0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0):
        out.append(f'<text x="{sx(x):.1f}" y="{h - b + 16}" class="tick" text-anchor="middle">{x:g}</text>')
    out.append(f'<text x="{l + pw / 2:.1f}" y="{h - 6}" class="lab" text-anchor="middle">'
               f'fraction of breakpoint reads carrying the insertion</text>')
    bw = pw * 0.05
    for lo, n in vals:
        x = sx(lo) + 1
        if n.isdigit():
            v = int(n)
            cls = "bar-c" if lo < 0.2 else ("bar-b" if lo >= 0.8 else "bar-a")
            out.append(f'<rect x="{x:.1f}" y="{sy(v):.1f}" width="{bw - 2:.1f}" height="{sy(0) - sy(v):.1f}" class="{cls}"/>')
            out.append(f'<text x="{x + bw / 2 - 1:.1f}" y="{sy(v) - 4:.1f}" class="tiny" text-anchor="middle">{v:,}</text>')
        else:
            out.append(f'<text x="{x + bw / 2 - 1:.1f}" y="{sy(0) - 4:.1f}" class="tiny" text-anchor="middle">{esc(n)}</text>')
    for x, lab in ((0.2, "cut 0.2"), (0.8, "cut 0.8")):
        out.append(f'<line x1="{sx(x):.1f}" y1="{t}" x2="{sx(x):.1f}" y2="{sy(0):.1f}" class="refline"/>')
        out.append(f'<text x="{sx(x) + 4:.1f}" y="{t + 10}" class="reflabel">{lab}</text>')
    out.append(f'<line x1="{sx(0.5):.1f}" y1="{t + 40}" x2="{sx(0.5):.1f}" y2="{sy(0):.1f}" class="arr"/>')
    out.append(f'<text x="{sx(0.5) - 4:.1f}" y="{t + 36}" class="tiny" text-anchor="end">balanced het = 0.5</text>')
    out.append(f'<text x="{sx(0.1):.1f}" y="{sy(3647) + 40:.1f}" class="tiny strong" text-anchor="middle">0 copies</text>')
    out.append(f'<text x="{sx(0.62):.1f}" y="{sy(811) - 22:.1f}" class="tiny strong" text-anchor="middle">1 copy (mode ≈ 0.6)</text>')
    out.append(f'<text x="{sx(0.93):.1f}" y="{sy(2237) - 20:.1f}" class="tiny strong" text-anchor="middle">2 copies</text>')
    out.append(f'<line x1="{l}" y1="{sy(0):.1f}" x2="{w - 14}" y2="{sy(0):.1f}" class="axis"/>')
    out.append("</svg>")
    return "".join(out)


# ----------------------------------------------------------------- 67: where errors happen

def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def err_class(v):
    if v is None:
        return "h0"
    for cut, c in ((12, "h6"), (9, "h5"), (6.5, "h4"), (4.5, "h3"), (3, "h2"), (1.5, "h1")):
        if v >= cut:
            return c
    return "h0"


def _range_lab(r):
    a, b = r.split("..")
    if a == "-inf":
        return f"< {float(b):.2f}"
    if b == "inf":
        return f"≥ {float(a):.2f}"
    return f"{float(a):.2f}–{float(b):.2f}"


def chart_67_pair_grid(rel):
    """67 A: discordance of relative pairs, binned on (lower depth of the pair) x (higher depth of the pair),
    tertiles of all pair-member depths. Upper triangle only (lower <= higher by construction)."""
    d = pd.read_csv(find(rel), sep="\t", dtype=str, keep_default_na=False)
    nb = int(d["x_bin"].astype(int).max())
    labs = {int(r.x_bin): _range_lab(r.x_range) for r in d.itertuples()}
    cw, ch, l, t = 150, 74, 120, 30
    w, h = l + cw * nb + 20, t + ch * nb + 70
    out = [f'<svg viewBox="0 0 {w} {h}" class="chart heat heatwide" role="img" '
           f'aria-label="Relative-pair discordance by the depth of both members">']
    for k in range(1, nb + 1):
        x = l + cw * (k - 1)
        out.append(f'<text x="{x + cw / 2}" y="{t + ch * nb + 18}" class="tick" text-anchor="middle">{esc(labs[k])}</text>')
        y = t + ch * (nb - k)
        out.append(f'<text x="{l - 8}" y="{y + ch / 2 + 4}" class="tick" text-anchor="end">{esc(labs[k])}</text>')
    out.append(f'<text x="{l + cw * nb / 2}" y="{t + ch * nb + 40}" class="lab" text-anchor="middle">'
               f'depth of the LOWER-depth relative at this gene (÷ platform median)</text>')
    out.append(f'<text x="14" y="{t + ch * nb / 2}" class="lab" text-anchor="middle" '
               f'transform="rotate(-90 14 {t + ch * nb / 2})">depth of the HIGHER-depth relative</text>')
    for r in d.itertuples():
        ix, iy = int(r.x_bin), int(r.y_bin)
        x, y = l + cw * (ix - 1), t + ch * (nb - iy)
        v = _num(r.error_pct)
        out.append(f'<rect x="{x + 2}" y="{y + 2}" width="{cw - 4}" height="{ch - 4}" class="{err_class(v)}"/>')
        n = _num(r.n_records)
        cap = f"< {100 * 19 / n:.1f}%" if v is None and n else ""
        main = f"{v:.1f}%" if v is not None else cap
        ci = f"{r.lo}–{r.hi}%" if v is not None and r.lo else "fewer than 20 errors"
        if v is None:
            v = 100 * 19 / n / 2 if n else None
        out.append(f'<text x="{x + cw / 2}" y="{y + ch / 2 - 4}" class="cell big {err_class(v)}t" text-anchor="middle">{main}</text>')
        out.append(f'<text x="{x + cw / 2}" y="{y + ch / 2 + 12}" class="cell {err_class(v)}t" text-anchor="middle">{esc(ci)}</text>')
        out.append(f'<text x="{x + cw / 2}" y="{y + ch / 2 + 26}" class="cell {err_class(v)}t" text-anchor="middle">n≈{esc(r.n_records)}</text>')
    out.append(f'<text x="{l}" y="{h - 8}" class="tiny">Cells below the diagonal are empty by construction (the lower depth cannot exceed the higher).</text>')
    out.append("</svg>")
    return "".join(out)


def _ci_series_chart(series, xlabels, ymax, xlab, ylab, aria, width=700, height=300, yticks=None):
    """series: list of (legend, css s1/s2/s3, [(value, lo, hi) or ('cap', bound) or None per x]).
    A ('cap', bound) point is a masked rate (fewer than 20 events): drawn as a hollow marker at the bound."""
    l, t, b = 56, 20, 74
    pw, ph = width - l - 16, height - t - b
    n = len(xlabels)
    sx = lambda i: l + pw * (i + 0.5) / n  # noqa: E731
    sy = lambda v: t + ph * (1 - v / ymax)  # noqa: E731
    out = [f'<svg viewBox="0 0 {width} {height}" class="chart" role="img" aria-label="{esc(aria)}">']
    for v in yticks or []:
        out.append(f'<line x1="{l}" y1="{sy(v):.1f}" x2="{l + pw}" y2="{sy(v):.1f}" class="grid"/>')
        out.append(f'<text x="{l - 6}" y="{sy(v) + 4:.1f}" class="tick" text-anchor="end">{v:g}%</text>')
    for i, xl in enumerate(xlabels):
        out.append(f'<text x="{sx(i):.1f}" y="{t + ph + 16}" class="tick" text-anchor="middle">{esc(xl)}</text>')
    out.append(f'<text x="{l + pw / 2}" y="{t + ph + 36}" class="lab" text-anchor="middle">{esc(xlab)}</text>')
    out.append(f'<text x="12" y="{t + ph / 2}" class="lab" text-anchor="middle" transform="rotate(-90 12 {t + ph / 2})">{esc(ylab)}</text>')
    for j, (leg, cls, vals) in enumerate(series):
        off = (j - (len(series) - 1) / 2) * 6
        pts = []
        for i, p in enumerate(vals):
            if p is None:
                continue
            if p[0] == "cap":
                x = sx(i) + off
                out.append(f'<path d="M{x - 4:.1f},{sy(p[1]) - 3:.1f} L{x + 4:.1f},{sy(p[1]) - 3:.1f} L{x:.1f},{sy(p[1]) + 4:.1f} Z" class="cap {cls}"/>')
                continue
            v, lo, hi = p
            x = sx(i) + off
            pts.append((x, sy(v)))
            if lo is not None and hi is not None:
                out.append(f'<line x1="{x:.1f}" y1="{sy(min(hi, ymax)):.1f}" x2="{x:.1f}" y2="{sy(lo):.1f}" class="ci {cls}"/>')
        if len(pts) > 1:
            out.append('<polyline points="' + " ".join(f"{x:.1f},{y:.1f}" for x, y in pts) + f'" class="ser {cls}"/>')
        for x, y in pts:
            out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.5" class="pt {cls}"/>')
        lx, ly = l + j * 250, height - 12
        out.append(f'<circle cx="{lx + 4}" cy="{ly - 4}" r="4" class="pt {cls}"/>')
        out.append(f'<text x="{lx + 12}" y="{ly}" class="tiny">{esc(leg)}</text>')
    out.append("</svg>")
    return "".join(out)


def _trip(r, col="error_pct", ncol="n_records"):
    v = _num(getattr(r, col))
    if v is None:
        n = _num(getattr(r, ncol, None))
        return ("cap", 100 * 19 / n) if n else None
    return (v, _num(r.lo), _num(r.hi))


def chart_67_pair_deciles(rel):
    d = pd.read_csv(find(rel), sep="\t", dtype=str, keep_default_na=False)
    lo = d[d["axis"] == "lower_depth_of_pair"].sort_values("bin", key=lambda s: s.astype(int))
    hi = d[d["axis"] == "higher_depth_of_pair"].sort_values("bin", key=lambda s: s.astype(int))
    xl = [f"{int(b)}" for b in lo["bin"]]
    vals = [x[0] for s in (lo, hi) for x in [_trip(r) for r in s.itertuples()] if x and x[0] != "cap"]
    ymax = max(16, (max(vals + [0]) // 4 + 1) * 4)
    return _ci_series_chart(
        [("ranked by the LOWER depth of the two relatives", "s1", [_trip(r) for r in lo.itertuples()]),
         ("ranked by the HIGHER depth", "s2", [_trip(r) for r in hi.itertuples()])],
        xl, ymax, "tenth of pair records (1 = lowest depth, 10 = highest)", "pairs disagreeing",
        "Pair discordance by depth decile", yticks=list(range(0, int(ymax) + 1, 4)))


def chart_67_trio(rel):
    d = pd.read_csv(find(rel), sep="\t", dtype=str, keep_default_na=False)
    c = d[d["axis"] == "child_depth"].sort_values("bin")
    p = d[d["axis"] == "transmitting_parent_depth"].sort_values("bin")
    xl = ["lowest third", "middle third", "highest third"][:len(c)]
    return _ci_series_chart(
        [("ranked by the child's depth", "s1", [_trip(r) for r in c.itertuples()]),
         ("ranked by the transmitting parent's depth", "s2", [_trip(r) for r in p.itertuples()])],
        xl, 12, "depth at the gene (thirds of all trio calls)", "Mendel errors",
        "Trio Mendel error by depth", width=620, height=280, yticks=[0, 4, 8, 12])


TOPO_LAB = {"classI": "within class I (A, C, B)", "classIII": "within class III",
            "classII_DRDQ": "within DR–DQ, TAP", "classII_DP": "within DP",
            "classI->classIII": "HLA-B → MICA", "classI->classII_DRDQ": "class I → DR–DQ",
            "classIII->classII_DRDQ": "MICA → DRB (~1.1 Mb)",
            "classII_DRDQ->classII_DP": "TAP → DPA1 (~220 kb)", "ALL": "all steps"}
TOPO_ORDER = ["classI", "classI->classIII", "classIII->classII_DRDQ", "classII_DRDQ",
              "classII_DRDQ->classII_DP", "classII_DP", "classI->classII_DRDQ", "ALL"]


def chart_67_topology(rel):
    d = pd.read_csv(find(rel), sep="\t", dtype=str, keep_default_na=False)
    d = d[d["view"] == "interval_class"]
    d = d.assign(o=d["level"].map(lambda x: TOPO_ORDER.index(x) if x in TOPO_ORDER else 99)).sort_values("o")
    rows = []
    for r in d.itertuples():
        v = _num(r.switch_pct)
        lab = TOPO_LAB.get(r.level, r.level.replace("->", " → ").replace("_", " "))
        cls = "bar-b" if "->" in r.level else ("bar-c" if r.level == "ALL" else "bar-a")
        n = _num(r.n_intervals)
        if v is not None:
            txt = f"{v:.1f}%  ({r.lo}–{r.hi}),  {r.n_intervals} steps"
        elif n:
            txt = f"< {100 * 19 / n:.0f}% (fewer than 20 switches in {r.n_intervals} steps)"
        else:
            txt = "fewer than 20 steps"
        rows.append((lab, v, cls, txt))
    xmax = max([x[1] for x in rows if x[1] is not None] + [10])
    xmax = (int(xmax) // 10 + 1) * 10 + 20
    return svg_hbars(rows, xmax, ticks=list(range(0, xmax - 19, 10)), width=740, label_w=170,
                     title="Switch rate between neighbouring genes, by where in the MHC the step falls")


def chart_67_prune(rel):
    """Prune curve: outcome among the people kept, after removing the worst x% by each metric, against
    removing the same number at random."""
    d = pd.read_csv(find(rel), sep="\t", dtype=str, keep_default_na=False)
    panels = [(o, lab) for o, lab in (("kir_incomplete", "KIR haplotype incomplete"),
                                      ("mhc_fragmented", "MHC on ≥ 7 contigs"),
                                      ("hla_absent", "an HLA gene absent"))
              if f"{o}_pct" in d.columns]
    metrics = list(dict.fromkeys(d["metric"]))
    cls = {m: f"s{i + 1}" for i, m in enumerate(metrics)}
    xs = sorted({int(x) for x in d["prune_pct"]})
    pw, ph, l, t, gap = 200, 170, 40, 30, 36
    w = len(panels) * (l + pw + gap)
    h = t + ph + 60 + 18 * (len(metrics) + 1)
    out = [f'<svg viewBox="0 0 {w} {h}" class="chart" role="img" aria-label="Outcome after pruning the worst people">']
    for k, (o, ptitle) in enumerate(panels):
        x0 = k * (l + pw + gap) + l
        allv = [_num(v) for c in (f"{o}_pct", f"{o}_pct_random") for v in d[c]]
        allv = [v for v in allv if v is not None]
        top = (int(max(allv)) // 10 + 1) * 10
        sy = lambda v, top=top: t + ph * (1 - v / top)  # noqa: E731
        sx = lambda x, x0=x0: x0 + pw * x / max(xs)  # noqa: E731
        out.append(f'<text x="{x0}" y="16" class="lab strong">{ptitle}</text>')
        for v in range(0, top + 1, 10 if top > 20 else 5):
            out.append(f'<line x1="{x0}" y1="{sy(v):.1f}" x2="{x0 + pw}" y2="{sy(v):.1f}" class="grid"/>')
            out.append(f'<text x="{x0 - 5}" y="{sy(v) + 4:.1f}" class="tick" text-anchor="end">{v}%</text>')
        for x in (0, 10, 20, 30):
            if x <= max(xs):
                out.append(f'<text x="{sx(x):.1f}" y="{t + ph + 15}" class="tick" text-anchor="middle">{x}%</text>')
        out.append(f'<text x="{x0 + pw / 2}" y="{t + ph + 32}" class="tiny" text-anchor="middle">share of people removed</text>')
        for m in metrics:
            s = d[d["metric"] == m].assign(x=lambda z: z["prune_pct"].astype(int)).sort_values("x")
            for col, extra in ((f"{o}_pct", ""), (f"{o}_pct_random", " dashed")):
                pts = [(sx(r.x), sy(_num(getattr(r, col)))) for r in s.itertuples() if _num(getattr(r, col)) is not None]
                out.append('<polyline points="' + " ".join(f"{a:.1f},{b:.1f}" for a, b in pts)
                           + f'" class="ser {cls[m]}{extra}"/>')
                if not extra:
                    for a, b in pts:
                        out.append(f'<circle cx="{a:.1f}" cy="{b:.1f}" r="2.8" class="pt {cls[m]}"/>')
    ly = t + ph + 52
    for i, m in enumerate(metrics):
        y = ly + 18 * i
        out.append(f'<line x1="{l}" y1="{y - 4}" x2="{l + 22}" y2="{y - 4}" class="ser {cls[m]}"/>')
        out.append(f'<text x="{l + 28}" y="{y}" class="tiny">remove the worst by {esc(m)}</text>')
    y = ly + 18 * len(metrics)
    out.append(f'<line x1="{l}" y1="{y - 4}" x2="{l + 22}" y2="{y - 4}" class="ser s0 dashed"/>')
    out.append(f'<text x="{l + 28}" y="{y}" class="tiny">dashed, same colour: remove the same number of people at random</text>')
    out.append("</svg>")
    return "".join(out)


R67 = "reports/hla_popgen/67_phasing_error_vs_coverage/"
CHARTS = {
    "65_zero_copy": chart_65_zero_copy,
    "59_rate": lambda: chart_59(
        "reports/hla_popgen/59_depth_outcome_join/vm_2026-10-12/outcome_by_depth_bin.tsv"),
    "trb_map": chart_trb_map,
    "65_af": lambda: chart_65_af(
        "reports/hla_popgen/65_trb_assembly_genotype/reads_unrelated/allele_fraction_histogram.tsv"),
    "67_pair_grid": lambda: chart_67_pair_grid(R67 + "67_pair_depth_grid.tsv"),
    "67_pair_deciles": lambda: chart_67_pair_deciles(R67 + "67_pair_error_by_depth.tsv"),
    "67_trio": lambda: chart_67_trio(R67 + "67_trio_error_by_depth.tsv"),
    "67_topology": lambda: chart_67_topology(R67 + "67_switch_topology.tsv"),
    "67_prune": lambda: chart_67_prune(R67 + "67_prune_curve.tsv"),
    "61_heat": lambda: chart_61_heatmap(
        "reports/hla_popgen/61_kir_by_ancestry_panel/vm_2026-10-12/gene_freq_by_ancestry.tsv"),
}


def build():
    src = "\n".join(open(os.path.join(SRC_DIR, f)).read()
                    for f in sorted(os.listdir(SRC_DIR)) if f.endswith(".html"))
    src = re.sub(r'<img data-src="([^"]+)"', lambda m: f'<img src="{img_uri(m.group(1))}"', src)
    src = re.sub(r"<!--SVG:(\w+)-->", lambda m: CHARTS[m.group(1)](), src)
    ids = re.findall(r'<figure class="plate[^"]*" id="([^"]+)"', src)
    num = {fid: i + 1 for i, fid in enumerate(ids)}
    src = re.sub(r'(<figure class="plate[^"]*" id="([^"]+)"[^>]*>)', lambda m: m.group(1), src)
    src = re.sub(r'<span class="fnum" data-for="([^"]+)"></span>',
                 lambda m: f'<span class="fnum">Figure {num[m.group(1)]}</span>', src)
    src = re.sub(r'<a class="fref" href="#([^"]+)"></a>',
                 lambda m: f'<a class="fref" href="#{m.group(1)}">Figure {num[m.group(1)]}</a>', src)
    missing = re.findall(r'class="fref" href="#([^"]+)"></a>', src)
    assert not missing, missing
    open(OUT, "w").write(src)
    print(f"wrote {os.path.relpath(OUT, REPO)}: {len(src) / 1e6:.2f} MB, {len(ids)} figures")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--extra-root", action="append", default=[])
    a = ap.parse_args()
    ROOTS.extend(a.extra_root)
    build()
