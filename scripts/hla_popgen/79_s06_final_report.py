#!/usr/bin/env python3
"""S06 79: the final, standalone S06 report (HLA + KIR), built from AGGREGATE tables only (laptop).

A NEW report (not an edit of 72/74's): what a call is and what "no call" means, the three outcome streams (A missing call, B wrong allele among fully
called genotypes, C novel call), why trios and pairs behave differently, the pruning evidence with its transfer / robustness numbers, and the same
for KIR. Every number in the prose is computed from the tables; a section whose tables are not there yet shows a clearly labelled "pending" box
(nothing is ever invented). Small cells (<20) are already masked upstream and are never unmasked here.

Inputs (aggregate_tables/<folder>/*.tsv; whichever exists is used):
  76   call outcomes (A_*, B_*, C_*, definition_card)        77  KIR call table aggregates        77b  KIR yak origin summaries
  77c  KIR trio / pair classes                               71f call-status audit                 76b  HLA pruning evidence
  78   KIR depth aggregates                                  78b  KIR pruning evidence

Output: <out>/S06_FINAL_REPORT.html (single file, figures embedded as base64) and <out>/figs/*.png|pdf (layout-linted by `save_fig`).
The section id="think-tank" is left empty on purpose.

    python3 scripts/hla_popgen/79_s06_final_report.py \
        --tables reports/hla_popgen/79_s06_final/aggregate_tables --out reports/hla_popgen/79_s06_final
Tests: scripts/hla_popgen/tests/test_79_s06_final_report.py
"""
import argparse
import base64
import html
import io
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _viz_common as V  # noqa: E402

REPO = os.path.dirname(os.path.dirname(HERE))
CLASSICAL = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
C_BLUE, C_RED, C_GREEN, C_GREY, C_ORANGE, C_PURPLE = "#3C5488", "#E64B35", "#00A087", "#8A8F98", "#E18727", "#7E6148"
PALETTE = [C_BLUE, C_RED, C_GREEN, C_ORANGE, C_PURPLE, "#4DBBD5", "#B09C85", "#1B1919"]
COST_RATIOS = (10, 25, 50, 100)      # "calls pruned per error removed" at which pruning stops being worth it (slider-free sensitivity)


def log(m):
    print("[79] " + m, file=sys.stderr, flush=True)


def style():
    """Repo Nature style with Arial first (Helvetica.ttc cannot be loaded by this matplotlib build on the laptop)."""
    return matplotlib.rc_context(rc=dict(V.NATURE_RC, **{"font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                                                         "axes.unicode_minus": False}))


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------
class Tables:
    def __init__(self, root):
        self.root = root
        self._c = {}

    def rd(self, sub, name):
        key = (sub, name)
        if key not in self._c:
            p = os.path.join(self.root, sub, name)
            self._c[key] = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False) if os.path.exists(p) else None
            if self._c[key] is not None:
                self._c[key] = typed(self._c[key])
        return self._c[key]

    def all_files(self):
        out = []
        for sub in sorted(os.listdir(self.root)) if os.path.isdir(self.root) else []:
            d = os.path.join(self.root, sub)
            if os.path.isdir(d):
                out += [(sub, f) for f in sorted(os.listdir(d)) if f.endswith(".tsv")]
        return out


def typed(df):
    """Numeric columns become float; masked cells ('<20') and blanks become NaN in numeric columns; text columns stay text."""
    out = df.copy()
    for c in out.columns:
        s = out[c].astype(str).str.strip()
        conv = pd.to_numeric(s.replace({"": np.nan, "<20": np.nan, "NA": np.nan, "nan": np.nan, "inf": "inf"}), errors="coerce")
        nonblank = s[~s.isin(["", "<20", "NA", "nan"])]
        if len(nonblank) and conv[nonblank.index].notna().all():
            out[c] = conv
    return out


def esc(s):
    return html.escape(str(s))


def fmt_n(x):
    return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else "{:,.0f}".format(x)


def fmt_p(x, d=1):
    return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else ("%." + str(d) + "f%%") % x


def wilson(k, n, z=1.96):
    if not n:
        return (np.nan, np.nan)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (100 * max(c - h, 0), 100 * min(c + h, 1))


# ---------------------------------------------------------------------------
# HTML pieces
# ---------------------------------------------------------------------------
CSS = """
:root{--bg:#fff;--fg:#1b1b1b;--mut:#666;--line:#ddd;--acc:#3C5488;--soft:#f4f6fa;--warn:#fff7e6;--pend:#fdeeee}
@media (prefers-color-scheme: dark){:root{--bg:#15171c;--fg:#e8e8e8;--mut:#9aa0aa;--line:#2d313a;--acc:#8fa6d9;--soft:#1d2027;--warn:#2a2418;--pend:#2c1c1c}}
body{background:var(--bg);color:var(--fg);font:15px/1.55 -apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:0}
main{max-width:1080px;margin:0 auto;padding:24px 20px 90px}
nav{position:sticky;top:0;background:var(--bg);border-bottom:1px solid var(--line);padding:6px 20px;font-size:13px;z-index:5}
nav a{margin-right:12px;color:var(--acc);text-decoration:none}
h1{font-size:27px;margin:.2em 0}h2{font-size:21px;margin:2.2em 0 .4em;border-bottom:1px solid var(--line);padding-bottom:4px}
h3{font-size:16px;margin:1.6em 0 .3em;color:var(--acc)}.sub{color:var(--mut)}
figure{margin:12px 0 4px}figure img{width:100%;height:auto;border:1px solid var(--line);background:#fff}
figcaption{color:var(--mut);font-size:13px;margin-top:3px}
table{border-collapse:collapse;width:100%;font-size:12.5px;margin:8px 0 14px}
th,td{border:1px solid var(--line);padding:3px 6px;text-align:left;vertical-align:top}th{background:var(--soft)}
td.num{text-align:right}
.box{background:var(--soft);border-left:3px solid var(--acc);padding:8px 12px;margin:12px 0}
.mind{background:var(--warn);border-left:3px solid #E18727;padding:8px 12px;margin:12px 0;font-size:14px}
.pending{background:var(--pend);border-left:3px solid #E64B35;padding:8px 12px;margin:12px 0}
.read{color:var(--mut);font-size:13.5px}code{background:var(--soft);padding:0 3px;border-radius:3px}
.tbl-wrap{overflow-x:auto;max-height:420px;overflow-y:auto}
details{margin:8px 0}summary{cursor:pointer;color:var(--acc)}
svg text{font-family:Arial,Helvetica,sans-serif}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:12px}@media(max-width:800px){.grid2{grid-template-columns:1fr}}
"""

JS = """
function flt(id,q){q=q.toLowerCase();document.querySelectorAll('#'+id+' tbody tr').forEach(function(r){
 r.style.display=r.textContent.toLowerCase().indexOf(q)>=0?'':'none';});}
"""


def table_html(df, cls="", max_rows=None, fmt=None, searchable=None):
    fmt = fmt or {}
    d = df if max_rows is None else df.head(max_rows)
    h = "".join("<th>%s</th>" % esc(c) for c in d.columns)
    rows = []
    for _, r in d.iterrows():
        cells = []
        for c in d.columns:
            v = r[c]
            numeric = isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool)
            if c in fmt:
                v = fmt[c](v)
            elif isinstance(v, (float, np.floating)):
                v = "" if np.isnan(v) else ("%.4g" % v)
            cells.append('<td class="%s">%s</td>' % ("num" if numeric else "", esc(v)))
        rows.append("<tr>%s</tr>" % "".join(cells))
    box = ""
    tid = searchable or ""
    if searchable:
        box = '<input type="search" placeholder="filter rows" oninput="flt(\'%s\',this.value)" style="margin:4px 0;padding:3px 6px">' % tid
    return '%s<div class="tbl-wrap"%s><table class="%s"%s><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>' % (
        box, "", cls, (' id="%s"' % tid) if tid else "", h, "".join(rows))


def pending(what, needs):
    return ('<div class="pending"><b>Pending:</b> %s. Needs <code>%s</code>, which is not in the aggregate folder yet; nothing is shown in its place.</div>'
            % (esc(what), esc(needs)))


def mind(text):
    return '<div class="mind"><b>What would change my mind.</b> %s</div>' % text


def img_tag(path, alt, caption=""):
    if not path or not os.path.exists(path):
        return '<div class="pending">[figure not produced: %s]</div>' % esc(alt)
    from PIL import Image
    im = Image.open(path).convert("RGBA")
    bg = Image.new("RGB", im.size, "white")
    bg.paste(im, mask=im.split()[3])
    if bg.width > 1800:
        bg = bg.resize((1800, int(bg.height * 1800 / bg.width)), Image.LANCZOS)
    buf = io.BytesIO()
    bg.save(buf, "PNG", optimize=True)
    b64 = base64.b64encode(buf.getvalue()).decode()
    cap = "<figcaption>%s <span class='sub'>(file: %s, with PDF)</span></figcaption>" % (caption, esc(os.path.basename(path)))
    return '<figure><img alt="%s" src="data:image/png;base64,%s"/>%s</figure>' % (esc(alt), b64, cap)


# ---------------------------------------------------------------------------
# inline SVG schematics
# ---------------------------------------------------------------------------
def svg_flow(boxes, arrows, width=1000, height=220, title=""):
    """boxes: list of dict(x, y, w, h, title, lines, fill); arrows: list of (i, j, label)."""
    parts = ['<svg viewBox="0 0 %d %d" width="100%%" role="img" aria-label="%s" xmlns="http://www.w3.org/2000/svg">' % (width, height, esc(title)),
             '<defs><marker id="ah" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#555"/></marker></defs>']
    for b in boxes:
        parts.append('<rect x="%d" y="%d" width="%d" height="%d" rx="6" fill="%s" stroke="#3C5488" stroke-width="1.2"/>' % (
            b["x"], b["y"], b["w"], b["h"], b.get("fill", "#eef2fa")))
        parts.append('<text x="%d" y="%d" font-size="13" font-weight="bold" fill="#1b1b1b" text-anchor="middle">%s</text>' % (
            b["x"] + b["w"] / 2, b["y"] + 18, esc(b["title"])))
        for k, ln in enumerate(b["lines"]):
            parts.append('<text x="%d" y="%d" font-size="11" fill="#333" text-anchor="middle">%s</text>' % (b["x"] + b["w"] / 2, b["y"] + 36 + 14 * k, esc(ln)))
    for i, j, lab in arrows:
        a, c = boxes[i], boxes[j]
        if c["x"] > a["x"] + a["w"] - 5:     # left to right
            x1, y1, x2, y2 = a["x"] + a["w"], a["y"] + a["h"] / 2, c["x"], c["y"] + c["h"] / 2
        else:                                 # downward
            x1, y1, x2, y2 = a["x"] + a["w"] / 2, a["y"] + a["h"], c["x"] + c["w"] / 2, c["y"]
        parts.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#555" stroke-width="1.4" marker-end="url(#ah)"/>' % (x1, y1, x2, y2))
        if lab:
            parts.append('<text x="%d" y="%d" font-size="10.5" fill="#B22" text-anchor="middle">%s</text>' % ((x1 + x2) / 2, (y1 + y2) / 2 - 5, esc(lab)))
    parts.append("</svg>")
    return "".join(parts)


def schematic_call():
    boxes = [dict(x=10, y=60, w=150, h=84, title="Assembly", lines=["one FASTA per haplotype", "(hap1, hap2), many contigs"]),
             dict(x=190, y=60, w=150, h=84, title="Trim to MHC", lines=["chr6:29.5-33.5 Mb", "+100 kb pad"], fill="#f3f3f3"),
             dict(x=370, y=60, w=170, h=84, title="Immuannot", lines=["best template per gene", "keep if adjNM<=3", "or adjNM/match<0.03"]),
             dict(x=580, y=10, w=190, h=70, title="CALL", lines=["gene row: allele, distances,", "warning, ties, novel tag"], fill="#e6f5ee"),
             dict(x=580, y=110, w=190, h=70, title="NO CALL", lines=["template failed the filter", "or sequence not in window"], fill="#fdeeee"),
             dict(x=810, y=60, w=170, h=84, title="Table 1", lines=["one row per person x hap", "x contig x gene x copy"])]
    return svg_flow(boxes, [(0, 1, ""), (1, 2, ""), (2, 3, "passes"), (2, 4, "fails"), (3, 5, "row"), ], 1000, 200,
                    "How a call is made and what no call means")


def schematic_trio():
    boxes = [dict(x=10, y=10, w=190, h=84, title="Child haplotype call", lines=["gene on a child contig", "window = gene +-10 kb"]),
             dict(x=240, y=10, w=190, h=84, title="yak trioeval", lines=["k-mers unique to father /", "unique to mother in window"]),
             dict(x=470, y=10, w=200, h=84, title="Origin rule", lines=[">=5 k-mers and >=3x the", "other parent -> that parent;", "else ambiguous (dropped)"]),
             dict(x=710, y=10, w=270, h=84, title="Transmitting parent's calls", lines=["both haplotypes at the gene", "2-field name of the child call", "is looked up there"]),
             dict(x=10, y=150, w=300, h=66, title="MATCH", lines=["name in transmitting parent"], fill="#e6f5ee"),
             dict(x=350, y=150, w=300, h=66, title="OTHER PARENT ONLY", lines=["origin error / recombination,", "not a typing error"], fill="#fff7e6"),
             dict(x=690, y=150, w=290, h=66, title="NEITHER", lines=["typing-error candidate"], fill="#fdeeee")]
    return svg_flow(boxes, [(0, 1, ""), (1, 2, ""), (2, 3, ""), (3, 4, ""), (3, 5, ""), (3, 6, "")], 1000, 230, "How the trio label works")


def schematic_streams():
    boxes = [dict(x=10, y=10, w=310, h=130, title="A  Missing call", lines=["unit: person x gene x haplotype", "outcome: no gene annotation", "mostly assembly structure", "(contig break / region end)", "studied cohort-wide"]),
             dict(x=345, y=10, w=310, h=130, title="B  Wrong allele", lines=["unit: child call, fully called", "genotype, fully called parent", "outcome: allele in neither parent", "class MATCH / OTHER ONLY / NEITHER", "tens of events: rates, not models"], fill="#fdeeee"),
             dict(x=680, y=10, w=310, h=130, title="C  Novel call", lines=["alleles outside IPD-IMGT naming", "(new at field 1-2)", "no 2-field name, so no label", "checked by depth and twins", "never counted as an error"], fill="#fff7e6")]
    return svg_flow(boxes, [], 1000, 160, "The three outcome streams")


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def save(fig, figdir, name):
    stem = os.path.join(figdir, name)
    V.save_fig(fig, stem, dpi=300)
    return stem + ".png"


def dotci(ax, y, x, lo, hi, color=C_BLUE, label=None, marker="o"):
    ax.errorbar(x, y, xerr=[np.maximum(np.asarray(x) - np.asarray(lo), 0), np.maximum(np.asarray(hi) - np.asarray(x), 0)], fmt=marker, ms=3.2,
                color=color, ecolor=color, elinewidth=0.8, capsize=0, label=label, lw=0)


def fig_missing(T, figdir):
    g = T.rd("76", "A_missing_by_gene.tsv")
    hc = T.rd("76", "A_missing_by_hap_class.tsv")
    ex = T.rd("76", "A_missing_explained.tsv")
    if g is None or hc is None or ex is None:
        return None
    with style():
        fig, axs = plt.subplots(1, 3, figsize=(V.mm(183), V.mm(70)), gridspec_kw=dict(width_ratios=[1.1, 1, 1.3]), constrained_layout=True)
        ax = axs[0]
        g = g[g["gene"] != "ALL"]
        order = list(g[g["scope"] == "classical8"].sort_values("pct")["gene"])
        for sc, col, off in (("classical8", C_BLUE, -0.15), ("all_fixed_genes", C_RED, 0.15)):
            sub = g[g["scope"] == sc].set_index("gene").reindex(order).dropna(subset=["pct"])
            dotci(ax, [order.index(i) + off for i in sub.index], sub["pct"], sub["lo"], sub["hi"], col, sc.replace("_", " "))
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels(order)
        ax.set_xlabel("slots with no gene annotation (%)")
        ax.legend(frameon=False, loc="lower right")
        V.panel_letter(ax, "a")
        ax = axs[1]
        h = hc[hc["scope"] == "classical8"].copy()
        h["lab"] = h["hap_class"].str.replace("_", " ")
        h = h.sort_values("pct")
        ax.barh(h["lab"], h["pct"], color=C_BLUE)
        ax.set_xscale("symlog", linthresh=0.1)
        ax.set_xlabel("missing (%), log-like scale")
        for i, (p, n) in enumerate(zip(h["pct"], h["n"])):
            ax.text(min(p * 1.15 + 0.02, 150), i, "n=%s" % fmt_n(n), va="center", fontsize=5.5)
        ax.set_xlim(0, 1000)
        V.panel_letter(ax, "b")
        ax = axs[2]
        e = ex[ex["scope"].isin(["classical8", "all_fixed_genes"])]
        cats = list(e.groupby("explanation")["n_missing_slots"].sum().sort_values(ascending=False).index)
        left = {s: 0.0 for s in ("classical8", "all_fixed_genes")}
        for k, cat in enumerate(cats):
            for yi, sc in enumerate(("classical8", "all_fixed_genes")):
                v = e[(e["scope"] == sc) & (e["explanation"] == cat)]["pct_of_missing"]
                v = float(v.iloc[0]) if len(v) else 0.0
                ax.barh(yi, v, left=left[sc], color=PALETTE[k % len(PALETTE)], label=cat.replace("_", " ") if yi == 0 else None)
                left[sc] += v
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["classical 8", "all fixed genes"])
        ax.set_xlabel("share of missing slots (%)")
        ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=2, fontsize=5.5)
        V.panel_letter(ax, "c")
        return save(fig, figdir, "fig_missing_overview")


def fig_missing_features(T, figdir):
    b = T.rd("76", "A_missing_by_feature_bins.tsv")
    if b is None:
        return None
    feats = ["hap_n_genes", "hap_n_contigs", "x_rel_cds", "x_rel_mhc", "x_cds_min_rel", "x_cds_share_lt10"]
    b = b[(b["scope"] == "classical8") & b["feature"].isin(feats)]
    with style():
        fig, axs = plt.subplots(2, 3, figsize=(V.mm(183), V.mm(95)), constrained_layout=True)
        for ax, f, L in zip(axs.ravel(), feats, "abcdef"):
            d = b[(b["feature"] == f) & (b["n"] >= 20)].sort_values("bin_lo")
            if len(d):
                x = (d["bin_lo"] + d["bin_hi"]) / 2
                ax.fill_between(x, d["lo"], d["hi"], color=C_BLUE, alpha=0.2, lw=0)
                ax.plot(x, d["pct"], color=C_BLUE, lw=1, marker=".", ms=2)
            ax.set_xlabel(f[2:].replace("_", " ") if f.startswith("x_") else f.replace("_", " "))
            ax.set_ylabel("slots missing (%)")
            V.panel_letter(ax, L)
        return save(fig, figdir, "fig_missing_by_feature")


def fig_missing_auc(T, figdir):
    a = T.rd("76", "A_missing_feature_auc.tsv")
    if a is None:
        return None
    a = a[a["scope"].isin(["classical8", "all_fixed_genes"])]
    feats = list(a[a["scope"] == "classical8"].sort_values("auc")["feature"])
    with style():
        fig, ax = plt.subplots(figsize=(V.mm(89), V.mm(65)), constrained_layout=True)
        for sc, col, off in (("classical8", C_BLUE, -0.15), ("all_fixed_genes", C_RED, 0.15)):
            d = a[a["scope"] == sc].set_index("feature").reindex(feats).dropna(subset=["auc"])
            dotci(ax, [feats.index(i) + off for i in d.index], d["auc"], d["lo"], d["hi"], col, sc.replace("_", " "))
        ax.axvline(0.5, color=C_GREY, lw=0.6, ls="--")
        ax.set_yticks(range(len(feats)))
        ax.set_yticklabels([f.replace("_", " ") for f in feats])
        ax.set_xlabel("AUC for a missing slot")
        ax.legend(frameon=False, loc="lower right")
        return save(fig, figdir, "fig_missing_auc")


def fig_parent_status(T, figdir):
    ps = T.rd("71f", "trio_parent_status.tsv")
    if ps is None:
        return None
    ps = ps[ps["child_call_named"].astype(str).str.lower() == "true"]
    with style():
        fig, ax = plt.subplots(figsize=(V.mm(89), V.mm(65)), constrained_layout=True)
        for k, (sc, col) in enumerate((("CLASSICAL8", C_BLUE), ("ALL", C_RED))):
            d = ps[ps["scope"] == sc].sort_values("parent_calls_at_gene")
            x = np.arange(len(d)) + (k - 0.5) * 0.38
            ax.bar(x, d["missing_pct"], width=0.36, color=col, label=sc.lower())
            for xi, (p, n, kk) in zip(x, zip(d["missing_pct"], d["child_calls"], d["child_allele_missing_from_parent"])):
                ax.text(xi, p + 1, "%d/%d" % (kk, n), ha="center", fontsize=5.2)
        ticks = sorted(ps["parent_calls_at_gene"].unique())
        ax.set_xticks(range(len(ticks)))
        ax.set_xticklabels([("%d call(s)" % t) for t in ticks])
        ax.set_xlabel("transmitting parent's calls at the gene")
        ax.set_ylabel("child allele absent from parent (%)")
        ax.set_ylim(0, 112)
        ax.legend(frameon=False, loc="upper right")
        return save(fig, figdir, "fig_parent_status")


def fig_status_by_gene(T, figdir):
    t = T.rd("71f", "call_status_taxonomy.tsv")
    if t is None:
        return None
    d = t[t["scope"].isin(["gene", "CLASSICAL8"]) | (t["gene"] == "ALL")]
    genes = [g for g in CLASSICAL if g in set(d["gene"])]
    if not genes:
        return None
    cats = ["valid + valid", "valid + no_gene_row", "no_gene_row + no_gene_row", "valid + no_haplotype", "no_gene_row + no_haplotype"]
    with style():
        fig, ax = plt.subplots(figsize=(V.mm(120), V.mm(70)), constrained_layout=True)
        left = np.zeros(len(genes))
        for k, c in enumerate(cats):
            v = np.array([float(d[(d["gene"] == g) & (d["status_pair"] == c)]["pct"].sum()) for g in genes])
            ax.barh(genes, v, left=left, color=PALETTE[k], label=c)
            left += v
        ax.set_xlim(0, 100)
        ax.set_xlabel("person-genes (%)")
        ax.invert_yaxis()
        ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=5.5)
        return save(fig, figdir, "fig_call_status_by_gene")


def fig_stream_b(T, figdir):
    cc = T.rd("76", "B_class_counts.tsv")
    if cc is None:
        return None
    d0 = cc[(cc["trio_set"] == "primary") & (cc["origin_subset"] == "default")]
    with style():
        fig, axs = plt.subplots(1, 2, figsize=(V.mm(183), V.mm(70)), gridspec_kw=dict(width_ratios=[1, 1.2]), constrained_layout=True)
        ax = axs[0]
        sub = cc[(cc["trio_set"] == "primary") & (cc["scope"] == "classical8") & (cc["gene"] == "ALL")]
        order = ["default", "strong_origin", "no_switch_suspect", "strong_and_no_switch"]
        sub = sub.set_index("origin_subset").reindex([o for o in order if o in set(sub["origin_subset"])])
        y = np.arange(len(sub))
        dotci(ax, y - 0.12, sub["pct"], sub["lo"], sub["hi"], C_RED, "neither (typing-error candidate)")
        oo = [(100 * r["n_other_parent_only"] / r["n_eligible"], *wilson(r["n_other_parent_only"], r["n_eligible"])) for _, r in sub.iterrows()]
        oo = np.array(oo)
        dotci(ax, y + 0.12, oo[:, 0], oo[:, 1], oo[:, 2], C_ORANGE, "other parent only (origin error)")
        ax.set_yticks(y)
        ax.set_yticklabels([s.replace("_", " ") + " (n=%s)" % fmt_n(n) for s, n in zip(sub.index, sub["n_eligible"])])
        ax.invert_yaxis()
        ax.set_xlim(left=0)
        ax.set_xlabel("share of eligible child calls (%)")
        ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=5.5)
        V.panel_letter(ax, "a")
        ax = axs[1]
        g = d0[(d0["scope"] == "classical8") & (d0["gene"] != "ALL")].sort_values("gene")
        y = np.arange(len(g))
        dotci(ax, y, g["pct"], g["lo"], g["hi"], C_RED)
        ax.set_yticks(y)
        ax.set_yticklabels(["%s (n=%s, %s event%s)" % (a, fmt_n(n), fmt_n(k), "" if k == 1 else "s") for a, n, k in zip(g["gene"], g["n_eligible"], g["n_neither"])])
        ax.invert_yaxis()
        ax.set_xlim(left=0)
        ax.set_xlabel("neither-parent rate (%), Wilson 95% interval")
        V.panel_letter(ax, "b")
        return save(fig, figdir, "fig_stream_b")


def fig_stream_b_support(T, figdir):
    fh = T.rd("76", "B_family_event_hist.tsv")
    om = T.rd("76", "B_origin_margin.tsv")
    if fh is None or om is None:
        return None
    with style():
        fig, axs = plt.subplots(1, 3, figsize=(V.mm(183), V.mm(60)), constrained_layout=True)
        ax = axs[0]
        ax.bar(fh["neither_events_in_family"].astype(int).astype(str), fh["families"], color=C_BLUE)
        ax.set_yscale("log")
        ax.set_xlabel("neither events in one family")
        ax.set_ylabel("families (log)")
        V.panel_letter(ax, "a")
        for ax, feat, L, xl in ((axs[1], "log2_kmer_ratio", "b", "log2 (winner / loser k-mers)"), (axs[2], "winner_kmers", "c", "k-mers for the winning parent")):
            for k, cls in enumerate(("match", "other_parent_only", "neither")):
                d = om[(om["cls"] == cls) & (om["feature"] == feat)].sort_values("bin_lo")
                tot = d["n"].sum()
                if tot > 0:
                    ax.step((d["bin_lo"] + d["bin_hi"]) / 2, 100 * d["n"] / tot, where="mid", color=[C_BLUE, C_ORANGE, C_RED][k], label=cls.replace("_", " "))
            ax.set_xlabel(xl)
            ax.set_ylabel("share of class (%)")
            if feat == "winner_kmers":
                ax.set_xscale("symlog", linthresh=10)
            ax.legend(frameon=False, fontsize=5.5)
            V.panel_letter(ax, L)
        return save(fig, figdir, "fig_stream_b_support")


def fig_stream_c(T, figdir):
    r = T.rd("76", "C_novel_rates_by_gene.tsv")
    bd = T.rd("76", "C_novel_by_depth.tsv")
    if r is None or bd is None:
        return None
    with style():
        fig, axs = plt.subplots(1, 2, figsize=(V.mm(183), V.mm(70)), constrained_layout=True)
        ax = axs[0]
        g = r.sort_values("novel_any_pct")
        y = np.arange(len(g))
        dotci(ax, y - 0.12, g["novel_any_pct"], g["novel_any_lo"], g["novel_any_hi"], C_BLUE, "any novel tag")
        dotci(ax, y + 0.12, g["unnamed_pct"], g["unnamed_lo"], g["unnamed_hi"], C_RED, "novel at field 1-2 (no 2-field name)")
        ax.set_yticks(y)
        ax.set_yticklabels(g["gene"])
        ax.set_xlabel("calls (%)")
        ax.legend(frameon=False, loc="lower right", fontsize=5.5)
        V.panel_letter(ax, "a")
        ax = axs[1]
        for k, o in enumerate(("novel_any", "novel_unnamed")):
            d = bd[(bd["outcome"] == o) & (bd["n"] >= 20)].sort_values("bin_lo")
            x = (d["bin_lo"] + d["bin_hi"]) / 2
            ax.fill_between(x, d["lo"], d["hi"], color=[C_BLUE, C_RED][k], alpha=0.2, lw=0)
            ax.plot(x, d["pct"], color=[C_BLUE, C_RED][k], lw=1, label=o.replace("_", " "))
        ax.set_xlabel("CDS depth relative to the person's median")
        ax.set_ylabel("calls with a novel tag (%)")
        ax.legend(frameon=False, fontsize=5.5)
        V.panel_letter(ax, "b")
        return save(fig, figdir, "fig_stream_c")


# ---- pruning (76b / 78b share their table layout) ----
def top_metrics(auc, view, k=5):
    d = auc[auc["view"] == view].copy()
    d["dev"] = (d["auc"] - 0.5).abs()
    return list(d.sort_values("auc", ascending=False).head(k)["metric"])


def fig_prune(T, sub, figdir, name, views_wanted):
    cur = T.rd(sub, "prune_curves.tsv")
    auc = T.rd(sub, "auc.tsv")
    if cur is None or auc is None:
        return None
    views = [v for v in views_wanted if v in set(cur["view"])]
    if not views:
        return None
    nc = 2 if len(views) > 1 else 1
    nr = int(np.ceil(len(views) / nc))
    with style():
        fig, axs = plt.subplots(nr, nc, figsize=(V.mm(183), V.mm(62 * nr)), squeeze=False, constrained_layout=True)
        for ax, v, L in zip(axs.ravel(), views, "abcdefgh"):
            ms = top_metrics(auc, v, 5)
            first = True
            for k, m in enumerate(ms):
                d = cur[(cur["view"] == v) & (cur["metric"] == m)].sort_values("prune_pct")
                if d.empty:
                    continue
                K = d["n_errors"].iloc[0]
                if first and K > 0:
                    ax.fill_between(d["prune_pct"], 100 * d["rand_removed_lo"].fillna(0) / K, 100 * d["rand_removed_hi"].fillna(0) / K, color=C_GREY, alpha=0.25, lw=0)
                    first = False
                ax.plot(d["prune_pct"], d["captured_pct"], color=PALETTE[k], lw=1, label=m.replace("_", " "))
            mx = float(cur[cur["view"] == v]["prune_pct"].max())
            ax.plot([0, mx], [0, mx], color=C_GREY, lw=0.6, ls="--")
            e = cur[cur["view"] == v]["n_errors"]
            ax.set_title("%s (%s events)" % (v.replace("_", " "), fmt_n(e.iloc[0]) if len(e) else "?"), fontsize=6.5)
            ax.set_xlabel("worst x% of units removed")
            ax.set_ylabel("errors removed (%)")
            ax.legend(frameon=False, fontsize=5, loc="lower right")
            V.panel_letter(ax, L)
        for ax in axs.ravel()[len(views):]:
            ax.axis("off")
        return save(fig, figdir, name)


def fig_prune_residual(T, sub, figdir, name, views_wanted):
    cur = T.rd(sub, "prune_curves.tsv")
    auc = T.rd(sub, "auc.tsv")
    if cur is None or auc is None or "residual_error_pct" not in cur:
        return None
    views = [v for v in views_wanted if v in set(cur["view"])]
    if not views:
        return None
    nc = 2 if len(views) > 1 else 1
    nr = int(np.ceil(len(views) / nc))
    with style():
        fig, axs = plt.subplots(nr, nc, figsize=(V.mm(183), V.mm(62 * nr)), squeeze=False, constrained_layout=True)
        for ax, v, L in zip(axs.ravel(), views, "abcdefgh"):
            ms = top_metrics(auc, v, 1)
            if not ms:
                continue
            d = cur[(cur["view"] == v) & (cur["metric"] == ms[0])].sort_values("prune_pct")
            if "residual_lo" in d:
                ax.fill_between(d["prune_pct"], d["residual_lo"], d["residual_hi"], color=C_BLUE, alpha=0.2, lw=0)
            ax.plot(d["prune_pct"], d["residual_error_pct"], color=C_BLUE, lw=1, label="residual error, " + ms[0].replace("_", " "))
            ax.axhline(d["base_error_pct"].iloc[0], color=C_RED, lw=0.7, ls="--", label="no pruning")
            ax.set_title(v.replace("_", " "), fontsize=6.5)
            ax.set_xlabel("worst x% of units removed")
            ax.set_ylabel("error rate among the rest (%)")
            ax.legend(frameon=False, fontsize=5.5)
            V.panel_letter(ax, L)
        for ax in axs.ravel()[len(views):]:
            ax.axis("off")
        return save(fig, figdir, name)


def fig_auc_views(T, sub, figdir, name):
    auc = T.rd(sub, "auc.tsv")
    if auc is None:
        return None
    views = list(dict.fromkeys(auc["view"]))
    nc = 2 if len(views) > 1 else 1
    nr = int(np.ceil(len(views) / nc))
    with style():
        fig, axs = plt.subplots(nr, nc, figsize=(V.mm(183), V.mm(max(55, 8 + 3.2 * len(auc["metric"].unique())) * nr)), squeeze=False, constrained_layout=True)
        for ax, v, L in zip(axs.ravel(), views, "abcdefgh"):
            d = auc[auc["view"] == v].dropna(subset=["auc"]).sort_values("auc")
            y = np.arange(len(d))
            dotci(ax, y, d["auc"], d["auc_lo"].fillna(d["auc"]), d["auc_hi"].fillna(d["auc"]), C_BLUE)
            ax.axvline(0.5, color=C_GREY, lw=0.6, ls="--")
            ax.set_yticks(y)
            ax.set_yticklabels([m.replace("_", " ") for m in d["metric"]], fontsize=5.5)
            ax.set_title("%s (%s events)" % (v.replace("_", " "), fmt_n(d["n_errors"].max() if len(d) else 0)), fontsize=6.5)
            ax.set_xlabel("AUC")
            V.panel_letter(ax, L)
        for ax in axs.ravel()[len(views):]:
            ax.axis("off")
        return save(fig, figdir, name)


def fig_cv(T, sub, figdir, name):
    cv = T.rd(sub, "cv_by_metric.tsv")
    if cv is None or cv.empty:
        return None
    views = list(dict.fromkeys(cv["view"]))
    with style():
        fig, axs = plt.subplots(1, len(views), figsize=(V.mm(183), V.mm(80)), squeeze=False, constrained_layout=True)
        for ax, v, L in zip(axs.ravel(), views, "abcdefgh"):
            d = cv[cv["view"] == v].sort_values("enrichment")
            y = np.arange(len(d))
            dotci(ax, y, d["enrichment"], d["enrichment_lo"], d["enrichment_hi"], C_BLUE)
            ax.axvline(1, color=C_GREY, lw=0.6, ls="--")
            ax.set_yticks(y)
            ax.set_yticklabels([m.replace("_", " ") for m in d["metric"]], fontsize=5)
            ax.set_title(v.replace("_", " "), fontsize=6.5)
            ax.set_xlabel("held-out enrichment")
            V.panel_letter(ax, L)
        return save(fig, figdir, name)


# ---- KIR ----
def fig_kir_structure(T, figdir):
    hc = T.rd("77", "hap_class_counts.tsv")
    nc = T.rd("77", "gene_n_called.tsv")
    nt = T.rd("77", "gene_by_novelty_tier.tsv")
    if hc is None or nc is None or nt is None:
        return None
    with style():
        fig, axs = plt.subplots(1, 3, figsize=(V.mm(183), V.mm(80)), gridspec_kw=dict(width_ratios=[0.8, 1.2, 1.2]), constrained_layout=True)
        ax = axs[0]
        h = hc.sort_values("n_haps")
        ax.barh(h["hap_class"].str.replace("_", " "), h["n_haps"], color=C_BLUE)
        ax.set_xlabel("haplotypes")
        V.panel_letter(ax, "a")
        ax = axs[1]
        genes = sorted(nc["gene"].unique())
        one = nc[nc["n_called"] == 1].set_index("gene")["n_people"].reindex(genes).fillna(0)
        two = nc[nc["n_called"] == 2].set_index("gene")["n_people"].reindex(genes).fillna(0)
        tot = (one + two).replace(0, np.nan)
        ax.barh(genes, 100 * two / tot, color=C_BLUE, label="2 calls")
        ax.barh(genes, 100 * one / tot, left=100 * two / tot, color=C_ORANGE, label="1 call")
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels([g.replace("KIR", "") for g in genes], fontsize=5.5)
        ax.invert_yaxis()
        ax.set_xlabel("people with the gene called (%)")
        ax.legend(frameon=False, fontsize=5.5, loc="lower right")
        V.panel_letter(ax, "b")
        ax = axs[2]
        piv = nt.pivot_table(index="gene", columns="novelty_tier", values="n", aggfunc="sum").reindex(genes).fillna(0)
        piv = piv.div(piv.sum(axis=1).replace(0, np.nan), axis=0) * 100
        left = np.zeros(len(piv))
        for k, c in enumerate(piv.columns):
            ax.barh(range(len(piv)), piv[c].fillna(0), left=left, color=PALETTE[k % len(PALETTE)], label=c.replace("_", " "))
            left += piv[c].fillna(0).to_numpy()
        ax.set_yticks(range(len(piv)))
        ax.set_yticklabels([g.replace("KIR", "") for g in piv.index], fontsize=5.5)
        ax.invert_yaxis()
        ax.set_xlabel("calls (%)")
        ax.legend(frameon=False, fontsize=4.8, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2)
        V.panel_letter(ax, "c")
        return save(fig, figdir, "fig_kir_structure")


def fig_kir_origin(T, figdir):
    g = T.rd("77b", "gene_origin_summary_kir.tsv")
    k = T.rd("77b", "kmers_per_window_kir.tsv")
    if g is None or k is None:
        return None
    g = g[g["gene"] != "ALL"].sort_values("gene")
    with style():
        fig, axs = plt.subplots(1, 2, figsize=(V.mm(183), V.mm(75)), constrained_layout=True)
        ax = axs[0]
        y = np.arange(len(g))
        ax.barh(y - 0.2, g["ambiguous_pct"].fillna(0), height=0.38, color=C_RED, label="ambiguous windows")
        ax.barh(y + 0.2, g["strong_of_resolved_pct"].fillna(0), height=0.38, color=C_GREEN, label="strong among resolved")
        ax.set_yticks(y)
        ax.set_yticklabels([x.replace("KIR", "") for x in g["gene"]], fontsize=5.5)
        ax.invert_yaxis()
        ax.set_xlabel("windows (%); blank bars are masked (<20)")
        ax.legend(frameon=False, fontsize=5.5, loc="lower right")
        V.panel_letter(ax, "a")
        ax = axs[1]
        order = list(dict.fromkeys(k["bin"]))
        tot = k.groupby("bin")["n"].sum().reindex(order)
        ax.bar(range(len(order)), tot, color=C_BLUE)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, fontsize=5.5, rotation=45, ha="right")
        ax.set_xlabel("parent-specific k-mers in the window")
        ax.set_ylabel("windows (unmasked genes only)")
        V.panel_letter(ax, "b")
        return save(fig, figdir, "fig_kir_origin")


def fig_kir_trio(T, figdir):
    tg = T.rd("77c", "trio_classes_group.tsv")
    pg = T.rd("77c", "pair_classes_group.tsv")
    if tg is None:
        return None
    with style():
        fig, axs = plt.subplots(1, 2, figsize=(V.mm(183), V.mm(75)), constrained_layout=True)
        ax = axs[0]
        d = tg[tg["gene"] == "ALL"]
        order = ["all", "tx_named", "strong_origin", "complete_only"]
        d = d.set_index("variant").reindex([o for o in order if o in set(d["variant"])])
        y = np.arange(len(d))
        ok = d["neither_pct"].notna()
        dotci(ax, y[ok.to_numpy()], d["neither_pct"][ok], d["lo"][ok], d["hi"][ok], C_RED)
        ax.set_yticks(y)
        ax.set_yticklabels([("%s (n=%s)" % (v.replace("_", " "), fmt_n(n))) for v, n in zip(d.index, d["n_informative"])])
        ax.invert_yaxis()
        ax.set_xlim(left=0)
        ax.set_xlabel("neither-parent rate (%), protein level")
        V.panel_letter(ax, "a")
        ax = axs[1]
        if pg is not None:
            q = pg[(pg["kind"] == "pair") & (pg["variant"].isin(["all", "complete_only"])) & (pg["locus"] != "ALL")]
            loci = list(q[q["variant"] == "all"].sort_values("no_share_pct")["locus"])
            for k, (v, col, off) in enumerate((("all", C_BLUE, -0.15), ("complete_only", C_GREEN, 0.15))):
                dd = q[q["variant"] == v].set_index("locus").reindex(loci).dropna(subset=["no_share_pct"])
                dotci(ax, [loci.index(i) + off for i in dd.index], dd["no_share_pct"], dd["lo"], dd["hi"], col, v.replace("_", " "))
            ax.set_yticks(range(len(loci)))
            ax.set_yticklabels([x.replace("KIR", "") for x in loci], fontsize=5.5)
            ax.set_xlabel("first-degree pairs sharing no label (%)")
            ax.legend(frameon=False, fontsize=5.5, loc="lower right")
        V.panel_letter(ax, "b")
        return save(fig, figdir, "fig_kir_trio_pair")


# ---------------------------------------------------------------------------
# computed numbers for the prose
# ---------------------------------------------------------------------------
class N:
    """All numbers the prose uses, computed once from the tables (None when the table is absent)."""

    def __init__(self, T):
        self.T = T
        g = T.rd("76", "A_missing_by_gene.tsv")
        self.miss_pct = self.miss_n = self.miss_k = None
        if g is not None:
            d = g[g["scope"] == "classical8"]
            if len(d) and d["n"].sum() > 0:
                self.miss_k, self.miss_n = d["k"].sum(), d["n"].sum()
                self.miss_pct = 100 * self.miss_k / self.miss_n
        ex = T.rd("76", "A_missing_explained.tsv")
        self.break_pct = self.end_pct = self.inside_pct = None
        if ex is not None:
            e = ex[ex["scope"] == "classical8"].set_index("explanation")["pct_of_missing"]
            self.break_pct = e.get("assembly_break_at_gene")
            self.end_pct = e.get("region_end")
            self.inside_pct = (e.get("inside_contig_depth_ok", 0) or 0) + (e.get("inside_contig_low_depth", 0) or 0) if len(e) else None
        hc = T.rd("76", "A_missing_by_hap_class.tsv")
        self.inside_rate = self.between_rate = None
        if hc is not None:
            h = hc[hc["scope"] == "classical8"].set_index("hap_class")["pct"]
            self.inside_rate, self.between_rate = h.get("inside_contig"), h.get("between_contigs")
        cc = T.rd("76", "B_class_counts.tsv")
        self.b = None
        if cc is not None:
            d = cc[(cc["trio_set"] == "primary") & (cc["scope"] == "classical8") & (cc["origin_subset"] == "default") & (cc["gene"] == "ALL")]
            if len(d):
                self.b = d.iloc[0]
        ps = T.rd("71f", "trio_parent_status.tsv")
        self.p1 = self.p2 = None
        if ps is not None:
            d = ps[(ps["scope"] == "CLASSICAL8") & (ps["child_call_named"].astype(str).str.lower() == "true")].set_index("parent_calls_at_gene")
            if 1 in d.index and 2 in d.index:
                self.p1 = (d.loc[1, "child_allele_missing_from_parent"], d.loc[1, "child_calls"], d.loc[1, "missing_pct"])
                self.p2 = (d.loc[2, "child_allele_missing_from_parent"], d.loc[2, "child_calls"], d.loc[2, "missing_pct"])

    def sentence_missing(self):
        if self.miss_pct is None:
            return pending("missing-call rate", "76/A_missing_by_gene.tsv")
        s = "Across the classical 8 genes, %s of %s person-haplotype slots (%s) carry no gene annotation." % (fmt_n(self.miss_k), fmt_n(self.miss_n), fmt_p(self.miss_pct, 2))
        if self.break_pct is not None:
            s += " Of the missing slots, %s sit where the haplotype assembly is broken between contigs at the gene and %s where the gene lies beyond the end of the assembled region; only %s are inside a contig." % (
                fmt_p(self.break_pct), fmt_p(self.end_pct), fmt_p(self.inside_pct))
        if self.inside_rate is not None and self.between_rate is not None:
            s += " A gene inside an unbroken contig is missing in %s of slots, one between contigs in %s." % (fmt_p(self.inside_rate, 3), fmt_p(self.between_rate))
        return "<p>%s</p>" % s

    def sentence_b(self):
        if self.b is None:
            return pending("wrong-allele rate among fully called trios", "76/B_class_counts.tsv")
        b = self.b
        return ("<p>Among %s trio child calls with a fully called genotype, a fully called transmitting parent and a resolved parent of origin, "
                "%s carry an allele found in neither parent (%s, Wilson 95%% %s-%s) and %s carry one found only in the other parent (an origin error, not a typing error).</p>" % (
                    fmt_n(b["n_eligible"]), fmt_n(b["n_neither"]), fmt_p(b["pct"], 2), fmt_p(b["lo"], 2), fmt_p(b["hi"], 2), fmt_n(b["n_other_parent_only"])))

    def sentence_parent(self):
        if self.p1 is None:
            return pending("trio mismatch by parent completeness", "71f/trio_parent_status.tsv")
        return ("<p>The same child calls look very different depending on how many calls the transmitting parent has at the gene: with one call the child's allele is absent "
                "from the parent in %s/%s records (%s); with two calls in %s/%s (%s). A parent with one call simply cannot show the allele the child received from the "
                "missing haplotype, so the old pooled mismatch rate measured missing calls, not wrong alleles.</p>" % (
                    fmt_n(self.p1[0]), fmt_n(self.p1[1]), fmt_p(self.p1[2]), fmt_n(self.p2[0]), fmt_n(self.p2[1]), fmt_p(self.p2[2])))


# ---------------------------------------------------------------------------
# pruning decision table
# ---------------------------------------------------------------------------
def xopt_at_cost(d, cost):
    """Largest prune % before the windowed marginal NNP exceeds `cost` (the 76b rule applied to the stored curve)."""
    d = d.sort_values("prune_pct")
    m = d["marginal_nnp"].to_numpy(float)
    x = d["prune_pct"].to_numpy(float)
    best = 0.0
    for i in range(len(d)):
        if np.isfinite(m[i]) and m[i] > cost:
            return best
        if np.isfinite(m[i]) or i == 0:
            best = float(x[i])
    return best


def evidence_label(p5, share, n_err, auc, auc_min):
    if n_err is None or not np.isfinite(n_err) or n_err < 5:
        return "none (fewer than 5 events)"
    if auc_min is not None and np.isfinite(auc_min) and np.isfinite(auc) and auc < auc_min and n_err < 30:
        return "underpowered (AUC below the minimal detectable %.2f)" % auc_min
    if np.isfinite(p5) and p5 < 0.001 and np.isfinite(share) and share >= 0.9:
        return "strong"
    if np.isfinite(p5) and p5 < 0.05 and np.isfinite(share) and share >= 0.7:
        return "moderate"
    return "weak / not distinguishable from random"


def decision_table(T, sub, views, top=3):
    cur, auc, summ = T.rd(sub, "prune_curves.tsv"), T.rd(sub, "auc.tsv"), T.rd(sub, "prune_summary.tsv")
    if cur is None or auc is None:
        return None
    cvm = T.rd(sub, "cv_by_metric.tsv")
    pw = T.rd(sub, "power.tsv")
    rows = []
    for v in views:
        if v not in set(auc["view"]):
            continue
        auc_min = np.nan
        if pw is not None and v in set(pw["view"]) and "auc_min_detectable" in pw:
            auc_min = float(pw[pw["view"] == v]["auc_min_detectable"].iloc[0])
        for m in top_metrics(auc, v, top):
            c = cur[(cur["view"] == v) & (cur["metric"] == m)]
            if c.empty:
                continue
            a = auc[(auc["view"] == v) & (auc["metric"] == m)].iloc[0]
            r = dict(view=v, metric=m, AUC=a["auc"], events=a["n_errors"])
            for cost in COST_RATIOS:
                x = xopt_at_cost(c, cost)
                row = c.iloc[(c["prune_pct"] - x).abs().argmin()]
                r["x%% @%d" % cost] = x
                r["captured%% @%d" % cost] = row["captured_pct"]
            p5 = np.nan
            if summ is not None:
                s = summ[(summ["view"] == v) & (summ["metric"] == m)]
                if len(s) and "p_tail_at_5" in s:
                    p5 = float(s["p_tail_at_5"].iloc[0])
            share = np.nan
            if cvm is not None:
                q = cvm[(cvm["view"] == v) & (cvm["metric"] == m)]
                if len(q):
                    share = float(q["share_reps_enrichment_gt1"].iloc[0])
                    r["held-out enrichment"] = float(q["enrichment"].iloc[0])
            r["p(random>=obs) @5%"] = p5
            r["evidence"] = evidence_label(p5, share, a["n_errors"], a["auc"], auc_min)
            rows.append(r)
    return pd.DataFrame(rows) if rows else None


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
TABLE_DOC = {
    "76/A_missing_by_gene.tsv": "stream A: share of person-haplotype slots with no gene annotation, per gene (Wilson 95%)",
    "76/A_missing_by_hap_class.tsv": "stream A: missing rate by where the gene lies relative to the assembly contigs",
    "76/A_missing_explained.tsv": "stream A: how the missing slots break down (assembly break, region end, inside contig)",
    "76/A_missing_feature_auc.tsv": "stream A: AUC of haplotype features for a missing slot",
    "76/A_missing_by_feature_bins.tsv": "stream A: missing rate in fine bins of each feature",
    "76/B_class_counts.tsv": "stream B: match / other-parent-only / neither counts per gene, scope, origin subset",
    "76/B_family_clustering.tsv": "stream B: events per family (effective independent events)",
    "76/B_family_event_hist.tsv": "stream B: histogram of neither events per family",
    "76/B_origin_margin.tsv": "stream B: yak k-mer margin distribution by class",
    "76/B_by_child_homozygosity.tsv": "stream B: events by child homozygosity at the gene",
    "76/B_status_counts.tsv": "stream B: why child calls were or were not labelled",
    "76/B_name_vs_cds.tsv": "stream B: 2-field name equal but CDS differs",
    "76/C_novel_rates_by_gene.tsv": "stream C: novel-call rates per gene",
    "76/C_novel_by_depth.tsv": "stream C: novel-call rate by CDS depth bin",
    "76/C_novel_vs_incomplete.tsv": "stream C: novel calls versus incomplete genotypes",
    "76/C_twin_concordance.tsv": "stream C: twin concordance with and without a novel call",
    "76/definition_card.tsv": "definitions and denominators of every label used",
    "71f/trio_parent_status.tsv": "child allele absent from the transmitting parent, by the parent's number of calls",
    "71f/call_status_taxonomy.tsv": "what each person-gene's two haplotype slots look like (valid / no gene row / no haplotype)",
    "71f/trio_candidates.tsv": "children with more than one candidate trio",
    "77/hap_class_counts.tsv": "KIR haplotypes by completeness class",
    "77/gene_n_called.tsv": "KIR people with 1 or 2 calls per gene",
    "77/gene_by_novelty_tier.tsv": "KIR calls per novelty tier per gene",
    "77b/gene_origin_summary_kir.tsv": "KIR yak origin: windows, ambiguity, strong origin, k-mers per window",
    "77b/kmers_per_window_kir.tsv": "KIR yak origin: k-mers per window histogram",
    "77c/trio_classes_group.tsv": "KIR trio classes at paralog-group level, per variant",
    "77c/trio_classes_gene.tsv": "KIR trio classes per gene, per variant",
    "77c/pair_classes_group.tsv": "KIR first-degree pairs: share of pairs with no shared label",
    "76b/prune_curves.tsv": "HLA pruning: error removal along the worst-first ordering (0.25% steps)",
    "76b/prune_summary.tsv": "HLA pruning: x_opt and fixed-point summaries with bootstrap spread",
    "76b/auc.tsv": "HLA pruning: AUC per view and metric",
    "76b/cv_by_metric.tsv": "HLA pruning: family-grouped cross-validation, thresholds learned on training families",
    "76b/cv_selected.tsv": "HLA pruning: cross-validation with the best training metric chosen per split",
    "76b/transfer.tsv": "HLA pruning: thresholds learned in one view applied to another",
    "76b/power.tsv": "HLA pruning: events and minimal detectable AUC per view",
    "76b/threshold_stability.tsv": "HLA pruning: spread of the learned thresholds across folds",
    "78/depth_coverage_by_gene.tsv": "KIR depth: which genes have a depth metric",
    "78b/prune_curves.tsv": "KIR pruning: error removal curves", "78b/auc.tsv": "KIR pruning: AUC per view and metric",
    "78b/per_stratum_auc.tsv": "KIR pruning: AUC by gene, paralog group, gene class, haplotype content",
}

GLOSSARY = [
    ("Call", "one gene annotation on one haplotype assembly produced by Immuannot; carries allele name, distances, warning flag, ties and a novel tag."),
    ("No call", "no gene annotation on that haplotype: Immuannot found no template passing its filter, or the sequence is not in the window."),
    ("Slot", "person x gene x haplotype (two per person-gene)."),
    ("Fully called", "both haplotype slots of a person-gene carry a call."),
    ("Origin (pat/mat)", "the parent whose parent-specific k-mers (yak trioeval) dominate the child's haplotype window at the gene (gene +-10 kb, >=5 k-mers and >=3x the other); the labels are the two parent ids, not sexes."),
    ("MATCH / OTHER ONLY / NEITHER", "child allele present in the transmitting parent / only in the other parent (origin error or recombination) / in neither (typing-error candidate)."),
    ("Pruning", "removing the worst x% of units ordered by a metric; 'worst' direction is fixed per metric (oriented columns o_*, higher = worse)."),
    ("Enrichment", "share of errors removed divided by share of units removed (1 = random)."),
    ("NNP", "calls pruned per error removed; marginal NNP is the same over the last 2 percentage points."),
    ("Cost ratio", "the marginal NNP above which pruning stops being worth it: a choice, so four values are shown."),
    ("Contig co-location", "two genes on one assembled contig; never called 'physical phasing' in this project."),
]


def build(T, figs, out_dir):
    n = N(T)
    parts = []
    A = parts.append
    A("<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
      "<title>S06 call quality final report</title><style>%s</style><script>%s</script></head><body>" % (CSS, JS))
    A("<nav>" + "".join('<a href="#%s">%s</a>' % (i, t) for i, t in (("sum", "Summary"), ("mech", "Mechanics"), ("hla", "HLA results"), ("prune", "Pruning"), ("kir", "KIR"),
                                                                      ("raw", "Raw evidence"), ("next", "Open questions"), ("think-tank-section", "Think tank"))) + "</nav><main>")
    A("<h1>Call-quality metrics and pruning: HLA and KIR</h1><p class='sub'>Sprint S06. Built from aggregate tables only; every number below is computed from them. "
      "Cells below 20 are masked upstream.</p>")

    # (0) summary
    A("<h2 id='sum'>0. Summary</h2><div class='box'><b>Question.</b> Can measurable per-call quality metrics (Immuannot distances and flags, exon / CDS depth, assembly "
      "context) replace the subjective 0-7 score, and if we remove the worst x% of calls by a metric, which percent and which metric remove the most error for the least loss, "
      "in a way that holds on other families and other kinds of relatives?</div>")
    A("<h3>What we found</h3>")
    A("<ol><li><b>Most of what looked like typing error is missing calls.</b></li></ol>")
    A(n.sentence_missing())
    A(n.sentence_parent())
    A(n.sentence_b())
    A("<p>So there are two different problems that were mixed in the old label: a call that is absent (common, structural, predicted by assembly context) and a call that is present "
      "but wrong (rare, tens of events, so only rates and pre-specified tail tests are defensible).</p>")
    summ_dec = None
    if T.rd("76b", "prune_curves.tsv") is not None:
        summ_dec = decision_table(T, "76b", ["A_classical8", "B_child", "B_worst", "pair_complete"], top=1)
    if summ_dec is not None and len(summ_dec):
        A("<h3>What to choose (HLA)</h3><p>Best single metric per view with the percent where marginal cost exceeds the stated cost ratio; the full table, all cost ratios and the raw curves are in section 3.</p>")
        A(table_html(summ_dec, fmt={c: (lambda v: "" if (isinstance(v, float) and np.isnan(v)) else (("%.3g" % v) if isinstance(v, (float, np.floating)) else v)) for c in summ_dec.columns}))
    else:
        A(pending("the pruning recommendation (percent and metric, with evidence strength)", "76b/prune_curves.tsv, auc.tsv, prune_summary.tsv, cv_by_metric.tsv"))
    A("<h3>Caveats</h3><ul><li>Parent of origin is decided per gene call (gene +-10 kb), not per haplotype; windows of close genes overlap, so neighbouring genes share evidence. "
      "Switches along a contig are allowed but not checked at haplotype level; the 'no switch suspect' subset is the check.</li>"
      "<li>The trio set is small (about 170 children); events are clustered by family; power is stated next to every curve.</li>"
      "<li>Depth is measured on GRCh38-aligned reads, reference-biased at divergent exons (DRB1 exon 2 about 0.65 of other exons); KIR genes absent from GRCh38 have no depth metric.</li></ul>")

    # (1) mechanics
    A("<h2 id='mech'>1. Mechanics in plain language</h2>")
    A("<h3>1.1 What a call is, and what no call means</h3><p>Immuannot runs once per haplotype assembly on a trimmed MHC window (KIR: chr19 window). For every gene it picks the best "
      "template allele; the gene is annotated only if that template passes the filter (adjusted mismatches &le;3, or mismatches/matches &lt;0.03). Otherwise the gene is silently absent. "
      "A 'no call' is therefore not a measurement of absence: it is the caller declining, usually because the assembly breaks near the gene or the region ends.</p>")
    A(schematic_call())
    A("<h3>1.2 How the trio label works</h3><p>Each child haplotype call gets its own yak origin: parent-unique k-mers are counted on the child's contig over the gene +-10 kb; a parent wins with "
      "&ge;5 k-mers and &ge;3x the other, else the call is ambiguous and is left out. The child's 2-field allele name is then looked up among the transmitting parent's calls.</p>")
    A(schematic_trio())
    A("<h3>1.3 Why pairs and trios behave differently</h3>")
    A(table_html(pd.DataFrame([
        ["unit", "child haplotype call with a resolved parent of origin", "person pair x gene"],
        ["error means", "child allele absent from the transmitting parent", "no shared 2-field allele between the two people"],
        ["a missing call in...", "...the parent: scored as an error; the child's own missing call: no record", "...either person: scored as an error"],
        ["biology scored as error", "rarely", "DRB3/4/5 presence differences, 33 of the earlier 142 pair errors"],
        ["therefore", "the pruning metric of the child's own call barely matters; the parent's completeness does", "any missing call anywhere in the pair matters, so n_called predicts it well"],
    ], columns=["", "trios", "first-degree pairs"])))
    A("<h3>1.4 The three outcome streams</h3>")
    A(schematic_streams())
    A("<p class='read'>Definitions and denominators (data-driven definition card):</p>")
    dc = T.rd("76", "definition_card.tsv")
    A(table_html(dc, max_rows=60) if dc is not None else pending("definition card", "76/definition_card.tsv"))

    # (2) HLA results
    A("<h2 id='hla'>2. HLA results</h2><h3>2.1 Stream A: the missing call</h3>")
    A(n.sentence_missing())
    A(img_tag(figs.get("missing"), "missing call overview", "a, rate per gene (Wilson 95%); b, rate by where the gene lies relative to the contigs; c, what the missing slots are.") if figs.get("missing")
      else pending("stream A overview figure", "76/A_missing_by_gene.tsv, A_missing_by_hap_class.tsv, A_missing_explained.tsv"))
    A(img_tag(figs.get("missing_features"), "missing by feature", "Missing rate in fine bins of haplotype features; bins with fewer than 20 slots are not drawn.") if figs.get("missing_features")
      else pending("missing rate by feature", "76/A_missing_by_feature_bins.tsv"))
    A(img_tag(figs.get("missing_auc"), "missing AUC", "AUC of each haplotype feature for a missing slot, with bootstrap intervals.") if figs.get("missing_auc") else "")
    A(mind("If missing slots inside unbroken contigs turned out to be common (they are %s here), or if short-read data showed the gene is present on those haplotypes (the earlier HLA-A check found most "
           "'deletions' contradicted by short reads), the missing call would be a calling failure rather than assembly structure, and the cohort-wide fix would be recalling rather than pruning." % fmt_p(n.inside_rate, 3)))
    A("<h3>2.2 Stream B: the wrong allele, among fully called genotypes</h3>")
    A(n.sentence_b())
    A(img_tag(figs.get("stream_b"), "stream B", "a, neither and other-parent-only rates under four origin rules; b, per-gene neither rate with the number of events.") if figs.get("stream_b")
      else pending("stream B figure", "76/B_class_counts.tsv"))
    A(img_tag(figs.get("stream_b_support"), "stream B support", "a, events per family; b-c, yak k-mer margins by class.") if figs.get("stream_b_support") else "")
    for nm, why in (("B_family_clustering.tsv", "events per family"), ("B_by_child_homozygosity.tsv", "events by child homozygosity"), ("B_status_counts.tsv", "why calls were not labelled"),
                    ("B_name_vs_cds.tsv", "name versus CDS")):
        d = T.rd("76", nm)
        if d is not None:
            A("<p class='read'>%s</p>%s" % (esc(why), table_html(d, max_rows=30)))
    A(mind("A neither event can still be a yak origin error with the other parent also carrying a similar allele, or a recombination inside the +-10 kb window; the 'no switch suspect' and 'strong origin' rows "
           "show how much the count moves. If several events came from one family, the effective number of independent events is lower than the count (see family table)."))
    A("<h3>2.3 Stream C: novel calls</h3>")
    A(img_tag(figs.get("stream_c"), "stream C", "a, novel-call rates per gene; b, novel-call rate by CDS depth.") if figs.get("stream_c") else pending("stream C figure", "76/C_novel_rates_by_gene.tsv, C_novel_by_depth.tsv"))
    for nm, why in (("C_novel_vs_incomplete.tsv", "novel unnamed calls versus incomplete genotypes"), ("C_twin_concordance.tsv", "twin concordance")):
        d = T.rd("76", nm)
        if d is not None:
            A("<p class='read'>%s</p>%s" % (esc(why), table_html(d)))
    A("<h3>2.4 Why trios and pairs behave so differently, reconciled</h3>")
    A(n.sentence_parent())
    A(img_tag(figs.get("parent_status"), "parent status", "Child allele absent from the transmitting parent, by the parent's number of calls at the gene (labels: events/records).") if figs.get("parent_status")
      else pending("parent-completeness figure", "71f/trio_parent_status.tsv"))
    A(img_tag(figs.get("status_gene"), "call status by gene", "What the two haplotype slots of every person-gene look like, per gene.") if figs.get("status_gene") else "")

    # (3) pruning
    A("<h2 id='prune'>3. Pruning evidence (HLA)</h2><p>Nothing is chosen for you: curves, bands, cost-ratio sensitivity and held-out results are all shown. The gray band is where random pruning of the same size would land (exact 95%).</p>")
    if T.rd("76b", "prune_curves.tsv") is None:
        A(pending("all HLA pruning evidence", "76b/prune_curves.tsv, prune_summary.tsv, auc.tsv, cv_by_metric.tsv, cv_selected.tsv, transfer.tsv, power.tsv, threshold_stability.tsv"))
    else:
        views = ["A_classical8", "A_all", "B_child", "B_worst", "B_child_all", "combined_child", "pair_complete", "pair_all"]
        A("<h3>3.1 Errors removed against units removed</h3>")
        A(img_tag(figs.get("prune76"), "pruning curves", "Top five metrics per view by AUC; dashed = random expectation, gray = random 95% band."))
        A(img_tag(figs.get("prune76_res"), "residual error", "Error rate among the units that remain after pruning by the best metric (bootstrap 95% band)."))
        A("<h3>3.2 AUC per metric</h3>" + img_tag(figs.get("auc76"), "AUC", "AUC with cluster bootstrap intervals; 0.5 = no information."))
        A("<h3>3.3 Decision table and cost-ratio sensitivity</h3><p class='read'>x = largest percent pruned before the marginal calls-pruned-per-error-removed exceeds the cost ratio (10, 25, 50, 100). "
          "Evidence: strong needs p(random) &lt;0.001 at 5% pruning and enrichment &gt;1 in &ge;90% of held-out repeats; moderate p&lt;0.05 and &ge;70%.</p>")
        dt = decision_table(T, "76b", views, top=3)
        A(table_html(dt, searchable="dectab") if dt is not None else pending("decision table", "76b/prune_curves.tsv + auc.tsv"))
        A("<h3>3.4 Do the thresholds transfer?</h3>" + (img_tag(figs.get("cv76"), "cross-validation", "Family-grouped held-out enrichment per metric (thresholds learned on training families).") if figs.get("cv76") else ""))
        for nm, why in (("cv_selected.tsv", "metric chosen on training families, scored on held-out families"), ("transfer.tsv", "threshold learned in one view, applied to another"),
                        ("power.tsv", "events and minimal detectable AUC"), ("threshold_stability.tsv", "threshold spread across folds")):
            d = T.rd("76b", nm)
            A("<p class='read'>%s</p>%s" % (esc(why), table_html(d, max_rows=80, searchable=nm.replace(".tsv", "").replace("_", "")) if d is not None else pending(why, "76b/" + nm)))
        A(mind("If the held-out enrichment intervals include 1, or the thresholds learned on missing-call units do not transfer to wrong-allele units or pairs, a single pruning percent is not supported "
               "and the right unit of action is the missing call itself (structural), not a metric cut-off."))

    # (4) KIR
    A("<h2 id='kir'>4. KIR</h2><div class='box'><b>Nuances that shape every KIR number.</b> Gene content differs between haplotypes (presence/absence is biology, not error). "
      "Framework genes (3DL3, 2DL4, 3DP1, 3DL2) are expected on a complete haplotype; B-content genes (2DL2, 2DL5A/B, 2DS1/2/3/5, 3DS1) only on B-type content. 2DL2/2DL3, 3DL1/3DS1 are alleles of one locus, so trio labels are also given per paralog group. "
      "Alleles are 3/5/7-digit; the label here is the protein (3-digit) field and novelty is common. A haplotype is complete only if all KIR calls sit on one contig with both 3DL3 and 3DL2. "
      "GRCh38 carries depth only for 10 KIR genes; for the rest depth is NA (never 0), and MAPQ0 fraction is the quantitative mismapping indicator across the segmental duplications.</div>")
    A("<h3>4.1 Structure of the KIR calls</h3>" + (img_tag(figs.get("kir_struct"), "KIR structure", "a, haplotype completeness; b, people with 1 or 2 calls per gene; c, novelty tier per gene.") if figs.get("kir_struct")
                                                    else pending("KIR structure", "77/hap_class_counts.tsv, gene_n_called.tsv, gene_by_novelty_tier.tsv")))
    A("<h3>4.2 yak parent of origin for KIR</h3>")
    g = T.rd("77b", "gene_origin_summary_kir.tsv")
    if g is not None and (g["gene"] == "ALL").any():
        r = g[g["gene"] == "ALL"].iloc[0]
        A("<p>%s windows in %s children; %s ambiguous; %s of resolved windows are strong (&ge;10 k-mers and &ge;5x); median %s parent-specific k-mers per window.</p>" % (
            fmt_n(r["n_windows"]), fmt_n(r["n_children"]), fmt_p(r["ambiguous_pct"]), fmt_p(r["strong_of_resolved_pct"]), fmt_n(r["kmers_median"])))
    A(img_tag(figs.get("kir_origin"), "KIR origin", "a, ambiguous and strong windows per gene; b, k-mers per window.") if figs.get("kir_origin") else pending("KIR origin figure", "77b/*.tsv"))
    for nm, why in (("region_summary_kir.tsv", "region-level yak switch / hamming error by platform"), ("both_haps_same_parent_kir.tsv", "children whose two haplotypes at a gene got the same parent")):
        d = T.rd("77b", nm)
        if d is not None:
            A("<p class='read'>%s</p>%s" % (esc(why), table_html(d, max_rows=40)))
    A("<h3>4.3 KIR trio and pair classes</h3>")
    tg = T.rd("77c", "trio_classes_group.tsv")
    if tg is not None:
        d = tg[(tg["gene"] == "ALL") & (tg["variant"] == "all")]
        if len(d):
            r = d.iloc[0]
            A("<p>At the protein field, %s of %s informative KIR trio records have a child label in neither parent (%s, %s-%s). Restricting to records whose transmitting parent's own label is named "
              "(not novel) changes this, as the figure shows: unnamed labels are not comparable.</p>" % (fmt_n(r["neither_n"]), fmt_n(r["n_informative"]), fmt_p(r["neither_pct"], 2), fmt_p(r["lo"], 2), fmt_p(r["hi"], 2)))
    A(img_tag(figs.get("kir_trio"), "KIR trio and pair", "a, neither-parent rate by variant; b, share of first-degree pairs with no shared label per locus.") if figs.get("kir_trio") else pending("KIR trio figure", "77c/trio_classes_group.tsv"))
    A(mind("In KIR the 'neither' rate mixes wrong alleles with novel-name instability and with presence/absence of paralogs; if the rate collapses when only named labels or only complete haplotypes are used, "
           "the signal is naming, not typing error."))
    A("<h3>4.4 KIR pruning evidence</h3>")
    if T.rd("78b", "prune_curves.tsv") is None:
        A(pending("KIR pruning evidence", "78b/prune_curves.tsv, auc.tsv, power.tsv, per_stratum_auc.tsv, definition_card_kir.tsv, depth_metric_coverage.tsv"))
    else:
        kv = ["B_child", "B_child_strong", "B_worst", "H_hap", "P_pair", "P_pair_complete"]
        A(img_tag(figs.get("prune78"), "KIR pruning", "Errors removed against units removed for the top metrics per KIR view."))
        A(img_tag(figs.get("prune78_res"), "KIR residual", "Residual error after pruning by the best metric."))
        A(img_tag(figs.get("auc78"), "KIR AUC", "AUC per KIR view and metric."))
        dt = decision_table(T, "78b", kv, top=3)
        A(table_html(dt, searchable="dectab78") if dt is not None else "")
        for nm in ("per_stratum_auc.tsv", "power.tsv", "definition_card_kir.tsv", "depth_metric_coverage.tsv"):
            d = T.rd("78b", nm)
            if d is not None:
                A("<p class='read'>%s</p>%s" % (esc(nm), table_html(d, max_rows=80, searchable=nm.replace(".tsv", "").replace("_", ""))))
    d = T.rd("78", "depth_coverage_by_gene.tsv")
    if d is not None:
        A("<h3>4.5 KIR depth coverage by gene</h3>" + table_html(d, max_rows=40))

    # (5) raw evidence
    A("<h2 id='raw'>5. Raw evidence</h2><h3>5.1 Tables in this report</h3>")
    inv = []
    for sub, f in T.all_files():
        d = T.rd(sub, f)
        inv.append(dict(table="%s/%s" % (sub, f), rows=len(d), columns=", ".join(list(d.columns)[:8]) + (" ..." if d.shape[1] > 8 else ""),
                        what=TABLE_DOC.get("%s/%s" % (sub, f), "see the producing script's docstring")))
    A(table_html(pd.DataFrame(inv), searchable="inv") if inv else pending("table inventory", "any aggregate table"))
    A("<h3>5.2 Key tables, searchable</h3>")
    for sub, nm in (("76", "B_class_counts.tsv"), ("76", "A_missing_by_gene.tsv"), ("71f", "trio_parent_status.tsv"), ("77c", "trio_classes_group.tsv"), ("77b", "gene_origin_summary_kir.tsv")):
        d = T.rd(sub, nm)
        if d is not None and len(d) <= 300:
            A("<details><summary>%s/%s (%d rows)</summary>%s</details>" % (sub, nm, len(d), table_html(d, searchable=("t%s%s" % (sub, nm)).replace(".tsv", "").replace("_", ""))))
    A("<h3>5.3 Literature context</h3><ul><li>Pollock, Harrison, Norman 2022 (KIR and HLA class I co-evolution and health; PMC10038757): ligand specificities (KIR2DL1-C2, KIR2DL2/3-C1, KIR3DL1-Bw4, KIR3DL2-A3/11), "
      "KIR-A / KIR-B haplotype structure and the 2DS3-3DS1 segmental duplication (their Figure 2 is a template for a gene-map figure); KIR studies mostly use gene presence/absence; non-European populations are under-studied.</li>"
      "<li>Ahn et al. 2021 (KIR*IMP copy-number imputation in psoriasis; PMC8231283): imputed vs sequencing-based typing concordance 91.2% over six genes, KIR2DL1 weakest at 76%; all psoriasis samples were European. "
      "A reference point for how an allele-level long-read call set compares with imputation.</li>"
      "<li>Not read (access blocked): the Advanced Science paper (10.1002/advs.202521531, Figure 4) and the Immunological Reviews piece (10.1111/imr.70145); neither is cited here.</li></ul>")
    A("<h3>5.4 Glossary</h3>" + table_html(pd.DataFrame(GLOSSARY, columns=["term", "meaning"])))

    # (6) open questions
    A("<h2 id='next'>6. Open questions and next steps</h2><ul>"
      "<li>Verify the trio parent-child relationship by IBD0/IBS0 (the relatedness table in use may lack it; the 'trio unverifiable' row of the definition card says how many).</li>"
      "<li>Check missing slots against short-read heterozygosity at the gene (separates assembly structure from a true absence).</li>"
      "<li>A haplotype-aware origin check using neighbouring genes' origins and yak switch intervals; the 'no switch suspect' subset is the first step.</li>"
      "<li>Sequence-level comparison of child and parent calls (CDS) instead of names, in particular for KIR and for novel calls.</li>"
      "<li>KIR depth on the child's own haplotype contigs for the genes absent from GRCh38; the 991 sequel2 people are not in the KIR call set.</li>"
      "<li>Small-cell disclosure decision (open in DECISIONS.md) before anything leaves the private folder.</li></ul>")

    # (7) think tank
    A("<h2 id='think-tank-section'>7. Think tank</h2><div id='think-tank'></div>")
    A("</main></body></html>")
    return "\n".join(parts)


FIGS = [("missing", fig_missing), ("missing_features", fig_missing_features), ("missing_auc", fig_missing_auc), ("parent_status", fig_parent_status),
        ("status_gene", fig_status_by_gene), ("stream_b", fig_stream_b), ("stream_b_support", fig_stream_b_support), ("stream_c", fig_stream_c),
        ("kir_struct", fig_kir_structure), ("kir_origin", fig_kir_origin), ("kir_trio", fig_kir_trio)]
HLA_VIEWS = ["A_classical8", "A_all", "B_child", "B_worst", "B_child_all", "combined_child", "pair_complete", "pair_all"]
KIR_VIEWS = ["B_child", "B_child_strong", "B_worst", "H_hap", "P_pair", "P_pair_complete"]


def make_figures(T, figdir):
    os.makedirs(figdir, exist_ok=True)
    out = {}
    jobs = list(FIGS) + [
        ("prune76", lambda T, d: fig_prune(T, "76b", d, "fig_prune_hla", HLA_VIEWS)),
        ("prune76_res", lambda T, d: fig_prune_residual(T, "76b", d, "fig_prune_hla_residual", HLA_VIEWS)),
        ("auc76", lambda T, d: fig_auc_views(T, "76b", d, "fig_auc_hla")),
        ("cv76", lambda T, d: fig_cv(T, "76b", d, "fig_cv_hla")),
        ("prune78", lambda T, d: fig_prune(T, "78b", d, "fig_prune_kir", KIR_VIEWS)),
        ("prune78_res", lambda T, d: fig_prune_residual(T, "78b", d, "fig_prune_kir_residual", KIR_VIEWS)),
        ("auc78", lambda T, d: fig_auc_views(T, "78b", d, "fig_auc_kir"))]
    for key, fn in jobs:
        try:
            out[key] = fn(T, figdir)
            log("figure %-12s %s" % (key, "ok" if out[key] else "skipped (tables absent)"))
        except Exception as ex:      # one broken figure must not lose the report; it shows as a labelled placeholder
            out[key] = None
            log("figure %s FAILED: %s: %s" % (key, type(ex).__name__, ex))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", default=os.path.join(REPO, "reports", "hla_popgen", "79_s06_final", "aggregate_tables"))
    ap.add_argument("--out", default=os.path.join(REPO, "reports", "hla_popgen", "79_s06_final"))
    ap.add_argument("--no-figs", action="store_true", help="reuse figures already in <out>/figs")
    a = ap.parse_args(argv)
    T = Tables(a.tables)
    figdir = os.path.join(a.out, "figs")
    if a.no_figs:
        figs = {}
        for k, nm in (("missing", "fig_missing_overview"), ("missing_features", "fig_missing_by_feature"), ("missing_auc", "fig_missing_auc"), ("parent_status", "fig_parent_status"),
                      ("status_gene", "fig_call_status_by_gene"), ("stream_b", "fig_stream_b"), ("stream_b_support", "fig_stream_b_support"), ("stream_c", "fig_stream_c"),
                      ("kir_struct", "fig_kir_structure"), ("kir_origin", "fig_kir_origin"), ("kir_trio", "fig_kir_trio_pair"), ("prune76", "fig_prune_hla"),
                      ("prune76_res", "fig_prune_hla_residual"), ("auc76", "fig_auc_hla"), ("cv76", "fig_cv_hla"), ("prune78", "fig_prune_kir"),
                      ("prune78_res", "fig_prune_kir_residual"), ("auc78", "fig_auc_kir")):
            p = os.path.join(figdir, nm + ".png")
            figs[k] = p if os.path.exists(p) else None
    else:
        figs = make_figures(T, figdir)
    doc = build(T, figs, a.out)
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, "S06_FINAL_REPORT.html")
    with open(path, "w") as f:
        f.write(doc)
    log("wrote %s (%.1f MB)" % (path, os.path.getsize(path) / 1e6))
    return 0


if __name__ == "__main__":
    sys.exit(main())
