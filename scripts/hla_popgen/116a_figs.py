"""116a Figure and layout helpers for the full S08+S09 report (inline SVG, light/dark safe). Imported by 116_full_report.py.

All charts take plain Python lists (the same lists are printed in the report as data tables), so every figure has its numbers next to it.
Colour keys map to CSS variables defined in the report stylesheet: a (blue), b (orange), c (green), d (purple), warn (red), mut (grey).
"""
import math
from html import escape as esc

W = 880
KEY = {"main": "--a", "pos": "--c", "test": "--b", "neg": "--mut", "rev": "--warn", "ref": "--d", "a": "--a", "b": "--b", "c": "--c", "d": "--d", "warn": "--warn", "mut": "--mut"}


def _svg(h, body, w=W):
    return '<svg viewBox="0 0 %d %d" role="img" xmlns="http://www.w3.org/2000/svg">%s</svg>' % (w, h, body)


def _t(x, y, s, size=12, anchor="start", bold=False, muted=False, extra=""):
    return '<text x="%.1f" y="%.1f" font-size="%s" text-anchor="%s"%s%s%s>%s</text>' % (x, y, size, anchor, ' font-weight="600"' if bold else "", ' class="m"' if muted else "", extra, esc(str(s)))


def _line(x1, y1, x2, y2, col="--mut", sw=1, op=0.5, dash=False):
    return '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" style="stroke:var(%s);stroke-width:%s;stroke-opacity:%s%s"/>' % (x1, y1, x2, y2, col, sw, op, ";stroke-dasharray:4 3" if dash else "")


def _rect(x, y, w, h, col, fill=0.85, rx=3, stroke=False):
    return '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="%d" style="fill:var(%s);fill-opacity:%.2f%s"/>' % (x, y, max(w, 0), h, rx, col, fill, ";stroke:var(%s);stroke-width:1.5" % col if stroke else "")


def _ticks(lo, hi, n=5):
    span = hi - lo
    raw = span / n
    mag = 10 ** math.floor(math.log10(raw))
    step = min((m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= raw), default=raw)
    t = math.ceil(lo / step) * step
    out = []
    while t <= hi + 1e-9:
        out.append(0.0 if abs(t) < 1e-12 else round(t, 10))
        t += step
    return out


def forest(rows, lim, title=None, label_w=330, row_h=25, right_w=170, xlabel="", fmt="{:+.3f}"):
    """rows: (label, beta, se, kind) or ("#", header text). lim = (lo, hi). Draws beta with a 95% CI (1.96 se); right column prints the numbers."""
    lo, hi = lim
    x0, x1 = label_w, W - right_w - 12
    sx = lambda v: x0 + (min(max(v, lo), hi) - lo) / (hi - lo) * (x1 - x0)
    top = 34 if title else 10
    h = top + len(rows) * row_h + 44
    p = []
    if title:
        p.append(_t(10, 20, title, 14, bold=True))
    for t in _ticks(lo, hi):
        p.append(_line(sx(t), top - 4, sx(t), h - 36, op=0.12))
        p.append(_t(sx(t), h - 20, ("%g" % t), 11, "middle", muted=True))
    if lo < 0 < hi:
        p.append(_line(sx(0), top - 4, sx(0), h - 36, sw=1.5, op=0.55))
    if xlabel:
        p.append(_t((x0 + x1) / 2, h - 4, xlabel, 11, "middle", muted=True))
    for i, r in enumerate(rows):
        y = top + i * row_h + row_h / 2
        if r[0] == "#":
            p.append(_t(10, y + 4, r[1], 12, bold=True))
            continue
        label, b, se, kind = r[:4]
        col = KEY.get(kind, "--a")
        lab = label if len(label) < 72 else label[:70] + "…"
        p.append(_t(18, y + 4, lab, 12))
        if b is None:
            p.append(_t(x0, y + 4, "n/a", 11, muted=True))
            continue
        l, u = b - 1.96 * se, b + 1.96 * se
        p.append(_line(sx(l), y, sx(u), y, col, sw=2.5, op=0.9))
        p.append('<circle cx="%.1f" cy="%.1f" r="4.5" style="fill:var(%s)"/>' % (sx(b), y, col))
        if l < lo:
            p.append(_t(sx(lo) - 2, y + 4, "◂", 10, "end"))
        if u > hi:
            p.append(_t(sx(hi) + 2, y + 4, "▸", 10))
        p.append(_t(W - right_w, y + 4, (fmt + "  [" + fmt + ", " + fmt + "]").format(b, l, u), 11, muted=True))
    return _svg(h, "".join(p))


def hbars(rows, vmax=None, title=None, label_w=330, bar_h=22, gap=7, fmt="{:,.0f}", unit="", right=150):
    """rows: (label, value, kind[, note]). Horizontal bars with the value printed at the bar end."""
    vmax = vmax or max(r[1] for r in rows) * 1.12
    x0, x1 = label_w, W - right
    top = 34 if title else 8
    h = top + len(rows) * (bar_h + gap) + 8
    p = [_t(10, 20, title, 14, bold=True)] if title else []
    for i, r in enumerate(rows):
        y = top + i * (bar_h + gap)
        col = KEY.get(r[2], "--a")
        w = (r[1] / vmax) * (x1 - x0)
        p.append(_t(x0 - 8, y + bar_h * 0.7, r[0], 12, "end"))
        p.append(_rect(x0, y, w, bar_h, col))
        txt = fmt.format(r[1]) + unit + (("  " + r[3]) if len(r) > 3 else "")
        p.append(_t(x0 + w + 6, y + bar_h * 0.7, txt, 11.5, muted=False))
    return _svg(h, "".join(p))


def grouped(cats, series, ymax, title=None, height=300, fmt="{:.2f}", ylabel="", legend=True, ymin=0.0):
    """cats: category labels; series: [(name, values, kind)]."""
    left, bottom, top = 56, 46, (62 if len(series) > 1 and legend else 40) if title else 20
    pw, ph = W - left - 12, height - bottom - top
    g = pw / len(cats)
    bw = g * 0.78 / len(series)
    p = [_t(10, 20, title, 14, bold=True)] if title else []
    yv = lambda v: top + ph * (ymax - v) / (ymax - ymin)
    for t in range(5):
        val = ymin + (ymax - ymin) * t / 4
        p.append(_line(left, yv(val), W - 8, yv(val), op=0.12))
        p.append(_t(left - 6, yv(val) + 4, fmt.format(val), 11, "end", muted=True))
    p.append(_line(left, yv(0), W - 8, yv(0), sw=1.5, op=0.5))
    if ylabel:
        p.append('<text x="14" y="%.1f" font-size="11" class="m" transform="rotate(-90 14 %.1f)" text-anchor="middle">%s</text>' % (top + ph / 2, top + ph / 2, esc(ylabel)))
    for i, c in enumerate(cats):
        gx = left + i * g + g * 0.11
        for j, (name, vals, kind) in enumerate(series):
            v = vals[i]
            if v is None:
                continue
            col = KEY.get(kind, "--a")
            x = gx + j * bw
            y0, y1 = yv(0), yv(min(max(v, ymin), ymax))
            p.append(_rect(x, min(y0, y1), bw - 2, max(abs(y1 - y0), 1.5), col))
            p.append(_t(x + (bw - 2) / 2, (y1 - 4) if v >= 0 else (y1 + 12), fmt.format(v), 10, "middle"))
        p.append(_t(left + i * g + g / 2, height - bottom + 18 + (0 if ymin >= 0 else 0), c, 12, "middle"))
    if legend and len(series) > 1:
        lx = 14
        ly = 32 if title else 8
        for name, vals, kind in series:
            p.append(_rect(lx, ly, 12, 12, KEY.get(kind, "--a")))
            p.append(_t(lx + 17, ly + 10, name, 11.5))
            lx += 7.2 * len(name) + 40
    return _svg(height, "".join(p))


def dumbbell(rows, xlim, title=None, labels=("before", "after"), kinds=("warn", "c"), label_w=300, row_h=26, vlines=(), xlabel=""):
    """rows: (label, v1, v2). Two points per row joined by a line; vlines: [(x, text)]."""
    lo, hi = xlim
    x0, x1 = label_w, W - 20
    sx = lambda v: x0 + (min(max(v, lo), hi) - lo) / (hi - lo) * (x1 - x0)
    top = 66 if title else 26
    h = top + len(rows) * row_h + 44
    p = [_t(10, 20, title, 14, bold=True)] if title else []
    for t in _ticks(lo, hi):
        p.append(_line(sx(t), top - 6, sx(t), h - 36, op=0.12))
        p.append(_t(sx(t), h - 20, "%g" % t, 11, "middle", muted=True))
    p.append(_line(sx(0), top - 6, sx(0), h - 36, sw=1.5, op=0.5))
    for v, text in vlines:
        p.append(_line(sx(v), top - 6, sx(v), h - 36, "--warn", sw=1.2, op=0.7, dash=True))
        p.append(_t(sx(v), top - 12, text, 10.5, "middle", muted=True))
    if xlabel:
        p.append(_t((x0 + x1) / 2, h - 4, xlabel, 11, "middle", muted=True))
    for i, (lab, a, b) in enumerate(rows):
        y = top + i * row_h + row_h / 2
        p.append(_t(10, y + 4, lab, 12))
        if a is None or b is None:
            continue
        p.append(_line(sx(a), y, sx(b), y, "--mut", sw=2, op=0.6))
        p.append('<circle cx="%.1f" cy="%.1f" r="5" style="fill:var(%s)"/>' % (sx(a), y, KEY[kinds[0]]))
        p.append('<circle cx="%.1f" cy="%.1f" r="5" style="fill:var(%s)"/>' % (sx(b), y, KEY[kinds[1]]))
    lx = 16
    ly = 38 if title else 14
    for name, kind in zip(labels, kinds):
        p.append('<circle cx="%.1f" cy="%d" r="5" style="fill:var(%s)"/>' % (lx, ly, KEY[kind]))
        p.append(_t(lx + 10, ly + 4, name, 11.5))
        lx += 7.2 * len(name) + 44
    return _svg(h, "".join(p))


def stack100(rows, segs, title=None, label_w=300, bar_h=26, gap=10):
    """rows: (label, [v1, v2, ...]) normalised to 100%; segs: [(name, kind)]."""
    x0, x1 = label_w, W - 20
    top = 50 if title else 28
    h = top + len(rows) * (bar_h + gap) + 10
    p = [_t(10, 20, title, 14, bold=True)] if title else []
    lx = x0
    for name, kind in segs:
        p.append(_rect(lx, top - 24, 12, 12, KEY[kind]))
        p.append(_t(lx + 16, top - 14, name, 11.5))
        lx += 7 * len(name) + 36
    for i, (lab, vals) in enumerate(rows):
        y = top + i * (bar_h + gap)
        tot = sum(vals) or 1
        p.append(_t(x0 - 8, y + bar_h * 0.68, lab, 12, "end"))
        x = x0
        for (name, kind), v in zip(segs, vals):
            w = v / tot * (x1 - x0)
            p.append(_rect(x, y, w, bar_h, KEY[kind], 0.85, rx=0))
            if v / tot >= 0.06:
                p.append(_t(x + w / 2, y + bar_h * 0.68, "%.0f%%" % (100 * v / tot), 11.5, "middle", bold=True, extra=' style="fill:#fff"'))
            x += w
    return _svg(h, "".join(p))


def heat(rows, cols, matrix, vmax=8.0, title=None, label_w=300, cell_w=96, cell_h=28, fmt="{:+.1f}", note=""):
    """Diverging heat table. matrix[r][c] is a float or None; positive = orange, negative = blue; cell text is the value."""
    top = 56 if title else 36
    h = top + len(rows) * cell_h + (26 if note else 10)
    p = [_t(10, 20, title, 14, bold=True)] if title else []
    for j, c in enumerate(cols):
        p.append(_t(label_w + j * cell_w + cell_w / 2, top - 8, c, 12, "middle", bold=True))
    for i, lab in enumerate(rows):
        y = top + i * cell_h
        p.append(_t(10, y + cell_h * 0.66, lab, 12))
        for j in range(len(cols)):
            v = matrix[i][j]
            x = label_w + j * cell_w
            if v is None:
                p.append(_rect(x + 2, y + 2, cell_w - 4, cell_h - 4, "--mut", 0.06, rx=4))
                p.append(_t(x + cell_w / 2, y + cell_h * 0.66, "n/a", 11, "middle", muted=True))
                continue
            a = min(abs(v) / vmax, 1.0)
            p.append(_rect(x + 2, y + 2, cell_w - 4, cell_h - 4, "--b" if v > 0 else "--a", 0.12 + 0.78 * a, rx=4))
            p.append(_t(x + cell_w / 2, y + cell_h * 0.66, fmt.format(v), 12, "middle", bold=abs(v) >= 4))
    if note:
        p.append(_t(10, h - 6, note, 11, muted=True))
    return _svg(h, "".join(p))


def flow(boxes, title=None, box_w=None, h=96):
    """boxes: [(title, subtitle, kind)] left to right."""
    n = len(boxes)
    box_w = box_w or (W - 20 - (n - 1) * 22) / n
    top = 34 if title else 10
    p = ['<defs><marker id="fa" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto"><path d="M0 0 L6 3 L0 6z" fill="currentColor" opacity=".6"/></marker></defs>']
    if title:
        p.append(_t(10, 20, title, 14, bold=True))
    for i, (t, s, k) in enumerate(boxes):
        x = 10 + i * (box_w + 22)
        p.append('<rect x="%.1f" y="%d" width="%.1f" height="%d" rx="10" style="fill:var(%s);fill-opacity:.10;stroke:var(%s);stroke-width:2"/>' % (x, top, box_w, h - 20, KEY[k], KEY[k]))
        p.append(_t(x + box_w / 2, top + 28, t, 11.5 if box_w < 150 else 12.5, "middle", bold=True))
        for li, line in enumerate(s.split("|")):
            p.append(_t(x + box_w / 2, top + 48 + li * 14, line, 10.5, "middle", muted=True))
        if i < n - 1:
            p.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" style="stroke:currentColor;stroke-opacity:.5;stroke-width:1.5" marker-end="url(#fa)"/>' % (x + box_w + 1, top + (h - 20) / 2, x + box_w + 20, top + (h - 20) / 2))
    return _svg(top + h, "".join(p))


def sbars(rows, lim, title=None, label_w=300, bar_h=18, gap=6, fmt="{:+.1f}", xlabel="", marks=()):
    """Signed horizontal bars from zero. rows: (label, value, kind[, note]) or ("#", header). lim=(lo, hi); values beyond are clipped with an arrow but their number is printed.
    marks: [(x, text)] dashed vertical reference lines."""
    lo, hi = lim
    x0, x1 = label_w, W - 130
    sx = lambda v: x0 + (min(max(v, lo), hi) - lo) / (hi - lo) * (x1 - x0)
    top = (46 if marks else 34) if title else 14
    h = top + len(rows) * (bar_h + gap) + 40
    p = [_t(10, 20, title, 14, bold=True)] if title else []
    for t in _ticks(lo, hi):
        p.append(_line(sx(t), top - 4, sx(t), h - 30, op=0.12))
        p.append(_t(sx(t), h - 14, "%g" % t, 11, "middle", muted=True))
    p.append(_line(sx(0), top - 4, sx(0), h - 30, sw=1.5, op=0.55))
    for v, text in marks:
        p.append(_line(sx(v), top - 4, sx(v), h - 30, "--warn", sw=1.2, op=0.7, dash=True))
        p.append(_t(sx(v), top - 8, text, 10, "middle", muted=True))
    if xlabel:
        p.append(_t((x0 + x1) / 2, h - 1, xlabel, 11, "middle", muted=True))
    for i, r in enumerate(rows):
        y = top + i * (bar_h + gap)
        if r[0] == "#":
            p.append(_t(10, y + bar_h * 0.7, r[1], 12, bold=True))
            continue
        label, v, kind = r[:3]
        col = KEY.get(kind, "--a")
        p.append(_t(x0 - 8, y + bar_h * 0.7, label, 12, "end"))
        a, b = sorted((sx(0), sx(v)))
        p.append(_rect(a, y, max(b - a, 1.5), bar_h, col))
        clip = "▸" if v > hi else ("◂" if v < lo else "")
        txt = fmt.format(v) + clip + (("  " + r[3]) if len(r) > 3 else "")
        if v >= 0:
            p.append(_t(max(sx(v), sx(0)) + 6, y + bar_h * 0.7, txt, 11.5))
        else:
            p.append(_t(min(sx(v), sx(0)) - 6, y + bar_h * 0.7, txt, 11.5, "end"))
    return _svg(h, "".join(p))


# ---------------------------------------------------------------- HTML helpers
def table(headers, rows, cls=""):
    th = "".join("<th>%s</th>" % esc(str(h)) for h in headers)
    tr = "".join("<tr>%s</tr>" % "".join("<td>%s</td>" % (c if str(c).startswith("<") else esc(str(c))) for c in r) for r in rows)
    return '<div class="tw"><table%s><tr>%s</tr>%s</table></div>' % (' class="%s"' % cls if cls else "", th, tr)


def fig(svg, caption, source=None, data=None):
    """Figure with caption, a one-line source tag and an optional collapsed data table."""
    d = ""
    if data:
        d = '<details class="dt"><summary>numbers behind this figure</summary>%s</details>' % data
    s = '<p class="src">Source: %s</p>' % source if source else ""
    return '<figure>%s<figcaption>%s</figcaption>%s%s</figure>' % (svg, caption, s, d)


def ev(title, inner, kind="ev", open_=False):
    """Collapsible block. kind: ev (evidence), bg (background), how (method), cav (caveats)."""
    title = title[len("Background: "):] if title.startswith("Background: ") else title
    label = {"ev": "Evidence", "bg": "Background", "how": "How it was done", "cav": "Caveats and red-team"}[kind]
    return '<details class="%s"%s><summary><span class="k">%s</span> %s</summary><div class="in">%s</div></details>' % (kind, " open" if open_ else "", label, esc(title), inner)


def badge(kind, text=None):
    names = {"strong": "Strong", "lead": "Lead", "open": "Open question", "known": "Known / validation", "null": "Null result", "retired": "Retired", "notrun": "Not run"}
    return '<span class="badge %s">%s</span>' % (kind, esc(text or names[kind]))


def cards(items):
    return '<div class="cards">%s</div>' % "".join('<div class="card"><b>%s</b><span>%s</span></div>' % (esc(a), b) for a, b in items)


# ---------------------------------------------------------------- page shell
CSS = """
:root{--bg:#fbfaf7;--fg:#1d2330;--mut:#6b7280;--a:#2f6f9f;--b:#c8663a;--c:#5f8f3a;--d:#7a5ea8;--warn:#c43d3d;--card:#fff;--line:rgba(120,120,120,.28)}
@media (prefers-color-scheme:dark){:root{--bg:#14171c;--fg:#e6e8ec;--mut:#9aa3b2;--a:#6fb0e0;--b:#e88a5c;--c:#8cc467;--d:#b69be0;--warn:#ef7b7b;--card:#1b1f26}}
*{box-sizing:border-box}html{scroll-behavior:smooth}
body{margin:0;background:var(--bg);color:var(--fg);font:16px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}
main{max-width:920px;margin:0 auto;padding:22px 18px 80px}
h1{font-size:30px;line-height:1.2;margin:.2em 0 .3em}h2{font-size:23px;margin:2.2em 0 .5em;padding-top:.6em;border-top:1px solid var(--line)}h3{font-size:17px;margin:1.6em 0 .4em}
p{margin:.55em 0}ul,ol{margin:.4em 0 .6em;padding-left:1.3em}li{margin:.25em 0}code{font-size:.88em;background:rgba(128,128,128,.14);padding:1px 5px;border-radius:4px}
.sub{color:var(--mut);font-size:14px;margin:0}.m{fill:var(--mut)}svg text{fill:var(--fg)}svg text.m{fill:var(--mut)}
.toc{display:flex;flex-wrap:wrap;gap:6px 14px;font-size:14px;margin:12px 0 4px}.toc a{color:var(--a);text-decoration:none}.toc a:hover{text-decoration:underline}
.glance{display:grid;grid-template-columns:1fr;gap:14px;margin:16px 0}
.box{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:14px 18px}
.box h3{margin:.1em 0 .5em;font-size:16px;letter-spacing:.02em}.box.imp{border-left:5px solid var(--c)}.box.prom{border-left:5px solid var(--a)}.box.under{border-left:5px solid var(--mut)}
.box ol,.box ul{padding-left:1.15em}.box li{margin:.5em 0}.why{display:block;color:var(--mut);font-size:14px}
.badge{display:inline-block;font-size:12px;font-weight:600;padding:1px 9px;border-radius:99px;border:1px solid;margin-right:4px;vertical-align:1px;white-space:nowrap}
.badge.strong{color:var(--c);border-color:var(--c)}.badge.lead{color:var(--a);border-color:var(--a)}.badge.open{color:var(--b);border-color:var(--b)}.badge.known{color:var(--d);border-color:var(--d)}
.badge.null{color:var(--mut);border-color:var(--mut)}.badge.retired{color:var(--warn);border-color:var(--warn)}.badge.notrun{color:var(--mut);border-color:var(--mut);border-style:dashed}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:10px;margin:12px 0}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px 12px}.card b{display:block;font-size:22px;line-height:1.2;color:var(--a)}.card span{font-size:13.5px;color:var(--mut)}
.bl{border-left:5px solid var(--a);background:rgba(120,120,120,.07);padding:10px 16px;border-radius:0 10px 10px 0;margin:10px 0 14px}
.bl b.lbl{font-size:12.5px;text-transform:uppercase;letter-spacing:.06em;color:var(--mut);display:block}
.cheat{font-size:14.5px}.cheat dt{font-weight:600;float:left;clear:left;width:110px}.cheat dd{margin:0 0 .35em 118px}
figure{margin:14px 0;background:var(--card);border:1px solid var(--line);border-radius:12px;padding:10px 12px}figure svg{width:100%;height:auto;display:block}
figcaption{font-size:14.5px;margin-top:8px}.src{font-size:12.5px;color:var(--mut);margin:.3em 0 0}
details{margin:10px 0;border:1px solid var(--line);border-radius:10px;background:rgba(120,120,120,.04)}details>summary{cursor:pointer;padding:9px 14px;font-weight:600;font-size:15px;list-style-position:inside}
details>.in{padding:2px 16px 12px}details[open]>summary{border-bottom:1px solid var(--line)}
details.ev{border-left:5px solid var(--a)}details.bg{border-left:5px solid var(--d)}details.how{border-left:5px solid var(--c)}details.cav{border-left:5px solid var(--b)}
details .k{display:inline-block;font-size:11.5px;text-transform:uppercase;letter-spacing:.06em;color:var(--mut);min-width:78px;margin-right:6px}
details.dt{margin:8px 0 0;border-style:dashed}details.dt>summary{font-weight:500;font-size:13.5px;color:var(--mut)}
.tw{overflow-x:auto;margin:8px 0}table{border-collapse:collapse;width:100%;font-size:14px}th,td{border-bottom:1px solid var(--line);padding:6px 9px;text-align:left;vertical-align:top}th{font-size:12.5px;text-transform:uppercase;letter-spacing:.04em;color:var(--mut)}
.note{border-left:5px solid var(--b);background:rgba(200,102,58,.08);padding:9px 14px;border-radius:0 10px 10px 0;margin:10px 0}
.kicker{font-size:12.5px;text-transform:uppercase;letter-spacing:.07em;color:var(--mut)}
@media (min-width:860px){.glance.three{grid-template-columns:1fr}}
"""


def page(title, kicker, h1, lede, toc, body):
    t = "".join('<a href="#%s">%s</a>' % (i, esc(n)) for i, n in toc)
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>%s</title><style>%s</style></head>'
            '<body><main><p class="kicker">%s</p><h1>%s</h1>%s<nav class="toc">%s</nav>%s</main></body></html>') % (esc(title), CSS, esc(kicker), esc(h1), lede, t, body)


def section(sid, title, badge_html, bottom, cards_html, blocks):
    return '<h2 id="%s">%s %s</h2><div class="bl"><b class="lbl">Bottom line</b>%s</div>%s%s' % (sid, esc(title), badge_html, bottom, cards_html, blocks)
