"""94 S08 final report (HTML, inline SVG figures; no external assets).

Built from aggregate numbers only (every count shown is >= 20 or written '<20'); no participant-level data. The think-tank section is read
from an optional JSON file (--thinktank) so the report can be regenerated after more discussion rounds.
  python3 scripts/hla_popgen/94_s08_report.py --out reports/hla_popgen/91_novel_recurrence/S08_report.html --thinktank thinktank.json
"""
import argparse
import html
import json
import os

ANC = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
# mean substitution-only / indel-bearing novel-protein calls per person, classical genes (92_novel_burden, 12,233 persons; groups >= 20 persons)
REVIO_INDEL = [.443, .381, .343, .363, .303, .304]
SEQ2E_INDEL = [.705, .875, .746, .621, .864, .979]
REVIO_SNV = [.023, .013, .020, .013, .019, .013]
SEQ2E_SNV = [.037, .024, .008, .017, .000, .021]
ALL_SNV = [.022, .013, .019, .013, .020, .014]
ALL_NOVEL = [3.42, 2.86, 3.52, 2.50, 3.06, 3.31]
PERSONS = [3108, 2777, 1494, 3078, 495, 1258]
# 95: effective number of flanking HLA-A + (G or F) backgrounds among carrier haplotypes vs median of 500 matched non-carrier draws
HAP_NAMES = ["F T237I", "F P205T", "F Y106C\n(= F*01:25, known)", "G syn.\n(= G*01:01:33, known)"]
HAP_CARRIER = [1.12, 1.41, 1.60, 1.00]
HAP_CONTROL = [8.26, 7.20, 5.33, 7.45]
# 95: mean substitution-only novel signatures per person by template-distance tercile (distance from the reference), within ancestry
TD_LOW = [.0331, .0326, .0599, .0184, .0301, .0307]
TD_MID = [.0545, .0344, .1135, .0210, .0663, .0570]
TD_HIGH = [.0912, .0717, .1556, .0530, .0798, .0652]

CSS = """
:root{--bg:#fbfaf7;--fg:#1d1d1f;--mut:#666a73;--card:#fff;--line:#e2e0d8;--a:#2f6f8f;--b:#c8693a;--c:#6a8f3a;--d:#8a5aa8;--warn:#b3402a}
@media (prefers-color-scheme:dark){:root{--bg:#16181c;--fg:#ecebe6;--mut:#9aa0a8;--card:#1f2228;--line:#333842;--a:#6fb3d2;--b:#e8956a;--c:#9cc46a;--d:#c19adb;--warn:#f08a70}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:16px/1.55 -apple-system,Segoe UI,Helvetica,Arial,sans-serif}
main{max-width:980px;margin:0 auto;padding:32px 16px 80px}h1{font-size:2rem;line-height:1.2;margin:.2em 0}h2{margin:2.2em 0 .4em;font-size:1.35rem;border-top:1px solid var(--line);padding-top:1em}
h3{margin:1.4em 0 .3em;font-size:1.05rem}p,li{max-width:70ch}.sub{color:var(--mut)}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:12px;margin:18px 0}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px}.card b{display:block;font-size:1.7rem;color:var(--a)}.card span{color:var(--mut);font-size:.9rem}
figure{margin:18px 0;background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px}figcaption{color:var(--mut);font-size:.9rem;margin-top:6px}svg{width:100%;height:auto;display:block}
table{border-collapse:collapse;width:100%;font-size:.92rem;background:var(--card)}th,td{border:1px solid var(--line);padding:6px 8px;text-align:left;vertical-align:top}th{background:rgba(127,127,127,.1)}
.tag{display:inline-block;padding:1px 8px;border-radius:99px;font-size:.78rem;border:1px solid var(--line)}.known{color:var(--c)}.new{color:var(--b)}.art{color:var(--warn)}
.note{border-left:4px solid var(--b);padding:6px 12px;background:rgba(200,105,58,.08);margin:14px 0}.tt{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:12px}
code{background:rgba(127,127,127,.15);padding:0 4px;border-radius:4px;font-size:.9em}svg text{fill:var(--fg);font-family:inherit}svg .m{fill:var(--mut)}
"""


def esc(s):
    return html.escape(str(s))


def hbars(rows, width=900, label_w=210, bar_h=26, gap=8, color_var="--a", fmt="{:,}", vmax=None):
    vmax = vmax or max(v for _, v, *_ in rows)
    h = len(rows) * (bar_h + gap) + 10
    out = ['<svg viewBox="0 0 %d %d" role="img">' % (width, h)]
    for i, r in enumerate(rows):
        lab, v = r[0], r[1]
        col = r[2] if len(r) > 2 else color_var
        txt = r[3] if len(r) > 3 else fmt.format(v)
        y = 5 + i * (bar_h + gap)
        w = max(2, (width - label_w - 90) * v / vmax)
        out.append('<text x="0" y="%d" font-size="14">%s</text>' % (y + bar_h * 0.7, esc(lab)))
        out.append('<rect x="%d" y="%d" width="%.1f" height="%d" rx="4" style="fill:var(%s)"/>' % (label_w, y, w, bar_h, col))
        out.append('<text x="%.1f" y="%d" font-size="14">%s</text>' % (label_w + w + 8, y + bar_h * 0.7, esc(txt)))
    out.append("</svg>")
    return "".join(out)


def grouped(series, names, cols, title, width=900, height=300, ymax=None, fmt="{:.2f}"):
    ymax = ymax or max(max(s) for s in series) * 1.15
    left, bottom, top = 50, 40, 28
    pw, ph = width - left - 10, height - bottom - top
    g = pw / len(names)
    bw = g * 0.7 / len(series)
    out = ['<svg viewBox="0 0 %d %d" role="img"><text x="%d" y="18" font-size="14" font-weight="600">%s</text>' % (width, height, left, esc(title))]
    for t in range(5):
        y = top + ph - ph * t / 4
        out.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="currentColor" opacity=".12"/>' % (left, width - 10, y, y))
        out.append('<text class="m" x="%d" y="%.1f" font-size="11" text-anchor="end">%s</text>' % (left - 6, y + 4, fmt.format(ymax * t / 4)))
    for i, n in enumerate(names):
        gx = left + i * g + g * 0.15
        for j, s in enumerate(series):
            h = ph * s[i] / ymax
            out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="3" style="fill:var(%s)"/>' % (gx + j * bw, top + ph - h, bw - 2, h, cols[j][1]))
        out.append('<text x="%.1f" y="%d" font-size="13" text-anchor="middle">%s</text>' % (left + i * g + g / 2, height - 16, esc(n)))
    lx = width - 220
    for j, (lab, cv) in enumerate(cols):
        out.append('<rect x="%d" y="%d" width="12" height="12" rx="2" style="fill:var(%s)"/><text x="%d" y="%d" font-size="12">%s</text>' % (lx, 6 + j * 16, cv, lx + 18, 17 + j * 16, esc(lab)))
    out.append("</svg>")
    return "".join(out)


def flow():
    boxes = [
        ("Long-read assemblies", "~12.2k people, 14 HLA genes", "--a"),
        ("Immuannot calls vs IPD 3.55", "'new' = no matching CDS", "--a"),
        ("Rebuild each alignment", "public IPD CDS + cs string", "--b"),
        ("Classify the edits", "indel / substitution / synonymous", "--b"),
        ("Independent checks", "short reads, depth, platform, ancestry", "--c"),
        ("Re-test vs IPD 3.66", "47,419 named CDSs", "--c"),
        ("Recurrence + biology", "unrelated carriers, literature", "--d"),
        ("Think tank", "next experiments", "--d"),
    ]
    w, h, bw, bh = 960, 250, 215, 70
    out = ['<svg viewBox="0 0 %d %d" role="img">' % (w, h)]
    for i, (t, s, c) in enumerate(boxes):
        x = 8 + (i % 4) * (bw + 26)
        y = 12 + (i // 4) * 120
        out.append('<rect x="%d" y="%d" width="%d" height="%d" rx="10" style="fill:none;stroke:var(%s);stroke-width:2"/>' % (x, y, bw, bh, c))
        out.append('<text x="%d" y="%d" font-size="14" font-weight="600" text-anchor="middle">%s</text>' % (x + bw / 2, y + 28, esc(t)))
        out.append('<text class="m" x="%d" y="%d" font-size="12" text-anchor="middle">%s</text>' % (x + bw / 2, y + 50, esc(s)))
        if i % 4 < 3:
            out.append('<path d="M%d %d h22" stroke="currentColor" opacity=".5" marker-end="url(#ar)"/>' % (x + bw + 2, y + bh / 2))
    out.append('<path d="M%d %d C %d %d, 120 %d, 120 %d" fill="none" stroke="currentColor" opacity=".35" stroke-dasharray="4 4" marker-end="url(#ar)"/>' % (8 + 3 * (bw + 26) + bw / 2, 12 + bh, 8 + 3 * (bw + 26) + bw / 2, 118, 112, 132))
    out.append('<defs><marker id="ar" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto"><path d="M0 0 L6 3 L0 6z" fill="currentColor" opacity=".6"/></marker></defs></svg>')
    return "".join(out)


def experiments_html(ase):
    hap = grouped([HAP_CARRIER, HAP_CONTROL], [n.split("\n")[0] for n in HAP_NAMES], [("carrier haplotypes", "--b"), ("matched non-carriers (median)", "--a")],
                  "Effective number of flanking backgrounds (HLA-A + HLA-G/F) - lower = fewer origins", ymax=12, height=280, fmt="{:.0f}")
    td = grouped([TD_LOW, TD_MID, TD_HIGH], ANC, [("low divergence", "--a"), ("middle", "--c"), ("high divergence", "--b")],
                 "Substitution-only novel signatures per person by divergence-from-reference tercile", ymax=0.17, height=280, fmt="{:.2f}")
    ase_html = ("<p>%s</p>" % esc(ase)) if ase else "<p class='sub'>RNA-seq allele-specific expression is running on the VM (personal two-haplotype reference, bias-controlled); result pending.</p>"
    return """
<p>The think tank’s consensus order was run as far as the data allow without new sequencing. Each test includes <b>known alleles as positive controls</b> (F*01:25 and G*01:01:33 have the same recurrence pattern but are named in IPD), so “the two unnamed alleles look like a real allele” can be judged against alleles we <i>know</i> are real.</p>
<h3>A. One haplotype origin or many?</h3>
<figure>%s<figcaption>A recurrent <i>mutation</i> or a recurrent <i>artefact</i> would appear on many different backgrounds; an inherited allele sits on one. Both unnamed HLA-F signatures are on essentially a single flanking haplotype (2 and 3 distinct backgrounds among carrier haplotypes), exactly like the two known controls, whereas matched non-carrier haplotypes of the same HLA-F lineage and ancestry carry 5–8 effective backgrounds (carriers fall below every one of 500 control draws).</figcaption></figure>
<h3>B. Hidden batch effects: platform and sequencing centre</h3>
<p>Observed carriers per platform and per centre match the ancestry-matched expectation for all four signatures (e.g. Revio 22 observed vs 22.1 expected for HLA-F T237I; the Broad centre 86 vs 85.7 for the G control). No platform- or centre-specific pattern, so the run-fingerprint concern is not supported for these alleles. (One centre is mildly above expectation for T237I; not significant at these counts.)</p>
<h3>C. Do relatives agree?</h3>
<p>Only a few carriers have a relative in the long-read cohort (counts &lt;20 are withheld), but wherever one exists the relative carries the same signature in almost every case — for the HLA-F T237I and the G control alike. Consistent with inheritance, too few for a formal transmission rate. Carriers with blood RNA-seq: fewer than 20 for each HLA-F signature (exact counts withheld).</p>
<h3>D. Reference-bias proxy (inconclusive by design)</h3>
<figure>%s<figcaption>Within every ancestry, the novel-substitution burden rises with how far a person’s HLA haplotypes sit from the reference templates. That is expected if novel alleles live on poorly-represented haplotypes — but it is <i>also</i> what a mapping artefact would do, so this proxy cannot separate them. The real discriminator is the simulated-read null (reads simulated from planted alleles through the same pipeline), still to do.</figcaption></figure>
<h3>E. RNA-seq allele-specific expression</h3>%s""" % (hap, td, ase_html)


def thinktank_html(tt):
    if not tt:
        return "<p class='sub'>No think-tank output attached.</p>"
    out = ["<div class='tt'>"]
    for r in tt.get("roles", []):
        out.append("<div class='card'><b style='font-size:1.05rem'>%s</b><ul>%s</ul><span>Favourite: %s</span></div>" % (
            esc(r["role"]), "".join("<li>%s</li>" % esc(i) for i in r["ideas"]), esc(r.get("favourite", ""))))
    out.append("</div>")
    if tt.get("round2"):
        out.append("<h3>Round 2 — they read each other</h3><ul>%s</ul>" % "".join("<li>%s</li>" % esc(x) for x in tt["round2"]))
    if tt.get("synthesis"):
        out.append("<h3>Synthesis: what we would do next</h3><ol>%s</ol>" % "".join("<li>%s</li>" % x for x in tt["synthesis"]))
    return "".join(out)


def build(tt, ase=None):
    funnel = hbars([
        ("Novel HLA calls (vs IPD 3.55)", 8015, "--a"),
        ("  novel-protein", 7635, "--a"),
        ("  substitution-only, rebuildable", 675, "--b"),
        ("  distinct signatures", 441, "--b"),
        ("  seen in >= 2 people", 49, "--c"),
        ("  unnamed in 3.66 and recurrent", 2, "--d", "2 (both HLA-F)"),
    ], label_w=300)
    short = hbars([
        ("Shortlist (protein-altering, subst.-only)", 337, "--a"),
        ("  not an exact 3.66 match", 279, "--b"),
        ("  and no 3.66 partial hit", 270, "--b"),
        ("  + tier-A, <=15 edits, plain ref", 98, "--c"),
    ], label_w=330)
    plat = (grouped([REVIO_INDEL, SEQ2E_INDEL], ANC, [("Revio", "--a"), ("Sequel II/IIe", "--b")], "Indel-bearing novel-protein calls per person (classical genes)", ymax=1.2) +
            grouped([REVIO_SNV, SEQ2E_SNV], ANC, [("Revio", "--a"), ("Sequel II/IIe", "--b")], "Substitution-only novel-protein calls per person (classical genes)", ymax=0.04, fmt="{:.3f}"))
    anc = (grouped([ALL_SNV], ANC, [("all platforms", "--c")], "Substitution-only novel-protein calls per person, by ancestry", ymax=0.024, height=250, fmt="{:.3f}") +
           grouped([ALL_NOVEL], ANC, [("all platforms", "--d")], "Novel calls of any kind per person (incl. synonymous/beyond-CDS)", ymax=4.0, height=250, fmt="{:.1f}"))
    rows = [
        ("HLA-G synonymous (on G*01:01:03:01)", "88", "EAS 82, others <20", "known", "= G*01:01:33"),
        ("HLA-F T237I (on F*01:01:01:01)", "27", "AFR-predominant, some AMR", "new", "no 3.66 allele has this CDS"),
        ("HLA-F P205T (on F*01:01:02:01)", "<20", "AFR-predominant", "new", "no 3.66 allele has this CDS"),
        ("HLA-F Y106C", "<20", "AMR only", "known", "= F*01:25"),
        ("DQB1 A57D (Asp57)", "<20", "AFR + AMR", "known", "= DQB1*02:180:02 (T1D-relevant residue)"),
        ("HLA-E K195R", "<20", "AFR + AMR", "known", "= E*01:152"),
        ("HLA-G D130H / G107S / K170E", "<20 each", "AFR+AMR / EAS / EAS", "known", "= G*01:58 / 01:64 / 01:66"),
        ("DRB1 R262T; B/C many-edit substitutions", "<20", "n/a", "art", "3 different references; truncated or suffix-flagged templates"),
        ("DPB1 T16M", "<20", "EAS", "known", "CpG transition; base already in older IPD"),
    ]
    trs = "".join("<tr><td>%s</td><td>%s</td><td>%s</td><td><span class='tag %s'>%s</span></td><td>%s</td></tr>" % (
        esc(a), esc(b), esc(c), d, {"known": "named in IPD", "new": "unnamed", "art": "artefact"}[d], esc(e)) for a, b, c, d, e in rows)
    return """<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>S08 novel HLA alleles report</title><style>%s</style></head><body><main>
<p class="sub">Omni-HLA · sprint S08 · 23 Oct 2026 · aggregate results only (counts of 1–19 are written &lt;20)</p>
<h1>Are the “novel” HLA alleles real, and do they matter?</h1>
<p>We took 8,015 HLA calls that Immuannot labelled <i>new</i> in ~12,200 long-read genomes and tried to break them. The short answer: <b>most of the apparent novelty is the reference database, not biology</b> — and what survives is small, concentrated in non-classical genes, and organised by ancestry.</p>
<div class="cards">
<div class="card"><b>86%%</b><span>of novel-protein calls are indel-only: mostly 1-bp homopolymer consensus errors</span></div>
<div class="card"><b>~2×</b><span>more indel “novelty” on Sequel II/IIe than Revio, at every ancestry</span></div>
<div class="card"><b>62 / 441</b><span>substitution signatures are already named alleles in IPD 3.66 (54 exact, 8 partial)</span></div>
<div class="card"><b>2</b><span>recurrent signatures still unnamed: HLA-F T237I and P205T</span></div></div>

<h2>1. Method in one picture</h2>
<figure>%s<figcaption>Every step uses only Immuannot’s own alignments plus public IPD sequence. Participant-level files stay in the Workbench; only suppressed aggregates leave it.</figcaption></figure>

<h2>2. From 8,015 calls to a defensible set</h2>
<figure>%s<figcaption>Funnel over the whole cohort (91_novel_recurrence.py on the VM, one core). “Rebuildable” = the alignment can be reconstructed from the public reference allele and Immuannot’s cs string.</figcaption></figure>
<figure>%s<figcaption>The 337-case shortlist (local 93_shortlist_ipd_release.py): 58 are exact matches to a named allele in IPD 3.66 and 9 contain a partial entry. A further filter (tier-A short-read support, ≤15 edits, reference allele not suffix-flagged) leaves 98, dominated by HLA-F, -E and DQA1.</figcaption></figure>

<h2>3. Skeptic test: sequencing platform</h2>
<p>If “novel” meant error, it would follow the sequencing platform; if real, it would follow ancestry. Indels do the first, substitutions the second.</p>
<figure>%s<figcaption>12,233 people (Revio ≈10.4k, Sequel II/IIe ≈0.9k, groups ≥20 persons shown). Indel-bearing novel-protein calls roughly double on Sequel II/IIe. Substitution-only calls do not inflate comparably (the Sequel MID bar rests on 22 people — noise).</figcaption></figure>

<h2>4. Ancestry structure — the hallmark of real variants</h2>
<figure>%s<figcaption>Per-person burden of substitution-only novel-protein calls is higher in AFR, EAS and MID than EUR/AMR/SAS, and novel calls of any kind per person order AFR ≈ EAS &gt; SAS &gt; MID &gt; AMR &gt; EUR — the pattern expected when the reference is built mostly from European donors.</figcaption></figure>

<h2>5. The recurrent signatures, one by one</h2>
<table><tr><th>Signature</th><th>People</th><th>Unrelated carriers by ancestry</th><th>Status vs IPD 3.66</th><th>Note</th></tr>%s</table>
<div class="note"><b>Why T237I and P205T matter.</b> Both sit in exon 4 (α3 domain, outside the peptide groove), each on a single reference allele, with balanced strands, mapping quality ≈60, no CpG context and clean short-read support, and both follow ancestry. HLA-F has few characterised variants (open conformers bind KIR3DS1), so effects are unknown — the first test is whether they travel on a conserved haplotype and are transmitted in trios.</div>

<h2>6. What it would take to claim a new allele</h2>
<ul><li>IPD-IMGT/HLA requires phased, full-length exons (2–3 class I; 2 class II), the new allele sequenced apart from its partner, and names are assigned by the WHO Nomenclature Committee.</li>
<li>Transmission in trios is the missing evidence: only 827 of 12,233 long-read people have a relative in the long-read cohort, so it must come from S06’s yak trio tables.</li>
<li>Immuannot’s own paper found ~41%% of assembled HLA/KIR sequences novel against v3.52, concentrated in non-classical genes — consistent with what we see.</li></ul>

<h2>7. Honest limits</h2>
<ul><li>Novelty was defined against IPD 3.55; 3.66 changes the answer, 3.67 will again.</li>
<li>Substitution signatures are keyed on one reference allele; ancestry labels are AoU predictions; relatedness data is sparse.</li>
<li>No phenotype comparison yet: the EHR/condition query path is unprobed, and rare alleles have little power.</li>
<li>Platform effects were tested only for classical genes; DRB3/4/5 have no short-read check.</li></ul>

<h2>8. Think-tank experiments, run</h2>%s

<h2>9. Think tank</h2>%s

<h2>10. Reproduce</h2>
<p><code>91_novel_recurrence.py</code> (VM, <code>--ipd-nuc hla_nuc.fasta</code>), <code>92_novel_burden.py</code> (VM), <code>93_shortlist_ipd_release.py</code> (local), <code>95_signature_context.py</code> and <code>96_f_signature_ase.py</code> (VM), <code>94_s08_report.py</code> (this page) on branch <code>s08-case-browser</code>. VM runs: one core, <code>nice -n 10</code>, throwaway clone, removed afterwards; raw tables and carrier lists are mode 600 on the VM and are not in the repo.</p>
</main></body></html>""" % (CSS, flow(), funnel, short, plat, anc, trs, experiments_html(ase), thinktank_html(tt))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--thinktank", default=None)
    ap.add_argument("--ase", default=None, help="text file with the public ASE summary (optional)")
    a = ap.parse_args(argv)
    tt = json.load(open(a.thinktank)) if a.thinktank and os.path.exists(a.thinktank) else None
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as f:
        f.write(build(tt, open(a.ase).read() if a.ase and os.path.exists(a.ase) else None))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
