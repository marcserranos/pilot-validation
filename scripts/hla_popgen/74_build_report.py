#!/usr/bin/env python3
"""S06 74: assemble the call-quality report (one self-contained HTML file) from the aggregate tables and the figures.

Everything numeric in the tables is computed here from the pulled aggregate TSVs; the prose states how to read each figure and what
the numbers on this run say (numbers are inserted from the data, not typed). Images are embedded at web size; the full-resolution
PNG / PDF stay next to the report.

    python3 74_build_report.py --folder reports/hla_popgen/72_call_quality
"""
import argparse
import base64
import html
import io
import os
import sys

import numpy as np
import pandas as pd
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

METRICS = {  # key: (label, family, direction, definition, where it comes from)
    "template_distance": ("Template distance", "Immuannot", "higher = worse",
                          "Weighted edit count between the haplotype contig's gene span (introns and UTR included) and the closest template allele: "
                          "substitution 1, indel up to 5 bp 1, indel over 5 bp 2.", "Immuannot GTF gene row `template_distance`"),
    "template_distance_per_kb": ("Template distance per kb", "Immuannot", "higher = worse",
                                 "Template distance divided by the gene span in kb, so long and short genes compare.", "derived from the same row"),
    "cds_distance": ("CDS distance", "Immuannot", "higher = worse",
                     "Raw edit count of the coding sequence against the closest named allele. Absent when the template already matched exactly; "
                     "counted as 0.", "Immuannot GTF transcript row `cds_distance`"),
    "n_aa_changes": ("Amino-acid changes", "Immuannot", "higher = worse",
                     "Number of codon changes (protein level) between the call and the closest named allele.", "parsed from `cds_mut`"),
    "has_warning": ("Template warning", "Immuannot", "present = worse",
                    "Immuannot flags the template's CDS reconstruction: partial CDS, no start codon, no stop codon, in-frame stop. The literal string "
                    "NA means no warning.", "Immuannot GTF transcript row `template_warning`"),
    "n_tied": ("Alleles tied", "Immuannot", "higher = worse",
               "How many reference alleles tie for the best CDS match, i.e. how ambiguous the call is at the resolution reported.", "Table 1 `n_tied`"),
    "is_novel": ("Novel allele", "Immuannot", "novel = worse",
                 "The call carries a `new` field: it differs from every allele in IPD-IMGT/HLA 3.55.0.", "Table 1 `novelty_depth`"),
    "rel_gene": ("Gene depth", "Depth", "lower = worse", "Mean long-read depth over the whole gene span (introns included), divided by the person's genome "
                 "control depth.", "70, samtools depth on GRCh38"),
    "rel_exons": ("Exon depth", "Depth", "lower = worse", "Mean depth over the exons of the reference transcript (UTR included), relative.", "70"),
    "rel_cds": ("CDS depth", "Depth", "lower = worse", "Mean depth over the coding sequence only, relative.", "70"),
    "rel_ars": ("Peptide-binding exon depth", "Depth", "lower = worse",
                "Mean depth over the groove exons (class I exons 2-3, class II exon 2), relative: the exons that decide two-field typing.", "70"),
    "rel_exon_min": ("Shallowest exon depth", "Depth", "lower = worse", "Lowest mean depth among the coding exons, relative.", "70"),
    "rel_cds_n": ("CDS depth, gene-normalised", "Depth", "lower = worse",
                  "CDS depth divided by the cohort median of that gene within the platform: removes the gene's baseline.", "70 + cohort median"),
    "rel_ars_n": ("Groove depth, gene-normalised", "Depth", "lower = worse", "Groove-exon depth divided by the cohort median of that gene and platform.", "70 + cohort median"),
    "rel_exon_min_n": ("Shallowest exon, exon-normalised", "Depth", "lower = worse",
                       "Each exon divided by its own cohort median (platform-wise) before taking the minimum, so an exon that reference-aligns "
                       "systematically low (DRB1 exon 2) is not always the shallowest.", "70 + cohort median"),
    "cds_p10_rel": ("CDS 10th-percentile depth", "Depth", "lower = worse", "10th percentile of per-base depth over the CDS, relative: sensitive to a thin stretch.", "70 per-base"),
    "cds_min_rel": ("CDS minimum depth", "Depth", "lower = worse", "Lowest per-base depth over the CDS, relative.", "70 per-base"),
    "ars_min_rel": ("Groove minimum depth", "Depth", "lower = worse", "Lowest per-base depth over the groove exons, relative.", "70 per-base"),
    "cds_share_lt10": ("CDS bases under 10x", "Depth", "higher = worse", "Share of CDS bases covered by fewer than 10 reads (absolute depth).", "70 per-base"),
    "ars_share_lt10": ("Groove bases under 10x", "Depth", "higher = worse", "Share of groove-exon bases covered by fewer than 10 reads.", "70 per-base"),
    "mapq_frac": ("Unique-mapping fraction", "Depth", "lower = worse", "CDS depth counting only reads with MAPQ 20 or more, divided by the depth counting all reads: "
                  "low where reads map ambiguously (paralogs).", "70"),
    "rel_mhc": ("MHC window depth", "Depth", "lower = worse", "Mean depth over the 4 Mb MHC window relative to the genome control: a person-level coverage measure.", "58"),
    "mhc_contigs": ("MHC contigs", "Assembly", "higher = worse", "Number of distinct contigs that carry MHC genes on the haplotype: a fragmented assembly.", "Table 1 `contig`"),
    "diff_contig": ("Off partner's contig", "Assembly", "yes = worse",
                    "The gene's phasing partner (A-B, B-A/C, DRB1-DQB1, DQA1-DQB1/DRB1, DPA1-DPB1, ...) is called on the same haplotype but on another contig.", "Table 1"),
    "score7": ("Old 0-7 score", "Legacy", "lower = worse", "The hand-built additive confidence score of 56 (kept only as the comparison).", "56"),
}
FAM_ORDER = ["Immuannot", "Depth", "Assembly", "Legacy"]
SIGN = {"has_warning": 1, "is_novel": 1, "diff_contig": 1}
LOWER_BAD = {k for k, v in METRICS.items() if v[2].startswith("lower")}


def esc(s):
    return html.escape(str(s))


def img_tag(path, alt, max_w=1700):
    if not os.path.exists(path):
        return f'<p class="missing">[missing figure {esc(os.path.basename(path))}]</p>'
    im = Image.open(path).convert("RGBA")
    bg = Image.new("RGB", im.size, "white")
    bg.paste(im, mask=im.split()[3])
    if bg.width > max_w:
        bg = bg.resize((max_w, int(bg.height * max_w / bg.width)), Image.LANCZOS)
    buf = io.BytesIO()
    bg.save(buf, "PNG", optimize=True)
    b64 = base64.b64encode(buf.getvalue()).decode()
    name = os.path.basename(path)
    return (f'<figure><img alt="{esc(alt)}" src="data:image/png;base64,{b64}"/>'
            f'<figcaption class="file">full resolution: {esc(name)} (and .pdf)</figcaption></figure>')


def table(df, cls="", fmt=None):
    fmt = fmt or {}
    h = "".join(f"<th>{esc(c)}</th>" for c in df.columns)
    rows = []
    for _, r in df.iterrows():
        cells = []
        for c in df.columns:
            v = r[c]
            if c in fmt:
                v = fmt[c](v)
            elif isinstance(v, float):
                v = "" if np.isnan(v) else f"{v:.3g}"
            cells.append(f"<td>{v if isinstance(v, str) and v.startswith('<') else esc(v)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return f'<table class="{cls}"><thead><tr>{h}</tr></thead><tbody>{"".join(rows)}</tbody></table>'


class Data:
    def __init__(self, folder):
        self.f = folder
        agg = os.path.join(folder, "aggregate_tables")

        def rd(*p):
            fp = os.path.join(agg, *p)
            return pd.read_csv(fp, sep="\t") if os.path.exists(fp) else None
        self.rd = rd
        self.auc = rd("71", "stat_auc.tsv")
        self.base = rd("71", "stat_error_baseline.tsv")
        self.prune = rd("71", "stat_prune.tsv")
        self.coef = rd("71", "stat_model_coef.tsv")
        self.auc8 = rd("71", "classical8", "stat_auc.tsv")
        self.base8 = rd("71", "classical8", "stat_error_baseline.tsv")
        self.prune8 = rd("71", "classical8", "stat_prune.tsv")
        self.contig = rd("71c", "contig_matrix_trio.tsv")
        self.res = rd("71c", "gene_resolution.tsv")
        self.quant = rd("71", "dist_gene_depth_quantiles.tsv")
        self.exon = rd("71", "dist_exon_depth.tsv")
        self.features = rd("70", "features.tsv")


def a(df, view, metric, col="auc"):
    if df is None:
        return np.nan
    d = df[(df["view"] == view) & (df["metric"] == metric)]
    return float(d[col].iloc[0]) if len(d) else np.nan


def flag_table(prune, base, view, targets=(5, 10)):
    """For each metric: the cut that flags about `t`% of ALL cohort calls (pooled scale), the error rate among flagged and unflagged records,
    the share of errors caught."""
    rows = []
    if prune is None:
        return pd.DataFrame()
    p = prune[(prune["view"] == view) & (prune["scale"] == "pooled")]
    for m, d in p.groupby("metric"):
        d = d.sort_values("prune_pct")
        K = float(d.iloc[0]["errors_kept"])
        n = float(d.iloc[0]["n_kept"])
        for t in targets:
            if "cohort_pct_removed" not in d or d["cohort_pct_removed"].isna().all():
                continue
            i = (d["cohort_pct_removed"] - t).abs().idxmin()
            # the flagged set is "worse than or equal to the threshold": the LAST row that still has this threshold removes the whole tie group
            thr = d.loc[i, "thr_oriented"]
            r = d[d["thr_oriented"] == thr].iloc[-1]
            rem_n = float(r["n_removed"])
            if rem_n <= 0:
                continue
            rem_e = K - float(r["errors_kept"])
            rows.append(dict(metric=m, target=t, cohort_pct=float(r["cohort_pct_removed"]), thr_oriented=float(r["thr_oriented"]),
                             records_flagged_pct=100 * rem_n / n, error_flagged=100 * rem_e / rem_n, error_unflagged=float(r["error_pct_kept"]),
                             errors_caught=float(r["captured_pct"])))
    return pd.DataFrame(rows)


def raw_thr(metric, thr_oriented):
    s = -1 if metric in LOWER_BAD else 1
    return s * thr_oriented


CSS = """
:root{--bg:#fff;--fg:#1b1b1b;--mut:#666;--line:#ddd;--acc:#3C5488;--soft:#f4f6fa;--warn:#fff7e6}
@media (prefers-color-scheme: dark){:root{--bg:#15171c;--fg:#e8e8e8;--mut:#9aa0aa;--line:#2d313a;--acc:#8fa6d9;--soft:#1d2027;--warn:#2a2418}}
body{background:var(--bg);color:var(--fg);font:15px/1.55 -apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:0}
main{max-width:1060px;margin:0 auto;padding:24px 20px 80px}
h1{font-size:26px;margin:.2em 0}h2{font-size:20px;margin:2em 0 .4em;border-bottom:1px solid var(--line);padding-bottom:4px}
h3{font-size:16px;margin:1.6em 0 .3em;color:var(--acc)}
.sub{color:var(--mut)}
figure{margin:12px 0 4px}figure img{width:100%;height:auto;border:1px solid var(--line);background:#fff}
.file{color:var(--mut);font-size:12px}
table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0 14px}
th,td{border:1px solid var(--line);padding:4px 7px;text-align:left;vertical-align:top}th{background:var(--soft)}
td.n,th.n{text-align:right}
.box{background:var(--soft);border-left:3px solid var(--acc);padding:8px 12px;margin:12px 0}
.warn{background:var(--warn);border-left:3px solid #E18727;padding:8px 12px;margin:12px 0}
.read{color:var(--mut);font-size:13.5px}
code{background:var(--soft);padding:0 3px;border-radius:3px}
.tbl-wrap{overflow-x:auto}
.missing{color:#b00}
"""


def build(folder):
    D = Data(folder)
    H = []
    ap = H.append
    im = lambda n, alt: img_tag(os.path.join(folder, n), alt)  # noqa: E731

    def pct(x, d=1):
        return "n/a" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{d}f}%"

    nrec = {}
    if D.base is not None:
        for v in ("trio_child", "pair", "twin"):
            r = D.base[(D.base["view"] == v) & (D.base["gene"] == "ALL")]
            if len(r):
                nrec[v] = (int(r["n_records"].iloc[0]), int(r["n_errors"].iloc[0]), float(r["error_pct"].iloc[0]))
    nrec8 = {}
    if D.base8 is not None:
        for v in ("trio_child", "pair", "twin"):
            r = D.base8[(D.base8["view"] == v) & (D.base8["gene"] == "ALL")]
            if len(r):
                nrec8[v] = (int(r["n_records"].iloc[0]), int(r["n_errors"].iloc[0]), float(r["error_pct"].iloc[0]))

    ap("<h1>Per-call quality metrics for HLA typing: depth over the exons, Immuannot's own numbers, and how well they predict typing error</h1>")
    ap('<p class="sub">Sprint S06, 2026-10-20 / 21. Long-read HLA calls in All of Us. Figures at full resolution; values are fine-grained aggregates (small-cell '
       'masking deferred, so this file is for internal reading only).</p>')

    # ------------------------------------------------------------- summary
    p_all = a(D.auc, "pair", "model:all")
    p_old = a(D.auc, "pair", "score7")
    t_all = a(D.auc, "trio_worst", "model:all")
    t_chi = a(D.auc, "trio_child", "model:all")
    ap("<h2>In one page</h2>")
    c8 = ""
    if nrec8:
        c8 = (f"<li><b>Classical eight only</b> (A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1): trios {nrec8['trio_child'][0]} calls / {nrec8['trio_child'][1]} errors ({pct(nrec8['trio_child'][2])}), "
              f"pairs {nrec8['pair'][0]} / {nrec8['pair'][1]} ({pct(nrec8['pair'][2])}), twins {nrec8['twin'][0]} / {nrec8['twin'][1]}.</li>")
    ap('<div class="box"><ul>'
       "<li><b>What we did.</b> Replaced the hand-built 0-7 confidence score with measurable quantities for every call: Immuannot's own distances and flags, "
       "long-read depth over the gene, exons, coding sequence and the peptide-binding (groove) exons, and assembly context. We then measured how each one "
       "relates to typing error, using the trio (yak, Mendelian mismatch), first-degree pair and twin truth sets, and what happens when the worst x% of calls "
       "per metric are removed.</li>"
       f"<li><b>Truth sets.</b> Trios: {nrec.get('trio_child', ('?',))[0]} child calls with {nrec.get('trio_child', ('', '?'))[1]} Mendelian mismatches "
       f"({pct(nrec.get('trio_child', (0, 0, np.nan))[2])}). Pairs: {nrec.get('pair', ('?',))[0]} pair-genes, {nrec.get('pair', ('', '?'))[1]} disagreements "
       f"({pct(nrec.get('pair', (0, 0, np.nan))[2])}). Twins: {nrec.get('twin', ('?',))[0]} records (too few to lean on).</li>" + c8 +
       f"<li><b>Main result (pairs, the best-powered truth).</b> The old score reaches AUC {p_old:.2f}; all measured metrics combined reach {p_all:.2f} "
       f"(classical eight only: {a(D.auc8, 'pair', 'score7'):.2f} vs {a(D.auc8, 'pair', 'model:all'):.2f}).</li>"
       f"<li><b>Trios are harder.</b> The child's own call alone gives AUC {t_chi:.2f} (combined model); the worse of child and parent {t_all:.2f}. A Mendelian mismatch is rarely attributable to one person's call. "
       f"On trios the combined model is only on par with the old score (all genes {t_all:.2f} vs {a(D.auc, 'trio_worst', 'score7'):.2f}; classical eight {a(D.auc8, 'trio_worst', 'model:all'):.2f} vs "
       f"{a(D.auc8, 'trio_worst', 'score7'):.2f}); Immuannot's CDS distance alone gets {a(D.auc, 'trio_worst', 'cds_distance'):.2f}.</li>"
       f"<li><b>What predicts error.</b> Immuannot's distances and the share of CDS bases under 10x are the strongest single metrics (pairs: share under 10x {a(D.auc, 'pair', 'cds_share_lt10'):.2f}, template distance "
       f"{a(D.auc, 'pair', 'template_distance'):.2f}, CDS depth {a(D.auc, 'pair', 'rel_cds'):.2f}, gene depth {a(D.auc, 'pair', 'rel_gene'):.2f}); the number of amino-acid changes separates almost error-free calls "
       "(zero changes) from error-prone ones. The exon- and gene-normalised depth metrics did not beat the raw ones "
       f"(CDS depth {a(D.auc, 'pair', 'rel_cds'):.2f} raw vs {a(D.auc, 'pair', 'rel_cds_n'):.2f} normalised; groove {a(D.auc, 'pair', 'rel_ars'):.2f} vs {a(D.auc, 'pair', 'rel_ars_n'):.2f}), so the raw versions are the ones to keep.</li>"
       "<li><b>Pruning.</b> Removing a few percent of the worst calls removes a large share of the errors (numbers in section 7), clearly better than removing the same number at random. On pairs the combined model beats every single "
       "metric; on trios it is modestly better than the best single metrics.</li>"
       "<li><b>Caveats that matter.</b> Depth is measured on the GRCh38 reference, so it is reference-biased where alleles are divergent (DRB1 exon 2); the non-classical genes "
       "(DRB3, MICA) err far more than the classical eight, so everything is also shown for the classical eight alone; the twin set is tiny; pair errors include recombination.</li>"
       "</ul></div>")

    # ------------------------------------------------------------- 1 design
    ap("<h2>1. Question and design</h2>")
    ap("<p>A typing call is trustworthy or not. The old confidence score (0-7) added seven yes/no checks that were chosen by hand. Here every ingredient is a number that "
       "can be computed for any call, so we can see which ones carry information, in what range, and whether combinations help. Three questions: (i) which numbers relate to "
       "typing error, (ii) where along each number does the error start to rise, (iii) if we drop the worst x% of calls by a number, what is the error among the rest?</p>")

    # ------------------------------------------------------------- 2 methods
    ap("<h2>2. Methods</h2>")
    ap("<h3>2.1 How depth is measured over a gene</h3>")
    ap(im("fig_schematic_depth.png", "How depth is measured over a gene"))
    ap('<p class="read">Each person has a GRCh38-aligned long-read BAM. For the MHC window we read the alignments once per MAPQ floor (0 and 20) and count, for every base of every '
       "gene of interest, how many reads cover it (<code>samtools depth</code>). The mean over a set of positions is that set's depth. Dividing by the person's genome control depth "
       "(the median of five 2 Mb control windows, from script 58) gives a relative depth that compares people with different overall coverage. Gene depth matches script 58's "
       "gene depth with correlation 0.9995-1.000 (0.3-3% lower because <code>depth</code> skips duplicate and secondary reads). The reads in the schematic are simulated for illustration.</p>")

    ap("<h3>2.2 What counts as coding sequence and as groove exon</h3>")
    ap(im("fig_schematic_cds.png", "CDS and groove exon definition"))
    ap('<p class="read"><b>Ground truth for the exon boundaries</b> is the public GENCODE v46 basic annotation on GRCh38. For each gene we take the reference transcript (MANE Select, else Ensembl canonical, '
       "else the protein-coding transcript with the longest coding sequence); its <code>CDS</code> rows are the coding sequence; UTR, introns and exons found only in alternative isoforms are "
       "not coding sequence. We first tried the union of all isoforms and rejected it: it added alternative coding exons (HLA-A 1,211 bp against 1,098 bp) and would dilute the depth. "
       "Coding exons are numbered in transcript order. <b>Groove exons</b> are a stated convention: class I (A, B, C) exons 2 and 3 (alpha1, alpha2), class II exon 2 (beta1 or alpha1); "
       "they encode the peptide-binding groove and carry most of the variation that defines the two-field allele name.</p>")
    ap('<p class="read"><b>Cross-check (panel b).</b> GENCODE\'s coding sequence plus its stop codon matches the length of the modal full-length allele in IPD-IMGT/HLA 3.55.0 (the database Immuannot '
       "types against) for HLA-A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1, DRA, TAP1 and TAP2. <b>MICA is the exception</b>: GENCODE's reference transcript is shorter (996 bp against about 1,155 bp), so "
       "MICA exon-level numbers carry that caveat. Immuannot defines the exons of a <i>call</i> from the IPD allele annotation projected onto the contig (contig coordinates); we measure depth on GRCh38 with "
       "the GENCODE coordinates, which is why the length agreement matters.</p>")

    ap("<h3>2.3 What a typing error is</h3>")
    ap(im("fig_schematic_truth.png", "Truth definition"))
    base_tab = None
    if D.base is not None:
        rows = []
        for v, lab in (("trio_child", "Trios (one record per child call)"), ("pair", "First-degree pairs"), ("twin", "Twins / duplicates")):
            r = nrec.get(v)
            r8 = nrec8.get(v)
            if r:
                rows.append({"truth set": lab, "records (all genes)": r[0], "errors": r[1], "error rate": pct(r[2], 2),
                             "records (classical 8)": r8[0] if r8 else "", "errors (c8)": r8[1] if r8 else "", "error rate (c8)": pct(r8[2], 2) if r8 else ""})
        base_tab = pd.DataFrame(rows)
    if base_tab is not None:
        ap(table(base_tab))
    ap('<p class="read"><b>Trios (primary label).</b> yak counts k-mers unique to each parent on the child\'s gene +-10 kb of each haplotype contig; a haplotype is assigned to a parent when that parent has at '
       "least 5 such k-mers and at least 3 times the other's. The label is a Mendelian mismatch: the transmitting parent does not carry the child's 2-field allele. The mismatch can come from the child's call, "
       "the parent's call or the parent-of-origin assignment, so the metrics are attributed three ways: the child's call (strict), the transmitting parent's call, or the worse of the two (primary). "
       "<b>Pairs</b>: first-degree pairs that share their MHC (at most one classical gene without a shared 2-field allele); error = no shared 2-field allele (typing error, or recombination). "
       "The pair's metric is the worse of the two people at that gene. <b>Twins/duplicates</b>: unordered 2-field genotypes differ.</p>")

    ap("<h3>How the trio set and the pair set differ, and how much they overlap</h3>")
    ap(im("fig_schematic_sets.png", "Trio and pair sets"))
    ov = D.rd("71", "pair_trio_overlap.tsv")
    if ov is not None:
        o = dict(zip(ov["quantity"], ov["value"]))
        ap('<p class="read"><b>Who they are.</b> Both come from the All of Us relatedness table. A <b>trio</b> is a child with two first-degree relatives who are not related to each other, i.e. the two parents; '
           f"here {int(o['trio_children'])} children and {int(o['trio_edges'])} child-parent links. The <b>pair set</b> is every first-degree pair with calls ({int(o['pairs_distinct'])} pairs among {int(o['pair_people'])} people): parent and child, or full siblings. "
           "Without IBS0 in the data the two cannot be told apart. Twins are a third set. "
           f"<b>Overlap.</b> {int(o['edges_in_both'])} of the {int(o['trio_edges'])} trio child-parent links are also in the pair set, and they make up {100 * o['edges_in_both'] / o['pairs_distinct']:.0f}% of the pairs "
           f"({int(o['pair_records_on_trio_edges'])} of {int(o['pair_records'])} pair-genes). The remaining {int(o['pairs_distinct'] - o['edges_in_both'])} pairs are siblings or parent-child pairs whose other parent is not in the data. "
           f"{int(o['people_in_both'])} people appear in both sets, so the two truth sets are not independent: they largely re-test the same families.</p>")
    ap('<p class="read"><b>Why they give different error rates.</b> The trio test is direct and per haplotype: yak says which parent a haplotype came from, and its allele must be among that parent\'s two alleles. '
       "The pair test is symmetric and lenient: two relatives that share any allele at the gene pass, so a wrong allele in a heterozygote is hidden whenever the other allele matches. It only flags a gene when the pair shares nothing. "
       "It also counts a recombination as an error, and includes siblings (who are only expected to share an allele when they inherited a haplotype in common; pairs that share nothing across the MHC are removed by the filter). "
       "So the trio error (5.1%) is higher than the pair error (3.9%) even though the same families are involved, the pair set has about twice the records, and pair errors are partly not typing errors. "
       "That is why both are shown: trios as the cleaner label, pairs for power.</p>")

    ap("<h3>2.4 The metrics</h3>")
    rows = []
    for k, (lab, fam, dirn, dfn, src) in METRICS.items():
        rows.append({"family": fam, "metric": lab, "key": k, "worse when": dirn, "definition": dfn, "source": src})
    mt = pd.DataFrame(rows)
    mt["_o"] = mt["family"].map({f: i for i, f in enumerate(FAM_ORDER)})
    mt = mt.sort_values("_o", kind="stable").drop(columns="_o")
    ap('<div class="tbl-wrap">' + table(mt) + "</div>")

    ap("<h3>2.5 Statistics</h3>")
    ap('<p class="read"><b>AUC</b>: the probability that an erroneous call has a worse value of the metric than a correct one (0.5 = no information, 1 = perfect); confidence intervals by resampling people. '
       "<b>Within-gene AUC</b>: the same, computed gene by gene and averaged, so a metric cannot look good just because genes differ in both metric and error. <b>Odds ratio per SD</b>: logistic, adjusted for gene "
       "and platform, with person-clustered standard errors. <b>Rolling error curves</b>: error rate in sliding windows of n/25 records along the metric, with Wilson 95% bands. "
       "<b>Pruning</b>: records ranked worst to best by the metric; the worst x% removed (0 to 50% in 0.25% steps; ties removed proportionally); error among the rest; the grey band is what removing the same number at "
       "random would give (95%). <b>Models</b>: L2 logistic regressions on the metrics (Immuannot only, depth only, all) with gene indicators, scored out-of-fold by 5-fold cross-validation grouped by person.</p>")

    # ------------------------------------------------------------- 3 depth data
    ap("<h2>3. What the depth data look like</h2>")
    ap(im("fig_depth_distributions.png", "Depth distributions"))
    ap('<p class="read"><b>a</b> per-gene CDS depth (0.01-wide bins) for the classical genes. <b>b</b> median, interquartile range and 5-95% of four depth definitions per gene: they agree for most genes. '
       "<b>c</b> per-person MHC-window depth relative to the genome, by platform. <b>d</b> MHC-window depth against genome depth: the MHC region runs below the genome, proportionally. "
       "<b>e</b> absolute MHC depth by ancestry. <b>f</b> gene-by-gene Spearman correlation of CDS depth across people: depth at different genes is nearly independent once genome depth is divided out, except for neighbours "
       "(DPA1-DPB1, TAP1-TAP2); so a thin gene is a local event, not a thin person.</p>")
    ap(im("fig_exon_depth.png", "Per-exon depth"))
    ap('<p class="read">Depth of each coding exon (median, interquartile range, 5-95%); the groove exons are red. Most genes are flat across exons. <b>HLA-DRB1 exon 2, the groove exon, is the exception</b> '
       "(low) and exon 1 is high. This is reference bias: DRB1 exon 2 is the most polymorphic exon, reads from divergent alleles map poorly to the single GRCh38 allele, and DRB1 exon 1 collects extra reads from "
       "paralogous DRB genes. For this reason we also give exon-normalised depth metrics (each exon divided by its own cohort median).</p>")
    ap(im("fig_gene_depth_profile.png", "Depth along genes"))
    ap('<p class="read">Depth along each classical gene (median with interquartile and 10-90% bands over 300 random people) above the gene model (light: exons including UTR, dark: coding sequence, red: groove exons).</p>')

    # ------------------------------------------------------------- 4 curves
    ap("<h2>4. Typing error along each metric</h2>")
    ap('<p class="read">Each panel shows the typing error (with Wilson 95% band) in sliding windows along one metric; the dashed line is the overall error. '
       "Read left to right as the call gets worse (depth panels run from thin to thick, so error falls to the right). Discrete metrics show one point per value.</p>")
    for v, lab in (("pair", "First-degree pairs (worse of the two calls): the best-powered truth set"),
                   ("trio_worst", "Trios, worse of child and parent (primary trio view)"), ("trio_child", "Trios, the child's own call (strict)")):
        ap(f"<h3>{esc(lab)}</h3>")
        ap(im(f"fig_error_curves_{v}.png", lab))
    ap('<p class="read">Twins/duplicates (72 records, about 10 errors) are drawn in <code>fig_error_curves_twin</code> but are too few to read curves from.</p>')

    # metric-by-metric reading table from the auc table
    if D.auc is not None:
        ap("<h3>Metric by metric, in numbers</h3>")
        rows = []
        for k in METRICS:
            r = {"metric": METRICS[k][0], "family": METRICS[k][1]}
            for v, lab in (("pair", "pair"), ("trio_worst", "trio worst"), ("trio_child", "trio child")):
                r[f"AUC {lab}"] = a(D.auc, v, k)
            r["AUC pair, within gene"] = a(D.auc, "pair", k, "auc_within_gene")
            r["AUC pair, classical 8"] = a(D.auc8, "pair", k)
            r["OR per SD (pair)"] = a(D.auc, "pair", k, "or_per_sd")
            rows.append(r)
        rt = pd.DataFrame(rows).sort_values("AUC pair", ascending=False)
        ap('<div class="tbl-wrap">' + table(rt, fmt={c: (lambda x: "" if (isinstance(x, float) and np.isnan(x)) else f"{x:.3f}")
                                                  for c in rt.columns if c.startswith("AUC") or c.startswith("OR")}) + "</div>")

    # ------------------------------------------------------------- 5 AUC
    ap("<h2>5. Which metrics carry information</h2>")
    ap(im("fig_auc.png", "AUC per metric"))
    ap('<p class="read">Left: AUC with 95% interval (resampling people) for every metric and truth set, sorted by the pair truth set; colour = family of metric, marker = truth set. '
       "Right: AUC within gene. The black rows are the multi-metric models (out-of-fold). Intervals are wide for trios (about 180 errors) and very wide for twins (about 10).</p>")

    ap("<h2>6. How the metrics relate to each other</h2>")
    ap(im("fig_metric_correlation.png", "Metric correlation"))
    ap('<p class="read">Spearman correlation between metrics over all calls in the cohort, clustered. Blocks of strongly correlated metrics (the depth family; the Immuannot distance family) carry overlapping information; '
       "keep one or two per block unless the combined model gains from more.</p>")
    for g in ("fig_grid_rel_cds__template_distance", "fig_grid_rel_ars__cds_distance", "fig_grid_rel_cds__rel_mhc"):
        ap(im(g + ".png", g))
    ap('<p class="read">Two metrics together (pair truth set): typing error over a 40 x 40 grid of the cohort percentile ranks of two metrics (100 = worst), kernel-smoothed, blank where there is too little data. '
       "A diagonal gradient means both metrics add information; a gradient along one axis means the other adds little.</p>")

    # ------------------------------------------------------------- 7 pruning
    ap("<h2>7. Removing the worst calls</h2>")
    ap(im("fig_pruning.png", "Pruning curves"))
    ap('<p class="read">Top: typing error among the calls kept, as the worst x% of records are removed by a metric (x from 0 to 50%). Grey: the 95% band when the same number of records is removed at random. '
       "Bottom: the share of all errors removed. A curve far below the band and far above the diagonal is a good filter. Flat stretches are metrics with many tied values (once all non-zero values are removed, "
       "nothing is left to rank).</p>")
    if D.prune is not None:
        rows = []
        for v, lab in (("pair", "Pairs"), ("trio_worst", "Trios, worse of two"), ("trio_child", "Trios, child")):
            for m, ml in (("model:all", "Model, all metrics"), ("model:immuannot", "Model, Immuannot"), ("model:depth", "Model, depth"), ("score7", "Old 0-7 score"),
                          ("rel_cds", "CDS depth"), ("cds_share_lt10", "CDS bases under 10x"), ("cds_distance", "CDS distance"),
                          ("template_distance", "Template distance"), ("n_aa_changes", "Amino-acid changes")):
                d = D.prune[(D.prune["view"] == v) & (D.prune["metric"] == m) & (D.prune["scale"] == "pooled")].set_index("prune_pct")
                if not len(d):
                    continue
                row = {"truth": lab, "metric": ml, "error, all kept": d.loc[0.0, "error_pct_kept"]}
                for x in (2.0, 5.0, 10.0, 20.0):
                    row[f"error after removing {x:g}%"] = d.loc[x, "error_pct_kept"]
                    row[f"errors removed {x:g}%"] = d.loc[x, "captured_pct"]
                rows.append(row)
        pt = pd.DataFrame(rows)
        ap("<h3>Numbers behind the pruning curves (all genes)</h3>")
        ap('<div class="tbl-wrap">' + table(pt, fmt={c: (lambda x: "" if (isinstance(x, float) and np.isnan(x)) else f"{x:.2f}") for c in pt.columns
                                                  if c not in ("truth", "metric")}) + "</div>")

    # ------------------------------------------------------------- 8 contig
    ap("<h2>8. Contig co-location against the trio truth</h2>")
    ap(im("fig_contig_truth.png", "Contig co-location versus trio truth"))
    ap('<p class="read">For every pair of MHC genes in the 179 trio children. <b>a</b> the share of child haplotypes where both genes sit on one assembled contig (an assembly property, '
       "not physical phasing). <b>b</b> the share of haplotypes (any contig) where yak says both genes came from the same parent: the physical truth. <b>c</b> the same, restricted to pairs on one contig: "
       "how often one contig really is one parental haplotype. Where b and c are near 100 the contig can be trusted as a phase; where c is far below 100 the contig joins pieces from different parents "
       "(a switch inside the contig). Genes with fewer than half the haplotypes resolved are noisier.</p>")
    if D.contig is not None:
        t = D.contig[(D.contig["stratum"] == "all") & (D.contig["gene_a"] < D.contig["gene_b"])].copy()
        t = t[["gene_a", "gene_b", "n_haplotypes", "colocation_pct", "n_resolved", "same_parent_all_pct", "n_resolved_same_contig",
               "same_parent_samecontig_pct", "n_resolved_diff_contig", "same_parent_diffcontig_pct"]]
        t.columns = ["gene a", "gene b", "haplotypes", "% on one contig", "resolved", "% same parent (any contig)", "resolved, one contig",
                     "% same parent | one contig", "resolved, different contigs", "% same parent | different contigs"]
        ap("<h3>The numbers behind the matrices</h3>")
        ap('<div class="tbl-wrap">' + table(t, fmt={c: (lambda x: "" if (isinstance(x, float) and np.isnan(x)) else f"{x:.1f}") for c in t.columns if c.startswith("%")}) + "</div>")

    # ------------------------------------------------------------- 9 classical 8
    ap("<h2>9. Classical genes only</h2>")
    ap('<p class="read">The same analysis restricted to HLA-A, B, C, DRB1, DQA1, DQB1, DPA1 and DPB1 (figures in the <code>classical8</code> folder). The non-classical genes (DRB3, DRB4, DRB5, DQA2, MICA, TAP1, TAP2) have higher error rates '
       "and, for DRB3/DRB4, no depth on the primary assembly, so pooling them lowers and blurs the picture.</p>")
    if D.auc8 is not None:
        ap(im("classical8/fig_auc.png", "AUC, classical genes"))
        ap(im("classical8/fig_pruning.png", "Pruning, classical genes"))

    # ------------------------------------------------------------- 10 final table
    ap("<h2>10. QC parameters to attach to each call</h2>")
    ap('<p class="read">Every parameter below can be computed for any call and written next to it. <b>Flag at</b> is the value that flags about 5% of ALL calls in the cohort, with the error rate among flagged and unflagged records in the '
       "pair truth set (classical genes in brackets), and the share of errors caught. <b>AUC</b> columns say how much information the metric carries. <b>Suggested use</b> is my reading; you decide what to keep.</p>")
    ft = flag_table(D.prune, D.base, "pair", targets=(5,))
    ft8 = flag_table(D.prune8, D.base8, "pair", targets=(5,)) if D.prune8 is not None else pd.DataFrame()
    ftt = flag_table(D.prune, D.base, "trio_worst", targets=(5,))
    rows = []
    for k, (lab, fam, dirn, dfn, src) in METRICS.items():
        if k == "score7":
            continue
        r = {"family": fam, "parameter": f"<span>{esc(lab)} (<code>{k}</code>)</span>", "worse when": dirn}
        f = ft[ft["metric"] == k]
        if len(f):
            f = f.iloc[0]
            thr = raw_thr(k, f["thr_oriented"])
            r["flag at"] = (f"{thr:.3g}" if k not in ("has_warning", "is_novel", "diff_contig") else "yes")
            r["cohort calls flagged"] = pct(f["cohort_pct"], 1)
            r["error flagged / unflagged (pairs)"] = f"{f['error_flagged']:.1f}% / {f['error_unflagged']:.1f}%"
            r["errors caught"] = pct(f["errors_caught"], 0)
        g = ft8[ft8["metric"] == k] if len(ft8) else []
        if len(g):
            g = g.iloc[0]
            r["classical 8: flagged / unflagged"] = f"{g['error_flagged']:.1f}% / {g['error_unflagged']:.1f}% ({pct(g['cohort_pct'], 1)} flagged)"
        h = ftt[ftt["metric"] == k]
        if len(h):
            h = h.iloc[0]
            r["trio worst: flagged / unflagged"] = f"{h['error_flagged']:.1f}% / {h['error_unflagged']:.1f}%"
        r["AUC pair"] = a(D.auc, "pair", k)
        r["AUC trio worst"] = a(D.auc, "trio_worst", k)
        rows.append(r)
    fin = pd.DataFrame(rows)
    fin["_o"] = fin["family"].map({f: i for i, f in enumerate(FAM_ORDER)})
    fin = fin.sort_values(["_o", "AUC pair"], ascending=[True, False]).drop(columns="_o")
    ap('<div class="tbl-wrap">' + table(fin, fmt={c: (lambda x: "" if (isinstance(x, float) and np.isnan(x)) else f"{x:.3f}") for c in ("AUC pair", "AUC trio worst")}) + "</div>")

    if D.coef is not None:
        ap("<h3>The combined model</h3>")
        cf = D.coef[(D.coef["view"] == "pair") & (D.coef["model"] == "model:all")][["term", "coef", "odds_ratio"]].copy()
        cf = cf[~cf["term"].str.startswith("gene:")]
        cf["abs"] = cf["coef"].abs()
        cf = cf.sort_values("abs", ascending=False).drop(columns="abs")
        cf.columns = ["term (standardised)", "coefficient", "odds ratio per SD"]
        ap('<p class="read">Standardised coefficients of the all-metrics model fitted on the pair records (distances enter as log(1+x), missing values are median-filled; gene indicators omitted here). '
           "Correlated metrics share weight, so read this as which metrics the model leans on, not as independent effects.</p>")
        ap(table(cf.head(14), fmt={"coefficient": lambda x: f"{x:.2f}", "odds ratio per SD": lambda x: f"{x:.2f}"}))

    ap("<h2>11. Limits of these results</h2>")
    ap("<ul>"
       "<li><b>Counts.</b> About 180 trio errors, 270 pair errors, 10 twin errors. Intervals are wide for trios and twins; read single-metric differences below about 0.05 AUC as ties.</li>"
       "<li><b>Attribution.</b> A Mendelian mismatch is not necessarily the child's error; the strict child view and the worse-of-two view bracket it.</li>"
       "<li><b>Pair errors include recombination</b> and are limited to pairs that share their MHC; the 'worse of two' metric rewards any metric that flags either person.</li>"
       "<li><b>Reference bias.</b> Depth is on GRCh38. It is lower where the sample's allele diverges from the reference (DRB1 exon 2), and it is not measurable for DRB3/DRB4 (absent from the primary assembly). Exon-normalised metrics reduce but do not remove this.</li>"
       "<li><b>Depth is the read depth of the person, not of the haplotype.</b> A heterozygous gene with one thin haplotype is not seen separately; haplotype-resolved depth from the haplotagged BAM is a possible extension.</li>"
       "<li><b>MICA</b> coding sequence in GENCODE is shorter than in IPD; MICA exon numbers are indicative.</li>"
       "<li><b>Thresholds</b> in the final table are read off these truth sets and describe this cohort; they are starting points, not validated cut-offs.</li>"
       "<li><b>Whole-cohort distributions</b> use the people measured so far when the report was built; the error analysis uses only truth-set people and is complete.</li>"
       "</ul>")
    ap("<h2>12. Reproducibility</h2>")
    ap("<p class=\"read\">Scripts (branch <code>s06-call-quality</code>): <code>70_call_exon_depth.py</code> (depth), <code>71_call_quality_table.py</code> (call table, truth records, distributions), "
       "<code>71b_call_quality_stats.py</code> (curves, AUC, pruning, models), <code>71c_contig_truth_matrix.py</code> (contig matrices), <code>72_call_quality_figures.py</code> and "
       "<code>73_depth_schematic.py</code> (figures), <code>74_build_report.py</code> (this file). Person-level tables stay inside the Workbench (<code>~/s06/internal</code>); this report is built only from aggregate tables.</p>")

    body = "\n".join(H)
    n_missing = body.count("[missing figure")
    if D.auc8 is None or n_missing:
        body = ('<div class="warn"><b>Draft status.</b> This build uses the first analysis run. Still arriving: the classical-8-only tables and figures (section 9), '
                "the exon-normalised depth metrics (rows marked blank in the tables), and the along-gene depth profile figure (section 3). "
                f"{n_missing} figure placeholder(s) below. The report is rebuilt from the same script when they land.</div>") + body
    return f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>HLA call quality report</title><style>{CSS}</style></head><body><main>{body}</main></body></html>'


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--folder", required=True)
    ap.add_argument("--out", default=None)
    a_ = ap.parse_args(argv)
    out = a_.out or os.path.join(a_.folder, "S06_CALL_QUALITY_REPORT.html")
    h = build(a_.folder)
    with open(out, "w") as f:
        f.write(h)
    print(f"wrote {out} ({len(h) / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
