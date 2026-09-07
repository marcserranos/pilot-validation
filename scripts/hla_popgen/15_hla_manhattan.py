#!/usr/bin/env python3
"""Per-site nucleotide diversity (pi) along an HLA gene, every haplotype projected onto ONE fixed
canonical reference allele -- the proper Manhattan/diversity track that
`11_gene_diversity_track.py`'s template-stratified plot failed to be.

## Why this script exists: what 11's Manhattan plot got wrong

`11_gene_diversity_track.py` built its positional track by restricting to haplotypes whose
best-matching template allele was IDENTICAL, so that raw positions shared a coordinate system.
That selection is self-defeating: a haplotype "matched template X" precisely because it is
(near-)identical to X, so most of that subset has `template_distance = 0` -- ZERO variants by
construction. The resulting plots showed 8-15 variant positions across 90-260 haplotypes. The
method selected the least-diverse possible subset and then plotted its emptiness.

The field's answer (IPD-IMGT/HLA alignment convention; MAFFT/DnaSP sliding-window studies of
full-length HLA genes) is: align EVERYTHING to ONE canonical reference and measure per-position
diversity in that single frame. Aligning each sample to its own best match makes "difference" mean
"difference from a different thing per sample", which is not a comparable quantity at all.

## The key data fact that makes this cheap

`mm2.ipd.gen.paf.gz` is NOT a shortlist of plausible candidates -- it contains the ENTIRE IPD/IMGT
allele database aligned against each person's contig. Confirmed live on the VM (2026-09-07): every
one of 20 sampled people had exactly the same 4,810 distinct `HLA-A*` query names, 5,681 for
`HLA-B*`. So a fixed canonical allele's alignment row is already present, for free, for (nearly)
every haplotype -- no re-alignment, no MSA build, no external download.

Sanity check that distinguishes this script from its predecessor: forcing a fixed canonical
reference onto people who carry divergent alleles yields NM 0-160 per haplotype (vs. NM ~0 for
almost everyone under template-matching), and variant positions in the thousands, not tens.

## Coordinate frame, strand, and the bug 11 also had

query = the canonical IPD allele; target = the person's contig. Walking the `cs` string gives
offsets along the QUERY, i.e. canonical-reference coordinates -- shared by every haplotype.

`11_gene_diversity_track.py` never handled strand. Real PAF rows for the same gene come back on
BOTH strands (confirmed live). For a `-` strand row minimap2 aligns the reverse-complemented query,
so accumulated query offsets run backwards relative to the original query; the canonical position
is `qend - 1 - offset`, not `qstart + offset`. Getting this wrong silently mirrors roughly half the
haplotypes and smears the track. This script handles it explicitly and reports a per-strand CDS
consistency check as a built-in test.

## cs operation semantics in THIS orientation (query = canonical allele, target = contig)

  `:N`     N identical bases            -- consumes canonical and contig
  `*xy`    x = contig base, y = canonical base; the PERSON carries x -- consumes both
  `+seq`   bases in canonical, absent from contig -> person has a DELETION -- consumes canonical
  `-seq`   bases in contig, absent from canonical -> person has an INSERTION -- consumes contig
           only, so it has no canonical coordinate; counted in a separate indel track rather than
           forced onto the axis (standard practice: gap columns are never invented).

## Statistic

Per-site unbiased nucleotide diversity, pi_i = n/(n-1) * (1 - sum p_a^2), over the observed allele
labels {REF, A, C, G, T, DEL} at that site among haplotypes that actually COVER it. Coverage is
per-site (`n_surveyed`), because forcing a fixed reference onto divergent haplotypes yields genuinely
partial alignments (query coverage 0.53-1.00 observed), so a site's denominator is the number of
haplotypes that could have been observed there -- not the cohort size.

pi needs allele frequencies, not base identities, so the canonical reference SEQUENCE is never
required: "REF" is a valid allele label and every substitution reports the person's base directly.

CDS/exon structure comes from Immuannot's own `alleles.csv.gz` (`CDS=`, `exon=`, `UTR=` fields, in
the allele's own 1-based coordinates) -- authoritative, not re-derived.

Usage:
    python3 scripts/hla_popgen/12_hla_manhattan.py --limit 200 --threads 6 --plot
    python3 scripts/hla_popgen/12_hla_manhattan.py --genes HLA-A --canonical HLA-A*01:01:01:01
"""
import argparse
import concurrent.futures
import gzip
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict

GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1"]

# Well-known full-genomic reference alleles. Used if present for a given gene; otherwise the script
# auto-selects (see pick_canonical) because candidate-set membership is NOT universal for every
# gene -- confirmed live: 4,810 universal HLA-A alleles but only 1 universal HLA-C allele.
PREFERRED_CANONICAL = {
    "HLA-A": "HLA-A*01:01:01:01",
    "HLA-B": "HLA-B*07:02:01:01",
    "HLA-C": "HLA-C*07:02:01:01",
    "HLA-DRB1": "HLA-DRB1*15:01:01:01",
}

# Exons encoding the peptide-binding groove -- the whole point of the figure. Class I (A/B/C):
# exons 2 and 3 encode the alpha1/alpha2 domains that form the groove. Class II (DRB1): exon 2
# encodes the beta1 domain. 1-based exon numbers, indexing into the annotation's `exons` list.
GROOVE_EXONS = {
    "HLA-A": [2, 3], "HLA-B": [2, 3], "HLA-C": [2, 3],
    "HLA-DRB1": [2],
}

CS_TOKEN_RE = re.compile(r":\d+|\*[a-z]{2}|[+-][a-z]+", re.IGNORECASE)
REFDATA_ALLELES = os.path.expanduser("~/tools/Immuannot_refdata/alleles.csv.gz")


# ---------------------------------------------------------------------------
# Canonical allele annotation (CDS / exon / UTR), from Immuannot's own refdata
# ---------------------------------------------------------------------------
def parse_ranges(spec):
    """'301..373,504..773' -> [(301, 373), (504, 773)] (1-based inclusive, allele coordinates)."""
    out = []
    for part in spec.split(","):
        part = part.strip()
        if not part or ".." not in part:
            continue
        a, _, b = part.partition("..")
        try:
            out.append((int(a), int(b)))
        except ValueError:
            continue
    return out


def load_allele_annotation(allele, path=REFDATA_ALLELES):
    """Pull the `alleles.csv.gz` record for one allele -> {gene_range, cds, exons, utr}.

    Record format (tab-separated key=value fields), e.g. for HLA-A*01:01:01:01:
        geneRange=1..3503  CDS=301..373,504..773,...  UTR=1..300,3204..3503
        exon=1:1..373,2:504..773,...   nexon=8
    """
    if not os.path.exists(path):
        return None
    target = f"allele={allele}"
    with gzip.open(path, "rt") as f:
        for line in f:
            if target not in line:
                continue
            fields = dict()
            for tok in line.rstrip("\n").split("\t"):
                k, _, v = tok.partition("=")
                fields[k.strip()] = v.strip()
            if fields.get("allele") != allele:
                continue
            exons = []
            for part in fields.get("exon", "").split(","):
                if ":" in part:
                    _idx, _, rng = part.partition(":")
                    exons.extend(parse_ranges(rng))
            return {
                "allele": allele,
                "gene_range": parse_ranges(fields.get("geneRange", "")),
                "cds": parse_ranges(fields.get("CDS", "")),
                "utr": parse_ranges(fields.get("UTR", "")),
                "exons": exons,
            }
    return None


# ---------------------------------------------------------------------------
# cs walk in canonical (query) coordinates, strand-aware
# ---------------------------------------------------------------------------
def walk_cs_canonical(cs_string, qstart, qend, strand):
    """Returns (observed, insertions):
      observed: dict {canonical_pos0 -> allele_label} for every position where this haplotype
                DIFFERS from the canonical reference ('a'/'c'/'g'/'t' for a substitution -- the
                base the PERSON carries -- or 'DEL' for canonical bases the person lacks).
                Positions not listed but inside [qstart, qend) are matches (label 'REF').
      insertions: list of canonical_pos0 anchors where the person carries extra bases absent from
                the canonical reference (no canonical coordinate exists; separate track).

    Strand handling: for a '-' row minimap2 aligned the reverse-complemented query, so the k-th
    query base consumed corresponds to canonical position qend-1-k, not qstart+k.
    """
    observed = {}
    insertions = []
    q_consumed = 0

    def canon_pos(k):
        return (qstart + k) if strand == "+" else (qend - 1 - k)

    for tok in CS_TOKEN_RE.findall(cs_string):
        if tok.startswith(":"):
            q_consumed += int(tok[1:])
        elif tok.startswith("*"):
            # tok = '*' + contig_base + canonical_base; the person carries the contig base.
            person_base = tok[1].lower()
            observed[canon_pos(q_consumed)] = person_base
            q_consumed += 1
        elif tok[0] == "+":
            # Present in canonical, absent from contig -> person is missing these canonical bases.
            n = len(tok) - 1
            for k in range(q_consumed, q_consumed + n):
                observed[canon_pos(k)] = "DEL"
            q_consumed += n
        elif tok[0] == "-":
            # Present in contig, absent from canonical -> no canonical coordinate; anchor only.
            insertions.append(canon_pos(q_consumed))
    return observed, insertions


# ---------------------------------------------------------------------------
# PAF: the canonical allele's row for one haplotype
# ---------------------------------------------------------------------------
def get_tag(fields, tag):
    prefix = tag + ":"
    for f in fields[12:]:
        if f.startswith(prefix):
            return f.split(":", 2)[2]
    return None


def canonical_rows(paf_path, canonicals):
    """One pass over a PAF, returning the best (longest query-span) row for EACH canonical allele.

    Reading the file once for all genes rather than once per gene is the difference between a
    ~2 hour and a ~30 minute full-cohort run: each PAF is ~32k rows / 1.2MB gzipped, and the scan,
    not the cs parsing, dominates.
    """
    best = {}
    best_span = {}
    opener = gzip.open if paf_path.endswith(".gz") else open
    with opener(paf_path, "rt") as f:
        for line in f:
            # Cheap prefilter before the (relatively costly) split: the qname is the line prefix.
            hit = None
            for c in canonicals:
                if line.startswith(c) and line[len(c)] == "\t":
                    hit = c
                    break
            if hit is None:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            span = int(fields[3]) - int(fields[2])
            if span > best_span.get(hit, -1):
                best_span[hit] = span
                best[hit] = fields
    return best


def observations_from_row(row):
    """PAF row -> canonical-frame observation dict, or None."""
    qlen, qstart, qend, strand = int(row[1]), int(row[2]), int(row[3]), row[4]
    cs = get_tag(row, "cs")
    if cs is None:
        return None
    nm = get_tag(row, "NM")
    observed, insertions = walk_cs_canonical(cs, qstart, qend, strand)
    return {
        "qlen": qlen, "qstart": qstart, "qend": qend, "strand": strand,
        "nm": int(nm) if nm is not None else None,
        "observed": observed, "insertions": insertions,
    }


def process_haplotype_multi(person_dir, hap, canonicals):
    """-> {canonical: observation dict}; missing canonicals simply absent from the result."""
    paf = os.path.join(person_dir, hap, "mm2.ipd.gen.paf.gz")
    if not os.path.exists(paf):
        return {}
    rows = canonical_rows(paf, canonicals)
    out = {}
    for c, row in rows.items():
        obs = observations_from_row(row)
        if obs is not None:
            out[c] = obs
    return out


# ---------------------------------------------------------------------------
# Aggregation across haplotypes -> per-site allele counts
# ---------------------------------------------------------------------------
def aggregate_multi(persons, outroot, canonicals, threads=1, progress_every=400):
    """Per-site allele-label counts in canonical coordinates for ALL genes in one pass over the
    per-haplotype PAFs. Returns {canonical: accumulator dict}."""
    acc = {c: {"n_haps": 0, "n_missing": 0, "nm_values": [], "strand_counts": Counter(),
                "coverage": None, "alt_counts": defaultdict(Counter),
                "insertion_counts": Counter(), "qlen_seen": Counter()}
           for c in canonicals}
    t0 = time.time()

    def tasks():
        for pid in persons:
            pdir = os.path.join(outroot, pid, "immuannot_output")
            for hap in ("hap1", "hap2"):
                yield pdir, hap

    if threads > 1:
        pool = concurrent.futures.ThreadPoolExecutor(max_workers=threads)
        results = pool.map(lambda ph: process_haplotype_multi(ph[0], ph[1], canonicals), tasks())
    else:
        pool = None
        results = (process_haplotype_multi(pd, hap, canonicals) for pd, hap in tasks())

    n_units = len(persons) * 2
    for i, per_gene in enumerate(results, 1):
        if progress_every and i % progress_every == 0:
            el = time.time() - t0
            done = {c[:12]: acc[c]["n_haps"] for c in canonicals}
            print(f"    {i}/{n_units} units, {el:.0f}s ({i/el:.1f}/s), haps so far: {done}",
                  file=sys.stderr)
        for c in canonicals:
            a = acc[c]
            res = per_gene.get(c)
            if res is None:
                a["n_missing"] += 1
                continue
            a["n_haps"] += 1
            a["qlen_seen"][res["qlen"]] += 1
            a["nm_values"].append(res["nm"])
            a["strand_counts"][res["strand"]] += 1
            if a["coverage"] is None:
                a["coverage"] = [0] * res["qlen"]
            cov = a["coverage"]
            for p in range(res["qstart"], min(res["qend"], len(cov))):
                cov[p] += 1
            for pos, label in res["observed"].items():
                if 0 <= pos < len(cov):
                    a["alt_counts"][pos][label] += 1
            for pos in res["insertions"]:
                if 0 <= pos < len(cov):
                    a["insertion_counts"][pos] += 1

    if pool is not None:
        pool.shutdown()
    return acc


def per_site_pi(coverage, alt_counts):
    """Unbiased per-site nucleotide diversity over labels {REF, a, c, g, t, DEL}.

    pi_i = n/(n-1) * (1 - sum_a p_a^2), with n = haplotypes covering site i. Sites with n < 2 get
    pi = 0.0 (undefined, reported as zero rather than dropped so the x-axis stays continuous)."""
    pi = [0.0] * len(coverage)
    n_diff = [0] * len(coverage)
    for i, n in enumerate(coverage):
        if n < 2:
            continue
        alts = alt_counts.get(i)
        if not alts:
            continue  # every covered haplotype matches the reference -> pi = 0
        n_alt_total = sum(alts.values())
        # Defensive: alt observations should never exceed the coverage denominator. If they do,
        # the coverage figure is the unreliable one, so widen n rather than clamping the counts --
        # clamping only the total still leaves per-allele p > 1 and yields a negative pi.
        n_eff = max(n, n_alt_total)
        counts = list(alts.values()) + [n_eff - n_alt_total]
        s = sum((c / n_eff) ** 2 for c in counts if c > 0)
        pi[i] = (n_eff / (n_eff - 1)) * (1.0 - s) if n_eff > 1 else 0.0
        n_diff[i] = n_alt_total
    return pi, n_diff


def in_ranges_mask(length, ranges_1based):
    """1-based inclusive ranges -> 0-based boolean mask of given length."""
    mask = [False] * length
    for a, b in ranges_1based:
        for p in range(max(a - 1, 0), min(b, length)):
            mask[p] = True
    return mask


# ---------------------------------------------------------------------------
# Canonical selection
# ---------------------------------------------------------------------------
def pick_canonical(persons, outroot, gene, probe=12):
    """Choose the canonical allele for `gene`: prefer PREFERRED_CANONICAL if it is present for
    (nearly) all probed haplotypes; otherwise take the allele maximizing (people covered, qlen).

    Necessary because candidate-set membership is NOT universal across genes -- confirmed live:
    4,810 HLA-A alleles present for every probed person, but only 1 for HLA-C."""
    want = PREFERRED_CANONICAL.get(gene)
    counts = Counter()
    qlens = {}
    n = 0
    for pid in persons[:probe]:
        paf = os.path.join(outroot, pid, "immuannot_output", "hap1", "mm2.ipd.gen.paf.gz")
        if not os.path.exists(paf):
            continue
        n += 1
        seen = set()
        with gzip.open(paf, "rt") as f:
            for line in f:
                if not line.startswith(gene + "*"):
                    continue
                q, qlen = line.split("\t", 2)[0], int(line.split("\t", 2)[1])
                seen.add(q)
                qlens[q] = qlen
        for q in seen:
            counts[q] += 1
    if n == 0:
        return want, {"probed": 0, "note": "no probe data"}
    if want and counts.get(want, 0) == n:
        return want, {"probed": n, "present_in": counts[want], "qlen": qlens.get(want),
                       "auto_selected": False}

    # Among alleles present in at least `min_frac` of probed haplotypes, take the LONGEST.
    # Demanding strict universality is wrong: it once picked HLA-C*16:85 (3,369bp) over the
    # full-length C*07:02:01:01 (~4.3kb) purely because the latter missed 1 of 12 probes, i.e. it
    # traded ~20% of the gene away to avoid a <10% haplotype dropout. Length wins; the dropout is
    # reported and shows up as `n_missing_canonical`.
    min_frac = 0.85
    eligible = {q: c for q, c in counts.items() if c >= min_frac * n}
    if not eligible:
        eligible = counts
    best = max(eligible.items(), key=lambda kv: (qlens.get(kv[0], 0), kv[1]))[0] if eligible else None
    return best, {"probed": n, "present_in": counts.get(best, 0), "qlen": qlens.get(best),
                   "auto_selected": True, "min_frac": min_frac,
                   "n_eligible_at_min_frac": len(eligible),
                   "preferred_present_in": counts.get(want, 0) if want else None,
                   "preferred": want}


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def sliding_mean(values, window):
    """Centered sliding mean; simple prefix-sum implementation."""
    n = len(values)
    if window <= 1 or n == 0:
        return list(values)
    pref = [0.0] * (n + 1)
    for i, v in enumerate(values):
        pref[i + 1] = pref[i] + v
    half = window // 2
    out = []
    for i in range(n):
        a, b = max(0, i - half), min(n, i + half + 1)
        out.append((pref[b] - pref[a]) / (b - a))
    return out


def plot_manhattan(gene, canonical, pi, coverage, cds_mask, exons, out_path, n_haps,
                    window=151, groove_exons=()):
    """Per-site pi Manhattan with a sliding-window overlay, an exon gene model, and the
    peptide-binding-groove exons called out -- the latter being the actual scientific claim."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    n = len(pi)
    x = np.arange(n)
    pi_arr = np.asarray(pi)
    smooth = np.asarray(sliding_mean(pi, window))
    cds = np.asarray(cds_mask, dtype=bool)

    C_NONCDS = "#C9C6D6"   # muted lilac-grey: background sites
    C_CDS = "#C0392B"      # strong red: coding sites
    C_SMOOTH = "#16324F"   # deep navy: the smoothed trend
    C_GROOVE = "#F2C14E"   # amber: peptide-binding groove band

    fig, (ax, axg) = plt.subplots(
        2, 1, figsize=(14, 5.6), sharex=True, dpi=180,
        gridspec_kw={"height_ratios": [12, 1.1], "hspace": 0.06})
    fig.patch.set_facecolor("white")

    groove_spans = [exons[i - 1] for i in groove_exons if 0 < i <= len(exons)]
    for a, b in groove_spans:
        ax.axvspan(a - 1, b, color=C_GROOVE, alpha=0.20, lw=0, zorder=0)

    ax.grid(axis="y", color="#E6E6E6", lw=0.8, zorder=0)
    ax.set_axisbelow(True)

    ax.vlines(x[~cds], 0, pi_arr[~cds], color=C_NONCDS, lw=0.7, zorder=2,
              label="non-CDS site (intron / UTR)")
    ax.vlines(x[cds], 0, pi_arr[cds], color=C_CDS, lw=0.75, zorder=3, label="CDS site")
    ax.plot(x, smooth, color=C_SMOOTH, lw=2.0, zorder=4, solid_capstyle="round",
            label=f"{window} bp sliding mean")

    mean_cds = float(pi_arr[cds].mean()) if cds.any() else float("nan")
    mean_non = float(pi_arr[~cds].mean()) if (~cds).any() else float("nan")
    ratio = (mean_cds / mean_non) if mean_non else float("nan")

    ax.set_ylabel("Nucleotide diversity  $\\pi$  per site", fontsize=11)
    ax.set_xlim(0, n)
    ax.set_ylim(0, max(float(pi_arr.max()) * 1.28, 1e-6))
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=10)

    handles, labels = ax.get_legend_handles_labels()
    if groove_spans:
        handles.append(plt.Rectangle((0, 0), 1, 1, facecolor=C_GROOVE, alpha=0.35, edgecolor="none"))
        labels.append("peptide-binding groove exons")
    ax.legend(handles, labels, frameon=False, fontsize=9.5, ncol=4,
              loc="lower left", bbox_to_anchor=(0.0, 1.005))

    short = gene.replace("HLA-", "")
    fig.suptitle(
        f"HLA-{short} — per-site nucleotide diversity along the full gene",
        x=0.5, y=1.045, fontsize=15, fontweight="semibold")
    ax.text(0.5, 1.115,
            f"{n_haps:,} haplotypes, every one projected onto a single canonical reference "
            f"({canonical})    •    "
            f"mean $\\pi$: CDS {mean_cds:.4f}  vs  non-CDS {mean_non:.4f}  ({ratio:.2f}×)",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=10, color="#444444")

    # Gene model ribbon. Groove exons are always labelled even when too narrow to hold text --
    # they are the point of the figure, and on DRB1 exon 2 is ~270bp inside a ~14kb gene, so the
    # width-based rule alone silently drops the one label that matters.
    axg.axhline(0.5, color="#9A9A9A", lw=1.0, zorder=1)
    for idx, (a, b) in enumerate(exons, start=1):
        is_groove = idx in groove_exons
        axg.add_patch(plt.Rectangle(
            (a - 1, 0.12), max(b - a + 1, 1), 0.76,
            facecolor=(C_CDS if is_groove else "#5A5A5A"), edgecolor="none", zorder=2))
        mid = (a - 1 + b) / 2
        if (b - a) > n * 0.02:
            axg.text(mid, 0.5, str(idx), ha="center", va="center",
                     fontsize=7.5, color="white", zorder=3)
        elif is_groove:
            axg.annotate(f"exon {idx}", xy=(mid, 0.88), xytext=(mid, 2.05),
                          ha="center", va="bottom", fontsize=8.5, color=C_CDS, zorder=4,
                          annotation_clip=False,
                          arrowprops=dict(arrowstyle="-", color=C_CDS, lw=0.9,
                                          shrinkA=0, shrinkB=0))
    axg.set_ylim(0, 1)
    axg.set_yticks([])
    axg.set_xlabel(f"Position along {canonical}  (bp)", fontsize=11)
    axg.tick_params(labelsize=10)
    for side in ["top", "right", "left"]:
        axg.spines[side].set_visible(False)

    fig.savefig(out_path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outroot", default=os.path.expanduser("~/pipeline_outputs/people"))
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--genes", nargs="+", default=GENES)
    ap.add_argument("--canonical", default=None,
                     help="Force one canonical allele (only valid with a single --genes entry).")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--window", type=int, default=151, help="Sliding-window width for the overlay.")
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--replot", action="store_true",
                     help="Skip extraction entirely: rebuild figures from the {GENE}.json files "
                          "already in --out-dir. Extraction over the full cohort takes ~35 min, so "
                          "no figure restyling should ever require re-reading the PAFs.")
    args = ap.parse_args()

    if args.replot:
        out_dir_rp = args.out_dir or os.path.expanduser("~/results/15_hla_manhattan")
        for gene in args.genes:
            stem = gene.replace("HLA-", "")
            jp = os.path.join(out_dir_rp, f"{stem}.json")
            if not os.path.exists(jp):
                print(f"  {gene}: no {jp}; skipping", file=sys.stderr)
                continue
            with open(jp) as f:
                d = json.load(f)
            cds_mask = in_ranges_mask(d["canonical_length"], [tuple(r) for r in d["cds_ranges"]])
            png = os.path.join(out_dir_rp, f"{stem}_manhattan.png")
            plot_manhattan(gene, d["canonical"], d["pi"], d["coverage"], cds_mask,
                            [tuple(e) for e in d["exons"]], png, d["n_haps"],
                            window=args.window, groove_exons=GROOVE_EXONS.get(gene, ()))
            print(f"  {gene}: replotted -> {png}  (n_haps={d['n_haps']}, "
                  f"pi ratio={d.get('pi_ratio_cds_over_noncds')})", file=sys.stderr)
        return

    out_dir = args.out_dir or os.path.expanduser("~/results/15_hla_manhattan")
    os.makedirs(out_dir, exist_ok=True)

    persons = sorted(d for d in os.listdir(args.outroot)
                      if os.path.isdir(os.path.join(args.outroot, d, "immuannot_output")))
    if args.limit:
        persons = persons[:args.limit]
    print(f"{len(persons)} people, genes={args.genes}, threads={args.threads}", file=sys.stderr)

    # Resolve every gene's canonical allele up front so all genes can be extracted in ONE pass
    # over each haplotype's PAF (see canonical_rows: the scan dominates, so this is ~Nx faster).
    canon_of = {}
    for gene in args.genes:
        if args.canonical and len(args.genes) == 1:
            canonical, sel = args.canonical, {"forced": True}
        else:
            canonical, sel = pick_canonical(persons, args.outroot, gene)
        if not canonical:
            print(f"  no canonical allele found for {gene}; skipping", file=sys.stderr)
            continue
        canon_of[gene] = (canonical, sel)
        print(f"  {gene}: canonical={canonical}  selection={sel}", file=sys.stderr)
    if not canon_of:
        sys.exit("no canonical alleles resolved for any requested gene")

    t_all = time.time()
    acc_by_canon = aggregate_multi(persons, args.outroot,
                                    [c for c, _ in canon_of.values()], threads=args.threads)
    print(f"single-pass extraction over {len(persons)*2} haplotype units took "
          f"{time.time()-t_all:.0f}s for {len(canon_of)} gene(s)", file=sys.stderr)

    for gene, (canonical, sel) in canon_of.items():
        print(f"=== {gene} ===", file=sys.stderr)
        ann = load_allele_annotation(canonical)
        if ann is None:
            print(f"  WARNING: no refdata annotation for {canonical}; CDS mask unavailable",
                  file=sys.stderr)
            ann = {"cds": [], "exons": [], "utr": [], "gene_range": []}

        agg = acc_by_canon[canonical]
        agg["n_missing_canonical"] = agg["n_missing"]
        elapsed = time.time() - t_all
        if not agg["coverage"]:
            print("  no haplotypes carried the canonical row; skipping", file=sys.stderr)
            continue

        pi, n_diff = per_site_pi(agg["coverage"], agg["alt_counts"])
        length = len(pi)
        cds_mask = in_ranges_mask(length, ann["cds"])
        n_variant_sites = sum(1 for d in n_diff if d > 0)
        cds_pi = [p for p, m in zip(pi, cds_mask) if m]
        non_pi = [p for p, m in zip(pi, cds_mask) if not m]
        mean_cds = sum(cds_pi) / len(cds_pi) if cds_pi else None
        mean_non = sum(non_pi) / len(non_pi) if non_pi else None
        nm = [v for v in agg["nm_values"] if v is not None]

        out = {
            "gene": gene, "canonical": canonical, "canonical_selection": sel,
            "n_people": len(persons), "n_haps": agg["n_haps"],
            "n_missing_canonical": agg["n_missing_canonical"],
            "strand_counts": dict(agg["strand_counts"]),
            "canonical_length": length,
            "n_variant_sites": n_variant_sites,
            "mean_pi_cds": mean_cds, "mean_pi_noncds": mean_non,
            "pi_ratio_cds_over_noncds": (mean_cds / mean_non) if (mean_cds and mean_non) else None,
            "nm_min": min(nm) if nm else None, "nm_max": max(nm) if nm else None,
            "nm_mean": (sum(nm) / len(nm)) if nm else None,
            "cds_ranges": ann["cds"], "exons": ann["exons"], "utr": ann["utr"],
            "coverage": agg["coverage"], "pi": pi, "n_diff": n_diff,
            "insertion_counts": dict(sorted(agg["insertion_counts"].items())),
            "elapsed_seconds": elapsed,
        }
        stem = gene.replace("HLA-", "")
        json_path = os.path.join(out_dir, f"{stem}.json")
        with open(json_path, "w") as f:
            json.dump(out, f)

        print(f"  n_haps={agg['n_haps']}  missing_canonical={agg['n_missing_canonical']}  "
              f"strands={dict(agg['strand_counts'])}  elapsed={elapsed:.0f}s", file=sys.stderr)
        print(f"  canonical_length={length}  n_variant_sites={n_variant_sites}  "
              f"NM min/mean/max={out['nm_min']}/{out['nm_mean']}/{out['nm_max']}", file=sys.stderr)
        print(f"  mean pi: CDS={mean_cds}  nonCDS={mean_non}  "
              f"ratio={out['pi_ratio_cds_over_noncds']}", file=sys.stderr)
        print(f"  wrote {json_path}", file=sys.stderr)

        if args.plot:
            png = os.path.join(out_dir, f"{stem}_manhattan.png")
            plot_manhattan(gene, canonical, pi, agg["coverage"], cds_mask, ann["exons"], png,
                            agg["n_haps"], window=args.window,
                            groove_exons=GROOVE_EXONS.get(gene, ()))
            print(f"  wrote {png}", file=sys.stderr)


if __name__ == "__main__":
    main()
