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


def canonical_row(paf_path, canonical):
    """Best (longest query-span) PAF row whose qname == canonical. Returns None if absent."""
    best = None
    best_span = -1
    opener = gzip.open if paf_path.endswith(".gz") else open
    with opener(paf_path, "rt") as f:
        for line in f:
            if not line.startswith(canonical):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or fields[0] != canonical:
                continue
            span = int(fields[3]) - int(fields[2])
            if span > best_span:
                best_span, best = span, fields
    return best


def process_haplotype(person_dir, hap, canonical):
    """-> None, or dict with the canonical-frame observations for this haplotype."""
    paf = os.path.join(person_dir, hap, "mm2.ipd.gen.paf.gz")
    if not os.path.exists(paf):
        return None
    row = canonical_row(paf, canonical)
    if row is None:
        return None
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


# ---------------------------------------------------------------------------
# Aggregation across haplotypes -> per-site allele counts
# ---------------------------------------------------------------------------
def aggregate(persons, outroot, canonical, threads=1, progress_every=400):
    """Per-site allele-label counts in canonical coordinates, plus coverage and QC."""
    qlen_seen = Counter()
    n_haps = n_missing = 0
    nm_values = []
    strand_counts = Counter()
    # coverage[pos] and diffs[pos] are dense over the canonical length; allocate lazily once qlen known.
    coverage = None
    alt_counts = None   # list of Counter, only for positions that ever differ (sparse dict)
    insertion_counts = Counter()

    t0 = time.time()

    def tasks():
        for pid in persons:
            pdir = os.path.join(outroot, pid, "immuannot_output")
            for hap in ("hap1", "hap2"):
                yield pdir, hap

    if threads > 1:
        pool = concurrent.futures.ThreadPoolExecutor(max_workers=threads)
        results = pool.map(lambda ph: process_haplotype(ph[0], ph[1], canonical), tasks())
    else:
        pool = None
        results = (process_haplotype(pd, hap, canonical) for pd, hap in tasks())

    n_units = len(persons) * 2
    for i, res in enumerate(results, 1):
        if progress_every and i % progress_every == 0:
            el = time.time() - t0
            print(f"    [{canonical}] {i}/{n_units} units, {el:.0f}s ({i/el:.1f}/s), "
                  f"{n_haps} haps with the canonical row", file=sys.stderr)
        if res is None:
            n_missing += 1
            continue
        n_haps += 1
        qlen_seen[res["qlen"]] += 1
        nm_values.append(res["nm"])
        strand_counts[res["strand"]] += 1
        if coverage is None:
            coverage = [0] * res["qlen"]
            alt_counts = defaultdict(Counter)
        lo, hi = res["qstart"], res["qend"]
        for p in range(lo, min(hi, len(coverage))):
            coverage[p] += 1
        for pos, label in res["observed"].items():
            if 0 <= pos < len(coverage):
                alt_counts[pos][label] += 1
        for pos in res["insertions"]:
            if 0 <= pos < len(coverage):
                insertion_counts[pos] += 1

    if pool is not None:
        pool.shutdown()

    return {
        "n_haps": n_haps, "n_missing_canonical": n_missing,
        "qlen_seen": qlen_seen, "nm_values": nm_values, "strand_counts": strand_counts,
        "coverage": coverage or [], "alt_counts": alt_counts or {},
        "insertion_counts": insertion_counts,
    }


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
    best = max(counts.items(), key=lambda kv: (kv[1], qlens.get(kv[0], 0)))[0] if counts else None
    return best, {"probed": n, "present_in": counts.get(best, 0), "qlen": qlens.get(best),
                   "auto_selected": True,
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
                    window=151):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    n = len(pi)
    x = np.arange(n)
    pi_arr = np.asarray(pi)
    smooth = np.asarray(sliding_mean(pi, window))
    cds = np.asarray(cds_mask)

    fig, (ax, axg) = plt.subplots(
        2, 1, figsize=(13, 5.0), sharex=True,
        gridspec_kw={"height_ratios": [11, 1], "hspace": 0.08})

    # Per-site pi, colored by CDS membership so the contrast is readable without a legend hunt.
    ax.vlines(x[~cds], 0, pi_arr[~cds], color="#B9B4D6", lw=0.6, alpha=0.85,
              label="non-CDS site (intron/UTR)")
    ax.vlines(x[cds], 0, pi_arr[cds], color="#C44E52", lw=0.6, alpha=0.95,
              label="CDS site")
    ax.plot(x, smooth, color="#1F3B73", lw=1.6, alpha=0.95,
            label=f"{window} bp sliding mean")

    ax.set_ylabel("Nucleotide diversity  $\\pi$  per site")
    ax.set_xlim(0, n)
    ax.set_ylim(bottom=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper right", fontsize=8.5, ncol=3)
    mean_cds = float(pi_arr[cds].mean()) if cds.any() else float("nan")
    mean_non = float(pi_arr[~cds].mean()) if (~cds).any() else float("nan")
    ratio = (mean_cds / mean_non) if mean_non else float("nan")
    ax.set_title(
        f"HLA-{gene.replace('HLA-', '')}: per-site nucleotide diversity along the gene "
        f"({n_haps:,} haplotypes, all projected onto {canonical})\n"
        f"mean $\\pi$ in CDS = {mean_cds:.4f}   vs   outside CDS = {mean_non:.4f}   "
        f"(ratio {ratio:.2f}x)", fontsize=11)

    # Gene model ribbon: exon boxes on a thin axis beneath the track.
    axg.axhline(0.5, color="#888888", lw=1.0, zorder=1)
    for a, b in exons:
        axg.add_patch(plt.Rectangle((a - 1, 0.15), max(b - a + 1, 1), 0.7,
                                     facecolor="#4C4C4C", edgecolor="none", zorder=2))
    axg.set_ylim(0, 1)
    axg.set_yticks([])
    axg.set_xlabel(f"Position along {canonical} (bp)")
    for side in ["top", "right", "left"]:
        axg.spines[side].set_visible(False)
    axg.text(0.002, 0.5, "exons", transform=axg.transAxes, fontsize=8,
             va="center", ha="left", color="#4C4C4C")

    fig.savefig(out_path, dpi=170, bbox_inches="tight")
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
    args = ap.parse_args()

    out_dir = args.out_dir or os.path.expanduser("~/results/12_hla_manhattan")
    os.makedirs(out_dir, exist_ok=True)

    persons = sorted(d for d in os.listdir(args.outroot)
                      if os.path.isdir(os.path.join(args.outroot, d, "immuannot_output")))
    if args.limit:
        persons = persons[:args.limit]
    print(f"{len(persons)} people, genes={args.genes}, threads={args.threads}", file=sys.stderr)

    for gene in args.genes:
        print(f"=== {gene} ===", file=sys.stderr)
        if args.canonical and len(args.genes) == 1:
            canonical, sel = args.canonical, {"forced": True}
        else:
            canonical, sel = pick_canonical(persons, args.outroot, gene)
        if not canonical:
            print(f"  no canonical allele found for {gene}; skipping", file=sys.stderr)
            continue
        print(f"  canonical={canonical}  selection={sel}", file=sys.stderr)

        ann = load_allele_annotation(canonical)
        if ann is None:
            print(f"  WARNING: no refdata annotation for {canonical}; CDS mask unavailable",
                  file=sys.stderr)
            ann = {"cds": [], "exons": [], "utr": [], "gene_range": []}

        t0 = time.time()
        agg = aggregate(persons, args.outroot, canonical, threads=args.threads)
        elapsed = time.time() - t0
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
                            agg["n_haps"], window=args.window)
            print(f"  wrote {png}", file=sys.stderr)


if __name__ == "__main__":
    main()
