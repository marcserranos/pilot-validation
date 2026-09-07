#!/usr/bin/env python3
"""Per-position variant density across the FULL gene body (introns+UTR+CDS), compared inside vs.
outside the CDS -- a selection-pressure track, not just "which codon is this one variant at"
(that's `09_mutation_topology.py`, CDS-only). See scripts/hla_popgen/NEEDLE_VIEW_BRIEF.md for the
full design rationale and the live-VM validation this script's join logic is built on.

## Hypothesis (read before "fixing" any CDS<flanking asymmetry you see)

The obvious framing -- "CDS under purifying selection -> LOWER diversity inside CDS" -- is
backwards for classical HLA class I/II peptide-binding-groove exons: these are the textbook case of
DIVERSIFYING/BALANCING selection (Hughes & Nei; research/NOVEL_LIT.md Section 4's Ewens-Watterson
neutrality-rejection results). Expect CDS diversity to be >= flanking, concentrated at the groove-
encoding exons (exon 2-3 for class I, exon 2 for class II) -- a non-monotonic, position-dependent
signal, not a uniform CDS-vs-flanking mean shift.

## Data source and the join that makes cross-person aggregation possible

Per haplotype: `mm2.ipd.gen.paf.gz` is minimap2 `-cx asm5 --cs` output, **query = the IPD/IMGT
reference gene FASTA, target = the person's contig** (confirmed live on the VM, 2026-09-07 --
NEEDLE_VIEW_BRIEF.md). Walking the `cs` string's QUERY-side offset gives a position in the shared
reference-gene coordinate frame -- the only way ~12,000 people's differently-coordinate contigs can
be aggregated onto one track at all.

This PAF file holds the FULL CANDIDATE UNIVERSE, not just the winning call (confirmed live: one
hap1 file had 32,616 candidate rows across every HLA/immune gene tested, e.g. 5 different
`HLA-A*01:01:01:0x` candidates mapped to the same contig region). The join that isolates the
accepted call: the same `hap{N}.gtf.gz`'s `gene` row carries `template_allele` (the winning
reference allele name) and `template_distance` (0 = perfect match). Filtering the PAF for
`qname == template_allele` returns exactly one row -- confirmed live (`HLA-A*26:01:01:01`,
`NM:i:0`, `cs:Z::3517`, matching `template_distance 0` in the GTF).

## CDS-boundary translation (the second, less obvious join)

The GTF's `CDS` rows give exon boundaries in CONTIG (target) coordinates -- not the reference-gene
(query) coordinate frame the diversity track lives in. This script builds a target->query
breakpoint map by walking the SAME cs string that produced the variant track (every `:N`/`*xy` op
advances query and target in lockstep; only indels break that lockstep), then looks up each CDS
row's contig-space boundary in that map to get its reference-gene-space equivalent. Boundary
estimates are collected per-haplotype and reduced to a per-gene consensus (median) with the spread
reported as a QC diagnostic -- large spread would flag alignment-quality-dependent boundary
instability worth investigating before trusting the CDS shading.

Usage (prototype, small sample):
    python3 scripts/hla_popgen/11_gene_diversity_track.py --outroot ~/pipeline_outputs/people \\
        --limit 50 --out-dir /tmp/diversity_track_prototype

Real run (VM, full cohort): python3 scripts/hla_popgen/11_gene_diversity_track.py
"""
import argparse
import gzip
import json
import os
import re
import sys
import time
from bisect import bisect_right
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module  # noqa: E402
extract_rich = import_module("01_extract_rich")  # noqa: E402  (parse_attrs, clean_gene_name)

GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1"]  # first-pass set: single-copy, best-annotated,
# strongest expected diversifying-selection signal at the groove (NEEDLE_VIEW_BRIEF.md).

CS_TOKEN_RE = re.compile(r":\d+|\*[a-z]{2}|[+-][a-z]+", re.IGNORECASE)
PAF_QNAME_COL, PAF_TSTART_COL, PAF_TEND_COL, PAF_TNAME_COL = 0, 7, 8, 5


# ---------------------------------------------------------------------------
# GTF: template_allele + CDS boundaries, per gene, contig-space
# ---------------------------------------------------------------------------
def parse_gtf_for_genes(path, genes):
    """One hap{N}.gtf.gz -> {gene: {"template_allele", "gene_start", "gene_end", "contig",
    "template_distance", "cds_ranges": [(start, end), ...]}} (contig-space, 1-based inclusive,
    GTF convention). Only classical, single-copy genes in `genes` with copy_index==1 are kept --
    multi-copy disambiguation is out of scope for this first pass (NEEDLE_VIEW_BRIEF.md)."""
    out = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            contig, _src, feature, start, end, _score, _strand, _frame, attr_str = fields[:9]
            if feature not in ("gene", "CDS"):
                continue
            attrs = extract_rich.parse_attrs(attr_str)
            gene_id = attrs.get("gene_id")
            if gene_id is None:
                continue
            base_id, copy_index = extract_rich.gene_id_copy_index(gene_id)
            gene_name_raw = attrs.get("gene_name")
            if gene_name_raw is None:
                continue
            gene, _c4 = extract_rich.clean_gene_name(gene_name_raw)
            if gene not in genes or copy_index != 1:
                continue
            if feature == "gene":
                out.setdefault(gene, {})["gene_start"] = int(start)
                out[gene]["gene_end"] = int(end)
                out[gene]["contig"] = contig
                out[gene]["template_allele"] = attrs.get("template_allele")
                td = attrs.get("template_distance")
                out[gene]["template_distance"] = int(td) if td is not None else None
            elif feature == "CDS":
                out.setdefault(gene, {}).setdefault("cds_ranges", []).append((int(start), int(end)))
    for gene in out:
        out[gene].setdefault("cds_ranges", []).sort()
    return out


# ---------------------------------------------------------------------------
# PAF: find the one accepted row for a gene, per hap
# ---------------------------------------------------------------------------
def find_accepted_paf_row(paf_path, template_allele, contig, gene_start, gene_end):
    """Streams the (large, ~30k-row) gene-level PAF once and returns the single row whose qname
    matches template_allele and whose target region overlaps the GTF gene span -- the
    candidate-vs-accepted join confirmed live on the VM (NEEDLE_VIEW_BRIEF.md). Returns None if
    template_allele is missing/unmatched (fails toward "no data" rather than a guessed row)."""
    if not template_allele:
        return None
    opener = gzip.open if paf_path.endswith(".gz") else open
    with opener(paf_path, "rt") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or fields[PAF_QNAME_COL] != template_allele:
                continue
            if fields[PAF_TNAME_COL] != contig:
                continue
            tstart, tend = int(fields[PAF_TSTART_COL]), int(fields[PAF_TEND_COL])
            # PAF target coords are 0-based half-open; GTF is 1-based inclusive.
            if tend < gene_start - 1 or tstart > gene_end:
                continue
            return fields
    return None


def get_tag(fields, tag):
    prefix = tag + ":"
    for f in fields[12:]:
        if f.startswith(prefix):
            return f.split(":", 2)[2]
    return None


# ---------------------------------------------------------------------------
# cs-string walk: variant events (query-space) + target<->query breakpoint map
# ---------------------------------------------------------------------------
def walk_cs(cs_string):
    """One pass over the cs string. Returns (variant_positions, breakpoints):
      variant_positions: sorted list of 0-based query offsets touched by a substitution or indel
        (an indel's query-side span is fully included: every reference position it disturbs, not
        just its start).
      breakpoints: list of (t_off, q_off) at the START of every op, 0-based, target-space and
        query-space respectively -- used to translate CDS boundaries (target-space) into query-
        space by nearest-preceding-breakpoint lookup.
    query = reference gene (advances on ':', '*', '+'), target = contig (advances on ':', '*', '-').
    """
    t_off = q_off = 0
    variants = []
    breakpoints = [(0, 0)]
    for tok in CS_TOKEN_RE.findall(cs_string):
        if tok.startswith(":"):
            n = int(tok[1:])
            t_off += n
            q_off += n
        elif tok.startswith("*"):
            variants.append(q_off)
            t_off += 1
            q_off += 1
        elif tok[0] == "+":
            # insertion: query (reference) has extra bases the contig doesn't -- a deletion from
            # the person's perspective. Query advances; target doesn't.
            n = len(tok) - 1
            variants.extend(range(q_off, q_off + n))
            q_off += n
        elif tok[0] == "-":
            # deletion: target (contig) has extra bases the query (reference) doesn't -- an
            # insertion from the person's perspective. Target advances; query doesn't -- record
            # one tick at the current query offset (the anchor point) since it's still evidence of
            # variation at that reference position, per 09_mutation_topology.py's convention.
            n = len(tok) - 1
            variants.append(q_off)
            t_off += n
        breakpoints.append((t_off, q_off))
    return variants, breakpoints


def target_to_query(breakpoints, t_target):
    """Nearest-preceding-breakpoint lookup: target offset (0-based) -> query offset (0-based)."""
    t_list = [b[0] for b in breakpoints]
    i = bisect_right(t_list, t_target) - 1
    i = max(0, min(i, len(breakpoints) - 1))
    t_b, q_b = breakpoints[i]
    return q_b + (t_target - t_b)  # within a matched run this is exact; across an indel it's an
    # approximation anchored to the nearest run, which is the best a breakpoint map can do.


# ---------------------------------------------------------------------------
# Per-person, per-gene, per-hap processing
# ---------------------------------------------------------------------------
def process_person_gene_hap(person_dir, hap, gene):
    """Returns None (no data / no accepted row) or a dict: variant_counter (Counter of query-space
    0-based positions touched by a variant in THIS haplotype -- caller sums across haplotypes),
    q_span (0, qlen) covered by this haplotype's alignment, cds_query_ranges (this haplotype's own
    CDS boundaries translated into query space), nm (alignment edit distance, for QC)."""
    gtf_path = os.path.join(person_dir, f"{hap}.gtf.gz")
    paf_path = os.path.join(person_dir, hap, "mm2.ipd.gen.paf.gz")
    if not (os.path.exists(gtf_path) and os.path.exists(paf_path)):
        return None
    gtf_info = parse_gtf_for_genes(gtf_path, [gene]).get(gene)
    if not gtf_info or not gtf_info.get("template_allele"):
        return None
    row = find_accepted_paf_row(paf_path, gtf_info["template_allele"], gtf_info["contig"],
                                 gtf_info["gene_start"], gtf_info["gene_end"])
    if row is None:
        return None
    qlen = int(row[1])
    cs = get_tag(row, "cs")
    nm = get_tag(row, "NM")
    if cs is None:
        return None
    variants, breakpoints = walk_cs(cs)
    tstart = int(row[PAF_TSTART_COL])
    cds_query_ranges = []
    for c_start, c_end in gtf_info["cds_ranges"]:
        # GTF 1-based inclusive contig coords -> 0-based offset from this alignment's tstart.
        t0, t1 = (c_start - 1) - tstart, c_end - tstart
        if t1 < 0 or t0 > (int(row[PAF_TEND_COL]) - tstart):
            continue  # CDS segment outside this alignment's covered target range
        q0, q1 = target_to_query(breakpoints, max(t0, 0)), target_to_query(breakpoints, t1)
        cds_query_ranges.append((q0, q1))
    # cds_bp is the EXACT GTF-derived CDS length (c_end - c_start + 1, summed) -- not derived from
    # the target->query translation, so it carries none of that lookup's boundary imprecision. This
    # is what makes the CDS-vs-flanking DENSITY comparison (see aggregate_gene / main) robust even
    # though the raw query-offset positional track is not directly comparable across haplotypes
    # that matched different-length template alleles (see module docstring "Hypothesis" section --
    # confirmed live, 2026-09-07: HLA-DRB1 haplotypes in a 50-person prototype matched templates
    # ranging 10,850-16,110bp, a >5kb spread, entirely from intron-length polymorphism between
    # reference alleles). Density only needs per-haplotype (count, bp) pairs, never a shared axis.
    cds_bp = sum(c_end - c_start + 1 for c_start, c_end in gtf_info["cds_ranges"])
    n_in_cds = sum(1 for v in set(variants) if in_any_range(v, cds_query_ranges))
    return {
        "variants": variants, "qlen": qlen, "cds_query_ranges": cds_query_ranges,
        "nm": int(nm) if nm is not None else None, "cds_bp": cds_bp,
        "n_variants_in_cds": n_in_cds, "n_variants_total": len(set(variants)),
    }


def in_any_range(pos, ranges):
    return any(lo <= pos <= hi for lo, hi in ranges)


def aggregate_gene(persons, outroot, gene, progress_every=200):
    """Streams all persons for one gene; returns {position: {"n_obs": int, "n_haps": int}} plus
    diagnostics: n_haps_used, qlen (reference gene length, should be ~constant across haps since
    they align to the same/similar reference alleles), cds_query_ranges (list of per-hap boundary
    estimates, for the consensus + spread QC)."""
    variant_counts = Counter()
    n_haps = 0
    qlens = Counter()
    all_cds_ranges = []
    nm_values = []
    total_cds_bp = total_qlen_bp = total_variants_in_cds = total_variants_all = 0
    t0 = time.time()
    for i, person_id in enumerate(persons, 1):
        person_dir = os.path.join(outroot, person_id, "immuannot_output")
        for hap in ("hap1", "hap2"):
            result = process_person_gene_hap(person_dir, hap, gene)
            if result is None:
                continue
            n_haps += 1
            qlens[result["qlen"]] += 1
            nm_values.append(result["nm"])
            for pos in set(result["variants"]):  # dedupe within-haplotype double-hits at one pos
                variant_counts[pos] += 1
            if result["cds_query_ranges"]:
                all_cds_ranges.append(result["cds_query_ranges"])
            # Density accounting (robust to the cross-haplotype coordinate-frame problem -- see
            # process_person_gene_hap's comment -- because it never compares raw positions across
            # haplotypes, only sums per-haplotype (count, bp) pairs).
            total_cds_bp += result["cds_bp"]
            total_qlen_bp += result["qlen"]
            total_variants_in_cds += result["n_variants_in_cds"]
            total_variants_all += result["n_variants_total"]
        if progress_every and i % progress_every == 0:
            elapsed = time.time() - t0
            print(f"    [{gene}] {i}/{len(persons)} people, {elapsed:.0f}s elapsed "
                  f"({i/elapsed:.1f} people/sec), {n_haps} haps with data so far", file=sys.stderr)
    return {
        "variant_counts": variant_counts, "n_haps": n_haps, "qlens": qlens,
        "all_cds_ranges": all_cds_ranges, "nm_values": nm_values,
        "total_cds_bp": total_cds_bp, "total_qlen_bp": total_qlen_bp,
        "total_variants_in_cds": total_variants_in_cds, "total_variants_all": total_variants_all,
    }


def cds_vs_flanking_density(total_cds_bp, total_qlen_bp, total_variants_in_cds, total_variants_all):
    """The primary, coordinate-frame-robust result: variant density (per kb) inside vs. outside the
    CDS, plus an exact binomial test of whether variants are distributed disproportionately to CDS
    vs. its share of total gene length (the null: variants land uniformly at random along the gene
    regardless of CDS membership -- expected P(in CDS) = total_cds_bp / total_qlen_bp). This tests
    diversity CONCENTRATION, not raw magnitude, and needs no cross-haplotype position alignment."""
    total_noncds_bp = total_qlen_bp - total_cds_bp
    n_out = total_variants_all - total_variants_in_cds
    density_in = 1000.0 * total_variants_in_cds / total_cds_bp if total_cds_bp else None
    density_out = 1000.0 * n_out / total_noncds_bp if total_noncds_bp else None
    p_cds = total_cds_bp / total_qlen_bp if total_qlen_bp else None
    # Exact two-sided binomial test via scipy (which works in log-space internally) rather than a
    # hand-rolled sum of math.comb(n, i) * p**i * (1-p)**(n-i) terms -- that direct approach
    # OverflowErrors (int too large to convert to float) once n reaches the low thousands, which
    # HLA-DRB1's variant count does at cohort scale (confirmed live, 2026-09-07: crashed a 500-
    # person run). scipy.stats.binomtest is exact and doesn't have this failure mode.
    p_value = None
    if p_cds is not None and total_variants_all > 0:
        from scipy.stats import binomtest
        p_value = float(binomtest(total_variants_in_cds, total_variants_all, p_cds,
                                   alternative="two-sided").pvalue)
    return {
        "total_cds_bp": total_cds_bp, "total_noncds_bp": total_noncds_bp,
        "n_variants_in_cds": total_variants_in_cds, "n_variants_outside_cds": n_out,
        "density_per_kb_in_cds": density_in, "density_per_kb_outside_cds": density_out,
        "density_ratio_in_over_out": (density_in / density_out)
        if (density_in is not None and density_out) else None,
        "expected_p_in_cds_if_uniform": p_cds, "binomial_p_value": p_value,
    }


def consensus_cds_ranges(all_cds_ranges):
    """Per-boundary (start/end of each CDS segment, matched by rank) median across haplotypes, plus
    spread (max-min) as a QC diagnostic. Returns (consensus_ranges, spread_per_boundary)."""
    if not all_cds_ranges:
        return [], []
    n_segments = Counter(len(r) for r in all_cds_ranges).most_common(1)[0][0]
    usable = [r for r in all_cds_ranges if len(r) == n_segments]
    consensus, spreads = [], []
    for seg_idx in range(n_segments):
        starts = sorted(r[seg_idx][0] for r in usable)
        ends = sorted(r[seg_idx][1] for r in usable)
        mid = len(starts) // 2
        consensus.append((starts[mid], ends[mid]))
        spreads.append((starts[-1] - starts[0], ends[-1] - ends[0]))
    return consensus, spreads


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outroot", default=os.path.expanduser("~/pipeline_outputs/people"))
    ap.add_argument("--limit", type=int, default=None, help="Cap number of people (prototype runs).")
    ap.add_argument("--persons-file", default=None,
                     help="Optional file, one person_id per line, to use an exact sample (e.g. the "
                          "existing 50-person recon set) instead of the first --limit by sort order.")
    ap.add_argument("--genes", nargs="+", default=GENES)
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args()

    out_dir = args.out_dir or os.path.expanduser("~/results/11_gene_diversity_track")
    os.makedirs(out_dir, exist_ok=True)

    if args.persons_file:
        with open(args.persons_file) as f:
            persons = [l.strip() for l in f if l.strip()]
    else:
        persons = sorted(
            d for d in os.listdir(args.outroot)
            if os.path.isdir(os.path.join(args.outroot, d, "immuannot_output"))
        )
    if args.limit:
        persons = persons[:args.limit]
    print(f"{len(persons)} people to process, genes={args.genes}", file=sys.stderr)

    for gene in args.genes:
        print(f"=== {gene} ===", file=sys.stderr)
        result = aggregate_gene(persons, args.outroot, gene)
        consensus, spreads = consensus_cds_ranges(result["all_cds_ranges"])
        qlen_summary = result["qlens"].most_common()
        nm = result["nm_values"]
        out = {
            "gene": gene, "n_people": len(persons), "n_haps_with_data": result["n_haps"],
            "qlen_distribution": qlen_summary,
            "cds_query_ranges_consensus": consensus,
            "cds_query_ranges_spread": spreads,
            "n_cds_boundary_estimates": len(result["all_cds_ranges"]),
            "nm_min": min(nm) if nm else None, "nm_max": max(nm) if nm else None,
            "nm_mean": (sum(nm) / len(nm)) if nm else None,
            "variant_counts": dict(sorted(result["variant_counts"].items())),
            "cds_vs_flanking_density": cds_vs_flanking_density(
                result["total_cds_bp"], result["total_qlen_bp"],
                result["total_variants_in_cds"], result["total_variants_all"]),
        }
        out_path = os.path.join(out_dir, f"{gene.replace('HLA-', '')}.json")
        with open(out_path, "w") as f:
            json.dump(out, f, indent=1)
        d = out["cds_vs_flanking_density"]
        print(f"  n_haps_with_data={result['n_haps']}  qlen_top={qlen_summary[:3]}  "
              f"n_variant_positions={len(result['variant_counts'])}  "
              f"cds_consensus={consensus}  cds_spread={spreads}  "
              f"NM min/mean/max={out['nm_min']}/{out['nm_mean']}/{out['nm_max']}", file=sys.stderr)
        print(f"  CDS-vs-flanking density (per kb): in_cds={d['density_per_kb_in_cds']}  "
              f"outside_cds={d['density_per_kb_outside_cds']}  "
              f"ratio_in_over_out={d['density_ratio_in_over_out']}  "
              f"n_in/n_out={d['n_variants_in_cds']}/{d['n_variants_outside_cds']}  "
              f"expected_p_in_cds_if_uniform={d['expected_p_in_cds_if_uniform']}  "
              f"binomial_p={d['binomial_p_value']}", file=sys.stderr)
        print(f"  wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
