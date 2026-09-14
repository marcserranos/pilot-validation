#!/usr/bin/env python3
"""Bypass HLA nomenclature entirely and compare the ACTUAL assembled DNA sequence (Immuannot's own
`cds.fa.gz` per person/hap) between the 545 high-sharing relative pairs from
12_phasing_mendelian_validation.py, at every gene both people have a sequence for.

## Why this script exists

12_phasing_mendelian_validation.py's ~4.9%/5.6% "error rate" is measured on the *named* 2-field
allele (Immuannot's `consensus` column) -- and Immuannot's naming step (`consensusCall()`) collapses
several equally-good reference candidates down to one representative name when there's a tie. Marc's
critique, 2026-09-08: a real difference of 1-2 nucleotides can flip which named allele two people's
otherwise-near-identical sequences get collapsed to, making a trivial technical difference look like
a full categorical "mismatch" -- and conversely, a real single-base difference SHOULD register as its
own distinct identity, not get silently absorbed into "the same named allele." Named-allele comparison
answers neither question correctly. Comparing the raw assembled sequence content directly does.

A worked example (2026-09-08, real data) already showed both directions of this failure: HLA-DPB1
was called a "mismatch" by name but the raw sequences were BYTE-IDENTICAL (pure naming-collapse
artifact); HLA-B, W, DQB2, H, T, TAP2 remained genuinely different at the raw sequence level despite
looking name-adjacent in some cases. This script makes that check systematic across the whole cohort
instead of one example pair, and quantifies not just "same or different" but exactly how many bases
differ and whether those differences are scattered points or one contiguous block.

## Method

For each gene, per person, collect the SET of distinct observed sequences across both haplotypes (a
person carries 1 sequence if homozygous, up to 2 if heterozygous, occasionally more for a real
segmental-duplication copy number). Two people are called `raw_shared` at that gene if any of person
A's sequences exactly equals any of person B's sequences (a real parent-child transmission event
implies this must be true for the actual transmitted copy, exact byte-for-byte, modulo real technical
error). When no exact match exists, this script finds the closest pairing across all A-vs-B sequence
combinations using Python's difflib.SequenceMatcher (a standard, dependency-free approximate
alignment) and reports:

  - `n_diff_bases`: total bases involved in a difference (substituted, inserted, or deleted) in the
    best-aligned pairing -- the direct, literal answer to "how many bases actually differ."
  - `n_diff_blocks`: how many separate contiguous stretches those differing bases fall into. This is
    the "shape" of the difference: `n_diff_blocks == 1` with `n_diff_bases` small means one localized
    change (e.g. a single point substitution or a small indel) -- consistent with real, tightly
    localized biological variation or a single base-calling slip. `n_diff_blocks > 1` (differences
    scattered at multiple separate positions along the same sequence) is a different signature --
    more consistent with comparing two genuinely different alleles/haplotypes, or an assembly issue
    affecting an extended stretch, than with a single clean transmission event carrying one new
    variant.
  - `same_length`: whether the closest-matching pair has identical length (a pure substitution-style
    difference) or not (an indel/structural-style difference).

Usage (fixtures):
    python3 scripts/hla_popgen/17_raw_sequence_divergence.py \\
        --outroot /tmp/hla_fixtures/people \\
        --pair-classification /tmp/hla_fixtures/phasing_pair_classification.sample.tsv \\
        --table1 /tmp/hla_fixtures/hla_calls_rich.sample.tsv \\
        --out-dir /tmp/hla_fixtures/reports/17_raw_sequence_divergence

Real run (VM):
    python3 scripts/hla_popgen/17_raw_sequence_divergence.py
"""
import argparse
import difflib
import gzip
import hashlib
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _viz_common as vc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

DEFAULT_PAIR_CLASSIFICATION = os.path.expanduser("~/pipeline_outputs/phasing_pair_classification.tsv")
DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs/people")

_MIC_TAP = {"MICA", "MICB", "TAP1", "TAP2"}


# ---------------------------------------------------------------------------
# Raw CDS FASTA parsing
# ---------------------------------------------------------------------------
def parse_cds_fasta(path):
    """Header format (confirmed empirically 2026-09-08): '>{contig}_{gene}_{i} {splice_sites}' --
    the splice-site annotation after the first space is NOT part of the identifying key and must be
    dropped before applying the trailing-integer-index regex, or every header fails to match.
    Returns {gene_with_HLA_prefix: [(copy_index, seq), ...]}."""
    out = defaultdict(list)
    if not os.path.exists(path):
        return out
    key, chunks = None, []

    def flush():
        if key is None:
            return
        seq = "".join(chunks).upper()
        m = re.match(r"^(.*)_(\d+)$", key)
        if m:
            gene = m.group(1).split("_")[-1]
            out[gene].append((int(m.group(2)), seq))

    with gzip.open(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                key = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    flush()
    return out


def load_person_sequences(outroot, person_id, gene_bare_names):
    """{(hap, gene_bare): [seq, ...]} for one person, both haps, restricted to the requested genes."""
    target_names = {g: ("HLA-" + g if g not in _MIC_TAP else g) for g in gene_bare_names}
    result = {}
    for hap in ("hap1", "hap2"):
        path = os.path.join(outroot, person_id, "immuannot_output", hap, "cds.fa.gz")
        idx = parse_cds_fasta(path)
        for g, target in target_names.items():
            seqs = [seq for _, seq in idx.get(target, [])]
            if seqs:
                result[(hap, g)] = seqs
    return result


# ---------------------------------------------------------------------------
# Sequence comparison
# ---------------------------------------------------------------------------
def best_alignment_diff(seq_a, seq_b):
    """Returns (n_diff_bases, n_diff_blocks) describing how seq_a and seq_b differ.

    Equal-length sequences (the dominant case -- both are the same gene's CDS, so length differs
    only when there's a real indel) use a direct POSITION-ANCHORED comparison: count mismatched
    positions, and count contiguous runs of them as blocks. This is deliberate, not a simplification
    -- a real bug caught by this script's own fixture test: difflib.SequenceMatcher optimizes for the
    LONGEST matching block, which for repetitive sequences (homopolymer runs and tandem repeats are
    a documented, real feature of this project's own novel-allele QC, see NOVEL_LIT.md's
    "is_homopolymer_indel_only" flag) can align on the wrong block entirely -- a single true
    substitution in a run of repeated bases got reported as a full-length insert+delete (20 diff
    bases instead of 1) before this fix. Position-anchoring sidesteps that failure mode completely
    whenever no indel is actually possible (equal length).

    Only genuinely different-length pairs fall back to difflib.SequenceMatcher, since some kind of
    indel-aware alignment is unavoidable there; this is a real, acknowledged limitation for that
    minority case (repetitive regions can still misalign), not silently assumed to be exact."""
    if len(seq_a) == len(seq_b):
        mismatches = [i for i, (x, y) in enumerate(zip(seq_a, seq_b)) if x != y]
        if not mismatches:
            return 0, 0
        n_blocks = 1
        for prev, cur in zip(mismatches, mismatches[1:]):
            if cur != prev + 1:
                n_blocks += 1
        return len(mismatches), n_blocks

    sm = difflib.SequenceMatcher(None, seq_a, seq_b, autojunk=False)
    n_bases, n_blocks = 0, 0
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal":
            continue
        n_blocks += 1
        n_bases += max(i2 - i1, j2 - j1)
    return n_bases, n_blocks


def compare_gene(a_seqs, b_seqs):
    """a_seqs, b_seqs: lists of observed sequences (1-2+ per person). Returns a dict describing the
    closest possible pairing. `raw_shared=True` short-circuits (exact match exists, 0 diff)."""
    a_set, b_set = set(a_seqs), set(b_seqs)
    if a_set & b_set:
        return {"raw_shared": True, "n_diff_bases": 0, "n_diff_blocks": 0,
                "same_length": True, "len_a": len(next(iter(a_set & b_set))), "len_b": None}

    best = None
    for a in a_set:
        for b in b_set:
            n_bases, n_blocks = best_alignment_diff(a, b)
            if best is None or n_bases < best[0]:
                best = (n_bases, n_blocks, len(a), len(b))
    n_bases, n_blocks, len_a, len_b = best
    return {"raw_shared": False, "n_diff_bases": n_bases, "n_diff_blocks": n_blocks,
            "same_length": len_a == len_b, "len_a": len_a, "len_b": len_b}


# ---------------------------------------------------------------------------
# Aggregation across all pairs
# ---------------------------------------------------------------------------
def run_all(pairs_df, outroot, all_genes, gene_class):
    rows = []
    # Cache per-person sequences to avoid re-reading the same person's files across multiple pairs.
    seq_cache = {}

    def get_person_seqs(person_id):
        if person_id not in seq_cache:
            seq_cache[person_id] = load_person_sequences(outroot, person_id, all_genes)
        return seq_cache[person_id]

    for row in pairs_df.itertuples(index=False):
        a_seqs_by_gene = get_person_seqs(row.person_i)
        b_seqs_by_gene = get_person_seqs(row.person_j)
        for gene in all_genes:
            a_seqs = list(set(a_seqs_by_gene.get(("hap1", gene), []) +
                               a_seqs_by_gene.get(("hap2", gene), [])))
            b_seqs = list(set(b_seqs_by_gene.get(("hap1", gene), []) +
                               b_seqs_by_gene.get(("hap2", gene), [])))
            if not a_seqs or not b_seqs:
                continue
            cmp = compare_gene(a_seqs, b_seqs)
            rows.append({"gene": gene, "gene_class": gene_class.get(gene, "other"), **cmp})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_diff_bases_histogram(diffs, out_path):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    vals = diffs["n_diff_bases"].values
    bins = [0, 1, 2, 3, 4, 5, 10, 20, 50, 100, 250, 500, 1000, 2500, 5000]
    ax.hist(vals, bins=bins, color="#0072B2", edgecolor="white")
    ax.set_xscale("symlog")
    ax.set_xlabel("Bases differing from the closest possible match (log scale)")
    ax.set_ylabel("Number of gene-comparisons")
    ax.set_title("Raw-sequence divergence, among gene-comparisons with NO exact sequence match\n"
                 "(named-allele nomenclature bypassed entirely)")
    ax.spines[["top", "right"]].set_visible(False)
    vc.savefig(fig, out_path)


def plot_diff_blocks_histogram(diffs, out_path):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    counts = diffs["n_diff_blocks"].clip(upper=5)
    labels = ["1\n(one clustered\nchange)", "2", "3", "4", "5+\n(scattered)"]
    vals = [int((counts == i).sum()) for i in range(1, 5)] + [int((counts >= 5).sum())]
    ax.bar(labels, vals, color="#D55E00")
    for i, v in enumerate(vals):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Number of gene-comparisons")
    ax.set_title("Shape of raw-sequence differences: one clustered block vs scattered\n"
                 "(among gene-comparisons with no exact sequence match)")
    ax.spines[["top", "right"]].set_visible(False)
    vc.savefig(fig, out_path)


def plot_shared_rate_by_class(df, out_path):
    summary = df.groupby("gene_class")["raw_shared"].agg(["sum", "count"]).reset_index()
    summary["rate"] = 1 - summary["sum"] / summary["count"]
    p, lo, hi = vc.wilson_ci((summary["count"] - summary["sum"]).values, summary["count"].values)
    summary["p"], summary["lo"], summary["hi"] = p, lo, hi
    summary = summary.sort_values("p", ascending=False)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(summary))
    ax.bar(x, summary["p"] * 100, color="#009E73")
    yerr_lo = np.clip((summary["p"] - summary["lo"]) * 100, 0, None)
    yerr_hi = np.clip((summary["hi"] - summary["p"]) * 100, 0, None)
    ax.errorbar(x, summary["p"] * 100, yerr=[yerr_lo, yerr_hi], fmt="none", ecolor="black",
                elinewidth=0.8, capsize=2)
    ax.set_xticks(x)
    ax.set_xticklabels(summary["gene_class"], rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Raw-sequence divergence rate (%)")
    ax.set_title("Raw-sequence (nomenclature-bypassed) divergence rate by gene class")
    ax.spines[["top", "right"]].set_visible(False)
    vc.savefig(fig, out_path)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def write_report(df, out_dir):
    total = len(df)
    n_shared = int(df["raw_shared"].sum())
    n_diff = total - n_shared
    p, lo, hi = vc.wilson_ci(np.array([n_diff]), np.array([total]))

    lines = ["# Raw-sequence divergence -- nomenclature-bypassed phasing/typing validation\n"]
    lines.append(f"Total gene-comparisons (545 high-sharing pairs x genes both sides have a "
                 f"sequence for): **{total}**\n")
    lines.append(f"Exact sequence match (raw_shared): **{n_shared}** "
                 f"({100*n_shared/total:.3f}%)\n")
    lines.append(f"No exact match, closest pairing still differs: **{n_diff}** "
                 f"({100*p[0]:.3f}%, 95% CI {100*lo[0]:.3f}-{100*hi[0]:.3f}%)\n")

    diffs = df[~df["raw_shared"]]
    if not diffs.empty:
        lines.append("\n## How many bases actually differ (among the non-matching comparisons)\n")
        desc = diffs["n_diff_bases"].describe(percentiles=[.25, .5, .75, .9])
        lines.append("| stat | value |\n|---|---|\n" +
                     "\n".join(f"| {k} | {v:.1f} |" for k, v in desc.items()))
        lines.append(f"\n\n- Same length (pure substitution-style, no indel): "
                     f"**{int(diffs['same_length'].sum())}/{len(diffs)}** "
                     f"({100*diffs['same_length'].mean():.1f}%)\n")
        lines.append(f"- 1-3 differing bases specifically (near-miss / point-level): "
                     f"**{int((diffs['n_diff_bases']<=3).sum())}/{len(diffs)}** "
                     f"({100*(diffs['n_diff_bases']<=3).mean():.1f}%)\n")

        lines.append("\n## Shape: one clustered block vs scattered across the sequence\n")
        block_counts = diffs["n_diff_blocks"].value_counts().sort_index()
        lines.append("| n_diff_blocks | n_gene_comparisons | meaning |\n|---|---|---|")
        meanings = {1: "one localized change (point substitution or small indel)",
                    2: "two separate differing regions"}
        for k, v in block_counts.items():
            m = meanings.get(k, "scattered across >=3 separate regions")
            lines.append(f"| {k} | {v} | {m} |")

        lines.append("\n## By gene class\n")
        lines.append("| gene_class | n_compared | n_diverged | divergence_rate_% |\n"
                     "|---|---|---|---|")
        by_class = df.groupby("gene_class").agg(
            n_compared=("raw_shared", "size"), n_diverged=("raw_shared", lambda s: (~s).sum()))
        by_class["rate"] = 100 * by_class["n_diverged"] / by_class["n_compared"]
        for gc, row in by_class.sort_values("rate", ascending=False).iterrows():
            lines.append(f"| {gc} | {int(row.n_compared)} | {int(row.n_diverged)} | "
                         f"{row.rate:.3f} |")

    vc.write_report(os.path.join(out_dir, "raw_sequence_divergence_report.md"), lines)
    df.drop(columns=[]).to_csv(os.path.join(out_dir, "raw_sequence_divergence_detail.tsv"),
                                sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=os.path.join(vc.DEFAULT_OUTROOT, "hla_calls_rich.tsv"))
    ap.add_argument("--pair-classification", default=DEFAULT_PAIR_CLASSIFICATION)
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT,
                    help="Root of per-person immuannot_output directories (contains cds.fa.gz).")
    ap.add_argument("--out-dir", default=os.path.join(vc.DEFAULT_REPORT_ROOT,
                                                       "17_raw_sequence_divergence"))
    args = ap.parse_args()

    print("Loading Table 1 (for gene list + gene_class)...", file=sys.stderr)
    table1 = vc.load_table1(args.table1)
    gene_class = dict(zip(table1["gene_bare"], table1["gene_class"]))
    all_genes = sorted(table1["gene_bare"].dropna().unique())

    print("Loading pair classification (545 high-sharing pairs expected)...", file=sys.stderr)
    pairs_df = pd.read_csv(args.pair_classification, sep="\t", dtype={"person_i": str, "person_j": str})
    pairs_df = pairs_df[pairs_df["pair_class"] == "high_sharing"]
    print(f"  {len(pairs_df)} high-sharing pairs.", file=sys.stderr)

    print("Comparing raw sequences gene-by-gene for every pair (reading cds.fa.gz)...", file=sys.stderr)
    df = run_all(pairs_df, args.outroot, all_genes, gene_class)
    print(f"  {len(df)} gene-comparisons with a sequence on both sides.", file=sys.stderr)
    print(f"  raw_shared: {int(df['raw_shared'].sum())} / {len(df)}", file=sys.stderr)

    os.makedirs(args.out_dir, exist_ok=True)
    diffs = df[~df["raw_shared"]]
    if not diffs.empty:
        plot_diff_bases_histogram(diffs, os.path.join(args.out_dir, "diff_bases_histogram.png"))
        plot_diff_blocks_histogram(diffs, os.path.join(args.out_dir, "diff_blocks_histogram.png"))
    plot_shared_rate_by_class(df, os.path.join(args.out_dir, "shared_rate_by_class.png"))
    write_report(df, args.out_dir)
    print(f"Done. Wrote outputs to {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
