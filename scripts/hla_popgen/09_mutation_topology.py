#!/usr/bin/env python3
"""Where, along the CDS, do novel differences actually land? A gene-faceted lollipop/needle plot
(cBioPortal MutationMapper / ProteinPaint / Bioconductor trackViewer idiom, per the 2026-09
visualization-literature review in research/VIZ_LIT.md's addendum): one vertical stem per distinct
codon position, height = number of haplotype-level observations of a novel difference at that codon
(a direct recurrence-weighted count, not a per-cluster count -- a position hit by many people's
independently-assembled haplotypes is a stronger "real hotspot" signal than a raw candidate count),
colored by `novelty_class`.

## Where the position comes from

`cds_mut` (Table 1, reference/IMMUANNOT_GTF_SPEC.md part A) is
`[best-matching reference allele]|[minimap2 cs string]|[aa diff]`, comma-joined for tied candidates.
The aa-diff field (`REFaa(REFcodon)<OBSaa(OBScodon)`) gives the codon TRIPLET SEQUENCE, not a
position -- there is no position number in that field. The position comes from the middle field
instead: the minimap2 `cs` short-form string, whose `:N` tokens are runs of N identical bases --
summing them gives the running nucleotide offset at each substitution (`*xy`) or indel (`+seq`/
`-seq`) event, which converts to a 1-based codon position as `nt_offset // 3 + 1`. This script
takes the FIRST (best) tied candidate's cs-string and extracts every event's codon position this
way (the identical token regex 03_novel_alleles.py's `is_homopolymer_indel_only()` already uses on
the same field, reused here for position rather than homopolymer-shape).

## Why this scales where a raw per-cluster scatter wouldn't (VIZ_LIT.md addendum)

At hundreds-to-thousands of distinct novel-allele clusters per gene, an unaggregated per-cluster
scatter becomes an unreadable smear the moment multiple clusters share a codon -- which for HLA is
the norm, not the exception, since polymorphism concentrates at a small number of
peptide-binding-groove positions (exons 2-3 for class I, exon 2 for class II; Bjorkman & Parham
1990). Aggregating to one mark per codon position, height = count, is the standard fix.

Usage (fixtures): python3 scripts/hla_popgen/09_mutation_topology.py --outroot /tmp/hla_fixtures \\
    --table1 /tmp/hla_fixtures/hla_calls_rich.sample.tsv --sample

Real run (VM): python3 scripts/hla_popgen/09_mutation_topology.py
"""
import argparse
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _viz_common as vc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

GENES = vc.CLASSICAL_GENES_BARE
NOVELTY_COLORS = {"protein_altering": "#C44E52", "synonymous": "#55A868", "beyond_cds": "#8172B2",
                   "undetermined": "#999999"}
# `cds_mut`'s aa-diff field is "REFaa(REFcodon)<OBSaa(OBScodon)" where (REFcodon)/(OBScodon) are the
# CODON TRIPLET SEQUENCE (e.g. "GGC"), not a position -- there is no position number in that field
# at all. The actual position lives in the MIDDLE pipe field, the minimap2 `cs` short-form string
# (reference/IMMUANNOT_GTF_SPEC.md part A/C; same field 03_novel_alleles.py's
# `is_homopolymer_indel_only()` parses with this identical token regex): ':N' tokens are runs of N
# identical bases, so summing them gives the running nucleotide OFFSET into the CDS at each event.
CS_TOKEN_RE = re.compile(r":\d+|\*[a-z]{2}|[+-][a-z]+", re.IGNORECASE)


def parse_codon_positions(cds_mut):
    """First (best) tied candidate's cs-string -> list of 1-based codon positions (nt_offset // 3
    + 1) for every substitution/indel event. Returns [] for missing/unparseable input -- fails
    toward "no position data" rather than a guessed one."""
    if not isinstance(cds_mut, str) or not cds_mut.strip():
        return []
    first_candidate = cds_mut.split(",")[0]
    parts = first_candidate.split("|")
    if len(parts) < 2:
        return []
    cs = parts[1]
    tokens = CS_TOKEN_RE.findall(cs)
    nt_pos = 0
    codons = []
    for tok in tokens:
        if tok.startswith(":"):
            nt_pos += int(tok[1:])
        elif tok.startswith("*"):
            codons.append(nt_pos // 3 + 1)
            nt_pos += 1
        elif tok[0] in "+-":
            codons.append(nt_pos // 3 + 1)
            if tok[0] == "-":
                nt_pos += len(tok) - 1  # deletion consumes reference bases; insertion doesn't
    return codons


def build_position_counts(table1, genes):
    """Returns {gene: Counter({(codon_pos, novelty_class): n_observations})} and per-gene max
    codon seen (for axis scaling)."""
    sub = table1[table1["gene_bare"].isin(genes) & table1["is_novel"]]
    counts = defaultdict(Counter)
    for gene, cds_mut, novelty_class in zip(sub["gene_bare"], sub["cds_mut"], sub["novelty_class"]):
        positions = parse_codon_positions(cds_mut)
        nc = novelty_class if isinstance(novelty_class, str) and novelty_class in NOVELTY_COLORS \
            else "undetermined"
        for pos in positions:
            counts[gene][(pos, nc)] += 1
    return counts


def plot_lollipops(counts, genes, out_path, title):
    fig, axes = plt.subplots(2, 4, figsize=(22, 9))
    for ax, gene in zip(axes.flat, genes):
        gene_counts = counts.get(gene)
        if not gene_counts:
            ax.axis("off")
            continue
        by_pos = defaultdict(lambda: defaultdict(int))
        for (pos, nc), n in gene_counts.items():
            by_pos[pos][nc] += n
        positions = sorted(by_pos)
        bottoms = np.zeros(len(positions))
        for nc in ["synonymous", "beyond_cds", "undetermined", "protein_altering"]:
            heights = np.array([by_pos[p].get(nc, 0) for p in positions], dtype=float)
            if heights.sum() == 0:
                continue
            ax.vlines(positions, bottoms, bottoms + heights, color=NOVELTY_COLORS[nc], lw=1.2,
                      alpha=0.85, label=nc)
            bottoms += heights
        totals_by_pos = [(sum(by_pos[p].values()), p) for p in positions]
        top_pos = sorted(totals_by_pos, reverse=True)[:3]
        for total, pos in top_pos:
            ax.annotate(f"codon {pos}", xy=(pos, sum(by_pos[pos].values())), fontsize=6.5,
                        ha="center", va="bottom", rotation=90)
        ax.set_title(f"{vc.gene_display(gene)} (n_positions={len(positions)})", fontsize=10.5)
        ax.set_xlabel("Codon position (CDS)", fontsize=8)
        ax.set_ylabel("N haplotype observations", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = [], []
    for nc, color in NOVELTY_COLORS.items():
        handles.append(plt.Line2D([0], [0], color=color, lw=3))
        labels.append(nc)
    fig.legend(handles, labels, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.04),
               frameon=False)
    fig.suptitle(title, y=1.08, fontsize=13)
    fig.tight_layout()
    vc.savefig(fig, out_path, dpi=140)


def write_report(path, counts, genes):
    lines = ["# Mutation topology: where novel differences land along the CDS "
             "(`09_mutation_topology.py`)\n",
             "Methodology: scripts/hla_popgen/research/VIZ_LIT.md addendum (2026-09). Each stem = "
             "one codon position; height = number of haplotype-level novel-difference observations "
             "at that position (recurrence-weighted, not a per-cluster count -- see module "
             "docstring). Position comes from `cds_mut`'s aa-diff field, first tied candidate "
             "only.\n"]
    lines.append("\n## Top-5 hottest codons per gene\n")
    lines.append("| gene | codon | n observations | dominant novelty_class |")
    lines.append("|---|---|---|---|")
    for gene in genes:
        gene_counts = counts.get(gene)
        if not gene_counts:
            continue
        by_pos = defaultdict(Counter)
        for (pos, nc), n in gene_counts.items():
            by_pos[pos][nc] += n
        totals = {pos: sum(c.values()) for pos, c in by_pos.items()}
        for pos, total in sorted(totals.items(), key=lambda kv: -kv[1])[:5]:
            dominant = by_pos[pos].most_common(1)[0][0]
            lines.append(f"| {vc.gene_display(gene)} | {pos} | {total} | {dominant} |")
    lines.append(f"\nFigure: `{os.path.join(os.path.dirname(path), os.path.basename(path).replace('_report.md', '.png'))}`\n")
    vc.write_report(path, lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outroot", default=vc.DEFAULT_OUTROOT,
                    help="Directory with hla_calls_rich.tsv (real run: ~/pipeline_outputs).")
    ap.add_argument("--table1", default=None, help="Override path to hla_calls_rich.tsv.")
    ap.add_argument("--out-dir", default=None,
                    help="Where to write the figure + report. Default: "
                         "reports/hla_popgen/09_mutation_topology/lr/")
    ap.add_argument("--sample", action="store_true")
    args = ap.parse_args()

    table1_path = args.table1 or os.path.join(args.outroot, "hla_calls_rich.tsv")
    out_dir = args.out_dir or vc.default_out_dir("09_mutation_topology", "lr")
    suffix = ".sample" if args.sample else ""

    print(f"Loading Table 1 from {table1_path!r} ...", file=sys.stderr)
    table1 = vc.load_table1(table1_path)
    if "novelty_class" not in table1.columns:
        sys.exit("FATAL: Table 1 has no 'novelty_class' column -- run 01_extract_rich.py first.")

    print("Parsing cds_mut codon positions for novel classical-gene calls ...", file=sys.stderr)
    counts = build_position_counts(table1, GENES)
    n_positions = sum(len(c) for c in counts.values())
    print(f"  {n_positions} distinct (gene, codon, novelty_class) cells", file=sys.stderr)

    fig_path = os.path.join(out_dir, f"mutation_topology_classical{suffix}.png")
    plot_lollipops(counts, GENES, fig_path,
                    "Novel-difference codon positions, classical genes (lollipop, "
                    "recurrence-weighted)")
    md_path = os.path.join(out_dir, f"09_mutation_topology_report{suffix}.md")
    write_report(md_path, counts, GENES)
    print(f"Wrote figure to {fig_path!r} and report to {md_path!r}.", file=sys.stderr)


if __name__ == "__main__":
    main()
