#!/usr/bin/env python3
"""Where, along the CDS, do novel differences actually land? A gene-faceted lollipop/needle plot
(cBioPortal MutationMapper / ProteinPaint / Bioconductor trackViewer idiom, per the 2026-09
visualization-literature review in research/VIZ_LIT.md's addendum): one vertical stem per distinct
codon position, height = number of haplotype-level observations of a novel difference at that codon
(a direct recurrence-weighted count, not a per-cluster count -- a position hit by many people's
independently-assembled haplotypes is a stronger "real hotspot" signal than a raw candidate count),
colored by the call's `novelty_class`.

## Where the position comes from

`cds_mut` (Table 1, reference/IMMUANNOT_GTF_SPEC.md part A) is
`[best-matching reference allele]|[minimap2 cs string]|[aa diff]`, comma-joined for tied candidates.
Position comes from the middle field: the minimap2 `cs` short-form string, whose `:N` tokens are
runs of N identical bases -- summing them gives the running nucleotide offset at each substitution
(`*xy`) or indel (`+seq`/`-seq`) event, converted to a 1-based codon position as
`nt_offset // 3 + 1`, deduped to one entry per distinct codon (two nt-level events in the same
codon must not be double-counted as two positions).

## Why coloring is the row's `novelty_class`, not a per-codon synonymous/non-synonymous call
(REVISED 2026-09, after a real-data regression -- read this before "fixing" it again)

An earlier version of this script tried to derive a TRUE per-codon synonymous call by parsing the
aa-diff field's `REFaa(REFcodon)<OBSaa(OBScodon)` entries directly and zipping them, in order,
against the cs-string's deduped codon positions. That was built and validated against
`tests/make_fixtures.py`'s synthetic `cds_mut` strings, which use clean single-letter amino acid
codes (K, D, E, ...) always in `X(codon)<Y(codon)` pairs. **Real Immuannot output does not look
like that**, confirmed against the actual VM data this session:
  `HLA-W*05:02|:67*tc:41*ga:7*ac:149-a:365-g:...|Tre(ACG)<Met(ATG):Gln(CAG)<Rrg(CGG):
   Rrg(CGC)<Ser(AGC):ProLys(CCCAAA)<CCCAAAA(CCCAAA):Val(GTG)<GTGG(GTGG):Ala(GCG):...`
Three real-data properties the fixture never modeled: (1) amino acid labels are multi-letter,
apparently-abbreviated/occasionally-garbled tokens (`Tre`, `Rrg`), not single IUPAC letters: no
fixed-width regex is safe; (2) indel-affected entries can span MULTIPLE consolidated residues with
variable-length codon strings (`CCCAAAA`, or much longer runs for larger indels) rather than one
fixed 3-letter codon per entry, so the entry COUNT does not reliably equal the cs-string's
substitution/indel TOKEN count; (3) some entries carry no `<OBS(...)` half at all (`Ala(GCG)` alone)
-- a single-sided annotation, not a diff pair. Attempting the position<->aa-diff zip against this
real format gave a 100% `n_position_aa_mismatch` rate on the actual cohort (0 of ~5,000 calls
usable) -- silently worse than the four-way row-level label it was meant to improve on. Reverted.

**What "I don't see synonymous events" actually reflects**: per 03_novel_alleles.py's real report,
resolved novel clusters ARE genuinely dominated by `protein_altering` (90.8%), with `synonymous`
only 3.2% and `beyond_cds` 6.1% -- a thin green sliver in a stacked lollipop plot is the correct
picture of the real class balance, not evidence of a mislabeling bug. A trustworthy true PER-CODON
synonymous call would need either Immuannot's own `algntools.findCodingDiff()` source (to learn its
exact abbreviation table and entry-consolidation rule for indels) or a from-scratch codon-by-codon
re-translation of the observed vs. reference CDS sequences (`cds.fa.gz` vs the matched reference's
own CDS) -- both bigger asks than this pass, flagged here rather than guessed at again.

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
# minimap2 short-form `cs` tag tokens: ':N' (N identical bases), '*xy' (substitution), '+seq'/
# '-seq' (indel). Standard, well-defined format -- unlike the aa-diff field (see module docstring),
# this part of cds_mut has held up against real data without surprises.
CS_TOKEN_RE = re.compile(r":\d+|\*[a-z]{2}|[+-][a-z]+", re.IGNORECASE)


def parse_codon_positions(cds_mut):
    """First (best) tied candidate's cs-string -> deduped, ordered list of 1-based codon positions
    (nt_offset // 3 + 1) for every substitution/indel event. Returns [] for missing/unparseable
    input -- fails toward "no position data" rather than a guessed one."""
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
            pos = nt_pos // 3 + 1
            if not codons or codons[-1] != pos:
                codons.append(pos)
            nt_pos += 1
        elif tok[0] in "+-":
            pos = nt_pos // 3 + 1
            if not codons or codons[-1] != pos:
                codons.append(pos)
            if tok[0] == "-":
                nt_pos += len(tok) - 1  # deletion consumes reference bases; insertion doesn't
    return codons


def build_position_counts(table1, genes):
    """Returns {gene: Counter({(codon_pos, novelty_class): n_observations})} -- novelty_class is
    the ROW's own whole-call label (see module docstring for why not per-codon)."""
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
        # Draw the rare classes first (bottom of stack) so protein_altering, the large majority
        # class, doesn't visually bury them.
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
        ax.set_ylabel("N haplotypes with a diff here\n(recurrence, not distinct AA types)",
                       fontsize=7.5)
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
             "one codon position; height = number of haplotype-level novel-difference OBSERVATIONS "
             "at that position -- a recurrence count across the cohort, NOT the number of distinct "
             "amino-acid substitution types (capped near 19-20 by definition; recurrence is not, "
             "since many different haplotypes can independently carry the same substitution). Color "
             "is the call's own `novelty_class` (protein_altering/synonymous/beyond_cds/"
             "undetermined) -- a true per-codon synonymous call was attempted and reverted after a "
             "real-data regression; see module docstring for exactly why (short version: real "
             "`cds_mut` amino-acid labels are multi-letter and occasionally-garbled, entries can be "
             "indel-consolidated across several residues, and some entries have no diff pair at "
             "all -- none of which the synthetic test fixture modeled).\n"]
    lines.append("\n## Why the plot is mostly protein_altering (red)\n"
                 "This is the real class balance, not a labeling bug: of resolved novel clusters "
                 "cohort-wide (03_novel_alleles.py's report), 90.8% are `protein_altering`, 3.2% "
                 "`synonymous`, 6.1% `beyond_cds`. A thin green sliver here is the correct picture, "
                 "not evidence something is hidden.\n")
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
