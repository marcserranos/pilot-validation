#!/usr/bin/env python3
"""Where, along the CDS, do novel differences actually land? A gene-faceted lollipop/needle plot
(cBioPortal MutationMapper / ProteinPaint / Bioconductor trackViewer idiom, per the 2026-09
visualization-literature review in research/VIZ_LIT.md's addendum): one vertical stem per distinct
codon position, height = number of haplotype-level observations of a novel difference at that codon
(a direct recurrence-weighted count, not a per-cluster count -- a position hit by many people's
independently-assembled haplotypes is a stronger "real hotspot" signal than a raw candidate count),
colored by `novelty_class`.

## Where the position AND the per-codon synonymous/non-synonymous call come from

`cds_mut` (Table 1, reference/IMMUANNOT_GTF_SPEC.md part A) is
`[best-matching reference allele]|[minimap2 cs string]|[aa diff]`, comma-joined for tied candidates.
Position comes from the middle field: the minimap2 `cs` short-form string, whose `:N` tokens are
runs of N identical bases -- summing them gives the running nucleotide offset at each substitution
(`*xy`) or indel (`+seq`/`-seq`) event, converted to a 1-based codon position as
`nt_offset // 3 + 1`, then DEDUPED to one entry per distinct codon (two nt-level events landing in
the same codon must not be double-counted as two positions).

`novelty_class` is deliberately NOT taken from the Table 1 row's own scalar `novelty_class` column
here -- that field describes the DEEPEST novelty across the WHOLE call (SCHEMA.md Table 1), so a
`protein_altering` row can still contain individually-silent codon diffs when it has multiple
diffs, and labelling every one of them "protein_altering" is exactly why an earlier version of this
plot showed almost no green (synonymous) even though real synonymous events are present -- see
NOVEL_LIT.md's discovery-certainty addendum for the incident writeup. The fix: the aa-diff field
(`REFaa(REFcodon)<OBSaa(OBScodon)`, colon-joined, one entry per distinct changed codon) already
carries the amino acid on both sides of that ONE codon -- comparing REFaa==OBSaa per entry gives a
TRUE per-codon synonymous/non-synonymous call, independent of the row's whole-call label. The
deduped cs-string codon positions and the aa-diff entries are zipped together IN ORDER (both are
written in left-to-right CDS order by the same upstream code); a count mismatch between the two
lists is a real parse failure, not something to guess through -- those calls are dropped from the
per-codon breakdown and counted in `n_position_aa_mismatch` in the report rather than silently
mislabeled.

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
# Two categories now, both TRUE per-codon calls (see parse_aa_diff_entries) -- not the four-way
# novelty_class, which describes the whole call rather than one codon within it.
SYN_COLORS = {False: "#C44E52", True: "#55A868"}  # False=non-synonymous, True=synonymous
SYN_LABELS = {False: "non-synonymous", True: "synonymous"}
# KNOWN FIXTURE GAP (does not affect real data): tests/make_fixtures.py's cds_mut generator always
# emits exactly ONE cs-string '*xy' token regardless of how many aa-diff entries it writes (1-4,
# via `cds_dist`) -- the two counts were never made consistent with each other in the fixture, so
# n_mismatch is near-100% against fixtures by construction. Real minimap2 cs output has no such
# gap: a '*xy' token is always exactly one base substitution, so a genuine N-codon diff call from
# Immuannot has N real substitution/indel tokens, not one. Verified the position/aa-diff zip logic
# directly against a hand-built matching example instead (same left-to-right semantics, just with
# a fixture that doesn't fake a single token for multiple diffs) -- see git history for the check.
# `cds_mut`'s aa-diff field is "REFaa(REFcodon)<OBSaa(OBScodon)" where (REFcodon)/(OBScodon) are the
# CODON TRIPLET SEQUENCE (e.g. "GGC"), not a position -- there is no position number in that field
# at all. The actual position lives in the MIDDLE pipe field, the minimap2 `cs` short-form string
# (reference/IMMUANNOT_GTF_SPEC.md part A/C; same field 03_novel_alleles.py's
# `is_homopolymer_indel_only()` parses with this identical token regex): ':N' tokens are runs of N
# identical bases, so summing them gives the running nucleotide OFFSET into the CDS at each event.
CS_TOKEN_RE = re.compile(r":\d+|\*[a-z]{2}|[+-][a-z]+", re.IGNORECASE)
AA_DIFF_ENTRY_RE = re.compile(r"^([A-Za-z*])\([ACGTacgt]+\)<([A-Za-z*])\([ACGTacgt]+\)$")


def parse_codon_positions(cds_mut):
    """First (best) tied candidate's cs-string -> deduped, ordered list of 1-based codon positions
    (nt_offset // 3 + 1) for every substitution/indel event -- deduped because two nt-level events
    in the same codon must produce one codon entry, matching how aa-diff itself is one-entry-per-
    changed-codon. Returns [] for missing/unparseable input -- fails toward "no position data"
    rather than a guessed one."""
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


def parse_aa_diff_entries(cds_mut):
    """First (best) tied candidate's aa-diff field -> ordered list of (ref_aa, obs_aa) tuples, one
    per colon-joined codon diff, same left-to-right order as parse_codon_positions()'s cs-string
    codons. Returns [] if the aa-diff field is absent/unparseable (e.g. an indel-only diff with no
    clean codon translation)."""
    if not isinstance(cds_mut, str) or not cds_mut.strip():
        return []
    first_candidate = cds_mut.split(",")[0]
    parts = first_candidate.split("|")
    if len(parts) < 3 or not parts[2]:
        return []
    out = []
    for entry in parts[2].split(":"):
        m = AA_DIFF_ENTRY_RE.match(entry.strip())
        if m is None:
            return []  # one unparseable entry invalidates the whole ordered match -- don't guess
        out.append((m.group(1), m.group(2)))
    return out


def build_position_counts(table1, genes):
    """Returns (counts, distinct_subs, n_mismatch, n_used).

    counts: {gene: Counter({(codon_pos, is_synonymous): n_observations})} -- is_synonymous is a
    TRUE per-codon call (see module docstring), not the row's whole-call novelty_class.
    distinct_subs: {gene: {codon_pos: set of (ref_aa, obs_aa) actually observed}} -- the answer to
    "how many DIFFERENT substitutions, not how many times", kept separate from the recurrence count
    since the two answer different questions (see write_report's table).
    n_mismatch: calls where the cs-string codon count and aa-diff entry count didn't line up --
    excluded from both structures above rather than mislabeled."""
    sub = table1[table1["gene_bare"].isin(genes) & table1["is_novel"]]
    counts = defaultdict(Counter)
    distinct_subs = defaultdict(lambda: defaultdict(set))
    n_mismatch = 0
    n_used = 0
    for gene, cds_mut in zip(sub["gene_bare"], sub["cds_mut"]):
        positions = parse_codon_positions(cds_mut)
        if not positions:
            continue
        aa_entries = parse_aa_diff_entries(cds_mut)
        if len(aa_entries) != len(positions):
            n_mismatch += 1
            continue
        n_used += 1
        for pos, (ref_aa, obs_aa) in zip(positions, aa_entries):
            counts[gene][(pos, ref_aa == obs_aa)] += 1
            distinct_subs[gene][pos].add((ref_aa, obs_aa))
    return counts, distinct_subs, n_mismatch, n_used


def plot_lollipops(counts, genes, out_path, title):
    fig, axes = plt.subplots(2, 4, figsize=(22, 9))
    for ax, gene in zip(axes.flat, genes):
        gene_counts = counts.get(gene)
        if not gene_counts:
            ax.axis("off")
            continue
        by_pos = defaultdict(lambda: defaultdict(int))
        for (pos, is_syn), n in gene_counts.items():
            by_pos[pos][is_syn] += n
        positions = sorted(by_pos)
        bottoms = np.zeros(len(positions))
        # Draw synonymous first (bottom of stack) so non-synonymous, the majority class, doesn't
        # visually bury it -- both are now real per-codon calls, not a whole-call label smeared
        # across every diff in a multi-diff call.
        for is_syn in [True, False]:
            heights = np.array([by_pos[p].get(is_syn, 0) for p in positions], dtype=float)
            if heights.sum() == 0:
                continue
            ax.vlines(positions, bottoms, bottoms + heights, color=SYN_COLORS[is_syn], lw=1.2,
                      alpha=0.85, label=SYN_LABELS[is_syn])
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
    for is_syn, color in SYN_COLORS.items():
        handles.append(plt.Line2D([0], [0], color=color, lw=3))
        labels.append(SYN_LABELS[is_syn])
    fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.04),
               frameon=False)
    fig.suptitle(title, y=1.08, fontsize=13)
    fig.tight_layout()
    vc.savefig(fig, out_path, dpi=140)


def write_report(path, counts, distinct_subs, n_mismatch, n_used, genes):
    lines = ["# Mutation topology: where novel differences land along the CDS "
             "(`09_mutation_topology.py`)\n",
             "Methodology: scripts/hla_popgen/research/VIZ_LIT.md addendum (2026-09). Each stem = "
             "one codon position; height = number of haplotype-level novel-difference OBSERVATIONS "
             "at that position -- a recurrence count across the cohort, NOT the number of distinct "
             "amino-acid substitution types (which is capped near 19 by definition; recurrence is "
             "not, since many different haplotypes can independently carry the same substitution). "
             "`n_distinct_subs` below is that other number, kept separate. `synonymous`/"
             "`non-synonymous` are now TRUE PER-CODON calls (REFaa vs OBSaa compared directly from "
             "the aa-diff field for that one codon), not the row's whole-call `novelty_class` -- a "
             "`protein_altering` call with several diffs can still contain individually-silent "
             "codons, and this script now shows them as such instead of painting the whole call one "
             "color.\n"]
    lines.append(f"\n{n_used} novel calls contributed a clean, position-matched breakdown; "
                 f"{n_mismatch} were excluded because the cs-string codon count and aa-diff entry "
                 f"count didn't line up (a real parse ambiguity, not guessed through).\n")
    lines.append("\n## Top-5 hottest codons per gene\n")
    lines.append("| gene | codon | n observations (recurrence) | n_distinct_subs | "
                 "dominant call |")
    lines.append("|---|---|---|---|---|")
    for gene in genes:
        gene_counts = counts.get(gene)
        if not gene_counts:
            continue
        by_pos = defaultdict(Counter)
        for (pos, is_syn), n in gene_counts.items():
            by_pos[pos][is_syn] += n
        totals = {pos: sum(c.values()) for pos, c in by_pos.items()}
        for pos, total in sorted(totals.items(), key=lambda kv: -kv[1])[:5]:
            dominant_syn = by_pos[pos].most_common(1)[0][0]
            dominant = "synonymous" if dominant_syn else "non-synonymous"
            n_distinct = len(distinct_subs.get(gene, {}).get(pos, ()))
            lines.append(f"| {vc.gene_display(gene)} | {pos} | {total} | {n_distinct} | "
                         f"{dominant} |")
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
    counts, distinct_subs, n_mismatch, n_used = build_position_counts(table1, GENES)
    n_positions = sum(len(c) for c in counts.values())
    print(f"  {n_positions} distinct (gene, codon, is_synonymous) cells; {n_used} calls used, "
          f"{n_mismatch} excluded on a position/aa-diff count mismatch", file=sys.stderr)

    fig_path = os.path.join(out_dir, f"mutation_topology_classical{suffix}.png")
    plot_lollipops(counts, GENES, fig_path,
                    "Novel-difference codon positions, classical genes (lollipop, "
                    "recurrence-weighted)")
    md_path = os.path.join(out_dir, f"09_mutation_topology_report{suffix}.md")
    write_report(md_path, counts, distinct_subs, n_mismatch, n_used, GENES)
    print(f"Wrote figure to {fig_path!r} and report to {md_path!r}.", file=sys.stderr)


if __name__ == "__main__":
    main()
