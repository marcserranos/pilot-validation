#!/usr/bin/env python3
"""Does enriching the allele space (known + novel, uncollapsed) change what the population
structure looks like -- and does collapsing it back down (2-field truncation, or SR's native
resolution) erase real signal? Four embeddings, same people where possible, same plotting/QC
conventions as `06_figures_structure.py` (manual sklearn-free PCA + optional UMAP, Okabe-Ito
ancestry palette, MIN_CELL_N_PEOPLE discipline).

## The four views (see README.md / chat discussion 2026-09 for the full request)

1. **`lr_full`** -- LR cohort, 8 classical genes, complete diploid cases, each allele call encoded
   at its OWN native resolution: the full-precision `consensus` string for known (non-novel) calls,
   `novel_id` (SCHEMA.md Table 3's `<gene>_nov_<sha1:8>`) for novel calls. This is the big-N,
   allele-enriched space and the PRIMARY fitted embedding -- both PCA and UMAP are fit ONCE here.
2. **`lr_known_only`** -- the exact SAME fitted `lr_full` embedding, subset down to just the people
   who carry ZERO novel calls across all 8 genes. Not a re-fit: literally the same coordinates,
   filtered. This answers "where do already-fully-catalogued people sit, and do novel-carriers
   occupy visibly different territory in the SAME coordinate frame" -- which is the only
   mathematically honest way to compare a subset against an enriched embedding (see module
   docstring section "Why not project a collapsed encoding through the same PCA loadings" below for
   why the fine-vs-collapsed comparison can NOT be done this way instead).
3. **`lr_collapsed`** -- a SEPARATELY FIT embedding: 2-field-truncated `consensus` (novel calls
   truncate to their known common prefix, exactly `cluster_hla_by_ancestry.py`/06's convention),
   same people as `lr_full` where a complete 2-field case exists.
4. **`sr_collapsed`** -- a separately fit embedding on the AoU-native short-read cohort, 2-field
   (SR's native resolution), restricted to people who ALSO have a complete `lr_full` case, so 3 and
   4 are the same individuals viewed through two different assays at the same nominal resolution.

## Why not project a collapsed encoding through the same PCA loadings

The natural-sounding ask ("fit once on the big space, then show collapsing in that same mapping")
is only literally possible when the coarser view's columns are a SUBSET of the fine view's columns
-- true for `lr_known_only` (dropping novel columns is a projection: zero out those columns, keep
the rest, multiply by the same loadings) but NOT true for `lr_collapsed`/`sr_collapsed`: a 2-field
bucket like `HLA-A*02` is a MERGER of dozens of distinct fine-grained `lr_full` columns
(`HLA-A*02:01:01:01`, `HLA-A*02:01:01:02`, `HLA-A_nov_...`, ...), and merging is not invertible --
there is no principled way to place `HLA-A*02` at a single point in the fine feature space without
guessing which fine allele each collapsed call "really" was. So `lr_collapsed`/`sr_collapsed` are
each their own PCA/UMAP fit (exactly 06's existing method), and the comparison is made by holding
the PEOPLE, the color coding, and the panel layout fixed across all four panels, not the coordinate
system -- structure that survives collapsing should still separate ancestry/novel-carrier groups in
its own fit; structure that was only visible at fine resolution should visibly flatten out.

## Encoding

Same dosage convention as 06_figures_structure.py throughout: one column per (gene, allele-identity)
pair, value = copy count (0/1/2) per person, unordered hap1/hap2, complete-case (both copies
resolved for all 8 classical genes) only -- incomplete cases are dropped and counted, never imputed
(SCHEMA.md discipline).

Usage (fixtures): python3 scripts/hla_popgen/tests/make_fixtures.py --outroot /tmp/hla_fixtures -n 300
    then: python3 scripts/hla_popgen/08_embedding_compare.py --outroot /tmp/hla_fixtures \\
        --table1 /tmp/hla_fixtures/hla_calls_rich.sample.tsv \\
        --cohort-membership /tmp/hla_fixtures/cohort_membership.sample.tsv \\
        --sr-genotypes /tmp/hla_fixtures/hla_genotypes.tsv --sample

Real run (VM): python3 scripts/hla_popgen/08_embedding_compare.py
"""
import argparse
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _viz_common as vc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

GENES = vc.CLASSICAL_GENES_BARE
NOVEL_COLORS = {False: "#4C72B0", True: "#C44E52"}
# Matches 03_novel_alleles.py's convention exactly: person_id dirs (and cds.fa.gz) live under
# people/, one level below where hla_calls_rich.tsv/cohort_membership.tsv actually sit.
DEFAULT_DATA_ROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_OUTROOT = os.path.join(DEFAULT_DATA_ROOT, "people")


# ---------------------------------------------------------------------------
# Identity assignment (mirrors 03_novel_alleles.py / 04_allele_saturation.py's matching logic)
# ---------------------------------------------------------------------------
def load_novel_matches(table1, outroot):
    """Re-derives the (person_id, hap, gene) -> novel_id map exactly as 03/04 do, by importing
    03_novel_alleles.py as a sibling module (never re-implement the matching logic twice)."""
    import importlib.util
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location(
        "hla_popgen_novel_alleles", os.path.join(here, "03_novel_alleles.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    matched_rows, _seqs, stats = mod.match_novel_rows(table1, outroot)
    lookup = {}
    for r in matched_rows:
        novel_id = f"{r['gene']}_nov_{r['cds_seq_sha1'][:8]}"
        lookup[(r["person_id"], r["hap"], r["gene"])] = novel_id
    return lookup, stats


def build_dosage_matrix(table1, person_ids, genes, resolution, novel_id_lookup=None):
    """resolution: 'full' (native precision: known consensus string, or novel_id for novel calls)
    or 'collapsed' (2-field truncation of consensus; novel calls fall back to whatever 2-field
    prefix survives truncation, same as 06's existing convention -- a novel call and its nearest
    documented allele share a 2-field bucket by construction of Immuannot's own naming rule,
    reference/IMMUANNOT_GTF_SPEC.md part B). Returns (dosage_df, novel_carrier: {person_id: bool},
    n_incomplete)."""
    sub = table1[table1["person_id"].isin(person_ids) & table1["gene_bare"].isin(genes)]
    per_person = {}
    novel_carrier = {}
    for pid, gene, hap, consensus, is_novel in zip(
            sub["person_id"], sub["gene_bare"], sub["hap"], sub["consensus"], sub["is_novel"]):
        if is_novel:
            novel_carrier[pid] = True  # tracked regardless of resolution -- "did this person carry
            # >=1 novel call in these genes" is a fact about the person, independent of whether the
            # identity space we're currently building can represent that call at fine resolution.
        else:
            novel_carrier.setdefault(pid, False)

        if resolution == "full" and is_novel:
            ident = (novel_id_lookup or {}).get((pid, hap, vc.gene_display(gene)))
            if ident is None:
                continue  # unresolved/ambiguous novel match -- never guess (SCHEMA.md)
        else:
            ident = vc.to_nfield(consensus, 2) if resolution == "collapsed" else (
                None if (isinstance(consensus, str) and consensus == "undetermined") else consensus)
            if ident is None:
                continue
        slots = per_person.setdefault(pid, {}).setdefault(gene, {})
        slots[len(slots) + 1] = ident

    complete, incomplete = [], 0
    rows = {}
    for pid, gene_map in per_person.items():
        if len(gene_map) < len(genes):
            incomplete += 1
            continue
        vals, ok = [], True
        for gene in genes:
            copies = gene_map.get(gene, {})
            v = [copies.get(1), copies.get(2)]
            if any(x is None for x in v) or len(copies) < 2:
                ok = False
                break
            vals.append((gene, v[0]))
            vals.append((gene, v[1]))
        if not ok:
            incomplete += 1
            continue
        rows[pid] = Counter(f"{g}:{a}" for g, a in vals)
        complete.append(pid)
        novel_carrier.setdefault(pid, False)

    feature_names = sorted({feat for c in rows.values() for feat in c})
    mat = pd.DataFrame(0.0, index=complete, columns=feature_names)
    for pid, counter in rows.items():
        for feat, n in counter.items():
            mat.at[pid, feat] = float(n)
    return mat, {pid: novel_carrier.get(pid, False) for pid in complete}, incomplete


def build_sr_dosage_matrix(sr_long, person_ids, genes):
    sub = sr_long[sr_long["person_id"].isin(person_ids) & sr_long["gene_bare"].isin(genes)]
    per_person = {}
    for pid, gene, copy, allele in zip(sub["person_id"], sub["gene_bare"], sub["copy"],
                                        sub["allele"]):
        a2 = vc.to_nfield(allele, 2)
        if a2 is None:
            continue
        per_person.setdefault(pid, {}).setdefault(gene, {})[copy] = a2
    complete, incomplete = [], 0
    rows = {}
    for pid, gene_map in per_person.items():
        if len(gene_map) < len(genes):
            incomplete += 1
            continue
        vals, ok = [], True
        for gene in genes:
            copies = gene_map.get(gene, {})
            v = [copies.get(1), copies.get(2)]
            if any(x is None for x in v):
                ok = False
                break
            vals.append((gene, v[0]))
            vals.append((gene, v[1]))
        if not ok:
            incomplete += 1
            continue
        rows[pid] = Counter(f"{g}:{a}" for g, a in vals)
        complete.append(pid)
    feature_names = sorted({feat for c in rows.values() for feat in c})
    mat = pd.DataFrame(0.0, index=complete, columns=feature_names)
    for pid, counter in rows.items():
        for feat, n in counter.items():
            mat.at[pid, feat] = float(n)
    return mat, incomplete


# ---------------------------------------------------------------------------
# Fitting + plotting
# ---------------------------------------------------------------------------
def fit_embedding(mat, seed, umap_neighbors, umap_min_dist):
    X = vc.standardize(mat.values)
    scores, var_ratio = vc.pca_fit_transform(X, n_components=10)
    umap_emb = None
    if vc.HAVE_UMAP:
        import umap
        umap_emb = umap.UMAP(n_neighbors=umap_neighbors, min_dist=umap_min_dist,
                              random_state=seed).fit_transform(X)
    return scores, var_ratio, umap_emb


def plot_two_color(coords, index, ancestry_by_person, novel_carrier, xlabel, ylabel, title,
                    out_path_prefix, var_note=""):
    """One embedding, two coloring panels side by side: ancestry, and novel-carrier status."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 6))
    anc_labels = [ancestry_by_person.get(pid) for pid in index]
    for anc in vc.ANCESTRY_ORDER:
        mask = np.array([l == anc for l in anc_labels])
        if not mask.any():
            continue
        ax1.scatter(coords[mask, 0], coords[mask, 1], s=10, alpha=0.55,
                    color=vc.ANCESTRY_COLORS[anc], label=f"{anc} (n={mask.sum()})",
                    edgecolors="none")
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title("Colored by ancestry", fontsize=10)
    ax1.legend(fontsize=7.5, frameon=False, markerscale=1.5)
    ax1.spines[["top", "right"]].set_visible(False)

    carrier = np.array([bool(novel_carrier.get(pid, False)) for pid in index])
    for flag, label in [(False, "known-only (this gene set)"), (True, "carries >=1 novel allele")]:
        mask = carrier == flag
        if not mask.any():
            continue
        ax2.scatter(coords[mask, 0], coords[mask, 1], s=10, alpha=0.55,
                    color=NOVEL_COLORS[flag], label=f"{label} (n={mask.sum()})", edgecolors="none")
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel)
    ax2.set_title("Colored by novel-carrier status", fontsize=10)
    ax2.legend(fontsize=7.5, frameon=False, markerscale=1.5)
    ax2.spines[["top", "right"]].set_visible(False)

    fig.suptitle(title + (f"\n{var_note}" if var_note else ""), fontsize=11)
    fig.tight_layout()
    vc.savefig(fig, out_path_prefix)


# ---------------------------------------------------------------------------
# Collapse-radius diagnostics (LR-collapsed vs SR-collapsed, same people, same nominal resolution)
# ---------------------------------------------------------------------------
def concentration_stats(mat, genes):
    """Per gene: distinct alleles observed, and the share of all copies covered by the top-10 most
    common alleles (a direct 'how concentrated on a few common alleles' statistic -- the quantity
    the SR-vs-LR 'collapse radius' hypothesis is actually about, per 2026-09 chat discussion)."""
    rows = []
    for gene in genes:
        cols = [c for c in mat.columns if c.startswith(f"{gene}:")]
        if not cols:
            continue
        totals = mat[cols].sum(axis=0).sort_values(ascending=False)
        total_copies = totals.sum()
        if total_copies == 0:
            continue
        top10_share = totals.head(10).sum() / total_copies
        freqs = (totals / total_copies).values
        simpson = float(np.sum(freqs ** 2))
        rows.append({"gene": gene, "n_distinct_alleles": int((totals > 0).sum()),
                     "total_copies": int(total_copies), "top10_share": round(top10_share, 4),
                     "simpson_concentration": round(simpson, 4)})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT,
                    help="Root holding <person_id>/immuannot_output/hap{1,2}/cds.fa.gz (needed to "
                         "re-derive novel_id per haplotype). Default: ~/pipeline_outputs/people "
                         "(NOT ~/pipeline_outputs itself, see RUNBOOK.md).")
    ap.add_argument("--table1", default=None, help="Path to hla_calls_rich.tsv. Default: "
                                                     "~/pipeline_outputs/hla_calls_rich.tsv")
    ap.add_argument("--cohort-membership", default=None,
                    help="Path to cohort_membership.tsv (Table 4). Default: "
                         "~/pipeline_outputs/cohort_membership.tsv")
    ap.add_argument("--sr-genotypes", default=None,
                    help="Path to hla_genotypes.tsv (AoU-native SR calls). Default: "
                         "~/pipeline_outputs/hla_genotypes.tsv")
    ap.add_argument("--out-dir", default=None,
                    help="Where to write PNGs + the markdown report. Default: "
                         "reports/hla_popgen/08_embedding_compare/lr/")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--umap-neighbors", type=int, default=15)
    ap.add_argument("--umap-min-dist", type=float, default=0.1)
    ap.add_argument("--sample", action="store_true")
    args = ap.parse_args()

    table1_path = args.table1 or os.path.join(DEFAULT_DATA_ROOT, "hla_calls_rich.tsv")
    cohort_path = args.cohort_membership or os.path.join(DEFAULT_DATA_ROOT, "cohort_membership.tsv")
    out_dir = args.out_dir or vc.default_out_dir("08_embedding_compare", "lr")
    suffix = ".sample" if args.sample else ""

    print(f"Loading Table 1 from {table1_path!r} ...", file=sys.stderr)
    table1 = vc.load_table1(table1_path)
    print(f"Loading Table 4 (ancestry) from {cohort_path!r} ...", file=sys.stderr)
    cohort_df = vc.load_cohort_membership(cohort_path)
    ancestry_by_person = dict(zip(cohort_df["person_id"], cohort_df["ancestry_pred"]))
    lr_people = set(cohort_df.loc[cohort_df["in_lr"], "person_id"]) if "in_lr" in \
        cohort_df.columns else set(table1["person_id"].unique())
    sr_people = set(cohort_df.loc[cohort_df["in_sr"], "person_id"]) if "in_sr" in \
        cohort_df.columns else set()

    print("Re-deriving novel-allele haplotype matches (03_novel_alleles.py's logic) ...",
          file=sys.stderr)
    novel_id_lookup, match_stats = load_novel_matches(table1, args.outroot)
    print(f"  {match_stats}", file=sys.stderr)

    report = ["# Embedding comparison: allele-space resolution vs population structure "
              "(`08_embedding_compare.py`)\n",
              "See this script's module docstring for why `lr_collapsed`/`sr_collapsed` are "
              "separately-fit embeddings (not projections of `lr_full`) while `lr_known_only` IS "
              "the same fitted `lr_full` coordinates, subset.\n"]

    # --- 1 & 2: lr_full (fit once) + lr_known_only (same coords, subset) ---------------------
    print("Building lr_full (known full-resolution + novel_id) dosage matrix ...", file=sys.stderr)
    mat_full, novel_carrier, n_incomplete_full = build_dosage_matrix(
        table1, lr_people, GENES, "full", novel_id_lookup)
    report.append(f"\n## lr_full\nComplete 8-classical-gene cases: {len(mat_full)} "
                   f"({n_incomplete_full} incomplete, dropped). Novel carriers: "
                   f"{sum(novel_carrier.values())} ({100*sum(novel_carrier.values())/max(len(mat_full),1):.1f}%).\n")
    scores_full, var_full, umap_full = fit_embedding(mat_full, args.seed, args.umap_neighbors,
                                                      args.umap_min_dist)
    pca_full_path = os.path.join(out_dir, f"pca_lr_full{suffix}.png")
    plot_two_color(scores_full[:, :2], mat_full.index, ancestry_by_person, novel_carrier,
                    f"PC1 ({100*var_full[0]:.1f}% var)", f"PC2 ({100*var_full[1]:.1f}% var)",
                    "PCA -- lr_full (known + novel, native resolution)", pca_full_path)
    report.append(f"Figure: `{pca_full_path}`\n")
    if umap_full is not None:
        umap_full_path = os.path.join(out_dir, f"umap_lr_full{suffix}.png")
        plot_two_color(umap_full, mat_full.index, ancestry_by_person, novel_carrier, "UMAP-1",
                        "UMAP-2", "UMAP -- lr_full (known + novel, native resolution)",
                        umap_full_path)
        report.append(f"Figure: `{umap_full_path}`\n")

    known_only_mask = np.array([not novel_carrier[pid] for pid in mat_full.index])
    report.append(f"\n## lr_known_only (same lr_full coordinates, subset)\n"
                   f"{known_only_mask.sum()} of {len(mat_full)} people carry zero novel calls "
                   f"across these 8 genes -- plotted as the SAME PCA/UMAP points as lr_full above, "
                   f"filtered, never a re-fit.\n")
    pca_known_path = os.path.join(out_dir, f"pca_lr_known_only_highlight{suffix}.png")
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    ax.scatter(scores_full[~known_only_mask, 0], scores_full[~known_only_mask, 1], s=10,
               alpha=0.25, color="#cccccc", label="novel-carrier (background)", edgecolors="none")
    ax.scatter(scores_full[known_only_mask, 0], scores_full[known_only_mask, 1], s=10, alpha=0.6,
               color="#4C72B0", label="known-only", edgecolors="none")
    ax.set_xlabel(f"PC1 ({100*var_full[0]:.1f}% var)")
    ax.set_ylabel(f"PC2 ({100*var_full[1]:.1f}% var)")
    ax.set_title("lr_known_only highlighted within the lr_full embedding", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    vc.savefig(fig, pca_known_path)
    report.append(f"Figure: `{pca_known_path}`\n")

    # --- 3: lr_collapsed (separately fit, same people as lr_full where possible) -------------
    print("Building lr_collapsed (2-field) dosage matrix ...", file=sys.stderr)
    mat_lr2, novel_carrier_lr2, n_incomplete_lr2 = build_dosage_matrix(
        table1, set(mat_full.index), GENES, "collapsed")
    report.append(
        f"\n## lr_collapsed\nComplete 2-field cases among lr_full's {len(mat_full)} people: "
        f"{len(mat_lr2)} ({n_incomplete_lr2} incomplete under 2-field truncation). This drop is "
        f"EXPECTED to be substantial, not a bug: a `protein_altering` novel call (depth-2 novelty, "
        f"reference/IMMUANNOT_GTF_SPEC.md part B) only has its FIRST field resolved before 'new' is "
        f"spliced in -- `to_nfield(consensus, 2)` needs 2 resolved fields and returns None for it. "
        f"Since protein-altering is the dominant novelty class (~91% of real novel calls, "
        f"03_novel_alleles.py's report), 2-field collapsing doesn't just blur these calls into a "
        f"coarser bucket -- it makes the PERSON'S ENTIRE 8-gene case incomplete and drops them from "
        f"this embedding outright. That is itself a quantitative answer to 'how much does "
        f"collapsing cost': compare this count directly against lr_full's zero incomplete cases.\n")
    scores_lr2, var_lr2, umap_lr2 = fit_embedding(mat_lr2, args.seed, args.umap_neighbors,
                                                   args.umap_min_dist)
    pca_lr2_path = os.path.join(out_dir, f"pca_lr_collapsed{suffix}.png")
    plot_two_color(scores_lr2[:, :2], mat_lr2.index, ancestry_by_person, novel_carrier_lr2,
                    f"PC1 ({100*var_lr2[0]:.1f}% var)", f"PC2 ({100*var_lr2[1]:.1f}% var)",
                    "PCA -- lr_collapsed (2-field)", pca_lr2_path)
    report.append(f"Figure: `{pca_lr2_path}`\n")
    if umap_lr2 is not None:
        umap_lr2_path = os.path.join(out_dir, f"umap_lr_collapsed{suffix}.png")
        plot_two_color(umap_lr2, mat_lr2.index, ancestry_by_person, novel_carrier_lr2, "UMAP-1",
                        "UMAP-2", "UMAP -- lr_collapsed (2-field)", umap_lr2_path)
        report.append(f"Figure: `{umap_lr2_path}`\n")

    # --- 4: sr_collapsed, same people as lr_collapsed where present in both cohorts ----------
    shared_people = set(mat_lr2.index) & sr_people
    sr_stats_note = ""
    if shared_people:
        print(f"Building sr_collapsed dosage matrix ({len(shared_people)} shared people) ...",
              file=sys.stderr)
        sr_long = vc.load_sr_genotypes(args.sr_genotypes or os.path.join(DEFAULT_DATA_ROOT,
                                                                          "hla_genotypes.tsv"))
        mat_sr2, n_incomplete_sr2 = build_sr_dosage_matrix(sr_long, shared_people, GENES)
        report.append(f"\n## sr_collapsed\nOf {len(shared_people)} people shared between lr and sr: "
                       f"{len(mat_sr2)} have a complete 2-field SR case "
                       f"({n_incomplete_sr2} incomplete/no-call, dropped).\n")
        scores_sr2, var_sr2, umap_sr2 = fit_embedding(mat_sr2, args.seed, args.umap_neighbors,
                                                       args.umap_min_dist)
        pca_sr2_path = os.path.join(out_dir, f"pca_sr_collapsed{suffix}.png")
        anc_only = {pid: ancestry_by_person.get(pid) for pid in mat_sr2.index}
        plot_two_color(scores_sr2[:, :2], mat_sr2.index, anc_only, {pid: False for pid in
                       mat_sr2.index}, f"PC1 ({100*var_sr2[0]:.1f}% var)",
                       f"PC2 ({100*var_sr2[1]:.1f}% var)", "PCA -- sr_collapsed (2-field, "
                       "AoU-native short-read)", pca_sr2_path)
        report.append(f"Figure: `{pca_sr2_path}`\n")

        # Collapse-radius: same people, LR-collapsed vs SR-collapsed concentration stats.
        shared_both = set(mat_lr2.index) & set(mat_sr2.index)
        lr2_shared = mat_lr2.loc[mat_lr2.index.isin(shared_both)]
        sr2_shared = mat_sr2.loc[mat_sr2.index.isin(shared_both)]
        stats_lr = concentration_stats(lr2_shared, GENES).set_index("gene")
        stats_sr = concentration_stats(sr2_shared, GENES).set_index("gene")
        report.append(
            f"\n## Collapse-radius comparison: LR-collapsed vs SR-collapsed, SAME "
            f"{len(shared_both)} people, both at 2-field resolution\n"
            "Hypothesis (2026-09 chat discussion): SR, having less intrinsic resolution/confidence, "
            "calls concentrate more heavily on a few common alleles than LR does even at the SAME "
            "nominal 2-field resolution. `top10_share` = fraction of all allele copies covered by "
            "the 10 most common alleles for that gene; `simpson_concentration` = sum(p_i^2) "
            "(higher = more concentrated on few alleles).\n")
        report.append("| gene | LR top10 share | SR top10 share | LR Simpson | SR Simpson | "
                       "LR n_distinct | SR n_distinct |")
        report.append("|---|---|---|---|---|---|---|")
        for gene in GENES:
            if gene not in stats_lr.index or gene not in stats_sr.index:
                continue
            lr_r, sr_r = stats_lr.loc[gene], stats_sr.loc[gene]
            report.append(f"| {gene} | {100*lr_r['top10_share']:.1f}% | "
                           f"{100*sr_r['top10_share']:.1f}% | {lr_r['simpson_concentration']:.3f} | "
                           f"{sr_r['simpson_concentration']:.3f} | {int(lr_r['n_distinct_alleles'])} |"
                           f" {int(sr_r['n_distinct_alleles'])} |")
        sr_stats_note = "computed"
    else:
        report.append("\n## sr_collapsed\nNo people shared between the lr_collapsed complete-case "
                       "set and the sr cohort membership -- skipped.\n")

    md_path = os.path.join(out_dir, f"08_embedding_compare_report{suffix}.md")
    vc.write_report(md_path, report)
    print(f"Wrote report to {md_path!r}.", file=sys.stderr)
    if sr_stats_note:
        print("Collapse-radius comparison computed -- see report.", file=sys.stderr)


if __name__ == "__main__":
    main()
