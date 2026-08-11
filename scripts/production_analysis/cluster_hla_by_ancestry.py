#!/usr/bin/env python3
"""PCA / UMAP on HLA genotype calls, colored by ancestry -- exploring whether the calling
resolution supports dimensionality-reduction clustering (Marc's own open question, 2026-08-10).

## Encoding choice (the resolution question, answered explicitly)

HLA calls are categorical star-allele strings, not continuous values -- they can't go into
PCA/UMAP directly. This script uses an **allele-dosage encoding**, the same family of method used
in published HLA-vs-ancestry population-structure work and analogous to how genome-wide SNP PCA
encodes genotype dosage (0/1/2 copies) per site, just restricted to the 8 hyper-polymorphic HLA
loci instead of genome-wide markers:
  - Each unique (gene, 2-field allele) combination observed in the cohort becomes one feature
    column (e.g. "A:02:01", "DRB1:15:01", ...).
  - Each person's value in that column is 0, 1, or 2 -- how many of their two (unordered) allele
    calls at that gene equal that allele. hap1/hap2 are NOT treated as ordered/parental -- Immuannot's
    haplotype labels aren't guaranteed parent-of-origin-consistent across genes, so a genotype is
    just an unordered 2-allele set, same convention as every comparison script in this project.
  - **2-field resolution only** (not 3/4-field) -- keeps the allele vocabulary a manageable,
    non-sparse size. This mirrors the project's existing "Field 2 is the headline metric"
    convention (context/DECISIONS.md) and sidesteps the still-open "does deeper resolution help or
    just add noise" question (DECISIONS.md, "Resolution vs. generalization") rather than
    presupposing an answer to it here.
  - **Complete cases only**: a person needs a real call at all 8 classical genes (both
    haplotypes) to be included -- no imputation. This is a real constraint on N, reported
    explicitly (ties back to analyze_completeness_and_demographics.py's completeness numbers).

**Honest caveat, stated not hidden:** 8 loci is a tiny feature set next to genome-wide SNP PCA
(hundreds of thousands of markers) -- don't expect textbook-clean ancestry separation. But HLA is
under strong selection with well-documented ancestry-correlated allele frequencies (this project's
own DPA1 homozygosity-by-ancestry finding, EXPERIMENTS.md 2026-07-10), so real structure is a
reasonable thing to look for, not a long shot.

Columns are standardized (mean-centered, unit variance) before PCA/UMAP -- a standard, defensible
prep, not a rigorous population-genetics eigen-decomposition (e.g. not Patterson/EIGENSOFT allele-
frequency normalization) -- adequate for exploratory clustering, flagged so it isn't overclaimed as
more rigorous than it is.

Reads: ~/pipeline_outputs/immuannot_calls.tsv, ~/pipeline_outputs/immuannot_cohort_full.tsv
Writes: markdown + PNGs under --out-dir. Point-level PCA/UMAP scatter plots are one dot per real
person (no identifiers, but structurally person-level, unlike this project's usual pure-aggregate
outputs) -- this is standard AoU Workbench research-figure practice (AoU itself publishes ancestry
PCA plots the same way), fine to build/view inside the Workbench; treat like any other research
figure for the sponsor-review question already open in context/DECISIONS.md if it's ever shared
externally.

Usage (via `pixi run -e spechla`; needs scikit-learn [already a dependency] and, for the UMAP
panel specifically, umap-learn [added 2026-08-10, see pixi.toml] -- PCA still runs fine without it):
  python3 scripts/production_analysis/cluster_hla_by_ancestry.py
"""
import argparse
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import umap
    HAVE_UMAP = True
except ImportError:
    HAVE_UMAP = False

GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]  # "NA" dropped from this analysis
ANCESTRY_COLOR = {
    "AFR": "#4C72B0", "AMR": "#DD8452", "EAS": "#55A868",
    "EUR": "#C44E52", "MID": "#8172B2", "SAS": "#937860",
}
NULL = {"", "NA", "nan", "None", ".", "-"}

DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")

# Confidence tiers for the "are the mini-clusters a confidence artifact?" experiments (Marc's own
# question, 2026-08-11) -- worst per-person template_distance (edit distance to nearest IMGT
# reference allele) across all 8 classical genes, both haplotypes.
DIST_TIER_ORDER = ["0 (exact)", "1", "2", "3-5", "6+", "no data"]
DIST_TIER_BOUNDS = [("0 (exact)", 0, 0), ("1", 1, 1), ("2", 2, 2), ("3-5", 3, 5), ("6+", 6, float("inf"))]


def to_2field(allele):
    if allele is None or (hasattr(pd, "isna") and pd.isna(allele)):
        return None
    s = str(allele).strip()
    if s in NULL:
        return None
    if "*" in s:
        s = s.split("*", 1)[1]
    fields = [f for f in s.split(":") if f != ""]
    fields = [f for f in fields if f.strip().lower() != "new"]
    if len(fields) < 2:
        return None
    return f"{fields[0]}:{fields[1]}"


def load_calls(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for c in ["person_id", "gene", "immuannot_1", "immuannot_2"]:
        if c not in df.columns:
            sys.exit(f"FATAL: {path} missing column '{c}'. Actual: {list(df.columns)}")
    df["gene_bare"] = df["gene"].str.replace("^HLA-", "", regex=True)
    filtered = df[df["gene_bare"].isin(GENES)]
    if filtered.empty:
        seen = sorted(df["gene_bare"].unique())[:20]
        sys.exit(f"FATAL: 0 of {len(df)} rows in {path} match any of the 8 classical genes "
                 f"{GENES}. Genes actually present (first 20): {seen}. This is the exact "
                 f"symptom of the 2026-08-10 merge_fragments() dedup bug -- run "
                 f"scripts/production_orchestrator/rebuild_immuannot_calls.py first.")
    return filtered


def load_ancestry(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    for c in ["person_id", "ancestry_pred"]:
        if c not in df.columns:
            sys.exit(f"FATAL: {path} missing column '{c}'. Actual: {list(df.columns)}")
    # Real production file uses lowercase ("afr", "amr", ...) -- found 2026-08-10.
    df["ancestry_pred"] = df["ancestry_pred"].str.upper()
    return dict(zip(df["person_id"], df["ancestry_pred"]))


def build_dosage_matrix(calls):
    """Returns (DataFrame[person_id x feature], list of person_ids kept, n_excluded_incomplete)."""
    per_person = {}
    for _, r in calls.iterrows():
        pid, gene = r["person_id"], r["gene_bare"]
        a1, a2 = to_2field(r["immuannot_1"]), to_2field(r["immuannot_2"])
        per_person.setdefault(pid, {})[gene] = (a1, a2)

    complete, incomplete = [], 0
    rows = {}
    for pid, gene_map in per_person.items():
        if len(gene_map) < len(GENES) or any(v[0] is None or v[1] is None for v in gene_map.values()):
            incomplete += 1
            continue
        counter = Counter()
        for gene in GENES:
            a1, a2 = gene_map[gene]
            counter[f"{gene}:{a1}"] += 1
            counter[f"{gene}:{a2}"] += 1
        complete.append(pid)
        rows[pid] = counter

    feature_names = sorted({feat for c in rows.values() for feat in c})
    mat = pd.DataFrame(0, index=complete, columns=feature_names, dtype=float)
    for pid, counter in rows.items():
        for feat, n in counter.items():
            mat.at[pid, feat] = n
    return mat, complete, incomplete


def load_confidence(path):
    """Returns dict[(person_id, gene_bare)] -> worst_template_distance, or None if the cache
    doesn't exist (confidence experiments are skipped in that case, everything else still runs)."""
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t", dtype=str)
    for c in ["person_id", "gene", "worst_template_distance"]:
        if c not in df.columns:
            sys.exit(f"FATAL: {path} missing column '{c}'. Actual: {list(df.columns)}")
    df["gene_bare"] = df["gene"].str.replace("^HLA-", "", regex=True)
    df["worst_template_distance"] = pd.to_numeric(df["worst_template_distance"], errors="coerce")
    out = {}
    for pid, gene, td in zip(df["person_id"], df["gene_bare"], df["worst_template_distance"]):
        if pd.isna(td):
            continue
        out[(pid, gene)] = td
    return out


def person_max_distance(pid, confidence):
    vals = [confidence[(pid, g)] for g in GENES if (pid, g) in confidence]
    return max(vals) if vals else None


def bin_distance(d):
    if d is None:
        return "no data"
    for label, lo, hi in DIST_TIER_BOUNDS:
        if lo <= d <= hi:
            return label
    return "6+"


def duplicate_vector_report(mat):
    """Tests the competing, simpler explanation for the mini-clusters: HLA allele frequencies are
    heavily skewed and only 8 loci are used as features, so unrelated people can land on the
    *exact same* dosage vector by chance -- a literal tie in high-dim space, before any embedding
    or confidence question even enters the picture."""
    counts = mat.groupby(list(mat.columns), sort=False).size()
    dup_counts = counts[counts > 1]
    n_dup = int(dup_counts.sum())
    if n_dup == 0:
        return f"0 of {len(mat)} people share an identical dosage vector with anyone else."
    return (f"{n_dup} of {len(mat)} people ({100 * n_dup / len(mat):.1f}%) share an identical "
            f"16-call dosage vector with at least one other person -- largest identical group: "
            f"{int(dup_counts.max())} people landing on the exact same point before any embedding "
            f"runs. This alone can produce dense point-stacks that read as 'mini clusters' in 2D, "
            f"independent of call confidence.")


def plot_embedding_by_distance_tier(coords, tier_labels, xlabel, ylabel, title, out_path, extra_note=""):
    """Same idea as plot_embedding but colored by confidence tier (sequential, one hue, light to
    dark -- magnitude, not identity) instead of ancestry (categorical). 'no data' is neutral gray,
    outside the ramp. Higher-distance (less confident) tiers are drawn on top (higher zorder) so
    the rare, more-interesting points aren't hidden under the majority distance-0 mass."""
    fig, ax = plt.subplots(figsize=(8, 7))
    cmap = plt.cm.Blues
    n_ramp = len(DIST_TIER_ORDER) - 1  # exclude "no data" from the ramp itself
    for i, tier in enumerate(DIST_TIER_ORDER):
        mask = [t == tier for t in tier_labels]
        if not any(mask):
            continue
        pts = coords[mask]
        color = "#B0B0B0" if tier == "no data" else cmap(0.35 + 0.55 * i / max(n_ramp - 1, 1))
        zorder = 1 if tier in ("0 (exact)", "no data") else 5 + i
        ax.scatter(pts[:, 0], pts[:, 1], s=10, alpha=0.6, color=color,
                   label=f"{tier} (n={sum(mask)})", edgecolors="none", zorder=zorder)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title + (f"\n{extra_note}" if extra_note else ""), fontsize=10)
    ax.legend(fontsize=8, frameon=False, markerscale=1.5,
              title="worst template_distance\n(any of the 8 genes)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def plot_embedding(coords, ancestry_labels, xlabel, ylabel, title, out_path, extra_note=""):
    fig, ax = plt.subplots(figsize=(8, 7))
    for anc in ANCESTRY_ORDER:
        mask = [a == anc for a in ancestry_labels]
        if not any(mask):
            continue
        pts = coords[mask]
        ax.scatter(pts[:, 0], pts[:, 1], s=10, alpha=0.55, color=ANCESTRY_COLOR[anc],
                   label=f"{anc} (n={sum(mask)})", edgecolors="none")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title + (f"\n{extra_note}" if extra_note else ""), fontsize=10)
    ax.legend(fontsize=8, frameon=False, markerscale=1.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", default=os.path.join(DEFAULT_OUTROOT, "immuannot_calls.tsv"))
    ap.add_argument("--cohort", default=os.path.join(DEFAULT_OUTROOT, "immuannot_cohort_full.tsv"))
    ap.add_argument("--confidence", default=os.path.join(DEFAULT_OUTROOT, "immuannot_confidence.tsv"),
                    help="Per-(person,gene) worst_template_distance cache written by "
                         "rebuild_immuannot_calls.py -- used for the confidence-collapsing "
                         "experiments below. If missing, those experiments are skipped and "
                         "everything else still runs.")
    ap.add_argument("--out-dir", default=os.path.join(DEFAULT_OUTROOT, "production_analysis", "clustering"))
    ap.add_argument("--umap-neighbors", type=int, default=30)
    ap.add_argument("--umap-min-dist", type=float, default=0.25)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading calls + ancestry...", file=sys.stderr)
    calls = load_calls(args.calls)
    ancestry = load_ancestry(args.cohort)

    print("Building allele-dosage matrix (complete cases only)...", file=sys.stderr)
    mat, kept_ids, n_incomplete = build_dosage_matrix(calls)
    print(f"  {len(kept_ids)} people with complete 8-gene calls kept "
          f"({n_incomplete} excluded for incomplete calls); "
          f"{mat.shape[1]} distinct (gene, 2-field allele) features", file=sys.stderr)

    labels = [ancestry.get(pid, "NA") for pid in kept_ids]
    keep_mask = [a in ANCESTRY_ORDER for a in labels]
    mat = mat.loc[[pid for pid, k in zip(kept_ids, keep_mask) if k]]
    labels = [a for a, k in zip(labels, keep_mask) if k]
    print(f"  {len(labels)} of those have a usable genetic-ancestry label (AFR/AMR/EAS/EUR/MID/SAS)",
          file=sys.stderr)

    if len(mat) < 20:
        sys.exit(f"FATAL: only {len(mat)} complete, ancestry-labeled people -- too few for a "
                 f"meaningful embedding. Check that immuannot_calls.tsv / immuannot_cohort_full.tsv "
                 f"point at the real production output.")

    X = StandardScaler().fit_transform(mat.values)

    print("Running PCA...", file=sys.stderr)
    pca = PCA(n_components=min(10, X.shape[0], X.shape[1]), random_state=args.seed)
    pcs = pca.fit_transform(X)
    var_expl = pca.explained_variance_ratio_

    pca_path = os.path.join(args.out_dir, "pca_pc1_pc2_by_ancestry.png")
    plot_embedding(pcs[:, :2], labels,
                   f"PC1 ({100 * var_expl[0]:.1f}% var)", f"PC2 ({100 * var_expl[1]:.1f}% var)",
                   "PCA of HLA allele-dosage, colored by genetic ancestry", pca_path)

    sil_pca = silhouette_score(pcs[:, :2], labels) if len(set(labels)) > 1 else float("nan")

    umap_path = None
    sil_umap = float("nan")
    if HAVE_UMAP:
        print("Running UMAP...", file=sys.stderr)
        reducer = umap.UMAP(n_neighbors=args.umap_neighbors, min_dist=args.umap_min_dist,
                            random_state=args.seed)
        emb = reducer.fit_transform(X)
        umap_path = os.path.join(args.out_dir, "umap_by_ancestry.png")
        plot_embedding(emb, labels, "UMAP-1", "UMAP-2",
                       "UMAP of HLA allele-dosage, colored by genetic ancestry",
                       umap_path, extra_note=f"n_neighbors={args.umap_neighbors}, min_dist={args.umap_min_dist}")
        sil_umap = silhouette_score(emb, labels) if len(set(labels)) > 1 else float("nan")
    else:
        print("  umap-learn not installed -- skipping UMAP panel (PCA above still valid). "
              "Add `umap-learn` to pixi.toml's spechla deps and `pixi install -e spechla` to enable it.",
              file=sys.stderr)

    # ---- Mini-cluster investigation (Marc's own question, 2026-08-11): is the sparsified,
    # blobby look of the embedding a confidence artifact (calls getting snapped to the nearest
    # IMGT reference allele), or just HLA's skewed allele frequencies + only 8 loci producing
    # literal duplicate genotypes? Test both explicitly instead of assuming either.
    print("\nChecking for literal duplicate dosage vectors (skewed-allele-frequency hypothesis)...",
          file=sys.stderr)
    dup_report = duplicate_vector_report(mat)
    print(f"  {dup_report}", file=sys.stderr)

    conf_section = [dup_report]
    confidence = load_confidence(args.confidence)
    if confidence is None:
        conf_section.append(f"`{args.confidence}` not found -- confidence-collapsing experiments "
                             f"skipped (run rebuild_immuannot_calls.py's real repair to generate it).")
        print(f"  {args.confidence} not found -- skipping the confidence-collapsing experiments.",
              file=sys.stderr)
    else:
        person_td = {pid: person_max_distance(pid, confidence) for pid in mat.index}
        td_tiers = [bin_distance(person_td[pid]) for pid in mat.index]

        # Experiment A ("map first, mask after"): keep the SAME full-cohort embedding coordinates
        # computed above, just recolor by confidence tier instead of ancestry. If low-confidence
        # people concentrate in specific blobs, that's evidence for the collapsing hypothesis.
        plot_embedding_by_distance_tier(
            pcs[:, :2], td_tiers, f"PC1 ({100 * var_expl[0]:.1f}% var)", f"PC2 ({100 * var_expl[1]:.1f}% var)",
            "Same PCA embedding as above, recolored by call confidence",
            os.path.join(args.out_dir, "pca_by_confidence_posthoc.png"))
        if HAVE_UMAP:
            plot_embedding_by_distance_tier(
                emb, td_tiers, "UMAP-1", "UMAP-2",
                "Same UMAP embedding as above, recolored by call confidence",
                os.path.join(args.out_dir, "umap_by_confidence_posthoc.png"),
                extra_note=f"n_neighbors={args.umap_neighbors}, min_dist={args.umap_min_dist}")

        # Experiment B ("filter first, remap"): keep only people who are distance-0 (exact
        # reference match) on all 8 genes, refit PCA/UMAP from scratch on that subset only.
        confident_mask = [t == "0 (exact)" for t in td_tiers]
        n_confident = sum(confident_mask)
        print(f"  {n_confident} of {len(mat)} people are distance-0 on all 8 genes -- refitting "
              f"PCA/UMAP on that subset only...", file=sys.stderr)
        conf_section.append(f"{n_confident} of {len(mat)} people ({100 * n_confident / len(mat):.1f}%) "
                             f"are distance-0 (exact IMGT reference match) on all 8 genes.")
        if n_confident >= 20:
            mat_c = mat.loc[confident_mask]
            labels_c = [l for l, k in zip(labels, confident_mask) if k]
            dup_report_c = duplicate_vector_report(mat_c)
            conf_section.append(f"Confident-only subset duplicate check: {dup_report_c}")

            X_c = StandardScaler().fit_transform(mat_c.values)
            pca_c = PCA(n_components=min(10, X_c.shape[0], X_c.shape[1]), random_state=args.seed)
            pcs_c = pca_c.fit_transform(X_c)
            var_c = pca_c.explained_variance_ratio_
            plot_embedding(pcs_c[:, :2], labels_c,
                           f"PC1 ({100 * var_c[0]:.1f}% var)", f"PC2 ({100 * var_c[1]:.1f}% var)",
                           "PCA, distance-0-only people, colored by ancestry",
                           os.path.join(args.out_dir, "pca_confident_only.png"))
            sil_pca_c = silhouette_score(pcs_c[:, :2], labels_c) if len(set(labels_c)) > 1 else float("nan")
            conf_section.append(f"PCA (confident-only, refit from scratch): PC1+PC2 explain "
                                 f"{100 * (var_c[0] + var_c[1]):.1f}% of variance, silhouette "
                                 f"{sil_pca_c:.3f} (full-cohort PCA silhouette was {sil_pca:.3f}).")
            if HAVE_UMAP:
                reducer_c = umap.UMAP(n_neighbors=args.umap_neighbors, min_dist=args.umap_min_dist,
                                      random_state=args.seed)
                emb_c = reducer_c.fit_transform(X_c)
                plot_embedding(emb_c, labels_c, "UMAP-1", "UMAP-2",
                               "UMAP, distance-0-only people, colored by ancestry",
                               os.path.join(args.out_dir, "umap_confident_only.png"),
                               extra_note=f"n_neighbors={args.umap_neighbors}, min_dist={args.umap_min_dist}")
                sil_umap_c = silhouette_score(emb_c, labels_c) if len(set(labels_c)) > 1 else float("nan")
                conf_section.append(f"UMAP (confident-only, refit from scratch) silhouette: "
                                     f"{sil_umap_c:.3f} (full-cohort UMAP silhouette was {sil_umap:.3f}).")
        else:
            conf_section.append(f"Too few distance-0 people ({n_confident}) to refit a meaningful "
                                 f"embedding -- skipping the filter-first experiment.")

    md = [
        "# HLA calls x ancestry -- dimensionality reduction\n",
        f"N = {len(mat)} people with complete 8-gene calls and a usable ancestry label "
        f"(of {len(kept_ids)} complete-case people, {n_incomplete} excluded for incomplete calls). "
        f"{mat.shape[1]} (gene, 2-field allele) features.\n",
        f"PCA: PC1+PC2 explain {100 * (var_expl[0] + var_expl[1]):.1f}% of variance "
        f"(PC1 {100 * var_expl[0]:.1f}%, PC2 {100 * var_expl[1]:.1f}%). "
        f"Silhouette score (ancestry labels on PC1-2): {sil_pca:.3f} "
        f"(range -1 to 1; higher = ancestry groups separate more cleanly in this embedding).\n",
    ]
    if HAVE_UMAP:
        md.append(f"UMAP silhouette score (ancestry labels): {sil_umap:.3f}.\n")
    md.append(f"Ancestry breakdown of the analyzed set: "
              f"{dict(Counter(labels))}\n")
    md.append(f"\nFigures: `{pca_path}`" + (f", `{umap_path}`" if umap_path else " (UMAP skipped, umap-learn not installed)") + "\n")

    md.append("\n## Mini-cluster investigation (2026-08-11)\n")
    md.append("Two competing explanations for the sparsified, blobby (not continuous) look of the "
              "embedding above -- tested rather than assumed:\n")
    for line in conf_section:
        md.append(f"- {line}\n")
    if confidence is not None:
        md.append("\nFigures: `pca_by_confidence_posthoc.png`, "
                  + ("`umap_by_confidence_posthoc.png`, " if HAVE_UMAP else "")
                  + "`pca_confident_only.png`"
                  + (", `umap_confident_only.png`" if HAVE_UMAP else "") + "\n")

    md_text = "\n".join(md)
    md_path = os.path.join(args.out_dir, "clustering_report.md")
    with open(md_path, "w") as f:
        f.write(md_text)
    print(md_text)
    print(f"\n(written to {md_path})", file=sys.stderr)


if __name__ == "__main__":
    main()
