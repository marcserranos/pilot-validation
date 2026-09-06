#!/usr/bin/env python3
"""Population-structure figures, reproducible across all three cohorts (`--cohort sr|lr|lr_td`).
Implements VIZ_LIT.md 1.6, 1.7, and Part 2 (the continuous-ancestry differentiator):

  1. PCA of HLA allele-dosage, colored by discrete ancestry (VIZ_LIT.md 1.6). UMAP degrades
     gracefully to PCA-only if umap-learn isn't installed -- matches
     scripts/production_analysis/cluster_hla_by_ancestry.py's existing behavior exactly.
  2. Fst-like genetic-distance matrix + dendrogram between ancestry groups (VIZ_LIT.md 1.7), using
     a Hudson-estimator-style Fst per classical gene, averaged.
  3. Ancestry-continuum scatter: allele carrier status plotted against two CONTINUOUS admixture
     proportions (`p_*` in Table 4) instead of binned ancestry groups (VIZ_LIT.md 2.3) -- the
     figure the task calls out explicitly ("carrier status against the ancestry continuum rather
     than binned groups").
  4. Ternary plot of a 3-component admixture subset (renormalized), colored by carrier status
     (VIZ_LIT.md 1.4/2.4).
  5. ADMIXTURE-style stacked "barcode" plot of the full continuous 6-way probabilities
     (VIZ_LIT.md 1.5), downsampled above 5,000 people per its own stated failure mode by default
     (override with --admixture-max-people; --barcode-only skips every other figure in this
     script, for cheaply re-running just this one against a large cohort like sr).

Encoding choice for PCA/UMAP mirrors cluster_hla_by_ancestry.py exactly (2-field allele-dosage,
complete 8-classical-gene cases only, unordered hap1/hap2, no imputation) so this script
supersedes rather than diverges from that precedent -- extended here to work across all three
cohorts and to use manual (sklearn-free) PCA per this sub-project's stated dependency list
(pandas/numpy/scipy/matplotlib/seaborn-optional/umap-learn-optional, no scikit-learn).

Per-person points in the PCA/UMAP/continuum scatter are UNLABELED aggregate research figures --
same compliance posture as cluster_hla_by_ancestry.py and AoU's own ancestry PCA plots (task
instructions explicitly carve this out as acceptable precedent). No bare person_id is ever
written to any file under reports/.

Usage (fixtures): see 05_figures_frequency.py's docstring for the fixture pipeline; then e.g.
    python3 scripts/hla_popgen/06_figures_structure.py --outroot /tmp/hla_fixtures \\
        --table1 /tmp/hla_fixtures/hla_calls_rich.sample.tsv \\
        --cohort-membership /tmp/hla_fixtures/cohort_membership.sample.tsv \\
        --sr-genotypes /tmp/hla_fixtures/hla_genotypes.tsv --cohort lr

Real run (VM): python3 scripts/hla_popgen/06_figures_structure.py --cohort lr
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


# ---------------------------------------------------------------------------
# Dosage matrix (same convention as cluster_hla_by_ancestry.py: unordered genotype, complete
# 8-gene cases only, 2-field)
# ---------------------------------------------------------------------------
def build_dosage_matrix(args, cohort_people):
    person_ids = set(cohort_people["person_id"])
    per_person = {}
    if args.cohort == "sr":
        sr_path = args.sr_genotypes or os.path.join(args.outroot, "hla_genotypes.tsv")
        sr_long = vc.load_sr_genotypes(sr_path)
        sub = sr_long[sr_long["person_id"].isin(person_ids) & sr_long["gene_bare"].isin(GENES)]
        for pid, gene, copy, allele in zip(sub["person_id"], sub["gene_bare"], sub["copy"],
                                            sub["allele"]):
            a2 = vc.to_nfield(allele, 2)
            per_person.setdefault(pid, {}).setdefault(gene, {})[copy] = a2
    else:
        table1_path = args.table1 or os.path.join(args.outroot, "hla_calls_rich.tsv")
        table1 = vc.load_table1(table1_path)
        sub = table1[table1["person_id"].isin(person_ids) & table1["gene_bare"].isin(GENES)]
        # Two hap rows per gene = the two copies; use an incrementing slot per (person, gene).
        for pid, gene, consensus in zip(sub["person_id"], sub["gene_bare"], sub["consensus"]):
            a2 = vc.to_nfield(consensus, 2)
            slots = per_person.setdefault(pid, {}).setdefault(gene, {})
            slots[len(slots) + 1] = a2

    complete, incomplete = [], 0
    rows = {}
    for pid, gene_map in per_person.items():
        if len(gene_map) < len(GENES):
            incomplete += 1
            continue
        vals = []
        ok = True
        for gene in GENES:
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
        counter = Counter(f"{g}:{a}" for g, a in vals)
        complete.append(pid)
        rows[pid] = counter

    feature_names = sorted({feat for c in rows.values() for feat in c})
    mat = pd.DataFrame(0.0, index=complete, columns=feature_names)
    for pid, counter in rows.items():
        for feat, n in counter.items():
            mat.at[pid, feat] = float(n)
    return mat, incomplete


def plot_embedding(coords, labels, xlabel, ylabel, title, out_path, extra_note=""):
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for anc in vc.ANCESTRY_ORDER:
        mask = [l == anc for l in labels]
        if not any(mask):
            continue
        pts = coords[mask]
        ax.scatter(pts[:, 0], pts[:, 1], s=10, alpha=0.55, color=vc.ANCESTRY_COLORS[anc],
                   label=f"{anc} (n={sum(mask)})", edgecolors="none")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title + (f"\n{extra_note}" if extra_note else ""), fontsize=10)
    ax.legend(fontsize=8, frameon=False, markerscale=1.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    vc.savefig(fig, out_path)


def run_pca_umap(mat, labels, out_dir, seed, umap_neighbors, umap_min_dist):
    X = vc.standardize(mat.values)
    scores, var_ratio = vc.pca_fit_transform(X, n_components=10)
    pca_path = os.path.join(out_dir, "pca_pc1_pc2_by_ancestry.png")
    plot_embedding(scores[:, :2], labels, f"PC1 ({100 * var_ratio[0]:.1f}% var)",
                   f"PC2 ({100 * var_ratio[1]:.1f}% var)",
                   "PCA of HLA allele-dosage, colored by ancestry", pca_path)
    sil_pca = vc.silhouette_like(scores[:, :2], labels)

    umap_path, sil_umap = None, np.nan
    if vc.HAVE_UMAP:
        import umap
        reducer = umap.UMAP(n_neighbors=umap_neighbors, min_dist=umap_min_dist, random_state=seed)
        emb = reducer.fit_transform(X)
        umap_path = os.path.join(out_dir, "umap_by_ancestry.png")
        plot_embedding(emb, labels, "UMAP-1", "UMAP-2",
                       "UMAP of HLA allele-dosage, colored by ancestry", umap_path,
                       extra_note=f"n_neighbors={umap_neighbors}, min_dist={umap_min_dist}")
        sil_umap = vc.silhouette_like(emb, labels)
    else:
        print("  umap-learn not installed -- skipping UMAP panel (PCA above still valid).",
              file=sys.stderr)
    return pca_path, umap_path, var_ratio, sil_pca, sil_umap


# ---------------------------------------------------------------------------
# Fst-like genetic distance + dendrogram (VIZ_LIT.md 1.7)
# ---------------------------------------------------------------------------
def hudson_fst(freq_a, freq_b, n_a, n_b):
    """Hudson's Fst estimator, averaged over alleles at one locus: 1 - Hs/Ht style but using the
    unbiased per-allele numerator/denominator (Hudson, Slatkin & Maddison 1992; Bhatia et al.
    2013 formulation), which handles unequal sample sizes across ancestry groups better than a
    naive Wright Fst -- flagged as necessary in VIZ_LIT.md 1.7 given AoU's very unequal ancestry
    group sizes. freq_a/freq_b: aligned arrays of per-allele frequency; n_a/n_b: number of
    chromosomes (2 x people) sampled in each group at this locus."""
    num_sum = 0.0
    den_sum = 0.0
    for pa, pb in zip(freq_a, freq_b):
        n_num = (pa - pb) ** 2 - pa * (1 - pa) / max(n_a - 1, 1) - pb * (1 - pb) / max(n_b - 1, 1)
        n_den = pa * (1 - pb) + pb * (1 - pa)
        num_sum += n_num
        den_sum += n_den
    if den_sum == 0:
        return np.nan
    return num_sum / den_sum


def compute_fst_matrix(long_df, genes, min_cell_n):
    """Per-locus Hudson Fst averaged across genes ('ratio of averages' -- standard practice),
    returned as a square ancestry x ancestry DataFrame. Cells involving a thin-N ancestry x gene
    combination are excluded from that gene's contribution rather than silently biasing the
    average."""
    ancs = vc.ANCESTRY_ORDER
    fst_by_pair = {(a, b): [] for a in ancs for b in ancs if a != b}
    for gene in genes:
        g = long_df[long_df["gene_bare"] == gene]
        if g.empty:
            continue
        counts = g.groupby(["allele_label", "ancestry_pred"]).size().unstack(fill_value=0)
        counts = counts.reindex(columns=ancs, fill_value=0)
        totals = counts.sum(axis=0)
        for i, a in enumerate(ancs):
            for b in ancs[i + 1:]:
                if totals.get(a, 0) < min_cell_n or totals.get(b, 0) < min_cell_n:
                    continue
                fa = (counts[a] / totals[a]).values
                fb = (counts[b] / totals[b]).values
                fst = hudson_fst(fa, fb, totals[a], totals[b])
                if not np.isnan(fst):
                    # Bug found against the fixtures: Hudson's estimator is UNBIASED but not
                    # bounded at 0 -- small samples (e.g. AMR/EAS/MID/SAS N in a td-filtered
                    # cohort) routinely push it slightly negative, and scipy.cluster.hierarchy.
                    # linkage() hard-rejects any negative distance ("Linkage 'Z' contains negative
                    # distances"), silently killing the dendrogram for exactly the small-N cohorts
                    # this figure most needs to work on. Floor at 0 -- standard population-genetics
                    # convention for a negative Fst point estimate (interpreted as "no detectable
                    # differentiation", not a real negative distance).
                    fst_by_pair[(a, b)].append(max(fst, 0.0))
                    fst_by_pair[(b, a)].append(max(fst, 0.0))
    mat = pd.DataFrame(np.nan, index=ancs, columns=ancs)
    for a in ancs:
        mat.loc[a, a] = 0.0
    for (a, b), vals in fst_by_pair.items():
        if vals:
            mat.loc[a, b] = float(np.mean(vals))
    return mat


def plot_fst(fst_mat, out_path):
    ancs = list(fst_mat.index)
    vals = fst_mat.fillna(0).values
    have_valid = np.isfinite(fst_mat.values.astype(float))
    fig, (ax_dendro, ax_heat) = plt.subplots(1, 2, figsize=(11, 5),
                                              gridspec_kw={"width_ratios": [1, 2]})
    order = list(range(len(ancs)))
    try:
        from scipy.cluster.hierarchy import linkage, dendrogram
        from scipy.spatial.distance import squareform
        d = np.nan_to_num(vals, nan=vals[have_valid].max() if have_valid.any() else 0)
        np.fill_diagonal(d, 0)
        d = (d + d.T) / 2
        Z = linkage(squareform(d, checks=False), method="average")
        dres = dendrogram(Z, labels=ancs, ax=ax_dendro, orientation="left",
                           color_threshold=0, above_threshold_color="#888888")
        order = [ancs.index(l) for l in dres["ivl"]][::-1]
    except Exception as e:  # pragma: no cover
        print(f"  WARNING: Fst dendrogram unavailable ({e}).", file=sys.stderr)
        ax_dendro.axis("off")

    ordered = [ancs[i] for i in order]
    disp = fst_mat.loc[ordered, ordered]
    im = ax_heat.imshow(disp.values.astype(float), cmap="magma")
    ax_heat.set_xticks(range(len(ordered)))
    ax_heat.set_xticklabels(ordered)
    ax_heat.set_yticks(range(len(ordered)))
    ax_heat.set_yticklabels(ordered)
    for i in range(len(ordered)):
        for j in range(len(ordered)):
            v = disp.values[i, j]
            if np.isfinite(v):
                ax_heat.text(j, i, f"{v:.3f}", ha="center", va="center", fontsize=8,
                             color="white" if v > np.nanmax(disp.values) * 0.5 else "black")
            else:
                ax_heat.text(j, i, "n/a", ha="center", va="center", fontsize=7, color="grey")
    ax_heat.set_title("Hudson Fst (averaged over classical genes)", fontsize=10)
    fig.colorbar(im, ax=ax_heat, fraction=0.04)
    ax_dendro.set_title("Clustering", fontsize=10)
    fig.tight_layout()
    vc.savefig(fig, out_path)


# ---------------------------------------------------------------------------
# Continuous-ancestry figures (VIZ_LIT.md Part 2)
# ---------------------------------------------------------------------------
def pick_target_allele(long_df, gene):
    g = long_df[long_df["gene_bare"] == gene]
    if g.empty:
        return None
    return g["allele_label"].value_counts().idxmax()


def carrier_status(long_df, cohort_people, gene, allele_label):
    g = long_df[(long_df["gene_bare"] == gene) & (long_df["allele_label"] == allele_label)]
    carriers = set(g["person_id"])
    out = cohort_people[["person_id"] + [f"p_{a.lower()}" for a in vc.ANCESTRY_ORDER]
                          + ["ancestry_pred"]].copy()
    out["carrier"] = out["person_id"].isin(carriers)
    return out.dropna(subset=[f"p_{a.lower()}" for a in vc.ANCESTRY_ORDER])


def plot_continuum_scatter(df, gene, allele_label, anc_x, anc_y, out_path):
    fig, ax = plt.subplots(figsize=(6.5, 6))
    xcol, ycol = f"p_{anc_x.lower()}", f"p_{anc_y.lower()}"
    non_carriers = df[~df["carrier"]]
    carriers = df[df["carrier"]]
    ax.scatter(non_carriers[xcol], non_carriers[ycol], s=8, alpha=0.25, color="#AAAAAA",
               label=f"non-carrier (n={len(non_carriers)})")
    ax.scatter(carriers[xcol], carriers[ycol], s=14, alpha=0.75, color="#D55E00",
               label=f"carrier (n={len(carriers)})")
    ax.set_xlabel(f"{anc_x} admixture proportion")
    ax.set_ylabel(f"{anc_y} admixture proportion")
    ax.set_title(f"{allele_label} carrier status vs. continuous ancestry "
                 f"({anc_x}/{anc_y} plane)\n(NOT binned ancestry groups)", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    vc.savefig(fig, out_path)


def plot_ternary(df, gene, allele_label, components, out_path):
    """Barycentric projection of a 3-component admixture subset (renormalized to sum 1),
    colored by carrier status. Manual implementation -- no python-ternary dependency."""
    cols = [f"p_{c.lower()}" for c in components]
    sub = df[cols].copy()
    sub["sum"] = sub[cols].sum(axis=1)
    sub = sub[sub["sum"] > 0]
    renorm = sub[cols].div(sub["sum"], axis=0)
    # Barycentric -> Cartesian: vertices at (0,0), (1,0), (0.5, sqrt(3)/2)
    verts = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]])
    xy = renorm.values @ verts
    carrier = df.loc[sub.index, "carrier"].values

    fig, ax = plt.subplots(figsize=(6.5, 6))
    tri = plt.Polygon(verts, fill=False, edgecolor="black", linewidth=1)
    ax.add_patch(tri)
    ax.scatter(xy[~carrier, 0], xy[~carrier, 1], s=8, alpha=0.25, color="#AAAAAA",
               label=f"non-carrier (n={int((~carrier).sum())})")
    ax.scatter(xy[carrier, 0], xy[carrier, 1], s=16, alpha=0.8, color="#D55E00",
               label=f"carrier (n={int(carrier.sum())})")
    for (vx, vy), name in zip(verts, components):
        ax.text(vx, vy - 0.05 if vy == 0 else vy + 0.03, name, ha="center", fontsize=10,
                fontweight="bold")
    ax.set_xlim(-0.15, 1.15)
    ax.set_ylim(-0.15, 1.0)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"{allele_label} carrier status in {'/'.join(components)} admixture simplex\n"
                 f"(renormalized 3-way subset of the full 6-way continuum)", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    fig.tight_layout()
    vc.savefig(fig, out_path)


def plot_admixture_barcode(cohort_people, out_path, max_people=5000, seed=0):
    """ADMIXTURE-style stacked barcode (VIZ_LIT.md 1.5). Downsampled above max_people per its own
    documented failure mode (overplotting). Sorted by predicted ancestry, then descending value
    of that component within each block, per VIZ_LIT.md's sorting-strategy warning."""
    df = cohort_people.dropna(subset=[f"p_{a.lower()}" for a in vc.ANCESTRY_ORDER] + ["ancestry_pred"])
    if len(df) > max_people:
        df = df.sample(n=max_people, random_state=seed)
    def sort_key(row):
        anc = row["ancestry_pred"]
        dominant = row.get(f"p_{anc.lower()}", 0) if anc in vc.ANCESTRY_ORDER else 0
        return (vc.ANCESTRY_ORDER.index(anc) if anc in vc.ANCESTRY_ORDER else 99, -dominant)
    df = df.assign(_key=df.apply(sort_key, axis=1)).sort_values("_key").drop(columns="_key")

    fig, ax = plt.subplots(figsize=(12, 3.5))
    bottom = np.zeros(len(df))
    x = np.arange(len(df))
    for anc in vc.ANCESTRY_ORDER:
        vals = df[f"p_{anc.lower()}"].values
        ax.bar(x, vals, bottom=bottom, width=1.0, color=vc.ANCESTRY_COLORS[anc], label=anc,
               linewidth=0)
        bottom += vals
    ax.set_xlim(0, len(df))
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_ylabel("admixture proportion")
    ax.set_title(f"Continuous 6-way admixture proportions per person (n={len(df)}"
                 f"{', downsampled' if len(cohort_people) > max_people else ''}), "
                 f"sorted by predicted ancestry then dominant fraction", fontsize=10)
    ax.legend(fontsize=8, ncol=6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.18))
    fig.tight_layout()
    vc.savefig(fig, out_path)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    vc.add_common_args(ap)
    ap.add_argument("--umap-neighbors", type=int, default=15)
    ap.add_argument("--umap-min-dist", type=float, default=0.25)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--continuum-gene", default="B",
                     help="Bare gene name used for the carrier-status-vs-continuum and ternary "
                          "figures (default HLA-B: highest-diversity classical locus).")
    ap.add_argument("--continuum-ancestries", default="AFR,EUR",
                     help="Comma-separated pair of ancestries for the 2D continuum scatter plane.")
    ap.add_argument("--ternary-ancestries", default="AFR,EUR,AMR",
                     help="Comma-separated triple of ancestries for the ternary plot.")
    ap.add_argument("--admixture-max-people", type=int, default=5000,
                     help="Downsample cap for the admixture barcode plot (VIZ_LIT.md 1.5's own "
                          "documented overplotting threshold). Pass a number >= cohort size "
                          "(e.g. a large value like 20000) to plot every person with no "
                          "downsampling.")
    ap.add_argument("--barcode-only", action="store_true",
                     help="Skip PCA/UMAP, Fst, and the continuum/ternary figures -- produce only "
                          "the admixture barcode. Use this for large cohorts (e.g. --cohort sr, "
                          "~500K people) where the other figures are expensive and not wanted.")
    args = ap.parse_args()

    cohort_df = vc.load_cohort_membership(
        args.cohort_membership or os.path.join(args.outroot, "cohort_membership.tsv"))
    cohort_people, cohort_label = vc.select_cohort_people(cohort_df, args.cohort, args.td_max)
    print(f"Cohort: {cohort_label} -- {len(cohort_people)} people.", file=sys.stderr)
    if len(cohort_people) < 20:
        sys.exit(f"FATAL: only {len(cohort_people)} people in cohort {cohort_label!r} -- too few.")

    out_dir = args.out_dir or vc.default_out_dir("06_figures_structure", args.cohort, args.td_max)
    vc.ensure_dir(out_dir)
    report = [f"# Population-structure figures -- cohort `{cohort_label}`\n",
              f"N people (cohort membership) = {len(cohort_people)}.\n"]
    fig_paths = []

    if args.barcode_only:
        barcode_path = os.path.join(out_dir, "admixture_barcode.png")
        plot_admixture_barcode(cohort_people, barcode_path, max_people=args.admixture_max_people,
                                seed=args.seed)
        report.append(f"\n## Figures (barcode-only run)\n- `{barcode_path}`")
        vc.write_report(os.path.join(out_dir, "structure_report.md"), report)
        return

    # ---- 1. PCA / UMAP on allele dosage --------------------------------------------------
    print("Building allele-dosage matrix (complete 8-gene cases only)...", file=sys.stderr)
    mat, n_incomplete = build_dosage_matrix(args, cohort_people)
    anc_map = dict(zip(cohort_people["person_id"], cohort_people["ancestry_pred"]))
    labels = [anc_map.get(pid) for pid in mat.index]
    # anc_map.get(pid) can be pd.NA/NaN for the small number of people with a missing/malformed
    # ancestry_pred (ENVIRONMENT.md quirk #30: ~24 blank/malformed rows in the real
    # immuannot_cohort_full.tsv, out of ~12,233 -- too rare to reliably appear in a 300-person
    # fixture, which is why this surfaced only at real scale). `pd.NA in list` raises TypeError
    # rather than returning False, because `in` evaluates `NA == x` for each x and that comparison
    # itself returns NA, not a bool. Guard with an explicit string check first so NA/NaN short-
    # circuits to "not kept" instead of reaching the `in` test at all.
    keep = [isinstance(l, str) and l in vc.ANCESTRY_ORDER for l in labels]
    mat_k = mat.loc[[pid for pid, k in zip(mat.index, keep) if k]]
    labels_k = [l for l, k in zip(labels, keep) if k]
    report.append(f"Allele-dosage matrix: {len(mat)} complete-case people ({n_incomplete} "
                   f"excluded for incomplete 8-gene calls), {len(mat_k)} with a usable ancestry "
                   f"label, {mat.shape[1]} (gene, 2-field allele) features.\n")

    if len(mat_k) >= 20 and mat_k.shape[1] >= 2:
        pca_path, umap_path, var_ratio, sil_pca, sil_umap = run_pca_umap(
            mat_k, labels_k, out_dir, args.seed, args.umap_neighbors, args.umap_min_dist)
        fig_paths += [p for p in (pca_path, umap_path) if p]
        report.append(f"PCA: PC1+PC2 explain {100 * (var_ratio[0] + var_ratio[1]):.1f}% of "
                       f"variance. Silhouette (ancestry labels, PC1-2): {sil_pca:.3f}.\n")
        if umap_path:
            report.append(f"UMAP silhouette (ancestry labels): {sil_umap:.3f}.\n")
        else:
            report.append("UMAP skipped -- umap-learn not installed; PCA above still valid.\n")
    else:
        report.append("Skipped PCA/UMAP: too few complete, ancestry-labeled people or too few "
                       "features for this cohort.\n")

    # ---- 2. Fst / genetic-distance dendrogram --------------------------------------------
    print("Computing Fst-like genetic distance between ancestry groups...", file=sys.stderr)
    long_df = _build_long_frame_local(args, cohort_people, GENES)
    if not long_df.empty:
        fst_mat = compute_fst_matrix(long_df, GENES, args.min_cell_n)
        fst_path = os.path.join(out_dir, "fst_dendrogram.png")
        plot_fst(fst_mat, fst_path)
        fig_paths.append(fst_path)
        report.append("\n## Pairwise Hudson Fst (averaged over classical genes)\n")
        report.append("| | " + " | ".join(vc.ANCESTRY_ORDER) + " |")
        report.append("|---|" + "---|" * len(vc.ANCESTRY_ORDER))
        for a in vc.ANCESTRY_ORDER:
            cells = [f"{fst_mat.loc[a, b]:.3f}" if np.isfinite(fst_mat.loc[a, b]) else "n/a"
                     for b in vc.ANCESTRY_ORDER]
            report.append(f"| {a} | " + " | ".join(cells) + " |")
    else:
        report.append("\nSkipped Fst matrix: no usable allele calls for this cohort.\n")

    # ---- 3/4/5. Continuous-ancestry figures ----------------------------------------------
    gene = args.continuum_gene
    allele_label = pick_target_allele(long_df, gene) if not long_df.empty else None
    if allele_label:
        df_carrier = carrier_status(long_df, cohort_people, gene, allele_label)
        ax_pair = args.continuum_ancestries.split(",")
        cont_path = os.path.join(out_dir, f"continuum_scatter_{gene}.png")
        plot_continuum_scatter(df_carrier, gene, allele_label, ax_pair[0], ax_pair[1], cont_path)
        fig_paths.append(cont_path)

        tri = args.ternary_ancestries.split(",")
        tern_path = os.path.join(out_dir, f"ternary_{gene}.png")
        plot_ternary(df_carrier, gene, allele_label, tri, tern_path)
        fig_paths.append(tern_path)

        report.append(f"\n## Continuous-ancestry figures (gene {gene}, allele {allele_label}, "
                       f"most common in this cohort)\n")
        report.append(f"- {int(df_carrier['carrier'].sum())} carriers of {len(df_carrier)} "
                       f"people with full admixture data.\n")
    else:
        report.append(f"\nSkipped continuum/ternary figures: no calls for gene {gene}.\n")

    barcode_path = os.path.join(out_dir, "admixture_barcode.png")
    plot_admixture_barcode(cohort_people, barcode_path, max_people=args.admixture_max_people,
                            seed=args.seed)
    fig_paths.append(barcode_path)

    report.append("\n## Figures\n")
    for p in fig_paths:
        report.append(f"- `{p}`")

    vc.write_report(os.path.join(out_dir, "structure_report.md"), report)


def _build_long_frame_local(args, cohort_people, genes):
    """Same shape as 05_figures_frequency.build_long_frame -- reimplemented locally (not
    imported cross-script) so 06 has no hard dependency on 05's internals."""
    person_ids = set(cohort_people["person_id"])
    anc_map = dict(zip(cohort_people["person_id"], cohort_people["ancestry_pred"]))
    if args.cohort == "sr":
        sr_path = args.sr_genotypes or os.path.join(args.outroot, "hla_genotypes.tsv")
        sr_long = vc.load_sr_genotypes(sr_path)
        sub = sr_long[sr_long["person_id"].isin(person_ids) & sr_long["gene_bare"].isin(genes)].copy()
        sub["allele_2field"] = sub["allele"].map(lambda a: vc.to_nfield(a, 2))
    else:
        table1_path = args.table1 or os.path.join(args.outroot, "hla_calls_rich.tsv")
        table1 = vc.load_table1(table1_path)
        sub = table1[table1["person_id"].isin(person_ids) & table1["gene_bare"].isin(genes)].copy()
        sub["allele_2field"] = sub["consensus"].map(lambda a: vc.to_nfield(a, 2))
    sub = sub.dropna(subset=["allele_2field"])
    sub["ancestry_pred"] = sub["person_id"].map(anc_map)
    sub = sub[sub["ancestry_pred"].isin(vc.ANCESTRY_ORDER)]
    sub["allele_label"] = sub["gene_bare"] + "*" + sub["allele_2field"]
    return sub[["person_id", "gene_bare", "allele_2field", "allele_label", "ancestry_pred"]]


if __name__ == "__main__":
    main()
