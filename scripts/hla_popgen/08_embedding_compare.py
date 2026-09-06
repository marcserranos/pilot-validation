#!/usr/bin/env python3
"""Does the enriched (known+novel) allele space look different once you collapse it back down --
and does it look different in the SAME coordinate frame, not just a similarly-shaped one fit
separately? Locked-mapping design from the 2026-09 chat discussion (git history has two earlier,
superseded versions of this script).

## The three columns (same people's coordinates, one fitted mapping, two mapped into it)

1. **`full_enriched`** -- LR cohort, 8 classical genes, complete diploid cases, each allele call at
   its own native resolution. The ONLY column that gets FIT (PCA and, per-variant, UMAP).
2. **`collapsed_nearest_ref`** -- ALL of the same people, every novel call replaced by its
   `nearest_allele` (Immuannot's own tied-best CDS-search candidate), standardized with
   `full_enriched`'s own mean/std and placed into the locked mapping via `.transform()`/retained
   PCA loadings -- a genuine projection, not a second fit.
3. **`known_only`** -- people with ZERO novel calls (whole-person drop). Already inside the
   column-1 fit; this column is just that fit's coordinates, filtered.

## Sanity check added 2026-09: is the locked mapping itself flattening real differences?

`collapsed_nearest_ref` looked almost indistinguishable from `full_enriched` in the raw-UMAP row of
the first real run -- worth asking whether that's because the underlying signal genuinely doesn't
move, or because forcing a `.transform()` into an already-fixed manifold can't help but look similar
regardless of what changed. `--sanity-check` fits `collapsed_nearest_ref` and `known_only`
INDEPENDENTLY (their own fresh PCA + UMAP fit, not a projection) and plots them next to the locked
versions -- if independent fits also look similar to `full_enriched`, that's evidence the
similarity is real, not a projection artifact; if they look different from each other while the
locked versions look artificially similar, that points the other way. See
`sanity_check/sanity_check_report.md` for the comparison once run.

## Why PC1/PC2 explain so little variance (worth having the number, not just the plot)

The dosage matrix has ~7,000+ columns (one per distinct allele identity actually observed across 8
genes) for ~9,000 people -- wide and extremely sparse. Standardizing every column to unit variance
before PCA (a deliberate, ordinary step) means a column only 3 people carry gets the same variance
weight as a column half the cohort carries, so total variance (~n_features units) is spread almost
evenly across thousands of columns -- any single PC capturing <1% of that total is the expected
arithmetic consequence of this encoding, not evidence the signal is weak. PCA-denoising before UMAP
still helps because CORRELATED signal (shared ancestry-driven allele co-occurrence across genes)
concentrates into the top components while independent per-column noise does not, even though each
top component is individually a small fraction of total variance -- see `matrix_report.md`'s "PCA
variance" section for the actual `n_samples x n_features` and per-component numbers from this run.

## Rows: coloring x whether UMAP ran on a PCA-denoised representation

PCA-then-UMAP (denoise to `--pca-denoise-k` components) is a PARALLEL EXPERIMENT, not a replacement
default -- both `raw` and `pca_denoised` are run and shown side by side. Colorings: ancestry
(categorical), %novel (continuous viridis gradient on fraction of 16 slots that are novel), and
(new, `--disease-alleles`) carrier status for a curated list of well-known HLA-disease-associated
alleles -- see `DISEASE_ALLELE_GROUPS` and its caveats below.

## Output layout

Everything lives under `<out-dir>/`: `matrix/` (the locked-mapping 3x4 grid + individual panels),
`sanity_check/` (independent-fit comparison), `disease/` (disease-allele coloring),
`collapse_radius/` (SR-vs-LR concentration diagnostic, unrelated question, kept from earlier).

## Caching

The expensive step (PCA + 2 UMAP fits on ~9,000 x ~7,000, plus 2 more for `--sanity-check`) is
cached to `<out-dir>/embedding_cache.pkl` after the first run. `--from-cache` skips straight to
plotting from that cache -- use it to iterate on figure layout or add a new coloring (like
`--disease-alleles`) without re-fitting UMAP.

Usage (fixtures): python3 scripts/hla_popgen/tests/make_fixtures.py --outroot /tmp/hla_fixtures -n 300
    then: python3 scripts/hla_popgen/08_embedding_compare.py --outroot /tmp/hla_fixtures \\
        --table1 /tmp/hla_fixtures/hla_calls_rich.sample.tsv \\
        --cohort-membership /tmp/hla_fixtures/cohort_membership.sample.tsv \\
        --sr-genotypes /tmp/hla_fixtures/hla_genotypes.tsv --sample

Real run (VM): python3 scripts/hla_popgen/08_embedding_compare.py --sanity-check --disease-alleles
"""
import argparse
import os
import pickle
import sys
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _viz_common as vc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

GENES = vc.CLASSICAL_GENES_BARE
N_SLOTS = len(GENES) * 2  # 8 genes x 2 haps
DEFAULT_DATA_ROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_OUTROOT = os.path.join(DEFAULT_DATA_ROOT, "people")

# Curated, NOT clinically authoritative: well-known HLA-disease associations from general
# immunogenetics literature/textbook knowledge, compiled for an exploratory population-genetics
# "do carriers cluster" check -- not verified against a citation per entry for this task, and NOT a
# diagnostic resource. (label, gene_bare, allele-group prefix to match against the 2-field-or-finer
# identity string, e.g. "B*27" matches "HLA-B*27:05:02:01"). Ordered by how well-established the
# association is; first MATCHING group wins if a person carries more than one (rare, logged).
DISEASE_ALLELE_GROUPS = [
    ("B*27 (spondyloarthritis)", "B", "B*27"),
    ("B*57:01 (abacavir hypersensitivity)", "B", "B*57:01"),
    ("B*58:01 (allopurinol SJS/TEN)", "B", "B*58:01"),
    ("B*15:02 (carbamazepine SJS/TEN)", "B", "B*15:02"),
    ("DRB1*15:01 (multiple sclerosis)", "DRB1", "DRB1*15:01"),
    ("DRB1*04 (rheumatoid arthritis)", "DRB1", "DRB1*04"),
    ("DQB1*06:02 (narcolepsy)", "DQB1", "DQB1*06:02"),
    ("DQB1*03:02 / DQ8 (T1D)", "DQB1", "DQB1*03:02"),
    ("DQA1*05:01 / DQ2 (celiac)", "DQA1", "DQA1*05:01"),
    ("B*51 (Behcet's disease)", "B", "B*51"),
    ("C*06:02 (psoriasis)", "C", "C*06:02"),
    ("A*29:02 (birdshot chorioretinopathy)", "A", "A*29:02"),
]
DISEASE_COLORS = ["#D55E00", "#E69F00", "#009E73", "#0072B2", "#CC79A7", "#56B4E9", "#8172B2",
                   "#C44E52", "#4C9F70", "#F2A541", "#2AACA4", "#8E6FCE"]
NONE_COLOR = "#d8dcd8"


# ---------------------------------------------------------------------------
# Identity assignment
# ---------------------------------------------------------------------------
def load_novel_matches(table1, outroot):
    """Re-derives the (person_id, hap, gene) -> novel_id map by importing 03_novel_alleles.py as a
    sibling module; also hands back the module so `nearest_allele_of()` can be reused verbatim."""
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
    return lookup, stats, mod


def build_full_and_collapsed(table1, person_ids, genes, novel_id_lookup, nearest_allele_mod):
    """One pass over Table 1 building BOTH the full_enriched identity per slot and its
    collapsed_nearest_ref fallback identity. Returns:
      per_person: {pid: {gene: {1: (full_ident, nearest_ident_or_None), 2: (...)}}}
      novel_carrier: {pid: bool}
      novel_slot_count: {pid: int}  -- out of N_SLOTS
    """
    sub = table1[table1["person_id"].isin(person_ids) & table1["gene_bare"].isin(genes)]
    per_person = {}
    novel_carrier = {}
    novel_slot_count = Counter()
    for _, row in sub.iterrows():
        pid, gene, hap = row["person_id"], row["gene_bare"], row["hap"]
        is_novel = bool(row["is_novel"])
        consensus = row["consensus"]
        if is_novel:
            novel_carrier[pid] = True
            novel_slot_count[pid] += 1
            full_ident = (novel_id_lookup or {}).get((pid, hap, vc.gene_display(gene)))
            nearest_ident = nearest_allele_mod.nearest_allele_of(row)
        else:
            novel_carrier.setdefault(pid, False)
            full_ident = None if (isinstance(consensus, str) and consensus == "undetermined") \
                else consensus
            nearest_ident = full_ident
        if full_ident is None:
            continue
        slots = per_person.setdefault(pid, {}).setdefault(gene, {})
        slots[len(slots) + 1] = (full_ident, nearest_ident)
    return per_person, novel_carrier, novel_slot_count


def matrices_from_slots(per_person, genes, novel_carrier):
    complete, incomplete = [], 0
    full_counters, collapsed_idents = {}, {}
    for pid, gene_map in per_person.items():
        if len(gene_map) < len(genes):
            incomplete += 1
            continue
        ok = True
        full_vals, collapsed_vals = [], []
        for gene in genes:
            copies = gene_map.get(gene, {})
            if len(copies) < 2 or any(copies.get(k) is None for k in (1, 2)):
                ok = False
                break
            for k in (1, 2):
                full_ident, nearest_ident = copies[k]
                full_vals.append((gene, full_ident))
                collapsed_vals.append((gene, nearest_ident))
        if not ok:
            incomplete += 1
            continue
        complete.append(pid)
        full_counters[pid] = Counter(f"{g}:{a}" for g, a in full_vals)
        collapsed_idents[pid] = collapsed_vals

    feature_names = sorted({feat for c in full_counters.values() for feat in c})
    mat_full = pd.DataFrame(0.0, index=complete, columns=feature_names)
    for pid, counter in full_counters.items():
        for feat, n in counter.items():
            mat_full.at[pid, feat] = float(n)
    return mat_full, collapsed_idents, complete, incomplete


def build_collapsed_matrix(mat_full_columns, collapsed_idents, complete):
    col_set = set(mat_full_columns)
    n_dropped_slots = 0
    n_total_slots = 0
    rows = {}
    for pid in complete:
        counter = Counter()
        for gene, ident in collapsed_idents[pid]:
            n_total_slots += 1
            feat = f"{gene}:{ident}"
            if feat in col_set:
                counter[feat] += 1
            else:
                n_dropped_slots += 1
        rows[pid] = counter
    mat = pd.DataFrame(0.0, index=complete, columns=list(mat_full_columns))
    for pid, counter in rows.items():
        for feat, n in counter.items():
            mat.at[pid, feat] = float(n)
    return mat, n_dropped_slots, n_total_slots


def load_real_disease_labels(path):
    """Reads 12_disease_phenotypes.py's person_disease_labels.tsv (VM-local, per-person, real
    EHR-confirmed diagnosis -- see that script's HLA_LINKED list, reused here verbatim for the
    color/legend order so this is the SAME 8 disease definitions, not a re-derivation). Returns
    ({pid: label_or_None}, {label: n}, n_multiple, label_order)."""
    df = pd.read_csv(path, sep="\t", dtype={"person_id": str})
    import importlib.util
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location(
        "hla_popgen_disease_phenotypes", os.path.join(here, "12_disease_phenotypes.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    label_order = [lbl for lbl, *_j in mod.HLA_LINKED]
    labels = {row.person_id: row.disease_label for row in df.itertuples()
              if isinstance(row.disease_label, str)}
    n_multiple = int((df["n_matches"] > 1).sum())
    counts = Counter(labels.values())
    return labels, counts, n_multiple, label_order


# ---------------------------------------------------------------------------
# Disease-allele carrier labels (PROXY, from allele identity -- see --real-disease-labels for the
# non-proxy version) -- computed directly from a dosage matrix's own columns, no need to touch
# Table 1 again (works identically on mat_full or mat_collapsed).
# ---------------------------------------------------------------------------
def disease_labels_from_matrix(mat):
    """Returns ({pid: group_label_or_None}, {group_label: n_carriers}, n_multiple)."""
    prefixes = [(label, f"{vc.gene_display(gene)}*{allele_group.split('*', 1)[1]}")
                for label, gene, allele_group in DISEASE_ALLELE_GROUPS]
    col_matches = {}  # label -> [columns]
    for label, prefix in prefixes:
        col_matches[label] = [c for c in mat.columns if c.split(":", 1)[-1].startswith(prefix)]
    labels = {}
    n_multiple = 0
    counts = Counter()
    for pid in mat.index:
        row = mat.loc[pid]
        matched = [label for label, cols in col_matches.items()
                   if cols and row[cols].sum() > 0]
        if not matched:
            continue
        if len(matched) > 1:
            n_multiple += 1
        labels[pid] = matched[0]
        counts[matched[0]] += 1
    return labels, counts, n_multiple


# ---------------------------------------------------------------------------
# PCA + UMAP fitting
# ---------------------------------------------------------------------------
def fit_pca_with_loadings(X, n_components=10):
    mu = X.mean(axis=0)
    Xc = X - mu
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    k = min(n_components, Vt.shape[0])
    scores = U[:, :k] * S[:k]
    var = (S ** 2) / max(X.shape[0] - 1, 1)
    var_ratio = var / var.sum()
    return scores, var_ratio[:k], mu, Vt[:k]


def transform_pca(X_new, mu, Vt_k):
    return (X_new - mu) @ Vt_k.T


def standardize_with(X, mu, sd):
    sd_safe = np.where(sd == 0, 1.0, sd)
    return (X - mu) / sd_safe


def fit_umap_pair(X_std, pca_scores, pca_denoise_k, seed, umap_neighbors, umap_min_dist):
    """One raw UMAP fit + one PCA-denoised UMAP fit, both FRESH (used for both the locked
    full_enriched fit and, under --sanity-check, the independent collapsed/known_only fits)."""
    if not vc.HAVE_UMAP:
        return None, None
    import umap
    raw = umap.UMAP(n_neighbors=umap_neighbors, min_dist=umap_min_dist,
                     random_state=seed).fit_transform(X_std)
    pca = umap.UMAP(n_neighbors=umap_neighbors, min_dist=umap_min_dist,
                     random_state=seed).fit_transform(pca_scores[:, :pca_denoise_k])
    return raw, pca


class LockedEmbedding:
    """PCA fit (with retained loadings) + raw UMAP fit + PCA-denoised UMAP fit on full_enriched,
    with `.transform()` available for projecting other matrices into the same frame."""
    def __init__(self, mat_full, seed, umap_neighbors, umap_min_dist, pca_denoise_k):
        self.mu = mat_full.values.mean(axis=0)
        self.sd = mat_full.values.std(axis=0)
        X = standardize_with(mat_full.values, self.mu, self.sd)
        self.pca_scores, self.pca_var, self.pca_mu, self.pca_Vt = fit_pca_with_loadings(
            X, n_components=max(pca_denoise_k, 10))
        self.pca_denoise_k = pca_denoise_k
        self.umap_raw_reducer = None
        self.umap_pca_reducer = None
        self.umap_raw_full = None
        self.umap_pca_full = None
        if vc.HAVE_UMAP:
            import umap
            self.umap_raw_reducer = umap.UMAP(n_neighbors=umap_neighbors, min_dist=umap_min_dist,
                                               random_state=seed)
            self.umap_raw_full = self.umap_raw_reducer.fit_transform(X)
            self.umap_pca_reducer = umap.UMAP(n_neighbors=umap_neighbors, min_dist=umap_min_dist,
                                               random_state=seed)
            self.umap_pca_full = self.umap_pca_reducer.fit_transform(
                self.pca_scores[:, :pca_denoise_k])

    def transform_raw(self, mat_new):
        X_new = standardize_with(mat_new.values, self.mu, self.sd)
        return None if self.umap_raw_reducer is None else self.umap_raw_reducer.transform(X_new)

    def transform_pca_denoised(self, mat_new):
        X_new = standardize_with(mat_new.values, self.mu, self.sd)
        pca_new = transform_pca(X_new, self.pca_mu, self.pca_Vt)[:, :self.pca_denoise_k]
        return None if self.umap_pca_reducer is None else self.umap_pca_reducer.transform(pca_new)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_panel(ax, coords, index, color_kind, ancestry_by_person=None, novel_frac=None,
               disease_labels=None, label_order=None, title=""):
    if coords is None:
        ax.text(0.5, 0.5, "umap-learn not installed", ha="center", va="center", fontsize=9,
                color="#999999", transform=ax.transAxes)
        ax.set_title(title, fontsize=9)
        ax.axis("off")
        return None
    sc = None
    if color_kind == "ancestry":
        anc_labels = [pid_anc if isinstance(pid_anc, str) else None
                      for pid_anc in (ancestry_by_person.get(pid) for pid in index)]
        for anc in vc.ANCESTRY_ORDER:
            mask = np.array([l == anc for l in anc_labels])
            if mask.any():
                ax.scatter(coords[mask, 0], coords[mask, 1], s=6, alpha=0.5,
                           color=vc.ANCESTRY_COLORS[anc], edgecolors="none")
    elif color_kind == "gradient":
        vals = np.array([novel_frac.get(pid, 0.0) for pid in index])
        sc = ax.scatter(coords[:, 0], coords[:, 1], s=6, alpha=0.6, c=vals, cmap="viridis",
                        vmin=0, vmax=1, edgecolors="none")
    elif color_kind == "disease":
        label_order = label_order if label_order is not None else \
            [lbl for lbl, *_j in DISEASE_ALLELE_GROUPS]
        labels = [disease_labels.get(pid) for pid in index]
        is_carrier = np.array([l is not None for l in labels])
        ax.scatter(coords[~is_carrier, 0], coords[~is_carrier, 1], s=5, alpha=0.25,
                   color=NONE_COLOR, edgecolors="none", zorder=1)
        color_map = {lbl: DISEASE_COLORS[i % len(DISEASE_COLORS)] for i, lbl in enumerate(label_order)}
        for lbl in label_order:
            mask = np.array([l == lbl for l in labels])
            if mask.any():
                ax.scatter(coords[mask, 0], coords[mask, 1], s=16, alpha=0.9,
                           color=color_map[lbl], edgecolors="none", zorder=2)
    ax.set_title(title, fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[["top", "right"]].set_visible(False)
    return sc


def build_matrix_figure(coords_by_cell, index_by_col, ancestry_by_person, novel_frac_by_col,
                         out_dir, suffix):
    cols = ["full_enriched", "collapsed_nearest_ref", "known_only"]
    col_titles = {"full_enriched": "full_enriched (fit)",
                  "collapsed_nearest_ref": "collapsed_nearest_ref (projected)",
                  "known_only": "known_only (subset of fit)"}
    rows = ["raw_ancestry", "raw_gradient", "pca_ancestry", "pca_gradient"]
    row_titles = {"raw_ancestry": "raw UMAP, ancestry", "raw_gradient": "raw UMAP, %novel",
                  "pca_ancestry": "PCA-denoised UMAP, ancestry",
                  "pca_gradient": "PCA-denoised UMAP, %novel"}
    matrix_dir = os.path.join(out_dir, "matrix")
    os.makedirs(matrix_dir, exist_ok=True)

    # Layout fix: reserve a dedicated top strip for legend/title (via height_ratios) and a
    # dedicated right column for the colorbar (via width_ratios), instead of letting matplotlib's
    # auto-placement steal space from a data axes -- that was the previous run's overlap bug.
    fig = plt.figure(figsize=(15, 18))
    gs = fig.add_gridspec(nrows=5, ncols=4, height_ratios=[0.35, 1, 1, 1, 1],
                          width_ratios=[1, 1, 1, 0.06], hspace=0.35, wspace=0.15)
    legend_ax = fig.add_subplot(gs[0, :])
    legend_ax.axis("off")
    handles = [plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=vc.ANCESTRY_COLORS[a],
                          markersize=7, label=a) for a in vc.ANCESTRY_ORDER]
    legend_ax.legend(handles=handles, loc="center", ncol=6, frameon=False, fontsize=9)
    legend_ax.set_title("Locked-mapping embedding comparison: 3 allele-space views x "
                        "raw/PCA-denoised UMAP x ancestry/%novel coloring", fontsize=12, pad=14)

    grad_scatter = None
    for ri, row in enumerate(rows):
        for ci, col in enumerate(cols):
            ax = fig.add_subplot(gs[ri + 1, ci])
            coords = coords_by_cell.get((col, row))
            index = index_by_col[col]
            color_kind = "gradient" if row.endswith("gradient") else "ancestry"
            title = f"{col_titles[col]}\n{row_titles[row]}" if ri == 0 else row_titles[row]
            sc = plot_panel(ax, coords, index, color_kind, ancestry_by_person=ancestry_by_person,
                             novel_frac=novel_frac_by_col[col], title=title)
            if sc is not None:
                grad_scatter = sc
            if coords is not None:
                fig_i, ax_i = plt.subplots(figsize=(6, 5.5))
                plot_panel(ax_i, coords, index, color_kind, ancestry_by_person=ancestry_by_person,
                           novel_frac=novel_frac_by_col[col],
                           title=f"{col_titles[col]} -- {row_titles[row]}")
                fig_i.tight_layout()
                vc.savefig(fig_i, os.path.join(matrix_dir, f"panel_{col}_{row}{suffix}.png"))

    if grad_scatter is not None:
        cbar_ax = fig.add_subplot(gs[1:3, 3])
        fig.colorbar(grad_scatter, cax=cbar_ax, label="fraction of 16 slots novel")

    grid_path = os.path.join(matrix_dir, f"embedding_matrix{suffix}.png")
    fig.savefig(grid_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {grid_path}", file=sys.stderr)
    return grid_path


def build_sanity_check_figure(full_raw, full_pca, coll_raw_indep, coll_pca_indep,
                              known_raw_indep, known_pca_indep, index_full, index_coll,
                              index_known, ancestry_by_person, out_dir, suffix):
    """Independent (not projected) fits for collapsed_nearest_ref and known_only, plotted next to
    the already-fit full_enriched -- compare against matrix/panel_*_raw_ancestry.png (the LOCKED
    projections of the same two views) to judge whether the locked mapping is itself flattening
    real differences, or the underlying signal genuinely doesn't move much."""
    sc_dir = os.path.join(out_dir, "sanity_check")
    os.makedirs(sc_dir, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    cells = [
        (axes[0, 0], full_raw, index_full, "full_enriched\n(fit, raw UMAP)"),
        (axes[0, 1], coll_raw_indep, index_coll, "collapsed_nearest_ref\n(INDEPENDENT fit, raw UMAP)"),
        (axes[0, 2], known_raw_indep, index_known, "known_only\n(INDEPENDENT fit, raw UMAP)"),
        (axes[1, 0], full_pca, index_full, "full_enriched\n(fit, PCA-denoised UMAP)"),
        (axes[1, 1], coll_pca_indep, index_coll, "collapsed_nearest_ref\n(INDEPENDENT fit, PCA-denoised)"),
        (axes[1, 2], known_pca_indep, index_known, "known_only\n(INDEPENDENT fit, PCA-denoised)"),
    ]
    for ax, coords, index, title in cells:
        plot_panel(ax, coords, index, "ancestry", ancestry_by_person=ancestry_by_person, title=title)
    handles = [plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=vc.ANCESTRY_COLORS[a],
                          markersize=7, label=a) for a in vc.ANCESTRY_ORDER]
    fig.legend(handles=handles, loc="upper center", ncol=6, bbox_to_anchor=(0.5, 1.06), frameon=False)
    fig.suptitle("Sanity check: independently-fit collapsed_nearest_ref/known_only vs the "
                "locked-projection versions in matrix/ -- do they still look similar?", y=1.12,
                fontsize=11)
    path = os.path.join(sc_dir, f"sanity_check{suffix}.png")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}", file=sys.stderr)
    return path


def build_disease_figure(coords_full_raw, coords_full_pca, index_full, disease_labels, counts,
                         n_multiple, out_dir, suffix, label_order=None, source_note="HLA-disease-"
                         "allele carrier status", filename="disease_alleles"):
    """label_order: ordered list of label strings (defaults to DISEASE_ALLELE_GROUPS' own labels
    for the allele-proxy coloring). Pass the real HLA_LINKED label order + a note saying so when
    coloring by actual EHR diagnosis instead -- see 12_disease_phenotypes.py."""
    label_order = label_order if label_order is not None else \
        [lbl for lbl, *_j in DISEASE_ALLELE_GROUPS]
    dis_dir = os.path.join(out_dir, "disease")
    os.makedirs(dis_dir, exist_ok=True)
    present = [lbl for lbl in label_order if counts.get(lbl, 0) > 0]
    fig, axes = plt.subplots(1, 2, figsize=(15, 7))
    for ax, coords, sub in [(axes[0], coords_full_raw, "raw UMAP"),
                            (axes[1], coords_full_pca, "PCA-denoised UMAP")]:
        plot_panel(ax, coords, index_full, "disease", disease_labels=disease_labels,
                  label_order=label_order,
                  title=f"full_enriched, {sub}, colored by {source_note}")
    color_map = {lbl: DISEASE_COLORS[i % len(DISEASE_COLORS)] for i, lbl in enumerate(label_order)}
    handles = [plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=NONE_COLOR,
                          markersize=7, label=f"none of these ({len(index_full) - sum(counts.values())})")]
    for lbl in present:
        handles.append(plt.Line2D([0], [0], marker="o", color="none",
                                  markerfacecolor=color_map[lbl], markersize=8,
                                  label=f"{lbl} (n={counts[lbl]})"))
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.25),
              frameon=False, fontsize=8)
    fig.suptitle(f"Colored by {source_note} ({n_multiple} people matched >1 group, shown by "
                "their first match)", y=1.03, fontsize=11)
    path = os.path.join(dis_dir, f"{filename}{suffix}.png")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}", file=sys.stderr)
    return path, present


# ---------------------------------------------------------------------------
# Collapse-radius diagnostics (LR-collapsed vs SR-collapsed, same people, same nominal resolution)
# -- unrelated question, kept from the previous version.
# ---------------------------------------------------------------------------
def build_2field_matrix(table1, person_ids, genes):
    sub = table1[table1["person_id"].isin(person_ids) & table1["gene_bare"].isin(genes)]
    per_person = {}
    for pid, gene, hap, consensus in zip(sub["person_id"], sub["gene_bare"], sub["hap"],
                                          sub["consensus"]):
        a2 = vc.to_nfield(consensus, 2)
        if a2 is None:
            continue
        per_person.setdefault(pid, {}).setdefault(gene, {})[hap] = a2
    complete, incomplete = [], 0
    rows = {}
    for pid, gene_map in per_person.items():
        if len(gene_map) < len(genes):
            incomplete += 1
            continue
        vals, ok = [], True
        for gene in genes:
            copies = gene_map.get(gene, {})
            if len(copies) < 2:
                ok = False
                break
            vals.extend((gene, a) for a in copies.values())
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


def concentration_stats(mat, genes):
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
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT)
    ap.add_argument("--table1", default=None)
    ap.add_argument("--cohort-membership", default=None)
    ap.add_argument("--sr-genotypes", default=None)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--umap-neighbors", type=int, default=15)
    ap.add_argument("--umap-min-dist", type=float, default=0.1)
    ap.add_argument("--pca-denoise-k", type=int, default=15)
    ap.add_argument("--sample", action="store_true")
    ap.add_argument("--sanity-check", action="store_true",
                    help="Also fit collapsed_nearest_ref and known_only INDEPENDENTLY (not "
                         "projected) and compare -- see module docstring.")
    ap.add_argument("--disease-alleles", action="store_true",
                    help="Also color full_enriched by curated HLA-disease-allele carrier status "
                         "(a proxy -- close to tautological on an allele-identity embedding).")
    ap.add_argument("--real-disease-labels", default=None,
                    help="Path to 12_disease_phenotypes.py's person_disease_labels.tsv (real "
                         "EHR-confirmed diagnosis, not an allele proxy) -- VM-local, per-person, "
                         "never commit this file's contents. Implies --disease-alleles's figure "
                         "output but colors by actual diagnosis instead.")
    ap.add_argument("--from-cache", action="store_true",
                    help="Skip loading Table 1 / fitting entirely; rebuild figures from "
                         "<out-dir>/embedding_cache.pkl (written by a prior non-cached run).")
    args = ap.parse_args()

    out_dir = args.out_dir or vc.default_out_dir("08_embedding_compare", "lr")
    suffix = ".sample" if args.sample else ""
    cache_path = os.path.join(out_dir, "embedding_cache.pkl")

    if args.from_cache:
        print(f"Loading cache from {cache_path!r} ...", file=sys.stderr)
        with open(cache_path, "rb") as f:
            cache = pickle.load(f)
        mat_full = pd.DataFrame(cache["mat_full_values"], index=cache["mat_full_index"],
                                columns=cache["mat_full_columns"])
        mat_collapsed = pd.DataFrame(cache["mat_collapsed_values"], index=cache["mat_collapsed_index"],
                                     columns=cache["mat_collapsed_columns"])
        coords_by_cell = cache["coords_by_cell"]
        index_by_col = cache["index_by_col"]
        novel_frac_by_col = cache["novel_frac_by_col"]
        ancestry_by_person = cache["ancestry_by_person"]
        known_only_ids = cache["known_only_ids"]
        sanity = cache.get("sanity")
        report = cache["report_so_far"]
    else:
        table1_path = args.table1 or os.path.join(DEFAULT_DATA_ROOT, "hla_calls_rich.tsv")
        cohort_path = args.cohort_membership or os.path.join(DEFAULT_DATA_ROOT,
                                                              "cohort_membership.tsv")
        print(f"Loading Table 1 from {table1_path!r} ...", file=sys.stderr)
        table1 = vc.load_table1(table1_path)
        print(f"Loading Table 4 (ancestry) from {cohort_path!r} ...", file=sys.stderr)
        cohort_df = vc.load_cohort_membership(cohort_path)
        ancestry_by_person = dict(zip(cohort_df["person_id"], cohort_df["ancestry_pred"]))
        lr_people = set(cohort_df.loc[cohort_df["in_lr"], "person_id"]) if "in_lr" in \
            cohort_df.columns else set(table1["person_id"].unique())

        print("Re-deriving novel-allele haplotype matches (03_novel_alleles.py's logic) ...",
              file=sys.stderr)
        novel_id_lookup, match_stats, novel_mod = load_novel_matches(table1, args.outroot)
        print(f"  {match_stats}", file=sys.stderr)

        report = ["# Embedding comparison: locked-mapping allele-space views "
                  "(`08_embedding_compare.py`)\n",
                  "3 columns (full_enriched, collapsed_nearest_ref, known_only) x 4 rows (raw/"
                  "PCA-denoised UMAP x ancestry/%novel-gradient coloring). Full 12-panel grid + "
                  "individual panels: `matrix/`.\n"]

        print("Building full_enriched + collapsed_nearest_ref slot data ...", file=sys.stderr)
        per_person, novel_carrier, novel_slot_count = build_full_and_collapsed(
            table1, lr_people, GENES, novel_id_lookup, novel_mod)
        mat_full, collapsed_idents, complete, n_incomplete = matrices_from_slots(
            per_person, GENES, novel_carrier)
        novel_frac = {pid: novel_slot_count.get(pid, 0) / N_SLOTS for pid in complete}
        report.append(f"\n## Cohort\n{len(complete)} complete 8-classical-gene cases "
                      f"({n_incomplete} incomplete, dropped). "
                      f"{sum(novel_carrier.get(p, False) for p in complete)} "
                      f"({100*sum(novel_carrier.get(p, False) for p in complete)/max(len(complete),1):.1f}%) "
                      f"carry >=1 novel allele. Dosage matrix: {mat_full.shape[0]} people x "
                      f"{mat_full.shape[1]} allele-identity columns.\n")

        mat_collapsed, n_dropped_slots, n_total_slots = build_collapsed_matrix(
            mat_full.columns, collapsed_idents, complete)
        report.append(f"`collapsed_nearest_ref` slot resolution: {n_total_slots - n_dropped_slots} "
                      f"of {n_total_slots} slots ({100*n_dropped_slots/max(n_total_slots,1):.1f}% "
                      f"dropped) had a `nearest_allele` that IS a column in the locked fit.\n")

        known_only_ids = [pid for pid in complete if not novel_carrier.get(pid, False)]
        report.append(f"`known_only`: {len(known_only_ids)} of {len(complete)} people "
                      f"(whole-person drop).\n")

        print(f"Fitting locked mapping (PCA + raw UMAP + PCA-denoised UMAP, "
              f"k={args.pca_denoise_k}) ...", file=sys.stderr)
        locked = LockedEmbedding(mat_full, args.seed, args.umap_neighbors, args.umap_min_dist,
                                  args.pca_denoise_k)
        report.append(f"\nPCA on `full_enriched` ({mat_full.shape[0]} x {mat_full.shape[1]}, "
                      f"standardized): PC1 {100*locked.pca_var[0]:.3f}% var, PC2 "
                      f"{100*locked.pca_var[1]:.3f}% var, cumulative top-{args.pca_denoise_k} "
                      f"{100*locked.pca_var[:args.pca_denoise_k].sum():.2f}% var. See module "
                      f"docstring for why this is expected to be small at this encoding's "
                      f"dimensionality, not evidence of weak signal.\n")

        print("Projecting collapsed_nearest_ref (raw + PCA-denoised) ...", file=sys.stderr)
        coords_by_cell = {
            ("full_enriched", "raw_ancestry"): locked.umap_raw_full,
            ("full_enriched", "raw_gradient"): locked.umap_raw_full,
            ("full_enriched", "pca_ancestry"): locked.umap_pca_full,
            ("full_enriched", "pca_gradient"): locked.umap_pca_full,
        }
        coords_collapsed_raw = locked.transform_raw(mat_collapsed)
        coords_collapsed_pca = locked.transform_pca_denoised(mat_collapsed)
        for row in ["raw_ancestry", "raw_gradient"]:
            coords_by_cell[("collapsed_nearest_ref", row)] = coords_collapsed_raw
        for row in ["pca_ancestry", "pca_gradient"]:
            coords_by_cell[("collapsed_nearest_ref", row)] = coords_collapsed_pca

        known_mask = np.array([pid in set(known_only_ids) for pid in mat_full.index])
        for row, full_coords in [("raw_ancestry", locked.umap_raw_full),
                                 ("raw_gradient", locked.umap_raw_full),
                                 ("pca_ancestry", locked.umap_pca_full),
                                 ("pca_gradient", locked.umap_pca_full)]:
            coords_by_cell[("known_only", row)] = full_coords[known_mask] \
                if full_coords is not None else None

        index_by_col = {"full_enriched": mat_full.index,
                        "collapsed_nearest_ref": mat_collapsed.index,
                        "known_only": pd.Index(known_only_ids)}
        novel_frac_by_col = {"full_enriched": novel_frac, "collapsed_nearest_ref": novel_frac,
                            "known_only": {pid: 0.0 for pid in known_only_ids}}

        sanity = None
        if args.sanity_check:
            print("Fitting collapsed_nearest_ref INDEPENDENTLY (sanity check) ...", file=sys.stderr)
            X_coll = standardize_with(mat_collapsed.values, mat_collapsed.values.mean(axis=0),
                                      mat_collapsed.values.std(axis=0))
            pca_coll, _v, _mu, _vt = fit_pca_with_loadings(X_coll, n_components=args.pca_denoise_k)
            coll_raw_indep, coll_pca_indep = fit_umap_pair(
                X_coll, pca_coll, args.pca_denoise_k, args.seed, args.umap_neighbors,
                args.umap_min_dist)

            print("Fitting known_only INDEPENDENTLY (sanity check) ...", file=sys.stderr)
            mat_known = mat_full.loc[known_only_ids]
            X_known = standardize_with(mat_known.values, mat_known.values.mean(axis=0),
                                       mat_known.values.std(axis=0))
            pca_known, _v, _mu, _vt = fit_pca_with_loadings(X_known, n_components=args.pca_denoise_k)
            known_raw_indep, known_pca_indep = fit_umap_pair(
                X_known, pca_known, args.pca_denoise_k, args.seed, args.umap_neighbors,
                args.umap_min_dist)
            sanity = {"coll_raw_indep": coll_raw_indep, "coll_pca_indep": coll_pca_indep,
                     "known_raw_indep": known_raw_indep, "known_pca_indep": known_pca_indep}

        # Cache everything expensive before doing anything else.
        os.makedirs(out_dir, exist_ok=True)
        with open(cache_path, "wb") as f:
            pickle.dump({
                "mat_full_values": mat_full.values.astype(np.int8),
                "mat_full_index": list(mat_full.index), "mat_full_columns": list(mat_full.columns),
                "mat_collapsed_values": mat_collapsed.values.astype(np.int8),
                "mat_collapsed_index": list(mat_collapsed.index),
                "mat_collapsed_columns": list(mat_collapsed.columns),
                "coords_by_cell": coords_by_cell, "index_by_col": index_by_col,
                "novel_frac_by_col": novel_frac_by_col, "ancestry_by_person": ancestry_by_person,
                "known_only_ids": known_only_ids, "sanity": sanity, "report_so_far": report,
            }, f)
        print(f"Cached fit results to {cache_path!r} (use --from-cache to replot without "
              f"re-fitting).", file=sys.stderr)

    grid_path = build_matrix_figure(coords_by_cell, index_by_col, ancestry_by_person,
                                     novel_frac_by_col, out_dir, suffix)
    report = report + [f"\nFigure (combined grid): `{grid_path}`\n"]

    if args.sanity_check and sanity is not None:
        sc_path = build_sanity_check_figure(
            coords_by_cell[("full_enriched", "raw_ancestry")],
            coords_by_cell[("full_enriched", "pca_ancestry")],
            sanity["coll_raw_indep"], sanity["coll_pca_indep"],
            sanity["known_raw_indep"], sanity["known_pca_indep"],
            index_by_col["full_enriched"], index_by_col["collapsed_nearest_ref"],
            index_by_col["known_only"], ancestry_by_person, out_dir, suffix)
        report.append(f"\n## Sanity check: independent fits vs locked projections\n"
                      f"Figure: `{sc_path}`\nCompare visually against "
                      f"`matrix/panel_collapsed_nearest_ref_raw_ancestry{suffix}.png` and "
                      f"`matrix/panel_known_only_raw_ancestry{suffix}.png` -- if these independent "
                      f"fits also resemble full_enriched, the locked-mapping similarity is real, "
                      f"not a projection artifact.\n")

    if args.real_disease_labels:
        print(f"Loading REAL diagnosis labels from {args.real_disease_labels!r} ...",
              file=sys.stderr)
        disease_labels, counts, n_multiple, label_order = load_real_disease_labels(
            args.real_disease_labels)
        dis_path, present = build_disease_figure(
            coords_by_cell[("full_enriched", "raw_ancestry")],
            coords_by_cell[("full_enriched", "pca_ancestry")],
            index_by_col["full_enriched"], disease_labels, counts, n_multiple, out_dir, suffix,
            label_order=label_order, source_note="ACTUAL EHR-confirmed diagnosis (not an allele "
            "proxy)", filename="disease_real_diagnosis")
        report.append(f"\n## Real EHR-confirmed disease diagnosis coloring\n"
                      f"Actual diagnosis, not allele-carrier proxy -- see "
                      f"`scripts/hla_popgen/12_disease_phenotypes.py`'s `HLA_LINKED` list "
                      f"(reused verbatim from `origin/aleix/hla-resolve-phase1`'s "
                      f"`deep_immune_breakdown.py`, ICD-10-3char + SNOMED-substring definitions). "
                      f"{sum(counts.values())} of {len(index_by_col['full_enriched'])} people "
                      f"matched >=1 disease definition ({n_multiple} matched more than one, shown "
                      f"by first match). Figure: `{dis_path}`\n\n| disease | n people |\n|---|---|\n")
        for lbl in label_order:
            report.append(f"| {lbl} | {counts.get(lbl, 0)} |")
    elif args.disease_alleles:
        # Rebuild from the cache's raw arrays directly (works whether or not we just fit this run
        # -- the cache was written above either way -- disease_labels_from_matrix only needs a
        # DataFrame with the right columns).
        with open(cache_path, "rb") as f:
            c = pickle.load(f)
        mat_full_for_disease = pd.DataFrame(c["mat_full_values"], index=c["mat_full_index"],
                                            columns=c["mat_full_columns"])
        disease_labels, counts, n_multiple = disease_labels_from_matrix(mat_full_for_disease)
        dis_path, present = build_disease_figure(
            coords_by_cell[("full_enriched", "raw_ancestry")],
            coords_by_cell[("full_enriched", "pca_ancestry")],
            index_by_col["full_enriched"], disease_labels, counts, n_multiple, out_dir, suffix)
        report.append(f"\n## HLA-disease-allele carrier coloring (PROXY -- see "
                      f"--real-disease-labels for actual diagnosis)\n"
                      f"Curated, exploratory list ({len(DISEASE_ALLELE_GROUPS)} groups) -- see "
                      f"`DISEASE_ALLELE_GROUPS` in the script for caveats. "
                      f"{sum(counts.values())} people matched at least one group "
                      f"({n_multiple} matched more than one, shown by first match). "
                      f"Figure: `{dis_path}`\n\n| disease/allele group | n carriers |\n|---|---|\n")
        for lbl, _gene, _grp in DISEASE_ALLELE_GROUPS:
            report.append(f"| {lbl} | {counts.get(lbl, 0)} |")

    # --- Collapse-radius: LR-2field vs SR-2field, same people, matched resolution -------------
    # Skipped entirely under --from-cache: it needs Table 1 fresh (no caching for this diagnostic,
    # it's cheap relative to the UMAP fits) and --from-cache runs are specifically for iterating on
    # figures without touching Table 1 again.
    if not args.from_cache:
        sr_people = set(cohort_df.loc[cohort_df["in_sr"], "person_id"]) if "in_sr" in \
            cohort_df.columns else set()
        print("Building 2-field LR + SR matrices for the collapse-radius diagnostic ...",
              file=sys.stderr)
        mat_lr2, _ = build_2field_matrix(table1, set(index_by_col["full_enriched"]), GENES)
        shared_people = set(mat_lr2.index) & sr_people
        cr_dir = os.path.join(out_dir, "collapse_radius")
        os.makedirs(cr_dir, exist_ok=True)
        if shared_people:
            sr_long = vc.load_sr_genotypes(args.sr_genotypes or os.path.join(DEFAULT_DATA_ROOT,
                                                                              "hla_genotypes.tsv"))
            mat_sr2, n_incomplete_sr2 = build_sr_dosage_matrix(sr_long, shared_people, GENES)
            shared_both = set(mat_lr2.index) & set(mat_sr2.index)
            stats_lr = concentration_stats(mat_lr2.loc[mat_lr2.index.isin(shared_both)],
                                           GENES).set_index("gene")
            stats_sr = concentration_stats(mat_sr2.loc[mat_sr2.index.isin(shared_both)],
                                           GENES).set_index("gene")
            report.append(
                f"\n## Collapse-radius: LR-2field vs SR-2field, SAME {len(shared_both)} people\n"
                "`top10_share` = fraction of allele copies in the 10 most common alleles per "
                "gene; `simpson_concentration` = sum(p_i^2).\n")
            report.append("| gene | LR top10 share | SR top10 share | LR Simpson | SR Simpson | "
                           "LR n_distinct | SR n_distinct |")
            report.append("|---|---|---|---|---|---|---|")
            for gene in GENES:
                if gene not in stats_lr.index or gene not in stats_sr.index:
                    continue
                lr_r, sr_r = stats_lr.loc[gene], stats_sr.loc[gene]
                report.append(f"| {gene} | {100*lr_r['top10_share']:.1f}% | "
                               f"{100*sr_r['top10_share']:.1f}% | "
                               f"{lr_r['simpson_concentration']:.3f} | "
                               f"{sr_r['simpson_concentration']:.3f} | "
                               f"{int(lr_r['n_distinct_alleles'])} | "
                               f"{int(sr_r['n_distinct_alleles'])} |")
        else:
            report.append("\n## Collapse-radius\nNo people shared -- skipped.\n")

    md_path = os.path.join(out_dir, f"08_embedding_compare_report{suffix}.md")
    vc.write_report(md_path, report)
    print(f"Wrote report to {md_path!r}.", file=sys.stderr)


if __name__ == "__main__":
    main()
