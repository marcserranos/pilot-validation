#!/usr/bin/env python3
"""Does the enriched (known+novel) allele space look different once you collapse it back down --
and does it look different in the SAME coordinate frame, not just a similarly-shaped one fit
separately? This is the locked-mapping design from the 2026-09 chat discussion, replacing an
earlier version of this script that fit three independent embeddings and only argued they were
comparable (see git history if you need that version back).

## The three columns (same people's coordinates, one fitted mapping, two mapped into it)

1. **`full_enriched`** -- LR cohort, 8 classical genes, complete diploid cases, each allele call at
   its own native resolution (full-precision `consensus` string for known calls, `novel_id` for
   novel calls). This is the ONLY column that gets FIT (PCA and, per-variant, UMAP) -- the mapping
   used everywhere else is locked from this fit and never re-fit.
2. **`collapsed_nearest_ref`** -- ALL of the same people, but every novel call is replaced by its
   `nearest_allele` (Immuannot's own tied-best CDS-search candidate, already computed per novel
   call -- no new distance metric invented here). Standardized with `full_enriched`'s own mean/std
   (not refit) and placed into the locked mapping via `.transform()` (UMAP) or the retained PCA
   loadings (PCA) -- a genuine projection, not a second fit, which is what makes this column
   comparable to column 1 in the same axes. A `nearest_allele` that isn't itself a column anyone in
   the cohort carries (rare) has nowhere to project and is dropped for that one slot; counted, not
   guessed.
3. **`known_only`** -- people with ZERO novel calls across these 8 genes (dropped as a whole
   person, not per-slot, per 2026-09 chat clarification). Needs no projection at all -- these
   people are already inside the column-1 fit; this column is that fit's coordinates, filtered.

## The four rows -- coloring x whether UMAP ran on a PCA-denoised representation

PCA-then-UMAP (denoise to `--pca-denoise-k` components before fitting/transforming UMAP) is a
PARALLEL EXPERIMENT here, not a replacement default -- both `raw` (UMAP direct on the standardized
dosage matrix) and `pca_denoised` are run and shown side by side, exactly so that question can be
answered by looking rather than by picking one. Each variant is shown twice: colored by ancestry
(categorical, Okabe-Ito), and colored by a CONTINUOUS gradient (viridis) on each person's fraction
of novel slots out of 16 (8 genes x 2 haps) -- not a binary carrier flag. Column 3 by construction
has 0 novel slots for everyone, so its gradient panel is a uniform-color sanity check that the
`known_only` filter is doing what it says, not a substantive result.

## Output layout (2026-09 debloat pass)

All 12 individual panels plus the combined 3x4 grid live under `<out-dir>/matrix/` -- nothing
scattered at the top level. The SR-vs-LR collapse-radius diagnostic (a genuinely different
question: does SR concentrate on fewer, more common alleles than LR at matched 2-field resolution)
is kept, under `<out-dir>/collapse_radius/`.

## Encoding

Same dosage convention as 06_figures_structure.py: one column per (gene, allele-identity) pair,
value = copy count (0/1/2) per person, unordered hap1/hap2, complete-case only for the fitted
column -- incomplete cases are dropped and counted, never imputed (SCHEMA.md discipline).

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
N_SLOTS = len(GENES) * 2  # 8 genes x 2 haps
DEFAULT_DATA_ROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_OUTROOT = os.path.join(DEFAULT_DATA_ROOT, "people")


# ---------------------------------------------------------------------------
# Identity assignment
# ---------------------------------------------------------------------------
def load_novel_matches(table1, outroot):
    """Re-derives the (person_id, hap, gene) -> novel_id map by importing 03_novel_alleles.py as a
    sibling module (never re-implement the matching logic twice); also hands back the module itself
    so `nearest_allele_of()` can be reused verbatim for the collapsed_nearest_ref column."""
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
    collapsed_nearest_ref fallback identity, so the two matrices are guaranteed to share a person
    set and slot structure. Returns:
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
            nearest_ident = full_ident  # known calls: collapsed == full, nothing to fall back to
        if full_ident is None:
            continue  # unresolved/ambiguous -- never guess (SCHEMA.md)
        slots = per_person.setdefault(pid, {}).setdefault(gene, {})
        slots[len(slots) + 1] = (full_ident, nearest_ident)
    return per_person, novel_carrier, novel_slot_count


def matrices_from_slots(per_person, genes, novel_carrier):
    """Turns the per-person slot dict into (mat_full, mat_collapsed_raw_idents, complete_person_ids,
    n_incomplete). mat_collapsed_raw_idents is a DataFrame of the SAME shape as mat_full but built
    from each slot's fallback identity string -- restricting it to mat_full's own columns (the only
    columns the locked fit knows about) happens later, once mat_full's column universe is fixed."""
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
    """Builds the collapsed_nearest_ref matrix IN mat_full's exact column space -- a nearest_allele
    that isn't one of those columns has nowhere to project and is dropped for that slot (counted,
    not guessed). This restriction is what makes `.transform()`/PCA-loading projection valid."""
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


# ---------------------------------------------------------------------------
# Locked-mapping fit + transform (PCA with retained loadings; UMAP with retained reducer)
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


class LockedEmbedding:
    """One PCA fit (with retained loadings) + optionally one UMAP fit (raw, on the standardized
    dosage matrix) + optionally one UMAP fit on the PCA-denoised representation -- everything a
    `collapsed_nearest_ref` matrix needs to be projected into the SAME coordinate frame as
    `full_enriched`, for both the raw and pca_denoised rows of the comparison matrix."""
    def __init__(self, mat_full, seed, umap_neighbors, umap_min_dist, pca_denoise_k):
        self.mu = mat_full.values.mean(axis=0)
        self.sd = mat_full.values.std(axis=0)
        X = standardize_with(mat_full.values, self.mu, self.sd)
        self.pca_scores, self.pca_var, self.pca_mu, self.pca_Vt = fit_pca_with_loadings(
            X, n_components=max(pca_denoise_k, 10))
        self.pca_denoise_k = pca_denoise_k
        self.umap_raw_reducer = None
        self.umap_raw_full = None
        self.umap_pca_reducer = None
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
        if self.umap_raw_reducer is None:
            return None
        return self.umap_raw_reducer.transform(X_new)

    def transform_pca_denoised(self, mat_new):
        X_new = standardize_with(mat_new.values, self.mu, self.sd)
        pca_new = transform_pca(X_new, self.pca_mu, self.pca_Vt)[:, :self.pca_denoise_k]
        if self.umap_pca_reducer is None:
            return None
        return self.umap_pca_reducer.transform(pca_new)


# ---------------------------------------------------------------------------
# Plotting: one panel, one coloring
# ---------------------------------------------------------------------------
def plot_panel(ax, coords, index, ancestry_by_person, novel_frac, color_mode, title):
    if coords is None:
        ax.text(0.5, 0.5, "umap-learn not installed", ha="center", va="center", fontsize=9,
                color="#999999", transform=ax.transAxes)
        ax.set_title(title, fontsize=9)
        ax.axis("off")
        return
    sc = None
    if color_mode == "ancestry":
        anc_labels = [pid_anc if isinstance(pid_anc, str) else None
                      for pid_anc in (ancestry_by_person.get(pid) for pid in index)]
        for anc in vc.ANCESTRY_ORDER:
            mask = np.array([l == anc for l in anc_labels])
            if not mask.any():
                continue
            ax.scatter(coords[mask, 0], coords[mask, 1], s=6, alpha=0.5,
                       color=vc.ANCESTRY_COLORS[anc], edgecolors="none")
    else:  # gradient: fraction of novel slots, viridis
        vals = np.array([novel_frac.get(pid, 0.0) for pid in index])
        sc = ax.scatter(coords[:, 0], coords[:, 1], s=6, alpha=0.6, c=vals, cmap="viridis",
                        vmin=0, vmax=1, edgecolors="none")
    ax.set_title(title, fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[["top", "right"]].set_visible(False)
    return sc  # None for ancestry panels; the scatter handle for gradient panels (shared colorbar)


def build_matrix_figure(coords_by_cell, index_by_col, ancestry_by_person, novel_frac_by_col,
                         out_dir, suffix):
    """coords_by_cell: {(col, row): coords_or_None} for the 3x4 grid.
    col in ['full_enriched','collapsed_nearest_ref','known_only'],
    row in ['raw_ancestry','raw_gradient','pca_ancestry','pca_gradient']."""
    cols = ["full_enriched", "collapsed_nearest_ref", "known_only"]
    col_titles = {"full_enriched": "full_enriched (fit)",
                  "collapsed_nearest_ref": "collapsed_nearest_ref (projected)",
                  "known_only": "known_only (subset of fit)"}
    rows = ["raw_ancestry", "raw_gradient", "pca_ancestry", "pca_gradient"]
    row_titles = {"raw_ancestry": "raw UMAP, ancestry", "raw_gradient": "raw UMAP, %novel",
                  "pca_ancestry": "PCA-denoised UMAP, ancestry",
                  "pca_gradient": "PCA-denoised UMAP, %novel"}

    fig, axes = plt.subplots(4, 3, figsize=(15, 19))
    grad_scatter = None
    os.makedirs(os.path.join(out_dir, "matrix"), exist_ok=True)
    for ri, row in enumerate(rows):
        for ci, col in enumerate(cols):
            ax = axes[ri, ci]
            coords = coords_by_cell.get((col, row))
            index = index_by_col[col]
            color_mode = "gradient" if row.endswith("gradient") else "ancestry"
            title = f"{col_titles[col]}\n{row_titles[row]}" if ri == 0 else row_titles[row]
            sc = plot_panel(ax, coords, index, ancestry_by_person,
                             novel_frac_by_col[col], color_mode, title)
            if sc is not None:
                grad_scatter = sc

            # Also save each panel standalone into matrix/ (the "put it all in a folder" ask).
            if coords is not None:
                fig_i, ax_i = plt.subplots(figsize=(6, 5.5))
                plot_panel(ax_i, coords, index, ancestry_by_person, novel_frac_by_col[col],
                           color_mode, f"{col_titles[col]} -- {row_titles[row]}")
                fig_i.tight_layout()
                vc.savefig(fig_i, os.path.join(out_dir, "matrix", f"panel_{col}_{row}{suffix}.png"))

    if grad_scatter is not None:
        fig.colorbar(grad_scatter, ax=axes[[1, 3], :].ravel().tolist(), shrink=0.5,
                     label="fraction of 16 slots novel")
    # Ancestry legend, once, off to the side.
    handles = [plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=vc.ANCESTRY_COLORS[a],
                          markersize=7, label=a) for a in vc.ANCESTRY_ORDER]
    fig.legend(handles=handles, loc="upper center", ncol=6, bbox_to_anchor=(0.5, 0.995),
              frameon=False, fontsize=9)
    fig.suptitle("Locked-mapping embedding comparison: 3 allele-space views x "
                 "raw/PCA-denoised UMAP x ancestry/%novel coloring", y=1.015, fontsize=12)
    grid_path = os.path.join(out_dir, "matrix", f"embedding_matrix{suffix}.png")
    fig.savefig(grid_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {grid_path}", file=sys.stderr)
    return grid_path


# ---------------------------------------------------------------------------
# Collapse-radius diagnostics (LR-collapsed vs SR-collapsed, same people, same nominal resolution)
# -- unchanged question from the previous version of this script: still 2-field truncation, still
# separately fit (SR has no fine resolution to lock a mapping from), just re-homed under its own
# subfolder per the 2026-09 debloat request.
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
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT,
                    help="Root holding <person_id>/immuannot_output/hap{1,2}/cds.fa.gz. Default: "
                         "~/pipeline_outputs/people (NOT ~/pipeline_outputs itself).")
    ap.add_argument("--table1", default=None, help="Path to hla_calls_rich.tsv.")
    ap.add_argument("--cohort-membership", default=None, help="Path to cohort_membership.tsv.")
    ap.add_argument("--sr-genotypes", default=None, help="Path to hla_genotypes.tsv (AoU-native SR).")
    ap.add_argument("--out-dir", default=None,
                    help="Where to write output. Default: reports/hla_popgen/08_embedding_compare/lr/")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--umap-neighbors", type=int, default=15)
    ap.add_argument("--umap-min-dist", type=float, default=0.1)
    ap.add_argument("--pca-denoise-k", type=int, default=15,
                    help="Number of PCA components to denoise to before the parallel "
                         "PCA-denoised UMAP variant (not the default path -- see module docstring).")
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
    novel_id_lookup, match_stats, novel_mod = load_novel_matches(table1, args.outroot)
    print(f"  {match_stats}", file=sys.stderr)

    report = ["# Embedding comparison: locked-mapping allele-space views "
              "(`08_embedding_compare.py`)\n",
              "3 columns (full_enriched, collapsed_nearest_ref, known_only) x 4 rows (raw/"
              "PCA-denoised UMAP x ancestry/%novel-gradient coloring) -- see module docstring for "
              "why `collapsed_nearest_ref` is a genuine `.transform()` projection into the "
              "`full_enriched` fit, not a second fit. Full 12-panel grid + individual panels: "
              "`matrix/`.\n"]

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
                  f"carry >=1 novel allele.\n")

    mat_collapsed, n_dropped_slots, n_total_slots = build_collapsed_matrix(
        mat_full.columns, collapsed_idents, complete)
    report.append(f"`collapsed_nearest_ref` slot resolution: {n_total_slots - n_dropped_slots} of "
                  f"{n_total_slots} slots ({100*n_dropped_slots/max(n_total_slots,1):.1f}% dropped) "
                  f"had a `nearest_allele` that IS a column in the locked fit and could be placed; "
                  f"the rest had nowhere to project (a nearest reference nobody in this cohort's "
                  f"own known calls happens to carry) and are left at 0 for that slot.\n")

    known_only_ids = [pid for pid in complete if not novel_carrier.get(pid, False)]
    report.append(f"`known_only`: {len(known_only_ids)} of {len(complete)} people "
                  f"(whole-person drop, per 2026-09 chat clarification).\n")

    print(f"Fitting locked mapping (PCA + raw UMAP + PCA-denoised UMAP, "
          f"k={args.pca_denoise_k}) ...", file=sys.stderr)
    locked = LockedEmbedding(mat_full, args.seed, args.umap_neighbors, args.umap_min_dist,
                              args.pca_denoise_k)
    report.append(f"\nPCA on `full_enriched`: PC1 {100*locked.pca_var[0]:.2f}% var, PC2 "
                  f"{100*locked.pca_var[1]:.2f}% var (10 components retained for the record; "
                  f"the PCA-denoised UMAP row uses the top {args.pca_denoise_k}).\n")

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
        coords_by_cell[("known_only", row)] = full_coords[known_mask] if full_coords is not None \
            else None

    index_by_col = {"full_enriched": mat_full.index, "collapsed_nearest_ref": mat_collapsed.index,
                    "known_only": pd.Index(known_only_ids)}
    novel_frac_by_col = {"full_enriched": novel_frac, "collapsed_nearest_ref": novel_frac,
                        "known_only": {pid: 0.0 for pid in known_only_ids}}

    grid_path = build_matrix_figure(coords_by_cell, index_by_col, ancestry_by_person,
                                     novel_frac_by_col, out_dir, suffix)
    report.append(f"\nFigure (combined grid): `{grid_path}`\n")

    # --- Collapse-radius: LR-2field vs SR-2field, same people, matched resolution -------------
    print("Building 2-field LR + SR matrices for the collapse-radius diagnostic ...",
          file=sys.stderr)
    mat_lr2, _ = build_2field_matrix(table1, set(complete), GENES)
    shared_people = set(mat_lr2.index) & sr_people
    cr_dir = os.path.join(out_dir, "collapse_radius")
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
            "Does SR concentrate on fewer, more common alleles than LR at the SAME matched "
            "2-field resolution? `top10_share` = fraction of allele copies in the 10 most common "
            "alleles per gene; `simpson_concentration` = sum(p_i^2).\n")
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
    else:
        report.append("\n## Collapse-radius\nNo people shared between the LR-2field complete-case "
                       "set and the SR cohort membership -- skipped.\n")

    os.makedirs(cr_dir, exist_ok=True)
    md_path = os.path.join(out_dir, f"08_embedding_compare_report{suffix}.md")
    vc.write_report(md_path, report)
    print(f"Wrote report to {md_path!r}.", file=sys.stderr)


if __name__ == "__main__":
    main()
