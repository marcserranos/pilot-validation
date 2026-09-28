#!/usr/bin/env python3
"""01 -- Visualize the SCEPTR TRB repertoire embeddings (report: reports/01_sceptr_embedding_viz/).

Two levels, two figures, plus a PCA supplement:

  Fig 1  CLONOTYPE map. Each dot = one TRB clonotype (unique TRBV + CDR3aa per person), a
         seeded random subsample of the embedded pool. UMAP of the SCEPTR vectors, coloured by
           A  TRBV gene family       -- sanity check; V is a b_sceptr input, so this SHOULD
                                         structure the map. It validates, it doesn't discover.
           B  CDR3 length            -- the dominant known axis of TCR sequence variation
           C  clonal expansion       -- log10 read support; expanded = likely antigen-driven
           D  publicness             -- how many people carry the identical clonotype
  Fig 2  PERSON map. Each dot = one person, repertoire vector = mean of their clonotype
         vectors. UMAP coloured by
           E  genetic ancestry
           F  clonotypes in the repertoire (capped at --max-per-person in embed_cdr3s.py),
              to see whether person-level structure is biology or just recovery depth
  Fig S1 the same colourings on PCA (linear, no tuning -- a check that UMAP isn't
         inventing the structure).

Runs ON THE VM (inputs have real research_ids). Figures + a de-identified summary.csv go to
--outdir (VM-local by default). Figure dots are individual clonotypes/people, so whether they
can be committed to the public repo is the open small-cell/disclosure question in
context/DECISIONS.md -- this script never pushes anything.

Usage (from aleix/RNA-seq/, after embed_cdr3s.py --models sceptr):
  pixi run python3 scripts/01_sceptr_embedding_viz.py [--tag cohort_full_vcdr3]
      [--n-clonotypes 200000] [--seed 0] [--outdir ~/pipeline_outputs/rnaseq/reports/01_sceptr_embedding_viz]

Needs: cnsplots umap-learn (pixi run pip install cnsplots umap-learn).
"""
import argparse
import os
import sys
import time

import matplotlib

matplotlib.use("Agg")

import cnsplots as cns  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402

ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
PUBLIC_BINS = [(1, 1, "1 (private)"), (2, 9, "2–9"), (10, 99, "10–99"), (100, np.inf, "≥100")]
PUBLIC_COLORS = ["#D9D9D9", "#4DBBD5", "#F39B7F", "#B2182B"]  # grey under, warm on top
OTHER_COLOR = "#D9D9D9"
TOP_V_FAMILIES = 10
PT_CLONOTYPE = 0.15
PT_PERSON = 2.0
PANEL = 150


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def v_family(v):
    """TRBV10-3 -> TRBV10; TRBV20-1 -> TRBV20; empty stays empty."""
    return v.split("-")[0] if v else ""


def run_umap(x, seed, cache_path):
    """Seeded (reproducible, so single-threaded). Cached, so restyling a figure doesn't
    recompute it. The cache holds per-clonotype/per-person coordinates: VM-local only."""
    if os.path.exists(cache_path):
        emb = np.load(cache_path)
        if len(emb) == len(x):
            print(f"  UMAP: cached {cache_path}", file=sys.stderr)
            return emb
    import umap
    t0 = time.time()
    emb = umap.UMAP(n_neighbors=15, min_dist=0.3, metric="euclidean",
                    random_state=seed).fit_transform(x)
    print(f"  UMAP on {len(x):,} x {x.shape[1]}: {time.time() - t0:.0f}s", file=sys.stderr)
    np.save(cache_path, emb)
    return emb


def categorical(ax, xy, labels, order, colors, title, legend_title, s, layered=False):
    """layered=False: random draw order, so no category systematically hides another.
    layered=True: draw categories in `order` (first = bottom), for a rare-on-top encoding."""
    labels = np.asarray(labels)
    rng = np.random.default_rng(0)
    idx = rng.permutation(len(xy))
    if layered:
        rank = pd.Series(range(len(order)), index=order)
        idx = idx[np.argsort(rank.reindex(labels[idx]).to_numpy(), kind="stable")]
    df = pd.DataFrame({"x": xy[idx, 0], "y": xy[idx, 1], "g": labels[idx]})
    cns.scatterplot(data=df, x="x", y="y", hue="g", hue_order=order,
                    palette=dict(zip(order, colors)), s=s, linewidth=0, rasterized=True, ax=ax)
    ax.set(title=title, xticks=[], yticks=[])
    cns.take_legend_out(title=legend_title, ax=ax)
    for h in ax.get_legend().legend_handles:  # dots this small make invisible legend keys
        if hasattr(h, "set_sizes"):
            h.set_sizes([14])
        elif hasattr(h, "set_markersize"):
            h.set_markersize(4)


def continuous(ax, xy, values, title, cbar_label, s, vmin=None, vmax=None, cmap="parula"):
    order = np.argsort(values, kind="stable")  # high values drawn last, on top
    sc = ax.scatter(xy[order, 0], xy[order, 1], c=np.asarray(values)[order], s=s,
                    linewidths=0, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
    ax.set(title=title, xticks=[], yticks=[])
    cax = ax.inset_axes([1.04, 0.15, 0.04, 0.7])
    cb = plt.colorbar(sc, cax=cax)
    cb.set_label(cbar_label)
    cb.outline.set_linewidth(0.5)


def save(stem, outdir):
    for ext in ("svg", "pdf", "png"):
        cns.savefig(os.path.join(outdir, f"{stem}.{ext}"))
    plt.close("all")
    print(f"  wrote {stem}.svg/.pdf/.png", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", default="cohort_full_vcdr3",
                    help="embed_cdr3s.py output tag (default matches the full-cohort v+cdr3 run)")
    ap.add_argument("--embeddings-dir", default="~/pipeline_outputs/rnaseq/embeddings")
    ap.add_argument("--n-clonotypes", type=int, default=200_000,
                    help="clonotypes subsampled for the clonotype-level map (default 200k)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", default="~/pipeline_outputs/rnaseq/reports/01_sceptr_embedding_viz")
    args = ap.parse_args()

    edir = os.path.expanduser(args.embeddings_dir)
    emb_path = os.path.join(edir, f"embeddings_sceptr_{args.tag}.npy")
    pool_path = os.path.join(edir, f"pool_sceptr_{args.tag}.tsv")
    for p in (emb_path, pool_path):
        if not os.path.exists(p):
            die(f"missing {p} -- run embed_cdr3s.py --models sceptr first (or check --tag)")
    outdir = os.path.expanduser(args.outdir)
    cache = os.path.join(outdir, "_cache")
    os.makedirs(cache, exist_ok=True)

    embs = np.load(emb_path)
    pool = pd.read_csv(pool_path, sep="\t", dtype={"research_id": str}, keep_default_na=False)
    if len(embs) != len(pool):
        die(f"embeddings ({len(embs):,}) and pool ({len(pool):,}) are not row-aligned")
    pool["trbv"] = pool["v_gene"].astype(str).str.split("*").str[0]
    pool["reads"] = pd.to_numeric(pool["reads"], errors="coerce")
    print(f"Loaded {len(pool):,} clonotypes x {embs.shape[1]} dims, "
          f"{pool['research_id'].nunique():,} people", file=sys.stderr)

    # Publicness over the FULL pool, not the subsample: number of distinct people carrying
    # the identical (TRBV, CDR3aa) clonotype.
    pool["n_people_sharing"] = (pool.groupby(["trbv", "cdr3aa"])["research_id"]
                                    .transform("nunique"))

    # ---------------- Fig 1: clonotype level ----------------
    rng = np.random.default_rng(args.seed)
    n = min(args.n_clonotypes, len(pool))
    sel = np.sort(rng.choice(len(pool), size=n, replace=False))
    sub = pool.iloc[sel].reset_index(drop=True)
    x_sub = embs[sel]
    print(f"\nFig 1: {n:,} clonotypes (seed {args.seed})", file=sys.stderr)
    xy = run_umap(x_sub, args.seed,
                  os.path.join(cache, f"umap_clonotype_{args.tag}_n{n}_s{args.seed}.npy"))
    pca_c = PCA(n_components=10, random_state=args.seed).fit(x_sub)
    pc_c = pca_c.transform(x_sub)[:, :2]

    fam = sub["trbv"].map(v_family)
    top = fam.value_counts().head(TOP_V_FAMILIES).index.tolist()
    fam_lab = np.where(fam.isin(top), fam, "other")
    fam_order = sorted(top, key=lambda f: int("".join(c for c in f[4:] if c.isdigit()) or 0)) + ["other"]
    cdr3_len = sub["cdr3aa"].str.len().to_numpy()
    log_reads = np.log10(sub["reads"].clip(lower=1).to_numpy())
    pub = sub["n_people_sharing"].to_numpy()
    pub_lab = np.select([(pub >= lo) & (pub <= hi) for lo, hi, _ in PUBLIC_BINS],
                        [lab for *_, lab in PUBLIC_BINS], default="")
    pub_order = [lab for *_, lab in PUBLIC_BINS]

    nature = cns.palettes("Nature")
    fam_colors = list(nature[:len(fam_order) - 1]) + [OTHER_COLOR]
    s = PT_CLONOTYPE

    def clonotype_fig(coords, axis_name, stem):
        mp = cns.multipanel(max_width=540)
        categorical(mp.panel("A", width=PANEL, height=PANEL, margin_right=80), coords, fam_lab,
                    fam_order, fam_colors, "TRBV family", "TRBV family", s)
        continuous(mp.panel("B", width=PANEL, height=PANEL, margin_right=60), coords, cdr3_len,
                   "CDR3 length", "CDR3 length (aa)", s, vmin=10, vmax=20)
        continuous(mp.panel("C", width=PANEL, height=PANEL, margin_right=80), coords, log_reads,
                   "Clonal expansion", "log$_{10}$ reads", s,
                   vmax=float(np.quantile(log_reads, 0.99)))
        categorical(mp.panel("D", width=PANEL, height=PANEL, margin_right=60), coords, pub_lab,
                    pub_order, PUBLIC_COLORS, "Publicness", "People sharing", s, layered=True)
        for ax in mp.axes:
            ax.set(xlabel=f"{axis_name} 1", ylabel=f"{axis_name} 2")
        save(stem, outdir)

    clonotype_fig(xy, "UMAP", "fig1_clonotype_umap")

    # ---------------- Fig 2: person level ----------------
    codes, people = pd.factorize(pool["research_id"])
    person_vec = np.zeros((len(people), embs.shape[1]))
    np.add.at(person_vec, codes, embs)
    n_clono = np.bincount(codes)
    person_vec /= n_clono[:, None]
    anc = (pool.drop_duplicates("research_id").set_index("research_id")
               .reindex(people)["ancestry"].fillna("").to_numpy()
           if "ancestry" in pool.columns else np.full(len(people), ""))
    print(f"\nFig 2: {len(people):,} people", file=sys.stderr)
    xy_p = run_umap(person_vec, args.seed,
                    os.path.join(cache, f"umap_person_{args.tag}_s{args.seed}.npy"))
    pca_p = PCA(n_components=10, random_state=args.seed).fit(person_vec)
    pc_p = pca_p.transform(person_vec)[:, :2]
    anc_order = [a for a in ANCESTRY_ORDER if a in set(anc)] + \
                sorted(set(anc) - set(ANCESTRY_ORDER) - {""})

    def person_fig(coords, axis_name, stem):
        mp = cns.multipanel(max_width=540)
        categorical(mp.panel("E", width=PANEL, height=PANEL, margin_right=70), coords, anc,
                    anc_order, list(nature[:len(anc_order)]), "Genetic ancestry", "Ancestry",
                    PT_PERSON)
        continuous(mp.panel("F", width=PANEL, height=PANEL, margin_right=80), coords, n_clono,
                   "Repertoire size", "TRB clonotypes", PT_PERSON)
        for ax in mp.axes:
            ax.set(xlabel=f"{axis_name} 1", ylabel=f"{axis_name} 2")
        save(stem, outdir)

    person_fig(xy_p, "UMAP", "fig2_person_umap")

    # ---------------- Fig S1: PCA versions ----------------
    clonotype_fig(pc_c, "PC", "figS1a_clonotype_pca")
    person_fig(pc_p, "PC", "figS1b_person_pca")

    # ---------------- de-identified summary ----------------
    rho_pc1, p_pc1 = spearmanr(pc_p[:, 0], n_clono)
    rho_pc2, p_pc2 = spearmanr(pc_p[:, 1], n_clono)
    rows = [
        ("clonotypes_embedded", len(pool)),
        ("people", len(people)),
        ("clonotypes_in_map", n),
        ("seed", args.seed),
        ("median_clonotypes_per_person", float(np.median(n_clono))),
        # The 500 cap is applied before V-usability exclusion, so "at cap" is ~490-500, not
        # exactly 500 -- report the low tail instead of an exact-equality count.
        ("frac_people_lt_100_clonotypes", float(np.mean(n_clono < 100))),
        ("frac_people_lt_400_clonotypes", float(np.mean(n_clono < 400))),
        ("frac_clonotypes_public_ge2", float(np.mean(pool["n_people_sharing"] >= 2))),
        ("clonotype_pca_var_PC1", pca_c.explained_variance_ratio_[0]),
        ("clonotype_pca_var_PC2", pca_c.explained_variance_ratio_[1]),
        ("person_pca_var_PC1", pca_p.explained_variance_ratio_[0]),
        ("person_pca_var_PC2", pca_p.explained_variance_ratio_[1]),
        ("spearman_person_PC1_vs_repertoire_size", rho_pc1),
        ("spearman_person_PC1_vs_repertoire_size_p", p_pc1),
        ("spearman_person_PC2_vs_repertoire_size", rho_pc2),
        ("spearman_person_PC2_vs_repertoire_size_p", p_pc2),
    ]
    for a in anc_order:
        rows.append((f"people_{a}", int(np.sum(anc == a))))
    summary = pd.DataFrame(rows, columns=["metric", "value"])
    summary.to_csv(os.path.join(outdir, "summary.csv"), index=False)
    print("\n" + summary.to_string(index=False))
    print(f"\nAll outputs in {outdir}")


if __name__ == "__main__":
    main()
