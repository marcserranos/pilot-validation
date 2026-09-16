#!/usr/bin/env python3
"""How much of the HLA allele space have we seen, and what would it take to see more?
Spec: sprints/S01_novelty_qc_coverage/WS3_allele_space_coverage.md ("script 27").

Reads 24's `allele_counts_by_resolution.tsv` (one row per gene x resolution x
exclude_flagged x ancestry_scheme x ancestry x unrelated_only x allele_id, with a
`n_haplotypes` count) and applies the Good-Turing / Chao & Jost sample-coverage
estimators in `_coverage.py` to answer, per gene x ancestry:

    "A catalogue built from the haplotypes we've already typed would type X% of the
    next haplotype we see. Reaching 99% would need about N more haplotypes."

Coverage (not richness) is the headline statistic -- see `_coverage.py`'s docstring and
WS3's brief for why: this project's allele-frequency spectra have a long singleton tail
that makes richness estimators (Chao1/Chao2) unstable, while coverage is well-behaved
because it only depends on the shape of the spectrum near k=1,2.

USAGE
-----
    pixi run -e spechla -- python3 scripts/hla_popgen/27_allele_space_coverage.py \\
        [--allele-counts PATH] [--out-dir DIR] [--exclude-flagged clean_only|all] \\
        [--ancestry-scheme strict|pred] [--unrelated-only / --no-unrelated-only] \\
        [--limit-genes N]

OUTPUTS (in --out-dir, default ~/results/27_allele_space_coverage)
--------------------------------------------------------------------
    coverage_by_gene_ancestry.tsv   -- per gene x resolution x ancestry: n, S_obs, f1, f2,
                                        coverage (+ bootstrap CI), unseen mass, Chao1
                                        (secondary/unstable), sizes for 95/99/99.9% coverage
    coverage_curves.tsv             -- <=200 points per (gene, resolution, ancestry) line
    coverage_curves_protein.png     -- small multiples per gene, one line per ancestry
    coverage_curves_cds.png
    cross_ancestry_matrix.tsv       -- per-gene + classical-gene-mean cross-coverage
    cross_ancestry_matrix.png
    imgt_coverage_today.tsv         -- fraction of haplotypes already in IPD-IMGT/HLA
    imgt_coverage_today.png
    design_scenarios.tsv            -- fixed-N sampling-design coverage scenarios
    design_scenarios.png
    summary.json
    allele_space_coverage_report.md -- plain-language report

This script is pure aggregation/plotting glue around `_coverage.py`'s estimators -- no
new coverage math lives here. All outputs are aggregate (per gene/ancestry); no
allele-level rows are ever written, and any count in [1, 20) is suppressed as "<20"
(`SUPPRESS_BELOW`, matching 24's convention).

Compatibility: must run under python 3.8 / pandas 2.0.3 / numpy 1.22 as well as newer
stacks -- no `DataFrame.map` (use `.apply`/`.applymap`/dict-based `.map` sparingly, see
note at `_series_map`), no `list[str]`-style annotations, no `match` statements, no
`include_groups=` kwarg to groupby.apply.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _viz_common as vc
import _coverage as cov

SUPPRESS_BELOW = 20
DEFAULT_ALLELE_COUNTS = os.path.expanduser(
    "~/pipeline_outputs/s01_local/24/allele_counts_by_resolution.tsv")
DEFAULT_OUT_DIR = os.path.expanduser("~/results/27_allele_space_coverage")

REQUIRED_COLUMNS = ["gene", "resolution", "exclude_flagged", "ancestry_scheme", "ancestry",
                    "unrelated_only", "allele_id", "n_haplotypes"]

TARGET_COVERAGES = [0.95, 0.99, 0.999]
N_CURVE_POINTS = 200
N_BOOT = 200


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def suppress(n, threshold=SUPPRESS_BELOW):
    """Matches 24_novelty_by_field.suppress: bare small counts become '<threshold'."""
    if pd.isna(n):
        return n
    n = int(n)
    if 1 <= n < threshold:
        return f"<{threshold}"
    return str(n)


def suppress_df(df, cols, threshold=SUPPRESS_BELOW):
    out = df.copy()
    for c in cols:
        if c in out.columns:
            out[c] = out[c].apply(lambda v: suppress(v, threshold))
    return out


def load_allele_counts(path):
    if not os.path.exists(path):
        sys.exit(f"FATAL: allele counts table (24's output) not found at {path!r}. "
                  f"Run 24_novelty_by_field.py first.")
    df = pd.read_csv(path, sep="\t", dtype=str)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        sys.exit(f"FATAL: {path} is missing column(s) {missing}. "
                  f"Available columns: {list(df.columns)}")
    df["n_haplotypes"] = pd.to_numeric(df["n_haplotypes"], errors="coerce")
    df["unrelated_only"] = df["unrelated_only"].astype(str).str.strip().str.lower().isin(
        {"true", "1", "t", "yes"})
    df["ancestry"] = vc.normalize_ancestry(df["ancestry"])
    return df


def is_known_allele(allele_id):
    """24's `allele_ids()` gives known (already-catalogued IPD-IMGT) alleles a name of the
    form '<gene_display>*<fields>' (contains '*'); novel-cluster / non-functional ids are
    '<gene_bare>_prot_<hash>', '<gene_bare>_cds_<hash>', or '<gene_bare>_nonfunc_<hash>'
    (no '*'). See 24_novelty_by_field.normalize_allele_name / allele_ids."""
    if not isinstance(allele_id, str):
        return False
    return "*" in allele_id


def filter_scope(df, resolution, exclude_flagged, ancestry_scheme, unrelated_only):
    d = df[(df["resolution"] == resolution) & (df["exclude_flagged"] == exclude_flagged)
           & (df["ancestry_scheme"] == ancestry_scheme)
           & (df["unrelated_only"] == unrelated_only)]
    return d[d["ancestry"] != "POOLED"].copy()


def counts_by_gene_ancestry(scoped_df):
    """scoped_df: rows already filtered to one (resolution, exclude_flagged, scheme,
    unrelated_only). Returns {(gene, ancestry): {allele_id: n_haplotypes}}."""
    out = {}
    for (gene, anc), g in scoped_df.groupby(["gene", "ancestry"], observed=True):
        out[(gene, anc)] = dict(zip(g["allele_id"], g["n_haplotypes"]))
    return out


# ---------------------------------------------------------------------------
# 1. coverage_by_gene_ancestry
# ---------------------------------------------------------------------------
def build_coverage_table(scoped_df, resolution, n_boot=N_BOOT, genes=None):
    rows = []
    by_ga = counts_by_gene_ancestry(scoped_df)
    for (gene, anc), counts in sorted(by_ga.items()):
        if genes is not None and gene not in genes:
            continue
        X = np.array(list(counts.values()), dtype=float)
        n = int(X.sum())
        s_obs = int(len(X))
        f1 = int(np.sum(X == 1))
        f2 = int(np.sum(X == 2))
        chat = cov.sample_coverage(X)
        lo, hi = cov.bootstrap_ci(X, cov.sample_coverage, n_boot=n_boot, seed=0)
        chao1_f0 = cov._chao1_f0(n, f1, f2) if n > 0 else np.nan
        chao1 = s_obs + chao1_f0 if n > 0 else np.nan
        rec = {"gene": gene, "resolution": resolution, "ancestry": anc,
               "n_haplotypes": n, "S_obs": s_obs, "f1": f1, "f2": f2,
               "coverage": chat, "coverage_ci_lo": lo, "coverage_ci_hi": hi,
               "unseen_mass": 1.0 - chat if not np.isnan(chat) else np.nan,
               "chao1_richness": chao1, "chao1_is_secondary_unstable": True}
        for target in TARGET_COVERAGES:
            m, reliable = cov.size_for_coverage(X, target) if n > 0 else (np.inf, False)
            pct = int(round(target * 1000)) / 10.0
            rec[f"n_for_{pct}pct_coverage"] = m
            rec[f"n_for_{pct}pct_coverage_reliable"] = reliable
            rec[f"n_more_for_{pct}pct_coverage"] = (
                (m - n) if np.isfinite(m) else np.inf)
        rows.append(rec)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 2. coverage curves (interpolate to n, extrapolate to 2n)
# ---------------------------------------------------------------------------
def build_coverage_curves(scoped_df, resolution, n_points=N_CURVE_POINTS, genes=None):
    rows = []
    by_ga = counts_by_gene_ancestry(scoped_df)
    for (gene, anc), counts in sorted(by_ga.items()):
        if genes is not None and gene not in genes:
            continue
        X = np.array(list(counts.values()), dtype=float)
        n = int(X.sum())
        if n == 0:
            continue
        m_max = 2 * n
        n_pts = min(n_points, m_max) if m_max > 0 else 1
        sizes = np.unique(np.round(np.linspace(1, m_max, max(n_pts, 1)))).astype(int)
        sizes = sizes[sizes >= 1]
        cvals = cov.coverage_at_size(X, sizes)
        cvals = np.atleast_1d(cvals)
        for m, c in zip(sizes, cvals):
            rows.append({"gene": gene, "resolution": resolution, "ancestry": anc,
                         "n_haplotypes": n, "m": int(m), "coverage": float(c),
                         "is_extrapolated": bool(m > n)})
    return pd.DataFrame(rows)


def fig_coverage_curves(curves_df, path, resolution):
    d = curves_df[curves_df["resolution"] == resolution]
    genes = sorted(d["gene"].unique())
    if not genes:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        vc.savefig(fig, path, dpi=200)
        return
    ncols = min(4, len(genes))
    nrows = int(np.ceil(len(genes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.0 * ncols, 2.4 * nrows), squeeze=False)
    for i, gene in enumerate(genes):
        ax = axes[i // ncols][i % ncols]
        dg = d[d["gene"] == gene]
        for anc in vc.ANCESTRY_ORDER:
            da = dg[dg["ancestry"] == anc].sort_values("m")
            if da.empty:
                continue
            color = vc.ANCESTRY_COLORS.get(anc, vc.MISSING_COLOR)
            obs = da[~da["is_extrapolated"]]
            ext = da[da["is_extrapolated"]]
            ax.plot(obs["m"], obs["coverage"], color=color, lw=1.2, label=anc)
            if len(ext):
                bridge_m = pd.concat([obs["m"].tail(1), ext["m"]])
                bridge_c = pd.concat([obs["coverage"].tail(1), ext["coverage"]])
                ax.plot(bridge_m, bridge_c, color=color, lw=1.0, linestyle="--")
        ax.set_title(gene, fontsize=8)
        ax.set_ylim(0, 1.02)
        ax.tick_params(labelsize=6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    for j in range(len(genes), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    handles = [plt.Line2D([0], [0], color=vc.ANCESTRY_COLORS[a], lw=1.5, label=a)
               for a in vc.ANCESTRY_ORDER]
    fig.legend(handles=handles, loc="upper center", ncol=len(vc.ANCESTRY_ORDER), fontsize=7,
               bbox_to_anchor=(0.5, 1.02), frameon=False)
    fig.suptitle(f"Coverage curves ({resolution}) -- solid = interpolated, dashed = "
                 f"extrapolated to 2n", fontsize=8, y=1.06)
    fig.tight_layout()
    vc.savefig(fig, path, dpi=200)


# ---------------------------------------------------------------------------
# 3. cross-ancestry catalogue matrix
# ---------------------------------------------------------------------------
def build_cross_ancestry_matrix(scoped_df, resolution, genes=None):
    rows = []
    by_ga = counts_by_gene_ancestry(scoped_df)
    genes_present = sorted(set(g for g, _ in by_ga.keys()))
    if genes is not None:
        genes_present = [g for g in genes_present if g in genes]
    for gene in genes_present:
        gene_counts = {anc: counts for (g, anc), counts in by_ga.items() if g == gene}
        ancs = sorted(gene_counts.keys())
        if len(ancs) < 2:
            continue
        m = min(int(sum(gene_counts[a].values())) for a in ancs)
        if m <= 0:
            continue
        for a_ref in ancs:
            for a_target in ancs:
                cc, unseen_b = cov.cross_coverage(gene_counts[a_ref], gene_counts[a_target],
                                                   m_ref=m)
                rows.append({"gene": gene, "resolution": resolution, "ancestry_ref": a_ref,
                             "ancestry_target": a_target, "m_ref": m,
                             "cross_coverage": cc, "unseen_in_target": unseen_b,
                             "is_diagonal": a_ref == a_target})
    df = pd.DataFrame(rows)
    if len(df):
        classical = df[df["gene"].isin(vc.CLASSICAL_GENES) | df["gene"].isin(
            [g.replace("HLA-", "") for g in vc.CLASSICAL_GENES])]
        if len(classical):
            mean_rows = (classical.groupby(["ancestry_ref", "ancestry_target", "is_diagonal"],
                                            observed=True)["cross_coverage"].mean()
                         .reset_index())
            mean_rows.insert(0, "gene", "classical_gene_mean")
            mean_rows.insert(1, "resolution", resolution)
            mean_rows["m_ref"] = np.nan
            mean_rows["unseen_in_target"] = np.nan
            df = pd.concat([df, mean_rows], ignore_index=True, sort=False)
    return df


def fig_cross_ancestry_matrix(matrix_df, path, resolution):
    d = matrix_df[(matrix_df["resolution"] == resolution)
                  & (matrix_df["gene"] == "classical_gene_mean")]
    ancs = [a for a in vc.ANCESTRY_ORDER if a in set(d["ancestry_ref"]) | set(d["ancestry_target"])]
    if not ancs:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        vc.savefig(fig, path, dpi=200)
        return
    mat = np.full((len(ancs), len(ancs)), np.nan)
    for i, a_ref in enumerate(ancs):
        for j, a_tgt in enumerate(ancs):
            row = d[(d["ancestry_ref"] == a_ref) & (d["ancestry_target"] == a_tgt)]
            if len(row):
                mat[i, j] = row["cross_coverage"].iloc[0]
    fig, ax = plt.subplots(figsize=(1.0 + 0.6 * len(ancs), 1.0 + 0.6 * len(ancs)))
    im = ax.imshow(mat, vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(range(len(ancs)))
    ax.set_xticklabels(ancs, fontsize=7)
    ax.set_yticks(range(len(ancs)))
    ax.set_yticklabels(ancs, fontsize=7)
    ax.set_xlabel("catalogued population (target)", fontsize=7)
    ax.set_ylabel("reference catalogue built from", fontsize=7)
    for i in range(len(ancs)):
        for j in range(len(ancs)):
            if np.isfinite(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if mat[i, j] < 0.6 else "black")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title(f"Cross-ancestry catalogue coverage, classical-gene mean ({resolution})",
                 fontsize=8)
    vc.savefig(fig, path, dpi=200)


# ---------------------------------------------------------------------------
# 4. IMGT coverage today
# ---------------------------------------------------------------------------
def build_imgt_coverage_today(scoped_df, resolution, genes=None):
    rows = []
    by_ga = counts_by_gene_ancestry(scoped_df)
    for (gene, anc), counts in sorted(by_ga.items()):
        if genes is not None and gene not in genes:
            continue
        n = sum(counts.values())
        if n == 0:
            continue
        known_n = sum(v for k, v in counts.items() if is_known_allele(k))
        rows.append({"gene": gene, "resolution": resolution, "ancestry": anc,
                     "n_haplotypes": int(n),
                     "frac_in_imgt_today": known_n / n})
    return pd.DataFrame(rows)


def fig_imgt_coverage_today(df, path, resolution):
    d = df[df["resolution"] == resolution]
    genes = sorted(d["gene"].unique())
    if not genes:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        vc.savefig(fig, path, dpi=200)
        return
    ancs = [a for a in vc.ANCESTRY_ORDER if a in set(d["ancestry"])]
    width = 0.8 / max(len(ancs), 1)
    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(genes)), 3.5))
    x = np.arange(len(genes))
    for ai, anc in enumerate(ancs):
        vals = []
        for g in genes:
            row = d[(d["gene"] == g) & (d["ancestry"] == anc)]
            vals.append(row["frac_in_imgt_today"].iloc[0] if len(row) else np.nan)
        ax.bar(x + (ai - (len(ancs) - 1) / 2) * width, vals, width=width,
               color=vc.ANCESTRY_COLORS.get(anc, vc.MISSING_COLOR), label=anc)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, fontsize=7, rotation=45, ha="right")
    ax.set_ylabel("fraction already in IPD-IMGT/HLA", fontsize=8)
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=6, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title(f"Reference-database coverage today ({resolution})", fontsize=9)
    vc.savefig(fig, path, dpi=200)


# ---------------------------------------------------------------------------
# 5. sampling-design scenarios at fixed total N
# ---------------------------------------------------------------------------
def greedy_allocation(n_by_pop, total_n, step=None):
    """Greedy allocation of `total_n` haplotypes across populations (capped at each
    population's observed n) that approximately maximizes the *unweighted mean* coverage
    achievable across populations. Not from `_coverage.py` (that module has no allocator);
    here we just decide sizes, then hand them to `pooled_design_coverage`."""
    pops = sorted(n_by_pop.keys())
    if step is None:
        step = max(1, int(total_n / 200))
    sizes = dict((p, 0) for p in pops)
    remaining = total_n
    caps = dict(n_by_pop)
    # start everyone at 0; greedily add `step` haplotypes to whichever population's
    # ancestry-mean coverage benefits most per haplotype (marginal gain), similar in
    # spirit to a knapsack greedy heuristic. Since exact marginal coverage gain per pop
    # needs its own catalogue (counts), the caller supplies coverage_fn via closure below;
    # this generic version equalizes proportional to remaining capacity as a simple,
    # deterministic fallback used only when no coverage function is given.
    while remaining > 0 and any(sizes[p] < caps[p] for p in pops):
        available = [p for p in pops if sizes[p] < caps[p]]
        take = min(step, remaining, min(caps[p] - sizes[p] for p in available))
        take = max(take, 1)
        for p in available:
            if remaining <= 0:
                break
            add = min(take, caps[p] - sizes[p], remaining)
            sizes[p] += add
            remaining -= add
    return sizes


def greedy_allocation_by_gain(counts_by_pop, total_n, target_weights, n_steps=40):
    """Greedy allocation that iteratively assigns chunks of `total_n` to whichever
    population currently has the largest marginal coverage gain per haplotype added,
    evaluated via `pooled_design_coverage`. Sizes never exceed observed per-pop n."""
    pops = sorted(counts_by_pop.keys())
    n_by_pop = {p: int(sum(counts_by_pop[p].values())) for p in pops}
    cap = sum(min(n_by_pop[p], total_n) for p in pops)
    total_n = min(total_n, cap)
    step = max(1, int(total_n / n_steps))
    sizes = dict((p, 0) for p in pops)
    remaining = total_n
    while remaining > 0:
        candidates = [p for p in pops if sizes[p] < n_by_pop[p]]
        if not candidates:
            break
        best_p, best_gain = None, -np.inf
        _, base_global = pooled_design_coverage_safe(counts_by_pop, sizes, target_weights)
        for p in candidates:
            trial = dict(sizes)
            add = min(step, remaining, n_by_pop[p] - sizes[p])
            trial[p] = trial[p] + add
            _, g = pooled_design_coverage_safe(counts_by_pop, trial, target_weights)
            gain = g - base_global
            if gain > best_gain:
                best_gain, best_p, best_add = gain, p, add
        if best_p is None:
            break
        sizes[best_p] += best_add
        remaining -= best_add
    return sizes


def pooled_design_coverage_safe(counts_by_pop, sizes_by_pop, target_weights):
    sizes_by_pop = dict((p, s) for p, s in sizes_by_pop.items() if s > 0)
    return cov.pooled_design_coverage(counts_by_pop, sizes_by_pop, target_weights)


def build_design_scenarios(scoped_df, resolution, genes=None):
    """One row per (gene, design, target ancestry) at fixed total N = the observed total
    unrelated LR haplotype count for that gene (sum over ancestries, at this resolution)."""
    rows = []
    by_ga = counts_by_gene_ancestry(scoped_df)
    genes_present = sorted(set(g for g, _ in by_ga.keys()))
    if genes is not None:
        genes_present = [g for g in genes_present if g in genes]
    for gene in genes_present:
        counts_by_pop = {anc: counts for (g, anc), counts in by_ga.items() if g == gene}
        ancs = sorted(counts_by_pop.keys())
        if len(ancs) < 2:
            continue
        n_by_pop = {a: int(sum(counts_by_pop[a].values())) for a in ancs}
        total_n = sum(n_by_pop.values())
        if total_n <= 0:
            continue
        weights = {a: n_by_pop[a] / total_n for a in ancs}

        actual_sizes = dict(n_by_pop)
        equal_n = total_n // len(ancs)
        equal_sizes = dict((a, min(equal_n, n_by_pop[a])) for a in ancs)
        eur_key = "EUR" if "EUR" in ancs else None
        if eur_key:
            eur_target = int(round(0.94 * total_n))
            eur_sizes = {eur_key: min(eur_target, n_by_pop[eur_key])}
            rest = total_n - eur_sizes[eur_key]
            other_ancs = [a for a in ancs if a != eur_key]
            per_other = rest // len(other_ancs) if other_ancs else 0
            for a in other_ancs:
                eur_sizes[a] = min(per_other, n_by_pop[a])
        else:
            eur_sizes = None
        greedy_sizes = greedy_allocation_by_gain(counts_by_pop, total_n, weights)

        designs = [("actual_aou_composition", actual_sizes),
                   ("equal_allocation", equal_sizes)]
        if eur_sizes is not None:
            designs.append(("94pct_eur_ukb_like", eur_sizes))
        else:
            designs.append(("94pct_eur_ukb_like", None))
        designs.append(("greedy_max_mean_coverage", greedy_sizes))

        for design_name, sizes in designs:
            if sizes is None:
                for a in ancs:
                    rows.append({"gene": gene, "resolution": resolution, "design": design_name,
                                 "target_ancestry": a, "total_n": total_n,
                                 "coverage": np.nan,
                                 "note": "not reachable with this cohort (no EUR group)"})
                continue
            needs_more = any(s > n_by_pop[a] for a, s in sizes.items())
            note = ""
            if needs_more:
                note = "not reachable with this cohort (design needs more than observed n)"
                sizes = dict((a, min(s, n_by_pop[a])) for a, s in sizes.items())
            cov_by_pop, global_cov = pooled_design_coverage_safe(counts_by_pop, sizes, weights)
            for a in ancs:
                rows.append({"gene": gene, "resolution": resolution, "design": design_name,
                             "target_ancestry": a, "total_n": total_n,
                             "size_used": sizes.get(a, 0),
                             "coverage": cov_by_pop.get(a, np.nan), "note": note})
            rows.append({"gene": gene, "resolution": resolution, "design": design_name,
                         "target_ancestry": "GLOBAL_MEAN", "total_n": total_n,
                         "size_used": sum(sizes.values()),
                         "coverage": global_cov, "note": note})
    return pd.DataFrame(rows)


def fig_design_scenarios(design_df, path, resolution):
    d = design_df[(design_df["resolution"] == resolution)
                  & (design_df["target_ancestry"] == "GLOBAL_MEAN")]
    if not len(d):
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        vc.savefig(fig, path, dpi=200)
        return
    genes = sorted(d["gene"].unique())
    designs = ["actual_aou_composition", "equal_allocation", "94pct_eur_ukb_like",
               "greedy_max_mean_coverage"]
    colors = {"actual_aou_composition": "#555555", "equal_allocation": "#0072B2",
              "94pct_eur_ukb_like": "#D55E00", "greedy_max_mean_coverage": "#009E73"}
    width = 0.8 / len(designs)
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(genes)), 3.5))
    x = np.arange(len(genes))
    for di, design in enumerate(designs):
        vals = []
        for g in genes:
            row = d[(d["gene"] == g) & (d["design"] == design)]
            vals.append(row["coverage"].iloc[0] if len(row) and pd.notna(
                row["coverage"].iloc[0]) else 0.0)
        ax.bar(x + (di - (len(designs) - 1) / 2) * width, vals, width=width,
               color=colors[design], label=design)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, fontsize=7, rotation=45, ha="right")
    ax.set_ylabel("mean coverage across ancestries", fontsize=8)
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=6, frameon=False, ncol=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title(f"Sampling-design scenarios at fixed total N ({resolution})", fontsize=9)
    vc.savefig(fig, path, dpi=200)


# ---------------------------------------------------------------------------
# 6. plain-language sentences
# ---------------------------------------------------------------------------
def coverage_sentences(coverage_df, resolution="protein"):
    d = coverage_df[coverage_df["resolution"] == resolution]
    lines = []
    for _, r in d.sort_values("ancestry").iterrows():
        if pd.isna(r["coverage"]):
            continue
        pct = 100 * r["coverage"]
        m99 = r.get("n_for_99.0pct_coverage", np.nan)
        n = r["n_haplotypes"]
        if np.isfinite(m99) and m99 > n:
            more = f"about {int(round(m99 - n))} more"
        elif np.isfinite(m99):
            more = "no more (already reached)"
        else:
            more = "an unknown (unreliable-to-estimate) number more"
        lines.append(
            f"{r['ancestry']}: a catalogue from these {int(n)} haplotypes types "
            f"{pct:.1f}% of haplotypes; reaching 99% needs {more}.")
    return lines


# ---------------------------------------------------------------------------
# report / summary
# ---------------------------------------------------------------------------
DOUBTS = [
    "Coverage estimates assume haplotypes are independent draws. Relatives are removed, "
    "but population structure within an \"ancestry\" label mildly violates this "
    "assumption; strict ancestry labels (not `ancestry_pred`) are used as the primary "
    "analysis to reduce that risk.",
    "Novel-protein/CDS cluster identity depends on script 24's artifact labels; this "
    "report defaults to `clean_only` (flagged calls excluded) -- rerun with "
    "`--exclude-flagged all` to see the sensitivity analysis.",
    "AoU is not a random sample of any real-world population; every coverage number here "
    "is \"coverage of the AoU-sampled population\", not of a general human population.",
    "Extrapolated sizes for 95/99/99.9% coverage are only trustworthy out to about 2x the "
    "observed sample size (Chao & Jost); rows beyond that are flagged 'not reliable' and "
    "should be read as a lower bound only.",
    "Chao1 richness is reported only as a secondary, explicitly unstable number -- do not "
    "use it to claim 'there are exactly N alleles'; use coverage for that question.",
]


def write_report(path, coverage_df, cross_df, imgt_df, design_df, args):
    L = []
    L.append("# How much of the HLA allele space have we seen?")
    L.append("")
    L.append("_Generated by `27_allele_space_coverage.py`. "
             "Spec: `sprints/S01_novelty_qc_coverage/WS3_allele_space_coverage.md`._")
    L.append("")
    L.append("## What \"coverage\" means, in plain language")
    L.append("")
    L.append("Every person's HLA gene copy (haplotype) carries one allele. If we had typed "
              "every haplotype in a population, we would have a complete catalogue. We "
              "haven't -- we've typed a sample. **Sample coverage** is the fraction of the "
              "population's haplotypes whose allele is already in our catalogue; "
              "equivalently, it is the chance that the *next* new haplotype we type carries "
              "an allele we've already seen. A catalogue with 99% coverage will correctly "
              "recognize 99 out of the next 100 haplotypes typed against it.")
    L.append("")
    L.append("## Why coverage instead of \"how many alleles exist\" (richness)")
    L.append("")
    L.append("A natural first question is \"how many distinct alleles are there, total?\" "
              "(richness). But HLA allele-frequency distributions have a long tail of very "
              "rare alleles seen only once or twice (singletons/doubletons); estimating the "
              "total number that exist requires extrapolating that tail, which earlier "
              "analyses in this project (04/18) found to be unstable ('low confidence "
              "Chao2' on nearly every row). Coverage sidesteps this: it only depends on how "
              "many alleles were seen exactly once or twice in the sample we already have, "
              "which is a much better-behaved quantity. Richness (Chao1) is still reported "
              "here, but only as a secondary, explicitly-flagged number.")
    L.append("")
    L.append("## Headline sentences (protein resolution)")
    L.append("")
    for line in coverage_sentences(coverage_df, "protein"):
        L.append(f"- {line}")
    L.append("")
    L.append("## What's in each output file")
    L.append("")
    L.append("- `coverage_by_gene_ancestry.tsv`: per gene x resolution x ancestry -- sample "
              "size, observed allele count, singleton/doubleton counts, estimated coverage "
              "with a bootstrap 95% interval, unseen frequency mass, secondary Chao1 "
              "richness, and how many haplotypes it would take to reach 95%/99%/99.9% "
              "coverage (with a reliability flag).")
    L.append("- `coverage_curves_protein.png` / `coverage_curves_cds.png` + "
              "`coverage_curves.tsv`: coverage as sample size grows, per ancestry, small "
              "multiples per gene. Solid = interpolated (within observed data), dashed = "
              "extrapolated (beyond it, out to double the sample).")
    L.append("- `cross_ancestry_matrix.tsv` / `.png`: if you build a catalogue from one "
              "ancestry's haplotypes, how well does it type another ancestry's haplotypes? "
              "The diagonal is \"self\" coverage; off-diagonal cells below the diagonal show "
              "reference-catalogue bias.")
    L.append("- `imgt_coverage_today.tsv` / `.png`: the fraction of each ancestry's "
              "haplotypes whose allele is already a known, named entry in IPD-IMGT/HLA "
              "today (versus a novel cluster this project's own sequencing discovered).")
    L.append("- `design_scenarios.tsv` / `.png`: at a fixed total sample size (this "
              "cohort's actual unrelated haplotype count), how would different ways of "
              "splitting that sample across ancestries change mean coverage? Compares the "
              "actual AoU mix, an equal split, a 94%-European mix (UK-Biobank-like), and a "
              "greedy allocation chosen to maximize mean coverage across ancestries.")
    L.append("")
    L.append("## Caveats")
    L.append("")
    for d in DOUBTS:
        L.append(f"- {d}")
    L.append("")
    L.append("## Run parameters")
    L.append("")
    L.append(f"- Input: `{args.allele_counts}`")
    L.append(f"- `exclude_flagged`: `{args.exclude_flagged}`")
    L.append(f"- `ancestry_scheme`: `{args.ancestry_scheme}`")
    L.append(f"- `unrelated_only`: `{args.unrelated_only}`")
    L.append(f"- Counts below {SUPPRESS_BELOW} are suppressed as `<{SUPPRESS_BELOW}` "
             f"throughout; no allele-level rows are written anywhere in this output.")
    vc.write_report(path, L)


def build_summary(coverage_df, cross_df, imgt_df, design_df, args):
    prot = coverage_df[coverage_df["resolution"] == "protein"]
    summary = {
        "params": {"allele_counts": args.allele_counts, "exclude_flagged": args.exclude_flagged,
                   "ancestry_scheme": args.ancestry_scheme,
                   "unrelated_only": bool(args.unrelated_only),
                   "limit_genes": args.limit_genes},
        "n_gene_ancestry_rows": {"protein": int(len(prot)),
                                  "cds": int(len(coverage_df[coverage_df["resolution"] == "cds"]))},
        "mean_protein_coverage_by_ancestry": (
            prot.groupby("ancestry")["coverage"].mean().round(4).to_dict() if len(prot) else {}),
        "sentences_protein": coverage_sentences(coverage_df, "protein"),
        "n_cross_ancestry_rows": int(len(cross_df)),
        "n_imgt_rows": int(len(imgt_df)),
        "n_design_rows": int(len(design_df)),
    }
    return summary


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--allele-counts", default=DEFAULT_ALLELE_COUNTS,
                     help="24's allele_counts_by_resolution.tsv")
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    ap.add_argument("--exclude-flagged", choices=["clean_only", "all"], default="clean_only",
                     help="clean_only (default): drop artifact-flagged calls; all: sensitivity.")
    ap.add_argument("--ancestry-scheme", choices=["strict", "pred"], default="strict",
                     help="strict (default, admixture proportion >= 0.9) or pred (ancestry_pred).")
    unrel = ap.add_mutually_exclusive_group()
    unrel.add_argument("--unrelated-only", dest="unrelated_only", action="store_true",
                        default=True, help="Unrelated people only (default).")
    unrel.add_argument("--no-unrelated-only", dest="unrelated_only", action="store_false",
                        help="Include related people too.")
    ap.add_argument("--limit-genes", type=int, default=None,
                     help="Only analyze the first N genes (sorted name) -- fast local runs.")
    ap.add_argument("--n-boot", type=int, default=N_BOOT,
                     help="Bootstrap replicates for coverage CIs.")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    vc.ensure_dir(args.out_dir)

    df = load_allele_counts(args.allele_counts)

    genes_all = sorted(df["gene"].unique())
    genes = set(genes_all[:args.limit_genes]) if args.limit_genes else None

    cov_parts, curve_parts, cross_parts, imgt_parts = [], [], [], []
    for resolution in ("protein", "cds"):
        scoped = filter_scope(df, resolution, args.exclude_flagged, args.ancestry_scheme,
                               args.unrelated_only)
        if not len(scoped):
            print(f"WARNING: no rows for resolution={resolution!r} with the given filters",
                  file=sys.stderr)
            continue
        cov_parts.append(build_coverage_table(scoped, resolution, n_boot=args.n_boot,
                                               genes=genes))
        curve_parts.append(build_coverage_curves(scoped, resolution, genes=genes))
        cross_parts.append(build_cross_ancestry_matrix(scoped, resolution, genes=genes))
        imgt_parts.append(build_imgt_coverage_today(scoped, resolution, genes=genes))

    coverage_df = (pd.concat(cov_parts, ignore_index=True) if cov_parts
                   else pd.DataFrame(columns=["gene", "resolution", "ancestry"]))
    curves_df = (pd.concat(curve_parts, ignore_index=True) if curve_parts
                 else pd.DataFrame(columns=["gene", "resolution", "ancestry", "m", "coverage"]))
    cross_df = (pd.concat(cross_parts, ignore_index=True) if cross_parts
                else pd.DataFrame(columns=["gene", "resolution"]))
    imgt_df = (pd.concat(imgt_parts, ignore_index=True) if imgt_parts
               else pd.DataFrame(columns=["gene", "resolution", "ancestry"]))

    # Design scenarios: protein resolution only (the classical clinical resolution), per WS3.
    scoped_protein = filter_scope(df, "protein", args.exclude_flagged, args.ancestry_scheme,
                                   args.unrelated_only)
    design_df = (build_design_scenarios(scoped_protein, "protein", genes=genes) if len(scoped_protein)
                 else pd.DataFrame(columns=["gene", "resolution", "design", "target_ancestry"]))

    # ---- write TSVs (suppress raw counts <20; no allele-level rows anywhere) ----
    suppress_df(coverage_df, ["n_haplotypes", "S_obs"]).to_csv(
        os.path.join(args.out_dir, "coverage_by_gene_ancestry.tsv"), sep="\t", index=False)
    curves_df.to_csv(os.path.join(args.out_dir, "coverage_curves.tsv"), sep="\t", index=False)
    cross_df.to_csv(os.path.join(args.out_dir, "cross_ancestry_matrix.tsv"), sep="\t", index=False)
    suppress_df(imgt_df, ["n_haplotypes"]).to_csv(
        os.path.join(args.out_dir, "imgt_coverage_today.tsv"), sep="\t", index=False)
    suppress_df(design_df, ["size_used"]).to_csv(
        os.path.join(args.out_dir, "design_scenarios.tsv"), sep="\t", index=False)

    # ---- figures ----
    fig_coverage_curves(curves_df, os.path.join(args.out_dir, "coverage_curves_protein.png"),
                         "protein")
    fig_coverage_curves(curves_df, os.path.join(args.out_dir, "coverage_curves_cds.png"), "cds")
    fig_cross_ancestry_matrix(cross_df, os.path.join(args.out_dir, "cross_ancestry_matrix.png"),
                               "protein")
    fig_imgt_coverage_today(imgt_df, os.path.join(args.out_dir, "imgt_coverage_today.png"),
                             "protein")
    fig_design_scenarios(design_df, os.path.join(args.out_dir, "design_scenarios.png"), "protein")

    # ---- report + summary ----
    write_report(os.path.join(args.out_dir, "allele_space_coverage_report.md"),
                 coverage_df, cross_df, imgt_df, design_df, args)
    summary = build_summary(coverage_df, cross_df, imgt_df, design_df, args)
    with open(os.path.join(args.out_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"Wrote outputs to {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
