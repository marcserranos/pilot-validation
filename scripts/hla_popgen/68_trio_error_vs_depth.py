#!/usr/bin/env python3
"""S05 68: does thin long-read coverage explain the errors trios reveal, and does pruning low-coverage people help?

Marc (2026-10-13): plot the depth of one relative near a call against the depth of the other relative, coloured by
whether the call agrees; give an error measure; prune the people with the worst coverage step by step and see whether
the quality metrics improve. 59 already showed depth is ASSOCIATED with KIR/MHC outcomes, but those outcomes have no
truth (a missing gene may be real). Trios do have truth, so this asks the question where it can be answered:

  truth 1  Mendelian typing (57b): the parent a child copy came from (parental k-mers, +-10 kb) must carry the
           child's allele. A mismatch is an error in the child's OR the parent's call (both are long-read calls).
  truth 2  physical phasing (57b): two genes on one child contig must come from the same parent.
  truth 3  switch / hamming error (57, 50 trios): parental k-mer flips along the child's MHC haplotype pieces.

Depth = 58's rel_depth (region mean depth / the person's own genome-wide control depth, MAPQ >= 0), normalised by
the cohort median of that region x platform (59's `norm` scale: 1.0 = typical for that region). "Near the call"
= the gene's own GENCODE span (58 regions); person-level = the 4 Mb MHC window.

Questions, and why each is the right form of Marc's (naive) ask:
  (a) error vs depth of the child at the gene and of the transmitting parent at the gene: 2D binned rates (a point
      per call would publish per-person depth and genotype status, so the scatter stays INTERNAL on the VM);
      logistic OR per 0.1 lower depth with both depths in one model (which side's coverage matters?).
  (b) pruning: drop the fraction f of trios whose child has the lowest MHC depth (and, separately, of calls whose
      gene depth is lowest in child or parent) and recompute each error rate. A drop from pruning means nothing
      alone: removing any people changes the rate by chance, so every step is compared with random removal of
      the same number (1,000 draws; 2.5-97.5% band). Pruning helps only if it beats random.
  (c) phasing errors vs child MHC depth and contig count (do switch-prone children have thin coverage, or
      fragmented assemblies?).

Inputs (VM, INTERNAL): ~/s05/internal/57b/gene_origin.INTERNAL.tsv, ~/s05/internal/57/per_trio_region.INTERNAL.tsv,
~/s05/internal/58/per_person_depth.tsv, Table 1 (~/pipeline_outputs/hla_calls_rich.tsv), the cohort file and the
relatedness pairs (cached LR list when the AoU table is gone).
PUBLIC (~/s05/results/68/): mendel_by_depth_2d.tsv, error_by_depth_bin.tsv, depth_logit.tsv, prune_curve.tsv,
phase_vs_depth.tsv, fig_68_*.png. Counts 1-19 masked; rates blanked when numerator or denominator is 1-19; prune
rows blanked where the removed set's error count (cumulative or step) is 1-19, so differences between rows cannot
give a small count. INTERNAL (~/s05/internal/68/): per-call table + the point scatter (VM only).

    cd ~/repos/<checkout>/scripts/hla_popgen && $PY 68_trio_error_vs_depth.py
Python 3.8 compatible.
"""
import argparse
import importlib.util
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _s05_common as C  # noqa: E402

OUT = os.path.expanduser("~/s05/results/68/")
INT = os.path.expanduser("~/s05/internal/68/")
GO = os.path.expanduser("~/s05/internal/57b/gene_origin.INTERNAL.tsv")
TRIOEVAL = os.path.expanduser("~/s05/internal/57/per_trio_region.INTERNAL.tsv")
DEPTH = os.path.expanduser("~/s05/internal/58/per_person_depth.tsv")
KIN_LO, KIN_HI = 2 ** -2.5, 2 ** -1.5
FRACS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]
BINS = [0.0, 0.8, 0.9, 1.0, 1.1, np.inf]   # normalised depth (1.0 = cohort median for that region x platform)
SMALL = 20


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def small(*v):
    return any(0 < int(x) < SMALL for x in v)


def rate(k, n, digits=1):
    """(n display, pct, lo, hi) with PUBLIC masking."""
    k, n = int(k), int(n)
    if n == 0 or small(k, n):
        return (C.mask_n(n), "", "", "")
    p, lo, hi = C.wilson(k, n)
    f = f"%.{digits}f"
    return (str(n), f % (100 * p), f % (100 * lo), f % (100 * hi))


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def find_trios(children, pairs):
    """child -> (parent1, parent2): two first-degree relatives of the child not related to each other (57's rule)."""
    nb, anyp = defaultdict(set), set()
    for i, j, k in pairs:
        anyp.add((i, j))
        anyp.add((j, i))
        if KIN_LO <= k <= KIN_HI:
            nb[i].add(j)
            nb[j].add(i)
    out = {}
    for c in children:
        ps = sorted(nb.get(c, ()))
        for a in range(len(ps)):
            for b in range(a + 1, len(ps)):
                if (ps[a], ps[b]) not in anyp and c not in out:
                    out[c] = (ps[a], ps[b])
    return out


def gene_region(g):
    return g if g in ("MICA", "TAP1", "TAP2") or g.startswith("HLA-") else "HLA-" + g


def depth_table(path, people, plat):
    d = pd.read_csv(path, sep="\t", dtype={"person_id": str}, comment="#")
    d = d[pd.to_numeric(d["mapq"], errors="coerce") == 0].copy()
    d["rel_depth"] = pd.to_numeric(d["rel_depth"], errors="coerce")
    d["platform"] = d["person_id"].map(plat).fillna("unknown")
    # 59's norm scale: per region x platform cohort median (all people), so 1.0 = typical for that region
    med = d.groupby(["region", "platform"])["rel_depth"].median().rename("med")
    d = d.join(med, on=["region", "platform"])
    d["norm"] = d["rel_depth"] / d["med"]
    d = d[d["person_id"].isin(people)]
    return {(p, r): v for p, r, v in zip(d["person_id"], d["region"], d["norm"])}


def per_call(o, calls, trio_of, dep):
    """One row per resolved child call: error flags + child / transmitting parent / other parent depth."""
    key = calls[calls["copy_index"] == 1].set_index(["person_id", "hap", "gene"])
    o = o.join(key[["tf", "consensus", "score"]], on=["child", "hap", "gene"], how="left")
    par = {}
    for (pid, g), d in calls.groupby(["person_id", "gene"]):
        par[(pid, g)] = (set(d["tf"]), set(d["consensus"]))
    rows = []
    for r in o.itertuples(index=False):
        if r.origin == "amb" or pd.isna(r.tf) or r.child not in trio_of:
            continue
        p1, p2 = trio_of[r.child]
        tp, op = (p1, p2) if r.origin == "pat" else (p2, p1)
        if (tp, r.gene) not in par:
            continue
        ptf, pcons = par[(tp, r.gene)]
        reg = gene_region(r.gene)
        rows.append(dict(child=r.child, plat=r.plat, hap=r.hap, gene=r.gene, score=int(r.score),
                         bad2=int(r.tf not in ptf), badfull=int(r.consensus not in pcons),
                         d_child=dep.get((r.child, reg), np.nan), d_tpar=dep.get((tp, reg), np.nan),
                         d_opar=dep.get((op, reg), np.nan),
                         mhc_child=dep.get((r.child, "window_MHC"), np.nan),
                         mhc_tpar=dep.get((tp, "window_MHC"), np.nan)))
    return pd.DataFrame(rows)


def phase_pairs(o):
    """Same-contig gene pairs in a child with both origins resolved: bad = different parents (57b's rule)."""
    pairs = [("A", "B"), ("A", "C"), ("B", "C"), ("DRB1", "DQA1"), ("DRB1", "DQB1"), ("DQA1", "DQB1"),
             ("DPA1", "DPB1"), ("DQB1", "DPB1")]
    idx = {(c, h, g): (ctg, og) for c, h, g, ctg, og in zip(o["child"], o["hap"], o["gene"], o["contig"], o["origin"])}
    rec = []
    for (c, h, g), (ctg, og) in idx.items():
        for ga, gb in pairs:
            if g == ga and (c, h, gb) in idx:
                ctg2, og2 = idx[(c, h, gb)]
                if ctg2 == ctg and "amb" not in (og, og2):
                    rec.append(dict(child=c, hap=h, pair=f"{ga}~{gb}", bad=int(og != og2)))
    return pd.DataFrame(rec, columns=["child", "hap", "pair", "bad"])


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------
def bin_of(x, edges=BINS):
    lab = [f"{a:g}-{b:g}" if np.isfinite(b) else f">={a:g}" for a, b in zip(edges[:-1], edges[1:])]
    i = np.searchsorted(edges, x, side="right") - 1
    return [lab[j] if np.isfinite(v) and 0 <= j < len(lab) else "NA" for v, j in zip(x, i)]


def mendel_2d(pc, metric):
    """Rates by child-gene-depth bin x transmitting-parent-gene-depth bin (3 x 3 tertiles of the pooled values).
    No marginal totals are published on these bins; a grid with exactly one masked cell masks one more."""
    d = pc.dropna(subset=["d_child", "d_tpar"]).copy()
    q = np.nanquantile(np.r_[d["d_child"], d["d_tpar"]], [1 / 3, 2 / 3])
    edges = [0.0, round(float(q[0]), 2), round(float(q[1]), 2), np.inf]
    d["bc"], d["bp"] = bin_of(d["d_child"].to_numpy(), edges), bin_of(d["d_tpar"].to_numpy(), edges)
    rows = []
    for (bc, bp), g in d.groupby(["bc", "bp"]):
        k, n = int(g[metric].sum()), len(g)
        rows.append(dict(metric=metric, child_depth_bin=bc, parent_depth_bin=bp, kk=k, nn=n))
    t = pd.DataFrame(rows)
    t["masked"] = [small(k, n) or n == 0 for k, n in zip(t["kk"], t["nn"])]
    if t["masked"].sum() == 1:   # complementary: one masked cell + the published ALL rate would recover it
        t.loc[t[~t["masked"]].sort_values("kk").index[0], "masked"] = True
    out = []
    for r in t.itertuples(index=False):
        if r.masked:
            out.append(dict(metric=r.metric, child_depth_bin=r.child_depth_bin, parent_depth_bin=r.parent_depth_bin,
                            n_calls="<20*" if r.nn >= SMALL else C.mask_n(r.nn), error_pct="", lo="", hi=""))
        else:
            n, p, lo, hi = rate(r.kk, r.nn)
            out.append(dict(metric=r.metric, child_depth_bin=r.child_depth_bin, parent_depth_bin=r.parent_depth_bin,
                            n_calls=str(int(round(r.nn / 10.0)) * 10), error_pct=p, lo=lo, hi=hi))
    return pd.DataFrame(out), edges


def logit_table(pc, M59):
    """Logistic: error ~ (1 - child depth)/0.1 + (1 - parent depth)/0.1 [+ score], per metric. OR per 0.1 lower depth."""
    rows = []
    for metric in ("bad2", "badfull"):
        for label, cols in (("gene_depth", ["d_child", "d_tpar"]), ("mhc_window_depth", ["mhc_child", "mhc_tpar"])):
            for adj in (False, True):
                d = pc.dropna(subset=cols)
                X = np.column_stack([np.ones(len(d))] + [(1 - d[c].clip(0, 2).to_numpy()) / 0.1 for c in cols]
                                    + ([d["score"].to_numpy(float)] if adj else []))
                y = d[metric].to_numpy(float)
                if len(d) < 50 or y.sum() < SMALL:
                    continue
                fit = M59.fit_logit(X, y)
                if not fit[2]:
                    continue
                beta, se = fit[0], fit[1]
                for i, who in ((1, "child"), (2, "transmitting_parent")):
                    z = beta[i] / se[i] if se[i] > 0 else np.nan
                    rows.append(dict(metric=metric, depth=label, adjusted_for_score=adj, term=who,
                                     n_calls=str(int(round(len(d) / 50.0)) * 50),
                                     or_per_0_1_lower=round(float(np.exp(beta[i])), 3),
                                     lo=round(float(np.exp(beta[i] - 1.96 * se[i])), 3),
                                     hi=round(float(np.exp(beta[i] + 1.96 * se[i])), 3),
                                     p=float(f"{2 * M59._norm_sf(abs(z)):.2g}") if np.isfinite(z) else ""))
    return pd.DataFrame(rows)


def error_by_bin(pc, ph, trio):
    rows = []
    for metric, d, col in (("mendel_2field", pc, "bad2"), ("mendel_full", pc, "badfull")):
        for who, dc in (("child_gene", "d_child"), ("transmitting_parent_gene", "d_tpar"),
                        ("min_child_parent_gene", None)):
            x = d[dc].to_numpy() if dc else np.fmin(d["d_child"].to_numpy(), d["d_tpar"].to_numpy())
            for b, g in d.assign(_b=bin_of(x)).groupby("_b"):
                n, p, lo, hi = rate(g[col].sum(), len(g))
                rows.append(dict(metric=metric, depth_of=who, depth_bin=b, n=n, error_pct=p, lo=lo, hi=hi))
    if len(ph):
        for b, g in ph.assign(_b=bin_of(ph["mhc_child"].to_numpy())).groupby("_b"):
            n, p, lo, hi = rate(g["bad"].sum(), len(g))
            rows.append(dict(metric="phase_pair_diff_parent", depth_of="child_mhc_window", depth_bin=b, n=n,
                             error_pct=p, lo=lo, hi=hi))
    t = pd.DataFrame(rows)
    # complementary: within one metric x depth_of, the bins sum to a published total (57b) -> a single masked
    # bin is recoverable; mask the smallest other bin's rate too
    for _, g in t.groupby(["metric", "depth_of"]):
        m = g[g["error_pct"] == ""]
        if len(m) == 1 and (g["error_pct"] != "").any():
            t.loc[g[g["error_pct"] != ""].index[0], ["error_pct", "lo", "hi"]] = "<20*"
    return t


def prune_curve(units, ycol, order_col, rng, label, n_rand=1000):
    """units: one row per unit (call / trio-hap pair / child) with error flag ycol and the depth to prune on.
    Lowest depth removed first. Returned rates are blanked where the removed set's error count, cumulative or
    since the previous shown row, is 1-19."""
    u = units.dropna(subset=[order_col]).sort_values(order_col).reset_index(drop=True)
    y = u[ycol].to_numpy(int)
    N, K = len(u), int(y.sum())
    rows, last_k = [], 0
    for f in FRACS:
        r = int(round(f * N))
        kept = y[r:]
        removed_k = K - int(kept.sum())
        sims = np.array([y[rng.choice(N, N - r, replace=False)].mean() for _ in range(n_rand)]) if r else None
        blank = small(removed_k) or small(removed_k - last_k)
        n, p, lo, hi = rate(kept.sum(), len(kept), digits=2)
        row = dict(analysis=label, prune_fraction=f, n_kept=str(int(round(len(kept) / 10.0)) * 10),
                   error_pct="" if blank else p, lo="" if blank else lo, hi="" if blank else hi)
        if sims is not None and not blank:
            row.update(random_median_pct=f"{100 * np.median(sims):.2f}",
                       random_lo_pct=f"{100 * np.quantile(sims, 0.025):.2f}",
                       random_hi_pct=f"{100 * np.quantile(sims, 0.975):.2f}",
                       beats_random=bool(kept.mean() < np.quantile(sims, 0.025)))
        rows.append(row)
        if not blank:
            last_k = removed_k
    return rows


def figures(m2d, edges, prune, ebb, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labs = bin_of(np.array([0.5 * (a + min(b, 2)) for a, b in zip(edges[:-1], edges[1:])]), edges)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, metric in zip(axes, ("badfull", "bad2")):
        g = m2d[m2d["metric"] == metric]
        M = np.full((3, 3), np.nan)
        for r in g.itertuples(index=False):
            if r.error_pct not in ("", "<20*") and r.child_depth_bin in labs and r.parent_depth_bin in labs:
                M[labs.index(r.parent_depth_bin), labs.index(r.child_depth_bin)] = float(r.error_pct)
        im = ax.imshow(M, origin="lower", cmap="viridis", vmin=0)
        for i in range(3):
            for j in range(3):
                ax.text(j, i, "masked" if np.isnan(M[i, j]) else f"{M[i, j]:.1f}%", ha="center", va="center",
                        color="w" if np.isnan(M[i, j]) or M[i, j] < np.nanmax(M) * 0.6 else "k", fontsize=9)
        ax.set_xticks(range(3))
        ax.set_xticklabels(labs, fontsize=8)
        ax.set_yticks(range(3))
        ax.set_yticklabels(labs, fontsize=8)
        ax.set_xlabel("child depth at the gene (normalised)")
        ax.set_ylabel("transmitting parent depth at the gene")
        ax.set_title({"badfull": "Mendelian mismatch, full resolution", "bad2": "Mendelian mismatch, 2-field"}[metric],
                     fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="% calls")
    fig.tight_layout()
    fig.savefig(out + "fig_68_mendel_depth_grid.png", dpi=200)
    plt.close(fig)
    p = pd.DataFrame(prune)
    an = list(dict.fromkeys(p["analysis"]))
    fig, axes = plt.subplots(1, len(an), figsize=(3.3 * len(an), 3.4), squeeze=False)
    for ax, a in zip(axes[0], an):
        g = p[(p["analysis"] == a) & (p["error_pct"] != "")]
        f = g["prune_fraction"].astype(float) * 100
        ax.plot(f, g["error_pct"].astype(float), "o-", color="C0", label="prune lowest depth")
        gr = g[g.get("random_median_pct", pd.Series(dtype=str)).notna() & (g.get("random_median_pct", "") != "")] \
            if "random_median_pct" in g else g.iloc[0:0]
        if len(gr):
            fr = gr["prune_fraction"].astype(float) * 100
            ax.fill_between(fr, gr["random_lo_pct"].astype(float), gr["random_hi_pct"].astype(float), color="0.8",
                            label="random removal 95%")
            ax.plot(fr, gr["random_median_pct"].astype(float), "--", color="0.4")
        ax.set_title(a.replace("_", " "), fontsize=9)
        ax.set_xlabel("% removed")
        ax.set_ylabel("error % in the rest")
    axes[0][0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(out + "fig_68_prune_curves.png", dpi=200)
    plt.close(fig)


def internal_scatter(pc, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, metric in zip(axes, ("badfull", "bad2")):
        d = pc.dropna(subset=["d_child", "d_tpar"])
        ok, bad = d[d[metric] == 0], d[d[metric] == 1]
        ax.scatter(ok["d_child"], ok["d_tpar"], s=4, c="0.7", alpha=0.5, label="parent carries allele")
        ax.scatter(bad["d_child"], bad["d_tpar"], s=8, c="C3", alpha=0.9, label="mismatch")
        ax.set_xlim(0, 2)
        ax.set_ylim(0, 2)
        ax.set_xlabel("child depth at gene (norm)")
        ax.set_ylabel("transmitting parent depth at gene (norm)")
        ax.set_title(metric)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=OUT)
    ap.add_argument("--internal-dir", default=INT)
    ap.add_argument("--seed", type=int, default=68)
    a = ap.parse_args(argv)
    out, intd = a.out_dir.rstrip("/") + "/", a.internal_dir.rstrip("/") + "/"
    for d in (out, intd):
        os.makedirs(d, exist_ok=True)
    B = load("b56", "56_phasing_typing_benchmark.py")
    M59 = load("m59", "59_depth_outcome_join.py")
    o = pd.read_csv(GO, sep="\t", dtype={"child": str, "contig": str})
    cohort = C.load_cohort(C.DEFAULT_COHORT)
    plat = dict(zip(cohort["person_id"], cohort["platform"].astype(str).str.lower()))
    trio_of = find_trios(sorted(set(o["child"])), C.load_relatedness_pairs(C.DEFAULT_RELATEDNESS))
    C.log(f"57b children {o['child'].nunique()}, trios re-found {len(trio_of)}")
    people = set(trio_of) | {p for t in trio_of.values() for p in t}
    C.TABLE1_USECOLS = C.TABLE1_USECOLS + ["n_tied"]
    calls = B.prepare(C.load_table1(C.DEFAULT_TABLE1), plat)
    calls = calls[calls["person_id"].isin(people)].copy()
    dep = depth_table(DEPTH, people, plat)
    pc = per_call(o, calls, trio_of, dep)
    pc.to_csv(intd + "per_call.INTERNAL.tsv", sep="\t", index=False)
    C.log(f"per-call rows {len(pc)}; with both gene depths {pc.dropna(subset=['d_child', 'd_tpar']).shape[0]}")
    ph = phase_pairs(o)
    ph["mhc_child"] = [dep.get((c, "window_MHC"), np.nan) for c in ph["child"]]

    m2d = pd.concat([mendel_2d(pc, m)[0] for m in ("badfull", "bad2")])
    edges = mendel_2d(pc, "badfull")[1]
    C.write_tsv(m2d, out + "mendel_by_depth_2d.tsv")
    C.write_tsv(logit_table(pc, M59), out + "depth_logit.tsv")
    ebb = error_by_bin(pc, ph, None)
    C.write_tsv(ebb, out + "error_by_depth_bin.tsv")

    rng = np.random.default_rng(a.seed)
    prune = []
    per_child = pc.groupby("child").agg(mhc_child=("mhc_child", "first")).reset_index()
    order = dict(zip(per_child["child"], per_child["mhc_child"]))
    # (1) trios pruned by the child's MHC depth: every call / phase pair of a removed child goes
    for metric, units, col in (("mendel_full_by_child_mhc_depth", pc, "badfull"),
                               ("mendel_2field_by_child_mhc_depth", pc, "bad2"),
                               ("phase_pairs_by_child_mhc_depth", ph, "bad")):
        u = units.copy()
        # prune whole children: rank children by MHC depth, then expand to their calls / pairs
        ch = sorted({c for c in u["child"] if np.isfinite(order.get(c, np.nan))}, key=lambda c: order[c])
        rank = {c: i for i, c in enumerate(ch)}
        u["_o"] = u["child"].map(rank)
        prune += prune_curve(u, col, "_o", rng, metric)
    # (2) calls pruned by the lower of child and transmitting-parent gene depth (a per-call filter)
    pcm = pc.assign(_o=np.fmin(pc["d_child"], pc["d_tpar"]))
    prune += prune_curve(pcm, "badfull", "_o", rng, "mendel_full_by_min_gene_depth")
    prune += prune_curve(pcm, "bad2", "_o", rng, "mendel_2field_by_min_gene_depth")
    C.write_tsv(pd.DataFrame(prune), out + "prune_curve.tsv")

    # (3) switch / hamming (57, MHC target) vs child MHC depth and number of pieces
    rows = []
    if os.path.exists(TRIOEVAL):
        te = pd.read_csv(TRIOEVAL, sep="\t", dtype={"child": str})
        for tgt, g in te.groupby("target"):
            g = g.groupby("child").agg(sw_err=("sw_err", "sum"), sw_tot=("sw_tot", "sum"), hm_err=("hm_err", "sum"),
                                       hm_tot=("hm_tot", "sum"), pieces=("n_pieces", "sum")).reset_index()
            g["mhc"] = [dep.get((c, "window_MHC"), np.nan) for c in g["child"]]
            g["sw"] = g["sw_err"] / g["sw_tot"].where(g["sw_tot"] > 0)
            g["hm"] = g["hm_err"] / g["hm_tot"].where(g["hm_tot"] > 0)
            g = g.dropna(subset=["mhc", "sw", "hm"])
            for y in ("sw", "hm"):
                for x in ("mhc", "pieces"):
                    if len(g) >= SMALL:
                        rho = pd.Series(g[x]).rank().corr(pd.Series(g[y]).rank())
                        rows.append(dict(target=tgt, error={"sw": "switch", "hm": "hamming"}[y],
                                         versus={"mhc": "child_mhc_depth", "pieces": "n_pieces"}[x],
                                         n_children=str(len(g)), spearman_rho=round(float(rho), 3)))
            if tgt == "MHC" and len(g) >= 2 * SMALL:
                lo = g["mhc"] <= g["mhc"].median()
                for y in ("sw", "hm"):
                    rows.append(dict(target=tgt, error={"sw": "switch", "hm": "hamming"}[y],
                                     versus="median_split_low_vs_high_depth", n_children=str(len(g)),
                                     low_depth_median_pct=round(100 * float(g.loc[lo, y].median()), 3),
                                     high_depth_median_pct=round(100 * float(g.loc[~lo, y].median()), 3)))
    C.write_tsv(pd.DataFrame(rows), out + "phase_vs_depth.tsv")
    figures(m2d, edges, prune, ebb, out)
    internal_scatter(pc, intd + "fig_68_scatter.INTERNAL.png")
    C.log("DONE 68")
    return 0


if __name__ == "__main__":
    sys.exit(main())
