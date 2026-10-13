#!/usr/bin/env python3
"""S05 57d: do the PARENTS' long-read data drive the trio error rates? (VM; run from the scripts dir next to
57, 57b, 58, 67 and _s05_common; reads INTERNAL tables and the long-read manifest only, no BAM is opened.)

Why (Marc, report comment 2026-10-13): 57's figure splits trios by the child's platform only, but yak's
truth is built from the parents' reads. A parent with thin or older (Sequel II) data misses some of the
k-mers he or she carries; those k-mers then look specific to the other parent, so a correct child contig
can be scored as a switch or hamming error. The error can come from either side.

What it does, per trio (the same trios 57 and 57b scored):
  - parent platform = the platform of the HiFi BAM 57 actually used (57.manifest_index), child platform
    = the platform of the child's assembly;
  - depth = absolute mean depth (MAPQ >= 0) from 58 over the window the k-mer databases came from:
    window_MHC (MHC targets), window_KIR (KIR), the person's genome-wide control depth (chr20);
  - parents' mix: both_revio / both_sequel2e / mixed_revio_sequel2e / any_sequel2.

Outputs (aggregate PUBLIC, ~/s05/results/57d/):
  57d_trio_platform_mix.tsv     trios per child platform x parents' mix (counts masked 1-19)
  57d_parent_depth.tsv          parents per platform, median MHC depth
  57d_region_error_vs_parents.tsv  per 57 target x metric: Spearman of per-trio error with the thinner
                                parent's depth (all trios and within child platform), and a regression
                                log10(error + 0.01%) ~ child Sequel IIe + number of non-Revio parents
                                + log2(thinner parent depth) + log2(child depth), with bootstrap 95% CIs
  57d_origin_vs_parents.tsv     57b gene-level: ambiguous-origin share by tertile of the thinner parent's
                                depth, and whether resolved calls lean towards the deeper parent
INTERNAL: ~/s05/internal/57d/per_trio.INTERNAL.tsv (one row per trio, never leaves the VM).
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

import _s05_common as C

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


M67 = _load("m67", "67_phasing_error_vs_coverage.py")
rnd, rate, read_int, spearman = M67.rnd, M67.rate, M67.read_int, M67.spearman
MIN = M67.MIN
OUT = os.path.expanduser("~/s05/results/57d/")
INT = os.path.expanduser("~/s05/internal/57d/")
WINDOW = {"MHC": "window_MHC", "MHC_classI_A_to_B": "window_MHC", "MHC_classII_DR_DQ": "window_MHC",
          "KIR": "window_KIR", "chr20_control": "control"}
EPS = 1e-4          # 0.01%: lets a zero-error trio enter the log model
N_BOOT = 2000
DEEPER_RATIO = 1.25  # deeper-parent lean test only where one parent has >= 25% more depth


def parents_mix(pa, pb):
    s = {pa, pb}
    if "sequel2" in s:
        return "any_sequel2"
    if s == {"revio"}:
        return "both_revio"
    if s == {"sequel2e"}:
        return "both_sequel2e"
    return "mixed_revio_sequel2e"


def depth_index(depth):
    """58 per_person_depth (MAPQ 0) -> {(person, 'window_MHC'|'window_KIR'|'control'): absolute mean depth}."""
    d = depth[pd.to_numeric(depth["mapq"], errors="coerce") == 0]
    out = {}
    for p, r, md, cd in zip(d["person_id"], d["region"], d["mean_depth"], d["control_depth"]):
        if r in ("window_MHC", "window_KIR"):
            out[(p, r)] = float(md)
        out[(p, "control")] = float(cd)
    return out


def per_trio_rates(tr):
    """57 checkpoint -> one row per (child, target): mean over the two haplotypes, as in 67."""
    tr = tr.copy()
    for c in ("sw_err", "sw_tot", "hm_err", "hm_tot"):
        tr[c] = pd.to_numeric(tr[c], errors="coerce")
    tr["sw"] = tr["sw_err"] / tr["sw_tot"].replace(0, np.nan)
    tr["hm"] = tr["hm_err"] / tr["hm_tot"].replace(0, np.nan)
    return tr.groupby(["child", "target"]).agg(child_platform=("child_platform", "first"),
                                               sw=("sw", "mean"), hm=("hm", "mean")).reset_index()


def ols_boot(y, X, names, seed=57):
    """OLS with bootstrap (over rows) 95% CIs. Returns {name: (beta, lo, hi)}; empty if too few rows."""
    y, X = np.asarray(y, float), np.asarray(X, float)
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    y, X = y[ok], X[ok]
    if len(y) < 2 * MIN:
        return {}
    X1 = np.column_stack([np.ones(len(y)), X])
    beta = np.linalg.lstsq(X1, y, rcond=None)[0]
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(y), len(y))
        if np.linalg.matrix_rank(X1[i]) < X1.shape[1]:
            continue
        bs.append(np.linalg.lstsq(X1[i], y[i], rcond=None)[0])
    bs = np.array(bs)
    lo, hi = np.percentile(bs, 2.5, axis=0), np.percentile(bs, 97.5, axis=0)
    return {n: (beta[k + 1], lo[k + 1], hi[k + 1]) for k, n in enumerate(names)}


def region_table(rates, trio, dep):
    """rates: per_trio_rates; trio: per-trio frame (child, pa, pb, pa_plat, pb_plat, child_platform)."""
    rows = []
    m = rates.merge(trio[["child", "pa", "pb", "n_nonrevio_parents", "parents_mix"]], on="child", how="inner")
    for target, d in m.groupby("target"):
        win = WINDOW[target]
        d = d.copy()
        d["dep_child"] = [dep.get((c, win), np.nan) for c in d["child"]]
        d["dep_parent_min"] = [np.fmin(dep.get((a, win), np.nan), dep.get((b, win), np.nan))
                               for a, b in zip(d["pa"], d["pb"])]
        for metric in ("sw", "hm"):
            row = {"target": target, "metric": "switch_error" if metric == "sw" else "hamming_error",
                   "n_trios": rnd(len(d))}
            row["rho_vs_thinner_parent_depth"], row["p_vs_thinner_parent_depth"] = spearman(d["dep_parent_min"], d[metric])
            row["rho_vs_child_depth"], row["p_vs_child_depth"] = spearman(d["dep_child"], d[metric])
            for pl in ("revio", "sequel2e"):
                s = d[d["child_platform"] == pl]
                row[f"rho_parent_depth_within_{pl}_child"], row[f"p_parent_depth_within_{pl}_child"] = \
                    spearman(s["dep_parent_min"], s[metric])
            for mix in ("both_revio", "mixed_revio_sequel2e", "both_sequel2e", "any_sequel2"):
                s = d[d["parents_mix"] == mix][metric].dropna()
                row[f"median_pct_parents_{mix}"] = f"{100 * s.median():.3f}" if len(s) >= MIN else ""
            X = np.column_stack([(d["child_platform"] == "sequel2e").astype(float), d["n_nonrevio_parents"],
                                 np.log2(d["dep_parent_min"]), np.log2(d["dep_child"])])
            fit = ols_boot(np.log10(d[metric] + EPS), X,
                           ["child_sequel2e", "per_nonrevio_parent", "per_doubling_thinner_parent_depth",
                            "per_doubling_child_depth"])
            for n, (b, lo, hi) in fit.items():
                row[f"beta_{n}"], row[f"lo_{n}"], row[f"hi_{n}"] = f"{b:.3f}", f"{lo:.3f}", f"{hi:.3f}"
            rows.append(row)
    return pd.DataFrame(rows)


def origin_table(o, trio, dep):
    """57b per-call origins -> ambiguous share by thinner-parent depth tertile; lean towards the deeper parent.
    In 57b 'pat' is the trio's first parent (pa) and 'mat' the second (pb)."""
    win = "window_MHC"
    t = trio.copy()
    t["da"] = [dep.get((p, win), np.nan) for p in t["pa"]]
    t["db"] = [dep.get((p, win), np.nan) for p in t["pb"]]
    t["dmin"] = t[["da", "db"]].min(axis=1)
    o = o.merge(t[["child", "da", "db", "dmin", "parents_mix"]], on="child", how="inner")
    o["amb"] = (o["origin"] == "amb").astype(int)
    o["k_inf"] = pd.to_numeric(o["pat_k"], errors="coerce") + pd.to_numeric(o["mat_k"], errors="coerce")
    rows = []
    per = o.groupby("child").agg(amb=("amb", "mean"), dmin=("dmin", "first"), k=("k_inf", "median")).reset_index()
    rho, p = spearman(per["dmin"], per["amb"])
    rho_k, p_k = spearman(per["dmin"], per["k"])
    rows.append({"test": "spearman_per_trio", "group": "ambiguous_share_vs_thinner_parent_depth",
                 "n_trios": rnd(len(per)), "value": rho, "p": p})
    rows.append({"test": "spearman_per_trio", "group": "informative_kmers_vs_thinner_parent_depth",
                 "n_trios": rnd(len(per)), "value": rho_k, "p": p_k})
    per = per[per["dmin"].notna()]
    if len(per) >= 3 * MIN:
        edges = np.quantile(per["dmin"], [1 / 3, 2 / 3])
        per["tert"] = np.digitize(per["dmin"], edges) + 1
        tert = dict(zip(per["child"], per["tert"]))
        o["tert"] = o["child"].map(tert)
        for k in (1, 2, 3):
            d = o[o["tert"] == k]
            pct, lo, hi = rate(int(d["amb"].sum()), len(d))
            rows.append({"test": "ambiguous_by_thinner_parent_depth_tertile", "group": f"tertile_{k}",
                         "n_trios": rnd(d["child"].nunique()), "n_calls": rnd(len(d)),
                         "depth_median": f"{per.loc[per['tert'] == k, 'dmin'].median():.1f}",
                         "value": pct, "lo": lo, "hi": hi,
                         "informative_kmers_median": f"{d['k_inf'].median():.0f}"})
    for mix, d in o.groupby("parents_mix"):
        if d["child"].nunique() < MIN:
            continue
        pct, lo, hi = rate(int(d["amb"].sum()), len(d))
        rows.append({"test": "ambiguous_by_parents_mix", "group": mix, "n_trios": rnd(d["child"].nunique()),
                     "n_calls": rnd(len(d)), "value": pct, "lo": lo, "hi": hi})
    # lean: among resolved calls in trios with unequal parent depth, share assigned to the deeper parent
    r = o[(o["origin"] != "amb") & o["da"].notna() & o["db"].notna()].copy()
    r = r[np.maximum(r["da"], r["db"]) >= DEEPER_RATIO * np.minimum(r["da"], r["db"])]
    deeper = np.where(r["da"] > r["db"], "pat", "mat")
    k = int((r["origin"].values == deeper).sum())
    pct, lo, hi = rate(k, len(r))
    rows.append({"test": "resolved_calls_assigned_to_deeper_parent", "group": f"depth_ratio_ge_{DEEPER_RATIO}",
                 "n_trios": rnd(r["child"].nunique()), "n_calls": rnd(len(r)), "value": pct, "lo": lo, "hi": hi})
    return pd.DataFrame(rows)


def build_trios(T, lr, kids57, kids57b, bam, asm):
    """57's trio per child = first trio with assembly + both parent BAMs; 57b's = first trio found."""
    first57, first57b = {}, {}
    for t in T.find_trios(lr):
        c = t[0]
        if c in kids57b and c not in first57b:
            first57b[c] = t
        if c in kids57 and c not in first57 and c in asm and t[1] in bam and t[2] in bam:
            first57[c] = t
    rows = []
    for name, d in (("57", first57), ("57b", first57b)):
        for c, (_, a, b) in d.items():
            if a not in bam or b not in bam or c not in asm:
                continue
            pa, pb = bam[a]["plat"], bam[b]["plat"]
            rows.append({"set": name, "child": c, "pa": a, "pb": b, "child_platform": asm[c]["plat"],
                         "pa_plat": pa, "pb_plat": pb, "parents_mix": parents_mix(pa, pb),
                         "n_nonrevio_parents": int(pa != "revio") + int(pb != "revio")})
    return pd.DataFrame(rows)


def mix_table(trio):
    rows = []
    for (s, cp), d in trio.groupby(["set", "child_platform"]):
        for mix in ("both_revio", "mixed_revio_sequel2e", "both_sequel2e", "any_sequel2"):
            rows.append({"set": s, "child_platform": cp, "parents_mix": mix,
                         "n_trios": C.mask_n(int((d["parents_mix"] == mix).sum()))})
        same = int(((d["pa_plat"] == cp) & (d["pb_plat"] == cp)).sum())
        q = C.rate_cells(len(d) - same, len(d))  # report the not-matching count, so a small one is masked
        rows.append({"set": s, "child_platform": cp, "parents_mix": "not_both_same_as_child",
                     "n_trios": q["n"], "pct": q["pct"], "lo": q["lo"], "hi": q["hi"]})
    return pd.DataFrame(rows)


def parent_depth_table(trio, dep):
    p = pd.concat([trio[["pa", "pa_plat"]].set_axis(["person", "plat"], axis=1),
                   trio[["pb", "pb_plat"]].set_axis(["person", "plat"], axis=1)]).drop_duplicates("person")
    p["mhc"] = [dep.get((x, "window_MHC"), np.nan) for x in p["person"]]
    p["ctrl"] = [dep.get((x, "control"), np.nan) for x in p["person"]]
    rows = []
    for pl, d in p.groupby("plat"):
        small = len(d) < MIN
        rows.append({"parent_platform": pl, "n_parents": C.mask_n(len(d)),
                     "mhc_depth_median": "" if small else f"{d['mhc'].median():.1f}",
                     "mhc_depth_p25": "" if small else f"{d['mhc'].quantile(.25):.1f}",
                     "mhc_depth_p75": "" if small else f"{d['mhc'].quantile(.75):.1f}",
                     "genome_depth_median": "" if small else f"{d['ctrl'].median():.1f}"})
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=OUT)
    a = ap.parse_args(argv)
    T = _load("t57", "57_yak_trioeval_pilot.py")
    if not os.path.isdir(T.MOUNT + "v9"):
        sys.exit("FATAL: AoU bucket mount missing")
    os.makedirs(a.out_dir, exist_ok=True)
    os.makedirs(INT, exist_ok=True)
    tr = read_int(M67.TRIO57, dtype={"child": str})
    o = read_int(M67.ORIGIN57B, dtype={"child": str, "contig": str})
    dep = depth_index(read_int(M67.DEPTH58, dtype={"person_id": str}))
    cohort = C.load_cohort(C.DEFAULT_COHORT)
    lr = set(cohort.loc[cohort["in_lr"].astype(str).str.lower() == "true", "person_id"])
    kids57, kids57b = set(tr["child"]), set(o["child"])
    cand = [t for t in T.find_trios(lr) if t[0] in kids57 | kids57b]
    bam, asm = T.manifest_index({x for t in cand for x in t})
    trio = build_trios(T, lr, kids57, kids57b, bam, asm)
    C.log(f"trios: 57 {int((trio['set'] == '57').sum())} of {len(kids57)}, "
          f"57b {int((trio['set'] == '57b').sum())} of {len(kids57b)}")
    trio.to_csv(INT + "per_trio.INTERNAL.tsv", sep="\t", index=False)
    w = lambda df, name: C.write_tsv(df, os.path.join(a.out_dir, name))  # noqa: E731
    w(mix_table(trio), "57d_trio_platform_mix.tsv")
    w(parent_depth_table(trio, dep), "57d_parent_depth.tsv")
    w(region_table(per_trio_rates(tr), trio[trio["set"] == "57"], dep), "57d_region_error_vs_parents.tsv")
    w(origin_table(o, trio[trio["set"] == "57b"], dep), "57d_origin_vs_parents.tsv")
    C.log("DONE 57d")


if __name__ == "__main__":
    main()
