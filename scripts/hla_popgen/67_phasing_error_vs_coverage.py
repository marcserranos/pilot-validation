#!/usr/bin/env python3
"""S05 67: where do typing and phasing errors happen, and does pruning low-coverage people help?
(VM; run from the scripts dir next to 56, 57, 58, 59 and _s05_common; reads INTERNAL tables only.)

Four questions, each answered with aggregate (PUBLIC) tables:

(A) Relative pairs (56/56b truth): for first-degree pairs that share their MHC (<= 1 classical gene with
    no shared 2-field allele, as in 56b), each gene called in both relatives is one record:
    bad = no shared 2-field allele (typing error or recombination), plus each relative's sequencing
    depth AT THAT GENE (58, MAPQ 0, normalised to the cohort median for that gene x platform).
      67_pair_depth_grid.tsv     error rate per (lower-depth relative, higher-depth relative) cell,
                                 depth tertiles over all records -> the binned version of a two-axis
                                 scatter (a point per pair would be participant-level data, which the
                                 All of Us dissemination policy does not allow outside the Workbench)
      67_pair_error_by_depth.tsv error rate by decile of the LOWER of the two depths (and of the higher)
      67_pair_logit.tsv          logistic OR per 0.1 lower depth (min of the pair), adjusted for gene
                                 and for the call-confidence score (0-7), and the same for the max
(B) Trios (57b truth): each child call with a resolved parent of origin is one record: bad = the
    transmitting parent lacks the child's 2-field allele (Mendelian typing error); depths of child and
    transmitting parent at that gene.
      67_trio_depth_grid.tsv     2 x 2 grid (below / above the median depth for each of the two people)
      67_trio_error_by_depth.tsv error by child-depth tertile
(C) Phasing switches.
      67_switch_vs_depth.tsv     57 yak trioeval per child x target (50 trios): Spearman rho of switch
                                 and hamming error vs the child's depth over the target window and vs
                                 the number of contig pieces; mean error below / above median depth
      67_switch_topology.tsv     57b: along each child contig, consecutive genes (hg38 order) with
                                 resolved but DIFFERENT parents = a switch inside that interval. Per
                                 interval class (within class I, class I to III, class III, class III to
                                 II, within class II, class II to DP): intervals seen, % with a switch;
                                 and the switch rate by child MHC-depth half and by contig count
(D) Pruning (59 per-person outcomes, unrelated set): remove the worst x% of people by a QC metric
    (normalised MHC or KIR window depth, or MHC contig count) and recompute the quality metrics on
    who is left, next to removing the same number of people at random (mean of 20 draws), plus the
    ancestry and platform make-up of the people removed.
      67_prune_curve.tsv

Disclosure: rates are blanked when numerator or denominator is 1-19; record and people counts are
shown rounded to the nearest 10 ('<20' for 1-19), so cells cannot be recovered by subtraction from
other tables. The INTERNAL per-record tables stay on the VM (~/s05/internal/67/).

VM:
    cd ~/repos/pilot-validation-main/scripts/hla_popgen
    $PY 67_phasing_error_vs_coverage.py            # all parts
    $PY 67_phasing_error_vs_coverage.py --parts D  # one part
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


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


OUT = os.path.expanduser("~/s05/results/67/")
INT = os.path.expanduser("~/s05/internal/67/")
DEPTH58 = os.path.expanduser("~/s05/internal/58/per_person_depth.tsv")
PP59 = os.path.expanduser("~/s05/internal/59/per_person_depth_outcome.tsv")
TRIO57 = os.path.expanduser("~/s05/internal/57/per_trio_region.INTERNAL.tsv")
ORIGIN57B = os.path.expanduser("~/s05/internal/57b/gene_origin.INTERNAL.tsv")
REGIONS_BED = os.path.join(HERE, "..", "..", "reports", "hla_popgen", "58_regional_depth", "regions.bed")
TARGET_WINDOW = {"MHC": "window_MHC", "MHC_classI_A_to_B": "window_MHC", "MHC_classII_DR_DQ": "window_MHC",
                 "KIR": "window_KIR", "chr20_control": None}
PRUNE_PCT = [0, 1, 2, 5, 10, 15, 20, 30]   # >= 1% steps, rates to 1 decimal: differencing
# neighbouring rows cannot pin the events among the people removed between two steps
PRUNE_OUTCOMES = ["kir_incomplete", "kir_fw_missing", "hla_absent", "hla_hap_empty", "dq_split",
                  "mhc_fragmented"]
MIN = 20


# --------------------------------------------------------------------------- disclosure helpers
def rnd(n):
    n = int(n)
    if n == 0:
        return "0"
    if n < MIN:
        return f"<{MIN}"
    return str(int(round(n / 10.0)) * 10)


def rate(k, n, digits=1):
    """(pct, lo, hi) strings; blank when k or n is 1-19 or n == 0."""
    k, n = int(k), int(n)
    if n == 0 or 0 < k < MIN or 0 < n < MIN:
        return "", "", ""
    p, lo, hi = C.wilson(k, n)
    return f"{100*p:.{digits}f}", f"{100*lo:.{digits}f}", f"{100*hi:.{digits}f}"


def read_int(path, **kw):
    """INTERNAL tables may start with a '# INTERNAL ...' banner line; skip it."""
    with open(path) as f:
        first = f.readline()
    return pd.read_csv(path, sep="\t", skiprows=1 if first.startswith("#") else 0, **kw)


# --------------------------------------------------------------------------- depth
def region_gene(region):
    r = str(region)
    return r[4:] if r.startswith("HLA-") else r


def gene_depth_index(depth, platform_of, mapq=0):
    """58 per_person_depth -> {(person, gene_or_window): rel_depth / median(rel_depth) over people of the
    same platform}. Windows keep their names (window_MHC, window_KIR)."""
    d = depth[pd.to_numeric(depth["mapq"], errors="coerce") == mapq].copy()
    d["rel_depth"] = pd.to_numeric(d["rel_depth"], errors="coerce")
    d = d[d["rel_depth"].notna()]
    d["plat"] = d["person_id"].map(platform_of).fillna("other")
    med = d.groupby(["region", "plat"])["rel_depth"].transform("median")
    d["norm"] = d["rel_depth"] / med
    return {(p, region_gene(r)): v for p, r, v in zip(d["person_id"], d["region"], d["norm"])}


def logit_or(x, y, covars=None):
    """OR per 0.1 LOWER x (x = normalised depth), logistic, optional dummy covariates (DataFrame)."""
    m59 = sys.modules.get("m59") or _load("m59", "59_depth_outcome_join.py")
    X = [np.ones(len(x)), -10.0 * np.asarray(x, float)]
    if covars is not None:
        for c in covars.columns:
            v = pd.to_numeric(covars[c], errors="coerce")
            if v.notna().all():
                X.append(v.to_numpy(float))
            else:   # categorical: dummies only for levels with >= MIN events, others pooled in the reference
                y_ = pd.Series(np.asarray(y, float), index=covars.index)
                ev = y_.groupby(covars[c].astype(str)).sum()
                for lev in sorted(ev[ev >= MIN].index)[1:]:
                    X.append((covars[c].astype(str) == lev).to_numpy(float))
    X = np.column_stack(X)
    b, s, conv, how = m59.fit_logit(X, np.asarray(y, float))
    if not np.isfinite(b[1]):
        return {"or_per_0.1_lower": "", "lo": "", "hi": "", "p": "", "converged": False}
    z = b[1] / s[1] if s[1] > 0 else np.nan
    p = float(2 * m59._norm_sf(abs(z))) if np.isfinite(z) else np.nan
    return {"or_per_0.1_lower": f"{np.exp(b[1]):.3f}", "lo": f"{np.exp(b[1] - 1.96 * s[1]):.3f}",
            "hi": f"{np.exp(b[1] + 1.96 * s[1]):.3f}", "p": f"{p:.2g}", "converged": bool(conv)}


# --------------------------------------------------------------------------- (A) relative pairs
def pair_records(fd, geno, dep, genes):
    """fd: pairs (i, j, pp); geno: 56 index_calls genotype dict. Returns DataFrame of records for pairs
    that share their MHC (56b filter)."""
    rec = []
    for i, j, pp in zip(fd["i"], fd["j"], fd["pp"]):
        fails = [g for g in genes["G8"] if geno.get((i, g)) and geno.get((j, g))
                 and not set(geno[(i, g)][0]) & set(geno[(j, g)][0])]
        if len(fails) > 1:
            continue
        for g in genes["ALL"]:
            a, b = geno.get((i, g)), geno.get((j, g))
            di, dj = dep.get((i, g)), dep.get((j, g))
            if a and b and di is not None and dj is not None:
                rec.append((i, j, g, pp, int(not set(a[0]) & set(b[0])), min(a[2], b[2]), min(di, dj), max(di, dj)))
    return pd.DataFrame(rec, columns=["i", "j", "gene", "pp", "bad", "score", "dep_lo", "dep_hi"])


def quantile_edges(x, q):
    e = np.unique(np.nanquantile(np.asarray(x, float), np.linspace(0, 1, q + 1)))
    e[0], e[-1] = -np.inf, np.inf
    return e


def depth_grid(r, xcol, ycol, edges, label):
    rows = []
    lab = [f"{'-inf' if not np.isfinite(a) else f'{a:.2f}'}..{'inf' if not np.isfinite(b) else f'{b:.2f}'}"
           for a, b in zip(edges[:-1], edges[1:])]
    bx = pd.cut(r[xcol], edges, labels=False, include_lowest=True)
    by = pd.cut(r[ycol], edges, labels=False, include_lowest=True)
    for ix in range(len(lab)):
        for iy in range(len(lab)):
            d = r[(bx == ix) & (by == iy)]
            if label == "pair" and iy < ix:   # lower <= higher: the other triangle is empty by construction
                continue
            pct, lo, hi = rate(int(d["bad"].sum()), len(d))
            rows.append({"grid": label, "x_bin": ix + 1, "x_range": lab[ix], "y_bin": iy + 1,
                         "y_range": lab[iy], "n_records": rnd(len(d)), "error_pct": pct, "lo": lo, "hi": hi})
    return rows


def by_quantile(r, col, q, label):
    rows = []
    e = quantile_edges(r[col], q)
    b = pd.cut(r[col], e, labels=False, include_lowest=True)
    for k in range(len(e) - 1):
        d = r[b == k]
        pct, lo, hi = rate(int(d["bad"].sum()), len(d))
        rows.append({"axis": label, "bin": k + 1, "depth_median": f"{d[col].median():.2f}" if len(d) else "",
                     "n_records": rnd(len(d)), "error_pct": pct, "lo": lo, "hi": hi})
    return rows


# --------------------------------------------------------------------------- (B) trios
def trio_records(o, calls, trio_of, dep):
    """o: 57b gene_origin INTERNAL; calls: 56-prepared calls of trio people. Records per resolved child call."""
    key = calls[calls["copy_index"] == 1].set_index(["person_id", "hap", "gene"])
    o = o.join(key[["tf", "score"]], on=["child", "hap", "gene"], how="left")
    par = {(p, g): set(d["tf"]) for (p, g), d in calls.groupby(["person_id", "gene"])}
    rec = []
    for c, g, og, tf, sc in zip(o["child"], o["gene"], o["origin"], o["tf"], o["score"]):
        if og not in ("pat", "mat") or pd.isna(tf) or c not in trio_of:
            continue
        pid = trio_of[c][0] if og == "pat" else trio_of[c][1]
        if (pid, g) not in par:
            continue
        dc, dp = dep.get((c, g)), dep.get((pid, g))
        if dc is None or dp is None:
            continue
        rec.append((c, pid, g, int(sc), int(tf not in par[(pid, g)]), dc, dp))
    return pd.DataFrame(rec, columns=["child", "parent", "gene", "score", "bad", "dep_child", "dep_parent"])


# --------------------------------------------------------------------------- (C) switches
def spearman(x, y):
    x, y = pd.Series(x, dtype=float), pd.Series(y, dtype=float)
    ok = x.notna() & y.notna()
    if ok.sum() < MIN:
        return "", ""
    rho = x[ok].rank().corr(y[ok].rank())
    n = int(ok.sum())
    t = rho * np.sqrt((n - 2) / max(1e-12, 1 - rho * rho))
    m59 = sys.modules.get("m59") or _load("m59", "59_depth_outcome_join.py")
    return f"{rho:.3f}", f"{float(2 * m59._norm_sf(abs(t))):.2g}"


def switch_vs_depth(tr, dep):
    rows = []
    tr = tr.copy()
    for c in ("sw_err", "sw_tot", "hm_err", "hm_tot", "n_pieces"):
        tr[c] = pd.to_numeric(tr[c], errors="coerce")
    tr["sw"] = tr["sw_err"] / tr["sw_tot"].replace(0, np.nan)
    tr["hm"] = tr["hm_err"] / tr["hm_tot"].replace(0, np.nan)
    for target, d in tr.groupby("target"):
        win = TARGET_WINDOW.get(target)
        # one value per child: mean over its two haplotypes
        g = d.groupby("child").agg(sw=("sw", "mean"), hm=("hm", "mean"), pieces=("n_pieces", "sum")).reset_index()
        g["depth"] = [dep.get((c, win)) if win else np.nan for c in g["child"]]
        for metric in ("sw", "hm"):
            rho_d, p_d = spearman(g["depth"], g[metric])
            rho_p, p_p = spearman(g["pieces"], g[metric])
            row = {"target": target, "metric": "switch_error" if metric == "sw" else "hamming_error",
                   "n_children": rnd(len(g)), "rho_vs_depth": rho_d, "p_vs_depth": p_d,
                   "rho_vs_contig_pieces": rho_p, "p_vs_contig_pieces": p_p}
            dd = g[g["depth"].notna()]
            if win and len(dd) >= 2 * MIN:
                med = dd["depth"].median()
                lo_, hi_ = dd[dd["depth"] < med], dd[dd["depth"] >= med]
                row.update(mean_pct_low_depth_half=f"{100*lo_[metric].mean():.2f}",
                           mean_pct_high_depth_half=f"{100*hi_[metric].mean():.2f}",
                           median_pct_low_depth_half=f"{100*lo_[metric].median():.2f}",
                           median_pct_high_depth_half=f"{100*hi_[metric].median():.2f}")
            rows.append(row)
    return pd.DataFrame(rows)


def gene_hg38_mid(bed_path):
    pos = {}
    for line in open(bed_path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4 and f[0] == "chr6":
            pos[region_gene(f[3])] = (int(f[1]) + int(f[2])) / 2
    return pos


# interval classes along chr6 (hg38): class I up to ~31.4 Mb (A..B, C), class III ~31.4-32.4 Mb, class II
# DR/DQ ~32.4-32.8 Mb, DP ~33.0-33.1 Mb
def mhc_zone(x):
    if x < 31_400_000:
        return "classI"
    if x < 32_400_000:
        return "classIII"
    if x < 32_950_000:
        return "classII_DRDQ"
    return "classII_DP"


def switch_intervals(o, calls, pos):
    """Consecutive resolved genes along each child contig (ordered by gene_start on the contig).
    Returns interval records: child, zone_a, zone_b, switch (0/1), contig gene count."""
    c1 = calls[calls["copy_index"] == 1][["person_id", "hap", "gene", "contig", "gene_start"]]
    c1 = c1.assign(ctg=c1["contig"].astype(str).str.partition(":")[0])
    m = o[o["origin"].isin(["pat", "mat"])].merge(
        c1, left_on=["child", "hap", "gene", "contig"], right_on=["person_id", "hap", "gene", "ctg"], how="inner")
    m["gene_start"] = pd.to_numeric(m["gene_start"], errors="coerce")
    m = m[m["gene_start"].notna() & m["gene"].isin(pos)]
    rec = []
    for (c, h, ctg), d in m.groupby(["child", "hap", "ctg"]):
        d = d.sort_values("gene_start")
        g, og = list(d["gene"]), list(d["origin"])
        for k in range(len(g) - 1):
            za, zb = sorted([mhc_zone(pos[g[k]]), mhc_zone(pos[g[k + 1]])],
                            key=["classI", "classIII", "classII_DRDQ", "classII_DP"].index)
            rec.append((c, h, za if za == zb else f"{za}->{zb}", int(og[k] != og[k + 1]), len(g)))
    return pd.DataFrame(rec, columns=["child", "hap", "interval", "switch", "genes_on_contig"])


def topology_table(iv, dep, nctg_of):
    rows = []
    for name, d in [("ALL", iv)] + list(iv.groupby("interval")):
        pct, lo, hi = rate(int(d["switch"].sum()), len(d))
        rows.append({"view": "interval_class", "level": name, "n_intervals": rnd(len(d)),
                     "switch_pct": pct, "lo": lo, "hi": hi})
    iv = iv.copy()
    iv["depth"] = [dep.get((c, "window_MHC")) for c in iv["child"]]
    dd = iv[iv["depth"].notna()]
    if len(dd):
        med = dd.groupby("child")["depth"].first().median()
        for lev, d in (("below_median_MHC_depth", dd[dd["depth"] < med]),
                       ("at_or_above_median_MHC_depth", dd[dd["depth"] >= med])):
            pct, lo, hi = rate(int(d["switch"].sum()), len(d))
            rows.append({"view": "child_depth_half", "level": lev, "n_intervals": rnd(len(d)),
                         "switch_pct": pct, "lo": lo, "hi": hi})
    iv["nctg"] = [nctg_of.get((c, h)) for c, h in zip(iv["child"], iv["hap"])]
    for lev, d in (("1-2 MHC contigs on the haplotype", iv[iv["nctg"] <= 2]),
                   ("3-5", iv[(iv["nctg"] >= 3) & (iv["nctg"] <= 5)]),
                   ("6+", iv[iv["nctg"] >= 6])):
        pct, lo, hi = rate(int(d["switch"].sum()), len(d))
        rows.append({"view": "haplotype_contig_count", "level": lev, "n_intervals": rnd(len(d)),
                     "switch_pct": pct, "lo": lo, "hi": hi})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- (D) pruning
def prune_curve(pp, metrics, pair_rec=None, trio_rec=None, everyone=None, draws=20, seed=67):
    """pp: per-person UNRELATED table with outcomes, ancestry, platform and the metric columns. The worst
    x% of pp is removed (ties broken by order). Relatives are by definition not in the unrelated set, so
    the pair and trio metrics use the same metric cut value applied to everyone (`everyone`: per-person
    table of all people with the metric columns): a record is kept when every person in it passes."""
    rng = np.random.default_rng(seed)
    rows = []
    n_all = len(pp)
    for mname, col, worse_high in metrics:
        v = pd.to_numeric(pp[col], errors="coerce")
        ok = v.notna()
        order = (v[ok].sort_values(ascending=not worse_high)).index  # worst first
        for x in PRUNE_PCT:
            k = int(round(n_all * x / 100.0))
            removed = set(order[:k])
            kept = pp.loc[~pp.index.isin(removed)]
            rem = pp.loc[pp.index.isin(removed)]
            row = {"metric": mname, "prune_pct": x, "n_removed": rnd(len(rem)), "n_kept": rnd(len(kept))}
            for o in PRUNE_OUTCOMES:
                if o not in pp:
                    continue
                y = pd.to_numeric(kept[o], errors="coerce").dropna()
                row[f"{o}_pct"] = rate(int(y.sum()), len(y), 1)[0]
                # random null: same number removed at random
                if k:
                    vals = []
                    for _ in range(draws):
                        idx = rng.choice(pp.index.to_numpy(), size=k, replace=False)
                        yy = pd.to_numeric(pp.loc[~pp.index.isin(set(idx)), o], errors="coerce").dropna()
                        vals.append(100.0 * yy.mean())
                    row[f"{o}_pct_random"] = f"{np.mean(vals):.1f}"
                else:
                    row[f"{o}_pct_random"] = row[f"{o}_pct"]
            for a in C.ANCESTRY_ORDER:
                # whole percents of a rounded n_removed only bound the count to a few people, so also blank shares
                # under 30 people (a margin over MIN) to keep a back-calculated count from landing below 20
                ka = int((rem["ancestry"] == a).sum())
                row[f"removed_pct_{a}"] = rate(ka, len(rem), 0)[0] if len(rem) and ka >= MIN + 10 else ""
            row["removed_pct_sequel2e"] = (rate(int((rem["platform"] == "sequel2e").sum()), len(rem), 0)[0]
                                           if len(rem) else "")
            if everyone is not None and col in everyone:
                ev = pd.to_numeric(everyone[col], errors="coerce")
                if k:
                    cut = v.loc[order[k - 1]]
                    bad_ids = set(everyone.loc[(ev <= cut) if not worse_high else (ev >= cut), "person_id"])
                else:
                    bad_ids = set()
                if pair_rec is not None and len(pair_rec):
                    pr = pair_rec[~pair_rec["i"].isin(bad_ids) & ~pair_rec["j"].isin(bad_ids)]
                    row["pair_records_kept"] = rnd(len(pr))
                    row["pair_error_pct"] = rate(int(pr["bad"].sum()), len(pr), 1)[0]
                if trio_rec is not None and len(trio_rec):
                    tr = trio_rec[~trio_rec["child"].isin(bad_ids) & ~trio_rec["parent"].isin(bad_ids)]
                    row["trio_calls_kept"] = rnd(len(tr))
                    row["trio_mendel_error_pct"] = rate(int(tr["bad"].sum()), len(tr), 1)[0]
            rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parts", default="ABCD")
    ap.add_argument("--out-dir", default=OUT)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    os.makedirs(INT, exist_ok=True)
    w = lambda df, name: (C.write_tsv(df, os.path.join(a.out_dir, name)),  # noqa: E731
                          C.log(f"wrote {name} ({len(df)} rows)"))
    B = _load("b56", "56_phasing_typing_benchmark.py")
    cohort = C.load_cohort(C.DEFAULT_COHORT)
    plat = dict(zip(cohort["person_id"], cohort["platform"]))
    anc = dict(zip(cohort["person_id"], cohort["ancestry_pred"]))
    C.log("loading 58 depth")
    dep = gene_depth_index(read_int(DEPTH58, dtype={"person_id": str}), plat)
    C.TABLE1_USECOLS = C.TABLE1_USECOLS + ["n_tied"]
    t1 = C.load_table1(C.DEFAULT_TABLE1)
    calls = B.prepare(t1, plat)
    genes = {"G8": B.G8, "ALL": B.ALLG}

    pair_rec = trio_rec = None
    if set("AD") & set(a.parts):
        pairs = B.build_pairs(B.load_rel(C.DEFAULT_RELATEDNESS), set(t1["person_id"]), plat)
        fd = pairs[pairs["cls"] != "twin"]
        geno, _hp = B.index_calls(calls[calls["person_id"].isin(set(fd["i"]) | set(fd["j"]))])
        pair_rec = pair_records(fd, geno, dep, genes)
        pair_rec.to_csv(INT + "pair_records.INTERNAL.tsv", sep="\t", index=False)
        C.log(f"pair records {len(pair_rec)} from {pair_rec[['i', 'j']].drop_duplicates().shape[0]} pairs")
    if "A" in a.parts and len(pair_rec):
        e = quantile_edges(pd.concat([pair_rec["dep_lo"], pair_rec["dep_hi"]]), 3)
        w(pd.DataFrame(depth_grid(pair_rec, "dep_lo", "dep_hi", e, "pair")), "67_pair_depth_grid.tsv")
        w(pd.DataFrame(by_quantile(pair_rec, "dep_lo", 10, "lower_depth_of_pair")
                       + by_quantile(pair_rec, "dep_hi", 10, "higher_depth_of_pair")), "67_pair_error_by_depth.tsv")
        rows = []
        for xcol in ("dep_lo", "dep_hi"):
            for name, d in (("all_scores", pair_rec), ("score_ge_6", pair_rec[pair_rec["score"] >= 6])):
                r = logit_or(d[xcol].clip(upper=2.0), d["bad"], d[["score"]])
                rows.append({"depth": xcol, "subset": name, "n_records": rnd(len(d)), **r})
        w(pd.DataFrame(rows), "67_pair_logit.tsv")

    o = None
    if set("BCD") & set(a.parts) and os.path.exists(ORIGIN57B):
        o = read_int(ORIGIN57B, dtype={"child": str, "contig": str})
        T = _load("t57", "57_yak_trioeval_pilot.py")
        lr = set(cohort.loc[cohort["in_lr"].astype(str).str.lower() == "true", "person_id"])
        kids = set(o["child"])
        trio_of = {t[0]: (t[1], t[2]) for t in T.find_trios(lr) if t[0] in kids}
        trio_rec = trio_records(o, calls[calls["person_id"].isin(set(kids) | {x for v in trio_of.values() for x in v})],
                                trio_of, dep)
        C.log(f"trio records {len(trio_rec)}")
    if "B" in a.parts and trio_rec is not None and len(trio_rec):
        rows = []
        mc, mp = trio_rec["dep_child"].median(), trio_rec["dep_parent"].median()
        for xi, xl in ((trio_rec["dep_child"] < mc, "child_below_median"), (trio_rec["dep_child"] >= mc, "child_at_or_above")):
            for yi, yl in ((trio_rec["dep_parent"] < mp, "parent_below_median"),
                           (trio_rec["dep_parent"] >= mp, "parent_at_or_above")):
                d = trio_rec[xi & yi]
                pct, lo, hi = rate(int(d["bad"].sum()), len(d))
                rows.append({"child_depth": xl, "parent_depth": yl, "n_calls": rnd(len(d)),
                             "mendel_error_pct": pct, "lo": lo, "hi": hi})
        w(pd.DataFrame(rows), "67_trio_depth_grid.tsv")
        w(pd.DataFrame(by_quantile(trio_rec, "dep_child", 3, "child_depth")
                       + by_quantile(trio_rec, "dep_parent", 3, "transmitting_parent_depth")),
          "67_trio_error_by_depth.tsv")

    if "C" in a.parts:
        if os.path.exists(TRIO57):
            tr = read_int(TRIO57, dtype={"child": str})
            w(switch_vs_depth(tr, dep), "67_switch_vs_depth.tsv")
        if o is not None:
            pos = gene_hg38_mid(REGIONS_BED)
            iv = switch_intervals(o, calls, pos)
            nctg = calls.groupby(["person_id", "hap"])["contig"].nunique().to_dict()
            w(topology_table(iv, dep, nctg), "67_switch_topology.tsv")

    if "D" in a.parts:
        everyone = read_int(PP59, dtype={"person_id": str})
        pp = everyone
        keep = C.greedy_unrelated(pp["person_id"], C.load_relatedness_pairs(C.DEFAULT_RELATEDNESS))
        pp = pp[pp["person_id"].isin(keep)].reset_index(drop=True)
        if "ancestry" not in pp:
            pp["ancestry"] = pp["person_id"].map(anc)
        if "platform" not in pp:
            pp["platform"] = pp["person_id"].map(plat)
        metrics = [("MHC window depth (normalised)", "norm_window_MHC_q0", False),
                   ("KIR window depth (normalised)", "norm_window_KIR_q0", False),
                   ("MHC contigs per haplotype (mean)", "mhc_contigs_mean", True)]
        metrics = [m for m in metrics if m[1] in pp]
        w(prune_curve(pp, metrics, pair_rec, trio_rec, everyone), "67_prune_curve.tsv")
    C.log("DONE 67")
    return 0


if __name__ == "__main__":
    sys.exit(main())
