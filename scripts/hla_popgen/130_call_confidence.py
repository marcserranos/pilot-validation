#!/usr/bin/env python3
"""S12 130: how reliable is a called HLA allele, and does a per-call number tell us when it is not (VM).

Framing (Marc, 2026-10-30): a missing call is visible in the call matrix, so confidence only matters for calls
that exist. The error is therefore measured on called alleles only: trio child calls whose genotype and whose
transmitting parent are fully called (allele in neither parent = error) and first-degree pairs where both people are
fully called (no shared 2-field allele = error). Calls fall into three cohorts:
  C1  template distance 0 (the assembly equals the named reference allele)
  C2  template distance > 0 and not one of the S10 novel alleles (collapsed to the nearest reference; the noisy middle)
  C3  carrier haplotypes of the S10 novel-ladder selection (SR level >= moderate, >= 1 unrelated carrier, 243 alleles)
Experiments (all aggregate, no person ids in the outputs)
  E1 cohort sizes over every called haplotype            cohort_sizes
  E2 error rate per cohort with Wilson intervals         error_by_cohort
  E3 error rate by template distance / CDS depth bins    error_by_metric_bin (all records, C1 only, C2 only)
  E4 pruning curves and fixed prune policies             prune_curves, prune_policies
  E5 does the number predict the error: AUC + bootstrap  predictor_auc
Reads (VM): ~/s06/internal/71/call_table.INTERNAL.tsv, ~/s06/internal/76/units_B_all.INTERNAL.tsv,
            ~/s06/internal/76/units_pair_complete.INTERNAL.tsv, ~/s10b/ladder_data.private.json, pool_carriers.private.tsv
Writes <out>/<name>.private.tsv (exact) and <name>.pub.tsv (counts 1-19 masked; rates dropped on masked rows).

VM:  python3 scripts/hla_popgen/130_call_confidence.py --out ~/s13/130
"""
import argparse
import json
import os
import re

import numpy as np
import pandas as pd

CLASSICAL8 = {"A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"}
SRL_MIN = 3          # short-read level moderate or better
THRESH = 20
N_PERM = 100


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def auc(y, s):
    y = np.asarray(y, bool)
    s = np.asarray(s, float)
    n1, n0 = y.sum(), (~y).sum()
    if n1 == 0 or n0 == 0:
        return np.nan
    r = pd.Series(s).rank(method="average").to_numpy()
    return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def boot_auc(y, s, groups, n_boot=500, seed=0):
    rng = np.random.default_rng(seed)
    y, s = np.asarray(y, bool), np.asarray(s, float)
    idx_by = pd.Series(np.arange(len(y))).groupby(np.asarray(groups)).apply(np.array).tolist()
    out = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(idx_by), len(idx_by))
        ii = np.concatenate([idx_by[i] for i in pick])
        a = auc(y[ii], s[ii])
        if not np.isnan(a):
            out.append(a)
    return (np.percentile(out, 2.5), np.percentile(out, 97.5)) if len(out) > 20 else (np.nan, np.nan)


def gname(g):
    return re.sub(r"^HLA-", "", str(g))


def hnum(h):
    return re.sub(r"\D", "", str(h))


def load_carriers(ladder_json, carriers_tsv):
    """Haplotype keys of the S10 selection: SR level >= moderate and >= 1 unrelated carrier (xu, or xuc lower bound)."""
    L = json.load(open(ladder_json))
    xu = lambda a: a["xu"] if a["xu"] is not None else (1 if a["xuc"] == 0 else 2)
    ids = {a["id"] for a in L["alleles"] if a["srl"] >= SRL_MIN and xu(a) >= 1}
    P = pd.read_csv(carriers_tsv, sep="\t")
    P = P[P.allele_id.isin(ids)]
    hap = {f"{p}|{gname(g)}|{hnum(h)}" for p, g, h in zip(P.person.astype(str), P.gene, P.hap)}
    pg = {f"{p}|{gname(g)}" for p, g in zip(P.person.astype(str), P.gene)}
    return len(ids), hap, pg


def cohort(dist, is_c3):
    return np.where(is_c3, "C3", np.where(np.nan_to_num(dist, nan=0.0) == 0, "C1", "C2"))


def edges(call_depth):
    pos = call_depth[call_depth > 0]
    return np.quantile(pos, [1 / 3, 2 / 3])


def bin_metric(metric, v, depth_edges):
    v = np.nan_to_num(np.asarray(v, float), nan=0.0)
    if metric == "template_distance":
        return pd.cut(v, [-1, 0, 1, 2, 5, 1e9], labels=["0", "1", "2", "3-5", "6+"]).astype(str)
    if metric == "n_aa_changes":
        return pd.cut(v, [-1, 0, 1, 3, 1e9], labels=["0", "1", "2-3", "4+"]).astype(str)
    e1, e2 = depth_edges
    return pd.cut(v, [-1, 0, e1, e2, 2], labels=["zero", "low>0", "mid>0", "high>0"]).astype(str)


def rate_row(k, n, **kw):
    lo, hi = wilson(k, n)
    return dict(kw, n=int(n), errors=int(k), rate=k / n if n else np.nan, lo=lo, hi=hi)


def curve(y, s, seed=0):
    """Expected share of errors removed when pruning the worst fraction f of records (random tie-break)."""
    rng = np.random.default_rng(seed)
    y, s = np.asarray(y, float), np.asarray(s, float)
    n = len(y)
    fr = np.round(np.arange(0, 0.3001, 0.005), 4)
    acc = np.zeros(len(fr))
    for _ in range(N_PERM):
        o = np.lexsort((rng.random(n), -s))
        cum = np.concatenate([[0], np.cumsum(y[o])])
        acc += cum[np.round(fr * n).astype(int)]
    return fr, acc / N_PERM / max(y.sum(), 1)


def save(df, out, name, count_cols):
    df.to_csv(os.path.join(out, name + ".private.tsv"), sep="\t", index=False)
    pub = df.copy().astype(object)
    hide = pd.Series(False, index=df.index)
    for c in count_cols:
        m = (df[c] > 0) & (df[c] < THRESH)
        pub.loc[m, c] = f"<{THRESH}"
        hide |= m
    for c in ("rate", "lo", "hi", "share", "frac_errors_removed", "remaining_rate", "auc", "auc_lo", "auc_hi"):
        if c in pub:
            pub.loc[hide, c] = ""
    pub.to_csv(os.path.join(out, name + ".pub.tsv"), sep="\t", index=False)


def main():
    h = os.path.expanduser("~")
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=h + "/s13/130")
    ap.add_argument("--call-table", default=h + "/s06/internal/71/call_table.INTERNAL.tsv")
    ap.add_argument("--trio", default=h + "/s06/internal/76/units_B_all.INTERNAL.tsv")
    ap.add_argument("--pair", default=h + "/s06/internal/76/units_pair_complete.INTERNAL.tsv")
    ap.add_argument("--ladder", default=h + "/s10b/ladder_data.private.json")
    ap.add_argument("--carriers", default=h + "/s10b/pool_carriers.private.tsv")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    n_alleles, hap_keys, pg_keys = load_carriers(a.ladder, a.carriers)

    # E1 cohort sizes over every called haplotype
    cols = ["person_id", "hap", "gene", "x_template_distance", "x_cds_share_lt10", "x_n_aa_changes"]
    C = pd.read_csv(a.call_table, sep="\t", usecols=cols, low_memory=False)
    C["gene"] = C.gene.map(gname)
    C["c3"] = [f"{p}|{g}|{hnum(x)}" in hap_keys for p, g, x in zip(C.person_id.astype(str), C.gene, C.hap)]
    C["cohort"] = cohort(C.x_template_distance.to_numpy(float), C.c3.to_numpy())
    depth_edges = edges(C.x_cds_share_lt10.dropna().to_numpy(float))
    rows = []
    for scope, sub in (("classical8", C[C.gene.isin(CLASSICAL8)]), ("all_genes", C)):
        for co, g in sub.groupby("cohort"):
            rows.append(dict(scope=scope, cohort=co, n_calls=len(g), n_persons=g.person_id.nunique(),
                             share=len(g) / len(sub)))
    save(pd.DataFrame(rows), a.out, "cohort_sizes", ["n_calls", "n_persons"])
    by_gene = C.groupby(["gene", "cohort"]).size().rename("n_calls").reset_index()
    save(by_gene, a.out, "cohort_sizes_by_gene", ["n_calls"])

    # records
    T = pd.read_csv(a.trio, sep="\t", low_memory=False)
    T["err"] = (T.cls == "neither").astype(int)
    T["oth"] = (T.cls == "other_parent_only").astype(int)
    T["c3"] = [f"{p}|{gname(g)}|{hnum(x)}" in hap_keys for p, g, x in zip(T.person.astype(str), T.gene, T.hap)]
    T["cohort"] = cohort(T.c_template_distance.to_numpy(float), T.c3.to_numpy())
    Q = pd.read_csv(a.pair, sep="\t", low_memory=False)
    Q["err"] = Q.bad.astype(int)
    Q["oth"] = 0
    Q["c3"] = [f"{p}|{gname(g)}" in pg_keys or f"{o}|{gname(g)}" in pg_keys
               for p, o, g in zip(Q.person.astype(str), Q.other.astype(str), Q.gene)]
    Q["cohort"] = cohort(Q.w_template_distance.to_numpy(float), Q.c3.to_numpy())
    sets = {"trio": (T, ["c", "w"]), "pair": (Q, ["w"])}

    e2, e3, e4, e5, pol = [], [], [], [], []
    for sname, (R, prefs) in sets.items():
        for scope, sub in (("classical8", R[R.gene.map(gname).isin(CLASSICAL8)]), ("all_genes", R)):
            sub = sub.reset_index(drop=True)
            # E2 error by cohort
            for co in ["C1", "C2", "C3", "all"]:
                g = sub if co == "all" else sub[sub.cohort == co]
                e2.append(dict(rate_row(g.err.sum(), len(g), record_set=sname, scope=scope, cohort=co),
                               other_parent_only=int(g.oth.sum())))
            fam = sub["family"] if "family" in sub else sub["person"]
            for pre in prefs:
                for metric in ("template_distance", "cds_share_lt10", "n_aa_changes"):
                    v = sub[f"{pre}_{metric}"].to_numpy(float)
                    b = bin_metric(metric, v, depth_edges)
                    for subset, m in (("all", np.ones(len(sub), bool)), ("C1_only", (sub.cohort == "C1").to_numpy()),
                                      ("C2_only", (sub.cohort == "C2").to_numpy())):
                        for lab in pd.unique(b[m]):
                            mm = m & (b == lab)
                            e3.append(rate_row(sub.err[mm].sum(), mm.sum(), record_set=sname, scope=scope,
                                               metric_of=pre, metric=metric, subset=subset, bin=lab))
                    # E5 AUC (cluster bootstrap by family) on all records and on C2 records
                    for subset, m in (("all", np.ones(len(sub), bool)), ("C2_only", (sub.cohort == "C2").to_numpy())):
                        y, s = sub.err[m].to_numpy(), np.nan_to_num(v[m], nan=0.0)
                        lo, hi = boot_auc(y, s, fam[m].to_numpy())
                        e5.append(dict(record_set=sname, scope=scope, metric_of=pre, metric=metric, subset=subset,
                                       n=int(m.sum()), errors=int(y.sum()), auc=auc(y, s), auc_lo=lo, auc_hi=hi))
                # E4 curves on the oriented 'worse of' prefix; combined = mean percentile rank of distance and depth
                d = np.nan_to_num(sub[f"{pre}_template_distance"].to_numpy(float), nan=0.0)
                p = np.nan_to_num(sub[f"{pre}_cds_share_lt10"].to_numpy(float), nan=0.0)
                comb = (pd.Series(d).rank(pct=True) + pd.Series(p).rank(pct=True)).to_numpy() / 2
                for mname, s in (("template_distance", d), ("cds_share_lt10", p), ("combined", comb)):
                    fr, cap = curve(sub.err.to_numpy(), s)
                    for f, c in zip(fr, cap):
                        e4.append(dict(record_set=sname, scope=scope, metric_of=pre, metric=mname,
                                       frac_pruned=f, frac_errors_removed=c, errors=int(sub.err.sum()), n=len(sub)))
                # fixed policies: prune a flagged set of calls; cost = share of all called haplotypes flagged
                Cs = C if scope == "all_genes" else C[C.gene.isin(CLASSICAL8)]
                hi_depth = depth_edges[1]
                for pname, rflag, cflag in (
                        ("prune d>0", d > 0, Cs.x_template_distance.fillna(0) > 0),
                        ("prune d>=3", d >= 3, Cs.x_template_distance.fillna(0) >= 3),
                        ("prune d>=6", d >= 6, Cs.x_template_distance.fillna(0) >= 6),
                        ("prune depth high tertile", p > hi_depth, Cs.x_cds_share_lt10.fillna(0) > hi_depth),
                        ("prune d>=3 or depth high", (d >= 3) | (p > hi_depth),
                         (Cs.x_template_distance.fillna(0) >= 3) | (Cs.x_cds_share_lt10.fillna(0) > hi_depth))):
                    k_out, n_out = sub.err[rflag].sum(), rflag.sum()
                    rem = rate_row(sub.err.sum() - k_out, len(sub) - n_out)
                    pol.append(dict(record_set=sname, scope=scope, metric_of=pre, policy=pname,
                                    calls_flagged_share=float(cflag.mean()), records_flagged=int(n_out),
                                    errors_removed=int(k_out), errors=int(sub.err.sum()), n=len(sub),
                                    remaining_rate=rem["rate"], lo=rem["lo"], hi=rem["hi"]))
    save(pd.DataFrame(e2), a.out, "error_by_cohort", ["n", "errors", "other_parent_only"])
    save(pd.DataFrame(e3), a.out, "error_by_metric_bin", ["n", "errors"])
    save(pd.DataFrame(e4), a.out, "prune_curves", ["errors"])
    save(pd.DataFrame(pol), a.out, "prune_policies", ["records_flagged", "errors_removed", "errors", "n"])
    save(pd.DataFrame(e5), a.out, "predictor_auc", ["n", "errors"])
    man = dict(s10_alleles_selected=n_alleles, carrier_haplotypes=len(hap_keys), depth_edges=[float(x) for x in depth_edges],
               calls=len(C), trio_records=len(T), pair_records=len(Q),
               trio_c3_records=int(T.c3.sum()), pair_c3_records=int(Q.c3.sum()))
    json.dump(man, open(os.path.join(a.out, "manifest.private.json"), "w"), indent=1)
    print(json.dumps(man))


if __name__ == "__main__":
    main()
