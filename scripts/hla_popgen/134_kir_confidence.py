#!/usr/bin/env python3
"""S06 v2 134: KIR call confidence, the counterpart of 130 (HLA) and 132 (missingness), on the same three cohorts (VM).

Calls (77's table, one row per person x hap x gene copy) fall into
  C1  template distance 0 with a protein-field label        C2  template distance > 0 with a label (collapsed to the nearest reference)
  C3  novel protein (no label)                               U   undetermined (no label)
Error is measured on called, labelled alleles only (77c):
  trio  child call with a resolved parent of origin, transmitting parent complete (both haplotypes `complete`, 51's rule), carrying >= 1 call at
        the locus and no unlabelled call there (tx_unnamed False); error = the protein-field label is in neither parent (cls3 == 'neither').
        Unlabelled transmitting-parent calls make the label unstable (S06): reported separately, never pooled.
  pair  first-degree pair at the gene, both people complete, both with labelled calls; error = no shared label (cls == 'no_share').
        Siblings can legitimately share no label, so pair rates are an upper bound.
Metrics per record: template distance, CDS distance (worse of the call and, for trios, the transmitting parent's worst call; for pairs the worse of the
two people), and for the 10 genes that exist on GRCh38 chr19 the depth summaries of 78 (share of CDS < 10x, median / 10th percentile / minimum
relative CDS depth, MAPQ share), oriented so higher = worse. The other KIR genes have no depth (NaN, never 0).
Missingness counterpart: KIR gene presence is biology; what is missing is assembly completeness, so haplotype completeness (complete / one_contig_open_end /
multi_contig / empty) by ancestry and platform is reported.
Reads (VM): ~/s06/internal/77/{kir_calls,kir_haps}.INTERNAL.tsv, ~/s06/internal/77c/records_{trio,pair}_kir.INTERNAL.tsv,
            ~/s06/internal/78/gene_depth_metrics.INTERNAL.tsv, ~/pipeline_outputs/cohort_membership.tsv. Helpers come from 130 (same folder).
Writes <out>/k_*.private.tsv (exact) and .pub.tsv (counts 1-19 masked).
VM:  python3 scripts/hla_popgen/134_kir_confidence.py --out ~/s13/134
"""
import argparse
import importlib.util
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("m130", os.path.join(HERE, "130_call_confidence.py"))
M = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M)

DEPTH = {"cds_share_lt10": 1, "rel_cds": -1, "cds_p10_rel": -1, "cds_min_rel": -1, "mapq_frac": -1}   # sign makes higher = worse
SCOPE_NOTE = "KIR genes with depth: the 10 on GRCh38 chr19"


def tf(x):
    return x.astype(str).str.lower().isin(["true", "1", "1.0"])


def cohort_calls(K):
    lab = K["label3"].notna() & K["label3"].astype(str).ne("")
    tier = K["novelty_tier"].astype(str)
    d = K["template_distance"].fillna(0)
    co = np.where(tier.eq("novel_protein"), "C3", np.where(~lab, "U", np.where(d == 0, "C1", "C2")))
    return pd.Series(co, index=K.index)


def worst_by(df, keys, cols):
    return df.groupby(keys)[cols].max().reset_index()


def prep_trio(T, D):
    T = T[T["mode"] == "gene"].copy()
    for c in [c for c in ("tx_complete", "tx_unnamed", "child_complete") if c in T]:
        T[c] = tf(T[c])
    T["tx_n_calls"] = T["tx_n_calls"].fillna(0)
    if "origin" in T:
        T = T[T["origin"].isin(["pat", "mat"])]
    nm = T[(T["stream"] == "named") & T["cls3"].isin(["match", "neither", "only_other_parent"]) & (T["tx_n_calls"] >= 1)].copy()
    nm["err"] = (nm["cls3"] == "neither").astype(int)
    nm["oth"] = (nm["cls3"] == "only_other_parent").astype(int)
    nm = nm.merge(D, left_on=["child", "gene"], right_on=["person_id", "gene"], how="left")
    nm["cohort"] = np.where(nm["template_distance"].fillna(0) == 0, "C1", "C2")
    return T, nm


def prep_pair(P, K, D):
    P = P[(P["mode"] == "gene") & (P["kind"] == "pair") & P["cls"].isin(["share", "no_share"])].copy()
    P["both_complete"] = tf(P["both_complete"])
    P["err"] = (P["cls"] == "no_share").astype(int)
    w = worst_by(K, ["person_id", "gene"], ["template_distance", "cds_distance"]).rename(columns={"person_id": "p"})
    for who in ("person", "other"):
        x = w.rename(columns={"p": who, "template_distance": f"td_{who}", "cds_distance": f"cd_{who}"})
        P = P.merge(x, left_on=[who, "locus"], right_on=[who, "gene"], how="left").drop(columns=["gene"])
        y = D.rename(columns={"person_id": who}).add_suffix(f"_{who}").rename(columns={f"{who}_{who}": who, f"gene_{who}": "gene"})
        P = P.merge(y, left_on=[who, "locus"], right_on=[who, "gene"], how="left").drop(columns=["gene"])
    P["template_distance"] = P[["td_person", "td_other"]].max(axis=1)
    P["cds_distance"] = P[["cd_person", "cd_other"]].max(axis=1)
    for m, sg in DEPTH.items():
        P[m] = np.where(sg > 0, P[[f"{m}_person", f"{m}_other"]].max(axis=1), -P[[f"{m}_person", f"{m}_other"]].min(axis=1))
    P["cohort"] = np.where(P["template_distance"].fillna(0) == 0, "C1", "C2")
    P["family"] = P["person"].astype(str)
    return P


def orient_depth(df):
    for m, sg in DEPTH.items():
        if m in df and sg < 0:
            df[m] = -df[m]
        # (shared with pair records, where max/min already produced the oriented worst value)
    return df


def depth_cols_for(df):
    return [m for m in DEPTH if m in df and df[m].notna().any()]


def bins(metric, v, edges):
    v = np.asarray(v, float)
    if metric == "template_distance":
        return np.where(np.isnan(v), "NA", pd.cut(np.nan_to_num(v), [-1, 0, 1, 2, 5, 1e9], labels=["0", "1", "2", "3-5", "6+"]).astype(str))
    if metric == "cds_distance":
        return np.where(np.isnan(v), "NA", pd.cut(np.nan_to_num(v), [-1, 0, 1, 3, 1e9], labels=["0", "1", "2-3", "4+"]).astype(str))
    out = np.full(len(v), "NA", dtype=object)
    ok = ~np.isnan(v)
    if ok.sum() >= 8:
        out[ok] = pd.qcut(pd.Series(v[ok]).rank(method="first"), 4, labels=["Q1 best", "Q2", "Q3", "Q4 worst"]).astype(str).to_numpy()
    return out


def analyse(S, out, K, D, T, P, H, coh):
    man = {}
    # K1 cohort sizes
    K = K.copy()
    K["cohort"] = cohort_calls(K)
    rows = [dict(scope="all", cohort=c, n_calls=len(g), n_persons=g["person_id"].nunique()) for c, g in K.groupby("cohort")]
    tot = len(K)
    for r in rows:
        r["share"] = r["n_calls"] / tot
    rows += [dict(scope=gn, cohort=c, n_calls=len(g), n_persons=g["person_id"].nunique(), share=len(g) / len(K[K.gene == gn])) for (gn, c), g in K.groupby(["gene", "cohort"])]
    M.save(pd.DataFrame(rows), out, "k_cohort_sizes", ["n_calls", "n_persons"])
    # trio / pair records
    Tall, Tn = prep_trio(T, D)
    Tn = orient_depth(Tn)
    P = prep_pair(P, K, D)
    sets = {"trio": (Tn[Tn["tx_complete"] & ~Tn["tx_unnamed"]].reset_index(drop=True), "child"), "pair": (P[P["both_complete"]].reset_index(drop=True), "family")}
    instab = []
    for lab, g in (("tx_unnamed False", Tn[~Tn["tx_unnamed"]]), ("tx_unnamed True", Tn[Tn["tx_unnamed"]])):
        instab.append(M.rate_row(g["err"].sum(), len(g), transmitting_parent=lab))
    M.save(pd.DataFrame(instab), out, "k_label_instability", ["n", "errors"])
    e2, e3, e4, e5, pol = [], [], [], [], []
    for sname, (R, clus) in sets.items():
        for co in ("C1", "C2", "all"):
            g = R if co == "all" else R[R.cohort == co]
            e2.append(M.rate_row(g["err"].sum(), len(g), record_set=sname, cohort=co))
        mets = ["template_distance", "cds_distance"] + depth_cols_for(R)
        grp = R[clus].astype(str).to_numpy()
        for m in mets:
            v = R[m].to_numpy(float)
            b = bins(m, v, None)
            for lab in pd.unique(b):
                mm = b == lab
                e3.append(M.rate_row(R["err"][mm].sum(), mm.sum(), record_set=sname, metric=m, bin=lab))
            ok = ~np.isnan(v)
            if ok.sum() > 50 and R["err"][ok].sum() > 0:
                lo, hi = M.boot_auc(R["err"][ok].to_numpy(), v[ok], grp[ok])
                e5.append(dict(record_set=sname, metric=m, n=int(ok.sum()), errors=int(R["err"][ok].sum()), auc=M.auc(R["err"][ok].to_numpy(), v[ok]), auc_lo=lo, auc_hi=hi))
                fr, cap = M.curve(R["err"][ok].to_numpy(), v[ok])
                for f, c in zip(fr, cap):
                    e4.append(dict(record_set=sname, metric=m, frac_pruned=f, frac_errors_removed=c, errors=int(R["err"][ok].sum()), n=int(ok.sum())))
        d = np.nan_to_num(R["template_distance"].to_numpy(float), nan=0.0)
        calls_d = K["template_distance"].fillna(0)
        for pname, rflag, cshare in (("prune d>0", d > 0, float((calls_d > 0).mean())), ("prune d>=3", d >= 3, float((calls_d >= 3).mean())),
                                     ("prune d>=6", d >= 6, float((calls_d >= 6).mean()))):
            k_out, n_out = R["err"][rflag].sum(), int(rflag.sum())
            rem = M.rate_row(R["err"].sum() - k_out, len(R) - n_out)
            pol.append(dict(record_set=sname, policy=pname, calls_flagged_share=cshare, records_flagged=n_out, errors_removed=int(k_out), errors=int(R["err"].sum()), n=len(R),
                            remaining_rate=rem["rate"], lo=rem["lo"], hi=rem["hi"]))
    M.save(pd.DataFrame(e2), out, "k_error_by_cohort", ["n", "errors"])
    M.save(pd.DataFrame(e3), out, "k_error_by_metric_bin", ["n", "errors"])
    M.save(pd.DataFrame(e4), out, "k_prune_curves", ["errors"])
    M.save(pd.DataFrame(pol), out, "k_prune_policies", ["records_flagged", "errors_removed", "errors", "n"])
    M.save(pd.DataFrame(e5), out, "k_predictor_auc", ["n", "errors"])
    # haplotype completeness (the missingness counterpart)
    H = H.copy()
    H["anc"] = H["person_id"].astype(str).map(coh["anc"]).fillna("NA")
    H["plat"] = H["person_id"].astype(str).map(coh["plat"]).fillna("NA")
    hc = []
    for key, g in [("all", H)] + [("ancestry:" + k, x) for k, x in H.groupby("anc")] + [("platform:" + k, x) for k, x in H.groupby("plat")]:
        for c, x in g.groupby("hap_class"):
            lo, hi = M.wilson(len(x), len(g))
            hc.append(dict(group=key, hap_class=c, n=len(x), total=len(g), rate=len(x) / len(g), lo=lo, hi=hi))
    M.save(pd.DataFrame(hc), out, "k_hap_completeness", ["n", "total"])
    man.update(calls=len(K), haps=len(H), trio_records=len(T), trio_named_informative=len(Tn), trio_eval=len(sets["trio"][0]), pair_eval=len(sets["pair"][0]),
               persons_depth=int(D["person_id"].nunique()), note=SCOPE_NOTE)
    json.dump(man, open(os.path.join(out, "manifest_134.private.json"), "w"), indent=1)
    return man


def load_all(a):
    rd = lambda p, **k: pd.read_csv(p, sep="\t", comment="#", low_memory=False, **k)
    K = rd(a.calls, dtype={"person_id": str})
    H = rd(a.haps, dtype={"person_id": str})
    D = rd(a.depth, dtype={"person_id": str})
    D = D[["person_id", "gene", "depth_defined"] + [m for m in DEPTH if m in D]].copy()
    T = rd(a.trio, dtype={"child": str, "tx": str, "other": str})
    P = rd(a.pair, dtype={"person": str, "other": str})
    C = rd(a.cohort, dtype={"person_id": str})
    anc = C["ancestry_pred"].astype(str).str.upper().where(~C["ancestry_pred"].isna(), "NA")
    coh = dict(anc=dict(zip(C["person_id"], anc)), plat=dict(zip(C["person_id"], C["platform"].astype(str))))
    return K, H, D, T, P, coh


def main():
    h = os.path.expanduser("~")
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=h + "/s13/134")
    ap.add_argument("--calls", default=h + "/s06/internal/77/kir_calls.INTERNAL.tsv")
    ap.add_argument("--haps", default=h + "/s06/internal/77/kir_haps.INTERNAL.tsv")
    ap.add_argument("--depth", default=h + "/s06/internal/78/gene_depth_metrics.INTERNAL.tsv")
    ap.add_argument("--trio", default=h + "/s06/internal/77c/records_trio_kir.INTERNAL.tsv")
    ap.add_argument("--pair", default=h + "/s06/internal/77c/records_pair_kir.INTERNAL.tsv")
    ap.add_argument("--cohort", default=h + "/pipeline_outputs/cohort_membership.tsv")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    K, H, D, T, P, coh = load_all(a)
    print(json.dumps(analyse(None, a.out, K, D, T, P, H, coh)))


if __name__ == "__main__":
    main()
