#!/usr/bin/env python3
"""S06 v2 132: where do missing HLA calls fall, how many are there, and who has them (VM).

A missing call is visible in the call matrix (the caller declined: no gene annotation on that haplotype), so the
matter is how much is missing, where, and whether it is unevenly spread over genes, ancestries and platforms.
Unit = person x gene x haplotype slot (2 per person-gene) from 76's slot table (status called | undetermined | missing,
hap_class = structure around the gene, also computed for called slots so rates P(missing | class) are meaningful).
Scopes: classical8 (A,B,C,DRB1,DQA1,DQB1,DPA1,DPB1) and all_fixed (the 12 anchor genes of 76).
Tables (each as <name>.private.tsv exact and <name>.pub.tsv counts 1-19 masked):
  m_overall         slots, missing, rate; persons with >= 1 missing slot
  m_per_person      persons by number of missing slots (classical8, 0..16)
  m_by_gene         rate per gene in genomic order, with the share of person-genes with 0 / 1 / 2 calls
  m_by_class        where missing slots fall (hap_class): share of missing and P(missing | class)
  m_by_ancestry     rate per ancestry, chi-square across ancestries
  m_ancestry_gene   rate per ancestry x gene
  m_by_platform     rate per platform; m_ancestry_platform ancestry x platform
  m_by_contigs      rate by number of contigs on the haplotype
  m_units_by_anc    person-gene units with 0 / 1 / 2 calls per ancestry (what the person x gene matrix looks like)
VM:  python3 scripts/hla_popgen/132_missingness.py --out ~/s13/132
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

ANCHORS = ["A", "C", "B", "MICA", "DRA", "DRB1", "DQA1", "DQB1", "TAP2", "TAP1", "DPA1", "DPB1"]   # genomic order (76)
CLASSICAL8 = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
THRESH = 20


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def rate(g, **kw):
    n, k = len(g), int(g["miss"].sum())
    lo, hi = wilson(k, n)
    return dict(kw, n=n, missing=k, rate=k / n if n else np.nan, lo=lo, hi=hi)


def save(df, out, name, count_cols=("n", "missing", "persons", "persons_any")):
    df.to_csv(os.path.join(out, name + ".private.tsv"), sep="\t", index=False)
    pub = df.copy().astype(object)
    hide = pd.Series(False, index=df.index)
    for c in count_cols:
        if c in df:
            m = (df[c] > 0) & (df[c] < THRESH)
            pub.loc[m, c] = f"<{THRESH}"
            hide |= m
    for c in df.columns:
        if c in ("rate", "lo", "hi", "share", "share_of_missing", "p_missing", "share_0", "share_1", "share_2", "pct"):
            pub.loc[hide, c] = ""
    pub.to_csv(os.path.join(out, name + ".pub.tsv"), sep="\t", index=False)


def load(slots_path, cohort_path):
    S = pd.read_csv(slots_path, sep="\t", dtype={"person": str, "hap": str}, low_memory=False)
    S["miss"] = (S["status"] == "missing").astype(int)
    C = pd.read_csv(cohort_path, sep="\t", dtype={"person_id": str}, low_memory=False)
    anc = C["ancestry_pred"].astype(str).str.upper().where(~C["ancestry_pred"].isna(), "NA")
    S["anc"] = S["person"].map(dict(zip(C["person_id"], anc))).fillna("NA")
    S["plat"] = S["person"].map(dict(zip(C["person_id"], C["platform"].astype(str)))).fillna("NA")
    return S


def chi2_p(tab):
    try:
        from scipy.stats import chi2_contingency
        return float(chi2_contingency(tab)[1])
    except Exception:
        return np.nan


def analyse(S, out):
    manifest = {}
    scopes = {"classical8": S[S.gene.isin(CLASSICAL8)], "all_fixed": S[S.gene.isin(ANCHORS)]}
    ov, pp = [], []
    for sc, d in scopes.items():
        pers = d.groupby("person")["miss"].sum()
        ov.append(dict(rate(d, scope=sc), persons=int(pers.size), persons_any=int((pers > 0).sum())))
        if sc == "classical8":
            pp = pers.value_counts().sort_index().rename_axis("n_missing_slots").reset_index(name="persons")
            pp["share"] = pp["persons"] / pp["persons"].sum()
    save(pd.DataFrame(ov), out, "m_overall")
    save(pp, out, "m_per_person", ("persons",))
    d8 = scopes["classical8"]
    # by gene in genomic order + person-gene units with 0/1/2 calls
    units = S[S.gene.isin(ANCHORS)].assign(called=lambda x: 1 - x["miss"]).groupby(["person", "gene"]).agg(called=("called", "sum"), anc=("anc", "first")).reset_index()
    rows = []
    for i, g in enumerate([x for x in ANCHORS]):
        dg = S[S.gene == g]
        u = units[units.gene == g]["called"].value_counts(normalize=True)
        rows.append(dict(rate(dg, gene=g, order=i, classical=g in CLASSICAL8), share_0=u.get(0, 0.0), share_1=u.get(1, 0.0), share_2=u.get(2, 0.0)))
    save(pd.DataFrame(rows), out, "m_by_gene")
    # where they fall
    cl = []
    for sc, d in scopes.items():
        tot = d.miss.sum()
        for c, g in d.groupby("hap_class"):
            cl.append(dict(rate(g, scope=sc, hap_class=c), share_of_missing=g.miss.sum() / tot if tot else np.nan))
    save(pd.DataFrame(cl), out, "m_by_class")
    # ancestry, platform
    a = [rate(g, ancestry=k) for k, g in d8.groupby("anc")]
    tab = pd.crosstab(d8["anc"], d8["miss"])
    for r in a:
        r["chi2_p_across_ancestries"] = chi2_p(tab)
    save(pd.DataFrame(a), out, "m_by_ancestry")
    save(pd.DataFrame([rate(g, ancestry=k[0], gene=k[1]) for k, g in d8.groupby(["anc", "gene"])]), out, "m_ancestry_gene")
    save(pd.DataFrame([rate(g, platform=k) for k, g in d8.groupby("plat")]), out, "m_by_platform")
    save(pd.DataFrame([rate(g, ancestry=k[0], platform=k[1]) for k, g in d8.groupby(["anc", "plat"])]), out, "m_ancestry_platform")
    bins = pd.cut(d8["hap_n_contigs"], [-1, 0, 1, 2, 5, 10, 1e9], labels=["0", "1", "2", "3-5", "6-10", "11+"])
    save(pd.DataFrame([rate(g, n_contigs=str(k)) for k, g in d8.groupby(bins, observed=True)]), out, "m_by_contigs")
    ua = units[units.gene.isin(CLASSICAL8)].groupby(["anc", "called"]).size().unstack(fill_value=0)
    ua = ua.div(ua.sum(1), axis=0).rename(columns={0: "share_0", 1: "share_1", 2: "share_2"}).reset_index()
    ua["n"] = units[units.gene.isin(CLASSICAL8)].groupby("anc").size().reindex(ua["anc"]).to_numpy()
    save(ua, out, "m_units_by_anc")
    manifest.update(slots=len(S), persons=int(S.person.nunique()), ancestries=sorted(S.anc.unique()), platforms=sorted(S.plat.unique()))
    json.dump(manifest, open(os.path.join(out, "manifest_132.private.json"), "w"), indent=1)
    return manifest


def main():
    h = os.path.expanduser("~")
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=h + "/s13/132")
    ap.add_argument("--slots", default=h + "/s06/internal/76/slots.INTERNAL.tsv")
    ap.add_argument("--cohort", default=h + "/pipeline_outputs/cohort_membership.tsv")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    print(json.dumps(analyse(load(a.slots, a.cohort), a.out)))


if __name__ == "__main__":
    main()
