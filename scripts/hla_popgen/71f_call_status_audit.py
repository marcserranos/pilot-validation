#!/usr/bin/env python3
"""S06 71f: what does "one called haplotype" mean, and where do the trio / pair "errors" come from? A read-only audit (VM, aggregate output).

Three tables, all counts (no person ids):
  call_status_taxonomy.tsv   per gene (and ALL / CLASSICAL8): for every person x gene, the status of each of the two haplotypes, from Table 1 BEFORE any filter:
                               valid         a call exists (named allele, or novel `new` call)
                               undetermined  Immuannot annotated the gene but could not decide the allele (consensus == "undetermined")
                               no_gene_row   the haplotype has other MHC rows for this person but none for this gene
                               no_haplotype  the haplotype has no MHC rows at all for this person (assembly or trimming gave nothing)
                             counted as an unordered pair of statuses, e.g. valid+no_gene_row = "one call, gene absent on the other haplotype"
  trio_candidates.tsv        how many children have 2+ candidate trios (57b used the first; 71 used the last before the fix of 2026-10-21)
  trio_parent_status.tsv     for every child haplotype call whose parent of origin is resolved (57b): how many calls the transmitting parent has at that gene
                             (0 = the record is SKIPPED by 71, 1, 2), the child's allele name status, and how often the child's 2-field allele is absent from
                             the parent's calls. This is the "what if the parent has no call" question.
  trio_pc_grid_<metric>.tsv  error rate of the trio records on a 40 x 40 grid of cohort percentile ranks, child's call on one axis and the transmitting
                             parent's call on the other (100 = worst), by stratum (all / parent has two calls / parent has one call): n and errors per cell.
Run on the named-allele records (71_named_v2) so novel calls are not counted as errors.

    python3 71f_call_status_audit.py --internal ~/s06/internal/71_named_v2/ --out ~/s06/results/71f/
Tests: scripts/hla_popgen/tests/test_71f_call_status_audit.py
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _s05_common as C  # noqa: E402

G8 = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def hap_status_table(t1, genes_of_interest=None):
    """Table 1 (person_id, hap, gene, consensus) -> one row per person x gene with the status of hap1 and hap2."""
    t = t1[["person_id", "hap", "gene", "consensus"]].copy()
    t["hap"] = t["hap"].astype(str)
    t["hap"] = t["hap"].str.replace("hap", "", regex=False)
    t = t[t["hap"].isin(["1", "2"])]
    cons = t["consensus"].astype(str)
    t["valid"] = t["consensus"].notna() & ~cons.isin(["undetermined", "nan", "None", ""])
    t["undet"] = cons.eq("undetermined")
    hap_has_rows = t.groupby(["person_id", "hap"]).size().unstack(fill_value=0)
    persons = hap_has_rows.index
    genes = sorted(t["gene"].unique()) if genes_of_interest is None else genes_of_interest
    g = t[t["gene"].isin(genes)].groupby(["person_id", "gene", "hap"]).agg(valid=("valid", "max"), undet=("undet", "max")).reset_index()
    base = pd.MultiIndex.from_product([persons, genes, ["1", "2"]], names=["person_id", "gene", "hap"]).to_frame(index=False)
    base = base.merge(g, how="left", on=["person_id", "gene", "hap"])
    hh = hap_has_rows.stack().rename("n_rows_on_hap").reset_index()
    hh = hh.rename(columns={"level_1": "hap"}) if "hap" not in hh.columns else hh
    base = base.merge(hh[["person_id", "hap", "n_rows_on_hap"]], how="left", on=["person_id", "hap"])
    base["n_rows_on_hap"] = base["n_rows_on_hap"].fillna(0)
    status = np.where(base["valid"].fillna(False).astype(bool), "valid",
                      np.where(base["undet"].fillna(False).astype(bool), "undetermined",
                               np.where(base["n_rows_on_hap"] > 0, "no_gene_row", "no_haplotype")))
    base["status"] = status
    w = base.pivot_table(index=["person_id", "gene"], columns="hap", values="status", aggfunc="first").reset_index()
    w.columns = ["person_id", "gene", "hap1", "hap2"]
    return w


def taxonomy(w):
    rows = []
    order = {"valid": 0, "undetermined": 1, "no_gene_row": 2, "no_haplotype": 3}
    w = w.copy()
    w["pair"] = [" + ".join(sorted([a, b], key=lambda s: order[s])) for a, b in zip(w["hap1"], w["hap2"])]
    for scope, sub in (("ALL", w), ("CLASSICAL8", w[w["gene"].isin(G8)])):
        for pair, n in sub["pair"].value_counts().items():
            rows.append(dict(scope=scope, gene="ALL", status_pair=pair, n_person_genes=int(n), pct=100.0 * n / len(sub)))
    for gn, sub in w.groupby("gene"):
        for pair, n in sub["pair"].value_counts().items():
            rows.append(dict(scope="gene", gene=gn, status_pair=pair, n_person_genes=int(n), pct=100.0 * n / len(sub)))
    return pd.DataFrame(rows)


def trio_parent_status(origin, ct, trio_of):
    """Every child haplotype call with a resolved parent of origin: how many calls does the transmitting parent have at the gene (0 / 1 / 2) and how often is the
    child's 2-field allele missing from the parent's calls."""
    ncalls = ct.groupby(["person_id", "gene"])["hap"].nunique().to_dict()
    alle = {k: set(d["tf"].dropna()) for k, d in ct.groupby(["person_id", "gene"])}
    named = ct.set_index(["person_id", "hap", "gene"])["named"].to_dict()
    tf = ct.set_index(["person_id", "hap", "gene"])["tf"].to_dict()
    rows = []
    for c, h, g, og in zip(origin["child"], origin["hap"], origin["gene"], origin["origin"]):
        if og not in ("pat", "mat") or c not in trio_of or (c, h, g) not in tf:
            continue
        par = trio_of[c][0] if og == "pat" else trio_of[c][1]
        n = ncalls.get((par, g), 0)
        other = trio_of[c][1] if og == "pat" else trio_of[c][0]
        absent = isinstance(tf[(c, h, g)], str) and tf[(c, h, g)] not in alle.get((par, g), set())
        in_other = isinstance(tf[(c, h, g)], str) and tf[(c, h, g)] in alle.get((other, g), set())
        rows.append((g in G8, n, bool(named.get((c, h, g), True)), absent, absent and in_other))
    d = pd.DataFrame(rows, columns=["classical8", "parent_calls", "child_named", "child_allele_missing_from_parent", "missing_but_in_other_parent"])
    out = []
    for scope, sub in (("ALL", d), ("CLASSICAL8", d[d["classical8"]])):
        for (pc, nm), s in sub.groupby(["parent_calls", "child_named"]):
            out.append(dict(scope=scope, parent_calls_at_gene=int(pc), child_call_named=bool(nm), child_calls=len(s),
                            child_allele_missing_from_parent=int(s["child_allele_missing_from_parent"].sum()),
                            missing_but_in_other_parent=int(s["missing_but_in_other_parent"].sum()),
                            missing_pct=100.0 * s["child_allele_missing_from_parent"].mean()))
    return pd.DataFrame(out)


def pc_grid(rec, metric, cohort_ct, strata, bins=40):
    """Error rate over (child's percentile, parent's percentile) for one metric; rec = records_trio with c_/p_ oriented columns."""
    S = sys.modules["m71b"]
    a = S.percentile_rank(rec["c_" + metric], cohort_ct["o_" + metric])
    b = S.percentile_rank(rec["p_" + metric], cohort_ct["o_" + metric])
    edges = np.linspace(0, 100, bins + 1)
    out = []
    for sname, m in strata.items():
        ok = np.isfinite(a) & np.isfinite(b) & m
        ia = np.clip(np.digitize(a[ok], edges) - 1, 0, bins - 1)
        ib = np.clip(np.digitize(b[ok], edges) - 1, 0, bins - 1)
        y = rec["bad"].to_numpy(int)[ok]
        n = np.zeros((bins, bins))
        k = np.zeros((bins, bins))
        np.add.at(n, (ia, ib), 1)
        np.add.at(k, (ia, ib), y)
        for i in range(bins):
            for j in range(bins):
                if n[i, j] > 0:
                    out.append(dict(stratum=sname, child_lo=edges[i], parent_lo=edges[j], n=int(n[i, j]), k=int(k[i, j])))
    return pd.DataFrame(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--internal", default=os.path.expanduser("~/s06/internal/71_named_v2/"))
    ap.add_argument("--out", default=os.path.expanduser("~/s06/results/71f/"))
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    Q = _load("m71", "71_call_quality_table.py")
    _load("m71b", "71b_call_quality_stats.py")
    w = lambda df, n: C.write_tsv(df, os.path.join(a.out, n))  # noqa: E731
    t1 = C.load_table1(C.DEFAULT_TABLE1)
    w(taxonomy(hap_status_table(t1)), "call_status_taxonomy.tsv")
    ct = Q.read_int(os.path.join(a.internal, "call_table.INTERNAL.tsv"), dtype={"person_id": str, "hap": str})
    T = _load("t57", "57_yak_trioeval_pilot.py")
    cohort = C.load_cohort(C.DEFAULT_COHORT)
    lr = set(cohort.loc[cohort["in_lr"].astype(str).str.lower() == "true", "person_id"])
    origin = Q.read_int(Q.ORIGIN57B, dtype={"child": str, "contig": str})
    kids = set(origin["child"])
    cand = [t for t in T.find_trios(lr) if t[0] in kids]
    trio_of = {}
    for t in cand:                      # first candidate trio per child, as 57b used
        trio_of.setdefault(t[0], (t[1], t[2]))
    multi = pd.Series([t[0] for t in cand]).value_counts()
    pd.DataFrame([dict(children=len(trio_of), children_with_2plus_candidate_trios=int((multi >= 2).sum()), max_candidates=int(multi.max()),
                       candidate_trios=len(cand))]).to_csv(os.path.join(a.out, "trio_candidates.tsv"), sep="\t", index=False)
    ct["hap"] = ct["hap"].astype(str)
    w(trio_parent_status(origin, ct, trio_of), "trio_parent_status.tsv")
    rec = pd.read_csv(os.path.join(a.internal, "records_trio.INTERNAL.tsv"), sep="\t", dtype={"person": str, "other": str})
    pn = ct.groupby(["person_id", "gene"])["hap"].nunique().to_dict()
    pcalls = np.array([pn.get((p, g), 0) for p, g in zip(rec["other"], rec["gene"])])
    strata = {"all": np.ones(len(rec), bool), "parent_two_calls": pcalls == 2, "parent_one_call": pcalls == 1}
    for metric in ("rel_ars", "rel_cds", "cds_share_lt10"):
        w(pc_grid(rec, metric, ct, strata), f"trio_pc_grid_{metric}.tsv")
    C.log("DONE 71f")
    return 0


if __name__ == "__main__":
    sys.exit(main())
