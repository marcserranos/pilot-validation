#!/usr/bin/env python3
"""S06 71e: the incomplete-genotype ("dropout") outcome, over the whole cohort (VM).

Finding that motivates it (named-allele analysis, 2026-10-20): most trio mismatches and pair disagreements are not a wrong allele among a fully
called genotype but a person with only ONE called haplotype at the gene (the second allele is missing, so the person looks homozygous or the
relative's allele is "not found"). In the classical eight, 59 of 67 pair errors and 29 of 50 trio mismatches sit in such records, against about 0.2%
and 1% error when everyone is fully called. That is observable for every call, not only for relatives, so it can be studied on the whole cohort:

    outcome   incomplete = the person has a call on only one haplotype at that gene (n_called < 2)
    unit      one row per person x gene (metrics = the mean of the person's calls at the gene)
    outputs   dropout_auc.tsv    AUC of every metric (and models / rules) for predicting incompleteness; within-gene AUC
              dropout_prune.tsv  remove the worst x% of person-genes by the metric: share still incomplete among the rest, share of all
                                 incomplete genotypes removed (0.25% steps), plus the random band
              dropout_curves.tsv rolling incomplete rate along each metric
              dropout_rate_by_gene.tsv  incomplete share per gene (genes that are absent on some haplotypes, DRB3/4/5, are biology, not dropout)

Scopes: `all` genes and `classical8`. Not a typing-error label by itself: whether an incomplete genotype is a failure of the call or a gene
that is truly absent on one haplotype has to be read gene by gene (see the rate table).

VM:  pixi run -e spechla -- python3 scripts/hla_popgen/71e_dropout_outcome.py [--scope classical8]
Tests: scripts/hla_popgen/tests/test_71e_dropout_outcome.py
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


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


S = sys.modules.get("m71b") or _load("m71b", "71b_call_quality_stats.py")
Q = S.Q
G8 = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
INT = os.path.expanduser("~/s06/internal/71/")
OUT = os.path.expanduser("~/s06/results/71e/")
METRIC_KEYS = [k for k in S.ALL_KEYS if k != "n_called"]


def log(m):
    C.log("[71e] " + m)


def person_gene_table(ct):
    """Call table -> one row per person x gene: incomplete flag and the MEAN over the person's calls of every oriented metric (a max over two calls
    would be larger than a single call by construction, which would make complete genotypes look worse than incomplete ones)."""
    ocols = ["o_" + k for k in Q.KEYS if k != "n_called"]
    g = ct.groupby(["person_id", "gene"])
    out = g[ocols].mean()
    out["n_called"] = g["hap"].nunique()
    out["plat"] = g["plat"].first()
    out = out.reset_index().rename(columns={"person_id": "person"})
    out["incomplete"] = (out["n_called"] < 2).astype(int)
    return out


def analyse_dropout(pg, out_dir, scope="all"):
    os.makedirs(out_dir, exist_ok=True)
    w = lambda df, name: C.write_tsv(df, os.path.join(out_dir, name))  # noqa: E731
    y = pg["incomplete"].to_numpy(int)
    v = pg.rename(columns={"incomplete": "bad"})
    rate = [dict(gene=gn, n_person_genes=len(d), n_incomplete=int(d["incomplete"].sum()), incomplete_pct=100 * d["incomplete"].mean())
            for gn, d in pg.groupby("gene")]
    rate.append(dict(gene="ALL", n_person_genes=len(pg), n_incomplete=int(y.sum()), incomplete_pct=100 * y.mean()))
    w(pd.DataFrame(rate).sort_values("incomplete_pct", ascending=False), "dropout_rate_by_gene.tsv")
    scores = {k: pg["o_" + k].to_numpy(float) for k in METRIC_KEYS if ("o_" + k) in pg}
    pcts = S.metric_percentiles(v, v, [k for k in S.ALL_KEYS if k != "n_called"])
    for rname, keys_ in (("rule:worst_rank_depth", [k for k in S.DEPTH_KEYS]), ("rule:worst_rank_all", [k for k in S.ALL_KEYS if k != "n_called"])):
        stack = np.vstack([pcts[k] for k in keys_ if np.isfinite(pcts[k]).any()])
        scores[rname] = np.nanmax(stack, axis=0)
    X, _ = S.design_matrix(v, [k for k in S.ALL_KEYS if k != "n_called"])
    scores["model:all"] = S.oof_risk(X, y, pg["person"])
    Xd, _ = S.design_matrix(v, S.DEPTH_KEYS)
    scores["model:depth"] = S.oof_risk(Xd, y, pg["person"])
    Xt, _ = S.tail_design(pcts, [k for k in S.ALL_KEYS if k != "n_called"], v)
    scores["model:tail"] = S.oof_risk(Xt, y, pg["person"])
    try:
        Xg = np.column_stack([pg["o_" + k].to_numpy(float) for k in S.ALL_KEYS if k != "n_called" and ("o_" + k) in pg])
        scores["model:gbm"] = S.gbm_oof(Xg, y, pg["person"])
    except ImportError:
        pass
    aucs, prunes, curves = [], [], []
    for k, o in scores.items():
        fam = "model" if k.startswith("model:") else ("rule" if k.startswith("rule:") else Q.FAMILY[k])
        aucs.append(dict(scope=scope, metric=k, family=fam, n=int(np.isfinite(o).sum()), n_incomplete=int(y[np.isfinite(o)].sum()),
                         auc=S.auc_score(o, y), auc_within_gene=S.auc_within_gene(o, y, pg["gene"])))
        pc = S.prune_curve(o, y, cohort_o=None if k.startswith(("model:", "rule:")) else o)
        if len(pc):
            pc.insert(0, "metric", k)
            pc.insert(0, "scope", scope)
            prunes.append(pc)
        if not k.startswith(("model:", "rule:")):
            c = S.rolling_curve(o, y, k)
            if len(c):
                c.insert(0, "metric", k)
                c.insert(0, "scope", scope)
                curves.append(c)
    w(pd.DataFrame(aucs).sort_values("auc", ascending=False), "dropout_auc.tsv")
    w(pd.concat(prunes, ignore_index=True), "dropout_prune.tsv")
    w(pd.concat(curves, ignore_index=True), "dropout_curves.tsv")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--internal-dir", default=INT)
    ap.add_argument("--out-dir", default=OUT)
    ap.add_argument("--scope", choices=["all", "classical8"], default="all")
    a = ap.parse_args(argv)
    ct = Q.read_int(os.path.join(a.internal_dir, "call_table.INTERNAL.tsv"), dtype={"person_id": str, "hap": str})
    pg = person_gene_table(ct)
    out = a.out_dir
    if a.scope == "classical8":
        pg = pg[pg["gene"].isin(G8)].reset_index(drop=True)
        out = os.path.join(out, "classical8")
    log(f"{len(pg)} person-genes, {int(pg['incomplete'].sum())} incomplete ({100 * pg['incomplete'].mean():.1f}%)")
    analyse_dropout(pg, out, a.scope)
    log("DONE 71e")
    return 0


if __name__ == "__main__":
    sys.exit(main())
