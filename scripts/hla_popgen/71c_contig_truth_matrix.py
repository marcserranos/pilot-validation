#!/usr/bin/env python3
"""S06 71c: contig co-location versus trio truth, for every pair of MHC genes (VM).

The figure "% haplotypes with both genes on one contig" (Marc, 10-20) measures an ASSEMBLY property: two genes were called on
the same contig. It is not physical phasing. yak trioeval gives the physical truth for the 179 trio children: for each
child gene call, which parent the surrounding k-mers came from (57b `origin`, pat or mat when one side has >= 5 k-mers and
>= 3x the other). Two genes on one child haplotype are truly in phase when they came from the same parent. So, per gene pair
(a, b), over child haplotypes where both genes are called:

  colocation_pct           % of haplotypes with both genes on one contig                        (assembly, as in 50/56)
  same_parent_all_pct      among haplotypes with both origins resolved, % from the same parent  (truth, ANY contig)
  same_parent_samecontig   the same, restricted to pairs on one contig: how often one-contig really is one parent
  same_parent_diffcontig   the same, restricted to pairs on different contigs: how much true phase the assembly split
  fraction_all/same/diff   the counts behind each

Outputs (AGGREGATE, ~/s06/results/71c):
  contig_matrix_cohort.tsv  all unrelated people with calls (the figure of 10-20): colocation only
  contig_matrix_trio.tsv    trio children: colocation + the four truth columns, per gene pair, plus platform strata
  gene_resolution.tsv       per gene: child calls, % with a resolved parent of origin

Genes are those in 57b's origin table; DRA and C4A need the extra run
  57b_trio_gene_origin.py --extra-genes DRA,C4A   (writes gene_origin_extra.INTERNAL.tsv, read here when present)

VM:  pixi run -e spechla -- python3 scripts/hla_popgen/71c_contig_truth_matrix.py
Tests: scripts/hla_popgen/tests/test_71c_contig_truth_matrix.py
"""
import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _s05_common as C  # noqa: E402

OUT = os.path.expanduser("~/s06/results/71c/")
ORIGIN = os.path.expanduser("~/s05/internal/57b/gene_origin.INTERNAL.tsv")
ORIGIN_EXTRA = os.path.expanduser("~/s05/internal/57b/gene_origin_extra.INTERNAL.tsv")
# row / column order of the 10-20 figure
GENE_ORDER = ["A", "C", "B", "MICA", "C4A", "DRA", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1", "TAP1"]


def log(m):
    C.log("[71c] " + m)


def hap_gene_contigs(calls):
    """calls with person_id, hap, gene, contig, copy_index, consensus -> {(person, hap): {gene: contig}} (copy 1, called)."""
    t = calls[(calls["copy_index"] == 1) & calls["consensus"].notna() & (calls["consensus"] != "undetermined")]
    out = {}
    for p, h, g, c in zip(t["person_id"], t["hap"].astype(str), t["gene"], t["contig"]):
        out.setdefault((p, h), {})[g] = c
    return out


def colocation_table(hg, genes=GENE_ORDER, strata=None):
    """hg: {(person, hap): {gene: contig}}; strata: {person: label}. -> long DataFrame per ordered gene pair."""
    rows = []
    labels = ["all"] + (sorted(set(strata.values())) if strata else [])
    cnt = {(lab, a, b): [0, 0] for lab in labels for a, b in itertools.permutations(genes, 2)}
    for (p, h), gc in hg.items():
        labs = ["all"] + ([strata[p]] if strata and p in strata else [])
        present = [g for g in genes if g in gc]
        for a, b in itertools.permutations(present, 2):
            for lab in labs:
                c = cnt[(lab, a, b)]
                c[0] += 1
                c[1] += int(gc[a] == gc[b])
    for (lab, a, b), (n, k) in cnt.items():
        if n:
            rows.append(dict(stratum=lab, gene_a=a, gene_b=b, n_haplotypes=n, n_same_contig=k,
                             colocation_pct=100.0 * k / n))
    return pd.DataFrame(rows)


def truth_table(origin, genes=GENE_ORDER, strata=None):
    """origin: 57b rows (child, hap, gene, contig, origin). -> per unordered-then-mirrored gene pair counts and rates."""
    o = origin[origin["gene"].isin(genes)]
    by = {}
    for ch, h, g, c, og in zip(o["child"], o["hap"], o["gene"], o["contig"], o["origin"]):
        by.setdefault((ch, h), {})[g] = (c, og)
    labels = ["all"] + (sorted(set(strata.values())) if strata else [])
    keys = ("called", "same_contig", "resolved", "resolved_same_parent", "res_same_contig",
            "res_same_contig_same_parent", "res_diff_contig", "res_diff_contig_same_parent")
    cnt = {(lab, a, b): dict.fromkeys(keys, 0) for lab in labels for a, b in itertools.permutations(genes, 2)}
    for (ch, h), gc in by.items():
        labs = ["all"] + ([strata[ch]] if strata and ch in strata else [])
        present = [g for g in genes if g in gc]
        for a, b in itertools.permutations(present, 2):
            (ca, oa), (cb, ob) = gc[a], gc[b]
            same_c = ca == cb
            res = oa in ("pat", "mat") and ob in ("pat", "mat")
            same_p = res and oa == ob
            for lab in labs:
                c = cnt[(lab, a, b)]
                c["called"] += 1
                c["same_contig"] += int(same_c)
                if res:
                    c["resolved"] += 1
                    c["resolved_same_parent"] += int(same_p)
                    if same_c:
                        c["res_same_contig"] += 1
                        c["res_same_contig_same_parent"] += int(same_p)
                    else:
                        c["res_diff_contig"] += 1
                        c["res_diff_contig_same_parent"] += int(same_p)
    rows = []
    for (lab, a, b), c in cnt.items():
        if not c["called"]:
            continue
        pct = lambda k, n: 100.0 * k / n if n else np.nan  # noqa: E731
        rows.append(dict(stratum=lab, gene_a=a, gene_b=b, n_haplotypes=c["called"],
                         colocation_pct=pct(c["same_contig"], c["called"]),
                         n_resolved=c["resolved"], same_parent_all_pct=pct(c["resolved_same_parent"], c["resolved"]),
                         n_resolved_same_contig=c["res_same_contig"],
                         same_parent_samecontig_pct=pct(c["res_same_contig_same_parent"], c["res_same_contig"]),
                         n_resolved_diff_contig=c["res_diff_contig"],
                         same_parent_diffcontig_pct=pct(c["res_diff_contig_same_parent"], c["res_diff_contig"]),
                         n_same_parent_all=c["resolved_same_parent"],
                         n_same_parent_samecontig=c["res_same_contig_same_parent"],
                         n_same_parent_diffcontig=c["res_diff_contig_same_parent"]))
    return pd.DataFrame(rows)


def gene_resolution(origin):
    rows = []
    for g, d in origin.groupby("gene"):
        res = d["origin"].isin(["pat", "mat"])
        rows.append(dict(gene=g, n_child_calls=len(d), n_resolved=int(res.sum()), resolved_pct=100.0 * res.mean()))
    return pd.DataFrame(rows)


def load_origin(paths=(ORIGIN, ORIGIN_EXTRA)):
    parts = []
    for p in paths:
        if os.path.exists(p):
            with open(p) as f:
                first = f.readline()
            parts.append(pd.read_csv(p, sep="\t", skiprows=1 if first.startswith("#") else 0,
                                     dtype={"child": str, "contig": str}))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=OUT)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    cohort = C.load_cohort(C.DEFAULT_COHORT)
    plat = dict(zip(cohort["person_id"], cohort["platform"]))
    C.TABLE1_USECOLS = C.TABLE1_USECOLS + [c for c in ("n_tied",) if c not in C.TABLE1_USECOLS]
    t1 = C.load_table1(C.DEFAULT_TABLE1)
    # cohort-wide co-location among unrelated people (the figure of 10-20)
    keep = C.greedy_unrelated(t1["person_id"].unique(), C.load_relatedness_pairs(C.DEFAULT_RELATEDNESS))
    hg = hap_gene_contigs(t1[t1["person_id"].isin(keep)])
    w = lambda df, n: C.write_tsv(df, os.path.join(a.out_dir, n))  # noqa: E731
    w(colocation_table(hg, strata=None), "contig_matrix_cohort.tsv")
    o = load_origin()
    if len(o):
        strata = {c: plat.get(c, "NA") for c in o["child"].unique()}
        w(truth_table(o, strata=strata), "contig_matrix_trio.tsv")
        w(gene_resolution(o), "gene_resolution.tsv")
    log("DONE 71c")
    return 0


if __name__ == "__main__":
    sys.exit(main())
