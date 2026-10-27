"""114 Does HLA-B*14:02 shift blood RNA of LRPAP1 / LRP1 / SGSH (VM side)?

Tests B*14:02 (and C*08:02, B*14:01 for resolution) dosage against log2(TPM+1) of the genes in ~/s09/b1402_rna_tpm.private.tsv (extracted with 101 --genes ... --out-name b1402_rna_tpm.private.tsv)
in unrelated long-read people with blood RNA-seq: age, age^2, sex, ancestry dummies, 10 ancestry PCs; HC1 robust; per-ancestry where >= 8 carriers.
Interpretation: an RNA shift implies a transcript-level genetic effect (regulatory); no RNA shift with a plasma shift implies protein-level handling or an assay (epitope/interference) effect.
Output: ~/s09/olink/b1402_rna_pub.txt (aggregate).
"""
import argparse
import csv
import importlib.util
import math
import os
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def run(a):
    r91 = load("r91", "91_novel_recurrence.py")
    e102 = load("e102", "102_expression_dosage_tests.py")
    e107 = load("e107", "107_phewas_engine.py")
    e108 = load("e108", "108_mechanistic_rna_tests.py")
    tpm = {}
    with open(os.path.expanduser(a.tpm)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            tpm[r["person_id"]] = {k: (math.log2(float(v) + 1) if v not in ("", None) else np.nan) for k, v in r.items() if k != "person_id"}
    genes = [g for g in next(iter(tpm.values()))]
    anc, lr = r91.read_membership(os.path.expanduser(a.membership))
    first, anyk = r91.read_relatives(os.path.expanduser(a.relatedness), lr)
    unrel = set(r91.greedy_unrelated(lr, anyk))
    info = {}
    for ln in open(os.path.expanduser(a.ehr)):
        x = ln.rstrip("\n").split("\t")
        if x[0] == "P":
            info[x[1]] = (float(x[2]), float(x[3]))
    persons = sorted(p for p in tpm if p in unrel and p in info)

    class A:
        kir_root = a.kir_root
    calls, _, _ = e107.build_exposures(A, set(persons), os.path.expanduser(a.rich), os.path.expanduser(a.kir_root), r91)
    ancv = np.array([anc[p] for p in persons])
    ag = np.array([info[p][1] for p in persons], float)
    ag = (ag - ag.mean()) / ag.std()
    sx = np.array([info[p][0] for p in persons], float)
    pcmap = e108.read_pcs(os.path.expanduser(a.ancestry_preds), a.n_pcs) or {}
    pcs = np.array([pcmap.get(p, [0.0] * a.n_pcs) for p in persons], float)
    use = sorted(k for k, n in Counter(ancv).items() if n >= 20 and k != "EUR")
    C = np.column_stack([np.ones(len(persons))] + [(ancv == k).astype(float) for k in use] + [ag, ag ** 2, sx, pcs])

    def dose(al):
        g = al.split("*")[0]
        return np.array([calls[p][g].count(al) for p in persons], float)

    L = ["persons with RNA, unrelated, EHR: %d" % len(persons)]
    for al in ("B*14:02", "C*08:02", "B*14:01"):
        x = dose(al)
        L.append("%s carriers %s" % (al, r91.supp(int((x > 0).sum()))))
    L.append("gene\texposure\tbeta(log2TPM)\tse\tp\tn")
    for al in ("B*14:02", "C*08:02", "B*14:01"):
        x = dose(al)
        for g in genes:
            y = np.array([tpm[p].get(g, np.nan) for p in persons], float)
            ok = ~np.isnan(y)
            if ok.sum() < 200 or (x[ok] > 0).sum() < 8:
                continue
            X = np.column_stack([C[ok], x[ok]])
            if np.linalg.matrix_rank(X) < X.shape[1]:
                continue
            b, se, pv = e102.ols(y[ok] - y[ok].mean(), X)
            L.append("%s\t%s\t%.3f\t%.3f\t%.2g\t%d" % (g, al, b, se, pv, int(ok.sum())))
    open(os.path.expanduser(a.out), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    base = "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux"
    p.add_argument("--tpm", default="~/s09/b1402_rna_tpm.private.tsv")
    p.add_argument("--out", default="~/s09/olink/b1402_rna_pub.txt")
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--relatedness", default=base + "/relatedness/samples_relatedness.tsv")
    p.add_argument("--ancestry-preds", default=base + "/ancestry/ancestry_preds.tsv")
    p.add_argument("--kir-root", default="~/pipeline_outputs_kir")
    p.add_argument("--ehr", default="~/s09/phewas/ehr_icd3.private.tsv")
    p.add_argument("--n-pcs", type=int, default=10)
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
