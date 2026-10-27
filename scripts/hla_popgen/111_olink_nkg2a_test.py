"""111 Plasma protein test of the HLA-E 107G -> NKG2A (KLRC1) association (VM side; All of Us Olink Explore HT, v9 batch-normalised).

Reads the per-assay subset extracted from Olink_10k_aou_v9_normalized.tsv.gz (long format; SAMPLE rows, SampleQC and AssayQC PASS), joins unrelated
long-read people with a two-copy HLA-E call, and regresses plasma NPX on E*01:03 (107G) dose with age, age^2, sex, ancestry dummies and ancestry PCs;
per ancestry and as an inverse-variance meta. Specificity: other NK-receptor proteins, HLA-E (soluble), B2M, and a B*27 dose negative control.
Validity: correlation of KLRC1 NPX with KLRC1 blood-RNA in people who have both. Caveat printed with the result: plasma KLRC1 is a soluble/shed fraction, not cell-surface NKG2A.
Outputs: ~/s09/olink/olink_tests_pub.tsv (counts < 20 masked) and olink_tests_raw.private.tsv (mode 600).
"""
import argparse
import csv
import importlib.util
import math
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def log(m):
    print("[111] " + m, file=sys.stderr, flush=True)


def read_olink(path, value_col):
    """person -> {assay: mean NPX} over SAMPLE rows with passing sample and assay QC."""
    acc = defaultdict(lambda: defaultdict(list))
    with open(path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            if r["SampleType"] != "SAMPLE" or r.get("SampleQC") != "PASS" or r.get("AssayQC") != "PASS":
                continue
            v = r.get(value_col)
            if v in (None, "", "NA"):
                continue
            acc[r["ResearchID"]][r["Assay"]].append(float(v))
    return {p: {a: float(np.mean(v)) for a, v in d.items()} for p, d in acc.items()}


def meta(bs):
    """Inverse-variance fixed-effect meta of [(beta, se, n)] -> (beta, se, p, Q)."""
    w = [1 / se ** 2 for _, se, _ in bs]
    bm = sum(wi * b for wi, (b, _, _) in zip(w, bs)) / sum(w)
    sem = 1 / math.sqrt(sum(w))
    q = sum(wi * (b - bm) ** 2 for wi, (b, _, _) in zip(w, bs))
    return bm, sem, math.erfc(abs(bm / sem) / math.sqrt(2)), q


def run(a):
    r91 = load("r91", "91_novel_recurrence.py")
    e102 = load("e102", "102_expression_dosage_tests.py")
    e108 = load("e108", "108_mechanistic_rna_tests.py")
    out = os.path.expanduser(a.out_dir)
    olink = read_olink(os.path.join(out, "subset.private.tsv"), a.value_col)
    log("Olink persons with passing samples: %d" % len(olink))
    anc, lr = r91.read_membership(os.path.expanduser(a.membership))
    first, anyk = r91.read_relatives(os.path.expanduser(a.relatedness), lr)
    unrel = set(r91.greedy_unrelated(lr, anyk))
    calls = defaultdict(lambda: defaultdict(list))
    with open(os.path.expanduser(a.rich)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            p = r["person_id"]
            g = r["gene"][4:] if r["gene"].startswith("HLA-") else r["gene"]
            if p in olink and p in unrel and g in ("B", "E") and str(r.get("copy_index", "1")) in ("1", "1.0"):
                c = r["consensus"]
                if c and "new" not in c and "*" in c:
                    cc = c[4:] if c.startswith("HLA-") else c
                    gg, rest = cc.split("*", 1)
                    f_ = rest.split(":")
                    if len(f_) >= 2:
                        calls[p][g].append("%s*%s:%s" % (gg, f_[0], re.sub(r"[A-Za-z]+$", "", f_[1])))
    info = {}
    with open(os.path.expanduser(a.ehr)) as f:
        for ln in f:
            x = ln.rstrip("\n").split("\t")
            if x[0] == "P" and x[1] in olink:
                info[x[1]] = (float(x[2]), float(x[3]))
    pcmap = e108.read_pcs(os.path.expanduser(a.ancestry_preds), a.n_pcs) or {}
    persons = sorted(p for p in olink if p in unrel and len(calls[p]["E"]) == 2 and p in info)
    log("unrelated LR persons with Olink + two-copy E call + EHR: %d" % len(persons))
    e107g = np.array([sum(1 for al in calls[p]["E"] if al == "E*01:03") for p in persons], float)
    d27 = np.array([sum(1 for al in calls[p]["B"] if al.startswith("B*27:")) for p in persons], float)
    ancv = np.array([anc[p] for p in persons])
    ag = np.array([info[p][1] for p in persons], float)
    ag = (ag - ag.mean()) / ag.std()
    sx = np.array([info[p][0] for p in persons], float)
    zero = [0.0] * a.n_pcs
    pcs = np.array([pcmap.get(p, zero) for p in persons], float)
    rows = [("n persons (unrelated, LR, Olink, E called)", len(persons), None, None, None, "")]
    rows += [("n by ancestry %s" % k, n, None, None, None, "") for k, n in sorted(Counter(ancv).items())]
    rows.append(("E*01:03 dose distribution 0/1/2", len(persons), None, None, None, "/".join(str(int((e107g == v).sum())) for v in (0, 1, 2))))

    def test(label, xvec, assay, only=None):
        y = np.array([olink[p].get(assay, np.nan) for p in persons], float)
        ok = ~np.isnan(y)
        if only:
            ok &= ancv == only
        if ok.sum() < 60 or xvec[ok].std() == 0 or (xvec[ok] > 0).sum() < 8:
            rows.append((label, int(ok.sum()), None, None, None, "insufficient"))
            return None
        use = [] if only else sorted(k for k, n in Counter(ancv[ok]).items() if n >= 20 and k != "EUR")
        cols = [np.ones(ok.sum())] + [(ancv[ok] == k).astype(float) for k in use] + [ag[ok], ag[ok] ** 2, sx[ok]] + [pcs[ok][:, j] for j in range(pcs.shape[1])]
        X = np.column_stack(cols + [xvec[ok]])
        if np.linalg.matrix_rank(X) < X.shape[1]:
            rows.append((label, int(ok.sum()), None, None, None, "rank deficient"))
            return None
        b, se, pv = e102.ols(y[ok] - np.mean(y[ok]), X)
        sd = float(np.std(y[ok]))
        rows.append((label, int(ok.sum()), b, se, pv, "NPX SD %.2f" % sd))
        return b, se, int(ok.sum())

    for assay in a.assays.split(","):
        res = []
        test("107G dose -> %s NPX [all ancestries; age+sex+PCs]" % assay, e107g, assay)
        for k in ("EUR", "AFR", "AMR", "EAS", "SAS"):
            t = test("107G dose -> %s NPX [%s]" % (assay, k), e107g, assay, only=k)
            if t:
                res.append(t)
        if len(res) >= 2:
            bm, sem, pv, q = meta(res)
            rows.append(("META over %d ancestries: 107G -> %s NPX; Q=%.2f" % (len(res), assay, q), sum(n for _, _, n in res), bm, sem, pv, ""))
    test("NEG B*27 dose -> KLRC1 NPX [all ancestries]", d27, "KLRC1")
    # protein vs RNA validity
    tpm_path = os.path.expanduser(a.rna)
    if os.path.exists(tpm_path):
        rna = {}
        with open(tpm_path) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["person_id"] in olink and r.get("KLRC1") not in (None, ""):
                    rna[r["person_id"]] = math.log2(float(r["KLRC1"]) + 1)
        both = [p for p in olink if p in rna and "KLRC1" in olink[p]]
        if len(both) >= 30:
            xs = np.array([rna[p] for p in both])
            ys = np.array([olink[p]["KLRC1"] for p in both])
            rows.append(("validity: corr(KLRC1 RNA log2TPM, KLRC1 plasma NPX)", len(both), float(np.corrcoef(xs, ys)[0, 1]), None, None, "all persons with both"))
    for suffix, pub in (("pub", True), ("raw.private", False)):
        pth = os.path.join(out, "olink_tests_%s.tsv" % suffix)
        with open(pth, "w") as f:
            f.write("test\tn\tbeta\tse\tp\tnote\n")
            for name, n, b, se, pv, note in rows:
                f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (name, r91.supp(n) if pub else n, "" if b is None else "%.4f" % b, "" if se is None else "%.4f" % se, "" if pv is None else "%.3g" % pv, note))
        if not pub:
            os.chmod(pth, 0o600)
    log("wrote %d rows" % len(rows))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    base = "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux"
    p.add_argument("--out-dir", default="~/s09/olink")
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--relatedness", default=base + "/relatedness/samples_relatedness.tsv")
    p.add_argument("--ancestry-preds", default=base + "/ancestry/ancestry_preds.tsv")
    p.add_argument("--ehr", default="~/s09/phewas/ehr_icd3.private.tsv")
    p.add_argument("--rna", default="~/s09/rna_panel_tpm.private.tsv")
    p.add_argument("--n-pcs", type=int, default=10)
    p.add_argument("--value-col", default="NPX")
    p.add_argument("--assays", default="KLRC1,KLRD1,KLRK1,KLRB1,KLRF1,NCR1,NCR3,NCAM1,HLA-E,B2M")
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
