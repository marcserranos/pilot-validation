"""108 mechanistic hypotheses tested against blood RNA (S09-15, S09-16; VM side, aggregate output only).

M1  HLA-B*27 misfolding -> unfolded-protein-response (UPR) module in carriers' blood RNA (hypothesised since 1999, human evidence small and conflicting).
    exposure: B*27 dose (0/1/2) from long-read 2-field calls; sub-contrast B*27:05 (AS-associated) vs other B*27 subtypes;
    outcome: UPR module score (mean z of HSPA5, DDIT3, XBP1, ATF6, ATF4, SEL1L, SYVN1, HERPUD1, EDEM1, DNAJB9, PDIA4, PDIA6, CALR, CANX), ERAP1, ERAP2, IL23 axis;
    controls: B*07 dose and B*57 dose (same locus, no misfolding), housekeeping module (negative), T-cell score (negative).
M2  HLA-B/A/C/G leader residue -21 (Met vs Thr) supplies the HLA-E-binding peptide -> HLA-E expression and NKG2A (KLRC1) tone; HLA-E R107G (E*01:01 vs E*01:03).
    exposure: diploid count of leader -21M alleles over A, B, C, G; HLA-E 107G dosage;
    outcomes: HLA-E TPM, KLRC1, KLRC1 - KLRC2 (NKG2A:NKG2C), NK score; negative: GAPDH, T-cell score.
OLS (HC1) on log2(TPM+1) with ancestry dummies (>= 20 people) and T/B/myeloid marker scores. Output: mechanistic_tests_pub.tsv / ..._raw.private.tsv (carrier groups < 20 skipped).
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
UPR = ["HSPA5", "DDIT3", "XBP1", "ATF6", "ATF4", "SEL1L", "SYVN1", "HERPUD1", "EDEM1", "DNAJB9", "PDIA4", "PDIA6", "CALR", "CANX"]
HK = ["PPIA", "RPL13A", "HPRT1", "UBC", "GUSB", "PGK1", "RPLP0", "TFRC", "RPS18", "YWHAZ"]
NK = ["NKG7", "GNLY", "KLRF1", "NCAM1", "GZMB", "PRF1"]
T_SET, B_SET, MYE_SET = ["CD3E", "CD3D", "CD8A"], ["CD19", "MS4A1", "CD79A"], ["CD14", "LYZ", "FCGR3A"]


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def leader_m21_map(prot_dir):
    """allele 2-field -> True if leader residue -21 is Met (alignment string index 3 of the leader), False if Thr/other."""
    import urllib.request
    out = {}
    for g in ("A", "B", "C", "G"):
        path = os.path.join(prot_dir, "%s_prot.txt" % g)
        if not os.path.exists(path):
            urllib.request.urlretrieve("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/%s_prot.txt" % g, path)
        seqs, order = defaultdict(str), []
        for ln in open(path):
            m = re.match(r"^ (\w+\*[\w:]+)\s+(.*)$", ln.rstrip("\n"))
            if m:
                if m.group(1) not in seqs:
                    order.append(m.group(1))
                seqs[m.group(1)] += m.group(2).replace(" ", "")
        ref = seqs[order[0]]
        votes = defaultdict(Counter)
        for n in order:
            s = "".join(r if c == "-" else c for r, c in zip(ref, seqs[n]))
            if len(s) > 3 and s[3] not in "*.":
                f = n.split(":")
                votes["%s:%s" % (f[0], f[1])][s[3]] += 1
        for k, c in votes.items():
            out[k] = c.most_common(1)[0][0] == "M"
    return out


def run(a):
    r91 = load("r91", "91_novel_recurrence.py")
    e102 = load("e102", "102_expression_dosage_tests.py")
    out = os.path.expanduser(a.out_dir)
    tpm = {}
    with open(os.path.join(out, "rna_panel_tpm.private.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            tpm[r["person_id"]] = {k: (math.log2(float(v) + 1) if v not in ("", None) else None) for k, v in r.items() if k != "person_id"}
    anc, lr = r91.read_membership(os.path.expanduser(a.membership))
    first, anyk = r91.read_relatives(os.path.expanduser(a.relatedness), lr)
    unrel = set(r91.greedy_unrelated(lr, anyk))
    persons = sorted(p for p in tpm if p in unrel)
    calls = defaultdict(lambda: defaultdict(list))
    with open(os.path.expanduser(a.rich)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            p = r["person_id"]
            g = r["gene"][4:] if r["gene"].startswith("HLA-") else r["gene"]
            if p in tpm and p in unrel and g in ("A", "B", "C", "G", "E") and str(r.get("copy_index", "1")) in ("1", "1.0"):
                c = r["consensus"]
                if c and "new" not in c and "*" in c:
                    cc = c[4:] if c.startswith("HLA-") else c
                    gg, rest = cc.split("*", 1)
                    f_ = rest.split(":")
                    if len(f_) >= 2:
                        calls[p][g].append("%s*%s:%s" % (gg, f_[0], re.sub(r"[A-Za-z]+$", "", f_[1])))
    m21 = leader_m21_map(out)
    ancs_all = sorted(set(anc[p] for p in persons) - {"EUR"})
    pers_ok = [p for p in persons if len(calls[p]["B"]) == 2]
    log = lambda m: print("[108] " + m, file=sys.stderr, flush=True)
    log("persons with RNA, unrelated, B called: %d" % len(pers_ok))

    def mean_z(genes, plist):
        mat = np.array([[tpm[p].get(g) if tpm[p].get(g) is not None else np.nan for g in genes] for p in plist], float)
        mat = (mat - np.nanmean(mat, 0)) / np.where(np.nanstd(mat, 0) > 0, np.nanstd(mat, 0), 1)
        return np.nanmean(mat, axis=1)

    rows = []

    def test(name, plist, xfun, yarr_or_gene, group_fn=None, extra=None, only=None):
        if isinstance(yarr_or_gene, str):
            y = np.array([tpm[p].get(yarr_or_gene) if tpm[p].get(yarr_or_gene) is not None else np.nan for p in plist], float)
        else:
            y = np.asarray(yarr_or_gene, float)
        x = np.array([xfun(p) for p in plist], float)
        ok = ~np.isnan(y)
        if group_fn is not None:
            gc = Counter(group_fn(p) for p, o in zip(plist, ok) if o)
            if sum(1 for g_, n in gc.items() if g_ != 0 and n >= 20) == 0 or any(0 < n < 20 for g_, n in gc.items() if g_ != 0 and n < 20 and False):
                rows.append((name, int(ok.sum()), None, None, None, "carrier group < 20")); return
            small = {g_ for g_, n in gc.items() if n < 20}
            ok &= np.array([group_fn(p) not in small for p in plist])
        if only:
            ok &= np.array([anc[p] == only for p in plist])
        pl = [p for p, o in zip(plist, ok) if o]
        ac = Counter(anc[p] for p in pl)
        use = [] if only else sorted(k for k, n in ac.items() if n >= 20 and k != "EUR")
        comp = np.array([[np.mean([tpm[p][g] for g in S if tpm[p].get(g) is not None] or [0.0]) for S in (T_SET, B_SET, MYE_SET)] for p in pl])
        cols = [np.ones(len(pl)), np.array([[1.0 if anc[p] == k else 0.0 for k in use] for p in pl]).reshape(len(pl), -1), comp]
        if extra:
            cols.append(np.array([extra(p) for p in pl], float).reshape(len(pl), -1))
        X = np.column_stack(cols + [x[ok]])
        if len(pl) < 60 or np.linalg.matrix_rank(X) < X.shape[1]:
            rows.append((name, len(pl), None, None, None, "insufficient/rank deficient")); return
        b, se, pv = e102.ols(y[ok], X)
        rows.append((name, len(pl), b, se, pv, ""))

    B27 = lambda p: calls[p]["B"].count
    d27 = lambda p: sum(1 for al in calls[p]["B"] if al.startswith("B*27:"))
    upr = mean_z(UPR, pers_ok)
    hk = mean_z(HK, pers_ok)
    nk = mean_z(NK, pers_ok)
    test("M1 B*27 dose -> UPR module score", pers_ok, d27, upr, d27)
    test("M1 B*27 dose -> ERAP1", pers_ok, d27, "ERAP1", d27)
    test("M1 B*27 dose -> ERAP2", pers_ok, d27, "ERAP2", d27)
    test("M1 B*27 dose -> IL23A", pers_ok, d27, "IL23A", d27)
    test("M1 B*27 dose -> HSPA5 (BiP)", pers_ok, d27, "HSPA5", d27)
    test("M1 B*27 dose -> XBP1", pers_ok, d27, "XBP1", d27)
    d2705 = lambda p: sum(1 for al in calls[p]["B"] if al == "B*27:05")
    test("M1 B*27:05 dose -> UPR module score", pers_ok, d2705, upr, d2705)
    d07 = lambda p: sum(1 for al in calls[p]["B"] if al == "B*07:02")
    d57 = lambda p: sum(1 for al in calls[p]["B"] if al.startswith("B*57:"))
    test("M1 NEG B*07:02 dose -> UPR module score", pers_ok, d07, upr, d07)
    test("M1 NEG B*57 dose -> UPR module score", pers_ok, d57, upr, d57)
    test("M1 NEG B*27 dose -> housekeeping module", pers_ok, d27, hk, d27)
    test("M1 NEG B*27 dose -> T-cell score", pers_ok, d27, mean_z(T_SET, pers_ok), d27)
    # M2
    pers2 = [p for p in persons if all(len(calls[p][g]) == 2 for g in ("A", "B", "C", "E"))]
    log("M2 persons: %d" % len(pers2))

    def load21(p):
        return sum(1 for g in ("A", "B", "C", "G") for al in calls[p][g] if m21.get(al) is True)

    e107g = lambda p: sum(1 for al in calls[p]["E"] if al == "E*01:03")
    nkgr = np.array([(tpm[p].get("KLRC1") or 0) - (tpm[p].get("KLRC2") or 0) for p in pers2])
    nk2 = mean_z(NK, pers2)
    test("M2 E*01:03 (107G) dose -> HLA-E TPM [positive control]", pers2, e107g, "HLA-E", e107g)
    test("M2 leader -21M load (A,B,C,G; 0-8) -> HLA-E TPM", pers2, load21, "HLA-E")
    test("M2 leader -21M load -> KLRC1 (NKG2A)", pers2, load21, "KLRC1")
    test("M2 leader -21M load -> NKG2A:NKG2C (KLRC1-KLRC2)", pers2, load21, nkgr)
    test("M2 leader -21M load -> NK score", pers2, load21, nk2)
    test("M2 E*01:03 (107G) dose -> KLRC1", pers2, e107g, "KLRC1", e107g)
    test("M2 NEG leader -21M load -> GAPDH", pers2, load21, "GAPDH")
    test("M2 NEG leader -21M load -> T-cell score", pers2, load21, mean_z(T_SET, pers2))
    # M2 robustness: age + sex (+age^2), within-ancestry, per-gene leaders, E107 adjustment
    info = {}
    ep = os.path.join(out, "phewas", "ehr_icd3.private.tsv")
    if os.path.exists(ep):
        for ln in open(ep):
            x_ = ln.rstrip("\n").split("\t")
            if x_[0] == "P":
                info[x_[1]] = (float(x_[2]), float(x_[3]))
    p2 = [p for p in pers2 if p in info]
    if len(p2) > 500:
        am = np.mean([info[p][1] for p in p2]); asd = np.std([info[p][1] for p in p2])
        ag = lambda p: ((info[p][1] - am) / asd)
        cov = lambda p: [info[p][0], ag(p), ag(p) ** 2]
        cov_e = lambda p: [info[p][0], ag(p), ag(p) ** 2, e107g(p)]
        nkgr2 = np.array([(tpm[p].get("KLRC1") or 0) - (tpm[p].get("KLRC2") or 0) for p in p2])
        nk3 = mean_z(NK, p2)
        nkidx = {p: v for p, v in zip(p2, nk3)}
        cov_nk = lambda p: [info[p][0], ag(p), ag(p) ** 2, nkidx[p]]
        for xl, xf in (("leader -21M load", load21), ("E*01:03 (107G) dose", e107g)):
            for lab, only in (("all ancestries", None), ("EUR", "EUR"), ("AFR", "AFR"), ("AMR", "AMR"), ("EAS", "EAS"), ("SAS", "SAS")):
                test("R2 %s -> KLRC1 [age+sex+age2; %s]" % (xl, lab), p2, xf, "KLRC1", extra=cov, only=only)
                test("R2 %s -> NK score [age+sex+age2; %s]" % (xl, lab), p2, xf, nk3, extra=cov, only=only)
                test("R2 %s -> KLRC1 given NK score [age+sex+age2+NK; %s]" % (xl, lab), p2, xf, "KLRC1", extra=cov_nk, only=only)
                test("R2 %s -> HLA-E TPM [age+sex+age2; %s]" % (xl, lab), p2, xf, "HLA-E", extra=cov, only=only)
            for kind in ("-> KLRC1 [age+sex+age2;", "-> NK score [age+sex+age2;", "-> KLRC1 given NK score [", "-> HLA-E TPM [age+sex+age2;"):
                bs = [(r[2], r[3], r[1]) for r in rows if r[0].startswith("R2 %s %s" % (xl, kind)) and r[2] is not None and not r[0].endswith("all ancestries]")]
                if len(bs) >= 2:
                    w = [1 / se ** 2 for _, se, _ in bs]
                    bm = sum(wi * b for wi, (b, _, _) in zip(w, bs)) / sum(w)
                    sem = 1 / math.sqrt(sum(w))
                    het = sum(wi * (b - bm) ** 2 for wi, (b, _, _) in zip(w, bs))
                    rows.append(("R2 META over %d ancestries: %s %s heterogeneity Q=%.2f)" % (len(bs), xl, kind, het), sum(n for _, _, n in bs), bm, sem, math.erfc(abs(bm / sem) / math.sqrt(2)), ""))
        test("R2 leader -21M load -> KLRC1 [age+sex+age2 + E107G]", p2, load21, "KLRC1", extra=cov_e)
        for g_ in ("A", "B", "C", "G"):
            fx = lambda p, g_=g_: sum(1 for al in calls[p][g_] if m21.get(al) is True)
            test("R2 leader -21M load (%s only) -> KLRC1 [age+sex+age2]" % g_, p2, fx, "KLRC1", extra=cov)
        test("R2 NEG leader -21M load -> NCAM1-independent T-cell score [age+sex+age2]", p2, load21, mean_z(T_SET, p2), extra=cov)
    for suffix, pub in (("pub", True), ("raw.private", False)):
        pth = os.path.join(out, "mechanistic_tests_%s.tsv" % suffix)
        with open(pth, "w") as f:
            f.write("test\tn\tbeta\tse\tp\tnote\n")
            for name, n, b, se, pv, note in rows:
                f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (name, r91.supp(n) if pub else n, "" if b is None else "%.4f" % b, "" if se is None else "%.4f" % se, "" if pv is None else "%.3g" % pv, note))
        if not pub:
            os.chmod(pth, 0o600)
    log("wrote %d tests" % len(rows))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    base = "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux"
    p.add_argument("--out-dir", default="~/s09")
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--relatedness", default=base + "/relatedness/samples_relatedness.tsv")
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
