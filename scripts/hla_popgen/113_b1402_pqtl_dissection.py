"""113 Dissect the B*14:02 / C*08:02 -> LRPAP1 / LRP1 / SGSH plasma-protein signal (VM side; Olink Explore HT x long-read HLA calls).

A  haplotype tagging: B*14:02 / C*08:02 / DRB1*01:02 / DQB1*05:01 co-occurrence (counts >= 20), single-allele vs joint models, forward selection among all alleles.
B  sub-allele resolution: effect of every B*14:xx and C*08:xx 2-field allele on LRPAP1.
C  amino-acid resolution: per-residue dosage scan over mature HLA-B / HLA-C (IPD alignments) vs LRPAP1, before and after conditioning on B*14:02.
D  specificity and pathway: B*14:02 z over all assays (top list, number past the FWER threshold), a hand-picked set of LRP1-ligand / receptor / lysosomal assays,
   dose-response (0/1/2 copies, mean residual NPX in SD), correlation of LRPAP1 / LRP1 / SGSH residuals in carriers vs non-carriers.
Reads the cached person x assay matrix (112). Output: ~/s09/olink/b1402_pub.txt (aggregate only).
"""
import argparse
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
    print("[113] " + m, file=sys.stderr, flush=True)


LIGAND_ASSAYS = ["LRPAP1", "LRP1", "LRP1B", "LRP2", "LRP8", "LRP10", "LRP11", "LRP12", "LDLR", "VLDLR", "SORT1", "SORL1", "A2M", "APOE", "SERPINE1", "PLAU", "PLAUR", "MMP9", "MMP13",
                 "THBS1", "THBS2", "CCN2", "PCSK9", "APP", "F8", "F9", "F10", "IGFBP3", "TFPI", "PSAP", "GAA", "GLA", "IDUA", "IDS", "SGSH", "HEXA", "HEXB", "GBA", "CTSD", "CTSB", "NAGLU", "GUSB", "LAMP1", "LAMP2", "TPP1"]


def aligned_alleles(path):
    """2-field allele -> aligned residue string over the reference coordinate system; (ref, lead) where lead is the alignment index of mature residue 1."""
    seqs, order = defaultdict(str), []
    for ln in open(path):
        m = re.match(r"^ (\w+\*[\w:]+)\s+(.*)$", ln.rstrip("\n"))
        if m:
            if m.group(1) not in seqs:
                order.append(m.group(1))
            seqs[m.group(1)] += m.group(2).replace(" ", "")
    ref = seqs[order[0]]
    lead = ref.find("SHSMRYF") - 1
    votes = defaultdict(lambda: defaultdict(Counter))
    for n in order:
        s = "".join(r if c == "-" else c for r, c in zip(ref, seqs[n]))
        f = n.split(":")
        key = "%s:%s" % (f[0], f[1])
        for i, ch in enumerate(s):
            votes[key][i][ch] += 1
    out = {}
    for k, d in votes.items():
        out[k] = {i: c.most_common(1)[0][0] for i, c in d.items()}
    return out, lead


def run(a):
    e112 = load("e112", "112_hla_kir_c4_pqtl_scan.py")
    e107 = load("e107", "107_phewas_engine.py")
    e108 = load("e108", "108_mechanistic_rna_tests.py")
    e102 = load("e102", "102_expression_dosage_tests.py")
    r91 = load("r91", "91_novel_recurrence.py")
    out = os.path.expanduser(a.out_dir)
    persons_o, assays, M = e112.build_matrix(os.path.expanduser(a.npx), os.path.join(out, "npx_matrix.private.tsv.gz"))
    anc, lr = r91.read_membership(os.path.expanduser(a.membership))
    first, anyk = r91.read_relatives(os.path.expanduser(a.relatedness), lr)
    unrel = set(r91.greedy_unrelated(lr, anyk))
    info = {}
    for ln in open(os.path.expanduser(a.ehr)):
        x = ln.rstrip("\n").split("\t")
        if x[0] == "P":
            info[x[1]] = (float(x[2]), float(x[3]))
    pidx = [i for i, p in enumerate(persons_o) if p in unrel and p in info]
    persons = [persons_o[i] for i in pidx]
    M = M[pidx]

    class A:
        kir_root = a.kir_root
    calls, c4a, c4b = e107.build_exposures(A, set(persons), os.path.expanduser(a.rich), os.path.expanduser(a.kir_root), r91)
    ancv = np.array([anc[p] for p in persons])
    ag = np.array([info[p][1] for p in persons], float)
    ag = (ag - ag.mean()) / ag.std()
    sx = np.array([info[p][0] for p in persons], float)
    pcmap = e108.read_pcs(os.path.expanduser(a.ancestry_preds), a.n_pcs) or {}
    pcs = np.array([pcmap.get(p, [0.0] * a.n_pcs) for p in persons], float)
    use = sorted(k for k, n in Counter(ancv).items() if n >= 20 and k != "EUR")
    C = np.column_stack([np.ones(len(persons))] + [(ancv == k).astype(float) for k in use] + [ag, ag ** 2, sx, pcs])
    ai = {x: i for i, x in enumerate(assays)}
    L = ["persons %d, assays %d" % (len(persons), len(assays))]

    def dose(allele):
        g = allele.split("*")[0]
        return np.array([calls[p][g].count(allele) for p in persons], float)

    def yvec(assay):
        y = M[:, ai[assay]]
        return y

    def fit(assay, xs):
        y = yvec(assay)
        ok = ~np.isnan(y)
        X = np.column_stack([C[ok]] + [x[ok] for x in xs])
        b, se, pv = e102.ols(y[ok] - y[ok].mean(), X)
        sd = float(np.std(y[ok]))
        return b / sd, b / se, int(ok.sum())

    b14, c08, d12, q51 = dose("B*14:02"), dose("C*08:02"), dose("DRB1*01:02"), dose("DQB1*05:01")
    # A haplotype tagging
    L.append("\n== A. haplotype tagging")
    L.append("carriers: B*14:02 %s, C*08:02 %s, DRB1*01:02 %s, DQB1*05:01 %s" % tuple(r91.supp(int((x > 0).sum())) for x in (b14, c08, d12, q51)))
    for nm, (u, v) in {"B*14:02 & C*08:02": (b14, c08), "B*14:02 & DRB1*01:02": (b14, d12), "B*14:02 & DQB1*05:01": (b14, q51), "C*08:02 & DRB1*01:02": (c08, d12)}.items():
        L.append("%s: both %s, first only %s, second only %s" % (nm, r91.supp(int(((u > 0) & (v > 0)).sum())), r91.supp(int(((u > 0) & (v == 0)).sum())), r91.supp(int(((u == 0) & (v > 0)).sum()))))
    for assay in ("LRPAP1", "LRP1", "SGSH"):
        if assay not in ai:
            continue
        res = []
        for lab, xs in (("B*14:02 alone", [b14]), ("C*08:02 alone", [c08]), ("B*14:02 + C*08:02", [b14, c08]), ("B*14:02 + C*08:02 + DRB1*01:02 + DQB1*05:01", [b14, c08, d12, q51])):
            sdb, z, n = fit(assay, xs)
            res.append("%s: beta(SD) %.2f z %.1f" % (lab, sdb, z))
        L.append("%s | %s" % (assay, " | ".join(res)))
        # joint-model coefficients for the 2-allele model
        y = yvec(assay)
        ok = ~np.isnan(y)
        X = np.column_stack([C[ok], b14[ok], c08[ok]])
        bb = np.linalg.lstsq(X, y[ok] - y[ok].mean(), rcond=None)[0]
        L.append("   joint model coefficients (SD): B*14:02 %.2f, C*08:02 %.2f" % (bb[-2] / np.std(y[ok]), bb[-1] / np.std(y[ok])))
    # forward selection among all alleles with >= 30 carriers
    names, cols = [], []
    cnt = Counter(al for p in persons for g in e107.GENES for al in set(calls[p][g]))
    for al, n in sorted(cnt.items()):
        if n >= 30:
            names.append(al)
            cols.append(dose(al))
    Xall = np.array(cols).T
    L.append("forward selection of independent HLA alleles for LRPAP1 (%d candidate alleles):" % len(names))
    y = yvec("LRPAP1")
    ok = ~np.isnan(y)
    cur = C[ok]
    yv = y[ok] - y[ok].mean()
    chosen = []
    for _ in range(4):
        Q, _ = np.linalg.qr(cur)
        Xr = Xall[ok] - Q @ (Q.T @ Xall[ok])
        yr = yv - Q @ (Q.T @ yv)
        z = e112.zscores(Xr - 0, yr[:, None], cur.shape[1] - 0)[:, 0]
        j = int(np.nanargmax(np.abs(z)))
        chosen.append("%s (z %.1f)" % (names[j], z[j]))
        cur = np.column_stack([cur, Xall[ok][:, j]])
    L.append("   " + " -> ".join(chosen))
    # B sub-allele resolution
    L.append("\n== B. sub-allele resolution (LRPAP1)")
    for gene, pre in (("B", "B*14:"), ("C", "C*08:"), ("B", "B*15:"), ("C", "C*07:")):
        for al, n in sorted(cnt.items()):
            if al.startswith(pre) and n >= 15:
                sdb, z, nn = fit("LRPAP1", [dose(al)])
                L.append("%s carriers %s: beta(SD) %.2f z %.1f" % (al, r91.supp(n), sdb, z))
    # C amino acid resolution
    L.append("\n== C. amino-acid resolution (LRPAP1, mature positions; z before / after conditioning on B*14:02)")
    Q, _ = np.linalg.qr(C[ok])
    for gene in ("B", "C"):
        path = os.path.expanduser(os.path.join(a.prot_dir, "%s_prot.txt" % gene))
        if not os.path.exists(path):
            L.append("%s: alignment file missing (%s)" % (gene, path))
            continue
        al_map, lead = aligned_alleles(path)
        resid = []
        names_r, cols_r = [], []
        npos = max(max(d) for d in al_map.values())
        for pos in range(1, 183):
            idx = lead + pos - 1
            cnt_res = Counter()
            vecs = defaultdict(lambda: np.zeros(len(persons)))
            for k_i, p in enumerate(persons):
                for al in calls[p][gene]:
                    k = al.split("*")[1]
                    key = "%s*%s" % (gene, k)
                    d = al_map.get(key) or al_map.get(k)
                    if d is None:
                        continue
                    r = d.get(idx)
                    if r and r not in "*.-":
                        vecs[r][k_i] += 1
            for r, v in vecs.items():
                if 40 <= (v > 0).sum() <= len(persons) - 40:
                    names_r.append("%s%d%s" % (gene, pos, r))
                    cols_r.append(v)
        if not cols_r:
            L.append("%s: no residue dosages (alignment keys not matched)" % gene)
            continue
        X = np.array(cols_r).T
        Xr = X[ok] - Q @ (Q.T @ X[ok])
        yr = yv - Q @ (Q.T @ yv)
        z0 = e112.zscores(Xr, yr[:, None], C.shape[1])[:, 0]
        Cc = np.column_stack([C[ok], b14[ok]])
        Qc, _ = np.linalg.qr(Cc)
        Xrc = X[ok] - Qc @ (Qc.T @ X[ok])
        yrc = yv - Qc @ (Qc.T @ yv)
        z1 = e112.zscores(Xrc, yrc[:, None], Cc.shape[1])[:, 0]
        top = np.argsort(-np.abs(z0))[:8]
        L.append("HLA-%s: %d residue-dosage variables; top by |z|: " % (gene, len(names_r)) + "; ".join("%s z %.1f -> %.1f after B*14:02" % (names_r[j], z0[j], z1[j]) for j in top))
    # D specificity / pathway
    L.append("\n== D. specificity and pathway")
    Mi = np.where(np.isnan(M), np.nanmean(M, axis=0, keepdims=True), M)
    Yr, k = e112.residualise(Mi, C)
    xr, _ = e112.residualise((b14 - b14.mean())[:, None], C)
    z = e112.zscores(xr, Yr, k)[0]
    thr = a.fwer_z
    L.append("B*14:02: assays past |z| >= %.1f: %d of %d; top 25: " % (thr, int((np.abs(z) >= thr).sum()), len(z)) + "; ".join("%s %.1f" % (assays[i], z[i]) for i in np.argsort(-np.abs(z))[:25]))
    xr2, _ = e112.residualise((c08 - c08.mean())[:, None], C)
    z2 = e112.zscores(xr2, Yr, k)[0]
    L.append("C*08:02: assays past %.1f: %d; top 12: " % (thr, int((np.abs(z2) >= thr).sum())) + "; ".join("%s %.1f" % (assays[i], z2[i]) for i in np.argsort(-np.abs(z2))[:12]))
    L.append("LRP1-ligand / receptor / lysosomal panel (B*14:02 z): " + "; ".join("%s %.1f" % (x, z[ai[x]]) for x in LIGAND_ASSAYS if x in ai))
    L.append("ligand assays not on panel: " + ",".join(x for x in LIGAND_ASSAYS if x not in ai))
    # dose response in SD of residual NPX
    for assay in ("LRPAP1", "LRP1", "SGSH"):
        if assay not in ai:
            continue
        yy = Yr[:, ai[assay]]
        sd = yy.std()
        parts = []
        for dv in (0, 1, 2):
            m = b14 == dv
            parts.append("%d copies: n %s mean %.2f SD" % (dv, r91.supp(int(m.sum())), yy[m].mean() / sd if m.sum() >= 20 else float("nan")))
        L.append("dose-response %s: " % assay + "; ".join(parts))
    carr = b14 > 0
    for pair in (("LRPAP1", "LRP1"), ("LRPAP1", "SGSH"), ("LRP1", "SGSH")):
        if pair[0] in ai and pair[1] in ai:
            c1 = float(np.corrcoef(Yr[carr, ai[pair[0]]], Yr[carr, ai[pair[1]]])[0, 1])
            c0 = float(np.corrcoef(Yr[~carr, ai[pair[0]]], Yr[~carr, ai[pair[1]]])[0, 1])
            L.append("residual correlation %s~%s: carriers %.2f vs non-carriers %.2f" % (pair[0], pair[1], c1, c0))
    # per-ancestry carriers and effect sizes
    for k_ in ("AFR", "AMR", "EAS", "EUR", "SAS"):
        m = ancv == k_
        L.append("%s: n %d, B*14:02 carriers %s, mean LRPAP1 residual (SD) carriers %.2f vs non %.2f" % (k_, m.sum(), r91.supp(int((b14[m] > 0).sum())),
                 Yr[m & (b14 > 0), ai["LRPAP1"]].mean() / Yr[:, ai["LRPAP1"]].std() if (m & (b14 > 0)).sum() >= 20 else float("nan"), Yr[m & (b14 == 0), ai["LRPAP1"]].mean() / Yr[:, ai["LRPAP1"]].std()))
    pth = os.path.join(out, "b1402_pub.txt")
    open(pth, "w").write("\n".join(L) + "\n")
    log("wrote " + pth)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    base = "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux"
    p.add_argument("--out-dir", default="~/s09/olink")
    p.add_argument("--npx", default="~/workspace/vwb-aou-datasets-controlled-v9/v9/multiomics/proteomics/npx/batch_normalized/Olink_10k_aou_v9_normalized.tsv.gz")
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--relatedness", default=base + "/relatedness/samples_relatedness.tsv")
    p.add_argument("--ancestry-preds", default=base + "/ancestry/ancestry_preds.tsv")
    p.add_argument("--kir-root", default="~/pipeline_outputs_kir")
    p.add_argument("--ehr", default="~/s09/phewas/ehr_icd3.private.tsv")
    p.add_argument("--prot-dir", default="~/s09/phewas")
    p.add_argument("--n-pcs", type=int, default=10)
    p.add_argument("--fwer-z", type=float, default=5.91)
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
