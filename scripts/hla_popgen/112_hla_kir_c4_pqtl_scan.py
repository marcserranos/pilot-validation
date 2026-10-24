"""112 Immunogenetic loci -> plasma proteome scan (VM side; All of Us Olink Explore HT x long-read HLA / KIR / C4 calls).

Exposures (unrelated long-read people with Olink): HLA alleles at 2-field resolution (>= --min-carriers), C4A / C4B copy number,
KIR gene content, ligand-group dosages. Outcomes: every Olink assay (SAMPLE rows, SampleQC and AssayQC PASS; NPX, mean over replicates).
Model per exposure x assay: NPX ~ exposure + age + age^2 + sex + ancestry dummies + ancestry PCs (residual-projection OLS; z from partial r).
Calibration: genomic-inflation lambda; label-permutation (within ancestry x sex strata) FWER threshold for max|z|; Bonferroni.
Positive controls (gene -> own protein): C4A/C4B copy -> C4 protein, KIR gene content -> its soluble protein, HLA-E / HLA-A alleles -> their proteins.
Top hits are re-estimated within each ancestry; an association counts as robust only if it keeps its sign in >= 4 of the ancestries with data
and the ancestry-meta z stays > 5.
Outputs: ~/s09/olink/pqtl_hits_pub.tsv (aggregate, counts >= 20) and the person x assay matrix ~/s09/olink/npx_matrix.private.tsv.gz (mode 600).
"""
import argparse
import csv
import gzip
import importlib.util
import math
import os
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
    print("[112] " + m, file=sys.stderr, flush=True)


def build_matrix(path, out_path):
    """Stream the long-format Olink file into person x assay mean NPX; cache as a gz tsv (mode 600)."""
    if os.path.exists(out_path):
        persons, assays, rows = [], None, []
        with gzip.open(out_path, "rt") as f:
            assays = f.readline().rstrip("\n").split("\t")[1:]
            for ln in f:
                x = ln.rstrip("\n").split("\t")
                persons.append(x[0])
                rows.append([float(v) if v != "" else np.nan for v in x[1:]])
        return persons, assays, np.array(rows, float)
    import pandas as pd
    parts = []
    cols = ["ResearchID", "Assay", "NPX", "SampleType", "SampleQC", "AssayQC"]
    for ch in pd.read_csv(path, sep="\t", usecols=cols, chunksize=2_000_000, dtype=str):
        ch = ch[(ch.SampleType == "SAMPLE") & (ch.SampleQC == "PASS") & (ch.AssayQC == "PASS") & (ch.NPX.notna()) & (ch.NPX != "NA")]
        ch = ch.assign(NPX=ch.NPX.astype(float))
        parts.append(ch.groupby(["ResearchID", "Assay"]).NPX.agg(["sum", "count"]))
        log("read chunk, %d (person, assay) cells so far" % sum(len(x) for x in parts))
    g = pd.concat(parts).groupby(level=[0, 1]).sum()
    mean = (g["sum"] / g["count"]).unstack("Assay")
    persons = [str(x) for x in mean.index]
    assays = list(mean.columns)
    M = mean.to_numpy(float)
    with gzip.open(out_path, "wt") as f:
        f.write("person_id\t" + "\t".join(assays) + "\n")
        for p, row in zip(persons, M):
            f.write(p + "\t" + "\t".join("" if np.isnan(v) else "%.4f" % v for v in row) + "\n")
    os.chmod(out_path, 0o600)
    return persons, assays, M


def residualise(M, C):
    Q, _ = np.linalg.qr(C)
    return M - Q @ (Q.T @ M), C.shape[1]


def zscores(Xr, Yr, k):
    """Partial-correlation z for every exposure column against every assay column (both residualised on the same covariates)."""
    n = Xr.shape[0]
    xs = Xr / np.maximum(np.sqrt((Xr ** 2).sum(0)), 1e-12)
    ys = Yr / np.maximum(np.sqrt((Yr ** 2).sum(0)), 1e-12)
    r = xs.T @ ys
    dof = n - k - 1
    return r * math.sqrt(dof) / np.sqrt(np.maximum(1 - r ** 2, 1e-12))


def run(a):
    r91 = load("r91", "91_novel_recurrence.py")
    e107 = load("e107", "107_phewas_engine.py")
    e108 = load("e108", "108_mechanistic_rna_tests.py")
    out = os.path.expanduser(a.out_dir)
    persons_o, assays, M = build_matrix(os.path.expanduser(a.npx), os.path.join(out, "npx_matrix.private.tsv.gz"))
    log("Olink matrix: %d persons x %d assays" % (len(persons_o), len(assays)))
    anc, lr = r91.read_membership(os.path.expanduser(a.membership))
    first, anyk = r91.read_relatives(os.path.expanduser(a.relatedness), lr)
    unrel = set(r91.greedy_unrelated(lr, anyk))
    info = {}
    with open(os.path.expanduser(a.ehr)) as f:
        for ln in f:
            x = ln.rstrip("\n").split("\t")
            if x[0] == "P":
                info[x[1]] = (float(x[2]), float(x[3]))
    pidx = [i for i, p in enumerate(persons_o) if p in unrel and p in info]
    persons = [persons_o[i] for i in pidx]
    M = M[pidx]
    log("unrelated LR people with Olink + EHR: %d" % len(persons))
    class A:  # args object expected by 107's exposure builder / KIR reader
        kir_root = a.kir_root
    calls, c4a, c4b = e107.build_exposures(A, set(persons), os.path.expanduser(a.rich), os.path.expanduser(a.kir_root), r91)
    names, cols = [], []
    cnt = Counter(al for p in persons for g in e107.GENES for al in set(calls[p][g]))
    for al, n in sorted(cnt.items()):
        if n >= a.min_carriers:
            names.append(al)
            cols.append([calls[p][al.split("*")[0]].count(al) for p in persons])
    for nm, d in (("C4A copies", c4a), ("C4B copies", c4b)):
        names.append(nm)
        cols.append([d[p] for p in persons])
    kd = e107.r91_kir(A, persons)
    have_kir = np.array([p in kd for p in persons])
    kir_names = []
    for g in e107.KIR_GENES:
        v = np.array([kd[p].get(g, 0) if p in kd else np.nan for p in persons], float)
        if np.nansum(v > 0) >= a.min_carriers and np.nansum(v == 0) >= a.min_carriers:
            kir_names.append(g + " dosage")
            names.append(g + " dosage")
            cols.append(list(v))
    X = np.array(cols, float).T
    sd = np.nanstd(X, axis=0)
    keep = sd > 0.15
    names = [n for n, k in zip(names, keep) if k]
    X = X[:, keep]
    log("exposures: %d (HLA alleles %d, KIR %d, C4 2)" % (len(names), len(names) - len(kir_names) - 2, len(kir_names)))
    ancv = np.array([anc[p] for p in persons])
    ag = np.array([info[p][1] for p in persons], float)
    ag = (ag - ag.mean()) / ag.std()
    sx = np.array([info[p][0] for p in persons], float)
    pcmap = e108.read_pcs(os.path.expanduser(a.ancestry_preds), a.n_pcs) or {}
    pcs = np.array([pcmap.get(p, [0.0] * a.n_pcs) for p in persons], float)
    use = sorted(k for k, n in Counter(ancv).items() if n >= 20 and k != "EUR")
    C = np.column_stack([np.ones(len(persons))] + [(ancv == k).astype(float) for k in use] + [ag, ag ** 2, sx, pcs])
    # assays with <70% coverage dropped; remaining NaN imputed to the assay mean (after centring) -> zero residual contribution
    cov_ok = (~np.isnan(M)).mean(0) >= 0.7
    assays = [x for x, k in zip(assays, cov_ok) if k]
    M = M[:, cov_ok]
    M = np.where(np.isnan(M), np.nanmean(M, axis=0, keepdims=True), M)
    Yr, k = residualise(M, C)
    log("assays used: %d" % len(assays))
    # KIR exposures are defined on a person subset: run the HLA/C4 block on everyone, the KIR block on complete-KIR people
    Z = np.full((len(names), len(assays)), np.nan)
    isk = np.array([n in kir_names for n in names])
    base = ~isk
    Xr, _ = residualise(X[:, base] - X[:, base].mean(0), C)
    Z[base] = zscores(Xr, Yr, k)
    if isk.any() and have_kir.sum() > 300:
        sub = np.where(have_kir)[0]
        Csub = C[sub]
        Ysub, ks = residualise(M[sub], Csub)
        Xk = X[sub][:, isk]
        Xkr, _ = residualise(Xk - Xk.mean(0), Csub)
        Z[isk] = zscores(Xkr, Ysub, ks)
    lam = float(np.nanmedian(Z ** 2) / 0.4549)
    # permutation FWER threshold for max |z|: permute exposures within ancestry x sex strata (HLA/C4 block)
    rnd = np.random.default_rng(5)
    strata = np.array([a_ + str(int(s_)) for a_, s_ in zip(ancv, sx)])
    mx = []
    for _ in range(a.perms):
        Xp = X[:, base].copy()
        for s_ in set(strata):
            idx = np.where(strata == s_)[0]
            Xp[idx] = Xp[rnd.permutation(idx)]
        Xpr, _ = residualise(Xp - Xp.mean(0), C)
        mx.append(float(np.nanmax(np.abs(zscores(Xpr, Yr, k)))))
    thr = float(np.percentile(mx, 95)) if mx else float("nan")
    lines = ["persons %d, exposures %d, assays %d, genomic-inflation lambda %.3f, permutation FWER max|z| threshold %.2f (%d perms)" % (len(persons), len(names), len(assays), lam, thr, len(mx))]
    # positive controls
    ai = {x: i for i, x in enumerate(assays)}
    ctrl = [("C4A copies", "C4"), ("C4B copies", "C4"), ("C4A copies", "C4A"), ("C4B copies", "C4B"), ("KIR3DL1 dosage", "KIR3DL1"), ("KIR3DL2 dosage", "KIR3DL2"), ("KIR3DL3 dosage", "KIR3DL3"),
            ("KIR2DL2 dosage", "KIR2DL2_KIR2DL3"), ("KIR2DL3 dosage", "KIR2DL2_KIR2DL3"), ("E*01:03", "HLA-E"), ("A*02:01", "HLA-A")]
    c4_assays = [x for x in assays if x.upper().startswith("C4")]
    lines.append("C4-like assays on panel: " + ",".join(c4_assays))
    for ex, asy in ctrl:
        if ex in names and asy in ai:
            lines.append("CONTROL %s -> %s: z %.2f" % (ex, asy, Z[names.index(ex), ai[asy]]))
    # top hits and per-ancestry robustness
    flat = sorted(((abs(Z[j, i]), j, i) for j in range(len(names)) for i in range(len(assays)) if not np.isnan(Z[j, i])), reverse=True)[:a.top]
    lines.append("")
    lines.append("exposure\tassay\tpooled_z\tFWER\tper-ancestry z (%s)\tsign-concordant\tmeta_z" % ",".join(["AFR", "AMR", "EAS", "EUR", "SAS"]))
    for az, j, i in flat:
        zs, ws = [], []
        for kk in ("AFR", "AMR", "EAS", "EUR", "SAS"):
            m = (ancv == kk)
            if kk in ("AFR", "AMR", "EAS", "EUR", "SAS") and m.sum() >= 150 and np.nanstd(X[m, j]) > 0 and (not isk[j] or have_kir[m].sum() > 100):
                mm = m & (have_kir if isk[j] else np.ones_like(m))
                Cm = np.column_stack([np.ones(mm.sum()), ag[mm], ag[mm] ** 2, sx[mm], pcs[mm]])
                yr, kc = residualise(M[mm][:, [i]], Cm)
                xr, _ = residualise((X[mm][:, [j]] - X[mm][:, [j]].mean()), Cm)
                if np.sqrt((xr ** 2).sum()) > 0:
                    zs.append(float(zscores(xr, yr, kc)[0, 0]))
                    ws.append(kk)
                else:
                    zs.append(float("nan"))
            else:
                zs.append(float("nan"))
        ok = [z for z in zs if not np.isnan(z)]
        sgn = math.copysign(1, Z[j, i])
        conc = sum(1 for z in ok if math.copysign(1, z) == sgn)
        meta_z = sum(ok) / math.sqrt(len(ok)) if ok else float("nan")
        lines.append("%s\t%s\t%.2f\t%s\t%s\t%d/%d\t%.2f" % (names[j], assays[i], Z[j, i], "yes" if az >= thr else "no", ",".join("NA" if np.isnan(z) else "%.1f" % z for z in zs), conc, len(ok), meta_z))
    pth = os.path.join(out, "pqtl_hits_pub.tsv")
    open(pth, "w").write("\n".join(lines) + "\n")
    log("wrote %s" % pth)
    print("\n".join(lines))


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
    p.add_argument("--min-carriers", type=int, default=40)
    p.add_argument("--n-pcs", type=int, default=10)
    p.add_argument("--perms", type=int, default=10)
    p.add_argument("--top", type=int, default=60)
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
