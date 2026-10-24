"""107 hypothesis-free HLA/KIR/C4 x ICD-10 phenome scan with controls (S09-14; VM side; person-level tables stay on the VM).

Exposures (long-read resolved, unrelated people, one analysis per set):
  hla   dosage (0/1/2) of every 2-field HLA allele with >= --min-carriers carriers (A,B,C,E,F,G,DRB1,DQA1,DQB1,DPA1,DPB1), + C1/C2 and Bw4 ligand dosages
  kir   (people with two complete KIR haplotypes) presence dosage of 15 KIR genes, activating-gene count, B-content score and KIR x ligand compound genotypes
        (KIR2DL1 x C2, KIR2DL2/3 x C1, KIR3DL1 x Bw4, KIR3DS1 x Bw4)
  c4    C4A / C4B / total copy number
Outcomes: every ICD-10-CM 3-character category with >= 2 distinct dated records in >= --min-cases people (BigQuery condition_occurrence, source concept).
Test: logistic SCORE test per exposure x outcome, null model = ancestry + sex + age (+ age^2); effect = U/I (approx. log odds ratio), z = U/sqrt(I).
Controls built in: positive-control pairs (B*27-M45, DRB1*15:01-G35, DQB1*02:01/DRB1*03:01-E10, C*06:02-L40, ...) are reported; negative-control outcomes
(chapters S, T, V, W, X, Y, Z) must show no FDR hits and genomic inflation near 1; max|z| permutation (rows of X shuffled within ancestry x sex) gives FWER thresholds.
Outputs (--out-dir): phewas_hits_<set>.tsv (exposure, outcome, carriers, cases, OR, z, FDR; all counts >= 20 by construction), phewas_controls_<set>.txt, phewas_lambda.txt.
Person-level EHR pulls are cached in <out>/ehr_icd3.private.tsv (mode 600, never committed).
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
GENES = ["A", "B", "C", "E", "F", "G", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
KIR_GENES = ["KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DS1"]
ACTIVATING = ["KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DS1"]
B_SPECIFIC = ["KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS5", "KIR2DL2", "KIR2DL5A", "KIR2DL5B", "KIR3DS1"]
NEG_CHAPTERS = set("STVWXYZ")
POSITIVE = [("B*27:05", "M45"), ("DRB1*15:01", "G35"), ("DQB1*02:01", "E10"), ("DRB1*03:01", "E10"), ("DRB1*04:01", "E10"), ("C*06:02", "L40"), ("DQB1*03:02", "E10"), ("A*31:01", "M05")]


def log(m):
    print("[107] " + m, file=sys.stderr, flush=True)


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ---------------------------------------------------------------- statistics
def null_logistic(y, Z, iters=25):
    beta = np.zeros(Z.shape[1])
    beta[0] = math.log(max(y.mean(), 1e-6) / max(1 - y.mean(), 1e-6))
    for _ in range(iters):
        eta = Z @ beta
        p = 1 / (1 + np.exp(-eta))
        w = np.clip(p * (1 - p), 1e-9, None)
        z = eta + (y - p) / w
        A = Z.T @ (Z * w[:, None]) + 1e-6 * np.eye(Z.shape[1])
        nb = np.linalg.solve(A, Z.T @ (w * z))
        if np.max(np.abs(nb - beta)) < 1e-8:
            beta = nb
            break
        beta = nb
    p = 1 / (1 + np.exp(-(Z @ beta)))
    return p


def score_tests(X, Y, Z):
    """X: n x J exposures, Y: n x K binary outcomes, Z: n x q null covariates (with intercept). Returns U, I (J x K)."""
    X = X - X.mean(axis=0, keepdims=True)           # centre: avoids cancellation when an exposure is nearly constant (common alleles)
    n, J = X.shape
    K = Y.shape[1]
    U = np.zeros((J, K))
    I = np.zeros((J, K))
    X2 = X ** 2
    for k in range(K):
        y = Y[:, k]
        p = null_logistic(y, Z)
        w = np.clip(p * (1 - p), 1e-9, None)
        U[:, k] = X.T @ (y - p)
        A_inv = np.linalg.inv(Z.T @ (Z * w[:, None]) + 1e-9 * np.eye(Z.shape[1]))
        B = Z.T @ (X * w[:, None])                       # q x J
        I[:, k] = X2.T @ w - np.einsum("qj,qr,rj->j", B, A_inv, B)
    return U, I


def zmat(U, I):
    return U / np.sqrt(np.clip(I, 1e-12, None))


def bh(p):
    p = np.asarray(p)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1)
    return out


def norm_sf2(z):
    return np.array([math.erfc(abs(v) / math.sqrt(2)) for v in np.ravel(z)]).reshape(np.shape(z))


# ---------------------------------------------------------------- data
def ligand_maps(prot_dir):
    """allele 2-field -> group from IPD protein alignments: C1/C2 (C residue 80 N/K), Bw4/Bw6 (A,B residue 80 I,T vs N)."""
    import urllib.request
    maps = {}
    for g in ("A", "B", "C"):
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
        lead = ref.find("SHSMRYF") - 1
        idx = lead + 79                                    # mature residue 80 (0-based)
        votes = defaultdict(Counter)
        for n in order:
            s = "".join(r if c == "-" else c for r, c in zip(ref, seqs[n]))
            if idx < len(s) and s[idx] not in "*.":
                f = n.split(":")
                votes["%s:%s" % (f[0], f[1])][s[idx]] += 1
        for k, c in votes.items():
            r80 = c.most_common(1)[0][0]
            if g == "C":
                maps[k] = "C1" if r80 == "N" else ("C2" if r80 == "K" else None)
            else:
                maps[k] = "Bw4" if r80 in "IT" else ("Bw6" if r80 == "N" else None)
    return maps


def build_exposures(a, persons, rich_path, kir_root, r91, hla_only=False):
    calls = defaultdict(lambda: defaultdict(list))
    c4a, c4b = Counter(), Counter()
    with open(rich_path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            p = r["person_id"]
            if p not in persons:
                continue
            g = r["gene"][4:] if r["gene"].startswith("HLA-") else r["gene"]
            g2 = re.sub(r"^(C4[AB]).*$", r"\1", g)
            if g2 == "C4A":
                c4a[p] += 1
            elif g2 == "C4B":
                c4b[p] += 1
            elif g in GENES and str(r.get("copy_index", "1")) in ("1", "1.0"):
                c = r["consensus"]
                if c and "new" not in c and "*" in c:
                    cc = c[4:] if c.startswith("HLA-") else c
                    gg, rest = cc.split("*", 1)
                    f_ = rest.split(":")
                    if len(f_) >= 2:
                        calls[p][g].append("%s*%s:%s" % (gg, f_[0], re.sub(r"[A-Za-z]+$", "", f_[1])))
    return calls, c4a, c4b


def run(a):
    out = os.path.expanduser(a.out_dir)
    os.makedirs(out, exist_ok=True)
    r91 = load("r91", "91_novel_recurrence.py")
    anc, lr = r91.read_membership(os.path.expanduser(a.membership))
    first, anyk = r91.read_relatives(os.path.expanduser(a.relatedness), lr)
    unrel = sorted(r91.greedy_unrelated(lr, anyk))
    log("unrelated LR persons: %d" % len(unrel))
    ehr = load_ehr(a, unrel, out)
    persons = [p for p in unrel if p in ehr["info"]]
    log("with EHR person info: %d" % len(persons))
    calls, c4a, c4b = build_exposures(a, set(persons), os.path.expanduser(a.rich), os.path.expanduser(a.kir_root), r91)
    lig = ligand_maps(out)
    # covariates
    sex = np.array([ehr["info"][p][0] for p in persons], float)
    age = np.array([ehr["info"][p][1] for p in persons], float)
    age = (age - age.mean()) / age.std()
    ancs = sorted(set(anc[p] for p in persons) - {"EUR"})
    A = np.array([[1.0 if anc[p] == k else 0.0 for k in ancs] for p in persons])
    Z = np.column_stack([np.ones(len(persons)), A, sex, age, age ** 2])
    strata = np.array([anc[p] + str(int(s)) for p, s in zip(persons, sex)])
    # outcomes
    cnt = Counter(c for p in persons for c in ehr["cats"].get(p, ()))
    cats = sorted(c for c, n in cnt.items() if n >= a.min_cases)
    Y = np.array([[1.0 if c in ehr["cats"].get(p, ()) else 0.0 for c in cats] for p in persons])
    log("outcome categories: %d" % len(cats))
    sets = {}
    # HLA
    alleles = Counter(al for p in persons for g in GENES for al in set(calls[p][g]))
    names = sorted(al for al, n in alleles.items() if n >= a.min_carriers)
    Xh = np.array([[calls[p][al.split("*")[0]].count(al) for al in names] for p in persons], float)
    lig_names, lig_cols = [], []
    for grp in ("C1", "C2", "Bw4"):
        genes_ = ("C",) if grp in ("C1", "C2") else ("A", "B")
        lig_cols.append([sum(1 for g in genes_ for al in calls[p][g] if lig.get(al) == grp) for p in persons])
        lig_names.append("LIGAND " + grp)
    sets["hla"] = (names + lig_names, np.column_stack([Xh] + [np.array(c, float) for c in lig_cols]))
    # C4
    tot = np.array([c4a[p] + c4b[p] for p in persons], float)
    sets["c4"] = (["C4A copies", "C4B copies", "C4 total"], np.column_stack([[c4a[p] for p in persons], [c4b[p] for p in persons], tot]).astype(float))
    # KIR
    kd = r91_kir(a, persons)
    ksub = [i for i, p in enumerate(persons) if p in kd]
    log("complete-KIR persons: %d" % len(ksub))
    if len(ksub) > 500:
        Pk = [persons[i] for i in ksub]
        cols, nm = [], []
        for g in KIR_GENES:
            d = np.array([kd[p].get(g, 0) for p in Pk], float)
            if (d > 0).sum() >= a.min_carriers and (d == 0).sum() >= a.min_carriers:
                cols.append(d)
                nm.append(g + " dosage")
        act = np.array([sum(1 for g in ACTIVATING if kd[p].get(g, 0) > 0) for p in Pk], float)
        bc = np.array([sum(1 for g in B_SPECIFIC if kd[p].get(g, 0) > 0) for p in Pk], float)
        cols += [act, bc]
        nm += ["KIR activating-gene count", "KIR B-content score"]
        ligc = {grp: np.array([sum(1 for g in (("C",) if grp in ("C1", "C2") else ("A", "B")) for al in calls[p][g] if lig.get(al) == grp) for p in Pk], float) for grp in ("C1", "C2", "Bw4")}
        for kg, grp in (("KIR2DL1", "C2"), ("KIR2DL2", "C1"), ("KIR2DL3", "C1"), ("KIR3DL1", "Bw4"), ("KIR3DS1", "Bw4")):
            pres = np.array([1.0 if kd[p].get(kg, 0) > 0 else 0.0 for p in Pk])
            cols.append(pres * ligc[grp])
            nm.append("%s x %s dose" % (kg, grp))
        sets["kir"] = (nm, np.column_stack(cols), np.array(ksub))
    rnd = np.random.default_rng(13)
    lam_lines = []
    for sname, spec in sets.items():
        nm, X = spec[0], spec[1]
        sd = X.std(axis=0)
        okc = sd > 0.15
        nm, X = [n_ for n_, o in zip(nm, okc) if o], X[:, okc]
        rows = np.arange(len(persons)) if len(spec) == 2 else spec[2]
        Xs, Ys, Zs, ss = X, Y[rows], Z[rows], strata[rows]
        keep = Ys.sum(0) >= a.min_cases
        cats_s = [c for c, k in zip(cats, keep) if k]
        Ys = Ys[:, keep]
        log("%s: %d exposures x %d outcomes (n=%d)" % (sname, Xs.shape[1], Ys.shape[1], Xs.shape[0]))
        U, I = score_tests(Xs, Ys, Zs)
        zz = zmat(U, I)
        pv = norm_sf2(zz)
        neg = np.array([c[0] in NEG_CHAPTERS for c in cats_s])
        lam = np.median(zz ** 2) / 0.4549
        lam_neg = np.median(zz[:, neg] ** 2) / 0.4549 if neg.any() else float("nan")
        lam_lines.append("%s: genomic inflation lambda (all outcomes) %.3f; (negative-control chapters S,T,V,W,X,Y,Z) %.3f" % (sname, lam, lam_neg))
        # permutation FWER threshold
        mx = []
        for _ in range(a.perms):
            Xp = Xs.copy()
            for s in set(ss):
                idx = np.where(ss == s)[0]
                Xp[idx] = Xs[rnd.permutation(idx)]
            Up, Ip = score_tests(Xp, Ys, Zs)
            zp = np.abs(zmat(Up, Ip))
            zp[Ip <= 0.5] = 0.0
            mx.append(zp.max())
        thr = float(np.percentile(mx, 95)) if mx else float("nan")
        lam_lines.append("%s: permutation 95%% max|z| threshold (FWER 0.05) = %.2f (%d perms)" % (sname, thr, len(mx)))
        fdr = bh(pv.ravel()).reshape(pv.shape)
        flat = sorted(((abs(zz[j, k]), j, k) for j in range(zz.shape[0]) for k in range(zz.shape[1]) if I[j, k] > 0.5), reverse=True)
        with open(os.path.join(out, "phewas_hits_%s.tsv" % sname), "w") as f:
            f.write("exposure\toutcome_icd10_3\tcarriers\tcases\tapprox_OR\tz\tp\tFDR\tpasses_perm_FWER\tnegative_control_chapter\n")
            for az, j, k in flat[:300]:
                carriers = int((Xs[:, j] > 0).sum())
                cases = int(Ys[:, k].sum())
                beta = float(np.clip(U[j, k] / max(I[j, k], 1e-12), -20, 20))
                f.write("%s\t%s\t%d\t%d\t%.3f\t%.2f\t%.3g\t%.3g\t%s\t%s\n" % (nm[j], cats_s[k], carriers, cases, math.exp(beta), zz[j, k], pv[j, k], fdr[j, k],
                                                                           "yes" if az >= thr else "no", "yes" if cats_s[k][0] in NEG_CHAPTERS else "no"))
        with open(os.path.join(out, "phewas_controls_%s.txt" % sname), "w") as f:
            for ex, oc in POSITIVE:
                if ex in nm and oc in cats_s:
                    j, k = nm.index(ex), cats_s.index(oc)
                    f.write("POS %s -> %s: z=%.2f, approx OR %.2f, carriers %d, cases %d\n" % (ex, oc, zz[j, k], math.exp(float(np.clip(U[j, k] / max(I[j, k], 1e-12), -20, 20))), int((Xs[:, j] > 0).sum()), int(Ys[:, k].sum())))
            f.write("NEG-chapter outcomes with FDR<0.05: %d of %d tests\n" % (int(((fdr < 0.05) & neg[None, :]).sum()), int(neg.sum() * zz.shape[0])))
    with open(os.path.join(out, "phewas_lambda.txt"), "w") as f:
        f.write("\n".join(lam_lines) + "\n")
    log("done")


def load_ehr(a, persons, out):
    pth = os.path.join(out, "ehr_icd3.private.tsv")
    if os.path.exists(pth):
        info, cats = {}, defaultdict(set)
        for ln in open(pth):
            x = ln.rstrip("\n").split("\t")
            if x[0] == "P":
                info[x[1]] = (float(x[2]), float(x[3]))
            else:
                cats[x[1]].add(x[2])
        return {"info": info, "cats": cats}
    from google.cloud import bigquery
    cdr = os.environ["WORKSPACE_CDR"]
    client = bigquery.Client()
    ids = [int(p) for p in persons]
    cfg = bigquery.QueryJobConfig(query_parameters=[bigquery.ArrayQueryParameter("ids", "INT64", ids)])
    people = client.query("SELECT person_id, CASE WHEN sex_at_birth_concept_id=45878463 THEN 1 ELSE 0 END AS female, "
                          "2026 - year_of_birth AS age FROM `%s.person` WHERE person_id IN UNNEST(@ids)" % cdr, job_config=cfg).to_dataframe()
    cond = client.query("SELECT co.person_id, SUBSTR(cs.concept_code,1,3) AS icd3 FROM `%s.condition_occurrence` co JOIN `%s.concept` cs ON cs.concept_id=co.condition_source_concept_id "
                        "WHERE cs.vocabulary_id='ICD10CM' AND co.person_id IN UNNEST(@ids) GROUP BY co.person_id, icd3 HAVING COUNT(DISTINCT co.condition_start_date) >= 2" % (cdr, cdr),
                        job_config=cfg).to_dataframe()
    with open(pth, "w") as f:
        for r in people.itertuples():
            f.write("P\t%d\t%d\t%d\n" % (r.person_id, r.female, r.age))
        for r in cond.itertuples():
            f.write("C\t%d\t%s\n" % (r.person_id, r.icd3))
    os.chmod(pth, 0o600)
    return load_ehr(a, persons, out)


def r91_kir(a, persons):
    e = load("e102", "102_expression_dosage_tests.py")
    return e.read_kir(os.path.expanduser(a.kir_root), persons)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    base = "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux"
    p.add_argument("--out-dir", default="~/s09/phewas")
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--relatedness", default=base + "/relatedness/samples_relatedness.tsv")
    p.add_argument("--kir-root", default="~/pipeline_outputs_kir")
    p.add_argument("--min-carriers", type=int, default=60)
    p.add_argument("--min-cases", type=int, default=30)
    p.add_argument("--perms", type=int, default=10)
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
