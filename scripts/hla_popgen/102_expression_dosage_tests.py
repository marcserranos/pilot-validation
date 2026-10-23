"""102 long-read-only features -> blood RNA expression (S09; VM side, aggregate output only).

Pre-registered tests (see sprints/S09_discovery_loop/LEDGER.md) of what short reads cannot see, using the private TPM panel from 101:
  positive controls  gene copy number / presence (KIR genes, C4A/C4B, DRB3/DRB4/DRB5, HLA-A hemizygosity) -> that gene's own log2(TPM+1)
  negative controls  the same predictors -> ACTB / GAPDH and a T-cell marker
  exploratory        activating-KIR count and total KIR gene count -> NK transcript score (covariates: ancestry, T-cell/B-cell/myeloid marker scores);
                     total C4 copies -> interferon-stimulated-gene score; HLA-F T237I / P205T carriage -> HLA-F expression
OLS with HC1 robust SEs; ancestry dummies; person-level counts 1-19 per non-reference group are reported '<20' with no estimate.
Features: KIR gene dosage from ~/pipeline_outputs_kir/<pid>/immuannot_output/hap{1,2}.gtf.gz restricted to people whose BOTH haplotypes are 'complete'
(KIR3DL3 and KIR3DL2 on the same contig); C4A/C4B, DRB3/4/5 and HLA-A copies from hla_calls_rich.tsv; HLA-F carriers from 91's carriers.private.tsv.
Outputs: <out>/expression_tests_pub.tsv (masked) and expression_tests_raw.private.tsv (exact, mode 600).
"""
import argparse
import csv
import gzip
import math
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

KIR_GENES = ["KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DS1"]
ACTIVATING = ["KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DS1"]
NK_SET = ["NKG7", "GNLY", "KLRF1", "NCAM1", "GZMB", "PRF1"]
T_SET = ["CD3E", "CD3D", "CD8A"]
B_SET = ["CD19", "MS4A1", "CD79A"]
MYE_SET = ["CD14", "LYZ", "FCGR3A"]
ISG_SET = ["ISG15", "IFI44L", "IFI44", "MX1", "MX2", "IFIT1", "IFIT3", "OAS1", "OAS3", "RSAD2"]


def log(m):
    print("[102] " + m, file=sys.stderr, flush=True)


def supp(n, floor=20):
    return "<%d" % floor if 0 < n < floor else str(n)


def ols(y, X):
    """Return (beta, se, p) for the LAST column of X (HC1 robust)."""
    n, k = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    b = XtX_inv @ X.T @ y
    r = y - X @ b
    meat = (X * r[:, None]).T @ (X * r[:, None]) * n / max(1, n - k)
    V = XtX_inv @ meat @ XtX_inv
    se = math.sqrt(max(V[-1, -1], 1e-300))
    z = b[-1] / se
    p = math.erfc(abs(z) / math.sqrt(2))
    return float(b[-1]), se, p


def read_kir(outroot, persons):
    """person -> {gene: dosage 0/1/2} for people with two complete KIR haplotypes."""
    res = {}
    for p in persons:
        dos, ok = Counter(), True
        for hap in ("hap1", "hap2"):
            path = os.path.join(outroot, p, "immuannot_output", hap + ".gtf.gz")
            if not os.path.exists(path):
                ok = False
                break
            genes_by_contig = defaultdict(set)
            with gzip.open(path, "rt") as f:
                for ln in f:
                    x = ln.split("\t")
                    if len(x) < 9 or x[2] != "gene":
                        continue
                    m = re.search(r'gene_name "([^"]+)"', x[8])
                    if m and m.group(1) in KIR_GENES:
                        genes_by_contig[x[0]].add(m.group(1))
            if not any({"KIR3DL3", "KIR3DL2"} <= g for g in genes_by_contig.values()):
                ok = False
                break
            for g in set().union(*genes_by_contig.values()):
                dos[g] += 1
        if ok:
            res[p] = dos
    return res


def read_rich_features(path, persons):
    c4a, c4b, hlaa = Counter(), Counter(), Counter()
    drb = {"DRB3": Counter(), "DRB4": Counter(), "DRB5": Counter()}
    seen_a = defaultdict(set)
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            p = r["person_id"]
            if p not in persons:
                continue
            g = r["gene"][4:] if r["gene"].startswith("HLA-") else r["gene"]
            g = re.sub(r"^(C4[AB]).*$", r"\1", g)
            if g == "C4A":
                c4a[p] += 1
            elif g == "C4B":
                c4b[p] += 1
            elif g in drb and str(r.get("copy_index", "1")) in ("1", "1.0"):
                drb[g][p] += 1
            elif g == "A" and str(r.get("copy_index", "1")) in ("1", "1.0"):
                seen_a[p].add(r["hap"])
    for p, s in seen_a.items():
        hlaa[p] = len(s)
    return c4a, c4b, hlaa, drb


def zscore(v):
    v = np.asarray(v, float)
    sd = v.std()
    return (v - v.mean()) / sd if sd > 0 else v * 0


def run(a):
    out = os.path.expanduser(a.out_dir)
    os.makedirs(out, exist_ok=True)
    tpm = {}
    with open(os.path.join(out, "rna_panel_tpm.private.tsv")) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            tpm[r["person_id"]] = {k: (math.log2(float(v) + 1) if v not in ("", None) else None) for k, v in r.items() if k != "person_id"}
    anc = {}
    with open(os.path.expanduser(a.membership)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("in_lr", "").lower() in ("true", "1"):
                anc[r["person_id"]] = (r.get("ancestry_pred") or "NA").upper()
    persons = sorted(p for p in tpm if p in anc)
    log("persons with RNA panel and LR: %d" % len(persons))
    kir = read_kir(os.path.expanduser(a.kir_root), persons)
    c4a, c4b, hlaa, drb = read_rich_features(os.path.expanduser(a.rich), set(persons))
    log("complete-KIR persons: %d" % len(kir))
    # HLA-F carriers from 91
    fcar = defaultdict(set)
    raw = list(csv.DictReader(open(os.path.join(os.path.expanduser(a.r366), "recurrence_raw.private.tsv")), delimiter="\t"))
    carriers = defaultdict(set)
    for r in csv.DictReader(open(os.path.join(os.path.expanduser(a.r366), "carriers.private.tsv")), delimiter="\t"):
        carriers[(r["gene"], r["ref_allele"], r["cds_events"])].add(r["person_id"])
    for r in raw:
        if r["gene"] == "F" and r["protein"] in ("T237I", "P205T"):
            fcar[r["protein"]] |= carriers[(r["gene"], r["ref_allele"], r["cds_events"])]
    ancs = sorted({anc[p] for p in persons} - {"EUR"})

    def score(p, genes):
        v = [tpm[p].get(g) for g in genes]
        return [x for x in v if x is not None]

    def build(plist, predictor, outcome):
        """design matrix with intercept, ancestry dummies, composition covariates (unless outcome is a lineage score) and predictor LAST."""
        X, y = [], []
        for p in plist:
            xv = predictor(p)
            yv = outcome(p)
            if xv is None or yv is None:
                continue
            row = [1.0] + [1.0 if anc[p] == k else 0.0 for k in ancs]
            row += a_cov(p)
            X.append(row + [xv])
            y.append(yv)
        return np.array(X, float), np.array(y, float)

    covmode = {"on": True}

    def a_cov(p):
        if not covmode["on"]:
            return []
        t, b, m = np.mean(score(p, T_SET)), np.mean(score(p, B_SET)), np.mean(score(p, MYE_SET))
        return [t, b, m]

    rows = []

    def test(name, plist, predictor, outcome, group_fn=None, cov=True):
        covmode["on"] = cov
        groups = Counter(group_fn(p) for p in plist) if group_fn else None
        if groups is not None:
            small = [g for g, n in groups.items() if g != 0 and 0 < n < 20]
            if small or sum(1 for g in groups if g != 0 and groups[g] >= 20) == 0:
                rows.append({"test": name, "n": len(plist), "groups": dict(groups), "beta": None, "se": None, "p": None, "note": "a non-reference group has <20 people"})
                return
        X, y = build(plist, predictor, outcome)
        if len(y) < 40 or np.linalg.matrix_rank(X) < X.shape[1]:
            rows.append({"test": name, "n": len(y), "groups": dict(groups) if groups else None, "beta": None, "se": None, "p": None, "note": "insufficient/rank deficient"})
            return
        b, se, pv = ols(y, X)
        rows.append({"test": name, "n": len(y), "groups": dict(groups) if groups else None, "beta": b, "se": se, "p": pv, "note": ""})

    kp = sorted(kir)
    for g in a.kir_controls.split(","):
        test("POS KIR dosage %s -> own TPM" % g, kp, lambda p, g=g: kir[p].get(g, 0), lambda p, g=g: tpm[p].get(g), lambda p, g=g: kir[p].get(g, 0))
        test("NEG KIR dosage %s -> CD3E (T-cell marker)" % g, kp, lambda p, g=g: kir[p].get(g, 0), lambda p: tpm[p].get("CD3E"), lambda p, g=g: kir[p].get(g, 0), cov=False)
    nk = lambda p: float(np.mean(score(p, NK_SET)))
    test("EXPL activating-KIR count -> NK score (covariates: ancestry, T/B/myeloid markers)", kp, lambda p: sum(1 for g in ACTIVATING if kir[p].get(g, 0) > 0), nk)
    test("EXPL total KIR gene copies -> NK score", kp, lambda p: sum(kir[p].values()), nk)
    test("NEG activating-KIR count -> T-cell score", kp, lambda p: sum(1 for g in ACTIVATING if kir[p].get(g, 0) > 0), lambda p: float(np.mean(score(p, T_SET))), cov=False)
    allp = persons
    c4tot = lambda p: c4a[p] + c4b[p]
    test("POS C4A copies -> C4A TPM", allp, lambda p: c4a[p], lambda p: tpm[p].get("C4A"), lambda p: c4a[p])
    test("POS C4B copies -> C4B TPM", allp, lambda p: c4b[p], lambda p: tpm[p].get("C4B"), lambda p: c4b[p])
    test("NEG C4 total copies -> GAPDH", allp, c4tot, lambda p: tpm[p].get("GAPDH"), lambda p: c4tot(p))
    test("EXPL C4 total copies -> interferon-stimulated-gene score", allp, c4tot, lambda p: float(np.mean(score(p, ISG_SET))), lambda p: c4tot(p))
    for g in ("DRB3", "DRB4", "DRB5"):
        test("POS %s copies -> own TPM" % g, allp, lambda p, g=g: drb[g][p], lambda p, g=g: tpm[p].get("HLA-" + g), lambda p, g=g: drb[g][p])
    test("POS HLA-A haplotype count (1 = hemizygous) -> HLA-A TPM", allp, lambda p: hlaa[p], lambda p: tpm[p].get("HLA-A"), lambda p: hlaa[p])
    for s, cs in sorted(fcar.items()):
        test("EXPL HLA-F %s carrier -> HLA-F TPM" % s, allp, lambda p, cs=cs: 1.0 if p in cs else 0.0, lambda p: tpm[p].get("HLA-F"), lambda p, cs=cs: 1 if p in cs else 0)
    cols = ["test", "n", "groups", "beta", "se", "p", "note"]
    for suffix, pub in (("pub", True), ("raw.private", False)):
        pth = os.path.join(out, "expression_tests_%s.tsv" % suffix)
        with open(pth, "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows:
                grp = r["groups"]
                gs = ";".join("%s=%s" % (k, supp(v) if pub else v) for k, v in sorted(grp.items())) if grp else ""
                f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (r["test"], supp(r["n"]) if pub else r["n"], gs,
                                                          "" if r["beta"] is None else "%.4f" % r["beta"], "" if r["se"] is None else "%.4f" % r["se"],
                                                          "" if r["p"] is None else "%.3g" % r["p"], r["note"]))
        if suffix != "pub":
            os.chmod(pth, 0o600)
    log("wrote %d tests to %s" % (len(rows), out))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    base = "~/workspace/vwb-aou-datasets-controlled-v9/v9"
    p.add_argument("--out-dir", default="~/s09")
    p.add_argument("--kir-root", default="~/pipeline_outputs_kir")
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--r366", default="~/s08/r366")
    p.add_argument("--kir-controls", default="KIR2DL1,KIR2DL3,KIR2DS4,KIR3DL1,KIR3DS1,KIR2DS2")
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
