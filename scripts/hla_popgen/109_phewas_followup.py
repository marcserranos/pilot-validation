"""109 Within-ancestry follow-up of PheWAS leads from 107 (VM side).

For each exposure -> outcome pair: per-ancestry carrier/case counts (small cells suppressed), per-ancestry score z, score-meta z across
ancestries, and the pooled test under added adjustments: (a) healthcare utilisation (log number of distinct ICD-10 3-char categories),
(b) the most correlated co-carried allele / the sibling C4 copy number. Also lists the alleles most correlated with the exposure.
A real association should survive within ancestry (same sign, meta z similar) and utilisation adjustment; a confound will collapse.
Output: <out>/followup_pub.txt (aggregate, counts <20 masked).
"""
import argparse
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def meta_score(U, I):
    """Score-test meta over strata: z = sum U / sqrt(sum I)."""
    I = np.asarray(I)
    return float(np.sum(U) / np.sqrt(max(np.sum(I), 1e-12)))


def run(a):
    e = load("e107", "107_phewas_engine.py")
    r91 = load("r91", "91_novel_recurrence.py")
    out = os.path.expanduser(a.out_dir)
    anc, lr = r91.read_membership(os.path.expanduser(a.membership))
    first, anyk = r91.read_relatives(os.path.expanduser(a.relatedness), lr)
    unrel = sorted(r91.greedy_unrelated(lr, anyk))
    ehr = e.load_ehr(a, unrel, out)
    persons = [p for p in unrel if p in ehr["info"]]
    calls, c4a, c4b = e.build_exposures(a, set(persons), os.path.expanduser(a.rich), os.path.expanduser(a.kir_root), r91)
    sex = np.array([ehr["info"][p][0] for p in persons], float)
    age = np.array([ehr["info"][p][1] for p in persons], float)
    age = (age - age.mean()) / age.std()
    util = np.log1p(np.array([len(ehr["cats"].get(p, ())) for p in persons], float))
    util = (util - util.mean()) / util.std()
    ancv = np.array([anc[p] for p in persons])
    pcs = read_pcs(os.path.expanduser(a.ancestry_preds), persons, a.n_pcs)
    if pcs is not None:
        lines.append("ancestry PCs available for %d persons (first %d used)" % (int(pcs[1].sum()), a.n_pcs))
        pcs = pcs[0]
    ancs = sorted(set(ancv))
    A = np.array([[1.0 if anc[p] == k else 0.0 for k in ancs if k != "EUR"] for p in persons])
    lines = ["n unrelated persons with EHR: %d" % len(persons)]

    def expo(name):
        if name.startswith("C4"):
            d = {"C4A": c4a, "C4B": c4b}[name[:3]]
            return np.array([d[p] for p in persons], float)
        g = name.split("*")[0]
        return np.array([calls[p][g].count(name) for p in persons], float)

    allele_names = sorted({al for p in persons for g in e.GENES for al in set(calls[p][g])})
    for exp_name, outcomes in [x.rsplit(":", 1) for x in a.pairs.split(",")]:
        x = expo(exp_name)
        lines.append("\n=== %s (carriers/alleles per ancestry; dosage mean %.3f)" % (exp_name, x.mean()))
        # correlated alleles
        if not exp_name.startswith("C4"):
            cors = []
            for al in allele_names:
                if al == exp_name:
                    continue
                v = expo(al)
                if v.std() > 0 and (v > 0).sum() >= 30:
                    cors.append((abs(np.corrcoef(x, v)[0, 1]), np.corrcoef(x, v)[0, 1], al))
            cors.sort(reverse=True)
            lines.append("most correlated alleles (pooled): " + ", ".join("%s r=%.2f" % (c[2], c[1]) for c in cors[:4]))
            partner = cors[0][2] if cors else None
        else:
            partner = None
        for oc in outcomes.split("+"):
            y = np.array([1.0 if oc in ehr["cats"].get(p, ()) else 0.0 for p in persons])
            lines.append("-- outcome %s: cases %s of %d" % (oc, r91.supp(int(y.sum())), len(y)))
            base = np.column_stack([np.ones(len(persons)), A, sex, age, age ** 2])
            variants = [("pooled (ancestry+sex+age)", base), ("+ utilisation", np.column_stack([base, util]))]
            if partner:
                variants.append(("+ utilisation + partner %s" % partner, np.column_stack([base, util, expo(partner)])))
            if pcs is not None:
                variants.append(("+ utilisation + %d PCs (pooled)" % a.n_pcs, np.column_stack([base, util, pcs])))
            if exp_name.startswith("C4A"):
                variants.append(("+ utilisation + C4B copies", np.column_stack([base, util, c4b_arr(persons, c4b)])))
            for nm, Z in variants:
                U, I = e.score_tests(x[:, None], y[:, None], Z)
                lines.append("  %-38s z %6.2f  (beta approx %.3f)" % (nm, U[0, 0] / np.sqrt(max(I[0, 0], 1e-12)), U[0, 0] / max(I[0, 0], 1e-12)))
            Us, Is, Ups, Ips = [], [], [], []
            for k in ancs:
                m = ancv == k
                if m.sum() < 200:
                    continue
                Zk = np.column_stack([np.ones(m.sum()), sex[m], age[m], age[m] ** 2, util[m]])
                Zp = np.column_stack([Zk, pcs[m]]) if pcs is not None else None
                car = int((x[m] > 0).sum())
                cs = int(y[m].sum())
                cc = int(((x[m] > 0) & (y[m] > 0)).sum())
                if car < 5 or cs < 5 or x[m].std() == 0:
                    lines.append("  %-4s n %5d carriers %s cases %s: too sparse" % (k, m.sum(), r91.supp(car), r91.supp(cs)))
                    continue
                U, I = e.score_tests(x[m][:, None], y[m][:, None], Zk)
                Us.append(U[0, 0])
                Is.append(I[0, 0])
                if Zp is not None:
                    Up, Ip = e.score_tests(x[m][:, None], y[m][:, None], Zp)
                    Ups.append(Up[0, 0])
                    Ips.append(Ip[0, 0])
                    lines.append("       with %d within-ancestry PCs: z %6.2f" % (a.n_pcs, Up[0, 0] / np.sqrt(max(Ip[0, 0], 1e-12))))
                lines.append("  %-4s n %5d carriers %5s cases %5s carrier-cases %5s  z %6.2f  beta %.3f" % (k, m.sum(), r91.supp(car), r91.supp(cs), r91.supp(cc), U[0, 0] / np.sqrt(max(I[0, 0], 1e-12)), U[0, 0] / max(I[0, 0], 1e-12)))
            if Us:
                lines.append("  within-ancestry score meta (utilisation-adjusted): z %.2f over %d ancestries" % (meta_score(Us, Is), len(Us)))
                if Ups:
                    lines.append("  within-ancestry score meta with PCs: z %.2f" % meta_score(Ups, Ips))
                if len(Us) > 1:
                    zs = [u / np.sqrt(max(i, 1e-12)) for u, i in zip(Us, Is)]
                    for j in range(len(Us)):
                        keep = [t for t in range(len(Us)) if t != j]
                        lines.append("    leave-out #%d: z %.2f" % (j + 1, meta_score([Us[t] for t in keep], [Is[t] for t in keep])))
    with open(os.path.join(out, a.out_name + "_pub.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


def read_pcs(path, persons, k):
    """First k ancestry PCs (AoU ancestry_preds.tsv pca_features), standardised; returns (matrix, has_pc mask) or None."""
    import ast
    if not os.path.exists(path):
        return None
    d = {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ci, cj = hdr.index("research_id"), hdr.index("pca_features")
        for ln in f:
            x = ln.rstrip("\n").split("\t")
            d[x[ci]] = ast.literal_eval(x[cj])[:k]
    has = np.array([p in d for p in persons])
    M = np.array([d.get(p, [0.0] * k) for p in persons], float)
    mu, sd = M[has].mean(0), M[has].std(0)
    M = np.where(has[:, None], (M - mu) / np.where(sd > 0, sd, 1), 0.0)
    return M, has


def c4b_arr(persons, c4b):
    return np.array([c4b[p] for p in persons], float)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    base = "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux"
    p.add_argument("--out-dir", default="~/s09/phewas")
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--relatedness", default=base + "/relatedness/samples_relatedness.tsv")
    p.add_argument("--kir-root", default="~/pipeline_outputs_kir")
    p.add_argument("--ancestry-preds", default=base + "/ancestry/ancestry_preds.tsv")
    p.add_argument("--n-pcs", type=int, default=10)
    p.add_argument("--out-name", default="followup_pcs")
    p.add_argument("--pairs", default="DQA1*05:03:E87+K72+K70+A41,C4A copies:I77,B*27:05:H20")
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
