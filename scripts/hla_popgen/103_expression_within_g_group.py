"""103 expression differences between alleles identical over the peptide-binding exons (S09-11; VM side, aggregate output only).

Question only long reads can ask: among people whose two HLA alleles lie in the SAME G group (identical exon 2+3 class I / exon 2 class II, so every
exon-based typer sees one allele), does the 3-field allele identity (differences in exon 4+, introns, UTRs) still predict that gene's blood RNA level?
Per gene: restrict to G-homozygous persons; person class = sorted pair of 3-field names; keep classes with >= 20 persons inside a G group with >= 2 such classes.
y = log2(TPM+1) residualised on ancestry + T/B/myeloid marker scores. Statistics:
  between-G   sum over people of (G-group mean - grand mean)^2 / total SS        (positive control: known allele-level expression differences)
  within-G    sum over G groups of between-class SS / total SS, with a label-permutation p inside each G group (B permutations)
Caveat stated up front: RNA is mapped to GRCh38, so haplotypes that diverge from the reference at UTR/exon-4 sites may look lower for mapping reasons;
a positive within-G signal is an association with mapping-bias risk, not yet a regulatory claim (follow-up: personal-reference quantification, 53/96).
Output: <out>/within_g_expression_pub.tsv (counts < 20 suppressed) and ..._raw.private.tsv (mode 600).
"""
import argparse
import csv
import importlib.util
import math
import os
import random
import sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
T_SET, B_SET, MYE_SET = ["CD3E", "CD3D", "CD8A"], ["CD19", "MS4A1", "CD79A"], ["CD14", "LYZ", "FCGR3A"]


def log(m):
    print("[103] " + m, file=sys.stderr, flush=True)


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def three_field(a):
    f = a.split(":")
    return ":".join(f[:3])


def residualise(y, X):
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ b


def within_stat(res, groups):
    """groups: dict G -> dict class -> list of idx. Returns (between-class SS inside each G, summed)."""
    tot = 0.0
    for G, cls in groups.items():
        allidx = [i for v in cls.values() for i in v]
        gm = res[allidx].mean()
        tot += sum(len(v) * (res[v].mean() - gm) ** 2 for v in cls.values())
    return tot


def run(a):
    n100 = load("n100", "100_native_vs_lr.py")
    pm = n100.load_gmap(os.path.expanduser(a.gfile))
    out = os.path.expanduser(a.out_dir)
    tpm = {}
    with open(os.path.join(out, "rna_panel_tpm.private.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            tpm[r["person_id"]] = {k: (math.log2(float(v) + 1) if v not in ("", None) else None) for k, v in r.items() if k != "person_id"}
    anc = {}
    with open(os.path.expanduser(a.membership)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("in_lr", "").lower() in ("true", "1"):
                anc[r["person_id"]] = (r.get("ancestry_pred") or "NA").upper()
    calls = defaultdict(lambda: defaultdict(dict))        # gene -> person -> hap -> full name
    with open(os.path.expanduser(a.rich)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            g = r["gene"][4:] if r["gene"].startswith("HLA-") else r["gene"]
            p = r["person_id"]
            if g in GENES and p in tpm and p in anc and str(r.get("copy_index", "1")) in ("1", "1.0"):
                nm = n100.full_name(r["consensus"])
                if nm:
                    calls[g][p][r["hap"]] = nm
    rows = []
    rnd = random.Random(5)

    def comp(p):
        t = tpm[p]
        v = lambda S: np.mean([t[g] for g in S if t.get(g) is not None] or [0.0])
        return [v(T_SET), v(B_SET), v(MYE_SET)]

    for g in GENES:
        col = "HLA-" + g
        people = []
        for p, d in calls[g].items():
            if len(d) != 2 or tpm[p].get(col) is None:
                continue
            Gs = []
            for nm in d.values():
                s = pm.get(nm)
                Gs.append(next(iter(s)) if s and len(s) == 1 else None)
            if None in Gs or Gs[0] != Gs[1]:
                continue
            people.append((p, Gs[0], "|".join(sorted(three_field(n) for n in d.values()))))
        if len(people) < 200:
            rows.append({"gene": g, "n_people": len(people), "note": "too few G-homozygous people"})
            continue
        ancl = sorted(k for k, n in Counter(anc[p] for p, _, _ in people).items() if n >= 20 and k != "EUR")
        X = np.array([[1.0] + [1.0 if anc[p] == k else 0.0 for k in ancl] + comp(p) for p, _, _ in people])
        y = np.array([tpm[p][col] for p, _, _ in people])
        res = residualise(y, X)
        tss = float((res ** 2).sum())
        # between-G (positive control)
        byG = defaultdict(list)
        for i, (_, G, _) in enumerate(people):
            byG[G].append(i)
        okG = {G: v for G, v in byG.items() if len(v) >= 20}
        gm = res.mean()
        between = sum(len(v) * (res[v].mean() - gm) ** 2 for v in okG.values())
        # within-G classes
        groups = {}
        for G, idx in okG.items():
            cls = defaultdict(list)
            for i in idx:
                cls[people[i][2]].append(i)
            cls = {c: v for c, v in cls.items() if len(v) >= 20}
            if len(cls) >= 2:
                groups[G] = cls
        if not groups:
            rows.append({"gene": g, "n_people": len(people), "n_G_groups": len(okG), "between_G_R2": between / tss, "note": "no G group with >= 2 three-field classes of >= 20 people"})
            continue
        obs = within_stat(res, groups)
        B, ge = a.perms, 0
        for _ in range(B):
            perm = {}
            for G, cls in groups.items():
                allidx = [i for v in cls.values() for i in v]
                labels = [c for c, v in cls.items() for _ in v]
                rnd.shuffle(labels)
                nc = defaultdict(list)
                for i, c in zip(allidx, labels):
                    nc[c].append(i)
                perm[G] = nc
            if within_stat(res, perm) >= obs:
                ge += 1
        n_in = sum(len(v) for cls in groups.values() for v in cls.values())
        rows.append({"gene": g, "n_people": len(people), "n_G_groups": len(okG), "n_testable_G": len(groups), "n_people_testable": n_in,
                     "between_G_R2": between / tss, "within_G_R2": obs / tss, "perm_p": (1 + ge) / (1 + B), "note": ""})
        log("%s: %d G-hom people, within-G R2=%.4f perm p=%.4f" % (g, len(people), obs / tss, (1 + ge) / (1 + B)))
    cols = ["gene", "n_people", "n_G_groups", "n_testable_G", "n_people_testable", "between_G_R2", "within_G_R2", "perm_p", "note"]
    for suffix, fm in (("pub", n100.supp), ("raw.private", str)):
        pth = os.path.join(out, "within_g_expression_%s.tsv" % suffix)
        with open(pth, "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows:
                cells = []
                for c in cols:
                    v = r.get(c, "")
                    if c in ("n_people", "n_G_groups", "n_testable_G", "n_people_testable") and v != "":
                        v = fm(int(v))
                    elif isinstance(v, float):
                        v = "%.4f" % v
                    cells.append(str(v))
                f.write("\t".join(cells) + "\n")
        if suffix != "pub":
            os.chmod(pth, 0o600)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--out-dir", default="~/s09")
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--gfile", default="~/s09/hla_nom_g.txt")
    p.add_argument("--perms", type=int, default=500)
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
