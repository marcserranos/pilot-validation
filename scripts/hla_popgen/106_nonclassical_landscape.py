"""106 diversity landscape of HLA-E/F/G vs the classical genes, by ancestry (S09-13; VM side, aggregate output only).

Per gene x ancestry, on UNRELATED long-read people, rarefied to the same number of people so ancestries and genes are comparable:
  distinct 2-field alleles (mean of --reps random subsamples), expected heterozygosity at 2-field (1 - sum p^2), top-allele frequency,
  share of haplotypes whose call contains 'new' (protein-altering novelty against the IPD release Immuannot used), named null/questionable alleles (suffix N/Q/S/L),
  and the haplotype share carrying each recurrent unnamed HLA-F/E/G signature from 91 (T237I, P205T ...).
Counts < 20 -> '<20' in the public table. Output: <out>/landscape_pub.tsv, landscape_raw.private.tsv (mode 600).
"""
import argparse
import csv
import importlib.util
import os
import random
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
GENES = ["A", "B", "C", "E", "F", "G", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]


def supp(n, floor=20):
    return "<%d" % floor if 0 < n < floor else str(n)


def two_field(c):
    if not c or "*" not in c:
        return None
    c = c[4:] if c.startswith("HLA-") else c
    g, rest = c.split("*", 1)
    f = rest.split(":")
    if len(f) < 2 or "new" in f[:2]:
        return None
    return "%s*%s:%s" % (g, f[0], re.sub(r"[A-Za-z]+$", "", f[1]))


def null_like(c):
    return bool(c) and bool(re.search(r"\d[NQSL]$", c.split("*")[-1] if "*" in c else c))


def heterozygosity(alleles):
    n = len(alleles)
    if n == 0:
        return float("nan")
    return 1.0 - sum((v / n) ** 2 for v in Counter(alleles).values())


def rarefy(persons, calls, gene, n, reps, rnd):
    """persons: list of ids with both haps called for `gene`; returns mean distinct alleles, mean het, mean top freq over reps."""
    d, h, t = [], [], []
    for _ in range(reps):
        sub = rnd.sample(persons, n)
        al = [a for p in sub for a in calls[gene][p]]
        c = Counter(al)
        d.append(len(c))
        h.append(heterozygosity(al))
        t.append(max(c.values()) / len(al))
    m = lambda x: sum(x) / len(x)
    return m(d), m(h), m(t)


def run(a):
    r91 = importlib.util.spec_from_file_location("r91", os.path.join(HERE, "91_novel_recurrence.py"))
    m91 = importlib.util.module_from_spec(r91)
    r91.loader.exec_module(m91)
    out = os.path.expanduser(a.out_dir)
    os.makedirs(out, exist_ok=True)
    anc, lr = m91.read_membership(os.path.expanduser(a.membership))
    first, anyk = m91.read_relatives(os.path.expanduser(a.relatedness), lr)
    unrel = set(m91.greedy_unrelated(lr, anyk))
    calls = defaultdict(lambda: defaultdict(list))   # gene -> person -> [2-field alleles]
    news = defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))   # gene -> person -> [n haps, n new protein-level, n null-like]
    with open(os.path.expanduser(a.rich)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            g = r["gene"][4:] if r["gene"].startswith("HLA-") else r["gene"]
            p = r["person_id"]
            if g not in GENES or p not in unrel or str(r.get("copy_index", "1")) not in ("1", "1.0"):
                continue
            c = r["consensus"]
            nw = news[g][p]
            nw[0] += 1
            tf = two_field(c)
            nw[1] += 1 if (c and "new" in c.split("*")[-1].split(":")[:2]) else 0
            nw[2] += 1 if null_like(c) else 0
            if tf:
                calls[g][p].append(tf)
    sig = defaultdict(set)                           # (gene, protein) -> persons (recurrent unnamed signatures, from 91)
    raw = list(csv.DictReader(open(os.path.join(os.path.expanduser(a.tables), "recurrence_raw.private.tsv")), delimiter="\t"))
    car = defaultdict(set)
    for r in csv.DictReader(open(os.path.join(os.path.expanduser(a.tables), "carriers.private.tsv")), delimiter="\t"):
        car[(r["gene"], r["ref_allele"], r["cds_events"])].add(r["person_id"])
    for r in raw:
        if r["gene"] in ("E", "F", "G") and r["ipd_match"] == "none" and int(r["persons"]) >= 2:
            sig[(r["gene"], r["protein"], r["ref_allele"])] |= car[(r["gene"], r["ref_allele"], r["cds_events"])]
    rnd = random.Random(3)
    rows = []
    byanc = defaultdict(list)
    for p in unrel:
        byanc[anc.get(p, "NA")].append(p)
    for an in sorted(byanc):
        if len(byanc[an]) < a.n:
            continue
        for g in GENES:
            ok = [p for p in byanc[an] if len(calls[g].get(p, [])) == 2]
            if len(ok) < a.n:
                continue
            d, h, t = rarefy(ok, calls, g, a.n, a.reps, rnd)
            hap = sum(news[g][p][0] for p in byanc[an] if p in news[g])
            newh = sum(news[g][p][1] for p in byanc[an] if p in news[g])
            nullh = sum(news[g][p][2] for p in byanc[an] if p in news[g])
            rows.append({"gene": g, "ancestry": an, "unrelated_persons": len(byanc[an]), "rarefied_n": a.n, "distinct_2field": d, "exp_het_2field": h, "top_allele_freq": t,
                         "haplotypes": hap, "new_protein_haps": newh, "null_like_haps": nullh,
                         "new_share": newh / max(1, hap)})
    sigrows = []
    for (g, prot, ref), ps in sorted(sig.items()):
        for an in sorted(byanc):
            c = len([p for p in ps if p in unrel and anc.get(p) == an])
            if c:
                sigrows.append((g, prot, an, c, len(byanc[an])))
    for suffix, pub in (("pub", True), ("raw.private", False)):
        pth = os.path.join(out, "landscape_%s.tsv" % suffix)
        with open(pth, "w") as f:
            f.write("gene\tancestry\tunrelated_persons\trarefied_n\tdistinct_2field\texp_het_2field\ttop_allele_freq\thaplotypes\tnew_protein_haps\tnull_like_haps\tnew_share\n")
            for r in rows:
                fm = (lambda v: supp(int(v))) if pub else (lambda v: str(int(v)))
                f.write("%s\t%s\t%s\t%d\t%.1f\t%.4f\t%.4f\t%s\t%s\t%s\t%.4f\n" % (r["gene"], r["ancestry"], fm(r["unrelated_persons"]), r["rarefied_n"], r["distinct_2field"],
                                                                                 r["exp_het_2field"], r["top_allele_freq"], fm(r["haplotypes"]), fm(r["new_protein_haps"]), fm(r["null_like_haps"]), r["new_share"]))
        with open(os.path.join(out, "signature_frequency_%s.tsv" % suffix), "w") as f:
            f.write("gene\tprotein\tancestry\tcarriers_unrelated\tunrelated_in_ancestry\n")
            for g, prot, an, c, n in sigrows:
                f.write("%s\t%s\t%s\t%s\t%s\n" % (g, prot, an, supp(c) if pub else c, supp(n) if pub else n))
        if not pub:
            for fn in ("landscape_raw.private.tsv", "signature_frequency_raw.private.tsv"):
                os.chmod(os.path.join(out, fn), 0o600)
    print("wrote landscape tables to", out, file=sys.stderr)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    base = "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux"
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--relatedness", default=base + "/relatedness/samples_relatedness.tsv")
    p.add_argument("--tables", default="~/s08/r366")
    p.add_argument("--out-dir", default="~/s09/landscape")
    p.add_argument("--n", type=int, default=400, help="rarefied number of unrelated persons per ancestry")
    p.add_argument("--reps", type=int, default=30)
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
