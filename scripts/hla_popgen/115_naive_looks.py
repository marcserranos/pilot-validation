"""115 Naive looks (VM side): cheap aggregate histograms of the long-read HLA/KIR calls, to see whether anything stands out.

1 template distance (how far a person's HLA haplotypes are from the nearest reference template) by platform
2 number of HLA genes called per person by platform and trim tier
3 number of KIR genes present per person by ancestry (and activating-gene count)
4 DRB1 family vs presence of DRB3 / DRB4 / DRB5 (known linkage: DR3/11/12/13/14 -> DRB3, DR4/7/9 -> DRB4, DR15/16 -> DRB5, DR1/8/10 -> none)
(phase-block lengths are not in the call tables and are not computed)

Output: <out>/naive_looks_pub.tsv, long format  panel \t group \t bin \t count  (counts 1-19 written <20). No person-level output.
"""
import argparse
import csv
import importlib.util
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def supp(n, floor=20):
    return "<%d" % floor if 0 < n < floor else str(n)


def hist(values, edges):
    """Counter of bin label -> count; bins are [edges[i], edges[i+1]) with a final open bin."""
    c = Counter()
    for v in values:
        for i in range(len(edges) - 1, -1, -1):
            if v >= edges[i]:
                c["%g+" % edges[i] if i == len(edges) - 1 else "%g-%g" % (edges[i], edges[i + 1])] += 1
                break
    return c


DRB345 = {"DRB3": {"03", "11", "12", "13", "14"}, "DRB4": {"04", "07", "09"}, "DRB5": {"15", "16"}}


def expected_drb345(drb1_group):
    return {g for g, s in DRB345.items() if drb1_group in s}


def run(a):
    rows = []
    mem = {}
    with open(os.path.expanduser(a.membership)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("in_lr") in ("True", "true", "1"):
                mem[r["person_id"]] = r
    byp = defaultdict(list)
    for p, r in mem.items():
        try:
            byp[r["platform"]].append(float(r["mean_template_distance"]))
        except (KeyError, ValueError):
            pass
    for plat, vals in byp.items():
        for b, n in sorted(hist(vals, [0, 1, 2, 5, 10, 20]).items()):
            rows.append(("template_distance_mean", plat, b, supp(n)))
    for key in ("platform", "trim_tier"):
        g = defaultdict(list)
        for p, r in mem.items():
            try:
                g[r[key]].append(float(r["n_genes_called"]))
            except (KeyError, ValueError):
                pass
        for grp, vals in g.items():
            for b, n in sorted(hist(vals, [0, 8, 10, 12, 13, 14]).items()):
                rows.append(("genes_called_by_" + key, grp, b, supp(n)))
    if a.kir_root and os.path.isdir(os.path.expanduser(a.kir_root)):
        e102 = load("e102", "102_expression_dosage_tests.py")
        kd = e102.read_kir(os.path.expanduser(a.kir_root), list(mem))
        act = {"KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DS1"}
        cnt, cact = defaultdict(list), defaultdict(list)
        for p, d in kd.items():
            anc = mem[p].get("ancestry_pred", "NA")
            cnt[anc].append(sum(1 for g, v in d.items() if v > 0))
            cact[anc].append(sum(1 for g in act if d.get(g, 0) > 0))
        for anc in cnt:
            for b, n in sorted(hist(cnt[anc], [0, 8, 10, 12, 14, 16]).items()):
                rows.append(("kir_genes_present", anc, b, supp(n)))
            for b, n in sorted(hist(cact[anc], [0, 1, 2, 3, 4, 5]).items()):
                rows.append(("kir_activating_present", anc, b, supp(n)))
    calls = defaultdict(lambda: defaultdict(list))
    with open(os.path.expanduser(a.rich)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            p = r["person_id"]
            if p not in mem:
                continue
            g = r["gene"][4:] if r["gene"].startswith("HLA-") else r["gene"]
            if g in ("DRB1", "DRB3", "DRB4", "DRB5"):
                c = r.get("consensus", "")
                if c and "*" in c:
                    calls[p][g].append(c.split("*")[1].split(":")[0])
    tab = Counter()
    for p, d in calls.items():
        if len(d["DRB1"]) != 2:
            continue
        exp = set()
        for grp in d["DRB1"]:
            exp |= expected_drb345(grp)
        have = {g for g in ("DRB3", "DRB4", "DRB5") if d[g]}
        tab[("expected=" + (",".join(sorted(exp)) or "none"), "observed=" + (",".join(sorted(have)) or "none"))] += 1
    for (e, o), n in sorted(tab.items()):
        rows.append(("drb1_vs_drb345", e, o, supp(n)))
    os.makedirs(os.path.expanduser(a.out_dir), exist_ok=True)
    pth = os.path.join(os.path.expanduser(a.out_dir), "naive_looks_pub.tsv")
    with open(pth, "w") as f:
        f.write("panel\tgroup\tbin\tcount\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print("wrote", pth, len(rows), "rows", file=sys.stderr)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--kir-root", default="~/pipeline_outputs_kir")
    p.add_argument("--out-dir", default="~/s09/naive")
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
