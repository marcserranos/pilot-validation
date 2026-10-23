"""95 context tests for recurrent novel signatures (VM side; think-tank experiments that need no new sequencing data).

Inputs are the outputs of 91 (--out-dir/recurrence_raw.private.tsv + carriers.private.tsv), hla_calls_rich.tsv, cohort_membership.tsv, the lrWGS
manifest (platform, centre), the RNA-seq manifest and the relatedness table. For each chosen signature (the two unnamed HLA-F signatures plus KNOWN
alleles as positive controls) it reports, as aggregates only:
  carriers, ancestry mix, carriers in the RNA-seq cohort, carriers with a first-degree relative in the LR cohort (and how many relatives carry it)
  platform/centre: observed carriers per platform vs expectation from the ancestry-matched cohort (a hidden batch effect shows as a skew)
  flanking haplotype (HLA-A + HLA-G two-field on the same contig): effective number of backgrounds among carrier haplotypes vs matched non-carrier
    haplotypes of the same F lineage and ancestry (single origin = few backgrounds, recurrent mutation/artefact = many)
and cohort-wide: per-person substitution-only novel burden by ancestry x template-distance tercile (reference-bias proxy: burden should NOT track
divergence from the reference within an ancestry if the calls are not mapping artefacts).
Output: <out>/ctx95_pub.txt (counts < 20 -> '<20') and <out>/ctx95_raw.private.txt (exact, mode 600, VM only).
"""
import argparse
import csv
import importlib.util
import os
import random
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
CHOSEN = [("F", "T237I"), ("F", "P205T"), ("F", "Y106C"), ("G", "synonymous", "01:01:03:01")]


def load91():
    spec = importlib.util.spec_from_file_location("r91", os.path.join(HERE, "91_novel_recurrence.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def two(consensus):
    """'HLA-A*02:01:01:01' -> 'A*02:01'; novel or empty -> None"""
    if not consensus or "new" in consensus or "undetermined" in consensus:
        return None
    c = consensus[4:] if consensus.startswith("HLA-") else consensus
    g, _, rest = c.partition("*")
    return "%s*%s" % (g, ":".join(rest.split(":")[:2]))


def eff_number(items):
    n = len(items)
    if not n:
        return 0.0
    return 1.0 / sum((v / n) ** 2 for v in Counter(items).values())


def tersile_bounds(vals):
    s = sorted(vals)
    return s[len(s) // 3], s[2 * len(s) // 3]


def build(a, fmt):
    r91 = load91()
    out = os.path.expanduser(a.out_dir)
    raw = list(csv.DictReader(open(os.path.join(out, "recurrence_raw.private.tsv")), delimiter="\t"))
    car = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(out, "carriers.private.tsv")), delimiter="\t"):
        car[(r["gene"], r["ref_allele"], r["cds_events"])][r["person_id"]] = r["haps"].split(",")
    anc, lr = r91.read_membership(os.path.expanduser(a.membership))
    td = {}
    with open(os.path.expanduser(a.membership)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                td[r["person_id"]] = float(r["mean_template_distance"])
            except (KeyError, ValueError):
                pass
    plat, centre = {}, {}
    with open(os.path.expanduser(a.manifest)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            plat.setdefault(r["research_id"], r.get("platform") or "NA")
            centre.setdefault(r["research_id"], r.get("center") or "NA")
    rna = set()
    if os.path.exists(os.path.expanduser(a.rna_manifest)):
        with open(os.path.expanduser(a.rna_manifest)) as f:
            rna = {r["research_id"] for r in csv.DictReader(f, delimiter="\t")}
    first, anyk = r91.read_relatives(os.path.expanduser(a.relatedness), lr)
    # haplotype background from the rich table
    hapinfo = defaultdict(dict)
    with open(os.path.expanduser(a.rich)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            g = r["gene"][4:] if r["gene"].startswith("HLA-") else r["gene"]
            if g in ("A", "G", "F") and r["person_id"] in lr:
                hapinfo[(r["person_id"], r["hap"], r["contig"])][g] = r["consensus"]
    lines = []
    for sel in CHOSEN:
        gene, prot = sel[0], sel[1]
        rows = [r for r in raw if r["gene"] == gene and r["protein"] == prot and (len(sel) < 3 or r["ref_allele"].endswith(sel[2]) or sel[2] in r["ref_allele"])]
        if not rows:
            continue
        r = max(rows, key=lambda x: int(x["persons"]))
        key = (r["gene"], r["ref_allele"], r["cds_events"])
        carriers = car[key]
        cs = sorted(carriers)
        lines.append("== %s %s on %s: carriers %s (known-allele control: %s)" % (gene, prot, r["ref_allele"], fmt(len(cs)), r["ipd_match"]))
        lines.append("ancestry: " + ", ".join("%s=%s" % (k, fmt(v)) for k, v in sorted(Counter(anc.get(p, "NA") for p in cs).items())))
        lines.append("in RNA-seq cohort: %s | with a first-degree relative in the LR cohort: %s | relatives who also carry it: %s" % (
            fmt(sum(1 for p in cs if p in rna)), fmt(sum(1 for p in cs if first.get(p))),
            fmt(sum(1 for p in cs if any(q in carriers for q in first.get(p, ()))))))
        # platform / centre skew vs ancestry-matched expectation
        for label, attr in (("platform", plat), ("centre", centre)):
            by_anc = defaultdict(Counter)
            for q in lr:
                by_anc[anc.get(q, "NA")][attr.get(q, "NA")] += 1
            exp = Counter()
            for q in cs:
                tot = sum(by_anc[anc.get(q, "NA")].values())
                for k, n in by_anc[anc.get(q, "NA")].items():
                    exp[k] += n / tot
            obs = Counter(attr.get(q, "NA") for q in cs)
            lines.append("%s observed/expected: " % label + ", ".join("%s=%s/%.1f" % (k, fmt(obs.get(k, 0)), exp[k]) for k in sorted(exp)))
        # flanking background: HLA-A plus the other of F/G on the same contig as the carrier gene call
        other = "G" if gene == "F" else "F"
        cb, ctl_pool = [], []
        ancs = {anc.get(q) for q in cs}
        lineage = two(r["ref_allele"] if r["ref_allele"].startswith("HLA-") else "HLA-" + r["ref_allele"])
        for (q, h, c), d in hapinfo.items():
            bg = (two(d.get("A")), two(d.get(other)))
            if None in bg or gene not in d:
                continue
            if q in carriers and h in carriers[q] and "new" in d[gene]:
                cb.append(bg)
            elif q not in carriers and anc.get(q) in ancs and two(d[gene]) == lineage:
                ctl_pool.append(bg)
        if len(cb) >= 2 and len(ctl_pool) >= len(cb):
            rnd = random.Random(7)
            obs_eff = eff_number(cb)
            draws = sorted(eff_number(rnd.sample(ctl_pool, len(cb))) for _ in range(500))
            pct = sum(1 for d_ in draws if d_ <= obs_eff) / len(draws)
            lines.append("flanking A+%s backgrounds: carrier haplotypes %s, distinct %d, effective number %.2f; matched controls (same %s lineage+ancestry) median %.2f; carriers at percentile %.0f%% of control draws" % (
                other, fmt(len(cb)), len(set(cb)), obs_eff, gene, draws[len(draws) // 2], 100 * pct))
        else:
            lines.append("flanking background: too few haplotypes (%s carriers, %d controls)" % (fmt(len(cb)), len(ctl_pool)))
        lines.append("")
    # reference-bias proxy: burden vs template distance within ancestry
    burden = Counter()
    for key, byp in car.items():
        for p in byp:
            burden[p] += 1
    lines.append("== per-person substitution-only novel signatures (any) by ancestry x template-distance tercile (low/mid/high); persons, mean burden")
    for an in sorted({anc[p] for p in lr}):
        ps = [p for p in lr if anc[p] == an and p in td]
        if len(ps) < 60:
            continue
        lo, hi = tersile_bounds([td[p] for p in ps])
        cells = []
        for name, f in (("low", lambda x: x <= lo), ("mid", lambda x: lo < x <= hi), ("high", lambda x: x > hi)):
            g = [p for p in ps if f(td[p])]
            cells.append("%s n=%s mean=%.4f" % (name, fmt(len(g)), sum(burden[p] for p in g) / max(1, len(g))))
        lines.append("%s: %s" % (an, " | ".join(cells)))
    return "\n".join(lines) + "\n"


def run(a):
    out = os.path.expanduser(a.out_dir)
    r91 = load91()
    pub = build(a, r91.supp)
    rawtxt = build(a, str)
    with open(os.path.join(out, "ctx95_pub.txt"), "w") as f:
        f.write(pub)
    p = os.path.join(out, "ctx95_raw.private.txt")
    with open(p, "w") as f:
        f.write(rawtxt)
    os.chmod(p, 0o600)
    print("wrote ctx95 files to", out, file=sys.stderr)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    base = "~/workspace/vwb-aou-datasets-controlled-v9/v9"
    p.add_argument("--out-dir", default="~/s08/r366", help="where 91 wrote its tables")
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--manifest", default=base + "/wgs/long_read/manifest.tsv")
    p.add_argument("--rna-manifest", default=base + "/multiomics/rnaseq/manifest.tsv")
    p.add_argument("--relatedness", default=base + "/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv")
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
