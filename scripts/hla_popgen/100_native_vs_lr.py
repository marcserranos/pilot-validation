"""100 native short-read HLA calls vs long-read truth (S09-01; VM side, aggregate output).

AoU-native (HLA-HD + Polysolver + OptiType ensemble, hla_genotypes.tsv, 535k people) against full-length long-read Immuannot calls (hla_calls_rich.tsv,
~12k people), per gene x ancestry at 2-field. Truth tiers use the long-read per-call template_distance (0 = assembled gene identical to a known allele):
  td0      both haplotypes distance 0 and no warning        (ground-truth-like)
  td1-2    max distance 1-2     td3+  max distance >= 3     (divergence dose-response: does native accuracy fall as a person's haplotype diverges?)
New relative to the August baseline (analyze_confidence_vs_aou_native.py): by ancestry, by divergence, and the ERROR STRUCTURE — what the native caller
reports instead (allele confusions), hidden homozygosity (LR het called hom, LR hom called het), and uncallable rates.
Output (--out-dir, default ~/s09): native_vs_lr_pub.tsv (n and counts < 20 -> '<20'), native_vs_lr_raw.private.tsv (exact, mode 600).
"""
import argparse
import csv
import os
import sys
from collections import Counter, defaultdict

GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]


def supp(n, floor=20):
    return "<%d" % floor if 0 < n < floor else str(n)


def two(a):
    """'A*01:01:01' / 'HLA-A*01:01:01:01' / '01:01' -> '01:01'; uncallable -> None."""
    if a is None:
        return None
    a = str(a).strip()
    if a == "" or a.lower() in {"na", "nan", "none", "-", "."} or "new" in a or "undetermined" in a:
        return None
    if "*" in a:
        a = a.split("*", 1)[1]
    f = a.split(":")
    return "%s:%s" % (f[0], f[1]) if len(f) >= 2 and f[0] and f[1] else None


def compare(native, lr):
    """native, lr: 2-tuples of 2-field alleles (None allowed). -> dict of outcome flags."""
    n = [x for x in native if x]
    out = {"callable": len(n) == 2, "geno_ok": False, "alleles_ok": 0, "lr_hom": lr[0] == lr[1], "nat_hom": None, "kind": "uncallable"}
    if len(n) < 2:
        return out
    out["nat_hom"] = n[0] == n[1]
    out["geno_ok"] = sorted(n) == sorted(lr)
    c = Counter(n) & Counter(lr)
    out["alleles_ok"] = sum(c.values())
    if out["geno_ok"]:
        out["kind"] = "concordant"
    elif out["lr_hom"] and not out["nat_hom"]:
        out["kind"] = "lr_hom_native_het"
    elif (not out["lr_hom"]) and out["nat_hom"]:
        out["kind"] = "lr_het_native_hom"
    else:
        out["kind"] = "allele_swap"
    out["conf"] = [(a, b) for a, b in zip(sorted((Counter(n) - Counter(lr)).elements()), sorted((Counter(lr) - Counter(n)).elements()))]
    return out


def td_stratum(tds, warn):
    if warn:
        return "warned"
    m = max(tds)
    return "td0" if m == 0 else ("td1-2" if m <= 2 else "td3+")


def read_membership(path):
    anc = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("in_lr", "").lower() in ("true", "1"):
                anc[r["person_id"]] = (r.get("ancestry_pred") or "NA").upper()
    return anc


def run(a):
    out = os.path.expanduser(a.out_dir)
    os.makedirs(out, exist_ok=True)
    anc = read_membership(os.path.expanduser(a.membership))
    lr = defaultdict(lambda: defaultdict(dict))    # person -> gene -> hap -> (allele, td, warn)
    with open(os.path.expanduser(a.rich)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            g = r["gene"][4:] if r["gene"].startswith("HLA-") else r["gene"]
            if g not in GENES or r["person_id"] not in anc or str(r.get("copy_index", "1")) not in ("1", "1.0"):
                continue
            try:
                td = int(float(r["template_distance"]))
            except (ValueError, KeyError):
                td = 99
            w = r.get("template_warning", "") not in ("", "NA", "nan")
            lr[r["person_id"]][g][r["hap"]] = (two(r["consensus"]), td, w)
    rows = defaultdict(Counter)       # (gene, ancestry, stratum) -> counters
    conf = defaultdict(Counter)       # gene -> (native, lr) counts (td0 only)
    n_native = 0
    with open(os.path.expanduser(a.native)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            p = r["research_id"]
            if p not in lr:
                continue
            n_native += 1
            for g in GENES:
                d = lr[p].get(g, {})
                if len(d) < 2 or any(v[0] is None for v in d.values()):
                    continue
                haps = [d[h] for h in sorted(d)]
                st = td_stratum([h[1] for h in haps], any(h[2] for h in haps))
                res = compare((two(r.get(g + "_1")), two(r.get(g + "_2"))), (haps[0][0], haps[1][0]))
                for key in ((g, anc[p], st), (g, "ALL", st), (g, anc[p], "ALLTIERS"), (g, "ALL", "ALLTIERS")):
                    c = rows[key]
                    c["n"] += 1
                    c[res["kind"]] += 1
                    c["alleles_ok"] += res["alleles_ok"]
                if st == "td0":
                    for x in res.get("conf", []) if res["kind"] == "allele_swap" else []:
                        conf[g][x] += 1
    cols = ["n", "concordant", "allele_swap", "lr_hom_native_het", "lr_het_native_hom", "uncallable", "alleles_ok"]
    for suffix, fm in (("pub", supp), ("raw.private", str)):
        p = os.path.join(out, "native_vs_lr_%s.tsv" % suffix)
        with open(p, "w") as f:
            f.write("gene\tancestry\ttier\t" + "\t".join(cols) + "\tgenotype_concordance\tallele_concordance\n")
            for (g, an, st), c in sorted(rows.items()):
                if c["n"] < 20 and fm is supp:
                    f.write("%s\t%s\t%s\t<20\n" % (g, an, st))
                    continue
                cc = c["n"] - c["uncallable"]
                f.write("%s\t%s\t%s\t%s\t%.4f\t%.4f\n" % (g, an, st, "\t".join(fm(c[k]) for k in cols),
                                                          c["concordant"] / max(1, cc), c["alleles_ok"] / max(1, 2 * cc)))
        if suffix != "pub":
            os.chmod(p, 0o600)
    with open(os.path.join(out, "native_confusions_pub.tsv"), "w") as f:
        f.write("gene\tnative_called\tlr_truth\tn\n")
        for g, cc in sorted(conf.items()):
            for (nat, tru), n in cc.most_common(15):
                if n >= 20:
                    f.write("%s\t%s\t%s\t%d\n" % (g, nat, tru, n))
    print("native rows joined: %d; tables in %s" % (n_native, out), file=sys.stderr)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    base = "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux"
    p.add_argument("--rich", default="~/pipeline_outputs/hla_calls_rich.tsv")
    p.add_argument("--membership", default="~/pipeline_outputs/cohort_membership.tsv")
    p.add_argument("--native", default=base + "/hla_variants/hla_genotypes.tsv")
    p.add_argument("--out-dir", default="~/s09")
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
