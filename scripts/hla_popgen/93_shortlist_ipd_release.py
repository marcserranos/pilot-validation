"""93 re-test the S08 shortlist against a CURRENT IPD-IMGT/HLA release (local; reads the gitignored case_browser/ files).

Immuannot called 'new' against IPD 3.55. For every shortlist case the observed CDS (from the stored Immuannot alignment) is compared with every
CDS of the same gene in `hla_nuc.fasta` of the release you pass: exact match -> already a named allele; observed CDS contains a shorter
(partial) entry -> 'partial'; otherwise 'still_novel'. Output: shortlist_ipd.tsv (adds ipd_status, ipd_alleles) + a printed count table.
  python3 scripts/hla_popgen/93_shortlist_ipd_release.py --nuc hla_nuc.fasta --cases case_browser/cases.json \\
      --shortlist case_browser/shortlist.tsv --out case_browser/shortlist_ipd.tsv
"""
import argparse
import csv
import json
from collections import Counter, defaultdict


def read_nuc(path):
    nuc, name, buf = defaultdict(dict), None, []

    def flush():
        if name:
            g, al = name.split("*", 1)
            nuc[g][al] = "".join(buf).upper()
    with open(path) as f:
        for ln in f:
            if ln.startswith(">"):
                flush()
                name, buf = ln.split()[1], []
            else:
                buf.append(ln.strip())
    flush()
    return nuc


def status(nuc_gene, obs):
    hit = [al for al, s in nuc_gene.items() if s == obs]
    if hit:
        return "exact", hit[:3]
    part = [al for al, s in nuc_gene.items() if len(s) < len(obs) and s in obs]
    return ("partial", part[:3]) if part else ("still_novel", [])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--nuc", required=True)
    ap.add_argument("--cases", default="case_browser/cases.json")
    ap.add_argument("--shortlist", default="case_browser/shortlist.tsv")
    ap.add_argument("--out", default="case_browser/shortlist_ipd.tsv")
    a = ap.parse_args(argv)
    nuc = read_nuc(a.nuc)
    by = {c["id"]: c for c in json.load(open(a.cases))["cases"]}
    rows = list(csv.DictReader(open(a.shortlist), delimiter="\t"))
    tot = Counter()
    for r in rows:
        obs = by[r["case"]]["alignments"][0]["aobs"].replace("-", "").upper()
        r["ipd_status"], al = status(nuc.get(r["gene"], {}), obs)
        r["ipd_alleles"] = ";".join(al)
        tot[r["ipd_status"]] += 1
    with open(a.out, "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(dict(tot))
    return tot


if __name__ == "__main__":
    main()
