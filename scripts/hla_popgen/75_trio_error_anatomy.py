#!/usr/bin/env python3
"""S06 75: what is a trio "typing error", exactly? Anatomy of the Mendelian mismatches and four alternative definitions (VM).

71's label: a child call is an error when the transmitting parent carries no call with the same TWO-FIELD name (A*01:01).
Immuannot names a call by the closest IPD allele; when the sequence matches none, the name carries a `new` field at the depth where
it differs (63: depth 1 undetermined, 2 protein-altering incl. frameshifts, 3 synonymous, 4 non-coding). 56's `tf` is the 2-field name,
or, when `new` sits in the first two fields (depth <= 2), the raw string (`HLA-A*01:new`). So:
  - novel only at depth 3-4: the 2-field name equals the named allele's -> NOT an error (the caller already collapsed it);
  - novel at depth <= 2 (protein-altering, incl. a homopolymer indel that shifts the frame): no named allele, the raw string
    cannot equal any parent name -> counted as an ERROR by construction.
This script separates those cases and tests definitions that treat them differently.

Per trio record (71's `records_trio`, same set, same child haplotype and transmitting parent) the child's and the parent's calls
are re-read from the call table and Table 1 (`consensus`, `cds_mut`). `cds_mut` starts with the closest reference allele
(`ALLELE|cs|aa`), which is what Immuannot would have named the call had the indel not been there.

Definitions (n errors reported for each):
  P0  current label (2-field name, `new` at depth <= 2 is its own string)
  P1  homopolymer-indel-only novel calls (child and parent) replaced by their closest reference allele's 2-field name
  P2  ALL protein-novel calls (child and parent) replaced by the closest reference allele's 2-field name
  P3  records where neither the child's call nor any parent call at the gene is protein-novel (named alleles only)
Anatomy of the P0 errors: both calls named and different; child protein-novel (by mutation class: homopolymer indel only, other
indel, substitution, mixed, no cds_mut); child named but the parent's call is protein-novel.

Output (aggregate, ~/s06/results/75): definitions.tsv, anatomy.tsv, anatomy_by_gene.tsv. Reads INTERNAL tables, writes counts only.

    cd ~/repos/pilot-validation-main && pixi run -e spechla -- python3 ~/s06/75_trio_error_anatomy.py
"""
import argparse
import importlib.util
import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _s05_common as C  # noqa: E402

I71 = os.path.expanduser("~/s06/internal/71/")
OUT = os.path.expanduser("~/s06/results/75/")
CS_TOKEN_RE = re.compile(r"(:\d+|\*[a-z]{2}|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})", re.I)


def read_int(path, **kw):
    with open(path) as f:
        first = f.readline()
    return pd.read_csv(path, sep="\t", skiprows=1 if first.startswith("#") else 0, **kw)


def is_homopolymer_only(cds_mut):
    """Same rule as 56/50/03: every difference is an indel of a single repeated base, no substitution."""
    if not isinstance(cds_mut, str) or not cds_mut:
        return False
    any_indel = False
    for entry in cds_mut.split(","):
        parts = entry.split("|")
        if len(parts) < 2:
            continue
        for tok in CS_TOKEN_RE.findall(parts[1]):
            if tok.startswith("*"):
                return False
            if tok[0] in "+-":
                any_indel = True
                if len(set(tok[1:].upper())) != 1:
                    return False
    return any_indel


def mut_class(cds_mut):
    """homopolymer_indel / other_indel / substitution / mixed / none, from the cs string of the first candidate allele."""
    if not isinstance(cds_mut, str) or not cds_mut:
        return "none"
    toks = []
    for entry in cds_mut.split(","):
        parts = entry.split("|")
        if len(parts) >= 2:
            toks += CS_TOKEN_RE.findall(parts[1])
    sub = any(t.startswith("*") for t in toks)
    ind = any(t[0] in "+-" for t in toks)
    if sub and ind:
        return "mixed"
    if sub:
        return "substitution"
    if ind:
        return "homopolymer_indel" if is_homopolymer_only(cds_mut) else "other_indel"
    return "none"


def closest_two_field(cds_mut):
    """Two-field name of the closest reference allele (first entry of cds_mut), or None."""
    if not isinstance(cds_mut, str) or not cds_mut:
        return None
    return C.two_field(cds_mut.split(",")[0].split("|")[0])


def protein_novel(consensus):
    """True when Immuannot could not give a named 2-field allele: `new` within the first two fields, or undetermined."""
    return C.two_field(consensus) is None


def call_frame(ct, t1):
    """One row per (person, hap, gene): consensus, cds_mut, name policies."""
    k = ["person_id", "hap", "gene"]
    t1 = t1[t1["copy_index"] == 1].drop_duplicates(k)[k + ["cds_mut"]]
    ct = ct.drop_duplicates(k)
    d = ct.merge(t1, on=k, how="left")
    d["hap"] = d["hap"].astype(str)
    d["named"] = [C.two_field(c) for c in d["consensus"]]
    d["pnovel"] = d["named"].isna()
    d["mclass"] = [mut_class(m) if pn else "" for m, pn in zip(d["cds_mut"], d["pnovel"])]
    d["closest"] = [closest_two_field(m) for m in d["cds_mut"]]
    raw = d["consensus"].astype(str)
    d["n0"] = d["named"].where(~d["pnovel"], raw)                       # P0: raw string for protein-novel
    homo = (d["mclass"] == "homopolymer_indel") & d["closest"].notna()
    d["n1"] = d["n0"].where(~homo, d["closest"])                          # P1
    allp = d["pnovel"] & d["closest"].notna()
    d["n2"] = d["n0"].where(~allp, d["closest"])                          # P2
    return d


def evaluate(records, calls):
    """records: person (child), other (parent), gene, hap, bad. -> (definitions DataFrame, anatomy DataFrame, by-gene DataFrame)."""
    key = {(p, h, g): i for i, (p, h, g) in enumerate(zip(calls["person_id"], calls["hap"], calls["gene"]))}
    pset = {}
    for (p, g), d in calls.groupby(["person_id", "gene"]):
        pset[(p, g)] = dict(n0=set(d["n0"]), n1=set(d["n1"]), n2=set(d["n2"]), pnovel=bool(d["pnovel"].any()))
    rows = []
    for c, par, g, h, bad in zip(records["person"], records["other"], records["gene"], records["hap"], records["bad"]):
        i = key.get((c, str(h), g))
        ps = pset.get((par, g))
        if i is None or ps is None:
            continue
        r = calls.iloc[i]
        rows.append(dict(gene=g, bad0_file=int(bad), pn_child=bool(r["pnovel"]), pn_parent=ps["pnovel"], mclass=r["mclass"],
                         e0=int(r["n0"] not in ps["n0"]), e1=int(r["n1"] not in ps["n1"]), e2=int(r["n2"] not in ps["n2"])))
    R = pd.DataFrame(rows)
    named_only = R[~R["pn_child"] & ~R["pn_parent"]]
    defs = pd.DataFrame([
        dict(definition="P0 current label (2-field name; protein-novel calls keep their `new` string)", n_records=len(R), n_errors=int(R["e0"].sum())),
        dict(definition="P1 homopolymer-indel-only novel calls collapsed to their closest reference allele", n_records=len(R), n_errors=int(R["e1"].sum())),
        dict(definition="P2 all protein-novel calls collapsed to their closest reference allele", n_records=len(R), n_errors=int(R["e2"].sum())),
        dict(definition="P3 named alleles only (no protein-novel call by the child or at the parent's gene)", n_records=len(named_only),
             n_errors=int(named_only["e0"].sum())),
    ])
    defs["error_pct"] = 100 * defs["n_errors"] / defs["n_records"]
    err = R[R["e0"] == 1]

    def cat(r):
        if r["pn_child"]:
            return "child protein-novel: " + (r["mclass"] or "none")
        if r["pn_parent"]:
            return "child named, parent has a protein-novel call"
        return "both named, different 2-field allele"
    err = err.assign(category=[cat(r) for _, r in err.iterrows()])
    an = err.groupby("category").size().rename("n_errors").reset_index()
    an["share_of_errors_pct"] = 100 * an["n_errors"] / an["n_errors"].sum()
    denom = R.assign(category=["child protein-novel" if r["pn_child"] else ("parent protein-novel" if r["pn_parent"] else "both named")
                               for _, r in R.iterrows()])
    dn = denom.groupby("category").agg(n_records=("e0", "size"), n_errors=("e0", "sum")).reset_index()
    dn["error_pct"] = 100 * dn["n_errors"] / dn["n_records"]
    bg = R.groupby("gene").agg(n_records=("e0", "size"), n_errors_P0=("e0", "sum"), n_errors_P1=("e1", "sum"), n_errors_P2=("e2", "sum")).reset_index()
    return defs, an, dn, bg, int((R["bad0_file"] != R["e0"]).sum())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--internal-dir", default=I71)
    ap.add_argument("--out-dir", default=OUT)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    ct = read_int(os.path.join(a.internal_dir, "call_table.INTERNAL.tsv"), dtype={"person_id": str, "hap": str},
                  usecols=["person_id", "hap", "gene", "consensus", "tf"])
    rec = read_int(os.path.join(a.internal_dir, "records_trio.INTERNAL.tsv"), dtype={"person": str, "other": str, "hap": str},
                   usecols=["person", "other", "gene", "hap", "bad"])
    t1 = C.load_table1(C.DEFAULT_TABLE1)
    calls = call_frame(ct, t1)
    defs, an, dn, bg, disagree = evaluate(rec, calls)
    for name, df in (("definitions.tsv", defs), ("anatomy.tsv", an), ("denominators.tsv", dn), ("anatomy_by_gene.tsv", bg)):
        df.to_csv(os.path.join(a.out_dir, name), sep="\t", index=False)
    print("records whose P0 label differs from 71's file:", disagree)
    for df in (defs, an, dn):
        print(df.round(2).to_string(index=False))
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
