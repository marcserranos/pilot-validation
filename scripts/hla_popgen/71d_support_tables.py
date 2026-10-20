#!/usr/bin/env python3
"""S06 71d: two small aggregate tables for the report, from 71's INTERNAL tables (VM only).

  pair_trio_overlap.tsv   how the pair truth set and the trio truth set overlap (people, relationships, records), for one 71 run
  named_effect.tsv        what the named-alleles-only fix changed: how many calls are novel in field 1-2 and how many records / errors of the first
                          analysis involve such a call (needs the first analysis' tables: --main-internal)

    python3 71d_support_tables.py --internal ~/s06/internal/71_named/ --out ~/s06/results/71_named/71 \\
        --main-internal ~/s06/internal/71/
"""
import argparse
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _s05_common as C  # noqa: E402


def rd(d, name, cols=None):
    p = os.path.join(d, name)
    with open(p) as f:
        first = f.readline()
    return pd.read_csv(p, sep="\t", dtype=str, usecols=cols, skiprows=1 if first.startswith("#") else 0)


def key(a, b):
    return tuple(sorted((a, b)))


def overlap_table(pair, trio):
    pp = set(key(a, b) for a, b in zip(pair["person"], pair["other"]))
    tp = set(key(a, b) for a, b in zip(trio["person"], trio["other"]))
    kids, parents = set(trio["person"]), set(trio["other"])
    pair_people = set(pair["person"]) | set(pair["other"])
    trio_people = kids | parents
    on_edge = sum(key(a, b) in tp for a, b in zip(pair["person"], pair["other"]))
    rows = [("pairs_distinct", len(pp), "first-degree pairs in the pair truth set"), ("pair_people", len(pair_people), "people in those pairs"),
            ("pair_records", len(pair), "pair-genes"), ("trio_children", len(kids), "children with a resolved trio record"),
            ("trio_edges", len(tp), "distinct child-parent links in the trio set"), ("trio_people", len(trio_people), "children plus parents"),
            ("trio_records", len(trio), "child haplotype calls with a resolved parent of origin"),
            ("edges_in_both", len(pp & tp), "trio child-parent links that are also pairs"),
            ("pair_records_on_trio_edges", on_edge, "pair-genes that belong to a trio child-parent link"),
            ("people_in_both", len(pair_people & trio_people), "people in both sets")]
    return pd.DataFrame(rows, columns=["quantity", "value", "note"])


def named_effect(ct, trio, pair, named_trio, named_pair):
    ct = ct.copy()
    un = set(zip(ct.loc[ct["named"] == "False", "person_id"], ct.loc[ct["named"] == "False", "gene"]))
    unh = set(zip(ct.loc[ct["named"] == "False", "person_id"], ct.loc[ct["named"] == "False", "hap"], ct.loc[ct["named"] == "False", "gene"]))
    t = trio.copy()
    t["bad"] = t["bad"].astype(int)
    t["any_un"] = [((p, h, g) in unh) or ((o, g) in un) for p, h, g, o in zip(t["person"], t["hap"], t["gene"], t["other"])]
    p = pair.copy()
    p["bad"] = p["bad"].astype(int)
    p["any_un"] = [((a, g) in un) or ((b, g) in un) for a, b, g in zip(p["person"], p["other"], p["gene"])]
    rows = [("calls_total", len(ct)), ("calls_unnamed", int((ct["named"] == "False").sum())),
            ("trio_records", len(t)), ("trio_errors", int(t["bad"].sum())),
            ("trio_any_unnamed_records", int(t["any_un"].sum())), ("trio_any_unnamed_errors", int(t.loc[t["any_un"], "bad"].sum())),
            ("trio_named_records", len(named_trio)), ("trio_named_errors", int(named_trio["bad"].astype(int).sum())),
            ("pair_records", len(p)), ("pair_errors", int(p["bad"].sum())),
            ("pair_any_unnamed_records", int(p["any_un"].sum())), ("pair_any_unnamed_errors", int(p.loc[p["any_un"], "bad"].sum())),
            ("pair_named_records", len(named_pair)), ("pair_named_errors", int(named_pair["bad"].astype(int).sum()))]
    return pd.DataFrame(rows, columns=["quantity", "value"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--internal", required=True, help="71 internal dir of the run to describe")
    ap.add_argument("--out", required=True)
    ap.add_argument("--main-internal", help="internal dir of the first (all-alleles) run, for named_effect.tsv")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    pair = rd(a.internal, "records_pair.INTERNAL.tsv", ["person", "other", "gene", "bad"])
    trio = rd(a.internal, "records_trio.INTERNAL.tsv", ["person", "other", "gene", "hap", "bad"])
    overlap_table(pair, trio).to_csv(os.path.join(a.out, "pair_trio_overlap.tsv"), sep="\t", index=False)
    if a.main_internal:
        ct = rd(a.main_internal, "call_table.INTERNAL.tsv", ["person_id", "hap", "gene", "named"])
        mp = rd(a.main_internal, "records_pair.INTERNAL.tsv", ["person", "other", "gene", "bad"])
        mt = rd(a.main_internal, "records_trio.INTERNAL.tsv", ["person", "other", "gene", "hap", "bad"])
        named_effect(ct, mt, mp, trio, pair).to_csv(os.path.join(a.out, "named_effect.tsv"), sep="\t", index=False)
    C.log("DONE 71d")
    return 0


if __name__ == "__main__":
    sys.exit(main())
