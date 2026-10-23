#!/usr/bin/env python3
"""S06 77c: Mendelian-mismatch records for KIR (trio / pair / twin), the KIR analogue of 71's records (VM).

Allele level: label3 (protein field, 'KIR2DL1*003') is the main level, label5 ('KIR2DL1*00302') the second one; calls without a label
(novel_protein, undetermined) form a separate NOVEL stream and are never scored as errors.
KIR presence/absence is biology (gene content varies between haplotypes), so a missing gene is not an error by default. For a child call
with a resolved yak origin (77b) the transmitting parent T and the other parent O are compared at the locus (pooled over both haps,
copy 1 and further copies):
  match                  child's label is among T's labels
  only_other_parent      T has the locus but not that label, O has it   -> yak origin error (or a call inconsistency), not a typing error
  neither                T has the locus, neither parent carries the label   -> true-error candidate (the only 'bad')
  absent_in_parents      neither parent has any call at the locus           -> biology / dropout, not an error
  tx_absent_other_has    T has no call at the locus, O has                  -> origin error or dropout of T's call
Locus mode 'gene' compares per gene; 'group' merges allelic / paralogous partners (2DL2/2DL3, 3DL1/3DS1, 2DL5A/2DL5B) into one locus.
Variants are columns, not separate files: child_complete / tx_complete / other_complete (both haps of the person are 'complete' by 51's
rule), child_hap_class, strong (>=10 k-mers and >=5x), tx_unnamed (T carries an unlabelled call at the locus), tx_n_calls.
Pairs (first degree, kin 2^-2.5..2^-1.5) and twins (kin > 2^-1.5) among 77's people: per locus 'share' (>=1 shared label) / 'no_share'
(both have calls and none shared = bad) / 'presence_discordant' (one lacks the locus) / 'absent_both'; twins also record exact label-set
equality.  Aggregates (counts masked 1-19): ~/s06/results/77c/; person-level records: ~/s06/internal/77c/.

    cd ~/repos/pilot-validation-main && pixi run -e spechla -- python3 ~/repos/pv-s06/scripts/hla_popgen/77c_kir_trio_records.py
"""
import argparse
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _s05_common as C  # noqa: E402

OUT = os.path.expanduser("~/s06/results/77c/")
INT = os.path.expanduser("~/s06/internal/77c/")
I77 = os.path.expanduser("~/s06/internal/77/")
I77B = os.path.expanduser("~/s06/internal/77b/")
GROUP = {"KIR2DL2": "KIR2DL2/3", "KIR2DL3": "KIR2DL2/3", "KIR3DL1": "KIR3DL1/S1", "KIR3DS1": "KIR3DL1/S1",
         "KIR2DL5A": "KIR2DL5", "KIR2DL5B": "KIR2DL5"}
CLASSES = ["match", "only_other_parent", "neither", "absent_in_parents", "tx_absent_other_has"]


def locus_of(gene, mode):
    return GROUP.get(gene, gene) if mode == "group" else gene


def read_int(path, **kw):
    with open(path) as f:
        first = f.readline()
    return pd.read_csv(path, sep="\t", skiprows=1 if first.startswith("#") else 0, **kw)


def person_index(calls, mode, col):
    """(person, locus) -> dict(labels=set, n=calls at the locus, unnamed=any call without a label at this level)."""
    idx = {}
    for p, g, lab in zip(calls["person_id"], calls["gene"], calls[col]):
        d = idx.setdefault((p, locus_of(g, mode)), {"labels": set(), "n": 0, "unnamed": False})
        d["n"] += 1
        if isinstance(lab, str):
            d["labels"].add(lab)
        else:
            d["unnamed"] = True
    return idx


def complete_map(haps):
    """person -> True when both hap1 and hap2 are 'complete' (51's rule)."""
    out = {}
    for p, d in haps.groupby("person_id"):
        out[p] = bool(set(d["hap"]) == {"hap1", "hap2"} and (d["hap_class"] == "complete").all())
    return out


def classify(label, tx, ot):
    """label: child label; tx / ot: index entries or None. Returns one of CLASSES."""
    tl = tx["labels"] if tx else set()
    ol = ot["labels"] if ot else set()
    t_has, o_has = bool(tx and tx["n"]), bool(ot and ot["n"])
    if t_has:
        if label in tl:
            return "match"
        return "only_other_parent" if label in ol else "neither"
    if not o_has:
        return "absent_in_parents"
    return "tx_absent_other_has"


def trio_records(calls, origin, trio_of, comp, mode="gene"):
    """origin: 77b table (child, hap, gene, copy_index, pat_k, mat_k, origin, strong); trio_of: child -> (pat, mat)."""
    i3, i5 = person_index(calls, mode, "label3"), person_index(calls, mode, "label5")
    ck = calls.set_index(["person_id", "hap", "gene", "copy_index"])
    ck = ck[~ck.index.duplicated()]
    rows = []
    for c, h, g, cp, og, st, pk, mk in zip(origin["child"], origin["hap"], origin["gene"], origin["copy_index"], origin["origin"],
                                           origin["strong"], origin["pat_k"], origin["mat_k"]):
        if og not in ("pat", "mat") or c not in trio_of or (c, h, g, cp) not in ck.index:
            continue
        tx_id, ot_id = (trio_of[c][0], trio_of[c][1]) if og == "pat" else (trio_of[c][1], trio_of[c][0])
        r = ck.loc[(c, h, g, cp)]
        loc = locus_of(g, mode)
        tx3, ot3 = i3.get((tx_id, loc)), i3.get((ot_id, loc))
        tx5, ot5 = i5.get((tx_id, loc)), i5.get((ot_id, loc))
        stream = "named" if isinstance(r["label3"], str) else "novel"
        row = dict(child=c, tx=tx_id, other=ot_id, gene=g, locus=loc, mode=mode, hap=h, copy_index=cp, origin=og,
                   strong=int(st), pat_k=pk, mat_k=mk, tier=r["novelty_tier"], stream=stream,
                   label3=r["label3"], label5=r["label5"], child_hap_class=r["hap_class"],
                   child_complete=bool(comp.get(c, False)), tx_complete=bool(comp.get(tx_id, False)),
                   other_complete=bool(comp.get(ot_id, False)), tx_n_calls=tx3["n"] if tx3 else 0,
                   tx_unnamed=int(bool(tx3 and tx3["unnamed"])), template_distance=r.get("template_distance"),
                   cds_distance=r.get("cds_distance"), n_called=r.get("n_called"))
        if stream == "named":
            row["cls3"] = classify(r["label3"], tx3, ot3)
            row["cls5"] = classify(r["label5"], tx5, ot5) if isinstance(r["label5"], str) else None
            row["bad3"] = int(row["cls3"] == "neither")
            row["origin_err3"] = int(row["cls3"] in ("only_other_parent", "tx_absent_other_has"))
        else:
            row["cls3"] = row["cls5"] = None
            row["bad3"] = row["origin_err3"] = None
        rows.append(row)
    return pd.DataFrame(rows)


def pair_records(calls, pairs, comp, mode="gene", kind="pair"):
    """pairs: iterable of (i, j). One record per pair and locus present in either person."""
    i3 = person_index(calls, mode, "label3")
    loci = {}
    for (p, loc) in i3:
        loci.setdefault(p, set()).add(loc)
    rows = []
    for i, j in pairs:
        for loc in sorted(loci.get(i, set()) | loci.get(j, set())):
            a, b = i3.get((i, loc)), i3.get((j, loc))
            ha, hb = bool(a and a["n"]), bool(b and b["n"])
            if not ha and not hb:
                continue
            if ha and hb:
                la, lb = a["labels"], b["labels"]
                if not la or not lb:
                    cls = "unnamed"
                else:
                    cls = "share" if la & lb else "no_share"
                same = int(la == lb) if la and lb else None
            else:
                cls, same = "presence_discordant", None
            rows.append(dict(person=i, other=j, locus=loc, mode=mode, kind=kind, cls=cls,
                             bad=int(cls == "no_share"), same_set=same,
                             both_complete=bool(comp.get(i, False) and comp.get(j, False))))
    return pd.DataFrame(rows)


def summarize_trio(rec):
    rows = []
    variants = [("all", rec), ("complete_only", rec[rec["child_complete"] & rec["tx_complete"] & rec["other_complete"]]),
                ("strong_origin", rec[rec["strong"] == 1]),
                ("tx_named", rec[rec["tx_unnamed"] == 0])]
    for vn, d in variants:
        d = d[d["stream"] == "named"]
        for gene, g in [("ALL", d)] + list(d.groupby("gene")):
            vc = g["cls3"].value_counts()
            n_inf = int(vc.get("match", 0) + vc.get("only_other_parent", 0) + vc.get("neither", 0))
            r = C.rate_cells(int(vc.get("neither", 0)), n_inf)
            ro = C.rate_cells(int(vc.get("only_other_parent", 0)), n_inf)
            row = {"variant": vn, "gene": gene, "n_records": C.mask_n(len(g)), "n_informative": C.mask_n(n_inf),
                   "neither_n": r["n"], "neither_pct": r["pct"], "lo": r["lo"], "hi": r["hi"],
                   "only_other_pct": ro["pct"]}
            for k in CLASSES:
                row["n_" + k] = C.mask_n(int(vc.get(k, 0)))
            rows.append(row)
    return pd.DataFrame(rows)


def summarize_pair(rec):
    rows = []
    for kind, d0 in rec.groupby("kind"):
        for vn, d in (("all", d0), ("complete_only", d0[d0["both_complete"]])):
            for loc, g in [("ALL", d)] + list(d.groupby("locus")):
                inf = g[g["cls"].isin(["share", "no_share"])]
                r = C.rate_cells(int((inf["cls"] == "no_share").sum()), len(inf))
                pd_ = C.rate_cells(int((g["cls"] == "presence_discordant").sum()), len(g))
                rows.append({"kind": kind, "variant": vn, "locus": loc, "n_records": C.mask_n(len(g)),
                             "n_informative": r["d"], "no_share_n": r["n"], "no_share_pct": r["pct"], "lo": r["lo"], "hi": r["hi"],
                             "presence_discordant_pct": pd_["pct"]})
    return pd.DataFrame(rows)


def select_pairs(people, rel):
    """rel: iterable of (i, j, kin) -> (pairs, twins) lists among `people`."""
    pairs, twins = [], []
    for i, j, k in rel:
        if i not in people or j not in people or not (k == k):
            continue
        if k > 2 ** -1.5:
            twins.append((i, j))
        elif k >= 2 ** -2.5:
            pairs.append((i, j))
    return pairs, twins


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", default=I77 + "kir_calls.INTERNAL.tsv")
    ap.add_argument("--haps", default=I77 + "kir_haps.INTERNAL.tsv")
    ap.add_argument("--origin", default=I77B + "gene_origin_kir.INTERNAL.tsv")
    ap.add_argument("--trios", default=I77B + "trios.INTERNAL.tsv")
    ap.add_argument("--out-dir", default=OUT)
    ap.add_argument("--int-dir", default=INT)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    os.makedirs(a.int_dir, exist_ok=True)
    dt = {"person_id": str, "contig": str}
    calls = read_int(a.calls, dtype=dt)
    haps = read_int(a.haps, dtype=dt)
    comp = complete_map(haps)
    origin = pd.read_csv(a.origin, sep="\t", dtype={"child": str, "contig": str})
    tr = pd.read_csv(a.trios, sep="\t", dtype=str)
    trio_of = {c: (p, m) for c, p, m in zip(tr["child"], tr["pat"], tr["mat"])}
    recs = []
    for mode in ("gene", "group"):
        r = trio_records(calls, origin, trio_of, comp, mode)
        recs.append(r)
    trio = pd.concat(recs, ignore_index=True)
    trio.to_csv(os.path.join(a.int_dir, "records_trio_kir.INTERNAL.tsv"), sep="\t", index=False)
    for mode, d in trio.groupby("mode"):
        s = summarize_trio(d)
        C.write_tsv(s, os.path.join(a.out_dir, "trio_classes_%s.tsv" % mode))
    nov = trio[trio["mode"] == "gene"].groupby(["gene", "stream"]).size().rename("n").reset_index()
    nov["n"] = [C.mask_n(v) for v in nov["n"]]
    C.write_tsv(nov, os.path.join(a.out_dir, "trio_streams.tsv"))
    rel = C.load_relatedness_pairs(C.DEFAULT_RELATEDNESS)
    people = set(calls["person_id"])
    pairs, twins = select_pairs(people, rel)
    out = []
    for mode in ("gene", "group"):
        out.append(pair_records(calls, pairs, comp, mode, "pair"))
        out.append(pair_records(calls, twins, comp, mode, "twin"))
    pr = pd.concat(out, ignore_index=True)
    pr.to_csv(os.path.join(a.int_dir, "records_pair_kir.INTERNAL.tsv"), sep="\t", index=False)
    for mode, d in pr.groupby("mode"):
        C.write_tsv(summarize_pair(d), os.path.join(a.out_dir, "pair_classes_%s.tsv" % mode))
    C.log("DONE 77c: trio records %d, pair records %d" % (len(trio), len(pr)))


if __name__ == "__main__":
    main()
