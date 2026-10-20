#!/usr/bin/env python3
"""S05 57b: per-gene parental origin in trio children -> direct truth for physical phasing and typing
(VM; run from ~/s05 next to 57, 56 and _s05_common; after 57 so the same trios are reused).

For every child HLA call (copy 1) take the gene +-10 kb on its hap contig, count paternal-only vs
maternal-only k-mers (yak trioeval S lines, parental MHC databases from 57's DB window) and assign
origin = pat / mat when one side has >= 5 k-mers and >= 3x the other, else ambiguous.
  (1) gene_origin_summary.tsv   per gene: resolved %, hap1/hap2 copies from different parents %
  (2) pair_phase_truth.tsv      per gene pair on ONE contig: same-parent % (phasing truth for 56's A~B 54%)
  (3) mendel_typing.tsv         per gene x min score: assigned parent carries the child's allele (2-field / full)
Aggregate output only: ~/s05/results/57b/*.tsv; per-call INTERNAL: ~/s05/internal/57b/ (checkpoint/resume).
"""
import argparse
import importlib.util
import os
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

import _s05_common as C


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


T = load("t57", "57_yak_trioeval_pilot.py")
B = load("b56", "56_phasing_typing_benchmark.py")
OUT = os.path.expanduser("~/s05/results/57b/")
INT = os.path.expanduser("~/s05/internal/57b/")
WORK = os.path.expanduser("~/s05/work57b/")
FLANK, MIN_K, RATIO = 10_000, 5, 3
GENES = B.G8 + ["DRB3", "DRB4", "DRB5", "DQA2", "MICA", "TAP1", "TAP2"]
CK_FIELDS = ["child", "plat", "hap", "gene", "contig", "pat_k", "mat_k", "origin"]


def origin(pk, mk):
    if pk >= MIN_K and pk >= RATIO * mk:
        return "pat"
    if mk >= MIN_K and mk >= RATIO * pk:
        return "mat"
    return "amb"


def run_child(trio, bam, asm, calls):
    c, p, m = trio
    wd = WORK + c + "/"
    os.makedirs(wd, exist_ok=True)
    try:
        for who, tag in ((p, "pat"), (m, "mat")):
            T.yak_db(bam[who]["bam"], T.DB["MHC"], wd + f"{tag}.yak", 2)
        rows = []
        for h in ("1", "2"):
            sub = calls[calls["hap"].str.endswith(h)]
            if sub.empty:
                continue
            fa = wd + f"h{h}.fa.gz"
            T.sh(f"cp '{asm[c]['h' + h]}' {fa}.part && mv {fa}.part {fa} && {T.SAM} faidx {fa}")
            ln = dict(l.split("\t")[:2] for l in open(fa + ".fai"))
            reg, names = [], []
            for g, ctg_r, s, e in zip(sub["gene"], sub["contig"], sub["gene_start"], sub["gene_end"]):
                # call-table contigs are extraction ranges "ctg:S-E"; gene coords are relative to S
                ctg, _, rng = str(ctg_r).partition(":")
                off = int(rng.split("-")[0]) - 1 if rng else 0
                if ctg not in ln or pd.isna(s) or pd.isna(e):
                    continue
                s, e = off + int(min(s, e)), off + int(max(s, e))
                a, b = max(1, s - FLANK), min(int(ln[ctg]), e + FLANK)
                reg.append(f"{ctg}:{a}-{b}")
                names.append((g, ctg))
            if not reg:
                continue
            open(wd + "w.reg", "w").write("\n".join(reg) + "\n")
            T.sh(f"{T.SAM} faidx {fa} -r {wd}w.reg > {wd}w.fa")
            te = T.sh(f"{T.YAK} trioeval -t2 {wd}pat.yak {wd}mat.yak {wd}w.fa 2>/dev/null")
            s_lines = {x[1]: (int(x[2]), int(x[3])) for x in (l.split("\t") for l in te.splitlines())
                       if x[0] == "S" and len(x) >= 4}
            for (g, ctg), r in zip(names, reg):
                pk, mk = s_lines.get(r, (0, 0))
                rows.append([c, asm[c]["plat"], "hap" + h, g, ctg, pk, mk, origin(pk, mk)])
            T.sh(f"rm -f {fa} {fa}.fai {fa}.gzi")
        return rows, None
    except Exception as ex:
        return [], f"{c}: {ex}"
    finally:
        T.sh(f"rm -rf {wd}")


def summaries(o, calls_all, trio_of):
    """o: per-call origins (INTERNAL); calls_all: prepared calls (tf, consensus, score) for trio people."""
    key = calls_all[calls_all["copy_index"] == 1].set_index(["person_id", "hap", "gene"])
    o = o.join(key[["tf", "consensus", "score"]], on=["child", "hap", "gene"], how="left")
    # (1) per gene: resolved, complementary hap1/hap2 origins
    g1 = []
    for g, d in o.groupby("gene"):
        res = d[d["origin"] != "amb"]
        piv = res.pivot_table(index="child", columns="hap", values="origin", aggfunc="first")
        both = piv.dropna() if {"hap1", "hap2"} <= set(piv.columns) else piv.iloc[0:0]
        # report the failure (both copies from one parent) so a small count is masked, not implied by its complement
        same = int((both.get("hap1") == both.get("hap2")).sum()) if len(both) else 0
        r1, r2 = C.rate_cells(len(d) - len(res), len(d)), C.rate_cells(same, len(both))
        g1.append({"gene": g, "n_calls": r1["d"], "n_ambiguous": r1["n"], "ambiguous_pct": r1["pct"],
                   "n_children_both_haps_resolved": r2["d"], "n_both_copies_same_parent": r2["n"],
                   "both_copies_same_parent_pct": r2["pct"], "lo": r2["lo"], "hi": r2["hi"]})
    C.write_tsv(pd.DataFrame(g1), OUT + "gene_origin_summary.tsv")
    # (2) gene pairs on one contig: same parent?
    rows = []
    idx = {(c, h, g): (ctg, og) for c, h, g, ctg, og in zip(o["child"], o["hap"], o["gene"], o["contig"], o["origin"])}
    recs = []
    for (c, h, g), (ctg, og) in idx.items():
        for ga, gb in B.HAPPAIRS + [("A", "C"), ("DRB1", "DQA1"), ("DQB1", "DPB1")]:
            if g != ga or (c, h, gb) not in idx:
                continue
            ctg2, og2 = idx[(c, h, gb)]
            if ctg2 == ctg and og != "amb" and og2 != "amb":
                recs.append((ga + "~" + gb, o.loc[o["child"] == c, "plat"].iloc[0], int(og != og2)))
    r = pd.DataFrame(recs, columns=["pair", "plat", "bad"])
    for (pair, plat), d in [((p, "all"), d) for p, d in r.groupby("pair")] + list(r.groupby(["pair", "plat"])):
        q = C.rate_cells(int(d["bad"].sum()), len(d))
        rows.append({"pair": pair, "child_platform": plat, "n_same_contig_haps": q["d"],
                     "n_parent_mismatch": q["n"], "mismatch_pct": q["pct"], "lo": q["lo"], "hi": q["hi"]})
    pt = pd.DataFrame(rows)
    # complementary suppression: if one platform's mismatch count is masked, blank the other platform's count and
    # rate too, otherwise "all" minus the shown platform recovers the masked count
    for pair, g in pt[pt["child_platform"] != "all"].groupby("pair"):
        if (g["n_parent_mismatch"] == "<20").any():
            ix = g.index[g["n_parent_mismatch"] != "<20"]
            pt.loc[ix, ["n_parent_mismatch", "mismatch_pct", "lo", "hi"]] = "<20*"
    C.write_tsv(pt, OUT + "pair_phase_truth.tsv")
    # (3) Mendelian typing: does the assigned parent carry the child's allele?
    par = {}
    for (pid, g), d in calls_all.groupby(["person_id", "gene"]):
        par[(pid, g)] = (set(d["tf"]), set(d["consensus"]))
    recs = []
    for c, g, og, tf, cons, sc in zip(o["child"], o["gene"], o["origin"], o["tf"], o["consensus"], o["score"]):
        if og == "amb" or pd.isna(tf):
            continue
        pid = trio_of[c][0] if og == "pat" else trio_of[c][1]
        if (pid, g) not in par:
            continue
        ptf, pcons = par[(pid, g)]
        recs.append((g, int(sc), int(tf not in ptf), int(cons not in pcons)))
    r = pd.DataFrame(recs, columns=["gene", "score", "bad2", "badfull"])
    rows = []
    for gene, d in [("ALL", r), ("classical8", r[r["gene"].isin(B.G8)])] + list(r.groupby("gene")):
        for t in (0, 4, 6):
            s = d[d["score"] >= t]
            q2, qf = C.rate_cells(int(s["bad2"].sum()), len(s)), C.rate_cells(int(s["badfull"].sum()), len(s))
            rows.append({"gene": gene, "min_score": t, "n_calls": q2["d"], "mismatch_2field_pct": q2["pct"],
                         "lo": q2["lo"], "hi": q2["hi"], "mismatch_full_pct": qf["pct"]})
    C.write_tsv(pd.DataFrame(rows), OUT + "mendel_typing.tsv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--summarise-only", action="store_true")
    ap.add_argument("--all-trios", action="store_true", help="every trio with assembly + HiFi parents, not only 57's 50")
    ap.add_argument("--extra-genes", help="comma list (e.g. DRA,C4A): measure ONLY these genes, into "
                    "gene_origin_extra.INTERNAL.tsv (S06 71c contig-vs-truth matrix); 57b's own tables are not rebuilt")
    a = ap.parse_args()
    global GENES
    for d in (OUT, INT, WORK):
        os.makedirs(d, exist_ok=True)
    kids = set(pd.read_csv(os.path.expanduser("~/s05/internal/57/per_trio_region.INTERNAL.tsv"), sep="\t",
                           dtype={"child": str})["child"])
    cohort = C.load_cohort(C.DEFAULT_COHORT)
    lr = set(cohort.loc[cohort["in_lr"].astype(str).str.lower() == "true", "person_id"])
    plat = dict(zip(cohort["person_id"], cohort["platform"]))
    trios, seen = [], set()
    for t in T.find_trios(lr):
        if (a.all_trios or t[0] in kids) and t[0] not in seen:
            seen.add(t[0])
            trios.append(t)
    trio_of = {t[0]: (t[1], t[2]) for t in trios}
    people = {x for t in trios for x in t}
    C.TABLE1_USECOLS = C.TABLE1_USECOLS + ["n_tied"]
    t1 = C.load_table1(C.DEFAULT_TABLE1)
    extra = [g for g in (a.extra_genes or "").split(",") if g]
    if extra:       # in-process only: 56's gene universe and published numbers are untouched
        GENES = list(extra)
        B.ALLG.extend(g for g in extra if g not in B.ALLG)
    calls = B.prepare(t1, plat)  # score q75 on the whole cohort, as in 56
    calls = calls[calls["person_id"].isin(people) & calls["gene"].isin(GENES)].copy()
    ck = INT + ("gene_origin_extra.INTERNAL.tsv" if extra else "gene_origin.INTERNAL.tsv")
    if not a.summarise_only:
        bam, asm = T.manifest_index(people)
        done = set(pd.read_csv(ck, sep="\t", dtype={"child": str})["child"]) if os.path.exists(ck) else set()
        if not os.path.exists(ck):
            open(ck, "w").write("\t".join(CK_FIELDS) + "\n")
        todo = [t for t in trios if t[0] not in done and t[0] in asm and t[1] in bam and t[2] in bam]
        C.log(f"trios {len(trios)}, to run {len(todo)}")
        c1 = calls[calls["copy_index"] == 1]
        with ThreadPoolExecutor(a.workers) as ex:
            futs = [ex.submit(run_child, t, bam, asm, c1[c1["person_id"] == t[0]]) for t in todo]
            for k, fu in enumerate(futs, 1):
                rows, err = fu.result()
                with open(ck, "a") as f:
                    f.write("".join("\t".join(map(str, r)) + "\n" for r in rows))
                    f.flush()
                    os.fsync(f.fileno())
                C.log(f"child {k}/{len(todo)} calls {len(rows)}" + (f" ERROR {err}" if err else ""))
    o = pd.read_csv(ck, sep="\t", dtype={"child": str, "contig": str})
    if not extra:
        summaries(o, calls, trio_of)
    C.log("DONE 57b")


if __name__ == "__main__":
    main()
