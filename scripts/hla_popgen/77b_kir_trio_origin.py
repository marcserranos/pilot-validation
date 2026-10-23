#!/usr/bin/env python3
"""S06 77b: per-KIR-gene parental origin in trio children (KIR analogue of 57b; VM).

For each child KIR call (every copy_index) the gene +-10 kb is cut from the child's hap{1,2}.trimmed.fa (the Immuannot-trimmed KIR
window, chr19:54.6-54.92 Mb +100 kb pad; contig names = GTF column 1 and gene coordinates are contig-relative, so no offset) and scored
with `yak trioeval` against the parents' KIR k-mer databases (57's DB['KIR'] = chr19:54.5-55.1 Mb, built from each parent's GRCh38 BAM).
origin = pat / mat when one side has >= 5 k-mers and >= 3x the other, else amb (57b's rule); the raw counts are kept and `strong` marks
>= 10 k-mers and >= 5x. 'nowin' = the contig was not found in the trimmed fasta (loud in the log). pat/mat are just the two parent ids of
the trio tuple (sorted neighbour ids), not sexes.
Trios = 57's find_trios, the FIRST candidate per child (as 57b and 71 do), restricted to long-read children with assembly + two HiFi parents
and with KIR calls in 77's table. Each parent's DB is built ONCE (cached under work/db, used by every trio the person belongs to, deleted
when its last trio is done); per-trio work dirs are removed.
Region-level summary per trio and hap (as 57): yak trioeval W (switch) and H (hamming) lines over the whole trimmed KIR window.
Caveat (state in the report): parental k-mers whose reads map elsewhere than chr19:54.5-55.1 Mb in GRCh38 (KIR genes missing from the
primary assembly: 2DL2, 2DL5, 2DS1-3,5, 3DS1) are absent from the DBs -> more ambiguous windows for those genes.

INTERNAL (VM only): ~/s06/internal/77b/{gene_origin_kir,region_kir,trios}.INTERNAL.tsv.  Aggregates: ~/s06/results/77b/*.tsv.

    cd ~/repos/pilot-validation-main && pixi run -e spechla -- python3 ~/repos/pv-s06/scripts/hla_popgen/77b_kir_trio_origin.py --limit 3 --workers 1   # smoke
    ... --workers 3        # all trios, resumable (finished children are skipped)
    ... --summarise-only
"""
import argparse
import importlib.util
import os
import sys
import threading
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _s05_common as C  # noqa: E402

OUT = os.path.expanduser("~/s06/results/77b/")
INT = os.path.expanduser("~/s06/internal/77b/")
WORK = os.path.expanduser("~/s06/work77b/")
CALLS = os.path.expanduser("~/s06/internal/77/kir_calls.INTERNAL.tsv")
FLANK, MIN_K, RATIO = 10_000, 5, 3
STRONG_K, STRONG_RATIO = 10, 5
GENE_FIELDS = ["child", "plat", "hap", "gene", "copy_index", "contig", "win_start", "win_end", "pat_k", "mat_k", "origin", "strong"]
REGION_FIELDS = ["child", "plat", "hap", "bp", "sw_err", "sw_tot", "hm_err", "hm_tot"]
_T = []


def T():
    if not _T:
        spec = importlib.util.spec_from_file_location("t57_for_77b", os.path.join(HERE, "57_yak_trioeval_pilot.py"))
        m = importlib.util.module_from_spec(spec)
        sys.modules["t57_for_77b"] = m
        spec.loader.exec_module(m)
        _T.append(m)
    return _T[0]


def origin(pk, mk, min_k=MIN_K, ratio=RATIO):
    if pk >= min_k and pk >= ratio * mk:
        return "pat"
    if mk >= min_k and mk >= ratio * pk:
        return "mat"
    return "amb"


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for ln in f:
            ln = ln.rstrip("\r\n")
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = ln[1:].split()[0], []
            elif ln:
                buf.append(ln)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def make_windows(seqs, calls, flank=FLANK):
    """calls: iterable of (gene, copy_index, contig, gene_start, gene_end). -> (fasta text, meta list aligned to window names w0, w1, ...);
    meta rows (gene, copy, contig, a, b) with 1-based inclusive window a-b; a call whose contig is missing gets origin 'nowin'."""
    fa, meta, nowin = [], [], []
    for g, cp, ctg, s, e in calls:
        if ctg not in seqs or pd.isna(s) or pd.isna(e):
            nowin.append((g, cp, ctg))
            continue
        L = len(seqs[ctg])
        s, e = int(min(s, e)), int(max(s, e))
        a, b = max(1, s - flank), min(L, e + flank)
        name = "w%d" % len(meta)
        fa.append(">%s\n%s\n" % (name, seqs[ctg][a - 1:b]))
        meta.append((g, cp, ctg, a, b))
    return "".join(fa), meta, nowin


def strong(pk, mk):
    return int(origin(pk, mk, STRONG_K, STRONG_RATIO) != "amb")


def run_child(trio, bam, plat, calls, kir_root, dbpath, dblock):
    c, p, m = trio
    wd = WORK + c + "/"
    os.makedirs(wd, exist_ok=True)
    t = T()
    rows, reg, notes = [], [], []
    try:
        for who in (p, m):
            with dblock[who]:
                t.yak_db(bam[who]["bam"], t.DB["KIR"], dbpath(who), 2)
        pat, mat = dbpath(p), dbpath(m)
        for h in ("1", "2"):
            sub = calls[calls["hap"] == "hap" + h]
            fa = os.path.join(kir_root, c, "immuannot_output", "hap%s.trimmed.fa" % h)
            if sub.empty or not os.path.exists(fa):
                continue
            seqs = read_fasta(fa)
            txt, meta, nowin = make_windows(
                seqs, zip(sub["gene"], sub["copy_index"], sub["contig"], sub["gene_start"], sub["gene_end"]))
            for g, cp, ctg in nowin:
                rows.append([c, plat, "hap" + h, g, cp, ctg, 0, 0, 0, 0, "nowin", 0])
            if nowin:
                notes.append("hap%s: %d calls without contig in trimmed.fa" % (h, len(nowin)))
            if meta:
                open(wd + "w.fa", "w").write(txt)
                te = t.sh("%s trioeval -t2 %s %s %sw.fa 2>/dev/null" % (t.YAK, pat, mat, wd))
                sl = {x[1]: (int(x[2]), int(x[3])) for x in (l.split("\t") for l in te.splitlines())
                      if x[0] == "S" and len(x) >= 4}
                for i, (g, cp, ctg, a, b) in enumerate(meta):
                    pk, mk = sl.get("w%d" % i, (0, 0))
                    rows.append([c, plat, "hap" + h, g, cp, ctg, a, b, pk, mk, origin(pk, mk), strong(pk, mk)])
            # region level: the whole trimmed KIR window
            te = t.sh("%s trioeval -t2 %s %s %s 2>/dev/null" % (t.YAK, pat, mat, fa))
            (sw, hm) = t.parse_trioeval(te)
            reg.append([c, plat, "hap" + h, sum(len(v) for v in seqs.values()), sw[0], sw[1], hm[0], hm[1]])
        return rows, reg, notes, None
    except Exception as ex:
        return [], [], notes, "%s: %s" % (c, ex)
    finally:
        t.sh("rm -rf %s" % wd)


def summaries():
    g = pd.read_csv(INT + "gene_origin_kir.INTERNAL.tsv", sep="\t", dtype={"child": str, "contig": str})
    rows = []
    for gene, d in [("ALL", g)] + list(g.groupby("gene")):
        n = len(d)
        q = {k: C.rate_cells(int((d["origin"] == k).sum()), n) for k in ("amb", "nowin")}
        res = d[d["origin"].isin(["pat", "mat"])]
        st = C.rate_cells(int(res["strong"].sum()), len(res))
        k = (d["pat_k"] + d["mat_k"])[d["origin"] != "nowin"]
        rows.append({"gene": gene, "n_windows": C.mask_n(n), "n_children": C.mask_n(d["child"].nunique()),
                     "ambiguous_pct": q["amb"]["pct"], "ambiguous_n": q["amb"]["n"], "nowin_pct": q["nowin"]["pct"],
                     "strong_of_resolved_pct": st["pct"],
                     "kmers_median": float(k.median()) if len(k) else np.nan,
                     "kmers_p10": float(k.quantile(.1)) if len(k) else np.nan,
                     "kmers_p90": float(k.quantile(.9)) if len(k) else np.nan})
    C.write_tsv(pd.DataFrame(rows), OUT + "gene_origin_summary_kir.tsv")
    # k-mer count distribution per gene (fine bins, for the report's "how much evidence per window" figure)
    k = g[g["origin"] != "nowin"].copy()
    k["tot"] = k["pat_k"] + k["mat_k"]
    k["bin"] = pd.cut(k["tot"], [-1, 0, 4, 9, 19, 49, 99, 199, 499, 10 ** 9],
                      labels=["0", "1-4", "5-9", "10-19", "20-49", "50-99", "100-199", "200-499", "500+"])
    kt = k.groupby(["gene", "bin"], observed=True).size().rename("n").reset_index()
    kt["n"] = [C.mask_n(v) for v in kt["n"]]
    C.write_tsv(kt, OUT + "kmers_per_window_kir.tsv")
    # both haps of a child resolved to the same parent (per gene), as in 57b's summary
    s1 = []
    for gene, d in g[g["origin"].isin(["pat", "mat"]) & (g["copy_index"] == 1)].groupby("gene"):
        piv = d.pivot_table(index="child", columns="hap", values="origin", aggfunc="first")
        both = piv.dropna() if {"hap1", "hap2"} <= set(piv.columns) else piv.iloc[0:0]
        same = int((both["hap1"] == both["hap2"]).sum()) if len(both) else 0
        r = C.rate_cells(same, len(both))
        s1.append({"gene": gene, "n_children_both_resolved": r["d"], "n_same_parent": r["n"], "same_parent_pct": r["pct"]})
    C.write_tsv(pd.DataFrame(s1), OUT + "both_haps_same_parent_kir.tsv")
    rg = pd.read_csv(INT + "region_kir.INTERNAL.tsv", sep="\t", dtype={"child": str})
    rg["sw_pct"] = 100 * rg["sw_err"] / rg["sw_tot"].replace(0, np.nan)
    rg["hm_pct"] = 100 * rg["hm_err"] / rg["hm_tot"].replace(0, np.nan)
    out = rg.groupby("plat").agg(n_hap_regions=("child", "size"), sw_pct_median=("sw_pct", "median"),
                                 hm_pct_median=("hm_pct", "median"), bp_median=("bp", "median")).reset_index()
    out["n_hap_regions"] = [C.mask_n(v) for v in out["n_hap_regions"]]
    C.write_tsv(out, OUT + "region_summary_kir.tsv")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--limit", type=int, help="first N trios (smoke)")
    ap.add_argument("--calls", default=CALLS)
    ap.add_argument("--kir-root", default=C.DEFAULT_KIR_OUTROOT)
    ap.add_argument("--summarise-only", action="store_true")
    a = ap.parse_args(argv)
    for d in (OUT, INT, WORK, WORK + "db"):
        os.makedirs(d, exist_ok=True)
    ck, rk = INT + "gene_origin_kir.INTERNAL.tsv", INT + "region_kir.INTERNAL.tsv"
    if not a.summarise_only:
        t = T()
        calls = pd.read_csv(a.calls, sep="\t", skiprows=1, dtype={"person_id": str, "contig": str})
        have = set(calls["person_id"])
        cohort = C.load_cohort(C.DEFAULT_COHORT)
        lr = set(cohort.loc[cohort["in_lr"].astype(str).str.lower() == "true", "person_id"])
        plat = dict(zip(cohort["person_id"], cohort["platform"]))
        trios, seen = [], set()
        for tr in t.find_trios(lr):
            if tr[0] not in seen and tr[0] in have:     # first candidate per child
                seen.add(tr[0])
                trios.append(tr)
        people = {x for tr in trios for x in tr}
        bam, asm = t.manifest_index(people)
        trios = [tr for tr in trios if tr[0] in asm and tr[1] in bam and tr[2] in bam]
        pd.DataFrame(trios, columns=["child", "pat", "mat"]).to_csv(INT + "trios.INTERNAL.tsv", sep="\t", index=False)
        done = set(pd.read_csv(ck, sep="\t", dtype={"child": str})["child"]) if os.path.exists(ck) else set()
        for f, fields in ((ck, GENE_FIELDS), (rk, REGION_FIELDS)):
            if not os.path.exists(f):
                open(f, "w").write("\t".join(fields) + "\n")
        todo = [tr for tr in trios if tr[0] not in done]
        if a.limit:
            todo = todo[:a.limit]
        remaining = {}
        for tr in todo:
            for who in tr[1:]:
                remaining[who] = remaining.get(who, 0) + 1
        dblock = {w: threading.Lock() for w in remaining}
        dbpath = lambda w: WORK + "db/%s.KIR.yak" % w  # noqa: E731
        C.log("trios %d, to run %d, parent DBs %d" % (len(trios), len(todo), len(remaining)))
        with ThreadPoolExecutor(a.workers) as ex:
            futs = [(tr, ex.submit(run_child, tr, bam, plat.get(tr[0], asm[tr[0]]["plat"]),
                                   calls[calls["person_id"] == tr[0]], a.kir_root, dbpath, dblock)) for tr in todo]
            for k, (tr, fu) in enumerate(futs, 1):
                rows, reg, notes, err = fu.result()
                if not err:    # a failed child is retried by the next --resume
                    for f, rr in ((ck, rows), (rk, reg)):
                        with open(f, "a") as fh:
                            fh.write("".join("\t".join(map(str, r)) + "\n" for r in rr))
                            fh.flush()
                            os.fsync(fh.fileno())
                for who in tr[1:]:
                    remaining[who] -= 1
                    if remaining[who] <= 0 and os.path.exists(dbpath(who)):
                        os.remove(dbpath(who))
                C.log("child %d/%d windows %d" % (k, len(todo), len(rows)) + (" NOTE " + "; ".join(notes) if notes else "")
                      + (" ERROR " + err if err else ""))
    if os.path.exists(ck) and os.path.getsize(ck) > 60:
        summaries()
    C.log("DONE 77b")


if __name__ == "__main__":
    main()
