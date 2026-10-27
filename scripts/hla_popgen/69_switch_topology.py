#!/usr/bin/env python3
"""S05 69: where along the MHC do phase switches happen, in whom, and does local coverage explain them?

57 counted switch and hamming errors per trio and target (one number per child haplotype); 57b gave the parent of
origin of each gene. Neither says WHERE on the haplotype the parent flips. Here every child haplotype piece that
aligns to the MHC (chr6:29.5-33.5 Mb, hifiasm -> hg38 PAF, primary blocks) is cut into --win (10 kb) windows along
the contig; each window gets an origin from the parents' k-mers (yak trioeval S line: paternal-only vs maternal-only
k-mer counts; pat/mat when one side has >= 5 and >= 3x the other, else ambiguous; same rule as 57b) and an hg38
position (linear within the PAF block). Along a contig, two consecutive resolved windows with different origins are
a switch at the midpoint between them.

Switch shapes:
  isolated  a run of exactly one window flanked on both sides by the other parent (a short flip: gene-conversion-
            like tract, a collapsed/duplicated segment, or k-mer noise);
  block     the haplotype continues with the other parent for >= 2 windows (a real phase switch: hifiasm joined the
            two parental haplotypes across a homozygous stretch).
Local depth: the child's own HiFi depth (samtools bedcov, MAPQ >= 0) in --bin (100 kb) hg38 bins over the MHC,
divided by the child's genome-wide control depth from 58 (rel depth), compared between bins with and without a
block switch in that child (paired within child).

Trios: 57b's children (internal/57b/gene_origin.INTERNAL.tsv), parents re-found with 68's rule.
PUBLIC (~/s05/results/69/, aggregate): switch_rate_by_bin.tsv (per 250 kb hg38 bin: resolved adjacencies,
block / isolated switches, rate per 1,000 adjacencies; counts 1-19 masked), switch_by_child_summary.tsv
(distribution of block switches per child haplotype: shares, masked), switch_vs_depth.tsv (paired local-depth
contrast; per-child Spearman of block switches vs child/parent MHC depth and number of pieces),
fig_69_switch_topology.png.
INTERNAL (~/s05/internal/69/): per-window origins (checkpoint, resume per child), per-child local depth.
    cd ~/repos/<checkout>/scripts/hla_popgen && setsid nohup $PY -u 69_switch_topology.py --workers 3 \\
        > ~/s05/logs/69.log 2>&1 < /dev/null & disown
    $PY 69_switch_topology.py --summarise-only
Python 3.8 compatible.
"""
import argparse
import importlib.util
import os
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _s05_common as C  # noqa: E402


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


T = load("t57", "57_yak_trioeval_pilot.py")
OUT = os.path.expanduser("~/s05/results/69/")
INT = os.path.expanduser("~/s05/internal/69/")
WORK = os.path.expanduser("~/s05/work69/")
GO = os.path.expanduser("~/s05/internal/57b/gene_origin.INTERNAL.tsv")
DEPTH = os.path.expanduser("~/s05/internal/58/per_person_depth.tsv")
CHROM, S, E = "chr6", 29_500_000, 33_500_000
MIN_K, RATIO, SMALL = 5, 3, 20
WIN_FIELDS = ["child", "plat", "hap", "ctg", "qa", "qb", "ref_pos", "pat_k", "mat_k", "origin"]
DEP_FIELDS = ["child", "bin_start", "mean_depth"]


def origin(pk, mk):
    if pk >= MIN_K and pk >= RATIO * mk:
        return "pat"
    if mk >= MIN_K and mk >= RATIO * pk:
        return "mat"
    return "amb"


def windows_from_paf(paf, win, minlen=1000):
    """[(ctg, qa, qb, ref_pos)] 1-based inclusive contig windows of `win` bp inside primary PAF blocks over the MHC,
    with the hg38 position of the window midpoint (linear within the block, strand-aware)."""
    out = []
    with open(paf) as f:
        for ln in f:
            x = ln.split("\t")
            if x[5] != CHROM or "tp:A:P" not in ln:
                continue
            ql, qs, qe, strand, ts, te = int(x[1]), int(x[2]), int(x[3]), x[4], int(x[7]), int(x[8])
            a, b = max(ts, S), min(te, E)
            if b - a < minlen:
                continue
            f0, f1 = (a - ts) / (te - ts), (b - ts) / (te - ts)
            if strand == "+":
                q0, q1 = qs + f0 * (qe - qs), qs + f1 * (qe - qs)
            else:
                q0, q1 = qe - f1 * (qe - qs), qe - f0 * (qe - qs)
            q0, q1 = max(0, int(q0)), min(ql, int(q1))
            for wa in range(q0, q1 - win // 2, win):
                wb = min(q1, wa + win)
                mid = (wa + wb) / 2.0
                fr = (mid - qs) / (qe - qs)
                ref = ts + (fr if strand == "+" else 1 - fr) * (te - ts)
                out.append((x[0], wa + 1, wb, int(ref)))
    return out


def run_child(trio, bam, asm, win, dbin):
    c, p, m = trio
    wd = WORK + c + "/"
    os.makedirs(wd, exist_ok=True)
    try:
        for who, tag in ((p, "pat"), (m, "mat")):
            T.yak_db(bam[who]["bam"], T.DB["MHC"], wd + f"{tag}.yak", 2)
        wrows = []
        for h in ("1", "2"):
            ws = windows_from_paf(asm[c]["p" + h], win)
            if not ws:
                continue
            fa = wd + f"h{h}.fa.gz"
            T.sh(f"cp '{asm[c]['h' + h]}' {fa}.part && mv {fa}.part {fa} && {T.SAM} faidx {fa}")
            ln = dict(l.split("\t")[:2] for l in open(fa + ".fai"))
            ws = [w for w in ws if w[0] in ln and w[2] <= int(ln[w[0]])]
            open(wd + "w.reg", "w").write("".join(f"{a}:{b}-{e}\n" for a, b, e, _ in ws))
            T.sh(f"{T.SAM} faidx {fa} -r {wd}w.reg > {wd}w.fa")
            te = T.sh(f"{T.YAK} trioeval -t2 {wd}pat.yak {wd}mat.yak {wd}w.fa 2>/dev/null")
            sl = {x[1]: (int(x[2]), int(x[3])) for x in (l.split("\t") for l in te.splitlines())
                  if x[0] == "S" and len(x) >= 4}
            for ctg, a, b, ref in ws:
                pk, mk = sl.get(f"{ctg}:{a}-{b}", (0, 0))
                wrows.append([c, asm[c]["plat"], "hap" + h, ctg, a, b, ref, pk, mk, origin(pk, mk)])
            T.sh(f"rm -f {fa} {fa}.fai {fa}.gzi {wd}w.fa")
        # child local depth in hg38 bins (MAPQ >= 0, primary alignments, like 58)
        bed = wd + "bins.bed"
        open(bed, "w").write("".join(f"{CHROM}\t{s}\t{s + dbin}\n" for s in range(S, E, dbin)))
        bc = T.sh(f"{T.SAM} bedcov -c -Q 0 {bed} '{bam[c]['bam']}'")
        drows = [[c, int(x[1]), int(x[3]) / dbin] for x in (l.split("\t") for l in bc.splitlines()) if len(x) >= 4]
        return wrows, drows, None
    except Exception as ex:  # one bad trio must not stop the run
        return [], [], f"{c}: {ex}"
    finally:
        T.sh(f"rm -rf {wd}")


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------
def call_switches(w):
    """Per child haplotype contig, ordered along the contig, resolved windows only, as runs of one origin.
    An interior run of exactly one window between two runs of the other parent = one isolated flip (at that
    window); every other change of run = one block switch (at the midpoint between the two windows)."""
    ev, adj = [], []
    for (c, h, ctg), g in w[w["origin"] != "amb"].groupby(["child", "hap", "ctg"]):
        g = g.sort_values("qa")
        o, ref = g["origin"].to_numpy(), g["ref_pos"].to_numpy(float)
        for i in range(1, len(o)):
            adj.append((c, h, (ref[i - 1] + ref[i]) / 2.0))
        runs = []  # [origin, first index, last index]
        for i, x in enumerate(o):
            if runs and runs[-1][0] == x:
                runs[-1][2] = i
            else:
                runs.append([x, i, i])
        iso = {j for j in range(1, len(runs) - 1) if runs[j][1] == runs[j][2] and runs[j - 1][0] == runs[j + 1][0]}
        for j in sorted(iso):
            ev.append(dict(child=c, hap=h, ctg=ctg, pos=ref[runs[j][1]], kind="isolated"))
        for j in range(1, len(runs)):
            if j in iso or j - 1 in iso:
                continue
            i = runs[j][1]
            ev.append(dict(child=c, hap=h, ctg=ctg, pos=(ref[i - 1] + ref[i]) / 2.0, kind="block"))
    return (pd.DataFrame(ev, columns=["child", "hap", "ctg", "pos", "kind"]),
            pd.DataFrame(adj, columns=["child", "hap", "pos"]))


def m(n):
    return C.mask_n(int(n))


def rate_by_bin(ev, adj, binw=250_000):
    rows = []
    for s in range(S, E, binw):
        a = int(((adj["pos"] >= s) & (adj["pos"] < s + binw)).sum())
        k = ev[(ev["pos"] >= s) & (ev["pos"] < s + binw)]
        b, i = int((k["kind"] == "block").sum()), int((k["kind"] == "isolated").sum())
        nb_child = k[k["kind"] == "block"]["child"].nunique()
        rows.append(dict(bin_start=s, bin_end=s + binw, n_adjacencies=m(a), n_block_switches=m(b),
                         n_isolated_flips=m(i), n_children_with_block_switch=m(nb_child),
                         block_per_1000=round(1000.0 * b / a, 2) if a >= SMALL and (b == 0 or b >= SMALL) else "",
                         isolated_per_1000=round(1000.0 * i / a, 2) if a >= SMALL and (i == 0 or i >= SMALL) else ""))
    return pd.DataFrame(rows)


def child_summary(ev, w):
    haps = w.groupby(["child", "hap"]).size().reset_index()[["child", "hap"]]
    nb = ev[ev["kind"] == "block"].groupby(["child", "hap"]).size().rename("nb")
    haps = haps.join(nb, on=["child", "hap"]).fillna({"nb": 0})
    tot = len(haps)
    rows = []
    for lab, sel in (("0", haps["nb"] == 0), ("1", haps["nb"] == 1), ("2-3", haps["nb"].between(2, 3)),
                     (">=4", haps["nb"] >= 4)):
        k = int(sel.sum())
        q = C.rate_cells(k, tot)
        rows.append(dict(block_switches_per_haplotype=lab, n_haplotypes=q["n"], pct=q["pct"]))
    t = pd.DataFrame(rows)
    if (t["n_haplotypes"] == "<20").sum() == 1:  # complementary: categories sum to the total
        i = t[t["n_haplotypes"] != "<20"].index[-1]
        t.loc[i, ["n_haplotypes", "pct"]] = ["<20*", ""]
    t.loc[len(t)] = dict(block_switches_per_haplotype="total", n_haplotypes=m(tot), pct="100")
    return t


def vs_depth(ev, dep, ctrl, mhcdep, pieces, trio_of):
    rows = []
    # paired local contrast: per child, mean rel depth of bins holding a block switch vs the other bins
    d = dep.copy()
    d["rel"] = d["mean_depth"] / d["child"].map(ctrl)
    blk = ev[ev["kind"] == "block"]
    diffs = []
    for c, g in d.groupby("child"):
        pos = blk[blk["child"] == c]["pos"].to_numpy()
        if not len(pos) or g["rel"].isna().all():
            continue
        hit = np.zeros(len(g), bool)
        for p in pos:
            hit |= (g["bin_start"].to_numpy() <= p) & (p < g["bin_start"].to_numpy() + (g["bin_start"].diff().median()
                                                                                          or 100_000))
        if hit.any() and (~hit).any():
            diffs.append(g.loc[hit, "rel"].mean() - g.loc[~hit, "rel"].mean())
    diffs = np.array(diffs)
    if len(diffs) >= SMALL:
        from math import erfc, sqrt
        z = diffs.mean() / (diffs.std(ddof=1) / sqrt(len(diffs))) if diffs.std(ddof=1) > 0 else 0.0
        rows.append(dict(test="local_depth_switch_bins_minus_other_bins", n_children=str(len(diffs)),
                         estimate=round(float(diffs.mean()), 4), median=round(float(np.median(diffs)), 4),
                         p=float(f"{erfc(abs(z) / sqrt(2)):.2g}")))
    # per child: number of block switches vs depth / pieces
    nb = blk.groupby("child").size()
    kids = sorted(set(dep["child"]))
    x = pd.DataFrame(dict(child=kids))
    x["nb"] = x["child"].map(nb).fillna(0)
    x["mhc_child"] = x["child"].map(mhcdep)
    x["mhc_parent_min"] = [min(mhcdep.get(trio_of[c][0], np.nan), mhcdep.get(trio_of[c][1], np.nan))
                           if c in trio_of else np.nan for c in x["child"]]
    x["pieces"] = x["child"].map(pieces)
    for col in ("mhc_child", "mhc_parent_min", "pieces"):
        g = x.dropna(subset=[col])
        if len(g) >= SMALL:
            rows.append(dict(test=f"spearman_block_switches_vs_{col}", n_children=str(len(g)),
                             estimate=round(float(g["nb"].rank().corr(g[col].rank())), 3), median="", p=""))
    return pd.DataFrame(rows)


def figure(rb, genes, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 3.6))
    x = (rb["bin_start"] + rb["bin_end"]) / 2e6
    for col, lab, colr in (("block_per_1000", "block switch", "C3"), ("isolated_per_1000", "isolated flip", "C0")):
        y = pd.to_numeric(rb[col], errors="coerce")
        ax.step(x, y, where="mid", color=colr, label=lab)
        ax.plot(x[y.isna()], np.zeros(int(y.isna().sum())), "x", color=colr, ms=4)
    ymax = ax.get_ylim()[1]
    for g, s0 in genes:
        ax.axvline(s0 / 1e6, color="0.85", lw=0.6, zorder=0)
        ax.text(s0 / 1e6, ymax * 1.01, g, rotation=90, fontsize=6, ha="center", va="bottom")
    ax.set_xlabel("chr6 position (Mb, GRCh38)")
    ax.set_ylabel("per 1,000 resolved adjacencies")
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.set_title("x = masked (fewer than 20 events) or no data", fontsize=8, loc="right")
    fig.tight_layout()
    fig.savefig(out + "fig_69_switch_topology.png", dpi=200)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--win", type=int, default=10_000)
    ap.add_argument("--depth-bin", type=int, default=100_000)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--summarise-only", action="store_true")
    ap.add_argument("--out-dir", default=OUT)
    a = ap.parse_args(argv)
    out = a.out_dir.rstrip("/") + "/"
    for d in (out, INT, WORK):
        os.makedirs(d, exist_ok=True)
    M68 = load("m68", "68_trio_error_vs_depth.py")
    kids = sorted(set(pd.read_csv(GO, sep="\t", dtype={"child": str})["child"]))
    trio_of = M68.find_trios(kids, C.load_relatedness_pairs(C.DEFAULT_RELATEDNESS))
    trios = [(c,) + trio_of[c] for c in kids if c in trio_of][: a.limit]
    wck, dck = INT + "windows.INTERNAL.tsv", INT + "local_depth.INTERNAL.tsv"
    if not a.summarise_only:
        people = {x for t in trios for x in t}
        bam, asm = T.manifest_index(people)
        done = set(pd.read_csv(dck, sep="\t", dtype={"child": str})["child"]) if os.path.exists(dck) else set()
        for path, cols in ((wck, WIN_FIELDS), (dck, DEP_FIELDS)):
            if not os.path.exists(path):
                open(path, "w").write("\t".join(cols) + "\n")
        todo = [t for t in trios if t[0] not in done and t[0] in asm and t[0] in bam and t[1] in bam and t[2] in bam]
        C.log(f"trios {len(trios)}, to run {len(todo)}, {a.workers} workers")
        with ThreadPoolExecutor(a.workers) as ex:
            futs = [ex.submit(run_child, t, bam, asm, a.win, a.depth_bin) for t in todo]
            for k, fu in enumerate(futs, 1):
                wr, dr, err = fu.result()
                if wr and dr:
                    for path, rows in ((wck, wr), (dck, dr)):  # depth last = the 'done' marker
                        with open(path, "a") as f:
                            f.write("".join("\t".join(map(str, r)) + "\n" for r in rows))
                            f.flush()
                            os.fsync(f.fileno())
                C.log(f"child {k}/{len(todo)} windows {len(wr)}" + (f" ERROR {err}" if err else ""))
    w = pd.read_csv(wck, sep="\t", dtype={"child": str, "ctg": str})
    dep = pd.read_csv(dck, sep="\t", dtype={"child": str})
    w = w[w["child"].isin(set(dep["child"]))]
    ev, adj = call_switches(w)
    ev.to_csv(INT + "switch_events.INTERNAL.tsv", sep="\t", index=False)
    rb = rate_by_bin(ev, adj)
    C.write_tsv(rb, out + "switch_rate_by_bin.tsv")
    C.write_tsv(child_summary(ev, w), out + "switch_by_child_summary.tsv")
    d58 = pd.read_csv(DEPTH, sep="\t", dtype={"person_id": str}, comment="#")
    d58 = d58[pd.to_numeric(d58["mapq"], errors="coerce") == 0]
    ctrl = d58.groupby("person_id")["control_depth"].first().astype(float).to_dict()
    mhc = d58[d58["region"] == "window_MHC"].set_index("person_id")["rel_depth"].astype(float).to_dict()
    pieces = w.groupby("child")["ctg"].nunique().to_dict()
    C.write_tsv(vs_depth(ev, dep, ctrl, mhc, pieces, trio_of), out + "switch_vs_depth.tsv")
    genes = []
    bed = os.path.join(os.path.dirname(HERE), "..", "reports", "hla_popgen", "58_regional_depth", "regions.bed")
    for p in (bed, os.path.expanduser("~/s05/results/58/regions.bed")):
        if os.path.exists(p):
            r = pd.read_csv(p, sep="\t", header=None, comment="#")
            r = r[(r[0] == CHROM) & r[3].astype(str).isin(["HLA-" + g for g in ("A", "C", "B", "DRB1", "DQA1", "DQB1",
                                                                                 "DPA1", "DPB1", "E", "G", "F")]
                                                          + ["MICA", "TAP1"])]
            genes = list(zip(r[3].astype(str).str.replace("HLA-", ""), r[1].astype(int)))
            break
    figure(rb, genes, out)
    C.log(f"DONE 69: {dep['child'].nunique()} children, {len(ev)} switch events")
    return 0


if __name__ == "__main__":
    sys.exit(main())
