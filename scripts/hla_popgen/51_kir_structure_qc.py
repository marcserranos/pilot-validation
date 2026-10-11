#!/usr/bin/env python3
"""S05 KIR structural QC: the three KIR objections a reviewer will raise first.

  a. framework_by_completeness.tsv  Framework genes (3DL3, 2DL4, 3DP1, 3DL2) should be on ~100% of
     haplotypes; we see 93-96.5%. Split haplotypes into
       complete   = every KIR call of the hap on ONE contig AND both ends (3DL3, 3DL2) present
       one_contig_open_end = one contig but an end gene missing (contig may stop inside the cluster)
       multi_contig = KIR calls spread over >1 contig (assembly fragmentation)
     and report framework presence in each class. If presence is ~100% on complete haps, the
     deficit is assembly incompleteness, not biology.
  b. allelic_pair_events.tsv  2DL2+2DL3 and 3DL1+3DS1 called on one hap. Each event is typed
       overlap        = the two calls overlap in coordinates on one contig (one locus annotated twice
                        = reference cross-mapping, an annotation artefact)
       same_contig    = distinct, non-overlapping loci on one contig (candidate real duplication
                        haplotype; reported with the bp distance between the two loci)
       diff_contig    = on different contigs of the same hap (assembly/hap-assignment artefact)
     plus the same split restricted to complete haplotypes.
  c. gene_content_haplotypes.tsv  gene-content strings (presence/absence of the 17 genes) on
     complete haps; top content haplotypes by ancestry; cA/cB share on complete vs all haps.
  d. kir_same_contig_matrix.tsv  KIR gene x gene same-contig fraction (physical phasing level).

Reads only ~/pipeline_outputs_kir/<pid>/immuannot_output/hap{1,2}.gtf.gz (no re-annotation).
Aggregate-only; participant/haplotype counts 1-19 -> "<20".

    python3 51_kir_structure_qc.py --out-dir ~/s05/results/51 --workers 4
    python3 51_kir_structure_qc.py --synthetic --out-dir /tmp/s05_51
"""
import argparse
import gzip
import os
import re
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations

import numpy as np
import pandas as pd

import _s05_common as C

KIR_GENES = ["KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DP1", "KIR2DL1", "KIR3DP1",
             "KIR2DL4", "KIR3DL1", "KIR3DS1", "KIR2DL5A", "KIR2DS3", "KIR2DS5", "KIR2DS1",
             "KIR2DS4", "KIR2DL5B", "KIR3DL2"]  # approx. cen -> tel order, for readable strings
FRAMEWORK = ["KIR3DL3", "KIR2DL4", "KIR3DP1", "KIR3DL2"]
B_CONTENT = {"KIR2DL2", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS5", "KIR3DS1"}
ALLELIC_PAIRS = [("KIR2DL2", "KIR2DL3"), ("KIR3DL1", "KIR3DS1")]
GENE_RE = re.compile(r'gene_name "([^"]+)"')


def parse_gtf(path):
    """[(contig, gene, start, end)] for transcript rows (one per gene-copy call)."""
    out = []
    if not os.path.exists(path):
        return None
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) < 9 or x[2] != "transcript":
                continue
            m = GENE_RE.search(x[8])
            if m:
                out.append((x[0], C.norm_gene(m.group(1)), int(x[3]), int(x[4])))
    return out


def person_records(args):
    pid, root = args
    recs = []
    for hap in ("hap1", "hap2"):
        calls = parse_gtf(os.path.join(root, pid, "immuannot_output", f"{hap}.gtf.gz"))
        if calls is None:
            continue
        recs.append((pid, hap, calls))
    return recs


def classify_hap(calls):
    genes = {g for _, g, _, _ in calls}
    contigs = {c for c, _, _, _ in calls}
    if not calls:
        return "empty"
    if len(contigs) > 1:
        return "multi_contig"
    if "KIR3DL3" in genes and "KIR3DL2" in genes:
        return "complete"
    return "one_contig_open_end"


def pair_event(calls, ga, gb):
    A = [(c, s, e) for c, g, s, e in calls if g == ga]
    B = [(c, s, e) for c, g, s, e in calls if g == gb]
    if not A or not B:
        return None
    best = None
    for ca, sa, ea in A:
        for cb, sb, eb in B:
            if ca != cb:
                kind, dist = "diff_contig", None
            elif sa <= eb and sb <= ea:
                kind, dist = "overlap", 0
            else:
                kind, dist = "same_contig", max(sb - ea, sa - eb)
            rank = {"overlap": 0, "same_contig": 1, "diff_contig": 2}[kind]
            if best is None or rank < best[0]:
                best = (rank, kind, dist)
    return best[1], best[2]


def run(haps, ancestry_of, unrelated, mode):
    tabs = {}
    # a. framework by completeness
    cls_n = Counter()
    fw = Counter()
    for pid, hap, calls in haps:
        k = classify_hap(calls)
        cls_n[k] += 1
        genes = {g for _, g, _, _ in calls}
        for g in FRAMEWORK:
            fw[(k, g)] += g in genes
        fw[("ALL", "n")] += 1
        for g in FRAMEWORK:
            fw[("ALL", g)] += g in genes
    rows = []
    for k in ["ALL", "complete", "one_contig_open_end", "multi_contig", "empty"]:
        n = cls_n[k] if k != "ALL" else fw[("ALL", "n")]
        row = {"hap_class": k, "n_haps": C.mask_n(n, mode)}
        for g in FRAMEWORK:
            r = C.rate_cells(fw[(k, g)], n, mode, digits=2)
            row[f"{g}_pct"] = r["pct"]
        rows.append(row)
    tabs["framework_by_completeness"] = pd.DataFrame(rows)

    # b. allelic-pair events
    rows = []
    dists = defaultdict(list)
    for ga, gb in ALLELIC_PAIRS:
        for scope in ["all_haps", "complete_haps"]:
            cnt = Counter()
            denom = 0
            for pid, hap, calls in haps:
                if scope == "complete_haps" and classify_hap(calls) != "complete":
                    continue
                genes = {g for _, g, _, _ in calls}
                if ga in genes or gb in genes:
                    denom += 1
                ev = pair_event(calls, ga, gb)
                if ev:
                    cnt[ev[0]] += 1
                    if ev[0] == "same_contig" and scope == "all_haps":
                        dists[(ga, gb)].append(ev[1])
            tot = sum(cnt.values())
            row = {"pair": f"{ga}+{gb}", "scope": scope, "haps_with_either": C.mask_n(denom, mode),
                   "co_called": C.mask_n(tot, mode),
                   "co_called_pct": C.rate_cells(tot, denom, mode, 2)["pct"]}
            for kind in ["overlap", "same_contig", "diff_contig"]:
                row[kind] = C.mask_n(cnt[kind], mode)
                row[f"{kind}_share_pct"] = C.rate_cells(cnt[kind], tot, mode, 1)["pct"]
            if scope == "all_haps" and len(dists[(ga, gb)]) >= 20:
                d = np.array(dists[(ga, gb)])
                row["same_contig_dist_bp_q25_q50_q75"] = "/".join(
                    str(int(v)) for v in np.percentile(d, [25, 50, 75]))
            rows.append(row)
    tabs["allelic_pair_events"] = pd.DataFrame(rows)

    # c. gene-content haplotypes on complete haps (unrelated people)
    content = []
    for pid, hap, calls in haps:
        if pid not in unrelated:
            continue
        genes = {g for _, g, _, _ in calls}
        comp = classify_hap(calls) == "complete"
        s = "".join("1" if g in genes else "0" for g in KIR_GENES)
        content.append((ancestry_of.get(pid, "NA"), comp, s, bool(genes & B_CONTENT)))
    cdf = pd.DataFrame(content, columns=["ancestry", "complete", "content", "is_cB"])
    comp = cdf[cdf["complete"]]
    top = comp["content"].value_counts()
    rows = []
    for rank, (s, n) in enumerate(top.items(), 1):
        if n < 20 and mode == "public":
            break
        row = {"rank": rank, "content_string": s,
               "genes": ",".join(g.replace("KIR", "") for g, b in zip(KIR_GENES, s) if b == "1"),
               "is_cB": any(b == "1" and g in B_CONTENT for g, b in zip(KIR_GENES, s)),
               "n_haps_pooled": C.mask_n(n, mode)}
        for a in C.ANCESTRY_ORDER:
            ca = comp[comp["ancestry"] == a]
            r = C.rate_cells(int((ca["content"] == s).sum()), len(ca), mode, 1)
            row[f"{a}_pct"] = r["pct"]
        rows.append(row)
    rest = int(top[top < 20].sum()) if mode == "public" else 0
    rows.append({"rank": "other(<20 each)", "content_string": "", "genes": "",
                 "n_haps_pooled": C.mask_n(rest, mode),
                 "n_distinct": int((top < 20).sum())})
    tabs["gene_content_haplotypes"] = pd.DataFrame(rows)
    rows = []
    for scope, d in [("all_haps", cdf), ("complete_haps", comp)]:
        for a, g in d.groupby("ancestry"):
            r = C.rate_cells(int(g["is_cB"].sum()), len(g), mode, 1)
            rows.append({"scope": scope, "ancestry": a, "n_haps": r["d"], "cB_pct": r["pct"],
                         "lo": r["lo"], "hi": r["hi"], "n_distinct_content": int(g["content"].nunique())})
    tabs["cB_share_by_completeness"] = pd.DataFrame(rows)

    # d. same-contig matrix
    both, same = Counter(), Counter()
    for pid, hap, calls in haps:
        gc = {}
        for c, g, _, _ in calls:
            gc.setdefault(g, c)
        for g1, g2 in combinations([g for g in KIR_GENES if g in gc], 2):
            both[(g1, g2)] += 1
            same[(g1, g2)] += gc[g1] == gc[g2]
    rows = []
    for (g1, g2), d in both.items():
        r = C.rate_cells(same[(g1, g2)], d, mode, 1)
        rows.append({"gene_a": g1, "gene_b": g2, "n_haps_both": r["d"], "pct_same_contig": r["pct"]})
    tabs["kir_same_contig_matrix"] = pd.DataFrame(rows)
    tabs["hap_class_counts"] = pd.DataFrame(
        [{"hap_class": k, "n_haps": C.mask_n(v, mode)} for k, v in cls_n.items()])
    return tabs


def make_synthetic(n=400, seed=3):
    rng = np.random.default_rng(seed)
    haps = []
    for i in range(n):
        pid = f"P{i:05d}"
        for hap in ("hap1", "hap2"):
            genes = list(FRAMEWORK) + list(rng.choice(KIR_GENES, 5, replace=False))
            if rng.random() < .05:
                genes += ["KIR2DL2", "KIR2DL3"]
            split = rng.random() < .2
            calls = []
            for j, g in enumerate(sorted(set(genes), key=KIR_GENES.index)):
                ctg = "c2" if split and j > 4 else "c1"
                if rng.random() < .03 and g in FRAMEWORK:
                    continue
                calls.append((ctg, g, j * 15000, j * 15000 + 14000))
            haps.append((pid, hap, calls))
    anc = {f"P{i:05d}": rng.choice(C.ANCESTRY_ORDER) for i in range(n)}
    return haps, anc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kir-outroot", default=C.DEFAULT_KIR_OUTROOT)
    ap.add_argument("--cohort-membership", default=C.DEFAULT_COHORT)
    ap.add_argument("--relatedness-table", default=C.DEFAULT_RELATEDNESS)
    ap.add_argument("--out-dir", default=os.path.expanduser("~/s05/results/51"))
    ap.add_argument("--disclosure", default="public", choices=["public", "internal"])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--synthetic", action="store_true")
    args = ap.parse_args()
    if args.synthetic:
        haps, ancestry_of = make_synthetic()
        unrelated = {p for p, _, _ in haps}
    else:
        root = os.path.expanduser(args.kir_outroot)
        pids = sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d)) and d.isdigit())
        if args.limit:
            pids = pids[:args.limit]
        C.log(f"KIR people dirs: {len(pids)}")
        haps = []
        with ProcessPoolExecutor(args.workers) as ex:
            for i, recs in enumerate(ex.map(person_records, [(p, root) for p in pids], chunksize=64)):
                haps.extend(recs)
                if i % 2000 == 0:
                    C.log(f"  parsed {i}/{len(pids)}")
        cohort = C.load_cohort(args.cohort_membership)
        ancestry_of = dict(zip(cohort["person_id"], cohort["ancestry_pred"]))
        unrelated = C.greedy_unrelated({p for p, _, _ in haps}, C.load_relatedness_pairs(args.relatedness_table))
    C.log(f"haplotypes parsed: {len(haps)}; unrelated people: {len(unrelated)}")
    tabs = run(haps, ancestry_of, unrelated, args.disclosure)
    out = os.path.expanduser(args.out_dir)
    for name, df in tabs.items():
        C.write_tsv(df, f"{out}/{name}.tsv")
    C.log("DONE 51")


if __name__ == "__main__":
    main()
