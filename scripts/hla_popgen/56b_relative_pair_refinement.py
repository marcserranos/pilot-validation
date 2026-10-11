#!/usr/bin/env python3
"""S05 56b: de-confound 56's first-degree checks (VM; run from ~/s05 next to 56 and _s05_common).
(1) per pair, number of the 8 classical genes with no shared 2-field allele -> bimodal: MHC-IBD0 sibs
    fail at most genes, typing errors at one; pairs with <=1 failing gene = "MHC-sharing".
(2) among MHC-sharing pairs: per-gene non-sharing rate by min confidence score (error + recombination).
(3) haplotype transmission only where BOTH people carry the gene pair on one contig.
Aggregate output: ~/s05/results/56b/*.tsv (counts 1-19 -> <20)."""
import importlib.util
import sys
from collections import Counter

import pandas as pd

import _s05_common as C

spec = importlib.util.spec_from_file_location("b56", "56_phasing_typing_benchmark.py")
B = importlib.util.module_from_spec(spec)
sys.modules["b56"] = B
spec.loader.exec_module(B)
OUT = "/home/jupyter/s05/results/56b/"
C.TABLE1_USECOLS = C.TABLE1_USECOLS + ["n_tied"]
plat = dict(zip(*[C.load_cohort(C.DEFAULT_COHORT)[c] for c in ("person_id", "platform")]))
t1 = C.load_table1(C.DEFAULT_TABLE1)
pairs = B.build_pairs(B.load_rel(C.DEFAULT_RELATEDNESS), set(t1["person_id"]), plat)
fd = pairs[pairs["cls"] != "twin"]
cs = B.prepare(t1[t1["person_id"].isin(set(fd["i"]) | set(fd["j"]))], plat)
geno, hp = B.index_calls(cs)
dist, recs, haps = Counter(), [], []
for i, j, pp in zip(fd["i"], fd["j"], fd["pp"]):
    fails = [g for g in B.G8 if geno.get((i, g)) and geno.get((j, g))
             and not set(geno[(i, g)][0]) & set(geno[(j, g)][0])]
    dist[min(len(fails), 5)] += 1
    if len(fails) > 1:
        continue
    for g in B.ALLG:
        a, b = geno.get((i, g)), geno.get((j, g))
        if a and b:
            recs.append((g, pp, int(not set(a[0]) & set(b[0])), min(a[2], b[2])))
    for ga, gb in B.HAPPAIRS:
        xh, yh = B.hapset(hp, i, ga, gb, True), B.hapset(hp, j, ga, gb, True)
        if len(xh) == 2 and len(yh) == 2:
            sc = min(min(geno[(p, g)][2] for g in (ga, gb)) for p in (i, j))
            haps.append((ga + "~" + gb, pp, int(not any(x == y for x in xh for y in yh)), sc))
rows = [{"n_classical_genes_no_shared_allele": ("5+" if k == 5 else k), "n_pairs": C.mask_n(v)}
        for k, v in sorted(dist.items())]
C.write_tsv(pd.DataFrame(rows), OUT + "pair_sharing_distribution.tsv")
r = pd.DataFrame(recs, columns=["gene", "pp", "bad", "score"])
h = pd.DataFrame(haps, columns=["pair", "pp", "bad", "score"])
out = []
for name, d, key in (("gene_nonsharing", r, "gene"), ("hap_not_transmitted", h, "pair")):
    for item, g in [("ALL", d)] + list(d.groupby(key)):
        for t in range(8):
            s = g[g["score"] >= t]
            q = C.rate_cells(int(s["bad"].sum()), len(s), "public")
            out.append({"check": name, "item": item, "min_score": t, "n": q["d"], "k": q["n"],
                        "pct": q["pct"], "lo": q["lo"], "hi": q["hi"]})
    for pp, g in d.groupby("pp"):
        q = C.rate_cells(int(g["bad"].sum()), len(g), "public")
        out.append({"check": name, "item": "ALL", "min_score": "platform:" + pp, "n": q["d"], "k": q["n"],
                    "pct": q["pct"], "lo": q["lo"], "hi": q["hi"]})
C.write_tsv(pd.DataFrame(out), OUT + "refined_sweep.tsv")
C.log("DONE 56b")
