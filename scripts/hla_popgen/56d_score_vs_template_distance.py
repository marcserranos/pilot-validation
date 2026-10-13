#!/usr/bin/env python3
"""S05 56d: does the 7-point confidence score beat filtering on template distance alone? (Marc, 2026-10-13)

Same records as 56b's typing sweep: first-degree pairs that share an MHC haplotype (<= 1 classical gene with no
shared 2-field allele), one record per pair x gene with both people typed; "error" = the two people share no
2-field allele. Each record gets
  score     = min confidence score over both people's calls (56.prepare, 0-7), as in 56b;
  dist      = max template distance (Immuannot edit distance to the closest known allele) over both people's
              calls of the gene; a missing distance counts as infinite.
Two filters are swept: score >= t (0..7) and dist <= t (thresholds at the 5%..100% quantiles of dist). The fair
comparison is error at the same share of records kept, so the figure plots error against % retained for both.
(A within-gene percentile of distance was tried and dropped: about half of all calls sit at distance 0, so the
tied percentile jumps straight past the low thresholds and the curve is empty there.)

Disclosure: the sweep is cumulative, so two neighbouring exact rows give the records (and errors) between them
by subtraction. A row is dropped when its record or error step from the previous kept row is 1-19.

VM (run from the scripts dir next to 56 and _s05_common):  python3 56d_score_vs_template_distance.py
    -> ~/s05/results/56d/score_vs_distance_sweep.tsv (aggregate; counts 1-19 masked by rate_cells)
Laptop/cloud, from the committed PUBLIC table:           python3 56d_score_vs_template_distance.py --figure
    -> reports/hla_popgen/56d_score_vs_template_distance/fig_score_vs_distance.png
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
PUB = os.path.join(REPO, "reports", "hla_popgen", "56d_score_vs_template_distance")
OUT = os.path.expanduser("~/s05/results/56d/")
QUANTS = [round(q, 2) for q in np.arange(0.05, 1.0001, 0.05)]


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def call_distance(cs):
    """Prepared calls -> {(person, gene): max template distance}; a missing distance counts as infinite."""
    d = cs[["person_id", "gene", "template_distance"]].copy()
    d["td"] = pd.to_numeric(d["template_distance"], errors="coerce").fillna(np.inf)
    return d.groupby(["person_id", "gene"])["td"].max().to_dict()


def records(fd, geno, dist, G8, ALLG):
    """56b's gene records (MHC-sharing pairs), with score and dist per pair x gene."""
    rec = []
    for i, j in zip(fd["i"], fd["j"]):
        fails = [g for g in G8 if geno.get((i, g)) and geno.get((j, g))
                 and not set(geno[(i, g)][0]) & set(geno[(j, g)][0])]
        if len(fails) > 1:
            continue
        for g in ALLG:
            a, b = geno.get((i, g)), geno.get((j, g))
            if a and b and (i, g) in dist and (j, g) in dist:
                rec.append((g, int(not set(a[0]) & set(b[0])), min(a[2], b[2]), max(dist[(i, g)], dist[(j, g)])))
    return pd.DataFrame(rec, columns=["gene", "bad", "score", "dist"])


def drop_differencing(rows, small=20):
    """Cumulative rows of one filter -> drop any row whose record or error step from the previous kept row
    (in order of growing n) is 1..small-1, so no subtraction of neighbours yields a small count."""
    kept, prev = [], (0, 0)
    for r in sorted(rows, key=lambda r: (int(r["_n"]), int(r["_k"]))):
        dn, dk = r["_n"] - prev[0], r["_k"] - prev[1]
        if 0 < dn < small or 0 < dk < small:
            continue
        kept.append(r)
        prev = (r["_n"], r["_k"])
    return kept


def sweep(r, rate_cells):
    """Error and retention for both filters. rate_cells = _s05_common.rate_cells (masks 1-19)."""
    n0 = len(r)
    raw = []

    def add(filt, thr, keep):
        s = r[keep]
        q = rate_cells(int(s["bad"].sum()), len(s))
        small = q["pct"] == "" and len(s) > 0
        raw.append({"filter": filt, "threshold": thr, "n": q["d"], "k": q["n"],
                    "retained_pct": "" if small else (f"{100 * len(s) / n0:.1f}" if n0 else ""),
                    "pct": q["pct"], "lo": q["lo"], "hi": q["hi"], "_n": len(s), "_k": int(s["bad"].sum())})

    for t in range(8):
        add("score_min", t, r["score"] >= t)
    finite = r.loc[np.isfinite(r["dist"]), "dist"]
    thr = sorted({float(finite.quantile(q, interpolation="lower")) for q in QUANTS}) if len(finite) else []
    for t in thr:
        add("template_distance_max", f"{t:g}", r["dist"] <= t)
    add("template_distance_max", "all", np.ones(n0, bool))
    rows = []
    for f in ("score_min", "template_distance_max"):
        rows += drop_differencing([x for x in raw if x["filter"] == f])
    return pd.DataFrame(rows).drop(columns=["_n", "_k"])


def figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    m56c = _load("m56c", "56c_benchmark_figures.py")  # shared style
    s = pd.read_csv(os.path.join(PUB, "score_vs_distance_sweep.tsv"), sep="\t", dtype=str)
    for c in ("retained_pct", "pct", "lo", "hi"):
        s[c] = pd.to_numeric(s[c], errors="coerce")
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    styles = [("score_min", "confidence score ≥ t (7 points)", m56c.BLUE, "o-"),
              ("template_distance_max", "template distance ≤ t (edits to the closest known allele)", m56c.ORANGE, "s-")]
    for f, lab, col, fmt in styles:
        # score 7 keeps 0.6% (CI to 8.8%) and would only stretch the axis; it is in the table and the caption
        d = s[(s["filter"] == f) & s["pct"].notna() & (s["retained_pct"] >= 1)].sort_values("retained_pct")
        ax.errorbar(d["retained_pct"], d["pct"], yerr=[d["pct"] - d["lo"], d["hi"] - d["pct"]], fmt=fmt, color=col,
                    ms=4, lw=1.4, capsize=2, label=lab)
        for x, y, t in zip(d["retained_pct"], d["pct"], d["threshold"]):
            if f == "score_min":
                ax.annotate(f"≥{t}", (x, y), xytext=(-14, 6), textcoords="offset points", fontsize=7, color=col)
            elif t not in ("all", d["threshold"].iloc[-2]):
                ax.annotate(f"≤{t}", (x, y), xytext=(5, -11), textcoords="offset points", fontsize=7, color=col)
    ax.set_xlim(30, 103)
    ax.set_ylim(0, 5)
    ax.set_xlabel("pair-genes kept (%)")
    ax.set_ylabel("pair-genes with no shared allele (%)")
    ax.set_title("Same data kept, which filter removes more errors?")
    ax.legend(loc="upper left", fontsize=8)
    ax.text(102, 0.15, "right = keep everything", ha="right", fontsize=7, color=m56c.INK2)
    os.makedirs(PUB, exist_ok=True)
    fig.savefig(os.path.join(PUB, "fig_score_vs_distance.png"))
    print("wrote", os.path.join(PUB, "fig_score_vs_distance.png"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figure", action="store_true", help="draw the figure from the committed PUBLIC table")
    a = ap.parse_args(argv)
    if a.figure:
        return figure()
    import _s05_common as C
    B = _load("b56", "56_phasing_typing_benchmark.py")
    C.TABLE1_USECOLS = C.TABLE1_USECOLS + ["n_tied"]
    plat = dict(zip(*[C.load_cohort(C.DEFAULT_COHORT)[c] for c in ("person_id", "platform")]))
    t1 = C.load_table1(C.DEFAULT_TABLE1)
    pairs = B.build_pairs(B.load_rel(C.DEFAULT_RELATEDNESS), set(t1["person_id"]), plat)
    fd = pairs[pairs["cls"] != "twin"]
    cs = B.prepare(t1[t1["person_id"].isin(set(fd["i"]) | set(fd["j"]))], plat)
    geno, _ = B.index_calls(cs)
    r = records(fd, geno, call_distance(cs), B.G8, B.ALLG)
    C.log(f"records {len(r)} (56b had 6,455 at score 0); distance missing in {int((~np.isfinite(r['dist'])).sum())}")
    C.write_tsv(sweep(r, C.rate_cells), OUT + "score_vs_distance_sweep.tsv")
    C.log("DONE 56d")


if __name__ == "__main__":
    main()
