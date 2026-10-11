#!/usr/bin/env python3
"""S05 56: orthogonal benchmark of HLA TYPING and PHASING quality, with confidence-threshold sweeps.

Pair universe is kinship-only (as in 26): both people have HLA calls in Table 1; kin > 2^-1.5 =
twin/duplicate; 2^-2.5 <= kin <= 2^-1.5 = first degree, split into parent-offspring (PO) vs full
sib by IBS0 (PO: ibs0 < 0.002) or, if only IBD0 exists, ibd0 < 0.1. With neither column the split
is impossible and all first-degree pairs are labelled PO_or_sib (and used as PO, flagged).

  pairs_summary.tsv            pair counts by class x platform combination
  twin_concordance_by_gene.tsv unordered genotype discordance in twin/duplicate pairs (2-field, full)
  novel_replication.tsv        novel calls (depth 2/3/4) whose identical consensus+cds_mut is absent
                               in the co-twin = false-novel-rate estimate
  transmission_consistency.tsv PO: >=1 shared allele per gene (sib contrast); PO haplotype check
                               (same-contig a~b in one person must exist as same-contig a~b in the
                               other when that person carries both alleles) = switch-error proxy;
                               cross-contig control shows how arbitrary hap labels are across contigs
  threshold_sweep.tsv          additive per-call confidence score (0-7) vs error among passing calls

Score = same-contig with phasing partner + MHC contigs on hap <= 2 + template_distance <= 75th pct
        + not homopolymer-indel-only + no template_warning + n_tied <= 1 + platform revio.
A pair-gene (or pair-hap) passes threshold t when the minimum score over its calls is >= t.
Aggregate only; counts 1-19 -> "<20", rates blanked when either side is 1-19 (PUBLIC).

    python3 56_phasing_typing_benchmark.py --out-dir ~/s05/results/56
    python3 56_phasing_typing_benchmark.py --synthetic --out-dir /tmp/s05_56   # self-test
"""
import argparse
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

import _s05_common as C

G8 = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
ALLG = G8 + ["DRB3", "DRB4", "DRB5", "DQA2", "MICA", "TAP1", "TAP2", "C4A", "C4B"]
HAPPAIRS = [("DQA1", "DQB1"), ("DRB1", "DQB1"), ("B", "C"), ("DPA1", "DPB1"), ("A", "B")]
PART = {"A": ["B"], "B": ["A", "C"], "C": ["B"], "DRB1": ["DQB1"], "DQB1": ["DQA1", "DRB1"],
        "DQA1": ["DQB1"], "DPA1": ["DPB1"], "DPB1": ["DPA1"], "DRB3": ["DRB1"],
        "DRB4": ["DRB1"], "DRB5": ["DRB1"]}
KIN_FIRST, KIN_TWIN = 2 ** -2.5, 2 ** -1.5
IBS0_PO_MAX, IBD0_PO_MAX = 0.002, 0.1
CS_TOKEN_RE = re.compile(r"(:\d+|\*[a-z]{2}|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})", re.I)


def is_homopolymer_indel_only(cds_mut):
    """Copied from 50_qc_foundation (same rule as 03_novel_alleles) to keep deploy flat."""
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


# --------------------------------------------------------------------------- pairs
def load_rel(path):
    """Relatedness table -> DataFrame(i, j, kin, ibs0 (or NaN), ibd0 (or NaN)).
    Assumed columns: i.s / j.s (or person_i / person_j, or i / j), kin, and optionally ibs0 / ibd0
    (case-insensitive)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    low = {c.lower(): c for c in df.columns}
    for a, b in [("i.s", "j.s"), ("person_i", "person_j"), ("i", "j")]:
        if a in low and b in low:
            break
    else:
        raise SystemExit("relatedness table: no recognised id columns in %s" % list(df.columns))
    out = pd.DataFrame({"i": df[low[a]], "j": df[low[b]],
                        "kin": pd.to_numeric(df[low["kin"]], errors="coerce")})
    for c in ("ibs0", "ibd0"):
        out[c] = pd.to_numeric(df[low[c]], errors="coerce") if c in low else np.nan
    return out


def plat_pair(a, b):
    return "-".join(sorted([str(a), str(b)]))


def build_pairs(rel, people, platform_of):
    d = rel[rel["i"].isin(people) & rel["j"].isin(people) & (rel["kin"] >= KIN_FIRST)].copy()
    d = d.drop_duplicates(["i", "j"])
    if d["ibs0"].notna().any():
        split, how = d["ibs0"] < IBS0_PO_MAX, "ibs0"
    elif d["ibd0"].notna().any():
        split, how = d["ibd0"] < IBD0_PO_MAX, "ibd0"
    else:
        split, how = None, "none"
        C.log("WARNING: no ibs0/ibd0 column; PO and sibs cannot be separated (class PO_or_sib).")
    cls = np.where(d["kin"] > KIN_TWIN, "twin", "x")
    if split is None:
        cls = np.where(cls == "x", "PO_or_sib", cls)
    else:
        cls = np.where(cls == "x", np.where(split, "PO", "sib"), cls)
    d["cls"] = cls
    d["pp"] = [plat_pair(platform_of.get(i), platform_of.get(j)) for i, j in zip(d["i"], d["j"])]
    d.attrs["split"] = how
    return d


def pairs_summary(pairs, mode):
    rows = []
    for cls, g in pairs.groupby("cls"):
        for pp in ["all"] + sorted(g["pp"].unique()):
            s = g if pp == "all" else g[g["pp"] == pp]
            rows.append({"pair_class": cls, "platform_pair": pp, "n_pairs": C.mask_n(len(s), mode),
                         "median_kin": "%.3f" % s["kin"].median() if len(s) >= 20 else "",
                         "po_sib_split_by": pairs.attrs.get("split", "")})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- call features
def prepare(t1, platform_of):
    t = t1[t1["consensus"].notna() & (t1["consensus"] != "undetermined")].copy()
    t["hap"] = t["hap"].astype(str)
    t["tf"] = [C.two_field(c) or c for c in t["consensus"]]
    t["plat"] = t["person_id"].map(platform_of)
    nct = t.groupby(["person_id", "hap"])["contig"].nunique()
    t["nctg"] = nct.reindex(pd.MultiIndex.from_frame(t[["person_id", "hap"]])).to_numpy()
    c1 = t[t["copy_index"] == 1]
    cof = {(p, h, g): c for p, h, g, c in zip(c1["person_id"], c1["hap"], c1["gene"], c1["contig"])}
    t = t[t["gene"].isin(ALLG)].copy()
    sc = []
    for p, h, g, c in zip(t["person_id"], t["hap"], t["gene"], t["contig"]):
        sc.append(True if g not in PART else any(cof.get((p, h, q)) == c for q in PART[g]))
    td = pd.to_numeric(t["template_distance"], errors="coerce")
    q75 = np.nanquantile(td, 0.75) if td.notna().any() else np.nan
    warn = t["template_warning"].notna() & ~t["template_warning"].astype(str).isin(["", "nan", "None"])
    tied = (pd.to_numeric(t["n_tied"], errors="coerce") > 1) if "n_tied" in t else False
    homo = t["cds_mut"].map(is_homopolymer_indel_only)
    t["score"] = (np.array(sc, int) + (t["nctg"] <= 2).astype(int) + (~(td > q75)).astype(int)
                  + (~homo).astype(int) + (~warn).astype(int) + (~tied).astype(int)
                  + (t["plat"] == "revio").astype(int)).astype(int)
    return t


def index_calls(t):
    geno = {}
    for (p, g), d in t.groupby(["person_id", "gene"]):
        geno[(p, g)] = (tuple(sorted(d["tf"])), tuple(sorted(d["consensus"])), int(d["score"].min()))
    hp = defaultdict(lambda: defaultdict(dict))
    c1 = t[t["copy_index"] == 1]
    for p, h, g, tf, c, s in zip(c1["person_id"], c1["hap"], c1["gene"], c1["tf"], c1["contig"],
                                 c1["score"]):
        hp[p][h][g] = (tf, c, int(s))
    return geno, hp


# --------------------------------------------------------------------------- records
def twin_records(tw, geno):
    rec = []
    for i, j, pp in zip(tw["i"], tw["j"], tw["pp"]):
        for g in ALLG:
            a, b = geno.get((i, g)), geno.get((j, g))
            if a and b:
                rec.append((g, pp, int(a[0] != b[0]), int(a[1] != b[1]), min(a[2], b[2])))
    return pd.DataFrame(rec, columns=["gene", "pp", "d2", "dfull", "score"])


def novel_records(tw, t):
    ks = defaultdict(set)
    for p, g, c, m in zip(t["person_id"], t["gene"], t["consensus"], t["cds_mut"]):
        ks[(p, g)].add((c, str(m)))
    nov = t[t["novelty_depth"].isin([2, 3, 4])]
    by_p = defaultdict(list)
    for r in nov.itertuples(index=False):
        by_p[r.person_id].append(r)
    rec = []
    for i, j, pp in zip(tw["i"], tw["j"], tw["pp"]):
        for a, b in ((i, j), (j, i)):
            for r in by_p.get(a, []):
                if (b, r.gene) in ks:
                    rec.append((r.gene, int(r.novelty_depth), pp,
                                int((r.consensus, str(r.cds_mut)) not in ks[(b, r.gene)]), int(r.score)))
    return pd.DataFrame(rec, columns=["gene", "depth", "pp", "nonrep", "score"])


def hapset(hp, p, a, b, same):
    return [(gd[a][0], gd[b][0]) for gd in hp[p].values()
            if a in gd and b in gd and ((gd[a][1] == gd[b][1]) == same)]


def transmission_records(po, geno, hp):
    rec = []  # kind, item, cls, pp, bad, score
    for i, j, cls, pp in zip(po["i"], po["j"], po["cls"], po["pp"]):
        for g in ALLG:
            a, b = geno.get((i, g)), geno.get((j, g))
            if a and b:
                rec.append(("gene_share", g, cls, pp, int(not (set(a[0]) & set(b[0]))), min(a[2], b[2])))
        if cls == "sib":
            continue
        for ga, gb in HAPPAIRS:
            gs = [geno.get((p, g)) for p in (i, j) for g in (ga, gb)]
            if not all(gs):
                continue
            sc = min(x[2] for x in gs)
            for same, kind in ((True, "hap_same_contig"), (False, "hap_cross_contig")):
                xh, yh = hapset(hp, i, ga, gb, same), hapset(hp, j, ga, gb, same)
                inf = (any(x[0] in gs[2][0] and x[1] in gs[3][0] for x in xh)
                       or any(y[0] in gs[0][0] and y[1] in gs[1][0] for y in yh))
                if not inf:
                    continue
                if same:
                    ok = any(x == y for x in xh for y in yh)
                else:
                    xa = xh + hapset(hp, i, ga, gb, True)
                    ya = yh + hapset(hp, j, ga, gb, True)
                    ok = any(x in ya for x in xh) or any(y in xa for y in yh)
                rec.append((kind, ga + "~" + gb, cls, pp, int(not ok), sc))
    return pd.DataFrame(rec, columns=["kind", "item", "cls", "pp", "bad", "score"])


# --------------------------------------------------------------------------- aggregation
def strata(df):
    yield "all", df
    for pp in sorted(df["pp"].unique()):
        yield pp, df[df["pp"] == pp]


def rrow(prefix, k, n, mode):
    r = C.rate_cells(k, n, mode)
    return {prefix + "_n": r["d"], prefix + "_k": r["n"], prefix + "_pct": r["pct"],
            prefix + "_lo": r["lo"], prefix + "_hi": r["hi"]}


def twin_table(tr, mode):
    rows = []
    for sname, s in strata(tr):
        pooled = s.assign(gene="ALL_GENES")
        for g, d in pd.concat([s, pooled]).groupby("gene"):
            row = {"platform_pair": sname, "gene": g}
            row.update(rrow("disc_2field", int(d["d2"].sum()), len(d), mode))
            r = C.rate_cells(int(d["dfull"].sum()), len(d), mode)
            row.update({"disc_full_k": r["n"], "disc_full_pct": r["pct"]})
            rows.append(row)
    return pd.DataFrame(rows)


def novel_table(nr, mode):
    rows = []
    for sname, s in strata(nr):
        for g in ["ALL"] + sorted(s["gene"].unique()):
            sg = s if g == "ALL" else s[s["gene"] == g]
            for dep in ["all", 2, 3, 4]:
                d = sg if dep == "all" else sg[sg["depth"] == dep]
                if len(d):
                    row = {"platform_pair": sname, "gene": g, "novelty_depth": dep}
                    row.update(rrow("nonreplicating", int(d["nonrep"].sum()), len(d), mode))
                    rows.append(row)
    return pd.DataFrame(rows)


def transmission_table(tr, mode):
    rows = []
    for (kind, cls), g in tr.groupby(["kind", "cls"]):
        for sname, s in strata(g):
            for item in ["ALL"] + sorted(s["item"].unique()):
                d = s if item == "ALL" else s[s["item"] == item]
                row = {"check": kind, "item": item, "pair_class": cls, "platform_pair": sname}
                row.update(rrow("bad", int(d["bad"].sum()), len(d), mode))
                rows.append(row)
    return pd.DataFrame(rows)


def sweep(cs, tw, nr, tr, mode):
    po = tr[tr["cls"].isin(["PO", "PO_or_sib"])]
    gs, hs = po[po["kind"] == "gene_share"], po[po["kind"] == "hap_same_contig"]
    rows, raw = [], []
    for t in range(int(cs["score"].max()) + 1):
        parts = [("calls", cs["score"] >= t, None), ("twin", tw["score"] >= t, tw["d2"]),
                 ("twin_full", tw["score"] >= t, tw["dfull"]),
                 ("novel", nr["score"] >= t, nr["nonrep"]), ("po_gene", gs["score"] >= t, gs["bad"]),
                 ("po_hap", hs["score"] >= t, hs["bad"])]
        row, rw = {"threshold": t}, {"threshold": t}
        for name, m, flag in parts:
            n = int(m.sum())
            k = n if flag is None else int(flag[m].sum())
            if name == "calls":
                r = C.rate_cells(n, len(cs), mode)
                row.update({"calls_retained_n": r["n"], "calls_total": r["d"], "calls_retained_pct": r["pct"]})
            else:
                row.update(rrow(name, k, n, mode))
            rw[name + "_n"], rw[name + "_k"] = n, k
        rows.append(row)
        raw.append(rw)
    return pd.DataFrame(rows), pd.DataFrame(raw)


# --------------------------------------------------------------------------- synthetic
def make_synth(seed):
    """Truth genomes drawn from a haplotype pool; twins = same genome sequenced twice; PO child gets
    one parental hap. Errors (allele flips, spurious novel calls) depend on a latent quality
    (n bad features); hap labels are random across contig segments (as in non-trio hifiasm)."""
    rng = np.random.default_rng(seed)
    gl = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
    pool = [({g: (int(rng.integers(1, 9)), int(rng.integers(1, 4))) for g in gl},
             {g: "n%d_%s" % (k, g) for g in gl if rng.random() < 0.03}) for k in range(40)]
    w = rng.dirichlet(np.ones(40))
    rh = lambda: pool[rng.choice(40, p=w)]
    fnd = [[rh(), rh()] for _ in range(500)]
    people, relrows = [], []

    def add(pid, genome):
        people.append((pid, genome))

    for k, gnm in enumerate(fnd):
        add("F%d" % k, gnm)
    for k in range(80):
        add("T%d" % k, fnd[k])
        relrows.append(("F%d" % k, "T%d" % k, 0.49, 0.0002))
    for k in range(80, 220):
        add("C%d" % k, [fnd[k][rng.integers(2)], rh()])
        relrows.append(("F%d" % k, "C%d" % k, 0.25, rng.uniform(0, 0.0008)))
    for k in range(220, 290):
        for s in "ab":
            add("S%d%s" % (k, s), [fnd[k][rng.integers(2)], rh()])
        relrows.append(("S%da" % k, "S%db" % k, 0.25, rng.uniform(0.01, 0.03)))
    plat, rows, ctr = {}, [], [0]
    for pid, genome in people:
        pl = "revio" if rng.random() < 0.7 else "sequel2e"
        plat[pid] = pl
        nseg = rng.choice([1, 2, 3], p=[.5, .35, .15] if pl == "revio" else [.2, .5, .3])
        cuts = sorted(rng.choice(np.arange(1, 8), nseg - 1, replace=False)) if nseg > 1 else []
        orient = rng.integers(0, 2, nseg)
        for i, g in enumerate(gl):
            seg = int(sum(i >= c for c in cuts))
            for tr_ in (0, 1):
                hh = 1 + (tr_ ^ int(orient[seg]))
                f1, f2 = genome[tr_][0][g]
                novk = genome[tr_][1].get(g)
                warn, tied, tdh = rng.random() < .08, rng.random() < .08, rng.random() < .1
                perr = .005 + .05 * (warn + tied + tdh) + (.02 if pl == "sequel2e" else 0)
                cons, dep, cds = "HLA-%s*%02d:%02d:01:01" % (g, f1, f2), np.nan, np.nan
                if novk:
                    cons, dep, cds = "HLA-%s*%02d:new" % (g, f1), 2, "%s*01:01|:10*ag:20|%s" % (g, novk)
                if rng.random() < perr:
                    ctr[0] += 1
                    if rng.random() < .7:
                        cons, dep, cds = "HLA-%s*%02d:%02d:01:01" % (g, rng.integers(1, 9), rng.integers(1, 4)), np.nan, np.nan
                    else:
                        dep = int(rng.choice([2, 3, 4]))
                        mut = "+aaa" if rng.random() < .5 else "*ag"
                        cons, cds = "HLA-%s*%02d:new" % (g, f1), "%s*01:01|:10%s:20|e%d" % (g, mut, ctr[0])
                rows.append({"person_id": pid, "hap": "hap%d" % hh, "contig": "%s_h%d_c%d" % (pid, hh, seg),
                             "gene": g, "copy_index": 1, "consensus": cons, "novelty_depth": dep,
                             "template_distance": rng.uniform(5, 20) if tdh else rng.uniform(0, 3),
                             "cds_mut": cds, "template_warning": "w" if warn else np.nan,
                             "n_tied": 2 if tied else 1})
    rel = pd.DataFrame(relrows, columns=["i", "j", "kin", "ibs0"])
    rel["ibd0"] = np.nan
    return pd.DataFrame(rows), plat, rel


def synth_check(raw):
    from scipy.stats import spearmanr
    ok = True
    for name, k, n in [("twin", "twin_k", "twin_n"), ("novel", "novel_k", "novel_n"),
                       ("po_hap", "po_hap_k", "po_hap_n")]:
        r = raw[raw[n] >= 30]
        rate = (r[k] / r[n]).to_numpy()
        rho = spearmanr(r["threshold"], rate)[0]
        good = rate[-1] < rate[0] and rho < -0.5
        ok &= bool(good)
        print("SYNTH %-7s rate t=%d: %.3f -> t=%d: %.3f  spearman=%.2f  %s" % (
            name, r["threshold"].iloc[0], rate[0], r["threshold"].iloc[-1], rate[-1], rho,
            "OK" if good else "FAIL"), file=sys.stderr)
    print("SYNTH RESULT:", "PASS" if ok else "FAIL", file=sys.stderr)
    return ok


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=C.DEFAULT_TABLE1)
    ap.add_argument("--cohort-membership", default=C.DEFAULT_COHORT)
    ap.add_argument("--relatedness-table", default=C.DEFAULT_RELATEDNESS)
    ap.add_argument("--out-dir", default=os.path.expanduser("~/s05/results/56"))
    ap.add_argument("--disclosure", default="public", choices=["public", "internal"])
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--seed", type=int, default=3)
    args = ap.parse_args()
    mode = args.disclosure
    if args.synthetic:
        t1, platform_of, rel = make_synth(args.seed)
    else:
        C.TABLE1_USECOLS = C.TABLE1_USECOLS + ["n_tied"]  # loader filters on this module global
        platform_of = dict(zip(*[C.load_cohort(args.cohort_membership)[c] for c in ("person_id", "platform")]))
        t1 = C.load_table1(args.table1, args.limit)
        rel = load_rel(args.relatedness_table)
    C.log("table1 rows=%d people=%d" % (len(t1), t1["person_id"].nunique()))
    pairs = build_pairs(rel, set(t1["person_id"]), platform_of)
    C.log("pairs: %s (split by %s)" % (pairs["cls"].value_counts().to_dict(), pairs.attrs["split"]))
    out = os.path.expanduser(args.out_dir)
    C.write_tsv(pairs_summary(pairs, mode), out + "/pairs_summary.tsv")
    keep = set(pairs["i"]) | set(pairs["j"])
    cs = prepare(t1[t1["person_id"].isin(keep)], platform_of)
    geno, hp = index_calls(cs)
    tw = pairs[pairs["cls"] == "twin"]
    po = pairs[pairs["cls"].isin(["PO", "PO_or_sib", "sib"])]
    twr, nvr, trr = twin_records(tw, geno), novel_records(tw, cs), transmission_records(po, geno, hp)
    C.write_tsv(twin_table(twr, mode), out + "/twin_concordance_by_gene.tsv")
    C.write_tsv(novel_table(nvr, mode), out + "/novel_replication.tsv")
    C.write_tsv(transmission_table(trr, mode), out + "/transmission_consistency.tsv")
    sw, raw = sweep(cs, twr, nvr, trr, mode)
    C.write_tsv(sw, out + "/threshold_sweep.tsv")
    print(raw.to_string(index=False), file=sys.stderr)
    if args.synthetic and not synth_check(raw):
        sys.exit(1)
    C.log("DONE 56")


if __name__ == "__main__":
    main()
