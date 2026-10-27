#!/usr/bin/env python3
"""S06 76: the three call outcomes, kept apart (VM; aggregate output + INTERNAL unit tables for 76b).

71/71e showed that most of what the trio and pair tests call a "typing error" is a MISSING call (one haplotype without a gene annotation), not a wrong
allele among a fully called genotype. This script separates them:

  STREAM A  missing call.   unit = person x gene x haplotype slot; label = no gene annotation on that haplotype. Whole cohort. Features of the
            haplotype that lacks the call: its contig structure (gene between two contigs / at the end of the region / inside one contig), how many
            genes it carries, depth at the gene, the other haplotype's call. Answers "why is the call missing" before any pruning is asked.
  STREAM B  wrong allele.   unit = one child haplotype call, with the transmitting parent resolved by yak (57b) AND the child fully called at the gene
            AND the transmitting parent fully called. Three-way class: `match` (2-field name found in the transmitting parent), `other_parent_only`
            (not in the transmitting parent, found in the other: yak origin error or recombination, NOT a typing error), `neither` (a typing-error
            candidate). DRB3/4/5 are reported apart (biological absence); calls novel in field 1-2 are not labelled (stream C). Origin robustness
            (k-mer margin, strong-origin subset, neighbouring-gene origin switch) and trio verification (parent-child IBD0, number of candidate trios)
            are columns, so the subset is chosen at analysis time and every subset's count is in the tables.
  STREAM C  novel / undetermined calls. Rates by gene and by depth / missingness, twin concordance; no error label.

Outputs
  INTERNAL (~/s06/internal/76/, never leaves the Workbench): slots, units_A, trio_calls, units_B, units_B_all, units_combined, units_pair_all,
           units_pair_complete (all `.INTERNAL.tsv`; 76b reads these)
  AGGREGATE (~/s06/results/76/): definition_card.tsv, A_*.tsv, B_*.tsv, C_*.tsv

VM:
    cd ~/repos/pilot-validation-main && pixi run -e spechla -- python3 scripts/hla_popgen/76_call_outcomes.py \
        --internal-71 ~/s06/internal/71_named_v4/ --out-dir ~/s06/results/76/ --internal-dir ~/s06/internal/76/
~10-20 min (Table 1 read ~ 5 min, per-slot loops dominate). Tests: scripts/hla_popgen/tests/test_76_call_outcomes.py (synthetic).
"""
import argparse
import importlib.util
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _s05_common as C  # noqa: E402


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


Q = sys.modules.get("m71") or _load("m71", "71_call_quality_table.py")
S = sys.modules.get("m71b") or _load("m71b", "71b_call_quality_stats.py")
E = sys.modules.get("m71e") or _load("m71e", "71e_dropout_outcome.py")

G8 = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
BIO = ["DRB3", "DRB4", "DRB5"]            # absent on many haplotypes by biology, never read as missing calls
ANCHORS = ["A", "C", "B", "MICA", "DRA", "DRB1", "DQA1", "DQB1", "TAP2", "TAP1", "DPA1", "DPB1"]   # genomic order, fixed-content genes
INT71 = os.path.expanduser("~/s06/internal/71_named_v4/")
INT = os.path.expanduser("~/s06/internal/76/")
OUT = os.path.expanduser("~/s06/results/76/")
OFFSET_WINDOW = 100_000           # neighbouring gene = within 100 kb on the same contig (switch-suspect check)
MIN_K, RATIO = 5, 3               # 57b's own rule (default origin)
STRONG_K, STRONG_RATIO = 10, 5    # strong-origin subset
IBD0_MAX = 0.1                    # parent-offspring share ~0 IBD0 segments; siblings ~0.25


def log(m):
    C.log("[76] " + m)


def hapn(h):
    """'hap1' / '1' / 1 -> '1'."""
    return str(h).replace("hap", "")


def wilson_cols(k, n):
    p, lo, hi = S.wilson_arr(np.asarray([k], float), np.asarray([n], float))
    return 100 * float(p[0]), 100 * float(lo[0]), 100 * float(hi[0])


def rate_row(k, n, **kw):
    p, lo, hi = wilson_cols(k, n) if n else (np.nan, np.nan, np.nan)
    return dict(n=int(n), k=int(k), pct=p, lo=lo, hi=hi, **kw)


# ===========================================================================
# STREAM A: missing calls (haplotype slots)
# ===========================================================================
def slot_table(t1, genes):
    """Table 1 -> one row per person x gene x haplotype slot (hap '1'/'2') for the persons' genes in `genes`.

    status per slot: called | undetermined | missing. For a missing slot the haplotype's structure around the gene (anchor genes in genomic order):
      hap_class   no_haplotype      the haplotype has no rows at all
                  no_flank          the haplotype carries no anchor gene on either side
                  region_end        anchors only on one side: the gene lies beyond the called part of the region (contig or region end)
                  between_contigs   the nearest anchors on each side sit on DIFFERENT contigs: an assembly break at the gene
                  inside_contig     both neighbours on one contig: gene silently absent inside a contiguous stretch (filtered by Immuannot, divergent, or deleted)
                  not_anchored      gene without a fixed position (DRB3/4/5): no class
    The same class is computed for CALLED slots too (as if the gene were absent), so rates P(missing | class) are meaningful."""
    t = t1.copy()
    t["hap"] = t["hap"].map(hapn)
    t = t[t["hap"].isin(["1", "2"])]
    ctg = t["contig"].astype(str).str.partition(":")[0]
    t["ctg"] = ctg
    cons = t["consensus"].astype(str)
    t["undet"] = cons.eq("undetermined")
    persons = sorted(t["person_id"].unique())
    on_hap = defaultdict(dict)          # (person, hap) -> gene -> ctg
    undet = set()
    for p, h, g, c, u in zip(t["person_id"], t["hap"], t["gene"], t["ctg"], t["undet"]):
        d = on_hap[(p, h)]
        if g not in d:
            d[g] = c
        if u:
            undet.add((p, h, g))
    rows = []
    n_ctg = {k: len(set(d.values())) for k, d in on_hap.items()}
    for p in persons:
        for h in ("1", "2"):
            d = on_hap.get((p, h), {})
            anchors_called = [(i, d[a]) for i, a in enumerate(ANCHORS) if a in d]
            for g in genes:
                if (p, h) not in on_hap:
                    cls = "no_haplotype"
                elif g not in ANCHORS:
                    cls = "not_anchored"
                else:
                    gi = ANCHORS.index(g)
                    left = [c for i, c in anchors_called if i < gi]
                    right = [c for i, c in anchors_called if i > gi]
                    if not left and not right:
                        cls = "no_flank"
                    elif not left or not right:
                        cls = "region_end"
                    else:
                        cls = "inside_contig" if left[-1] == right[0] else "between_contigs"
                if g in d:
                    status = "undetermined" if (p, h, g) in undet else "called"
                else:
                    status = "missing"
                rows.append((p, g, h, status, cls, len(d), n_ctg.get((p, h), 0)))
    return pd.DataFrame(rows, columns=["person", "gene", "hap", "status", "hap_class", "hap_n_genes", "hap_n_contigs"])


def attach_slot_features(slots, ct):
    """Depth / MHC / other-haplotype features. Gene depth is a person x gene quantity, taken from any call the person has at the gene (a person with
    no call at the gene on either haplotype has no depth here: that case is counted, not imputed)."""
    c = ct.copy()
    c["hap"] = c["hap"].map(hapn)
    pg = c.groupby(["person_id", "gene"]).agg(x_rel_cds=("x_rel_cds", "max"), x_cds_share_lt10=("x_cds_share_lt10", "max"),
                                              x_rel_mhc=("x_rel_mhc", "max"), x_cds_min_rel=("x_cds_min_rel", "max")).reset_index()
    pg = pg.rename(columns={"person_id": "person"})
    s = slots.merge(pg, on=["person", "gene"], how="left")
    called = slots[slots["status"] != "missing"][["person", "gene", "hap"]].assign(c=1)
    other = called.assign(hap=called["hap"].map({"1": "2", "2": "1"})).rename(columns={"c": "other_called"})
    s = s.merge(other, on=["person", "gene", "hap"], how="left")
    s["other_called"] = s["other_called"].fillna(0).astype(int)
    s["missing"] = (s["status"] == "missing").astype(int)
    return s


def bin_table(s, feature, edges, scope):
    v = s[feature].to_numpy(float)
    ok = np.isfinite(v)
    ix = np.digitize(v[ok], edges) - 1
    ix = np.clip(ix, 0, len(edges) - 2)
    m = s["missing"].to_numpy(int)[ok]
    rows = []
    n = np.bincount(ix, minlength=len(edges) - 1)
    k = np.bincount(ix, weights=m, minlength=len(edges) - 1)
    for j in range(len(edges) - 1):
        if n[j]:
            rows.append(dict(scope=scope, feature=feature, bin_lo=float(edges[j]), bin_hi=float(edges[j + 1]), **rate_row(k[j], n[j])))
    return rows


def analyse_slots(s, out_dir, boot=200):
    """Stream A tables. `s` is attach_slot_features output."""
    w = lambda df, name: C.write_tsv(df, os.path.join(out_dir, name))  # noqa: E731
    scopes = {"classical8": s[s["gene"].isin(G8)], "all_fixed_genes": s[~s["gene"].isin(BIO)], "drb345_biology": s[s["gene"].isin(BIO)]}
    by_gene, by_class, bins, auc, explained = [], [], [], [], []
    for sc, d in scopes.items():
        for g, dg in d.groupby("gene"):
            by_gene.append(dict(scope=sc, gene=g, **rate_row(dg["missing"].sum(), len(dg))))
        by_gene.append(dict(scope=sc, gene="ALL", **rate_row(d["missing"].sum(), len(d))))
        for cls, dc in d.groupby("hap_class"):
            by_class.append(dict(scope=sc, hap_class=cls, **rate_row(dc["missing"].sum(), len(dc))))
        bins += bin_table(d, "hap_n_genes", np.arange(0, 20, 1), sc)
        bins += bin_table(d, "hap_n_contigs", np.arange(0, 15, 1), sc)
        for f in ("x_rel_cds", "x_rel_mhc", "x_cds_min_rel", "x_cds_share_lt10"):
            bins += bin_table(d, f, np.round(np.arange(0, 2.0001, 0.05), 2), sc)
        # AUC of each feature for "slot missing" (oriented so that higher = more likely missing)
        feats = {"hap_n_genes": -d["hap_n_genes"], "hap_n_contigs": d["hap_n_contigs"], "x_rel_cds": -d["x_rel_cds"],
                 "x_rel_mhc": -d["x_rel_mhc"], "x_cds_min_rel": -d["x_cds_min_rel"], "x_cds_share_lt10": d["x_cds_share_lt10"],
                 "other_called": -d["other_called"], "class_between_or_end": d["hap_class"].isin(["between_contigs", "region_end", "no_flank"]).astype(int)}
        y = d["missing"].to_numpy(int)
        for f, v in feats.items():
            v = np.asarray(v, float)
            a = S.auc_score(v, y)
            lo, hi = S.cluster_boot_auc(v, y, d["person"].to_numpy(), B=boot) if np.isfinite(a) else (np.nan, np.nan)
            auc.append(dict(scope=sc, feature=f, n=int(np.isfinite(v).sum()), n_missing=int(y[np.isfinite(v)].sum()), auc=a, lo=lo, hi=hi))
        # what explains the missing slots (mutually exclusive, in order): no haplotype > assembly structure > low depth > unexplained
        m = d[d["missing"] == 1]
        low_depth = (m["x_rel_cds"] < 0.3) | (m["x_cds_share_lt10"] > 0.5)
        expl = np.select([m["hap_class"] == "no_haplotype", m["hap_class"].isin(["between_contigs"]), m["hap_class"].isin(["region_end", "no_flank"]),
                          (m["hap_class"] == "inside_contig") & low_depth, m["hap_class"] == "inside_contig"],
                         ["haplotype_absent", "assembly_break_at_gene", "region_end", "inside_contig_low_depth", "inside_contig_depth_ok"],
                         default="not_anchored_gene")
        for cat, n in pd.Series(expl).value_counts().items():
            explained.append(dict(scope=sc, explanation=cat, n_missing_slots=int(n), pct_of_missing=100 * n / max(len(m), 1)))
    w(pd.DataFrame(by_gene), "A_missing_by_gene.tsv")
    w(pd.DataFrame(by_class), "A_missing_by_hap_class.tsv")
    w(pd.DataFrame(bins), "A_missing_by_feature_bins.tsv")
    w(pd.DataFrame(auc), "A_missing_feature_auc.tsv")
    w(pd.DataFrame(explained), "A_missing_explained.tsv")


def person_gene_units(ct, genes_scope):
    """Cohort person x gene units for pruning (71e's unit: metrics of one random call, so the metric cannot reveal the call count)."""
    pg = E.person_gene_table(ct)
    pg["scope_classical8"] = pg["gene"].isin(G8)
    pg = pg[~pg["gene"].isin(BIO)].reset_index(drop=True)
    pg["bad"] = pg["incomplete"]
    return pg


# ===========================================================================
# STREAM B: wrong allele among fully called genotypes
# ===========================================================================
def three_field(c):
    """'HLA-A*01:01:01:02N' -> 'A*01:01:01' (protein + synonymous = identical CDS in IPD naming); None if not determinable."""
    if not isinstance(c, str) or "*" not in c:
        return None
    g, rest = C.norm_gene(c).split("*", 1)
    f = rest.split(":")
    if len(f) < 3 or "new" in f[:3]:
        return None
    return "%s*%s:%s:%s" % (g, f[0], f[1], "".join(ch for ch in f[2] if ch.isdigit()) or f[2])


def trio_meta(trios, rel=None, ibd0_max=IBD0_MAX):
    """trios: list of (child, a, b) from 57 find_trios (candidates, in order). -> trio_of (first candidate per child, as 57b) and a per-child
    DataFrame: n_candidates, verified ('yes' | 'no' | 'unverifiable') from the parent-child IBD0 (or IBS0 when only that exists), family id.
    Verified = both child-parent links have IBD0 <= ibd0_max (parent-offspring ~0; full sibling ~0.25, grandparent/half-sibling ~0.5)."""
    ncand = defaultdict(int)
    trio_of = {}
    for t in trios:
        ncand[t[0]] += 1
        trio_of.setdefault(t[0], (t[1], t[2]))
    link = {}
    if rel is not None and len(rel):
        col = "ibd0" if rel["ibd0"].notna().any() else ("ibs0" if rel["ibs0"].notna().any() else None)
        if col:
            for i, j, v in zip(rel["i"], rel["j"], rel[col]):
                link[(i, j)] = link[(j, i)] = float(v)
    rows = []
    for c, (a, b) in trio_of.items():
        va, vb = link.get((c, a)), link.get((c, b))
        if va is None or vb is None or not np.isfinite(va) or not np.isfinite(vb):
            ver = "unverifiable"
        else:
            ver = "yes" if max(va, vb) <= ibd0_max else "no"
        rows.append((c, ncand[c], ver))
    meta = pd.DataFrame(rows, columns=["child", "n_candidates", "verified"])
    # family = connected component over everyone who appears in a (first-candidate) trio: siblings share parents
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for c, (a, b) in trio_of.items():
        for x in (a, b):
            parent[find(x)] = find(c)
    roots = {r: i for i, r in enumerate(sorted({find(c) for c in trio_of}))}
    meta["family"] = [roots[find(c)] for c in meta["child"]]
    meta["primary"] = (meta["n_candidates"] == 1) & (meta["verified"] != "no")
    return trio_of, meta


def load_positions(t1):
    """(person, 'hapN', gene) -> (contig name, absolute midpoint) for copy 1; gene_start/end are relative to the extraction range `ctg:S-E`."""
    t = t1[t1["copy_index"] == 1]
    ctg = t["contig"].astype(str).str.partition(":")
    rng = ctg[2].str.partition("-")[0]
    off = pd.to_numeric(rng, errors="coerce").fillna(1) - 1
    mid = off + (pd.to_numeric(t["gene_start"], errors="coerce") + pd.to_numeric(t["gene_end"], errors="coerce")) / 2
    return {(p, "hap" + hapn(h), g): (c, m) for p, h, g, c, m in zip(t["person_id"], t["hap"], t["gene"], ctg[0], mid)}


def switch_suspects(origin, positions=None):
    """origin rows -> set of (child, hap, gene) whose window origin disagrees with a resolved neighbouring gene on the same child contig:
    within OFFSET_WINDOW bp when positions are known, else the adjacent anchor gene(s) in genomic order. A suspect window may straddle a
    parent switch, so its assigned parent is unreliable."""
    grp = defaultdict(list)
    for c, h, g, ctg, og in zip(origin["child"], origin["hap"], origin["gene"], origin["contig"], origin["origin"]):
        grp[(c, h, str(ctg))].append((g, og))
    sus = set()
    for (c, h, ctg), lst in grp.items():
        for g, og in lst:
            if og not in ("pat", "mat"):
                continue
            for g2, og2 in lst:
                if g2 == g or og2 not in ("pat", "mat") or og2 == og:
                    continue
                if positions is not None and (c, h, g) in positions and (c, h, g2) in positions:
                    near = abs(positions[(c, h, g)][1] - positions[(c, h, g2)][1]) <= OFFSET_WINDOW
                else:
                    near = g in ANCHORS and g2 in ANCHORS and abs(ANCHORS.index(g) - ANCHORS.index(g2)) == 1
                if near:
                    sus.add((c, h, g))
                    break
    return sus


def classify_trio_calls(ct, origin, trio_of, meta, positions=None):
    """One row per child haplotype call with a resolved parent of origin (57b origin pat/mat).

    status   eligible            child fully called at the gene, transmitting parent fully called, child call named  -> labelled
             child_incomplete    child has one call at the gene (the missing second call is stream A)
             parent_incomplete   transmitting parent has 0 or 1 call at the gene (a missing parent call cannot be told from a wrong allele)
             novel_unlabelled    child's call, or a call of the parent at the gene, is novel in field 1-2 (stream C)
    DRB3/4/5: completeness is not required (status eligible when both people have >=1 call and the parent has one) and `bio` = True; kept out of the main tables.
    cls      match | other_parent_only | neither  (eligible rows only)"""
    c0 = ct.copy()
    c0["hap"] = c0["hap"].map(lambda h: "hap" + hapn(h))
    ncalls = c0.groupby(["person_id", "gene"])["hap"].nunique().to_dict()
    two = {k: set(d["tf"].dropna()) for k, d in c0.groupby(["person_id", "gene"])}
    cons = {k: set(map(three_field, d["consensus"])) - {None} for k, d in c0.groupby(["person_id", "gene"])}
    unnamed = set(zip(c0.loc[~c0["named"], "person_id"], c0.loc[~c0["named"], "gene"]))
    key = c0.set_index(["person_id", "hap", "gene"])
    tf = key["tf"].to_dict()
    cs = key["consensus"].to_dict()
    named = key["named"].to_dict()
    plat = key["plat"].to_dict()
    ocols = ["o_" + k for k in Q.KEYS]
    pg = Q.person_gene_worst(c0)
    pgi = {k: i for i, k in enumerate(pg.index)}
    pgv = pg.to_numpy(float)
    cv = key[ocols]
    cvi = {k: i for i, k in enumerate(cv.index)}
    cvv = cv.to_numpy(float)
    sus = switch_suspects(origin, positions)
    mm = meta.set_index("child")
    rows = []
    for c, h, g, og, pk, mk in zip(origin["child"], origin["hap"], origin["gene"], origin["origin"], origin["pat_k"], origin["mat_k"]):
        if og not in ("pat", "mat") or c not in trio_of or (c, h, g) not in tf:
            continue
        a, b = trio_of[c]
        par, oth = (a, b) if og == "pat" else (b, a)
        bio = g in BIO
        n_c, n_p, n_o = ncalls.get((c, g), 0), ncalls.get((par, g), 0), ncalls.get((oth, g), 0)
        if n_p == 0:
            status = "parent_incomplete"
        elif (not named.get((c, h, g), False)) or (par, g) in unnamed or (oth, g) in unnamed:
            status = "novel_unlabelled"
        elif not bio and n_c < 2:
            status = "child_incomplete"
        elif not bio and n_p < 2:
            status = "parent_incomplete"
        else:
            status = "eligible"
        t = tf[(c, h, g)]
        cls = ""
        if status == "eligible":
            cls = "match" if t in two.get((par, g), set()) else ("other_parent_only" if t in two.get((oth, g), set()) else "neither")
        win, lose = max(pk, mk), min(pk, mk)
        ratio = win / max(lose, 1)
        row = dict(person=c, other=par, other2=oth, gene=g, hap=h, origin=og, plat=plat.get((c, h, g), "NA"), status=status, cls=cls, bio=bio,
                   n_child_calls=n_c, n_parent_calls=n_p, n_other_calls=n_o,
                   child_homozygous=int(len(two.get((c, g), set())) == 1 and n_c == 2),
                   cds_differs_given_name_match=int(cls == "match" and three_field(cs[(c, h, g)]) not in cons.get((par, g), set())
                                                    and three_field(cs[(c, h, g)]) is not None and bool(cons.get((par, g)))),
                   k_win=int(win), k_lose=int(lose), ratio=float(ratio),
                   strong=int(win >= STRONG_K and ratio >= STRONG_RATIO), switch_suspect=int((c, h, g) in sus),
                   family=int(mm.loc[c, "family"]), n_candidates=int(mm.loc[c, "n_candidates"]), verified=mm.loc[c, "verified"],
                   primary=bool(mm.loc[c, "primary"]))
        i = cvi[(c, h, g)]
        pv = pgv[pgi[(par, g)]] if (par, g) in pgi else np.full(len(ocols), np.nan)
        for j, k in enumerate(Q.KEYS):
            row["c_" + k], row["p_" + k], row["w_" + k] = cvv[i, j], pv[j], np.fmax(cvv[i, j], pv[j])
        rows.append(row)
    d = pd.DataFrame(rows)
    d["bad"] = (d["cls"] == "neither").astype(int)
    return d


def _subsets(d):
    """Origin-robustness subsets as boolean masks."""
    return {"default": np.ones(len(d), bool), "strong_origin": d["strong"].to_numpy() == 1,
            "no_switch_suspect": d["switch_suspect"].to_numpy() == 0,
            "strong_and_no_switch": (d["strong"].to_numpy() == 1) & (d["switch_suspect"].to_numpy() == 0)}


def summarise_trio(rec, out_dir):
    w = lambda df, name: C.write_tsv(df, os.path.join(out_dir, name))  # noqa: E731
    main = rec[~rec["bio"]]
    # status accounting: where does every resolved child call go?
    st = []
    for sc, d in (("classical8", main[main["gene"].isin(G8)]), ("all_fixed_genes", main), ("drb345_biology", rec[rec["bio"]])):
        for s_, n in d["status"].value_counts().items():
            st.append(dict(scope=sc, status=s_, n_child_calls=int(n), pct=100 * n / len(d)))
    w(pd.DataFrame(st), "B_status_counts.tsv")
    el = rec[rec["status"] == "eligible"]
    rows = []
    for tset, dd in (("primary", el[el["primary"]]), ("all_trios", el)):
        masks = _subsets(dd)
        for sc, d0 in (("classical8", dd[dd["gene"].isin(G8) & ~dd["bio"]]), ("all_fixed_genes", dd[~dd["bio"]]), ("drb345_biology", dd[dd["bio"]])):
            for sname, mk in _subsets(d0).items():
                d = d0[mk]
                n = len(d)
                rows.append(dict(trio_set=tset, scope=sc, origin_subset=sname, gene="ALL", n_eligible=n, n_match=int((d["cls"] == "match").sum()),
                                 n_other_parent_only=int((d["cls"] == "other_parent_only").sum()), n_neither=int((d["cls"] == "neither").sum()),
                                 **{k: v for k, v in rate_row((d["cls"] == "neither").sum(), n).items() if k in ("pct", "lo", "hi")}))
            if sc != "drb345_biology":
                for g, d in d0.groupby("gene"):
                    rows.append(dict(trio_set=tset, scope=sc, origin_subset="default", gene=g, n_eligible=len(d), n_match=int((d["cls"] == "match").sum()),
                                     n_other_parent_only=int((d["cls"] == "other_parent_only").sum()), n_neither=int((d["cls"] == "neither").sum()),
                                     **{k: v for k, v in rate_row((d["cls"] == "neither").sum(), len(d)).items() if k in ("pct", "lo", "hi")}))
    w(pd.DataFrame(rows), "B_class_counts.tsv")
    # name vs sequence: among 2-field matches, how many differ at the 3-field (synonymous) level
    ex = el[~el["bio"]]
    w(pd.DataFrame([dict(scope="all_fixed_genes", n_name_matches=int((ex["cls"] == "match").sum()),
                         n_name_match_cds_differs=int(ex["cds_differs_given_name_match"].sum()),
                         note="2-field name equal but 3-field (synonymous) name differs: same protein, different CDS; not an error at the protein level")]),
        "B_name_vs_cds.tsv")
    # origin robustness: k-mer margin by class (log2 ratio bins; winner k-mers bins)
    mrows = []
    for cl, d in ex.groupby("cls"):
        lr = np.log2(np.maximum(d["ratio"].to_numpy(float), 1.0))
        cnt, edges = np.histogram(np.clip(lr, 0, 8), bins=np.arange(0, 8.01, 0.5))
        for lo, n in zip(edges[:-1], cnt):
            mrows.append(dict(cls=cl, feature="log2_kmer_ratio", bin_lo=float(lo), bin_hi=float(lo + 0.5), n=int(n)))
        cnt, edges = np.histogram(np.clip(d["k_win"].to_numpy(float), 0, 200), bins=np.arange(0, 201, 5))
        for lo, n in zip(edges[:-1], cnt):
            mrows.append(dict(cls=cl, feature="winner_kmers", bin_lo=float(lo), bin_hi=float(lo + 5), n=int(n)))
    w(pd.DataFrame(mrows), "B_origin_margin.tsv")
    # families: events per family (family ids are integers local to this run; per-family rows stay INTERNAL)
    fam = ex.groupby("family").agg(n_records=("cls", "size"), n_neither=("bad", "sum"),
                                   n_other_only=("cls", lambda x: int((x == "other_parent_only").sum())))
    fh = fam["n_neither"].value_counts().sort_index()
    w(pd.DataFrame([dict(families=len(fam), families_with_neither=int((fam["n_neither"] > 0).sum()),
                         families_with_other_only=int((fam["n_other_only"] > 0).sum()),
                         max_neither_in_family=int(fam["n_neither"].max()) if len(fam) else 0, n_neither=int(fam["n_neither"].sum()),
                         effective_events=int((fam["n_neither"] > 0).sum()),
                         records_per_family_median=float(fam["n_records"].median()) if len(fam) else np.nan)]), "B_family_clustering.tsv")
    w(pd.DataFrame({"neither_events_in_family": fh.index, "families": fh.values}), "B_family_event_hist.tsv")
    # homozygous child: a homozygous child call that is "neither" is more likely a genuine call inconsistency
    hom = ex.groupby("child_homozygous").agg(n=("bad", "size"), n_neither=("bad", "sum"), n_other_only=("cls", lambda x: int((x == "other_parent_only").sum()))).reset_index()
    w(hom, "B_by_child_homozygosity.tsv")
    return fam


# ===========================================================================
# combined unit (child person x gene: missing OR wrong) and pair units
# ===========================================================================
def combined_units(rec, pgu):
    """Child person x gene units for the 'any problem' endpoint. missing = the child has one call at the gene; wrong = a `neither` call in a fully
    called child with a fully called transmitting parent. Child complete but a parent call incomplete or novel: excluded (cannot be labelled)."""
    r = rec[~rec["bio"]]
    g = r.groupby(["person", "gene"])
    inc = g["status"].apply(lambda s: bool((s == "child_incomplete").any())).rename("missing")
    wrong = g["bad"].max().rename("wrong")
    anyel = g["status"].apply(lambda s: bool((s == "eligible").any() and not s.isin(["parent_incomplete", "novel_unlabelled"]).any())).rename("eligible")
    fam = g["family"].first()
    prim = g["primary"].first()
    u = pd.concat([inc, wrong, anyel, fam, prim], axis=1).reset_index()
    u = u[u["missing"] | u["eligible"]].copy()
    u["wrong"] = np.where(u["missing"], 0, u["wrong"]).astype(int)
    u["missing"] = u["missing"].astype(int)
    u["bad"] = ((u["missing"] == 1) | (u["wrong"] == 1)).astype(int)
    m = pgu.drop(columns=["bad"], errors="ignore").merge(u[["person", "gene", "missing", "wrong", "bad", "family", "primary"]], on=["person", "gene"])
    return m


def pair_units(pair_rec, genes_exclude=BIO):
    """71's pair records -> two views: all (as in 71) and `complete` (both people fully called at the gene); the gene is restricted to fixed-content genes.
    A family id comes from the connected components of the pair graph."""
    d = pair_rec[~pair_rec["gene"].isin(genes_exclude)].copy()
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i, j in zip(d["person"], d["other"]):
        parent[find(i)] = find(j)
    roots = {r: n for n, r in enumerate(sorted({find(p) for p in d["person"]}))}
    d["family"] = [roots[find(p)] for p in d["person"]]
    d["complete"] = (d["w_n_called"] == -2).astype(int) if "w_n_called" in d else 0
    return d, d[d["complete"] == 1].reset_index(drop=True)


# ===========================================================================
# STREAM C: novel / undetermined
# ===========================================================================
def analyse_novel(ct, slots, twin_rec, out_dir):
    w = lambda df, name: C.write_tsv(df, os.path.join(out_dir, name))  # noqa: E731
    c = ct.copy()
    c["novel_any"] = c["x_is_novel"].fillna(0).astype(int)
    c["novel_unnamed"] = (~c["named"]).astype(int)
    rows = []
    for g, d in list(c.groupby("gene")) + [("ALL", c)]:
        rows.append(dict(gene=g, n_calls=len(d), novel_any=int(d["novel_any"].sum()), novel_unnamed_field12=int(d["novel_unnamed"].sum()),
                         **{"novel_any_" + k: v for k, v in rate_row(d["novel_any"].sum(), len(d)).items() if k in ("pct", "lo", "hi")},
                         **{"unnamed_" + k: v for k, v in rate_row(d["novel_unnamed"].sum(), len(d)).items() if k in ("pct", "lo", "hi")}))
    w(pd.DataFrame(rows), "C_novel_rates_by_gene.tsv")
    # novelty vs depth and vs completeness (fine bins of gene depth)
    b = []
    edges = np.round(np.arange(0, 2.0001, 0.05), 2)
    v = c["x_rel_cds"].to_numpy(float)
    ok = np.isfinite(v)
    ix = np.clip(np.digitize(v[ok], edges) - 1, 0, len(edges) - 2)
    for col in ("novel_any", "novel_unnamed"):
        n = np.bincount(ix, minlength=len(edges) - 1)
        k = np.bincount(ix, weights=c[col].to_numpy(float)[ok], minlength=len(edges) - 1)
        for j in range(len(edges) - 1):
            if n[j]:
                b.append(dict(outcome=col, feature="x_rel_cds", bin_lo=float(edges[j]), bin_hi=float(edges[j + 1]), **rate_row(k[j], n[j])))
    w(pd.DataFrame(b), "C_novel_by_depth.tsv")
    nc = c.groupby(["person_id", "gene"]).agg(n=("hap", "nunique"), nov=("novel_unnamed", "max")).reset_index()
    w(pd.DataFrame([dict(has_unnamed_novel=int(k), person_genes=len(d), incomplete=int((d["n"] < 2).sum()),
                         **{x: v for x, v in rate_row((d["n"] < 2).sum(), len(d)).items() if x in ("pct", "lo", "hi")})
                    for k, d in nc.groupby("nov")]), "C_novel_vs_incomplete.tsv")
    # twin / duplicate concordance of the full consensus string for person-genes where one person has a novel call
    if twin_rec is not None and len(twin_rec):
        pairs = set(zip(twin_rec["person"], twin_rec["other"]))
        cons = {k: set(d["consensus"].dropna()) for k, d in c.groupby(["person_id", "gene"])}
        nov = {k: bool(d["novel_any"].max()) for k, d in c.groupby(["person_id", "gene"])}
        r = []
        for i, j in pairs:
            for g in set(c.loc[c["person_id"] == i, "gene"]) & set(c.loc[c["person_id"] == j, "gene"]):
                if (i, g) in cons and (j, g) in cons:
                    r.append((nov[(i, g)] or nov[(j, g)], cons[(i, g)] == cons[(j, g)]))
        d = pd.DataFrame(r, columns=["novel", "same"])
        w(pd.DataFrame([dict(novel_in_pair=bool(k), person_genes=len(s), identical_consensus=int(s["same"].sum()),
                             **{x: v for x, v in rate_row(s["same"].sum(), len(s)).items() if x in ("pct", "lo", "hi")}) for k, s in d.groupby("novel")]),
          "C_twin_concordance.tsv")


# ===========================================================================
def definition_card(out_dir, slots, rec, units_a, pairs_all, pairs_complete, meta):
    r = rec[~rec["bio"]]
    rows = [
        ("A", "slot", "person x gene x haplotype; 2 per person-gene", len(slots)),
        ("A", "missing slot", "no gene annotation on that haplotype (Immuannot found no template passing its filter, or no sequence there)", int((slots["status"] == "missing").sum())),
        ("A", "incomplete person-gene unit", "a person-gene with a call on only one haplotype (units_A; 71e unit)", int(units_a["incomplete"].sum())),
        ("B", "trio child call with resolved origin", "yak origin pat/mat (>=5 k-mers and >=3x the other parent, gene +-10 kb)", len(rec)),
        ("B", "eligible", "child fully called, transmitting parent fully called, child call named: labelled", int((r["status"] == "eligible").sum())),
        ("B", "child_incomplete", "child has one call at the gene: counted in stream A, not labelled", int((r["status"] == "child_incomplete").sum())),
        ("B", "parent_incomplete", "transmitting parent has <2 calls: absent parent allele cannot be told from a typing error", int((r["status"] == "parent_incomplete").sum())),
        ("B", "novel_unlabelled", "child or parent call novel in field 1-2: stream C", int((r["status"] == "novel_unlabelled").sum())),
        ("B", "match", "child's 2-field name is in the transmitting parent's calls", int((r["cls"] == "match").sum())),
        ("B", "other_parent_only", "name absent in the transmitting parent but present in the other parent: yak origin error / recombination, not a typing error", int((r["cls"] == "other_parent_only").sum())),
        ("B", "neither", "name absent from both parents: typing-error candidate (the error label of 76b)", int((r["cls"] == "neither").sum())),
        ("B", "primary trio", "child with exactly one candidate trio and parent-child IBD0 not failing", int(meta["primary"].sum())),
        ("B", "trio verified", "both child-parent IBD0 <= %.2f" % IBD0_MAX, int((meta["verified"] == "yes").sum())),
        ("B", "trio verification failed", "a child-parent IBD0 above the limit", int((meta["verified"] == "no").sum())),
        ("B", "trio unverifiable", "no IBD0/IBS0 column in the relatedness table", int((meta["verified"] == "unverifiable").sum())),
        ("pair", "pair record (all)", "first-degree pair x gene (71, named only); bad = no shared 2-field allele", len(pairs_all)),
        ("pair", "pair record (both fully called)", "both people have two calls at the gene", len(pairs_complete)),
    ]
    C.write_tsv(pd.DataFrame(rows, columns=["stream", "term", "definition", "n"]), os.path.join(out_dir, "definition_card.tsv"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--internal-71", default=INT71)
    ap.add_argument("--origin", default=Q.ORIGIN57B)
    ap.add_argument("--out-dir", default=OUT)
    ap.add_argument("--internal-dir", default=INT)
    ap.add_argument("--ibd0-max", type=float, default=IBD0_MAX)
    ap.add_argument("--boot", type=int, default=200)
    ap.add_argument("--limit-people", type=int, default=None, help="debug")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    os.makedirs(a.internal_dir, exist_ok=True)
    T = _load("t57", "57_yak_trioeval_pilot.py")
    B = _load("b56", "56_phasing_typing_benchmark.py")
    ct = Q.read_int(os.path.join(a.internal_71, "call_table.INTERNAL.tsv"), dtype={"person_id": str, "hap": str})
    pair_rec = Q.read_int(os.path.join(a.internal_71, "records_pair.INTERNAL.tsv"), dtype={"person": str, "other": str})
    twin_rec = Q.read_int(os.path.join(a.internal_71, "records_twin.INTERNAL.tsv"), dtype={"person": str, "other": str})
    C.TABLE1_USECOLS = C.TABLE1_USECOLS + [c for c in Q.EXTRA_T1 if c not in C.TABLE1_USECOLS]
    t1 = C.load_table1(C.DEFAULT_TABLE1, a.limit_people)
    log("Table 1 %d rows, call table %d" % (len(t1), len(ct)))
    # ---- stream A
    genes = sorted(set(t1["gene"]) & set(G8 + BIO + ANCHORS))
    slots = attach_slot_features(slot_table(t1, genes), ct)
    C.write_tsv(slots, os.path.join(a.internal_dir, "slots.INTERNAL.tsv"))
    analyse_slots(slots, a.out_dir, boot=a.boot)
    units_a = person_gene_units(ct, genes)
    units_a.to_csv(os.path.join(a.internal_dir, "units_A.INTERNAL.tsv"), sep="\t", index=False)
    # ---- stream B
    cohort = C.load_cohort(C.DEFAULT_COHORT)
    lr = set(cohort.loc[cohort["in_lr"].astype(str).str.lower() == "true", "person_id"])
    trios = list(T.find_trios(lr))
    rel = B.load_rel(C.DEFAULT_RELATEDNESS)
    trio_of, meta = trio_meta(trios, rel, a.ibd0_max)
    origin = Q.read_int(a.origin, dtype={"child": str, "contig": str})
    rec = classify_trio_calls(ct, origin, trio_of, meta, load_positions(t1))
    rec.to_csv(os.path.join(a.internal_dir, "trio_calls.INTERNAL.tsv"), sep="\t", index=False)
    fam = summarise_trio(rec, a.out_dir)
    fam.to_csv(os.path.join(a.internal_dir, "family_events.INTERNAL.tsv"), sep="\t")
    el = rec[(rec["status"] == "eligible") & ~rec["bio"]]
    el.to_csv(os.path.join(a.internal_dir, "units_B_all.INTERNAL.tsv"), sep="\t", index=False)
    el[el["primary"]].to_csv(os.path.join(a.internal_dir, "units_B.INTERNAL.tsv"), sep="\t", index=False)
    comb = combined_units(rec, units_a)
    comb.to_csv(os.path.join(a.internal_dir, "units_combined.INTERNAL.tsv"), sep="\t", index=False)
    pa, pc = pair_units(pair_rec)
    pa.to_csv(os.path.join(a.internal_dir, "units_pair_all.INTERNAL.tsv"), sep="\t", index=False)
    pc.to_csv(os.path.join(a.internal_dir, "units_pair_complete.INTERNAL.tsv"), sep="\t", index=False)
    # ---- stream C and the card
    analyse_novel(ct, slots, twin_rec, a.out_dir)
    definition_card(a.out_dir, slots, rec, units_a, pa, pc, meta)
    log("DONE 76")
    return 0


if __name__ == "__main__":
    sys.exit(main())
