#!/usr/bin/env python3
"""S06 71: one table of per-call quality metrics, the typing-error truth sets, and the depth distributions (VM).

Replaces the hand-built 0-7 score of 56 with measurable quantities. For every HLA call (person x hap x gene,
copy 1) this joins
  Immuannot's own numbers   template_distance (+ per kb), cds_distance, n_aa_changes, template_warning, n_tied, novelty
  read depth over the exons 70's gene / exon-union / CDS / ARS / minimum-exon depth, per-base minimum, 10th percentile,
                            share of CDS bases under 10x, MAPQ20/MAPQ0 ratio; all relative to the person's genome
                            control depth (58), plus the person's MHC-window depth
  assembly context          MHC contigs on the haplotype, gene on the same contig as its phasing partner
  the old 0-7 score         56's `score`, kept only as a comparison.
Every metric is stored twice: raw (`x_<name>`) and oriented (`o_<name>`, sign flipped where LOW is bad, so a higher
`o_` always means a worse call). Records for the error analysis (71b) are built from the oriented columns.

Typing-error truth sets (one row per record; `bad` = 1 is an error):
  trio   57b/yak: the parent that yak says transmitted the child haplotype does not carry the child's 2-field allele
         (`bad`) or full allele (`bad_full`). Metrics of the child call, of the transmitting parent (worst of its calls at
         that gene) and the worse of the two. A mismatch can come from the parent's call, hence the parent columns.
  pair   56b: first-degree pair sharing their MHC (at most one classical gene without a shared 2-field allele); bad = no
         shared 2-field allele (typing error or recombination). Metrics = worse of the two people at that gene.
  twin   56: twin / duplicate pair; bad = unordered 2-field genotypes differ.

Outputs
  INTERNAL (~/s06/internal/71, never leaves the Workbench): call_table, records_trio, records_pair, records_twin,
           person_depth (person-level), priority_ids.txt (people in any truth set, for 70's run order)
  AGGREGATE (~/s06/results/71): fine-grained distribution tables (counts per bin; no person ids)
           dist_gene_depth_hist.tsv      per gene x feature x platform: histogram of relative depth, 0.01 bins
           dist_gene_depth_quantiles.tsv per gene x feature x platform: quantiles of relative and absolute depth
           dist_exon_depth.tsv           per gene x exon of the reference transcript: quantiles of relative depth
           dist_person_mhc_hist.tsv      per-person MHC-window depth, relative (0.01 bins) and absolute (0.5x bins)
           dist_mhc_vs_genome_grid.tsv   2-D counts: genome control depth x MHC-window depth
           corr_gene_depth.tsv           gene x gene Spearman correlation of CDS relative depth across people
           corr_metric_metric.tsv        metric x metric Spearman correlation over all calls
           truth_set_sizes.tsv           records and errors per truth set x gene

VM:
    cd ~/repos/pilot-validation
    pixi run -e spechla -- python3 scripts/hla_popgen/71_call_quality_table.py --emit-priority   # before 70
    pixi run -e spechla -- python3 scripts/hla_popgen/71_call_quality_table.py                   # after 70 (partial is fine)

Tests: scripts/hla_popgen/tests/test_71_call_quality_table.py (synthetic).
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


OUT = os.path.expanduser("~/s06/results/71/")
INT = os.path.expanduser("~/s06/internal/71/")
DEPTH70 = os.path.expanduser("~/s06/internal/70/feature_depth.tsv")
PERBASE70 = os.path.expanduser("~/s06/internal/70/cds_perbase.tsv")
FEATURES70 = os.path.expanduser("~/s06/results/70/features.tsv")
DEPTH58 = os.path.expanduser("~/s05/internal/58/per_person_depth.tsv")
ORIGIN57B = os.path.expanduser("~/s05/internal/57b/gene_origin.INTERNAL.tsv")
EXTRA_T1 = ["n_tied", "cds_distance", "n_aa_changes", "template_distance_per_kb"]

# key, label, family, sign (+1: higher is worse, -1: lower is worse)
METRICS = [
    ("template_distance", "Immuannot template distance (gene span)", "immuannot", +1),
    ("template_distance_per_kb", "Immuannot template distance per kb", "immuannot", +1),
    ("cds_distance", "Immuannot CDS distance", "immuannot", +1),
    ("n_aa_changes", "Amino-acid changes vs closest allele", "immuannot", +1),
    ("has_warning", "Immuannot template warning", "immuannot", +1),
    ("n_tied", "Alleles tied for best CDS match", "immuannot", +1),
    ("is_novel", "Novel allele call", "immuannot", +1),
    ("rel_gene", "Gene depth (relative)", "depth", -1),
    ("rel_exons", "Exon depth (relative)", "depth", -1),
    ("rel_cds", "CDS depth (relative)", "depth", -1),
    ("rel_ars", "Peptide-binding exons depth (relative)", "depth", -1),
    ("rel_exon_min", "Shallowest coding exon depth (relative)", "depth", -1),
    ("rel_cds_n", "CDS depth, normalised to its gene and platform", "depth", -1),
    ("rel_ars_n", "Peptide-binding exon depth, normalised to its gene and platform", "depth", -1),
    ("rel_exon_min_n", "Shallowest exon depth, each exon normalised to its own median", "depth", -1),
    ("cds_p10_rel", "CDS 10th-percentile base depth (relative)", "depth", -1),
    ("cds_min_rel", "CDS minimum base depth (relative)", "depth", -1),
    ("ars_min_rel", "Peptide-binding exons minimum base depth (relative)", "depth", -1),
    ("cds_share_lt10", "Share of CDS bases under 10x", "depth", +1),
    ("ars_share_lt10", "Share of peptide-binding bases under 10x", "depth", +1),
    ("mapq_frac", "Unique-mapping fraction of CDS depth (MAPQ20/MAPQ0)", "depth", -1),
    ("rel_mhc", "MHC-window depth (relative, person)", "depth", -1),
    ("n_called", "Haplotypes with a call at the gene (2 = complete genotype)", "assembly", -1),
    ("mhc_contigs", "MHC contigs on the haplotype", "assembly", +1),
    ("diff_contig", "Gene off its partner's contig", "assembly", +1),
    ("score7", "Old 0-7 confidence score", "legacy", -1),
]
SIGN = {k: s for k, _, _, s in METRICS}
LABEL = {k: l for k, l, _, _ in METRICS}
FAMILY = {k: f for k, _, f, _ in METRICS}
KEYS = [k for k, *_ in METRICS]


def log(msg):
    C.log("[71] " + msg)


def read_int(path, **kw):
    with open(path) as f:
        first = f.readline()
    return pd.read_csv(path, sep="\t", skiprows=1 if first.startswith("#") else 0, **kw)


# ---------------------------------------------------------------------------
# depth -> per person x gene table
# ---------------------------------------------------------------------------
def control_depth(per58):
    """58 per_person_depth -> {(person, mapq): control depth} and {person: MHC window relative depth (mapq 0)}."""
    d = per58.drop_duplicates(["person_id", "mapq"])
    ctrl = {(p, int(q)): float(c) for p, q, c in zip(d["person_id"], d["mapq"], d["control_depth"])}
    w = per58[(per58["region"] == "window_MHC") & (per58["mapq"] == 0)]
    return ctrl, dict(zip(w["person_id"], w["rel_depth"])), dict(zip(w["person_id"], w["mean_depth"]))


def gene_depth_table(fd, features, ctrl, perbase, plat_of=None):
    """fd: 70 feature_depth (person_id, feature, mapq, sum_depth, mean_depth); features: 70 features.tsv;
    perbase: 70 cds_perbase. -> DataFrame person_id, gene with absolute and relative depth columns."""
    f = features.set_index("name")
    fd = fd.copy()
    fd["gene"] = fd["feature"].map(f["gene"])
    fd["kind"] = fd["feature"].map(f["kind"])
    fd["ars"] = fd["feature"].map(f["ars"]).astype(bool)
    fd["length"] = fd["feature"].map(f["length"]).astype(float)
    fd = fd[fd["gene"].notna()]

    def wmean(sub, label):
        g = sub.groupby(["person_id", "gene", "mapq"]).agg(s=("sum_depth", "sum"), L=("length", "sum")).reset_index()
        g[label] = g["s"] / g["L"]
        return g[["person_id", "gene", "mapq", label]]

    parts = {
        "gene": wmean(fd[fd["kind"] == "gene"], "abs_gene"),
        "exons": wmean(fd[fd["kind"] == "exons"], "abs_exons"),
        "cds": wmean(fd[fd["kind"] == "cds"], "abs_cds"),
        "ars": wmean(fd[(fd["kind"] == "ex") & fd["ars"]], "abs_ars"),
    }
    ex = fd[fd["kind"] == "ex"]
    exmin = ex.groupby(["person_id", "gene", "mapq"])["mean_depth"].min().reset_index().rename(
        columns={"mean_depth": "abs_exon_min"})
    out = None
    for tab in list(parts.values()) + [exmin]:
        out = tab if out is None else out.merge(tab, on=["person_id", "gene", "mapq"], how="outer")
    q0 = out[out["mapq"] == 0].drop(columns="mapq")
    q20 = out[out["mapq"] == 20][["person_id", "gene", "abs_cds"]].rename(columns={"abs_cds": "abs_cds_q20"})
    t = q0.merge(q20, on=["person_id", "gene"], how="left")
    c0 = t["person_id"].map(lambda p: ctrl.get((p, 0), np.nan)).to_numpy(float)
    for k in ("gene", "exons", "cds", "ars", "exon_min"):
        t["rel_" + k] = t["abs_" + k] / c0
    t["mapq_frac"] = np.where(t["abs_cds"] > 0, t["abs_cds_q20"] / t["abs_cds"], np.nan)
    pb = perbase[perbase["mapq"] == 0].copy()
    pc = pb["person_id"].map(lambda p: ctrl.get((p, 0), np.nan)).to_numpy(float)
    for src, dst in (("cds_p10", "cds_p10_rel"), ("cds_min", "cds_min_rel"), ("ars_min", "ars_min_rel")):
        pb[dst] = pb[src] / pc
    pb = pb[["person_id", "gene", "cds_p10_rel", "cds_min_rel", "ars_min_rel", "cds_share_lt10", "ars_share_lt10",
             "cds_median", "cds_mean"]]
    t = t.merge(pb, on=["person_id", "gene"], how="left")
    t["control_depth"] = c0
    # Normalised versions: divide by the cohort median of the same gene (and exon) within a platform, so that a gene or an exon
    # that reference-aligns systematically low (HLA-DRB1 exon 2: allelic divergence from GRCh38) is not read as "thin".
    plat = t["person_id"].map(plat_of) if plat_of else pd.Series("all", index=t.index)
    for src, dst in (("rel_cds", "rel_cds_n"), ("rel_ars", "rel_ars_n")):
        med = t.groupby([t["gene"], plat])[src].transform("median")
        t[dst] = t[src] / med
    exd = fd[(fd["kind"] == "ex") & (fd["mapq"] == 0)].copy()
    exd["rel"] = exd["mean_depth"] / exd["person_id"].map(lambda p: ctrl.get((p, 0), np.nan))
    exd["plat"] = exd["person_id"].map(plat_of) if plat_of else "all"
    exd["reln"] = exd["rel"] / exd.groupby(["feature", "plat"])["rel"].transform("median")
    mn = exd.groupby(["person_id", "gene"])["reln"].min().rename("rel_exon_min_n").reset_index()
    t = t.merge(mn, on=["person_id", "gene"], how="left")
    return t


# ---------------------------------------------------------------------------
# call table
# ---------------------------------------------------------------------------
def build_call_table(calls, gdepth, rel_mhc):
    """calls: 56 `prepare` output (copy_index kept); gdepth: gene_depth_table; rel_mhc: {person: rel}.
    -> one row per copy-1 call with x_<metric> (raw) and o_<metric> (oriented, higher = worse)."""
    t = calls[calls["copy_index"] == 1].copy()
    t["hap"] = t["hap"].astype(str)
    td = pd.to_numeric(t["template_distance"], errors="coerce")
    glen = (pd.to_numeric(t.get("gene_end"), errors="coerce") - pd.to_numeric(t.get("gene_start"), errors="coerce") + 1)
    t["x_template_distance"] = td
    t["x_template_distance_per_kb"] = (pd.to_numeric(t["template_distance_per_kb"], errors="coerce")
                                       if "template_distance_per_kb" in t else 1000 * td / glen)
    # cds_distance / n_aa_changes are absent when Immuannot found a perfect template: that is "no difference", i.e. 0
    t["x_cds_distance"] = pd.to_numeric(t["cds_distance"], errors="coerce").fillna(0) if "cds_distance" in t else np.nan
    if "n_aa_changes" in t:
        t["x_n_aa_changes"] = pd.to_numeric(t["n_aa_changes"], errors="coerce").fillna(0)
    else:
        t["x_n_aa_changes"] = t["cds_mut"].map(_count_aa)
    w = t["template_warning"].astype(str)
    t["x_has_warning"] = (t["template_warning"].notna() & ~w.isin(["", "nan", "None", "NA"])).astype(int)
    t["x_n_tied"] = pd.to_numeric(t["n_tied"], errors="coerce") if "n_tied" in t else np.nan
    t["x_is_novel"] = t["novelty_depth"].notna().astype(int)
    t = t.merge(gdepth.drop(columns=["control_depth"], errors="ignore"), on=["person_id", "gene"], how="left")
    for k in ("rel_gene", "rel_exons", "rel_cds", "rel_ars", "rel_exon_min", "rel_cds_n", "rel_ars_n", "rel_exon_min_n",
              "cds_p10_rel", "cds_min_rel",
              "ars_min_rel", "cds_share_lt10", "ars_share_lt10", "mapq_frac"):
        t["x_" + k] = t[k] if k in t else np.nan
    t["x_rel_mhc"] = t["person_id"].map(rel_mhc)
    t["x_n_called"] = t.groupby(["person_id", "gene"])["hap"].transform("nunique")
    t["x_mhc_contigs"] = t["nctg"]
    t["x_diff_contig"] = _off_partner_contig(t)
    t["x_score7"] = t["score"]
    for k in KEYS:
        t["o_" + k] = SIGN[k] * pd.to_numeric(t["x_" + k], errors="coerce")
    t["named"] = t["consensus"].map(lambda c: C.two_field(c) is not None)
    keep = ["person_id", "hap", "contig", "gene", "consensus", "tf", "plat", "named"] + ["x_" + k for k in KEYS] \
           + ["o_" + k for k in KEYS]
    return t[keep].reset_index(drop=True)


def _count_aa(cds_mut):
    if not isinstance(cds_mut, str) or not cds_mut:
        return 0
    parts = cds_mut.split(",")[0].split("|")
    return len(parts[2].split(":")) if len(parts) > 2 and parts[2] else 0


def _off_partner_contig(t):
    """1 when the gene's phasing partner (56's PART) is called on the same haplotype but on another contig, 0 when on
    the same contig, NaN when the partner is not called."""
    B = sys.modules.get("b56") or _load("b56", "56_phasing_typing_benchmark.py")
    where = {(p, h, g): c for p, h, g, c in zip(t["person_id"], t["hap"], t["gene"], t["contig"])}
    out = []
    for p, h, g, c in zip(t["person_id"], t["hap"], t["gene"], t["contig"]):
        vals = [where.get((p, h, q)) for q in B.PART.get(g, [])]
        vals = [v for v in vals if v is not None]
        out.append(np.nan if not vals else float(any(v != c for v in vals)))
    return out


def person_gene_worst(call_table):
    """Per (person, gene): the worse (max oriented) of the person's calls, for pair / twin / parent records."""
    ocols = ["o_" + k for k in KEYS]
    return call_table.groupby(["person_id", "gene"])[ocols].max()


# ---------------------------------------------------------------------------
# truth records
# ---------------------------------------------------------------------------
def _fill_worst(rec, own, other, prefix_other="p"):
    for k in KEYS:
        o = "o_" + k
        rec["c_" + k] = own[o].to_numpy(float) if own is not None else np.nan
        rec[prefix_other + "_" + k] = other[o].to_numpy(float)
        rec["w_" + k] = np.fmax(rec["c_" + k], rec[prefix_other + "_" + k])


def trio_records(origin, call_table, trio_of, named_only=False):
    """origin: 57b per-call table (child, hap, gene, origin); trio_of: child -> (parent_a, parent_b) with a the 'pat' side
    as 57b defined it."""
    ct = call_table.set_index(["person_id", "hap", "gene"])
    pg = person_gene_worst(call_table)
    pidx = {k: i for i, k in enumerate(pg.index)}
    parr = pg.to_numpy(float)
    ocols = ["o_" + k for k in KEYS]
    carr = ct[ocols].to_numpy(float)
    cidx = {k: i for i, k in enumerate(ct.index)}
    tfs, cons, plats = ct["tf"].to_numpy(object), ct["consensus"].to_numpy(object), ct["plat"].to_numpy(object)
    allele = {(p, g): (set(d["tf"].dropna()), set(d["consensus"].dropna())) for (p, g), d in call_table.groupby(["person_id", "gene"])}
    # named_only: a call whose 2-field name is undefined (novel in field 1 or 2, "A*01:new") cannot be compared with a parent's collapsed
    # name; without this it is counted as a mismatch by construction. Skip a record when the child's call, or ANY call of the parent at
    # that gene, is such a call.
    unnamed_pg = set(zip(call_table.loc[~call_table["named"], "person_id"], call_table.loc[~call_table["named"], "gene"])) if named_only else set()
    cnamed = ct["named"].to_numpy() if named_only else None
    rows = []
    for c, h, g, og in zip(origin["child"], origin["hap"], origin["gene"], origin["origin"]):
        if og not in ("pat", "mat") or c not in trio_of or (c, h, g) not in cidx:
            continue
        par = trio_of[c][0] if og == "pat" else trio_of[c][1]
        if (par, g) not in allele or (par, g) not in pidx:
            continue
        i = cidx[(c, h, g)]
        if not isinstance(tfs[i], str):
            continue
        if named_only and (not cnamed[i] or (par, g) in unnamed_pg):
            continue
        ptf, pcons = allele[(par, g)]
        row = dict(person=c, other=par, gene=g, hap=h, origin=og, plat=plats[i], bad=int(tfs[i] not in ptf),
                   bad_full=int(cons[i] not in pcons))
        pv = parr[pidx[(par, g)]]
        for j, k in enumerate(KEYS):
            row["c_" + k], row["p_" + k], row["w_" + k] = carr[i, j], pv[j], np.fmax(carr[i, j], pv[j])
        rows.append(row)
    return pd.DataFrame(rows)


def pair_records(pairs, geno, call_table, genes, kind):
    """pairs: DataFrame i, j, pp (and cls); geno: 56 index_calls dict; kind 'pair' or 'twin'. For pairs the 56b
    filter applies: at most one classical gene without a shared 2-field allele."""
    pg = person_gene_worst(call_table)
    pidx = {k: i for i, k in enumerate(pg.index)}
    parr = pg.to_numpy(float)
    rows = []
    for i, j in zip(pairs["i"], pairs["j"]):
        if kind == "pair":
            fails = [g for g in genes["G8"] if geno.get((i, g)) and geno.get((j, g))
                     and not set(geno[(i, g)][0]) & set(geno[(j, g)][0])]
            if len(fails) > 1:
                continue
        for g in genes["ALL"]:
            a, b = geno.get((i, g)), geno.get((j, g))
            if not a or not b or (i, g) not in pidx or (j, g) not in pidx:
                continue
            bad = int(not set(a[0]) & set(b[0])) if kind == "pair" else int(a[0] != b[0])
            row = dict(person=i, other=j, gene=g, bad=bad, bad_full=int(a[1] != b[1]))
            vi, vj = parr[pidx[(i, g)]], parr[pidx[(j, g)]]
            for n, k in enumerate(KEYS):
                row["c_" + k], row["p_" + k], row["w_" + k] = vi[n], vj[n], np.fmax(vi[n], vj[n])
            rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# distributions (aggregate)
# ---------------------------------------------------------------------------
REL_BINS = np.round(np.arange(0, 2.0001, 0.01), 2)
DEPTH_FEATURES = ["rel_gene", "rel_exons", "rel_cds", "rel_ars", "rel_exon_min", "cds_p10_rel"]


def hist_table(gd, plat_of):
    d = gd.copy()
    d["platform"] = d["person_id"].map(plat_of)
    rows = []
    for (gene, plat), g in list(d.groupby(["gene", "platform"])) + [((gn, "all"), gg) for gn, gg in d.groupby("gene")]:
        for f in DEPTH_FEATURES:
            v = g[f].dropna().to_numpy(float)
            if not len(v):
                continue
            cnt, _ = np.histogram(np.clip(v, 0, 2 - 1e-9), bins=REL_BINS)
            for lo, n in zip(REL_BINS[:-1], cnt):
                rows.append(dict(gene=gene, platform=plat, feature=f, bin_lo=lo, bin_hi=round(lo + 0.01, 2), n=int(n)))
    return pd.DataFrame(rows)


def quantile_table(gd, plat_of):
    d = gd.copy()
    d["platform"] = d["person_id"].map(plat_of)
    qs = [1, 5, 10, 25, 50, 75, 90, 95, 99]
    rows = []
    for (gene, plat), g in list(d.groupby(["gene", "platform"])) + [((gn, "all"), gg) for gn, gg in d.groupby("gene")]:
        for f in DEPTH_FEATURES + ["abs_gene", "abs_exons", "abs_cds", "abs_ars", "cds_share_lt10", "mapq_frac"]:
            v = g[f].dropna().to_numpy(float) if f in g else np.array([])
            if len(v):
                rows.append(dict(gene=gene, platform=plat, feature=f, n_people=len(v), mean=v.mean(),
                                 **{f"p{q:02d}": float(np.percentile(v, q)) for q in qs}))
    return pd.DataFrame(rows)


def exon_table(fd, features, ctrl):
    ex = features[features["kind"] == "ex"].set_index("name")
    d = fd[(fd["mapq"] == 0) & fd["feature"].isin(ex.index)].copy()
    d["rel"] = d["mean_depth"] / d["person_id"].map(lambda p: ctrl.get((p, 0), np.nan))
    d["gene"] = d["feature"].map(ex["gene"])
    d["exon_number"] = d["feature"].map(ex["exon_number"])
    d["cds_ord"] = d["feature"].map(ex["cds_ord"])
    d["ars"] = d["feature"].map(ex["ars"])
    d["length"] = d["feature"].map(ex["length"])
    rows = []
    for (g, n, co, ars, L), s in d.groupby(["gene", "exon_number", "cds_ord", "ars", "length"]):
        v = s["rel"].dropna().to_numpy(float)
        rows.append(dict(gene=g, exon_number=n, cds_ord=co, ars=ars, length=L, n_people=len(v),
                         **{f"p{q:02d}": float(np.percentile(v, q)) for q in (5, 10, 25, 50, 75, 90, 95)}))
    return pd.DataFrame(rows)


def person_mhc_tables(per58, plat_of, anc_of):
    ctrl = per58[per58["mapq"] == 0].drop_duplicates("person_id").set_index("person_id")["control_depth"]
    w = per58[(per58["region"] == "window_MHC") & (per58["mapq"] == 0)].set_index("person_id")
    d = pd.DataFrame({"rel": w["rel_depth"], "abs": w["mean_depth"], "genome": ctrl.reindex(w.index)})
    d["platform"] = d.index.map(plat_of)
    d["ancestry"] = d.index.map(anc_of)
    rows = []
    for (kind, edges) in (("relative", np.round(np.arange(0, 2.0001, 0.01), 2)), ("absolute", np.arange(0, 60.01, 0.5))):
        col = "rel" if kind == "relative" else "abs"
        for key, sub in [(("all", "all"), d)] + [((c, v), g) for c in ("platform", "ancestry") for v, g in d.groupby(c)]:
            cnt, _ = np.histogram(np.clip(sub[col].dropna(), edges[0], edges[-1] - 1e-9), bins=edges)
            for lo, hi, n in zip(edges[:-1], edges[1:], cnt):
                rows.append(dict(scale=kind, stratum_type=key[0], stratum=key[1], bin_lo=lo, bin_hi=hi, n=int(n)))
    ex = np.arange(0, 60.01, 1.0)
    ey = np.arange(0, 60.01, 1.0)
    h, _, _ = np.histogram2d(d["genome"].dropna().clip(0, 59.99), d.loc[d["genome"].notna(), "abs"].clip(0, 59.99),
                             bins=[ex, ey])
    grid = [dict(genome_lo=ex[i], mhc_lo=ey[j], n=int(h[i, j])) for i in range(len(ex) - 1)
            for j in range(len(ey) - 1) if h[i, j] > 0]
    return pd.DataFrame(rows), pd.DataFrame(grid)


def gene_corr(gd, col="rel_cds", min_people=200):
    w = gd.pivot_table(index="person_id", columns="gene", values=col)
    w = w.loc[:, w.notna().sum() >= min_people]
    r = w.corr(method="spearman")
    r.index.name, r.columns.name = "gene_a", "gene_b"
    return r.stack().rename("spearman").reset_index()


def metric_corr(call_table, cap=300_000, seed=71):
    t = call_table
    if len(t) > cap:
        t = t.sample(cap, random_state=seed)
    cols = ["o_" + k for k in KEYS if t["o_" + k].notna().sum() > 100]
    r = t[cols].corr(method="spearman")
    out = r.stack().reset_index()
    out.columns = ["metric_a", "metric_b", "spearman"]
    out["metric_a"], out["metric_b"] = out["metric_a"].str[2:], out["metric_b"].str[2:]
    out["n_calls"] = len(t)
    return out


def truth_sizes(recs):
    rows = []
    for name, df in recs.items():
        if df is None or not len(df):
            continue
        for g, d in [("ALL", df)] + list(df.groupby("gene")):
            rows.append(dict(truth=name, gene=g, n_records=len(d), n_errors=int(d["bad"].sum()),
                             error_pct=round(100 * d["bad"].mean(), 3)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
def priority_ids(pairs, trios):
    ids = set(pairs["i"]) | set(pairs["j"])
    for t in trios:
        ids |= set(t)
    return sorted(ids)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--emit-priority", action="store_true", help="only write the list of truth-set people for 70")
    ap.add_argument("--named-only", action="store_true",
                    help="error labels only from calls with a collapsed (2-field) allele name: records where a child, parent or relative has a novel "
                         "call in field 1-2 at that gene are dropped. Defaults the output dirs to ~/s06/{internal,results}/71_named/")
    ap.add_argument("--out-dir", default=OUT)
    ap.add_argument("--internal-dir", default=INT)
    a = ap.parse_args(argv)
    if a.named_only:
        if a.out_dir == OUT:
            a.out_dir = os.path.expanduser("~/s06/results/71_named/")
        if a.internal_dir == INT:
            a.internal_dir = os.path.expanduser("~/s06/internal/71_named/")
    os.makedirs(a.out_dir, exist_ok=True)
    os.makedirs(a.internal_dir, exist_ok=True)
    B = _load("b56", "56_phasing_typing_benchmark.py")
    T = _load("t57", "57_yak_trioeval_pilot.py")
    cohort = C.load_cohort(C.DEFAULT_COHORT)
    plat = dict(zip(cohort["person_id"], cohort["platform"]))
    anc = dict(zip(cohort["person_id"], cohort["ancestry_pred"]))
    lr = set(cohort.loc[cohort["in_lr"].astype(str).str.lower() == "true", "person_id"])
    C.TABLE1_USECOLS = C.TABLE1_USECOLS + [c for c in EXTRA_T1 if c not in C.TABLE1_USECOLS]
    t1 = C.load_table1(C.DEFAULT_TABLE1)
    rel = B.load_rel(C.DEFAULT_RELATEDNESS)
    pairs = B.build_pairs(rel, set(t1["person_id"]), plat)
    trios = [t for t in T.find_trios(lr)]
    if a.emit_priority:
        ids = priority_ids(pairs, trios)
        with open(os.path.join(a.internal_dir, "priority_ids.txt"), "w") as f:
            f.write("\n".join(ids) + "\n")
        log(f"priority ids: {len(ids)}")
        return 0

    log("loading depth tables")
    per58 = read_int(DEPTH58, dtype={"person_id": str})
    ctrl, rel_mhc, abs_mhc = control_depth(per58)
    features = pd.read_csv(FEATURES70, sep="\t")
    fd = read_int(DEPTH70, dtype={"person_id": str})
    pb = read_int(PERBASE70, dtype={"person_id": str})
    gd = gene_depth_table(fd, features, ctrl, pb, plat)
    log(f"depth: {gd['person_id'].nunique()} people x {gd['gene'].nunique()} genes")

    extra = [g for g in ("DRA",) if g not in B.ALLG]
    B.ALLG.extend(extra)        # in-process only: 56's file and its published numbers are untouched
    calls = B.prepare(t1, plat)
    ct = build_call_table(calls, gd, rel_mhc)
    ct.to_csv(os.path.join(a.internal_dir, "call_table.INTERNAL.tsv"), sep="\t", index=False)
    log(f"call table {len(ct)} calls, {ct['person_id'].nunique()} people")

    genes = {"G8": B.G8, "ALL": B.ALLG}
    if a.named_only:
        # a person's gene with ANY call whose 2-field name is undefined (novel in field 1 or 2) is removed whole, so relatives are never
        # compared on half a genotype
        cn = calls["consensus"].map(lambda c: C.two_field(c) is not None)
        bad_pg = set(zip(calls.loc[~cn, "person_id"], calls.loc[~cn, "gene"]))
        drop = pd.Series([(p, g) in bad_pg for p, g in zip(calls["person_id"], calls["gene"])], index=calls.index)
        log(f"named-only: dropping {int(drop.sum())} call rows of {len(calls)} ({len(bad_pg)} person-genes with a novel call in field 1-2)")
        calls = calls[~drop]
    fd1 = pairs[pairs["cls"] != "twin"]
    tw = pairs[pairs["cls"] == "twin"]
    geno, _ = B.index_calls(calls[calls["person_id"].isin(set(fd1["i"]) | set(fd1["j"]) | set(tw["i"]) | set(tw["j"]))])
    rec_pair = pair_records(fd1, geno, ct, genes, "pair")
    rec_twin = pair_records(tw, geno, ct, genes, "twin")
    rec_trio = pd.DataFrame()
    if os.path.exists(ORIGIN57B):
        o = read_int(ORIGIN57B, dtype={"child": str, "contig": str})
        kids = set(o["child"])
        trio_of = {t[0]: (t[1], t[2]) for t in trios if t[0] in kids}
        rec_trio = trio_records(o, ct, trio_of, named_only=a.named_only)
    for name, df in (("trio", rec_trio), ("pair", rec_pair), ("twin", rec_twin)):
        df.to_csv(os.path.join(a.internal_dir, f"records_{name}.INTERNAL.tsv"), sep="\t", index=False)
        log(f"records {name}: {len(df)}, errors {int(df['bad'].sum()) if len(df) else 0}")

    w = lambda df, name: C.write_tsv(df, os.path.join(a.out_dir, name))  # noqa: E731
    w(truth_sizes({"trio": rec_trio, "pair": rec_pair, "twin": rec_twin}), "truth_set_sizes.tsv")
    w(hist_table(gd, plat), "dist_gene_depth_hist.tsv")
    w(quantile_table(gd, plat), "dist_gene_depth_quantiles.tsv")
    w(exon_table(fd, features, ctrl), "dist_exon_depth.tsv")
    hist, grid = person_mhc_tables(per58, plat, anc)
    w(hist, "dist_person_mhc_hist.tsv")
    w(grid, "dist_mhc_vs_genome_grid.tsv")
    w(gene_corr(gd), "corr_gene_depth.tsv")
    w(metric_corr(ct), "corr_metric_metric.tsv")
    log("DONE 71")
    return 0


if __name__ == "__main__":
    sys.exit(main())
