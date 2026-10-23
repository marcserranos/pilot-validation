#!/usr/bin/env python3
"""S06 78b: pruning evidence for the KIR genes (VM; the KIR counterpart of 76b, which supplies all curve / CV / transfer machinery).

Outcomes (a view = a table of units with `bad`, a family `group`, and oriented o_<metric> columns, higher = worse):
  B_child, B_child_strong, B_child_anytx, B_worst, B_worst_strong
        trio child calls (77c records, gene level, stream 'named'); bad = class `neither` at the protein (3-digit) field
        (the transmitting parent has the locus and neither parent carries the child's label). `only_other_parent` / `tx_absent_other_has`
        are yak-origin or dropout cases and are NOT errors; `absent_in_parents` is biology. B_child / B_worst require the transmitting
        parent to be complete (both haplotypes `complete`, 51's rule) and to carry >=1 call at the locus (informative); *_anytx drops the
        completeness condition; *_strong keeps only origins with >=10 k-mers and >=5x. Metrics of B_child = the child's call; of B_worst =
        the worse of the child's call and the transmitting parent's worst call at the locus.
  H_hap  person x haplotype with calls; bad = the haplotype is not `complete` (multi_contig or one_contig_open_end; `empty` excluded).
        Metrics are haplotype aggregates of the Immuannot fields and of the depth of the haplotype's genes. Contig count, contig-end distance and
        gene count are NOT used here (they define the outcome).
  P_pair, P_pair_complete  first-degree pair x locus (77c), only loci where both people have named calls; bad = no shared label; metrics =
        worse of the two people's worst call at the locus; _complete keeps pairs where both people are complete. Family groups are connected
        components over pair and trio links.
Metrics: Immuannot (template_distance, cds_distance, n_aa_changes, has_warning, n_tied, tier_rank), assembly (n_called, hap_incomplete,
n_contigs_hap, end_proximity = log10(1 + distance to nearer contig end)) and depth from 78 (rel_gene, rel_cds, cds_p10_rel, cds_min_rel,
cds_share_lt10, mapq_frac, dev_cds_q0). Depth exists ONLY for the 10 genes on GRCh38 chr19; for KIR2DL2, 2DL5A/B, 2DS1/2/3/5, 3DS1 it is NaN
(never 0) and such rows simply do not enter that metric's curve. mapq_frac is the quantitative mismapping indicator across the segmental
duplications.

Strata (per_stratum_auc.tsv): gene, paralog group (2DL2/3, 3DL1/S1, 2DL5), gene class (framework / B-content / other), haplotype content
(A: no B-content gene on the haplotype; B: any). Only strata with >=3 events and >=50 units are scored.
Curves, bootstrap bands, family-grouped CV, transfer (trio <-> pair) and power are 76b's (`run`), applied unchanged. Power: with a few tens of events
at most, only AUC >~0.65-0.7 is detectable; power.tsv and POWER_STATEMENT.txt say so explicitly. No "optimal x" is declared.
Small-cell masking: counts 1-19 in the PUBLIC aggregates (this script's out-dir) become "<20", as in 77c; unmasked curves stay in
<internal-dir>/unmasked/ (INTERNAL, never pulled or committed).

    cd ~/repos/pilot-validation-main && pixi run -e spechla -- python3 ~/repos/pv-s06/scripts/hla_popgen/78b_kir_pruning_evidence.py
    ... --boot 50 --reps 5      # quick pass
"""
import argparse
import importlib.util
import os
import sys

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


B76 = sys.modules.get("m76b") or _load("m76b", "76b_pruning_evidence.py")
S = B76.S
OUT = os.path.expanduser("~/s06/results/78b/")
INT = os.path.expanduser("~/s06/internal/78b/")
I77 = os.path.expanduser("~/s06/internal/77/")
I77C = os.path.expanduser("~/s06/internal/77c/")
I78 = os.path.expanduser("~/s06/internal/78/")
GROUP = {"KIR2DL2": "KIR2DL2/3", "KIR2DL3": "KIR2DL2/3", "KIR3DL1": "KIR3DL1/S1", "KIR3DS1": "KIR3DL1/S1",
         "KIR2DL5A": "KIR2DL5", "KIR2DL5B": "KIR2DL5"}
FRAMEWORK = {"KIR3DL3", "KIR2DL4", "KIR3DP1", "KIR3DL2"}
B_CONTENT = {"KIR2DL2", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS5", "KIR3DS1"}
TIER_RANK = {"known": 0, "novel_genomic_known_cds": 1, "novel_cds_synonymous": 2, "novel_protein": 3, "undetermined": 4}
NO_WARN = {"", "0", ".", "false", "nan", "none", "no"}
# key -> sign (+1: higher is worse, -1: lower is worse)
CALL_METRICS = {"template_distance": 1, "cds_distance": 1, "n_aa_changes": 1, "has_warning": 1, "n_tied": 1, "tier_rank": 1,
                "n_called": -1, "hap_incomplete": 1, "n_contigs_hap": 1, "end_proximity": -1}
DEPTH_METRICS = {"rel_gene": -1, "rel_cds": -1, "cds_p10_rel": -1, "cds_min_rel": -1, "cds_share_lt10": 1, "mapq_frac": -1, "dev_cds_q0": 1}
HAP_METRICS = {"h_td_mean": 1, "h_td_max": 1, "h_cds_max": 1, "h_aa_sum": 1, "h_warn_n": 1, "h_novel_share": 1, "h_rel_cds_min": -1,
               "h_rel_cds_mean": -1, "h_mapq_frac_min": -1, "h_share_lt10_max": 1, "h_dev_q0_max": 1}
ALL_SIGN = dict(CALL_METRICS, **DEPTH_METRICS, **HAP_METRICS)
TRANSFERS = [("B_child", "P_pair_complete"), ("P_pair_complete", "B_child"), ("B_child_strong", "P_pair_complete"),
             ("B_worst", "P_pair_complete"), ("P_pair_complete", "B_worst"), ("B_child_anytx", "P_pair_all"), ("P_pair_all", "B_child_anytx")]
COUNT_COLS = ("n_units", "n_errors", "n_pruned", "errors_removed", "test_units", "test_errors", "pruned", "rand_removed_lo", "rand_removed_hi",
              "target_n", "target_n_errors", "target_n_pruned", "target_errors_removed", "events")


def log(m):
    C.log("[78b] " + m)


def read_int(path, **kw):
    with open(path) as f:
        first = f.readline()
    return pd.read_csv(path, sep="\t", skiprows=1 if first.startswith("#") else 0, **kw)


# ---------------------------------------------------------------------------
# call / person-gene metric tables
# ---------------------------------------------------------------------------
def count_aa(cds_mut):
    if not isinstance(cds_mut, str) or not cds_mut or cds_mut.lower() in NO_WARN:
        return 0
    return sum(1 for e in cds_mut.split(",") if e.strip())


def has_warning(v):
    return 0 if (v is None or (isinstance(v, float) and np.isnan(v)) or str(v).strip().lower() in NO_WARN) else 1


def call_table(calls, depth=None):
    """One row per call with x_<key> and oriented o_<key> for the call and depth metrics (depth is person x gene; NaN where undefined)."""
    c = calls.copy()
    c["person_id"] = c["person_id"].astype(str)
    x = pd.DataFrame({"person_id": c["person_id"], "hap": c["hap"], "gene": c["gene"], "copy_index": c["copy_index"]})
    x["x_template_distance"] = pd.to_numeric(c["template_distance"], errors="coerce")
    x["x_cds_distance"] = pd.to_numeric(c["cds_distance"], errors="coerce")
    x["x_n_aa_changes"] = c["cds_mut"].map(count_aa)
    x["x_has_warning"] = c["template_warning"].map(has_warning)
    x["x_n_tied"] = pd.to_numeric(c["n_tied"], errors="coerce")
    x["x_tier_rank"] = c["novelty_tier"].map(TIER_RANK)
    x["x_n_called"] = pd.to_numeric(c["n_called"], errors="coerce")
    x["x_hap_incomplete"] = (c["hap_class"] != "complete").astype(int)
    x["x_n_contigs_hap"] = pd.to_numeric(c["n_contigs_hap"], errors="coerce")
    prox = np.fmin(pd.to_numeric(c["dist_start"], errors="coerce"), pd.to_numeric(c["dist_end"], errors="coerce"))
    x["x_end_proximity"] = np.log10(1 + prox.clip(lower=0))
    x["hap_class"] = c["hap_class"]
    if depth is not None and len(depth):
        d = depth.copy()
        d["person_id"] = d["person_id"].astype(str)
        x = x.merge(d[["person_id", "gene"] + [k for k in DEPTH_METRICS if k in d]].rename(columns={k: "x_" + k for k in DEPTH_METRICS}),
                    on=["person_id", "gene"], how="left")
    for k in DEPTH_METRICS:
        if "x_" + k not in x:
            x["x_" + k] = np.nan
    for k, s in dict(CALL_METRICS, **DEPTH_METRICS).items():
        x["o_" + k] = s * x["x_" + k]
    return x


def person_gene_worst(ct):
    ocols = [c for c in ct if c.startswith("o_")]
    return ct.groupby(["person_id", "gene"])[ocols].max().reset_index()


def hap_content(calls):
    """(person, hap) -> 'B' if any B-content gene is called on the haplotype else 'A'."""
    b = calls[calls["gene"].isin(B_CONTENT)].groupby(["person_id", "hap"]).size()
    keys = calls[["person_id", "hap"]].drop_duplicates()
    return {(p, h): ("B" if (p, h) in b.index else "A") for p, h in zip(keys["person_id"].astype(str), keys["hap"])}


def hap_table(ct, haps, calls):
    """person x hap aggregates. Outcome = not complete (empty haps excluded)."""
    h = haps.copy()
    h["person_id"] = h["person_id"].astype(str)
    h = h[h["hap_class"] != "empty"]
    g = ct.groupby(["person_id", "hap"])
    agg = pd.DataFrame({
        "h_td_mean": g["x_template_distance"].mean(), "h_td_max": g["x_template_distance"].max(), "h_cds_max": g["x_cds_distance"].max(),
        "h_aa_sum": g["x_n_aa_changes"].sum(), "h_warn_n": g["x_has_warning"].sum(),
        "h_novel_share": g["x_tier_rank"].apply(lambda v: float((v >= 3).mean())),
        "h_rel_cds_min": g["x_rel_cds"].min(), "h_rel_cds_mean": g["x_rel_cds"].mean(), "h_mapq_frac_min": g["x_mapq_frac"].min(),
        "h_share_lt10_max": g["x_cds_share_lt10"].max(), "h_dev_q0_max": g["x_dev_cds_q0"].max()}).reset_index()
    out = h[["person_id", "hap", "hap_class"]].merge(agg, on=["person_id", "hap"], how="left")
    out["bad"] = (out["hap_class"] != "complete").astype(int)
    for k, s in HAP_METRICS.items():
        out["o_" + k] = s * out[k]
    cont = hap_content(calls)
    out["content"] = [cont.get((p, hp), "A") for p, hp in zip(out["person_id"], out["hap"])]
    out["group"] = out["person_id"]
    out["person"], out["gene"] = out["person_id"], "ALL"
    return out


# ---------------------------------------------------------------------------
# families: connected components over pair and trio links
# ---------------------------------------------------------------------------
def components(edges):
    parent = {}

    def find(a):
        parent.setdefault(a, a)
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    return {k: find(k) for k in list(parent)}


# ---------------------------------------------------------------------------
# views
# ---------------------------------------------------------------------------
def trio_views(trio, ct, pgw, cont):
    t = trio[(trio["mode"] == "gene") & (trio["stream"] == "named")].copy()
    t["child"], t["tx"], t["other"] = t["child"].astype(str), t["tx"].astype(str), t["other"].astype(str)
    t = t[t["cls3"].isin(["match", "only_other_parent", "neither"])]
    t["bad"] = (t["cls3"] == "neither").astype(int)
    key = ct.drop_duplicates(["person_id", "hap", "gene", "copy_index"]).set_index(["person_id", "hap", "gene", "copy_index"])
    ocols = [c for c in ct if c.startswith("o_")]
    kk = list(zip(t["child"], t["hap"], t["gene"], t["copy_index"]))
    co = key.reindex(kk)[ocols].reset_index(drop=True)
    pg = pgw.set_index(["person_id", "gene"])
    po = pg.reindex(list(zip(t["tx"], t["gene"])))[ocols].reset_index(drop=True)
    t = t.reset_index(drop=True)
    t["group"] = ["|".join(sorted([a, b])) for a, b in zip(t["tx"], t["other"])]
    t["person"] = t["child"]
    t["content"] = [cont.get((c, h), "A") for c, h in zip(t["child"], t["hap"])]
    base = t[["person", "gene", "bad", "group", "content", "strong", "tx_complete"]].copy()
    base["strong"] = base["strong"].astype(int)
    vi = {}
    child = pd.concat([base, co], axis=1)
    worst = pd.concat([base, pd.DataFrame(np.fmax(co.to_numpy(float), po.to_numpy(float)), columns=ocols)], axis=1)
    comp = child["tx_complete"].astype(str).str.lower() == "true"
    vi["B_child"] = child[comp]
    vi["B_child_strong"] = child[comp & (child["strong"] == 1)]
    vi["B_child_anytx"] = child
    wcomp = worst["tx_complete"].astype(str).str.lower() == "true"
    vi["B_worst"] = worst[wcomp]
    vi["B_worst_strong"] = worst[wcomp & (worst["strong"] == 1)]
    return {k: v.reset_index(drop=True) for k, v in vi.items()}


def pair_views(pair, pgw, trio, comp_of):
    p = pair[(pair["mode"] == "gene") & (pair["kind"] == "pair") & pair["cls"].isin(["share", "no_share"])].copy()
    p["person"], p["other"] = p["person"].astype(str), p["other"].astype(str)
    edges = list(zip(p["person"], p["other"]))
    if trio is not None and len(trio):
        edges += list(zip(trio["child"].astype(str), trio["tx"].astype(str))) + list(zip(trio["child"].astype(str), trio["other"].astype(str)))
    comp = components(edges)
    ocols = [c for c in pgw if c.startswith("o_")]
    pg = pgw.set_index(["person_id", "gene"])
    a = pg.reindex(list(zip(p["person"], p["locus"])))[ocols].reset_index(drop=True)
    b = pg.reindex(list(zip(p["other"], p["locus"])))[ocols].reset_index(drop=True)
    p = p.reset_index(drop=True)
    d = pd.DataFrame({"person": p["person"], "gene": p["locus"], "bad": p["bad"].astype(int), "group": [comp.get(x, x) for x in p["person"]],
                      "both_complete": p["both_complete"].astype(str).str.lower() == "true"})
    d = pd.concat([d, pd.DataFrame(np.fmax(a.to_numpy(float), b.to_numpy(float)), columns=ocols)], axis=1)
    d["content"] = "NA"
    return {"P_pair_all": d.drop(columns="both_complete"), "P_pair_complete": d[d["both_complete"]].drop(columns="both_complete").reset_index(drop=True)}


def usable_metrics(d, min_finite=50):
    return [k[2:] for k in d.columns if k.startswith("o_") and np.isfinite(d[k].to_numpy(float)).sum() >= min_finite
            and d[k].nunique(dropna=True) > 1]


# ---------------------------------------------------------------------------
# strata, definition card, masking
# ---------------------------------------------------------------------------
def gene_class(g):
    return "framework" if g in FRAMEWORK else ("B_content" if g in B_CONTENT else "other")


def stratum_auc(views, min_events=3, min_units=50):
    rows = []
    for vn, d in views.items():
        if "gene" not in d:
            continue
        kinds = {"gene": d["gene"], "paralog_group": d["gene"].map(lambda g: GROUP.get(g, g)), "gene_class": d["gene"].map(gene_class)}
        if "content" in d and (d["content"] != "NA").any():
            kinds["hap_content"] = d["content"]
        for kind, lab in kinds.items():
            for s, idx in d.groupby(lab).groups.items():
                sub = d.loc[idx]
                y = sub["bad"].to_numpy(int)
                if y.sum() < min_events or len(sub) < min_units:
                    continue
                for k in usable_metrics(sub, min_finite=min_units):
                    o = sub["o_" + k].to_numpy(float)
                    a = S.auc_score(o, y)
                    rows.append(dict(view=vn, stratum_kind=kind, stratum=s, metric=k, n_units=len(sub), n_errors=int(y.sum()),
                                     base_error_pct=100 * y.mean(), auc=a))
    return pd.DataFrame(rows)


def definition_card(views):
    txt = {"B_child": "trio child call; bad = `neither` (protein field); informative locus, transmitting parent complete; origin >=5 k-mers and >=3x",
           "B_child_strong": "as B_child, origin >=10 k-mers and >=5x", "B_child_anytx": "as B_child, transmitting parent need not be complete",
           "B_worst": "as B_child, metrics = worse of child call and transmitting parent's worst call", "B_worst_strong": "B_worst, strong origin",
           "H_hap": "person x haplotype; bad = haplotype not complete (multi_contig or one_contig_open_end)",
           "P_pair_all": "first-degree pair x locus, both named; bad = no shared protein-field label",
           "P_pair_complete": "as P_pair_all, both people complete"}
    return pd.DataFrame([dict(view=k, definition=txt.get(k, ""), n_units=len(d), events=int(d["bad"].sum()),
                              base_error_pct=100 * d["bad"].mean() if len(d) else np.nan) for k, d in views.items()])


def mask_frame(df):
    """Counts 1-19 -> '<20' (as 77c); 0 and >=20 exact."""
    out = df.copy()
    for c in out.columns:
        if c in COUNT_COLS:
            out[c] = [C.mask_n(v) if pd.notna(v) and float(v) == int(float(v)) else (("<20" if 0 < float(v) < 20 else v) if pd.notna(v) else v)
                      for v in out[c]]
    return out


def mask_dir(src, dst):
    os.makedirs(dst, exist_ok=True)
    for fn in sorted(os.listdir(src)):
        if fn.endswith(".tsv"):
            df = pd.read_csv(os.path.join(src, fn), sep="\t", low_memory=False)
            C.write_tsv(mask_frame(df), os.path.join(dst, fn))


def power_text(power):
    lines = ["POWER STATEMENT (KIR pruning evidence)",
             "Only the events below carry information about which calls are wrong. With K events among n units, the smallest AUC detectable at",
             "one-sided alpha 0.05 and power 0.8 is 0.5 + 2.49 * sqrt((n + 1) / (12 K (n - K))); independent units assumed, so real power is lower.",
             "A flat curve or an AUC near 0.5 in a view with few events means 'not detectable', not 'does not predict'.", ""]
    for r in power.itertuples():
        lines.append("%s: n=%s events=%s min detectable AUC=%s" % (r.view, r.n_units, r.n_errors, ("%.3f" % r.auc_min_detectable)
                                                                     if np.isfinite(r.auc_min_detectable) else "NA"))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
def build_views(i77, i77c, i78):
    calls = read_int(os.path.join(i77, "kir_calls.INTERNAL.tsv"), dtype={"person_id": str, "contig": str})
    haps = read_int(os.path.join(i77, "kir_haps.INTERNAL.tsv"), dtype={"person_id": str})
    dpath = os.path.join(i78, "gene_depth_metrics.INTERNAL.tsv")
    depth = read_int(dpath, dtype={"person_id": str}) if os.path.exists(dpath) else None
    if depth is None:
        log("no 78 depth table at %s: depth metrics are NaN everywhere" % dpath)
    ct = call_table(calls, depth)
    pgw = person_gene_worst(ct)
    cont = hap_content(calls)
    trio = read_int(os.path.join(i77c, "records_trio_kir.INTERNAL.tsv"), dtype={"child": str, "tx": str, "other": str})
    pair = read_int(os.path.join(i77c, "records_pair_kir.INTERNAL.tsv"), dtype={"person": str, "other": str})
    views = trio_views(trio, ct, pgw, cont)
    tg = trio[trio["mode"] == "gene"]
    views.update(pair_views(pair, pgw, tg, None))
    views["H_hap"] = hap_table(ct, haps, calls)
    return views, ct


def run(i77, i77c, i78, internal_dir, out_dir, boot=200, folds=5, reps=20, grid_max=20.0, step=0.25, cost=25.0, window=2.0):
    views, ct = build_views(i77, i77c, i78)
    os.makedirs(internal_dir, exist_ok=True)
    for vn, d in views.items():
        log("%s: %d units, %d events" % (vn, len(d), int(d["bad"].sum())))
    unm = os.path.join(internal_dir, "unmasked")
    os.makedirs(unm, exist_ok=True)
    B76.make_views = lambda internal_dir_: views
    B76.usable_metrics = usable_metrics
    B76.TRANSFERS = TRANSFERS
    B76.run(internal_dir, unm, grid_max=grid_max, step=step, cost=cost, window=window, boot=boot, folds=folds, reps=reps)
    mask_dir(unm, out_dir)
    sa = stratum_auc(views)
    C.write_tsv(mask_frame(sa), os.path.join(out_dir, "per_stratum_auc.tsv"))
    C.write_tsv(mask_frame(definition_card(views)), os.path.join(out_dir, "definition_card_kir.tsv"))
    pw = pd.read_csv(os.path.join(unm, "power.tsv"), sep="\t")
    with open(os.path.join(out_dir, "POWER_STATEMENT.txt"), "w") as fh:
        fh.write(power_text(pw))
    cov = pd.DataFrame([dict(metric=k, has_depth_for_genes="only genes on GRCh38 chr19 (10 of 18); NaN for 2DL2, 2DL5A/B, 2DS1/2/3/5, 3DS1")
                        for k in DEPTH_METRICS])
    C.write_tsv(cov, os.path.join(out_dir, "depth_metric_coverage.tsv"))
    log("DONE 78b")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--int77", default=I77)
    ap.add_argument("--int77c", default=I77C)
    ap.add_argument("--int78", default=I78)
    ap.add_argument("--internal-dir", default=INT)
    ap.add_argument("--out-dir", default=OUT)
    ap.add_argument("--max-prune", type=float, default=20.0)
    ap.add_argument("--step", type=float, default=0.25)
    ap.add_argument("--cost-ratio", type=float, default=25.0)
    ap.add_argument("--nnp-window", type=float, default=2.0)
    ap.add_argument("--boot", type=int, default=200)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--reps", type=int, default=20)
    a = ap.parse_args(argv)
    run(a.int77, a.int77c, a.int78, a.internal_dir, a.out_dir, a.boot, a.folds, a.reps, a.max_prune, a.step, a.cost_ratio, a.nnp_window)
    return 0


if __name__ == "__main__":
    sys.exit(main())
