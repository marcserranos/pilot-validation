#!/usr/bin/env python3
"""09 -- W0: one analysis table, one row per person (report: reports/09_analysis_table/).

Every disease and ancestry analysis from here on (W1 V/J usage by ancestry, W4 disease
association, W5 classification) starts from this table, so inclusion rules, covariates and
depth-robust repertoire metrics are defined once and identically everywhere.

ROWS: the 7,922 people of cohort_full.tsv. Nobody is dropped here. Exclusions are boolean
flags, so every downstream analysis picks its own subset and the cost of each rule is visible
in the waterfall.

COLUMNS (VM-local file; see data_dictionary.csv for every column):
  identity      ancestry label, max ancestry probability, purity flags, PC1-10 (AoU)
  demographics  age (2020 - year of birth), sex at birth
  EHR           observation window in years, distinct conditions, conditions per EHR-year
                (an access-to-care proxy; it differs ~2x across ancestry groups)
  technical     TRB reads, TRB clonotypes, all-chain CDR3 counts (TRUST4), RNA-SeQC2 depth
                and quality metrics, rnaseq_metadata.tsv fields (discovered, see below)
  repertoire    raw Shannon / inverse Simpson / top-10 share (as Cole's PheWAS uses), and the
                same measures at a FIXED read depth (--depth, default 823 = report 02's
                rarefaction depth): exact rarefied richness (Hurlbert) plus Shannon, effective
                diversity e^H, inverse Simpson and public fraction averaged over --draws seeded
                subsamples. Fixed-depth metrics are NaN below the depth.
  relatedness   family_id (connected components of kinship >= --kin within this cohort) and
                unrelated (a maximal set with no pair at kinship >= --kin)
  flags         pass_read_floor, has_ehr_window, adult, has_sex, analysis_main (floor, EHR
                window, adult, unrelated). Sex is NOT a filter: sex_model is Female / Male /
                Unknown, because requiring a recorded sex removed 13% of AFR vs 1.5% of EAS
                (2026-10-08 run). has_sex stays available for sex-specific diseases.

NON-OBVIOUS CHOICES (full reasoning in DECISIONS.md, 2026-10-08):
  - Relatedness is resolved INSIDE this cohort, from AoU's pairwise samples_relatedness.tsv,
    not from AoU's global flagged-sample list: the global list removes people whose relative
    is not in our 7,922, which loses data and protects against nothing. Which relative to
    keep: the one with more TRB reads (better data), ties by id.
  - Ancestry probabilities are parsed from AoU's `probabilities` list, assumed in
    alphabetical order (AFR, AMR, EAS, EUR, MID, SAS). The script checks the assumption: the
    argmax must reproduce ancestry_pred for >= 99% of people, or it stops.
  - Public = (TRBV, CDR3 aa) carried by >= 2 of the 7,922. Public fraction is computed on
    the fixed-depth subsample, because the raw fraction falls with depth.
  - rnaseq_metadata.tsv is joined as-is (prefixed meta_) and its columns are profiled
    (name, % filled, n distinct) in metadata_columns.csv. Which one is the blood-draw date is
    decided after reading that profile, not guessed here.

INPUTS (VM): cohort_full.tsv, pheno/{person,observation_period,conditions}.tsv,
  atlas/trb_clonotypes_full.pkl, batch_summary_detail.tsv, qc/<RNA-SeQC2 metrics>.txt.gz,
  ~/aou_local/... copies of ancestry_preds.tsv, samples_relatedness.tsv, rnaseq_metadata.tsv
  and the RNA-seq manifest (run_w0.sh fetches them).

OUTPUTS: VM-local ~/pipeline_outputs/rnaseq/analysis/person_table.{tsv.gz,pkl};
  aggregate-only to --outdir: waterfall.csv, waterfall_by_floor.csv, column_summary.csv,
  data_dictionary.csv, metadata_columns.csv, summary.csv, fig_depth_and_diversity.

Usage (from aleix/RNA-seq/):
  pixi run python3 -u scripts/09_analysis_table.py [--depth 823] [--draws 10] [--kin 0.0884]
"""
import argparse
import ast
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import repfig as R  # noqa: E402

ANC = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]   # AoU probability order (checked below)
FLOORS = [250, 500, 823, 1000, 1500, 2000]
QC_KEEP = ["Total Reads", "Mapped Reads", "Non-Globin Reads", "Genes Detected",
           "Estimated Library Complexity", "Mapping Rate", "rRNA Rate", "Exonic Rate",
           "Duplicate Rate of Mapped", "Median 3' bias"]


def p(path):
    return os.path.expanduser(path)


def note(msg):
    print(msg, file=sys.stderr, flush=True)


# ------------------------------------------------------------------ ancestry

def load_ancestry(path, ids):
    a = pd.read_csv(p(path), sep="\t", dtype={"research_id": str})
    a = a[a["research_id"].isin(ids)].copy()
    out = pd.DataFrame({"research_id": a["research_id"].to_numpy()})
    if "probabilities" in a.columns:
        P = np.array([ast.literal_eval(str(v)) for v in a["probabilities"]], float)
        if P.shape[1] != len(ANC):
            R.die(f"ancestry probabilities have {P.shape[1]} entries, expected {len(ANC)}")
        pred = a["ancestry_pred"].astype(str).str.upper().to_numpy()
        agree = np.mean(np.array(ANC)[P.argmax(1)] == pred)
        note(f"Ancestry probabilities: argmax reproduces ancestry_pred for {100 * agree:.2f}%")
        if agree < 0.99:
            R.die("probability order assumption (AFR, AMR, EAS, EUR, MID, SAS) fails; "
                  "inspect the column before continuing")
        for i, g in enumerate(ANC):
            out[f"p_{g}"] = P[:, i]
        out["anc_maxprob"] = P.max(1)
    else:
        note("!! no `probabilities` column: purity flags will be missing")
    if "pca_features" in a.columns:
        PC = np.array([ast.literal_eval(str(v))[:10] for v in a["pca_features"]], float)
        for i in range(PC.shape[1]):
            out[f"PC{i + 1}"] = PC[:, i]
    return out, (float(agree) if "probabilities" in a.columns else np.nan)


# ------------------------------------------------------------------ relatedness

def resolve_relatedness(path, people, kin, eligible):
    """family_id = connected component; unrelated = greedy maximal set within the cohort.
    Greedy: first drop relatives who fail the other main filters (`eligible` False), so a
    pair never loses its only usable member; then repeatedly drop the person with most
    remaining relatives (ties: fewer TRB reads, then larger id)."""
    ids = set(people["research_id"])
    rel = pd.read_csv(p(path), sep="\t", dtype={"i.s": str, "j.s": str})
    for c in ("i.s", "j.s", "kin"):
        if c not in rel.columns:
            R.die(f"relatedness table lacks column {c}: {list(rel.columns)}")
    rel = rel[(rel["kin"] >= kin) & rel["i.s"].isin(ids) & rel["j.s"].isin(ids)]
    rel = rel[rel["i.s"] != rel["j.s"]]
    adj = {}
    for i, j in zip(rel["i.s"], rel["j.s"]):
        adj.setdefault(i, set()).add(j)
        adj.setdefault(j, set()).add(i)
    fam, fid = {}, 0
    for s in adj:
        if s in fam:
            continue
        fid += 1
        stack = [s]
        while stack:
            x = stack.pop()
            if x in fam:
                continue
            fam[x] = fid
            stack.extend(adj[x] - fam.keys())
    reads = people.set_index("research_id")["trb_reads"].fillna(0).to_dict()
    elig = dict(zip(people["research_id"], np.asarray(eligible, bool)))
    live = {k: set(v) for k, v in adj.items()}
    dropped = set()
    while any(live.values()):
        x = max((k for k, v in live.items() if v),
                key=lambda k: (not elig.get(k, False), len(live[k]), -reads.get(k, 0), k))
        dropped.add(x)
        for y in live.pop(x):
            live[y].discard(x)
    people["family_id"] = people["research_id"].map(fam).fillna(0).astype(int)
    people["family_id"] = np.where(people["family_id"] > 0, people["family_id"],
                                   -np.arange(1, len(people) + 1))   # singletons: own group
    people["unrelated"] = ~people["research_id"].isin(dropped)
    return {"related_pairs_in_cohort": len(rel), "people_with_relative": len(adj),
            "families": fid, "dropped_for_relatedness": len(dropped)}


# ------------------------------------------------------------------ repertoire metrics

def repertoire_metrics(clono, depth, draws, seed):
    """Raw and fixed-depth diversity per person from the full TRB clonotype cache."""
    clono = clono[["research_id", "trbv", "cdr3aa", "reads"]].copy()
    clono["research_id"] = clono["research_id"].astype(str)
    r = pd.to_numeric(clono["reads"], errors="coerce").fillna(0)
    frac_nonint = float(np.mean(r != np.round(r)))
    clono["reads_int"] = np.maximum(np.round(r), 1).astype(np.int64)
    key = clono["trbv"].astype(str) + "|" + clono["cdr3aa"].astype(str)
    clono["public"] = key.map(clono.groupby(key)["research_id"].nunique()) >= 2
    rng = np.random.default_rng(seed)
    rows = []
    t0 = time.time()
    for i, (rid, g) in enumerate(clono.groupby("research_id", sort=True), 1):
        n = g["reads_int"].to_numpy()
        N = int(n.sum())
        pr = n / N
        pub = g["public"].to_numpy()
        row = {"research_id": rid, "trb_clonotypes": len(n), "trb_reads": N,
               "shannon_raw": float(-np.sum(pr * np.log(pr))),
               "inv_simpson_raw": float(1 / np.sum(pr ** 2)),
               "top10_frac": float(np.sort(n)[::-1][:10].sum() / N),
               "public_frac_raw": float(pub.mean())}
        if N >= depth:
            row["richness_at_depth"] = R.rarefied_richness(n, depth)
            H, D2, PF = [], [], []
            for _ in range(draws):
                s = rng.multivariate_hypergeometric(n, depth)
                k = s > 0
                q = s[k] / depth
                H.append(-np.sum(q * np.log(q)))
                D2.append(1 / np.sum(q ** 2))
                PF.append(pub[k].mean())
            row.update(shannon_at_depth=float(np.mean(H)),
                       eff_diversity_at_depth=float(np.mean(np.exp(H))),
                       inv_simpson_at_depth=float(np.mean(D2)),
                       public_frac_at_depth=float(np.mean(PF)))
        rows.append(row)
        if i % 1000 == 0:
            note(f"  repertoire metrics {i:,} people ({time.time() - t0:.0f}s)")
    return pd.DataFrame(rows), frac_nonint


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    base = "~/pipeline_outputs/rnaseq"
    ap.add_argument("--cohort", default=f"{base}/cohort_full.tsv")
    ap.add_argument("--pheno-dir", default=f"{base}/pheno")
    ap.add_argument("--cache-dir", default=f"{base}/atlas")
    ap.add_argument("--batch-detail", default=f"{base}/batch_summary_detail.tsv")
    ap.add_argument("--qc", default=f"{base}/qc/aou_rnaseq_20260413.metrics.txt.gz")
    ap.add_argument("--aou-local", default="~/aou_local/v9",
                    help="local copies of AoU files, same relative layout as the bucket")
    ap.add_argument("--depth", type=int, default=823,
                    help="fixed read depth for rarefied metrics and the main read floor")
    ap.add_argument("--draws", type=int, default=10)
    ap.add_argument("--kin", type=float, default=0.0884,
                    help="kinship threshold for 'related' (0.0884 = second degree, KING)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--local-out", default=f"{base}/analysis")
    ap.add_argument("--outdir", default=f"{base}/reports/09_analysis_table")
    args = ap.parse_args()
    outdir, local = p(args.outdir), p(args.local_out)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(local, exist_ok=True)
    aou = p(args.aou_local)
    summary = [("depth", args.depth), ("draws", args.draws), ("kin_threshold", args.kin)]

    # ---- people, demographics
    people = R.load_people(args.cohort, args.pheno_dir)[["research_id", "ancestry", "sex", "age"]]
    ids = set(people["research_id"])
    note(f"Cohort: {len(people):,} people")

    # ---- repertoire
    t0 = time.time()
    clono = pd.read_pickle(os.path.join(p(args.cache_dir), "trb_clonotypes_full.pkl"))
    rep, frac_nonint = repertoire_metrics(clono, args.depth, args.draws, args.seed)
    del clono
    people = people.merge(rep, on="research_id", how="left")
    note(f"Repertoire metrics: {time.time() - t0:.0f}s; non-integer read counts {100 * frac_nonint:.2f}% "
         f"(rounded, min 1)")
    summary.append(("frac_nonint_reads_rounded", frac_nonint))

    # ---- all-chain TRUST4 counts
    if os.path.exists(p(args.batch_detail)):
        bd = pd.read_csv(p(args.batch_detail), sep="\t", dtype={"research_id": str})
        keep = ["research_id", "n_cdr3"] + [c for c in bd.columns
                                            if c.startswith("n_") and not c.endswith("_vonly")
                                            and c != "n_cdr3"]
        people = people.merge(bd[keep], on="research_id", how="left")

    # ---- ancestry probabilities and PCs
    anc, agree = load_ancestry(os.path.join(aou, "wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"), ids)
    people = people.merge(anc, on="research_id", how="left")
    summary.append(("ancestry_prob_argmax_agreement", agree))
    for t in (0.90, 0.95):
        people[f"pure{int(t * 100)}"] = people.get("anc_maxprob", pd.Series(np.nan)) >= t

    # ---- EHR window and burden
    obs = pd.read_csv(os.path.join(p(args.pheno_dir), "observation_period.tsv"), sep="\t",
                      dtype={"research_id": str})
    cond = pd.read_csv(os.path.join(p(args.pheno_dir), "conditions.tsv"), sep="\t",
                       dtype={"research_id": str}, usecols=["research_id", "condition_concept_id"])
    burden = cond.groupby("research_id")["condition_concept_id"].nunique().rename("n_conditions")
    people = (people.merge(obs[["research_id", "ehr_years"]], on="research_id", how="left")
                    .merge(burden, left_on="research_id", right_index=True, how="left"))
    people["n_conditions"] = people["n_conditions"].fillna(0).astype(int)
    people["conditions_per_ehr_year"] = np.where(people["ehr_years"] > 0,
                                                 people["n_conditions"] / people["ehr_years"], np.nan)

    # ---- RNA-SeQC2 metrics (QC sample_id -> research_id through the manifest if needed)
    if os.path.exists(p(args.qc)):
        qc = pd.read_csv(p(args.qc), sep="\t", dtype=str, compression="gzip")
        if len(ids & set(qc["sample_id"])) < 0.5 * len(ids):
            man = pd.read_csv(os.path.join(aou, "multiomics/rnaseq/manifest.tsv"), sep="\t", dtype=str)
            qc = qc.merge(man[["sampleid", "research_id"]], left_on="sample_id", right_on="sampleid")
        else:
            qc = qc.rename(columns={"sample_id": "research_id"})
        cols = [c for c in QC_KEEP if c in qc.columns]
        q = qc[["research_id"] + cols].drop_duplicates("research_id").copy()
        for c in cols:
            q[c] = pd.to_numeric(q[c], errors="coerce")
        q.columns = ["research_id"] + ["qc_" + c.lower().replace(" ", "_").replace("'", "")
                                       for c in cols]
        if {"qc_non-globin_reads", "qc_mapped_reads"} <= set(q.columns):
            q["qc_globin_fraction"] = 1 - q["qc_non-globin_reads"] / q["qc_mapped_reads"]
        people = people.merge(q, on="research_id", how="left")
    else:
        note(f"!! no RNA-SeQC2 table at {args.qc}")

    # ---- rnaseq_metadata.tsv: joined as-is, profiled, not interpreted
    meta_path = os.path.join(aou, "multiomics/rnaseq/rnaseq_metadata.tsv")
    if os.path.exists(meta_path):
        meta = pd.read_csv(meta_path, sep="\t", dtype=str)
        idcol = "research_id" if "research_id" in meta.columns else None
        sid = next((c for c in ("sampleid", "sample_id") if c in meta.columns), None)
        if idcol is None and sid:
            man = pd.read_csv(os.path.join(aou, "multiomics/rnaseq/manifest.tsv"), sep="\t", dtype=str)
            meta = meta.merge(man[["sampleid", "research_id"]], left_on=sid, right_on="sampleid",
                              how="left", suffixes=("", "_manifest"))
            idcol = "research_id"
        if idcol is None:
            note(f"!! rnaseq_metadata.tsv has no id column to join on: {list(meta.columns)}")
        if idcol:
            meta = meta[meta["research_id"].isin(ids)].drop_duplicates("research_id")
            prof = pd.DataFrame({
                "column": [c for c in meta.columns
                           if c not in ("research_id", "sampleid", "sample_id", "sampleid_manifest")]})
            prof["pct_filled"] = [100 * meta[c].notna().mean() for c in prof["column"]]
            prof["n_distinct"] = [meta[c].nunique() for c in prof["column"]]
            prof["looks_like_date"] = [bool(pd.to_datetime(meta[c].dropna().head(200),
                                                           errors="coerce", format="mixed").notna().mean() > 0.8)
                                       if meta[c].notna().any() else False for c in prof["column"]]
            prof.to_csv(os.path.join(outdir, "metadata_columns.csv"), index=False)
            meta = meta.rename(columns={c: f"meta_{c}" for c in meta.columns if c != "research_id"})
            people = people.merge(meta, on="research_id", how="left")
            note("rnaseq_metadata.tsv columns:\n" + prof.to_string(index=False))
    else:
        note(f"!! no rnaseq_metadata.tsv at {meta_path}")

    # ---- flags (before relatedness, which prefers to keep the eligible relative)
    people["pass_read_floor"] = people["trb_reads"] >= args.depth
    people["has_ehr_window"] = people["ehr_years"] > 0
    people["adult"] = people["age"] >= 18
    people["has_sex"] = people["sex"].isin(["Female", "Male"])
    people["sex_model"] = people["sex"].where(people["has_sex"], "Unknown")
    eligible = people["pass_read_floor"] & people["has_ehr_window"] & people["adult"]

    # ---- relatedness
    rel_path = os.path.join(aou, "wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv")
    if os.path.exists(rel_path):
        rs = resolve_relatedness(rel_path, people, args.kin, eligible)
        summary += list(rs.items())
        note(f"Relatedness: {rs}")
    else:
        note(f"!! no relatedness table at {rel_path}: unrelated set to True for everyone")
        people["family_id"] = -np.arange(1, len(people) + 1)
        people["unrelated"] = True
        summary.append(("relatedness", "missing"))

    people["analysis_main"] = (people["pass_read_floor"] & people["has_ehr_window"]
                               & people["adult"] & people["unrelated"])

    # ---- write the person table (VM-local only)
    people.to_pickle(os.path.join(local, "person_table.pkl"))
    people.to_csv(os.path.join(local, "person_table.tsv.gz"), sep="\t", index=False)
    note(f"Person table: {people.shape[0]:,} x {people.shape[1]} -> {local} (VM-local)")

    # ---- aggregate outputs
    steps = [("cohort", np.ones(len(people), bool)),
             (f"TRB reads >= {args.depth}", people["pass_read_floor"]),
             ("EHR window > 0", people["has_ehr_window"]),
             ("adult (>= 18)", people["adult"]),
             ("unrelated", people["unrelated"])]
    m = np.ones(len(people), bool)
    wf = []
    for name, f in steps:
        m = m & np.asarray(f, bool)
        row = {"step": name, "all": int(m.sum())}
        for a in ANC:
            k = int((m & (people["ancestry"] == a)).sum())
            row[a] = k if k >= R.MIN_PEOPLE or k == 0 else "<20"
        wf.append(row)
    for t in (90, 95):
        mm = m & people[f"pure{t}"].fillna(False).to_numpy()
        row = {"step": f"(stratified figures only) ancestry prob >= 0.{t}", "all": int(mm.sum())}
        for a in ANC:
            k = int((mm & (people["ancestry"] == a)).sum())
            row[a] = k if k >= R.MIN_PEOPLE or k == 0 else "<20"
        wf.append(row)
    pd.DataFrame(wf).to_csv(os.path.join(outdir, "waterfall.csv"), index=False)
    note("\n" + pd.DataFrame(wf).to_string(index=False))

    fl = []
    for f in FLOORS:
        k = people["trb_reads"] >= f
        row = {"trb_read_floor": f, "kept": int(k.sum()), "pct_kept": round(100 * k.mean(), 2)}
        for a in ANC:
            g = people["ancestry"] == a
            row[f"pct_kept_{a}"] = round(100 * float(k[g].mean()), 2)
        fl.append(row)
    pd.DataFrame(fl).to_csv(os.path.join(outdir, "waterfall_by_floor.csv"), index=False)

    num = [c for c in people.columns if c not in ("research_id", "family_id")
           and not c.startswith("meta_") and pd.api.types.is_numeric_dtype(people[c])
           and not pd.api.types.is_bool_dtype(people[c])]
    cs = []
    for c in num:
        for a in ["ALL"] + ANC:
            v = people[c] if a == "ALL" else people.loc[people["ancestry"] == a, c]
            v = v.dropna()
            if len(v) >= R.MIN_PEOPLE:
                cs.append((c, a, len(v), v.quantile(.05), v.quantile(.25), v.median(),
                           v.quantile(.75), v.quantile(.95)))
    pd.DataFrame(cs, columns=["column", "ancestry", "n", "p05", "p25", "median", "p75", "p95"]
                 ).round(5).to_csv(os.path.join(outdir, "column_summary.csv"), index=False)

    dd = pd.DataFrame({"column": people.columns,
                       "dtype": [str(people[c].dtype) for c in people.columns],
                       "pct_filled": [round(100 * people[c].notna().mean(), 2) for c in people.columns]})
    dd.to_csv(os.path.join(outdir, "data_dictionary.csv"), index=False)

    # ---- does fixed depth remove the depth dependence? (Spearman with log reads)
    from scipy.stats import spearmanr
    lr = np.log10(people["trb_reads"])
    for c in ("shannon_raw", "shannon_at_depth", "inv_simpson_raw", "inv_simpson_at_depth",
              "trb_clonotypes", "richness_at_depth", "public_frac_raw", "public_frac_at_depth"):
        ok = people[c].notna() & lr.notna()
        summary.append((f"spearman_{c}_vs_log_reads", float(spearmanr(people.loc[ok, c], lr[ok])[0])))
    summary += [("people", len(people)), ("analysis_main", int(people["analysis_main"].sum()))]
    pd.DataFrame(summary, columns=["metric", "value"]).to_csv(os.path.join(outdir, "summary.csv"),
                                                              index=False)

    # ---- figure: diversity vs depth, raw vs fixed depth; read-floor cost by ancestry
    import cnsplots as cns
    mp = cns.multipanel(max_width=540)
    col = R.ancestry_colors()
    bins = np.quantile(lr.dropna(), np.linspace(0, 1, 21))
    people["_bin"] = pd.cut(lr, np.unique(bins), include_lowest=True)

    def binned(ax, c, title, ylab):
        g = people.dropna(subset=[c]).groupby("_bin", observed=True)
        s = pd.DataFrame({"x": g["trb_reads"].median(), "n": g.size(), "med": g[c].median(),
                          "lo": g[c].quantile(.25), "hi": g[c].quantile(.75)})
        s = s[s["n"] >= R.MIN_PEOPLE]
        ax.fill_between(np.log10(s["x"]), s["lo"], s["hi"], color=R.nature()[3], alpha=0.25, lw=0)
        ax.plot(np.log10(s["x"]), s["med"], "-o", ms=2.5, lw=0.9, color=R.nature()[3])
        ax.axvline(np.log10(args.depth), color="0.5", lw=0.5, ls="--")
        R.style(ax, title, "log$_{10}$ TRB reads", ylab)

    ax = mp.panel("a", width=120, height=105, margin_right=40, margin_bottom=34)
    binned(ax, "shannon_raw", "Shannon, raw", "Shannon (nats)")
    ax = mp.panel("b", width=120, height=105, margin_right=40, margin_bottom=34)
    binned(ax, "shannon_at_depth", f"Shannon at {args.depth} reads", "Shannon (nats)")
    ax = mp.panel("c", width=150, height=105, margin_right=12, margin_bottom=34)
    for a in ANC:
        ax.plot(FLOORS, [r[f"pct_kept_{a}"] for r in fl], "-o", ms=2.5, lw=0.9, color=col[a], label=a)
    ax.axvline(args.depth, color="0.5", lw=0.5, ls="--")
    ax.set_xscale("log")
    ax.set_xticks(FLOORS)
    ax.set_xticklabels([str(f) for f in FLOORS], rotation=45)
    ax.minorticks_off()
    ax.legend(frameon=False, fontsize=5.5, ncol=2)
    R.style(ax, "Cost of a read floor", "Minimum TRB reads", "People kept (%)")
    R.save("fig_depth_and_diversity", outdir)
    note(f"\nAll aggregate outputs in {outdir}")


if __name__ == "__main__":
    main()
