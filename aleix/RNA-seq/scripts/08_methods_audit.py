#!/usr/bin/env python3
"""08 -- Methods audit: the numbers a methods section needs that no report produced yet
(report: reports/08_methods_audit/).

  1  FILTER STEPS. Rows surviving each step of the shared clonotype filter
     (embed_cdr3s.load_person_cdr3s): raw cdr3.out rows -> CDR3_score >= 0.02 -> canonical
     junction -> <= 30 aa -> chain called -> unique clonotypes -> TRB clonotypes.
     Cohort totals and per-person medians.

  2  TIES AT THE TOP-500 CUTOFF. The embedding pool keeps each person's 500 TRB clonotypes
     with most reads. Until 2026-10-04 ties were broken alphabetically by (chain, V, CDR3):
     groupby sorted the keys and the read sort was stable. With most clonotypes at 1-2 reads
     the cutoff almost always falls inside a tie, so the old pool favoured alphabetically
     early TRBV genes. This measures it: share of people whose cutoff is inside a tie, how
     much of the top 500 changes under the seeded random tie-break, and per-TRBV share of
     the full repertoire vs the old pool vs the new pool. If the old pool file is still in
     place, it is also checked that the old selection is reproduced exactly.

  3  THE 12 GB BAM CAP. People in the LR x RNA-seq overlap left out of the cohort for BAM
     size, compared with those kept on AoU's RNA-SeQC2 metrics (depth, complexity, genes
     detected...). TRUST4 never ran on them, so only QC can be compared.

Run BEFORE re-embedding (run_polish.sh does), so check 2 can compare against the old pool.
All outputs are aggregate; every reported group pools >= 20 people.

Usage (from aleix/RNA-seq/):
  pixi run python3 -u scripts/08_methods_audit.py [--workers 16]
"""
import argparse
import os
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import mannwhitneyu  # noqa: E402

from embed_cdr3s import load_person_cdr3s  # noqa: E402

MIN_PEOPLE = 20
CAP = 500
STEPS = ["rows_raw", "rows_score", "rows_canonical", "rows_len_ok", "rows_chain_called",
         "clonotypes", "clonotypes_trb"]
STEP_LABELS = {"rows_raw": "cdr3.out rows", "rows_score": "CDR3_score >= 0.02",
               "rows_canonical": "canonical junction C...F/W",
               "rows_len_ok": "length <= 30 aa", "rows_chain_called": "chain called (V/J/C)",
               "clonotypes": "unique clonotypes (chain, V, CDR3 aa)",
               "clonotypes_trb": "TRB clonotypes"}
QC_METRICS = ["Mapped Reads", "Total Reads", "Non-Globin Reads", "Genes Detected",
              "Estimated Library Complexity", "Mapping Rate", "rRNA Rate", "Exonic Rate",
              "Duplicate Rate of Mapped", "Median 3' bias"]


def strip(v):
    return v.split("*")[0] if v else "(no V call)"


def audit_one(job):
    trust4_dir, rid = job
    st = {}
    df, _ = load_person_cdr3s(trust4_dir, rid, 0.02, False, False, 30, stats=st)
    if df is None:
        return rid, None
    d = df[df["chain"] == "TRB"]
    old = (d.sort_values(["chain", "v_gene", "cdr3aa"])
            .sort_values("reads", ascending=False, kind="stable"))
    new_top, old_top = d.head(CAP), old.head(CAP)
    out = {"stats": st, "n_trb": len(d),
           "v_full": Counter(strip(v) for v in d["v_gene"]),
           "v_old": Counter(strip(v) for v in old_top["v_gene"]),
           "v_new": Counter(strip(v) for v in new_top["v_gene"])}
    if len(d) > CAP:
        r = d["reads"].to_numpy()
        r_cut = r[CAP - 1]
        out["reads_at_cutoff"] = float(r_cut)
        out["tied_total"] = int((r == r_cut).sum())
        out["tied_kept"] = int((r[:CAP] == r_cut).sum())
        k_old = set(zip(old_top["v_gene"], old_top["cdr3aa"]))
        k_new = set(zip(new_top["v_gene"], new_top["cdr3aa"]))
        out["frac_changed"] = 1 - len(k_old & k_new) / CAP
    return rid, out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", default="~/pipeline_outputs/rnaseq/cohort_full.tsv")
    ap.add_argument("--overlap", default="~/pipeline_outputs/rnaseq/lr_rnaseq_overlap_cohort.tsv")
    ap.add_argument("--trust4-dir", default="~/pipeline_outputs/rnaseq")
    ap.add_argument("--embeddings-dir", default="~/pipeline_outputs/rnaseq/embeddings")
    ap.add_argument("--old-pool", default=None,
                    help="old top-500 pool TSV to check the reconstruction against (default: "
                         "pool_sceptr_cohort_full_vcdr3.tsv, or its v1_alphabetical_ties backup)")
    ap.add_argument("--qc", default="~/pipeline_outputs/rnaseq/qc/aou_rnaseq_20260413.metrics.txt.gz")
    ap.add_argument("--mount", default="~/mnt/aou-controlled")
    ap.add_argument("--outdir", default="~/pipeline_outputs/rnaseq/reports/08_methods_audit")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    args = ap.parse_args()
    outdir = os.path.expanduser(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    summary = []

    cohort = pd.read_csv(os.path.expanduser(args.cohort), sep="\t", dtype=str)
    jobs = [(os.path.expanduser(args.trust4_dir), r) for r in cohort["research_id"]]

    # ---------------- 1 + 2: one pass over every cdr3.out ----------------
    t0 = time.time()
    res = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for i, (rid, o) in enumerate(ex.map(audit_one, jobs, chunksize=16), 1):
            if o is not None:
                res.append(o)
            if i % 1000 == 0:
                print(f"  audited {i:,}/{len(jobs):,} ({time.time() - t0:.0f}s)", file=sys.stderr)
    n = len(res)
    print(f"Audited {n:,} people in {time.time() - t0:.0f}s", file=sys.stderr)
    summary.append(("people_audited", n))

    S = pd.DataFrame([o["stats"] for o in res])[STEPS]
    tot = S.sum()
    steps = pd.DataFrame({"step": [STEP_LABELS[s] for s in STEPS],
                          "rows_total": tot.to_numpy(),
                          "pct_of_raw": (100 * tot / tot["rows_raw"]).round(2).to_numpy(),
                          "pct_of_previous": (100 * tot / tot.shift(1)).round(2).to_numpy(),
                          "median_per_person": S.median().to_numpy()})
    steps.to_csv(os.path.join(outdir, "filter_steps.csv"), index=False)
    print("\n" + steps.to_string(index=False), file=sys.stderr)

    capped = [o for o in res if "reads_at_cutoff" in o]
    rc = np.array([o["reads_at_cutoff"] for o in capped])
    tied_total = np.array([o["tied_total"] for o in capped])
    tied_kept = np.array([o["tied_kept"] for o in capped])
    changed = np.array([o["frac_changed"] for o in capped])
    inside = tied_total > tied_kept
    ties = [("people_with_more_than_500_trb", len(capped)),
            ("frac_people_cutoff_inside_tie", float(inside.mean())),
            ("frac_cutoff_at_1_read", float(np.mean(rc == 1))),
            ("frac_cutoff_at_2_reads", float(np.mean(rc == 2))),
            ("frac_cutoff_at_3plus_reads", float(np.mean(rc >= 3))),
            ("median_tied_block_size", float(np.median(tied_total))),
            ("median_tied_kept_of_block", float(np.median(tied_kept))),
            ("mean_frac_top500_changed", float(changed.mean())),
            ("median_frac_top500_changed", float(np.median(changed)))]

    def shares(key):
        c = Counter()
        for o in res:
            c.update(o[key])
        s = pd.Series(c, dtype=float)
        return s / s.sum(), s
    full_sh, _ = shares("v_full")
    old_sh, old_ct = shares("v_old")
    new_sh, new_ct = shares("v_new")
    carriers = Counter()
    for o in res:
        carriers.update(set(o["v_full"]))
    shift = pd.DataFrame({"people_carrying": pd.Series(carriers),
                          "pct_full_repertoire": 100 * full_sh,
                          "pct_old_pool": 100 * old_sh,
                          "pct_new_pool": 100 * new_sh}).fillna(0)
    shift = shift[shift["people_carrying"] >= MIN_PEOPLE]
    shift["old_minus_new_pp"] = shift["pct_old_pool"] - shift["pct_new_pool"]
    shift["old_over_new"] = shift["pct_old_pool"] / shift["pct_new_pool"].replace(0, np.nan)
    shift = shift.sort_values("old_minus_new_pp", ascending=False).round(4)
    shift.index.name = "trbv"
    shift.to_csv(os.path.join(outdir, "trbv_pool_shift.csv"))
    ties.append(("total_variation_old_vs_new_pool_pct",
                 float(0.5 * (shift["pct_old_pool"] - shift["pct_new_pool"]).abs().sum())))
    print("\nLargest TRBV shifts, old (alphabetical ties) vs new (random ties) pool:\n"
          + shift.head(8).to_string() + "\n...\n" + shift.tail(5).to_string(), file=sys.stderr)

    # does the old rule reproduce the pool that was actually embedded?
    edir = os.path.expanduser(args.embeddings_dir)
    cands = [args.old_pool] if args.old_pool else [
        os.path.join(edir, "v1_alphabetical_ties", "pool_sceptr_cohort_full_vcdr3.tsv"),
        os.path.join(edir, "pool_sceptr_cohort_full_vcdr3.tsv")]
    old_pool = next((os.path.expanduser(c) for c in cands if c and os.path.exists(os.path.expanduser(c))), None)
    if old_pool:
        op = pd.read_csv(old_pool, sep="\t", usecols=["v_gene"], keep_default_na=False)
        oc = op["v_gene"].map(strip).value_counts()
        diff = (oc - old_ct.reindex(oc.index).fillna(0)).abs()
        ties += [("old_pool_file", os.path.basename(os.path.dirname(old_pool)) + "/" +
                  os.path.basename(old_pool)),
                 ("old_pool_rows", len(op)),
                 ("old_pool_reconstruction_max_abs_diff_per_gene", float(diff.max())),
                 ("old_pool_reconstruction_total_abs_diff", float(diff.sum()))]
    else:
        ties.append(("old_pool_file", "not found -- reconstruction not checked"))
    pd.DataFrame(ties, columns=["metric", "value"]).to_csv(
        os.path.join(outdir, "tie_summary.csv"), index=False)
    print("\n" + pd.DataFrame(ties, columns=["metric", "value"]).to_string(index=False),
          file=sys.stderr)
    summary += ties

    # ---------------- 3: the BAM-size cap ----------------
    ov_path, qc_path = os.path.expanduser(args.overlap), os.path.expanduser(args.qc)
    if os.path.exists(ov_path) and os.path.exists(qc_path):
        ov = pd.read_csv(ov_path, sep="\t", dtype=str)
        if "ancestry" not in ov.columns:
            ov = ov.merge(cohort[["research_id", "ancestry"]], on="research_id", how="left")
        ov = ov[ov["ancestry"].isin(["AFR", "AMR", "EAS", "EUR", "MID", "SAS"])].copy()
        ov["group"] = np.where(ov["research_id"].isin(cohort["research_id"]), "kept", "excluded")
        qc = pd.read_csv(qc_path, sep="\t", dtype=str, compression="gzip")
        if len(set(ov["research_id"]) & set(qc["sample_id"])) >= 0.5 * len(ov):
            qc = qc.rename(columns={"sample_id": "research_id"})
        else:
            man = os.path.join(os.path.expanduser(args.mount), "v9/multiomics/rnaseq/manifest.tsv")
            if os.path.exists(man):
                m = pd.read_csv(man, sep="\t", dtype=str)[["sampleid", "research_id"]]
                qc = qc.merge(m, left_on="sample_id", right_on="sampleid")
            else:
                qc = None
                summary.append(("bam_cap", "QC ids need the manifest; mount not up -- skipped"))
        if qc is not None:
            mets = [c for c in QC_METRICS if c in qc.columns]
            for c in mets:
                qc[c] = pd.to_numeric(qc[c], errors="coerce")
            j = ov.merge(qc[["research_id"] + mets], on="research_id", how="inner")
            rows = []
            for c in mets:
                k_, e_ = j.loc[j["group"] == "kept", c].dropna(), j.loc[j["group"] == "excluded", c].dropna()
                if len(k_) >= MIN_PEOPLE and len(e_) >= MIN_PEOPLE:
                    rows.append((c, len(k_), len(e_), k_.median(), e_.median(),
                                 e_.median() / k_.median() if k_.median() else np.nan,
                                 mannwhitneyu(k_, e_).pvalue))
            cap = pd.DataFrame(rows, columns=["metric", "n_kept", "n_excluded", "median_kept",
                                              "median_excluded", "ratio_excluded_over_kept", "p_mwu"])
            cap.to_csv(os.path.join(outdir, "bam_cap_qc.csv"), index=False)
            by = j.groupby(["ancestry", "group"]).size().unstack(fill_value=0)
            by["pct_excluded"] = 100 * by.get("excluded", 0) / by.sum(axis=1)
            small = (by[["excluded", "kept"]] > 0) & (by[["excluded", "kept"]] < MIN_PEOPLE)
            by = by.astype(object)
            by[small.any(axis=1).to_numpy()] = "<20 suppressed"
            by.to_csv(os.path.join(outdir, "bam_cap_by_ancestry.csv"))
            dcol = next((c for c in ("Mapped Reads", "Total Reads", "Non-Globin Reads") if c in mets), None)
            if dcol:
                top = j[j[dcol] >= j[dcol].quantile(0.9)]
                summary += [("bam_cap_depth_metric", dcol),
                            ("bam_cap_pct_excluded_overall", 100 * float((j["group"] == "excluded").mean())),
                            ("bam_cap_pct_excluded_in_top_depth_decile",
                             100 * float((top["group"] == "excluded").mean()))]
            print("\n" + cap.to_string(index=False) + "\n" + by.to_string(), file=sys.stderr)
    else:
        summary.append(("bam_cap", f"skipped: need {ov_path} and {qc_path}"))

    pd.DataFrame(summary, columns=["metric", "value"]).to_csv(
        os.path.join(outdir, "summary.csv"), index=False)
    print(f"\nAll outputs in {outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
