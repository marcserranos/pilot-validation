#!/usr/bin/env python3
"""The five (plus a few) untested QC confounds for the ancestry recovery gap -- N4.

Background: ../results/rnaseq_100person_ancestry_writeup.md found AFR recovering 1.63x
more CDR3s than EAS, of which depth explains ~40%, with RQS and a cell-composition proxy
ruled out and Kruskal-Wallis p = 0.18 (not significant at n=16-17/group). AoU's RNA-SeQC2
metrics file ships ~82 per-sample technical measurements we have never looked at. This
checks all the plausible ones at once. Zero new compute -- every number was computed by
AoU in April 2026 and sits in a 2.9 MB file.

THE TEST, and it is the same one used throughout this workstream: a variable is only a
credible confound if it does BOTH
  (a) correlate with recovery (CDR3 count), AND
  (b) vary across the ancestry groups being compared.
Something that correlates with recovery but is flat across groups cannot produce a
between-group gap. Something that varies across groups but doesn't predict recovery
cannot either. This script reports both, side by side, so neither can be quietly assumed.

METRICS CHECKED, and why each is a real hypothesis rather than a fishing expedition:
  Median/Mean 3' bias   -- how badly RNA degraded from the 3' end. THE mechanistic
                           favourite: CDR3 sits mid-transcript, so 3' degradation should
                           hurt repertoire recovery specifically, more than it hurts
                           ordinary expression.
  Estimated Library Complexity -- distinct molecules before PCR. Arguably a better
                           denominator than raw read count: 100M reads off a
                           low-complexity library is not 100M reads' worth of information.
  Expression Profiling Efficiency -- AoU's own "ultimate benchmark" metric.
  rRNA Rate, globin_fraction   -- how well the two depletion chemistries worked. Depletion
                           is a wet-lab step with per-sample variation; poor depletion
                           spends the read budget on junk.
  Exonic/Intronic/Intergenic Rate, Duplicate rates, Unique Rate, Mapping Rate,
  Genes Detected, Fragment GC, Transcript Coverage CV -- standard technical axes.
  Read Length           -- *** checked separately and first ***. TRUST4's unmapped-read
                           k-mer screen requires an exact match of readLen/5 bases
                           (TRUST4_DEEP_DIVE.md section 2.2). If read length is NOT
                           constant across samples, then extraction stringency literally
                           differs per person, which would be a mechanistically direct,
                           previously unsuspected recovery confound. If it is constant at
                           151, that hypothesis is dead and we can stop worrying about it.

Outputs, same privacy split as every other script here:
  1. VM-local per-person joined detail -- not committed.
  2. DE-IDENTIFIED per-ancestry metric summary + a ranked candidate-confound table --
     safe to commit, written to ../results/.

Usage:
  pixi run python3 check_qc_confounds.py <cohort.tsv> [<cohort2.tsv> ...] \\
      [--qc ~/pipeline_outputs/rnaseq/qc/aou_rnaseq_20260413.metrics.txt.gz] \\
      [--detail ~/pipeline_outputs/rnaseq/batch_summary_detail.tsv] \\
      [--mount ~/mnt/aou-controlled]

Multiple cohort files may be passed -- pass both 100-person cohorts to analyse them
together and roughly double the per-group n.
"""
import argparse
import os
import sys

import pandas as pd
from scipy.stats import kruskal

RNASEQ_MANIFEST = "v9/multiomics/rnaseq/manifest.tsv"
QC_DEFAULT = os.path.expanduser(
    "~/pipeline_outputs/rnaseq/qc/aou_rnaseq_20260413.metrics.txt.gz")
DETAIL_DEFAULT = os.path.expanduser("~/pipeline_outputs/rnaseq/batch_summary_detail.tsv")
ANCESTRY_GROUPS = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]

# Exact column names as they appear in the file -- spaces and apostrophes included.
# Verified live 2026-08-17 against the real header; do not "tidy" these.
METRICS = [
    "Median 3' bias",
    "Mean 3' bias",
    "Estimated Library Complexity",
    "Expression Profiling Efficiency",
    "rRNA Rate",
    "Exonic Rate",
    "Intronic Rate",
    "Intergenic Rate",
    "Duplicate Rate of Mapped",
    "Duplicate Rate of Mapped, excluding Globins",
    "Unique Rate of Mapped",
    "Mapping Rate",
    "Genes Detected",
    "Median of Transcript Coverage CV",
    "Fragment GC Content Mean",
    "Base Mismatch",
]
DERIVED = ["globin_fraction"]


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cohorts", nargs="+", help="one or more cohort.tsv files")
    ap.add_argument("--qc", default=QC_DEFAULT)
    ap.add_argument("--detail", default=DETAIL_DEFAULT)
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    ap.add_argument("--repo-results", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "results"))
    args = ap.parse_args()

    qc_path = os.path.expanduser(args.qc)
    detail_path = os.path.expanduser(args.detail)
    for p in list(args.cohorts) + [qc_path, detail_path]:
        if not os.path.exists(p):
            die(f"not found: {p}")

    # ---------------- cohorts ----------------
    cohort = pd.concat([pd.read_csv(p, sep="\t", dtype=str) for p in args.cohorts],
                       ignore_index=True)
    if not {"research_id", "ancestry"} <= set(cohort.columns):
        die(f"cohort files need research_id + ancestry. Got: {list(cohort.columns)}")
    cohort["ancestry"] = cohort["ancestry"].astype(str).str.strip().str.upper()
    before = len(cohort)
    cohort = cohort.drop_duplicates(subset="research_id")
    if len(cohort) < before:
        print(f"NOTE: dropped {before - len(cohort)} duplicate research_ids across the "
              f"{len(args.cohorts)} cohort files.", file=sys.stderr)

    # ---------------- TRUST4 results ----------------
    detail = pd.read_csv(detail_path, sep="\t", dtype=str)
    if not {"research_id", "status", "n_cdr3"} <= set(detail.columns):
        die(f"detail file needs research_id, status, n_cdr3. Got: {list(detail.columns)}")
    detail = detail[detail["status"] == "ok"].copy()
    detail["n_cdr3"] = pd.to_numeric(detail["n_cdr3"], errors="coerce")

    # ---------------- QC metrics ----------------
    print(f"Reading QC metrics: {qc_path}", file=sys.stderr)
    qc = pd.read_csv(qc_path, sep="\t", dtype=str, compression="gzip")
    if "sample_id" not in qc.columns:
        die(f"QC file has no sample_id column. Got: {list(qc.columns)[:10]}")
    print(f"  {len(qc):,} samples x {len(qc.columns)} metrics", file=sys.stderr)

    missing = [m for m in METRICS if m not in qc.columns]
    if missing:
        die(f"QC file is missing expected metric column(s): {missing}\n"
            f"Full column list:\n  " + "\n  ".join(qc.columns))

    # ---------------- resolve the join key ----------------
    # The QC file keys on sample_id; the RNA-seq manifest carries BOTH sampleid and
    # research_id. Try the direct match first, fall back to routing through the manifest.
    # Never assume -- print which path was taken.
    cohort_ids = set(cohort["research_id"])
    direct_hits = len(cohort_ids & set(qc["sample_id"]))
    if direct_hits >= 0.5 * len(cohort_ids):
        print(f"Join: QC sample_id matches research_id directly "
              f"({direct_hits}/{len(cohort_ids)} of cohort).", file=sys.stderr)
        qc = qc.rename(columns={"sample_id": "research_id"})
    else:
        man_path = os.path.join(os.path.expanduser(args.mount), RNASEQ_MANIFEST)
        if not os.path.exists(man_path):
            die(f"QC sample_id does not look like research_id ({direct_hits} direct hits) "
                f"and the manifest is not reachable at {man_path} to translate. Remount?")
        print(f"Join: QC sample_id is NOT research_id ({direct_hits} direct hits) -- "
              f"routing through the manifest.", file=sys.stderr)
        man = pd.read_csv(man_path, sep="\t", dtype=str)
        if not {"sampleid", "research_id"} <= set(man.columns):
            die(f"manifest needs sampleid + research_id. Got: {list(man.columns)}")
        qc = qc.merge(man[["sampleid", "research_id"]],
                      left_on="sample_id", right_on="sampleid", how="inner")
        if qc.empty:
            die("manifest join produced zero rows -- sample_id matches neither "
                "research_id nor sampleid. Print a few of each and compare formats.")
        print(f"  translated {len(qc):,} QC rows to research_id.", file=sys.stderr)

    for m in METRICS:
        qc[m] = pd.to_numeric(qc[m], errors="coerce")
    for c in ("Non-Globin Reads", "Mapped Reads", "Read Length"):
        if c in qc.columns:
            qc[c] = pd.to_numeric(qc[c], errors="coerce")

    # Fraction of mapped reads that ARE globin -- how much of the budget the depletion
    # chemistry failed to reclaim. Higher = worse depletion.
    if {"Non-Globin Reads", "Mapped Reads"} <= set(qc.columns):
        qc["globin_fraction"] = 1 - (qc["Non-Globin Reads"] / qc["Mapped Reads"])
    else:
        qc["globin_fraction"] = pd.NA

    # ---------------- Read Length: the mechanistic hypothesis, checked first ----------
    print("\n" + "=" * 68)
    print("READ LENGTH -- does TRUST4's extraction stringency differ between people?")
    print("=" * 68)
    if "Read Length" in qc.columns:
        rl = qc["Read Length"].dropna()
        print(f"  distinct values across all {len(rl):,} samples: "
              f"{sorted(rl.unique().tolist())[:10]}")
        if rl.nunique() == 1:
            v = int(rl.iloc[0])
            print(f"  CONSTANT at {v} bp -> TRUST4's unmapped-read k-mer threshold is "
                  f"{v // 5} bp for everyone.")
            print("  Hypothesis DEAD: extraction stringency is identical across samples "
                  "and cannot contribute to the ancestry gap.")
        else:
            print("  *** READ LENGTH VARIES ACROSS SAMPLES ***")
            print("  TRUST4 sets its unmapped-read k-mer match requirement to readLen/5 "
                  "(capped at 101). Different read lengths therefore mean literally "
                  "different extraction stringency per person -- a direct, mechanistic "
                  "recovery confound. Check whether it varies BY ANCESTRY below.")
    else:
        print("  'Read Length' column absent -- skipped.")

    # ---------------- join everything ----------------
    keep = ["research_id"] + METRICS + DERIVED + \
           (["Read Length"] if "Read Length" in qc.columns else [])
    merged = cohort.merge(detail[["research_id", "n_cdr3"]], on="research_id", how="inner") \
                   .merge(qc[keep].drop_duplicates(subset="research_id"),
                          on="research_id", how="inner")
    print(f"\nJoined: {len(merged)} people with cohort + TRUST4 result + QC metrics")
    if merged.empty:
        die("zero overlap -- check that research_id formats match across the three inputs.")
    print("Per-ancestry n: " +
          ", ".join(f"{g}={int((merged['ancestry'] == g).sum())}" for g in ANCESTRY_GROUPS))

    os.makedirs(os.path.dirname(detail_path), exist_ok=True)
    jp = os.path.join(os.path.dirname(detail_path), "qc_confounds_detail.tsv")
    merged.to_csv(jp, sep="\t", index=False)
    print(f"Per-person joined detail (VM-local, NOT committed): {jp}")

    # ---------------- the two-part confound test ----------------
    all_metrics = METRICS + DERIVED + (["Read Length"] if "Read Length" in merged else [])
    rows = []
    for m in all_metrics:
        s = pd.to_numeric(merged[m], errors="coerce")
        if s.notna().sum() < 10 or s.nunique() < 2:
            continue
        r = s.corr(merged["n_cdr3"])
        gmeans = merged.assign(_m=s).groupby("ancestry")["_m"].mean().dropna()
        spread = (gmeans.max() / gmeans.min()) if (len(gmeans) > 1 and gmeans.min() > 0) \
            else float("nan")
        groups = [merged.loc[merged["ancestry"] == g, m].astype(float).dropna().values
                  for g in ANCESTRY_GROUPS if (merged["ancestry"] == g).sum() > 1]
        try:
            _, p_anc = kruskal(*groups) if len(groups) > 1 else (None, float("nan"))
        except ValueError:
            p_anc = float("nan")
        rows.append({
            "metric": m,
            "r_with_cdr3": round(r, 3) if pd.notna(r) else None,
            "ancestry_spread_ratio": round(spread, 3) if pd.notna(spread) else None,
            "p_varies_by_ancestry": round(p_anc, 4) if pd.notna(p_anc) else None,
            # Both conditions, explicitly. Neither alone is enough.
            "candidate_confound": bool(
                pd.notna(r) and abs(r) >= 0.2 and pd.notna(p_anc) and p_anc < 0.05),
        })

    res = pd.DataFrame(rows)
    res["_rank"] = res["r_with_cdr3"].abs().fillna(0)
    res = res.sort_values(["candidate_confound", "_rank"], ascending=[False, False]) \
             .drop(columns="_rank")

    os.makedirs(args.repo_results, exist_ok=True)
    out = os.path.join(args.repo_results, "rnaseq_qc_confounds_check.csv")
    res.to_csv(out, index=False)

    print("\n" + "=" * 68)
    print("CANDIDATE CONFOUNDS -- needs BOTH |r| >= 0.2 with recovery AND p < 0.05 across ancestry")
    print("=" * 68)
    print(res.to_string(index=False))
    print(f"\nDe-identified summary (safe to commit): {out}")

    flagged = res[res["candidate_confound"]]
    if flagged.empty:
        print("\nNOTHING FLAGGED. None of these technical metrics both predicts recovery "
              "and varies by ancestry. That does not prove the gap is biological -- it "
              "means the mundane explanations AoU measured for us are exhausted, and the "
              "remaining candidates are depth (already known, ~40%) and things nobody "
              "measured. Combined with p = 0.18 on the gap itself, the honest read stays "
              "'suggestive, underpowered, not explained away'.")
    else:
        print(f"\n{len(flagged)} METRIC(S) FLAGGED: "
              f"{', '.join(flagged['metric'])}. Each of these both predicts recovery and "
              f"differs across ancestry groups -- so each could manufacture the gap "
              f"without any biology. Next step is to re-run the ancestry comparison "
              f"controlling for them, exactly as check_depth_confound.py did for depth.")

    # Per-ancestry means for the flagged/leading metrics, so the direction is visible.
    lead = list(flagged["metric"]) if not flagged.empty else list(res["metric"].head(5))
    summ = merged.groupby("ancestry").agg(
        n_people=("research_id", "count"), mean_cdr3=("n_cdr3", "mean"),
        **{m: (m, "mean") for m in lead}).round(4)
    sp = os.path.join(args.repo_results, "rnaseq_qc_confounds_by_ancestry.csv")
    summ.to_csv(sp)
    print(f"\nPer-ancestry means for the leading metrics: {sp}")
    print(summ.to_string())


if __name__ == "__main__":
    main()
