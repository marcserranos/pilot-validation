#!/usr/bin/env python3
"""Experiment D -- resolution-tiered SR-vs-truth comparison (AoU-native and SpecHLA vs SpecImmune).

analyze_experiment_d.py answers "do all 3 methods say the same genotype" (binary exact match at
2-field). It cannot distinguish a near-miss (right allele family/group, wrong sub-type -- often a
resolution or reference-panel artifact) from a real miss (completely different allele group) --
and that distinction is exactly what determines whether a discordance matters clinically.

This script treats SpecImmune-LR as ground truth (per request) and, for each of AoU-native and
SpecHLA (both short-read), classifies every per-person-per-gene genotype call into one of 4 tiers
by aligning the two called alleles against the two truth alleles (best of the 2 possible pairings,
maximizing exact matches then family matches) and counting, among the 2 aligned allele-pairs:

  EXACT        both alleles match the truth at 2-field (e.g. 01:01 == 01:01)
  CLOSE        no allele is a different group, but not both exact (e.g. one 01:02 vs true 01:01 --
               same first-field group, wrong sub-type)
  PARTIAL_MISS exactly one allele is from a completely different group than its aligned truth allele
  FULL_MISS    both alleles are from different groups than truth (SpecImmune and the SR method are
               not just mis-resolved, they called different things)

Reads the same per-person comparison_log.csv files as analyze_experiment_d.py (run_label ==
"experiment_d"), only among people/genes where all 3 methods gave a callable 2-field genotype
(same denominator as the existing report, so the two are directly comparable).

Usage (via `pixi run -e spechla` for pandas):
  python3 analyze_experiment_d_resolution.py ~/pipeline_outputs/experiment_d/cohort.tsv
"""
import argparse
import os
import sys

import pandas as pd

GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
NULL = {"", "NA", "-", "nan", "None", ".", "na"}


def norm2(allele):
    if allele is None:
        return None
    a = str(allele).strip()
    if a in NULL:
        return None
    if "*" in a:
        a = a.split("*", 1)[1]
    fields = a.split(":")
    if len(fields) < 2:
        return None
    return f"{fields[0]}:{fields[1]}"


def group(allele2):
    """First field of an already-normalized 2-field allele, e.g. '01:01' -> '01'."""
    return allele2.split(":", 1)[0]


def geno(pair):
    n1, n2 = norm2(pair[0]), norm2(pair[1])
    if n1 is None or n2 is None:
        return None
    return (n1, n2)


def classify(sr, truth):
    """sr, truth: 2-tuples of normalized 2-field alleles (order not yet resolved).
    Returns one of EXACT / CLOSE / PARTIAL_MISS / FULL_MISS, aligning sr<->truth by the pairing
    that maximizes exact matches, then family (group) matches."""
    pairings = [(sr[0], truth[0], sr[1], truth[1]), (sr[0], truth[1], sr[1], truth[0])]

    def score(p):
        a1, t1, a2, t2 = p
        exact = (a1 == t1) + (a2 == t2)
        fam = (group(a1) == group(t1)) + (group(a2) == group(t2))
        return (exact, fam)

    best = max(pairings, key=score)
    a1, t1, a2, t2 = best
    kinds = []
    for a, t in ((a1, t1), (a2, t2)):
        if a == t:
            kinds.append("EXACT")
        elif group(a) == group(t):
            kinds.append("FAMILY")
        else:
            kinds.append("DIFFERENT")
    n_exact = kinds.count("EXACT")
    n_diff = kinds.count("DIFFERENT")
    if n_exact == 2:
        return "EXACT"
    if n_diff == 0:
        return "CLOSE"
    if n_diff == 1:
        return "PARTIAL_MISS"
    return "FULL_MISS"


def load_rows(cohort, outroot):
    anc = dict(zip(cohort["person_id"].astype(str), cohort["ancestry"].astype(str)))
    frames, missing = [], []
    for pid in cohort["person_id"].astype(str):
        csv = os.path.join(outroot, pid, "comparison_log.csv")
        if not os.path.exists(csv):
            missing.append(pid)
            continue
        df = pd.read_csv(csv, dtype=str)
        df = df[df["run_label"] == "experiment_d"]
        if df.empty:
            missing.append(pid)
            continue
        frames.append(df)
    if not frames:
        print("FATAL: no experiment_d rows found yet.", file=sys.stderr)
        sys.exit(1)
    allrows = pd.concat(frames, ignore_index=True)
    allrows = allrows.sort_values("timestamp").drop_duplicates(["person_id", "gene"], keep="last")
    allrows["ancestry"] = allrows["person_id"].astype(str).map(anc)
    if missing:
        print(f"NOTE: {len(missing)} cohort people have no experiment_d result yet: {missing}",
              file=sys.stderr)
    return allrows


TIERS = ["EXACT", "CLOSE", "PARTIAL_MISS", "FULL_MISS"]


def analyze(allrows):
    recs = []
    for gene in GENES:
        g = allrows[allrows["gene"] == gene]
        counts = {"aou": {t: 0 for t in TIERS}, "sh": {t: 0 for t in TIERS}}
        callable_n = 0
        for _, r in g.iterrows():
            aou = geno((r.get("aou_1"), r.get("aou_2")))
            sh = geno((r.get("spechla_1"), r.get("spechla_2")))
            si = geno((r.get("specimmune_1"), r.get("specimmune_2")))
            if aou is None or sh is None or si is None:
                continue
            callable_n += 1
            counts["aou"][classify(aou, si)] += 1
            counts["sh"][classify(sh, si)] += 1
        rec = {"gene": gene, "callable_n": callable_n}
        for method in ("aou", "sh"):
            for t in TIERS:
                rec[f"{method}_{t}"] = counts[method][t]
        recs.append(rec)
    return pd.DataFrame(recs)


def pct(n, d):
    return f"{100*n/d:.0f}%" if d else "—"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cohort")
    ap.add_argument("--outroot", default=os.path.expanduser("~/pipeline_outputs"))
    args = ap.parse_args()

    cohort = pd.read_csv(args.cohort, sep="\t", dtype=str)
    allrows = load_rows(cohort, args.outroot)
    df = analyze(allrows)

    print("# Experiment D -- resolution-tiered SR-vs-truth (SpecImmune-LR = ground truth)\n")
    print("Among people/genes where all 3 methods gave a callable 2-field genotype "
          "(same denominator as analyze_experiment_d.py's tables).\n")
    print("| Gene | n | AoU exact | AoU close | AoU partial-miss | AoU full-miss | "
          "SpecHLA exact | SpecHLA close | SpecHLA partial-miss | SpecHLA full-miss |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    totals = {f"{m}_{t}": 0 for m in ("aou", "sh") for t in TIERS}
    total_n = 0
    for _, r in df.iterrows():
        n = r["callable_n"]
        total_n += n
        cells = [r["gene"], str(n)]
        for m in ("aou", "sh"):
            for t in TIERS:
                totals[f"{m}_{t}"] += r[f"{m}_{t}"]
                cells.append(f"{r[f'{m}_{t}']} ({pct(r[f'{m}_{t}'], n)})")
        print("| " + " | ".join(cells) + " |")

    cells = ["**ALL**", str(total_n)]
    for m in ("aou", "sh"):
        for t in TIERS:
            cells.append(f"{totals[f'{m}_{t}']} ({pct(totals[f'{m}_{t}'], total_n)})")
    print("| " + " | ".join(cells) + " |")

    df.to_csv(os.path.join(os.path.dirname(args.cohort), "experiment_d_resolution_counts.csv")
              if os.path.dirname(args.cohort) else "experiment_d_resolution_counts.csv",
              index=False)
    print("\n(aggregate counts only -- paste this table back)", file=sys.stderr)


if __name__ == "__main__":
    main()
