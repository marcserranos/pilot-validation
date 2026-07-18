#!/usr/bin/env python3
"""HLA x disease sanity-check pilot (Celiac / Narcolepsy) on AoU-native calls.

First-principles check: does a KNOWN, huge-effect HLA-disease association show up using
nothing but AoU-native's pre-computed HLA calls + EHR-derived case/control status? Two
independent checks, run together so they cross-validate each other:

  (A) POSITIVE-CONTROL CHECK: hand-curated carrier flag for the disease's known risk marker
      (Celiac: DQ2.5 = DQA1*05:01 + DQB1*02:01 (unphased -- presence, not cis/trans, per
       standard clinical-typing convention) or DQ8 = DQA1*03:01 + DQB1*03:02.
       Narcolepsy: DQB1*06:02 dosage 0/1/2.)
      vs case/control -> 2x2 table, odds ratio, Fisher's exact p, pooled + per-ancestry.
      This answers "does the expected signal exist at all."

  (B) UNTARGETED CHECK: every 2-field allele at all 8 classical genes, dosage-encoded
      (0/1/2 copies), PLUS ancestry dummies (population-stratification adjustment) ->
      L1-penalized (LASSO) logistic regression -> top +/- coefficients (as odds ratios).
      Answers the stronger question: does a model with ZERO disease-specific hints
      independently rediscover the known risk allele(s) as its top signal?

Both checks use AoU-native calls -- the only callset with population-scale coverage (our
own SpecHLA/SpecImmune calling has only run on the n=60 pilot cohort) -- so this result is
itself further evidence about AoU-native's reliability (see reports/aou_callset_validation.md).

AGGREGATE-ONLY OUTPUT: contingency counts and fitted model coefficients only -- never a
participant-level row. Same egress discipline as every other script in this repo.

Inputs:
  --tsv               hla_genotypes.tsv (mounted path)
  --ancestry-tsv       AoU genetic-ancestry TSV (mounted path)
  --phenotype-csv      person_id,case  (0/1) -- built by you from a Workbench BigQuery query
                       against condition_occurrence (see scripts/build_phenotype_csv_template.sql
                       for the query shape; concept_id must be looked up via Athena first).
  --disease            celiac | narcolepsy

Usage (via `pixi run -e spechla`):
  python3 hla_disease_sanity_check.py --tsv <hla_genotypes.tsv> --ancestry-tsv <ancestry_preds.tsv> \
      --phenotype-csv <celiac_cases.csv> --disease celiac
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

CLASSICAL_GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]

# Known risk markers per disease -- the positive-control definition (2-field, unphased).
RISK_MARKERS = {
    "celiac": {
        "description": "DQ2.5 (DQA1*05:01+DQB1*02:01) or DQ8 (DQA1*03:01+DQB1*03:02) carrier, "
                        "unphased presence (standard clinical-typing convention -- both cis- and "
                        "trans-configured DQ2.5 confer risk, and AoU-native calls are not phased "
                        "across genes anyway)",
        "components": [("DQA1", "05:01"), ("DQB1", "02:01"), ("DQA1", "03:01"), ("DQB1", "03:02")],
    },
    "narcolepsy": {
        "description": "DQB1*06:02 dosage (0/1/2 copies)",
        "components": [("DQB1", "06:02")],
    },
}


def to2field(allele):
    if pd.isna(allele):
        return None
    a = str(allele).strip()
    if a == "" or a.lower() in {"na", "nan", "none", "-", "."}:
        return None
    if "*" in a:
        a = a.split("*", 1)[1]
    fields = a.split(":")
    if len(fields) < 2:
        return None
    return f"{fields[0]}:{fields[1]}"


def dosage(df, gene, allele2):
    """Count of copies of a specific 2-field allele at a gene, per person. NaN if uncallable."""
    c1, c2 = f"{gene}_1", f"{gene}_2"
    a1 = df[c1].map(to2field)
    a2 = df[c2].map(to2field)
    callable_mask = a1.notna() & a2.notna()
    d = ((a1 == allele2).astype(float) + (a2 == allele2).astype(float))
    d[~callable_mask] = np.nan
    return d


def positive_control_check(df, disease, anc_col):
    print(f"## (A) Positive-control check -- {disease}\n")
    print(f"Risk marker: {RISK_MARKERS[disease]['description']}\n")

    if disease == "celiac":
        dq25 = (dosage(df, "DQA1", "05:01") > 0) & (dosage(df, "DQB1", "02:01") > 0)
        dq8 = (dosage(df, "DQA1", "03:01") > 0) & (dosage(df, "DQB1", "03:02") > 0)
        carrier = (dq25.fillna(False) | dq8.fillna(False))
        callable_mask = dq25.notna() & dq8.notna()
    else:  # narcolepsy
        dose = dosage(df, "DQB1", "06:02")
        carrier = dose.fillna(0) > 0
        callable_mask = dose.notna()

    sub = df.loc[callable_mask].copy()
    sub["_carrier"] = carrier.loc[callable_mask]

    def table_and_or(g, label):
        n = len(g)
        if n == 0 or g["case"].nunique() < 2 or g["_carrier"].nunique() < 2:
            print(f"| {label} | {n} | -- | -- | -- | -- | skipped (degenerate) |")
            return
        a = ((g["_carrier"]) & (g["case"] == 1)).sum()
        b = ((g["_carrier"]) & (g["case"] == 0)).sum()
        c = ((~g["_carrier"]) & (g["case"] == 1)).sum()
        d = ((~g["_carrier"]) & (g["case"] == 0)).sum()
        odds_ratio, p = fisher_exact([[a, b], [c, d]])
        print(f"| {label} | {n} | {a} | {b} | {c} | {d} | {odds_ratio:.2f} | {p:.3e} |")

    print("| Group | n | carrier+case | carrier+ctrl | non-carrier+case | non-carrier+ctrl "
          "| OR | Fisher p |")
    print("|---|---|---|---|---|---|---|---|")
    table_and_or(sub, "ALL (pooled)")
    if anc_col:
        for label, g in sub.groupby(anc_col):
            table_and_or(g, label)
    print()


def untargeted_ml_check(df, anc_col, top_n=15):
    print("## (B) Untargeted check -- L1 logistic regression over ALL alleles, all 8 genes\n")
    print("No disease-specific hints given to the model -- features are every 2-field allele "
          "seen at each of the 8 classical genes (dosage 0/1/2), plus ancestry dummies for "
          "population-stratification adjustment. If the known risk allele(s) surface near the "
          "top of the coefficient ranking below, that's independent confirmation the calling + "
          "join + phenotype pipeline is sound -- not just a re-check of something we told it to "
          "find.\n")

    feature_cols = {}
    for gene in CLASSICAL_GENES:
        c1, c2 = f"{gene}_1", f"{gene}_2"
        a1 = df[c1].map(to2field)
        a2 = df[c2].map(to2field)
        alleles = pd.concat([a1, a2]).dropna().unique()
        for al in alleles:
            colname = f"{gene}*{al}"
            feature_cols[colname] = ((a1 == al).astype(float) + (a2 == al).astype(float))

    X = pd.DataFrame(feature_cols, index=df.index)
    callable_mask = pd.Series(True, index=df.index)
    for gene in CLASSICAL_GENES:
        c1, c2 = f"{gene}_1", f"{gene}_2"
        callable_mask &= df[c1].map(to2field).notna() & df[c2].map(to2field).notna()

    if anc_col:
        anc_dummies = pd.get_dummies(df[anc_col], prefix="ancestry")
        X = pd.concat([X, anc_dummies], axis=1)

    mask = callable_mask & df["case"].notna()
    X, y = X.loc[mask].fillna(0), df.loc[mask, "case"].astype(int)

    print(f"Rows used: {len(X)} ({int(y.sum())} cases / {int((1-y).sum())} controls). "
          f"Features: {X.shape[1]} (allele dosages + ancestry dummies).\n")

    if y.nunique() < 2 or int(y.sum()) < 10:
        print("SKIPPED -- fewer than 10 cases after filtering; not enough signal for a "
              "regression. Report the case count above and reconsider phenotype extraction.\n")
        return

    Xs = StandardScaler(with_mean=False).fit_transform(X)  # sparse-friendly, dosage stays >=0
    model = LogisticRegression(penalty="l1", solver="liblinear", C=0.1, max_iter=2000)
    model.fit(Xs, y)

    coefs = pd.Series(model.coef_[0], index=X.columns).sort_values()
    nonzero = coefs[coefs != 0]
    print(f"Non-zero coefficients after L1 shrinkage: {len(nonzero)} / {len(coefs)}\n")

    print(f"### Top {top_n} risk-associated features (positive coefficient = higher odds)\n")
    print("| Feature | Coefficient | Odds ratio (approx, per copy) |")
    print("|---|---|---|")
    for feat, coef in coefs.sort_values(ascending=False).head(top_n).items():
        print(f"| {feat} | {coef:.3f} | {np.exp(coef):.2f} |")
    print()

    print(f"### Top {top_n} protective-associated features (negative coefficient)\n")
    print("| Feature | Coefficient | Odds ratio (approx, per copy) |")
    print("|---|---|---|")
    for feat, coef in coefs.sort_values().head(top_n).items():
        print(f"| {feat} | {coef:.3f} | {np.exp(coef):.2f} |")
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--ancestry-tsv", required=True)
    ap.add_argument("--phenotype-csv", required=True,
                    help="CSV with person_id,case (0/1) for every person in --tsv")
    ap.add_argument("--disease", required=True, choices=list(RISK_MARKERS.keys()))
    ap.add_argument("--id-col", default="research_id")
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep="\t", dtype=str)
    anc = pd.read_csv(args.ancestry_tsv, sep=None, engine="python", dtype=str)
    anc_id = next((c for c in anc.columns if c.lower() in (args.id_col, "research_id", "person_id", "s")), None)
    anc_col_name = next((c for c in anc.columns if "ancestry" in c.lower()), None)
    if anc_id is None or anc_col_name is None:
        print(f"WARNING: could not identify ancestry columns in {list(anc.columns)}; "
              f"proceeding without ancestry adjustment.", file=sys.stderr)
        anc_col = None
    else:
        df = df.merge(anc[[anc_id, anc_col_name]].rename(columns={anc_col_name: "_ancestry"}),
                      left_on=args.id_col, right_on=anc_id, how="left")
        anc_col = "_ancestry"

    pheno = pd.read_csv(args.phenotype_csv, dtype=str)
    if args.id_col not in pheno.columns:
        # tolerate a raw "person_id" phenotype file joined against a research_id HLA table
        id_candidates = [c for c in pheno.columns if c.lower() in ("person_id", "research_id")]
        if not id_candidates:
            print(f"FATAL: --phenotype-csv has no id column matching {args.id_col}/person_id/"
                  f"research_id. Columns found: {list(pheno.columns)}", file=sys.stderr)
            sys.exit(1)
        pheno = pheno.rename(columns={id_candidates[0]: args.id_col})
    pheno[args.id_col] = pheno[args.id_col].astype(str)
    pheno["case"] = pheno["case"].astype(int)
    df = df.merge(pheno[[args.id_col, "case"]], on=args.id_col, how="left")

    n_total = len(df)
    n_labeled = df["case"].notna().sum()
    print(f"# HLA x {args.disease} sanity check\n")
    print(f"Loaded {n_total} people from the HLA-calls TSV; {n_labeled} have a phenotype label "
          f"({int(df['case'].sum())} cases, {int(n_labeled - df['case'].sum())} controls).\n")

    positive_control_check(df, args.disease, anc_col)
    untargeted_ml_check(df, anc_col)


if __name__ == "__main__":
    main()
