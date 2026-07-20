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
      L1-penalized (LASSO) logistic regression. Two separate questions, not one:
        (B1) GENERALIZATION: 5-fold stratified cross-validation, out-of-fold predictions ->
             pooled AUROC/AUPRC, a risk-decile enrichment table, PPV/sensitivity/specificity/
             NPV at top-5%/top-10% thresholds, and the same broken out per ancestry (gated on
             a minimum case count -- never reported for a group too small to mean anything).
             This is the "does it actually predict on people it hasn't seen" question -- NOT
             answered by coefficient inspection alone (2026-07-20 correction: an earlier version
             of this check only reported in-sample coefficients, which shows what the model
             learned but says nothing about generalization or per-ancestry transfer).
        (B2) WHAT IT LEARNED: one final fit on all data, top +/- coefficients (as odds ratios).
             Does a model with ZERO disease-specific hints independently rediscover the known
             risk allele(s) as its top signal? Interpretation only -- not a performance claim.

Both checks use AoU-native calls -- the only callset with population-scale coverage (our
own SpecHLA/SpecImmune calling has only run on the n=60 pilot cohort) -- so this result is
itself further evidence about AoU-native's reliability (see reports/aou_callset_validation.md).

AGGREGATE-ONLY OUTPUT: contingency counts and fitted model coefficients only -- never a
participant-level row. Same egress discipline as every other script in this repo.

Inputs:
  --tsv               hla_genotypes.tsv (mounted path)
  --ancestry-tsv       AoU genetic-ancestry TSV (mounted path)
  --phenotype-csv      person_id,case (0/1) for every person in --tsv. Built directly in a
                       Workbench notebook, no separate template file: look up the disease's
                       standard SNOMED concept_id via Athena (athena.ohdsi.org, public, no
                       AoU login), then query condition_occurrence for
                       `SELECT DISTINCT person_id WHERE condition_concept_id = <id>` against
                       the CDR (project `wb-silky-artichoke-2408`, dataset `C2025Q4R6` --
                       no WORKSPACE_CDR env var in this Workbench, see ENVIRONMENT.md quirk #5),
                       mark case=1 for matches. Confirmed concept_ids in use: Celiac disease
                       = 194992 (SNOMED 396331005, maps from ICD-10CM K90.0); Narcolepsy with
                       cataplexy = 437854 (SNOMED 193042000, maps from ICD-10CM G47.411 --
                       deliberately distinct from "without cataplexy", concept 43531721).
  --disease            celiac | narcolepsy

Usage (via `pixi run -e spechla`):
  python3 hla_disease_sanity_check.py --tsv <hla_genotypes.tsv> --ancestry-tsv <ancestry_preds.tsv> \
      --phenotype-csv <celiac_cases.csv> --disease celiac
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import fisher_exact
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

CLASSICAL_GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
MIN_CASES_FOR_GROUP_METRIC = 10  # below this, a per-ancestry AUROC/PPV is noise, not a number

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


def build_feature_matrix(df, anc_col):
    """Sparse construction, not a dense DataFrame: at real cohort scale (535k people x ~2,373
    distinct alleles, per experiment_a3) a dense float64 matrix is ~10GB, and StandardScaler
    would make a second full copy on top of it -- enough to OOM-kill a 25GB VM (confirmed
    2026-07-20, "Killed" with no traceback = the kernel OOM killer, not a Python exception).
    Each person carries at most 2 non-zero alleles per gene, so this matrix is >99% zeros --
    build it as a scipy.sparse CSR matrix instead, which costs memory proportional to actual
    non-zero entries, not rows x columns.
    Returns (X, y, ancestry_arr, feature_names) restricted to rows callable at all 8 genes
    with a non-null case label."""
    n = len(df)
    row_chunks, col_chunks = [], []
    feature_names = []
    col_offset = 0
    callable_mask = pd.Series(True, index=df.index)
    positions = np.arange(n)

    for gene in CLASSICAL_GENES:
        c1, c2 = f"{gene}_1", f"{gene}_2"
        a1 = df[c1].map(to2field)
        a2 = df[c2].map(to2field)
        callable_mask &= a1.notna() & a2.notna()

        alleles = sorted(pd.concat([a1, a2]).dropna().unique())
        allele_to_col = {al: col_offset + i for i, al in enumerate(alleles)}
        feature_names.extend(f"{gene}*{al}" for al in alleles)

        for a in (a1, a2):
            valid = a.notna().values
            row_chunks.append(positions[valid])
            col_chunks.append(a[valid].map(allele_to_col).to_numpy())

        col_offset += len(alleles)

    row_idx = np.concatenate(row_chunks)
    col_idx = np.concatenate(col_chunks)
    data = np.ones(len(row_idx), dtype=np.float32)
    # COO -> CSR summing duplicate (row, col) entries is exactly what gives dosage=2 for a
    # homozygous person (both haplotype columns map to the same allele -> same coordinate twice).
    X_alleles = sparse.coo_matrix((data, (row_idx, col_idx)), shape=(n, col_offset)).tocsr()

    if anc_col:
        anc_dummies = pd.get_dummies(df[anc_col], prefix="ancestry", dummy_na=False)
        feature_names.extend(anc_dummies.columns)
        X_anc = sparse.csr_matrix(anc_dummies.to_numpy(dtype=np.float32))
        X = sparse.hstack([X_alleles, X_anc], format="csr")
    else:
        X = X_alleles

    mask = (callable_mask & df["case"].notna()).to_numpy()
    idx = np.where(mask)[0]
    X = X[idx]
    y = df["case"].iloc[idx].astype(int).to_numpy()
    ancestry_arr = df[anc_col].iloc[idx].astype(str).to_numpy() if anc_col else None
    return X, y, ancestry_arr, feature_names


def cross_validated_oof_predictions(X, y, n_splits=5, seed=42):
    """Out-of-fold predicted probabilities: every person is scored by a model that never saw
    them during training. Scaling is fit on the training fold only and applied to the held-out
    fold -- fitting the scaler on all data first (as the interpretation-only fit below does) would
    leak test-fold statistics into training and quietly inflate the reported performance."""
    oof = np.zeros(len(y), dtype=np.float64)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for train_idx, test_idx in skf.split(X, y):
        scaler = StandardScaler(with_mean=False)
        X_train = scaler.fit_transform(X[train_idx])
        X_test = scaler.transform(X[test_idx])
        model = LogisticRegression(penalty="l1", solver="liblinear", C=0.1, max_iter=2000)
        model.fit(X_train, y[train_idx])
        oof[test_idx] = model.predict_proba(X_test)[:, 1]
    return oof


def threshold_metrics(y, proba, pct):
    """PPV/sensitivity/specificity/NPV for "flag the top pct fraction by predicted risk"."""
    n = len(y)
    k = max(1, int(round(n * pct)))
    order = np.argsort(-proba)
    flag = np.zeros(n, dtype=bool)
    flag[order[:k]] = True
    tp = int((flag & (y == 1)).sum())
    fp = int((flag & (y == 0)).sum())
    fn = int((~flag & (y == 1)).sum())
    tn = int((~flag & (y == 0)).sum())
    ppv = tp / (tp + fp) if (tp + fp) else float("nan")
    sens = tp / (tp + fn) if (tp + fn) else float("nan")
    spec = tn / (tn + fp) if (tn + fp) else float("nan")
    npv = tn / (tn + fn) if (tn + fn) else float("nan")
    return k, ppv, sens, spec, npv


def print_generalization_check(y, oof_proba, ancestry_arr):
    print("### (B1) Generalization -- 5-fold cross-validated, out-of-fold predictions\n")
    print("Every number below comes from predictions on people NOT used to train that "
          "prediction -- this is the actual answer to \"does it generalize,\" not the "
          "coefficient inspection in (B2) below.\n")

    prevalence = y.mean()
    auroc = roc_auc_score(y, oof_proba)
    auprc = average_precision_score(y, oof_proba)
    print(f"Pooled AUROC: {auroc:.3f} | Pooled AUPRC: {auprc:.3f} "
          f"(baseline/prevalence = {100*prevalence:.3f}%) -- AUPRC matters more than AUROC "
          f"here given how rare the outcome is; plain accuracy is not reported because "
          f"predicting \"nobody has it\" already scores {100*(1-prevalence):.1f}%.\n")

    print("**Risk-decile enrichment** (rank-based bins -- decile 1 = highest predicted risk)\n")
    print("| Decile | n | observed case rate | enrichment vs baseline |")
    print("|---|---|---|---|")
    order = np.argsort(-oof_proba)
    for i, bin_idx in enumerate(np.array_split(order, 10), start=1):
        rate = y[bin_idx].mean()
        enrich = rate / prevalence if prevalence > 0 else float("nan")
        print(f"| {i} | {len(bin_idx)} | {100*rate:.3f}% | {enrich:.1f}x |")
    print()

    print("**Threshold metrics** (flagging the top X% by predicted risk as \"high-risk\")\n")
    print("| Threshold | n flagged | PPV | Sensitivity | Specificity | NPV |")
    print("|---|---|---|---|---|---|")
    for pct in (0.05, 0.10):
        k, ppv, sens, spec, npv = threshold_metrics(y, oof_proba, pct)
        print(f"| top {int(pct*100)}% | {k} | {100*ppv:.2f}% | {100*sens:.2f}% | "
              f"{100*spec:.2f}% | {100*npv:.2f}% |")
    print()

    if ancestry_arr is not None:
        print("**Per-ancestry generalization** (same OOF predictions, evaluated within each "
              f"group; groups with fewer than {MIN_CASES_FOR_GROUP_METRIC} cases or "
              f"{MIN_CASES_FOR_GROUP_METRIC} controls are flagged as insufficient rather than "
              "given a number that would just be noise)\n")
        print("| Ancestry | n | cases | AUROC | AUPRC | top-10%-within-group PPV | sensitivity |")
        print("|---|---|---|---|---|---|---|")
        for label in sorted(pd.unique(ancestry_arr)):
            m = ancestry_arr == label
            y_g, p_g = y[m], oof_proba[m]
            n_g, cases_g = len(y_g), int(y_g.sum())
            if cases_g < MIN_CASES_FOR_GROUP_METRIC or (n_g - cases_g) < MIN_CASES_FOR_GROUP_METRIC:
                print(f"| {label} | {n_g} | {cases_g} | insufficient cases | -- | -- | -- |")
                continue
            auroc_g = roc_auc_score(y_g, p_g)
            auprc_g = average_precision_score(y_g, p_g)
            _, ppv_g, sens_g, _, _ = threshold_metrics(y_g, p_g, 0.10)
            print(f"| {label} | {n_g} | {cases_g} | {auroc_g:.3f} | {auprc_g:.3f} | "
                  f"{100*ppv_g:.2f}% | {100*sens_g:.2f}% |")
        print()


def print_learned_coefficients(X, y, feature_names, top_n=15):
    print("### (B2) What the model learned -- one fit on all data, for interpretation only\n")
    print("Not a performance claim (see (B1) above for that) -- this just asks whether the "
          "model, given zero disease-specific hints, lands on real biology.\n")

    Xs = StandardScaler(with_mean=False).fit_transform(X)
    model = LogisticRegression(penalty="l1", solver="liblinear", C=0.1, max_iter=2000)
    model.fit(Xs, y)

    coefs = pd.Series(model.coef_[0], index=feature_names).sort_values()
    nonzero = coefs[coefs != 0]
    print(f"Non-zero coefficients after L1 shrinkage: {len(nonzero)} / {len(coefs)}\n")

    print(f"#### Top {top_n} risk-associated features (positive coefficient = higher odds)\n")
    print("| Feature | Coefficient | Odds ratio (approx, per copy) |")
    print("|---|---|---|")
    for feat, coef in coefs.sort_values(ascending=False).head(top_n).items():
        print(f"| {feat} | {coef:.3f} | {np.exp(coef):.2f} |")
    print()

    print(f"#### Top {top_n} protective-associated features (negative coefficient)\n")
    print("| Feature | Coefficient | Odds ratio (approx, per copy) |")
    print("|---|---|---|")
    for feat, coef in coefs.sort_values().head(top_n).items():
        print(f"| {feat} | {coef:.3f} | {np.exp(coef):.2f} |")
    print()

    return coefs


def untargeted_ml_check(df, anc_col, top_n=15):
    print("## (B) Untargeted check -- L1 logistic regression over ALL alleles, all 8 genes\n")
    print("No disease-specific hints given to the model -- features are every 2-field allele "
          "seen at each of the 8 classical genes (dosage 0/1/2), plus ancestry dummies for "
          "population-stratification adjustment.\n")

    X, y, ancestry_arr, feature_names = build_feature_matrix(df, anc_col)

    print(f"Rows used: {X.shape[0]} ({int(y.sum())} cases / {int((1-y).sum())} controls). "
          f"Features: {X.shape[1]} (allele dosages + ancestry dummies).\n")

    if len(np.unique(y)) < 2 or int(y.sum()) < 10:
        print("SKIPPED -- fewer than 10 cases after filtering; not enough signal for a "
              "regression. Report the case count above and reconsider phenotype extraction.\n")
        return None

    oof_proba = cross_validated_oof_predictions(X, y)
    print_generalization_check(y, oof_proba, ancestry_arr)
    return print_learned_coefficients(X, y, feature_names, top_n=top_n)


class Tee:
    """Duplicates writes to multiple streams -- lets every print() reach both the terminal
    (so you can watch it live) and a results file (so the actual record doesn't depend on
    anyone's memory of a chat transcript)."""
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)

    def flush(self):
        for s in self.streams:
            s.flush()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--ancestry-tsv", required=True)
    ap.add_argument("--phenotype-csv", required=True,
                    help="CSV with person_id,case (0/1) for every person in --tsv")
    ap.add_argument("--disease", required=True, choices=list(RISK_MARKERS.keys()))
    ap.add_argument("--id-col", default="research_id")
    ap.add_argument("--out-dir", default=os.path.expanduser("~/pipeline_outputs/disease_sanity_check"),
                    help="Where to write the results file + full coefficient CSV. "
                         "Aggregate-only output, safe to commit into the repo afterward.")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

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

    results_path = os.path.join(args.out_dir, f"{args.disease}_results.md")
    real_stdout = sys.stdout
    with open(results_path, "w") as f:
        sys.stdout = Tee(real_stdout, f)
        try:
            print(f"# HLA x {args.disease} sanity check\n")
            print(f"Loaded {n_total} people from the HLA-calls TSV; {n_labeled} have a "
                  f"phenotype label ({int(df['case'].sum())} cases, "
                  f"{int(n_labeled - df['case'].sum())} controls).\n")

            positive_control_check(df, args.disease, anc_col)
            coefs = untargeted_ml_check(df, anc_col)
        finally:
            sys.stdout = real_stdout

    print(f"\n(full results written to {results_path})", file=sys.stderr)

    if coefs is not None:
        coefs_path = os.path.join(args.out_dir, f"{args.disease}_full_coefficients.csv")
        coefs.rename("coefficient").rename_axis("feature").reset_index().to_csv(
            coefs_path, index=False)
        print(f"(full coefficient table -- all {len(coefs)} features, not just the top "
              f"15 -- written to {coefs_path})", file=sys.stderr)


if __name__ == "__main__":
    main()
