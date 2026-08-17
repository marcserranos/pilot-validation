#!/usr/bin/env python3
"""Disease burden + rankings for the LR x RNA-seq overlap cohort (Task 2, step 3).

Pure pandas over the VM-local files written by query_overlap_phenotypes.py. No BigQuery,
no mount, no BAMs. Produces DE-IDENTIFIED, small-cell-suppressed tables that are safe to
commit and to put in slides.

Five outputs, in increasing order of how much we actually care:

  1. lr_rnaseq_cohort_demographics.csv  -- who these people are (age, sex, ancestry)
  2. lr_rnaseq_burden_summary.csv       -- how sick they are, normalized by EHR-years
  3. lr_rnaseq_top_conditions.csv       -- THE RANKING: most common diagnoses
  4. lr_rnaseq_icd_chapters.csv         -- burden by broad disease category
  5. lr_rnaseq_immune_feasibility.csv   -- *** the one that decides the science ***
                                           case counts for immune-mediated diseases, i.e.
                                           "which phenotypes could a repertoire model ever
                                           be trained on with this cohort?"

SMALL-CELL SUPPRESSION: any count below --min-cell (default 20) is reported as "<20" rather
than an exact number. This is the standard AoU disclosure floor and it applies to anything
that leaves the Workbench -- including slides.

Usage:
  pixi run python3 analyze_disease_burden.py [--indir ~/pipeline_outputs/rnaseq/pheno]
      [--cohort ~/pipeline_outputs/rnaseq/lr_rnaseq_overlap_cohort.tsv]
      [--outdir ../results] [--min-cell 20]
"""
import argparse
import os
import sys
from datetime import date

import pandas as pd

# ---------------------------------------------------------------------------
# ICD-10 chapter rollup. First character of the ICD-10-CM code -> broad category.
# Deliberately dependency-light: no phecode mapping file to download, no SNOMED
# hierarchy walk. Coarse, but every code maps, and it is enough to answer "is this
# cohort's burden cardiovascular, metabolic, psychiatric, or immune?"
# ---------------------------------------------------------------------------
ICD10_CHAPTERS = {
    "A": "Infectious & parasitic", "B": "Infectious & parasitic",
    "C": "Neoplasms",             "D": "Neoplasms / blood & immune",
    "E": "Endocrine, nutritional & metabolic",
    "F": "Mental & behavioural",
    "G": "Nervous system",
    "H": "Eye & ear",
    "I": "Circulatory",
    "J": "Respiratory",
    "K": "Digestive",
    "L": "Skin & subcutaneous",
    "M": "Musculoskeletal & connective tissue",
    "N": "Genitourinary",
    "O": "Pregnancy & childbirth",
    "P": "Perinatal",
    "Q": "Congenital malformations",
    "R": "Symptoms & abnormal findings",
    "S": "Injury & poisoning", "T": "Injury & poisoning",
    "V": "External causes", "W": "External causes",
    "X": "External causes", "Y": "External causes",
    "Z": "Health status & contact with services",
}

# ---------------------------------------------------------------------------
# Immune-mediated / immune-relevant disease families.
#
# WHY THIS LIST AND NOT ANOTHER: these are the phenotypes where a T/B-cell receptor
# repertoire could plausibly carry signal that a genotype cannot -- because the disease
# involves the adaptive immune system actually having responded to something. Everything
# here is either autoimmune (immune attacks self), an infection (immune responded to a
# pathogen), an allergy/hypersensitivity (immune responded to a harmless antigen), or a
# lymphoid malignancy (the immune cells themselves are the disease).
#
# Matched case-insensitively against the SNOMED concept_name. Substring matching is crude
# and WILL over-match (e.g. "arthritis" catches osteoarthritis, which is not immune-mediated)
# -- so the output is a SCREENING table for feasibility, not a validated phenotype
# definition. Anything that looks promising here gets a proper concept-set definition later.
# ---------------------------------------------------------------------------
IMMUNE_TERMS = {
    "Autoimmune -- rheumatologic": [
        "rheumatoid arthritis", "systemic lupus", "sjogren", "sjögren",
        "systemic sclerosis", "scleroderma", "ankylosing spondylitis",
        "psoriatic arthritis", "vasculitis", "polymyalgia", "antiphospholipid",
    ],
    "Autoimmune -- endocrine": [
        "type 1 diabetes", "hashimoto", "graves disease", "autoimmune thyroiditis",
        "addison",
    ],
    "Autoimmune -- gastrointestinal": [
        "crohn", "ulcerative colitis", "inflammatory bowel", "celiac", "coeliac",
        "autoimmune hepatitis", "primary biliary",
    ],
    "Autoimmune -- neurologic": [
        "multiple sclerosis", "myasthenia gravis", "guillain", "neuromyelitis",
    ],
    "Autoimmune -- dermatologic": [
        "psoriasis", "vitiligo", "pemphigus", "alopecia areata", "lichen planus",
    ],
    "Allergy & hypersensitivity": [
        "asthma", "allergic rhinitis", "atopic dermatitis", "eczema", "anaphyla",
        "food allergy", "drug allergy", "urticaria", "hypersensitivity",
    ],
    "Chronic infection": [
        "hiv", "hepatitis b", "hepatitis c", "tuberculosis", "cytomegalovirus",
        "epstein-barr", "herpes zoster", "chronic viral",
    ],
    "Acute / recurrent infection": [
        "pneumonia", "sepsis", "influenza", "covid", "sars-cov-2", "urinary tract infection",
    ],
    "Lymphoid & haematologic malignancy": [
        "lymphoma", "leukemia", "leukaemia", "myeloma", "myelodysplastic",
    ],
    "Immunodeficiency & immunosuppression": [
        "immunodeficiency", "agammaglobulin", "hypogammaglobulin", "transplant",
        "graft versus host", "immunosuppress",
    ],
}


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def suppress(n, min_cell):
    """AoU disclosure floor: never emit an exact count below min_cell."""
    return f"<{min_cell}" if 0 < n < min_cell else str(int(n))


def blank_below(df, count_col, cols, min_cell):
    """Blank derived stats for rows whose underlying count is suppressed.

    Suppressing the count alone is not enough: a percentage plus a known denominator
    reconstructs the exact count, and a median age over a 5-person group is close to
    individual-level. Anything derived from a below-floor cell has to go too.
    """
    mask = (df[count_col] > 0) & (df[count_col] < min_cell)
    for c in cols:
        if c in df.columns:
            df.loc[mask, c] = pd.NA
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default=os.path.expanduser("~/pipeline_outputs/rnaseq/pheno"))
    ap.add_argument("--cohort", default=os.path.expanduser(
                        "~/pipeline_outputs/rnaseq/lr_rnaseq_overlap_cohort.tsv"))
    ap.add_argument("--outdir", default=os.path.join(
                        os.path.dirname(os.path.abspath(__file__)), "..", "results"))
    ap.add_argument("--min-cell", type=int, default=20,
                    help="Disclosure floor. Counts below this are reported as '<N'.")
    ap.add_argument("--top-n", type=int, default=50, help="How many conditions to rank")
    args = ap.parse_args()

    indir = os.path.expanduser(args.indir)
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    os.makedirs(outdir, exist_ok=True)

    paths = {k: os.path.join(indir, f"{k}.tsv")
             for k in ("person", "observation_period", "conditions")}
    for k, p in paths.items():
        if not os.path.exists(p):
            die(f"Missing {p}. Run query_overlap_phenotypes.py first (in a notebook).")

    person = pd.read_csv(paths["person"], sep="\t", dtype={"research_id": str})
    obs = pd.read_csv(paths["observation_period"], sep="\t", dtype={"research_id": str})
    cond = pd.read_csv(paths["conditions"], sep="\t", dtype={"research_id": str})
    cohort = pd.read_csv(os.path.expanduser(args.cohort), sep="\t", dtype=str)

    base = cohort[["research_id"] + [c for c in ("ancestry",) if c in cohort.columns]]
    base = base.merge(person, on="research_id", how="left") \
               .merge(obs[["research_id", "ehr_years"]], on="research_id", how="left")
    if "ancestry" not in base.columns:
        base["ancestry"] = "unlabeled"
    base["age"] = date.today().year - pd.to_numeric(base["year_of_birth"], errors="coerce")

    n_cohort = len(base)
    n_with_ehr = int(base["ehr_years"].notna().sum())
    print("=" * 68)
    print(f"LR x RNA-seq overlap cohort: {n_cohort:,} people")
    print(f"  with any EHR record: {n_with_ehr:,} ({n_with_ehr / max(n_cohort,1):.1%})")
    print("=" * 68)
    if n_with_ehr < n_cohort:
        print(f"  NOTE: {n_cohort - n_with_ehr:,} people have NO EHR data. They are in the")
        print("  genomic cohort but contribute nothing to a disease study. All prevalence")
        print("  denominators below use the EHR-having subset only.")

    # ---------------- 1. demographics ----------------
    demo = base.groupby("ancestry").agg(
        n_people=("research_id", "count"),
        n_with_ehr=("ehr_years", lambda s: int(s.notna().sum())),
        median_age=("age", "median"),
        median_ehr_years=("ehr_years", "median"),
    ).reset_index()
    demo = blank_below(demo, "n_people",
                       ["n_with_ehr", "median_age", "median_ehr_years"], args.min_cell)
    demo["n_people"] = demo["n_people"].apply(lambda n: suppress(n, args.min_cell))
    demo["n_with_ehr"] = demo["n_with_ehr"].apply(
        lambda n: "" if pd.isna(n) else suppress(n, args.min_cell))
    demo.to_csv(os.path.join(outdir, "lr_rnaseq_cohort_demographics.csv"), index=False)

    # ---------------- 2. per-person burden ----------------
    ehr = base[base["ehr_years"].notna() & (base["ehr_years"] > 0)].copy()
    per_person = cond.groupby("research_id").agg(
        n_distinct_conditions=("condition_concept_id", "nunique"),
        n_total_occurrences=("n_occurrences", "sum"),
    ).reset_index()
    ehr = ehr.merge(per_person, on="research_id", how="left")
    ehr[["n_distinct_conditions", "n_total_occurrences"]] = \
        ehr[["n_distinct_conditions", "n_total_occurrences"]].fillna(0)
    # The normalized burden metric. Raw counts are confounded by chart length; this is not.
    ehr["conditions_per_ehr_year"] = ehr["n_distinct_conditions"] / ehr["ehr_years"]

    burden = ehr.groupby("ancestry").agg(
        n_people=("research_id", "count"),
        median_distinct_conditions=("n_distinct_conditions", "median"),
        median_ehr_years=("ehr_years", "median"),
        median_conditions_per_ehr_year=("conditions_per_ehr_year", "median"),
    ).reset_index()
    burden = blank_below(burden, "n_people",
                         ["median_distinct_conditions", "median_ehr_years",
                          "median_conditions_per_ehr_year"], args.min_cell)
    burden["n_people"] = burden["n_people"].apply(lambda n: suppress(n, args.min_cell))
    burden = burden.round(2)
    burden.to_csv(os.path.join(outdir, "lr_rnaseq_burden_summary.csv"), index=False)
    print("\nPer-person burden (median distinct conditions per EHR-year), by ancestry:")
    print(burden.to_string(index=False))

    # ---------------- 3. THE RANKING ----------------
    denom = len(ehr)
    cond_ehr = cond[cond["research_id"].isin(set(ehr["research_id"]))]
    ranking = cond_ehr.groupby(["condition_concept_id", "condition_name"]).agg(
        n_people=("research_id", "nunique"),
        total_occurrences=("n_occurrences", "sum"),
    ).reset_index().sort_values("n_people", ascending=False)
    ranking["prevalence_pct"] = (ranking["n_people"] / max(denom, 1) * 100).round(1)
    # Chronic-vs-incidental hint: a condition coded many times per affected person is more
    # likely a real, followed diagnosis than a single-visit code.
    ranking["mean_codes_per_affected"] = \
        (ranking["total_occurrences"] / ranking["n_people"]).round(1)
    n_all_concepts = len(ranking)
    ranking = ranking[ranking["n_people"] >= args.min_cell]
    cols = ["condition_name", "n_people", "prevalence_pct", "mean_codes_per_affected"]

    # FULL table -- every condition that clears the disclosure floor, not just the top N.
    # "Quantify the diseases fully" means this file; the top-N file below is only a
    # presentation convenience carved out of it.
    ranking[cols].to_csv(os.path.join(outdir, "lr_rnaseq_all_conditions.csv"), index=False)
    ranking[cols].head(args.top_n).to_csv(
        os.path.join(outdir, "lr_rnaseq_top_conditions.csv"), index=False)

    print(f"\nCONDITION COVERAGE (denominator = {denom:,} people with EHR):")
    print(f"  distinct conditions recorded in this cohort : {n_all_concepts:,}")
    print(f"  reportable at n >= {args.min_cell:<3}                        : {len(ranking):,}")
    print(f"  suppressed as below the disclosure floor     : "
          f"{n_all_concepts - len(ranking):,}")
    print("  -> full ranked list: lr_rnaseq_all_conditions.csv")
    print(f"\nTop 20 of {len(ranking):,}:")
    print(ranking[["condition_name", "n_people", "prevalence_pct"]].head(20).to_string(index=False))

    # ---------------- 4. ICD-10 chapters ----------------
    icd = cond_ehr[cond_ehr["source_vocabulary"].astype(str).str.startswith("ICD10")].copy()
    if len(icd):
        icd["chapter"] = icd["source_code"].astype(str).str[0].map(ICD10_CHAPTERS) \
                                                            .fillna("Unmapped")
        chap = icd.groupby("chapter").agg(
            n_people=("research_id", "nunique"),
            n_distinct_codes=("source_code", "nunique"),
        ).reset_index().sort_values("n_people", ascending=False)
        chap["pct_of_cohort"] = (chap["n_people"] / max(denom, 1) * 100).round(1)
        chap = blank_below(chap, "n_people", ["pct_of_cohort", "n_distinct_codes"],
                           args.min_cell)
        chap["n_people"] = chap["n_people"].apply(lambda n: suppress(n, args.min_cell))
        chap.to_csv(os.path.join(outdir, "lr_rnaseq_icd_chapters.csv"), index=False)
        print("\nBurden by ICD-10 chapter:")
        print(chap.head(12).to_string(index=False))

        # Middle granularity: the 3-character ICD-10 category (E11 = type 2 diabetes,
        # J45 = asthma, M06 = rheumatoid arthritis). Chapters are too coarse to be a
        # disease list and individual SNOMED concepts are too fine (dozens of near-identical
        # variants per disease). This is the level most people mean by "a disease".
        icd["icd3"] = icd["source_code"].astype(str).str[:3]
        icd3 = icd.groupby("icd3").agg(
            n_people=("research_id", "nunique"),
            n_distinct_concepts=("condition_concept_id", "nunique"),
            commonest_name=("condition_name",
                            lambda s: s.value_counts().idxmax() if len(s) else ""),
        ).reset_index().sort_values("n_people", ascending=False)
        icd3["prevalence_pct"] = (icd3["n_people"] / max(denom, 1) * 100).round(1)
        n_icd3_all = len(icd3)
        icd3 = icd3[icd3["n_people"] >= args.min_cell]
        icd3.to_csv(os.path.join(outdir, "lr_rnaseq_icd3_categories.csv"), index=False)
        print(f"\nICD-10 3-character categories: {n_icd3_all:,} present, "
              f"{len(icd3):,} reportable at n >= {args.min_cell} "
              f"-> lr_rnaseq_icd3_categories.csv")
        print(icd3[["icd3", "commonest_name", "n_people", "prevalence_pct"]]
              .head(15).to_string(index=False))
    else:
        print("\nNo ICD-10 source codes found -- chapter rollup skipped. Check whether this "
              "CDR populates condition_source_concept_id, or fall back to a SNOMED "
              "hierarchy rollup via concept_ancestor.")

    # ---------------- 5. THE FEASIBILITY TABLE ----------------
    names = cond_ehr["condition_name"].astype(str).str.lower()
    rows = []
    for family, terms in IMMUNE_TERMS.items():
        mask = pd.Series(False, index=cond_ehr.index)
        for t in terms:
            mask |= names.str.contains(t, regex=False, na=False)
        sub = cond_ehr[mask]
        n_cases = sub["research_id"].nunique()
        rows.append({
            "immune_family": family,
            "n_cases": n_cases,
            "pct_of_cohort": round(n_cases / max(denom, 1) * 100, 1),
            "n_distinct_concepts": sub["condition_concept_id"].nunique(),
            # Rough, honest feasibility read. These thresholds are conventions for
            # case/control work, not guarantees.
            "verdict": ("enough for a real case/control model" if n_cases >= 200 else
                        "descriptive only -- underpowered" if n_cases >= 50 else
                        "too few -- not viable in this cohort"),
        })
    feas = pd.DataFrame(rows).sort_values("n_cases", ascending=False)
    feas = blank_below(feas, "n_cases", ["pct_of_cohort", "n_distinct_concepts"],
                       args.min_cell)
    feas["n_cases"] = feas["n_cases"].apply(lambda n: suppress(n, args.min_cell))
    feas.to_csv(os.path.join(outdir, "lr_rnaseq_immune_feasibility.csv"), index=False)
    print("\n" + "=" * 68)
    print("IMMUNE-MEDIATED DISEASE FEASIBILITY -- the table that decides the science")
    print("=" * 68)
    print(feas.to_string(index=False))
    print("\nSubstring matching over-matches (e.g. 'arthritis' catches osteoarthritis).")
    print("Treat this as a SCREEN: anything that clears the bar gets a proper concept-set")
    print("definition before it is used as a real phenotype.")

    print(f"\nFive de-identified CSVs written to {outdir} -- safe to commit.")


if __name__ == "__main__":
    main()
