#!/usr/bin/env python3
"""Pull the EHR disease record for the LR x RNA-seq overlap cohort (Task 2, step 2).

WHERE THIS RUNS: **a Jupyter notebook on the Workbench, not the plain shell.** The CDR
(Curated Data Repository -- AoU's versioned clinical database) lives in a separate GCP
project and is only reachable via BigQuery with a fully-qualified table path;
ENVIRONMENT.md quirk #5 confirms it is NOT exposed via shell env vars. In a notebook the
variable `WORKSPACE_CDR` is normally set for you.

  %run query_overlap_phenotypes.py --cohort ~/pipeline_outputs/rnaseq/lr_rnaseq_overlap_cohort.tsv

WHAT IT PULLS, per person in the cohort:
  1. demographics          -- year of birth, sex at birth  (person)
  2. EHR observation window -- first and last date any record exists (observation_period)
  3. every condition       -- SNOMED standard concept + the raw ICD-10/9 source code
                              (condition_occurrence joined to concept)

WHY (2) MATTERS AND IS NOT OPTIONAL: raw disease counts are worthless without it. Someone
with 12 years of records will show more diagnoses than someone with 1 year of records,
regardless of how healthy they are. Every prevalence number must be interpreted against
how long we were actually watching. This is the single most common way EHR disease-burden
analyses go wrong.

PRIVACY: output contains real research_ids and individual-level clinical data -- the most
sensitive thing this workstream has ever touched. It stays VM-local under
~/pipeline_outputs, is never committed, and only the aggregated output of
analyze_disease_burden.py (group-level counts, small cells suppressed) goes anywhere near
git. Do not print per-person rows into notebook output that gets saved.
"""
import argparse
import os
import sys

import pandas as pd

MIN_CELL = 20  # cells smaller than this are suppressed downstream; see analyze_disease_burden.py


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def get_cdr(explicit):
    """Resolve the fully-qualified CDR dataset, e.g. wb-silky-artichoke-2408.C2025Q4R6."""
    if explicit:
        return explicit
    for var in ("WORKSPACE_CDR", "CDR_STORAGE_PATH"):
        v = os.environ.get(var)
        if v:
            return v
    die("Could not resolve the CDR dataset. Pass --cdr explicitly, e.g.\n"
        "  --cdr wb-silky-artichoke-2408.C2025Q4R6\n"
        "(ENVIRONMENT.md quirk #5: the CDR is not exposed via shell env vars on this VM; "
        "in a notebook, WORKSPACE_CDR is usually set. Confirm the current value from a "
        "Dataset Builder-generated snippet rather than trusting this default -- the "
        "dataset id changes with every CDR release.)")


def run_query(sql, project):
    try:
        import pandas_gbq  # noqa: F401
        return pd.read_gbq(sql, dialect="standard", progress_bar_type=None)
    except ImportError:
        pass
    try:
        from google.cloud import bigquery
    except ImportError:
        die("Neither pandas_gbq nor google-cloud-bigquery is importable. Both normally "
            "ship in the Workbench notebook image -- if this fires, you are probably in a "
            "pixi shell rather than the notebook kernel. Run this from a notebook.")
    client = bigquery.Client(project=project) if project else bigquery.Client()
    return client.query(sql).to_dataframe()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", default=os.path.expanduser(
                        "~/pipeline_outputs/rnaseq/lr_rnaseq_overlap_cohort.tsv"),
                    help="Output of check_lr_rnaseq_overlap.py (research_id, ancestry, ...)")
    ap.add_argument("--cdr", default=None,
                    help="Fully-qualified CDR dataset, e.g. wb-silky-artichoke-2408.C2025Q4R6")
    ap.add_argument("--project", default=os.environ.get("GOOGLE_PROJECT"),
                    help="Billing project for the BigQuery job")
    ap.add_argument("--outdir", default=os.path.expanduser("~/pipeline_outputs/rnaseq/pheno"),
                    help="VM-local output directory (NEVER commit its contents)")
    args = ap.parse_args()

    cdr = get_cdr(args.cdr)
    cohort_path = os.path.expanduser(args.cohort)
    if not os.path.exists(cohort_path):
        die(f"Cohort file not found: {cohort_path}. Run check_lr_rnaseq_overlap.py first.")

    cohort = pd.read_csv(cohort_path, sep="\t", dtype=str)
    if "research_id" not in cohort.columns:
        die(f"Cohort file has no research_id column. Columns: {list(cohort.columns)}")
    ids = sorted(set(cohort["research_id"].dropna()))
    if not ids:
        die("Cohort file has zero research_ids.")
    print(f"Cohort: {len(ids):,} people. CDR: {cdr}")

    # Inline the id list. Fine for a few hundred / few thousand people; if the overlap ever
    # runs to tens of thousands, materialize a temp table instead.
    id_list = ",".join(ids)

    os.makedirs(args.outdir, exist_ok=True)

    # ---------------- 1. demographics ----------------
    sql_person = f"""
    SELECT
      p.person_id            AS research_id,
      p.year_of_birth,
      c_sex.concept_name     AS sex_at_birth
    FROM `{cdr}.person` p
    LEFT JOIN `{cdr}.concept` c_sex
           ON c_sex.concept_id = p.sex_at_birth_concept_id
    WHERE p.person_id IN ({id_list})
    """

    # ---------------- 2. EHR observation window ----------------
    # The denominator for every prevalence number. Without this you cannot tell "healthy"
    # from "we only had 6 months of their chart".
    sql_obs = f"""
    SELECT
      person_id                                              AS research_id,
      MIN(observation_period_start_date)                     AS ehr_start,
      MAX(observation_period_end_date)                       AS ehr_end,
      DATE_DIFF(MAX(observation_period_end_date),
                MIN(observation_period_start_date), DAY)/365.25 AS ehr_years
    FROM `{cdr}.observation_period`
    WHERE person_id IN ({id_list})
    GROUP BY person_id
    """

    # ---------------- 3. conditions ----------------
    # Two views of the same diagnosis, deliberately:
    #   condition_concept_id  -> SNOMED standard concept (harmonized across all sites)
    #   condition_source_concept_id -> the raw ICD-10-CM/ICD-9-CM code the site actually
    #                                  recorded. Kept because the ICD code's first letter
    #                                  gives a free, dependency-light chapter rollup
    #                                  (A/B = infectious, C/D = neoplasm, M = musculoskeletal...)
    #                                  with no external phecode mapping file needed.
    # Aggregated to one row per person per condition, with occurrence count and first/last
    # date -- so a single stray code can be told apart from a chronic, repeatedly-coded one.
    sql_cond = f"""
    SELECT
      co.person_id                       AS research_id,
      co.condition_concept_id,
      c_std.concept_name                 AS condition_name,
      c_src.vocabulary_id                AS source_vocabulary,
      c_src.concept_code                 AS source_code,
      COUNT(*)                           AS n_occurrences,
      MIN(co.condition_start_date)       AS first_date,
      MAX(co.condition_start_date)       AS last_date
    FROM `{cdr}.condition_occurrence` co
    LEFT JOIN `{cdr}.concept` c_std ON c_std.concept_id = co.condition_concept_id
    LEFT JOIN `{cdr}.concept` c_src ON c_src.concept_id = co.condition_source_concept_id
    WHERE co.person_id IN ({id_list})
      AND co.condition_concept_id != 0
    GROUP BY 1,2,3,4,5
    """

    for name, sql in (("person", sql_person),
                      ("observation_period", sql_obs),
                      ("conditions", sql_cond)):
        print(f"  querying {name} ...", flush=True)
        df = run_query(sql, args.project)
        out = os.path.join(args.outdir, f"{name}.tsv")
        df.to_csv(out, sep="\t", index=False)
        print(f"    {len(df):,} rows -> {out}")

    print("\nDone. All three files are VM-LOCAL and contain individual-level clinical data.")
    print("DO NOT commit them, do not paste per-person rows into chat.")
    print("Next: pixi run python3 analyze_disease_burden.py")


if __name__ == "__main__":
    main()
