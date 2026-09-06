#!/usr/bin/env python3
"""Pull REAL clinical diagnosis labels (not allele-proxy labels) for the hla_popgen embedding
cohort, so `08_embedding_compare.py --real-disease-labels <path>` can color the locked UMAP by
actual EHR-confirmed disease rather than by carrying a disease-associated allele.

## Why this exists (2026-09 chat discussion)

The first disease-coloring pass in 08_embedding_compare.py colored by ALLELE carrier status
(e.g. "carries B*27"), which is close to tautological on an allele-identity embedding -- carriers
of one allele are, by construction, near each other in that space. This script instead pulls
whether each person has an ACTUAL diagnosis on record in the CDR, using the exact HLA-linked
disease definitions (ICD-10 3-char code + SNOMED name-substring fallback) already built and
reviewed on `origin/aleix/hla-resolve-phase1`'s `deep_immune_breakdown.py` -- reused verbatim here
(see HLA_LINKED below) rather than re-derived, since that work was already vetted.

## Where this runs

A plain shell works on this VM (confirmed 2026-09: `WORKSPACE_CDR` is set even outside a notebook
kernel here, unlike the general case Aleix's own scripts warn about -- verify this is still true
before assuming it on a different VM). Needs `google-cloud-bigquery` or `pandas-gbq` (both ship in
the Workbench image).

## Privacy (same posture as aleix/RNA-seq/scripts/query_overlap_phenotypes.py)

Output contains real per-person diagnosis labels keyed by `person_id` -- stays VM-local under
`~/pipeline_outputs/rnaseq/pheno/`, is NEVER committed, and is not printed to stdout. Only a
small-cell-suppressed (`--min-cell`, default 20) aggregate count table is safe to view/commit, and
even that is written to a separate file this script does not auto-commit.

Usage:
  python3 scripts/hla_popgen/12_disease_phenotypes.py \\
      --person-ids-from-cache ~/repos/pilot-validation/reports/hla_popgen/08_embedding_compare/lr/embedding_cache.pkl
"""
import argparse
import os
import pickle
import sys

import pandas as pd

MIN_CELL = 20

# ---------------------------------------------------------------------------
# Verbatim from origin/aleix/hla-resolve-phase1's deep_immune_breakdown.py HLA_LINKED list
# (2026-09, commit dc9563b) -- each entry: (label, [icd3_codes], [name_substrings]). A person
# matches if EITHER an ICD-10 3-char code matches OR a SNOMED concept-name substring matches.
# Reused rather than re-derived: this list was already reviewed (ICD3-primary with a documented
# reason for every substring-only fallback, e.g. "K90 dropped: includes other malabsorption
# syndromes" for celiac). Four of our original 12 allele-based groups (B*57:01/B*58:01/B*15:02
# drug-hypersensitivity reactions, A*29:02 birdshot chorioretinopathy) have no entry here --
# they're rare/differently-coded phenotypes Aleix's screen didn't cover, not omitted by us; left
# out rather than guessed at with an unreviewed definition.
# ---------------------------------------------------------------------------
HLA_LINKED = [
    ("Ankylosing spondylitis (HLA-B27)", ["M45"], ["ankylosing spondylitis"]),
    ("Celiac disease (HLA-DQ2/DQ8)", [], ["celiac", "coeliac"]),
    ("Type 1 diabetes (HLA-DR3/DR4)", ["E10"], ["type 1 diabetes"]),
    ("Psoriasis (HLA-Cw6)", ["L40"], ["psoriasis"]),
    ("Psoriatic arthritis (HLA-Cw6/B27)", [], ["psoriatic arthritis"]),
    ("Multiple sclerosis (HLA-DRB1*15:01)", ["G35"], ["multiple sclerosis"]),
    ("Rheumatoid arthritis (HLA shared epitope)", ["M05", "M06"], ["rheumatoid arthritis"]),
    ("Systemic lupus erythematosus (HLA-DR2/DR3)", ["M32"], ["systemic lupus", "lupus erythematosus"]),
    ("Graves disease (HLA-DR3)", [], ["graves"]),
    ("Narcolepsy (HLA-DQB1*06:02)", [], ["narcolepsy"]),
    ("Behcet disease (HLA-B51)", [], ["behcet"]),
]


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def get_cdr(explicit):
    if explicit:
        return explicit
    for var in ("WORKSPACE_CDR", "CDR_STORAGE_PATH"):
        v = os.environ.get(var)
        if v:
            return v
    die("Could not resolve the CDR dataset. Pass --cdr explicitly.")


def run_query(sql, project):
    try:
        import pandas_gbq
        return pandas_gbq.read_gbq(sql, dialect="standard", progress_bar_type=None,
                                   project_id=project)
    except ImportError:
        pass
    from google.cloud import bigquery
    client = bigquery.Client(project=project) if project else bigquery.Client()
    return client.query(sql).to_dataframe()


def load_person_ids(args):
    if args.person_ids_from_cache:
        with open(os.path.expanduser(args.person_ids_from_cache), "rb") as f:
            cache = pickle.load(f)
        ids = sorted(set(cache["mat_full_index"]))
        print(f"Loaded {len(ids)} person_ids from embedding cache "
              f"{args.person_ids_from_cache!r}.", file=sys.stderr)
        return ids
    if args.cohort_membership:
        df = pd.read_csv(os.path.expanduser(args.cohort_membership), sep="\t", dtype=str)
        col = "person_id" if "person_id" in df.columns else "research_id"
        sub = df[df["in_lr"].astype(str).str.lower() == "true"] if "in_lr" in df.columns else df
        ids = sorted(set(sub[col].dropna()))
        print(f"Loaded {len(ids)} person_ids from {args.cohort_membership!r} (in_lr).",
              file=sys.stderr)
        return ids
    die("Pass either --person-ids-from-cache or --cohort-membership.")


def assign_labels(cond_df):
    """cond_df: one row per (person_id, condition_source_value [ICD], condition_name [SNOMED]).
    Returns {person_id: [matched_label, ...]} -- ALL matches kept (not just first), so the caller
    can report both the "first match" simplification used for a single-color plot AND the true
    multi-match rate honestly."""
    icd3 = cond_df["condition_source_value"].astype(str).str.upper().str[:3]
    name = cond_df["condition_name"].astype(str).str.lower()
    matches = {}
    for label, icd_codes, substrings in HLA_LINKED:
        mask = icd3.isin(icd_codes) if icd_codes else pd.Series(False, index=cond_df.index)
        for s in substrings:
            mask = mask | name.str.contains(s, na=False, regex=False)
        for pid in cond_df.loc[mask, "person_id"].unique():
            matches.setdefault(pid, []).append(label)
    return matches


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--person-ids-from-cache", default=None,
                    help="Path to 08_embedding_compare.py's embedding_cache.pkl -- pulls exactly "
                         "the people being colored, guaranteeing a match.")
    ap.add_argument("--cohort-membership", default=None,
                    help="Alternative: path to cohort_membership.tsv, filtered to in_lr.")
    ap.add_argument("--cdr", default=None)
    ap.add_argument("--project", default=os.environ.get("GOOGLE_PROJECT"))
    ap.add_argument("--outdir", default=os.path.expanduser("~/pipeline_outputs/rnaseq/pheno"),
                    help="VM-local output directory (NEVER commit its contents).")
    ap.add_argument("--min-cell", type=int, default=MIN_CELL)
    args = ap.parse_args()

    cdr = get_cdr(args.cdr)
    ids = load_person_ids(args)
    if not ids:
        die("Zero person_ids resolved.")
    os.makedirs(args.outdir, exist_ok=True)

    id_list = ",".join(ids)
    print(f"Cohort: {len(ids):,} people. CDR: {cdr}", file=sys.stderr)
    print("Querying condition_occurrence + concept ...", file=sys.stderr)
    sql_cond = f"""
    SELECT DISTINCT
      co.person_id                       AS person_id,
      c_std.concept_name                 AS condition_name,
      c_src.concept_code                 AS condition_source_value
    FROM `{cdr}.condition_occurrence` co
    LEFT JOIN `{cdr}.concept` c_std
           ON c_std.concept_id = co.condition_concept_id
    LEFT JOIN `{cdr}.concept` c_src
           ON c_src.concept_id = co.condition_source_concept_id
    WHERE co.person_id IN ({id_list})
      AND co.condition_concept_id != 0
    """
    cond = run_query(sql_cond, args.project)
    print(f"  {len(cond):,} condition rows for {cond['person_id'].nunique():,} people with "
          f">=1 condition on record.", file=sys.stderr)

    matches = assign_labels(cond)
    n_any = len(matches)
    n_multi = sum(1 for v in matches.values() if len(v) > 1)
    print(f"{n_any} of {len(ids)} people matched >=1 HLA-linked disease definition "
          f"({n_multi} matched more than one).", file=sys.stderr)

    # Per-person file: FIRST match only, for a single-color plot -- VM-local, never committed.
    rows = [{"person_id": pid, "disease_label": labels[0], "n_matches": len(labels)}
            for pid, labels in matches.items()]
    for pid in ids:
        if pid not in matches:
            rows.append({"person_id": pid, "disease_label": None, "n_matches": 0})
    out_df = pd.DataFrame(rows)
    person_path = os.path.join(args.outdir, "person_disease_labels.tsv")
    out_df.to_csv(person_path, sep="\t", index=False)
    print(f"Wrote {person_path!r} ({len(out_df)} rows, person-level, VM-local, DO NOT COMMIT).",
          file=sys.stderr)

    # Wide multi-hot file: one boolean column per disease, ALL matches kept independently -- a
    # person diagnosed with both T1D and celiac counts toward BOTH tests. Needed for
    # 13_disease_allele_association.py's per-disease carrier-rate comparison; the single-label
    # file above (first match only) would silently under-count co-morbid people for every disease
    # except their first-listed one. Also VM-local, never committed.
    labels_all = [lbl for lbl, *_j in HLA_LINKED]
    wide_rows = []
    for pid in ids:
        hit = set(matches.get(pid, []))
        wide_rows.append({"person_id": pid, **{lbl: (lbl in hit) for lbl in labels_all}})
    wide_df = pd.DataFrame(wide_rows)
    wide_path = os.path.join(args.outdir, "person_disease_labels_wide.tsv")
    wide_df.to_csv(wide_path, sep="\t", index=False)
    print(f"Wrote {wide_path!r} ({len(wide_df)} rows x {len(labels_all)} disease columns, "
          f"person-level, VM-local, DO NOT COMMIT).", file=sys.stderr)

    # Aggregate, small-cell-suppressed -- safe to look at / eventually commit if wanted.
    counts = out_df["disease_label"].value_counts()
    agg_rows = []
    for label, _icd, _sub in HLA_LINKED:
        n = int(counts.get(label, 0))
        agg_rows.append({"disease": label,
                         "n_people": n if n >= args.min_cell else f"<{args.min_cell}"})
    agg_df = pd.DataFrame(agg_rows)
    agg_path = os.path.join(args.outdir, "disease_label_counts_aggregate.tsv")
    agg_df.to_csv(agg_path, sep="\t", index=False)
    print(f"Wrote {agg_path!r} (aggregate, small-cell-suppressed, safe to view):", file=sys.stderr)
    print(agg_df.to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
