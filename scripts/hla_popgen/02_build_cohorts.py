#!/usr/bin/env python3
"""Build Table 4 (`cohort_membership.tsv`) per scripts/hla_popgen/SCHEMA.md -- one row per person,
joining the long-read Immuannot calls (Table 1, from 01_extract_rich.py) with cohort/platform
metadata, AoU's continuous ancestry admixture, and the AoU-native short-read HLA calls.

## The three-cohort definition lives here and nowhere else (SCHEMA.md Table 4)

- **Cohort "in_lr"**: has Immuannot long-read calls (Table 1 non-empty for this person). ~12K.
- **Cohort "in_sr"**: has AoU-native short-read HLA calls (`hla_genotypes.tsv`). ~500K, a strict
  superset in practice (short-read WGS is run on everyone; long-read on a subset).
- **Cohort 3 -- `td_stratum`**: NOT a single hardcoded threshold. It's a sweep-friendly set of
  flags (`exact` / `near` / `loose` / `distant`) computed from `max_template_distance` over a
  person's classical genes, so every downstream frequency/structure figure can be reproduced at
  each stratum and the *movement* of a statistic as the filter tightens is itself the result
  (SCHEMA.md: "Cohort 3 is a sweep, not a fixed threshold").

## Correctness details this script exists to get right

- **`ancestry_pred` is lowercase in the real file** (ENVIRONMENT.md quirk #30) -- normalized with
  `normalize_ancestry()` (uppercase) immediately on load from EVERY source that carries it
  (`immuannot_cohort_full.tsv` AND `ancestry_preds.tsv` independently -- quirk #30's own incident
  was caused by checking only one of two ancestry-bearing files).
- **`probabilities` is a stringified array** in AFR/AMR/EAS/EUR/MID/SAS order (ENVIRONMENT.md
  "Confirmed data locations" table) -- parsed into `p_afr`..`p_sas` columns, the continuous
  admixture signal this project's non-linear PRS work actually wants (a hard ancestry label alone
  throws that away).
- **Hard-fails on a missing gcsfuse mount** (ENVIRONMENT.md quirk #26) as its very first action,
  with the exact remount command in the error -- this script is the one script in this sub-project
  that touches AoU-bucket-only files (`ancestry_preds.tsv`, `hla_genotypes.tsv`), unlike
  00_recon_vm.py/01_extract_rich.py which only ever touch the local `~/pipeline_outputs/` disk.

Usage:
    # Against fixtures (make_fixtures.py already lays these files out at --outroot):
    python3 scripts/hla_popgen/tests/make_fixtures.py --outroot /tmp/hla_fixtures -n 200
    python3 scripts/hla_popgen/01_extract_rich.py --outroot /tmp/hla_fixtures --sample
    python3 scripts/hla_popgen/02_build_cohorts.py --outroot /tmp/hla_fixtures \\
        --table1 /tmp/hla_fixtures/hla_calls_rich.sample.tsv \\
        --cohort-full /tmp/hla_fixtures/immuannot_cohort_full.tsv \\
        --ancestry-preds /tmp/hla_fixtures/ancestry_preds.tsv \\
        --hla-genotypes /tmp/hla_fixtures/hla_genotypes.tsv \\
        --skip-mount-check --sample

    # Real run on the VM (paths below are what ENVIRONMENT.md documents; override with flags if
    # they move -- never hand-build a `pooled/` path, always resolve via the v9 manifests):
    python3 scripts/hla_popgen/02_build_cohorts.py
"""
import argparse
import ast
import json
import os
import sys

import numpy as np
import pandas as pd

DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_MOUNT = os.path.expanduser("~/mnt/aou-controlled")
DEFAULT_BILLING_PROJECT = "wb-cordial-leechee-9743"
DEFAULT_ANCESTRY_PREDS = "v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
DEFAULT_HLA_GENOTYPES = "v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv"

ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]

# The 8 classical genes this project's tools call (ENVIRONMENT.md: "AoU-native HLA typing ...
# Restrict comparisons to the 8 classical genes our tools call").
CLASSICAL_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
                   "HLA-DRB1"]


def normalize_ancestry(val):
    """Uppercase and strip. Handles NaN/None/blank gracefully. ENVIRONMENT.md quirk #30: the real
    `ancestry_pred` column is lowercase; always normalize immediately on load, from every file that
    carries it independently -- never assume one normalized file means all of them are."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return None
    s = str(val).strip().upper()
    return s if s and s != "NAN" else None


def check_mount(mount_path, billing_project):
    """Hard-fail as the very first action if the gcsfuse mount isn't there (ENVIRONMENT.md quirk
    #26). This is the one script in this sub-project that needs it."""
    if os.path.isdir(mount_path) and os.listdir(mount_path):
        return
    sys.exit(
        f"FATAL: gcsfuse mount not found or empty at {mount_path!r}.\n"
        f"Remount with:\n"
        f"  mkdir -p {mount_path}\n"
        f"  gcsfuse --billing-project {billing_project} --implicit-dirs "
        f"vwb-aou-datasets-controlled {mount_path}\n"
        f"Then verify with `ls {mount_path}` BEFORE rerunning this script (ENVIRONMENT.md quirk "
        f"#14: gcsfuse needs a beat to become ready after mounting -- don't chain the mount and "
        f"the consumer in one paste)."
    )


def parse_probabilities(raw):
    """Parse the stringified probabilities array, ordered AFR/AMR/EAS/EUR/MID/SAS, into a dict of
    p_afr..p_sas. Returns a dict of NaN for every key if raw is missing/malformed (fail loud via
    the printed warning, not a silent wrong join)."""
    keys = [f"p_{a.lower()}" for a in ANCESTRY_ORDER]
    if raw is None or (isinstance(raw, float) and np.isnan(raw)):
        return {k: np.nan for k in keys}
    try:
        vals = ast.literal_eval(raw) if isinstance(raw, str) else list(raw)
        if len(vals) != 6:
            raise ValueError(f"expected 6 probabilities, got {len(vals)}")
        return dict(zip(keys, [float(v) for v in vals]))
    except (ValueError, SyntaxError) as e:
        print(f"  WARNING: could not parse probabilities {raw!r}: {e}", file=sys.stderr)
        return {k: np.nan for k in keys}


def compute_td_stratum(max_td):
    """Sweep-friendly stratum, NOT a single hardcoded threshold (SCHEMA.md: "Cohort 3 is a sweep").
    This function is the one place the thresholds live -- every downstream figure should call this
    rather than re-deriving its own cutoff."""
    if pd.isna(max_td):
        return None
    if max_td == 0:
        return "exact"
    if max_td <= 1:
        return "near"
    if max_td <= 5:
        return "loose"
    return "distant"


def load_table1(path):
    if not os.path.exists(path):
        sys.exit(f"FATAL: Table 1 not found at {path!r}. Run 01_extract_rich.py first.")
    df = pd.read_csv(path, sep="\t", dtype={"person_id": str})
    return df


def load_cohort_full(path):
    if not os.path.exists(path):
        sys.exit(f"FATAL: cohort metadata file not found at {path!r} (expected columns: "
                 f"person_id, platform, trim_tier, ancestry_pred).")
    df = pd.read_csv(path, sep="\t", dtype={"person_id": str})
    if "ancestry_pred" in df.columns:
        df["ancestry_pred"] = df["ancestry_pred"].apply(normalize_ancestry)
    return df


def load_ancestry_preds(path):
    if not os.path.exists(path):
        sys.exit(f"FATAL: ancestry_preds.tsv not found at {path!r}.")
    df = pd.read_csv(path, sep="\t", dtype={"research_id": str})
    df = df.rename(columns={"research_id": "person_id"})
    df["ancestry_pred_sr"] = df["ancestry_pred"].apply(normalize_ancestry)
    prob_cols = df["probabilities"].apply(parse_probabilities).apply(pd.Series)
    df = pd.concat([df[["person_id", "ancestry_pred_sr"]], prob_cols], axis=1)
    return df


def load_hla_genotypes(path):
    if not os.path.exists(path):
        sys.exit(f"FATAL: hla_genotypes.tsv (AoU-native short-read HLA calls) not found at "
                 f"{path!r}.")
    df = pd.read_csv(path, sep="\t", dtype={"research_id": str})
    df = df.rename(columns={"research_id": "person_id"})
    return df[["person_id"]].drop_duplicates()


def build_per_person_lr_stats(table1_df):
    """From Table 1's long grain, derive per-person: max/mean template_distance, n_genes_called,
    n_genes_exact -- restricted to the 8 classical genes (ENVIRONMENT.md: comparisons are
    restricted to the classical genes AoU also types)."""
    classical = table1_df[table1_df["gene"].isin(CLASSICAL_GENES)].copy()
    if classical.empty:
        return pd.DataFrame(columns=["person_id", "max_template_distance",
                                      "mean_template_distance", "n_genes_called",
                                      "n_genes_exact"])
    grouped = classical.groupby("person_id")
    out = grouped.agg(
        max_template_distance=("template_distance", "max"),
        mean_template_distance=("template_distance", "mean"),
    ).reset_index()
    n_called = grouped["gene"].nunique().rename("n_genes_called").reset_index()
    n_exact = classical[classical["template_distance"] == 0].groupby("person_id")["gene"] \
        .nunique().rename("n_genes_exact").reset_index()
    out = out.merge(n_called, on="person_id", how="left")
    out = out.merge(n_exact, on="person_id", how="left")
    out["n_genes_exact"] = out["n_genes_exact"].fillna(0).astype(int)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT,
                    help="Where per-person data lives / where cohort_membership.tsv is written "
                         "(SCHEMA.md hard rule #5: per-person rows stay in ~/pipeline_outputs/, "
                         "not reports/).")
    ap.add_argument("--table1", default=None,
                    help="Path to hla_calls_rich.tsv (Table 1). Default: <outroot>/"
                         "hla_calls_rich.tsv")
    ap.add_argument("--cohort-full", default=None,
                    help="Path to immuannot_cohort_full.tsv (platform/trim_tier/ancestry_pred). "
                         "Default: <outroot>/immuannot_cohort_full.tsv")
    ap.add_argument("--mount", default=DEFAULT_MOUNT, help="gcsfuse mount root.")
    ap.add_argument("--billing-project", default=DEFAULT_BILLING_PROJECT)
    ap.add_argument("--ancestry-preds", default=None,
                    help=f"Path to ancestry_preds.tsv. Default: <mount>/{DEFAULT_ANCESTRY_PREDS}")
    ap.add_argument("--hla-genotypes", default=None,
                    help=f"Path to hla_genotypes.tsv. Default: <mount>/{DEFAULT_HLA_GENOTYPES}")
    ap.add_argument("--skip-mount-check", action="store_true",
                    help="Skip the gcsfuse mount hard-check (for fixture/local testing only -- "
                         "NEVER pass this on the real VM).")
    ap.add_argument("--sample", action="store_true",
                    help="Write to a sample-suffixed output path (quirk #22b).")
    ap.add_argument("--out", default=None, help="Override the full output path.")
    args = ap.parse_args()

    table1_path = args.table1 or os.path.join(args.outroot, "hla_calls_rich.tsv")
    cohort_full_path = args.cohort_full or os.path.join(args.outroot, "immuannot_cohort_full.tsv")
    ancestry_preds_path = args.ancestry_preds or os.path.join(args.mount, DEFAULT_ANCESTRY_PREDS)
    hla_genotypes_path = args.hla_genotypes or os.path.join(args.mount, DEFAULT_HLA_GENOTYPES)

    if not args.skip_mount_check:
        # Only actually required when we're about to read files that live under the mount.
        needs_mount = (args.ancestry_preds is None or args.hla_genotypes is None)
        if needs_mount:
            check_mount(args.mount, args.billing_project)

    print(f"Loading Table 1 from {table1_path!r} ...", file=sys.stderr)
    table1 = load_table1(table1_path)

    print(f"Loading cohort metadata from {cohort_full_path!r} ...", file=sys.stderr)
    cohort_full = load_cohort_full(cohort_full_path)

    print(f"Loading ancestry predictions from {ancestry_preds_path!r} ...", file=sys.stderr)
    ancestry = load_ancestry_preds(ancestry_preds_path)

    print(f"Loading AoU-native short-read HLA calls from {hla_genotypes_path!r} ...",
          file=sys.stderr)
    sr_ids = load_hla_genotypes(hla_genotypes_path)
    sr_id_set = set(sr_ids["person_id"])

    lr_stats = build_per_person_lr_stats(table1)
    lr_id_set = set(table1["person_id"].unique())

    # Universe of people: everyone who appears in any of the LR pipeline's own bookkeeping
    # (cohort_full) union everyone who has LR calls at all -- cohort_full is expected to be a
    # superset (it also lists people with zero output), but union defensively rather than assuming.
    all_ids = sorted(set(cohort_full["person_id"]) | lr_id_set)
    base = pd.DataFrame({"person_id": all_ids})

    base["in_lr"] = base["person_id"].isin(lr_id_set)
    base["in_sr"] = base["person_id"].isin(sr_id_set)

    base = base.merge(lr_stats, on="person_id", how="left")
    base = base.merge(
        cohort_full[[c for c in ["person_id", "platform", "trim_tier"] if c in cohort_full.columns]],
        on="person_id", how="left")

    # Ancestry: prefer the AoU short-read prediction (continuous admixture available), fall back
    # to the LR pipeline's own recorded ancestry_pred if a person has no SR ancestry row at all.
    base = base.merge(ancestry, on="person_id", how="left")
    if "ancestry_pred" in cohort_full.columns:
        fallback = cohort_full.set_index("person_id")["ancestry_pred"]
        base["ancestry_pred"] = base["ancestry_pred_sr"].where(
            base["ancestry_pred_sr"].notna(), base["person_id"].map(fallback))
    else:
        base["ancestry_pred"] = base["ancestry_pred_sr"]
    base = base.drop(columns=["ancestry_pred_sr"])

    base["n_genes_called"] = base["n_genes_called"].fillna(0).astype(int)
    base["n_genes_exact"] = base["n_genes_exact"].fillna(0).astype(int)
    base["td_stratum"] = base["max_template_distance"].apply(compute_td_stratum)

    prob_cols = [f"p_{a.lower()}" for a in ANCESTRY_ORDER]
    columns = (["person_id", "in_sr", "in_lr", "max_template_distance",
                "mean_template_distance", "n_genes_called", "n_genes_exact", "td_stratum",
                "ancestry_pred"] + prob_cols + ["platform", "trim_tier"])
    for c in columns:
        if c not in base.columns:
            base[c] = pd.NA
    base = base[columns]

    suffix = ".sample" if args.sample else ""
    out_path = args.out or os.path.join(args.outroot, f"cohort_membership{suffix}.tsv")
    base.to_csv(out_path, sep="\t", index=False)

    n = len(base)
    print(f"\nWrote {n} rows to {out_path!r}.", file=sys.stderr)
    print(f"  in_lr: {base['in_lr'].sum()} ({base['in_lr'].mean():.1%})", file=sys.stderr)
    print(f"  in_sr: {base['in_sr'].sum()} ({base['in_sr'].mean():.1%})", file=sys.stderr)
    print(f"  in both: {(base['in_lr'] & base['in_sr']).sum()}", file=sys.stderr)
    print(f"  td_stratum distribution:\n{base['td_stratum'].value_counts(dropna=False)}",
          file=sys.stderr)
    n_missing_ancestry = base["ancestry_pred"].isna().sum()
    if n_missing_ancestry:
        print(f"  WARNING: {n_missing_ancestry} people have no ancestry_pred from either source.",
              file=sys.stderr)


if __name__ == "__main__":
    main()
