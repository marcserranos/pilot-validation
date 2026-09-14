#!/usr/bin/env python3
"""Real assertions against 10_allele_ancestry_geometry.py, per SCHEMA.md hard rule #6 ("every
script here must run against the synthetic fixtures before it is handed to Marc").

Two groups of tests:
  1. Unit tests locking in the two invariants explicitly required for this script: (a) centroid
     computation has NO minimum-carrier-count floor -- an n=1 allele's centroid must exactly equal
     that one carrier's own admixture proportions, and (b) the bootstrap enrichment p-value is
     well-calibrated -- an allele whose carriers are a uniform random draw from the cohort must NOT
     come back significant.
  2. An end-to-end smoke test against tests/make_fixtures.py's output, run through
     01/02/03_novel_alleles.py (which now also writes novel_alleles_carriers.tsv), asserting the
     script exits 0 and produces every figure + the report for both --renorm variants, plus a
     structure-recovery check: novel alleles should skew toward AFR more than known alleles do
     (the fixture's own documented novel-rate skew, afr 0.22 vs eur 0.05).

Run:
    python3 scripts/hla_popgen/tests/make_fixtures.py --outroot /tmp/hla_fixtures_test -n 300
    python3 scripts/hla_popgen/tests/test_allele_geometry.py --fixtures /tmp/hla_fixtures_test
"""
import argparse
import importlib.util
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)
sys.path.insert(0, HLA_POPGEN_DIR)


def _load_module(filename, modname):
    path = os.path.join(HLA_POPGEN_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


geo = _load_module("10_allele_ancestry_geometry.py", "hla_popgen_allele_geometry")

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  PASS  {name}")
    else:
        msg = f"  FAIL  {name}" + (f" -- {detail}" if detail else "")
        print(msg)
        FAILURES.append(name)


ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]


def _toy_cohort(n=60, seed=0):
    rng = np.random.default_rng(seed)
    probs = rng.dirichlet(np.ones(6), size=n)
    data = {f"p_{a.lower()}": probs[:, i] for i, a in enumerate(ANCESTRY_ORDER)}
    data["person_id"] = [f"p{i}" for i in range(n)]
    data["ancestry_pred"] = [ANCESTRY_ORDER[i] for i in probs.argmax(axis=1)]
    return pd.DataFrame(data)


# ---------------------------------------------------------------------------
# 1a. No minimum-carrier-count floor -- n=1 centroid == that person's own proportions
# ---------------------------------------------------------------------------
def test_single_carrier_centroid_equals_own_row():
    cohort = _toy_cohort(n=40, seed=1)
    components = ["AFR", "EUR", "AMR"]
    one_person = cohort.iloc[[5]]["person_id"].tolist()  # exactly one carrier
    centroid, n_used = geo.compute_centroid(set(one_person), cohort, components,
                                             renorm="per_person")
    check("single-carrier centroid computed (not skipped/gated)", centroid is not None)
    check("single-carrier n_used == 1", n_used == 1, detail=str(n_used))
    if centroid is not None:
        row = cohort[cohort["person_id"] == one_person[0]].iloc[0]
        raw = np.array([row[f"p_{c.lower()}"] for c in components])
        expected = raw / raw.sum()
        actual = np.array([centroid[c] for c in components])
        check("single-carrier centroid == that person's own renormalized proportions "
              "(the 'no floor' requirement)",
              np.allclose(actual, expected, atol=1e-9),
              detail=f"actual={actual.tolist()} expected={expected.tolist()}")

    # per_centroid variant: with exactly one carrier, averaging-then-renormalizing degenerates to
    # the same thing as renormalize-then-average (only one term in the average either way).
    centroid_alt, n_used_alt = geo.compute_centroid(set(one_person), cohort, components,
                                                     renorm="per_centroid")
    check("single-carrier centroid (per_centroid variant) also computed with no floor",
          centroid_alt is not None and n_used_alt == 1)
    if centroid is not None and centroid_alt is not None:
        check("both renorm variants agree exactly at n=1",
              np.allclose([centroid[c] for c in components],
                          [centroid_alt[c] for c in components], atol=1e-9))


def test_empty_carrier_set_returns_none_not_a_crash():
    cohort = _toy_cohort(n=20, seed=2)
    centroid, n_used = geo.compute_centroid(set(), cohort, ["AFR", "EUR", "AMR"])
    check("empty carrier set -> centroid is None (only case where it's ever skipped)",
          centroid is None)
    check("empty carrier set -> n_used == 0", n_used == 0)


def test_nan_centroid_gives_na_pvalue_not_spurious_zero():
    """Regression test for a real bug found against real cohort data (2026-09-06): a carrier
    whose admixture mass falls entirely OUTSIDE the chosen component subset (e.g. a 3-way AFR/EUR/
    AMR ternary for someone who is 100% EAS/MID/SAS) correctly gets a NaN centroid (0/0 in the
    per_person renormalization) -- but `null_dists >= nan` is elementwise False everywhere, so a
    naive `.mean()` silently returned 0.0, a spuriously 'significant' p-value on an undefined
    observation. bootstrap_enrichment must return NA for both z_score and p_value here, never 0."""
    components = ["AFR", "EUR", "AMR"]
    cohort = _toy_cohort(n=100, seed=7)
    # Force one person's AFR/EUR/AMR to sum to exactly 0 (all their mass in EAS/MID/SAS).
    cohort.loc[0, ["p_afr", "p_eur", "p_amr"]] = 0.0
    remaining = 1.0
    cohort.loc[0, ["p_eas", "p_mid", "p_sas"]] = [remaining / 3] * 3
    only_that_person = {cohort.loc[0, "person_id"]}

    centroid, n_used = geo.compute_centroid(only_that_person, cohort, components)
    check("degenerate-subset carrier still gets n_used == 1 (no floor)", n_used == 1)
    check("degenerate-subset centroid is NaN (0/0 in this component subset -- correct, not a bug)",
          all(np.isnan(centroid[c]) for c in components), detail=str(centroid))

    cohort_mean, _ = geo.compute_centroid(set(cohort["person_id"]), cohort, components)
    null = geo.bootstrap_null(cohort, n_used, components, n_boot=200, seed=0)
    result = geo.bootstrap_enrichment(centroid, null, cohort_mean, components)
    check("bootstrap p_value is NA (not spuriously 0.0) when the observed centroid is NaN",
          result["p_value"] is None or (isinstance(result["p_value"], float)
                                         and np.isnan(result["p_value"])),
          detail=str(result))
    check("bootstrap z_score is also NA in this case", np.isnan(result["z_score"]), detail=str(result))


# ---------------------------------------------------------------------------
# 1b. Bootstrap calibration -- a uniform-random-draw allele should NOT come back significant
# ---------------------------------------------------------------------------
def test_bootstrap_null_well_calibrated_on_random_allele():
    cohort = _toy_cohort(n=400, seed=3)
    components = ["AFR", "EUR", "AMR"]
    cohort_mean, _ = geo.compute_centroid(set(cohort["person_id"]), cohort, components)

    rng = np.random.default_rng(42)
    p_values = []
    for trial in range(30):
        carrier_ids = set(rng.choice(cohort["person_id"], size=15, replace=False))
        centroid, n_used = geo.compute_centroid(carrier_ids, cohort, components)
        null = geo.bootstrap_null(cohort, n_used, components, n_boot=500, seed=trial)
        result = geo.bootstrap_enrichment(centroid, null, cohort_mean, components)
        p_values.append(result["p_value"])

    p_values = np.array(p_values)
    frac_significant = float((p_values < 0.05).mean())
    check("bootstrap p-value calibration: a uniform-random-draw allele is flagged 'significant' "
          "(p<0.05) at roughly the nominal 5% rate across repeated trials, not systematically",
          frac_significant <= 0.25,
          detail=f"{frac_significant:.2%} of {len(p_values)} random-draw trials had p<0.05")


def test_bootstrap_flags_genuinely_skewed_allele():
    cohort = _toy_cohort(n=400, seed=4)
    components = ["AFR", "EUR", "AMR"]
    cohort_mean, _ = geo.compute_centroid(set(cohort["person_id"]), cohort, components)

    afr_heavy = cohort.sort_values("p_afr", ascending=False).head(15)["person_id"]
    centroid, n_used = geo.compute_centroid(set(afr_heavy), cohort, components)
    null = geo.bootstrap_null(cohort, n_used, components, n_boot=1000, seed=0)
    result = geo.bootstrap_enrichment(centroid, null, cohort_mean, components)
    check("bootstrap enrichment DOES flag a genuinely AFR-skewed carrier set as significant",
          result["p_value"] < 0.05, detail=str(result))


# ---------------------------------------------------------------------------
# 2. End-to-end smoke test against the real fixture pipeline
# ---------------------------------------------------------------------------
def test_end_to_end_smoke(fixtures_dir):
    py = sys.executable
    people_root = os.path.join(fixtures_dir, "people")
    t1 = os.path.join(fixtures_dir, "hla_calls_rich.sample.tsv")
    t4 = os.path.join(fixtures_dir, "cohort_membership.sample.tsv")

    with tempfile.TemporaryDirectory() as tmp:
        reports_dir = os.path.join(tmp, "reports")
        seqs_path = os.path.join(reports_dir, "novel_alleles_seqs.sample.fa")
        table3_path = os.path.join(reports_dir, "novel_alleles.sample.tsv")

        # NOTE: 03_novel_alleles.py is intentionally called exactly as it exists on main (no
        # --carriers-path, no changes) -- 00-04 is a different concurrent agent's ownership
        # boundary (see 10_allele_ancestry_geometry.py's module docstring); this test only ever
        # reads novel_alleles.tsv as a data file and independently recomputes carrier identity from
        # Table 1 + raw cds.fa.gz, exactly like the script under test does.
        r3 = subprocess.run([py, os.path.join(HLA_POPGEN_DIR, "03_novel_alleles.py"),
                              "--outroot", people_root, "--table1", t1,
                              "--cohort-membership", t4, "--out-dir", reports_dir,
                              "--seqs-path", seqs_path, "--sample"],
                             capture_output=True, text=True)
        check("03_novel_alleles.py exits 0 (prerequisite for the smoke test)", r3.returncode == 0,
              detail=r3.stderr[-500:])
        if r3.returncode != 0:
            return

        for renorm in ("per_person", "per_centroid"):
            out_dir = os.path.join(tmp, f"geo_{renorm}")
            r = subprocess.run([py, os.path.join(HLA_POPGEN_DIR, "10_allele_ancestry_geometry.py"),
                                 "--outroot", fixtures_dir, "--table1", t1,
                                 "--cohort-membership", t4, "--novel-alleles", table3_path,
                                 "--people-root", people_root, "--out-dir", out_dir,
                                 "--cohort", "lr", "--gene", "B", "--renorm", renorm],
                                capture_output=True, text=True)
            check(f"10_allele_ancestry_geometry.py exits 0 (--renorm {renorm})", r.returncode == 0,
                  detail=r.stderr[-800:])
            if r.returncode != 0:
                continue
            for fname in ("ternary_B.png", "tetrahedron_B.png", "pca_centroids_B.png",
                          "allele_ancestry_report.md"):
                check(f"({renorm}) {fname} written", os.path.exists(os.path.join(out_dir, fname)))

        # Structure-recovery check, run once (renorm choice doesn't matter for this direction):
        # find a gene with at least one novel allele in this fixture set and confirm novel alleles
        # skew more AFR than known alleles do, per the fixture's own documented novel-rate skew
        # (afr 0.22 vs eur 0.05, tests/make_fixtures.py's own docstring).
        table3 = pd.read_csv(table3_path, sep="\t")
        genes_with_novel = table3["gene"].unique().tolist()
        recovered = False
        for gene_prefixed in genes_with_novel:
            gene_bare = gene_prefixed.replace("HLA-", "")
            out_dir = os.path.join(tmp, f"geo_struct_{gene_bare}")
            r = subprocess.run([py, os.path.join(HLA_POPGEN_DIR, "10_allele_ancestry_geometry.py"),
                                 "--outroot", fixtures_dir, "--table1", t1,
                                 "--cohort-membership", t4, "--novel-alleles", table3_path,
                                 "--people-root", people_root, "--out-dir", out_dir,
                                 "--cohort", "lr", "--gene", gene_bare],
                                capture_output=True, text=True)
            if r.returncode != 0:
                continue
            report = os.path.join(out_dir, "allele_ancestry_report.md")
            if not os.path.exists(report):
                continue
            df = _parse_ternary_table(report)
            if df is None or "novel" not in df["kind"].values or "known" not in df["kind"].values:
                continue
            afr_novel = df.loc[df["kind"] == "novel", "AFR"].mean()
            afr_known = df.loc[df["kind"] == "known", "AFR"].mean()
            if pd.isna(afr_novel) or pd.isna(afr_known):
                continue
            recovered = True
            check(f"structure recovery (gene {gene_bare}): novel-allele centroids skew more AFR "
                  f"than known-allele centroids (fixture's afr 0.22 vs eur 0.05 novel-rate skew)",
                  afr_novel > afr_known,
                  detail=f"mean AFR share: novel={afr_novel:.3f} known={afr_known:.3f}")
            break
        check("structure-recovery check ran against at least one gene with both known and "
              "novel alleles present", recovered)

        # --mode multi: combined all-genes-colored-by-gene ternary + ancestry-colored small
        # multiples grid, both from ONE shared data pass (see build_records_for_genes).
        out_dir = os.path.join(tmp, "geo_multi")
        r_multi = subprocess.run(
            [py, os.path.join(HLA_POPGEN_DIR, "10_allele_ancestry_geometry.py"),
             "--outroot", fixtures_dir, "--table1", t1, "--cohort-membership", t4,
             "--novel-alleles", table3_path, "--people-root", people_root,
             "--out-dir", out_dir, "--cohort", "lr", "--mode", "multi",
             "--multi-genes", "A,B,C,DRB1,DQB1,DPB1,DQA1,DPA1",
             "--grid-genes", "A,B,C,DRB1,DQB1,DPB1"],
            capture_output=True, text=True)
        check("10_allele_ancestry_geometry.py --mode multi exits 0", r_multi.returncode == 0,
              detail=r_multi.stderr[-800:])
        if r_multi.returncode == 0:
            for fname in ("ternary_all_genes.png", "ternary_gene_grid.png"):
                check(f"(multi) {fname} written", os.path.exists(os.path.join(out_dir, fname)))

        # --grid-genes not a subset of --multi-genes must fail loudly, not silently drop genes.
        r_bad = subprocess.run(
            [py, os.path.join(HLA_POPGEN_DIR, "10_allele_ancestry_geometry.py"),
             "--outroot", fixtures_dir, "--table1", t1, "--cohort-membership", t4,
             "--novel-alleles", table3_path, "--people-root", people_root,
             "--out-dir", os.path.join(tmp, "geo_multi_bad"), "--cohort", "lr", "--mode", "multi",
             "--multi-genes", "A,B", "--grid-genes", "A,B,C"],
            capture_output=True, text=True)
        check("--mode multi with --grid-genes not a subset of --multi-genes fails loudly",
              r_bad.returncode != 0, detail=r_bad.stderr[-300:])


def _parse_ternary_table(report_path):
    """Pull the 'Ternary-space alleles' markdown table back into a DataFrame with a plain AFR
    column, for the structure-recovery check above."""
    with open(report_path) as f:
        lines = f.readlines()
    start = None
    for i, line in enumerate(lines):
        if line.startswith("## Ternary-space alleles"):
            start = i
            break
    if start is None:
        return None
    rows = []
    for line in lines[start:]:
        if line.startswith("| allele") or line.startswith("|---"):
            continue
        if not line.startswith("|"):
            if rows:
                break
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 4:
            continue
        allele_id, kind, n_carriers, centroid = cells[0], cells[1], cells[2], cells[3]
        afr = None
        for part in centroid.split(","):
            part = part.strip()
            if part.startswith("AFR="):
                afr = float(part.split("=", 1)[1])
        rows.append({"allele_id": allele_id, "kind": kind, "AFR": afr})
    return pd.DataFrame(rows) if rows else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixtures", default="/tmp/hla_fixtures_test",
                    help="Path to fixtures built by tests/make_fixtures.py. Built fresh if missing.")
    args = ap.parse_args()

    if not os.path.isdir(args.fixtures):
        print(f"Building fixtures at {args.fixtures!r} ...", file=sys.stderr)
        spec = importlib.util.spec_from_file_location(
            "hla_popgen_make_fixtures", os.path.join(HERE, "make_fixtures.py"))
        make_fixtures = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(make_fixtures)
        make_fixtures.build(args.fixtures, 300)

    print("\n-- Unit tests: no-minimum-carrier-count-floor invariant --")
    test_single_carrier_centroid_equals_own_row()
    test_empty_carrier_set_returns_none_not_a_crash()

    print("\n-- Unit tests: bootstrap enrichment calibration --")
    test_bootstrap_null_well_calibrated_on_random_allele()
    test_bootstrap_flags_genuinely_skewed_allele()
    test_nan_centroid_gives_na_pvalue_not_spurious_zero()

    print("\n-- End-to-end smoke test against fixtures --")
    test_end_to_end_smoke(args.fixtures)

    print(f"\n{'ALL PASS' if not FAILURES else f'{len(FAILURES)} FAILURE(S)'}")
    sys.exit(1 if FAILURES else 0)


if __name__ == "__main__":
    main()
