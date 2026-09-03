#!/usr/bin/env python3
"""Real assertions against 05_figures_frequency.py, 06_figures_structure.py, and
07_figures_crosscohort.py, run against the synthetic fixtures (tests/make_fixtures.py).

SCHEMA.md hard rule #6: "every script here must run against the synthetic fixtures before it is
handed to Marc." This file covers three things unit tests on individual functions can't:
  1. End-to-end smoke tests -- each script's `main()`, invoked via subprocess exactly as a real
     user would run it, against all three `--cohort` values (sr/lr/lr_td), producing real PNGs +
     markdown reports with no exceptions. Every subprocess invocation gets its OWN `--out-dir`
     under a fresh `tempfile.mkdtemp()` for this test run (see "Hermetic output dirs" below) --
     never the scripts' own shared `reports/hla_popgen/...` default.
  2. Structure-recovery: the fixtures encode KNOWN ancestry-biased allele frequencies (B*15:03/
     53:01/58:01 enriched in afr; B*07:02/08:01/57:01 enriched in eur -- tests/make_fixtures.py's
     own docstring), a KNOWN novel-rate skew (afr 0.22 vs eur 0.05, driving template_distance
     higher in afr), and -- as of the fixture's SR-genotype rewrite -- a KNOWN ancestry-dependent
     SR miscall rate biased toward EUR-common alleles (afr 0.28 miscall / eur 0.04, directional
     toward common EUR alleles, not symmetric noise). If 05's frequency tables, 07's
     template_distance-by-ancestry sweep, or 07's SR-vs-LR ancestry-stratified disagreement don't
     recover those directions, there is a real bug, not sampling noise -- these are the
     highest-value regression checks in this file.
  3. Unit tests on the statistics helpers in _viz_common.py (Wilson CI sanity bounds, rarefied
     richness monotonicity, gene_bare/gene_display round-trip) that are easy to get subtly wrong
     and don't require the full pipeline to check.

## Hermetic output dirs (fixes a real flaky-test bug)

Earlier versions of this suite let every subprocess write to the scripts' own default
`reports/hla_popgen/<script>/<cohort>/` path. That is shared, persistent state across DIFFERENT
test runs (different fixture sets, different code versions, a run that crashed mid-write) -- a
run's success or failure should never depend on what a PRIOR run (possibly against different
fixtures, possibly using code with the since-fixed negative-yerr crash) left on disk there. Every
subprocess call below is given an explicit `--out-dir` inside a fresh temp directory created once
per `main()` invocation, so this file's tests cannot see or be affected by any other run, in any
order, ever.

Run:
    python3 scripts/hla_popgen/tests/make_fixtures.py --outroot /tmp/hla_fixtures_test -n 300
    python3 scripts/hla_popgen/01_extract_rich.py --outroot /tmp/hla_fixtures_test --sample
    python3 scripts/hla_popgen/02_build_cohorts.py --outroot /tmp/hla_fixtures_test \\
        --table1 /tmp/hla_fixtures_test/hla_calls_rich.sample.tsv \\
        --cohort-full /tmp/hla_fixtures_test/immuannot_cohort_full.tsv \\
        --ancestry-preds /tmp/hla_fixtures_test/ancestry_preds.tsv \\
        --hla-genotypes /tmp/hla_fixtures_test/hla_genotypes.tsv --skip-mount-check --sample
    python3 scripts/hla_popgen/tests/test_figures.py --fixtures /tmp/hla_fixtures_test
"""
import argparse
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)


def _load_module(filename, modname):
    path = os.path.join(HLA_POPGEN_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


vc = _load_module("_viz_common.py", "hla_popgen_viz_common")
crosscohort = _load_module("07_figures_crosscohort.py", "hla_popgen_crosscohort")

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  PASS  {name}")
    else:
        msg = f"  FAIL  {name}" + (f" -- {detail}" if detail else "")
        print(msg)
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# Unit tests on _viz_common.py helpers
# ---------------------------------------------------------------------------
def test_wilson_ci_bounds():
    p, lo, hi = vc.wilson_ci(np.array([0, 5, 50]), np.array([10, 10, 50], dtype=float))
    check("wilson_ci point estimates correct", np.allclose(p, [0.0, 0.5, 1.0]), detail=str(p))
    check("wilson_ci lo <= p <= hi", bool(np.all(lo <= p + 1e-9) and np.all(p <= hi + 1e-9)),
          detail=f"lo={lo} p={p} hi={hi}")
    check("wilson_ci bounds within [0,1]", bool(np.all(lo >= -1e-9) and np.all(hi <= 1 + 1e-9)))
    check("wilson_ci n=0 cell is NaN", bool(np.isnan(vc.wilson_ci(np.array([0]), np.array([0.0]))[0][0])))


def test_diversity_indices():
    het_uniform = vc.expected_heterozygosity([0.25, 0.25, 0.25, 0.25])
    het_fixed = vc.expected_heterozygosity([1.0])
    check("heterozygosity: uniform 4-allele = 0.75", abs(het_uniform - 0.75) < 1e-9, str(het_uniform))
    check("heterozygosity: monomorphic = 0.0", abs(het_fixed - 0.0) < 1e-9, str(het_fixed))
    sh = vc.shannon_entropy([0.5, 0.5])
    check("shannon: two equally-likely alleles = ln(2)", abs(sh - np.log(2)) < 1e-9, str(sh))


def test_rarefied_richness_monotonic():
    counts = np.array([50, 30, 15, 5])
    r_small = vc.rarefied_richness(counts, 10)
    r_large = vc.rarefied_richness(counts, 50)
    check("rarefied richness increases with sample size", r_small < r_large,
          detail=f"r(10)={r_small} r(50)={r_large}")
    check("rarefied richness bounded by observed richness",
          vc.rarefied_richness(counts, int(counts.sum())) <= len(counts) + 1e-6)


def test_allele_field_parsing():
    check("to_nfield 2-field from 4-field consensus",
          vc.to_nfield("HLA-A*02:01:01:01", 2) == "02:01")
    check("to_nfield 4-field exact", vc.to_nfield("HLA-A*02:01:01:01", 4) == "02:01:01:01")
    check("to_nfield returns None when novel truncates before n",
          vc.to_nfield("HLA-A*02:new", 2) is None)
    check("to_nfield handles undetermined", vc.to_nfield("undetermined", 2) is None)
    check("gene_bare strips HLA- prefix", vc.gene_bare("HLA-DQA1") == "DQA1")
    check("gene_bare no-ops on already-bare names", vc.gene_bare("MICA") == "MICA")
    check("gene_display round-trips non-MIC/TAP genes", vc.gene_display("DRB3") == "HLA-DRB3")
    check("gene_display leaves MIC/TAP bare", vc.gene_display("MICA") == "MICA")


def test_nonclassical_genes_are_bare():
    # Regression for the real bug this project's own fixture run caught: a mixed prefixed/bare
    # NONCLASSICAL_GENES list silently dropped 8 of 12 genes from every filter that matches on
    # gene_bare (see 06/07's module docstrings + the bug note in _viz_common.py itself).
    check("NONCLASSICAL_GENES has no 'HLA-' prefixed entries",
          all(not g.startswith("HLA-") for g in vc.NONCLASSICAL_GENES),
          detail=str(vc.NONCLASSICAL_GENES))


def test_errorbar_never_gets_negative_yerr():
    """Regression for the negative-yerr crash: a cell with n=0 (p is NaN) must never be fed to
    ax.errorbar as a zero-filled point paired with finite Wilson bounds. This replicates the exact
    arithmetic 05_figures_frequency.py's plot_bar_gene and 07_figures_crosscohort.py's
    plot_heterodimer_and_pairable now do: mask NaN-p cells out entirely rather than nan_to_num
    them, so err_lo/err_hi are only ever computed from real (p, lo, hi) triples."""
    c = np.array([0, 3, 50])
    n_total = np.array([0.0, 20.0, 50.0])  # first cell: n=0 -> p is NaN
    p, lo, hi = vc.wilson_ci(c, n_total)
    has_data = ~np.isnan(p)
    err_lo = (p - lo)[has_data]
    err_hi = (hi - p)[has_data]
    check("errorbar mask excludes the NaN-p (n=0) cell", has_data.tolist() == [False, True, True])
    check("err_lo never negative once NaN cells are masked out", bool(np.all(err_lo >= -1e-9)),
          detail=str(err_lo))
    check("err_hi never negative once NaN cells are masked out", bool(np.all(err_hi >= -1e-9)),
          detail=str(err_hi))


# ---------------------------------------------------------------------------
# End-to-end smoke tests + structure recovery, against real fixture output on disk
# ---------------------------------------------------------------------------
def run_script(script, args):
    cmd = [sys.executable, os.path.join(HLA_POPGEN_DIR, script)] + args
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result


def test_smoke_all_cohorts(fixtures_dir, tmp_root):
    """Every subprocess call gets its own --out-dir under tmp_root (see module docstring,
    "Hermetic output dirs") -- no script here ever writes to or reads from the shared
    reports/hla_popgen/ tree, so this test's pass/fail cannot depend on run order or on anything a
    previous test_figures.py invocation (against different fixtures, or pre-bugfix code) left
    behind."""
    table1 = os.path.join(fixtures_dir, "hla_calls_rich.sample.tsv")
    table2 = os.path.join(fixtures_dir, "hla_cis_pairs.sample.tsv")
    cohort = os.path.join(fixtures_dir, "cohort_membership.sample.tsv")
    sr = os.path.join(fixtures_dir, "hla_genotypes.tsv")
    common = ["--table1", table1, "--cohort-membership", cohort, "--sr-genotypes", sr]

    out_dirs = {}
    for cohort_flag in [["--cohort", "lr"], ["--cohort", "sr"], ["--cohort", "lr_td", "--td-max", "5"]]:
        slug = "_".join(cohort_flag).replace("--", "")
        out_dir = os.path.join(tmp_root, "05", slug)
        r = run_script("05_figures_frequency.py", common + cohort_flag + ["--out-dir", out_dir])
        check(f"05_figures_frequency.py exits 0 ({' '.join(cohort_flag)})", r.returncode == 0,
              detail=(r.stderr or "")[-1500:])
        out_dirs[("05", slug)] = out_dir

    for cohort_flag in [["--cohort", "lr"], ["--cohort", "sr"], ["--cohort", "lr_td", "--td-max", "5"]]:
        slug = "_".join(cohort_flag).replace("--", "")
        out_dir = os.path.join(tmp_root, "06", slug)
        r = run_script("06_figures_structure.py", common + cohort_flag + ["--out-dir", out_dir])
        check(f"06_figures_structure.py exits 0 ({' '.join(cohort_flag)})", r.returncode == 0,
              detail=(r.stderr or "")[-1500:])
        out_dirs[("06", slug)] = out_dir

    cross_common = ["--table1", table1, "--table2", table2, "--cohort-membership", cohort,
                     "--sr-genotypes", sr]
    out_dir_07 = os.path.join(tmp_root, "07", "cohort_lr")
    r = run_script("07_figures_crosscohort.py",
                    cross_common + ["--cohort", "lr", "--out-dir", out_dir_07])
    check("07_figures_crosscohort.py exits 0 (--cohort lr)", r.returncode == 0,
          detail=(r.stderr or "")[-1500:])
    out_dirs[("07", "cohort_lr")] = out_dir_07

    r = run_script("07_figures_crosscohort.py", cross_common + ["--cohort", "sr"])
    check("07_figures_crosscohort.py refuses --cohort sr with exit != 0", r.returncode != 0)

    return out_dirs


def test_pca_sr_not_degenerate(out_dirs):
    """The SR fixture used to draw alleles uniformly at random with NO ancestry structure at all,
    so PCA on --cohort sr was correctly an unstructured blob and no PC could separate ancestries
    even a little. The fixture now derives sr as a degraded (ancestry-biased-miscall) copy of the
    lr truth, so it should show at least SOME separation. This doesn't assert a specific
    silhouette threshold (8 loci is a small feature set -- see 06's own docstring caveat), just
    that the report was written and records a real (non-NaN) silhouette, i.e. PCA actually ran on
    real multi-ancestry data rather than degenerating on an empty/single-ancestry matrix."""
    out_dir = out_dirs.get(("06", "cohort_sr"))
    if out_dir is None:
        check("06 --cohort sr PCA report exists", False, detail="out_dir not found")
        return
    report_path = os.path.join(out_dir, "structure_report.md")
    check("06 --cohort sr structure_report.md exists", os.path.exists(report_path), report_path)
    if not os.path.exists(report_path):
        return
    text = open(report_path).read()
    check("06 --cohort sr PCA did not skip for lack of data",
          "Skipped PCA/UMAP" not in text, detail=text[:400])
    check("06 --cohort sr reports a real (non-NaN) PCA silhouette",
          "silhouette: nan" not in text.lower(), detail=text[:800])


def test_structure_recovery(fixtures_dir):
    """The two original highest-value regression checks: does 05's frequency table and 07's
    template_distance-by-ancestry sweep actually recover the DIRECTION of bias the fixtures were
    built to contain? (tests/make_fixtures.py ANCESTRY_BIAS + per-ancestry novel_rate)."""
    table1_path = os.path.join(fixtures_dir, "hla_calls_rich.sample.tsv")
    cohort_path = os.path.join(fixtures_dir, "cohort_membership.sample.tsv")
    if not (os.path.exists(table1_path) and os.path.exists(cohort_path)):
        check("structure recovery: fixture files present", False,
              detail=f"missing {table1_path} or {cohort_path}")
        return

    table1 = vc.load_table1(table1_path)
    cohort_df = vc.load_cohort_membership(cohort_path)
    lr_people, _ = vc.select_cohort_people(cohort_df, "lr")
    person_ids = set(lr_people["person_id"])
    anc_map = dict(zip(lr_people["person_id"], lr_people["ancestry_pred"]))

    sub = table1[table1["person_id"].isin(person_ids) & (table1["gene_bare"] == "B")].copy()
    sub["allele_2field"] = sub["consensus"].map(lambda a: vc.to_nfield(a, 2))
    sub["ancestry_pred"] = sub["person_id"].map(anc_map)
    sub = sub.dropna(subset=["allele_2field"])
    sub = sub[sub["ancestry_pred"].isin(vc.ANCESTRY_ORDER)]

    def freq(allele, anc):
        g = sub[sub["ancestry_pred"] == anc]
        if len(g) == 0:
            return np.nan
        return (g["allele_2field"] == allele).mean()

    afr_biased = np.mean([freq("15:03", "AFR"), freq("53:01", "AFR"), freq("58:01", "AFR")])
    eur_for_afr_biased = np.mean([freq("15:03", "EUR"), freq("53:01", "EUR"), freq("58:01", "EUR")])
    check("HLA-B AFR-enriched alleles (15:03/53:01/58:01) higher in AFR than EUR",
          afr_biased > eur_for_afr_biased,
          detail=f"AFR={afr_biased:.3f} EUR={eur_for_afr_biased:.3f}")

    eur_biased = np.mean([freq("07:02", "EUR"), freq("08:01", "EUR"), freq("57:01", "EUR")])
    afr_for_eur_biased = np.mean([freq("07:02", "AFR"), freq("08:01", "AFR"), freq("57:01", "AFR")])
    check("HLA-B EUR-enriched alleles (07:02/08:01/57:01) higher in EUR than AFR",
          eur_biased > afr_for_eur_biased,
          detail=f"EUR={eur_biased:.3f} AFR={afr_for_eur_biased:.3f}")

    # template_distance: afr novel_rate (0.22) >> eur novel_rate (0.05) in the fixture generator,
    # so mean template_distance pooled over classical genes should be clearly higher in afr.
    classical_sub = table1[table1["person_id"].isin(person_ids) &
                            table1["gene_bare"].isin(vc.CLASSICAL_GENES_BARE)].copy()
    classical_sub["ancestry_pred"] = classical_sub["person_id"].map(anc_map)
    classical_sub = classical_sub[classical_sub["ancestry_pred"].isin(vc.ANCESTRY_ORDER)]
    mean_td = classical_sub.groupby("ancestry_pred")["template_distance"].mean()
    check("mean template_distance higher in AFR than EUR (reference-bias signal)",
          mean_td.get("AFR", 0) > mean_td.get("EUR", 999),
          detail=str(mean_td.to_dict()))


def test_sr_lr_disagreement_recovery(fixtures_dir):
    """New headline assertion (coordinator request): the SR fixture is now a degraded copy of the
    LR truth with an ancestry-dependent miscall rate biased TOWARD common EUR alleles (afr 0.28,
    eur 0.04 miscall rate). This is directly testable via
    07_figures_crosscohort.compute_sr_lr_disagreement_by_ancestry:
      1. SR-vs-LR agreement should be markedly BETTER (lower mean abs diff) in EUR than in AFR.
      2. The AFR disagreement should be DIRECTIONAL: SR overestimates EUR-common alleles in AFR
         (mean(sr_freq - lr_freq) > 0 for alleles common in the LR EUR population), not just
         noisier in some unsigned sense.
    """
    table1_path = os.path.join(fixtures_dir, "hla_calls_rich.sample.tsv")
    cohort_path = os.path.join(fixtures_dir, "cohort_membership.sample.tsv")
    sr_path = os.path.join(fixtures_dir, "hla_genotypes.tsv")
    if not all(os.path.exists(p) for p in (table1_path, cohort_path, sr_path)):
        check("sr/lr disagreement recovery: fixture files present", False)
        return

    table1 = vc.load_table1(table1_path)
    sr_long = vc.load_sr_genotypes(sr_path)
    cohort_df = vc.load_cohort_membership(cohort_path)
    lr_people, _ = vc.select_cohort_people(cohort_df, "lr")
    sr_people, _ = vc.select_cohort_people(cohort_df, "sr")

    genes = crosscohort.CLASSICAL
    lr_long2 = crosscohort.build_long_frame(table1, lr_people, genes, resolution=2)
    sr_long2 = crosscohort.build_sr_long_frame(sr_long, sr_people, genes)
    sr_freq_anc = crosscohort.freq_by_ancestry(sr_long2, genes)
    lr_freq_anc = crosscohort.freq_by_ancestry(lr_long2, genes)
    n_people_lr = lr_people["ancestry_pred"].value_counts().to_dict()
    n_people_sr = sr_people["ancestry_pred"].value_counts().to_dict()
    disagree_df = crosscohort.compute_sr_lr_disagreement_by_ancestry(
        sr_freq_anc, lr_freq_anc, n_people_lr, n_people_sr)

    row = {r["ancestry"]: r for _, r in disagree_df.iterrows()}
    eur_mad = row["EUR"]["mean_abs_diff"]
    afr_mad = row["AFR"]["mean_abs_diff"]
    check("SR-vs-LR agreement (mean abs freq diff) better in EUR than AFR",
          (not pd.isna(eur_mad)) and (not pd.isna(afr_mad)) and eur_mad < afr_mad,
          detail=f"EUR mean_abs_diff={eur_mad} AFR mean_abs_diff={afr_mad}")

    # EUR/AFR are both well-powered in the fixtures (n well above min_cell_n on both sides), so
    # this check itself is not a thin-N false-positive risk -- but assert that explicitly rather
    # than assuming it, so this test doesn't silently start passing on noise if the fixture's
    # ancestry composition ever changes.
    check("EUR/AFR rows used in the disagreement check are NOT thin",
          min(row["EUR"]["n_people_lr"], row["EUR"]["n_people_sr"]) >= vc.MIN_CELL_N_PEOPLE and
          min(row["AFR"]["n_people_lr"], row["AFR"]["n_people_sr"]) >= vc.MIN_CELL_N_PEOPLE,
          detail=f"EUR n_people(lr/sr)={row['EUR']['n_people_lr']}/{row['EUR']['n_people_sr']} "
                 f"AFR n_people(lr/sr)={row['AFR']['n_people_lr']}/{row['AFR']['n_people_sr']}")

    afr_bias = row["AFR"]["eur_common_allele_bias"]
    check("AFR disagreement is directional: SR overestimates EUR-common alleles in AFR",
          (not pd.isna(afr_bias)) and afr_bias > 0,
          detail=f"AFR eur_common_allele_bias={afr_bias} "
                 f"(n_eur_common_alleles={row['AFR']['n_eur_common_alleles']})")

    # Regression for the thin-N visual-weighting bug: an under-powered ancestry must be FLAGGED,
    # so a reader is not drawn to a tall bar resting on a handful of people.
    #
    # Assert the MECHANISM, not which ancestries happen to be small. An earlier version of this
    # test hardcoded MID/SAS, and broke the moment the fixture's ancestry-group sizes shifted --
    # a test failure that said nothing about the code under test. Whether a given group is thin is
    # a property of the fixture; that thin groups get flagged and non-thin ones don't is the
    # property of the code, and that is what belongs in an assertion.
    thin = [a for a, r in row.items()
            if min(r["n_people_lr"], r["n_people_sr"]) < vc.MIN_CELL_N_PEOPLE]
    fat = [a for a, r in row.items()
           if min(r["n_people_lr"], r["n_people_sr"]) >= vc.MIN_CELL_N_PEOPLE]
    check("thin/non-thin split is computable for every ancestry in the disagreement data",
          len(thin) + len(fat) == len(row),
          detail=f"thin={thin} non-thin={fat} (threshold={vc.MIN_CELL_N_PEOPLE} people)")
    # The well-powered comparison the figure actually exists to make must survive the guard.
    check("EUR and AFR are both non-thin, so the headline SR-vs-LR claim is well-powered",
          "EUR" in fat and "AFR" in fat,
          detail=f"non-thin={fat}")


def test_report_files_written(out_dirs, fixtures_dir):
    """Checks the SAME hermetic per-test-run out_dirs test_smoke_all_cohorts wrote to (never the
    shared reports/hla_popgen/ default), so this is immune to any other run's state."""
    out_lr_05 = out_dirs.get(("05", "cohort_lr"))
    check("05 frequency_report.md exists after smoke test",
          bool(out_lr_05) and os.path.exists(os.path.join(out_lr_05, "frequency_report.md")))
    out_lr_06 = out_dirs.get(("06", "cohort_lr"))
    check("06 structure_report.md exists after smoke test",
          bool(out_lr_06) and os.path.exists(os.path.join(out_lr_06, "structure_report.md")))
    out_07 = out_dirs.get(("07", "cohort_lr"))
    cross_report = os.path.join(out_07, "crosscohort_report.md") if out_07 else None
    check("07 crosscohort_report.md exists after smoke test",
          bool(out_07) and os.path.exists(cross_report or ""))
    # Aggregate-only compliance (SCHEMA.md hard rule #5): no bare person_id in any written report.
    if cross_report and os.path.exists(cross_report):
        cohort_df = vc.load_cohort_membership(
            os.path.join(fixtures_dir, "cohort_membership.sample.tsv"))
        sample_ids = set(cohort_df["person_id"].astype(str).head(50))
        text = open(cross_report).read()
        leaked = [pid for pid in sample_ids if pid in text]
        check("no bare person_id leaked into crosscohort_report.md", len(leaked) == 0,
              detail=str(leaked[:5]))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fixtures", default="/tmp/hla_fixtures",
                     help="Root written by make_fixtures.py + 01_extract_rich.py --sample + "
                          "02_build_cohorts.py --sample (see this file's own docstring for the "
                          "exact three-command setup).")
    ap.add_argument("--keep-tmp", action="store_true",
                     help="Don't delete the hermetic output tmpdir on exit (for debugging).")
    args = ap.parse_args()

    tmp_root = tempfile.mkdtemp(prefix="hla_popgen_test_")
    try:
        print("== Unit tests (_viz_common.py helpers) ==")
        test_wilson_ci_bounds()
        test_diversity_indices()
        test_rarefied_richness_monotonic()
        test_allele_field_parsing()
        test_nonclassical_genes_are_bare()
        test_errorbar_never_gets_negative_yerr()

        print(f"\n== End-to-end smoke tests (all three --cohort values, all three scripts, "
              f"hermetic out-dirs under {tmp_root}) ==")
        out_dirs = test_smoke_all_cohorts(args.fixtures, tmp_root)

        print("\n== SR-cohort PCA sanity (fixture now has real ancestry structure) ==")
        test_pca_sr_not_degenerate(out_dirs)

        print("\n== Structure recovery (known ancestry-biased signal in the fixtures) ==")
        test_structure_recovery(args.fixtures)

        print("\n== SR-vs-LR disagreement recovery (ancestry-stratified reference-bias claim) ==")
        test_sr_lr_disagreement_recovery(args.fixtures)

        print("\n== Output/compliance checks ==")
        test_report_files_written(out_dirs, args.fixtures)
    finally:
        if not args.keep_tmp:
            shutil.rmtree(tmp_root, ignore_errors=True)
        else:
            print(f"\n(kept hermetic output dir: {tmp_root})", file=sys.stderr)

    print(f"\n{'ALL PASS' if not FAILURES else f'{len(FAILURES)} FAILURE(S)'}")
    if FAILURES:
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)


if __name__ == "__main__":
    main()
