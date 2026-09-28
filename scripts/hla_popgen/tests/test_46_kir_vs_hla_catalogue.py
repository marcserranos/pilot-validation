#!/usr/bin/env python3
"""Unit tests for 46_kir_vs_hla_catalogue.py, synthetic data only.

Covers:
  1. `build_gene_metrics()` -- pct_novel_any/pct_novel_protein/completeness arithmetic on a small
     hand-computable synthetic `coverage_chao2.tsv`-shaped frame, including the merge against
     `equal_n_slope.tsv`-shaped slope data.
  2. `two_proportion_ztest()` reuse via the lazily-imported 45 module -- a known textbook z-test
     value (Marc's own worked hypothesis-test numbers from the S04 brief) reproduces to 2 decimal
     places.
  3. Figure rendering on synthetic data through `_viz_common.check_layout(strict=True)` (via
     `save_fig`) -- both figures must render with zero layout violations.
  4. **Row label <-> plotted value <-> TSV binding test** (2026-09-28 Cleveland-dot-plot rebuild):
     `TestBuildRowOrder` checks `build_row_order()`'s gene ordering against an order independently
     recomputed straight from the metrics table; `TestRowLabelBindingMatchesData` renders panel a
     and, for every row, confirms the shared y-axis tick label, the metrics table's own
     `gene_display` for that (gene, species), and the actual protein-level marker's (x, y) all
     agree -- a real regression guard against a mislabeled row (e.g. a sort that disagrees between
     the row builder and the plotting code, or a gene/species mismatch upstream).

Run: python3 scripts/hla_popgen/tests/test_46_kir_vs_hla_catalogue.py
"""
import importlib.util
import math
import os
import sys
import tempfile
import unittest

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)


def _load_module(filename, modname):
    spec = importlib.util.spec_from_file_location(modname, os.path.join(HLA_POPGEN_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


m46 = _load_module("46_kir_vs_hla_catalogue.py", "kir_vs_hla_catalogue_test_target")


def _synthetic_cov():
    """Two genes per species, genomic + any_novel + protein_novel levels, ancestry ALL only --
    hand-computable pct_novel/completeness."""
    rows = []

    def add(gene, species, level, s_obs, gt=0.99, chao2=None):
        rows.append({"gene": gene, "ancestry": "ALL", "level": level, "species": species,
                     "n_people": 1000, "s_obs": s_obs, "good_turing_coverage": gt,
                     "chao2": chao2 if chao2 is not None else s_obs * 1.5,
                     "chao2_se": 1.0, "chao2_undetected_f0hat": (chao2 or s_obs * 1.5) - s_obs,
                     "q1": "", "q2": "", "chao_new_by_2n": 0.0})

    # HLA-A: 100 genomic, 20 any-level novel (20%), 0 protein-novel (of 90 protein), chao2=200 ->
    # completeness 0.5. `protein` s_obs (90) is a SEPARATE granularity from `genomic` (100) on
    # purpose -- pct_novel_protein's fixed (2026-09-27) denominator is `protein`, not `genomic`.
    add("A", "hla", "genomic", 100, chao2=200.0)
    add("A", "hla", "any_novel", 20, chao2=40.0)
    add("A", "hla", "protein", 90, chao2=95.0)
    add("A", "hla", "protein_novel", 0, chao2=0.0)
    # HLA-B: 200 genomic, 100 any-level novel (50%), chao2=250 -> completeness 0.8
    add("B", "hla", "genomic", 200, chao2=250.0)
    add("B", "hla", "any_novel", 100, chao2=180.0)
    add("B", "hla", "protein", 150, chao2=170.0)
    add("B", "hla", "protein_novel", 0, chao2=0.0)
    # KIR2DL1: 50 genomic, 40 any-level novel (80%), chao2=55.55... -> completeness 0.9.
    # protein=50, protein_novel=10 -> pct_novel_protein = 20% (of `protein`, not `genomic`).
    add("KIR2DL1", "kir", "genomic", 50, chao2=50 / 0.9)
    add("KIR2DL1", "kir", "any_novel", 40, chao2=44.0)
    add("KIR2DL1", "kir", "protein", 50, chao2=60.0)
    add("KIR2DL1", "kir", "protein_novel", 10, chao2=12.0)
    # KIR2DL2: 40 genomic, 10 any-level novel (25%), chao2=80 -> completeness 0.5.
    # protein=20, protein_novel=2 -> pct_novel_protein = 10%.
    add("KIR2DL2", "kir", "genomic", 40, chao2=80.0)
    add("KIR2DL2", "kir", "any_novel", 10, chao2=15.0)
    add("KIR2DL2", "kir", "protein", 20, chao2=25.0)
    add("KIR2DL2", "kir", "protein_novel", 2, chao2=3.0)
    return pd.DataFrame(rows)


def _synthetic_slope():
    rows = []
    for gene, species, slope in [("A", "hla", 5.0), ("B", "hla", 8.0),
                                  ("KIR2DL1", "kir", 1.0), ("KIR2DL2", "kir", 0.5)]:
        rows.append({"gene": gene, "ancestry": "ALL", "level": "genomic", "species": species,
                     "curve": "distinct", "n_star": 1000, "slope_per_1000": slope,
                     "mean_at_n_star": 0.0})
    return pd.DataFrame(rows)


class TestBuildGeneMetrics(unittest.TestCase):
    def setUp(self):
        self.cov = _synthetic_cov()
        self.slope = _synthetic_slope()
        self.metrics = m46.build_gene_metrics(self.cov, self.slope)

    def test_pct_novel_any(self):
        row = self.metrics[(self.metrics.gene == "A") & (self.metrics.species == "hla")].iloc[0]
        self.assertAlmostEqual(row["pct_novel_any"], 20.0)
        row = self.metrics[(self.metrics.gene == "B") & (self.metrics.species == "hla")].iloc[0]
        self.assertAlmostEqual(row["pct_novel_any"], 50.0)
        row = self.metrics[(self.metrics.gene == "KIR2DL1") & (self.metrics.species == "kir")].iloc[0]
        self.assertAlmostEqual(row["pct_novel_any"], 80.0)

    def test_pct_novel_protein_zero_for_hla(self):
        for gene in ("A", "B"):
            row = self.metrics[(self.metrics.gene == gene) & (self.metrics.species == "hla")].iloc[0]
            self.assertAlmostEqual(row["pct_novel_protein"], 0.0)

    def test_completeness(self):
        row = self.metrics[(self.metrics.gene == "A") & (self.metrics.species == "hla")].iloc[0]
        self.assertAlmostEqual(row["completeness"], 0.5)
        row = self.metrics[(self.metrics.gene == "B") & (self.metrics.species == "hla")].iloc[0]
        self.assertAlmostEqual(row["completeness"], 0.8)
        row = self.metrics[(self.metrics.gene == "KIR2DL1") & (self.metrics.species == "kir")].iloc[0]
        self.assertAlmostEqual(row["completeness"], 0.9, places=6)
        row = self.metrics[(self.metrics.gene == "KIR2DL2") & (self.metrics.species == "kir")].iloc[0]
        self.assertAlmostEqual(row["completeness"], 0.5)

    def test_slope_merged(self):
        row = self.metrics[(self.metrics.gene == "A") & (self.metrics.species == "hla")].iloc[0]
        self.assertAlmostEqual(row["slope_per_1000"], 5.0)

    def test_gene_display(self):
        row = self.metrics[(self.metrics.gene == "A") & (self.metrics.species == "hla")].iloc[0]
        self.assertEqual(row["gene_display"], "HLA-A")
        row = self.metrics[(self.metrics.gene == "KIR2DL1") & (self.metrics.species == "kir")].iloc[0]
        self.assertEqual(row["gene_display"], "KIR2DL1")


class TestPctNovelProteinDenominator(unittest.TestCase):
    """2026-09-27 fix (S04 coordinator follow-up on 4a75657): pct_novel_protein's denominator must
    be `protein` s_obs (the SAME granularity as its own `protein_novel` numerator), not `genomic`
    s_obs -- dividing across two different identity granularities was itself part of why
    pct_novel_protein could exceed 100% even with 44's per-call identity fix in place."""

    def setUp(self):
        self.cov = _synthetic_cov()
        self.slope = _synthetic_slope()
        self.metrics = m46.build_gene_metrics(self.cov, self.slope)

    def test_pct_novel_protein_uses_protein_not_genomic_denominator(self):
        row = self.metrics[(self.metrics.gene == "KIR2DL1") & (self.metrics.species == "kir")].iloc[0]
        # protein_novel=10, protein=50 -> 20%. Dividing by genomic (50) would coincidentally also
        # give 20% here (genomic==protein by fixture construction) -- KIR2DL2 below disambiguates.
        self.assertAlmostEqual(row["pct_novel_protein"], 20.0)
        row = self.metrics[(self.metrics.gene == "KIR2DL2") & (self.metrics.species == "kir")].iloc[0]
        # protein_novel=2, protein=20 (NOT genomic=40) -> 10%, not 5% (which dividing by genomic=40
        # would give) -- this is the case that actually distinguishes the two denominators.
        self.assertAlmostEqual(row["pct_novel_protein"], 10.0)

    def test_pct_novel_protein_never_exceeds_100_when_protein_novel_is_a_true_subset(self):
        # protein_novel <= protein by construction in every fixture row here (44's own invariant:
        # a call is only ever added to protein_novel after being added to protein) -> the fixed
        # ratio must never exceed 100%, unlike the old genomic-denominator version.
        self.assertTrue((self.metrics["pct_novel_protein"] <= 100.0).all())

    def test_missing_protein_level_yields_nan_not_a_fabricated_value(self):
        """A gene whose `protein`/`protein_novel` rows are 44's literal "NA" (catalogue-uncovered,
        see kir_protein_catalogue_status()) must come through pd.to_numeric as NaN, so
        pct_novel_protein is NaN for that gene -- never silently 0% or divide-by-zero-as-100%."""
        rows = []

        def add(gene, species, level, s_obs, chao2):
            rows.append({"gene": gene, "ancestry": "ALL", "level": level, "species": species,
                         "n_people": 1000, "s_obs": s_obs, "good_turing_coverage": 0.99,
                         "chao2": chao2, "chao2_se": 1.0, "chao2_undetected_f0hat": chao2 - (
                             0 if s_obs == "NA" else s_obs),
                         "q1": "", "q2": "", "chao_new_by_2n": 0.0})

        add("KIR2DP1", "kir", "genomic", 30, 35.0)
        add("KIR2DP1", "kir", "any_novel", 5, 6.0)
        add("KIR2DP1", "kir", "protein", "NA", float("nan"))
        add("KIR2DP1", "kir", "protein_novel", "NA", float("nan"))
        cov = pd.DataFrame(rows)
        slope = pd.DataFrame([{"gene": "KIR2DP1", "ancestry": "ALL", "level": "genomic",
                               "species": "kir", "curve": "distinct", "n_star": 1000,
                               "slope_per_1000": 0.1, "mean_at_n_star": 0.0}])
        metrics = m46.build_gene_metrics(cov, slope)
        row = metrics[metrics.gene == "KIR2DP1"].iloc[0]
        self.assertTrue(math.isnan(row["pct_novel_protein"]))


class TestTwoProportionZTestReuse(unittest.TestCase):
    def test_reuses_45s_implementation_and_matches_worked_example(self):
        m45 = _load_module("45_kir_recurrence_figure.py", "kir_recurrence_figure_test_target_46")
        # Marc's brief worked example: HLA singleton share ~34%, KIR ~18% (approximate figures
        # given in the task) -- confirm the z-test direction and rough magnitude on the ACTUAL
        # pooled genomic-level numbers 44 produced (962/2823 vs 236/1460), not just synthetic
        # toy numbers, so this test would catch a sign/formula regression against real values.
        p1, p2, z, pval = m45.two_proportion_ztest(962, 2823, 236, 1460)
        self.assertAlmostEqual(p1, 0.3408, places=3)
        self.assertAlmostEqual(p2, 0.1616, places=3)
        self.assertGreater(z, 10)  # highly significant, HLA >> KIR singleton share
        self.assertLess(pval, 1e-30)


def _ancestry_expand(cov, ancestries=("AFR", "AMR", "EAS", "EUR", "SAS")):
    """Copies every `ancestry == "ALL"` row of a synthetic `coverage_chao2.tsv`-shaped frame to
    each of `ancestries` too (same values -- a smoke-test fixture, not meant to be realistic),
    so `build_ancestry_simpson_rows()`/panel c has real (non-NaN) numbers to plot when rendering
    the full figure end to end in a test."""
    all_rows = cov[cov["ancestry"] == "ALL"]
    extra = []
    for anc in ancestries:
        sub = all_rows.copy()
        sub["ancestry"] = anc
        extra.append(sub)
    return pd.concat([cov] + extra, ignore_index=True)


class TestFiguresRenderAndLayout(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.cov = _ancestry_expand(_synthetic_cov_with_cds())
        self.slope = _synthetic_slope()
        self.qc = pd.DataFrame([{"gene": "KIR2DL2", "protein_catalogue_covered": False,
                                 "reason": "no_catalogued_protein_entries"}])
        self.metrics = m46.build_gene_metrics(self.cov, self.slope, self.qc)

    def _synthetic_curve(self):
        """Minimal saturation_curves.tsv-shaped frame so pooled_recurrence_from_curves() (used by
        fig_recurrence_composition) has something to sum -- one gene per species, one level."""
        rows = []
        for gene, species in [("A", "hla"), ("KIR2DL1", "kir")]:
            for level in ("genomic", "any_novel", "protein_novel"):
                for n in (10, 20):
                    rows.append({"gene": gene, "ancestry": "ALL", "level": level,
                                 "species": species, "n": n, "mean_distinct": 10.0,
                                 "lo2_5": 8.0, "hi97_5": 12.0, "mean_eq1": 4.0, "mean_eq2": 2.0,
                                 "mean_gt2": 4.0, "mean_ge20": 1.0})
        return pd.DataFrame(rows)

    def test_fig_catalogue_completeness_renders(self):
        # `fig_catalogue_completeness` calls `vc.save_fig(..., strict=True)` internally, which
        # raises RuntimeError on any check_layout() error-severity violation -- a clean return
        # here already means "0 layout violations", not just "no Python exception".
        out_stem = os.path.join(self.tmp.name, "fig_catalogue_completeness")
        m46.fig_catalogue_completeness(self.metrics, self.cov, out_stem)
        self.assertTrue(os.path.exists(out_stem + ".png"))
        self.assertTrue(os.path.exists(out_stem + ".pdf"))

    def test_fig_recurrence_composition_renders(self):
        out_stem = os.path.join(self.tmp.name, "fig_recurrence_composition")
        m46.fig_recurrence_composition(self._synthetic_curve(), out_stem)
        self.assertTrue(os.path.exists(out_stem + ".png"))
        self.assertTrue(os.path.exists(out_stem + ".pdf"))


class TestBuildRowOrder(unittest.TestCase):
    """`build_row_order()` is the SINGLE source of truth for which gene sits in which row, shared
    by the plotting code and this test -- but that only guards against the plotting code
    disagreeing with `build_row_order()`. This class checks `build_row_order()` ITSELF against an
    order independently recomputed straight from the metrics table, so a bug inside
    `build_row_order()` (e.g. an ascending/descending flip, or NaN sorted first instead of last)
    cannot pass just because the plotting code faithfully reproduces whatever it returns."""

    def setUp(self):
        self.cov = _synthetic_cov_with_cds()
        self.slope = _synthetic_slope()
        self.qc = pd.DataFrame([{"gene": "KIR2DL2", "protein_catalogue_covered": False,
                                 "reason": "no_catalogued_protein_entries"}])
        self.metrics = m46.build_gene_metrics(self.cov, self.slope, self.qc)
        self.rows = m46.build_row_order(self.metrics)

    def test_two_block_headers_present_first_in_each_block(self):
        headers = [r for r in self.rows if r["kind"] == "header"]
        self.assertEqual(len(headers), 2)
        self.assertEqual(headers[0]["species"], "hla")
        self.assertEqual(headers[1]["species"], "kir")
        self.assertEqual(self.rows[0]["kind"], "header")

    def test_gene_rows_sorted_by_protein_completeness_descending_nan_last(self):
        for species in ("hla", "kir"):
            got = [r["gene_display"] for r in self.rows
                   if r["kind"] == "gene" and r["species"] == species]
            expected = self.metrics[self.metrics.species == species].sort_values(
                "completeness_protein", ascending=False, na_position="last"
            )["gene_display"].tolist()
            self.assertEqual(got, expected,
                              f"{species} row order does not match an independent "
                              "sort of the metrics table")

    def test_pseudogene_with_nan_protein_completeness_sorts_last_in_its_block(self):
        # KIR2DL2 is flagged catalogue-uncovered by the qc fixture -> completeness_protein is NaN
        # -> must be the LAST kir gene row, never dropped and never placed as if complete.
        kir_genes = [r["gene_display"] for r in self.rows
                     if r["kind"] == "gene" and r["species"] == "kir"]
        self.assertEqual(kir_genes[-1], "KIR2DL2")
        self.assertIn("KIR2DL2", kir_genes)


class TestRowLabelBindingMatchesData(unittest.TestCase):
    """The task's explicit ask: a binding test tying row label <-> plotted value <-> the metrics
    table (this script's "TSV") together, for the new Cleveland-dot-plot design -- there are no
    leader-lined point labels left to check (the whole point of the redesign), so what must be
    verified instead is that the SHARED y-axis tick label at a given row is the same gene as the
    marker plotted at that row's y position, and that marker's x is the metrics table's own value
    for that gene, not a stale or mismatched one."""

    def setUp(self):
        self.cov = _synthetic_cov_with_cds()
        self.slope = _synthetic_slope()
        self.qc = pd.DataFrame([{"gene": "KIR2DL2", "protein_catalogue_covered": False,
                                 "reason": "no_catalogued_protein_entries"}])
        self.metrics = m46.build_gene_metrics(self.cov, self.slope, self.qc)

    def test_tick_label_and_protein_marker_match_metrics_table_row_by_row(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rows = m46.build_row_order(self.metrics)
        n = len(rows)
        ys = [n - 1 - i for i in range(n)]

        fig, ax = plt.subplots(figsize=(4, 6), constrained_layout=True)
        m46._style_row_axis(ax, rows, ys)
        m46._plot_completeness_panel(ax, rows, ys)
        fig.canvas.draw()

        tick_text_by_y = {y: t.get_text() for t, y in zip(ax.get_yticklabels(), ys)}
        protein_lines = [ln for ln in ax.lines if ln.get_marker() == "o"]

        n_checked = 0
        for r, y in zip(rows, ys):
            if r["kind"] == "header":
                self.assertEqual(tick_text_by_y[y], r["label"])
                continue
            # row label <-> the row's own identity (a catalogue-uncovered gene's tick label is
            # its gene_display PLUS the "(no protein ref.)" suffix -- see `row_label()`).
            self.assertEqual(tick_text_by_y[y], m46.row_label(r))
            self.assertTrue(tick_text_by_y[y].startswith(r["gene_display"]))
            # row label <-> the metrics table ("TSV"), independent of build_row_order's own
            # bookkeeping -- re-fetch the row fresh from `self.metrics` by (gene, species).
            tsv_row = self.metrics[(self.metrics.gene == r["gene"]) &
                                    (self.metrics.species == r["species"])].iloc[0]
            self.assertEqual(r["gene_display"], tsv_row["gene_display"])
            # plotted value <-> the metrics table, at the row's OWN y (not just "somewhere").
            matches = [ln for ln in protein_lines if abs(ln.get_ydata()[0] - y) < 1e-9]
            if bool(tsv_row["protein_catalogue_covered"]):
                self.assertEqual(len(matches), 1,
                                  f"expected exactly one protein-level marker at row {y!r} "
                                  f"({r['gene_display']}), found {len(matches)}")
                self.assertAlmostEqual(matches[0].get_xdata()[0],
                                        tsv_row["completeness_protein"], places=6,
                                        msg=f"{r['gene_display']}: plotted protein completeness "
                                            "does not match the metrics table")
            else:
                # A catalogue-uncovered gene gets NO protein-level marker at all (2026-09-28 fix:
                # an earlier hollow-marker-at-x=0 design still read as "completeness 0") -- the
                # row's own tick-label suffix carries the caveat instead.
                self.assertEqual(len(matches), 0,
                                  f"{r['gene_display']} is catalogue-UNcovered and must get NO "
                                  f"protein-level marker, found {len(matches)}")
                self.assertIn(m46.NO_PROTEIN_REF_SUFFIX, tick_text_by_y[y])
            n_checked += 1
        self.assertEqual(n_checked, len(self.metrics))
        plt.close(fig)


class TestProteinCatalogueQCMerge(unittest.TestCase):
    """2026-09-27 coordinator ask: NA (catalogue-uncovered) genes must be rendered explicitly,
    never as 0 -- covers build_gene_metrics()'s merge of kir_protein_catalogue_qc.tsv and the
    completeness_protein column it enables."""

    def _qc(self):
        return pd.DataFrame([
            {"gene": "KIR2DL1", "protein_catalogue_covered": True, "reason": "ok"},
            {"gene": "KIR2DL2", "protein_catalogue_covered": False,
             "reason": "no_catalogued_protein_entries"},
        ])

    def test_completeness_protein_arithmetic(self):
        cov = _synthetic_cov()
        slope = _synthetic_slope()
        metrics = m46.build_gene_metrics(cov, slope, self._qc())
        row = metrics[(metrics.gene == "KIR2DL1") & (metrics.species == "kir")].iloc[0]
        # protein=50, chao2_protein=60 -> completeness_protein = 50/60.
        self.assertAlmostEqual(row["completeness_protein"], 50.0 / 60.0)

    def test_uncovered_gene_flagged_not_fabricated(self):
        cov = _synthetic_cov()
        slope = _synthetic_slope()
        metrics = m46.build_gene_metrics(cov, slope, self._qc())
        row = metrics[(metrics.gene == "KIR2DL2") & (metrics.species == "kir")].iloc[0]
        self.assertFalse(row["protein_catalogue_covered"])
        self.assertEqual(row["protein_catalogue_reason"], "no_catalogued_protein_entries")
        # Fixture KIR2DL2 has real (non-NaN) protein/protein_novel values -- the covered=False
        # flag from the QC table is independent of whether THIS fixture's numbers are NaN; the
        # actual v3 data's KIR2DP1/KIR3DP1 rows are simultaneously covered=False AND NaN. Confirm
        # covered=True genes are never accidentally flagged False by the merge.
        row2 = metrics[(metrics.gene == "KIR2DL1") & (metrics.species == "kir")].iloc[0]
        self.assertTrue(row2["protein_catalogue_covered"])

    def test_hla_defaults_covered_true_without_qc_table(self):
        # HLA genes never appear in kir_protein_catalogue_qc.tsv (KIR-only) -- must default to
        # covered=True/"ok", never inherit a stray False from the merge's NaN-fill.
        cov = _synthetic_cov()
        slope = _synthetic_slope()
        metrics = m46.build_gene_metrics(cov, slope, self._qc())
        for _, row in metrics[metrics.species == "hla"].iterrows():
            self.assertTrue(row["protein_catalogue_covered"])
            self.assertEqual(row["protein_catalogue_reason"], "ok")

    def test_no_qc_table_defaults_everyone_covered(self):
        cov = _synthetic_cov()
        slope = _synthetic_slope()
        metrics = m46.build_gene_metrics(cov, slope, None)
        self.assertTrue((metrics["protein_catalogue_covered"]).all())


class TestNaGenesRenderedExplicitly(unittest.TestCase):
    """Coordinator ask: a gene whose protein-level metrics are NaN (catalogue-uncovered, e.g.
    KIR2DP1/KIR3DP1) must be rendered explicitly -- never a fabricated 0 or a normal filled
    marker. 2026-09-28 fix: an earlier hollow-marker-at-x=0 design still read as "completeness 0"
    at a glance despite the hollow styling, so it now gets NO protein-level marker at all; the
    caveat lives in the row's own y-tick label instead (`NO_PROTEIN_REF_SUFFIX`)."""

    def _metrics_with_na_gene(self):
        cov = _synthetic_cov()
        rows = cov.to_dict("records")
        rows.append({"gene": "KIR2DP1", "ancestry": "ALL", "level": "genomic", "species": "kir",
                     "n_people": 1000, "s_obs": 30, "good_turing_coverage": 0.99, "chao2": 35.0,
                     "chao2_se": 1.0, "chao2_undetected_f0hat": 5.0, "q1": "", "q2": "",
                     "chao_new_by_2n": 0.0})
        rows.append({"gene": "KIR2DP1", "ancestry": "ALL", "level": "any_novel", "species": "kir",
                     "n_people": 1000, "s_obs": 5, "good_turing_coverage": 0.99, "chao2": 6.0,
                     "chao2_se": 1.0, "chao2_undetected_f0hat": 1.0, "q1": "", "q2": "",
                     "chao_new_by_2n": 0.0})
        rows.append({"gene": "KIR2DP1", "ancestry": "ALL", "level": "protein", "species": "kir",
                     "n_people": 1000, "s_obs": "NA", "good_turing_coverage": "NA",
                     "chao2": "NA", "chao2_se": "NA", "chao2_undetected_f0hat": "NA", "q1": "",
                     "q2": "", "chao_new_by_2n": "NA"})
        rows.append({"gene": "KIR2DP1", "ancestry": "ALL", "level": "protein_novel",
                     "species": "kir", "n_people": 1000, "s_obs": "NA",
                     "good_turing_coverage": "NA", "chao2": "NA", "chao2_se": "NA",
                     "chao2_undetected_f0hat": "NA", "q1": "", "q2": "", "chao_new_by_2n": "NA"})
        cov2 = pd.DataFrame(rows)
        slope = _synthetic_slope()
        qc = pd.DataFrame([{"gene": "KIR2DP1", "protein_catalogue_covered": False,
                            "reason": "no_catalogued_protein_entries"}])
        return m46.build_gene_metrics(cov2, slope, qc)

    def test_na_gene_gets_no_marker_and_a_labeled_row_not_a_fabricated_zero(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        metrics = self._metrics_with_na_gene()
        row = metrics[metrics.gene == "KIR2DP1"].iloc[0]
        self.assertTrue(math.isnan(row["pct_novel_protein"]))
        self.assertTrue(math.isnan(row["completeness_protein"]))

        rows = m46.build_row_order(metrics)
        n = len(rows)
        ys = [n - 1 - i for i in range(n)]
        fig, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)
        m46._style_row_axis(ax, rows, ys)
        m46._plot_completeness_panel(ax, rows, ys)

        kir2dp1_row = next(r for r in rows if r.get("gene") == "KIR2DP1")
        kir2dp1_y = ys[rows.index(kir2dp1_row)]
        protein_lines = [ln for ln in ax.lines if ln.get_marker() == "o"
                          and abs(ln.get_ydata()[0] - kir2dp1_y) < 1e-9]
        self.assertEqual(len(protein_lines), 0,
                          "a catalogue-uncovered gene must get NO protein-level marker at all "
                          "(never a fabricated 0, never even a hollow marker at 0)")

        tick_texts = [t.get_text() for t in ax.get_yticklabels()]
        self.assertTrue(any("KIR2DP1" in t and m46.NO_PROTEIN_REF_SUFFIX in t for t in tick_texts),
                        f"expected KIR2DP1's own row label to carry the caveat, got: {tick_texts}")
        plt.close(fig)



def _synthetic_cov_with_cds(base=None):
    """`_synthetic_cov()` (or `base`, if given) plus a `cds` row per gene (v4b headline
    granularity) -- hand-computable completeness_cds."""
    cov = base if base is not None else _synthetic_cov()
    rows = cov.to_dict("records")

    def add_cds(gene, species, s_obs, chao2):
        rows.append({"gene": gene, "ancestry": "ALL", "level": "cds", "species": species,
                     "n_people": 1000, "s_obs": s_obs, "good_turing_coverage": 0.98,
                     "chao2": chao2, "chao2_se": 1.0, "chao2_undetected_f0hat": chao2 - s_obs,
                     "q1": "", "q2": "", "chao_new_by_2n": 0.0})

    add_cds("A", "hla", 60, 120.0)      # completeness_cds = 0.5
    add_cds("B", "hla", 120, 150.0)     # completeness_cds = 0.8
    add_cds("KIR2DL1", "kir", 30, 40.0)  # completeness_cds = 0.75
    add_cds("KIR2DL2", "kir", 25, 50.0)  # completeness_cds = 0.5
    return pd.DataFrame(rows)


class TestCdsLevelMetrics(unittest.TestCase):
    """v4b (2026-09-27): CDS-level metrics added to `build_gene_metrics` -- Marc's ask to base the
    headline KIR-vs-HLA comparison on CDS and protein, not the genomic (span-level) upper bound."""

    def setUp(self):
        self.cov = _synthetic_cov_with_cds()
        self.slope = _synthetic_slope()
        self.metrics = m46.build_gene_metrics(self.cov, self.slope)

    def test_completeness_cds_arithmetic(self):
        row = self.metrics[(self.metrics.gene == "A") & (self.metrics.species == "hla")].iloc[0]
        self.assertAlmostEqual(row["completeness_cds"], 0.5)
        row = self.metrics[(self.metrics.gene == "B") & (self.metrics.species == "hla")].iloc[0]
        self.assertAlmostEqual(row["completeness_cds"], 0.8)
        row = self.metrics[(self.metrics.gene == "KIR2DL1") & (self.metrics.species == "kir")].iloc[0]
        self.assertAlmostEqual(row["completeness_cds"], 0.75)
        row = self.metrics[(self.metrics.gene == "KIR2DL2") & (self.metrics.species == "kir")].iloc[0]
        self.assertAlmostEqual(row["completeness_cds"], 0.5)

    def test_completeness_cds_nan_when_no_cds_row(self):
        # Plain _synthetic_cov() has no `cds` level rows at all -- every gene's completeness_cds
        # must come through as NaN, never a fabricated 0 or a silent KeyError.
        metrics = m46.build_gene_metrics(_synthetic_cov(), self.slope)
        self.assertTrue(metrics["completeness_cds"].isna().all())

    def test_existing_genomic_completeness_column_unchanged(self):
        # Adding the cds columns must not disturb the pre-existing genomic-level `completeness`.
        row = self.metrics[(self.metrics.gene == "A") & (self.metrics.species == "hla")].iloc[0]
        self.assertAlmostEqual(row["completeness"], 0.5)


class TestAncestrySimpsonPanel(unittest.TestCase):
    """Panel c (the optional per-ancestry strip): HLA-minus-KIR protein-level Chao2 completeness,
    aggregated straight from `coverage_chao2.tsv` (sum(s_obs)/sum(chao2) per species x ancestry,
    NOT an average of per-gene ratios -- see `_species_chao2_completeness()` docstring), for the
    5 well-powered ancestries plus a pooled row."""

    def test_rows_cover_five_ancestries_plus_pooled_in_order(self):
        cov = _ancestry_expand(_synthetic_cov())
        rows = m46.build_ancestry_simpson_rows(cov)
        self.assertEqual([r["label"] for r in rows],
                          ["AFR", "AMR", "EAS", "EUR", "SAS", "Pooled (ALL)"])
        self.assertFalse(rows[0]["pooled"])
        self.assertTrue(rows[-1]["pooled"])

    def test_diff_is_hla_minus_kir_aggregated_from_coverage_table(self):
        # Hand-computable: one HLA gene (genomic=100, protein s_obs=80, chao2=100 -> 0.8) and one
        # KIR gene (protein s_obs=40, chao2=100 -> 0.4) in a single ancestry -> diff = 0.4.
        rows = []

        def add(gene, species, level, s_obs, chao2, ancestry="AFR"):
            rows.append({"gene": gene, "ancestry": ancestry, "level": level, "species": species,
                         "n_people": 500, "s_obs": s_obs, "good_turing_coverage": 0.95,
                         "chao2": chao2, "chao2_se": 1.0, "chao2_undetected_f0hat": chao2 - s_obs,
                         "q1": "", "q2": "", "chao_new_by_2n": 0.0})

        add("A", "hla", "protein", 80, 100.0)
        add("KIR2DL1", "kir", "protein", 40, 100.0)
        cov = pd.DataFrame(rows)
        result = m46._species_chao2_completeness(cov, "AFR", "hla")
        self.assertAlmostEqual(result, 0.8)
        result_kir = m46._species_chao2_completeness(cov, "AFR", "kir")
        self.assertAlmostEqual(result_kir, 0.4)

    def test_missing_ancestry_level_yields_nan_not_a_fabricated_zero(self):
        cov = pd.DataFrame([{"gene": "A", "ancestry": "AFR", "level": "protein", "species": "hla",
                             "n_people": 500, "s_obs": 10, "good_turing_coverage": 0.9,
                             "chao2": 20.0, "chao2_se": 1.0, "chao2_undetected_f0hat": 10.0,
                             "q1": "", "q2": "", "chao_new_by_2n": 0.0}])
        # No KIR rows at all for AFR -> sum(chao2)==0 -> must be NaN, not a divide-by-zero 0/0->0.
        result = m46._species_chao2_completeness(cov, "AFR", "kir")
        self.assertTrue(math.isnan(result))


class TestFigCatalogueCompletenessFullFigure(unittest.TestCase):
    """End-to-end render of the Cleveland dot plot (all three panels together) on a small but
    ancestry-aware synthetic fixture -- confirms `check_layout(strict=True)` passes (via
    `vc.save_fig`'s own internal call) for the full figure, not just an isolated panel."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.cov = _ancestry_expand(_synthetic_cov_with_cds())
        self.slope = _synthetic_slope()
        self.qc = pd.DataFrame([{"gene": "KIR2DL2", "protein_catalogue_covered": False,
                                 "reason": "no_catalogued_protein_entries"}])
        self.metrics = m46.build_gene_metrics(self.cov, self.slope, self.qc)

    def test_renders_without_layout_violation(self):
        out_stem = os.path.join(self.tmp.name, "fig_catalogue_completeness")
        m46.fig_catalogue_completeness(self.metrics, self.cov, out_stem)
        self.assertTrue(os.path.exists(out_stem + ".png"))
        self.assertTrue(os.path.exists(out_stem + ".pdf"))


class TestMainWritesPerAncestryMetrics(unittest.TestCase):
    """Marc's ask: 'check it per gene and per ancestry' -- `main()` now also writes
    `46_catalogue_metrics_by_ancestry.tsv`. Tests the underlying loop logic directly (build a
    per-ancestry frame the same way `main()` does) rather than invoking `main()` and its argparse/
    file-path plumbing."""

    def test_per_ancestry_metrics_concat_has_ancestry_column_and_all_ancestries(self):
        rows = []

        def add(gene, species, level, ancestry, s_obs, chao2=None):
            rows.append({"gene": gene, "ancestry": ancestry, "level": level, "species": species,
                         "n_people": 500, "s_obs": s_obs, "good_turing_coverage": 0.95,
                         "chao2": chao2 if chao2 is not None else s_obs * 1.5, "chao2_se": 1.0,
                         "chao2_undetected_f0hat": (chao2 or s_obs * 1.5) - s_obs,
                         "q1": "", "q2": "", "chao_new_by_2n": 0.0})

        for anc in ("AFR", "EUR"):
            add("A", "hla", "genomic", anc, 100, chao2=200.0)
            add("A", "hla", "any_novel", anc, 20, chao2=40.0)
            add("A", "hla", "protein", anc, 90, chao2=95.0)
            add("A", "hla", "protein_novel", anc, 0, chao2=0.0)
        cov = pd.DataFrame(rows)
        slope = pd.DataFrame([{"gene": "A", "ancestry": a, "level": "genomic", "species": "hla",
                               "curve": "distinct", "n_star": 500, "slope_per_1000": 1.0,
                               "mean_at_n_star": 0.0} for a in ("AFR", "EUR")])
        per_anc = []
        for anc in ("AFR", "EUR"):
            m = m46.build_gene_metrics(cov, slope, ancestry=anc)
            m.insert(0, "ancestry", anc)
            per_anc.append(m)
        combined = pd.concat(per_anc, ignore_index=True)
        self.assertEqual(set(combined["ancestry"]), {"AFR", "EUR"})
        self.assertEqual(len(combined), 2)  # one HLA-A row per ancestry


if __name__ == "__main__":
    unittest.main()
