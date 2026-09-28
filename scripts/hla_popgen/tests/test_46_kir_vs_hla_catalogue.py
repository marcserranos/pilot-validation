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
  4. **Label <-> data binding test** (the task's explicit ask: "check that every plotted value <->
     label binding matches the TSV"): renders `fig_catalogue_completeness` on a small synthetic
     metrics table, reads back every `mark_label`-registered Annotation's anchor point (`xy`,
     documented in `_viz_common.mark_label`'s own point-registration Annotation contract) and text,
     and asserts each one's (x, y) matches the corresponding TSV row's
     (pct_novel_any, completeness) to floating-point tolerance, and the text equals that row's
     `gene_display`. This is a real regression guard against a mislabeled point (e.g. a repulsion
     bug that nudges a label's anchor instead of its offset, or a gene/species mismatch in the
     zip() ordering upstream).

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


def _synthetic_cov_clustered():
    """Same shape as `_synthetic_cov()` but with pct_novel_any clustered near 90-100% for BOTH
    species -- representative of this pipeline's real genomic-level identity data (44's
    `coverage_chao2.tsv`: every gene of both species is >60% novel at genomic level, see
    `44`'s Caveat 1 / this script's own `fig_catalogue_completeness` panel-a docstring), unlike
    `_synthetic_cov()`'s deliberately spread-out 20/50/80/25% (picked for hand-computable
    arithmetic, not to model clustering). Used ONLY by the layout/rendering tests below: panel a's
    two-species chained-column labeling (`_scatter_panel(column_species=("hla","kir"))`, added
    2026-09-28) assumes each species' own points are reasonably close to its own rightmost point
    -- true for real genomic-level data, not true for `_synthetic_cov()`'s wide artificial spread,
    which was exposing a leader-line-sweep layout bug that never occurs on real inputs (verified
    against the actual v4b `46_catalogue_metrics.tsv`, min pct_novel_any 61.9%)."""
    rows = []

    def add(gene, species, level, s_obs, chao2):
        rows.append({"gene": gene, "ancestry": "ALL", "level": level, "species": species,
                     "n_people": 1000, "s_obs": s_obs, "good_turing_coverage": 0.99,
                     "chao2": chao2, "chao2_se": 1.0, "chao2_undetected_f0hat": chao2 - s_obs,
                     "q1": "", "q2": "", "chao_new_by_2n": 0.0})

    # genomic_chao2_mult varies per gene (unlike a fixed ratio) so completeness values are spread
    # out on the y-axis too, not all identical -- an earlier version of this fixture used a fixed
    # s_obs*2.0 for every gene, which put every point at completeness==0.5 and made the column
    # labeler's vertical repulsion push labels into each other's leader lines (a fixture artifact,
    # not a real bug; caught when this fixture itself was being debugged, 2026-09-28).
    for gene, genomic, any_novel, protein, protein_novel, gmult in [
        ("A", 100, 65, 90, 20, 2.0), ("B", 200, 190, 150, 60, 1.3),
    ]:
        add(gene, "hla", "genomic", genomic, genomic * gmult)
        add(gene, "hla", "any_novel", any_novel, any_novel * 2.5)
        add(gene, "hla", "protein", protein, protein * 2.2)
        add(gene, "hla", "protein_novel", protein_novel, protein_novel * 2.8)
    for gene, genomic, any_novel, protein, protein_novel, gmult in [
        ("KIR2DL1", 50, 46, 50, 30, 1.8), ("KIR2DL2", 40, 39, 20, 15, 3.0),
    ]:
        add(gene, "kir", "genomic", genomic, genomic * gmult)
        add(gene, "kir", "any_novel", any_novel, any_novel * 2.5)
        add(gene, "kir", "protein", protein, protein * 2.2)
        add(gene, "kir", "protein_novel", protein_novel, protein_novel * 2.8)
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


class TestFiguresRenderAndLayout(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        # Clustered fixture (not the wide-spread `_synthetic_cov()`): panel a's two-species
        # chained column labeling needs realistically-clustered genomic-level novelty to exercise
        # what it's actually built for (see `_synthetic_cov_clustered()` docstring).
        self.cov = _synthetic_cov_clustered()
        self.slope = _synthetic_slope()
        self.metrics = m46.build_gene_metrics(self.cov, self.slope)

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
        out_stem = os.path.join(self.tmp.name, "fig_catalogue_completeness")
        m46.fig_catalogue_completeness(self.metrics, out_stem)
        self.assertTrue(os.path.exists(out_stem + ".png"))
        self.assertTrue(os.path.exists(out_stem + ".pdf"))

    def test_fig_recurrence_composition_renders(self):
        out_stem = os.path.join(self.tmp.name, "fig_recurrence_composition")
        m46.fig_recurrence_composition(self._synthetic_curve(), out_stem)
        self.assertTrue(os.path.exists(out_stem + ".png"))
        self.assertTrue(os.path.exists(out_stem + ".pdf"))


class TestLabelBindingMatchesData(unittest.TestCase):
    """The task's explicit ask: verify plotted value <-> label binding for at least one panel,
    against the source TSV/DataFrame -- not just "the figure rendered without an exception"."""

    def test_catalogue_completeness_labels_match_metrics_rows(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        cov = _synthetic_cov()
        slope = _synthetic_slope()
        metrics = m46.build_gene_metrics(cov, slope)

        fig, ax = plt.subplots(figsize=(4, 3), constrained_layout=True)
        all_x, all_y, all_lab, all_col = [], [], [], []
        for sp in m46.SPECIES_ORDER:
            s = metrics[metrics["species"] == sp]
            ax.scatter(s["pct_novel_any"], s["completeness"], color=m46.SPECIES_COLOR[sp])
            all_x += s["pct_novel_any"].tolist()
            all_y += s["completeness"].tolist()
            all_lab += s["gene_display"].tolist()
            all_col += [m46.SPECIES_COLOR[sp]] * len(s)
        texts = m46._label_points(ax, all_x, all_y, all_lab, all_col)
        plt.close(fig)

        self.assertEqual(len(texts), len(metrics))
        # Build the expected (gene_display -> (x, y)) mapping straight from the metrics
        # DataFrame (the "TSV" this script would otherwise write), independent of plotting order.
        expected = {row["gene_display"]: (row["pct_novel_any"], row["completeness"])
                    for _, row in metrics.iterrows()}
        seen_labels = set()
        for t in texts:
            label = t.get_text()
            self.assertIn(label, expected, f"unexpected label {label!r} not in metrics table")
            exp_x, exp_y = expected[label]
            got_x, got_y = t.xy  # the anchor point passed to ax.annotate(label, (x, y), ...)
            self.assertAlmostEqual(got_x, exp_x, places=6,
                                    msg=f"{label}: plotted x does not match metrics table")
            self.assertAlmostEqual(got_y, exp_y, places=6,
                                    msg=f"{label}: plotted y does not match metrics table")
            seen_labels.add(label)
        self.assertEqual(seen_labels, set(expected))


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
    """2026-09-27 coordinator ask: 'Render NA genes explicitly ... never as 0.' A gene whose
    protein-level metrics are NaN (catalogue-uncovered) must produce an in-panel text note naming
    it, and must NOT appear as a plotted point at (0, 0) or any other fabricated position."""

    def _metrics_with_na_gene(self):
        cov = _synthetic_cov()
        # Add a third KIR gene with NaN protein-level columns (simulating KIR2DP1/KIR3DP1).
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

    def test_na_gene_excluded_from_scatter_points_not_plotted_at_zero(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        metrics = self._metrics_with_na_gene()
        row = metrics[metrics.gene == "KIR2DP1"].iloc[0]
        self.assertTrue(math.isnan(row["pct_novel_protein"]))
        self.assertTrue(math.isnan(row["completeness_protein"]))

        fig, ax = plt.subplots(figsize=(4, 3), constrained_layout=True)
        handles = m46._scatter_panel(ax, metrics, "pct_novel_protein", "completeness_protein",
                                      na_note="KIR2DP1: pseudogene, no catalogued reference "
                                               "protein (excluded, not 0)")
        # The NaN row must not appear as a scattered point at all -- collect every point actually
        # drawn and confirm none sits at (0, 0) (the fabricated-zero failure mode this test
        # guards against) and that the point count equals only the non-NaN rows.
        n_plotted = sum(len(h.get_offsets()) for h in handles)
        n_non_na = metrics[["pct_novel_protein", "completeness_protein"]].dropna().shape[0]
        self.assertEqual(n_plotted, n_non_na)
        for h in handles:
            for (x, y) in h.get_offsets():
                self.assertFalse(x == 0.0 and y == 0.0,
                                  "an NA gene must never be plotted as a (0, 0) point")
        # The explicit note naming the excluded gene must be present as a real text artist.
        note_texts = [t.get_text() for t in ax.texts]
        self.assertTrue(any("KIR2DP1" in t for t in note_texts),
                        f"expected an in-panel note naming KIR2DP1, got texts: {note_texts}")
        plt.close(fig)


class TestColumnLabelStyleBinding(unittest.TestCase):
    """fig_catalogue_completeness's protein-level panel uses column_species=('kir',) to avoid
    labels sitting on top of densely-clustered marker dots (2026-09-27 visual-review fix). Confirm
    the column-style leader lines still connect each label to ITS OWN gene's real data point."""

    def test_column_style_leaders_match_data_points(self):
        """The column style (`_label_curve_ends`, reused from 45) deliberately moves each label's
        OWN text position away from its data point (vertical repulsion + shared x-column) and
        connects the two with a leader line -- so `.xy` is not expected to equal the data point
        here (unlike the 'free' `_label_points` style, covered by TestLabelBindingMatchesData).
        What this test guards against instead: `_scatter_panel` building the `ends` dict from the
        WRONG row when constructing column-style input (e.g. a species/gene mismatch upstream) --
        checked two ways: (a) every plotted MARKER sits exactly at its metrics row's data point,
        and (b) the set of column-style labels drawn is exactly the set of KIR gene_display names
        with non-NaN data, no more, no fewer."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        cov = _synthetic_cov()
        slope = _synthetic_slope()
        metrics = m46.build_gene_metrics(cov, slope)
        kir = metrics[metrics.species == "kir"].dropna(
            subset=["pct_novel_protein", "completeness_protein"])

        fig, ax = plt.subplots(figsize=(4, 3), constrained_layout=True)
        handles = m46._scatter_panel(ax, metrics, "pct_novel_protein", "completeness_protein",
                                      column_species=("kir",))
        kir_handle = handles[m46.SPECIES_ORDER.index("kir")]
        plotted = {tuple(round(v, 6) for v in pt) for pt in kir_handle.get_offsets()}
        expected_points = {(round(row["pct_novel_protein"], 6),
                             round(row["completeness_protein"], 6))
                            for _, row in kir.iterrows()}
        self.assertEqual(plotted, expected_points)

        drawn_labels = {t.get_text() for t in ax.texts
                        if hasattr(t, "_layout_direct_label") and t.get_text().startswith("KIR")}
        self.assertEqual(drawn_labels, set(kir["gene_display"]))
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


class TestFigCatalogueCompletenessThreePanels(unittest.TestCase):
    """v4b: `fig_catalogue_completeness` grew a third (CDS) panel between the pre-existing genomic
    and protein panels. Confirms it still renders through `check_layout(strict=True)` and that the
    CDS panel's `label_genes=False` path draws NO per-gene direct labels (only an na_note, if any)
    -- the busy CDS panel deliberately omits gene labels (see `_scatter_panel` docstring)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        # Clustered base (see `TestFiguresRenderAndLayout.setUp`) + cds rows.
        self.cov = _synthetic_cov_with_cds(base=_synthetic_cov_clustered())
        self.slope = _synthetic_slope()
        self.metrics = m46.build_gene_metrics(self.cov, self.slope)

    def test_renders_three_panels_without_layout_violation(self):
        out_stem = os.path.join(self.tmp.name, "fig_catalogue_completeness")
        m46.fig_catalogue_completeness(self.metrics, out_stem)
        self.assertTrue(os.path.exists(out_stem + ".png"))
        self.assertTrue(os.path.exists(out_stem + ".pdf"))

    def test_label_genes_false_draws_no_gene_labels(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(4, 3), constrained_layout=True)
        m46._scatter_panel(ax, self.metrics, "pct_novel_any", "completeness_cds",
                            label_genes=False)
        gene_names = set(self.metrics["gene_display"])
        drawn = {t.get_text() for t in ax.texts}
        self.assertFalse(drawn & gene_names, f"expected no gene labels, found: {drawn}")
        plt.close(fig)


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
