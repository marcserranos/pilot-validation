#!/usr/bin/env python3
"""Tests for 59_depth_outcome_join.py -- synthetic data only, no AoU data.

A synthetic cohort where people whose window depth is below a planted relative depth (0.55) get each
technical outcome at a high rate, and everyone else at a low background rate:
  KIR window : hap1 split over two contigs (kir_incomplete); hap2 loses KIR2DL4 (kir_fw_missing)
  MHC window : hap1 loses HLA-B (hla_absent) and pairs DQA1*01 with DQB1*02 on one contig
               (dq_cross_group_physical); hap2 is fragmented into one contig per gene
               (mhc_fragmented) and has DQB1 150 kb from DQA1 (dq_anomalous).
Every person's depth is also shifted by a platform offset, so the cohort-normalised scale must
remove it. Covers: exact per-person outcome recovery through 50/51's rules, the bin gradient on both
scales, the cutoff rule recovering the planted threshold, the OR direction, the zero-events guard,
the kinship filter, PUBLIC masking, no person ids in PUBLIC output, INTERNAL refusing repo paths,
the lock and the KIR resume cache.

Run: python3 -m pytest scripts/hla_popgen/tests/test_59_depth_outcome_join.py -q
"""
import gzip
import importlib.util
import os
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
spec = importlib.util.spec_from_file_location("dj59", os.path.join(SRC, "59_depth_outcome_join.py"))
dj = importlib.util.module_from_spec(spec)
sys.modules["dj59"] = dj
spec.loader.exec_module(dj)

import _disclosure as _disc  # noqa: E402

PLANTED = 0.55
PLATFORM_SHIFT = {"revio": 0.85, "sequel2e": 0.80}  # region depth under genome depth, by platform
KIR_GENES_HAP = ["KIR3DL3", "KIR2DL3", "KIR2DL4", "KIR3DP1", "KIR3DL1", "KIR3DL2"]
HLA_ORDER = ["F", "A", "E", "C", "B", "DRA", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
DQ_OK = [("HLA-DQA1*01:01:01", "HLA-DQB1*05:01:01"), ("HLA-DQA1*05:01:01", "HLA-DQB1*02:01:01")]
DQ_BAD = ("HLA-DQA1*01:02:01", "HLA-DQB1*02:01:01")
N_RELATED_PAIRS = 10


def rel_draw(rng, n):
    """Normalised depth: ~8% thin people over 0.2-0.7, the rest around 1.0 (continuous)."""
    thin = rng.random(n) < 0.08
    return np.clip(np.where(thin, rng.uniform(0.2, 0.7, n), rng.normal(1.0, 0.08, n)), 0.05, 1.6)


def write_kir_gtf(path, genes, split=False):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write("## synthetic\n")
        for i, g in enumerate(genes):
            s = 1000 + i * 20000
            ctg = "ctg2" if split and i >= 3 else "ctg1"
            fh.write(f'{ctg}\tImmuannot\tgene\t{s}\t{s+15000}\t.\t+\t.\tgene_id "{g}"; gene_name "{g}";\n')
            fh.write(f'{ctg}\tImmuannot\ttranscript\t{s}\t{s+15000}\t.\t+\t.\tgene_id "{g}"; '
                     f'gene_name "{g}"; consensus "{g}*001";\n')


class Cohort:
    def __init__(self, tmp, n=3000, seed=7, p_low=0.6, p_high=0.012):
        rng = np.random.default_rng(seed)
        self.tmp = tmp
        pids = [str(2_000_000 + i) for i in range(n)]
        self.pids = pids
        anc = np.array(["eur", "afr", "amr", "sas"])[rng.integers(0, 4, n)]
        plat = np.array(["revio", "sequel2e"])[(rng.random(n) < 0.3).astype(int)]
        self.norm_mhc, self.norm_kir = rel_draw(rng, n), rel_draw(rng, n)
        shift = np.array([PLATFORM_SHIFT[p] for p in plat])

        def hit(x):
            return rng.random(n) < np.where(x < PLANTED, p_low, p_high)
        self.kir_split, self.kir_drop = hit(self.norm_kir), hit(self.norm_kir)
        self.b_absent, self.dq_cross = hit(self.norm_mhc), hit(self.norm_mhc)
        self.frag, self.dq_far = hit(self.norm_mhc), hit(self.norm_mhc)

        rows = []
        for i, pid in enumerate(pids):
            for q in (0, 20):
                f = 1.0 if q == 0 else 0.97
                for region, rel in (("control_chr2", 1.0), ("window_MHC", self.norm_mhc[i] * shift[i]),
                                    ("window_KIR", self.norm_kir[i] * shift[i])):
                    r = round(rel * (f if region != "control_chr2" else 1.0), 4)
                    rows.append(dict(person_id=pid, region=region, mapq=q, mean_depth=r * 20,
                                     n_reads=100, control_depth=20.0, rel_depth=r))
        self.depth = os.path.join(tmp, "per_person_depth.tsv")
        with open(self.depth, "w") as fh:
            fh.write(_disc.INTERNAL_HEADER_COMMENT + "\n")
            pd.DataFrame(rows)[dj.m58().PER_PERSON_COLS].to_csv(fh, sep="\t", index=False)

        self.cohort = os.path.join(tmp, "cohort_membership.tsv")
        pd.DataFrame(dict(person_id=pids, in_lr=True, ancestry_pred=anc, platform=plat)).to_csv(
            self.cohort, sep="\t", index=False)

        # related pairs among the last 2*N people: greedy removal drops exactly one per pair
        rel = [dict(**{"i.s": pids[-2 * k - 1], "j.s": pids[-2 * k - 2]}, kin=0.25)
               for k in range(N_RELATED_PAIRS)]
        self.relatedness = os.path.join(tmp, "rel.tsv")
        pd.DataFrame(rel).to_csv(self.relatedness, sep="\t", index=False)

        t1 = []
        for i, pid in enumerate(pids):
            for hap in ("hap1", "hap2"):
                dqa, dqb = DQ_OK[rng.integers(0, 2)]
                if hap == "hap1" and self.dq_cross[i]:
                    dqa, dqb = DQ_BAD
                for k, g in enumerate(HLA_ORDER):
                    if hap == "hap1" and self.b_absent[i] and g == "B":
                        continue
                    cons = {"DQA1": dqa, "DQB1": dqb}.get(g, f"HLA-{g}*01:01:01")
                    ctg = f"{hap}_ctg1"
                    if hap == "hap2" and self.frag[i]:  # one contig per gene, DQA1/DQB1 kept together
                        ctg = f"{hap}_frag{HLA_ORDER.index('DQA1') if g == 'DQB1' else k}"
                    start = 1000 + 10000 * k
                    if g == "DQB1" and hap == "hap2" and self.dq_far[i]:
                        start += 150_000
                    t1.append(dict(person_id=pid, hap=hap, contig=ctg, gene=f"HLA-{g}", copy_index=1,
                                   gene_start=start, gene_end=start + 5000, consensus=cons,
                                   novelty_depth="", template_distance=0, cds_mut="",
                                   template_warning=""))
        self.table1 = os.path.join(tmp, "hla_calls_rich.tsv")
        pd.DataFrame(t1).to_csv(self.table1, sep="\t", index=False)

        # KIR GTFs (the first 5 people have no KIR directory: not assessable)
        self.kir_root = os.path.join(tmp, "pipeline_outputs_kir")
        os.makedirs(self.kir_root, exist_ok=True)
        for i, pid in enumerate(pids[5:], 5):
            base = os.path.join(self.kir_root, pid, "immuannot_output")
            write_kir_gtf(os.path.join(base, "hap1.gtf.gz"), KIR_GENES_HAP, split=self.kir_split[i])
            genes = [g for g in KIR_GENES_HAP if not (self.kir_drop[i] and g == "KIR2DL4")]
            write_kir_gtf(os.path.join(base, "hap2.gtf.gz"), genes)

    def run(self, out, internal, *extra):
        return dj.main(["--depth-table", self.depth, "--table1", self.table1, "--cohort", self.cohort,
                        "--relatedness-table", self.relatedness, "--kir-outroot", self.kir_root,
                        "--out-dir", out, "--internal-dir", internal, "--workers", "2", *extra])


class Test59(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.c = Cohort(cls.tmp)
        cls.out = os.path.join(cls.tmp, "out")
        cls.internal = os.path.join(cls.tmp, "internal")
        cls.c.run(cls.out, cls.internal)
        cls.pp = pd.read_csv(os.path.join(cls.internal, "per_person_depth_outcome.tsv"), sep="\t",
                             comment="#", dtype={"person_id": str}).set_index("person_id")

    def read(self, name):
        return pd.read_csv(os.path.join(self.out, name), sep="\t", dtype=str, keep_default_na=False)

    def test_internal_header_and_kinship(self):
        first = open(os.path.join(self.internal, "per_person_depth_outcome.tsv")).readline()
        self.assertEqual(first.strip(), _disc.INTERNAL_HEADER_COMMENT)
        self.assertEqual(len(self.pp), len(self.c.pids) - N_RELATED_PAIRS)

    def test_per_person_outcomes_match_plant(self):
        idx = [i for i, p in enumerate(self.c.pids) if p in self.pp.index]
        pp = self.pp.loc[[self.c.pids[i] for i in idx]]
        kir_ok = np.array([i >= 5 for i in idx])
        self.assertTrue(np.isnan(pp["kir_incomplete"].to_numpy()[~kir_ok]).all())
        for col, plant in (("kir_incomplete", self.c.kir_split), ("kir_fw_missing", self.c.kir_drop),
                           ("kir_fw_missing_complete", self.c.kir_drop)):
            np.testing.assert_array_equal(pp[col].to_numpy()[kir_ok], plant[idx][kir_ok].astype(float), col)
        for col, plant in (("hla_absent", self.c.b_absent), ("dq_cross_group_physical", self.c.dq_cross),
                           ("dq_cross_group_any", self.c.dq_cross), ("mhc_fragmented", self.c.frag),
                           ("dq_anomalous", self.c.dq_far)):
            np.testing.assert_array_equal(pp[col].to_numpy(), plant[idx].astype(float), col)
        self.assertEqual(pp["dq_split"].sum(), 0)
        self.assertEqual(pp["hla_hap_empty"].sum(), 0)
        self.assertEqual(int(pp["absent_B"].sum()), int(self.c.b_absent[idx].sum()))

    def test_normalised_scale_removes_platform_offset(self):
        idx = [self.c.pids.index(p) for p in self.pp.index]
        norm = self.pp["norm_window_MHC_q0"].to_numpy()
        self.assertLess(np.median(np.abs(norm - self.c.norm_mhc[idx])), 0.01)
        rs = dict(pd.read_csv(os.path.join(self.out, "run_summary.tsv"), sep="\t").values)
        self.assertAlmostEqual(float(rs["cohort_median_rel_window_MHC_q0_revio"]), 0.85, delta=0.02)
        self.assertAlmostEqual(float(rs["cohort_median_rel_window_MHC_q0_sequel2e"]), 0.80, delta=0.02)

    def test_bin_gradient(self):
        b = self.read("outcome_by_depth_bin.tsv")
        for scale, lo_bin, hi_bin in (("norm", "<0.5", "1-1.1"), ("abs", "<0.4", ">=1")):
            for outcome in ("kir_fw_missing", "dq_cross_group_physical", "hla_absent", "mhc_fragmented"):
                s = b[(b["outcome"] == outcome) & (b["mapq"] == "0") & (b["stratum"] == "ALL")
                      & (b["scale"] == scale)]
                rate = dict(zip(s["depth_bin"], s["pct_events"]))
                nev = dict(zip(s["depth_bin"], s["n_events"]))
                self.assertGreater(float(rate[lo_bin]), 30, (scale, outcome))
                if rate.get(hi_bin, "") == "":
                    self.assertIn(nev.get(hi_bin, "0"), ("<20", "0"), (scale, outcome))
                else:
                    self.assertLess(float(rate[hi_bin]), 3, (scale, outcome))
        self.assertIn("platform", set(b["stratum_type"]))
        cont = b[(b["outcome"] == "mhc_contigs_mean") & (b["scale"] == "norm") & (b["stratum"] == "ALL")
                 & (b["mapq"] == "0")]
        m = dict(zip(cont["depth_bin"], cont["mean_value"]))
        self.assertGreater(float(m["<0.5"]), float(m["1-1.1"]))

    def test_cutoff_rule_recovers_planted_threshold(self):
        p = self.read("proposed_cutoff.tsv")
        for outcome in ("kir_incomplete", "kir_fw_missing", "hla_absent", "dq_cross_group_physical",
                        "dq_anomalous", "mhc_fragmented"):
            for q in ("0", "20"):
                r = p[(p["outcome"] == outcome) & (p["mapq"] == q) & (p["scale"] == "norm")].iloc[0]
                self.assertEqual(r["status"], "proposed", (outcome, q))
                self.assertEqual(r["primary"], "True")
                # smallest c removing >= 50% of events; planted step at 0.55 (MAPQ20 is 0.97x, but the
                # normalisation divides that out)
                self.assertLessEqual(float(r["proposed_cutoff"]), PLANTED + 0.05 + 1e-9, (outcome, q))
                self.assertGreaterEqual(float(r["proposed_cutoff"]), 0.45 - 1e-9, (outcome, q))
                self.assertIn("smallest c", r["rule"])
        reg = p[(p["role"] == "region_proposal") & (p["scale"] == "norm")]
        self.assertEqual(set(reg["status"]), {"proposed"})
        neg = p[p["outcome"] == "kir_fw_missing_complete"]
        self.assertEqual(set(neg["role"]), {"negative_control"})
        # abs scale is also reported, cutoffs shifted by the ~0.84 platform offset
        ab = p[(p["outcome"] == "hla_absent") & (p["scale"] == "abs") & (p["mapq"] == "0")].iloc[0]
        self.assertEqual(ab["status"], "proposed")
        self.assertLessEqual(float(ab["proposed_cutoff"]), PLANTED * 0.85 + 0.05 + 1e-9)

    def test_scan_and_complementary_suppression(self):
        scan = self.read("cutoff_scan.tsv")
        s = scan[(scan["outcome"] == "kir_fw_missing") & (scan["mapq"] == "0") & (scan["ancestry"] == "ALL")
                 & (scan["scale"] == "norm")]
        removed = pd.to_numeric(s["pct_of_events_removed"], errors="coerce").dropna().to_numpy()
        self.assertGreater(len(removed), 2)
        self.assertTrue((np.diff(removed) >= 0).all())
        m = scan["n_events_kept"] == "<20"
        self.assertTrue(m.any())
        for col in ("n_events_removed", "pct_of_events_removed", "J", "outcome_rate_removed_pct"):
            self.assertTrue(scan.loc[m, col].isin(["", "<20", "0"]).all(), col)
        for anc in ("AFR", "AMR", "EUR", "SAS"):
            self.assertIn(anc, set(scan["ancestry"]))
        self.assertEqual(set(scan["scale"]), {"norm", "abs"})

    def test_logistic_or(self):
        lt = self.read("logistic_or.tsv")
        for outcome in ("kir_fw_missing", "dq_cross_group_physical", "hla_absent"):
            for scale in ("norm", "abs"):
                r = lt[(lt["outcome"] == outcome) & (lt["mapq"] == "0") & (lt["scale"] == scale)].iloc[0]
                self.assertEqual(r["status"], "ok", (outcome, scale))
                self.assertGreater(float(r["or_lo"]), 1.0)
                self.assertIn("platform[", r["covariates"])
                self.assertIn("ancestry[", r["covariates"])
        self.assertEqual(lt[lt["outcome"] == "dq_split"]["status"].unique().tolist(), ["too_few_events"])
        cl = self.read("contigs_linear.tsv")
        r = cl[(cl["scale"] == "norm") & (cl["mapq"] == "0")].iloc[0]
        self.assertGreater(float(r["lo"]), 0)

    def test_numpy_engine_matches_statsmodels(self):
        rng = np.random.default_rng(3)
        n = 3000
        x = rng.uniform(0.2, 1.4, n)
        z = rng.integers(0, 2, n).astype(float)
        y = (rng.random(n) < 1 / (1 + np.exp(-(1.0 - 4.0 * x + 0.5 * z)))).astype(float)
        X = np.column_stack([np.ones(n), x, z])
        b, se, conv = dj.logit_irls(X, y)
        self.assertTrue(conv)
        self.assertAlmostEqual(b[1], -4.0, delta=0.6)
        try:
            import statsmodels.api as sm
        except ImportError:
            return
        ref = sm.Logit(y, X).fit(disp=0)
        np.testing.assert_allclose(b, ref.params, rtol=1e-5, atol=1e-6)
        np.testing.assert_allclose(se, ref.bse, rtol=1e-4)

    def test_public_has_no_person_ids(self):
        for f in os.listdir(self.out):
            txt = open(os.path.join(self.out, f)).read()
            for pid in (self.c.pids[0], self.c.pids[1234], self.c.pids[-1]):
                self.assertNotIn(pid, txt, f)

    def test_rule_guards(self):
        x = np.linspace(0.2, 1.2, 1000)
        anc = np.array(["EUR", "AFR"] * 500)
        r = dj.propose_cutoff(np.zeros(1000), x, anc, dj.DEFAULT_NORM_CUTOFFS, 0.5, 0.1, 100)
        self.assertEqual(r["status"], "no_events")
        self.assertIsNone(r["cutoff"])
        rng = np.random.default_rng(1)
        x = rng.uniform(0.2, 1.4, 4000)
        y = (rng.random(4000) < 0.05).astype(float)
        anc = rng.choice(["EUR", "AFR"], 4000)
        r = dj.propose_cutoff(y, x, anc, dj.DEFAULT_NORM_CUTOFFS, 0.5, 0.1, 100)
        self.assertEqual(r["status"], "every_cutoff_removes_too_many_people")
        x = np.where(rng.random(4000) < 0.02, 0.3, 1.0)
        r = dj.propose_cutoff(y, x, anc, dj.DEFAULT_NORM_CUTOFFS, 0.5, 0.1, 100)
        self.assertEqual(r["status"], "depth_does_not_explain")
        self.assertIsNotNone(r["largest_admissible_cutoff"])
        # one ancestry thinner than the rest: the people limit binds on it
        x = np.where((anc == "AFR") & (rng.random(4000) < 0.3), 0.5, 1.0)
        y = (x < 0.6).astype(float)
        r = dj.propose_cutoff(y, x, anc, dj.DEFAULT_NORM_CUTOFFS, 0.5, 0.1, 100)
        self.assertEqual(r["binding_ancestry"], "AFR")
        self.assertIsNone(r["cutoff"])

    def test_masking_small_cohort_and_resume(self):
        out = os.path.join(self.tmp, "out_small")
        internal = os.path.join(self.tmp, "internal_small")
        self.c.run(out, internal, "--people-limit", "60")
        pp = pd.read_csv(os.path.join(internal, "per_person_depth_outcome.tsv"), sep="\t", comment="#")
        self.assertEqual(len(pp), 60)
        b = pd.read_csv(os.path.join(out, "outcome_by_depth_bin.tsv"), sep="\t", dtype=str,
                        keep_default_na=False)
        for c in ("n_people", "n_events"):
            ok = b[c].isin(["<20", "0", ""]) | b[c].str.isdigit()
            self.assertTrue(ok.all())
            nums = pd.to_numeric(b[c], errors="coerce")
            self.assertFalse(((nums > 0) & (nums < 20)).any(), c)
        self.assertTrue((b.loc[b["n_events"] == "<20", "pct_events"] == "").all())
        self.assertTrue((b.loc[b["n_people"] == "<20", "pct_events"] == "").all())
        self.assertTrue((b.loc[b["n_people"] == "<20", "mean_value"] == "").all())
        scan = pd.read_csv(os.path.join(out, "cutoff_scan.tsv"), sep="\t", dtype=str, keep_default_na=False)
        for col in scan.columns:
            if col.startswith("n_"):
                nums = pd.to_numeric(scan[col], errors="coerce")
                self.assertFalse(((nums > 0) & (nums < 20)).any(), col)
        lt = pd.read_csv(os.path.join(out, "logistic_or.tsv"), sep="\t", dtype=str, keep_default_na=False)
        self.assertFalse((lt["status"] == "ok").any())
        # resume: the KIR cache is reused, no duplicate rows
        cache = os.path.join(internal, "kir_outcomes.tsv")
        n1 = len(pd.read_csv(cache, sep="\t", comment="#"))
        self.c.run(out, internal, "--people-limit", "60")
        self.assertEqual(n1, len(pd.read_csv(cache, sep="\t", comment="#")))
        self.assertFalse(os.path.exists(os.path.join(internal, "59.lock")))

    def test_exclude_people_file(self):
        path = os.path.join(self.internal, "exclude_people.tsv")
        self.assertEqual(open(path).readline().strip(), _disc.INTERNAL_HEADER_COMMENT)
        ex = pd.read_csv(path, sep="\t", comment="#", dtype={"person_id": str})
        self.assertEqual(list(ex.columns), dj.EXCLUDE_COLS)
        # every person with depth, relatives included, one row per region
        self.assertEqual(len(ex), len(self.c.pids) * 2)
        self.assertEqual(set(ex["region"]), set(dj.REGIONS))
        self.assertEqual(set(ex["scale"]), {"norm"})
        self.assertEqual(set(ex["mapq"]), {0})
        p = self.read("proposed_cutoff.tsv")
        for region in dj.REGIONS:
            r = p[(p["role"] == "region_proposal") & (p["scale"] == "norm") & (p["mapq"] == "0")
                  & (p["region"] == region)].iloc[0]
            sub = ex[ex["region"] == region]
            self.assertEqual(set(sub["cutoff"]), {float(r["proposed_cutoff"])})
            np.testing.assert_array_equal(sub["below_cutoff"].to_numpy(),
                                          (sub["value"] < sub["cutoff"]).to_numpy())
        w = ex.pivot(index="person_id", columns="region", values="below_cutoff")
        e = ex.drop_duplicates("person_id").set_index("person_id")["exclude"]
        np.testing.assert_array_equal(e.loc[w.index].to_numpy(), w.any(axis=1).to_numpy())
        # planted thin people go, well-covered people stay
        idx = {p: i for i, p in enumerate(self.c.pids)}
        thin = [p for p in e.index if self.c.norm_mhc[idx[p]] < 0.4]
        fine = [p for p in e.index if min(self.c.norm_mhc[idx[p]], self.c.norm_kir[idx[p]]) > 0.8]
        self.assertTrue(e.loc[thin].all())
        self.assertFalse(e.loc[fine].any())
        # value is the normalised depth
        mhc = ex[ex["region"] == "window_MHC"].set_index("person_id")["value"]
        self.assertLess(np.median(np.abs(mhc.to_numpy() - self.c.norm_mhc[[idx[p] for p in mhc.index]])), 0.01)

    def test_exclude_region_option(self):
        out = os.path.join(self.tmp, "out_kironly")
        internal = os.path.join(self.tmp, "internal_kironly")
        self.c.run(out, internal, "--exclude-region", "window_KIR")
        ex = pd.read_csv(os.path.join(internal, "exclude_people.tsv"), sep="\t", comment="#",
                         dtype={"person_id": str})
        kir = ex[ex["region"] == "window_KIR"].set_index("person_id")
        mhc = ex[ex["region"] == "window_MHC"].set_index("person_id")
        np.testing.assert_array_equal(kir["exclude"].to_numpy(), kir["below_cutoff"].to_numpy())
        self.assertTrue((mhc["below_cutoff"] & ~mhc["exclude"]).any())  # MHC-only thin people kept
        with self.assertRaises(SystemExit):
            self.c.run(out, internal, "--exclude-region", "window_XYZ")

    def test_lock_refuses_second_run(self):
        internal = os.path.join(self.tmp, "internal_locked")
        os.makedirs(internal)
        with open(os.path.join(internal, "59.lock"), "w") as fh:
            fh.write(str(os.getpid()))  # a live pid
        with self.assertRaises(SystemExit):
            self.c.run(os.path.join(self.tmp, "out_locked"), internal, "--people-limit", "20")

    def test_internal_refuses_repo(self):
        with self.assertRaises(ValueError):
            self.c.run(os.path.join(self.tmp, "out_x"), os.path.join(SRC, "x_internal_59"))
        self.assertFalse(os.path.exists(os.path.join(SRC, "x_internal_59")))


if __name__ == "__main__":
    unittest.main()
