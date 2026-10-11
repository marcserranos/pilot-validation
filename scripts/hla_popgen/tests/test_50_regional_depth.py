#!/usr/bin/env python3
"""Tests for 50_regional_depth.py -- synthetic data only, no samtools, no AoU data.

Covers: region building from a GTF, bedcov parsing, relative depth against the control median,
TRB copy-number calls on a planted deletion, assessability of a paralog-confounded gene, the HWE
exact test, PUBLIC masking, failure accounting, and resume after a partial run.

Run: python3 -m pytest scripts/hla_popgen/tests/test_50_regional_depth.py -q
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
spec = importlib.util.spec_from_file_location("rd50", os.path.join(SRC, "50_regional_depth.py"))
rd = importlib.util.module_from_spec(spec)
sys.modules["rd50"] = rd
spec.loader.exec_module(rd)


def write_gtf(path):
    genes = [
        ("chr6", 29_941_000, 29_945_000, "HLA-A"),
        ("chr6", 32_637_000, 32_655_000, "HLA-DQB1"),
        ("chr7", 142_300_000, 142_300_600, "TRBV4-3"),
        ("chr7", 142_310_000, 142_310_600, "TRBV6-2"),
        ("chr7", 142_320_000, 142_320_600, "TRBV7-2"),
        ("chr7", 142_790_000, 142_790_050, "TRBJ1-1"),
        ("chr19", 54_724_000, 54_737_000, "KIR3DL3"),
        ("chr19", 54_850_000, 54_866_000, "KIR3DL2"),
        ("chr1", 1_000, 2_000, "OTHER1"),
        ("chr6", 100_000, 200_000, "NOTHLA"),
    ]
    with gzip.open(path, "wt") as fh:
        fh.write("##header\n")
        for c, s, e, n in genes:
            fh.write(f'{c}\tHAVANA\tgene\t{s}\t{e}\t.\t+\t.\tgene_id "X"; gene_name "{n}";\n')
            fh.write(f'{c}\tHAVANA\ttranscript\t{s}\t{e}\t.\t+\t.\tgene_name "{n}";\n')


class FakeBedcov:
    """Synthesises bedcov output: control depth = base; gene depth = base*rel (diploid=1).
    TRBV6-2 is paralog-confounded: MAPQ20 depth collapses to 10% of MAPQ0."""

    def __init__(self, regions, people_spec):
        self.regions = regions
        self.spec = people_spec  # bam -> dict(base=..., rel={name: rel})

    def __call__(self, bam, bed_path, mapq, timeout, samtools):
        s = self.spec[bam]
        lines = []
        for r in self.regions.itertuples():
            rel = s["rel"].get(r.name, 1.0)
            if r.name == "TRBV6-2" and mapq == 20:
                rel *= 0.1
            L = r.end - r.start
            sumd = s["base"] * rel * L
            lines.append(f"{r.chrom}\t{r.start}\t{r.end}\t{r.name}\t{r.region_class}\t{sumd}\t{int(sumd/15000)+1}")
        return rd.parse_bedcov("\n".join(lines))


class Test50(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.gtf = os.path.join(self.tmp, "g.gtf.gz")
        write_gtf(self.gtf)
        self.regions = rd.build_regions_from_gtf(self.gtf)

    def test_regions(self):
        names = set(self.regions["name"])
        self.assertIn("HLA-A", names)
        self.assertIn("TRBV4-3", names)
        self.assertIn("KIR3DL3", names)
        self.assertNotIn("OTHER1", names)
        self.assertNotIn("NOTHLA", names)
        self.assertEqual((self.regions["region_class"] == "control").sum(), 5)
        trbw = self.regions[self.regions["name"] == "window_TRB"].iloc[0]
        self.assertEqual(trbw.start, 142_300_000 - 1 - rd.LOCUS_PAD)
        self.assertEqual(trbw.end, 142_790_050 + rd.LOCUS_PAD)
        # BED is 0-based: GTF start 29,941,000 -> 29,940,999
        a = self.regions[self.regions["name"] == "HLA-A"].iloc[0]
        self.assertEqual(a.start, 29_940_999)

    def test_hwe(self):
        self.assertGreater(rd.hwe_exact_p(25, 50, 25), 0.9)
        self.assertLess(rd.hwe_exact_p(50, 0, 50), 1e-10)
        self.assertTrue(np.isnan(rd.hwe_exact_p(0, 0, 0)))

    def _cohort(self, n=600, seed=1):
        rng = np.random.default_rng(seed)
        anc = np.array(["SAS", "EUR"])[rng.integers(0, 2, n)]
        q = np.where(anc == "SAS", 0.6, 0.3)  # deletion allele freq
        cn = 2 - rng.binomial(2, q)
        spec, rows, man = {}, [], []
        for i in range(n):
            pid = str(1_000_000 + i)
            bam = os.path.join(self.tmp, f"{pid}.bam")
            open(bam, "w").close()
            base = rng.uniform(8, 30)
            rel = {"TRBV4-3": cn[i] / 2 + rng.normal(0, 0.04), "window_MHC": rng.uniform(0.4, 1.1)}
            spec[bam] = dict(base=base, rel=rel)
            rows.append(dict(person_id=pid, in_lr=True, ancestry_pred=anc[i].lower(), platform="revio"))
            man.append(dict(research_id=pid, platform="revio",
                            grch38_bam=rd.BUCKET_PREFIX + f"x/{pid}.bam"))
        cohort = os.path.join(self.tmp, "cohort.tsv")
        pd.DataFrame(rows).to_csv(cohort, sep="\t", index=False)
        root = os.path.join(self.tmp, "root")
        os.makedirs(os.path.join(root, "v9/wgs/long_read"), exist_ok=True)
        os.makedirs(os.path.join(root, "x"), exist_ok=True)
        for m in man:
            open(os.path.join(root, "x", f"{m['research_id']}.bam"), "w").close()
        pd.DataFrame(man).to_csv(os.path.join(root, rd.LR_MANIFEST_REL), sep="\t", index=False)
        spec = {os.path.join(root, "x", os.path.basename(k)): v for k, v in spec.items()}
        return cohort, root, spec, cn, anc

    def _run(self, cohort, root, spec, extra=()):
        bed = os.path.join(self.tmp, "regions.bed")
        rd.write_regions_bed(self.regions, bed)
        rd.run_bedcov = FakeBedcov(self.regions, spec)
        out = os.path.join(self.tmp, "out")
        internal = os.path.join(self.tmp, "internal")
        rd.main(["--regions-bed", bed, "--cohort", cohort, "--aou-root", root, "--out-dir", out,
                 "--internal-dir", internal, "--workers", "4", *extra])
        return out, internal

    def test_end_to_end(self):
        cohort, root, spec, cn, anc = self._cohort()
        out, internal = self._run(cohort, root, spec)
        assess = pd.read_csv(os.path.join(out, "trb_gene_assessability.tsv"), sep="\t")
        a = dict(zip(assess["gene"], assess["assessable"]))
        self.assertTrue(a["TRBV4-3"])
        self.assertFalse(a["TRBV6-2"])  # paralog-confounded
        pp = pd.read_csv(os.path.join(internal, "trb_cn_per_person.tsv"), sep="\t", comment="#",
                         dtype={"person_id": str})
        v = pp[pp["gene"] == "TRBV4-3"].sort_values("person_id")
        self.assertGreater((v["copies"].to_numpy() == cn).mean(), 0.97)
        tab = pd.read_csv(os.path.join(out, "trb_cn_by_ancestry.tsv"), sep="\t")
        sas = tab[(tab["gene"] == "TRBV4-3") & (tab["ancestry"] == "SAS")].iloc[0]
        eur = tab[(tab["gene"] == "TRBV4-3") & (tab["ancestry"] == "EUR")].iloc[0]
        self.assertAlmostEqual(float(sas["deletion_allele_freq"]), 0.6, delta=0.07)
        self.assertAlmostEqual(float(eur["deletion_allele_freq"]), 0.3, delta=0.07)
        self.assertGreater(float(sas["hwe_p"]), 1e-3)
        summ = pd.read_csv(os.path.join(out, "depth_summary_by_region.tsv"), sep="\t")
        mhc = summ[(summ["region"] == "window_MHC") & (summ["mapq"] == 0) & (summ["stratum"] == "ALL")].iloc[0]
        self.assertAlmostEqual(mhc["rel_p50"], 0.75, delta=0.05)
        # PUBLIC outputs never carry person ids
        for f in os.listdir(out):
            txt = open(os.path.join(out, f)).read()
            self.assertNotIn("1000001", txt, f)

    def test_masking_and_failures_and_resume(self):
        cohort, root, spec, cn, anc = self._cohort(n=30)
        # one person's BAM disappears -> counted as failure
        victim = sorted(spec)[0]
        os.remove(victim)
        out, internal = self._run(cohort, root, spec)
        run = dict(pd.read_csv(os.path.join(out, "run_summary.tsv"), sep="\t").values)
        self.assertEqual(int(run["people_measured"]), 29)
        self.assertEqual(int(run["people_failed"]), 1)
        tab = pd.read_csv(os.path.join(out, "trb_cn_by_ancestry.tsv"), sep="\t", dtype=str)
        self.assertTrue((tab["n_people"].isin(["<20"]) | tab["n_people"].str.isdigit()).all())
        self.assertTrue(tab.loc[tab["n_people"] == "<20", "deletion_allele_freq"].isna().all())
        # resume: a second run measures nobody new
        n_before = len(pd.read_csv(os.path.join(internal, "per_person_depth.tsv"), sep="\t", comment="#"))
        self._run(cohort, root, spec)
        n_after = len(pd.read_csv(os.path.join(internal, "per_person_depth.tsv"), sep="\t", comment="#"))
        self.assertEqual(n_before, n_after)

    def test_internal_refuses_repo(self):
        with self.assertRaises(ValueError):
            rd.append_rows(os.path.join(SRC, "x_internal.tsv"), [], rd.PER_PERSON_COLS)


if __name__ == "__main__":
    unittest.main()
