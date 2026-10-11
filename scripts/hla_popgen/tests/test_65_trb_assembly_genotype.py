#!/usr/bin/env python3
"""Tests for 65_trb_assembly_genotype.py -- synthetic assemblies only, no AoU data, no network.

Synthetic world: random-sequence TRB probes; TRBV6-3 is a 1-substitution copy of TRBV6-2 (must merge
into one class), TRBV4-2 is a 95%-identity paralog of TRBV4-3, and every haplotype also carries a
second, probe-less 95% near-copy of TRBV4-3 (the decoy). Neither may ever be called TRBV4-3.
Insertion haplotypes carry TRBV4-3/3-2/6-3/7-8 together, at planted per-ancestry frequencies.
One haplotype is fragmented (anchors on two contigs), one lacks the TRB locus, one person has no PAF
(anchor contig selection), one person's assembly is missing (failure).

Run: python3 -m pytest scripts/hla_popgen/tests/test_65_trb_assembly_genotype.py -q
"""
import gzip
import importlib.util
import os
import random
import sys
import tempfile
import unittest

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
spec = importlib.util.spec_from_file_location("tg65", os.path.join(SRC, "65_trb_assembly_genotype.py"))
tg = importlib.util.module_from_spec(spec)
sys.modules["tg65"] = tg
spec.loader.exec_module(tg)

INS_FREQ = {"SAS": 0.45, "EUR": 0.6}
INS = ["TRBV4-3", "TRBV3-2", "TRBV6-3", "TRBV7-8"]
CLS62 = "TRBV6-2/TRBV6-3"


def rs(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mut(rng, s, frac=None, n=None):
    s = list(s)
    k = n if n is not None else int(round(frac * len(s)))
    for i in rng.sample(range(len(s)), k):
        s[i] = rng.choice([c for c in "ACGT" if c != s[i]])
    return "".join(s)


def make_probes(rng):
    p = {g: rs(rng, 300) for g in ["TRBV2", "TRBV3-1", "TRBV12-4", "TRBV4-3", "TRBV3-2", "TRBV7-8",
                                    "TRBV6-2", "TRBV5-1"]}
    p["TRBV6-3"] = mut(rng, p["TRBV6-2"], n=1)
    p["TRBV4-2"] = mut(rng, p["TRBV4-3"], frac=0.05)
    p["TRBC1"] = rs(rng, 1200)
    return p


def build_hap(rng, probes, ins, frag=False, absent=False):
    """-> (contigs [(name, seq)], paf lines). Gene copies carry ~0.3% allele noise."""
    def g(name):
        return mut(rng, probes[name], n=1)

    def sp():
        return rs(rng, rng.randint(300, 700))

    left = [rs(rng, 2000), g("TRBV2"), sp(), g("TRBV3-1"), sp(), g("TRBV4-2"), sp(), g("TRBV6-2"), sp()]
    if ins:
        left += [g("TRBV4-3"), sp(), g("TRBV3-2"), sp(), g("TRBV6-3"), sp(), g("TRBV7-8"), sp()]
    left += [mut(rng, probes["TRBV4-3"], frac=0.05), sp(), g("TRBV5-1"), sp(), g("TRBV12-4"), sp()]
    right = [g("TRBC1"), rs(rng, 2000)]
    others = [("ctg_other1", rs(rng, 3000)), ("ctg_other2", rs(rng, 2500))]
    if absent:
        trb = []
    elif frag:
        trb = [("ctg_trb_a", "".join(left)), ("ctg_trb_b", "".join(right))]
    else:
        trb = [("ctg_trb", "".join(left + right))]
    trb = [(n, tg.revcomp(s) if rng.random() < 0.5 else s) for n, s in trb]
    paf, pos = [], 142_300_000
    for n, s in trb:
        paf.append(f"{n}\t{len(s)}\t0\t{len(s)}\t+\tchr7\t159345973\t{pos}\t{pos + len(s)}\t{len(s)}\t{len(s)}\t60\ttp:A:P")
        pos += len(s)
    for n, s in others:
        paf.append(f"{n}\t{len(s)}\t0\t{len(s)}\t+\tchr1\t248956422\t1000\t{1000 + len(s)}\t{len(s)}\t{len(s)}\t60\ttp:A:P")
    return others[:1] + trb + others[1:], paf


def write_fa_gz(path, contigs):
    with gzip.open(path, "wt") as fh:
        for n, s in contigs:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")


class World:
    def __init__(self, tmp, n, seed=7):
        rng = random.Random(seed)
        self.tmp = tmp
        self.probes = make_probes(rng)
        self.probe_fa = os.path.join(tmp, "probes.fa")
        with open(self.probe_fa, "w") as fh:
            for gname, s in self.probes.items():
                fh.write(f">X{abs(hash(gname)) % 10000}|{gname}*01|Homo sapiens|F|V-REGION|\n{s}\n")
        self.root = os.path.join(tmp, "root")
        os.makedirs(os.path.join(self.root, "v9/wgs/long_read"), exist_ok=True)
        os.makedirs(os.path.join(self.root, "x"), exist_ok=True)
        coh, man = [], []
        self.truth = {}  # pid -> (ins_h1, ins_h2)
        self.anc = {}
        for i in range(n):
            pid = str(2_000_000 + i)
            anc = "SAS" if rng.random() < 0.5 else "EUR"
            q = INS_FREQ[anc]
            h = (rng.random() < q, rng.random() < q)
            self.truth[pid], self.anc[pid] = h, anc
            row = dict(research_id=pid, platform="revio")
            for k in (1, 2):
                ctgs, paf = build_hap(rng, self.probes, h[k - 1], frag=(i == 0 and k == 2), absent=(i == 1 and k == 1))
                fa = f"x/{pid}.hap{k}.fa.gz"
                if i != 3:  # person 3: assembly missing on the mount
                    write_fa_gz(os.path.join(self.root, fa), ctgs)
                with open(os.path.join(self.root, f"x/{pid}.hap{k}.paf"), "w") as fh:
                    fh.write("\n".join(paf) + "\n")
                row[f"assembly_hap{k}_fa"] = "gs://vwb-aou-datasets-controlled/" + fa
                row[f"assembly_hap{k}_aln2_hg38_paf"] = ("" if i == 2 else
                                                         f"gs://vwb-aou-datasets-controlled/x/{pid}.hap{k}.paf")
            man.append(row)
            coh.append(dict(person_id=pid, in_lr=True, ancestry_pred=anc.lower(), platform="revio"))
        self.cohort = os.path.join(tmp, "cohort.tsv")
        pd.DataFrame(coh).to_csv(self.cohort, sep="\t", index=False)
        pd.DataFrame(man).to_csv(os.path.join(self.root, tg.LR_MANIFEST_REL), sep="\t", index=False)

    def run(self, aligner="kmer", out="out", internal="internal", extra=()):
        out, internal = os.path.join(self.tmp, out), os.path.join(self.tmp, internal)
        tg.main(["--probes", self.probe_fa, "--cohort", self.cohort, "--aou-root", self.root,
                 "--out-dir", out, "--internal-dir", internal, "--work-dir", os.path.join(self.tmp, "work"),
                 "--aligner", aligner, "--workers", "4", *extra])
        return out, internal


def read_int(path):
    return pd.read_csv(path, sep="\t", comment="#", dtype={"person_id": str}, keep_default_na=False)


class TestPrimitives(unittest.TestCase):
    def test_semiglobal_and_classes(self):
        self.assertEqual(tg.semiglobal_ed("ACGT", "TTACGTTT"), 0)
        self.assertEqual(tg.semiglobal_ed("ACGT", "TTACTTTT"), 1)
        self.assertEqual(tg.semiglobal_ed("ACGGT", "TTACGTTT"), 1)  # one deletion
        rng = random.Random(3)
        p = make_probes(rng)
        probes = [(f"{g}*01", g, s) for g, s in p.items()]
        cls, tab = tg.gene_classes(probes, 0.99)
        self.assertEqual(cls["TRBV6-2"], CLS62)
        self.assertEqual(cls["TRBV6-3"], CLS62)
        self.assertEqual(cls["TRBV4-3"], "TRBV4-3")  # 95% paralog stays separate
        self.assertEqual(cls["TRBV4-2"], "TRBV4-2")
        self.assertEqual(tg.feature_of_gene("TRBV6-3", cls), CLS62 + tg.COPY2)

    def test_gene_of_headers(self):
        self.assertEqual(tg.gene_of("M1|TRBV4-3*01|Homo sapiens|F|V-REGION|"), ("TRBV4-3*01", "TRBV4-3"))
        self.assertEqual(tg.gene_of("TRBC1*genomic_ENSG1 x"), ("TRBC1*genomic_ENSG1", "TRBC1"))
        page = "<html><pre>>a|TRBV2*01|Homo sapiens|F|\nacg..t\nAC\n>b|TRBV2*02|H|F|\nAAA</pre></html>"
        self.assertEqual(tg.parse_imgt_page(page), [("TRBV2*01", "TRBV2", "ACGTAC"), ("TRBV2*02", "TRBV2", "AAA")])

    def test_minimap2_paf_parse(self):
        txt = "TRBV2*01\t300\t0\t300\t-\tctg\t5000\t100\t400\t297\t300\t60\tNM:i:3\n"
        h = tg.parse_paf_hits(txt)[0]
        self.assertEqual((h.contig, h.t_st, h.t_en, h.strand, h.q_len), ("ctg", 100, 400, -1, 300))
        self.assertAlmostEqual(h.ident, 0.99)

    def test_paralog_not_called_even_at_low_threshold(self):
        """Haplotype WITHOUT the insertion: the TRBV4-2 locus and the probe-less decoy both hit the
        TRBV4-3 probe at ~95%. 0.97 calls neither; at 0.90 best-hit still credits TRBV4-2's locus to
        TRBV4-2 (only the probe-less decoy becomes TRBV4-3 -- the reason the default is 0.97)."""
        rng = random.Random(11)
        p = make_probes(rng)
        probes = [(f"{g}*01", g, s) for g, s in p.items()]
        cls, _ = tg.gene_classes(probes, 0.99)
        tmp = tempfile.mkdtemp()
        ctgs, _ = build_hap(rng, p, ins=False)
        fa = os.path.join(tmp, "h.fa.gz")
        write_fa_gz(fa, ctgs)
        rep = sorted({cls[g] for g in p}, key=tg._natkey)
        anc = [cls["TRBV2"], cls["TRBC1"]]
        for aligner in ["kmer"] + (["mappy"] if tg.mappy is not None else []):
            ad = dict(select="anchor", window=tg.parse_window(tg.DEFAULT_WINDOW), pad=1000,
                      anchors_kmers=tg.anchor_kmers(probes, set(p)), aligner=aligner, preset="sr",
                      minimap2="minimap2", work_dir=tmp, min_ident=0.97, min_cov=0.9, tie_eps=0.001)
            r = tg.genotype_haplotype(fa, None, probes, cls, rep, anc, ad)
            self.assertTrue(r["complete"], aligner)
            self.assertFalse(r["calls"]["TRBV4-3"]["present"], aligner)
            self.assertTrue(r["calls"]["TRBV4-2"]["present"], aligner)
            self.assertEqual(r["calls"][CLS62]["n_loci"], 1, aligner)
            self.assertGreater(r["calls"]["TRBV4-3"]["best_ident"], 0.93, aligner)  # it was seen, then rejected
            ad["min_ident"] = 0.90
            r = tg.genotype_haplotype(fa, None, probes, cls, rep, anc, ad)
            self.assertTrue(r["calls"]["TRBV4-2"]["present"], aligner)
            self.assertEqual(r["calls"]["TRBV4-2"]["n_loci"], 1, aligner)
            self.assertEqual(r["calls"]["TRBV4-3"]["n_loci"], 1, aligner)  # the decoy only


class TestEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.w = World(cls.tmp, n=400)
        cls.out, cls.internal = cls.w.run("kmer")

    def test_recovers_planted_genotypes(self):
        g = read_int(os.path.join(self.internal, "per_person_gene.tsv"))
        ok = g[g["status"] == "ok"]
        for feat in ("TRBV4-3", "TRBV3-2", "TRBV7-8", CLS62 + tg.COPY2):
            v = ok[ok["feature"] == feat].set_index("person_id")
            truth = pd.Series({p: int(h[0]) + int(h[1]) for p, h in self.w.truth.items()})
            self.assertTrue((v["copies"].astype(int) == truth.loc[v.index]).all(), feat)
            self.assertGreater(len(v), 380)
        six = ok[ok["feature"] == CLS62].set_index("person_id")
        for pid, r in six.iterrows():
            h = self.w.truth[pid]
            self.assertEqual((int(r["hap1_n_loci"]), int(r["hap2_n_loci"])), (1 + h[0], 1 + h[1]))
        ctrl = ok[ok["feature"].isin(["TRBV2", "TRBV4-2", "TRBV12-4", "TRBC1", "TRBV3-1"])]
        self.assertTrue((ctrl["copies"].astype(int) == 2).all())

    def test_qc_flags_and_failures(self):
        h = read_int(os.path.join(self.internal, "per_haplotype.tsv"))
        frag = h[(h["person_id"] == "2000000") & (h["hap"].astype(int) == 2)].iloc[0]
        self.assertEqual(frag["incomplete_reason"], "anchors_split")
        absent = h[(h["person_id"] == "2000001") & (h["hap"].astype(int) == 1)].iloc[0]
        self.assertEqual(absent["incomplete_reason"], "no_trb_contig")
        nopaf = h[h["person_id"] == "2000002"]
        self.assertTrue((nopaf["select"] == "anchor").all())
        self.assertTrue((nopaf["complete"].astype(str) == "True").all())
        g = read_int(os.path.join(self.internal, "per_person_gene.tsv"))
        self.assertTrue((g[g["person_id"].isin(["2000000", "2000001"])]["status"] == "incomplete").all())
        f = read_int(os.path.join(self.internal, "failures.tsv"))
        self.assertEqual(list(f["person_id"]), ["2000003"])
        self.assertEqual(f["reason"].iloc[0], "assembly_missing")
        run = dict(pd.read_csv(os.path.join(self.out, "run_summary.tsv"), sep="\t").values)
        self.assertEqual(run["people_genotyped"], "399")
        self.assertEqual(run["people_complete_both_haplotypes"], "397")
        comp = pd.read_csv(os.path.join(self.out, "haplotype_completeness.tsv"), sep="\t", dtype=str)
        allr = comp[comp["stratum"] == "ALL"].iloc[0]
        self.assertEqual(allr["n_haplotypes"], "798")
        self.assertEqual(allr["n_anchors_split"], "<20")

    def test_frequency_hwe_copresence(self):
        tab = pd.read_csv(os.path.join(self.out, "genotype_by_stratum.tsv"), sep="\t")
        for anc, q in INS_FREQ.items():
            r = tab[(tab["feature"] == "TRBV4-3") & (tab["stratum"] == anc)].iloc[0]
            self.assertAlmostEqual(float(r["insertion_allele_freq"]), q, delta=0.08)
            self.assertGreater(float(r["hwe_p"]), 1e-3)
        c = tab[(tab["feature"] == "TRBC1") & (tab["stratum"] == "ALL")].iloc[0]
        self.assertEqual(float(c["insertion_allele_freq"]), 1.0)
        cp = pd.read_csv(os.path.join(self.out, "copresence.tsv"), sep="\t")
        self.assertEqual(set(cp["feature_a"]) | set(cp["feature_b"]), {"TRBV4-3", "TRBV3-2", CLS62 + tg.COPY2})
        allr = cp[cp["stratum"] == "ALL"]
        self.assertTrue((allr["concordance"].astype(float) == 1.0).all())
        self.assertTrue((allr["r2"].astype(float) > 0.99).all())
        pats = pd.read_csv(os.path.join(self.out, "haplotype_patterns.tsv"), sep="\t")
        self.assertEqual(len(pats[pats["stratum"] == "ALL"]), 2)  # all-plus or all-minus only

    def test_public_has_no_person_ids(self):
        for f in os.listdir(self.out):
            txt = open(os.path.join(self.out, f)).read()
            for pid in ("2000000", "2000003", "2000117", "2000399"):
                self.assertNotIn(pid, txt, f)

    def test_resume_adds_nothing(self):
        n0 = len(read_int(os.path.join(self.internal, "per_haplotype.tsv")))
        self.w.run("kmer")
        self.assertEqual(n0, len(read_int(os.path.join(self.internal, "per_haplotype.tsv"))))


class TestSmallAndMappy(unittest.TestCase):
    def test_masking_small_cohort(self):
        tmp = tempfile.mkdtemp()
        w = World(tmp, n=30, seed=5)
        out, _ = w.run("kmer")
        tab = pd.read_csv(os.path.join(out, "genotype_by_stratum.tsv"), sep="\t", dtype=str, keep_default_na=False)
        anc = tab[tab["stratum"].isin(["SAS", "EUR"])]
        self.assertTrue((anc["n_people_measured"] == "<20").all())
        self.assertTrue((anc["insertion_allele_freq"] == "").all())
        self.assertTrue((anc["hwe_p"] == "").all())
        cp = pd.read_csv(os.path.join(out, "copresence.tsv"), sep="\t", dtype=str, keep_default_na=False)
        for c in ("n_both", "n_a_only", "n_b_only", "n_neither"):
            self.assertTrue(cp[c].isin(["0", "<20"] + [str(k) for k in range(20, 100)]).all())

    @unittest.skipIf(tg.mappy is None, "mappy not installed")
    def test_mappy_matches_truth(self):
        tmp = tempfile.mkdtemp()
        w = World(tmp, n=40, seed=9)
        out, internal = w.run("mappy")
        g = read_int(os.path.join(internal, "per_person_gene.tsv"))
        v = g[(g["status"] == "ok") & (g["feature"] == "TRBV4-3")].set_index("person_id")
        truth = pd.Series({p: int(h[0]) + int(h[1]) for p, h in w.truth.items()})
        self.assertGreater(len(v), 30)
        self.assertTrue((v["copies"].astype(int) == truth.loc[v.index]).all())
        self.assertTrue(v["method"].str.startswith("mappy:sr").all())

    def test_internal_refuses_repo_paths(self):
        with self.assertRaises(ValueError):
            tg.append_rows(os.path.join(SRC, "x_internal.tsv"), [], tg.GENE_COLS)
        with self.assertRaises(ValueError):
            tg.main(["--probes", "x.fa", "--work-dir", os.path.join(SRC, "work65"), "--aligner", "kmer"])

    def test_build_probes_route_b_and_reference_check(self):
        tmp = tempfile.mkdtemp()
        rng = random.Random(21)
        p = make_probes(rng)
        pf = os.path.join(tmp, "user.fa")
        with open(pf, "w") as fh:
            for gname, s in p.items():
                fh.write(f">{gname}*01\n{s}\n")
        alt, ref = os.path.join(tmp, "alt.fa"), os.path.join(tmp, "ref.fa")
        write_fa_gz(alt + ".gz", build_hap(rng, p, ins=True)[0][1:2])
        write_fa_gz(ref + ".gz", build_hap(rng, p, ins=False)[0][1:2])
        out = os.path.join(tmp, "out")
        tg.main(["--build-probes", "--probe-fasta", pf, "--alt-fasta", alt + ".gz", "--ref-fasta", ref + ".gz",
                 "--out-dir", out, "--aligner", "kmer", "--work-dir", os.path.join(tmp, "w")])
        self.assertEqual(len(tg.read_probes(os.path.join(out, "probes.fa"))), len(p))
        rc = pd.read_csv(os.path.join(out, "reference_check.tsv"), sep="\t")
        get = {(r.source, r.gene_class): r for r in rc.itertuples()}
        self.assertTrue(get[("alt", "TRBV4-3")].present)
        self.assertFalse(get[("primary", "TRBV4-3")].present)
        self.assertEqual(get[("alt", CLS62)].n_loci, 2)
        self.assertEqual(get[("primary", CLS62)].n_loci, 1)


if __name__ == "__main__":
    unittest.main()
