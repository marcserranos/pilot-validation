#!/usr/bin/env python3
"""Tests for 64_novel_allele_validation.py -- synthetic data only (no samtools, no network, no AoU).

A synthetic HLA-A-like gene (4 exons, one 6-bp poly-A run in exon 2) is placed on per-person
haplotype contigs, with Immuannot-style files (`<hap>.trimmed.fa`, `<hap>.gtf.gz`, `<hap>/cds.fa.gz`)
and simulated HiFi reads returned by a patched read fetcher:
  - a TRUE novel allele (one exon-3 SNV present in the reads)              -> supported
  - a homopolymer ARTIFACT (assembly lost one A; reads carry the known run) -> contradicted
  - a novel allele on a minus-strand contig that a fake newer IPD release
    contains (accession absent from the 3550 history column)               -> now_known + supported
  - a carrier without a BAM, a relative of a carrier, and a third carrier
    of the true-novel key beyond --max-per-key.
Also unit tests of the alignment primitives, IPD parsing, verdict thresholds, masking, absence of
person ids in PUBLIC output, refusal of INTERNAL paths in the repo, and resume.

Run: python3 -m pytest scripts/hla_popgen/tests/test_64_novel_allele_validation.py -q
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
spec = importlib.util.spec_from_file_location("nav64", os.path.join(SRC, "64_novel_allele_validation.py"))
nav = importlib.util.module_from_spec(spec)
sys.modules["nav64"] = nav
spec.loader.exec_module(nav)

HAVE_MAPPY = nav.mappy is not None
EXONS = [(2601, 2700), (3001, 3270), (3601, 3876), (4501, 4700)]  # 1-based on the base haplotype
GENE = (2501, 5500)
RUN = (3101, 3106)          # AAAAAA in exon 2
SNV_TRUE = 3700             # exon 3
SNV_NEW = 3200              # exon 2
HETS = [1000, 1800, 2400, 2800, 3400, 4100, 4300, 5200, 6000, 6600, 7000, 7600, 8000]
L = 9000


def _rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _sub(seq, pos1, base=None):
    old = seq[pos1 - 1]
    new = base or {"A": "C", "C": "G", "G": "T", "T": "A"}[old]
    return seq[:pos1 - 1] + new + seq[pos1:]


def make_base(seed=7):
    rng = random.Random(seed)
    s = list(_rand_seq(rng, L))
    s[EXONS[0][0] - 1:EXONS[0][0] + 2] = list("ATG")
    s[EXONS[-1][1] - 3:EXONS[-1][1]] = list("TGA")
    s[RUN[0] - 2] = "C"
    for p in range(RUN[0], RUN[1] + 1):
        s[p - 1] = "A"
    s[RUN[1]] = "G"
    return "".join(s)


def cds_of(seq, exons):
    return "".join(seq[s - 1:e] for s, e in exons)


def shift_exons(exons, after, delta):
    return [(s + (delta if s > after else 0), e + (delta if e >= after else 0)) for s, e in exons]


def with_hets(seq):
    for p in HETS:
        seq = _sub(seq, p)
    return seq


def sim_reads(rng, seq, n, rlen=4000, err=0.001, tag="r"):
    out = []
    for i in range(n):
        st = rng.randint(0, len(seq) - rlen)
        r = list(seq[st:st + rlen])
        for k in range(len(r)):
            if rng.random() < err:
                r[k] = rng.choice([b for b in "ACGT" if b != r[k]])
        r = "".join(r)
        if rng.random() < 0.5:
            r = nav.revcomp(r)
        out.append((f"{tag}{i}", r))
    return out


class World(object):
    """Builds every synthetic input on disk; `reads[bam_path]` feeds the patched fetcher."""

    def __init__(self, tmp):
        self.tmp = tmp
        self.rng = random.Random(11)
        self.base = make_base()
        self.outroot = os.path.join(tmp, "people")
        self.root = os.path.join(tmp, "aou")
        os.makedirs(os.path.join(self.root, "v9/wgs/long_read"), exist_ok=True)
        os.makedirs(os.path.join(self.root, "lr"), exist_ok=True)
        self.t1, self.cohort, self.man, self.reads = [], [], [], {}
        self.near_cds = cds_of(self.base, EXONS)
        b = self.base
        true_seq = _sub(b, SNV_TRUE)
        new_seq = _sub(b, SNV_NEW)
        art_asm = b[:RUN[0] - 1] + b[RUN[0]:]          # one A lost in the assembly
        hap2 = with_hets(b)
        cm_true = "HLA-A*01:01:01:01|:600*ag:245|"
        cm_art = "HLA-A*01:01:01:01|:500-a:345|"
        cm_new = "HLA-A*01:01:01:01|:400*tc:445|"
        self.new_cds = cds_of(new_seq, EXONS)
        # person: (pid, ancestry, platform, hap1 assembly, hap1 truth, exons1, gene1, cds_mut, minus, bam?)
        spec_ = [
            ("1000000", "eur", "revio", true_seq, true_seq, EXONS, GENE, cm_true, False, True),  # relative
            ("1000001", "eur", "revio", true_seq, true_seq, EXONS, GENE, cm_true, False, True),
            ("1000002", "afr", "sequel2e", true_seq, true_seq, EXONS, GENE, cm_true, False, True),
            ("1000003", "afr", "revio", true_seq, true_seq, EXONS, GENE, cm_true, False, True),
            ("1000004", "eur", "revio", art_asm, b, shift_exons(EXONS, RUN[0], -1), (GENE[0], GENE[1] - 1),
             cm_art, False, True),
            ("1000005", "sas", "revio", art_asm, b, shift_exons(EXONS, RUN[0], -1), (GENE[0], GENE[1] - 1),
             cm_art, False, False),                                                          # no BAM
            ("1000006", "eas", "revio", new_seq, new_seq, EXONS, GENE, cm_new, True, True),   # minus strand
        ]
        for pid, anc, plat, asm1, truth1, ex1, g1, cm, minus, has_bam in spec_:
            self._person(pid, asm1, ex1, g1, cm, minus, hap2)
            self.cohort.append(dict(person_id=pid, in_lr=True, ancestry_pred=anc, platform=plat))
            uri = f"gs://vwb-aou-datasets-controlled/lr/{pid}.bam"
            self.man.append(dict(research_id=pid, platform=plat, grch38_bam=uri))
            bam = os.path.join(self.root, "lr", f"{pid}.bam")
            if has_bam:
                open(bam, "w").close()
            self.reads[bam] = (sim_reads(self.rng, truth1, 14, tag=f"{pid}h1_")
                               + sim_reads(self.rng, hap2, 14, tag=f"{pid}h2_"))
        self.paths = self._write_tables()
        self.ipd_dir = self._write_ipd()

    def _person(self, pid, asm1, ex1, g1, cm, minus, hap2):
        d = os.path.join(self.outroot, pid, "immuannot_output")
        os.makedirs(os.path.join(d, "hap1"), exist_ok=True)
        os.makedirs(os.path.join(d, "hap2"), exist_ok=True)
        for hap, seq, ex, g, novel in (("hap1", asm1, ex1, g1, True), ("hap2", hap2, EXONS, GENE, False)):
            ctg = f"{hap[0]}{hap[-1]}tg000001l"
            strand = "+"
            cds = cds_of(seq, ex)
            if hap == "hap1" and minus:
                n = len(seq)
                seq = nav.revcomp(seq)
                ex = [(n - e + 1, n - s + 1) for s, e in ex]
                g = (n - g[1] + 1, n - g[0] + 1)
                strand = "-"
            with open(os.path.join(d, f"{hap}.trimmed.fa"), "w") as fh:
                fh.write(f">{ctg}\n{seq}\n")
            with gzip.open(os.path.join(d, f"{hap}.gtf.gz"), "wt") as fh:
                fh.write("## gene (copy num = 1): HLA-A\n")
                fh.write(f'{ctg}\tIPD\tgene\t{g[0]}\t{g[1]}\t.\t{strand}\t.\tgene_id "IAG1"; '
                         f'template_allele "HLA-A*01:01:01:01"; template_distance 1; gene_name "HLA-A";\n')
                for s, e in ex:
                    fh.write(f'{ctg}\tIPD\tCDS\t{s}\t{e}\t.\t{strand}\t0\tgene_id "IAG1"; '
                             f'transcript_id "IAT1"; gene_name "HLA-A";\n')
            with gzip.open(os.path.join(d, hap, "cds.fa.gz"), "wt") as fh:
                fh.write(f">{ctg}_HLA-A_1 GT_AG\n{cds}\n")
            self.t1.append(dict(person_id=pid, hap=hap, contig=ctg, gene="HLA-A", copy_index=1,
                                consensus="HLA-A*01:new" if novel else "HLA-A*01:01:01:01",
                                novelty_depth=2 if novel else "", gene_start=g[0], gene_end=g[1],
                                strand=strand, cds_mut=cm if novel else "", template_warning="NA",
                                template_allele="HLA-A*01:01:01:01", alleles="HLA-A*01:01:01:01"))

    def _write_tables(self):
        p = {}
        p["table1"] = os.path.join(self.tmp, "t1.tsv")
        pd.DataFrame(self.t1).to_csv(p["table1"], sep="\t", index=False)
        p["cohort"] = os.path.join(self.tmp, "cohort.tsv")
        pd.DataFrame(self.cohort).to_csv(p["cohort"], sep="\t", index=False)
        p["rel"] = os.path.join(self.tmp, "rel.tsv")
        pd.DataFrame([{"i.s": "1000000", "j.s": "1000001", "kin": 0.25}]).to_csv(p["rel"], sep="\t", index=False)
        pd.DataFrame(self.man).to_csv(os.path.join(self.root, nav.LR_MANIFEST_REL), sep="\t", index=False)
        p["bed"] = os.path.join(self.tmp, "regions.bed")
        with open(p["bed"], "w") as fh:
            fh.write("chr6\t29940000\t29946000\tHLA-A\tgene_HLA\n")
        return p

    def _write_ipd(self):
        d = os.path.join(self.tmp, "ipd")
        os.makedirs(d, exist_ok=True)
        other = _rand_seq(random.Random(3), 846)
        with open(os.path.join(d, "hla_nuc.fasta"), "w") as fh:
            fh.write(f">HLA:HLA00001 A*01:01:01:01 846 bp\n{self.near_cds}\n")
            fh.write(f">HLA:HLA00002 A*02:01:01:01 846 bp\n{other}\n")
            fh.write(f">HLA:HLA99999 A*01:999 846 bp\n{self.new_cds}\n")
        with open(os.path.join(d, "Allelelist.txt"), "w") as fh:
            fh.write("# file: Allelelist.txt\n# version: IPD-IMGT/HLA 3.99.0\nAlleleID,Allele\n")
        with open(os.path.join(d, "Allelelist_history.txt"), "w") as fh:
            fh.write("# file: Allelelist_history.txt\nHLA_ID,3990,3550\n"
                     "HLA00001,A*01:01:01:01,A*01:01:01:01\nHLA00002,A*02:01:01:01,A*02:01:01:01\n"
                     "HLA99999,A*01:999,NA\n")
        return d

    def args(self, out, internal, extra=()):
        return ["--ipd-dir", self.ipd_dir, "--refdata", os.path.join(self.tmp, "none"),
                "--table1", self.paths["table1"], "--cohort", self.paths["cohort"],
                "--relatedness", self.paths["rel"], "--outroot", self.outroot,
                "--regions-bed", self.paths["bed"], "--aou-root", self.root, "--out-dir", out,
                "--internal-dir", internal, "--workers", "2", "--max-per-key", "2"] + list(extra)


class FakeFetch(object):
    def __init__(self, reads):
        self.reads = reads
        self.calls = 0

    def __call__(self, path, region, timeout, samtools="samtools", ref=None):
        self.calls += 1
        return list(self.reads[path])


# ---------------------------------------------------------------------------
class TestPrimitives(unittest.TestCase):
    def test_align_events_roundtrip(self):
        rng = random.Random(1)
        a = _rand_seq(rng, 600)
        b = a[:100] + "T" + a[101:300] + a[302:450] + "GGG" + a[450:]
        if a[100] == "T":
            b = a[:100] + "C" + b[101:]
        ev = nav.align_events(a, b)
        self.assertEqual(nav.apply_events(a, ev), b)
        self.assertLessEqual(len(ev), 3)

    def test_event_hp(self):
        seq = "CCGTAAAAAGTC"
        hp = nav.event_hp(seq, (5, 1, ""), 3)       # called has one extra A: run 5 vs nearest 4
        self.assertEqual((hp["called_len"], hp["nearest_len"], hp["base"]), (5, 4, "A"))
        hp2 = nav.event_hp(seq, (9, 0, "A"), 3)     # nearest has one more A
        self.assertEqual((hp2["called_len"], hp2["nearest_len"]), (5, 6))
        self.assertIsNone(nav.event_hp(seq, (2, 1, "T"), 3))       # SNV
        self.assertIsNone(nav.event_hp("CGTACG", (2, 1, ""), 3))   # lone base, no run

    def test_cswalk(self):
        w = nav.CsWalk(":10*ag:5+tt:5-c:4", 100, 0)
        self.assertEqual(w.r_en, 126)
        self.assertEqual(w.q_at(100), 0)
        self.assertEqual(w.q_at(116), 18)       # boundary after the insertion
        self.assertEqual(w.q_at(121), 23)       # deleted target base -> query boundary
        self.assertEqual(w.edits_in(0, 30), 4)  # 1 sub + 2 ins + 1 del

    def test_lev(self):
        self.assertEqual(nav.lev("ACGTACGT", "ACGTACGT"), 0)
        self.assertEqual(nav.lev("ACGTAAAACGT", "ACGTAAACGT"), 1)
        self.assertEqual(nav.lev("AAAA", ""), 4)

    def test_verdict(self):
        base = dict(n_spanning=10, this_called=9, this_nearest=0, pool_called=12, pool_nearest=9)
        self.assertEqual(nav.verdict(base, True, 3, 0.7, 0.3, 0.15)[4], "supported")
        art = dict(n_spanning=10, this_called=0, this_nearest=9, pool_called=0, pool_nearest=20)
        self.assertEqual(nav.verdict(art, True, 3, 0.7, 0.3, 0.15)[4], "contradicted")
        low = dict(n_spanning=2, this_called=1, this_nearest=0, pool_called=1, pool_nearest=1)
        self.assertEqual(nav.verdict(low, True, 3, 0.7, 0.3, 0.15)[4:], ("ambiguous", "low_reads"))
        het = dict(n_spanning=20, this_called=0, this_nearest=0, pool_called=8, pool_nearest=9)
        v = nav.verdict(het, False, 3, 0.7, 0.3, 0.15)
        self.assertEqual((v[0], v[4]), ("2t", "supported"))
        mixed = dict(n_spanning=20, this_called=4, this_nearest=6, pool_called=4, pool_nearest=6)
        self.assertEqual(nav.verdict(mixed, True, 3, 0.7, 0.3, 0.15)[4], "ambiguous")

    def test_pooled_key_read_class(self):
        import pandas as pd
        def calls(rows):
            return pd.DataFrame(rows, columns=["lr_mode", "lr_verdict", "n_reads_support_called",
                                               "n_reads_support_nearest"])
        # three carriers, each ambiguous alone (2 reads), together a clear supported key
        real = calls([("3t", "ambiguous", 2, 0)] * 3)
        self.assertEqual(nav.pooled_key_read_class(real)[3], "supported")
        noise = calls([("3t", "ambiguous", 1, 2), ("3t", "ambiguous", 1, 1), ("3t", "ambiguous", 0, 2)])
        self.assertEqual(nav.pooled_key_read_class(noise)[3], "ambiguous")
        art = calls([("3t", "ambiguous", 0, 2)] * 3)
        self.assertEqual(nav.pooled_key_read_class(art)[3], "contradicted")
        # 2t calls are not pooled
        self.assertEqual(nav.pooled_key_read_class(calls([("2t", "supported", 8, 9)]))[3], "untested")

    def test_headers_and_history(self):
        self.assertEqual(nav.parse_allele_header("HLA:HLA00001 A*01:01:01:01 1098 bp"),
                         ("A*01:01:01:01", "A", "HLA00001"))
        self.assertEqual(nav.parse_allele_header("HLA-B*07:02:01:01 HLA00132 frame=1 1089bp"),
                         ("B*07:02:01:01", "B", "HLA00132"))
        tmp = tempfile.mkdtemp()
        p = os.path.join(tmp, "h.txt")
        with open(p, "w") as fh:
            fh.write("# x\nHLA_ID\t3600\t3550\nHLA00001\tA*01\tA*01\nHLA00009\tA*09\t\n")
        self.assertEqual(nav.load_baseline_history(p, "3550"), {"HLA00001"})
        self.assertIsNone(nav.load_baseline_history(p, "3400"))

    def test_fetch_ipd_lfs_fallback(self):
        seen = []

        def fake(url, timeout=300):
            seen.append(url)
            if "raw.githubusercontent" in url and url.endswith("hla_gen.fasta"):
                return b"version https://git-lfs.github.com/spec/v1\noid sha256:x\n"
            return b">HLA:HLA00001 A*01:01 3 bp\nATG\n"

        orig = nav.fetch_url
        nav.fetch_url = fake
        try:
            d = tempfile.mkdtemp()
            nav.fetch_ipd("Latest", d, files=["hla_nuc.fasta", "hla_gen.fasta"])
        finally:
            nav.fetch_url = orig
        self.assertIn("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_nuc.fasta", seen)
        self.assertIn("https://media.githubusercontent.com/media/ANHIG/IMGTHLA/Latest/hla_gen.fasta", seen)
        self.assertTrue(open(os.path.join(d, "hla_gen.fasta")).read().startswith(">HLA"))


class TestDepth4(unittest.TestCase):
    @unittest.skipUnless(HAVE_MAPPY, "mappy not installed")
    def test_genomic_events(self):
        rng = random.Random(5)
        wc = _rand_seq(rng, 3000)
        gen = wc[500:1000] + ("A" if wc[1000] != "A" else "C") + wc[1001:1500] + wc[1501:2500]
        ev, why = nav.genomic_events(wc, gen)
        self.assertIsNone(why)
        self.assertEqual(len(ev), 2)
        self.assertEqual(nav.apply_events(wc, ev)[500:500 + len(gen)], gen)

    def test_depth4_classify(self):
        rng = random.Random(9)
        old, new = _rand_seq(rng, 3000), _rand_seq(rng, 3000)
        d = tempfile.mkdtemp()
        with open(os.path.join(d, "hla_nuc.fasta"), "w") as fh:
            fh.write(">HLA:HLA00001 A*01:01:01:01 9 bp\nATGAAATGA\n")
        with open(os.path.join(d, "hla_gen.fasta"), "w") as fh:
            fh.write(f">HLA:HLA00001 A*01:01:01:01 3000 bp\n{old}\n>HLA:HLA99999 A*01:01:01:99 3000 bp\n{new}\n")
        ipd = nav.IpdRelease(d, {"HLA00001"}, "test", load_gen=True)
        self.assertEqual(ipd.classify("A", 4, new[100:2900], None), ("now_known", ["A*01:01:01:99"]))
        self.assertEqual(ipd.classify("A", 4, "GG" + old + "TT", None)[0], "known_in_baseline")
        self.assertEqual(ipd.classify("A", 4, _rand_seq(rng, 2000), None)[0], "still_novel")
        self.assertEqual(ipd.nearest_gen("A", "HLA-A*01:01:01:01")[2], "ipd_current")

    def test_depth2_protein_and_partial_classes(self):
        d = tempfile.mkdtemp()
        near = "ATGGCTGCTAAACCCGGGTTTTGA"                 # M A A K P G F *
        syn_new = "ATGGCAGCTAAACCCGGGTTTTGA"              # same protein, new CDS (A*01:02 new)
        obs = "ATGGCAGCTAAACCCGGGTTCTGA"                  # F->F synonymous twice -> protein known
        with open(os.path.join(d, "hla_nuc.fasta"), "w") as fh:
            fh.write(f">HLA:HLA00001 A*01:01 24 bp\n{near}\n>HLA:HLA99998 A*01:02 24 bp\n{syn_new}\n")
        ipd = nav.IpdRelease(d, {"HLA00001"}, "test")
        self.assertEqual(ipd.classify("A", 2, nav.cds_core(obs), nav.cds_core(near))[0], "protein_known_in_baseline")
        self.assertEqual(ipd.classify("A", 3, nav.cds_core(syn_new), nav.cds_core(near))[0], "now_known")
        ipd2 = nav.IpdRelease(d, None, "none")
        self.assertEqual(ipd2.classify("A", 3, nav.cds_core(syn_new), None)[0], "catalogued_age_unknown")


@unittest.skipUnless(HAVE_MAPPY, "mappy not installed")
class TestEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.w = World(cls.tmp)
        cls.orig = nav.fetch_reads
        cls.fetch = FakeFetch(cls.w.reads)
        nav.fetch_reads = cls.fetch
        cls.out = os.path.join(cls.tmp, "out")
        cls.internal = os.path.join(cls.tmp, "internal")
        nav.main(cls.w.args(cls.out, cls.internal))
        cls.per = pd.read_csv(os.path.join(cls.internal, "per_call_evidence.tsv"), sep="\t", comment="#",
                              dtype=str, keep_default_na=False)

    @classmethod
    def tearDownClass(cls):
        nav.fetch_reads = cls.orig

    def row(self, pid):
        r = self.per[self.per["person_id"] == pid]
        self.assertEqual(len(r), 1, pid)
        return r.iloc[0]

    def test_selection(self):
        pids = set(self.per["person_id"])
        self.assertNotIn("1000000", pids)                      # relative removed (greedy unrelated)
        true_key = self.per[self.per["key"].str.contains(":600\\*ag")]
        self.assertEqual(len(true_key), 2)                     # 3 unrelated carriers, --max-per-key 2
        self.assertEqual(set(true_key["platform"]), {"revio", "sequel2e"})   # stratified pick

    def test_true_novel_supported(self):
        for pid in self.per.loc[self.per["key"].str.contains(":600\\*ag"), "person_id"]:
            r = self.row(pid)
            self.assertEqual(r["status"], "ok")
            self.assertEqual(r["cds_check"], "agree")
            self.assertEqual(r["ipd_class"], "still_novel")
            self.assertEqual(r["lr_verdict"], "supported", r.to_dict())
            self.assertEqual(r["lr_mode"], "3t")
            self.assertGreaterEqual(float(r["frac_called"]), 0.9)
            self.assertGreater(int(r["lr_n_other"]), 0)        # other-haplotype reads recognised
            self.assertEqual(r["hp_only_seq"], "False")

    def test_homopolymer_artifact_contradicted(self):
        r = self.row("1000004")
        self.assertEqual(r["hp_only_proxy"], "True")
        self.assertEqual(r["hp_only_seq"], "True")
        self.assertEqual(r["lr_verdict"], "contradicted", r.to_dict())
        self.assertLessEqual(float(r["frac_called"]), 0.1)
        self.assertGreater(int(r["n_reads_support_nearest"]), 3)
        self.assertGreater(int(r["hp_n_nearest_len"]), int(r["hp_n_called_len"]))

    def test_missing_bam(self):
        r = self.row("1000005")
        self.assertEqual(r["status"], "ok")
        self.assertEqual(r["lr_reason"], "bam_missing")
        self.assertEqual(r["lr_verdict"], "not_run")

    def test_now_known_minus_strand(self):
        r = self.row("1000006")
        self.assertEqual(r["ipd_class"], "now_known")
        self.assertIn("A*01:999", r["ipd_match"])
        self.assertEqual(r["cds_check"], "agree")
        self.assertEqual(r["lr_verdict"], "supported", r.to_dict())
        nk = pd.read_csv(os.path.join(self.out, "now_known_alleles.tsv"), sep="\t", dtype=str)
        self.assertEqual(list(nk["ipd_alleles"]), ["A*01:999"])
        self.assertEqual(nk["n_unrelated_carriers"].iloc[0], "<20")
        self.assertEqual(nk["review_flag"].iloc[0], "True")

    def test_public_tables(self):
        s = pd.read_csv(os.path.join(self.out, "summary_by_stratum.tsv"), sep="\t", dtype=str,
                        keep_default_na=False)
        allrow = s[(s["gene"] == "A") & (s["platform"] == "ALL") & (s["ancestry"] == "ALL") & (s["hp_split"] == "all")]
        self.assertEqual(len(allrow), 1)
        a = allrow.iloc[0]
        self.assertEqual(a["n_keys"], "3")                    # distinct-key counts exact
        self.assertEqual(a["n_calls"], "<20")                 # call counts masked
        self.assertEqual(a["pct_supported"], "")              # rate over a <20 cell blanked
        self.assertEqual(float(a["pct_keys_now_known"]), round(100 / 3, 2))
        self.assertTrue({"hp_only", "other"} <= set(s["hp_split"]))
        self.assertIn("classical_8", set(s["gene"]))
        f = pd.read_csv(os.path.join(self.out, "funnel_with_evidence.tsv"), sep="\t", dtype=str)
        fa = f[(f["gene"] == "A") & (f["novelty_class"] == "protein")].set_index("step")
        self.assertEqual(fa.loc["distinct", "n_keys"], "3")
        self.assertEqual(fa.loc["ge2_unrelated_carriers", "n_keys"], "2")   # true key (3) + artifact key (2)
        self.assertEqual(fa.loc["not_homopolymer_only", "n_keys"], "1")
        self.assertEqual(fa.loc["distinct", "n_read_contradicted"], "1")
        self.assertEqual(fa.loc["distinct", "n_now_known"], "1")
        self.assertEqual(fa.loc["ge2_unrelated_carriers", "n_read_supported"], "1")
        rel = pd.read_csv(os.path.join(self.out, "ipd_release.tsv"), sep="\t", dtype=str)
        self.assertIn("IPD-IMGT/HLA 3.99.0", set(rel["value"]))
        self.assertEqual(dict(zip(rel["item"], rel["value"]))["gene_A"], "3 alleles; 1 added since baseline")

    def test_no_person_ids_in_public(self):
        for fn in os.listdir(self.out):
            txt = open(os.path.join(self.out, fn)).read()
            for pid in ("1000000", "1000001", "1000002", "1000003", "1000004", "1000005", "1000006"):
                self.assertNotIn(pid, txt, fn)
            self.assertNotIn("tg000001l", txt, fn)              # no contig names either

    def test_internal_mode_keeps_exact(self):
        out2 = os.path.join(self.tmp, "out_internal_mode")
        nav.main(self.w.args(out2, self.internal, ["--aggregate-only", "--disclosure", "internal"]))
        s = pd.read_csv(os.path.join(out2, "summary_by_stratum.tsv"), sep="\t", dtype=str, comment="#",
                        keep_default_na=False)
        a = s[(s["gene"] == "A") & (s["platform"] == "ALL") & (s["ancestry"] == "ALL") & (s["hp_split"] == "all")].iloc[0]
        self.assertEqual(a["n_calls"], "6")          # 6 unrelated carriers' calls
        self.assertEqual(a["n_calls_tested"], "4")   # 5 selected, one without a BAM
        self.assertEqual(a["n_supported"], "3")
        self.assertEqual(a["n_contradicted"], "1")

    def test_resume(self):
        n0 = len(self.per)
        calls0 = self.fetch.calls
        nav.main(self.w.args(self.out, self.internal))
        per = pd.read_csv(os.path.join(self.internal, "per_call_evidence.tsv"), sep="\t", comment="#", dtype=str)
        self.assertEqual(len(per), n0)
        self.assertEqual(self.fetch.calls, calls0)

    def test_internal_refuses_repo(self):
        with self.assertRaises(ValueError):
            nav.main(self.w.args(os.path.join(self.tmp, "o3"), os.path.join(SRC, "internal_64_test")))
        with self.assertRaises(ValueError):
            nav.append_rows(os.path.join(SRC, "x_internal_64.tsv"), [], nav.PER_CALL_COLS)
        self.assertFalse(os.path.exists(os.path.join(SRC, "internal_64_test")))

    def test_no_reads_mode(self):
        out = os.path.join(self.tmp, "o_noreads")
        internal = os.path.join(self.tmp, "i_noreads")
        calls0 = self.fetch.calls
        nav.main(self.w.args(out, internal, ["--no-reads"]))
        self.assertEqual(self.fetch.calls, calls0)
        per = pd.read_csv(os.path.join(internal, "per_call_evidence.tsv"), sep="\t", comment="#", dtype=str,
                          keep_default_na=False)
        self.assertEqual(set(per["lr_verdict"]), {"not_run"})
        self.assertIn("now_known", set(per["ipd_class"]))


if __name__ == "__main__":
    unittest.main()
