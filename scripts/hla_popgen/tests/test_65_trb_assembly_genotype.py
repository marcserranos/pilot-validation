#!/usr/bin/env python3
"""Tests for 65_trb_assembly_genotype.py -- synthetic assemblies only, no AoU data, no network.

Synthetic world: random-sequence TRB probes; TRBV6-3 is a 1-substitution copy of TRBV6-2 (must merge
into one class), TRBV4-2 is a 95%-identity paralog of TRBV4-3, and every haplotype also carries a
second, probe-less 95% near-copy of TRBV4-3 (the decoy). Neither may ever be called TRBV4-3.
Insertion haplotypes carry TRBV4-3/3-2/6-3/7-8 together, at planted per-ancestry frequencies.
One haplotype is fragmented (anchors on two contigs), one lacks the TRB locus, one person has no PAF
(anchor contig selection), one person's assembly is missing (failure).

Reads mode: a fixed synthetic primary TRB window (deletion allele) and its insertion version (the
"ALT"); per-person long reads sampled from both haplotypes and "aligned" to the primary the way a
long-read aligner places them (contiguous on the deletion allele; clipped at the breakpoint or with
one large I op on the insertion allele; reads entirely inside the insertion unmapped), served by a
patched samtools runner. Also: breakpoint discovery at probe build, stratified sampling, 65b compare.

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

    def run(self, aligner="kmer", out="out", internal="internal", extra=(), mode="assembly"):
        out, internal = os.path.join(self.tmp, out), os.path.join(self.tmp, internal)
        tg.main(["--mode", mode, "--probes", self.probe_fa, "--cohort", self.cohort, "--aou-root", self.root,
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
        # estimate checks read rates that 66's cross-table pass blanks on these tiny synthetic strata
        # (1-19 people not called everywhere); TestReadsEndToEnd.test_cross_table_pass covers the pass
        cls.out, cls.internal = cls.w.run("kmer", extra=("--no-cross-table",))

    def test_recovers_planted_genotypes(self):
        g = read_int(os.path.join(self.internal, "per_person_gene.tsv"))
        ok = g[g["status"] == "ok"]
        for feat in ("TRBV4-3", "TRBV3-2", "TRBV7-8", CLS62 + tg.COPY2):
            v = ok[ok["feature"] == feat].set_index("person_id")
            truth = pd.Series({p: int(h[0]) + int(h[1]) for p, h in self.w.truth.items()})
            self.assertTrue((v["copies"].astype(int) == truth.loc[v.index]).all(), feat)
            self.assertGreater(len(v), 390)
        six = ok[ok["feature"] == CLS62].set_index("person_id")
        for pid, r in six.iterrows():
            h = self.w.truth[pid]
            self.assertEqual((int(r["hap1_n_loci"]), int(r["hap2_n_loci"])), (1 + h[0], 1 + h[1]))
        ctrl = ok[ok["feature"].isin(["TRBV2", "TRBV4-2", "TRBV12-4", "TRBC1", "TRBV3-1"])]
        self.assertTrue((ctrl["copies"].astype(int) == 2).all())

    def test_qc_flags_and_failures(self):
        h = read_int(os.path.join(self.internal, "per_haplotype.tsv"))
        frag = h[(h["person_id"] == "2000000") & (h["hap"].astype(int) == 2)].iloc[0]
        self.assertEqual(str(frag["complete"]), "True")  # anchors on two contigs: complete collectively
        self.assertEqual(str(frag["complete_single_contig"]), "False")
        absent = h[(h["person_id"] == "2000001") & (h["hap"].astype(int) == 1)].iloc[0]
        self.assertEqual(absent["incomplete_reason"], "no_trb_contig")
        nopaf = h[h["person_id"] == "2000002"]
        self.assertTrue((nopaf["select"] == "anchor").all())
        self.assertTrue((nopaf["complete"].astype(str) == "True").all())
        g = read_int(os.path.join(self.internal, "per_person_gene.tsv"))
        self.assertTrue((g[g["person_id"] == "2000001"]["status"] == "incomplete").all())
        self.assertTrue((g[g["person_id"] == "2000000"]["status"] == "ok").all())
        f = read_int(os.path.join(self.internal, "failures.tsv"))
        self.assertEqual(list(f["person_id"]), ["2000003"])
        self.assertEqual(f["reason"].iloc[0], "assembly_missing")
        run = dict(pd.read_csv(os.path.join(self.out, "run_summary.tsv"), sep="\t").values)
        self.assertEqual(run["people_genotyped"], "399")
        # 399 genotyped - 398 called = 1 person: the called count is masked complementarily (S05 10-12)
        self.assertEqual(run["people_called_any_feature"], "<20*")
        comp = pd.read_csv(os.path.join(self.out, "haplotype_completeness.tsv"), sep="\t", dtype=str)
        allr = comp[comp["stratum"] == "ALL"].iloc[0]
        self.assertEqual(allr["n_haplotypes"], "798")
        # 798 = complete + <20 incomplete: the complete count is masked too, or the <20 cell is a
        # subtraction away
        self.assertEqual(allr["n_complete"], "<20*")
        self.assertEqual(allr["pct_complete"] if isinstance(allr["pct_complete"], str) else "", "")
        self.assertEqual(allr["n_complete_single_contig"], "796")
        self.assertEqual(allr["n_no_trb_contig"], "<20")

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
        self.w.run("kmer", extra=("--no-cross-table",))   # rewrites the shared self.out
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
        self.assertTrue(v["method"].str.startswith("assembly:mappy:sr").all())

    def test_internal_refuses_repo_paths(self):
        with self.assertRaises(ValueError):
            tg.append_rows(os.path.join(SRC, "x_internal.tsv"), [], tg.GENE_COLS)
        with self.assertRaises(ValueError):
            tg.main(["--probes", "x.fa", "--work-dir", os.path.join(SRC, "work65"), "--aligner", "kmer"])
        with self.assertRaises(ValueError):
            tg.main(["--internal-dir", os.path.join(SRC, "int65"), "--people-file", os.path.join(SRC, "p.tsv"),
                     "--sample-stratified", "5", "--aligner", "kmer", "--work-dir", tempfile.mkdtemp()])

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


# ---------------------------------------------------------------------------
# Reads mode
# ---------------------------------------------------------------------------
REF_OFF0 = 142_200_000  # synthetic primary window starts at chr7:142,200,001


class Locus:
    """Fixed synthetic TRB window: primary = L + R (deletion allele); alt = L + INS + R."""

    def __init__(self, seed=31):
        rng = random.Random(seed)
        self.probes = make_probes(rng)
        p = self.probes

        def sp(n=None):
            return rs(rng, n or rng.randint(300, 600))
        self.L = "".join([rs(rng, 3000), p["TRBV2"], sp(), p["TRBV3-1"], sp(), p["TRBV4-2"], sp(), p["TRBV6-2"], sp(1500)])
        self.INS = "".join([sp(800), p["TRBV4-3"], sp(), p["TRBV3-2"], sp(), p["TRBV6-3"], sp(), p["TRBV7-8"], sp(800)])
        self.R = "".join([sp(1500), mut(rng, p["TRBV4-3"], frac=0.05), sp(), p["TRBV5-1"], sp(), p["TRBV12-4"], sp(),
                          p["TRBC1"], rs(rng, 3000)])
        self.primary = self.L + self.R
        self.alt = self.L + self.INS + self.R
        self.B, self.I = len(self.L), len(self.INS)
        self.bp_left = REF_OFF0 + self.B          # 1-based last base of the left flank
        self.breakpoint = f"chr7:{self.bp_left}-{self.bp_left + 1}"

    def reads(self, rng, ins, depth, rlen, hp, name):
        """SAM lines for one haplotype, placed on the primary as an aligner would."""
        hap = self.alt if ins else self.primary
        n = int(depth * len(hap) / rlen)
        out = []
        B, I = self.B, self.I
        for k in range(n):
            st = rng.randint(-rlen + 1500, len(hap) - 1500)
            s, e = max(0, st), min(len(hap), st + rlen)
            if not ins:
                pos, cig = s, f"{e - s}M"
            elif e <= B:
                pos, cig = s, f"{e - s}M"
            elif s >= B + I:
                pos, cig = s - I, f"{e - s}M"
            elif s < B and e > B + I:
                pos, cig = s, f"{B - s}M{I}I{e - B - I}M"
            elif s < B:  # ends inside the insertion: clipped on the left flank
                if B - s < 200:
                    continue
                pos, cig = s, f"{B - s}M{e - B}S"
            elif e > B + I:  # starts inside the insertion: clipped on the right flank
                if e - B - I < 200:
                    continue
                pos, cig = B, f"{B + I - s}S{e - B - I}M"
            else:
                continue  # entirely inside the insertion: unmapped
            tag = f"\tHP:i:{hp}" if rng.random() < 0.9 else ""
            out.append(f"{name}_{k}\t0\tchr7\t{REF_OFF0 + pos + 1}\t60\t{cig}\t*\t0\t0\t{hap[s:e]}\t*{tag}")
        return out


class ReadsWorld:
    def __init__(self, tmp, n, seed=13, depth=15, rlen=6000):
        rng = random.Random(seed)
        self.tmp, self.loc = tmp, Locus()
        self.probe_fa = os.path.join(tmp, "probes.fa")
        with open(self.probe_fa, "w") as fh:
            for gname, sq in self.loc.probes.items():
                fh.write(f">{gname}*01\n{sq}\n")
        self.root = os.path.join(tmp, "root")
        os.makedirs(os.path.join(self.root, "v9/wgs/long_read"), exist_ok=True)
        os.makedirs(os.path.join(self.root, "x"), exist_ok=True)
        self.sam, self.truth, coh, man = {}, {}, [], []
        for i in range(n):
            pid = str(3_000_000 + i)
            anc = "SAS" if rng.random() < 0.5 else "EUR"
            h = (rng.random() < INS_FREQ[anc], rng.random() < INS_FREQ[anc])
            self.truth[pid] = h
            bam = os.path.join(self.root, "x", f"{pid}.bam")
            if i != 0:  # person 0: BAM missing
                open(bam, "w").close()
            lines = []
            for k in (1, 2):
                lines += self.loc.reads(rng, h[k - 1], depth, rlen, k, f"{pid}h{k}")
            self.sam[bam] = "\n".join(lines) + "\n"
            man.append(dict(research_id=pid, platform="revio", grch38_bam=f"gs://vwb-aou-datasets-controlled/x/{pid}.bam"))
            coh.append(dict(person_id=pid, in_lr=True, ancestry_pred=anc.lower(), platform="revio"))
        self.cohort = os.path.join(tmp, "cohort.tsv")
        pd.DataFrame(coh).to_csv(self.cohort, sep="\t", index=False)
        pd.DataFrame(man).to_csv(os.path.join(self.root, tg.LR_MANIFEST_REL), sep="\t", index=False)

    def runner(self, bam, region, timeout, samtools):
        assert region == tg.DEFAULT_READS_REGION
        return self.sam[bam]

    def run(self, aligner, extra=()):
        out, internal = os.path.join(self.tmp, "out_r"), os.path.join(self.tmp, "int_r")
        orig = tg.run_samtools_view
        tg.run_samtools_view = self.runner
        try:
            tg.main(["--probes", self.probe_fa, "--cohort", self.cohort, "--aou-root", self.root,
                     "--out-dir", out, "--internal-dir", internal, "--work-dir", os.path.join(self.tmp, "work"),
                     "--aligner", aligner, "--workers", "4", "--breakpoint", self.loc.breakpoint, *extra])
        finally:
            tg.run_samtools_view = orig
        return out, internal


class TestInsertionBlock(unittest.TestCase):
    def test_separate_alt_segment_is_split_off(self):
        # GRCh38 KI270803v1 layout: TRBV3-2 / TRBV4-3 / 6-3 copy together, TRBV7-8 ~180 kb downstream
        loc = [dict(contig="alt", gene_class=g, start=s0, end=s0 + 290, called=True)
               for g, s0 in [("TRBV3-2", 321_274), ("TRBV4-3", 324_727), ("TRBV7-8", 515_798)]]
        block, away = tg.insertion_block(loc, "alt", ["TRBV4-3", "TRBV3-2", "TRBV7-8"])
        self.assertEqual(sorted(a["gene_class"] for a in block), ["TRBV3-2", "TRBV4-3"])
        self.assertEqual(away, {"TRBV7-8"})
        block, away = tg.insertion_block(loc[:2], "alt", ["TRBV4-3", "TRBV3-2"])
        self.assertEqual((len(block), away), (2, set()))


class TestReadsPrimitives(unittest.TestCase):
    def test_spans_contiguously(self):
        self.assertTrue(tg.spans_contiguously(100, "5000M", 1000, 3000, 50))
        self.assertFalse(tg.spans_contiguously(1500, "5000M", 1000, 3000, 50))       # starts inside
        self.assertFalse(tg.spans_contiguously(100, "2000M3000S", 1000, 3000, 50))   # clipped at the site
        self.assertFalse(tg.spans_contiguously(100, "1500M4000I3000M", 1000, 3000, 50))  # carries the insertion
        self.assertTrue(tg.spans_contiguously(100, "1500M10I3000M", 1000, 3000, 50))  # small indel ok
        self.assertFalse(tg.spans_contiguously(100, "1500M200D3000M", 1000, 3000, 50))
        self.assertTrue(tg.spans_contiguously(100, "500M200D3000M", 1000, 3000, 50))  # deletion outside window
        r = tg.parse_sam("r1\t0\tchr7\t10\t60\t5M\t*\t0\t0\tACGTA\t*\tHP:i:2\n@HD\tx\n")
        self.assertEqual((r[0]["name"], r[0]["hp"], r[0]["pos"]), ("r1", 2, 10))
        self.assertEqual([tg.call_copies(x, (0.2, 0.8)) for x in (0.0, 0.19, 0.2, 0.5, 0.79, 0.8, 1.0)],
                         [0, 0, 1, 1, 1, 2, 2])

    def test_breakpoint_found_at_probe_build(self):
        loc = Locus()
        tmp = tempfile.mkdtemp()
        pf = os.path.join(tmp, "user.fa")
        with open(pf, "w") as fh:
            for gname, sq in loc.probes.items():
                fh.write(f">{gname}*01\n{sq}\n")
        alt = os.path.join(tmp, "alt.fa")
        with open(alt, "w") as fh:
            fh.write(f">chr7_KI270803v1_alt\n{loc.alt}\n")
        pad = 5000
        chrom = os.path.join(tmp, "chr7.fa.gz")
        write_fa_gz(chrom, [("chr7", rs(random.Random(2), pad) + loc.primary + rs(random.Random(3), 2000))])
        region = f"chr7:{pad + 1}-{pad + len(loc.primary)}"
        for aligner in ["kmer"] + (["mappy"] if tg.mappy is not None else []):
            out = os.path.join(tmp, "out_" + aligner)
            tg.main(["--build-probes", "--probe-fasta", pf, "--alt-fasta", alt, "--ref-fasta", chrom,
                     "--ref-region", region, "--out-dir", out, "--aligner", aligner, "--work-dir", os.path.join(tmp, "w")])
            meta = tg.read_probe_meta(os.path.join(out, "probe_meta.tsv"))
            _c, lft, rgt = tg.parse_window(meta["breakpoint"])
            self.assertLessEqual(abs(lft - (pad + loc.B)), 300, (aligner, meta))
            self.assertLessEqual(abs(rgt - (pad + loc.B + 1)), 300, (aligner, meta))
            self.assertLessEqual(abs(int(meta["ins_len_estimate"]) - loc.I), 50, (aligner, meta))
            rc = pd.read_csv(os.path.join(out, "reference_check.tsv"), sep="\t")
            get = {(r.source, r.gene_class): r for r in rc.itertuples()}
            self.assertTrue(get[("alt", "TRBV4-3")].present)
            self.assertFalse(get[("primary", "TRBV4-3")].present)
            self.assertEqual(get[("primary", CLS62)].n_loci, 1)
        # the ALT carries a second, unrelated insertion far out in the left flank (as GRCh38's
        # KI270803v1 does): the breakpoint must still come from the block next to the insertion
        rng = random.Random(7)
        cut = 2000
        alt2 = loc.alt[:cut] + rs(rng, 4000) + loc.alt[cut:]
        with open(alt, "w") as fh:
            fh.write(f">chr7_KI270803v1_alt\n{alt2}\n")
        out = os.path.join(tmp, "out_extra_indel")
        tg.main(["--build-probes", "--probe-fasta", pf, "--alt-fasta", alt, "--ref-fasta", chrom,
                 "--ref-region", region, "--out-dir", out, "--aligner", "kmer", "--work-dir", os.path.join(tmp, "w2")])
        meta = tg.read_probe_meta(os.path.join(out, "probe_meta.tsv"))
        _c, lft, rgt = tg.parse_window(meta["breakpoint"])
        self.assertLessEqual(abs(lft - (pad + loc.B)), 300, meta)
        self.assertLessEqual(abs(int(meta["ins_len_estimate"]) - loc.I), 50, meta)
        # samtools-faidx style header gives the offset without --ref-region
        win = os.path.join(tmp, "win.fa")
        with open(win, "w") as fh:
            fh.write(f">chr7:142038120-142040000\nACGT\n")
        self.assertEqual(tg.load_ref(win)[:2], ("chr7", 142038119))


@unittest.skipIf(tg.mappy is None, "mappy not installed")
class TestReadsEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.w = ReadsWorld(cls.tmp, n=160)
        cls.out, cls.internal = cls.w.run("mappy", extra=("--no-cross-table",))   # see TestEndToEnd

    def test_genotypes_and_haplotypes(self):
        g = read_int(os.path.join(self.internal, "per_person_gene.tsv"))
        truth = pd.Series({p: int(h[0]) + int(h[1]) for p, h in self.w.truth.items()})
        for feat in ("TRBV4-3", "TRBV3-2", tg.INS_ANY):
            v = g[(g["feature"] == feat) & (g["status"] == "ok")].set_index("person_id")
            self.assertGreater(len(v), 150, feat)
            acc = (v["copies"].astype(int) == truth.loc[v.index]).mean()
            self.assertGreater(acc, 0.93, feat)
        self.assertFalse(any(f.startswith("TRBV6-2") for f in g["feature"].unique()))  # merged class not callable
        v = g[g["feature"] == "TRBV4-3"].set_index("person_id")
        hom_del = [p for p, h in self.w.truth.items() if not any(h) and p in v.index]
        self.assertTrue((v.loc[hom_del, "n_ins"].astype(int) == 0).all())  # the 95% decoy never counts
        # HP-based per-haplotype calls
        both = v[(v["hap1_complete"].astype(int) == 1) & (v["hap2_complete"].astype(int) == 1)]
        ok = [(int(r.hap1_present), int(r.hap2_present)) == tuple(int(x) for x in self.w.truth[p])
              for p, r in both.iterrows()]
        self.assertGreater(len(ok), 140)
        self.assertGreater(sum(ok) / len(ok), 0.9)
        self.assertTrue(v["method"].str.startswith("reads:mappy").all())

    def test_public_tables(self):
        tab = pd.read_csv(os.path.join(self.out, "genotype_by_stratum.tsv"), sep="\t")
        r = tab[(tab["feature"] == "TRBV4-3") & (tab["stratum"] == "ALL")].iloc[0]
        self.assertAlmostEqual(float(r["insertion_allele_freq"]), 0.525, delta=0.08)
        hist = pd.read_csv(os.path.join(self.out, "allele_fraction_histogram.tsv"), sep="\t", dtype=str)
        self.assertEqual(len(hist[(hist["feature"] == "TRBV4-3") & (hist["stratum"] == "ALL")]), 20)
        self.assertTrue(hist["n_people"].isin(["0", "<20", "<20*"] + [str(k) for k in range(20, 200)]).all())
        cp = pd.read_csv(os.path.join(self.out, "copresence.tsv"), sep="\t")
        allr = cp[cp["stratum"] == "ALL"]
        self.assertEqual(len(allr), 1)   # TRBV4-3 x TRBV3-2 (TRBV7-8 is not a reads-mode feature)
        self.assertTrue((allr["concordance"].astype(float) > 0.9).all())
        run = dict(pd.read_csv(os.path.join(self.out, "run_summary.tsv"), sep="\t").values)
        self.assertEqual(run["mode"], "reads")
        self.assertEqual(run["people_failed"], "<20")
        fs = pd.read_csv(os.path.join(self.out, "failures_summary.tsv"), sep="\t")
        self.assertEqual(list(fs["reason"]), ["bam_missing"])
        for f in os.listdir(self.out):
            txt = open(os.path.join(self.out, f)).read()
            for pid in ("3000000", "3000001", "3000159"):
                self.assertNotIn(pid, txt, f)

    def test_no_single_hidden_cell_in_any_sum(self):
        """Masked counts are not a subtraction away: per AF bin, the ancestries (and the platforms)
        sum to ALL; the bins of a stratum sum to its called people; the genotype counts sum to the
        called people. Every such set has 0 or >= 2 hidden cells."""
        def hidden(v):
            return not str(v).isdigit()

        hist = pd.read_csv(os.path.join(self.out, "allele_fraction_histogram.tsv"), sep="\t", dtype=str)
        self.assertTrue(hist["n_people"].isin(["<20", "<20*"]).any())
        anc = hist[hist["stratum_type"] == "ancestry"]
        plat = hist[(hist["stratum_type"] == "platform") | (hist["stratum"] == "ALL")]
        for part in (anc, plat):
            for _k, g in part.groupby(["feature", "af_lo"]):
                if len(g) > 1:
                    self.assertNotEqual(sum(hidden(v) for v in g["n_people"]), 1, _k)
        for _k, g in hist.groupby(["feature", "stratum_type", "stratum"]):
            self.assertNotEqual(sum(hidden(v) for v in g["n_people"]), 1, _k)
        tab = pd.read_csv(os.path.join(self.out, "genotype_by_stratum.tsv"), sep="\t", dtype=str)
        for _i, r in tab.iterrows():
            self.assertNotEqual(sum(hidden(r[c]) for c in ("n_called", "n_0copies", "n_1copy", "n_2copies")), 1)

    def test_cross_table_histogram_vs_genotypes(self):
        """In-script pass (no 66): the AF bins in [0, 0.2), [0.2, 0.8), [0.8, 1] sum to n_0copies /
        n_1copy / n_2copies of the same feature x stratum, and all bins to n_called; none of these sets
        may hold exactly one hidden cell (S05 fact-check 2026-10-12)."""
        hist = pd.read_csv(os.path.join(self.out, "allele_fraction_histogram.tsv"), sep="\t", dtype=str)
        tab = pd.read_csv(os.path.join(self.out, "genotype_by_stratum.tsv"), sep="\t", dtype=str)
        tab = tab.set_index(["feature", "stratum_type", "stratum"])
        n = 0
        for k, h in hist.groupby(["feature", "stratum_type", "stratum"]):
            lo = h["af_lo"].astype(float)
            for col, sel in (("n_called", lo >= 0), ("n_0copies", lo < 0.2), ("n_1copy", (lo >= 0.2) & (lo < 0.8)),
                             ("n_2copies", lo >= 0.8)):
                cells = [tab.at[k, col]] + list(h.loc[sel, "n_people"])
                self.assertNotEqual(sum(not str(x).isdigit() for x in cells), 1, (k, col))
                n += 1
        self.assertGreater(n, 8)

    def test_cross_table_pass(self):
        """Default PUBLIC run: 66 masks what the in-script pass cannot see (people not called, people
        with no listed ancestry) and its checker then finds no count of 1-19 recoverable."""
        out = os.path.join(self.tmp, "out_cross")
        self.w.run("mappy", extra=("--aggregate-only", "--out-dir", out))
        s66 = importlib.util.spec_from_file_location("cts66_65", os.path.join(SRC, "66_cross_table_suppression.py"))
        cts = importlib.util.module_from_spec(s66)
        s66.loader.exec_module(cts)
        P = cts.build(dir65=out)
        leaks, inconsistent = P.analyse()
        self.assertEqual(inconsistent, [])
        self.assertEqual(leaks, [])
        run = dict(pd.read_csv(os.path.join(out, "run_summary.tsv"), sep="\t").values)
        self.assertEqual(run["genotyper"], "cuts")

    def test_aggregate_only_new_floor(self):
        out = os.path.join(self.tmp, "out_floor")
        self.w.run("mappy", extra=("--aggregate-only", "--min-reads", "1000", "--out-dir", out))
        tab = pd.read_csv(os.path.join(out, "genotype_by_stratum.tsv"), sep="\t", dtype=str, keep_default_na=False)
        self.assertTrue((tab["n_called"] == "0").all())

    def test_mode_guard_and_resume(self):
        n0 = len(read_int(os.path.join(self.internal, "per_haplotype.tsv")))
        self.w.run("mappy", extra=("--no-cross-table",))  # rewrites the shared self.out
        self.assertEqual(n0, len(read_int(os.path.join(self.internal, "per_haplotype.tsv"))))
        with self.assertRaises(SystemExit):
            tg.main(["--mode", "assembly", "--probes", self.w.probe_fa, "--cohort", self.w.cohort,
                     "--aou-root", self.w.root, "--out-dir", os.path.join(self.tmp, "o2"),
                     "--internal-dir", self.internal, "--work-dir", os.path.join(self.tmp, "work"),
                     "--aligner", "kmer"])


class TestReadsKmerSmall(unittest.TestCase):
    def test_kmer_path(self):
        tmp = tempfile.mkdtemp()
        w = ReadsWorld(tmp, n=8, seed=4)
        out, internal = w.run("kmer")
        g = read_int(os.path.join(internal, "per_person_gene.tsv"))
        v = g[(g["feature"] == "TRBV4-3") & (g["status"] == "ok")].set_index("person_id")
        truth = pd.Series({p: int(h[0]) + int(h[1]) for p, h in w.truth.items()})
        self.assertGreaterEqual(len(v), 6)
        self.assertGreaterEqual((v["copies"].astype(int) == truth.loc[v.index]).mean(), 0.85)
        tab = pd.read_csv(os.path.join(out, "genotype_by_stratum.tsv"), sep="\t", dtype=str, keep_default_na=False)
        self.assertTrue((tab["n_people_measured"] == "<20").all())
        self.assertTrue((tab["insertion_allele_freq"] == "").all())

    def test_no_breakpoint_refuses(self):
        tmp = tempfile.mkdtemp()
        w = ReadsWorld(tmp, n=2)
        with self.assertRaises(SystemExit):
            tg.main(["--probes", w.probe_fa, "--cohort", w.cohort, "--aou-root", w.root, "--aligner", "kmer",
                     "--out-dir", os.path.join(tmp, "o"), "--internal-dir", os.path.join(tmp, "i"),
                     "--work-dir", os.path.join(tmp, "w")])


class TestStratifiedAndCompare(unittest.TestCase):
    def test_sample_stratified_then_assembly_subset(self):
        tmp = tempfile.mkdtemp()
        w = World(tmp, n=60, seed=3)
        rel = os.path.join(tmp, "rel.tsv")
        sas = sorted(p for p, a in w.anc.items() if a == "SAS")
        pd.DataFrame([dict(**{"i.s": sas[0], "j.s": sas[1], "kin": 0.25}),
                      dict(**{"i.s": sas[2], "j.s": sas[3], "kin": 0.5})]).to_csv(rel, sep="\t", index=False)
        pf = os.path.join(tmp, "internal_x", "strat.tsv")
        args = ["--cohort", w.cohort, "--aou-root", w.root, "--relatedness", rel, "--aligner", "kmer",
                "--work-dir", os.path.join(tmp, "w"), "--people-file", pf, "--seed", "1"]
        tg.main(["--sample-stratified", "20"] + args)
        sel = pd.read_csv(pf, sep="\t", comment="#", dtype={"person_id": str})
        self.assertEqual(len(sel), 20)
        self.assertEqual(sel["ancestry"].value_counts().to_dict(), {"SAS": 10, "EUR": 10})
        self.assertFalse({sas[0], sas[1]} <= set(sel["person_id"]))
        self.assertFalse({sas[2], sas[3]} <= set(sel["person_id"]))
        tg.main(["--sample-stratified", "20"] + args[:-1] + ["2"])
        sel2 = pd.read_csv(pf, sep="\t", comment="#", dtype={"person_id": str})
        self.assertNotEqual(set(sel["person_id"]), set(sel2["person_id"]))  # seed matters
        # water-filling: ask for more than the smaller group holds
        tg.main(["--sample-stratified", "50"] + args)
        sel3 = pd.read_csv(pf, sep="\t", comment="#", dtype={"person_id": str})
        counts = sel3["ancestry"].value_counts().to_dict()
        self.assertEqual(sum(counts.values()), 50)
        self.assertTrue(counts["SAS"] == len(sas) - 2 or counts["EUR"] == 60 - len(sas))
        out, internal = w.run("kmer", extra=("--people-file", pf), out="o_sub", internal="i_sub")
        h = read_int(os.path.join(internal, "per_haplotype.tsv"))
        f = tg._read(os.path.join(internal, "failures.tsv"), tg.FAIL_COLS)
        self.assertEqual(set(h["person_id"]) | set(f["person_id"]), set(sel3["person_id"]))

    def test_compare(self):
        tmp = tempfile.mkdtemp()
        rng = random.Random(8)
        ids = [str(4_000_000 + i) for i in range(80)]
        coh = pd.DataFrame(dict(person_id=ids, ancestry_pred=["sas" if i % 2 else "eur" for i in range(80)],
                                platform="revio", in_lr=True))
        cpath = os.path.join(tmp, "cohort.tsv")
        coh.to_csv(cpath, sep="\t", index=False)
        cop = {p: rng.choice([0, 1, 2]) for p in ids}
        rows_r, rows_a = [], []
        for i, p in enumerate(ids):
            ca = cop[p] if i >= 6 else (cop[p] + 1) % 3  # 6 discordant
            rows_r.append(tg._row(tg.GENE_COLS, person_id=p, feature="TRBV4-3", copies=cop[p], status="ok"))
            rows_a.append(tg._row(tg.GENE_COLS, person_id=p, feature="TRBV4-3", copies=ca,
                                  status="ok" if i < 75 else "incomplete"))
        rp, apath = os.path.join(tmp, "r.tsv"), os.path.join(tmp, "a.tsv")
        pd.DataFrame(rows_r).to_csv(rp, sep="\t", index=False)
        pd.DataFrame(rows_a).to_csv(apath, sep="\t", index=False)
        out = os.path.join(tmp, "out")
        tg.main(["--compare", rp, apath, "--cohort", cpath, "--out-dir", out])
        t = pd.read_csv(os.path.join(out, "65b_method_concordance.tsv"), sep="\t", dtype=str, keep_default_na=False)
        allr = t[t["stratum"] == "ALL"].iloc[0]
        self.assertEqual(allr["n_both_called"], "75")
        self.assertEqual(allr["n_concordant"], "69")
        self.assertEqual(allr["pct_concordant"], "92.0")
        offdiag = [allr[f"n_reads{i}_asm{j}"] for i in range(3) for j in range(3) if i != j]
        self.assertTrue(all(x in ("0", "<20") for x in offdiag))
        self.assertFalse(any(p in open(os.path.join(out, "65b_method_concordance.tsv")).read() for p in ids))

    def test_compare_rederives_reads_copies(self):
        tmp = tempfile.mkdtemp()
        ids = [str(5_000_000 + i) for i in range(60)]
        coh = pd.DataFrame(dict(person_id=ids, ancestry_pred="eur", platform="revio", in_lr=True))
        cpath = os.path.join(tmp, "cohort.tsv")
        coh.to_csv(cpath, sep="\t", index=False)
        # 20 people at each of af 0.0, 0.85 (stored as 2 copies under cuts 0.2,0.8) and 1.0; assembly
        # says 0 / 1 / 2. Stored and re-derived cuts agree with each other; cuts 0.2,0.9 move the 0.85s to 1.
        rows_r, rows_a = [], []
        for i, p in enumerate(ids):
            k = i // 20
            ni = (0, 17, 20)[k]
            rows_r.append(tg._row(tg.GENE_COLS, person_id=p, feature="TRBV4-3", copies=(0, 2, 2)[k],
                                  status="ok", n_ins=ni, n_ref=20 - ni))
            rows_a.append(tg._row(tg.GENE_COLS, person_id=p, feature="TRBV4-3", copies=k, status="ok"))
        rp, apath = os.path.join(tmp, "r.tsv"), os.path.join(tmp, "a.tsv")
        pd.DataFrame(rows_r).to_csv(rp, sep="\t", index=False)
        pd.DataFrame(rows_a).to_csv(apath, sep="\t", index=False)

        def conc(*extra):
            out = os.path.join(tmp, "o" + "_".join(extra).replace(",", "").replace("-", ""))
            tg.main(["--compare", rp, apath, "--cohort", cpath, "--out-dir", out] + list(extra))
            t = pd.read_csv(os.path.join(out, "65b_method_concordance.tsv"), sep="\t", dtype=str,
                            keep_default_na=False)
            return t[t["stratum"] == "ALL"].iloc[0]["n_concordant"]
        self.assertEqual(conc(), "40")
        self.assertEqual(conc("--compare-genotyper", "cuts"), "40")
        self.assertEqual(conc("--compare-genotyper", "cuts", "--af-cuts", "0.2,0.9"), "60")


if __name__ == "__main__":
    unittest.main()


class TestAfMixture(unittest.TestCase):
    def test_recovers_shifted_het_fraction_and_calls(self):
        import numpy as np
        rng = np.random.default_rng(65)
        truth = rng.choice(3, size=3000, p=[0.33, 0.48, 0.19])
        n = rng.integers(6, 25, size=truth.size)
        f_true = np.array([0.01, 0.64, 0.99])[truth]
        k = rng.binomial(n, f_true)
        f, w, post = tg.fit_af_mixture(k, n)
        self.assertAlmostEqual(f[1], 0.64, delta=0.02)
        self.assertAlmostEqual(w[2], 0.19, delta=0.03)
        call = post.argmax(axis=1)
        sure = post.max(axis=1) >= 0.9
        self.assertGreater((call[sure] == truth[sure]).mean(), 0.97)
        # fixed 0.2/0.8 cuts on the same data overcall 2 copies, the S05 heterozygote deficit
        af = k / n
        cut2 = (af >= 0.8).mean()
        self.assertGreater(cut2, 0.19 + 0.03)

    def test_mixture_reads_calls_per_platform(self):
        import numpy as np
        rng = np.random.default_rng(1)
        rows, people = [], []
        for i in range(600):
            pl = "revio" if i % 3 else "sequel2e"
            g = rng.choice(3, p=[0.3, 0.5, 0.2])
            n = int(rng.integers(8, 20))
            k = int(rng.binomial(n, [0.01, 0.62 if pl == "revio" else 0.7, 0.99][g]))
            rows.append(dict(person_id=f"p{i}", feature="TRBV4-3", n_ins=k, n_ref=n - k))
            people.append(dict(person_id=f"p{i}", ancestry="EUR", platform=pl))
        rows.append(dict(person_id="low", feature="TRBV4-3", n_ins=1, n_ref=1))
        people.append(dict(person_id="low", ancestry="EUR", platform="revio"))
        calls, fit = tg.mixture_reads_calls(pd.DataFrame(rows), pd.DataFrame(people), min_reads=6)
        self.assertEqual(set(fit["platform"]), {"revio", "sequel2e"})
        self.assertEqual(calls.set_index("person_id").loc["low", "status"], "low_reads")
        self.assertTrue(set(calls["status"]) <= {"ok", "ambiguous", "low_reads"})
        rev = fit.set_index("platform").loc["revio", "f_1copy"]
        self.assertAlmostEqual(rev, 0.62, delta=0.04)
