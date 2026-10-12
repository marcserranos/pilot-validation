#!/usr/bin/env python3
"""Tests for the sequel2 (Tier 3) fix in scripts/run_immuannot_person.py and for the
63_sequel2_recall.py driver. Synthetic data only: no AoU data, no minimap2/samtools/immuannot.sh
(each is replaced by a fake that writes what the real tool would).

The sequel2 case reproduced: a haplotype whose assembly FASTA is PLAIN gzip (not BGZF) and that has
no assembly-to-hg38 PAF or BAM (paf_rel = aln_rel = None). Its contigs:
  ctgMHC   2.0 Mb, the true MHC contig. Whole-genome alignment: one 250 kb primary block on
           chr6:30.0-30.25 Mb at q 100-350 kb, plus a stray 12 kb primary block that also lands in
           the HLA window at q 1.80-1.812 Mb (the "balloon" case: min/max would span ~1.9 Mb).
  ctgSD    300 kb, a chr1 contig sharing a 12 kb segmental duplication with the MHC. Against a
           chr6-only reference that block is primary; against the whole genome the contig is
           primary on chr1 and the chr6 hit is secondary with MAPQ 0.
  ctgREP   80 kb, a repeat-rich contig: only 2 kb primary blocks (plus a 15 kb secondary).
  ctgOTHER 50 kb, no hit.

Run: python3 -m pytest scripts/hla_popgen/tests/test_63_sequel2_recall.py -q
"""
import argparse
import gzip
import hashlib
import importlib.util
import os
import random
import subprocess
import sys
import tempfile
import unittest
import zlib
import struct

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
SCRIPTS = os.path.dirname(SRC)
sys.path.insert(0, SRC)
sys.path.insert(0, SCRIPTS)
import run_immuannot_person as rip  # noqa: E402

spec = importlib.util.spec_from_file_location("sr63", os.path.join(SRC, "63_sequel2_recall.py"))
sr = importlib.util.module_from_spec(spec)
sys.modules["sr63"] = sr
spec.loader.exec_module(sr)

HLA_REGION = "chr6:29500000-33500000"
SLICE_NAME = "chr6:28500001-34500000"


def rseq(rng, n):
    return "".join(rng.choices("ACGT", k=n))


def write_fasta(path, records, gz=True, width=80):
    opener = gzip.open if gz else open
    with opener(path, "wt") as f:
        for name, seq in records:
            f.write(f">{name} some description\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + "\n")


def write_bgzf(path, text):
    """One BGZF block (gzip member with the 'BC' extra subfield) + the BGZF EOF block."""
    data = text.encode()
    co = zlib.compressobj(9, zlib.DEFLATED, -15)
    cdata = co.compress(data) + co.flush()
    bsize = 18 + len(cdata) + 8 - 1
    hdr = b"\x1f\x8b\x08\x04" + b"\x00" * 4 + b"\x00\xff" + struct.pack("<H", 6) + b"BC" + \
        struct.pack("<H", 2) + struct.pack("<H", bsize)
    tail = struct.pack("<I", zlib.crc32(data) & 0xffffffff) + struct.pack("<I", len(data))
    eof = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
    with open(path, "wb") as f:
        f.write(hdr + cdata + tail + eof)


def paf(q, qlen, qs, qe, t, tlen, ts, te, mapq, tp="P", nmatch=None):
    blen = qe - qs
    nmatch = blen if nmatch is None else nmatch
    return "\t".join(map(str, [q, qlen, qs, qe, "+", t, tlen, ts, te, nmatch, blen, mapq,
                               f"tp:A:{tp}"])) + "\n"


class World:
    def __init__(self, root):
        rng = random.Random(7)
        self.seqs = {"ctgMHC": rseq(rng, 2_000_000), "ctgSD": rseq(rng, 300_000),
                     "ctgREP": rseq(rng, 80_000), "ctgOTHER": rseq(rng, 50_000)}
        self.root = root
        self.mount = os.path.join(root, "mount")
        os.makedirs(os.path.join(self.mount, "asm"), exist_ok=True)
        self.fa_rel = "asm/p1.hap1.fa.gz"
        write_fasta(os.path.join(self.mount, self.fa_rel), list(self.seqs.items()))
        self.slice = os.path.join(root, "hla_stage1.fa")
        with open(self.slice, "w") as f:
            f.write(f">{SLICE_NAME}\nACGT\n")
        self.kir_slice = os.path.join(root, "kir_stage1.fa")
        with open(self.kir_slice, "w") as f:
            f.write(">chr19:53600001-55920000\nACGT\n")
        self.mmi = os.path.join(root, "hg38_primary.asm5.mmi")
        open(self.mmi, "w").write("fake index")
        self.calls = []

    def stage1_paf(self):
        L = {k: len(v) for k, v in self.seqs.items()}
        return "".join([
            paf("ctgMHC", L["ctgMHC"], 100_000, 350_000, SLICE_NAME, 6_000_000, 1_500_000, 1_750_000, 60),
            paf("ctgSD", L["ctgSD"], 150_000, 162_000, SLICE_NAME, 6_000_000, 2_000_000, 2_012_000, 60),
            paf("ctgREP", L["ctgREP"], 10_000, 12_000, SLICE_NAME, 6_000_000, 3_000_000, 3_002_000, 60),
            paf("ctgREP", L["ctgREP"], 30_000, 45_000, SLICE_NAME, 6_000_000, 3_100_000, 3_115_000, 0, "S"),
        ])

    def stage2_paf(self):
        L = {k: len(v) for k, v in self.seqs.items()}
        return "".join([
            paf("ctgMHC", L["ctgMHC"], 100_000, 350_000, "chr6", 170_805_979, 30_000_000, 30_250_000, 60),
            paf("ctgMHC", L["ctgMHC"], 1_800_000, 1_812_000, "chr6", 170_805_979, 32_000_000, 32_012_000, 60),
            paf("ctgSD", L["ctgSD"], 0, 300_000, "chr1", 248_956_422, 10_000_000, 10_300_000, 60),
            paf("ctgSD", L["ctgSD"], 150_000, 162_000, "chr6", 170_805_979, 31_000_000, 31_012_000, 0, "S"),
        ])

    def fake_self_align_paf(self, fa_path, ref, out_paf, threads, minimap2="minimap2"):
        self.calls.append((os.path.basename(fa_path), os.path.basename(ref)))
        if ref.endswith(".mmi"):
            with open(fa_path) as f:
                self.stage2_query_names = [l[1:].strip() for l in f if l.startswith(">")]
            text = self.stage2_paf()
        else:
            text = self.stage1_paf()
        with open(out_paf, "w") as f:
            f.write(text)
        return True, ""

    def args(self, **kw):
        a = dict(mount=self.mount, region=HLA_REGION, pad=100_000, enable_self_align_fallback=True,
                 self_align_ref=self.slice, chr6_ref_cache=None,
                 self_align_ref_dir=os.path.join(self.root, "ref"), self_align_flank=1_000_000,
                 self_align_confirm_ref=self.mmi, hg38_ref=os.path.join(self.root, "nohg38.fa"),
                 self_align_min_block=10_000, self_align_min_mapq=20, self_align_min_ident=0.5,
                 self_align_stage1_min_ident=0.8,
                 self_align_cluster_gap=250_000, self_align_max_mb=25.0,
                 self_align_max_candidate_mb=600.0, minimap2="minimap2", refdir="refdir", threads=2)
        a.update(kw)
        return argparse.Namespace(**a)


class TestTier3Pieces(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.w = World(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_plain_gzip_is_not_bgzf_but_bgzf_is(self):
        self.assertFalse(rip.is_bgzf(os.path.join(self.w.mount, self.w.fa_rel)))
        p = os.path.join(self.tmp.name, "b.fa.gz")
        write_bgzf(p, ">c1\nACGTACGTAC\n>c2\nGGGG\n")
        self.assertTrue(rip.is_bgzf(p))
        out = os.path.join(self.tmp.name, "o.fa")
        ok, _ = rip.stream_extract_fasta(p, ["c2", "c1:3-6"], out)
        self.assertTrue(ok)
        self.assertEqual(open(out).read(), ">c2\nGGGG\n>c1:3-6\nGTAC\n")

    def test_stream_extract_plain_gzip_matches_samtools_format(self):
        out = os.path.join(self.tmp.name, "o.fa")
        src = os.path.join(self.w.mount, self.w.fa_rel)
        ok, msg = rip.stream_extract_fasta(src, ["ctgREP:101-230", "ctgOTHER", "ctgREP:79990-90000"], out)
        self.assertTrue(ok, msg)
        lines = open(out).read().splitlines()
        self.assertEqual(lines[0], ">ctgREP:101-230")
        self.assertEqual(lines[1], self.w.seqs["ctgREP"][100:160])          # 60-col wrap
        self.assertEqual("".join(lines[1:4]), self.w.seqs["ctgREP"][100:230])
        i = lines.index(">ctgOTHER")
        j = lines.index(">ctgREP:79990-90000")
        self.assertEqual("".join(lines[i + 1:j]), self.w.seqs["ctgOTHER"])
        self.assertEqual("".join(lines[j + 1:]), self.w.seqs["ctgREP"][79989:])  # clamped at end
        ok, msg = rip.stream_extract_fasta(src, ["nope"], out)
        self.assertFalse(ok)
        self.assertIn("not found", msg)

    def test_regions_from_paf_default_is_unchanged_and_shows_the_old_bug(self):
        p = os.path.join(self.tmp.name, "s2.paf")
        open(p, "w").write(self.w.stage2_paf())
        regions, qlens, _ = rip.regions_from_paf(p, "chr6", 29_500_000, 33_500_000, 100_000)
        # Old (Tier 1) semantics: min(qstart)-pad .. max(qend)+pad over every overlapping block, any
        # tp/mapq -> ctgMHC balloons to ~1.9 Mb and the secondary ctgSD hit is kept.
        self.assertEqual(regions["ctgMHC"], (0, 1_912_000))
        self.assertIn("ctgSD", regions)

    def test_regions_from_paf_filters_and_clustering(self):
        p = os.path.join(self.tmp.name, "s2.paf")
        open(p, "w").write(self.w.stage2_paf())
        regions, _, _ = rip.regions_from_paf(p, "chr6", 29_500_000, 33_500_000, 100_000,
                                             min_block=10_000, min_mapq=20, min_ident=0.5,
                                             primary_only=True, cluster_gap=250_000)
        self.assertEqual(set(regions), {"ctgMHC"})
        self.assertEqual(regions["ctgMHC"], (0, 450_000))

    def test_candidates_from_stage1(self):
        p = os.path.join(self.tmp.name, "s1.paf")
        open(p, "w").write(self.w.stage1_paf())
        c = rip.candidate_contigs_from_paf(p, 10_000)
        self.assertEqual(set(c), {"ctgMHC", "ctgSD"})   # ctgREP: short primary, long secondary only

    def test_stage1_drops_low_identity_repeat_chains(self):
        # 63 smoke: long primary blocks at col10/col11 < 0.8 pulled in the whole assembly.
        p = os.path.join(self.tmp.name, "s1_junk.paf")
        open(p, "w").write(self.w.stage1_paf() + "".join(
            paf(f"ctgJUNK{i}", 500_000, 0, 40_000, SLICE_NAME, 6_000_000, 100_000 * i,
                100_000 * i + 40_000, 0, nmatch=12_000) for i in range(5)))
        self.assertEqual(len(rip.candidate_contigs_from_paf(p, 10_000)), 7)
        c = rip.candidate_contigs_from_paf(p, 10_000, min_ident=0.8)
        self.assertEqual(set(c), {"ctgMHC", "ctgSD"})

    def test_ensure_ref_slice_uses_region_chrom_and_checks_cache(self):
        hg = os.path.join(self.tmp.name, "hg38.fa")
        open(hg, "w").write(">chr19\nACGT\n")
        seen = []

        def fake_run(cmd, stdout=None, stderr=None, text=None):
            seen.append(cmd)
            stdout.write(f">{cmd[-1]}\nACGT\n")
            return subprocess.CompletedProcess(cmd, 0, "", "")

        orig = rip.subprocess.run
        rip.subprocess.run = fake_run
        try:
            cache = os.path.join(self.tmp.name, "ref", "chr19.fa")
            self.assertEqual(rip.ensure_ref_slice(hg, cache, "chr19"), cache)
            self.assertEqual(seen[0][-1], "chr19")
            self.assertEqual([f for f in os.listdir(os.path.dirname(cache)) if ".tmp" in f], [])
            # a chr19 cache must never be served for a chr6 request (and vice versa: the KIR bug)
            self.assertIsNone(rip.ensure_ref_slice(hg, cache, "chr6"))
            self.assertEqual(rip.ensure_chr6_ref(hg, os.path.join(self.tmp.name, "ref", "c6.fa")),
                             os.path.join(self.tmp.name, "ref", "c6.fa"))
            self.assertEqual(seen[-1][-1], "chr6")
        finally:
            rip.subprocess.run = orig


class TestSequel2Haplotype(unittest.TestCase):
    """process_haplotype() end to end on the synthetic sequel2 haplotype."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.w = World(self.tmp.name)
        self.pdir = os.path.join(self.tmp.name, "out", "p1", "immuannot_output")
        os.makedirs(self.pdir)
        self.orig = (rip.self_align_paf, rip.run_immuannot)
        rip.self_align_paf = self.w.fake_self_align_paf
        self.immu_inputs = []

        def fake_immuannot(script, refdir, contig_path, outprefix, threads):
            self.immu_inputs.append(open(contig_path).read())
            gtf = outprefix + ".gtf.gz"
            with gzip.open(gtf, "wt") as f:
                f.write('ctgMHC\timmuannot\ttranscript\t1\t9\t.\t+\t.\tgene_name "HLA-A"; '
                        'consensus "HLA-A*01:01:01:01";\n')
            return gtf
        rip.run_immuannot = fake_immuannot

    def tearDown(self):
        rip.self_align_paf, rip.run_immuannot = self.orig
        self.tmp.cleanup()

    def run_hap(self, **kw):
        return rip.process_haplotype("p1", "hap1", self.w.fa_rel, None, None, self.pdir, "immuannot.sh",
                                     self.w.args(**kw), 0.0, 2)

    def test_sequel2_haplotype_is_called_from_a_small_confirmed_trim(self):
        row, calls = self.run_hap()
        self.assertEqual(calls, {"HLA-A": "HLA-A*01:01:01:01"})
        self.assertEqual(row["trim_method"], "self_align_confirmed")
        self.assertEqual(row["self_align_candidates"], 2)
        self.assertEqual(sorted(self.w.stage2_query_names), ["ctgMHC", "ctgSD"])
        self.assertEqual(row["n_contigs"], 1)
        self.assertAlmostEqual(row["padded_mb"], 0.45)
        trimmed = self.immu_inputs[0].splitlines()
        self.assertEqual(trimmed[0], ">ctgMHC:1-450000")
        self.assertEqual("".join(trimmed[1:]), self.w.seqs["ctgMHC"][:450_000])
        self.assertFalse(os.path.exists(os.path.join(self.pdir, "hap1.self_align_candidates.fa")))
        self.assertTrue(os.path.exists(os.path.join(self.pdir, "hap1.self_align_confirm.paf")))

    def test_without_confirm_ref_tier3_refuses(self):
        row, calls = self.run_hap(self_align_confirm_ref=None)
        self.assertEqual(calls, {})
        self.assertEqual(self.w.calls, [])

    def test_kir_region_with_chr6_reference_is_refused(self):
        row, calls = self.run_hap(region="chr19:54600000-54920000")   # self_align_ref = chr6 slice
        self.assertEqual(calls, {})
        self.assertEqual(self.w.calls, [])
        row, calls = self.run_hap(region="chr19:54600000-54920000", self_align_ref=self.w.kir_slice)
        self.assertEqual(self.w.calls[0][1], "kir_stage1.fa")         # chr19 slice accepted

    def test_size_guard(self):
        row, calls = self.run_hap(self_align_max_mb=0.1)
        self.assertEqual(calls, {})
        self.assertEqual(self.immu_inputs, [])
        self.assertFalse(os.path.exists(os.path.join(self.pdir, "hap1.self_align_candidates.fa")))

    def test_fallback_flag_still_required(self):
        row, calls = self.run_hap(enable_self_align_fallback=False)
        self.assertEqual(calls, {})
        self.assertEqual(self.w.calls, [])


# ---------------------------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------------------------
def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


class TestDriver(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        t = self.tmp.name
        self.internal = os.path.join(t, "s05", "internal", "63")
        self.out = os.path.join(t, "s05", "results", "63")
        self.mount = os.path.join(t, "mount")
        os.makedirs(os.path.join(self.mount, "v9", "wgs", "long_read"))
        open(os.path.join(self.mount, "v9", "wgs", "long_read", "manifest.tsv"), "w").write("research_id\n")
        self.cohort = os.path.join(t, "cohort.tsv")
        pd.DataFrame({"person_id": ["a1", "s1", "s2", "a2", "s3"],
                      "platform": ["revio", "sequel2", "sequel2", "sequel2e", "sequel2"],
                      "trim_tier": ["paf_region", "self_align_needed", "self_align_needed",
                                    "paf_region", "self_align_needed"],
                      "n_rows": "1", "ancestry_pred": "afr"}).to_csv(self.cohort, sep="\t", index=False)
        # existing production tables (read-only inputs)
        self.ex_hla = os.path.join(t, "prod", "immuannot_calls.tsv")
        self.ex_kir = os.path.join(t, "prod_kir", "immuannot_calls.kir.tsv")
        for p, g in ((self.ex_hla, "HLA-A"), (self.ex_kir, "KIR2DL1")):
            os.makedirs(os.path.dirname(p))
            pd.DataFrame({"person_id": ["a1", "a2"], "gene": g, "immuannot_1": "x",
                          "immuannot_2": "y"}).to_csv(p, sep="\t", index=False)
        self.md5s = {p: md5(p) for p in (self.ex_hla, self.ex_kir)}
        self.saved = {k: dict(v) for k, v in sr.LOCI.items()}
        sr.LOCI["hla"]["existing"] = self.ex_hla
        sr.LOCI["kir"]["existing"] = self.ex_kir
        refs = sr.ref_paths(self.internal)
        os.makedirs(refs["dir"])
        for k in ("hla", "kir", "mmi"):
            open(refs[k], "w").write(">x\nA\n")
        self.cmds = []
        self.fail = {("s2", "kir")}

    def tearDown(self):
        for k, v in self.saved.items():
            sr.LOCI[k] = v
        self.tmp.cleanup()

    def runner(self, cmd, **kw):
        self.cmds.append(cmd)
        pid = cmd[2]
        outroot = cmd[cmd.index("--outroot") + 1]
        region = cmd[cmd.index("--region") + 1]
        sfx = cmd[cmd.index("--out-suffix") + 1]
        locus = "kir" if region.startswith("chr19") else "hla"
        assert "--enable-self-align-fallback" in cmd and "--self-align-confirm-ref" in cmd
        assert cmd[cmd.index("--self-align-ref") + 1].endswith(f"{locus}_stage1.fa")
        if (pid, locus) in self.fail:
            return subprocess.CompletedProcess(cmd, 1, "", "  hap1: SKIP -- no candidate contig has a "
                                               "primary, MAPQ>=20 block\n")
        d = os.path.join(outroot, pid, "immuannot_output")
        os.makedirs(d, exist_ok=True)
        gene = "KIR2DL1" if locus == "kir" else "HLA-A"
        for h in ("hap1", "hap2"):
            with gzip.open(os.path.join(d, f"{h}.gtf.gz"), "wt") as f:
                f.write(f'c\tI\ttranscript\t1\t2\t.\t+\t.\tgene_name "{gene}"; consensus "{gene}*01";\n')
        pd.DataFrame({"person_id": [pid], "gene": [gene], "immuannot_1": ["n1"],
                      "immuannot_2": ["n2"]}).to_csv(os.path.join(outroot, f"immuannot_calls{sfx}.tsv"),
                                                     sep="\t", index=False)
        pd.DataFrame({"person_id": [pid, pid], "hap": ["hap1", "hap2"],
                      "trim_method": "self_align_confirmed", "trimmed_mb": ["4.3", "4.1"],
                      "self_align_seconds": "60", "confirm_seconds": "90",
                      "immuannot_seconds": "170"}).to_csv(
            os.path.join(outroot, f"immuannot_timing{sfx}.tsv"), sep="\t", index=False)
        return subprocess.CompletedProcess(cmd, 0, "", "")

    def main(self, *extra):
        return sr.main(["--cohort", self.cohort, "--internal-dir", self.internal, "--out-dir", self.out,
                        "--mount", self.mount, "--threads", "2", "--python", sys.executable] + list(extra),
                       runner=self.runner)

    def test_assert_safe(self):
        home = self.tmp.name
        os.makedirs(os.path.join(home, "pipeline_outputs_kir"))
        for bad in ("pipeline_outputs", "pipeline_outputs/x", "pipeline_outputs_kir/y", "pipeline_outputs_new"):
            with self.assertRaises(SystemExit):
                sr.assert_safe(os.path.join(home, bad), home=home)
        sr.assert_safe(os.path.join(home, "s05", "internal", "63"), home=home)

    def test_full_cycle_resume_merge_and_public(self):
        self.main()
        # only the 3 self_align_needed people, both loci, HLA before KIR per person
        pairs = [(c[2], "kir" if c[c.index("--region") + 1].startswith("chr19") else "hla") for c in self.cmds]
        self.assertEqual(sorted(pairs), sorted([(p, l) for p in ("s1", "s2", "s3") for l in ("hla", "kir")]))
        self.assertTrue(all(c[c.index("--out-suffix") + 1] in (f".{c[2]}", f".kir.{c[2]}") for c in self.cmds))
        # existing tables untouched; combined tables written to new paths
        for p, h in self.md5s.items():
            self.assertEqual(md5(p), h)
        m_hla = pd.read_csv(os.path.join(self.internal, "merged", "immuannot_calls.with_sequel2.tsv"), sep="\t")
        m_kir = pd.read_csv(os.path.join(self.internal, "merged", "immuannot_calls.kir.with_sequel2.tsv"), sep="\t")
        self.assertEqual(sorted(m_hla["person_id"]), ["a1", "a2", "s1", "s2", "s3"])
        self.assertEqual(sorted(m_kir["person_id"]), ["a1", "a2", "s1", "s3"])
        rs = pd.read_csv(os.path.join(self.out, "run_summary.tsv"), sep="\t", dtype=str).set_index("locus")
        self.assertEqual(rs.loc["kir", "people_done"], "<20")      # small counts masked
        fr = pd.read_csv(os.path.join(self.out, "failure_reasons.tsv"), sep="\t", dtype=str)
        self.assertEqual(fr["reason"].tolist(), ["stage2_no_confirmed_contig"])
        ts = pd.read_csv(os.path.join(self.out, "trim_summary.tsv"), sep="\t", dtype=str)
        self.assertEqual(set(ts["trim_method"]), {"self_align_confirmed"})
        self.assertFalse(any("s1" in open(os.path.join(self.out, f)).read() for f in os.listdir(self.out)))
        # resume: only the failed person-locus is retried (attempt 2), then it is given up
        self.cmds.clear()
        self.main()
        self.assertEqual([(c[2], c[c.index("--region") + 1][:5]) for c in self.cmds], [("s2", "chr19")])
        self.cmds.clear()
        self.main()
        self.assertEqual(self.cmds, [])
        st = pd.read_csv(os.path.join(self.internal, "kir", "status.tsv"), sep="\t", dtype=str)
        self.assertEqual(st[st["person_id"] == "s2"]["attempt"].tolist(), ["1", "2"])
        self.assertFalse(os.path.exists(os.path.join(self.internal, "63.lock")))

    def test_limit_and_smoke(self):
        self.main("--limit", "1", "--loci", "hla")
        self.assertEqual([c[2] for c in self.cmds], ["s1"])
        self.cmds.clear()
        self.main("--smoke")
        self.assertEqual(sorted({c[2] for c in self.cmds}), ["s1", "s2"])
        self.assertTrue(all(self.internal + "_smoke" in c[c.index("--outroot") + 1] for c in self.cmds))
        self.assertTrue(os.path.exists(os.path.join(self.out + "_smoke", "run_summary.tsv")))

    def test_lock_refuses_second_instance(self):
        os.makedirs(self.internal, exist_ok=True)
        open(os.path.join(self.internal, "63.lock"), "w").write(str(os.getpid()))
        with self.assertRaises(SystemExit):
            self.main()
        os.remove(os.path.join(self.internal, "63.lock"))

    def test_resolve_tool_prefers_env_bin(self):
        envbin = os.path.join(self.tmp.name, "envbin")
        os.makedirs(envbin)
        py = os.path.join(envbin, "python3")
        self.assertEqual(sr.resolve_tool("minimap2", py), "minimap2")   # not there: left for PATH
        mm2 = os.path.join(envbin, "minimap2")
        open(mm2, "w").write("#!/bin/sh\n")
        os.chmod(mm2, 0o755)
        self.assertEqual(sr.resolve_tool("minimap2", py), mm2)
        self.assertEqual(sr.resolve_tool("/opt/x/minimap2", py), "/opt/x/minimap2")  # explicit path kept

    def test_refs_required(self):
        os.remove(sr.ref_paths(self.internal)["mmi"])
        with self.assertRaises(SystemExit):
            self.main()


if __name__ == "__main__":
    unittest.main()
