"""Tests for the 90 sr step's read logic (synthetic reads; no samtools/minimap2 needed)."""
import importlib.util
import os
import random
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


cb = load("cb90", "90_case_browser.py")
sr = load("cb90sr", "90_case_browser_sr.py")
RNG = random.Random(4)
REF = "ATG" + "".join(RNG.choice("ACGT") for _ in range(294)) + "TGA"     # 300 bp "reference allele" CDS
REF_START = 1001                                                           # contig coordinate of the reference slice start
CALL = {"strand": "+", "cds": [(REF_START, REF_START + 299)]}


def sam(seq, pos, cigar=None, flag=0, name="r"):
    return "%s\t%d\tref\t%d\t60\t%s\t*\t0\t0\t%s\t*" % (name, flag, pos, cigar or "%dM" % len(seq), seq)


def reads_over(contig, center, n_obs, n_alt, alt_seq=None, length=60):
    """n reads carrying `contig` and n carrying `alt_seq`, all spanning index `center` (0-based in the slice)."""
    out = []
    for k in range(n_obs):
        s = max(0, center - 20 - k % 5)
        out.append(sam(contig[s:s + length], s + 1, name="o%d" % k))
    for k in range(n_alt):
        s = max(0, center - 20 - k % 5)
        out.append(sam(alt_seq[s:s + length], s + 1, name="a%d" % k))
    return out


def make_case(called, ref=REF):
    aref, aobs = cb.banded_align(ref, called)
    return {"alignments": [cb.make_alignment("ref", "obs", aref, aobs, "immuannot-paf")]}


class TestSites(unittest.TestCase):
    def sub(self, s, k):
        return s[:k] + ("A" if s[k] != "A" else "C") + s[k + 1:]

    def test_heterozygous_snv_supported(self):
        called = self.sub(REF, 150)
        other = REF                                    # other haplotype keeps the reference base
        case = make_case(called)
        blk = sr.evaluate_case(cb, case, CALL, called, other, called, reads_over(called, 150, 14, 15, other), 1001, 1300)
        s = blk["sites"][0]
        self.assertTrue(s["het"])
        self.assertEqual((s["obs"], s["alt"], s["oth"], s["verdict"]), (14, 15, 0, "supported"))
        self.assertTrue(any(r["a"] == "alt" for r in s["window"]["reads"]))

    def test_homozygous_snv_and_unsupported_snv(self):
        called = self.sub(REF, 150)
        case = make_case(called)
        hom = sr.evaluate_case(cb, case, CALL, called, called, called, reads_over(called, 150, 30, 0), 1001, 1300)["sites"][0]
        self.assertEqual((hom["het"], hom["verdict"]), (False, "supported"))
        bad = sr.evaluate_case(cb, case, CALL, called, REF, called, reads_over(called, 150, 0, 30, REF), 1001, 1300)["sites"][0]
        self.assertEqual(bad["verdict"], "not seen in short reads")

    def test_called_insertion_vs_reads_with_deletion(self):
        called = REF[:100] + "G" + REF[100:]           # observed CDS has one extra base (needs a 301-bp contig slice)
        other = REF
        case = make_case(called)
        call = {"strand": "+", "cds": [(1001, 1001 + len(called) - 1)]}
        ev = case["alignments"][0]["events"][0]
        self.assertEqual(ev["type"], "ins")
        deleted = called[:100] + called[101:]          # reads from the other haplotype lack the base
        reads = reads_over(called, 100, 12, 0) + [sam(deleted[80:140], 81, cigar="20M1D40M", name="d%d" % k) for k in range(12)]
        # the 'D' reads align to the called contig with a deletion at the inserted base
        reads = [l if "D" not in l.split("\t")[5] else sam(called[80:100] + called[101:141], 81, cigar="20M1D40M", name=l.split("\t")[0]) for l in reads]
        s = sr.evaluate_case(cb, case, call, called, other, called, reads, 1001, 1001 + len(called) - 1)["sites"][0]
        self.assertTrue(s["het"])
        self.assertEqual((s["obs"], s["alt"]), (12, 12))

    def test_called_deletion_vs_reads_with_insertion(self):
        called = REF[:100] + REF[101:]                 # observed lacks one base
        case = make_case(called)
        call = {"strand": "+", "cds": [(1001, 1001 + len(called) - 1)]}
        ev = case["alignments"][0]["events"][0]
        self.assertEqual(ev["type"], "del")
        ins_reads = [sam(called[80:100] + REF[100] + called[100:140], 81, cigar="20M1I40M", name="i%d" % k) for k in range(10)]
        s = sr.evaluate_case(cb, case, call, called, REF, called, reads_over(called, 100, 12, 0) + ins_reads, 1001, 1001 + len(called) - 1)["sites"][0]
        self.assertTrue(s["het"])
        self.assertEqual((s["obs"], s["alt"]), (12, 10))

    def test_indel_in_homopolymer_tolerates_placement_differences(self):
        ref = REF[:90] + "TTTTT" + REF[95:]            # 5-T run at 90..94 in the reference allele
        called = ref[:90] + "TTTTTT" + ref[95:]        # called haplotype has a 6th T
        case = make_case(called, ref)
        call = {"strand": "+", "cds": [(1001, 1001 + len(called) - 1)]}
        ev = case["alignments"][0]["events"][0]
        self.assertEqual((ev["type"], ev["homopolymer"]), ("ins", True))
        # reads of the shorter haplotype, aligned with the deletion at the FAR end of the run (a different placement)
        reads = [sam(called[70:95] + called[96:130], 71, cigar="25M1D34M", name="d%d" % k) for k in range(10)]
        reads += reads_over(called, 92, 10, 0)
        s = sr.evaluate_case(cb, case, call, called, ref, called, reads, 1001, 1001 + len(called) - 1)["sites"][0]
        self.assertEqual((s["obs"], s["alt"], s["verdict"]), (10, 10, "supported"))

    def test_minus_strand_gene_complements_read_bases(self):
        called = self.sub(REF, 150)                    # CDS sense; the contig stores its reverse complement
        contig = cb.revcomp(called)
        call = {"strand": "-", "cds": [(1001, 1300)]}
        case = make_case(called)
        ev = case["alignments"][0]["events"][0]
        pos = sr.cds_contig_positions(call, len(called))
        self.assertEqual((pos[0], pos[-1]), (1300, 1001))
        p = pos[ev["qpos"] - 1] - 1001                 # 0-based index of the site inside the contig slice
        reads = reads_over(contig, p, 20, 0)
        s = sr.evaluate_case(cb, case, call, called, called, contig, reads, 1001, 1300)["sites"][0]
        self.assertEqual((s["obs"], s["oth"], s["verdict"]), (20, 0, "supported"))

    def test_cds_positions_across_exons_and_stop_codon(self):
        call = {"strand": "+", "cds": [(100, 104), (200, 204)]}
        self.assertEqual(sr.cds_contig_positions(call, 13), [100, 101, 102, 103, 104, 200, 201, 202, 203, 204, 205, 206, 207])

    def test_sam_parser_skips_secondary_and_parses_cigar(self):
        self.assertIsNone(sr.parse_sam_line(sam("ACGT", 5, flag=256), 1))
        r = sr.parse_sam_line(sam("ACGTACGT", 5, cigar="2M1D3M2I1M"), 1000)
        self.assertEqual((r.start, r.end), (1004, 1010))
        self.assertIn(1006, r.dele)
        self.assertEqual(r.ins, {1009: 2})


if __name__ == "__main__":
    unittest.main()
