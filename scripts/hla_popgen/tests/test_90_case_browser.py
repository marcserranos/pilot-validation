"""Tests for 90_case_browser core: Immuannot cs replay, event classification, direct aligner."""
import importlib.util
import json
import os
import random
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("cb80", os.path.join(HERE, "..", "90_case_browser.py"))
cb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cb)

# 12 codons: ATG GCC AAA AAA GAG CTG TTT GGG CAG TAC GTT TGA
REF = "ATGGCCAAAAAAGAGCTGTTTGGGCAGTACGTTTGA"


def paf(ref, obs, cs, qs=0, qe=None, ts=0, te=None):
    return cb.alignment_from_paf(ref, obs, qs, len(ref) if qe is None else qe, ts,
                                 len(obs) if te is None else te, cs)


class TestReplay(unittest.TestCase):
    def test_substitution_semantics(self):
        # observed has G where reference has A (target=obs x, query=ref y: '*ga')
        obs = REF[:13] + "G" + REF[14:]
        aref, aobs, c5, c3, bad = paf(REF, obs, ":13*ga:22")
        self.assertEqual(bad, 0)
        self.assertEqual((aref, aobs), (REF, obs))

    def test_indel_semantics_and_problems(self):
        obs = REF[:8] + REF[9:]          # observed lacks one of the A's (a deletion in the observed CDS)
        aref, aobs, _, _, bad = paf(REF, obs, ":8+a:%d" % (len(REF) - 9))
        self.assertEqual(bad, 0)
        self.assertEqual(aref.replace("-", ""), REF)
        self.assertEqual(aobs.replace("-", ""), obs)
        self.assertEqual(len(aref), len(aobs))
        _, _, _, _, bad = paf(REF, obs, ":9+a:%d" % (len(REF) - 9))   # wrong base consumed -> flagged
        self.assertGreater(bad, 0)

    def test_clipped_ends_kept_as_flanks(self):
        obs = "GG" + REF
        aref, aobs, c5, c3, bad = paf(REF, obs, ":%d" % len(REF), ts=2, te=len(obs))
        self.assertEqual((c5, c3, bad), (2, 0, 0))
        self.assertEqual(aobs[:2], "GG")
        self.assertEqual(aref[:2], "--")


class TestAnalyse(unittest.TestCase):
    def test_missense_vs_synonymous(self):
        # codon 5 GAG(E) -> GAA (synonymous) and codon 8 GGG(G) -> GAG (E, missense)
        obs = list(REF)
        obs[14] = "A"
        obs[22] = "A"
        out = cb.analyse(REF, "".join(obs))
        eff = {e["codon"]: e["effect"] for e in out["events"]}
        self.assertEqual(eff, {5: "synonymous", 8: "missense"})
        self.assertEqual(out["stats"]["nm"], 2)
        self.assertEqual(out["stats"]["missense"], 1)

    def test_homopolymer_deletion_and_frameshift(self):
        aref = REF
        aobs = REF[:8] + "-" + REF[9:]   # one A of the AAAAAA run (positions 6-11) removed
        ev = cb.analyse(aref, aobs)["events"]
        self.assertEqual(len(ev), 1)
        self.assertEqual(ev[0]["type"], "del")
        self.assertTrue(ev[0]["homopolymer"])
        self.assertEqual(ev[0]["run"], 6)
        self.assertEqual(ev[0]["effect"], "frameshift")

    def test_inframe_codon_insertion_is_not_a_frameshift(self):
        aref = REF[:12] + "---" + REF[12:]
        aobs = REF[:12] + "AAG" + REF[12:]
        ev = cb.analyse(aref, aobs)["events"][0]
        self.assertEqual((ev["type"], ev["len"], ev["effect"]), ("ins", 3, "inframe indel"))

    def test_tandem_slip_flagged(self):
        ref = "ATGCACACAGTC"
        out = cb.analyse(ref[:5] + "--" + ref[5:], ref[:5] + "CA" + ref[5:])
        self.assertTrue(out["events"][0]["tandem"] or out["events"][0]["homopolymer"])

    def test_clip_columns_not_counted_and_exons(self):
        obs = "GG" + REF
        aref, aobs, c5, c3, _ = paf(REF, obs, ":%d" % len(REF), ts=2, te=len(obs))
        out = cb.analyse(aref, aobs, c5, c3, exons=[(1, 20), (2, len(obs) - 20)], ars_exons=(2,))
        self.assertEqual(out["stats"]["nm"], 0)
        self.assertEqual([e["n"] for e in out["exon_cols"]], [1, 2])

    def test_exon_and_ars_assigned(self):
        obs = list(REF)
        obs[30] = "A"
        out = cb.analyse(REF, "".join(obs), exons=[(1, 12), (2, 12), (3, 12)], ars_exons=(2, 3))
        self.assertEqual((out["events"][0]["exon"], out["events"][0]["ars"]), (3, True))

    def test_amino_acid_tracks(self):
        out = cb.analyse(REF, REF)
        self.assertEqual(out["aa_ref"].replace(" ", ""), cb.translate(REF))
        self.assertEqual(out["stats"]["pct_cds"], 0.0)


class TestAligner(unittest.TestCase):
    def test_invariants_random(self):
        rng = random.Random(3)
        for _ in range(25):
            a = "".join(rng.choice("ACGT") for _ in range(rng.randrange(60, 160)))
            b = list(a)
            for _ in range(rng.randrange(0, 8)):
                p = rng.randrange(len(b))
                t = rng.random()
                if t < .5:
                    b[p] = rng.choice("ACGT")
                elif t < .75:
                    b.insert(p, rng.choice("ACGT"))
                else:
                    del b[p]
            b = "".join(b)
            ra, rb = cb.banded_align(a, b)
            self.assertEqual(ra.replace("-", ""), a)
            self.assertEqual(rb.replace("-", ""), b)
            self.assertEqual(len(ra), len(rb))
            self.assertFalse(any(x == "-" and y == "-" for x, y in zip(ra, rb)))

    def test_known_cases(self):
        self.assertEqual(cb.banded_align("ACGTAAAAGT", "ACGTAAAGT")[1].count("-"), 1)
        ra, rb = cb.banded_align(REF, REF[:20] + "A" + REF[21:])
        self.assertEqual(sum(1 for x, y in zip(ra, rb) if x != y), 1 if REF[20] != "A" else 0)


class TestNormalise(unittest.TestCase):
    def test_exon_fill_and_mendel_summary_are_idempotent(self):
        direct = cb.make_alignment("p", "c", *cb.banded_align(REF, REF[:8] + REF[9:]), "direct")
        c = {"tab": "mendel", "facts": [["closest parent haplotype", "1 edits"]], "alignments": [direct]}
        direct["exon_cols"] = [{"n": None, "c0": 0, "c1": 10, "ars": False}, {"n": 2, "c0": 10, "c1": 30, "ars": False}]
        cb.normalise_cases([c])
        self.assertEqual(direct["exon_cols"][0]["n"], 1)
        self.assertIn("1 edit: 1 indel (1 in homopolymers, 1 frameshifting)", c["facts"][0][1])
        before = json.dumps(c)
        cb.normalise_cases([c])
        self.assertEqual(json.dumps(c), before)


class TestClassify(unittest.TestCase):
    def test_edit_classes_and_cs_only_classification(self):
        self.assertEqual(cb.cs_edit_class(":10*ga:20"), "snv-only")
        self.assertEqual(cb.cs_edit_class(":10+a:20"), "indel-only")
        self.assertEqual(cb.cs_edit_class(":10*ga:5-c:20"), "mixed")
        self.assertEqual(cb.cs_edit_class(":30"), "identical")

    def test_replay_from_reference_alone(self):
        obs = REF[:13] + "G" + REF[14:]
        aref, aobs = cb.replay_from_ref(REF, ":13*ga:22")
        self.assertEqual((aref, aobs), (REF, obs))
        self.assertIsNone(cb.replay_from_ref(REF, ":13*ga:10"), "a clipped alignment cannot be rebuilt")
        self.assertIsNone(cb.replay_from_ref(REF, ":13*gc:22"), "reference disagreement is refused")

    def test_short_read_verdicts(self):
        self.assertEqual(cb.sr_verdict(30, 1, 0, het=False), "supported")
        self.assertEqual(cb.sr_verdict(2, 28, 0, het=False), "contradicted")
        self.assertEqual(cb.sr_verdict(14, 16, 0, het=True), "supported")
        self.assertEqual(cb.sr_verdict(0, 30, 0, het=True), "not seen in short reads")
        self.assertEqual(cb.sr_verdict(22, 8, 0, het=True), "skewed")
        self.assertEqual(cb.sr_verdict(3, 2, 0, het=True), "low depth")
        self.assertEqual(cb.sr_verdict(15, 15, 0, het=False), "half the reads (other haplotype may differ)")
        self.assertEqual(cb.sr_verdict(25, 5, 0, het=False), "weak")
        self.assertEqual(cb.sr_verdict(15, 0, 15, het=None), "supported")
        self.assertEqual(cb.sr_verdict(0, 0, 30, het=None), "not seen in short reads")


class TestMisc(unittest.TestCase):
    def test_novelty_class(self):
        self.assertEqual(cb.novelty_class("HLA-A*01:01:01:01"), "known")
        self.assertEqual(cb.novelty_class("HLA-A*new"), "novel-protein")
        self.assertEqual(cb.novelty_class("HLA-A*01:new"), "novel-protein")
        self.assertEqual(cb.novelty_class("HLA-A*01:01:new"), "novel-synonymous")
        self.assertEqual(cb.novelty_class("HLA-A*01:01:01:new"), "novel-beyond-CDS")

    def test_translate_and_revcomp(self):
        self.assertEqual(cb.translate("ATGGCCTGA"), "MA*")
        self.assertEqual(cb.revcomp("AACG"), "CGTT")


if __name__ == "__main__":
    unittest.main()
