"""End-to-end test of 80 extract on a synthetic Immuannot tree (no real data)."""
import gzip
import importlib.util
import json
import os
import random
import shutil
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


cb = load("cb80", "80_case_browser.py")
ex = load("cb80x", "80_case_browser_extract.py")
ex._CORE = cb


BLANK_EXON = False


def cs_from(aref, aobs):
    """minimap2-style cs with target = observed, query = reference."""
    out, run, i = [], 0, 0
    while i < len(aref):
        r, o = aref[i], aobs[i]
        if r == o:
            run += 1
            i += 1
            continue
        if run:
            out.append(":%d" % run)
            run = 0
        if r != "-" and o != "-":
            out.append("*%s%s" % (o.lower(), r.lower()))
            i += 1
        else:
            ref_only = o == "-"      # reference has bases the observed CDS lacks -> '+' (query-only)
            j = i
            while j < len(aref) and aref[j] != aobs[j] and ((aobs[j] == "-") == ref_only) and "-" in (aref[j], aobs[j]):
                j += 1
            out.append(("+" + aref[i:j].lower()) if ref_only else ("-" + aobs[i:j].lower()))
            i = j
    if run:
        out.append(":%d" % run)
    return "".join(out)


def mkseq(rng, n):
    return "ATG" + "".join(rng.choice("ACGT") for _ in range(n - 6)) + "TGA"


def write_gz(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as f:
        f.write(text)


def person(root, pid, gene, ref_name, ref, obs, contig, start, hap="1"):
    od = os.path.join(root, pid, "immuannot_output")
    aref, aobs = cb.banded_align(ref, obs)
    cs = cs_from(aref, aobs)
    nm = sum(1 for x, y in zip(aref, aobs) if x != y)
    g = '%s\tIPD\tgene\t%d\t%d\t.\t+\t.\tgene_id "IAG1"; template_allele "%s"; template_distance 9; gene_name "%s";\n' % (
        contig, start, start + 5000, ref_name, gene)
    t = ('%s\tIPD\ttranscript\t%d\t%d\t.\t+\t.\tgene_id "IAG1"; transcript_id "IAT1"; gene_name "%s"; consensus "%s"; '
         'cds_distance %d; cds_mut "%s|%s|x";\n') % (contig, start, start + 5000, gene, ref_name, nm, ref_name, cs)
    mid = len(obs) // 2
    e = ""
    for n, (s, en) in enumerate([(start + 100, start + 100 + mid - 1), (start + 3000, start + 3000 + len(obs) - mid - 4)], 1):
        e += '%s\tIPD\texon\t%d\t%d\t.\t+\t.\tgene_id "IAG1"; transcript_id "IAT1"; exon_number "%s";\n' % (contig, s, en, n if n != 2 or not BLANK_EXON else "")
        e += '%s\tIPD\tCDS\t%d\t%d\t.\t+\t0\tgene_id "IAG1"; transcript_id "IAT1";\n' % (contig, s, en)
    write_gz(od + "/hap%s.gtf.gz" % hap, "## format: gtf\n" + g + t + e)
    rec = "%s_%s_1" % (contig, gene)
    write_gz(od + "/hap%s/cds.fa.gz" % hap, ">%s GT_AG\n%s\n>other_X_1\nAAAA\n" % (rec, obs))
    write_gz(od + "/hap%s/mm2.ipd.cds.paf.gz" % hap,
             "%s\t%d\t0\t%d\t+\t%s\t%d\t0\t%d\t0\t0\t60\tNM:i:%d\tAS:i:100\tcs:Z:%s\n" % (
                 ref_name, len(ref), len(ref), rec, len(obs), len(obs), nm, cs))
    return nm


class TestExtract(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d)

    def test_blank_exon_numbers_are_numbered_by_position(self):
        global BLANK_EXON
        BLANK_EXON = True
        try:
            self.test_pipeline()
        finally:
            BLANK_EXON = False

    def test_pipeline(self):
        rng = random.Random(5)
        d = self.d
        gene, ref_name = "HLA-A", "HLA-A*01:01:01:01"
        ref = mkseq(rng, 300)
        write_gz(d + "/ref/CDSseq/HLA-A.fa.gz", ">%s HLA00001 frame=1 300bp\n%s\n" % (ref_name, ref))
        rows = []
        for k, pid in enumerate(["pers1", "pers2", "pers3"]):
            obs = list(ref)
            obs[50 + k] = "A" if obs[50 + k] != "A" else "C"
            del obs[200]
            nm = person(d + "/people", pid, gene, ref_name, ref, "".join(obs), "ctgx%d" % k, 1000)
            rows.append([pid, "hap1", "ctgx%d" % k, gene, "1", ref_name, str(nm), "%s|x|y" % ref_name, "9", "1000", "6000"])
        hdr = ["person_id", "hap", "contig", "gene", "copy_index", "consensus", "cds_distance", "cds_mut",
               "template_distance", "gene_start", "gene_end"]
        with open(d + "/rich.tsv", "w") as f:
            f.write("\t".join(hdr) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")
        args = type("A", (), dict(out_dir=d + "/out", hla_root=d + "/people", kir_root=d + "/nokir", refdata=d + "/ref",
                                  rich=d + "/rich.tsv", origin=d + "/none", relatedness=d + "/none", per_cell=3,
                                  kir_people=0, kir_per_cell=2, mendel=0, mendel_controls=0, switch=0,
                                  switch_controls=0, seed=1))()
        cases = ex.run(args, cb)
        self.assertEqual(len(cases), 3)
        for c in cases:
            a = c["alignments"][0]
            self.assertEqual(a["problems"], 0, "cs replay must reproduce both sequences exactly")
            self.assertEqual(a["stats"]["nm"], int(c["cds_distance"]))
            self.assertEqual(c["gene"], "A")
            self.assertEqual([e["n"] for e in a["exon_cols"]], [1, 2])
        txt = open(d + "/out/cases.json").read()
        for secret in ("pers1", "pers2", "pers3", "ctgx0", "ctgx1"):
            self.assertNotIn(secret, txt, "no person ids or contig names in the shared JSON")
        self.assertEqual(oct(os.stat(d + "/out/case_map.private.tsv").st_mode & 0o777), "0o600")
        j = json.loads(txt)
        out = cb.build_html(j["cases"], d + "/html", j["meta"])
        self.assertIn("c001", open(out).read())


class TestTrioAndSwitch(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d)

    def test_mendel_and_switch(self):
        import csv
        rng = random.Random(9)
        d = self.d
        gene = "HLA-A"
        r1, r2 = "HLA-A*01:01:01:01", "HLA-A*02:01:01:01"
        ref1, ref2 = mkseq(rng, 300), mkseq(rng, 300)
        write_gz(d + "/ref/CDSseq/HLA-A.fa.gz", ">%s x\n%s\n>%s x\n%s\n" % (r1, ref1, r2, ref2))
        mut = lambda s, k: s[:k] + ("A" if s[k] != "A" else "C") + s[k + 1:]
        person(d + "/people", "kid", gene, r1, ref1, mut(ref1, 60), "ctgk", 1000, "1")      # child hap1: A*01 + 1 edit
        person(d + "/people", "parA", gene, r2, ref2, mut(ref2, 70), "ctga", 1000, "1")     # parent A carries only A*02
        person(d + "/people", "parA", gene, r2, ref2, ref2[:100] + ref2[101:], "ctga", 1000, "2")
        hdr = ["person_id", "hap", "contig", "gene", "copy_index", "consensus", "cds_distance", "cds_mut", "template_distance",
               "gene_start", "gene_end"]
        rows = [["kid", "hap1", "ctgk", gene, "1", r1, "1", "x", "9", "1000", "6000"],
                ["parA", "hap1", "ctga", gene, "1", r2, "1", "x", "9", "1000", "6000"],
                ["parA", "hap2", "ctga", gene, "1", r2, "1", "x", "9", "1000", "6000"]]
        # a 5-gene contig for the switch builder (genes other than A carry no sequences)
        for k, (g, pos) in enumerate([("B", 100000), ("C", 200000), ("DRB1", 3000000), ("DQB1", 3100000)]):
            rows.append(["kid", "hap1", "ctgk", "HLA-" + g, "1", "HLA-%s*01:01" % g, "", "", "1", str(pos), str(pos + 5000)])
        with open(d + "/rich.tsv", "w") as f:
            f.write("\t".join(hdr) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")
        with open(d + "/origin.tsv", "w") as f:
            f.write("child\tplat\thap\tgene\tcontig\tpat_k\tmat_k\torigin\n")
            for g, o in [("A", "pat"), ("B", "pat"), ("C", "mat"), ("DRB1", "mat"), ("DQB1", "amb")]:
                f.write("kid\trevio\thap1\t%s\tctgk\t%d\t%d\t%s\n" % (g, 90 if o == "pat" else 2, 90 if o == "mat" else 2, o))
        with open(d + "/rel.tsv", "w") as f:
            f.write("i.s\tj.s\tkin\nkid\tparA\t0.25\nkid\tparB\t0.25\n")
        args = type("A", (), dict(out_dir=d + "/out", hla_root=d + "/people", kir_root=d + "/nokir", refdata=d + "/ref",
                                  rich=d + "/rich.tsv", origin=d + "/origin.tsv", relatedness=d + "/rel.tsv", per_cell=0,
                                  kir_people=0, kir_per_cell=0, mendel=5, mendel_controls=0, switch=5, switch_controls=0, seed=1))()
        cases = ex.run(args, cb)
        by = {c["tab"]: c for c in cases}
        m = by["mendel"]
        self.assertEqual(m["verdict"], "child allele absent from assigned parent")
        titles = [a["title"] for a in m["alignments"]]
        self.assertEqual(len(titles), 5, titles)      # child + 2 parent haps (Immuannot) + 2 direct
        self.assertTrue(all(a["problems"] == 0 for a in m["alignments"] if a["source"] == "immuannot-paf"))
        self.assertGreater(m["alignments"][3]["stats"]["nm"], 0)
        s = by["switch"]
        self.assertEqual(s["n_switch"], 1)
        self.assertEqual([t["gene"] for t in s["track"]], ["A", "B", "C", "DRB1", "DQB1"])
        self.assertEqual(s["switches"][0]["gap_kb"], 100.0)     # B (100 kb) -> C (200 kb); A starts the track at 1 kb
        self.assertNotIn("kid", open(d + "/out/cases.json").read())


if __name__ == "__main__":
    unittest.main()
