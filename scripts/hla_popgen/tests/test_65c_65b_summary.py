"""65c: summary-only 65b tables on synthetic INTERNAL inputs (no person ids or small counts leak)."""
import importlib.util
import os
import random
import sys
import tempfile
import unittest

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
_spec = importlib.util.spec_from_file_location("_s65c", os.path.join(os.path.dirname(HERE), "65c_65b_summary.py"))
s = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s)
tg = s.tg


class Test65c(unittest.TestCase):
    def test_tables(self):
        tmp = tempfile.mkdtemp()
        rng = random.Random(65)
        ids = [str(6_000_000 + i) for i in range(300)]
        coh = pd.DataFrame(dict(person_id=ids, ancestry_pred=["afr" if i % 2 else "sas" for i in range(300)],
                                platform="revio"))
        cpath = os.path.join(tmp, "cohort.tsv")
        coh.to_csv(cpath, sep="\t", index=False)
        sample = ids[:200]
        pd.DataFrame(dict(person_id=sample)).to_csv(os.path.join(tmp, "people.tsv"), sep="\t", index=False)
        asm = os.path.join(tmp, "asm")
        os.makedirs(asm)
        reads, ag, hp, fl = [], [], [], []
        for i, p in enumerate(ids):
            c = rng.choice([0, 1, 1, 2])
            ni = {0: 0, 1: 13, 2: 20}[c]
            for feat in ("TRBV4-3", "TRBV3-2"):
                reads.append(tg._row(tg.GENE_COLS, person_id=p, feature=feat, copies=c, status="ok",
                                     n_ins=ni, n_ref=20 - ni))
            if p in sample and i < 150:  # 150 done, 5 discordant
                for feat in ("TRBV4-3", "TRBV3-2"):
                    ag.append(tg._row(tg.GENE_COLS, person_id=p, feature=feat, copies=(c + 1) % 3 if i < 5 else c,
                                      status="ok"))
                for h in (1, 2):
                    hp.append(tg._row(tg.HAP_COLS, person_id=p, hap=h, complete=True))
            if 140 <= i < 160:
                fl.append(dict(person_id=p, reason="timeout"))
        pd.DataFrame(reads).to_csv(os.path.join(tmp, "reads.tsv"), sep="\t", index=False)
        pd.DataFrame(ag).to_csv(os.path.join(asm, "per_person_gene.tsv"), sep="\t", index=False)
        pd.DataFrame(hp).to_csv(os.path.join(asm, "per_haplotype.tsv"), sep="\t", index=False)
        pd.DataFrame(fl).to_csv(os.path.join(asm, "failures.tsv"), sep="\t", index=False)
        out = os.path.join(tmp, "out")
        s.main(["--reads-gene", os.path.join(tmp, "reads.tsv"), "--asm-dir", asm, "--people-file",
                os.path.join(tmp, "people.tsv"), "--cohort", cpath, "--out-dir", out])
        rd = lambda n: pd.read_csv(os.path.join(out, n), sep="\t", dtype=str, keep_default_na=False)  # noqa
        comp = rd("65c_completion.tsv").set_index("stratum")
        self.assertEqual(comp.loc["ALL", "n_sampled"], "200")
        self.assertEqual(comp.loc["ALL", "pct_genotyped"], "75.0")
        f = rd("65c_failures.tsv").iloc[0]
        self.assertEqual((f["n_people_ever_failed"], f["n_people_still_failed"]), ("20", "<20"))
        c = rd("65c_concordance.tsv")
        a = c[(c["reads_genotyper"] == "stored") & (c["feature"] == "TRBV4-3") & (c["stratum"] == "ALL")].iloc[0]
        self.assertEqual(a["n_both_called"], "150")
        self.assertEqual(a["pct_concordant"], "96.7")
        h = rd("65c_hwe.tsv")
        self.assertEqual(set(h["caller"]), {"assembly", "reads_cuts", "reads_mixture"})
        for name in os.listdir(out):
            txt = open(os.path.join(out, name)).read()
            self.assertFalse(any(p in txt for p in ids), name)


if __name__ == "__main__":
    unittest.main()
