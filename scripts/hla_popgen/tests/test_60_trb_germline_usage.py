#!/usr/bin/env python3
"""Tests for 60_trb_germline_usage.py -- synthetic data only, no AoU data.

Planted truth: an insertion haplotype carrying TRBV4-3, TRBV3-2, TRBV6-3 and TRBV7-8 (insertion
allele frequency 0.4 in SAS, 0.7 in EUR, 0.55 in AFR), an independent TRBV12-4 deletion (allele
frequency 0.15), usage proportional to copies (dosage), and a residual: 4% of people who carry
TRBV4-3 never use it anyway (planted non-functional allele). Two copy-number sources:
  58 table : insertion genes absent (not on primary chr7); their paralogs TRBV4-2/6-2/3-1 read
             1 + copies/2 relative depth (the proxy must recover the insertion copies);
  65 table : direct per-haplotype calls for every gene.
Covers: CN0 explains never-used (PPV ~1), residual ~4% among CN>=1, Kendall dosage, haplotype
concordance (by construction under the proxy, real under 65), pickle and long-TSV usage formats,
the not-found message listing candidates, PUBLIC masking, no person ids in PUBLIC output, INTERNAL
refusing repo paths.

Run: python3 -m pytest scripts/hla_popgen/tests/test_60_trb_germline_usage.py -q
"""
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
spec = importlib.util.spec_from_file_location("tg60", os.path.join(SRC, "60_trb_germline_usage.py"))
tg = importlib.util.module_from_spec(spec)
sys.modules["tg60"] = tg
spec.loader.exec_module(tg)

import _disclosure as _disc  # noqa: E402

INS_FREQ = {"sas": 0.4, "eur": 0.7, "afr": 0.55}
BASE_P = {"TRBV4-3": 0.02, "TRBV3-2": 0.012, "TRBV6-3": 0.015, "TRBV7-8": 0.02, "TRBV12-4": 0.03,
          "TRBV20-1": 0.10, "TRBV7-2": 0.05, "TRBV4-2": 0.01, "TRBV6-2": 0.02, "TRBV3-1": 0.02}
RESIDUAL = 0.04


class Synth:
    def __init__(self, tmp, n=1500, seed=11):
        rng = np.random.default_rng(seed)
        self.tmp = tmp
        self.pids = [str(3_000_000 + i) for i in range(n)]
        anc = rng.choice(list(INS_FREQ), n)
        f = np.array([INS_FREQ[a] for a in anc])
        self.ins = rng.binomial(2, f)
        self.d124 = 2 - rng.binomial(2, 0.15, n)
        self.dead43 = (rng.random(n) < RESIDUAL) & (self.ins >= 1)
        cn = {g: np.full(n, 2) for g in BASE_P}
        for g in ("TRBV4-3", "TRBV3-2", "TRBV6-3", "TRBV7-8"):
            cn[g] = self.ins
        cn["TRBV12-4"] = self.d124
        self.cn = cn

        # 58 table: primary-assembly genes only
        rows = []
        for i, pid in enumerate(self.pids):
            for g in ("TRBV4-2", "TRBV6-2", "TRBV3-1", "TRBV12-4", "TRBV20-1", "TRBV7-2"):
                if g in ("TRBV4-2", "TRBV6-2", "TRBV3-1"):
                    rel = 1 + self.ins[i] / 2 + rng.normal(0, 0.06)
                else:
                    rel = cn[g][i] / 2 + rng.normal(0, 0.04)
                copies = 0 if rel < 0.25 else (1 if rel < 0.75 else 2)
                rows.append(dict(person_id=pid, gene=g, rel_depth=round(rel, 4), copies=copies,
                                 assessable=g not in ("TRBV6-2", "TRBV12-4")))
        self.cn58 = os.path.join(tmp, "trb_cn_per_person.tsv")
        with open(self.cn58, "w") as fh:
            fh.write(_disc.INTERNAL_HEADER_COMMENT + "\n")
            pd.DataFrame(rows).to_csv(fh, sep="\t", index=False)
        # 65 table (its per_person_gene.tsv layout): 'feature' column, status, TRBV6-2 and TRBV6-3
        # merged into one class plus the derived '[copy2]' insertion feature; 3% incomplete people
        rows = []
        self.incomplete = set(p for p in self.pids if rng.random() < 0.03)
        feats = [g for g in BASE_P if g not in ("TRBV6-2", "TRBV6-3")]
        for i, pid in enumerate(self.pids):
            fc = [(g, int(cn[g][i])) for g in feats]
            fc += [("TRBV6-2/6-3", 2), ("TRBV6-2/6-3[copy2]", int(self.ins[i]))]
            for g, c in fc:
                rows.append(dict(person_id=pid, feature=g, hap1_present=int(c >= 1), hap2_present=int(c == 2),
                                 copies=c, status="incomplete" if pid in self.incomplete else "ok",
                                 method="minimap2:sr"))
        self.cn65 = os.path.join(tmp, "trb_assembly_genotype.tsv")
        pd.DataFrame(rows).to_csv(self.cn65, sep="\t", index=False)

        # usage counts: multinomial over genes, weight = base p x copies/2, dead TRBV4-3 -> 0
        genes = list(BASE_P)
        V = np.zeros((n, len(genes)), int)
        for i in range(n):
            w = np.array([BASE_P[g] * cn[g][i] / 2 for g in genes])
            if self.dead43[i]:
                w[genes.index("TRBV4-3")] = 0
            w = np.append(w, 1 - w.sum())  # all other TRBV genes
            V[i] = rng.multinomial(int(rng.integers(2500, 5000)), w)[:-1]
        other = rng.integers(1000, 2000, n)
        self.V = pd.DataFrame(V, index=pd.Index(self.pids, name="research_id"), columns=genes)
        self.V["TRBV5-1*01"] = other
        self.pkl = os.path.join(tmp, "vj_usage.pkl")
        pd.to_pickle({"V": self.V, "J": pd.DataFrame(index=self.V.index)}, self.pkl)
        long = self.V.reset_index().melt(id_vars="research_id", var_name="v_gene", value_name="count")
        long["v_gene"] = long["v_gene"].map(lambda g: g if "*" in g else g + "*01")
        self.long_tsv = os.path.join(tmp, "usage_long.tsv")
        long.to_csv(self.long_tsv, sep="\t", index=False)

        self.cohort = os.path.join(tmp, "cohort.tsv")
        pd.DataFrame(dict(person_id=self.pids, in_lr=True, ancestry_pred=anc,
                          platform="revio")).to_csv(self.cohort, sep="\t", index=False)
        self.rel = os.path.join(tmp, "rel.tsv")
        pd.DataFrame({"i.s": [self.pids[0]], "j.s": [self.pids[1]], "kin": [0.25]}).to_csv(
            self.rel, sep="\t", index=False)

    def run(self, out, internal, cn, usage, *extra):
        return tg.main(["--cn-table", cn, "--usage", usage, "--cohort", self.cohort,
                        "--relatedness-table", self.rel, "--out-dir", out, "--internal-dir", internal,
                        *extra])


def rd(out, name):
    return pd.read_csv(os.path.join(out, name), sep="\t", dtype=str, keep_default_na=False)


class Test60(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.s = Synth(cls.tmp)
        cls.out58 = os.path.join(cls.tmp, "out58")
        cls.s.run(cls.out58, os.path.join(cls.tmp, "int58"), cls.s.cn58, cls.s.pkl)
        cls.out65 = os.path.join(cls.tmp, "out65")
        cls.int65 = os.path.join(cls.tmp, "int65")
        cls.s.run(cls.out65, cls.int65, cls.s.cn65, cls.s.long_tsv, "--haplotype-genes",
                  ",".join(tg.HAPLOTYPE_GENES + ("TRBV6-2/6-3[copy2]",)))

    def test_source_detection_and_proxy(self):
        rs = dict(rd(self.out58, "run_summary.tsv").values)
        self.assertEqual(rs["cn_source"], "58")
        self.assertIn("proxy", rs["proxy_label"])
        self.assertEqual(dict(rd(self.out65, "run_summary.tsv").values)["cn_source"], "65")
        j = pd.read_csv(os.path.join(self.tmp, "int58", "joined_per_person.tsv"), sep="\t", comment="#",
                        dtype={"person_id": str})
        p = j[j["gene"] == "TRBV4-3"].set_index("person_id")
        self.assertEqual(set(p["cn_method"]), {tg.PROXY_METHOD})
        truth = pd.Series(self.s.ins, index=self.s.pids).loc[p.index]
        self.assertGreater((p["copies"] == truth).mean(), 0.97)

    def test_65_layout(self):
        j = pd.read_csv(os.path.join(self.int65, "joined_per_person.tsv"), sep="\t", comment="#",
                        dtype={"person_id": str})
        self.assertFalse(set(j["person_id"]) & self.s.incomplete)
        m = j[j["gene"] == "TRBV6-2/6-3"].set_index("person_id")
        v = self.s.V.loc[m.index]
        np.testing.assert_array_equal(m["usage_value"].to_numpy(), (v["TRBV6-2"] + v["TRBV6-3"]).to_numpy())
        tot = v.sum(axis=1).to_numpy()
        np.testing.assert_array_equal(m["usage_total"].to_numpy(), tot)  # merged column not double counted
        self.assertEqual(set(j["cn_method"]), {"minimap2:sr"})
        h = rd(self.out65, "haplotype_genes.tsv")
        r = h[(h["gene_a"] == "TRBV4-3") & (h["gene_b"] == "TRBV6-2/6-3[copy2]") & (h["ancestry"] == "ALL")]
        self.assertEqual(float(r.iloc[0]["pct_identical_copies"]), 100.0)
        self.assertEqual(tg.class_members("TRBV6-2/6-3"), ["TRBV6-2", "TRBV6-3"])
        self.assertEqual(tg.class_members("TRBV6-2/6-3[copy2]"), [])

    def test_cn0_explains_never_used(self):
        for out in (self.out58, self.out65):
            d = rd(out, "diagnostic.tsv")
            r = d[(d["gene"] == "TRBV4-3") & (d["ancestry"] == "ALL")].iloc[0]
            self.assertGreater(float(r["ppv_pct"]), 97, out)
            self.assertAlmostEqual(float(r["residual_never_used_cn_ge1_pct"]), 100 * RESIDUAL, delta=1.6)
            # never-used = CN0 people + residual: sensitivity = share explained by CN0
            self.assertGreater(float(r["sensitivity_pct"]), 70)
            self.assertLess(float(r["sensitivity_pct"]), 95)
            self.assertLess(float(r["expected_chance_zero_cn_ge1_pct"]), 1.0)
            r2 = d[(d["gene"] == "TRBV3-2") & (d["ancestry"] == "ALL")].iloc[0]  # no residual planted
            self.assertGreater(float(r2["sensitivity_pct"]), 97)

    def test_contingency_by_ancestry(self):
        c = rd(self.out65, "contingency.tsv")
        for anc, f in (("SAS", 0.4), ("EUR", 0.7)):
            s = c[(c["gene"] == "TRBV3-2") & (c["ancestry"] == anc)]
            n = {int(k): int(v) for k, v in zip(s["copies"], s["n_people"]) if v.isdigit()}
            self.assertAlmostEqual(n[0] / sum(n.values()), (1 - f) ** 2, delta=0.05)
            pct0 = float(s[s["copies"] == "0"]["pct_never_used"].iloc[0])
            self.assertGreater(pct0, 97)

    def test_dosage(self):
        d = rd(self.out65, "dosage.tsv")
        r = d[(d["gene"] == "TRBV12-4") & (d["ancestry"] == "ALL") & (d["copies"] == "2")].iloc[0]
        self.assertGreater(float(r["kendall_tau_b_cn_ge1"]), 0.2)
        self.assertLess(float(r["kendall_p_cn_ge1"]), 1e-6)
        med = {k: float(v) for k, v in zip(
            d[(d["gene"] == "TRBV12-4") & (d["ancestry"] == "ALL")]["copies"],
            d[(d["gene"] == "TRBV12-4") & (d["ancestry"] == "ALL")]["usage_pct_median"]) if v != ""}
        self.assertAlmostEqual(med["1"] / med["2"], 0.5, delta=0.1)

    def test_haplotype_concordance(self):
        h = rd(self.out65, "haplotype_genes.tsv")
        a = h[(h["gene_a"] == "TRBV4-3") & (h["gene_b"] == "TRBV3-2") & (h["ancestry"] == "ALL")].iloc[0]
        self.assertEqual(float(a["pct_identical_copies"]), 100.0)
        self.assertEqual(a["identical_by_construction"], "False")
        self.assertGreater(float(a["cohen_kappa"]), 0.99)
        b = h[(h["gene_a"] == "TRBV4-3") & (h["gene_b"] == "TRBV12-4") & (h["ancestry"] == "ALL")].iloc[0]
        self.assertLess(float(b["cohen_kappa"]), 0.1)
        h58 = rd(self.out58, "haplotype_genes.tsv")
        c = h58[(h58["gene_a"] == "TRBV4-3") & (h58["gene_b"] == "TRBV3-2")].iloc[0]
        self.assertEqual(c["identical_by_construction"], "True")
        self.assertIn("TRBV6-2", set(h58["gene_a"]) | set(h58["gene_b"]))

    def test_public_masking_and_no_ids(self):
        for out in (self.out58, self.out65):
            for f in os.listdir(out):
                df = pd.read_csv(os.path.join(out, f), sep="\t", dtype=str, keep_default_na=False)
                txt = open(os.path.join(out, f)).read()
                for pid in (self.s.pids[2], self.s.pids[700], self.s.pids[-1]):
                    self.assertNotIn(pid, txt, f)
                for col in df.columns:
                    if col.startswith(("n_", "TP_", "FP_", "FN_", "TN_")):
                        nums = pd.to_numeric(df[col], errors="coerce")
                        self.assertFalse(((nums > 0) & (nums < 20)).any(), (f, col))
        out, internal = os.path.join(self.tmp, "o_small"), os.path.join(self.tmp, "i_small")
        self.s.run(out, internal, self.s.cn65, self.s.pkl, "--people-limit", "40")
        c = rd(out, "contingency.tsv")
        self.assertIn("<20", set(c["n_people"]))
        self.assertTrue((c.loc[c["n_people"] == "<20", "pct_never_used"] == "").all())
        d = rd(out, "dosage.tsv")
        self.assertTrue((d.loc[d["n_people"] == "<20", "usage_pct_median"] == "").all())
        # complementary rule: masked never-used count never sits next to an exact used count
        full = rd(self.out65, "contingency.tsv")
        m = (full["n_never_used"] == "<20") & full["n_people"].str.isdigit()
        self.assertTrue(m.any())
        self.assertFalse(full.loc[m, "n_used"].str.isdigit().any())
        m = (full["n_used"] == "<20") & full["n_people"].str.isdigit()
        self.assertFalse(full.loc[m, "n_never_used"].str.isdigit().any())

    def test_formats_agree(self):
        out = os.path.join(self.tmp, "o_long58")
        self.s.run(out, os.path.join(self.tmp, "i_long58"), self.s.cn58, self.s.long_tsv)
        a, b = rd(self.out58, "diagnostic.tsv"), rd(out, "diagnostic.tsv")
        pd.testing.assert_frame_equal(a, b)

    def test_usage_not_found_lists_candidates(self):
        root = os.path.join(self.tmp, "aleix_root")
        os.makedirs(os.path.join(root, "sub"))
        for f in ("a_usage.tsv", "sub/b_usage.csv", "sub/other.tsv"):
            open(os.path.join(root, f), "w").close()
        with self.assertRaises(SystemExit) as cm:
            tg.main(["--cn-table", self.s.cn65, "--usage-root", root, "--cohort", self.s.cohort,
                     "--out-dir", os.path.join(self.tmp, "x"), "--internal-dir", os.path.join(self.tmp, "y")])
        msg = str(cm.exception)
        self.assertIn("a_usage.tsv", msg)
        self.assertIn("b_usage.csv", msg)
        with self.assertRaises(SystemExit) as cm:
            tg.main(["--cn-table", self.s.cn65, "--usage-root", os.path.join(self.tmp, "nope"),
                     "--cohort", self.s.cohort, "--out-dir", os.path.join(self.tmp, "x"),
                     "--internal-dir", os.path.join(self.tmp, "y")])
        self.assertIn("not mounted", str(cm.exception))
        # exactly one match is used
        root1 = os.path.join(self.tmp, "aleix_one")
        os.makedirs(root1)
        pd.to_pickle({"V": self.s.V}, os.path.join(root1, "vj_usage.pkl"))
        tg.main(["--cn-table", self.s.cn65, "--usage-root", root1, "--usage-glob", "*usage*",
                 "--cohort", self.s.cohort, "--relatedness-table", self.s.rel,
                 "--out-dir", os.path.join(self.tmp, "o1"), "--internal-dir", os.path.join(self.tmp, "i1"),
                 "--people-limit", "100"])
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "o1", "diagnostic.tsv")))

    def test_internal_refuses_repo(self):
        with self.assertRaises(ValueError):
            self.s.run(os.path.join(self.tmp, "o_x"), os.path.join(SRC, "x_internal_60"),
                       self.s.cn65, self.s.pkl)
        self.assertFalse(os.path.exists(os.path.join(SRC, "x_internal_60")))


if __name__ == "__main__":
    unittest.main()
