import importlib.util
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("l118", os.path.join(HERE, "..", "118_novel_ladder.py"))
L = importlib.util.module_from_spec(spec)
spec.loader.exec_module(L)


class T(unittest.TestCase):
    def test_lod_het_supported_vs_absent(self):
        self.assertGreater(L.site_lod(15, 15, 0, True), 10)
        self.assertLess(L.site_lod(0, 30, 0, True), -5)
        self.assertIsNone(L.site_lod(0, 0, 0, True))

    def test_lod_hom(self):
        self.assertGreater(L.site_lod(30, 0, 0, False), 40)
        self.assertLess(L.site_lod(10, 20, 0, False), 0)          # a minority of reads on a "hom" site is not support

    def test_case_lod_weakest_site(self):
        s = [{"obs": 15, "alt": 15, "oth": 0, "het": True, "verdict": "supported"}, {"obs": 0, "alt": 20, "oth": 0, "het": True, "verdict": "x"},
             {"verdict": "error: boom"}]
        self.assertLess(L.case_lod(s), 0)
        self.assertIsNone(L.case_lod([{"verdict": "error: x"}]))
        self.assertEqual(L.case_lod([{"obs": 300, "alt": 0, "oth": 0, "het": False, "verdict": "s"}]), L.LOD_CAP)

    def test_levels(self):
        self.assertEqual([L.sr_level(x) for x in (None, -1, 1, 4, 9)], [0, 1, 2, 3, 4])

    def test_effect_and_ref(self):
        self.assertEqual([L.effect_of(p) for p in ("synonymous", "A57D", "A57D;W10*")], [0, 1, 2])
        self.assertEqual(L.ref_flags("B*56:01:01:05S", 1000, 1000), (False, "suffix"))
        self.assertEqual(L.ref_flags("C*07:402", 900, 1000), (False, "truncated"))
        self.assertEqual(L.ref_flags("C*07:01", 1000, 1000), (True, ""))

    def test_bucket(self):
        b = dict(rc=True, near=True, srl=4, srn=2, xu=3, tr=True, mp=False)
        self.assertEqual(L.bucket(b), 0)
        self.assertEqual(L.bucket(dict(b, tr=False)), 1)
        self.assertEqual(L.bucket(dict(b, srl=0)), 2)
        self.assertEqual(L.bucket(dict(b, srl=1)), 3)
        self.assertEqual(L.bucket(dict(b, rc=False)), 3)

    def test_cube_masking(self):
        a = [dict(gc=0, eff=1, ipd=0, rc=True, near=True, srl=4, srn=1, tr=False, mp=False, hm=False, xu=n, xa=n + 1) for n in (1, 1, 5, 25)]
        priv = L.build_cube(a, False)
        pub = L.build_cube(a, True)
        self.assertEqual(sum(c["w"] for c in priv), 4)
        self.assertTrue(all(c["xu"] is None for c in pub if c["xuc"] < 2))
        self.assertEqual([c["xu"] for c in pub if c["xuc"] == 2], [25])
        self.assertEqual({c["xuc"] for c in pub}, {0, 1, 2})

    def test_saturation_monotone_and_complete(self):
        persons = set(range(200))
        car = {i: {i, (i * 7) % 200} for i in range(50)}
        s = L.saturation(car, persons, 20, L.sat_points(200))
        self.assertEqual(s["med"], sorted(s["med"]))
        self.assertEqual(s["med"][-1], 50)
        self.assertTrue(all(lo <= m <= hi for lo, m, hi in zip(s["lo"], s["med"], s["hi"])))

    def test_xclass(self):
        self.assertEqual([L.xclass(n) for n in (1, 2, 19, 20)], [0, 1, 1, 2])


if __name__ == "__main__":
    unittest.main()
