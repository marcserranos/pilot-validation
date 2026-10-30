import importlib.util
import json
import pathlib
import sys

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("m130", pathlib.Path(__file__).parents[1] / "130_call_confidence.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_wilson_and_auc():
    lo, hi = m.wilson(3, 2106)
    assert 0.0004 < lo < 0.001 and 0.003 < hi < 0.005
    assert m.auc([1, 1, 0, 0], [3, 2, 1, 0]) == 1.0
    assert m.auc([1, 0, 1, 0], [1, 1, 1, 1]) == 0.5


def test_curve_perfect_metric_captures_errors_first():
    y = np.array([1] * 5 + [0] * 95)
    fr, cap = m.curve(y, y.astype(float))
    assert cap[fr == 0.05][0] == 1.0 and cap[0] == 0.0


def test_end_to_end(tmp_path):
    rng = np.random.default_rng(1)
    n = 400
    gene = rng.choice(["A", "B", "DRB1", "G"], n)
    ct = pd.DataFrame(dict(person_id=np.arange(n) % 100, hap=rng.choice(["hap1", "hap2"], n), gene=gene,
                           x_template_distance=rng.choice([0, 0, 0, 1, 4], n), x_cds_share_lt10=rng.random(n) * (rng.random(n) > .5),
                           x_n_aa_changes=rng.integers(0, 3, n)))
    ct.to_csv(tmp_path / "ct.tsv", sep="\t", index=False)
    def rec(k, extra):
        d = dict(person=np.arange(k), other=np.arange(k) + 1000, gene=rng.choice(["A", "B", "DRB1", "G"], k), hap=rng.choice(["hap1", "hap2"], k),
                 family=np.arange(k) % 30)
        for p in extra:
            d[f"{p}_template_distance"] = rng.choice([0, 0, 1, 4], k).astype(float)
            d[f"{p}_cds_share_lt10"] = rng.random(k)
            d[f"{p}_n_aa_changes"] = rng.integers(0, 3, k).astype(float)
            d[f"{p}_rel_cds"] = -rng.random(k)
        return pd.DataFrame(d)
    T = rec(300, ["c", "w"]); T["cls"] = rng.choice(["match", "neither", "other_parent_only"], 300, p=[.96, .02, .02])
    Q = rec(300, ["w"]); Q["bad"] = rng.random(300) < .03
    T.to_csv(tmp_path / "t.tsv", sep="\t", index=False)
    Q.to_csv(tmp_path / "q.tsv", sep="\t", index=False)
    json.dump(dict(alleles=[dict(id=0, srl=4, xu=1, xuc=0)]), open(tmp_path / "l.json", "w"))
    pd.DataFrame(dict(allele_id=[0], gene=["A"], person=[1], hap=[1], template=["x"])).to_csv(tmp_path / "c.tsv", sep="\t", index=False)
    sys.argv = ["x", "--out", str(tmp_path / "o"), "--call-table", str(tmp_path / "ct.tsv"), "--trio", str(tmp_path / "t.tsv"),
                "--pair", str(tmp_path / "q.tsv"), "--ladder", str(tmp_path / "l.json"), "--carriers", str(tmp_path / "c.tsv")]
    m.main()
    out = tmp_path / "o"
    for name in ["cohort_sizes", "error_by_cohort", "error_by_metric_bin", "prune_curves", "prune_policies", "predictor_auc"]:
        assert (out / f"{name}.private.tsv").exists() and (out / f"{name}.pub.tsv").exists()
    e = pd.read_csv(out / "error_by_cohort.private.tsv", sep="\t")
    assert e[(e.record_set == "trio") & (e.scope == "all_genes") & (e.cohort == "all")].n.iloc[0] == 300
    assert "rel_cds" in set(pd.read_csv(out / "predictor_auc.private.tsv", sep="\t").metric)
    pub = pd.read_csv(out / "error_by_cohort.pub.tsv", sep="\t", dtype=str)
    assert "<20" in set(pub.errors)
