import importlib.util
import pathlib

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("m132", pathlib.Path(__file__).parents[1] / "132_missingness.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_end_to_end(tmp_path):
    rng = np.random.default_rng(0)
    people = [str(i) for i in range(300)]
    rows = [(p, g, h, rng.choice(["called", "missing"], p=[.95, .05]), rng.choice(["inside_contig", "between_contigs", "region_end"]),
             10, int(rng.integers(1, 8))) for p in people for g in m.ANCHORS for h in "12"]
    pd.DataFrame(rows, columns=["person", "gene", "hap", "status", "hap_class", "hap_n_genes", "hap_n_contigs"]).to_csv(tmp_path / "s.tsv", sep="\t", index=False)
    pd.DataFrame(dict(person_id=people, ancestry_pred=rng.choice(["eur", "afr", "amr"], 300), platform=rng.choice(["a", "b"], 300))).to_csv(tmp_path / "c.tsv", sep="\t", index=False)
    S = m.load(tmp_path / "s.tsv", tmp_path / "c.tsv")
    assert set(S.anc) == {"EUR", "AFR", "AMR"}
    man = m.analyse(S, str(tmp_path))
    assert man["slots"] == 300 * 12 * 2
    o = pd.read_csv(tmp_path / "m_overall.private.tsv", sep="\t")
    assert o.loc[o.scope == "classical8", "n"].iloc[0] == 300 * 8 * 2
    g = pd.read_csv(tmp_path / "m_by_gene.private.tsv", sep="\t")
    assert np.allclose(g.share_0 + g.share_1 + g.share_2, 1) and list(g.gene) == m.ANCHORS
    assert (tmp_path / "m_by_ancestry.pub.tsv").exists()
