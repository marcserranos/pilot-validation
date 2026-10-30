import importlib.util
import pathlib
import types

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("m134", pathlib.Path(__file__).parents[1] / "134_kir_confidence.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_end_to_end(tmp_path):
    rng = np.random.default_rng(3)
    genes = ["KIR2DL1", "KIR3DL1", "KIR2DS2"]
    people = [str(i) for i in range(400)]
    K = pd.DataFrame([(p, g, 1, 1, rng.choice(["known", "novel_protein", "novel_cds_synonymous"], p=[.8, .05, .15]), "KIR*001", float(rng.choice([0, 0, 1, 4])), 0.0)
                      for p in people for g in genes], columns=["person_id", "gene", "hap", "copy_index", "novelty_tier", "label3", "template_distance", "cds_distance"])
    K.loc[K.novelty_tier == "novel_protein", "label3"] = np.nan
    H = pd.DataFrame([(p, h, 3, 3, 1, rng.choice(["complete", "multi_contig", "empty"])) for p in people for h in (1, 2)], columns=["person_id", "hap", "n_calls", "n_genes", "n_contigs", "hap_class"])
    D = pd.DataFrame([(p, g, True, 1.0, rng.random(), rng.random() * 0.2, 0.9, 0.8, 0.95, 0.0, 0.0) for p in people for g in genes[:2]],
                     columns=["person_id", "gene", "depth_defined", "rel_gene", "rel_cds", "cds_p10_rel", "cds_min_rel", "cds_share_lt10", "mapq_frac", "rel_cds_q0", "dev_cds_q0"])
    n = 600
    T = pd.DataFrame(dict(child=rng.choice(people, n), tx="t", other="o", gene=rng.choice(genes, n), locus="x", mode="gene", hap=1, copy_index=1, origin="pat", strong=True,
                          stream=rng.choice(["named", "novel"], n, p=[.9, .1]), tx_complete=rng.choice([True, False], n, p=[.8, .2]), tx_unnamed=rng.choice([1, 0], n, p=[.1, .9]),
                          tx_n_calls=2, template_distance=rng.choice([0., 0., 1., 4.], n), cds_distance=0.0, cls3=rng.choice(["match", "neither", "only_other_parent"], n, p=[.95, .03, .02])))
    P = pd.DataFrame(dict(person=rng.choice(people, n), other=rng.choice(people, n), locus=rng.choice(genes, n), mode="gene", kind="pair",
                          cls=rng.choice(["share", "no_share"], n, p=[.9, .1]), both_complete=rng.choice([True, False], n, p=[.8, .2])))
    coh = dict(anc={p: rng.choice(["AFR", "EUR"]) for p in people}, plat={p: rng.choice(["revio", "sequel2e"]) for p in people})
    man = m.analyse(None, str(tmp_path), K, D, T, P, H, coh)
    assert man["calls"] == 1200 and man["trio_eval"] > 0 and man["pair_eval"] > 0
    for name in ["k_cohort_sizes", "k_error_by_cohort", "k_error_by_metric_bin", "k_prune_curves", "k_prune_policies", "k_predictor_auc", "k_hap_completeness", "k_label_instability"]:
        assert (tmp_path / f"{name}.private.tsv").exists(), name
    e = pd.read_csv(tmp_path / "k_error_by_cohort.private.tsv", sep="\t")
    assert set(e.cohort) == {"C1", "C2", "all"}
    a = pd.read_csv(tmp_path / "k_predictor_auc.private.tsv", sep="\t")
    assert "template_distance" in set(a.metric)
