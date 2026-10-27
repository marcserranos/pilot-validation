"""Synthetic aggregate tables with the exact columns of 70 / 71 / 71b / 71c, for testing 72's figures without the VM."""
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(SRC, fn))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def build(d, n_people=600, seed=7):
    q = _load("m71", "71_call_quality_table.py")
    s = _load("m71b", "71b_call_quality_stats.py")
    rng = np.random.default_rng(seed)
    os.makedirs(os.path.join(d, "71"), exist_ok=True)
    os.makedirs(os.path.join(d, "70"), exist_ok=True)
    os.makedirs(os.path.join(d, "71c"), exist_ok=True)
    genes = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1", "MICA", "TAP1", "TAP2", "DRA"]
    # ---- records and call table
    def views(n):
        lat = rng.normal(size=n)
        y = (rng.random(n) < 1 / (1 + np.exp(-(-3.3 + 1.4 * lat)))).astype(int)
        v = pd.DataFrame({"person": [f"p{i // 3}" for i in range(n)], "gene": rng.choice(genes[:8], n),
                          "plat": rng.choice(["revio", "sequel2e"], n), "bad": y})
        for k in q.KEYS:
            v["o_" + k] = rng.normal(size=n)
        v["o_rel_cds"] = lat + rng.normal(scale=0.6, size=n)
        v["o_rel_ars"] = lat + rng.normal(scale=0.5, size=n)
        v["o_cds_p10_rel"] = lat + rng.normal(scale=0.7, size=n)
        v["o_template_distance"] = np.exp(0.5 * lat + rng.normal(size=n))
        v["o_cds_distance"] = np.round(np.exp(0.4 * lat + rng.normal(size=n)))
        v["o_n_tied"] = rng.integers(1, 5, n).astype(float)
        v["o_has_warning"] = (lat + rng.normal(size=n) > 1).astype(float)
        v["o_mhc_contigs"] = rng.integers(1, 5, n).astype(float)
        v["o_score7"] = -np.clip(np.round(5 - lat + rng.normal(size=n)), 0, 7)
        return v
    ct = views(4000)
    ct["person_id"] = ct["person"]
    V = {"trio_child": views(1200), "trio_parent": views(1200), "trio_worst": views(1200), "pair": views(5000), "twin": views(900)}
    s.analyse(V, ct, os.path.join(d, "71"), do_boot=False)
    # ---- 71e dropout tables (synthetic person-gene outcomes)
    e = _load("m71e", "71e_dropout_outcome.py")
    rows = []
    for i in range(1200):
        for g in ("A", "B", "C", "DRB1"):
            lat = rng.normal()
            two = rng.random() > 1 / (1 + np.exp(-(-2.4 + 1.3 * lat)))
            for hap in (("1", "2") if two else ("1",)):
                r = dict(person_id=f"p{i}", hap=hap, gene=g, plat="revio")
                for k in q.KEYS:
                    r["o_" + k] = rng.normal()
                r["o_rel_cds"] = lat + rng.normal(scale=0.4)
                rows.append(r)
    e.analyse_dropout(e.person_gene_table(pd.DataFrame(rows)), os.path.join(d, "71e"))
    # ---- 71 distributions
    gd = pd.DataFrame([dict(person_id=f"p{i}", gene=g, **{f: np.clip(rng.normal(0.85, 0.18), 0.05, 1.9) for f in q.DEPTH_FEATURES},
                            abs_gene=rng.normal(17, 4), abs_exons=rng.normal(17, 4), abs_cds=rng.normal(17, 4),
                            abs_ars=rng.normal(17, 4), cds_share_lt10=rng.random() * 0.2, mapq_frac=rng.random())
                       for i in range(n_people) for g in genes])
    plat = {f"p{i}": rng.choice(["revio", "sequel2e"]) for i in range(n_people)}
    q.hist_table(gd, plat).to_csv(os.path.join(d, "71", "dist_gene_depth_hist.tsv"), sep="\t", index=False)
    q.quantile_table(gd, plat).to_csv(os.path.join(d, "71", "dist_gene_depth_quantiles.tsv"), sep="\t", index=False)
    q.gene_corr(gd).to_csv(os.path.join(d, "71", "corr_gene_depth.tsv"), sep="\t", index=False)
    cc = q.metric_corr(ct.assign(**{"o_" + k: ct["o_" + k] for k in q.KEYS}))
    cc.to_csv(os.path.join(d, "71", "corr_metric_metric.tsv"), sep="\t", index=False)
    per58 = pd.concat([pd.DataFrame({"person_id": [f"p{i}"] * 2, "region": ["window_MHC", "control_chr2"], "mapq": [0, 0],
                                     "mean_depth": [rng.normal(15, 3), 0], "control_depth": [rng.normal(19, 3)] * 2,
                                     "rel_depth": [rng.normal(0.8, 0.07), 1.0]}) for i in range(n_people)])
    plat_of = plat
    h, g2 = q.person_mhc_tables(per58, plat_of, {p: rng.choice(["AFR", "EUR", "EAS"]) for p in plat})
    h.to_csv(os.path.join(d, "71", "dist_person_mhc_hist.tsv"), sep="\t", index=False)
    g2.to_csv(os.path.join(d, "71", "dist_mhc_vs_genome_grid.tsv"), sep="\t", index=False)
    # ---- 70 features and profile
    rows, prof = [], []
    for k, g in enumerate(["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]):
        s0 = 1000 + k * 40000
        L = 8000 + 1000 * k
        rows.append(dict(chrom="chr6", start=s0, end=s0 + L, name=f"{g}|gene", gene=g, kind="gene", exon_number=0, cds_ord=0,
                         ars=False, length=L, strand="+" if k % 2 == 0 else "-", transcript_id=""))
        pos = s0
        ords = 0
        for e in range(1, 6):
            es, ee = s0 + e * 1400, s0 + e * 1400 + 250 + 20 * e
            rows.append(dict(chrom="chr6", start=es, end=ee, name=f"{g}|exons{e}", gene=g, kind="exons", exon_number=e, cds_ord=0,
                             ars=False, length=ee - es, strand="+", transcript_id="T"))
            rows.append(dict(chrom="chr6", start=es + 10, end=ee - 10, name=f"{g}|cds{e}", gene=g, kind="cds", exon_number=e, cds_ord=0,
                             ars=False, length=ee - es - 20, strand="+", transcript_id="T"))
            rows.append(dict(chrom="chr6", start=es + 10, end=ee - 10, name=f"{g}|ex{e}", gene=g, kind="ex", exon_number=e, cds_ord=e,
                             ars=e in (2, 3), length=ee - es - 20, strand="+", transcript_id="T"))
        bins = np.arange(s0 - 500, s0 + L + 500, 10)
        base = 0.85 + 0.1 * np.sin(bins / 700.0)
        for b, mu in zip(bins, base):
            prof.append(dict(gene=g, bin_start=b, bin_end=b + 10, n_people=300, rel_p10=mu - 0.25, rel_p25=mu - 0.1,
                             rel_p50=mu, rel_p75=mu + 0.1, rel_p90=mu + 0.25))
    pd.DataFrame(rows).to_csv(os.path.join(d, "70", "features.tsv"), sep="\t", index=False)
    pd.DataFrame(prof).to_csv(os.path.join(d, "70", "profile_by_gene.tsv"), sep="\t", index=False)
    ex = []
    for g in ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]:
        for e in range(1, 6):
            mu = rng.normal(0.85, 0.05)
            ex.append(dict(gene=g, exon_number=e, cds_ord=e, ars=e in (2, 3), length=300, n_people=500,
                           p05=mu - 0.3, p10=mu - 0.22, p25=mu - 0.1, p50=mu, p75=mu + 0.1, p90=mu + 0.22, p95=mu + 0.3))
    pd.DataFrame(ex).to_csv(os.path.join(d, "71", "dist_exon_depth.tsv"), sep="\t", index=False)
    # ---- 71c matrices
    c = _load("m71c", "71c_contig_truth_matrix.py")
    og = []
    for ch in range(60):
        for hap in ("hap1", "hap2"):
            par = rng.choice(["pat", "mat"])
            ctg = 0
            for g in c.GENE_ORDER:
                if g == "C4A":
                    continue
                if rng.random() < 0.15:
                    ctg += 1
                    par = rng.choice(["pat", "mat"])
                og.append(dict(child=f"c{ch}", hap=hap, gene=g, contig=f"ctg{ch}{hap}{ctg}", origin=par if rng.random() > 0.15 else "amb"))
    o = pd.DataFrame(og)
    c.truth_table(o, strata={f"c{i}": "revio" for i in range(60)}).to_csv(os.path.join(d, "71c", "contig_matrix_trio.tsv"), sep="\t", index=False)
    return d
