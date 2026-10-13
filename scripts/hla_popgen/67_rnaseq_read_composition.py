#!/usr/bin/env python3
"""S05 67: what the AoU blood RNA-seq reads are made of (rRNA, mitochondrial, globin, exonic, intronic, intergenic).

Marc (2026-10-13): "what percent of the RNA reads are rRNA, tRNA, or normal RNA?" AoU ships two per-sample QC tables
next to the STAR BAMs (v9/multiomics/rnaseq/): rnaseq_metadata.tsv (DRAGEN/Picard-style: ribosomal_bases_pct,
alignment_rate_pct, ...) and rnaseqc2/*.metrics.txt.gz (RNA-SeQC 2: Exonic / Intronic / Intergenic / rRNA Rate,
Non-Globin Reads, Mitochondrial reads if present, ...). Both are small. This script reads every numeric column
of both, turns count columns that have a matching total into fractions, and writes cohort quantiles only.

tRNA is not a column in either table: standard poly-A or rRNA-depleted libraries are size-selected (~200+ nt
inserts) and mature tRNAs (~76 nt, heavily modified, which stalls reverse transcription) are mostly lost, so their
share is not measured here (see the report's RNA-seq module).

PUBLIC (aggregate over ~8,980 samples, no ids): ~/s05/results/67/composition_quantiles.tsv (metric, source, n,
mean, p5, p25, median, p75, p95), composition_by_ancestry.tsv (median per ancestry, groups >= 20),
fig_67_composition.png, columns_seen.tsv (the exact column names, for the report's data trail).
    cd ~/repos/<checkout>/scripts/hla_popgen && $PY 67_rnaseq_read_composition.py
Python 3.8 compatible.
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _s05_common as C  # noqa: E402

ROOT = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/v9/multiomics/rnaseq/")
OUT = os.path.expanduser("~/s05/results/67/")
# count column -> its denominator (fractions are what the report needs)
RATIO = {"Non-Globin Reads": "Total Reads", "Mapped Reads": "Total Reads", "Mitochondrial Reads": "Mapped Reads",
         "rRNA Reads": "Mapped Reads", "Exonic Reads": "Mapped Reads", "Intronic Reads": "Mapped Reads",
         "Intergenic Reads": "Mapped Reads", "Chimeric Reads": "Mapped Reads"}
HEADLINE = ["rRNA Rate", "ribosomal_bases_pct", "Mitochondrial Reads / Mapped Reads", "Non-Globin Reads / Total Reads",
            "Exonic Rate", "Intronic Rate", "Intergenic Rate", "Mapping Rate", "alignment_rate_pct",
            "mrna_bases_pct", "intronic_bases_pct", "intergenic_bases_pct", "coding_bases_pct", "utr_bases_pct"]


def quantiles(x):
    x = pd.to_numeric(x, errors="coerce").dropna()
    if len(x) < 20:
        return None
    q = np.quantile(x, [0.05, 0.25, 0.5, 0.75, 0.95])
    return dict(n=int(round(len(x) / 10.0)) * 10, mean=x.mean(), p5=q[0], p25=q[1], median=q[2], p75=q[3], p95=q[4])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=ROOT)
    ap.add_argument("--out-dir", default=OUT)
    a = ap.parse_args(argv)
    out = a.out_dir.rstrip("/") + "/"
    os.makedirs(out, exist_ok=True)
    meta = pd.read_csv(a.root + "rnaseq_metadata.tsv", sep="\t", dtype=str)
    mpath = sorted(glob.glob(a.root + "rnaseqc2/*.metrics.txt.gz"))
    if not mpath:
        raise SystemExit(f"FATAL: no rnaseqc2/*.metrics.txt.gz under {a.root} (mount up?)")
    seq = pd.read_csv(mpath[-1], sep="\t")
    seq.columns = [c.strip() for c in seq.columns]
    seen = [dict(source="rnaseq_metadata.tsv", column=c) for c in meta.columns] + \
           [dict(source=os.path.basename(mpath[-1]), column=c) for c in seq.columns]
    C.write_tsv(pd.DataFrame(seen), out + "columns_seen.tsv")
    for num, den in RATIO.items():
        if num in seq.columns and den in seq.columns:
            seq[f"{num} / {den}"] = pd.to_numeric(seq[num], errors="coerce") / pd.to_numeric(seq[den], errors="coerce")
    rows = []
    for src, df in (("rnaseq_metadata", meta), ("rnaseqc2", seq)):
        for c in df.columns:
            if c in ("sample_id", "research_id", "sampleid", "person_id"):
                continue
            v = pd.to_numeric(df[c], errors="coerce")
            if v.notna().mean() < 0.9:
                continue
            q = quantiles(v)
            if q:
                rows.append(dict(metric=c, source=src, headline=c in HEADLINE, **{k: (round(x, 5) if isinstance(x, float)
                                                                                 else x) for k, x in q.items()}))
    t = pd.DataFrame(rows).sort_values(["headline", "source", "metric"], ascending=[False, True, True])
    C.write_tsv(t, out + "composition_quantiles.tsv")
    # by ancestry (medians; groups of >= 20 samples only)
    coh = C.load_cohort(C.DEFAULT_COHORT)[["person_id", "ancestry_pred"]]
    # RNA-SeQC sample_id == research_id (checked 2026-07-25); the metadata columns join on research_id
    sidc = "sample_id" if "sample_id" in seq.columns else seq.columns[0]
    j = seq.copy()
    j[sidc] = j[sidc].astype(str)
    if "research_id" in meta.columns:
        j = j.merge(meta.rename(columns={"research_id": sidc}), on=sidc, how="left", suffixes=("", "_meta"))
    j = j.merge(coh.rename(columns={"person_id": sidc}), on=sidc, how="left")
    keep = [c for c in HEADLINE if c in j.columns]
    by = []
    for anc, g in j.groupby(j["ancestry_pred"].fillna("NA")):
        if len(g) < 20 or anc in ("NA", "NAN"):
            continue
        by.append(dict(ancestry=anc, n=int(round(len(g) / 10.0)) * 10,
                       **{c: round(float(pd.to_numeric(g[c], errors="coerce").median()), 5) for c in keep}))
    C.write_tsv(pd.DataFrame(by), out + "composition_by_ancestry.tsv")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    h = t[t["headline"]].copy()
    fig, ax = plt.subplots(figsize=(8, 0.42 * len(h) + 1.2))
    for i, r in enumerate(h.itertuples(index=False)):
        sc = 100 if r.median <= 1.0 and not r.metric.endswith("_pct") else 1
        ax.plot([r.p5 * sc, r.p95 * sc], [i, i], color="0.6", lw=1)
        ax.plot([r.p25 * sc, r.p75 * sc], [i, i], color="C0", lw=5, solid_capstyle="butt")
        ax.plot(r.median * sc, i, "o", color="k", ms=4)
        ax.text(101, i, f"{r.median * sc:.1f}%", va="center", fontsize=8)
    ax.set_yticks(range(len(h)))
    ax.set_yticklabels([f"{m} ({s})" for m, s in zip(h["metric"], h["source"])], fontsize=8)
    ax.set_xlim(0, 112)
    ax.set_xlabel("% (bar = IQR, line = 5th-95th percentile, dot = median)")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out + "fig_67_composition.png", dpi=200)
    C.log("DONE 67")
    return 0


if __name__ == "__main__":
    sys.exit(main())
