#!/usr/bin/env python3
"""S05 68: what is in the All of Us blood RNA-seq libraries? (VM; aggregate output only.)

Answers "what share of the reads are ribosomal, globin, mitochondrial, intronic, ..." from the files
All of Us already ships, without touching BAMs:
  1. rnaseq_metadata.tsv (DRAGEN QC per sample): mrna_bases_pct, ribosomal_bases_pct, alignment rate,
     read pairs, insert size.
  2. rnaseqc2/*.metrics.txt.gz (RNA-SeQC 2 per sample): every numeric column (Exonic / Intronic /
     Intergenic / rRNA Rate, mapping and duplicate rates, Total Reads, Non-Globin Reads -> globin share).
  3. rnaseqc2/*.gene_tpm.gct.gz streamed in chunks: the share of each sample's TPM that falls on
     globin genes (HBA1/2, HBB, HBD, ...), mitochondrial rRNA (MT-RNR1/2), mitochondrial tRNA (MT-T*),
     other mitochondrial genes, nuclear rRNA (RNA5-8SN*, RNA18SN*, RNA28SN*, RNA45SN*, RNA5S*), and the
     HLA and KIR genes; plus the per-sample TPM total (should be 1e6).
Output (PUBLIC, cohort-level distribution only: n, mean, p1, p5, p25, p50, p75, p95, p99 per metric):
  ~/s05/results/68/68_rnaseq_composition.tsv
Nothing per sample leaves the VM. n is the number of samples (thousands), so no small cells arise.

VM:
    cd ~/repos/pilot-validation-main/scripts/hla_popgen && $PY 68_rnaseq_composition.py
"""
import argparse
import glob
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd

RNA = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/v9/multiomics/rnaseq")
OUT = os.path.expanduser("~/s05/results/68/")
QS = [1, 5, 25, 50, 75, 95, 99]
GENE_SETS = {
    "globin": re.compile(r"^(HBA1|HBA2|HBB|HBD|HBG1|HBG2|HBQ1|HBM|HBZ|HBE1)$"),
    "mt_rRNA": re.compile(r"^MT-RNR[12]$"),
    "mt_tRNA": re.compile(r"^MT-T[A-Z0-9]+$"),
    "mt_other": re.compile(r"^MT-(?!RNR|T[A-Z0-9]+$)"),
    "nuclear_rRNA": re.compile(r"^(RNA5-8SN\d*|RNA18SN\d*|RNA28SN\d*|RNA45SN\d*|RNA5S\d+)$"),
    "HLA_classical": re.compile(r"^HLA-(A|B|C|DRA|DRB1|DRB3|DRB4|DRB5|DQA1|DQB1|DPA1|DPB1)$"),
    "KIR": re.compile(r"^KIR\d"),
}


def log(msg):
    print(f"[68] {msg}", file=sys.stderr, flush=True)


def describe(name, source, x, unit):
    x = pd.to_numeric(pd.Series(x), errors="coerce").dropna()
    row = {"source": source, "metric": name, "unit": unit, "n_samples": len(x),
           "mean": f"{x.mean():.4g}" if len(x) else ""}
    for q in QS:
        row[f"p{q}"] = f"{np.percentile(x, q):.4g}" if len(x) else ""
    return row


def tpm_shares(path, chunksize=2000):
    """Stream the GCT: per sample, TPM summed over each gene set and over all genes."""
    with gzip.open(path, "rt") as f:
        f.readline()
        f.readline()
        header = f.readline().rstrip("\n").split("\t")
    samples = header[2:]
    tot = np.zeros(len(samples))
    sets = {k: np.zeros(len(samples)) for k in GENE_SETS}
    n_genes = 0
    for ch in pd.read_csv(path, sep="\t", skiprows=2, chunksize=chunksize, compression="gzip"):
        desc = ch.iloc[:, 1].astype(str)
        vals = ch.iloc[:, 2:].to_numpy(dtype=float)
        tot += vals.sum(axis=0)
        n_genes += len(ch)
        for k, rx in GENE_SETS.items():
            m = desc.str.match(rx).to_numpy()
            if m.any():
                sets[k] += vals[m].sum(axis=0)
        log(f"tpm rows {n_genes}")
    return samples, tot, sets, n_genes


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rna-root", default=RNA)
    ap.add_argument("--out-dir", default=OUT)
    ap.add_argument("--skip-tpm", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    rows = []
    meta = pd.read_csv(os.path.join(a.rna_root, "rnaseq_metadata.tsv"), sep="\t")
    for c in meta.columns:
        if pd.to_numeric(meta[c], errors="coerce").notna().mean() > 0.9 and c.lower() not in ("research_id", "sampleid", "sample"):
            rows.append(describe(c, "rnaseq_metadata.tsv (DRAGEN)", meta[c], "pct" if c.endswith("_pct") else "value"))
    mp = sorted(glob.glob(os.path.join(a.rna_root, "rnaseqc2", "*.metrics.txt.gz")))
    if mp:
        met = pd.read_csv(mp[-1], sep="\t")
        for c in met.columns:
            v = pd.to_numeric(met[c], errors="coerce")
            if v.notna().mean() > 0.9 and c.lower() not in ("sample", "sample id", "research_id"):
                rows.append(describe(c, "RNA-SeQC 2 metrics", v, "fraction" if "Rate" in c else "value"))
        if {"Non-Globin Reads", "Total Reads"} <= set(met.columns):
            g = 1 - pd.to_numeric(met["Non-Globin Reads"], errors="coerce") / pd.to_numeric(met["Total Reads"], errors="coerce")
            rows.append(describe("globin read share (1 - Non-Globin / Total)", "RNA-SeQC 2 metrics", 100 * g, "pct"))
    tp = sorted(glob.glob(os.path.join(a.rna_root, "rnaseqc2", "*.gene_tpm.gct.gz")))
    if tp and not a.skip_tpm:
        samples, tot, sets, n_genes = tpm_shares(tp[-1])
        rows.append(describe("TPM total per sample", "gene_tpm.gct", tot, "TPM"))
        for k, v in sets.items():
            rows.append(describe(f"share of TPM: {k}", "gene_tpm.gct", 100 * v / tot, "pct"))
        rows.append(describe("genes in matrix", "gene_tpm.gct", [n_genes], "count"))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(a.out_dir, "68_rnaseq_composition.tsv"), sep="\t", index=False)
    log(f"wrote {len(out)} rows to {a.out_dir}68_rnaseq_composition.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
