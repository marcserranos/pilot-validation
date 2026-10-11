#!/usr/bin/env python3
"""S05: orthogonal validation of structural (presence/absence) calls with blood RNA expression.

If an assembly-based presence/absence call is right, people carrying 0 copies of a gene should
show ~no expression of it, 1 copy intermediate, 2 copies the most (dosage). Uses the AoU v9
RNA-SeQC gene TPM matrix (GRCh38 primary annotation), so only genes present on the primary
assembly can be tested:
  - KIR2DL4 (framework gene; absent on contracted KIR haplotypes, see 51)
  - KIR3DL1 (absent on haplotypes carrying 3DS1 instead or a contracted haplotype; 3DS1 reads may
    cross-map, so a non-zero floor is expected)
  - KIR2DL1, KIR2DL3, KIR2DS4 (variable-content genes on the reference cA01~tA01 haplotype)
  - HLA-DRB5 (present only on DR51 haplotypes, DRB1*15/*16)
Copy number per person = number of haplotypes on which the gene is called (0/1/2). KIR: people
whose two haplotypes are both "complete" (51's definition) in the strict set; all people in the
loose set. DRB5: haplotypes where DRB1 and DQA1 are on one contig (absence trustworthy).

Output (PUBLIC, aggregate): structure_expression.tsv -- per gene x set x copy number: n people
(masked <20), median / IQR TPM, % with TPM < 1; Spearman rho(copies, TPM).

    python3 55_structure_expression_validation.py --out-dir ~/s05/results/55
"""
import argparse
import gzip
import importlib.util
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import _s05_common as C

RNA = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/v9/multiomics/rnaseq")
DEFAULT_TPM = os.path.join(RNA, "rnaseqc2", "aou_rnaseq_20260413.gene_tpm.gct.gz")
DEFAULT_MANIFEST = os.path.join(RNA, "manifest.tsv")
KIR_TEST = ["KIR2DL4", "KIR3DL1", "KIR2DL1", "KIR2DL3", "KIR2DS4"]
HLA_TEST = ["DRB5"]


def load_51():
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "51_kir_structure_qc.py")
    spec = importlib.util.spec_from_file_location("k51", p)
    m = importlib.util.module_from_spec(spec)
    sys.modules["k51"] = m
    spec.loader.exec_module(m)
    return m


def read_tpm(path, symbols):
    """Stream the GCT, keep rows whose Description is in symbols. Returns DataFrame symbol x sample."""
    with gzip.open(path, "rt") as f:
        f.readline()
        f.readline()
        header = f.readline().rstrip("\n").split("\t")
        rows = {}
        for line in f:
            i = line.find("\t")
            j = line.find("\t", i + 1)
            sym = line[i + 1:j]
            if sym in symbols:
                rows[sym] = np.array(line.rstrip("\n").split("\t")[2:], dtype=float)
    return pd.DataFrame(rows, index=header[2:]).T


def summarize(cn, tpm_by_pid, gene, setname, mode):
    rows = []
    df = pd.DataFrame({"cn": cn}).join(tpm_by_pid.rename("tpm"), how="inner").dropna()
    for k, g in df.groupby("cn"):
        q = g["tpm"].quantile([.25, .5, .75]).values if len(g) else [np.nan] * 3
        r = C.rate_cells(int((g["tpm"] < 1).sum()), len(g), mode, 1)
        small = mode == "public" and 0 < len(g) < 20
        rows.append({"gene": gene, "set": setname, "copies": int(k), "n_people": C.mask_n(len(g), mode),
                     "tpm_q25": "" if small else round(q[0], 2), "tpm_median": "" if small else round(q[1], 2),
                     "tpm_q75": "" if small else round(q[2], 2), "pct_tpm_lt1": r["pct"]})
    if df["cn"].nunique() > 1:
        rho = spearmanr(df["cn"], df["tpm"]).correlation
        rows.append({"gene": gene, "set": setname, "copies": "rho", "n_people": C.mask_n(len(df), mode),
                     "tpm_median": round(rho, 3)})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tpm", default=DEFAULT_TPM)
    ap.add_argument("--manifest", default=DEFAULT_MANIFEST)
    ap.add_argument("--table1", default=C.DEFAULT_TABLE1)
    ap.add_argument("--kir-outroot", default=C.DEFAULT_KIR_OUTROOT)
    ap.add_argument("--out-dir", default=os.path.expanduser("~/s05/results/55"))
    ap.add_argument("--disclosure", default="public", choices=["public", "internal"])
    ap.add_argument("--workers", type=int, default=2)
    args = ap.parse_args()
    mode = args.disclosure

    man = pd.read_csv(args.manifest, sep="\t", dtype=str).drop_duplicates("sampleid")
    pid_of_sample = dict(zip(man["sampleid"], man["research_id"]))
    syms = set(KIR_TEST) | {f"HLA-{g}" for g in HLA_TEST}
    tpm = read_tpm(args.tpm, syms)
    tpm.columns = [pid_of_sample.get(c, c) for c in tpm.columns]
    tpm = tpm.loc[:, ~tpm.columns.duplicated()]
    C.log(f"TPM rows found: {sorted(tpm.index)}; samples: {tpm.shape[1]}")

    rows = []
    k51 = load_51()
    root = os.path.expanduser(args.kir_outroot)
    pids = sorted(d for d in os.listdir(root) if d.isdigit() and d in set(tpm.columns))
    with ProcessPoolExecutor(args.workers) as ex:
        recs = [r for rs in ex.map(k51.person_records, [(p, root) for p in pids], chunksize=64) for r in rs]
    by_pid = {}
    for pid, hap, calls in recs:
        by_pid.setdefault(pid, []).append((k51.classify_hap(calls) == "complete", {g for _, g, _, _ in calls}))
    for gene in KIR_TEST:
        if gene not in tpm.index:
            continue
        for setname, strict in [("kir_both_haps_complete", True), ("kir_all", False)]:
            cn = {p: sum(gene in g for _, g in v) for p, v in by_pid.items()
                  if len(v) == 2 and (not strict or all(c for c, _ in v))}
            rows += summarize(pd.Series(cn), tpm.loc[gene], gene, setname, mode)

    t1 = C.load_table1(args.table1)
    t1 = t1[t1["person_id"].isin(set(tpm.columns))]
    for gene in HLA_TEST:
        sym = f"HLA-{gene}"
        if sym not in tpm.index:
            continue
        cn = {}
        for pid, g in t1.groupby("person_id"):
            ok = True
            k = 0
            for hap, gh in g.groupby("hap"):
                ctg = dict(zip(gh["gene"], gh["contig"]))
                if not ("DRB1" in ctg and "DQA1" in ctg and ctg["DRB1"] == ctg["DQA1"]):
                    ok = False
                k += int(gene in ctg and ctg[gene] == ctg.get("DRB1"))
            if ok and g["hap"].nunique() == 2:
                cn[pid] = k
        rows += summarize(pd.Series(cn), tpm.loc[sym], sym, "drb1_dqa1_same_contig_both_haps", mode)
    C.write_tsv(pd.DataFrame(rows), os.path.join(os.path.expanduser(args.out_dir), "structure_expression.tsv"))
    C.log("DONE 55")


if __name__ == "__main__":
    main()
