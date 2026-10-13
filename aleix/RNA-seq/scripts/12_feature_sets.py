#!/usr/bin/env python3
"""12 -- W3: person-level feature sets for disease models (VM-local; aggregate log only).

One matrix per feature family, rows = analysis_main people (report 09), so W4 and W5 use
identical inputs:

  base       age, age^2, sex (Female / Unknown vs Male), AoU PC1-10, log10 TRB reads
             -- the covariates every disease model is judged against
  diversity  fixed-depth (823 reads) richness, Shannon, inverse Simpson, public fraction,
             and top-10 clonal share (report 09)
  vj         CLR usage of the common TRBV and TRBJ genes (report 10's rules)
  sceptr     mean SCEPTR vector over the person's whole repertoire (full-cache embedding)
             + CLR TRBV usage -- report 06's best person representation
  kmer       CDR3 3-mer composition (unique clonotypes, frequency per person), reduced to
             100 dimensions by truncated SVD fitted on all people

Germline note (report 10): TRBV4-3 and TRBV12-4 usage partly encode a germline genotype
(27% / 16% of people never use them). They stay in vj / sceptr; the base set's ancestry PCs
absorb the ancestry-correlated part, so a gain over base is not ancestry re-detected.

Usage (from aleix/RNA-seq/, after 09, 10):
  pixi run python3 -u scripts/12_feature_sets.py
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402

KMER_DIM = 100
MIN_MEAN_USAGE, MIN_PRESENT = 0.001, 0.5


def clr(C, pseudo=0.5):
    L = np.log(C + pseudo)
    return L - L.mean(axis=1, keepdims=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    base = "~/pipeline_outputs/rnaseq"
    ap.add_argument("--person-table", default=f"{base}/analysis/person_table.pkl")
    ap.add_argument("--vj", default=f"{base}/analysis/vj_usage.pkl")
    ap.add_argument("--cache-dir", default=f"{base}/atlas")
    ap.add_argument("--embeddings-dir", default=f"{base}/embeddings")
    ap.add_argument("--local-out", default=f"{base}/analysis")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    t_all = time.time()
    pt = pd.read_pickle(os.path.expanduser(args.person_table))
    pt = pt[pt["analysis_main"]].set_index("research_id")
    ids = pt.index
    F = {}

    # ---- base covariates
    b = pd.DataFrame(index=ids)
    b["age"] = pt["age"]
    b["age2"] = (pt["age"] - pt["age"].mean()) ** 2
    b["female"] = (pt["sex_model"] == "Female").astype(float)
    b["sex_unknown"] = (pt["sex_model"] == "Unknown").astype(float)
    for k in range(1, 11):
        b[f"PC{k}"] = pt[f"PC{k}"]
    b["log10_trb_reads"] = np.log10(pt["trb_reads"])
    F["base"] = b.astype(float)

    # ---- diversity at fixed depth
    F["diversity"] = pt[["richness_at_depth", "shannon_at_depth", "inv_simpson_at_depth",
                         "public_frac_at_depth", "top10_frac"]].astype(float)

    # ---- V/J usage
    uc = pd.read_pickle(os.path.expanduser(args.vj))
    blocks = []
    for comp in ("V", "J"):
        C = uc[comp].reindex(ids).fillna(0)
        P = C.div(C.sum(1).replace(0, np.nan), axis=0)
        keep = (P.mean() >= MIN_MEAN_USAGE) & ((C > 0).mean() >= MIN_PRESENT)
        Ck = C.loc[:, keep]
        blocks.append(pd.DataFrame(clr(Ck.to_numpy(float)), index=ids, columns=Ck.columns))
    F["vj"] = pd.concat(blocks, axis=1)
    trbv_clr = blocks[0]

    # ---- clonotypes (full cache)
    t0 = time.time()
    clono = pd.read_pickle(os.path.join(os.path.expanduser(args.cache_dir), "trb_clonotypes_full.pkl"))
    clono = clono[clono["research_id"].astype(str).isin(set(ids))]
    rid = clono["research_id"].astype(str).to_numpy()
    pidx = pd.Series(np.arange(len(ids)), index=ids).reindex(rid).to_numpy()
    print(f"{len(clono):,} clonotypes for {len(ids):,} people ({time.time() - t0:.0f}s)", file=sys.stderr)

    # ---- SCEPTR mean over the whole repertoire
    edir = os.path.expanduser(args.embeddings_dir)
    E = np.load(os.path.join(edir, "embeddings_sceptr_fullcache.npy"), mmap_mode="r")
    keys = pd.read_csv(os.path.join(edir, "keys_sceptr_fullcache.tsv"), sep="\t", keep_default_na=False)
    lut = pd.Series(np.arange(len(keys)), index=keys["trbv"] + "|" + keys["cdr3aa"])
    row = lut.reindex(clono["trbv"].astype(str) + "|" + clono["cdr3aa"].astype(str)).to_numpy()
    ok = ~pd.isna(row)
    row, pi = row[ok].astype(np.int64), pidx[ok].astype(np.int64)
    order = np.argsort(pi, kind="stable")
    row, pi = row[order], pi[order]
    S = np.zeros((len(ids), E.shape[1]))
    n = np.bincount(pi, minlength=len(ids)).astype(float)
    for lo in range(0, len(row), 2_000_000):
        np.add.at(S, pi[lo:lo + 2_000_000], np.asarray(E[row[lo:lo + 2_000_000]], dtype=np.float64))
    S /= np.maximum(n, 1)[:, None]
    sc = pd.DataFrame(S, index=ids, columns=[f"sceptr_{k}" for k in range(S.shape[1])])
    F["sceptr"] = pd.concat([sc, trbv_clr.add_prefix("clr_")], axis=1)
    print(f"SCEPTR mean: {100 * ok.mean():.1f}% of clonotypes have a vector ({time.time() - t0:.0f}s)",
          file=sys.stderr)

    # ---- CDR3 3-mer composition -> SVD
    t0 = time.time()
    seqs = clono["cdr3aa"].astype(str).to_numpy()
    uniq, inv = np.unique(seqs, return_inverse=True)
    vocab, r_, c_ = {}, [], []
    for i, s in enumerate(uniq):
        for j in range(len(s) - 2):
            c_.append(vocab.setdefault(s[j:j + 3], len(vocab)))
            r_.append(i)
    K = csr_matrix((np.ones(len(r_), np.float32), (r_, c_)), shape=(len(uniq), len(vocab)))
    A = csr_matrix((np.ones(len(inv), np.float32), (pidx, inv)), shape=(len(ids), len(uniq)))
    PK = A @ K
    PK = csr_matrix(PK.multiply(1 / np.maximum(PK.sum(1).A1, 1)[:, None]))
    from sklearn.decomposition import TruncatedSVD
    Z = TruncatedSVD(KMER_DIM, random_state=args.seed).fit_transform(PK)
    F["kmer"] = pd.DataFrame(Z, index=ids, columns=[f"kmer_svd{k}" for k in range(KMER_DIM)])
    print(f"3-mers: {len(vocab):,} distinct over {len(uniq):,} unique CDR3s -> SVD {KMER_DIM} "
          f"({time.time() - t0:.0f}s)", file=sys.stderr)

    out = os.path.join(os.path.expanduser(args.local_out), "feature_sets.pkl")
    pd.to_pickle(F, out)
    for k, v in F.items():
        print(f"  {k:10s} {v.shape[0]:,} x {v.shape[1]}  missing {100 * v.isna().mean().mean():.2f}%",
              file=sys.stderr)
    print(f"Feature sets -> {out} (VM-local) in {time.time() - t_all:.0f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
