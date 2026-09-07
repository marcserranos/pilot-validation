#!/usr/bin/env python3
"""Third embedding arm: catELMo (TCR-specific, 4-layer BiLSTM, 1,024-dim), run separately
from embed_cdr3s.py because it needs python 3.6 / allennlp 0.9.0 / torch 1.9.1 -- a
dependency set that cannot coexist with ESMC's modern transformers/torch in one env. Run
this INSIDE the `catELMo` conda env that setup_repertoire_and_embedding_envs.sh builds:

  conda run -n catELMo python3 embed_catelmo.py <cdr3_pool.tsv> --weights-dir ~/tools/catELMo_weights

<cdr3_pool.tsv> is embed_cdr3s.py's own output (research_id, chain, v_gene, cdr3aa, score)
-- run that script FIRST so both scripts embed the exact same filtered CDR3 set. This keeps
the three-way comparison honest: same input, three models, nothing but the model differs.

Recipe follows catELMo's own embedders/README.md exactly (AllenNLP's ElmoEmbedder,
cuda_device=-1 for CPU, sum-then-mean pooling over the per-layer representations -- that
pooling choice is catELMo's own documented example, not something invented here).

Writes cdr3_embeddings_catelmo.npy next to embed_cdr3s.py's other two .npy files, same
naming convention, so a downstream compare-all-three step just globs for
cdr3_embeddings_*.npy.
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pool_tsv", help="embed_cdr3s.py's cdr3_pool_used_for_embedding.tsv")
    ap.add_argument("--weights-dir", default=os.path.expanduser("~/tools/catELMo_weights"))
    ap.add_argument("--outdir", default=os.path.expanduser("~/pipeline_outputs/rnaseq/embeddings"))
    args = ap.parse_args()

    try:
        import torch
        from allennlp.commands.elmo import ElmoEmbedder
    except ImportError as e:
        die(f"missing catELMo deps ({e}) -- this must run inside the 'catELMo' conda env: "
            f"conda run -n catELMo python3 {sys.argv[0]} ...")

    weights = os.path.join(args.weights_dir, "weights.hdf5")
    options = os.path.join(args.weights_dir, "options.json")
    for p in (weights, options):
        if not os.path.exists(p):
            die(f"missing {p} -- run setup_repertoire_and_embedding_envs.sh first, or check "
                f"{args.weights_dir} by hand (Dropbox folder layout can vary).")

    pool = pd.read_csv(os.path.expanduser(args.pool_tsv), sep="\t")
    if "cdr3aa" not in pool.columns:
        die(f"{args.pool_tsv} has no cdr3aa column -- expected embed_cdr3s.py's output")
    seqs = pool["cdr3aa"].tolist()
    print(f"Embedding {len(seqs)} CDR3s with catELMo (CPU, cuda_device=-1)...")

    embedder = ElmoEmbedder(options, weights, cuda_device=-1)

    def embed_one(seq):
        # catELMo's own documented pooling: sum across the 3 BiLSTM layers, mean across
        # sequence position -- see embedders/README.md Example 1.
        return torch.tensor(embedder.embed_sentence(list(seq))).sum(dim=0).mean(dim=0).tolist()

    t0 = time.time()
    vecs = [embed_one(s) for s in seqs]
    elapsed = time.time() - t0
    embs = np.asarray(vecs)

    os.makedirs(os.path.expanduser(args.outdir), exist_ok=True)
    out_npy = os.path.join(os.path.expanduser(args.outdir), "cdr3_embeddings_catelmo.npy")
    np.save(out_npy, embs)

    print(f"  {elapsed:.1f}s wall, {embs.shape[1]}-dim vectors, {embs.shape[0]} CDR3s")
    print(f"  saved to {out_npy}")
    print("\nThis script only embeds -- run the same-V-gene-vs-diff-V-gene comparison by "
          "loading this .npy alongside embed_cdr3s.py's two (same v_gene column, same row "
          "order from the shared pool_tsv) and reusing its same_vs_diff_vgene_contrast().")


if __name__ == "__main__":
    main()
