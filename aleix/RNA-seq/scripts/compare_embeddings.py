#!/usr/bin/env python3
"""Fold SCEPTR / ESMC / catELMo embeddings into one comparison table.

Run this LAST, after embed_cdr3s.py (SCEPTR + ESMC) and, optionally, embed_catelmo.py
(catELMo, run separately in its own conda env -- see setup_repertoire_and_embedding_envs.sh).
Works with whichever of the three .npy files actually exist -- 2 models is a valid partial
result, not an error.

Same check as embed_cdr3s.py: same-V-gene CDR3s should embed closer together than
different-V-gene ones, if the embedding captures anything biologically real. Larger gap =
more real structure recovered. A floor check, not a full evaluation -- passing it doesn't
prove a model is good at the harder downstream task (HLA/disease prediction), but failing it
is disqualifying.

Usage:
  pixi run python3 compare_embeddings.py <cdr3_pool.tsv>
      [--embeddings-dir ~/pipeline_outputs/rnaseq/embeddings] [--outdir ../results]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

MODELS = {
    "SCEPTR":   ("cdr3_embeddings_sceptr.npy",  "TCR-specific, 153K params"),
    "ESMC-300M":("cdr3_embeddings_esmc.npy",    "general protein LM, 300M params"),
    "catELMo":  ("cdr3_embeddings_catelmo.npy", "TCR-specific, 4-layer BiLSTM, ~93M params"),
}


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def same_vs_diff_vgene_contrast(embs, v_genes):
    norm = embs / (np.linalg.norm(embs, axis=1, keepdims=True) + 1e-9)
    sim = norm @ norm.T
    v = np.asarray(v_genes)
    same_mask = (v[:, None] == v[None, :])
    np.fill_diagonal(same_mask, False)
    diff_mask = ~same_mask
    np.fill_diagonal(diff_mask, False)
    return float(sim[same_mask].mean()), float(sim[diff_mask].mean())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pool_tsv", help="embed_cdr3s.py's cdr3_pool_used_for_embedding.tsv "
                                      "(row order must match what each embed_*.py script used)")
    ap.add_argument("--embeddings-dir", default=os.path.expanduser(
                        "~/pipeline_outputs/rnaseq/embeddings"))
    ap.add_argument("--outdir", default=os.path.join(
                        os.path.dirname(os.path.abspath(__file__)), "..", "results"))
    args = ap.parse_args()

    pool = pd.read_csv(os.path.expanduser(args.pool_tsv), sep="\t")
    v_genes = pool["v_gene"].tolist()
    n = len(pool)

    rows = []
    found_any = False
    for model, (fname, desc) in MODELS.items():
        path = os.path.join(os.path.expanduser(args.embeddings_dir), fname)
        if not os.path.exists(path):
            print(f"  [{model}] not found ({fname}) -- skipping, not an error", file=sys.stderr)
            continue
        embs = np.load(path)
        if embs.shape[0] != n:
            print(f"  [{model}] !! {embs.shape[0]} vectors but pool has {n} rows -- "
                  f"row order mismatch, skipping this model (re-run it against the same pool_tsv)",
                  file=sys.stderr)
            continue
        found_any = True
        same, diff = same_vs_diff_vgene_contrast(embs, v_genes)
        rows.append({"model": model, "description": desc, "dim": embs.shape[1],
                     "n_cdr3s": embs.shape[0], "same_vgene_cos_sim": round(same, 4),
                     "diff_vgene_cos_sim": round(diff, 4), "gap": round(same - diff, 4)})
        print(f"  [{model}] {embs.shape[1]}-dim, same-V {same:.4f} vs diff-V {diff:.4f}, "
              f"gap {same - diff:+.4f}")

    if not found_any:
        die("no embedding .npy files found -- run embed_cdr3s.py (and optionally "
            "embed_catelmo.py) first")

    df = pd.DataFrame(rows).sort_values("gap", ascending=False)
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "cdr3_embedding_comparison.csv")
    df.to_csv(out, index=False)

    print(f"\n=== {len(df)}-model comparison (ranked by gap, larger = more real structure "
          f"recovered) ===")
    print(df.to_string(index=False))
    print(f"\nWritten to {out} (de-identified, safe to commit -- aggregate stats only, no "
          f"research_ids or raw sequences).")


if __name__ == "__main__":
    main()
