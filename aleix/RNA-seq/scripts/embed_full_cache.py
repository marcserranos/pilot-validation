#!/usr/bin/env python3
"""Embed every distinct (TRBV, CDR3) in 02's full clonotype cache with b_sceptr, once.

WHY: embed_cdr3s.py embeds each person's top 500 clonotypes. Report 06 needs a vector for
every clonotype of a person, and until now looked vectors up from the top-500 pool, which
only covers sequences that are in somebody's top 500 (an estimated half of the cache,
biased toward expanded and public sequences). SCEPTR is deterministic -- a sequence gets the
same vector whoever carries it -- so embedding each distinct (TRBV, CDR3) once covers the
whole cache.

Same model and input as the production run (b_sceptr, TRBV + CDR3B; TRBV symbols SCEPTR
cannot use are excluded and counted). Checkpointed per chunk under --local-outdir/_fullcache_chunks/
so an interrupted run resumes where it stopped.

Outputs (VM-local; no research_ids, but never committed):
  embeddings_sceptr_fullcache.npy   float32, one row per key
  keys_sceptr_fullcache.tsv         trbv, cdr3aa (row-aligned)

Usage (from aleix/RNA-seq/, after 02):
  pixi run python3 -u scripts/embed_full_cache.py [--chunk 250000]
"""
import argparse
import os
import shutil
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from embed_cdr3s import accepted_trbv, sceptr_model  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache-dir", default="~/pipeline_outputs/rnaseq/atlas")
    ap.add_argument("--local-outdir", default="~/pipeline_outputs/rnaseq/embeddings")
    ap.add_argument("--chunk", type=int, default=250_000)
    args = ap.parse_args()
    out = os.path.expanduser(args.local_outdir)
    ck_dir = os.path.join(out, "_fullcache_chunks")
    os.makedirs(ck_dir, exist_ok=True)

    clono = pd.read_pickle(os.path.join(os.path.expanduser(args.cache_dir),
                                        "trb_clonotypes_full.pkl"))
    keys = (clono[["trbv", "cdr3aa"]].astype(str).drop_duplicates()
            .sort_values(["trbv", "cdr3aa"]).reset_index(drop=True))
    n_all = len(keys)
    model = sceptr_model("v+cdr3")
    good = accepted_trbv(model, keys["trbv"])
    keys = keys[keys["trbv"].isin(good)].reset_index(drop=True)
    print(f"{len(clono):,} cache rows, {n_all:,} distinct (TRBV, CDR3); {len(keys):,} with a "
          f"TRBV SCEPTR accepts ({len(good)} genes)", file=sys.stderr)

    t0 = time.time()
    n_chunks = (len(keys) + args.chunk - 1) // args.chunk
    for i in range(n_chunks):
        f = os.path.join(ck_dir, f"chunk_{i:05d}.npy")
        if os.path.exists(f):
            continue
        part = keys.iloc[i * args.chunk:(i + 1) * args.chunk]
        X = np.asarray(model.calc_vector_representations(
            pd.DataFrame({"TRBV": part["trbv"].to_numpy(), "CDR3B": part["cdr3aa"].to_numpy()})),
            dtype=np.float32)
        np.save(f + ".tmp.npy", X)
        os.replace(f + ".tmp.npy", f)
        print(f"  chunk {i + 1}/{n_chunks} ({min((i + 1) * args.chunk, len(keys)):,} keys, "
              f"{time.time() - t0:.0f}s)", file=sys.stderr)

    embs = np.concatenate([np.load(os.path.join(ck_dir, f"chunk_{i:05d}.npy"))
                           for i in range(n_chunks)])
    if len(embs) != len(keys):
        sys.exit(f"FATAL: {len(embs)} vectors for {len(keys)} keys -- delete {ck_dir} and rerun")
    np.save(os.path.join(out, "embeddings_sceptr_fullcache.npy"), embs)
    keys.to_csv(os.path.join(out, "keys_sceptr_fullcache.tsv"), sep="\t", index=False)
    shutil.rmtree(ck_dir)
    print(f"Wrote {embs.shape[0]:,} x {embs.shape[1]} vectors in {time.time() - t0:.0f}s -> {out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
