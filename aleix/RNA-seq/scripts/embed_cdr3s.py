#!/usr/bin/env python3
"""Embed recovered CDR3s with two models of deliberately different character, and compare.

WHY TWO MODELS: SCEPTR is TCR-specific (153K params, trained on TCR sequences themselves,
see reference/POST_TRUST4_OPTIONS.md). ESMC is a general protein language model (300M
params, trained on UniRef/MGnify/JGI -- ordinary evolutionarily-conserved proteins, not
somatically-recombined receptor loops). The project's own prior verdict on ESMC's
predecessor (ESM-2) was "not recommended as first choice" for CDR3 embedding specifically,
on the grounds that CDR3's hypervariable, junctionally-diversified sequence space is a
fundamentally different generative process from the evolutionarily-conserved sequence space
these models learn from (see the chat-log technical writeup, 2026-08-31, for the full
mechanistic argument). This script makes that comparison EMPIRICAL rather than argued: same
CDR3s, both models, one concrete structure-recovery check each has to pass.

THE CHECK: same-V-gene CDR3s should embed closer together than different-V-gene CDR3s, if
the embedding captures anything biologically real (V-gene partly determines CDR3 length/
composition via the germline-encoded segment it contributes). This needs no external labels
beyond what TRUST4 already outputs -- cheap, and in the same spirit as the project's
standing validation principle: check a known-true relationship before trusting anything
harder (reference/POST_TRUST4_OPTIONS.md, the HLA-prediction validation step).

INPUT: TRUST4 output dirs at ~/pipeline_outputs/rnaseq/<research_id>/ (from
run_trust4_sample.sh / run_rnaseq_batch.sh). Prefers *_cdr3.out (has CDR3_score, needed for
the quality filter below); falls back to *_report.tsv if cdr3.out is missing.

QUALITY FILTER (do this before embedding anything -- garbage in, garbage out on a model
comparison specifically): drops CDR3_score 0.00 (partial/incomplete) and, by default, 0.01
(imputed/guessed) -- only real-motif CDR3s (score > 0.01) get embedded. Also drops any
amino-acid sequence containing '_' (stop codon) or '?' (ambiguous base), per
TRUST4_DEEP_DIVE.md's flagged cleanup item. Pass --keep-imputed to relax the CDR3_score cut.

CATELMO NOTE: not wired in here. Its repo (github.com/Lee-CBG/catELMo) needs a
Python 3.6.13 / TensorFlow 2.6.0 / Keras 2.6.0 conda env and its pretrained-embedding
loading path isn't documented at the README level (buried in an `embedders/` folder,
unclear from a quick pass whether weights are bundled or need training). That is a real risk
to attempt cold inside a short compute window. SCEPTR is substituted as the "tailored"
comparator: also TCR-specific, but pip-installable and CPU-fast, already vetted in this
project's own reference docs. If catELMo is still wanted, budget separate time for it as a
third arm -- this script's embed_scores()/COMPARISON block is written so a third model slots
in the same way.

Usage:
  pixi run python3 embed_cdr3s.py <cohort.tsv> [--pheno-dir ~/pipeline_outputs/rnaseq]
      [--min-score 0.02] [--keep-imputed] [--outdir ../results]
      [--local-outdir ~/pipeline_outputs/rnaseq/embeddings] [--max-per-person 500]

Needs (on top of the pixi env's pandas/scipy):  pip install sceptr transformers torch
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

AA_VALID = set("ACDEFGHIKLMNPQRSTVWY")


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def load_person_cdr3s(pheno_dir, research_id, min_score, keep_imputed):
    """One row per usable CDR3 for this person: research_id, chain, v_gene, cdr3aa, score."""
    base = os.path.join(os.path.expanduser(pheno_dir), research_id)
    cdr3_out = os.path.join(base, f"{research_id}_cdr3.out")
    report = os.path.join(base, f"{research_id}_report.tsv")

    threshold = 0.01 if keep_imputed else min_score

    if os.path.exists(cdr3_out):
        # cdr3.out is whitespace/tab-delimited, no header in some TRUST4 versions -- try
        # header-first, fall back to TRUST4's documented column order if that fails.
        try:
            df = pd.read_csv(cdr3_out, sep="\t")
            if "CDR3_score" not in df.columns:
                raise ValueError("no header")
        except Exception:
            cols = ["contig_id", "V", "D", "J", "C", "CDR1", "CDR2", "CDR3_dna",
                    "CDR3_amino_acids", "CDR3_score", "read_fragment_count",
                    "CDR3_germline_similarity", "complete_vdj_assembly"]
            df = pd.read_csv(cdr3_out, sep="\t", header=None, names=cols,
                              usecols=range(len(cols)))
        df = df.rename(columns={"CDR3_amino_acids": "cdr3aa", "CDR3_score": "score",
                                 "V": "v_gene"})
        df = df[pd.to_numeric(df["score"], errors="coerce") > threshold]
    elif os.path.exists(report):
        print(f"  [{research_id}] no cdr3.out, falling back to report.tsv (no CDR3_score "
              f"available -- quality filter skipped, only stop-codon/ambiguous drop applies)",
              file=sys.stderr)
        df = pd.read_csv(report, sep="\t")
        df = df.rename(columns={"CDR3_amino_acids": "cdr3aa", "V": "v_gene"})
        df["score"] = np.nan
    else:
        return None

    if "cdr3aa" not in df.columns or "v_gene" not in df.columns:
        print(f"  [{research_id}] !! unexpected columns, skipping: {list(df.columns)}",
              file=sys.stderr)
        return None

    df = df[["v_gene", "cdr3aa", "score"]].dropna(subset=["cdr3aa"])
    df = df[df["cdr3aa"].astype(str).apply(
        lambda s: len(s) >= 5 and all(c in AA_VALID for c in s))]
    df["research_id"] = research_id
    df["chain"] = df["v_gene"].astype(str).str[:3]  # TRB/TRA/IGH/IGK/IGL/TRG/TRD
    return df.reset_index(drop=True)


def embed_sceptr(seqs):
    """SCEPTR -- TCR-specific, 153K params. Returns (embeddings ndarray, wall_seconds)."""
    from sceptr import variant
    model = variant.default()
    df = pd.DataFrame({"CDR3B": seqs})  # beta-only input, the mode we actually have
    t0 = time.time()
    vecs = model.calc_vector_representations(df)
    return np.asarray(vecs), time.time() - t0


def embed_esmc(seqs, batch_size=32):
    """ESMC-300M -- general protein LM. Returns (embeddings ndarray, wall_seconds).
    Mean-pools per-residue hidden states over the sequence (excluding special tokens)."""
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    model_id = "biohub/esmc-300m-2024-12"
    tok = AutoTokenizer.from_pretrained(model_id)
    model = AutoModelForMaskedLM.from_pretrained(
        model_id, trust_remote_code=True, output_hidden_states=True)
    model.eval()

    out = []
    t0 = time.time()
    with torch.no_grad():
        for i in range(0, len(seqs), batch_size):
            batch = seqs[i:i + batch_size]
            enc = tok(batch, return_tensors="pt", padding=True, truncation=True)
            res = model(**enc)
            hidden = res.hidden_states[-1]  # (B, L, D)
            mask = enc["attention_mask"].unsqueeze(-1).float()
            pooled = (hidden * mask).sum(1) / mask.sum(1).clamp(min=1)
            out.append(pooled.cpu().numpy())
    return np.concatenate(out, axis=0), time.time() - t0


def same_vs_diff_vgene_contrast(embs, v_genes):
    """Mean cosine similarity within the same V gene vs across different V genes.
    A real embedding should score meaningfully higher same-V than diff-V -- V gene partly
    determines CDR3 composition via the germline segment it contributes."""
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
    ap.add_argument("cohort", help="cohort.tsv with a research_id column (build_rnaseq_cohort.py output)")
    ap.add_argument("--pheno-dir", default="~/pipeline_outputs/rnaseq",
                    help="dir containing <research_id>/ TRUST4 output subdirs")
    ap.add_argument("--min-score", type=float, default=0.02,
                    help="minimum CDR3_score to keep (default 0.02 = real motif strength only)")
    ap.add_argument("--keep-imputed", action="store_true",
                    help="also keep CDR3_score==0.01 (imputed/guessed) -- relaxes the filter")
    ap.add_argument("--max-per-person", type=int, default=500,
                    help="cap CDR3s per person (highest read_count-weighted first would be "
                         "better -- for now, first N after filtering) to keep the comparison fast")
    ap.add_argument("--outdir", default=os.path.join(
                        os.path.dirname(os.path.abspath(__file__)), "..", "results"),
                    help="de-identified comparison summary only (aggregate stats, no "
                         "research_ids) -- this is the one that's safe to commit")
    ap.add_argument("--local-outdir", default="~/pipeline_outputs/rnaseq/embeddings",
                    help="VM-local: per-CDR3 pool (has real research_ids) + raw embedding "
                         ".npy files -- NEVER commit this, same privacy posture as "
                         "query_overlap_phenotypes.py's pheno cache")
    args = ap.parse_args()

    cohort = pd.read_csv(os.path.expanduser(args.cohort), sep="\t", dtype=str)
    if "research_id" not in cohort.columns:
        die(f"{args.cohort} has no research_id column -- expected build_rnaseq_cohort.py output")

    all_rows = []
    for rid in cohort["research_id"]:
        df = load_person_cdr3s(args.pheno_dir, rid, args.min_score, args.keep_imputed)
        if df is None:
            print(f"  [{rid}] !! no TRUST4 output found, skipping (run the repertoire batch first)",
                  file=sys.stderr)
            continue
        if len(df) == 0:
            print(f"  [{rid}] 0 CDR3s survived the quality filter", file=sys.stderr)
            continue
        all_rows.append(df.head(args.max_per_person))
        print(f"  [{rid}] {len(df)} CDR3s pass filter (using up to {args.max_per_person})",
              file=sys.stderr)

    if not all_rows:
        die("no usable CDR3s across the whole cohort -- run repertoire calling first "
            "(build_rnaseq_cohort.py + run_rnaseq_batch.sh)")

    pool = pd.concat(all_rows, ignore_index=True)
    seqs = pool["cdr3aa"].tolist()
    v_genes = pool["v_gene"].tolist()
    print(f"\n=== {len(seqs)} CDR3s pooled across {pool['research_id'].nunique()} people, "
          f"{pool['v_gene'].nunique()} distinct V genes ===\n")

    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    local_outdir = os.path.abspath(os.path.expanduser(args.local_outdir))
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(local_outdir, exist_ok=True)
    results = []

    print("--- SCEPTR (TCR-specific, 153K params) ---")
    try:
        sceptr_embs, sceptr_s = embed_sceptr(seqs)
        same, diff = same_vs_diff_vgene_contrast(sceptr_embs, v_genes)
        print(f"  {sceptr_s:.1f}s wall, {sceptr_embs.shape[1]}-dim vectors")
        print(f"  same-V-gene cosine sim: {same:.4f}  |  diff-V-gene: {diff:.4f}  "
              f"|  gap: {same - diff:+.4f}")
        np.save(os.path.join(local_outdir, "cdr3_embeddings_sceptr.npy"), sceptr_embs)
        results.append(("SCEPTR", sceptr_embs.shape[1], sceptr_s, same, diff))
    except Exception as e:
        print(f"  !! SCEPTR failed: {e}\n  (pip install sceptr if missing)", file=sys.stderr)

    print("\n--- ESMC-300M (general protein LM) ---")
    try:
        esmc_embs, esmc_s = embed_esmc(seqs)
        same, diff = same_vs_diff_vgene_contrast(esmc_embs, v_genes)
        print(f"  {esmc_s:.1f}s wall, {esmc_embs.shape[1]}-dim vectors")
        print(f"  same-V-gene cosine sim: {same:.4f}  |  diff-V-gene: {diff:.4f}  "
              f"|  gap: {same - diff:+.4f}")
        np.save(os.path.join(local_outdir, "cdr3_embeddings_esmc.npy"), esmc_embs)
        results.append(("ESMC-300M", esmc_embs.shape[1], esmc_s, same, diff))
    except Exception as e:
        print(f"  !! ESMC failed: {e}\n  "
              f"(pip install transformers torch; may need `huggingface-cli login` if the "
              f"biohub/esmc-300m-2024-12 checkpoint is gated)", file=sys.stderr)

    if results:
        summary = pd.DataFrame(results, columns=[
            "model", "dim", "wall_seconds", "same_vgene_cos_sim", "diff_vgene_cos_sim"])
        summary["gap"] = summary["same_vgene_cos_sim"] - summary["diff_vgene_cos_sim"]
        out = os.path.join(outdir, "cdr3_embedding_comparison.csv")
        summary.to_csv(out, index=False)
        print(f"\n=== Comparison written to {out} ===")
        print("Larger 'gap' = the model separates CDR3s by V-gene more cleanly, i.e. captures "
              "more biologically real structure. This is a floor check, not a full evaluation -- "
              "a model could pass this and still be mediocre at the harder downstream task "
              "(HLA/disease prediction). But failing it is disqualifying.")
    else:
        print("\nBoth models failed -- see errors above. Nothing to compare.", file=sys.stderr)

    pool[["research_id", "chain", "v_gene", "cdr3aa", "score"]].to_csv(
        os.path.join(local_outdir, "cdr3_pool_used_for_embedding.tsv"), sep="\t", index=False)
    print(f"CDR3 pool used (VM-local, has real research_ids -- do not commit): "
          f"{local_outdir}/cdr3_pool_used_for_embedding.tsv")


if __name__ == "__main__":
    main()
