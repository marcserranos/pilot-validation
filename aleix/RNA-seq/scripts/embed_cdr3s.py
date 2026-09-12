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

CATELMO NOTE: run separately, as embed_catelmo.py, inside its own `catELMo` conda env
(python 3.6 / allennlp 0.9.0 / torch 1.9.1 -- verified live against catELMo's own
embedders/README.md, 2026-09-07; an earlier pass wrongly assumed a TF2.6/Keras stack, which
is actually for catELMo's separate downstream binding-affinity trainer, not embedding).
That dependency set can't coexist with ESMC's modern transformers/torch in one interpreter,
so it isn't imported here -- run this script first (it writes the shared CDR3 pool TSV
below), then `conda run -n catELMo python3 embed_catelmo.py <pool.tsv>`, then
compare_embeddings.py to fold all three into one table. See
setup_repertoire_and_embedding_envs.sh for the one-shot environment setup.

Usage:
  pixi run python3 embed_cdr3s.py <cohort.tsv> [--pheno-dir ~/pipeline_outputs/rnaseq]
      [--min-score 0.02] [--keep-imputed] [--outdir ../results]
      [--local-outdir ~/pipeline_outputs/rnaseq/embeddings] [--max-per-person 500]

Needs (on top of the pixi env's pandas/scipy):
  pip install sceptr torch accelerate
  pip install git+https://github.com/Biohub/esm.git@v3.4.1
(NOT `pip install esm` alone -- PyPI's `esm` is still 3.2.3 as of 2026-09-12 and its ESMC
loading is broken upstream (biohub/esmc-300m-2024-12's HF config.json is empty, confirmed
live). v3.4.1, not yet on PyPI, ships a different loading path that works -- see
embed_esmc()'s docstring below for the full story.)
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

AA_VALID = set("ACDEFGHIKLMNPQRSTVWY")

# Standard genetic code. TRUST4's own convention (README + report.tsv docs): stop codon ->
# "_", any codon touching an ambiguous base (N, or anything non-ACGT) -> "?". Matched here
# exactly so translating cdr3.out's CDR3 *nucleotide* field ourselves (see below -- cdr3.out
# has NO amino-acid column at all, verified live against TRUST4's real output 2026-09-07,
# correcting an earlier wrong assumption that it did) produces the same convention TRUST4
# itself uses in report.tsv's CDR3aa column.
_CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'_','TAG':'_','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'_','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}


def translate_cdr3_dna(dna):
    """DNA -> amino acid, TRUST4's own convention: '_' for stop, '?' for any codon with a
    non-ACGT base or a length not divisible by 3 (frameshift/truncated -- can't translate)."""
    dna = str(dna).upper()
    if len(dna) % 3 != 0:
        return "?"
    aa = []
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]
        aa.append(_CODON_TABLE.get(codon, "?"))
    return "".join(aa)


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
        # cdr3.out is ALWAYS headerless (verified against TRUST4's own bundled example
        # output, 2026-09-07) -- the real 13-field schema, straight from TRUST4's README:
        #   consensus_id  index_within_consensus  V  D  J  C  CDR1  CDR2  CDR3(dna)
        #   CDR3_score  read_fragment_count  CDR3_germline_similarity  complete_vdj_assembly
        # NOTE: that CDR3 field is nucleotide, not amino acid -- cdr3.out has no AA column at
        # all. We translate it ourselves below so both input paths end up with real cdr3aa.
        cols = ["consensus_id", "V", "D", "J", "C", "CDR1", "CDR2", "CDR3_dna",
                "CDR3_score", "read_fragment_count", "CDR3_germline_similarity",
                "complete_vdj_assembly"]
        df = pd.read_csv(cdr3_out, sep="\t", header=None, names=cols)
        df = df.rename(columns={"CDR3_score": "score", "V": "v_gene"})
        df["cdr3aa"] = df["CDR3_dna"].apply(translate_cdr3_dna)
        df = df[pd.to_numeric(df["score"], errors="coerce") >= threshold]
    elif os.path.exists(report):
        print(f"  [{research_id}] no cdr3.out, falling back to report.tsv (no CDR3_score "
              f"available -- quality filter skipped, only stop-codon/ambiguous drop applies)",
              file=sys.stderr)
        df = pd.read_csv(report, sep="\t")
        # Real header (verified live): #count, frequency, CDR3nt, CDR3aa, V, D, J, C, cid,
        # cid_full_length -- "CDR3aa", not "CDR3_amino_acids" (an earlier wrong assumption).
        df = df.rename(columns={"CDR3aa": "cdr3aa", "V": "v_gene"})
        df["score"] = np.nan
    else:
        return None

    if "cdr3aa" not in df.columns or "v_gene" not in df.columns:
        print(f"  [{research_id}] !! unexpected columns, skipping: {list(df.columns)}",
              file=sys.stderr)
        return None

    # TRUST4 lists up to 3 ranked V-gene candidates comma-separated (e.g.
    # "IGHV3-11*04,IGHV3-21*01,IGHV3-48*01" -- verified live in the real example output).
    # Keep only the top-ranked candidate -- otherwise the same-V-gene-vs-diff-V-gene
    # comparison silently breaks (two CDR3s sharing the top candidate but differing in
    # ranked-2nd/3rd wouldn't match as "same V gene" on a raw string-equality basis).
    df["v_gene"] = df["v_gene"].astype(str).str.split(",").str[0]

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


def embed_esmc(seqs, model_name="esmc_300m"):
    """ESMC-300M -- general protein LM. Returns (embeddings ndarray, wall_seconds).
    Mean-pools per-residue hidden states over the sequence (excluding BOS/EOS).

    NOTE (corrected 2026-09-12, third attempt -- this one actually works end to end):
    both the plain transformers.AutoModel path AND the esm package's own
    EsmcForMaskedLM.from_pretrained("biohub/esmc-300m-2024-12") are broken -- that HF
    repo's config.json is a live upstream bug, confirmed empty ({}) as recently as
    2026-09-12. Not fixable locally.

    The fix: `esm` v3.4.0/3.4.1 (2026-08-27/09-08, NOT yet on PyPI as of this writing --
    PyPI's `esm` is still 3.2.3, install straight from GitHub instead:
    `pip install git+https://github.com/Biohub/esm.git@v3.4.1`) shipped a genuinely
    different loading path and API: `esm.models.esmc.ESMC.from_pretrained("esmc_300m")`
    (the bare model name, NOT the HF repo id) plus an SDK-client-shaped interface
    (`ESMProtein` -> `model.encode()` -> `model.logits(..., LogitsConfig(return_embeddings=True))`)
    that does not depend on the broken repo's config.json at all. Verified live: real
    960-dim embeddings, correct shape (len(seq)+2 for BOS/EOS).

    Also needs `pip install accelerate` (a new dependency this esm version pulls in that
    the old PyPI 3.2.3 didn't require -- easy to miss if esm is upgraded with --no-deps).

    One ESMProtein per call (this is an SDK-client-shaped API, mirroring the hosted
    Forge/Biohub Platform interface) -- no native batch-tensor call found, so this loops
    per sequence. CDR3s are short (10-20 aa) so each call is fast, but per-call Python/SDK
    overhead dominates at scale; time a small sample before committing to a large pool."""
    import torch
    from esm.models.esmc import ESMC
    from esm.sdk.api import ESMProtein, LogitsConfig

    model = ESMC.from_pretrained(model_name).eval()
    cfg = LogitsConfig(sequence=True, return_embeddings=True)

    out = []
    t0 = time.time()
    with torch.inference_mode():
        for s in seqs:
            encoded = model.encode(ESMProtein(sequence=s))
            res = model.logits(encoded, cfg)
            emb = res.embeddings[0]         # (L, D), L = len(s) + 2 (BOS/EOS)
            pooled = emb[1:-1].mean(dim=0)  # drop BOS/EOS, mean over real residues
            out.append(pooled.cpu().numpy())
    return np.stack(out, axis=0), time.time() - t0


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
    ap.add_argument("cohort", nargs="?", default=None,
                    help="cohort.tsv with a research_id column (build_rnaseq_cohort.py output). "
                         "Omit if using --from-pool-tsv instead.")
    ap.add_argument("--from-pool-tsv", default=None,
                    help="skip TRUST4 loading entirely and embed a pre-built pool TSV instead "
                         "(columns: research_id, chain, v_gene, cdr3aa, score) -- e.g. from "
                         "prep_vdjdb_pool.py, for testing the same two models against an open "
                         "dataset while the AoU mount is unreachable. Same downstream code "
                         "either way, only the input source differs.")
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

    if args.from_pool_tsv:
        pool = pd.read_csv(os.path.expanduser(args.from_pool_tsv), sep="\t")
        need = {"research_id", "chain", "v_gene", "cdr3aa"}
        missing = need - set(pool.columns)
        if missing:
            die(f"{args.from_pool_tsv} is missing column(s) {missing} -- expected the same "
                f"schema this script itself writes (see prep_vdjdb_pool.py for an example "
                f"adapter)")
        print(f"Using pre-built pool from {args.from_pool_tsv} (public/open dataset mode -- "
              f"no TRUST4 output needed).", file=sys.stderr)
    else:
        if not args.cohort:
            die("need either a cohort.tsv (TRUST4 mode) or --from-pool-tsv (open-dataset mode)")
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
              f"(need: pip install accelerate && pip install "
              f"git+https://github.com/Biohub/esm.git@v3.4.1 -- PyPI's esm 3.2.3 and the "
              f"old EsmcForMaskedLM path are both broken upstream, see embed_esmc() docstring)",
              file=sys.stderr)

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
