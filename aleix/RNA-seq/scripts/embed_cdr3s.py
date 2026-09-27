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

Usage (from aleix/RNA-seq/, so pixi finds the env where sceptr is installed):
  pixi run python3 scripts/embed_cdr3s.py <cohort.tsv> [--models sceptr,esmc]
      [--sceptr-input v+cdr3|cdr3] [--max-per-person 500] [--min-score 0.02]
      [--keep-imputed] [--allow-report-fallback] [--tag NAME]

Per person: unique (chain, V, CDR3aa) clonotypes, canonical junctions only (C...F/W),
ranked by summed read support, top --max-per-person kept.

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


CHAINS = {"TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL"}


def _top_gene(s):
    """TRUST4 lists up to 3 ranked candidates comma-separated ("TRBV6-2*01,TRBV6-3*01");
    keep the top one. '*' / '.' / empty means no call."""
    s = str(s).split(",")[0].strip()
    return "" if s in ("", "*", ".", "nan") else s


def _chain_of(v, j, c):
    """Chain from the first gene that was actually called, V then J then C. V-only labelling
    silently dropped CDR3s whose V wasn't called even when J/C clearly say TRB."""
    for g in (v, j, c):
        if g[:3] in CHAINS:
            return g[:3]
    return ""


def load_person_cdr3s(pheno_dir, research_id, min_score, keep_imputed, allow_report_fallback):
    """One row per unique clonotype (chain, V gene, CDR3aa) for this person, ranked by
    read support, highest first. Returns (df, source) where source is 'cdr3.out',
    'report.tsv', or None if nothing usable."""
    base = os.path.join(os.path.expanduser(pheno_dir), research_id)
    cdr3_out = os.path.join(base, f"{research_id}_cdr3.out")
    report = os.path.join(base, f"{research_id}_report.tsv")

    threshold = 0.01 if keep_imputed else min_score

    if os.path.exists(cdr3_out):
        # Headerless, 13 fields (TRUST4 README; verified against its bundled example
        # 2026-09-07). CDR3 is nucleotide -- translated below. All 13 named explicitly: with
        # 12 names pandas silently turns the first field into the index.
        cols = ["consensus_id", "index_within_consensus", "V", "D", "J", "C", "CDR1",
                "CDR2", "CDR3_dna", "score", "reads", "CDR3_germline_similarity",
                "complete_vdj_assembly"]
        df = pd.read_csv(cdr3_out, sep="\t", header=None, names=cols, index_col=False)
        df["cdr3aa"] = df["CDR3_dna"].apply(translate_cdr3_dna)
        df = df[pd.to_numeric(df["score"], errors="coerce") >= threshold]
        source = "cdr3.out"
    elif os.path.exists(report) and allow_report_fallback:
        # No CDR3_score in report.tsv, so the quality filter can't be applied. Off by
        # default: mixing filtered and unfiltered people makes the pool inhomogeneous.
        df = pd.read_csv(report, sep="\t")
        df = df.rename(columns={"CDR3aa": "cdr3aa", "#count": "reads"})
        df["score"] = np.nan
        source = "report.tsv"
    else:
        return None, None

    for g in ("V", "J", "C"):
        df[g] = df[g].apply(_top_gene)
    df["chain"] = [_chain_of(v, j, c) for v, j, c in zip(df["V"], df["J"], df["C"])]
    df["reads"] = pd.to_numeric(df["reads"], errors="coerce").fillna(0)

    # Productive, canonical junctions only: valid residues, starts with the conserved C,
    # ends with the conserved F/W (the format SCEPTR is trained on; also excludes the
    # truncated CDR3s TRUST4 reports when an assembly doesn't span the whole junction).
    aa = df["cdr3aa"].astype(str)
    ok = aa.str.fullmatch(r"C[ACDEFGHIKLMNPQRSTVWY]{3,}[FW]")
    df = df[ok & (df["chain"] != "")]

    # Collapse to unique clonotypes: the same CDR3 can appear in several TRUST4 consensus
    # assemblies. Sum their read support, then rank so a per-person cap keeps the
    # dominant clonotypes rather than whatever happened to come first in the file.
    df = (df.groupby(["chain", "V", "cdr3aa"], as_index=False)
            .agg(reads=("reads", "sum"), score=("score", "max"))
            .sort_values("reads", ascending=False, kind="stable"))
    df = df.rename(columns={"V": "v_gene"})
    df["research_id"] = research_id
    return df.reset_index(drop=True), source


def _imgt_gene(v):
    """'TRBV10-3*01' -> 'TRBV10-3' (IMGT gene symbol, allele stripped)."""
    return v.split("*")[0] if v else ""


def sceptr_model(sceptr_input):
    from sceptr import variant
    # b_sceptr: SCEPTR's beta-chain-only variant, the right model for unpaired bulk data.
    # cdr3_only: CDR3 loop alone -- used for the V-gene benchmark, where feeding the V gene
    # in would make the same-V vs diff-V check circular.
    return variant.b_sceptr() if sceptr_input == "v+cdr3" else variant.cdr3_only()


def accepted_trbv(model, v_genes):
    """Which TRBV symbols SCEPTR accepts (IMGT-functional genes only, per its docs).
    Probed one gene at a time so a single pseudogene/ORF call can't crash the whole run."""
    good = set()
    for v in sorted(set(v_genes) - {""}):
        try:
            model.calc_vector_representations(
                pd.DataFrame({"TRBV": [v], "CDR3B": ["CASSLGQGAEAFF"]}))
            good.add(v)
        except Exception:
            pass
    return good


def embed_sceptr(model, df_in, chunk=200_000):
    """Returns (embeddings ndarray, wall_seconds). Chunked for progress + bounded memory."""
    t0 = time.time()
    out = []
    for i in range(0, len(df_in), chunk):
        out.append(np.asarray(model.calc_vector_representations(
            df_in.iloc[i:i + chunk].reset_index(drop=True))))
        print(f"    {min(i + chunk, len(df_in)):,}/{len(df_in):,} embedded "
              f"({time.time() - t0:.0f}s)", file=sys.stderr)
    return np.concatenate(out, axis=0), time.time() - t0


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
    """Mean cosine similarity over all same-V-gene pairs vs all different-V-gene pairs
    (self-pairs excluded). A real embedding should score higher same-V than diff-V.

    Exact, in O(n*d) memory: for unit vectors, the sum of pairwise cosine sims within a
    group is ||sum of its vectors||^2 - n_group. The old n x n matrix needed ~60 TB at
    the full-cohort pool size."""
    norm = embs / (np.linalg.norm(embs, axis=1, keepdims=True) + 1e-9)
    n = len(norm)
    total = float(np.sum(norm.sum(axis=0) ** 2)) - n
    codes, groups = pd.factorize(pd.Series(v_genes))
    same_sum, same_pairs = 0.0, 0
    for g in range(len(groups)):
        vecs = norm[codes == g]
        k = len(vecs)
        same_sum += float(np.sum(vecs.sum(axis=0) ** 2)) - k
        same_pairs += k * (k - 1)
    diff_pairs = n * (n - 1) - same_pairs
    return same_sum / same_pairs, (total - same_sum) / diff_pairs


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
                    help="cap unique clonotypes per person, keeping the highest read-support "
                         "ones (default 500)")
    ap.add_argument("--sceptr-input", choices=["v+cdr3", "cdr3"], default="v+cdr3",
                    help="v+cdr3 (default): b_sceptr, the beta-only variant, given TRBV + CDR3B "
                         "-- the full beta-chain input, for HLA/disease work; clonotypes whose V "
                         "SCEPTR can't use are excluded and reported. cdr3: cdr3_only variant, "
                         "CDR3B alone -- use for the V-gene benchmark, where V as input would "
                         "make the same-V vs diff-V check circular.")
    ap.add_argument("--allow-report-fallback", action="store_true",
                    help="use report.tsv for people with no cdr3.out. Off by default: report.tsv "
                         "has no CDR3_score, so those people skip the quality filter and the "
                         "pool becomes a mix of filtered and unfiltered people.")
    ap.add_argument("--tag", default=None,
                    help="label for output files (default: cohort file stem + sceptr input), so "
                         "runs never overwrite each other's summaries or embeddings")
    ap.add_argument("--chain", default="TRB",
                    help="restrict the pool to one chain before embedding anything (default "
                         "TRB). SCEPTR is beta-chain-specific -- feeding it TRA/TRG/TRD/IGH/"
                         "IGK/IGL (which is ~80%% of recovered CDR3s per the Aug chain "
                         "breakdown) is outside its documented scope and is exactly what "
                         "triggers its own 'doesn't look like a standardised junction' "
                         "warnings. It also confounds the diff-V-gene bucket for BOTH models "
                         "with trivial cross-chain pairs (e.g. IGH vs TRB), which differ for "
                         "reasons that have nothing to do with embedding quality. Pass 'all' "
                         "to disable the filter and use every chain (not recommended for a "
                         "model comparison).")
    ap.add_argument("--outdir", default=os.path.join(
                        os.path.dirname(os.path.abspath(__file__)), "..", "results"),
                    help="de-identified comparison summary only (aggregate stats, no "
                         "research_ids) -- this is the one that's safe to commit")
    ap.add_argument("--local-outdir", default="~/pipeline_outputs/rnaseq/embeddings",
                    help="VM-local: per-CDR3 pool (has real research_ids) + raw embedding "
                         ".npy files -- NEVER commit this, same privacy posture as "
                         "query_overlap_phenotypes.py's pheno cache")
    ap.add_argument("--models", default="sceptr,esmc",
                    help="comma-separated subset of {sceptr,esmc} to run (default: both, for "
                         "the head-to-head comparison). At full-cohort scale ESMC is ~156x "
                         "slower than SCEPTR for a worse same/diff-V-gene gap (2026-09-10 "
                         "500-person comparison, DECISIONS.md) -- pass --models sceptr to "
                         "skip the ~16hr ESMC pass once the comparison question is settled.")
    args = ap.parse_args()
    models = {m.strip().lower() for m in args.models.split(",") if m.strip()}
    bad = models - {"sceptr", "esmc"}
    if bad:
        die(f"--models: unknown model(s) {bad} -- choose from sceptr, esmc")
    want_chain = args.chain.upper()

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
        for col in ("reads", "score"):
            if col not in pool.columns:
                pool[col] = np.nan
    else:
        if not args.cohort:
            die("need either a cohort.tsv (TRUST4 mode) or --from-pool-tsv (open-dataset mode)")
        cohort = pd.read_csv(os.path.expanduser(args.cohort), sep="\t", dtype=str)
        if "research_id" not in cohort.columns:
            die(f"{args.cohort} has no research_id column -- expected build_rnaseq_cohort.py output")

        n_missing, n_no_chain, sources = [], 0, {"cdr3.out": 0, "report.tsv": 0}
        all_rows = []
        for i, rid in enumerate(cohort["research_id"], 1):
            df, source = load_person_cdr3s(args.pheno_dir, rid, args.min_score,
                                           args.keep_imputed, args.allow_report_fallback)
            if df is None:
                n_missing.append(rid)
                continue
            sources[source] += 1
            # Chain filter BEFORE the per-person cap: capping first spends the cap on chains
            # about to be thrown away (found live 2026-09-12: lost 83 of 500 people).
            if want_chain != "ALL":
                df = df[df["chain"] == want_chain]
            if len(df) == 0:
                n_no_chain += 1
                continue
            all_rows.append(df.head(args.max_per_person))
            if i % 500 == 0:
                print(f"  loaded {i:,}/{len(cohort):,} people", file=sys.stderr)

        print(f"\nInput sources: {sources['cdr3.out']:,} cdr3.out (quality-filtered), "
              f"{sources['report.tsv']:,} report.tsv (unfiltered fallback)", file=sys.stderr)
        if n_missing:
            print(f"!! {len(n_missing):,} people skipped: no cdr3.out"
                  f"{'' if args.allow_report_fallback else ' (report.tsv fallback is off)'} "
                  f"-- e.g. {n_missing[:3]}", file=sys.stderr)
        if n_no_chain:
            print(f"{n_no_chain:,} people had zero {want_chain} clonotypes after filtering",
                  file=sys.stderr)
        if not all_rows:
            die("no usable CDR3s across the whole cohort -- run repertoire calling first")

        pool = pd.concat(all_rows, ignore_index=True)
        if "ancestry" in cohort.columns:
            pool = pool.merge(cohort[["research_id", "ancestry"]], on="research_id", how="left")

    if args.chain.lower() != "all":
        if "chain" not in pool.columns:
            die(f"--chain {args.chain} requested but the pool has no 'chain' column")
        before = len(pool)
        pool = pool[pool["chain"] == args.chain.upper()].reset_index(drop=True)
        print(f"--chain {args.chain.upper()}: {before:,} -> {len(pool):,} CDR3s "
              f"(dropped {before - len(pool):,} from other chains)", file=sys.stderr)
        if len(pool) == 0:
            die(f"no CDR3s left after restricting to chain={args.chain.upper()} -- check the "
                f"'chain' values actually present, or pass --chain all")

    # V-gene grouping at gene level (allele stripped): allele calls from ~146 bp RNA reads
    # are noisy, and SCEPTR takes IMGT gene symbols anyway.
    pool["v_gene"] = pool["v_gene"].fillna("").astype(str)
    pool["trbv"] = pool["v_gene"].map(_imgt_gene)
    print(f"\n=== {len(pool):,} clonotypes pooled across {pool['research_id'].nunique():,} "
          f"people, {pool.loc[pool['trbv'] != '', 'trbv'].nunique()} distinct V genes ===\n")

    src = args.cohort or args.from_pool_tsv
    tag = args.tag or (f"{os.path.splitext(os.path.basename(src))[0]}"
                       f"_{args.sceptr_input.replace('+', '')}")
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    local_outdir = os.path.abspath(os.path.expanduser(args.local_outdir))
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(local_outdir, exist_ok=True)
    results = []

    def contrast(embs, frame):
        has_v = (frame["trbv"] != "").to_numpy()
        return same_vs_diff_vgene_contrast(embs[has_v], frame["trbv"].to_numpy()[has_v])

    def report(name, sub, embs, secs, v_is_input):
        same, diff = contrast(embs, sub)
        print(f"  {secs:.1f}s wall, {embs.shape[1]}-dim vectors")
        print(f"  same-V-gene cosine sim: {same:.4f}  |  diff-V-gene: {diff:.4f}  "
              f"|  gap: {same - diff:+.4f}"
              f"{'  (V gene is a model input: gap is circular, not a quality check)' if v_is_input else ''}")
        np.save(os.path.join(local_outdir, f"embeddings_{name.lower()}_{tag}.npy"), embs)
        sub[["research_id", "chain", "v_gene", "cdr3aa", "reads", "score"]
            + (["ancestry"] if "ancestry" in sub.columns else [])].to_csv(
            os.path.join(local_outdir, f"pool_{name.lower()}_{tag}.tsv"), sep="\t", index=False)
        results.append((name, "TRBV+CDR3B" if v_is_input else "CDR3B", len(sub),
                        sub["research_id"].nunique(), embs.shape[1], secs, same, diff, v_is_input))

    if "sceptr" in models:
        v_in = args.sceptr_input == "v+cdr3"
        print(f"--- SCEPTR ({'b_sceptr: TRBV + CDR3B' if v_in else 'cdr3_only: CDR3B'}) ---")
        try:
            model = sceptr_model(args.sceptr_input)
            sub = pool
            if v_in:
                if want_chain != "TRB":
                    die("--sceptr-input v+cdr3 needs --chain TRB (b_sceptr is beta-only)")
                good = accepted_trbv(model, pool["trbv"])
                keep = pool["trbv"].isin(good)
                print(f"  SCEPTR accepts {len(good)} of {pool['trbv'].nunique()} TRBV symbols; "
                      f"{(~keep).sum():,} of {len(pool):,} clonotypes "
                      f"({(~keep).mean():.1%}) excluded (no V call, or V not usable)")
                if "ancestry" in pool.columns:
                    by_anc = (~keep).groupby(pool["ancestry"]).mean()
                    print("  excluded fraction by ancestry: "
                          + ", ".join(f"{a} {f:.1%}" for a, f in by_anc.items()))
                sub = pool[keep].reset_index(drop=True)
                df_in = pd.DataFrame({"TRBV": sub["trbv"], "CDR3B": sub["cdr3aa"]})
            else:
                df_in = pd.DataFrame({"CDR3B": sub["cdr3aa"]})
            embs, secs = embed_sceptr(model, df_in)
            report("SCEPTR", sub, embs, secs, v_in)
        except Exception as e:
            print(f"  !! SCEPTR failed: {e}\n  (run with `pixi run python3` from "
                  f"aleix/RNA-seq/ -- sceptr lives in the pixi env)", file=sys.stderr)

    if "esmc" in models:
        print("\n--- ESMC-300M (general protein LM, CDR3 only) ---")
        try:
            embs, secs = embed_esmc(pool["cdr3aa"].tolist())
            report("ESMC-300M", pool, embs, secs, False)
        except Exception as e:
            print(f"  !! ESMC failed: {e}\n  "
                  f"(need: pip install accelerate && pip install "
                  f"git+https://github.com/Biohub/esm.git@v3.4.1 -- PyPI's esm 3.2.3 and the "
                  f"old EsmcForMaskedLM path are both broken upstream, see embed_esmc() docstring)",
                  file=sys.stderr)

    if not results:
        die("every requested model failed -- see errors above")
    summary = pd.DataFrame(results, columns=[
        "model", "input", "n_clonotypes", "n_people", "dim", "wall_seconds",
        "same_vgene_cos_sim", "diff_vgene_cos_sim", "v_is_input"])
    summary["gap"] = summary["same_vgene_cos_sim"] - summary["diff_vgene_cos_sim"]
    out = os.path.join(outdir, f"cdr3_embedding_summary_{tag}.csv")
    summary.to_csv(out, index=False)
    print(f"\n=== Summary (de-identified, safe to commit): {out} ===")
    print(f"Embeddings + aligned pools (VM-local, real research_ids -- never commit): "
          f"{local_outdir}/embeddings_*_{tag}.npy, pool_*_{tag}.tsv")


if __name__ == "__main__":
    main()
