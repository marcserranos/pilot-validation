#!/usr/bin/env python3
"""Adapter: turn public VDJdb TCR data into the exact pool-TSV schema embed_cdr3s.py /
embed_catelmo.py / compare_embeddings.py already expect (research_id, chain, v_gene, cdr3aa,
score) -- so the SAME three embedding scripts run unmodified against open data, structurally
identical to a TRUST4-derived pool. Built for right now: the AoU workspace's gcsfuse mount is
blocked (VPC-SC perimeter denial, unresolved), so this lets embedding-model work continue
without waiting on that.

SOURCE: VDJdb (github.com/antigenomics/vdjdb-db), a curated, antigen-specificity-annotated
TCR database -- public, no auth, no controlled-tier restrictions. Schema confirmed live
2026-09-07 against release 2026-06-03: `gene, cdr3, species, antigen.epitope, antigen.gene,
antigen.species, complex.id, v.segm, j.segm, mhc.a, mhc.b, mhc.class, reference.id,
vdjdb.score, TCR_hash, j.start, v.end` (the `.slim.txt` file inside the release zip).

HONEST CAVEAT, read before trusting any comparison built on this: SCEPTR and catELMo were
both very likely trained on data overlapping this exact kind of public TCR database
(catELMo's own README says its training data is "collected from ImmunoSEQ," 4M+ sequences;
SCEPTR draws from similar public sources per its paper). So this is NOT a clean held-out
benchmark for those two -- there is real train/test overlap risk. ESMC is the one model
where this genuinely is fair: it was trained on ordinary proteins (UniRef/MGnify/JGI), never
TCR data, so VDJdb is honestly out-of-distribution for it either way. The same-V-gene-vs-
diff-V-gene structure check downstream is a weak, generic signal (not exact-sequence
retrieval), so it isn't worthless for SCEPTR/catELMo even with overlap risk -- just don't
call it a rigorous benchmark for those two in anything written up.

Usage:
  python3 prep_vdjdb_pool.py [--url <release-zip-url>] [--species HomoSapiens] [--gene TRB]
      [--max-rows 2000] [--out ~/pipeline_outputs/rnaseq/embeddings/vdjdb_pool.tsv]

  Then:  pixi run python3 embed_cdr3s.py --from-pool-tsv <out>
         conda run -n catELMo python3 embed_catelmo.py <out>
         python3 compare_embeddings.py <out>
"""
import argparse
import io
import os
import sys
import zipfile

import pandas as pd

AA_VALID = set("ACDEFGHIKLMNPQRSTVWY")
DEFAULT_URL = ("https://github.com/antigenomics/vdjdb-db/releases/download/"
               "2026-06-03-ZENODO/vdjdb-2026-06-03.zip")


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--url", default=DEFAULT_URL,
                    help="VDJdb release zip URL -- check github.com/antigenomics/vdjdb-db/"
                         "releases/latest if this one 404s (releases are dated, not stable)")
    ap.add_argument("--local-zip", default=None,
                    help="use an already-downloaded zip instead of fetching --url again")
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--gene", default="TRB",
                    help="TRB (beta) matches SCEPTR's beta-only convention already used "
                         "elsewhere in this project. TRA also exists in VDJdb if wanted.")
    ap.add_argument("--max-rows", type=int, default=2000,
                    help="cap after dedup, to keep the comparison the same rough scale as a "
                         "~10-person TRUST4 pool rather than VDJdb's full tens of thousands")
    ap.add_argument("--out", default=os.path.expanduser(
                        "~/pipeline_outputs/rnaseq/embeddings/vdjdb_pool.tsv"))
    args = ap.parse_args()

    if args.local_zip:
        with open(os.path.expanduser(args.local_zip), "rb") as f:
            raw = f.read()
    else:
        import urllib.request
        print(f"Downloading {args.url} ...", file=sys.stderr)
        with urllib.request.urlopen(args.url) as resp:
            raw = resp.read()

    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        slim_names = [n for n in zf.namelist() if n.endswith("vdjdb.slim.txt")]
        if not slim_names:
            die(f"no *vdjdb.slim.txt inside the zip -- contents: {zf.namelist()[:15]}")
        with zf.open(slim_names[0]) as f:
            df = pd.read_csv(f, sep="\t")

    print(f"Loaded {len(df):,} raw VDJdb rows.", file=sys.stderr)
    need = {"gene", "cdr3", "species", "v.segm"}
    missing = need - set(df.columns)
    if missing:
        die(f"VDJdb schema changed -- missing column(s) {missing}. "
            f"Actual columns: {list(df.columns)}")

    df = df[(df["gene"] == args.gene) & (df["species"] == args.species)]
    df = df.dropna(subset=["cdr3", "v.segm"])
    df = df[df["cdr3"].astype(str).apply(
        lambda s: len(s) >= 5 and all(c in AA_VALID for c in s))]
    # v.segm carries allele suffixes like "TRBV7-2*01" -- keep the full allele string,
    # consistent with what TRUST4's own V column looks like (same convention, not stripped).
    df = df.drop_duplicates(subset=["v.segm", "cdr3"])
    print(f"{len(df):,} unique, valid {args.gene}/{args.species} CDR3s after filtering.",
          file=sys.stderr)

    if len(df) > args.max_rows:
        df = df.sample(n=args.max_rows, random_state=0)
        print(f"Sampled down to {args.max_rows:,} (--max-rows) for scale-matched comparison.",
              file=sys.stderr)

    pool = pd.DataFrame({
        "research_id": "vdjdb_public",       # not a real person -- reference data, one label
        "chain": args.gene,
        "v_gene": df["v.segm"].values,
        "cdr3aa": df["cdr3"].values,
        "score": pd.NA,                      # not TRUST4's CDR3_score -- not comparable, left blank
    })

    outdir = os.path.dirname(os.path.expanduser(args.out))
    os.makedirs(outdir, exist_ok=True)
    pool.to_csv(os.path.expanduser(args.out), sep="\t", index=False)
    print(f"\nWrote {len(pool):,} rows, {pool['v_gene'].nunique()} distinct V genes, to "
          f"{args.out}")
    print("\nThis file is public reference data, not participant data -- fine to keep, "
          "move around, or eventually commit if useful, unlike the TRUST4-derived pool.")
    print("\nNext:")
    print(f"  pixi run python3 embed_cdr3s.py --from-pool-tsv {args.out}")
    print(f"  conda run -n catELMo python3 embed_catelmo.py {args.out}")
    print(f"  python3 compare_embeddings.py {args.out}")


if __name__ == "__main__":
    main()
