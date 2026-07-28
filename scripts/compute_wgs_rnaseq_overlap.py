#!/usr/bin/env python3
"""Computes and visualizes participant overlap between AoU's srWGS, lrWGS, and RNA-seq
(v9 multiomics) cohorts -- the gating question for the RNA-seq immune-repertoire scoping
thread (scoping memo delivered 2026-07-25; not yet logged to this repo).

Three manifests, v9 CDR (paths confirmed LIVE on the VM 2026-07-25 by browsing the mount --
the official "How the All of Us Genomic Data are Organized v9" PDF documents no RNA-seq
manifest at all, despite one existing; see reference/AOU_DATA_ACCESS_NOTES.md section 7):

  - srWGS CRAM manifest    v9/wgs/cram/manifest.csv          (person_id, cram_uri, ...)
  - lrWGS manifest         v9/wgs/long_read/manifest.tsv     (research_id, ..., grch38_bam, ...)
  - RNA-seq manifest       v9/multiomics/rnaseq/manifest.tsv (sampleid, research_id,
                                                               markduplicates_bam_file_path, ...)
  (RNA-seq manifest row count, 8980, exactly matches the PDF's stated cohort size -- a real
  cross-check that this is the right file.)

srWGS is trusted at manifest-presence (535k rows -- an existence check at that scale isn't
affordable over gcsfuse, and no quirk has ever surfaced a dead srWGS CRAM path, unlike lrWGS).
lrWGS and RNA-seq are both existence-checked against the mount before counting a person as
"has this data": lrWGS because ENVIRONMENT.md quirk #13 already proved manifest presence !=
a real usable file there; RNA-seq gets the same treatment on general principle, since AoU's own
docs don't even acknowledge this manifest exists, and the whole point of this script is to not
repeat the "trusted an inferred/undocumented path" mistake logged elsewhere in this project.

Writes an id-level TSV (person_id + 3 booleans) to ~/pipeline_outputs/ -- not the repo, same
IDs-not-in-repo posture as build_experiment_d_cohort.py's cohort.tsv -- and a schematic
(NOT area-proportional) 3-circle Venn PNG, which is counts-only and safe to bring into the repo.

Usage (needs pandas + matplotlib -- the spechla env has both):
  pixi run -e spechla -- python3 scripts/compute_wgs_rnaseq_overlap.py
"""
import argparse
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

BUCKET_PREFIX = "gs://vwb-aou-datasets-controlled/"
CRAM_MANIFEST = "v9/wgs/cram/manifest.csv"
LR_MANIFEST = "v9/wgs/long_read/manifest.tsv"
RNASEQ_MANIFEST = "v9/multiomics/rnaseq/manifest.tsv"


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def require_cols(df, cols, label):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        die(f"{label}: expected column(s) {missing} not found. Actual columns: "
            f"{list(df.columns)}. AoU may have renamed something -- update the constants "
            f"at the top of this script.")


def strip_bucket(uri):
    """gs://vwb-aou-datasets-controlled/pooled/... -> pooled/... (mount-relative)."""
    if not isinstance(uri, str) or not uri:
        return None
    if "fc-aou-datasets-controlled" in uri:
        return None  # wrong-bucket trap, ENVIRONMENT.md
    if uri.startswith(BUCKET_PREFIX):
        return uri[len(BUCKET_PREFIX):]
    if uri.startswith("gs://"):
        return None
    return uri


def check_exists(mount, rel_path):
    if rel_path is None:
        return False
    return os.path.exists(os.path.join(mount, rel_path))


def existence_check_column(df, uri_col, mount, label, workers=16):
    """Threaded os.path.exists() over a column of gs:// URIs -- FUSE I/O-bound, so threading
    helps despite the GIL (same rationale as ENVIRONMENT.md quirk #22's fix)."""
    rel_paths = list(df[uri_col].map(strip_bucket))
    n = len(rel_paths)
    print(f"  Existence-checking {n} {label} paths against the mount "
          f"(a few minutes over FUSE -- not hung)...", file=sys.stderr)
    results = [False] * n
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(check_exists, mount, p): i for i, p in enumerate(rel_paths)}
        done = 0
        for fut in as_completed(futures):
            i = futures[fut]
            results[i] = fut.result()
            done += 1
            if done % 2000 == 0:
                print(f"    ...checked {done}/{n}", file=sys.stderr)
    n_ok = sum(results)
    print(f"  {label}: {n_ok}/{n} resolve to a real file", file=sys.stderr)
    return pd.Series(results, index=df.index)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    ap.add_argument("--out-dir", default=os.path.expanduser("~/pipeline_outputs/rnaseq_overlap"))
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    cram_path = os.path.join(args.mount, CRAM_MANIFEST)
    lr_path = os.path.join(args.mount, LR_MANIFEST)
    rna_path = os.path.join(args.mount, RNASEQ_MANIFEST)
    for p in (cram_path, lr_path, rna_path):
        if not os.path.exists(p):
            die(f"manifest not found: {p} -- is the gcsfuse mount up? Remount + `ls`-verify "
                f"before running this (ENVIRONMENT.md quirk #11/#14).")

    # --- srWGS: trust manifest presence -- 535k rows, too large to existence-check over FUSE,
    # and no quirk has ever surfaced a dead srWGS CRAM path (unlike lrWGS, see quirk #13) ---
    print(f"Reading srWGS CRAM manifest: {cram_path}", file=sys.stderr)
    cram = pd.read_csv(cram_path, dtype=str, usecols=lambda c: c in ("person_id", "cram_uri"))
    require_cols(cram, ["person_id", "cram_uri"], "srWGS CRAM manifest")
    srwgs_ids = set(cram["person_id"].dropna().unique())
    print(f"  srWGS: {len(srwgs_ids)} unique person_id (manifest presence, not existence-checked "
          f"-- see docstring)", file=sys.stderr)

    # --- lrWGS: existence-checked (ENVIRONMENT.md quirk #13 -- manifest presence != real file) ---
    print(f"Reading lrWGS manifest: {lr_path}", file=sys.stderr)
    lr = pd.read_csv(lr_path, sep="\t", dtype=str)
    require_cols(lr, ["research_id", "grch38_bam"], "lrWGS manifest")
    ok_primary = existence_check_column(lr, "grch38_bam", args.mount, "lrWGS grch38_bam", args.workers)
    if "grch38_haplotagged_bam" in lr.columns:
        ok_secondary = existence_check_column(lr, "grch38_haplotagged_bam", args.mount,
                                               "lrWGS grch38_haplotagged_bam", args.workers)
        lr_resolves = ok_primary | ok_secondary
    else:
        lr_resolves = ok_primary
    lrwgs_ids = set(lr.loc[lr_resolves, "research_id"].dropna().unique())
    print(f"  lrWGS: {len(lrwgs_ids)} unique research_id with a confirmed-usable BAM "
          f"(of {lr['research_id'].nunique()} in the manifest)", file=sys.stderr)

    # --- RNA-seq: existence-checked on general principle -- AoU's own docs don't even
    # acknowledge this manifest exists, so trust nothing about it beyond what's verified live ---
    print(f"Reading RNA-seq manifest: {rna_path}", file=sys.stderr)
    rna = pd.read_csv(rna_path, sep="\t", dtype=str)
    require_cols(rna, ["research_id", "markduplicates_bam_file_path"], "RNA-seq manifest")
    rna_resolves = existence_check_column(rna, "markduplicates_bam_file_path", args.mount,
                                           "RNA-seq markduplicates_bam", args.workers)
    rnaseq_ids = set(rna.loc[rna_resolves, "research_id"].dropna().unique())
    print(f"  RNA-seq: {len(rnaseq_ids)} unique research_id with a confirmed-usable BAM "
          f"(of {rna['research_id'].nunique()} in the manifest)", file=sys.stderr)

    # --- overlap ---
    all_ids = sorted(srwgs_ids | lrwgs_ids | rnaseq_ids)
    df = pd.DataFrame({"person_id": all_ids})
    df["has_srwgs"] = df["person_id"].isin(srwgs_ids)
    df["has_lrwgs"] = df["person_id"].isin(lrwgs_ids)
    df["has_rnaseq"] = df["person_id"].isin(rnaseq_ids)

    out_tsv = os.path.join(args.out_dir, "overlap_ids.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nWrote {len(df)} people (union of all three) to {out_tsv}", file=sys.stderr)

    def region(s, l, r):
        mask = (df["has_srwgs"] == s) & (df["has_lrwgs"] == l) & (df["has_rnaseq"] == r)
        return int(mask.sum())

    only_sr = region(True, False, False)
    only_lr = region(False, True, False)
    only_rna = region(False, False, True)
    sr_lr = region(True, True, False)
    sr_rna = region(True, False, True)
    lr_rna = region(False, True, True)
    all3 = region(True, True, True)

    print("\n=== Overlap (confirmed-usable-file people, not raw manifest rows) ===", file=sys.stderr)
    print(f"  srWGS total:   {len(srwgs_ids)}", file=sys.stderr)
    print(f"  lrWGS total:   {len(lrwgs_ids)}", file=sys.stderr)
    print(f"  RNA-seq total: {len(rnaseq_ids)}", file=sys.stderr)
    print(f"  srWGS only:      {only_sr}", file=sys.stderr)
    print(f"  lrWGS only:      {only_lr}", file=sys.stderr)
    print(f"  RNA-seq only:    {only_rna}", file=sys.stderr)
    print(f"  srWGS + lrWGS:   {sr_lr}", file=sys.stderr)
    print(f"  srWGS + RNA-seq: {sr_rna}", file=sys.stderr)
    print(f"  lrWGS + RNA-seq: {lr_rna}", file=sys.stderr)
    print(f"  all three:       {all3}", file=sys.stderr)
    rna_and_lr = lr_rna + all3
    rna_and_either_wgs = sr_rna + lr_rna + all3
    print(f"\n  RNA-seq people who ALSO have lrWGS (the number this scoping thread hinges on -- "
          f"can this project's own HLA calls be linked to RNA-seq at all): {rna_and_lr} of "
          f"{len(rnaseq_ids)} RNA-seq participants", file=sys.stderr)
    print(f"  RNA-seq people who have EITHER WGS type: {rna_and_either_wgs} of "
          f"{len(rnaseq_ids)} RNA-seq participants", file=sys.stderr)

    # --- schematic 3-circle Venn -- NOT area-proportional (no matplotlib-venn dependency
    # added; deliberately labeled as schematic so circle size is never read as meaningful) ---
    fig, ax = plt.subplots(figsize=(8, 7))
    radius = 0.32
    centers = {"srWGS": (0.35, 0.60), "lrWGS": (0.65, 0.60), "RNA-seq": (0.50, 0.32)}
    totals = {"srWGS": len(srwgs_ids), "lrWGS": len(lrwgs_ids), "RNA-seq": len(rnaseq_ids)}
    colors = {"srWGS": "#4C72B0", "lrWGS": "#55A868", "RNA-seq": "#C44E52"}
    label_offset = {"srWGS": (-0.22, 0.30), "lrWGS": (0.22, 0.30), "RNA-seq": (0.0, -0.40)}

    for name, (cx, cy) in centers.items():
        ax.add_patch(Circle((cx, cy), radius, alpha=0.32, color=colors[name],
                             ec=colors[name], lw=2))
        ox, oy = label_offset[name]
        ax.text(cx + ox, cy + oy, f"{name}\n(n={totals[name]})", ha="center", va="center",
                fontsize=11, fontweight="bold", color=colors[name])

    region_labels = {
        (0.24, 0.68): only_sr,
        (0.76, 0.68): only_lr,
        (0.50, 0.14): only_rna,
        (0.50, 0.72): sr_lr,
        (0.38, 0.40): sr_rna,
        (0.62, 0.40): lr_rna,
        (0.50, 0.50): all3,
    }
    for (x, y), n in region_labels.items():
        ax.text(x, y, str(n), ha="center", va="center", fontsize=13, fontweight="bold")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("AoU v9: srWGS x lrWGS x RNA-seq participant overlap\n"
                  "(schematic layout -- circle size is NOT proportional to count)",
                  fontsize=12)
    fig_path = os.path.join(args.out_dir, "venn_srwgs_lrwgs_rnaseq.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    print(f"\nWrote {fig_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
