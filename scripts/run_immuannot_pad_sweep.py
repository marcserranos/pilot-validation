#!/usr/bin/env python3
"""Immuannot contig-trim padding sweep (2026-07-21, Marc). Investigates whether the production
100,000bp pad-on-the-whole-4Mb-window trim (run_immuannot_person.py) is far more generous than
Immuannot actually needs. Marc's reasoning: Immuannot works on an already-assembled contig, not
reads -- the read-binning margin logic behind pad100k (EXPERIMENTS.md's SpecImmune-LR sweep) may
simply not apply here. The real constraint is just "does each gene's full genomic span survive the
cut", which is a much tighter bound than a whole 4Mb reference window. Tests that directly, the
same way the SpecHLA truncation sanity check (EXPERIMENTS.md, 2026-07-10) confirmed a padding floor
by deliberately cutting INTO gene clusters as a positive control, not assuming a lever is free
until it's shown to be.

A SEPARATE, throwaway script (Marc: "make a separate script") -- does not touch or replace
run_immuannot_person.py's proven-safe 100kb-on-whole-window default; this only investigates
whether that default can be tightened, informing a possible future change, not applying one.

GENE_CLUSTERS reused from scripts/spechla_pad_helpers.py (GRCh38, Ensembl-derived, already used
for the SpecHLA/SpecImmune padding sweeps) rather than re-deriving.

Levels tested, per haplotype, on ONE fixed person (for direct comparability with the two data
points already measured by run_immuannot_person.py on 2026-07-21 -- naive/whole-contig and
pad100k-on-the-whole-window were both ~6-6.4min/4.1-4.5MB and are NOT re-run here, just reused):
  gene_pad_100k  -- 100,000bp around each of the 4 gene clusters (not the whole 4Mb window --
                    much tighter, since the clusters themselves are only kb-tens_of_kb wide)
  gene_pad_20k   -- intermediate
  gene_pad_5k    -- intermediate, near the SpecHLA-SR floor's order of magnitude
  gene_pad_0     -- "literal gene cut": exactly each cluster's own mapped span, no margin
  gene_infra     -- POSITIVE CONTROL: truncates INTO each cluster around its midpoint (mirrors
                    spechla_pad_helpers.py's cmd_truncated_windows exactly), deliberately smaller
                    than the gene bodies -- expected to degrade, confirms the sweep can actually
                    detect degradation rather than reporting "safe" no matter what.

MECHANISM CAVEAT (read before trusting gene_pad_0/gene_infra numbers): a .paf line's qstart/qend/
tstart/tend describe ONE alignment block's REPORTED span; where exactly a given reference position
falls WITHIN that block is estimated by linear interpolation (interpolate_query_pos), not exact
CIGAR-walked coordinates -- internal indels inside a block aren't accounted for. This is a real
approximation, tightest exactly where it matters most (the narrowest configs). Fine for an
exploratory sweep deciding whether to invest in exact CIGAR-based trimming next; not fine to adopt
gene_pad_0/gene_infra as a new production default without that follow-up if the numbers look good.

Usage (from ~/repos/pilot-validation, inside the specimmune pixi env):
  pixi run -e specimmune -- python3 scripts/run_immuannot_pad_sweep.py <person_id>
      [--mount ~/mnt/aou-controlled] [--immuannot-dir ~/tools/Immuannot]
      [--refdir ~/tools/Immuannot_refdata] [--outroot ~/pipeline_outputs] [--threads 4]
"""
import argparse
import gzip
import os
import sys
import time

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from run_immuannot_person import (  # noqa: E402
    FA_COLS, PAF_COLS, resolve_cols, trim_assembly, run_immuannot,
    parse_gtf, find_immuannot_script, DEFAULT_REGION,
)
from spechla_pad_helpers import GENE_CLUSTERS  # noqa: E402

CLASSICAL_8 = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1"]
INFRA_TRUNC_WIDTH = 2000  # bp, centered on each cluster's midpoint -- deliberately smaller than
                          # any classical gene body (smallest is ~8kb) to positively confirm
                          # degradation, same role pad5k/truncated-windows played for SR/LR sweeps.

CONFIGS = [
    ("gene_pad_100k", "pad", 100_000),
    ("gene_pad_20k", "pad", 20_000),
    ("gene_pad_5k", "pad", 5_000),
    ("gene_pad_0", "pad", 0),
    ("gene_infra", "trunc", INFRA_TRUNC_WIDTH),
]


def interpolate_query_pos(qstart, qend, tstart, tend, strand, ref_pos):
    """Estimates the query (contig) coordinate for a given reference position, by linear
    interpolation within one .paf block's reported span. See module docstring's MECHANISM CAVEAT."""
    ref_pos = min(max(ref_pos, tstart), tend)
    frac = (ref_pos - tstart) / max(1, (tend - tstart))
    return qstart + frac * (qend - qstart) if strand == "+" else qend - frac * (qend - qstart)


def windows_for_config(kind, value):
    """Returns the list of (start, end) GRCh38 windows to test contig overlap/extraction against,
    for one sweep config. 'pad': each gene cluster's own extent +/- value. 'trunc': a `value`-bp
    window centered on each cluster's midpoint, ignoring the cluster's own extent entirely --
    mirrors spechla_pad_helpers.py's cmd_truncated_windows."""
    if kind == "pad":
        return [(max(0, s - value), e + value) for s, e in GENE_CLUSTERS]
    half = value // 2
    return [(max(0, (s + e) // 2 - half), (s + e) // 2 + half) for s, e in GENE_CLUSTERS]


def regions_from_paf_windows(paf_path, chrom, windows):
    """Like run_immuannot_person.regions_from_paf, but against MULTIPLE reference windows (one
    per gene cluster) and using interpolation (not whole-block extraction) so a cut can actually
    be tighter than any single .paf block's own reported span -- necessary for gene_pad_0/
    gene_infra to be meaningfully different from the more generous levels."""
    opener = gzip.open if paf_path.endswith(".gz") else open
    accept_chrom = {chrom, chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"}
    regions, qlens, seen_tnames = {}, {}, set()
    with opener(paf_path, "rt") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qname, qlen, qstart, qend = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
            strand, tname = fields[4], fields[5]
            tstart, tend = int(fields[7]), int(fields[8])
            seen_tnames.add(tname)
            if tname not in accept_chrom:
                continue
            for wstart, wend in windows:
                if tend <= wstart or tstart >= wend:
                    continue
                lo = interpolate_query_pos(qstart, qend, tstart, tend, strand, max(tstart, wstart))
                hi = interpolate_query_pos(qstart, qend, tstart, tend, strand, min(tend, wend))
                q_lo, q_hi = (min(lo, hi), max(lo, hi))
                q_lo, q_hi = max(0, int(q_lo)), min(qlen, int(round(q_hi)))
                qlens[qname] = qlen
                if qname in regions:
                    s, e = regions[qname]
                    regions[qname] = (min(s, q_lo), max(e, q_hi))
                else:
                    regions[qname] = (q_lo, q_hi)
    return regions, qlens, seen_tnames


def run_one_config(pid, hap, config_name, kind, value, fa_path, paf_path,
                    person_dir, immuannot_script, refdir, threads, region):
    t0 = time.perf_counter()
    chrom = region.split(":")[0]
    windows = windows_for_config(kind, value)
    regions, qlens, seen_tnames = regions_from_paf_windows(paf_path, chrom, windows)

    if not regions:
        print(f"    {config_name}: 0 contigs overlap any gene-cluster window -- "
              f"tnames seen: {sorted(seen_tnames)[:10]}", file=sys.stderr)
        return None

    targets = [f"{c}:{s + 1}-{e}" for c, (s, e) in regions.items()]
    trimmed_bp = sum(e - s for s, e in regions.values())
    lookup_t1 = time.perf_counter()

    trimmed_fa = os.path.join(person_dir, f"{hap}.{config_name}.trimmed.fa")
    ok, err = trim_assembly(fa_path, targets, trimmed_fa)
    trim_t2 = time.perf_counter()
    if not ok:
        print(f"    {config_name}: SKIP -- trim failed: {err}", file=sys.stderr)
        return None
    trimmed_mb = os.path.getsize(trimmed_fa) / 1e6

    outprefix = os.path.join(person_dir, f"{hap}.{config_name}")
    gtf = run_immuannot(immuannot_script, refdir, trimmed_fa, outprefix, threads)
    t3 = time.perf_counter()
    calls = parse_gtf(gtf, show_raw_sample=False) if gtf else {}
    classical_called = [g for g in CLASSICAL_8 if g in calls]
    classical_missing = [g for g in CLASSICAL_8 if g not in calls]

    row = {
        "person_id": pid, "hap": hap, "config": config_name, "kind": kind, "value": value,
        "n_contigs": len(regions), "trimmed_mb": round(trimmed_mb, 2),
        "lookup_seconds": round(lookup_t1 - t0, 2), "trim_seconds": round(trim_t2 - lookup_t1, 2),
        "immuannot_seconds": round(t3 - trim_t2, 2), "config_total_seconds": round(t3 - t0, 2),
        "n_genes_total": len(calls), "n_classical8_called": len(classical_called),
        "classical8_missing": ",".join(classical_missing) if classical_missing else "",
    }
    print(f"    {config_name}: {trimmed_mb:.2f} MB, {(t3-t0)/60:.1f} min total, "
          f"{len(calls)} genes ({len(classical_called)}/8 classical) "
          f"{'MISSING: ' + row['classical8_missing'] if classical_missing else ''}",
          file=sys.stderr)
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("person_id")
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    ap.add_argument("--immuannot-dir", default=os.path.expanduser("~/tools/Immuannot"))
    ap.add_argument("--refdir", default=os.path.expanduser("~/tools/Immuannot_refdata"))
    ap.add_argument("--outroot", default=os.path.expanduser("~/pipeline_outputs"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--region", default=DEFAULT_REGION)
    args = ap.parse_args()

    lr_path = os.path.join(args.mount, "v9/wgs/long_read/manifest.tsv")
    if not os.path.exists(lr_path):
        print(f"FATAL: manifest not found: {lr_path} -- remount first.", file=sys.stderr)
        sys.exit(1)
    immuannot_script = find_immuannot_script(args.immuannot_dir)
    if immuannot_script is None:
        print(f"FATAL: immuannot.sh not found under {args.immuannot_dir}.", file=sys.stderr)
        sys.exit(1)

    lr = pd.read_csv(lr_path, sep="\t", dtype=str)
    pid = args.person_id
    fa_paths = resolve_cols(lr, pid, args.mount, FA_COLS)
    paf_paths = resolve_cols(lr, pid, args.mount, PAF_COLS)
    if fa_paths is None or any(v is None for v in fa_paths.values()):
        print(f"FATAL: {pid} missing assembly FASTA -- pick a person confirmed to have one "
              f"(e.g. one already in immuannot_calls.tsv).", file=sys.stderr)
        sys.exit(1)
    if any(v is None for v in paf_paths.values()):
        print(f"FATAL: {pid} missing a .paf -- this sweep needs it (no whole-contig fallback "
              f"here, that's a different, already-measured data point).", file=sys.stderr)
        sys.exit(1)

    person_dir = os.path.join(args.outroot, str(pid), "immuannot_pad_sweep")
    os.makedirs(person_dir, exist_ok=True)

    print(f"=== Immuannot pad sweep: person {pid} ===", file=sys.stderr)
    print(f"(reusing already-measured naive/whole-contig and pad100k-on-whole-window numbers "
          f"from immuannot_timing.tsv -- not re-run here)", file=sys.stderr)

    all_rows = []
    run_t0 = time.perf_counter()
    for hap in ("hap1", "hap2"):
        fa_path = os.path.join(args.mount, fa_paths[hap])
        paf_path = os.path.join(args.mount, paf_paths[hap])
        print(f"\n--- {hap} ---", file=sys.stderr)
        for config_name, kind, value in CONFIGS:
            row = run_one_config(pid, hap, config_name, kind, value, fa_path, paf_path,
                                  person_dir, immuannot_script, args.refdir, args.threads,
                                  args.region)
            if row:
                all_rows.append(row)

    run_total = time.perf_counter() - run_t0
    print(f"\n=== Sweep total: {run_total/60:.1f} min for {len(all_rows)} config runs ===",
          file=sys.stderr)

    out_path = os.path.join(args.outroot, "immuannot_pad_sweep.tsv")
    pd.DataFrame(all_rows).to_csv(out_path, sep="\t", index=False)
    print(f"Results: {out_path}", file=sys.stderr)

    print("\n=== Comparison (this sweep's new configs only -- see immuannot_timing.tsv for the "
          "already-measured naive/pad100k-whole-window baseline) ===", file=sys.stderr)
    df = pd.DataFrame(all_rows)
    if not df.empty:
        print(df[["hap", "config", "trimmed_mb", "immuannot_seconds", "n_genes_total",
                   "n_classical8_called", "classical8_missing"]].to_string(index=False),
              file=sys.stderr)


if __name__ == "__main__":
    main()
