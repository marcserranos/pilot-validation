#!/usr/bin/env python3
"""Run Immuannot (github.com/YingZhou001/Immuannot) on one or more AoU people.

Immuannot types HLA alleles from a PHASED ASSEMBLY contig (one haplotype's collapsed consensus
FASTA), not reads/BAM -- see DECISIONS.md "Assembly-based HLA typing". Only people on the
revio/sequel2e/sequel2 platforms have this data (reports/lr_data_census/README.md); ont-r10.4.1/
ont-r9.4.1 do not and will be skipped with a clear reason, not silently dropped.

TRIM STEP (2026-07-21, Marc: "trim the zone we need, not ridiculous paddings" before the 60-person
pilot run). A whole diploid assembly (both haplotypes, genome-wide) is far more than Immuannot
needs to type 8 classical HLA genes on chr6 -- but assembly contigs are NOT reference-coordinate
strings (they're de novo sequence, arbitrary contig IDs), so we can't naively slice
"chr6:29,500,000-33,500,000" out of them the way slice_and_fastq.sh does for aligned BAMs. Instead
we use each haplotype's OWN assembly-to-hg38 alignment file (assembly_hap{1,2}_aln2_hg38_bam --
confirmed present alongside the raw assembly FASTA for revio/sequel2e per the census) to find
which contig(s) actually overlap the HLA region in GRCh38 coordinates, then samtools-faidx just
those WHOLE contigs out of the full assembly. This is deliberately the simple version: take the
whole matching contig(s), not a sub-range within them -- no CIGAR/indel coordinate math, no risk of
clipping a gene at a padding boundary. Still a large cut versus the full genome (typically one
chromosome-arm-scale contig instead of the full ~1.5-3GB diploid assembly). If per-person timing
(see below) shows this isn't enough, true sub-range trimming via the companion .paf file is the
next lever -- not built yet, deliberately, per Marc's "try one/two, evaluate, then decide" plan.

Per person: resolves assembly_hap1_fa/assembly_hap2_fa AND assembly_hap1_aln2_hg38_bam/
assembly_hap2_aln2_hg38_bam mount-relative paths from the v9 lrWGS manifest (existence-checked
directly against the mount -- same discipline as build_experiment_d_cohort.py's
find_existing_bam(), never a path-pattern guess -- Marc's "I'm sure they're all revio" is treated
as a hypothesis to verify, not a fact, per this project's own repeated lesson about that). Then,
per haplotype: finds overlapping contig(s), trims, runs immuannot.sh (the tool's own documented
mode -- it does not accept both haplotypes in one run), verifies the expected *.gtf.gz output
actually exists (never trust exit code alone -- ENVIRONMENT.md quirks #17/#18), times every stage,
and writes a per-person, per-gene 2-allele call table plus a timing log for the scale-up decision.

GTF-PARSING NOTE: verified against the real reference/README_Immuannot.md now in this repo (not
just a web summary) -- "consensus"/"allele" attribute (not "template_allele") is correct, and it
lives on the GTF's "transcript" feature row, which also carries the inherited gene_id/gene_name
attribute on the same line, so the same-line regex match here is correct as designed.

Usage (from ~/repos/pilot-validation, inside the specimmune pixi env -- it already has both
minimap2, for immuannot.sh, and samtools>=1.10, for the trim step -- see pixi.toml):
  pixi run -e specimmune -- python3 scripts/run_immuannot_person.py <person_id> [<person_id> ...] \\
      [--mount ~/mnt/aou-controlled] [--immuannot-dir ~/tools/Immuannot] \\
      [--refdir ~/tools/Immuannot_refdata] [--outroot ~/pipeline_outputs] [--threads 4] \\
      [--region chr6:29500000-33500000] [--time-budget-min 30]

First run: pass exactly ONE person_id and read the printed timing summary before deciding whether
to pass all 60 cohort ids in one invocation.
"""
import argparse
import gzip
import os
import re
import subprocess
import sys
import time

import pandas as pd

BUCKET_PREFIX = "gs://vwb-aou-datasets-controlled/"
LR_MANIFEST = "v9/wgs/long_read/manifest.tsv"
FA_COLS = {"hap1": "assembly_hap1_fa", "hap2": "assembly_hap2_fa"}
ALN_COLS = {"hap1": "assembly_hap1_aln2_hg38_bam", "hap2": "assembly_hap2_aln2_hg38_bam"}
DEFAULT_REGION = "chr6:29500000-33500000"  # the project's standing HLA window (ENVIRONMENT.md)


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def strip_bucket(uri):
    """Mirrors build_experiment_d_cohort.py's strip_bucket -- same wrong-bucket trap guard."""
    if not isinstance(uri, str) or not uri:
        return None
    if "fc-aou-datasets-controlled" in uri:
        return None
    if uri.startswith(BUCKET_PREFIX):
        return uri[len(BUCKET_PREFIX):]
    if uri.startswith("gs://"):
        return None
    return uri


def resolve_cols(lr, person_id, mount, cols):
    """Returns {'hap1': rel_path_or_None, 'hap2': rel_path_or_None} for the given column map --
    existence-checked directly against the mount, never assumed from a non-null manifest cell.
    Returns None (not a dict) if the person has no manifest row at all."""
    rows = lr[lr["research_id"] == str(person_id)]
    if rows.empty:
        return None
    out = {}
    for hap, col in cols.items():
        if col not in lr.columns:
            out[hap] = None
            continue
        found = None
        for v in rows[col].dropna():
            rel = strip_bucket(v)
            if rel and os.path.exists(os.path.join(mount, rel)):
                found = rel
                break
        out[hap] = found
    return out


def contigs_overlapping_region(bam_path, region):
    """Returns (sorted unique contig names, stderr) via samtools view on the hap's
    assembly-to-hg38 alignment BAM. contigs is None (not empty) on a real samtools failure --
    callers must tell that apart from "ran fine, zero contigs overlap" (unexpected either way,
    but a different problem: a broken index/region string vs. a genuinely empty result)."""
    proc = subprocess.run(["samtools", "view", bam_path, region], capture_output=True, text=True)
    if proc.returncode != 0:
        return None, proc.stderr
    contigs = sorted({line.split("\t")[0] for line in proc.stdout.splitlines() if line.strip()})
    return contigs, proc.stderr


def trim_assembly(fa_path, contigs, out_fa_path):
    """samtools faidx pulls the WHOLE named contig(s) (not a sub-range -- see module docstring)
    out of the full assembly FASTA. Builds a .fai alongside fa_path on first use if missing."""
    with open(out_fa_path, "w") as out:
        proc = subprocess.run(["samtools", "faidx", fa_path] + contigs,
                               stdout=out, stderr=subprocess.PIPE, text=True)
    ok = proc.returncode == 0 and os.path.getsize(out_fa_path) > 0
    return ok, proc.stderr


def run_immuannot(immuannot_dir, refdir, contig_path, outprefix, threads):
    script = os.path.join(immuannot_dir, "scripts", "immuannot.sh")
    cmd = ["bash", script, "-c", contig_path, "-r", refdir, "-o", outprefix, "-t", str(threads)]
    print(f"    Running: {' '.join(cmd)}", file=sys.stderr)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    expected = outprefix + ".gtf.gz"
    if not os.path.exists(expected):
        # Never trust exit code alone (ENVIRONMENT.md quirks #17/#18) -- surface the real logs.
        print(f"    WARNING: expected output {expected} missing after run "
              f"(exit code {proc.returncode}). stdout/stderr follow:", file=sys.stderr)
        print(proc.stdout[-2000:], file=sys.stderr)
        print(proc.stderr[-2000:], file=sys.stderr)
        return None
    return expected


def parse_gtf(gtf_gz_path, show_raw_sample=True):
    """Extracts (gene, allele) pairs from Immuannot's gtf.gz. Prints the first record's raw
    attributes so the field spelling can be eyeballed against reference/README_Immuannot.md."""
    calls = {}
    printed_sample = False
    with gzip.open(gtf_gz_path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = fields[8]
            if show_raw_sample and not printed_sample:
                print(f"    [sanity check] first GTF attribute string seen: {attrs}", file=sys.stderr)
                printed_sample = True
            gene_m = re.search(r'gene_name "([^"]+)"', attrs) or re.search(r'gene_id "([^"]+)"', attrs)
            allele_m = re.search(r'consensus "([^"]+)"', attrs) or re.search(r'allele "([^"]+)"', attrs)
            if gene_m and allele_m:
                calls[gene_m.group(1)] = allele_m.group(1)
    return calls


def process_person(pid, lr, args):
    """Returns (gene_rows, timing_rows) for one person. Every stage is timestamped separately
    (manifest resolution, contig lookup, trim, immuannot.sh) and printed as it completes -- Marc's
    2026-07-21 ask, to see exactly which stage the "BAM things" (contig lookup / trim) actually
    cost versus immuannot.sh itself, rather than one bundled number. A timing row is recorded for
    every haplotype attempt, including ones that get skipped partway through, so a slow-then-failing
    stage is still visible in immuannot_timing.tsv, not silently dropped."""
    print(f"\n=== Person {pid} ===", file=sys.stderr)
    person_t0 = time.perf_counter()

    resolve_t0 = time.perf_counter()
    fa_paths = resolve_cols(lr, pid, args.mount, FA_COLS)
    aln_paths = resolve_cols(lr, pid, args.mount, ALN_COLS)
    resolve_seconds = time.perf_counter() - resolve_t0
    print(f"  manifest path resolution: {resolve_seconds:.2f}s", file=sys.stderr)

    if fa_paths is None:
        print(f"  SKIP: no manifest row for {pid} at all.", file=sys.stderr)
        return [], []
    missing_fa = [h for h, p in fa_paths.items() if p is None]
    missing_aln = [h for h, p in aln_paths.items() if p is None]
    if missing_fa or missing_aln:
        print(f"  SKIP: missing assembly FASTA {missing_fa or 'none'} / "
              f"aln-to-hg38 BAM {missing_aln or 'none'} -- this person's platform likely lacks "
              f"assembly data (only revio/sequel2e/sequel2 do -- reports/lr_data_census/README.md), "
              f"contradicting the 'they're all revio' assumption for this specific person. "
              f"Not treated as a bug -- flag which person_ids hit this so we know what fraction "
              f"of the 60 are actually usable.", file=sys.stderr)
        return [], []

    person_dir = os.path.join(args.outroot, str(pid), "immuannot_output")
    os.makedirs(person_dir, exist_ok=True)

    gene_calls = {}
    timing_rows = []
    for hap in ("hap1", "hap2"):
        row = {"person_id": pid, "hap": hap, "resolve_seconds": round(resolve_seconds, 2),
               "n_contigs": None, "trimmed_mb": None,
               "contig_lookup_seconds": None, "trim_seconds": None, "immuannot_seconds": None,
               "hap_total_seconds": None}
        hap_t0 = time.perf_counter()
        aln_path = os.path.join(args.mount, aln_paths[hap])
        fa_path = os.path.join(args.mount, fa_paths[hap])

        # --- Stage 1: contig lookup (samtools view on the remote/mounted aln-to-hg38 BAM) ---
        contigs, err = contigs_overlapping_region(aln_path, args.region)
        lookup_t1 = time.perf_counter()
        lookup_seconds = lookup_t1 - hap_t0
        row["contig_lookup_seconds"] = round(lookup_seconds, 2)
        print(f"  {hap} [1/3] contig lookup (samtools view): {lookup_seconds:.2f}s", file=sys.stderr)

        if contigs is None:
            print(f"  {hap}: SKIP -- samtools view failed against the aln-to-hg38 BAM: {err}",
                  file=sys.stderr)
            gene_calls[hap] = {}
            timing_rows.append(row)
            continue
        if not contigs:
            print(f"  {hap}: SKIP -- 0 contigs overlap {args.region} in this haplotype's own "
                  f"hg38 alignment. Unexpected (every revio person should have HLA region "
                  f"covered) -- worth a closer look for this specific person, not assumed benign.",
                  file=sys.stderr)
            gene_calls[hap] = {}
            timing_rows.append(row)
            continue
        row["n_contigs"] = len(contigs)
        print(f"  {hap}: {len(contigs)} contig(s) overlap {args.region}: {contigs}", file=sys.stderr)

        # --- Stage 2: trim (samtools faidx -- files-to-ready-for-immuannot step Marc flagged) ---
        trimmed_fa = os.path.join(person_dir, f"{hap}.trimmed.fa")
        ok, err = trim_assembly(fa_path, contigs, trimmed_fa)
        trim_t2 = time.perf_counter()
        trim_seconds = trim_t2 - lookup_t1
        row["trim_seconds"] = round(trim_seconds, 2)
        print(f"  {hap} [2/3] trim (samtools faidx): {trim_seconds:.2f}s", file=sys.stderr)

        if not ok:
            print(f"  {hap}: SKIP -- samtools faidx trim failed: {err}", file=sys.stderr)
            gene_calls[hap] = {}
            timing_rows.append(row)
            continue
        trimmed_mb = os.path.getsize(trimmed_fa) / 1e6
        row["trimmed_mb"] = round(trimmed_mb, 1)
        print(f"  {hap}: trimmed FASTA is {trimmed_mb:.1f} MB, ready for immuannot.sh", file=sys.stderr)

        # --- Stage 3: immuannot.sh itself ---
        outprefix = os.path.join(person_dir, hap)
        gtf = run_immuannot(args.immuannot_dir, args.refdir, trimmed_fa, outprefix, args.threads)
        hap_t3 = time.perf_counter()
        immuannot_seconds = hap_t3 - trim_t2
        row["immuannot_seconds"] = round(immuannot_seconds, 2)
        row["hap_total_seconds"] = round(hap_t3 - hap_t0, 2)
        gene_calls[hap] = parse_gtf(gtf) if gtf else {}
        print(f"  {hap} [3/3] immuannot.sh: {immuannot_seconds/60:.1f} min "
              f"({len(gene_calls[hap])} genes called)", file=sys.stderr)
        print(f"  {hap} stage breakdown -- lookup {lookup_seconds:.1f}s / trim {trim_seconds:.1f}s / "
              f"immuannot {immuannot_seconds/60:.1f}min / hap total {row['hap_total_seconds']/60:.1f}min",
              file=sys.stderr)

        timing_rows.append(row)

    genes = sorted(set(gene_calls.get("hap1", {})) | set(gene_calls.get("hap2", {})))
    gene_rows = [{
        "person_id": pid, "gene": gene,
        "immuannot_1": gene_calls.get("hap1", {}).get(gene, "NA"),
        "immuannot_2": gene_calls.get("hap2", {}).get(gene, "NA"),
    } for gene in genes]

    person_total = time.perf_counter() - person_t0
    budget = args.time_budget_min * 60
    verdict = "UNDER budget" if person_total <= budget else "OVER budget"
    print(f"  --- Person {pid} total: {person_total/60:.1f} min "
          f"(budget {args.time_budget_min} min) -- {verdict}. {len(genes)} genes called. ---",
          file=sys.stderr)
    return gene_rows, timing_rows


def write_incremental(df_new, out_path, key_cols):
    """Appends/overwrites rows for the people just processed, keeping prior people's rows from
    earlier invocations intact -- lets the eventual 60-person run be done in batches/reruns
    without losing already-computed results."""
    if os.path.exists(out_path):
        old = pd.read_csv(out_path, sep="\t", dtype=str)
        keys = set(df_new[key_cols[0]].astype(str))
        old = old[~old[key_cols[0]].astype(str).isin(keys)]
        df_new = pd.concat([old, df_new], ignore_index=True)
    df_new.to_csv(out_path, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("person_ids", nargs="+", help="One or more research_id values (cohort.tsv).")
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    ap.add_argument("--immuannot-dir", default=os.path.expanduser("~/tools/Immuannot"))
    ap.add_argument("--refdir", default=os.path.expanduser("~/tools/Immuannot_refdata"))
    ap.add_argument("--outroot", default=os.path.expanduser("~/pipeline_outputs"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--region", default=DEFAULT_REGION,
                    help=f"GRCh38 region used to pick overlapping contigs (default {DEFAULT_REGION}).")
    ap.add_argument("--time-budget-min", type=float, default=30,
                    help="Per-person minutes to compare against when printing the verdict "
                         "(informational only -- does not abort a run).")
    args = ap.parse_args()

    lr_path = os.path.join(args.mount, LR_MANIFEST)
    if not os.path.exists(lr_path):
        die(f"manifest not found: {lr_path} -- is the gcsfuse mount up? "
            f"(ENVIRONMENT.md quirk #11/#14: remount and `ls`-verify first.)")
    if not os.path.exists(os.path.join(args.immuannot_dir, "scripts", "immuannot.sh")):
        die(f"immuannot.sh not found under {args.immuannot_dir} -- run setup_immuannot.sh first.")

    print(f"Reading LR manifest: {lr_path}", file=sys.stderr)
    lr = pd.read_csv(lr_path, sep="\t", dtype=str)
    if "research_id" not in lr.columns:
        die(f"LR manifest missing research_id column. Actual columns: {list(lr.columns)}")

    run_t0 = time.perf_counter()
    all_gene_rows, all_timing_rows = [], []
    for pid in args.person_ids:
        gene_rows, timing_rows = process_person(pid, lr, args)
        all_gene_rows.extend(gene_rows)
        all_timing_rows.extend(timing_rows)

    n_ok = len({r["person_id"] for r in all_gene_rows})
    n_skipped = len(args.person_ids) - n_ok
    run_total = time.perf_counter() - run_t0
    print(f"\n=== Run summary: {n_ok}/{len(args.person_ids)} people produced calls "
          f"({n_skipped} skipped), {run_total/60:.1f} min total "
          f"({run_total/60/max(n_ok,1):.1f} min/person average) ===", file=sys.stderr)

    if all_timing_rows:
        timing_path = os.path.join(args.outroot, "immuannot_timing.tsv")
        write_incremental(pd.DataFrame(all_timing_rows), timing_path, ["person_id"])
        print(f"Timing log: {timing_path}", file=sys.stderr)

    if not all_gene_rows:
        die("No people produced any calls -- check the SKIP reasons above.")

    out_path = os.path.join(args.outroot, "immuannot_calls.tsv")
    write_incremental(pd.DataFrame(all_gene_rows), out_path, ["person_id"])
    print(f"Calls: {out_path}", file=sys.stderr)
    print("Aggregate-only note doesn't apply here the way it does for comparison_log.csv -- these "
          "ARE per-person real allele calls (participant data). Keep on the VM; do not commit or "
          "download raw calls (same rule as SMOKE_TEST_PICKS.local.md / comparison_log.csv).",
          file=sys.stderr)


if __name__ == "__main__":
    main()
