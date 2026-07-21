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
"chr6:29,500,000-33,500,000" out of them the way slice_and_fastq.sh does for aligned BAMs.

Two-tier design, preferring the tighter one:
  1. PREFERRED -- true sub-range extraction via each haplotype's own assembly-to-hg38 .paf file
     (assembly_hap{1,2}_aln2_hg38_paf). A .paf line's qstart/qend/tstart/tend are the REAL aligned
     boundaries of that block, not an estimate, so regions_from_paf() takes the union of
     (min qstart, max qend) across every block overlapping the target region, per contig, padded
     by --pad (default 100,000bp -- reuses the exact number EXPERIMENTS.md's SpecImmune-LR sweep
     already proved safe: zero degradation, 3.9x smaller than full-pad) and clamped to the contig's
     own length. Confirmed 2026-07-21 on person 1001871: cut the trimmed input further than
     whole-contig extraction alone -- see immuannot_timing.tsv's whole_contig_mb/padded_mb columns
     for the actual per-person reduction once run.
  2. FALLBACK -- whole-contig extraction via the assembly-to-hg38 .bam (assembly_hap{1,2}_aln2_hg38_bam),
     used only if this person's .paf is missing or has nothing recognizable overlapping the region.
     No CIGAR/indel math, no risk of clipping a gene -- just less tight than tier 1.
Both confirmed present alongside the raw assembly FASTA for revio/sequel2e per the census.

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
import concurrent.futures
import glob
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
PAF_COLS = {"hap1": "assembly_hap1_aln2_hg38_paf", "hap2": "assembly_hap2_aln2_hg38_paf"}
DEFAULT_REGION = "chr6:29500000-33500000"  # the project's standing HLA window (ENVIRONMENT.md)
# Reuses the exact number ENVIRONMENT.md/EXPERIMENTS.md already proved safe for SpecImmune-LR
# (pad100k: zero degradation across 8 genes, 3.9x smaller than full-pad) -- Marc, 2026-07-21:
# "a really safe reduction" for the analogous concern here (not clipping a gene near a coordinate-
# mapping boundary), not a new untested number.
DEFAULT_PAD = 100_000


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


def parse_region(region_str):
    chrom, coords = region_str.split(":")
    start_str, end_str = coords.split("-")
    return chrom, int(start_str), int(end_str)


def regions_from_paf(paf_path, chrom, region_start, region_end, pad):
    """True sub-range trim (2026-07-21, Marc: explored the padding-reduction margin properly
    rather than dismissing it -- see EXPERIMENTS.md's SpecImmune-LR pad100k sweep for the
    precedent this reuses). A .paf line's own qstart/qend/tstart/tend are the REAL aligned
    boundaries of that block (not an extrapolated estimate), so padding here only has to cover
    minor internal-indel drift and the chance a gene sits just past a block's reported edge --
    not a coordinate system we're guessing at. For each contig with >=1 PAF block overlapping the
    target region, takes the union of (min qstart, max qend) across its overlapping blocks, pads
    by `pad` on each side, and clamps to [0, qlen].

    Returns (regions, qlens, seen_tnames):
      regions: {contig: (start0, end0)} -- padded/clamped sub-range, ready for a faidx region string
      qlens: {contig: qlen} -- full contig length, from the PAF itself, no extra samtools call
      seen_tnames: set of every target name actually seen (for diagnosing a naming mismatch,
                   e.g. "chr6" vs "6", rather than silently returning empty and looking like "no
                   overlap" when it's really "wrong name assumed").
    """
    opener = gzip.open if paf_path.endswith(".gz") else open
    accept_chrom = {chrom, chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"}
    regions, qlens, seen_tnames = {}, {}, set()
    with opener(paf_path, "rt") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qname, qlen, qstart, qend = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
            tname, tstart, tend = fields[5], int(fields[7]), int(fields[8])
            seen_tnames.add(tname)
            if tname not in accept_chrom:
                continue
            if tend <= region_start or tstart >= region_end:
                continue
            qlens[qname] = qlen
            new_start, new_end = max(0, qstart - pad), min(qlen, qend + pad)
            if qname in regions:
                s, e = regions[qname]
                regions[qname] = (min(s, new_start), max(e, new_end))
            else:
                regions[qname] = (new_start, new_end)
    return regions, qlens, seen_tnames


def trim_assembly(fa_path, targets, out_fa_path):
    """samtools faidx extracts each entry in `targets` -- either a bare contig name (whole contig,
    the fallback path) or a "contig:start-end" region string (the preferred sub-range path via
    regions_from_paf) -- samtools faidx accepts both forms interchangeably, so this function
    doesn't need to know which kind it got. fa_path lives on the read-only gcsfuse mount, but
    samtools faidx's default index locations are alongside fa_path itself -- confirmed 2026-07-21
    (person 1001871) this fails with "Permission denied" building the .gzi/.fai next to the
    read-only source. Fix: --fai-idx/--gzi-idx redirect the index files into out_fa_path's own
    (writable) directory instead -- samtools still just READS fa_path off the mount, only the
    small index files themselves move. Confirmed present in this VM's samtools via
    `samtools faidx --help` before relying on it, not assumed."""
    fai_idx = out_fa_path + ".src.fai"
    gzi_idx = out_fa_path + ".src.gzi"
    with open(out_fa_path, "w") as out:
        proc = subprocess.run(
            ["samtools", "faidx", "--fai-idx", fai_idx, "--gzi-idx", gzi_idx, fa_path] + targets,
            stdout=out, stderr=subprocess.PIPE, text=True)
    ok = proc.returncode == 0 and os.path.getsize(out_fa_path) > 0
    return ok, proc.stderr


def find_immuannot_script(immuannot_dir):
    """Never hardcode the in-repo path -- confirmed 2026-07-21 the upstream repo ships immuannot.sh
    under a VERSIONED folder (scripts.pub.v3, not plain scripts/ as its own README shows), and that
    name could drift again on a future re-clone. Same lesson as ENVIRONMENT.md quirk #13: find the
    real file, don't trust a documented path pattern."""
    matches = glob.glob(os.path.join(immuannot_dir, "**", "immuannot.sh"), recursive=True)
    if not matches:
        return None
    matches.sort(key=lambda p: p.count(os.sep))  # prefer the shallowest match
    return matches[0]


def run_immuannot(script, refdir, contig_path, outprefix, threads):
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


def process_haplotype(pid, hap, fa_rel, paf_rel, aln_rel, person_dir, immuannot_script,
                       args, resolve_seconds, threads):
    """One haplotype's full trim + immuannot.sh + parse, fully self-contained (returns its own
    row/calls rather than mutating shared state) so it's safe to run hap1 and hap2 concurrently
    in separate threads (Marc, 2026-07-21: why sum the two haps instead of running them at once?)
    without any locking. `threads` is passed in explicitly (not read from args) because it differs
    between sequential mode (args.threads, e.g. 4) and --parallel-haps mode (args.threads split
    across the two concurrent haps, e.g. 2 each)."""
    row = {"person_id": pid, "hap": hap, "resolve_seconds": round(resolve_seconds, 2),
           "trim_method": None, "n_contigs": None,
           "whole_contig_mb": None, "padded_mb": None, "trimmed_mb": None,
           "contig_lookup_seconds": None, "trim_seconds": None, "immuannot_seconds": None,
           "hap_total_seconds": None}
    hap_t0 = time.perf_counter()
    fa_path = os.path.join(args.mount, fa_rel)

    # --- Stage 1: find trim targets. Prefer exact sub-ranges from the .paf (Marc, 2026-07-21:
    # explored the padding-reduction margin properly -- reuses the proven-safe pad100k from
    # the SpecImmune-LR sweep). Falls back to whole-contig via the .bam if this person's .paf
    # is missing or has nothing recognizable overlapping the region, so nobody gets skipped
    # just because the (optional) .paf isn't there. ---
    targets, trim_method = None, None
    if paf_rel:
        paf_path = os.path.join(args.mount, paf_rel)
        chrom, rstart, rend = parse_region(args.region)
        regions, qlens, seen_tnames = regions_from_paf(paf_path, chrom, rstart, rend, args.pad)
        if regions:
            targets = [f"{c}:{s + 1}-{e}" for c, (s, e) in regions.items()]
            trim_method = "paf_region"
            whole_bp = sum(qlens[c] for c in regions)
            padded_bp = sum(e - s for s, e in regions.values())
            row["n_contigs"] = len(regions)
            row["whole_contig_mb"] = round(whole_bp / 1e6, 2)
            row["padded_mb"] = round(padded_bp / 1e6, 2)
            pct = 100 * (1 - padded_bp / max(whole_bp, 1))
            print(f"  {hap}: {len(regions)} contig(s) overlap {args.region} via .paf "
                  f"(pad={args.pad:,}bp) -- whole-contig would be {whole_bp/1e6:.1f} MB, "
                  f"padded sub-range is {padded_bp/1e6:.1f} MB ({pct:.0f}% smaller)",
                  file=sys.stderr)
        else:
            print(f"  {hap}: .paf present but 0 contigs overlap {args.region} -- target names "
                  f"seen in this .paf: {sorted(seen_tnames)[:10]} -- falling back to "
                  f"whole-contig via .bam.", file=sys.stderr)
    else:
        print(f"  {hap}: no .paf resolved for this person -- using whole-contig fallback.",
              file=sys.stderr)

    if targets is None:
        aln_path = os.path.join(args.mount, aln_rel)
        contigs, err = contigs_overlapping_region(aln_path, args.region)
        if contigs is None:
            print(f"  {hap}: SKIP -- samtools view failed against the aln-to-hg38 BAM: {err}",
                  file=sys.stderr)
            return row, {}
        if not contigs:
            print(f"  {hap}: SKIP -- 0 contigs overlap {args.region} (checked both .paf and "
                  f".bam). Unexpected -- worth a closer look for this specific person, not "
                  f"assumed benign.", file=sys.stderr)
            return row, {}
        targets, trim_method = contigs, "bam_whole_contig"
        row["n_contigs"] = len(contigs)
        print(f"  {hap}: {len(contigs)} contig(s) overlap {args.region} via .bam "
              f"(whole-contig): {contigs}", file=sys.stderr)

    row["trim_method"] = trim_method
    lookup_t1 = time.perf_counter()
    lookup_seconds = lookup_t1 - hap_t0
    row["contig_lookup_seconds"] = round(lookup_seconds, 2)
    print(f"  {hap} [1/3] find trim targets ({trim_method}): {lookup_seconds:.2f}s", file=sys.stderr)

    # --- Stage 2: trim (samtools faidx -- files-to-ready-for-immuannot step Marc flagged) ---
    trimmed_fa = os.path.join(person_dir, f"{hap}.trimmed.fa")
    ok, err = trim_assembly(fa_path, targets, trimmed_fa)
    trim_t2 = time.perf_counter()
    trim_seconds = trim_t2 - lookup_t1
    row["trim_seconds"] = round(trim_seconds, 2)
    print(f"  {hap} [2/3] trim (samtools faidx): {trim_seconds:.2f}s", file=sys.stderr)

    if not ok:
        print(f"  {hap}: SKIP -- samtools faidx trim failed: {err}", file=sys.stderr)
        return row, {}
    trimmed_mb = os.path.getsize(trimmed_fa) / 1e6
    row["trimmed_mb"] = round(trimmed_mb, 1)
    print(f"  {hap}: trimmed FASTA is {trimmed_mb:.1f} MB, ready for immuannot.sh", file=sys.stderr)

    # --- Stage 3: immuannot.sh itself ---
    outprefix = os.path.join(person_dir, hap)
    gtf = run_immuannot(immuannot_script, args.refdir, trimmed_fa, outprefix, threads)
    hap_t3 = time.perf_counter()
    immuannot_seconds = hap_t3 - trim_t2
    row["immuannot_seconds"] = round(immuannot_seconds, 2)
    row["hap_total_seconds"] = round(hap_t3 - hap_t0, 2)
    calls = parse_gtf(gtf) if gtf else {}
    print(f"  {hap} [3/3] immuannot.sh (threads={threads}): {immuannot_seconds/60:.1f} min "
          f"({len(calls)} genes called)", file=sys.stderr)
    print(f"  {hap} stage breakdown -- lookup {lookup_seconds:.1f}s / trim {trim_seconds:.1f}s / "
          f"immuannot {immuannot_seconds/60:.1f}min / hap total {row['hap_total_seconds']/60:.1f}min",
          file=sys.stderr)
    return row, calls


def process_person(pid, lr, args, immuannot_script):
    """Returns (gene_rows, timing_rows) for one person. Every stage is timestamped separately
    (manifest resolution, contig lookup, trim, immuannot.sh) and printed as it completes -- Marc's
    2026-07-21 ask, to see exactly which stage the "BAM things" (contig lookup / trim) actually
    cost versus immuannot.sh itself, rather than one bundled number. A timing row is recorded for
    every haplotype attempt, including ones that get skipped partway through, so a slow-then-failing
    stage is still visible in immuannot_timing.tsv, not silently dropped.

    hap1/hap2 run SEQUENTIALLY by default (each gets the full --threads) -- total person time is
    the SUM of both haps' immuannot.sh runtime, not the slower one alone. Pass --parallel-haps to
    run them concurrently instead (each then gets --threads/2, minimum 1, so the two together don't
    oversubscribe the VM's cores) -- Marc, 2026-07-21: tested on person 1001871 before deciding
    whether to use this for the full 60-person run, same "test then decide" discipline as the pad
    sweep. Concurrent stderr from both haps interleaves in the terminal -- expected, not a bug."""
    print(f"\n=== Person {pid} ===", file=sys.stderr)
    person_t0 = time.perf_counter()

    resolve_t0 = time.perf_counter()
    fa_paths = resolve_cols(lr, pid, args.mount, FA_COLS)
    aln_paths = resolve_cols(lr, pid, args.mount, ALN_COLS)
    resolve_seconds = time.perf_counter() - resolve_t0
    print(f"  manifest path resolution: {resolve_seconds:.2f}s", file=sys.stderr)

    paf_paths = resolve_cols(lr, pid, args.mount, PAF_COLS)  # optional -- not gated on below

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

    hap_args = {
        hap: (pid, hap, fa_paths[hap], paf_paths.get(hap), aln_paths[hap], person_dir,
              immuannot_script, args, resolve_seconds)
        for hap in ("hap1", "hap2")
    }
    results = {}
    if args.parallel_haps:
        per_hap_threads = max(1, args.threads // 2)
        print(f"  --parallel-haps: running hap1+hap2 concurrently, {per_hap_threads} threads "
              f"each (of --threads {args.threads} total)", file=sys.stderr)
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            futures = {
                pool.submit(process_haplotype, *hap_args[hap], per_hap_threads): hap
                for hap in ("hap1", "hap2")
            }
            for fut in concurrent.futures.as_completed(futures):
                results[futures[fut]] = fut.result()
    else:
        for hap in ("hap1", "hap2"):
            results[hap] = process_haplotype(*hap_args[hap], args.threads)

    timing_rows = [results[hap][0] for hap in ("hap1", "hap2")]
    gene_calls = {hap: results[hap][1] for hap in ("hap1", "hap2")}

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
    ap.add_argument("--parallel-haps", action="store_true",
                    help="Run hap1 and hap2 concurrently instead of sequentially -- total "
                         "person time trends toward the slower haplotype instead of the sum of "
                         "both. Each hap gets --threads/2 (min 1) to avoid oversubscribing the "
                         "VM's cores. Untested at scale as of 2026-07-21 -- verify on 1-2 people "
                         "before using for a full 60-person run.")
    ap.add_argument("--region", default=DEFAULT_REGION,
                    help=f"GRCh38 region used to pick overlapping contigs (default {DEFAULT_REGION}).")
    ap.add_argument("--pad", type=int, default=DEFAULT_PAD,
                    help=f"bp padding on each side of a .paf-mapped block before extraction "
                         f"(default {DEFAULT_PAD:,} -- reuses the SpecImmune-LR pad100k value "
                         f"already proven safe in EXPERIMENTS.md). Only affects the .paf-based "
                         f"sub-range path; the .bam whole-contig fallback ignores this.")
    ap.add_argument("--time-budget-min", type=float, default=30,
                    help="Per-person minutes to compare against when printing the verdict "
                         "(informational only -- does not abort a run).")
    ap.add_argument("--force", action="store_true",
                    help="Re-process a person even if already present in immuannot_calls.tsv "
                         "(default: skip -- makes a 60-person run safe to kill and re-launch "
                         "with the exact same command).")
    args = ap.parse_args()

    lr_path = os.path.join(args.mount, LR_MANIFEST)
    if not os.path.exists(lr_path):
        die(f"manifest not found: {lr_path} -- is the gcsfuse mount up? "
            f"(ENVIRONMENT.md quirk #11/#14: remount and `ls`-verify first.)")
    immuannot_script = find_immuannot_script(args.immuannot_dir)
    if immuannot_script is None:
        die(f"immuannot.sh not found anywhere under {args.immuannot_dir} -- run setup_immuannot.sh first.")
    print(f"Using immuannot.sh at: {immuannot_script}", file=sys.stderr)

    print(f"Reading LR manifest: {lr_path}", file=sys.stderr)
    lr = pd.read_csv(lr_path, sep="\t", dtype=str)
    if "research_id" not in lr.columns:
        die(f"LR manifest missing research_id column. Actual columns: {list(lr.columns)}")

    calls_path = os.path.join(args.outroot, "immuannot_calls.tsv")
    already_done = set()
    if os.path.exists(calls_path) and not args.force:
        already_done = set(pd.read_csv(calls_path, sep="\t", dtype=str)["person_id"].astype(str))
        print(f"{len(already_done)} people already have calls in {calls_path} -- will be skipped "
              f"(pass --force to redo them). This is what makes killing and re-launching a 60-"
              f"person run with the exact same command safe.", file=sys.stderr)

    run_t0 = time.perf_counter()
    all_gene_rows, all_timing_rows = [], []
    n_skipped_done = 0
    for pid in args.person_ids:
        if str(pid) in already_done:
            print(f"\n=== Person {pid}: SKIP -- already in {calls_path} (--force to redo) ===",
                  file=sys.stderr)
            n_skipped_done += 1
            continue
        gene_rows, timing_rows = process_person(pid, lr, args, immuannot_script)
        all_gene_rows.extend(gene_rows)
        all_timing_rows.extend(timing_rows)

    n_ok = len({r["person_id"] for r in all_gene_rows})
    n_failed = len(args.person_ids) - n_ok - n_skipped_done
    run_total = time.perf_counter() - run_t0
    print(f"\n=== Run summary: {n_ok} produced calls this invocation, {n_skipped_done} already "
          f"done (skipped), {n_failed} failed/skipped this time, {run_total/60:.1f} min elapsed "
          f"({run_total/60/max(n_ok,1):.1f} min/person average, this invocation only) ===",
          file=sys.stderr)

    if all_timing_rows:
        timing_path = os.path.join(args.outroot, "immuannot_timing.tsv")
        write_incremental(pd.DataFrame(all_timing_rows), timing_path, ["person_id"])
        print(f"Timing log: {timing_path}", file=sys.stderr)

    if not all_gene_rows and n_skipped_done == 0:
        die("No people produced any calls -- check the SKIP reasons above.")
    if not all_gene_rows:
        print("Nothing new to write this invocation (everyone requested was already done).",
              file=sys.stderr)
        return

    out_path = os.path.join(args.outroot, "immuannot_calls.tsv")
    write_incremental(pd.DataFrame(all_gene_rows), out_path, ["person_id"])
    print(f"Calls: {out_path}", file=sys.stderr)
    print("Aggregate-only note doesn't apply here the way it does for comparison_log.csv -- these "
          "ARE per-person real allele calls (participant data). Keep on the VM; do not commit or "
          "download raw calls (same rule as SMOKE_TEST_PICKS.local.md / comparison_log.csv).",
          file=sys.stderr)


if __name__ == "__main__":
    main()
