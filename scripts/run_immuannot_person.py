#!/usr/bin/env python3
"""Run Immuannot (github.com/YingZhou001/Immuannot) on one or more AoU people, once Aleix's
10-person list arrives. Companion to setup_immuannot.sh (run that first, once).

Immuannot types HLA alleles from a PHASED ASSEMBLY contig (one haplotype's collapsed consensus
FASTA), not reads/BAM -- see DECISIONS.md "Assembly-based HLA typing". Only people on the
revio/sequel2e/sequel2 platforms have this data (reports/lr_data_census/README.md); ont-r10.4.1/
ont-r9.4.1 do not and will be skipped with a clear reason, not silently dropped.

Per person: resolves assembly_hap1_fa / assembly_hap2_fa mount-relative paths from the v9 lrWGS
manifest (existence-checked directly against the mount -- same discipline as
build_experiment_d_cohort.py's find_existing_bam(), never a path-pattern guess), then runs
immuannot.sh once per haplotype (the tool's own documented mode -- it does not accept both
haplotypes in one run), verifies each expected *.gtf.gz output actually exists (never trust exit
code alone -- same lesson as ENVIRONMENT.md quirks #17/#18 for SpecImmune/SpecHLA), and writes a
per-person, per-gene 2-allele call table.

GTF-PARSING CAVEAT (read before trusting the output table): Immuannot's own docs say to use the
"consensus" or "allele" GTF attribute, not "template_allele" -- but the exact attribute key
spelling has NOT been verified against real output (no local test data available while writing
this). parse_gtf() prints the first record's raw attribute string per haplotype run so this can be
eyeballed and the regex adjusted on the very first real run, rather than trusting it blind.

Usage (from ~/repos/pilot-validation, inside the specimmune pixi env -- reuses its minimap2):
  pixi run -e specimmune -- python3 scripts/run_immuannot_person.py <person_id> [<person_id> ...] \\
      [--mount ~/mnt/aou-controlled] [--immuannot-dir ~/tools/Immuannot] \\
      [--refdir ~/tools/Immuannot_refdata] [--outroot ~/pipeline_outputs] [--threads 4]
"""
import argparse
import gzip
import os
import re
import subprocess
import sys

import pandas as pd

BUCKET_PREFIX = "gs://vwb-aou-datasets-controlled/"
LR_MANIFEST = "v9/wgs/long_read/manifest.tsv"
HAP_COLS = {"hap1": "assembly_hap1_fa", "hap2": "assembly_hap2_fa"}


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


def resolve_haplotype_paths(lr, person_id, mount):
    """Returns {'hap1': rel_path_or_None, 'hap2': rel_path_or_None} -- existence-checked directly
    against the mount, never assumed from the manifest cell being non-null."""
    rows = lr[lr["research_id"] == str(person_id)]
    if rows.empty:
        return None
    out = {}
    for hap, col in HAP_COLS.items():
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


def run_immuannot(immuannot_dir, refdir, contig_path, outprefix, threads):
    script = os.path.join(immuannot_dir, "scripts", "immuannot.sh")
    cmd = ["bash", script, "-c", contig_path, "-r", refdir, "-o", outprefix, "-t", str(threads)]
    print(f"  Running: {' '.join(cmd)}", file=sys.stderr)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    expected = outprefix + ".gtf.gz"
    if not os.path.exists(expected):
        # Never trust exit code alone (ENVIRONMENT.md quirks #17/#18) -- surface the real logs.
        print(f"  WARNING: expected output {expected} missing after run "
              f"(exit code {proc.returncode}). stdout/stderr follow:", file=sys.stderr)
        print(proc.stdout[-2000:], file=sys.stderr)
        print(proc.stderr[-2000:], file=sys.stderr)
        return None
    return expected


def parse_gtf(gtf_gz_path, show_raw_sample=True):
    """Extracts (gene, allele) pairs from Immuannot's gtf.gz. UNVERIFIED against real output --
    see module docstring caveat. Prints the first record's raw attributes for manual sanity check."""
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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("person_ids", nargs="+", help="One or more research_id values from Aleix's list.")
    ap.add_argument("--mount", default=os.path.expanduser("~/mnt/aou-controlled"))
    ap.add_argument("--immuannot-dir", default=os.path.expanduser("~/tools/Immuannot"))
    ap.add_argument("--refdir", default=os.path.expanduser("~/tools/Immuannot_refdata"))
    ap.add_argument("--outroot", default=os.path.expanduser("~/pipeline_outputs"))
    ap.add_argument("--threads", type=int, default=4)
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

    all_rows = []
    for pid in args.person_ids:
        print(f"\n=== Person {pid} ===", file=sys.stderr)
        haps = resolve_haplotype_paths(lr, pid, args.mount)
        if haps is None:
            print(f"  SKIP: no manifest row for {pid} at all.", file=sys.stderr)
            continue
        missing = [h for h, p in haps.items() if p is None]
        if missing:
            print(f"  SKIP: missing {missing} -- this person's platform likely has no assembly "
                  f"data (only revio/sequel2e/sequel2 do -- reports/lr_data_census/README.md). "
                  f"Flag to Aleix if this person was on his list.", file=sys.stderr)
            continue

        person_dir = os.path.join(args.outroot, str(pid), "immuannot_output")
        os.makedirs(person_dir, exist_ok=True)
        gene_calls = {}
        for hap in ("hap1", "hap2"):
            contig_path = os.path.join(args.mount, haps[hap])
            outprefix = os.path.join(person_dir, hap)
            gtf = run_immuannot(args.immuannot_dir, args.refdir, contig_path, outprefix, args.threads)
            if gtf is None:
                gene_calls[hap] = {}
                continue
            gene_calls[hap] = parse_gtf(gtf)

        genes = sorted(set(gene_calls["hap1"]) | set(gene_calls["hap2"]))
        for gene in genes:
            all_rows.append({
                "person_id": pid,
                "gene": gene,
                "immuannot_1": gene_calls["hap1"].get(gene, "NA"),
                "immuannot_2": gene_calls["hap2"].get(gene, "NA"),
            })
        print(f"  {len(genes)} genes called.", file=sys.stderr)

    if not all_rows:
        die("No people produced any calls -- check the SKIP reasons above.")

    out_df = pd.DataFrame(all_rows)
    out_path = os.path.join(args.outroot, "immuannot_calls.tsv")
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"\nWrote {len(out_df)} gene-rows to {out_path}", file=sys.stderr)
    print("Aggregate-only note doesn't apply here the way it does for comparison_log.csv -- these "
          "ARE per-person real allele calls (participant data). Keep on the VM; do not commit or "
          "download raw calls (same rule as SMOKE_TEST_PICKS.local.md / comparison_log.csv).",
          file=sys.stderr)


if __name__ == "__main__":
    main()
