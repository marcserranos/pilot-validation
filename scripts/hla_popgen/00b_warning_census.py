#!/usr/bin/env python3
"""Break down `template_warning` tokens per gene across the cohort.

Why this exists, as its own script: `00_recon_vm.py` measured `template_warning` present on
**95.3%** of transcript rows on the real production cohort (2026-09-03, 50-person sample). That
single number invalidates a filter this project already relies on elsewhere -- `context/DECISIONS.md`'s
confidence convention is "`template_distance == 0` AND no `template_warning`", which on real data
would reject roughly 95% of all calls rather than a small unreliable tail.

But "95% warn" is not actionable on its own, because `template_warning` is not one thing. It
describes whether the TEMPLATE's CDS could be cleanly reconstructed from the gene-level alignment
(`searchTemplate.py`'s `checkCDScompleteness()`) -- not whether the typing call is wrong. Some
tokens are expected and benign:
  - `partial_CDS`     -- fires whenever the trimmed contig truncates a gene's span
  - `no-start_codon` / `no-stop_codon` -- structurally NORMAL for the pseudogenes in the 42-gene
                         panel (HLA-H/J/K/L/...), which genuinely lack valid codons
  - `inframe_stop`    -- the one token that actually suggests a broken reconstruction
So the decision this script informs is: which tokens should disqualify a call, and for which genes.
`03_novel_alleles.py --disqualifying-warnings` consumes that decision.

Read-only, local disk only (no gcsfuse mount needed), cheap on a sample.

Usage:
    python3 scripts/hla_popgen/00b_warning_census.py --limit 200
"""

import argparse
import collections
import glob
import gzip
import os
import re
import sys

DEFAULT_OUTROOT = "~/pipeline_outputs"
CLASSICAL = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1",
             "HLA-DPA1", "HLA-DPB1"]
GENE_RE = re.compile(r'gene_name "([^"]+)"')
WARN_RE = re.compile(r'template_warning "([^"]*)"')
NONE = "<none>"


def census(outroot, limit):
    dirs = sorted(glob.glob(os.path.join(outroot, "*", "immuannot_output")))
    if not dirs:
        sys.exit(f"FATAL: no person directories under {outroot!r}.\n"
                 f"       Expected {outroot}/<person_id>/immuannot_output/hap1.gtf.gz\n"
                 f"       If this is the wrong path, pass --outroot.")
    dirs = dirs[:limit]
    tok = collections.Counter()
    per_gene = collections.defaultdict(collections.Counter)
    n_gtf = 0
    for d in dirs:
        for hap in ("hap1", "hap2"):
            path = os.path.join(d, hap + ".gtf.gz")
            if not os.path.exists(path):
                continue
            n_gtf += 1
            with gzip.open(path, "rt") as f:
                for line in f:
                    if line.startswith("#") or "\ttranscript\t" not in line:
                        continue
                    attrs = line.split("\t")[8]
                    gm = GENE_RE.search(attrs)
                    if not gm:
                        continue
                    gene = gm.group(1)
                    wm = WARN_RE.search(attrs)
                    toks = [t.strip() for t in wm.group(1).split(",") if t.strip()] if wm else []
                    if not toks:
                        toks = [NONE]
                    for t in toks:
                        tok[t] += 1
                        per_gene[gene][t] += 1
    return tok, per_gene, len(dirs), n_gtf


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT)
    ap.add_argument("--limit", type=int, default=200,
                    help="people to sample (default 200; the token mix stabilizes fast)")
    args = ap.parse_args()

    outroot = os.path.expanduser(args.outroot)
    tok, per_gene, n_people, n_gtf = census(outroot, args.limit)
    if not tok:
        sys.exit("FATAL: parsed no transcript rows. Are the GTFs readable?")

    total = sum(tok.values())
    print(f"Sampled {n_people} people / {n_gtf} haplotype GTFs / {total} transcript rows\n")
    print("== warning tokens, all genes ==")
    for t, c in tok.most_common():
        print(f"{c:9d}  {100.0 * c / total:5.1f}%  {t}")

    print("\n== per gene: how often is the call CLEAN (no warning at all)? ==")
    print(f"{'gene':16s} {'n':>6s} {'clean':>7s} {'clean%':>7s}  top warning tokens")
    rest = sorted(g for g in per_gene if g not in CLASSICAL)
    for gene in [g for g in CLASSICAL if g in per_gene] + rest:
        c = per_gene[gene]
        n = sum(c.values())
        clean = c.get(NONE, 0)
        top = ", ".join(f"{k}:{v}" for k, v in c.most_common(4) if k != NONE)
        print(f"{gene:16s} {n:6d} {clean:7d} {100.0 * clean / n:6.1f}%  {top}")

    inframe = tok.get("inframe_stop", 0)
    print(f"\nIf only `inframe_stop` disqualifies, {100.0 * (1 - inframe / total):.1f}% of calls "
          f"remain eligible.\nIf ANY warning disqualifies, only "
          f"{100.0 * tok.get(NONE, 0) / total:.1f}% remain -- which is why the blanket gate is "
          f"unusable.\nSet the policy via 03_novel_alleles.py --disqualifying-warnings.")


if __name__ == "__main__":
    main()
