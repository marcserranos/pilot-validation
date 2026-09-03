#!/usr/bin/env python3
"""Break down `template_warning` tokens per gene across the cohort.

Why this exists, as its own script: `00_recon_vm.py` originally measured `template_warning`
present on **95.3%** of transcript rows on the real production cohort (2026-09-03, 50-person
sample). That number was itself a bug -- it counted the attribute's mere PRESENCE, and Immuannot
writes the literal string `template_warning "NA"` to mean *no warning* on 57.4% of transcript
rows, far more often than it omits the attribute entirely (4.6%). The 200-person census this
script performs (2026-09-03) measured the TRUE breakdown:

| value | share | meaning |
|---|---|---|
| `NA` (literal) | 57.4% | clean |
| attribute absent | 4.6% | clean |
| `partial_CDS` | 24.8% | real warning |
| `no-start_codon` | 6.8% | real warning (mostly pseudogene biology) |
| `no-stop_codon` | 6.3% | real warning (mostly pseudogene biology) |
| `inframe_stop` | 0% | never observed |

So 62% of calls are clean and the true warning rate is ~38%, not ~95%. This still matters:
`context/DECISIONS.md`'s confidence convention is "`template_distance == 0` AND no
`template_warning`", and a naive `bool(template_warning)`/`"template_warning" in attrs` check
would reject the wrong ~95% of calls instead of the real ~38%.

But even 38% is not actionable as a single number, because `template_warning` is not one thing.
It describes whether the TEMPLATE's CDS could be cleanly reconstructed from the gene-level
alignment (`searchTemplate.py`'s `checkCDScompleteness()`) -- not whether the typing call is
wrong. Some tokens are expected and benign:
  - `partial_CDS`     -- fires whenever the trimmed contig truncates a gene's span (disqualifying
                         by default: a truncated CDS can't support a novel-allele claim, and it's
                         only ~2.4% of classical-gene calls)
  - `no-start_codon` / `no-stop_codon` -- structurally NORMAL for the pseudogenes in the 42-gene
                         panel (HLA-H/J/K/L/...), which genuinely lack valid codons -- NOT
                         disqualifying by default
  - `inframe_stop`    -- the one token that would genuinely suggest a broken reconstruction if it
                         appeared; disqualifying by default despite never being observed
So the decision this script informs is: which tokens should disqualify a call, and for which genes.
`03_novel_alleles.py --disqualifying-warnings` consumes that decision (default:
`partial_CDS,inframe_stop`).

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
                    # Immuannot writes the literal string template_warning "NA" to mean *no
                    # warning* on 57.4% of transcript rows -- far more often than it omits the
                    # attribute (4.6%). Both spellings of "clean" are folded into the same NONE
                    # bucket here; treating "NA" as a real token (the original bug) makes clean
                    # calls look like they carry a warning token called "NA".
                    raw = wm.group(1).strip() if wm else ""
                    if raw.upper() in ("", "NA"):
                        toks = [NONE]
                    else:
                        toks = [t.strip() for t in raw.split(",") if t.strip()]
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

    clean_pct = 100.0 * tok.get(NONE, 0) / total
    inframe = tok.get("inframe_stop", 0)
    print(f"\nTrue clean rate ({NONE}, which folds in Immuannot's literal \"NA\" token): "
          f"{clean_pct:.1f}% -- NOT the bare-presence figure a naive `bool(template_warning)` or "
          f"`\"template_warning\" in attrs` check would report.")
    print(f"If only `inframe_stop` disqualifies, {100.0 * (1 - inframe / total):.1f}% of calls "
          f"remain eligible (a near-no-op: inframe_stop is essentially never observed).\n"
          f"If ANY real warning disqualifies, {clean_pct:.1f}% remain -- still a substantial "
          f"cut, which is why disqualification must be token-aware, not a blanket gate.\n"
          f"Set the policy via 03_novel_alleles.py --disqualifying-warnings "
          f"(default: partial_CDS,inframe_stop).")


if __name__ == "__main__":
    main()
