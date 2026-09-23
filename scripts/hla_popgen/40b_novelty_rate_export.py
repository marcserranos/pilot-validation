#!/usr/bin/env python3
"""Panel (c) source data for 40_figure1_v5.py -- gene x ancestry novelty TOTALS, exported once,
suppressed once. MUST run on the VM (recomputes from Table 1 directly, never from an
already-censored intermediate).

Why this script exists (orchestrator review, 2026-09-23): panel c's earlier rendering read
33_figure1_v3_compose.py's `rate_by_gene_ancestry`, which sources 24_novelty_by_field's
COMMITTED, already-<20-suppressed `novelty_field_counts.tsv`. That table is split into many small
sub-cells (gene x ancestry x field_class x artifact_label), each independently censored, so a
cell's *combined* novelty count (summed across sub-cells, honestly hatched per the 2026-09-23
correctness fix) ends up hatched almost everywhere even when the combined count would itself
clear 20. The fix is not a display trick -- it is computing the combined numerator BEFORE any
censoring, on the VM, and applying the <20 rule exactly once, to the number that is actually
written out.

Definitions -- reused EXACTLY from 24_novelty_by_field.py, not re-derived:
    - `add_call_labels(t1)` gives `field_class` (known / f4_noncoding / f3_synonymous /
      f2_protein / f1_undetermined / uncalled) from the nomenclature-field depth at which
      Immuannot wrote 'new', and `artifact_label` (clean / homopolymer_indel / partial_cds /
      inframe_stop).
    - `build_people(...)` gives the unrelated set and both ancestry schemes (anc_pred, anc_strict
      at 24's own default STRICT_MIN=0.9 -- the SAME scheme panel c has always used; this script
      does not switch to the 0.98 threshold used for panels a/b, and says so in its own output).
    - Denominator (`n_total`) = called haplotype-gene copies (copy_index==1, field_class !=
      'uncalled') for that gene x ancestry, among unrelated people.
    - `n_novel_any_field` = clean calls with field_class in {f2_protein, f3_synonymous,
      f4_noncoding} (Cole's "the 800,000 one" / 33's panel c metric -- dominated by non-coding).
    - `n_novel_cds`       = clean calls with field_class in {f2_protein, f3_synonymous} (coding
      sequence differs from the catalogue, whether or not the translated protein is also new).
    - `n_novel_protein`   = clean calls with field_class == f2_protein (translated protein not
      in IPD-IMGT/HLA by nomenclature -- same definition panel d already counts).

Disclosure: the <20 rule is applied exactly once, on the way out, to n_total and each n_novel_*
independently. A Wilson 95% CI (vc.wilson_ci, computed from the TRUE counts before suppression) is
exported alongside a rate only when both the numerator and denominator clear 20; otherwise the
count is written as the literal string "<20" and the rate/CI columns are left blank.

Usage (VM):
    python3 scripts/hla_popgen/40b_novelty_rate_export.py \\
        --out reports/hla_popgen/40_figure1_v5/panel_c_novelty_totals.tsv
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)
import _viz_common as vc

ANC = vc.ANCESTRY_ORDER
CLASSICAL = ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"]
SUPPRESS_BELOW = 20
METRICS = [
    ("any_field", ("f2_protein", "f3_synonymous", "f4_noncoding")),
    ("cds", ("f2_protein", "f3_synonymous")),
    ("protein", ("f2_protein",)),
]


def _load(fn, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_THIS_DIR, fn))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def suppressed(n):
    n = int(n)
    return "<%d" % SUPPRESS_BELOW if n < SUPPRESS_BELOW else str(n)


def run(args):
    m24 = _load("24_novelty_by_field.py", "novelty_by_field")

    print("[40b] loading table1 ...", flush=True)
    t1 = m24.load_table1(args.table1, args.limit)
    t1 = m24.add_call_labels(t1)

    people, n_removed = m24.build_people(
        t1["person_id"].unique(), args.cohort_membership, args.relatedness_table, args.kin_min,
        args.strict_threshold, args.skip_relatedness)
    anc_col = "anc_strict" if args.ancestry_scheme == "strict" else "anc_pred"
    keep = people[people["unrelated"]]
    anc_of = dict(zip(keep["person_id"].astype(str), keep[anc_col]))

    c1 = t1[pd.to_numeric(t1["copy_index"], errors="coerce") == 1].copy()
    c1["anc"] = c1["person_id"].astype(str).map(anc_of)
    c1 = c1[c1["anc"].isin(ANC) & c1["gene"].isin(CLASSICAL)]
    print("[40b] %d haplotype-gene calls in scope (unrelated, strict-ancestry, classical genes)"
          % len(c1), flush=True)

    rows = []
    for (g, a), sub in c1.groupby(["gene", "anc"]):
        called = sub[sub["field_class"] != "uncalled"]
        n_total = len(called)
        clean = called[called["artifact_label"] == "clean"]
        row = {"gene": g, "ancestry": a, "ancestry_scheme": args.ancestry_scheme,
              "ancestry_threshold": args.strict_threshold if args.ancestry_scheme == "strict"
              else None, "n_total": suppressed(n_total),
              "n_total_censored": n_total < SUPPRESS_BELOW}
        for key, classes in METRICS:
            n = int(clean["field_class"].isin(classes).sum())
            row["n_novel_%s" % key] = suppressed(n)
            if n >= SUPPRESS_BELOW and n_total >= SUPPRESS_BELOW:
                p, lo, hi = vc.wilson_ci(n, n_total)
                row["rate_%s_pct" % key] = round(float(p) * 100, 3)
                row["ci_%s_lo" % key] = round(float(lo) * 100, 3)
                row["ci_%s_hi" % key] = round(float(hi) * 100, 3)
            else:
                row["rate_%s_pct" % key] = ""
                row["ci_%s_lo" % key] = ""
                row["ci_%s_hi" % key] = ""
        rows.append(row)

    out = pd.DataFrame(rows)
    # completeness check: every (gene, ancestry) pair should be present even if all-zero, so 40's
    # compose can tell "cell genuinely absent from the cohort" apart from "not yet exported".
    full_idx = pd.MultiIndex.from_product([CLASSICAL, ANC], names=["gene", "ancestry"])
    have = pd.MultiIndex.from_frame(out[["gene", "ancestry"]])
    missing = full_idx.difference(have)
    if len(missing):
        print("[40b] WARNING: no calls at all for %d (gene, ancestry) pairs: %s"
              % (len(missing), list(missing)), flush=True)

    n_protein_ge20 = int((out["n_novel_protein"].apply(lambda s: not str(s).startswith("<")))
                         .sum())
    print("[40b] %d/%d gene x ancestry cells have n_novel_protein >= 20 (protein-level heatmap "
          "would be %s)" % (n_protein_ge20, len(out),
                            "usable" if n_protein_ge20 >= 0.5 * len(out) else "still mostly <20"),
          flush=True)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print("[40b] wrote %s (%d rows)" % (args.out, len(out)), flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=os.path.expanduser("~/pipeline_outputs/hla_calls_rich.tsv"))
    ap.add_argument("--cohort-membership",
                    default=os.path.expanduser("~/pipeline_outputs/cohort_membership.tsv"))
    ap.add_argument("--relatedness-table", default=os.path.expanduser(
        "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/"
        "relatedness/samples_relatedness.tsv"))
    ap.add_argument("--skip-relatedness", action="store_true")
    ap.add_argument("--kin-min", type=float, default=0.0442)
    ap.add_argument("--ancestry-scheme", choices=["strict", "pred"], default="strict")
    ap.add_argument("--strict-threshold", type=float, default=0.9,
                    help="Same scheme/threshold panel c has always used (24's own default), NOT "
                         "the 0.98 used for panels a/b -- deliberately kept separate, see module "
                         "docstring.")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", default=os.path.expanduser(
        "~/s03/results/40/panel_c_novelty_totals.tsv"))
    run(ap.parse_args(argv))


if __name__ == "__main__":
    main()
