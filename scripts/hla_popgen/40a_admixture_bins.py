#!/usr/bin/env python3
"""Panel (a) source data for 40_figure1_v5.py -- binned, disclosure-safe admixture strip.

MUST run on the VM: loads per-person admixture proportions (participant-level, never allowed to
leave the VM) and aggregates them into consecutive bins of >=20 people each before writing
anything out. Only bin-level MEANS are exported -- no bin, however narrow, can be traced back to
an individual's admixture vector.

Method (orchestrator spec, 2026-09-23 review of figure1_v5):
  1. Take every unrelated person (reuses 24_novelty_by_field.build_people -- never re-derived).
  2. Sort by predicted ancestry (anc_pred), then by that ancestry's own dominant-component
     probability, descending (least admixed -> most admixed within each block) -- same ordering
     36_figure1_native.py's panel (a) used, just binned instead of drawn per person.
  3. Walk the sorted list in consecutive bins of exactly `--bin-size` people (default 20, i.e. the
     AoU small-cell floor itself -- a bin can never be smaller than the disclosure threshold). The
     last, short bin within an ancestry block is merged into the previous bin rather than emitted
     under-sized.
  4. Per bin: mean of each p_<ANC> column, the bin's dominant (predicted) ancestry, and n_people
     (always >=20). ~500 bins over ~11,800 people at bin size 20-25.

Output: panel_a_admixture_bins.tsv (bin_id, anc, n_people, p_afr..p_sas, order_within_anc) --
this is the ONLY thing that leaves the VM for panel (a).

Usage (VM):
    python3 scripts/hla_popgen/40a_admixture_bins.py --out reports/hla_popgen/40_figure1_v5/panel_a_admixture_bins.tsv
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

ANC = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
MIN_BIN = 20


def _load(fn, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_THIS_DIR, fn))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def bin_people(pa, bin_size):
    """pa: one row per person, columns anc (predicted ancestry, sorted block), dom (dominant
    ancestry's own probability), p_afr..p_sas -- already sorted by (anc, dom desc), matching
    36_figure1_native.py's panel_a_data() ordering. Returns one row per bin."""
    cols = ["p_" + a.lower() for a in ANC]
    rows = []
    for anc, grp in pa.groupby("anc", sort=False):
        grp = grp.reset_index(drop=True)
        n = len(grp)
        starts = list(range(0, n, bin_size))
        # merge a short trailing bin into the previous one rather than emit it under-sized
        if len(starts) > 1 and (n - starts[-1]) < MIN_BIN:
            starts = starts[:-1]
        for bi, s in enumerate(starts):
            e = starts[bi + 1] if bi + 1 < len(starts) else n
            sub = grp.iloc[s:e]
            if len(sub) < MIN_BIN:
                continue
            row = {"anc": anc, "order_within_anc": bi, "n_people": len(sub)}
            for c in cols:
                row[c] = round(float(sub[c].mean()), 6)
            rows.append(row)
    return pd.DataFrame(rows)


def run(args):
    m36 = _load("36_figure1_native.py", "fig1_native_v36")
    m24 = _load("24_novelty_by_field.py", "novelty_by_field")

    t1 = pd.read_csv(args.table1, sep="\t", dtype=str, low_memory=False)
    cohort = pd.read_csv(args.cohort_membership, sep="\t", dtype=str)
    people, n_removed = m24.build_people(
        sorted(set(t1["person_id"].astype(str))), args.cohort_membership,
        args.relatedness_table, args.kin_min, args.strict_threshold, args.skip_relatedness)
    keep = people[people["unrelated"]]
    all_unrelated = set(keep["person_id"].astype(str))

    pa, _cols = m36.panel_a_data(cohort, all_unrelated, args.strict_threshold)
    print("[40a] %d unrelated people, predicted ancestry, before binning" % len(pa), flush=True)

    bins = bin_people(pa, args.bin_size)
    print("[40a] %d bins (bin size >= %d, %d)" % (len(bins), MIN_BIN, args.bin_size), flush=True)
    assert (bins["n_people"] >= MIN_BIN).all(), "a bin under the disclosure floor slipped through"

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    bins.to_csv(args.out, sep="\t", index=False)
    print("[40a] wrote %s" % args.out, flush=True)


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
    ap.add_argument("--strict-threshold", type=float, default=0.98)
    ap.add_argument("--bin-size", type=int, default=20)
    ap.add_argument("--out", default=os.path.expanduser(
        "~/s03/results/40/panel_a_admixture_bins.tsv"))
    run(ap.parse_args(argv))


if __name__ == "__main__":
    main()
