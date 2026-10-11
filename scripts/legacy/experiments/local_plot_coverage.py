#!/usr/bin/env python3
"""LOCAL-ONLY: static coverage chart (catalogued vs observed alleles per gene, 2-field).

Data below is the printed output of experiment_a3_allele_space_coverage.py (2026-07-19 run,
535,658-person AoU-native cohort) -- reproduced here as a small static table since the coverage
summary itself isn't written to a CSV by that script (only the per-group tree data is). Matches
the table in reports/aou_callset_validation.md section 3 exactly.

Usage:
  python3 local_plot_coverage.py --out reports/figures/coverage_catalogued_vs_observed.png
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np

DATA = {
    # gene: (catalogued_2f, observed_2f)
    "A": (5873, 305), "B": (7069, 578), "C": (5625, 482), "DRB1": (2756, 382),
    "DQA1": (562, 82), "DQB1": (1975, 172), "DPA1": (474, 89), "DPB1": (1974, 283),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="coverage_catalogued_vs_observed.png")
    args = ap.parse_args()

    genes = list(DATA.keys())
    catalogued = [DATA[g][0] for g in genes]
    observed = [DATA[g][1] for g in genes]

    x = np.arange(len(genes))
    width = 0.38

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x - width / 2, catalogued, width, label="Catalogued (IMGT)", color="#888780")
    ax.bar(x + width / 2, observed, width, label="Observed (AoU, n=535,658)", color="#1D9E75")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(genes)
    ax.set_ylabel("Distinct 2-field alleles (log scale)")
    ax.set_title("AoU-native calls vs. the full IPD-IMGT/HLA catalogue")
    ax.legend(frameon=False)
    for i, (c, o) in enumerate(zip(catalogued, observed)):
        ax.text(i - width / 2, c * 1.15, str(c), ha="center", fontsize=8)
        ax.text(i + width / 2, o * 1.15, str(o), ha="center", fontsize=8)
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
