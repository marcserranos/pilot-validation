#!/usr/bin/env python3
"""LD supplement table + heatmap -- re-rendered from the COMMITTED 29_hla_ld tables only.

Cole (call #8, sprint S03): "this should be a table too -- a simple CSV with r-squared and D-
prime", with a distinct colour for "not observed / below n=20".

This script is offline-only: it reads `reports/hla_popgen/29_hla_ld/ld_pairwise.tsv`, which
`29_hla_ld_by_ancestry.py`'s `pairwise_ld()` already writes as a tidy per-allele-pair table (every
allele that individually clears the 20-haplotype floor in its own marginal, crossed against every
other such allele, per gene pair per ancestry) -- see that script's docstring and
`reports/hla_popgen/29_hla_ld/README.md`. Nothing here is re-computed from raw calls; every number
was already committed and small-cell suppressed by script 29.

**What the committed table already has, and what it doesn't:**
`ld_pairwise.tsv` has per-ancestry, per-allele-pair `r2`, `D`, `Dprime`, `freq_a`, `freq_b`,
`freq_hap`, and `n_hap_ij_disp` (the suppressed joint-haplotype count: `"0"`, `"<20"`, or the
real integer) for FIVE gene pairs: A~B, B~C, DPA1~DPB1, DQA1~DQB1, DRB1~DQB1. That is sufficient
for the full task: a tidy CSV with r2 and D-prime per pair, per ancestry, with a `status` column
distinguishing `not_observed` (joint count is exactly 0) from `suppressed_lt20` (joint count is
1-19, written `<20` and NEVER read as 0 or 20 -- see 32_novel_allele_callouts.py's
`parse_count()` for the same convention elsewhere in this project) from `estimated` (a real,
disclosable count). No data gap here: this script produces the full table for all five pairs, not
just DQA1~DQB1/DPA1~DPB1 -- those two are singled out only for the heatmap panel, per the task
(they're the ones a downstream HLA-DQ analysis cares about; A~B/B~C/DRB1~DQB1 are LD-decay
controls already summarized in 29's own report).

Usage (offline, no VM):
    python3 scripts/hla_popgen/37b_ld_supplement_table.py

Run: scripts/hla_popgen/tests/test_ld_supplement_table.py
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import _viz_common as vc

_REPO_ROOT = os.path.dirname(os.path.dirname(_THIS_DIR))
DEFAULT_LD_DIR = os.path.join(_REPO_ROOT, "reports", "hla_popgen", "29_hla_ld")
DEFAULT_OUT_DIR = os.path.join(_REPO_ROOT, "reports", "hla_popgen", "37_dq_g1g2_signed_ld")

# The two gene pairs the heatmap panel is scoped to (task instructions); the CSV covers every
# pair in the committed table regardless.
HEATMAP_PAIRS = ["DQA1~DQB1", "DPA1~DPB1"]

STATUS_NOT_OBSERVED = "not_observed"
STATUS_SUPPRESSED = "suppressed_lt20"
STATUS_ESTIMATED = "estimated"


def status_from_disp(v):
    """`n_hap_ij_disp` -> {not_observed, suppressed_lt20, estimated}. `0` is a disclosable,
    genuine zero (the allele pair truly never co-occurred on a phased haplotype in this
    ancestry); `<20` is 1-19 and MUST NOT be read as either 0 or 20 -- it is censored, not
    absent. See 32_novel_allele_callouts.py's parse_count()/is_censored() for the identical
    convention used elsewhere in this project."""
    s = str(v).strip()
    if s.startswith("<"):
        return STATUS_SUPPRESSED
    try:
        n = int(float(s))
    except ValueError:
        return STATUS_SUPPRESSED
    return STATUS_NOT_OBSERVED if n == 0 else STATUS_ESTIMATED


def n_hap_from_disp(v, status):
    """The real integer where disclosable (0 or >=20), else NaN -- never a guessed value for a
    suppressed cell."""
    if status == STATUS_SUPPRESSED:
        return np.nan
    try:
        return float(int(float(str(v).strip())))
    except ValueError:
        return np.nan


def load_pairwise(ld_dir):
    path = os.path.join(ld_dir, "ld_pairwise.tsv")
    if not os.path.exists(path):
        sys.exit(f"FATAL: {path} not found. This script reads the COMMITTED 29_hla_ld tables "
                 f"only -- run 29_hla_ld_by_ancestry.py on the VM first if it's missing.")
    df = pd.read_csv(path, sep="\t", dtype=str)
    need = {"pair", "ancestry", "allele_a", "allele_b", "freq_a", "freq_b", "r2", "Dprime",
            "n_hap_ij_disp"}
    missing = need - set(df.columns)
    if missing:
        sys.exit(f"FATAL: {path} is missing columns {sorted(missing)}. Committed table schema "
                 f"has changed; this script needs to be updated, not silently patched around.")
    for c in ("freq_a", "freq_b", "r2", "Dprime"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["status"] = df["n_hap_ij_disp"].map(status_from_disp)
    df["n_hap"] = df.apply(lambda r: n_hap_from_disp(r["n_hap_ij_disp"], r["status"]), axis=1)
    return df


def build_tidy_table(df):
    out = df.rename(columns={"pair": "gene_pair"})[
        ["gene_pair", "ancestry", "allele_a", "allele_b", "n_hap", "freq_a", "freq_b", "r2",
         "Dprime", "status"]
    ].copy()
    out = out.sort_values(["gene_pair", "ancestry", "allele_a", "allele_b"]).reset_index(drop=True)
    return out


# ---------------------------------------------------------------------------
# heatmap panel: DQA1~DQB1 and DPA1~DPB1, one subplot per (gene_pair, ancestry), signed D'
# colored on a diverging scale centered at 0, not_observed/suppressed cells hatched grey and kept
# OUT of the color scale (never folded in as if they were a real low value).
# ---------------------------------------------------------------------------
def draw_heatmap_grid(df, out_stem):
    sub = df[df["gene_pair"].isin(HEATMAP_PAIRS)].copy()
    if sub.empty:
        print("[37b] no rows for the heatmap gene pairs -- skipping the heatmap panel.",
              file=sys.stderr)
        return None

    ancestries = [a for a in vc.ANCESTRY_ORDER if a in set(sub["ancestry"])]
    n_rows, n_cols = len(HEATMAP_PAIRS), len(ancestries)

    with vc.nature_style():
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(
            n_rows, n_cols, figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(45 * n_rows + 10)),
            squeeze=False)
        cmap = vc.diverging_cmap()
        norm = vc.diverging_norm(-1.0, 1.0)
        im = None
        for ri, pair in enumerate(HEATMAP_PAIRS):
            pair_df = sub[sub["gene_pair"] == pair]
            alleles_a = sorted(pair_df["allele_a"].unique())
            alleles_b = sorted(pair_df["allele_b"].unique())
            for ci, anc in enumerate(ancestries):
                ax = axes[ri][ci]
                cell = pair_df[pair_df["ancestry"] == anc]
                mat = np.full((len(alleles_a), len(alleles_b)), np.nan)
                ai = {a: i for i, a in enumerate(alleles_a)}
                bi = {b: i for i, b in enumerate(alleles_b)}
                for _, r in cell.iterrows():
                    if r["status"] == STATUS_ESTIMATED:
                        mat[ai[r["allele_a"]], bi[r["allele_b"]]] = r["Dprime"]
                im = ax.imshow(mat, cmap=cmap, norm=norm, aspect="auto", zorder=2)
                # not_observed / suppressed cells: distinct grey hatch, never in the color scale
                for _, r in cell.iterrows():
                    if r["status"] != STATUS_ESTIMATED:
                        i0, j0 = ai[r["allele_a"]], bi[r["allele_b"]]
                        vc.hatch_suppressed(ax, j0 - 0.5, i0 - 0.5, 1, 1)
                ax.set_xticks([])
                ax.set_yticks([])
                if ri == 0:
                    ax.set_title(anc, fontsize=6)
                if ci == 0:
                    ax.set_ylabel(pair, fontsize=6)
        if im is not None:
            cbar = fig.colorbar(im, ax=axes, shrink=0.6, pad=0.01, label="D'  (signed)")
            cbar.ax.tick_params(labelsize=5)
            cbar.set_label("D' (signed)", fontsize=6)
        from matplotlib.patches import Patch
        legend_handles = [Patch(facecolor=vc.SUPPRESSED_COLOR, edgecolor="#999999", hatch="////",
                                label="not observed / suppressed (n<20)")]
        fig.legend(handles=legend_handles, loc="lower center", frameon=False, fontsize=5.5,
                  bbox_to_anchor=(0.5, 0.0), ncol=1)
        fig.suptitle("Pairwise D' (signed), by ancestry -- rows: alleles of the first gene, "
                     "columns: alleles of the second", fontsize=7, y=0.995)
    return vc.save_fig(fig, out_stem)


README_TEMPLATE = """# 37 -- LD supplement table (r2 / D') and DQ/DP heatmap

Source: script `37b_ld_supplement_table.py`, re-rendered from the already-committed
`reports/hla_popgen/29_hla_ld/ld_pairwise.tsv`. No VM access, no re-computation from raw calls --
every number here was already public.

Cole (call #8): "this should be a table too -- a simple CSV with r2 and D-prime", with a distinct
colour for "not observed / below n=20".

## Data-gap check (per the task)

The committed `ld_pairwise.tsv` already carries the **full per-allele-pair signed table**: for
every allele pair (per gene pair, per ancestry) where BOTH alleles individually clear the
20-haplotype floor in their own marginal, it has `r2`, `D`, `Dprime`, `freq_a`, `freq_b`,
`freq_hap`, and the suppressed joint-haplotype count `n_hap_ij_disp`. This is **not** a top-12-only
or multiallelic-summary-only table (compare `ld_multiallelic.tsv`, which IS a whole-table summary
and was not sufficient on its own) -- so this script produces the complete CSV for **all five**
committed gene pairs (A~B, B~C, DPA1~DPB1, DQA1~DQB1, DRB1~DQB1), not only the two the heatmap is
scoped to. **No gap; nothing here is deferred to the VM version.**

## `supp_table_ld_pairs.csv`

Tidy, one row per (gene_pair, ancestry, allele_a, allele_b): `gene_pair`, `ancestry`, `allele_a`,
`allele_b`, `n_hap` (the joint haplotype count, blank/NaN where suppressed), `freq_a`, `freq_b`,
`r2`, `Dprime`, `status`.

`status` is one of:
- **`estimated`** -- the joint haplotype count is disclosable (>=20 or a genuine 0) and `r2`/
  `Dprime` are estimated from real data.
- **`not_observed`** -- the joint count is exactly 0 (this allele pair never co-occurred on a
  phased haplotype in this ancestry). Disclosable; not suppressed.
- **`suppressed_lt20`** -- the joint count is 1-19 under the AoU small-cell rule, written `<20` in
  the source table. **Parsed as censored throughout -- never read as 0 and never as 20.**

## Heatmap panel (DQA1~DQB1, DPA1~DPB1)

One small heatmap per (gene pair, ancestry): rows are the gene pair's first-locus alleles,
columns the second-locus alleles, colour is signed D' on a blue-white-red diverging scale
centered at 0 (`_viz_common.diverging_cmap()`/`diverging_norm()`). Cells with status
`not_observed` or `suppressed_lt20` are drawn as a distinct light-grey hatched cell
(`_viz_common.hatch_suppressed()` / `SUPPRESSED_COLOR`) and are **excluded from the colour scale
entirely** -- they never appear as if they were a real, low-magnitude D' value.

**How to read it:** a solid block of one colour within a gene pair/ancestry panel indicates
strong, consistently-signed linkage; hatched cells mean the underlying allele pair either never
co-occurred in this ancestry's phased haplotypes or its count fell under the 20-haplotype
disclosure floor -- in both cases, no D' is estimated or shown for that cell.

## Caveats

- This is a restyle + tabulation of already-committed, already-suppressed numbers -- no new
  statistic is computed here.
- r2/D' at low allele frequency are noisy; see `29_hla_ld/README.md`'s own caveats on rare-allele
  LD estimates (which this script does not repeat or override).
- `A~B`/`B~C`/`DRB1~DQB1` are in the CSV but not the heatmap panel, which is scoped to
  DQA1~DQB1/DPA1~DPB1 per the task.
"""


def write_readme(out_dir):
    path = os.path.join(out_dir, "README.md")
    with open(path, "w") as f:
        f.write(README_TEMPLATE)
    print(f"  wrote {path}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ld-dir", default=DEFAULT_LD_DIR,
                     help="Directory with the committed 29_hla_ld tables.")
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    args = ap.parse_args()

    vc.ensure_dir(args.out_dir)

    raw = load_pairwise(args.ld_dir)
    tidy = build_tidy_table(raw)

    csv_path = os.path.join(args.out_dir, "supp_table_ld_pairs.csv")
    tidy.to_csv(csv_path, index=False)
    print(f"  wrote {csv_path} ({len(tidy)} rows)", file=sys.stderr)

    heatmap_stem = os.path.join(args.out_dir, "supp_heatmap_dq_dp")
    draw_heatmap_grid(tidy, heatmap_stem)

    write_readme(args.out_dir)
    print(f"[37b] done -> {args.out_dir}")


if __name__ == "__main__":
    main()
