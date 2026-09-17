#!/usr/bin/env python3
"""Figure 1 v3 — the panel set the supervisors asked for on 2026-09-17.

Cole specified four panels (ask A1):

  (a) the admixture bar plot;
  (b) a class I ternary, "the best one you can make", HLA-A or HLA-B;
  (c) novel-allele rate by ancestry **broken out by gene** — "where is the reference data bias
      worst?";
  (d) the allele-frequency spectrum over **all** alleles with novel and known as two colours,
      overlaid, rather than three separate plots.

This script composes **(c) and (d)**, which are the two that can be rebuilt from committed
aggregates at S01's corrected novelty definition. Panels (a) and (b) are per-allele,
per-ancestry frequency products of scripts 06 and 10, which committed figures but not the
underlying tables, so they cannot be recomposed offline — and both need regenerating anyway at the
stricter admixture threshold Cole asked for separately (ask A9). `WS4_figure1_v3.md` records
exactly what that rerun needs.

Panel (c) is the one that changed most since the call. In the meeting the novel-allele rate was
one number per ancestry over everything flagged `new`. S01 showed that number is dominated by
non-coding novelty and by artifacts, so this panel splits the rate into the nomenclature fields
that actually mean different things:

    field 2  new protein     |  field 3  new synonymous CDS  |  field 4  new non-coding only

and draws the artifact rate as a **flat negative control** underneath. The artifact rate is flat
across ancestries by construction of the assembly pipeline, so if the field-2/3/4 rates have an
ancestry gradient and the artifact rate does not, the gradient is not an artifact of uneven
assembly quality. That control is the panel's argument, not decoration.

Usage:

    python3 scripts/hla_popgen/33_figure1_v3_compose.py
"""
import argparse
import json
import os
import re
import sys

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_THIS_DIR))
DEFAULT_REPORTS = os.path.join(_REPO_ROOT, "reports", "hla_popgen")
DEFAULT_FIELD_COUNTS = os.path.join(DEFAULT_REPORTS, "24_novelty_by_field",
                                    "novelty_field_counts.tsv")
DEFAULT_SFS = os.path.join(DEFAULT_REPORTS, "23_fig1_draft_panels", "sfs_known_vs_novel.tsv")
DEFAULT_OUT_DIR = os.path.join(DEFAULT_REPORTS, "33_figure1_v3")

SUPPRESS_BELOW = 20
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
ANCESTRY_COLORS = {"AFR": "#E69F00", "AMR": "#56B4E9", "EAS": "#009E73",
                   "EUR": "#0072B2", "MID": "#D55E00", "SAS": "#CC79A7"}
CLASSICAL = ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"]

FIELD_LABELS = {
    "f2_protein": "new protein (field 2)",
    "f3_synonymous": "new synonymous CDS (field 3)",
    "f4_noncoding": "new non-coding only (field 4)",
}
FIELD_ORDER = ["f2_protein", "f3_synonymous", "f4_noncoding"]

# Bin order for the SFS x-axis. The table stores bins as strings, so sorting them
# lexicographically would put "129-256" before "2" -- a plot that looks fine and is wrong.
SFS_BIN_ORDER = ["1", "2", "3-4", "5-8", "9-16", "17-32", "33-64", "65-128",
                 "129-256", "257-512", "513-1024", ">1024"]
SFS_BIN_WIDTH = {"1": 1, "2": 1, "3-4": 2, "5-8": 4, "9-16": 8, "17-32": 16, "33-64": 32,
                 "65-128": 64, "129-256": 128, "257-512": 256, "513-1024": 512, ">1024": 1024}


def _bounds(series):
    """(lower, upper, n_censored) total for a column of counts that may contain `<20`.

    **This function is why panel c is trustworthy.** 3,675 of the 8,784 cells in script 24's
    table are written `<20` under the AoU small-cell rule. Reading them with a plain
    `to_numeric(...).fillna(0)` turns every one of them into a zero -- which is what the first
    version of this script did, and it silently reported a 0.00% artifact rate for Middle Eastern
    ancestry across every gene, because MID is the smallest group and nearly all its artifact
    cells are censored. The plot looked completely normal.

    So a censored cell contributes [0, 19], never a point value, and the panel carries the
    resulting interval.
    """
    lo = hi = 0.0
    n_cens = 0
    for v in series:
        s = str(v).strip()
        if s.startswith("<"):
            try:
                cap = int(s[1:]) - 1
            except ValueError:
                cap = SUPPRESS_BELOW - 1
            hi += cap
            n_cens += 1
            continue
        try:
            x = float(s)
        except ValueError:
            continue
        lo += x
        hi += x
    return lo, hi, n_cens


def rate_by_gene_ancestry(counts, ancestry_scheme="strict", unrelated_only=True,
                          genes=None):
    """-> tidy frame: gene, ancestry, field_class, pct_lower/pct_upper, plus the artifact rate.

    The denominator is every haplotype the pipeline *called* at that gene in that ancestry --
    known plus novel, excluding `uncalled`. Using all rows as the denominator instead would let a
    gene that simply fails to be called more often in one ancestry masquerade as lower novelty
    there.

    Censoring is propagated as an interval, not a point: the rate's lower bound pairs the smallest
    possible numerator with the largest possible denominator, and vice versa.
    """
    df = counts.copy()
    df = df[(df["ancestry_scheme"].astype(str) == ancestry_scheme)
            & (df["unrelated_only"].astype(str).str.lower().isin({"true", "1"}) == bool(unrelated_only))
            & (df["ancestry"].isin(ANCESTRY_ORDER))]
    if genes:
        df = df[df["gene"].isin(genes)]

    called = df[df["field_class"] != "uncalled"]
    rows = []
    for (g, a), sub in called.groupby(["gene", "ancestry"]):
        den_lo, den_hi, den_cens = _bounds(sub["n_haplotypes"])
        if den_hi <= 0:
            continue
        clean = sub[sub["artifact_label"] == "clean"]
        for fc in FIELD_ORDER:
            lo, hi, nc = _bounds(clean[clean["field_class"] == fc]["n_haplotypes"])
            rows.append({"gene": g, "ancestry": a, "field_class": fc,
                         "n_called_lower": den_lo, "n_called_upper": den_hi,
                         "n_censored_cells": nc,
                         "pct_lower": 100.0 * lo / den_hi if den_hi else float("nan"),
                         "pct_upper": 100.0 * hi / den_lo if den_lo else float("nan")})
        lo, hi, nc = _bounds(sub[sub["artifact_label"] != "clean"]["n_haplotypes"])
        rows.append({"gene": g, "ancestry": a, "field_class": "artifact_control",
                     "n_called_lower": den_lo, "n_called_upper": den_hi,
                     "n_censored_cells": nc,
                     "pct_lower": 100.0 * lo / den_hi if den_hi else float("nan"),
                     "pct_upper": 100.0 * hi / den_lo if den_lo else float("nan")})
    return pd.DataFrame(rows)


def sfs_known_vs_novel(sfs, gene="ALL_CLASSICAL"):
    """Collapse script 23's four series into the two colours Cole asked for, density-normalised.

    The bins are powers of two of unequal width, so raw counts per bin produce a sawtooth that
    looks like structure and is purely binning. Dividing by bin width gives alleles per unit
    occurrence, which is comparable across bins.
    """
    d = sfs[sfs["gene"] == gene].copy()
    d["n_alleles"] = pd.to_numeric(d["n_alleles"], errors="coerce").fillna(0)
    d["group"] = np.where(d["series"] == "known", "known",
                          np.where(d["series"] == "flagged_artifact", "flagged", "novel"))
    out = (d.groupby(["occurrence_bin", "group"])["n_alleles"].sum()
           .unstack(fill_value=0.0).reindex(SFS_BIN_ORDER).fillna(0.0))
    width = pd.Series({b: SFS_BIN_WIDTH[b] for b in out.index})
    dens = out.div(width, axis=0)
    return out, dens


def compose(rates, counts_raw, sfs_counts, sfs_dens, out_png, out_pdf, gene_for_sfs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(13.0, 8.4))
    gs = GridSpec(2, 1, height_ratios=[1.25, 1.0], hspace=0.42, figure=fig)

    # ---- panel c: novel rate by gene x ancestry, split by nomenclature field ----------------
    axc = fig.add_subplot(gs[0])
    genes = [g for g in CLASSICAL if g in set(rates["gene"])]
    ancs = [a for a in ANCESTRY_ORDER if a in set(rates["ancestry"])]
    width = 0.8 / max(1, len(ancs))
    hatch = {"f2_protein": "///", "f3_synonymous": "...", "f4_noncoding": None}

    for ai, anc in enumerate(ancs):
        for gi, g in enumerate(genes):
            sub = rates[(rates["gene"] == g) & (rates["ancestry"] == anc)]
            x = gi - 0.4 + width * (ai + 0.5)
            bottom = 0.0
            top_hi = 0.0
            for fc in FIELD_ORDER:
                r = sub[sub["field_class"] == fc]
                v = float(r["pct_lower"].sum())
                top_hi += float(r["pct_upper"].sum())
                axc.bar(x, v, bottom=bottom, width=width * 0.9,
                        color=ANCESTRY_COLORS[anc],
                        alpha={"f4_noncoding": 1.0, "f3_synonymous": 0.72,
                               "f2_protein": 0.45}[fc],
                        hatch=hatch[fc], edgecolor="white", linewidth=0.3, zorder=3)
                bottom += v
            # How high the stack could go if every suppressed cell sat at its ceiling.
            if top_hi > bottom + 0.05:
                axc.plot([x, x], [bottom, top_hi], color="#555555", lw=0.7, zorder=4)
                axc.plot([x - width * 0.22, x + width * 0.22], [top_hi, top_hi],
                         color="#555555", lw=0.7, zorder=4)
            ar = sub[sub["field_class"] == "artifact_control"]
            art = float(ar["pct_lower"].sum())
            axc.plot([x - width * 0.45, x + width * 0.45], [art, art], color="#111111",
                     lw=1.1, zorder=5, solid_capstyle="butt")

    axc.set_xticks(range(len(genes)))
    axc.set_xticklabels([g.replace("HLA-", "") for g in genes], fontsize=10)
    axc.set_ylabel("% of called haplotypes carrying\nsequence absent from IPD-IMGT/HLA")
    axc.set_title("c   Where the reference is incomplete: novelty by gene, ancestry and "
                  "nomenclature field\n"
                  "stack = non-coding (solid) · synonymous CDS (dotted) · new protein (hatched) · "
                  "black bar = artifact rate (flat negative control) · grey whisker = how high the "
                  "stack could reach if every suppressed count sat at its ceiling", fontsize=10,
                  loc="left")
    axc.grid(axis="y", lw=0.4, alpha=0.35, zorder=0)
    axc.spines[["top", "right"]].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=ANCESTRY_COLORS[a]) for a in ancs]
    axc.legend(handles, ancs, frameon=False, ncol=len(ancs), fontsize=9,
               loc="upper left", bbox_to_anchor=(0.0, 1.0))

    # ---- panel d: SFS, known vs novel, one plot, two colours -------------------------------
    axd = fig.add_subplot(gs[1])
    x = np.arange(len(sfs_dens.index))
    for grp, color, lab in [("known", "#4C72B0", "already in IPD-IMGT/HLA"),
                            ("novel", "#C0392B", "absent from IPD-IMGT/HLA")]:
        if grp in sfs_dens.columns:
            axd.bar(x + (0.2 if grp == "novel" else -0.2), sfs_dens[grp], width=0.38,
                    color=color, label=lab, zorder=3)
    if "flagged" in sfs_dens.columns and sfs_dens["flagged"].sum() > 0:
        axd.plot(x, sfs_dens["flagged"], color="#777777", lw=1.0, ls="--", marker="o", ms=3,
                 label="flagged artifact (excluded from 'novel')", zorder=4)
    axd.set_yscale("log")
    axd.set_xticks(x)
    axd.set_xticklabels(sfs_dens.index, rotation=45, ha="right", fontsize=8)
    axd.set_xlabel("haplotypes carrying the allele")
    axd.set_ylabel("alleles per unit occurrence\n(density, log scale)")
    axd.set_title("d   Allele-frequency spectrum, %s — novel alleles are overwhelmingly rare, "
                  "known alleles are not\ncounts divided by bin width, because the bins double in "
                  "size and raw counts would show a sawtooth that is pure binning"
                  % gene_for_sfs.replace("ALL_CLASSICAL", "all 8 classical genes"),
                  fontsize=10, loc="left")
    axd.legend(frameon=False, fontsize=9)
    axd.grid(axis="y", lw=0.4, alpha=0.35, zorder=0)
    axd.spines[["top", "right"]].set_visible(False)

    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def write_readme(path, rates, gene_for_sfs):
    L = []
    L.append("# 33 — Figure 1 v3 (panels c and d)\n")
    L.append("*Supervisor ask A1 (Cole, 2026-09-17). Built from committed aggregates only.*\n")
    L.append("\n## What is here and what is not\n")
    L.append("| panel | Cole's ask | status |\n|---|---|---|\n")
    L.append("| a | admixture bar plot | **not in this file** — needs a rerun of script 06 at the "
             "stricter admixture threshold (ask A9) |\n")
    L.append("| b | class I ternary, HLA-A or HLA-B | **not in this file** — script 10 committed "
             "figures but not the per-allele frequency table, so it cannot be recomposed offline |\n")
    L.append("| c | novel rate by ancestry, broken out by gene | **here**, and rebuilt at S01's "
             "corrected novelty definition |\n")
    L.append("| d | SFS over all alleles, novel vs known in two colours | **here** |\n")

    L.append("\n## Panel c — read this before quoting a number from it\n")
    L.append("The rate presented in the call was a single 'novel' rate per ancestry. S01 showed "
             "that number is dominated by non-coding novelty and by artifacts, so this panel "
             "splits it by nomenclature field — new protein, new synonymous CDS, new non-coding "
             "only — and draws the artifact rate as a flat black bar underneath.\n")
    L.append("\n**Suppressed counts.** 3,675 of the 8,784 cells in script 24's table are written "
             "`<20` under the AoU small-cell rule. Every rate here is therefore an interval: the "
             "tabulated number is the **lower bound** (suppressed cells counted as zero) and the "
             "grey whisker in the figure shows where the stack would reach if every suppressed "
             "cell sat at 19. This matters most for MID, the smallest strict-ancestry group, "
             "where a naive reading of the suppressed cells as zero produces an artifact rate of "
             "exactly 0.00% across every gene — which is what the first version of this script "
             "reported, and it looked entirely plausible.\n")
    L.append("\nThe artifact bar is the control. Assembly and annotation artifacts are not "
             "ancestry-structured, so a gradient in the coloured stack against a flat black bar "
             "is reference incompleteness, not uneven assembly quality.\n")
    if not rates.empty:
        piv = (rates[rates["field_class"] != "artifact_control"]
               .pivot_table(index="gene", columns="ancestry", values="pct_lower", aggfunc="sum")
               .round(2))
        L.append("\n### Total clean novelty rate (%% of called haplotypes), all fields\n\n")
        L.append(piv.to_markdown() + "\n")
        art = (rates[rates["field_class"] == "artifact_control"]
               .pivot_table(index="gene", columns="ancestry", values="pct_lower").round(2))
        L.append("\n### Artifact rate — the control, expected to be flat across ancestries\n\n")
        L.append(art.to_markdown() + "\n")

    L.append("\n## Panel d\n")
    L.append("Spectrum for %s. Counts are divided by bin width: the bins double in size, so raw "
             "counts per bin produce a sawtooth that looks like structure and is entirely "
             "binning. Flagged-artifact alleles are drawn separately and are *not* counted as "
             "novel.\n" % gene_for_sfs)
    L.append("\n## Caveat on panel d's novelty definition\n")
    L.append("Panel d reuses script 23's spectrum, which predates S01's corrected split. Its "
             "'novel' series therefore still includes calls whose CDS turned out to be already "
             "catalogued. The shape of the spectrum — novel alleles collapsing after a couple of "
             "occurrences while known alleles reach hundreds — is unaffected, but the absolute "
             "allele counts in the novel series are an over-estimate. Regenerating it at the "
             "corrected definition needs the VM.\n")
    text = "\n".join(L)
    # The per-row appends above each end in a newline and the join adds another, which puts a
    # blank line between table rows -- markdown then renders the table as plain text.
    text = re.sub(r"\|\n\n\|", "|\n|", text)
    with open(path, "w") as fh:
        fh.write(text)


def run(args):
    os.makedirs(args.out_dir, exist_ok=True)
    for p in (args.field_counts, args.sfs):
        if not os.path.exists(p):
            sys.exit("FATAL: missing input %s" % p)
    counts = pd.read_csv(args.field_counts, sep="\t")
    sfs = pd.read_csv(args.sfs, sep="\t")

    rates = rate_by_gene_ancestry(counts, args.ancestry_scheme, True, CLASSICAL)
    if rates.empty:
        sys.exit("FATAL: no rows after filtering (scheme=%s). Check ancestry_scheme values."
                 % args.ancestry_scheme)
    sfs_counts, sfs_dens = sfs_known_vs_novel(sfs, args.sfs_gene)

    rates.to_csv(os.path.join(args.out_dir, "panel_c_novelty_rate.tsv"), sep="\t", index=False,
                 na_rep="NA")
    sfs_dens.to_csv(os.path.join(args.out_dir, "panel_d_sfs_density.tsv"), sep="\t", na_rep="NA")

    compose(rates, counts, sfs_counts, sfs_dens,
            os.path.join(args.out_dir, "figure1_v3_panels_cd.png"),
            os.path.join(args.out_dir, "figure1_v3_panels_cd.pdf"), args.sfs_gene)
    write_readme(os.path.join(args.out_dir, "README.md"), rates, args.sfs_gene)
    with open(os.path.join(args.out_dir, "summary.json"), "w") as fh:
        json.dump({"panels_built": ["c", "d"],
                   "panels_pending_vm": ["a_admixture", "b_class_I_ternary"],
                   "ancestry_scheme": args.ancestry_scheme,
                   "sfs_gene": args.sfs_gene}, fh, indent=2)
    print("[33] wrote panels c+d -> %s" % args.out_dir, flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--field-counts", default=DEFAULT_FIELD_COUNTS)
    ap.add_argument("--sfs", default=DEFAULT_SFS)
    ap.add_argument("--ancestry-scheme", choices=["strict", "pred"], default="strict")
    ap.add_argument("--sfs-gene", default="ALL_CLASSICAL")
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    run(ap.parse_args(argv))


if __name__ == "__main__":
    main()
