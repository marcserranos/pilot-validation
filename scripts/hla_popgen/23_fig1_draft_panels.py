#!/usr/bin/env python3
"""Draft Figure 1 panels (supervisor action items 3, 4, 5) from committed aggregate tables only.

Reads `novel_alleles.tsv`, `allele_richness.tsv` and `frequency_spectrum.tsv` from
`reports/hla_popgen/`. No person-level data and no VM round-trip -- runs on a laptop.

WHY NOVELTY IS SPLIT INTO THREE BUCKETS
---------------------------------------
Pooled "novel-call rate" hides that, for the classical genes, most novel haplotypes carry a CDS that
is byte-identical to a known allele and differ only outside it (`novelty_class == beyond_cds`): 96% of
DRB1's novel haplotypes, ~86% across the 8 classical genes. IPD-IMGT/HLA's full-genomic (intron/UTR)
sequences are far less complete than its exon sequences, so non-coding novelty measures catalogue
gaps as much as unseen biology. Reporting one number would let DRB1's intron novelty dominate
Figure 1. Every panel here keeps three buckets apart:

  cds_changing     protein_altering or synonymous, not flagged as an artifact
  noncoding_only   beyond_cds, not flagged as an artifact
  flagged_artifact confidence_tier == flagged_artifact (any novelty_class)

APPROXIMATION IN PANEL A (stated, not hidden)
---------------------------------------------
`ancestry_counts` counts carrier PEOPLE per cluster, not haplotypes. A person homozygous for a novel
cluster contributes 1, not 2, so the per-ancestry numerator is a slight undercount. The pooled
per-gene haplotype totals are reported next to the richness table's novel haplotype count as a
reconciliation check. The exact version is one groupby over Table 1 on the VM.

DISCLOSURE FLOOR
----------------
The AoU Data and Statistics Dissemination Policy forbids publishing participant counts of 1-19. The
novel-allele list writes such counts as "<20". Rates and allele counts (not participant counts) are
written as-is.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _viz_common import (ANCESTRY_COLORS, ANCESTRY_ORDER, CLASSICAL_GENES,  # noqa: E402
                         DEFAULT_REPORT_ROOT, plt, savefig, write_report)

BUCKETS = ["cds_changing", "noncoding_only", "flagged_artifact"]
BUCKET_LABELS = {"cds_changing": "CDS-changing", "noncoding_only": "non-coding only",
                 "flagged_artifact": "flagged artifact"}
DISCLOSURE_FLOOR = 20
# Power-of-2 bins; plotted counts are divided by bin width (number of integer k values in the bin)
# so unequal widths don't draw a sawtooth. The last bin is open-ended and its width is taken as
# 1024 (k = 1024..2047), the next power of 2.
SFS_BIN_EDGES = [1, 2, 3, 5, 9, 17, 33, 65, 129, 257, 513, 1025, np.inf]
SFS_BIN_LABELS = ["1", "2", "3-4", "5-8", "9-16", "17-32", "33-64", "65-128", "129-256",
                  "257-512", "513-1024", ">1024"]
SFS_BIN_WIDTHS = [1, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]


def novelty_bucket(row):
    if row["confidence_tier"] == "flagged_artifact":
        return "flagged_artifact"
    if row["novelty_class"] == "beyond_cds":
        return "noncoding_only"
    if row["novelty_class"] in ("protein_altering", "synonymous"):
        return "cds_changing"
    return None  # undetermined: counted in the report, excluded from panels


def load_novel(path):
    novel = pd.read_csv(path, sep="\t")
    novel["bucket"] = novel.apply(novelty_bucket, axis=1)
    novel["ancestry_counts"] = novel["ancestry_counts"].apply(
        lambda s: json.loads(s) if isinstance(s, str) and s.strip() else {})
    return novel


def novel_rate_table(novel, richness, genes):
    """Carrier incidences per 100 haplotypes, per gene x ancestry x bucket."""
    denom = (richness[richness["category"] == "all"]
             .set_index(["gene", "ancestry"])["n_haplotypes"])
    rows = []
    sub = novel[novel["gene"].isin(genes) & novel["bucket"].notna()]
    for (gene, bucket), grp in sub.groupby(["gene", "bucket"]):
        per_anc = {a: 0 for a in ANCESTRY_ORDER}
        for counts in grp["ancestry_counts"]:
            for anc, n in counts.items():
                if anc in per_anc:
                    per_anc[anc] += n
        for anc, n in per_anc.items():
            d = denom.get((gene, anc), np.nan)
            rows.append({"gene": gene, "ancestry": anc, "bucket": bucket, "carrier_incidences": n,
                         "n_haplotypes": d, "per_100_haplotypes": 100 * n / d if d else np.nan})
    out = pd.DataFrame(rows)
    full = pd.MultiIndex.from_product([genes, ANCESTRY_ORDER, BUCKETS],
                                      names=["gene", "ancestry", "bucket"])
    out = out.set_index(["gene", "ancestry", "bucket"]).reindex(full).reset_index()
    out["carrier_incidences"] = out["carrier_incidences"].fillna(0).astype(int)
    out["n_haplotypes"] = [denom.get((g, a), np.nan) for g, a in zip(out["gene"], out["ancestry"])]
    out["per_100_haplotypes"] = 100 * out["carrier_incidences"] / out["n_haplotypes"]
    return out


def reconciliation_table(novel, richness, genes):
    pooled = (richness[(richness["category"] == "novel") & (richness["ancestry"] == "POOLED")]
              .set_index("gene")["n_haplotypes"])
    sub = novel[novel["gene"].isin(genes)]
    t = sub.pivot_table(index="gene", columns="bucket", values="n_haplotypes", aggfunc="sum",
                        fill_value=0).reindex(columns=BUCKETS, fill_value=0)
    t["undetermined"] = sub[sub["bucket"].isna()].groupby("gene")["n_haplotypes"].sum()
    t["undetermined"] = t["undetermined"].fillna(0).astype(int)
    t["sum_all_buckets"] = t[BUCKETS + ["undetermined"]].sum(axis=1)
    t["richness_novel_haplotypes"] = pooled
    clean = t["cds_changing"] + t["noncoding_only"]
    t["noncoding_share_of_clean_%"] = (100 * t["noncoding_only"] / clean).round(1)
    return t.reindex(genes)


def bin_sfs(counts_by_k):
    """counts_by_k: Series indexed by occurrence k -> n_alleles. Returns Series over SFS bins."""
    k = counts_by_k.index.to_numpy(dtype=float)
    idx = np.digitize(k, SFS_BIN_EDGES) - 1
    binned = pd.Series(0, index=SFS_BIN_LABELS, dtype=int)
    for i, n in zip(idx, counts_by_k.to_numpy()):
        if 0 <= i < len(SFS_BIN_LABELS):
            binned.iloc[i] += int(n)
    return binned


def sfs_table(novel, spectrum, genes):
    """Allele counts per occurrence bin: known (from 04's spectrum) vs novel by bucket."""
    rows = []
    for gene in genes + ["ALL_CLASSICAL"]:
        gsel = genes if gene == "ALL_CLASSICAL" else [gene]
        known = spectrum[(spectrum["gene"].isin(gsel)) & (spectrum["category"] == "known")
                         & (spectrum["ancestry"] == "POOLED")]
        series = {"known": bin_sfs(known.groupby("occurrence_k")["n_alleles"].sum())}
        for bucket in BUCKETS:
            nb = novel[novel["gene"].isin(gsel) & (novel["bucket"] == bucket)]
            series[bucket] = bin_sfs(nb.groupby("n_haplotypes").size())
        for label, s in series.items():
            for (b, n), w in zip(s.items(), SFS_BIN_WIDTHS):
                rows.append({"gene": gene, "series": label, "occurrence_bin": b, "n_alleles": n,
                             "alleles_per_k": n / w})
    return pd.DataFrame(rows)


def suppress(n):
    return f"<{DISCLOSURE_FLOOR}" if 0 < n < DISCLOSURE_FLOOR else str(int(n))


def top_novel_list(novel, genes, bucket, min_persons):
    sub = novel[novel["gene"].isin(genes) & (novel["bucket"] == bucket)
                & novel["confidence_tier"].isin(["high", "recurrent"])
                & (novel["n_persons"] >= min_persons)]
    sub = sub.sort_values("n_persons", ascending=False)
    out = pd.DataFrame({
        "novel_id": sub["novel_id"], "gene": sub["gene"], "nearest_allele": sub["nearest_allele"],
        "novelty_class": sub["novelty_class"], "cds_distance": sub["cds_distance"],
        "n_aa_changes": sub["n_aa_changes"], "n_persons": sub["n_persons"].map(suppress),
        "n_ancestries": sub["n_ancestries"],
    })
    for anc in ANCESTRY_ORDER:
        out[f"n_{anc}"] = sub["ancestry_counts"].map(lambda c, a=anc: suppress(c.get(a, 0)))
    return out


def plot_novel_rate(rates, genes, path):
    fig, axes = plt.subplots(2, 4, figsize=(14, 6.5), sharey=False)
    x = np.arange(len(ANCESTRY_ORDER))
    styles = {"cds_changing": dict(hatch=None, alpha=1.0),
              "noncoding_only": dict(hatch=None, alpha=0.45),
              "flagged_artifact": dict(hatch="///", alpha=0.15)}
    for ax, gene in zip(axes.flat, genes):
        g = rates[rates["gene"] == gene].set_index(["ancestry", "bucket"])["per_100_haplotypes"]
        bottom = np.zeros(len(ANCESTRY_ORDER))
        for bucket in BUCKETS:
            vals = np.array([g.get((a, bucket), 0.0) for a in ANCESTRY_ORDER], dtype=float)
            vals = np.nan_to_num(vals)
            ax.bar(x, vals, bottom=bottom, color=[ANCESTRY_COLORS[a] for a in ANCESTRY_ORDER],
                   edgecolor="black", linewidth=0.4, **styles[bucket])
            bottom += vals
        ax.set_title(gene, fontsize=10)
        ax.set_xticks(x)
        ax.set_xticklabels(ANCESTRY_ORDER, fontsize=7)
        ax.tick_params(axis="y", labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("novel carrier incidences\nper 100 haplotypes", fontsize=8)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor="#777777", edgecolor="black",
                             **styles[b]) for b in BUCKETS]
    fig.legend(handles, [BUCKET_LABELS[b] for b in BUCKETS], loc="upper center", ncol=3,
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    savefig(fig, path, dpi=200)


def plot_sfs(sfs, genes, path):
    panels = ["ALL_CLASSICAL"] + genes
    fig, axes = plt.subplots(3, 3, figsize=(13, 10))
    colors = {"known": "#222222", "cds_changing": "#C0392B", "noncoding_only": "#2E86C1",
              "flagged_artifact": "#AAAAAA"}
    labels = {"known": "known (IPD-IMGT)", **BUCKET_LABELS}
    x = np.arange(len(SFS_BIN_LABELS))
    for ax, gene in zip(axes.flat, panels):
        g = sfs[sfs["gene"] == gene]
        for series in ["known"] + BUCKETS:
            s = g[g["series"] == series].set_index("occurrence_bin").reindex(SFS_BIN_LABELS)
            y = s["alleles_per_k"].to_numpy(dtype=float, copy=True)
            y[y == 0] = np.nan
            ax.plot(x, y, marker="o", ms=3.5, lw=1.4, color=colors[series],
                    ls="--" if series == "flagged_artifact" else "-", label=labels[series])
        ax.set_yscale("log")
        ax.set_title("8 classical genes pooled" if gene == "ALL_CLASSICAL" else gene, fontsize=10)
        ax.set_xticks(x)
        ax.set_xticklabels(SFS_BIN_LABELS, fontsize=6.5, rotation=45)
        ax.tick_params(axis="y", labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[-1, :]:
        ax.set_xlabel("haplotypes carrying the allele", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("distinct alleles per k\n(bin-width normalized)", fontsize=8)
    handles, lbls = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, lbls, loc="upper center", ncol=4, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, 1.01))
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    savefig(fig, path, dpi=200)


def md_table(df, max_rows=None):
    df = df if max_rows is None else df.head(max_rows)
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--report-root", default=DEFAULT_REPORT_ROOT)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--min-persons", type=int, default=11,
                    help="list recurrent novel alleles carried by at least this many people")
    args = ap.parse_args()
    root = args.report_root
    out = args.out_dir or os.path.join(root, "23_fig1_draft_panels")
    genes = CLASSICAL_GENES

    novel = load_novel(os.path.join(root, "novel_alleles.tsv"))
    richness = pd.read_csv(os.path.join(root, "allele_richness.tsv"), sep="\t")
    spectrum = pd.read_csv(os.path.join(root, "frequency_spectrum.tsv"), sep="\t")

    rates = novel_rate_table(novel, richness, genes)
    recon = reconciliation_table(novel, richness, genes)
    sfs = sfs_table(novel, spectrum, genes)
    cds_list = top_novel_list(novel, genes, "cds_changing", min_persons=2)
    nc_list = top_novel_list(novel, genes, "noncoding_only", args.min_persons)

    os.makedirs(out, exist_ok=True)
    rates.to_csv(os.path.join(out, "novel_rate_by_gene_ancestry.tsv"), sep="\t", index=False)
    recon.to_csv(os.path.join(out, "novel_haplotype_reconciliation.tsv"), sep="\t")
    sfs.to_csv(os.path.join(out, "sfs_known_vs_novel.tsv"), sep="\t", index=False)
    cds_list.to_csv(os.path.join(out, "recurrent_novel_cds_changing.tsv"), sep="\t", index=False)
    nc_list.to_csv(os.path.join(out, "recurrent_novel_noncoding.tsv"), sep="\t", index=False)
    plot_novel_rate(rates, genes, os.path.join(out, "panel_novel_rate_by_gene_ancestry.png"))
    plot_sfs(sfs, genes, os.path.join(out, "panel_sfs_known_vs_novel.png"))

    wide = (rates.pivot_table(index=["gene", "bucket"], columns="ancestry",
                              values="per_100_haplotypes").reindex(columns=ANCESTRY_ORDER).round(2)
            .reset_index())
    lines = [
        "# Draft Figure 1 panels -- novel-allele rate, SFS, recurrent novel alleles",
        "",
        "Generated by `scripts/hla_popgen/23_fig1_draft_panels.py` from committed aggregate tables "
        "only. See the module docstring for the three-bucket split and the Panel A approximation.",
        "",
        "## Reconciliation: novel haplotypes by bucket (pooled, classical genes)",
        "",
        *md_table(recon.reset_index()),
        "",
        "`sum_all_buckets` should sit within ~1-2% of `richness_novel_haplotypes` (different "
        "join paths in 03 vs 04).",
        "",
        "## Panel A -- novel carrier incidences per 100 haplotypes",
        "",
        *md_table(wide),
        "",
        "![novel rate](panel_novel_rate_by_gene_ancestry.png)",
        "",
        "## Panel B -- frequency spectrum, known vs novel",
        "",
        "![sfs](panel_sfs_known_vs_novel.png)",
        "",
        f"## Recurrent CDS-changing novel alleles (>=2 people): {len(cds_list)}",
        "",
        *md_table(cds_list),
        "",
        f"## Recurrent non-coding-only novel alleles (>={args.min_persons} people): {len(nc_list)}",
        "",
        "Top 30 shown; full list in `recurrent_novel_noncoding.tsv`.",
        "",
        *md_table(nc_list, max_rows=30),
        "",
        f"Participant counts of 1-{DISCLOSURE_FLOOR - 1} are written as `<{DISCLOSURE_FLOOR}` "
        "(AoU dissemination policy).",
    ]
    write_report(os.path.join(out, "fig1_draft_panels_report.md"), lines)


if __name__ == "__main__":
    main()
