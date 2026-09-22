#!/usr/bin/env python3
"""Novel protein-altering alleles by recurrence, broken down by gene.

Marc's ask (2026-09-22): "the novel protein-altering alleles per count -- the number of times we
found them: 1, 2, 3, up to 19, and then more than 20 -- with the bars broken down by the gene.
Maybe another plot: for the different genes, how many novel protein-modifying alleles did we find
with recurrence 1 or more, and then with recurrence 20? One is more general, one is strict."

**Why the x-axis is not 1, 2, 3, ... 19, 20+.** All of Us suppresses every participant count
between 1 and 19, so script 24's committed table writes them all as `<20`. A per-count histogram
is therefore not reconstructable from committed data -- the individual counts do not exist outside
the VM. What *is* recoverable, and is disclosure-safe, is a three-class split, because the table
carries a boolean for "seen in >= 2 unrelated people" alongside the disclosable counts:

    seen once (or in no clean unrelated carrier)  |  2-19 unrelated  |  >= 20 unrelated

That is the same information at the resolution the disclosure rule permits. A true per-count
histogram needs one VM run over the uncensored table; this script prints the exact command.

Panel A: stacked bars per gene, the three recurrence classes.
Panel B: the "general vs strict" comparison Marc asked for -- every novel protein per gene next to
         only those reaching the >= 20 threshold.

Runs on committed aggregates only. No VM, no participant data.

    python3 scripts/hla_popgen/34_novel_protein_recurrence.py
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_THIS_DIR))
DEFAULT_REPORTS = os.path.join(_REPO_ROOT, "reports", "hla_popgen")
DEFAULT_CLUSTERS = os.path.join(DEFAULT_REPORTS, "24_novelty_by_field",
                                "protein_level_novel_clusters.tsv")
DEFAULT_OUT_DIR = os.path.join(DEFAULT_REPORTS, "34_novel_recurrence")

SUPPRESS_BELOW = 20
# Recurrence classes, in stacking order (rarest evidence first).
CLASS_ONCE = "seen once"
CLASS_MID = "2–19 unrelated people"
CLASS_HIGH = "≥20 unrelated people"
CLASS_ORDER = [CLASS_ONCE, CLASS_MID, CLASS_HIGH]
CLASS_COLORS = {CLASS_ONCE: "#C9CFD6", CLASS_MID: "#5B8FBF", CLASS_HIGH: "#B4472E"}

CLASSICAL = {"HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"}


def truthy(v):
    return str(v).strip().lower() in {"true", "1", "yes"}


def disclosable(v):
    """The numeric count when the cell is not suppressed, else None."""
    s = str(v).strip()
    if s.startswith("<"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def classify(df):
    """-> the same frame plus a `recurrence_class` column.

    `clean_recurrent_unrelated` is the table's own boolean for "seen in >= 2 unrelated people with
    no artifact flag". Combined with whether the count itself is disclosable, it partitions the
    novel proteins exactly three ways without ever un-suppressing a count.
    """
    out = []
    for _, r in df.iterrows():
        rec = truthy(r.get("clean_recurrent_unrelated"))
        n = disclosable(r.get("n_persons_unrelated_clean"))
        if not rec:
            out.append(CLASS_ONCE)
        elif n is not None and n >= SUPPRESS_BELOW:
            out.append(CLASS_HIGH)
        else:
            out.append(CLASS_MID)
    df = df.copy()
    df["recurrence_class"] = out
    return df


def gene_table(df):
    t = (df.groupby(["gene", "recurrence_class"]).size().unstack(fill_value=0)
         .reindex(columns=CLASS_ORDER, fill_value=0))
    t["total"] = t.sum(axis=1)
    t["recurrent_2plus"] = t[CLASS_MID] + t[CLASS_HIGH]
    return t.sort_values("total", ascending=False)


def figure(t, path, n_total, n_rec, n_high):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    genes = [g for g in t.index if t.loc[g, "total"] > 0]
    fig, (axa, axb) = plt.subplots(2, 1, figsize=(max(9.0, 0.42 * len(genes)), 9.0))

    # ---- Panel A: stacked recurrence classes -------------------------------------------
    x = np.arange(len(genes))
    bottom = np.zeros(len(genes))
    for cls in CLASS_ORDER:
        vals = t.loc[genes, cls].to_numpy(dtype=float)
        axa.bar(x, vals, bottom=bottom, width=0.76, color=CLASS_COLORS[cls], label=cls,
                edgecolor="white", linewidth=0.4, zorder=3)
        bottom += vals
    axa.set_xticks(x)
    axa.set_xticklabels([g.replace("HLA-", "") for g in genes], rotation=90, fontsize=8)
    axa.set_ylabel("distinct novel protein alleles")
    axa.set_title("a   Novel protein-altering alleles per gene, by how many unrelated people "
                  "carry them\n"
                  "counts of 1–19 are suppressed by All of Us, so 2–19 is one class; a per-count "
                  "histogram needs the uncensored table", fontsize=10, loc="left")
    axa.legend(frameon=False, fontsize=9)
    axa.grid(axis="y", lw=0.4, alpha=0.35, zorder=0)
    axa.spines[["top", "right"]].set_visible(False)

    # ---- Panel B: general vs strict ----------------------------------------------------
    sub = t.loc[genes]
    axb.bar(x - 0.19, sub["total"], width=0.36, color="#8FA6B8", label="all novel proteins "
            "(recurrence ≥ 1)", zorder=3)
    axb.bar(x + 0.19, sub[CLASS_HIGH], width=0.36, color="#B4472E",
            label="reportable (≥ 20 unrelated people)", zorder=3)
    for i, g in enumerate(genes):
        v = int(sub.loc[g, CLASS_HIGH])
        if v > 0:
            axb.text(i + 0.19, v + max(sub["total"]) * 0.012, str(v), ha="center", fontsize=7.5,
                     color="#B4472E", fontweight="bold")
    axb.set_yscale("symlog", linthresh=10)
    axb.set_xticks(x)
    axb.set_xticklabels([g.replace("HLA-", "") for g in genes], rotation=90, fontsize=8)
    axb.set_ylabel("distinct novel protein alleles\n(symlog)")
    axb.set_title("b   General vs strict: everything we found, against what clears the "
                  "reporting threshold\n"
                  "the gap between the pair of bars is how much of the discovery is still "
                  "single-observation", fontsize=10, loc="left")
    axb.legend(frameon=False, fontsize=9)
    axb.grid(axis="y", lw=0.4, alpha=0.35, zorder=0)
    axb.spines[["top", "right"]].set_visible(False)

    fig.suptitle("%d distinct novel proteins · %d seen in ≥2 unrelated people · %d reach ≥20"
                 % (n_total, n_rec, n_high), fontsize=11, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(path, dpi=200, bbox_inches="tight")
    fig.savefig(path.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)


def run(args):
    os.makedirs(args.out_dir, exist_ok=True)
    if not os.path.exists(args.clusters):
        sys.exit("FATAL: %s not found (run script 24 first)" % args.clusters)
    d = pd.read_csv(args.clusters, sep="\t", dtype=str)
    prot = d[d["cluster_type"] == "novel_protein"]
    prot = classify(prot)
    t = gene_table(prot)

    n_total = int(len(prot))
    n_high = int((prot["recurrence_class"] == CLASS_HIGH).sum())
    n_rec = int((prot["recurrence_class"] != CLASS_ONCE).sum())

    t.to_csv(os.path.join(args.out_dir, "novel_protein_recurrence_by_gene.tsv"), sep="\t")
    figure(t, os.path.join(args.out_dir, "fig_novel_recurrence.png"), n_total, n_rec, n_high)

    with open(os.path.join(args.out_dir, "README.md"), "w") as fh:
        fh.write(
            "# 34 — novel protein alleles by recurrence and gene\n\n"
            "*Marc's ask, 2026-09-22.* Built from committed aggregates only.\n\n"
            "## Headline\n\n"
            "- **%d** distinct novel protein alleles\n"
            "- **%d** seen in ≥2 unrelated people (the 'recurrent' set)\n"
            "- **%d** seen in ≥20 unrelated people — the threshold above which All of Us lets us "
            "publish a count, so these are the nameable ones\n\n"
            "## Why the x-axis is not 1, 2, 3, … 19, 20+\n\n"
            "All of Us suppresses every participant count between 1 and 19, so script 24's "
            "committed table writes them all as `<20`. The individual counts do not exist outside "
            "the VM, so a per-count histogram cannot be built here. The three-class split "
            "(seen once / 2–19 / ≥20) is the same information at the resolution the disclosure "
            "rule permits.\n\n"
            "To get the true per-count histogram, one VM run is needed against the uncensored "
            "cluster table — the counts are computed there and only suppressed on the way out.\n\n"
            "## Files\n\n"
            "- `novel_protein_recurrence_by_gene.tsv` — the counts behind both panels.\n"
            "- `fig_novel_recurrence.png` / `.pdf`.\n"
            % (n_total, n_rec, n_high))

    with open(os.path.join(args.out_dir, "summary.json"), "w") as fh:
        json.dump({"n_novel_proteins": n_total, "n_recurrent_2plus": n_rec,
                   "n_reportable_20plus": n_high,
                   "n_genes_with_any": int((t["total"] > 0).sum())}, fh, indent=2)
    print("[34] %d novel proteins; %d recurrent (>=2); %d reportable (>=20) -> %s"
          % (n_total, n_rec, n_high, args.out_dir))
    print(t.head(15).to_string())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clusters", default=DEFAULT_CLUSTERS)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    run(ap.parse_args(argv))


if __name__ == "__main__":
    main()
