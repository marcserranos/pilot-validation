#!/usr/bin/env python3
"""Experiment D -- field-level cascade comparison (AoU-native and SpecHLA vs SpecImmune-LR truth).

Supersedes the coarser EXACT/CLOSE/PARTIAL_MISS/FULL_MISS buckets in
analyze_experiment_d_resolution.py with the real HLA nomenclature fields, per Marc's
2026-07-20 feedback: supervisors want to know exactly which field a discrepancy occurs at,
not an invented category name, and Field 2 (the specific protein -- non-synonymous) is the
threshold that actually matters; Field 3 (synonymous, coding) and Field 4 (non-coding) matter
much less but should still be visible, not hidden.

For every SR-vs-truth allele pair (AoU-vs-SpecImmune, SpecHLA-vs-SpecImmune), and for each
field level 1..4, this classifies the comparison as:
  True   -- this field and all shallower fields match
  False  -- this field or a shallower one is a real mismatch (cascade already broken)
  None   -- not assessable: one of the two calls simply doesn't report this many fields
            (e.g. AoU-native is structurally capped at 3 fields; SpecImmune sometimes only
            resolves to 2 fields for a given call) -- NOT counted as a disagreement.

This uses each call's *actual* reported field depth (no truncation to 2-field like the other
Experiment D scripts) -- the raw calls already carry up to 4 fields natively (SpecHLA emits
e.g. A*24:02:01:01; SpecImmune varies 2-4 fields per call depending on its own confidence).

Reads the same per-person comparison_log.csv files as the other Experiment D analysis scripts
(run_label == "experiment_d"). Aggregate output only (counts, no genotypes) -- keep on the VM
per the egress caveat (DECISIONS.md), paste tables back / bring PNGs back as needed.

Usage (via `pixi run -e spechla` for pandas+matplotlib):
  python3 analyze_experiment_d_field_cascade.py ~/pipeline_outputs/experiment_d/cohort.tsv
"""
import argparse
import os
import sys

import pandas as pd

import matplotlib
matplotlib.use("Agg")  # no display server on the VM
import matplotlib.pyplot as plt

GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
GROUPS = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
METHODS = ["aou", "sh"]
METHOD_LABEL = {"aou": "AoU-native", "sh": "SpecHLA"}
NULL = {"", "NA", "-", "nan", "None", ".", "na"}
FIELD_LABEL = {
    1: "Field 1 (allele group)",
    2: "Field 2 (protein -- non-synonymous)",
    3: "Field 3 (synonymous, coding)",
    4: "Field 4 (non-coding)",
}


def parse_fields(raw):
    """'A*24:02:01:01' -> ['24','02','01','01']; None if null/uncallable."""
    if raw is None:
        return None
    s = str(raw).strip()
    if s in NULL or (hasattr(pd, "isna") and pd.isna(raw)):
        return None
    if "*" in s:
        s = s.split("*", 1)[1]
    fields = [f for f in s.split(":") if f != ""]
    return fields or None


def compare_allele(sr_fields, truth_fields):
    """Per-field cascade for one SR allele vs one truth allele. Returns a list of 4
    values in {True, False, None} for field levels 1..4 (see module docstring)."""
    result = []
    mismatched = False
    unassessable = False
    for lvl in range(1, 5):
        if mismatched:
            result.append(False)
            continue
        if unassessable:
            result.append(None)
            continue
        if lvl > len(sr_fields) or lvl > len(truth_fields):
            result.append(None)
            unassessable = True
            continue
        if sr_fields[lvl - 1] == truth_fields[lvl - 1]:
            result.append(True)
        else:
            result.append(False)
            mismatched = True
    return result


def pairing_score(res_a, res_b):
    """Lexicographic score for choosing the best sr<->truth alignment: maximize the
    count of True matches at field 1 first, then field 2, then 3, then 4."""
    return tuple(
        (1 if res_a[i] is True else 0) + (1 if res_b[i] is True else 0)
        for i in range(4)
    )


def compare_genotype(sr_pair, truth_pair):
    """sr_pair, truth_pair: (allele1_raw, allele2_raw). Returns a list of 2 per-allele
    cascade results (each a list of 4 True/False/None), using whichever of the 2 possible
    alignments scores best, or None if either genotype isn't callable at all."""
    sr1, sr2 = parse_fields(sr_pair[0]), parse_fields(sr_pair[1])
    t1, t2 = parse_fields(truth_pair[0]), parse_fields(truth_pair[1])
    if sr1 is None or sr2 is None or t1 is None or t2 is None:
        return None
    pairing_a = (compare_allele(sr1, t1), compare_allele(sr2, t2))
    pairing_b = (compare_allele(sr1, t2), compare_allele(sr2, t1))
    best = max([pairing_a, pairing_b], key=lambda p: pairing_score(p[0], p[1]))
    return list(best)


def load_rows(cohort, outroot):
    anc = dict(zip(cohort["person_id"].astype(str), cohort["ancestry"].astype(str)))
    frames, missing = [], []
    for pid in cohort["person_id"].astype(str):
        csv = os.path.join(outroot, pid, "comparison_log.csv")
        if not os.path.exists(csv):
            missing.append(pid)
            continue
        df = pd.read_csv(csv, dtype=str)
        df = df[df["run_label"] == "experiment_d"]
        if df.empty:
            missing.append(pid)
            continue
        frames.append(df)
    if not frames:
        print("FATAL: no experiment_d rows found yet.", file=sys.stderr)
        sys.exit(1)
    allrows = pd.concat(frames, ignore_index=True)
    allrows = allrows.sort_values("timestamp").drop_duplicates(["person_id", "gene"], keep="last")
    allrows["ancestry"] = allrows["person_id"].astype(str).map(anc)
    if missing:
        print(f"NOTE: {len(missing)} cohort people have no experiment_d result yet: {missing}",
              file=sys.stderr)
    return allrows


def analyze(allrows):
    """Long-format counts: one row per (gene, ancestry, method, field_level) with
    n_true / n_false / n_na, counted per-allele (2 alleles per callable person-gene)."""
    counts = {}

    def bump(gene, grp, method, lvl, val):
        key = (gene, grp, method, lvl)
        rec = counts.setdefault(key, {"n_true": 0, "n_false": 0, "n_na": 0})
        if val is True:
            rec["n_true"] += 1
        elif val is False:
            rec["n_false"] += 1
        else:
            rec["n_na"] += 1

    for _, r in allrows.iterrows():
        gene, grp = r["gene"], r["ancestry"]
        if gene not in GENES or grp not in GROUPS:
            continue
        truth_pair = (r.get("specimmune_1"), r.get("specimmune_2"))
        for method, cols in (("aou", ("aou_1", "aou_2")), ("sh", ("spechla_1", "spechla_2"))):
            sr_pair = (r.get(cols[0]), r.get(cols[1]))
            res = compare_genotype(sr_pair, truth_pair)
            if res is None:
                continue
            for allele_res in res:
                for lvl in range(1, 5):
                    bump(gene, grp, method, lvl, allele_res[lvl - 1])

    recs = [
        {"gene": g, "ancestry": grp, "method": m, "field_level": lvl, **v}
        for (g, grp, m, lvl), v in counts.items()
    ]
    return pd.DataFrame(recs)


def rate(row):
    d = row["n_true"] + row["n_false"]
    if d == 0:
        return float("nan"), 0
    return 100 * row["n_true"] / d, d


def pooled_table(df):
    """Gene x method x field-level, summed across ancestry."""
    g = df.groupby(["gene", "method", "field_level"], as_index=False)[
        ["n_true", "n_false", "n_na"]
    ].sum()
    lines = ["## Pooled cascade -- match rate through each field, all ancestries combined\n",
             "Match rate = n_true / (n_true + n_false) among *assessable* comparisons; "
             "N/A = comparisons where one caller simply doesn't report that many fields "
             "(e.g. AoU-native never reports Field 4).\n"]
    for method in METHODS:
        lines.append(f"### {METHOD_LABEL[method]} vs SpecImmune-LR truth\n")
        lines.append("| Gene | " + " | ".join(FIELD_LABEL[l] for l in range(1, 5)) + " |")
        lines.append("|---|" + "---|" * 4)
        for gene in GENES:
            cells = [gene]
            for lvl in range(1, 5):
                row = g[(g.gene == gene) & (g.method == method) & (g.field_level == lvl)]
                if row.empty:
                    cells.append("— (no data)")
                    continue
                r0 = row.iloc[0]
                pct, d = rate(r0)
                na_pct = 100 * r0.n_na / (d + r0.n_na) if (d + r0.n_na) else 0
                if d == 0:
                    cells.append("N/A (100%)" if r0.n_na else "— (no data)")
                else:
                    cells.append(f"{pct:.0f}% (N/A {na_pct:.0f}%)")
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
    return "\n".join(lines), g


def field2_ancestry_table(df):
    lines = ["## Field 2 (protein-level) match rate by gene x ancestry\n",
             "The headline correctness metric -- does the caller get the actual HLA protein "
             "right, not just the allele family.\n"]
    for method in METHODS:
        lines.append(f"### {METHOD_LABEL[method]}\n")
        lines.append("| Gene | " + " | ".join(GROUPS) + " |")
        lines.append("|---|" + "---|" * len(GROUPS))
        for gene in GENES:
            cells = [gene]
            for grp in GROUPS:
                row = df[(df.gene == gene) & (df.ancestry == grp) & (df.method == method)
                         & (df.field_level == 2)]
                if row.empty:
                    cells.append("— (no data)")
                    continue
                r0 = row.iloc[0]
                pct, d = rate(r0)
                if d == 0:
                    cells.append("N/A (100%)" if r0.n_na else "— (no data)")
                else:
                    cells.append(f"{pct:.0f}% (n={d})")
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
    return "\n".join(lines)


def plot_cascade(pooled, path):
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    colors = {"aou": "#4C72B0", "sh": "#DD8452"}
    for ax, gene in zip(axes.flat, GENES):
        gene_total = pooled.loc[
            pooled.gene == gene, ["n_true", "n_false", "n_na"]
        ].sum().sum()
        if gene_total == 0:
            ax.text(0.5, 0.5, "no data", ha="center", va="center",
                     transform=ax.transAxes, fontsize=10, color="#999")
            ax.set_title(gene, fontsize=11)
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        x = range(4)
        width = 0.35
        for i, method in enumerate(METHODS):
            heights, na_labels = [], []
            for lvl in range(1, 5):
                row = pooled[(pooled.gene == gene) & (pooled.method == method)
                             & (pooled.field_level == lvl)]
                if row.empty:
                    heights.append(0)
                    na_labels.append("")
                    continue
                pct, d = rate(row.iloc[0])
                r0 = row.iloc[0]
                na_pct = 100 * r0.n_na / (d + r0.n_na) if (d + r0.n_na) else 0
                heights.append(0 if pct != pct else pct)
                na_labels.append(f"N/A {na_pct:.0f}%" if na_pct >= 1 else "")
            offs = [xi + (i - 0.5) * width for xi in x]
            bars = ax.bar(offs, heights, width=width, color=colors[method],
                           label=METHOD_LABEL[method] if gene == GENES[0] else None)
            # emphasize Field 2 (index 1) with a bold outline -- the field that matters most
            bars[1].set_edgecolor("black")
            bars[1].set_linewidth(1.8)
            # stagger AoU/SpecHLA N/A labels vertically so they never collide at low bar heights
            for xoff, h, lab in zip(offs, heights, na_labels):
                if lab:
                    ax.text(xoff, 6 + i * 9, lab, ha="center", va="bottom",
                             fontsize=6, color="#555", rotation=90)
        ax.set_title(gene, fontsize=11)
        ax.set_xticks(list(x))
        ax.set_xticklabels(["F1", "F2", "F3", "F4"], fontsize=8)
        ax.set_ylim(0, 108)
        ax.set_yticks(range(0, 101, 20))
        ax.axhline(100, color="#ddd", linewidth=0.6)
    axes.flat[0].set_ylabel("Match rate among assessable (%)")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 0.99))
    fig.suptitle("Experiment D -- field cascade vs SpecImmune-LR truth\n"
                  "F2 (protein) outlined bold -- the field that actually matters; F3/F4 shown lighter",
                  fontsize=10, y=0.935)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def heatmap_field2(df, method, path):
    grid = []
    for gene in GENES:
        row = []
        for grp in GROUPS:
            m = df[(df.gene == gene) & (df.ancestry == grp) & (df.method == method)
                   & (df.field_level == 2)]
            if m.empty:
                row.append(float("nan"))
                continue
            pct, d = rate(m.iloc[0])
            row.append(pct if d else float("nan"))
        grid.append(row)
    fig, ax = plt.subplots(figsize=(1.1 * len(GROUPS) + 2, 0.6 * len(GENES) + 2))
    im = ax.imshow(grid, cmap="RdYlGn", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(GROUPS)), GROUPS)
    ax.set_yticks(range(len(GENES)), GENES)
    for i in range(len(GENES)):
        for j in range(len(GROUPS)):
            v = grid[i][j]
            if v == v:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=8)
    ax.set_title(f"{METHOD_LABEL[method]} -- Field 2 (protein) match rate (%)", fontsize=10)
    fig.colorbar(im, ax=ax, label="%")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cohort", help="cohort.tsv from build_experiment_d_cohort.py")
    ap.add_argument("--outroot", default=os.path.expanduser("~/pipeline_outputs"))
    ap.add_argument("--analysis-dir", default=None)
    args = ap.parse_args()

    cohort = pd.read_csv(args.cohort, sep="\t", dtype=str)
    allrows = load_rows(cohort, args.outroot)
    df = analyze(allrows)

    expdir = os.path.join(args.outroot, "experiment_d")
    adir = args.analysis_dir or os.path.join(expdir, "analysis")
    os.makedirs(adir, exist_ok=True)

    pooled_md, pooled_df = pooled_table(df)
    anc_md = field2_ancestry_table(df)

    md = "# Experiment D -- field-level cascade (AoU-native / SpecHLA vs SpecImmune-LR truth)\n\n" \
         "Each field level is cumulative (matching through Field N implies Field 1..N-1 also " \
         "matched). N/A means the comparison isn't assessable at that depth (one caller didn't " \
         "report that many fields) -- not counted as a disagreement.\n\n" + pooled_md + "\n" + anc_md
    md_path = os.path.join(adir, "experiment_d_field_cascade.md")
    with open(md_path, "w") as f:
        f.write(md + "\n")

    plot_cascade(pooled_df, os.path.join(adir, "experiment_d_field_cascade.png"))
    heatmap_field2(df, "aou", os.path.join(adir, "experiment_d_field2_heatmap_aou.png"))
    heatmap_field2(df, "sh", os.path.join(adir, "experiment_d_field2_heatmap_spechla.png"))

    df.to_csv(os.path.join(adir, "experiment_d_field_cascade_counts.csv"), index=False)

    print(md)
    print(f"\n(written to {md_path} + 3 PNGs + counts CSV in {adir} -- aggregate only, "
          f"but keep on the VM per the egress caveat; paste tables/describe PNGs back)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
