#!/usr/bin/env python3
"""Experiment D aggregate analysis -- ancestry-stratified 3-way concordance.

Reads the per-person comparison_log.csv files that run_experiment_d.sh accumulated (one under
each ~/pipeline_outputs/<pid>/, appended by compare_hla_results.py), keeps only the
run_label=="experiment_d" rows, joins the ancestry label from cohort.tsv, and computes
per-locus x per-ancestry concordance.

NOT a naive vote (DECISIONS.md: AoU-native and SpecHLA share short-read data, so they are not
independent corroboration -- a majority-vote "outlier" framing is wrong). Instead it reports the
pairwise structure explicitly and, for the two carried-over open findings, the specific patterns:
  - DQA1: how often is AoU-native the lone outlier (SpecHLA == SpecImmune, AoU differs)? Is that
    rate ancestry-correlated or flat? (smoke tests: AoU-outlier at DQA1 in 3/3 people.)
  - DPA1: how often does SpecImmune diverge from an agreeing AoU+SpecHLA, and is that ancestry-
    specific? Weighted by SpecImmune's own confidence (tie flag + identity), because "over-
    resolution" and a real error look different in the confidence signal (DECISIONS.md).

All outputs are aggregate rates/counts -- no genotypes -- but they are derived from real calls,
so keep them on the VM (egress caveat, DECISIONS.md) and paste the tables/summary back rather
than committing them. Writes markdown + two PNG heatmaps under <expdir>/analysis/.

Usage (via `pixi run -e spechla` for pandas+matplotlib):
  python3 analyze_experiment_d.py ~/pipeline_outputs/experiment_d/cohort.tsv
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
NULL = {"", "NA", "-", "nan", "None", ".", "na"}


def norm2(allele):
    """Normalize an allele string to unordered-comparable 2-field, e.g. 'A*01:01:01:01'->'01:01',
    '01:01'->'01:01'. Returns None if uncallable or below 2-field."""
    if allele is None:
        return None
    a = str(allele).strip()
    if a in NULL:
        return None
    if "*" in a:
        a = a.split("*", 1)[1]
    fields = a.split(":")
    if len(fields) < 2:
        return None
    return f"{fields[0]}:{fields[1]}"


def geno(pair):
    """(a1, a2) -> frozenset-like sorted tuple of 2-field alleles, or None if either uncallable."""
    n1, n2 = norm2(pair[0]), norm2(pair[1])
    if n1 is None or n2 is None:
        return None
    return tuple(sorted((n1, n2)))


def truthy(x):
    return str(x).strip().lower() in {"true", "1", "yes"}


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
        print("FATAL: no experiment_d rows found in any per-person comparison_log.csv yet -- "
              "has run_experiment_d.sh produced any completed people?", file=sys.stderr)
        sys.exit(1)
    allrows = pd.concat(frames, ignore_index=True)
    # A person may have been re-compared on a resume -- keep the latest (person, gene) by timestamp.
    allrows = allrows.sort_values("timestamp").drop_duplicates(["person_id", "gene"], keep="last")
    allrows["ancestry"] = allrows["person_id"].astype(str).map(anc)
    if missing:
        print(f"NOTE: {len(missing)} cohort people have no experiment_d result yet "
              f"(not run / failed): {missing}", file=sys.stderr)
    n_people = allrows["person_id"].nunique()
    print(f"Loaded {len(allrows)} gene-rows across {n_people} people.", file=sys.stderr)
    return allrows


def analyze(allrows):
    """Returns a per (gene, ancestry) records list with concordance counts/rates."""
    recs = []
    for gene in GENES:
        g = allrows[allrows["gene"] == gene]
        for grp in GROUPS:
            sub = g[g["ancestry"] == grp]
            n = agree3 = sh_aou = si_aou = si_sh = aou_out = si_out = si_out_conf = callable_n = 0
            for _, r in sub.iterrows():
                n += 1
                aou = geno((r.get("aou_1"), r.get("aou_2")))
                sh = geno((r.get("spechla_1"), r.get("spechla_2")))
                si = geno((r.get("specimmune_1"), r.get("specimmune_2")))
                if aou is None or sh is None or si is None:
                    continue
                callable_n += 1
                if aou == sh:
                    sh_aou += 1
                if si == aou:
                    si_aou += 1
                if si == sh:
                    si_sh += 1
                if aou == sh == si:
                    agree3 += 1
                if sh == si and aou != sh:            # AoU-native lone outlier (the DQA1 pattern)
                    aou_out += 1
                if aou == sh and si != sh:            # SpecImmune diverges from agreeing SR pair (DPA1)
                    si_out += 1
                    conf = (not truthy(r.get("si_tied_1"))) and (not truthy(r.get("si_tied_2")))
                    try:
                        conf = conf and float(r.get("si_identity_1")) >= 0.99 \
                                     and float(r.get("si_identity_2")) >= 0.99
                    except (TypeError, ValueError):
                        conf = False
                    if conf:
                        si_out_conf += 1
            recs.append(dict(gene=gene, ancestry=grp, n=n, callable_n=callable_n,
                             agree3=agree3, sh_aou=sh_aou, si_aou=si_aou, si_sh=si_sh,
                             aou_out=aou_out, si_out=si_out, si_out_conf=si_out_conf))
    return pd.DataFrame(recs)


def rate(num, den):
    return f"{num}/{den} ({100*num/den:.0f}%)" if den else "—"


def md_table(df, value_fn, title):
    lines = [f"### {title}", "", "| Gene | " + " | ".join(GROUPS) + " |",
             "|---|" + "---|" * len(GROUPS)]
    for gene in GENES:
        cells = []
        for grp in GROUPS:
            row = df[(df.gene == gene) & (df.ancestry == grp)]
            cells.append(value_fn(row.iloc[0]) if not row.empty else "—")
        lines.append(f"| {gene} | " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def heatmap(df, num_col, den_col, title, path):
    grid = []
    for gene in GENES:
        row = []
        for grp in GROUPS:
            m = df[(df.gene == gene) & (df.ancestry == grp)]
            if m.empty or m.iloc[0][den_col] == 0:
                row.append(float("nan"))
            else:
                row.append(100 * m.iloc[0][num_col] / m.iloc[0][den_col])
        grid.append(row)
    fig, ax = plt.subplots(figsize=(1.1 * len(GROUPS) + 2, 0.6 * len(GENES) + 2))
    im = ax.imshow(grid, cmap="RdYlGn" if "agree" in num_col else "Reds",
                   vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(GROUPS)), GROUPS)
    ax.set_yticks(range(len(GENES)), GENES)
    for i in range(len(GENES)):
        for j in range(len(GROUPS)):
            v = grid[i][j]
            if v == v:  # not NaN
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=8)
    ax.set_title(title, fontsize=10)
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

    parts = [
        "# Experiment D — ancestry-stratified 3-way concordance",
        "",
        f"People with results: {allrows['person_id'].nunique()}. "
        "Cells are among people where all 3 methods gave a callable 2-field genotype at that "
        "locus. Not a vote — see the DQA1/DPA1 sections for the specific outlier patterns.",
        "",
        md_table(df, lambda r: rate(r.agree3, r.callable_n), "Full 3-way agreement (AoU = SpecHLA = SpecImmune)"),
        md_table(df, lambda r: rate(r.sh_aou, r.callable_n), "SpecHLA-SR vs AoU-native (both short-read)"),
        md_table(df, lambda r: rate(r.si_aou, r.callable_n), "SpecImmune-LR vs AoU-native (cross-technology)"),
        md_table(df, lambda r: rate(r.si_sh, r.callable_n), "SpecImmune-LR vs SpecHLA-SR (cross-technology)"),
        "## Carried-over question 1 — is AoU-native's DQA1 discordance ancestry-correlated?",
        "Rate at which AoU-native is the LONE outlier (SpecHLA == SpecImmune, AoU differs). "
        "Look at the DQA1 row: flat across ancestries ⇒ method/database bias; a gradient ⇒ ancestry-linked.",
        "",
        md_table(df, lambda r: rate(r.aou_out, r.callable_n), "AoU-native lone-outlier rate"),
        "## Carried-over question 2 — is SpecImmune's DPA1 divergence ancestry-specific?",
        "Rate at which SpecImmune diverges from an agreeing AoU+SpecHLA pair (raw), and the subset "
        "that is confident (SpecImmune not tied, identity ≥ 0.99 on both haplotypes). Look at DPA1.",
        "",
        md_table(df, lambda r: rate(r.si_out, r.callable_n), "SpecImmune divergence-from-SR-pair rate (raw)"),
        md_table(df, lambda r: rate(r.si_out_conf, r.callable_n), "…of which CONFIDENT (not tied, identity ≥ 0.99)"),
    ]
    md = "\n".join(parts)
    md_path = os.path.join(adir, "experiment_d_analysis.md")
    with open(md_path, "w") as f:
        f.write(md + "\n")

    heatmap(df, "agree3", "callable_n", "Experiment D — 3-way agreement (%)",
            os.path.join(adir, "experiment_d_3way_agreement.png"))
    heatmap(df, "aou_out", "callable_n", "AoU-native lone-outlier rate (%) — watch DQA1",
            os.path.join(adir, "experiment_d_aou_outlier.png"))

    df.to_csv(os.path.join(adir, "experiment_d_concordance_counts.csv"), index=False)
    print(md)
    print(f"\n(written to {md_path} + 2 PNG heatmaps + counts CSV in {adir} — aggregate only, "
          f"but keep on the VM per the egress caveat; paste the tables back)", file=sys.stderr)


if __name__ == "__main__":
    main()
