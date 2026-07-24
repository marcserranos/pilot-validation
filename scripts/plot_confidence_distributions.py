#!/usr/bin/env python3
"""Uncertainty distribution curves for SpecImmune-LR and Immuannot, plus a threshold-sweep
table -- Marc, 2026-07-24 weekly objective: before picking a fixed confidence cutoff (or a
downweighting scheme), "see the uncertainty distribution curve and then be able to cap the
tail, depending on its size."

Two signals, NOT on a shared scale (same discipline as plot_confidence_comparison() in
diagnose_immuannot_pilot.py -- never combine them into one number):
  - SpecImmune-LR: si_identity_1/2 (float, 0-1; higher = more confident) + si_tied_1/2 (bool).
    Per genotype (person, gene), this script uses the WORST of the two haplotypes'
    identities (min(id1, id2)) as the genotype-level confidence scalar -- a genotype is only
    as confident as its less-confident haplotype. The project's existing "confident" bar
    (analyze_experiment_d.py, reused by si_confidence_label()) is: not tied AND both
    identities >= 0.99. That is a THRESHOLD already in use, at t=0.99 -- this script shows
    what fraction of calls per gene that threshold (and looser alternatives) would exclude.
  - Immuannot: template_distance (edit distance to the nearest IMGT allele; 0 = perfect
    match, higher = less confident), read per-haplotype from each person's hap{1,2}.gtf.gz.
    Per genotype, this script uses the WORST of the two haplotypes' template_distance
    (max(td1, td2)) -- same "as confident as the worse haplotype" logic as SpecImmune above.
    NOTE: template_distance is a TEMPLATE-SELECTION confidence signal (how well the contig
    matches the chosen reference template used for annotation) -- it is NOT the same thing
    as the "new" token (a separate, orthogonal marker for "couldn't resolve nomenclature this
    deep", parsed out of the allele string itself, see analyze_experiment_d_field_cascade_
    immuannot.py's parse_fields()). A call can have template_distance=0 and still carry a
    "new" tag at a deep field, or vice versa -- don't conflate the two when deciding what to
    exclude. This script also reports template_warning (structural-annotation warning,
    e.g. incomplete_CDS/inframe_stop -- see reference/README_Immuannot.md) as a separate,
    binary-exclude flag: a real warning should probably exclude a call regardless of distance.

No new VM run needed -- reads the same completed-pilot files as diagnose_immuannot_pilot.py
(comparison_log.csv, hap{1,2}.gtf.gz). Aggregate-only output (rates/counts, no genotypes) --
keep on the VM per the standing egress rule; bring back the markdown + 2 PNGs.

Usage (via `pixi run -e spechla` for pandas+matplotlib):
  python3 scripts/plot_confidence_distributions.py ~/pipeline_outputs/experiment_d/cohort.tsv
"""
import argparse
import gzip
import os
import re
import sys

import pandas as pd

import matplotlib
matplotlib.use("Agg")  # no display server on the VM
import matplotlib.pyplot as plt

GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
DRB1_HIGHLIGHT = "#C44E52"
DEFAULT_COLOR = "#4C72B0"

# GTF attributes mix quoted and unquoted values -- see diagnose_immuannot_pilot.py's own comment
# for why the quoted-only-first version of this regex silently dropped unquoted attributes.
ATTR_RE = re.compile(r'(\w+)\s+(?:"([^"]*)"|([^";]+))')

SI_IDENTITY_THRESHOLDS = [0.99, 0.95, 0.90, 0.80]
IMMUANNOT_TD_THRESHOLDS = [0, 1, 2, 3, 5]


def parse_gtf_rich(gtf_gz_path):
    """Copied (not imported) from diagnose_immuannot_pilot.py -- same reasoning: keep this
    script runnable in a minimal env without pulling in the rest of that file's dependencies."""
    per_gene, gid_to_name = {}, {}
    with gzip.open(gtf_gz_path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = {k: (q if q else b.strip()) for k, q, b in ATTR_RE.findall(fields[8])}
            gid, gname = attrs.get("gene_id"), attrs.get("gene_name")
            if gid and gname:
                gid_to_name[gid] = gname
            key = gname or (gid_to_name.get(gid) if gid else None)
            if not key:
                continue
            slot = per_gene.setdefault(key, {})
            for k, v in attrs.items():
                slot.setdefault(k, v)
    return per_gene


def load_specimmune_genotype_confidence(cohort, outroot):
    """{(person_id, gene_bare): (worst_identity_float_or_None, tied_bool)} -- one row per
    (person, gene) with a resolved SpecImmune call, from comparison_log.csv."""
    out = {}
    for pid in cohort["person_id"].astype(str):
        log_path = os.path.join(outroot, pid, "comparison_log.csv")
        if not os.path.exists(log_path):
            continue
        df = pd.read_csv(log_path, dtype=str)
        df = df[df["run_label"] == "experiment_d"]
        if df.empty:
            continue
        df = df.sort_values("timestamp").drop_duplicates(["person_id", "gene"], keep="last")
        for _, r in df.iterrows():
            gene = str(r["gene"])

            def truthy(x):
                return str(x).strip().lower() in {"true", "1", "yes"}
            tied = truthy(r.get("si_tied_1")) or truthy(r.get("si_tied_2"))
            try:
                id1, id2 = float(r.get("si_identity_1")), float(r.get("si_identity_2"))
                worst = min(id1, id2)
            except (TypeError, ValueError):
                worst = None
            out[(str(r["person_id"]), gene)] = (worst, tied)
    return out


def load_immuannot_genotype_confidence(cohort, outroot):
    """{(person_id, gene_bare): (worst_template_distance_float_or_None, any_real_warning_bool)}
    -- one row per (person, gene) where at least one haplotype's GTF has a gene-row entry."""
    out = {}
    for pid in cohort["person_id"].astype(str):
        per_hap = {}
        for hap in ("hap1", "hap2"):
            gz = os.path.join(outroot, pid, "immuannot_output", f"{hap}.gtf.gz")
            if not os.path.exists(gz):
                continue
            try:
                per_hap[hap] = parse_gtf_rich(gz)
            except OSError:
                continue
        genes_seen = set()
        for attrs in per_hap.values():
            for gname in attrs:
                genes_seen.add(gname.replace("HLA-", ""))
        for gene in genes_seen:
            tds, warnings = [], []
            for hap, attrs in per_hap.items():
                a = attrs.get(f"HLA-{gene}", attrs.get(gene))
                if a is None:
                    continue
                try:
                    tds.append(float(a.get("template_distance")))
                except (TypeError, ValueError):
                    pass
                w = a.get("template_warning")
                if w and str(w).strip().upper() not in {"NA", ""}:
                    warnings.append(w)
            worst_td = max(tds) if tds else None
            out[(str(pid), gene)] = (worst_td, bool(warnings))
    return out


def plot_strip(values_by_gene, threshold_lines, ylabel, title, path, higher_is_confident):
    """One point per (person, gene) resolved call, jittered on x within each gene's column --
    handles this cohort's small per-gene n (~40-57) far better than a binned histogram would."""
    import random
    rng = random.Random(0)  # deterministic jitter -- a rerun on the same data looks the same
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for i, gene in enumerate(GENES):
        vals = [v for v in values_by_gene.get(gene, []) if v is not None]
        if not vals:
            continue
        xs = [i + rng.uniform(-0.18, 0.18) for _ in vals]
        color = DRB1_HIGHLIGHT if gene == "DRB1" else DEFAULT_COLOR
        ax.scatter(xs, vals, s=14, alpha=0.55, color=color, edgecolors="none")
        med = sorted(vals)[len(vals) // 2]
        ax.hlines(med, i - 0.25, i + 0.25, color="black", linewidth=1.5, zorder=5)
    for t in threshold_lines:
        ax.axhline(t, color="#888", linestyle="--", linewidth=0.8, zorder=1)
        ax.text(len(GENES) - 0.4, t, f"{t:g}", fontsize=7, color="#666", va="bottom")
    ax.set_xticks(range(len(GENES)))
    ax.set_xticklabels(GENES)
    ax.set_ylabel(ylabel)
    note = "higher = more confident" if higher_is_confident else "lower = more confident"
    ax.set_title(f"{title}\n({note}; one dot per resolved (person, gene) call, "
                 "black tick = median)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def threshold_sweep_table(values_by_gene, thresholds, exclude_if, header):
    """Markdown table: % of resolved calls per gene EXCLUDED at each candidate threshold --
    the 'how big is the tail if I cut here' view Marc asked for."""
    lines = [header, "", "| Gene | n resolved | " +
             " | ".join(f"excl. @ {t:g}" for t in thresholds) + " |",
             "|---|---|" + "---|" * len(thresholds)]
    for gene in GENES:
        vals = [v for v in values_by_gene.get(gene, []) if v is not None]
        n = len(vals)
        if n == 0:
            lines.append(f"| {gene} | 0 | " + " | ".join("—" for _ in thresholds) + " |")
            continue
        cells = []
        for t in thresholds:
            n_excl = sum(1 for v in vals if exclude_if(v, t))
            cells.append(f"{100*n_excl/n:.0f}%")
        lines.append(f"| {gene} | {n} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cohort", help="cohort.tsv from build_experiment_d_cohort.py")
    ap.add_argument("--outroot", default=os.path.expanduser("~/pipeline_outputs"))
    ap.add_argument("--analysis-dir", default=None)
    args = ap.parse_args()

    cohort = pd.read_csv(args.cohort, sep="\t", dtype=str)
    si_conf = load_specimmune_genotype_confidence(cohort, args.outroot)
    imm_conf = load_immuannot_genotype_confidence(cohort, args.outroot)

    si_by_gene = {g: [] for g in GENES}
    si_tied_count = {g: 0 for g in GENES}
    for (pid, gene), (worst_id, tied) in si_conf.items():
        if gene not in GENES:
            continue
        if tied:
            si_tied_count[gene] += 1
            continue  # tied calls have no meaningful identity scalar to plot -- counted separately
        si_by_gene[gene].append(worst_id)

    imm_by_gene = {g: [] for g in GENES}
    imm_warning_count = {g: 0 for g in GENES}
    for (pid, gene), (worst_td, has_warning) in imm_conf.items():
        if gene not in GENES:
            continue
        if has_warning:
            imm_warning_count[gene] += 1
        imm_by_gene[gene].append(worst_td)  # kept even with a warning -- warning is reported
                                             # separately, distance distribution stays complete

    expdir = os.path.join(args.outroot, "experiment_d")
    adir = args.analysis_dir or os.path.join(expdir, "analysis_confidence")
    os.makedirs(adir, exist_ok=True)

    plot_strip(si_by_gene, [0.99, 0.95],
               "Worst-haplotype si_identity (min of the 2 haplotypes)",
               "SpecImmune-LR genotype confidence by gene (tied calls excluded, counted separately)",
               os.path.join(adir, "specimmune_identity_distribution.png"),
               higher_is_confident=True)
    plot_strip(imm_by_gene, [0, 1, 2],
               "Worst-haplotype template_distance (max of the 2 haplotypes)",
               "Immuannot genotype confidence by gene",
               os.path.join(adir, "immuannot_template_distance_distribution.png"),
               higher_is_confident=False)

    si_table = threshold_sweep_table(
        si_by_gene, SI_IDENTITY_THRESHOLDS, lambda v, t: v < t,
        "## SpecImmune-LR: % of (non-tied) resolved calls EXCLUDED if requiring identity >= threshold\n\n"
        "Tied calls are already excluded from this table's denominator (n resolved) -- they are "
        "unconditionally low-confidence regardless of identity value; see the tied-count line below.")
    imm_table = threshold_sweep_table(
        imm_by_gene, IMMUANNOT_TD_THRESHOLDS, lambda v, t: v > t,
        "## Immuannot: % of resolved calls EXCLUDED if requiring template_distance <= threshold\n\n"
        "template_warning is reported separately below (a real structural warning, e.g. "
        "incomplete_CDS/inframe_stop, is a separate reason to exclude regardless of distance).")

    tied_lines = ["", "### SpecImmune-LR tied-call counts (excluded from the table above, by construction)",
                  "", "| Gene | n tied (unconditionally low-confidence) |", "|---|---|"]
    for g in GENES:
        tied_lines.append(f"| {g} | {si_tied_count[g]} |")

    warn_lines = ["", "### Immuannot template_warning counts (a real annotation-quality flag, independent of distance)",
                  "", "| Gene | n with a real template_warning |", "|---|---|"]
    for g in GENES:
        warn_lines.append(f"| {g} | {imm_warning_count[g]} |")

    md = ("# Uncertainty distributions -- SpecImmune-LR vs Immuannot\n\n"
          "Not the same scale (identity % vs. edit-distance count) -- read side by side, "
          "don't combine into one number.\n\n"
          + si_table + "\n" + "\n".join(tied_lines) + "\n\n"
          + imm_table + "\n" + "\n".join(warn_lines) + "\n")
    md_path = os.path.join(adir, "confidence_distributions.md")
    with open(md_path, "w") as f:
        f.write(md)

    print(md)
    print(f"\n(written to {md_path} + 2 PNGs in {adir} -- aggregate only, but keep on the VM "
          f"per the egress caveat; paste table/describe PNGs back)", file=sys.stderr)


if __name__ == "__main__":
    main()
