#!/usr/bin/env python3
"""Full-cohort AoU-native vs. Immuannot comparison: baseline concordance, a confidence-threshold
sweep, and a confidence distribution plot -- the three things Marc asked for (2026-08-10) to show
supervisors how AoU-native's quality compares to a high-confidence non-native "ground truth" proxy.

## The threshold-sweep idea, made concrete

Immuannot's own per-call confidence signal is `template_distance` (edit distance to the nearest
IMGT reference allele; 0 = perfect match, higher = less confident -- see
scripts/plot_confidence_distributions.py for the original single-cohort version of this signal).
As we restrict the comparison to only Immuannot calls at or below a given distance threshold t,
those calls become progressively closer to a genuinely trustworthy "ground truth" proxy -- so any
remaining AoU-native/Immuannot mismatch is more likely a real AoU-native error, not Immuannot's
own uncertainty. The discrepancy rate (1 - Field 2 concordance) should fall and level off as t
shrinks toward 0: that floor is the estimate of AoU-native's real error rate against
high-confidence ground truth (same interpretive logic as the confidence-matched-truth work
already done at n=60, scripts/analyze_confidence_matched_truth.py -- this generalizes it from one
fixed threshold to a full sweep, at full production scale, per gene).

X-axis is drawn with the axis REVERSED (large/loose distance on the left, 0/strict on the right)
so "raising the certainty threshold" reads naturally left-to-right, matching how Marc described it.

Reads (all on the VM):
  ~/pipeline_outputs/immuannot_calls.tsv
  ~/pipeline_outputs/<person_id>/immuannot_output/hap{1,2}.gtf.gz   (template_distance/_warning)
  ~/mnt/aou-controlled/v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv  (AoU-native)

No comparison_log.csv exists at production scale (that was only built for the earlier n=60 pilot)
-- this script does the AoU-native join itself, directly against the production calls file.

Writes markdown + 3 PNGs under --out-dir. Aggregate-only (rates/counts/distributions derived from
real calls, never the calls themselves) -- keep on the VM per the standing egress discipline.

Usage (via `pixi run -e spechla`):
  python3 scripts/production_analysis/analyze_confidence_vs_aou_native.py
"""
import argparse
import gzip
import os
import random
import re
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
GENE_COLOR = dict(zip(GENES, [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2", "#937860", "#DA8BC3", "#8C8C8C",
]))
NULL = {"", "NA", "nan", "None", ".", "-"}
ATTR_RE = re.compile(r'(\w+)\s+(?:"([^"]*)"|([^";]+))')

DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_MOUNT = os.path.expanduser("~/mnt/aou-controlled")
AOU_TSV_REL = "v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv"


def human_count(n):
    if n >= 1000:
        return f"{n / 1000:.1f}K"
    if n >= 100:
        return "hundreds"
    return str(n)


def parse_fields(raw):
    if raw is None or (hasattr(pd, "isna") and pd.isna(raw)):
        return None
    s = str(raw).strip()
    if s in NULL:
        return None
    if "*" in s:
        s = s.split("*", 1)[1]
    fields = [f for f in s.split(":") if f != ""]
    for i, f in enumerate(fields):
        if f.strip().lower() == "new":
            fields = fields[:i]
            break
    return fields or None


def compare_allele(sr_fields, truth_fields):
    result, mismatched, unassessable = [], False, False
    for lvl in range(1, 3):  # Field 1-2 only, same headline metric as everywhere else here
        if mismatched:
            result.append(False); continue
        if unassessable:
            result.append(None); continue
        if lvl > len(sr_fields) or lvl > len(truth_fields):
            result.append(None); unassessable = True; continue
        if sr_fields[lvl - 1] == truth_fields[lvl - 1]:
            result.append(True)
        else:
            result.append(False); mismatched = True
    return result


def pairing_score(res_a, res_b):
    return tuple((1 if res_a[i] is True else 0) + (1 if res_b[i] is True else 0) for i in range(2))


def compare_genotype(pair_a, pair_b):
    """Best-pairing Field-2 concordance for one genotype (2 alleles, unordered). Returns
    (n_allele_matches_field2, n_allele_comparable_field2) or None if either side unparseable."""
    a1, a2 = parse_fields(pair_a[0]), parse_fields(pair_a[1])
    b1, b2 = parse_fields(pair_b[0]), parse_fields(pair_b[1])
    if a1 is None or a2 is None or b1 is None or b2 is None:
        return None
    ra = (compare_allele(a1, b1), compare_allele(a2, b2))
    rb = (compare_allele(a1, b2), compare_allele(a2, b1))
    best = max([ra, rb], key=lambda p: pairing_score(p[0], p[1]))
    n_match = n_comp = 0
    for allele_res in best:
        v = allele_res[1]  # Field 2
        if v is True:
            n_match += 1; n_comp += 1
        elif v is False:
            n_comp += 1
    return n_match, n_comp


def parse_gtf_rich(gtf_gz_path):
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


def load_immuannot_calls(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for c in ["person_id", "gene", "immuannot_1", "immuannot_2"]:
        if c not in df.columns:
            sys.exit(f"FATAL: {path} missing column '{c}'. Actual: {list(df.columns)}")
    df["gene_bare"] = df["gene"].str.replace("^HLA-", "", regex=True)
    all_genes_seen = sorted(df["gene_bare"].unique())
    df = df[df["gene_bare"].isin(GENES)]
    if df.empty:
        sys.exit(f"FATAL: 0 of the loaded rows in {path} match any of the 8 classical genes "
                 f"{GENES}. Genes actually present (first 20): {all_genes_seen[:20]}. This is "
                 f"the exact symptom of the 2026-08-10 merge_fragments() dedup bug -- run "
                 f"scripts/production_orchestrator/rebuild_immuannot_calls.py first.")
    return {(r["person_id"], r["gene_bare"]): (r["immuannot_1"], r["immuannot_2"])
            for _, r in df.iterrows()}


def load_aou_native(path, cohort_person_ids):
    """AoU-native's own TSV covers the WHOLE AoU cohort (535,658+ people) -- filtering to just the
    production cohort's person_ids, not just for speed/memory (4.28M rows -> ~12k), but because an
    unfiltered load makes the printed count actively misleading (looked like 4.28M AoU-native
    calls existed for a ~12k-person run -- found 2026-08-10 comparing against the real production
    output)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    if "research_id" not in df.columns:
        sys.exit(f"FATAL: {path} missing 'research_id'. Actual columns (first 10): {list(df.columns)[:10]}")
    df = df[df["research_id"].isin(cohort_person_ids)]
    out = {}
    for _, r in df.iterrows():
        pid = r["research_id"]
        for g in GENES:
            c1, c2 = f"{g}_1", f"{g}_2"
            if c1 not in df.columns or c2 not in df.columns:
                continue
            out[(pid, g)] = (r[c1], r[c2])
    return out


def load_confidence(person_ids, outroot):
    """{(person_id, gene): (worst_template_distance_or_None, any_warning_bool)}."""
    out = {}
    for pid in person_ids:
        per_hap = {}
        for hap in ("hap1", "hap2"):
            gz = os.path.join(outroot, str(pid), "immuannot_output", f"{hap}.gtf.gz")
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
        for gene in genes_seen & set(GENES):
            tds, warnings = [], []
            for attrs in per_hap.values():
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
            out[(str(pid), gene)] = (max(tds) if tds else None, bool(warnings))
    return out


def build_combined(imm_calls, aou_calls, confidence):
    keys = set(imm_calls) | set(aou_calls)
    combined = {}
    for key in keys:
        td, warn = confidence.get(key, (None, False))
        combined[key] = {
            "immuannot": imm_calls.get(key), "aou": aou_calls.get(key),
            "worst_td": td, "warning": warn,
        }
    return combined


def field2_concordance(combined, keys):
    n_match = n_comp = 0
    for key in keys:
        row = combined[key]
        if row["immuannot"] is None or row["aou"] is None:
            continue
        res = compare_genotype(row["immuannot"], row["aou"])
        if res is None:
            continue
        nm, nc = res
        n_match += nm; n_comp += nc
    return n_match, n_comp


def wilson_ci(n_match, n_total, z=1.96):
    if n_total == 0:
        return (float("nan"), float("nan"), float("nan"))
    phat = n_match / n_total
    denom = 1 + z * z / n_total
    center = (phat + z * z / (2 * n_total)) / denom
    margin = z * np.sqrt(phat * (1 - phat) / n_total + z * z / (4 * n_total ** 2)) / denom
    return phat, max(0, center - margin), min(1, center + margin)


def plot_baseline_forest(combined, out_path):
    """Per-locus Field-2 concordance, unfiltered, ordered class I -> class II (the project's
    standing gene order). Error bars = Wilson 95% CI. This is the full-production-scale version
    of the same headline metric this project has tracked at n=60 throughout."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ys = list(range(len(GENES)))
    pts, los, his, ns = [], [], [], []
    for g in GENES:
        keys = [k for k in combined if k[1] == g]
        nm, nc = field2_concordance(combined, keys)
        p, lo, hi = wilson_ci(nm, nc)
        pts.append(100 * p); los.append(100 * (p - lo)); his.append(100 * (hi - p)); ns.append(nc)
    ax.errorbar(pts, ys, xerr=[los, his], fmt="o", color="#2C5C8A", ecolor="#88AACC",
                elinewidth=2, capsize=4, markersize=7)
    for y, p, n in zip(ys, pts, ns):
        ax.text(p + 2.5, y, f"{p:.0f}% (n={human_count(n)} alleles)", va="center", fontsize=8)
    ax.set_yticks(ys)
    ax.set_yticklabels(GENES)
    ax.invert_yaxis()
    ax.set_xlim(0, 118)
    ax.set_xlabel("Field 2 (protein) concordance, AoU-native vs Immuannot (%)")
    ax.set_title("Per-locus concordance, full production cohort (unfiltered, 95% Wilson CI)")
    ax.axvline(100, color="#CCCCCC", linewidth=0.8, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def build_threshold_grid(pooled_td, max_threshold=5.0, n_points=24):
    """Most real template_distance values cluster near 0 with a long thin tail (confirmed
    2026-08-10 against real production data -- every gene's median sits near 0 while a few
    outliers reach 100-200). A grid spanning the full range (e.g. pooled 99th percentile, ~30)
    wastes most of its points on territory where the threshold barely excludes anything -- 14 of
    16 points were flat. Capping at max_threshold (default 5, Marc's suggestion) concentrates the
    grid where the data actually lives; bounded by the real observed max so a sparse/low-range
    dataset never gets an artificially wide, mostly-empty grid either."""
    if not pooled_td:
        return [0, 1, 2, 3, 5]
    cap = min(max_threshold, max(pooled_td))
    grid = sorted(set(round(x, 2) for x in np.linspace(0, max(cap, 0.1), n_points)))
    return grid


def plot_threshold_sweep(combined, out_path, max_threshold=5.0):
    linthresh = 0.5
    all_td = [row["worst_td"] for row in combined.values() if row["worst_td"] is not None]
    grid = build_threshold_grid(all_td, max_threshold=max_threshold)

    per_gene_curve = {g: {"x": [], "y": [], "n": []} for g in GENES}
    mean_curve = {"x": [], "y": []}
    pooled_n_curve = {"x": [], "n": []}
    for t in grid:
        gene_pcts = []
        pooled_n = 0
        for g in GENES:
            keys = [k for k, row in combined.items()
                    if k[1] == g and row["worst_td"] is not None
                    and row["worst_td"] <= t and not row["warning"]]
            nm, nc = field2_concordance(combined, keys)
            if nc == 0:
                continue
            disc = 100 * (1 - nm / nc)
            per_gene_curve[g]["x"].append(t)
            per_gene_curve[g]["y"].append(disc)
            per_gene_curve[g]["n"].append(nc)
            gene_pcts.append(disc)
            pooled_n += nc
        if gene_pcts:
            mean_curve["x"].append(t)
            mean_curve["y"].append(float(np.mean(gene_pcts)))
        pooled_n_curve["x"].append(t)
        pooled_n_curve["n"].append(pooled_n)

    fig, ax = plt.subplots(figsize=(11, 6))
    for g in GENES:
        c = per_gene_curve[g]
        if c["x"]:
            ax.plot(c["x"], c["y"], marker="o", markersize=3, linewidth=1.3,
                    color=GENE_COLOR[g], label=g, alpha=0.9)
    ax.plot(mean_curve["x"], mean_curve["y"], linestyle="--", linewidth=2.2, color="black",
            label="mean across genes")

    ax2 = ax.twinx()
    ax2.fill_between(pooled_n_curve["x"], pooled_n_curve["n"], color="#CCCCCC", alpha=0.25, zorder=0)
    # This is an ALLELE-level count (up to 2 per person x gene, one per haplotype), not a
    # person-gene pair count -- field2_concordance()/compare_genotype() score each haplotype
    # separately. Mislabeled as "(person, gene) pairs" before 2026-08-10; fixed after the real
    # production run showed ~172K here against only 97,724 actual (person, gene) pairs -- the
    # ~1.76x ratio is explained by most genotypes contributing 2 comparable alleles, not a bug.
    ax2.set_ylabel("n allele-level comparisons, pooled across genes", fontsize=9, color="#888")
    ax2.tick_params(axis="y", colors="#888")
    # Annotate N at 3 points spread by VALUE (loosest, linthresh, strictest), not by linear index
    # position -- picking indices from a linearly-spaced grid and placing them on a symlog axis
    # bunched two of the three annotations on top of each other (found 2026-08-10: "2.6K" appeared
    # twice, overlapping, near the loosest end).
    target_xs = [pooled_n_curve["x"][-1], linthresh, 0.0]
    xs_arr = pooled_n_curve["x"]
    seen_i = set()
    for target in target_xs:
        i = min(range(len(xs_arr)), key=lambda j: abs(xs_arr[j] - target))
        if i in seen_i:
            continue
        seen_i.add(i)
        x, n = xs_arr[i], pooled_n_curve["n"][i]
        ax2.text(x, n, human_count(n), fontsize=8, color="#666", ha="center", va="bottom")

    ax.invert_xaxis()  # 0/strict on the right -- "raising certainty" reads left-to-right
    # symlog (not plain log): 0 is a real, common, meaningful value here (most confident calls) --
    # plain log can't render it. symlog linearizes near 0 and logs beyond linthresh, which is
    # exactly what a "cluster at 0, thin tail beyond" distribution needs (Marc, 2026-08-10: "add a
    # log scale... or a much more restrictive window (perhaps 5 as the max)" -- did both).
    ax.set_xscale("symlog", linthresh=linthresh)
    # symlog's default tick locator/formatter renders "10^0"/"10^-1" scientific notation -- not
    # what a "max allowed distance" axis should show. Don't subsample the (linearly-spaced) grid
    # for ticks either -- linear-spaced values bunch up once placed on a log-spaced axis. Instead
    # generate ticks that are themselves evenly spaced in symlog space: a doubling sequence
    # from linthresh up to the cap, plus one linear-region point below linthresh.
    tick_vals = [0.0, linthresh / 2, linthresh]
    v = linthresh
    cap = grid[-1]
    while v < cap:
        v = round(v * 2, 2)
        tick_vals.append(min(v, cap))
    tick_vals = sorted(set(round(t, 2) for t in tick_vals))
    ax.set_xticks(tick_vals)
    ax.set_xticklabels([f"{t:g}" for t in tick_vals])
    ax.set_xlabel("Max allowed template_distance kept  (looser ←  → stricter / higher certainty)")
    ax.set_ylabel("Discrepancy vs AoU-native, Field 2 (%)")
    ax.set_title("Discrepancy vs AoU-native as the Immuannot confidence bar tightens, per gene\n"
                 "(should fall and level off toward the right -- that floor estimates AoU-native's\n"
                 "real error rate against high-confidence ground truth. Right edge often gets noisy\n"
                 "-- watch the shaded N curve: a thin tail there means read that end skeptically)")
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.02), fontsize=8,
               ncol=5, frameon=False)
    ax.spines["top"].set_visible(False)
    fig.tight_layout(rect=[0, 0.09, 1, 1])
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return grid


def plot_confidence_violin_scatter(combined, out_path, max_points_per_gene=2500):
    """Enhanced version of plot_confidence_distributions.py's plot_strip(), for the full
    production cohort: violin (clustering/shape) + jittered scatter (subsampled for render
    clarity -- mean/median are still computed from the FULL data) + explicit mean (diamond) and
    median (tick) markers."""
    rng = random.Random(0)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    violin_data, positions = [], []
    for i, g in enumerate(GENES):
        vals = [row["worst_td"] for k, row in combined.items()
                if k[1] == g and row["worst_td"] is not None]
        if not vals:
            continue
        violin_data.append(vals)
        positions.append(i)

    if violin_data:
        parts = ax.violinplot(violin_data, positions=positions, showmeans=False,
                               showmedians=False, showextrema=False, widths=0.7)
        for pc in parts["bodies"]:
            pc.set_facecolor("#CCCCCC")
            pc.set_alpha(0.5)

    for i, g in enumerate(GENES):
        vals = [row["worst_td"] for k, row in combined.items()
                if k[1] == g and row["worst_td"] is not None]
        if not vals:
            continue
        plot_vals = vals if len(vals) <= max_points_per_gene else rng.sample(vals, max_points_per_gene)
        xs = [i + rng.uniform(-0.15, 0.15) for _ in plot_vals]
        color = "#C44E52" if g == "DRB1" else "#4C72B0"
        ax.scatter(xs, plot_vals, s=6, alpha=0.25, color=color, edgecolors="none", zorder=2)
        med = float(np.median(vals))
        mean = float(np.mean(vals))
        ax.hlines(med, i - 0.3, i + 0.3, color="black", linewidth=2, zorder=5)
        ax.scatter([i], [mean], marker="D", s=45, color="gold", edgecolors="black",
                   linewidths=0.8, zorder=6)

    ax.scatter([], [], marker="D", s=45, color="gold", edgecolors="black", label="mean")
    ax.hlines([], [], [], color="black", linewidth=2, label="median")
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    ax.set_xticks(range(len(GENES)))
    ax.set_xticklabels(GENES)
    # symlog, not plain log: template_distance=0 is the single most common, most meaningful value
    # (perfect match) -- plain log excludes it. Real production data (2026-08-10) showed most
    # calls clustered at/near 0 with a long thin tail to 100-200, squashing the informative part
    # into a thin band at the bottom of a linear axis -- symlog fixes that (Marc's request).
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylim(bottom=0)  # template_distance is never negative -- symlog defaults to a
                            # symmetric range around 0 otherwise, showing a misleading negative tick
    ax.set_ylabel("Worst-haplotype template_distance (lower = more confident, symlog scale)")
    ax.set_title("Immuannot confidence distribution by gene, full production cohort\n"
                 f"(scatter subsampled to {max_points_per_gene}/gene for render clarity; "
                 "mean/median use all data)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", default=os.path.join(DEFAULT_OUTROOT, "immuannot_calls.tsv"))
    ap.add_argument("--aou-tsv", default=os.path.join(DEFAULT_MOUNT, AOU_TSV_REL))
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT,
                     help="Where <person_id>/immuannot_output/hap{1,2}.gtf.gz live")
    ap.add_argument("--out-dir", default=os.path.join(DEFAULT_OUTROOT, "production_analysis", "confidence"))
    ap.add_argument("--max-threshold", type=float, default=5.0,
                     help="Cap on the threshold-sweep x-axis (default 5, per Marc's request "
                          "2026-08-10 -- most template_distance values cluster near 0, a wider "
                          "range mostly plots flat lines). Combined with a symlog x-axis, not "
                          "either/or.")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading Immuannot calls...", file=sys.stderr)
    imm_calls = load_immuannot_calls(args.calls)
    print(f"  {len(imm_calls)} (person, gene) Immuannot calls", file=sys.stderr)

    person_ids = sorted(set(pid for pid, _ in imm_calls))

    print(f"Loading AoU-native calls (filtered to this cohort's {len(person_ids)} people)...", file=sys.stderr)
    aou_calls = load_aou_native(args.aou_tsv, set(person_ids))
    print(f"  {len(aou_calls)} (person, gene) AoU-native calls", file=sys.stderr)
    print(f"Loading Immuannot confidence signals for {len(person_ids)} people "
          f"(reads hap1/hap2 .gtf.gz -- this is the slow step, I/O not CPU bound)...", file=sys.stderr)
    confidence = load_confidence(person_ids, args.outroot)
    print(f"  confidence signal found for {len(confidence)} (person, gene) pairs", file=sys.stderr)

    combined = build_combined(imm_calls, aou_calls, confidence)

    forest_path = os.path.join(args.out_dir, "baseline_concordance_by_gene.png")
    sweep_path = os.path.join(args.out_dir, "confidence_threshold_sweep.png")
    scatter_path = os.path.join(args.out_dir, "confidence_distribution.png")

    plot_baseline_forest(combined, forest_path)
    grid = plot_threshold_sweep(combined, sweep_path, max_threshold=args.max_threshold)
    plot_confidence_violin_scatter(combined, scatter_path)

    md = [
        "# AoU-native vs Immuannot -- full production cohort\n",
        f"Threshold sweep used {len(grid)} points: {[round(float(t), 2) for t in grid]}\n",
        f"Figures: `{forest_path}`, `{sweep_path}`, `{scatter_path}`\n",
        "See script docstring for the interpretation of the threshold-sweep floor as an estimate "
        "of AoU-native's real error rate against high-confidence ground truth.\n",
    ]
    md_text = "\n".join(md)
    md_path = os.path.join(args.out_dir, "confidence_report.md")
    with open(md_path, "w") as f:
        f.write(md_text)
    print(md_text)
    print(f"\n(written to {md_path} + 3 PNGs in {args.out_dir})", file=sys.stderr)


if __name__ == "__main__":
    main()
