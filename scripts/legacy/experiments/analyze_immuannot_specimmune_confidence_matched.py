#!/usr/bin/env python3
"""Confidence-filtered direct concordance: SpecImmune-LR vs Immuannot -- Marc, 2026-07-24:
"we saw the parts with less confidence got the highest discrepancy... cut down the calls with
less than X confidence, same threshold, equalize so not biased, compare SpecImmune vs Immuannot
directly, before/after plot."

This is a DIRECT pairwise comparison (unlike analyze_confidence_matched_truth.py, which scores
AoU/SpecHLA against two SEPARATE truth pools and needs sample-size matching between two
different comparison populations). Here both callers are compared on the SAME (person, gene)
slots in both the before and after states, so restricting to "both callers confident at this
slot" is already an apples-to-apples join -- no artificial rank-trimming needed. What IS worth
checking (and reported below) is that neither threshold is wildly more restrictive than the
other -- i.e., the join isn't silently dominated by one caller's abstention.

Thresholds are the SAME mathematically-grounded global bars already established and used in
analyze_confidence_matched_truth.py (sequencing-error-rate math for SpecImmune-LR identity,
IMGT minimum inter-allele-spacing for Immuannot template_distance) -- not re-derived here.
Required CLI args, no default, so this can't silently run with a placeholder.

Reads (no new VM run needed):
  ~/pipeline_outputs/experiment_d/cohort.tsv
  ~/pipeline_outputs/<pid>/comparison_log.csv          (specimmune_1/2, si_identity_1/2, si_tied_1/2)
  ~/pipeline_outputs/immuannot_calls.tsv                (immuannot_1/2)
  ~/pipeline_outputs/<pid>/immuannot_output/hap{1,2}.gtf.gz   (template_distance, template_warning)

Usage (via `pixi run -e spechla` for pandas+matplotlib):
  python3 scripts/analyze_immuannot_specimmune_confidence_matched.py \\
      ~/pipeline_outputs/experiment_d/cohort.tsv \\
      --si-identity-threshold 0.9995 --immuannot-distance-threshold 0
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
NULL = {"", "NA", "-", "nan", "None", ".", "na"}
ATTR_RE = re.compile(r'(\w+)\s+(?:"([^"]*)"|([^";]+))')

# Same slots/roles as the project's established palette (diagnose_immuannot_pilot.py):
# Immuannot=blue, SpecImmune=orange, DRB1 gets its own highlight. Muted = unfiltered ("before"),
# saturated + black outline = confidence-matched ("after") -- same visual grammar as
# analyze_confidence_matched_truth.py's plot_merged_field2.
BEFORE_COLOR = "#B7C6DA"
AFTER_COLOR = "#4C72B0"
DRB1_HIGHLIGHT = "#C44E52"


def parse_fields(raw):
    if raw is None:
        return None
    s = str(raw).strip()
    if s in NULL or (hasattr(pd, "isna") and pd.isna(raw)):
        return None
    if "*" in s:
        s = s.split("*", 1)[1]
    fields = [f for f in s.split(":") if f != ""]
    for i, f in enumerate(fields):
        if f.strip().lower() == "new":
            fields = fields[:i]
            break
    return fields or None


def compare_allele(a_fields, b_fields):
    result, mismatched, unassessable = [], False, False
    for lvl in range(1, 5):
        if mismatched:
            result.append(False); continue
        if unassessable:
            result.append(None); continue
        if lvl > len(a_fields) or lvl > len(b_fields):
            result.append(None); unassessable = True; continue
        if a_fields[lvl - 1] == b_fields[lvl - 1]:
            result.append(True)
        else:
            result.append(False); mismatched = True
    return result


def pairing_score(res_a, res_b):
    return tuple((1 if res_a[i] is True else 0) + (1 if res_b[i] is True else 0) for i in range(4))


def compare_genotype(a_pair, b_pair):
    a1, a2 = parse_fields(a_pair[0]), parse_fields(a_pair[1])
    b1, b2 = parse_fields(b_pair[0]), parse_fields(b_pair[1])
    if a1 is None or a2 is None or b1 is None or b2 is None:
        return None
    pa = (compare_allele(a1, b1), compare_allele(a2, b2))
    pb = (compare_allele(a1, b2), compare_allele(a2, b1))
    return list(max([pa, pb], key=lambda p: pairing_score(p[0], p[1])))


def geno_status_at(cascade, lvl):
    a, b = cascade[0][lvl - 1], cascade[1][lvl - 1]
    if a is False or b is False:
        return "discordant"
    if a is None or b is None:
        return "na"
    return "concordant"


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


def load_combined(cohort, outroot):
    """{(person_id, gene): row-dict} -- specimmune calls + confidence, immuannot calls +
    confidence. Same fields/logic as analyze_confidence_matched_truth.py's load_combined, minus
    the aou/spechla columns (not needed for a direct SpecImmune-vs-Immuannot comparison)."""
    combined = {}
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
            if gene not in GENES:
                continue

            def truthy(x):
                return str(x).strip().lower() in {"true", "1", "yes"}
            tied = truthy(r.get("si_tied_1")) or truthy(r.get("si_tied_2"))
            try:
                id1, id2 = float(r.get("si_identity_1")), float(r.get("si_identity_2"))
                worst_id = min(id1, id2)
            except (TypeError, ValueError):
                worst_id = None
            combined[(pid, gene)] = {
                "specimmune_1": r.get("specimmune_1"), "specimmune_2": r.get("specimmune_2"),
                "si_worst_id": worst_id, "si_tied": tied,
                "immuannot_1": None, "immuannot_2": None,
                "imm_worst_td": None, "imm_warning": False,
            }

    calls_path = os.path.join(outroot, "immuannot_calls.tsv")
    if os.path.exists(calls_path):
        calls = pd.read_csv(calls_path, sep="\t", dtype=str, keep_default_na=False)
        calls["gene_bare"] = calls["gene"].str.replace("^HLA-", "", regex=True)
        for _, r in calls.iterrows():
            key = (str(r["person_id"]), str(r["gene_bare"]))
            if key not in combined:
                continue
            combined[key]["immuannot_1"] = r["immuannot_1"]
            combined[key]["immuannot_2"] = r["immuannot_2"]

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
        for gene in GENES:
            key = (pid, gene)
            if key not in combined:
                continue
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
            combined[key]["imm_worst_td"] = max(tds) if tds else None
            combined[key]["imm_warning"] = bool(warnings)

    return combined


def both_resolved_keys(combined):
    out = []
    for key, row in combined.items():
        si_ok = parse_fields(row.get("specimmune_1")) is not None and parse_fields(row.get("specimmune_2")) is not None
        imm_ok = parse_fields(row.get("immuannot_1")) is not None and parse_fields(row.get("immuannot_2")) is not None
        if si_ok and imm_ok:
            out.append(key)
    return out


def si_confident(row, si_thresh):
    if row["si_tied"] or row["si_worst_id"] is None:
        return False
    return row["si_worst_id"] >= si_thresh


def imm_confident(row, imm_thresh):
    if row["imm_warning"] or row["imm_worst_td"] is None:
        return False
    return row["imm_worst_td"] <= imm_thresh


def concordance_by_gene(combined, keys, lvl):
    """{gene: (pct, n)} + overall, Field `lvl` concordance among the given (both-resolved) keys."""
    out, tot_c, tot_d = {}, 0, 0
    for g in GENES:
        c = d = 0
        for key in keys:
            if key[1] != g:
                continue
            row = combined[key]
            casc = compare_genotype((row["immuannot_1"], row["immuannot_2"]),
                                     (row["specimmune_1"], row["specimmune_2"]))
            if casc is None:
                continue
            st = geno_status_at(casc, lvl)
            if st == "concordant":
                c += 1
            elif st == "discordant":
                d += 1
        n = c + d
        out[g] = (100 * c / n if n else None, n)
        tot_c += c; tot_d += d
    tot_n = tot_c + tot_d
    out["overall"] = (100 * tot_c / tot_n if tot_n else None, tot_n)
    return out


def plot_before_after(before, after, path):
    """Field 2 concordance per gene, before (muted) vs after (saturated + outline) the
    confidence join -- same visual grammar as the rest of this project's confidence plots."""
    genes_plus = GENES + ["overall"]
    b_vals = [before[g][0] if before[g][0] is not None else 0 for g in genes_plus]
    b_ns = [before[g][1] for g in genes_plus]
    a_vals = [after[g][0] if after[g][0] is not None else 0 for g in genes_plus]
    a_ns = [after[g][1] for g in genes_plus]

    fig, ax = plt.subplots(figsize=(10, 5))
    x = range(len(genes_plus))
    width = 0.35
    b_colors = [DRB1_HIGHLIGHT if g == "DRB1" else BEFORE_COLOR for g in genes_plus]
    a_colors = [DRB1_HIGHLIGHT if g == "DRB1" else AFTER_COLOR for g in genes_plus]
    bars_b = ax.bar([xi - width / 2 for xi in x], b_vals, width=width, color=b_colors,
                     label="unfiltered (all resolved)")
    bars_a = ax.bar([xi + width / 2 for xi in x], a_vals, width=width, color=a_colors,
                     edgecolor="black", linewidth=1.3, label="confidence-matched (both confident)")
    for bar, v, n in zip(bars_b, b_vals, b_ns):
        ax.text(bar.get_x() + bar.get_width() / 2, v + 2, f"{v:.0f}%\n(n={n})",
                ha="center", va="bottom", fontsize=7)
    for bar, v, n in zip(bars_a, a_vals, a_ns):
        ax.text(bar.get_x() + bar.get_width() / 2, v + 2, f"{v:.0f}%\n(n={n})",
                ha="center", va="bottom", fontsize=7)
    ax.axvline(len(GENES) - 0.5, color="gray", linestyle=":", linewidth=1)
    ax.set_xticks(list(x))
    ax.set_xticklabels(genes_plus)
    ax.set_ylim(0, 122)
    ax.set_yticks(range(0, 101, 20))
    ax.set_ylabel("Field 2 (protein-level) concordance, %")
    ax.set_title("SpecImmune-LR vs Immuannot -- direct concordance, before/after confidence "
                  "filtering\n(cross-method agreement, NOT accuracy -- no ground truth)",
                  fontsize=10, pad=12)
    ax.legend(fontsize=8, loc="lower left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cohort", help="cohort.tsv from build_experiment_d_cohort.py")
    ap.add_argument("--si-identity-threshold", type=float, required=True,
                     help="Global SpecImmune worst-haplotype si_identity cutoff. No default -- "
                          "reuse the mathematically-grounded value from "
                          "analyze_confidence_matched_truth.py (0.9995).")
    ap.add_argument("--immuannot-distance-threshold", type=float, required=True,
                     help="Global Immuannot worst-haplotype template_distance cutoff. No "
                          "default -- reuse the mathematically-grounded value (0).")
    ap.add_argument("--outroot", default=os.path.expanduser("~/pipeline_outputs"))
    ap.add_argument("--analysis-dir", default=None)
    args = ap.parse_args()

    cohort = pd.read_csv(args.cohort, sep="\t", dtype=str)
    combined = load_combined(cohort, args.outroot)
    total_possible = len(cohort) * len(GENES)

    all_keys = both_resolved_keys(combined)
    si_pass = {k for k in all_keys if si_confident(combined[k], args.si_identity_threshold)}
    imm_pass = {k for k in all_keys if imm_confident(combined[k], args.immuannot_distance_threshold)}
    both_pass = si_pass & imm_pass

    if not both_pass:
        print("FATAL: zero (person, gene) slots pass BOTH confidence bars -- loosen the "
              "threshold(s) and rerun.", file=sys.stderr)
        sys.exit(1)

    before = concordance_by_gene(combined, all_keys, 2)
    after = concordance_by_gene(combined, both_pass, 2)

    # Retention/bias check -- not a rank-trim (this is one paired population, not two separate
    # ones), just confirms neither threshold is silently dominating the join.
    lines = ["# SpecImmune-LR vs Immuannot -- confidence-filtered direct concordance\n",
              f"Global thresholds (same as `analyze_confidence_matched_truth.py`): SpecImmune "
              f"si_identity >= {args.si_identity_threshold} (not tied); Immuannot "
              f"template_distance <= {args.immuannot_distance_threshold} (no template_warning). "
              f"Applied to BOTH callers' own calls, then joined on 'both confident at this "
              f"(person, gene) slot' -- a direct pairwise comparison needs the same slot on both "
              f"sides, so this join is inherently apples-to-apples (unlike scoring two "
              f"SEPARATE truth pools, which needed the count-matching in "
              f"analyze_confidence_matched_truth.py).\n",
              "## Retention (bias check)\n",
              f"Of {len(all_keys)} slots where both callers produced a genotype: "
              f"SpecImmune confident at {len(si_pass)} ({100*len(si_pass)/len(all_keys):.1f}%), "
              f"Immuannot confident at {len(imm_pass)} ({100*len(imm_pass)/len(all_keys):.1f}%), "
              f"BOTH confident (the comparison set used below) at {len(both_pass)} "
              f"({100*len(both_pass)/len(all_keys):.1f}%). If the two per-caller percentages "
              f"above are far apart, the join is being pulled toward one caller's stricter "
              f"behavior -- check before trusting the after-numbers.\n",
              "| Gene | both resolved | SI confident | Immuannot confident | BOTH confident |",
              "|---|---|---|---|---|"]
    for g in GENES:
        n_all = sum(1 for k in all_keys if k[1] == g)
        n_si = sum(1 for k in si_pass if k[1] == g)
        n_imm = sum(1 for k in imm_pass if k[1] == g)
        n_both = sum(1 for k in both_pass if k[1] == g)
        lines.append(f"| {g} | {n_all} | {n_si} | {n_imm} | {n_both} |")

    lines += ["", "## Field 2 (protein-level) concordance, before vs after\n",
              "| Gene | unfiltered (all resolved) | confidence-matched (both confident) |",
              "|---|---|---|"]
    for g in GENES + ["overall"]:
        b_pct, b_n = before[g]
        a_pct, a_n = after[g]
        b_str = f"{b_pct:.0f}% (n={b_n})" if b_pct is not None else "- (no data)"
        a_str = f"{a_pct:.0f}% (n={a_n})" if a_pct is not None else "- (no data)"
        lines.append(f"| {g} | {b_str} | {a_str} |")

    b_overall, a_overall = before["overall"][0], after["overall"][0]
    if b_overall is not None and a_overall is not None:
        lines += ["", f"**Overall Field 2 concordance moves from {b_overall:.1f}% (unfiltered, "
                       f"n={before['overall'][1]}) to {a_overall:.1f}% (confidence-matched, "
                       f"n={after['overall'][1]}) once both callers are held to a strict, "
                       f"mathematically-grounded confidence bar.**"]

    adir = args.analysis_dir or os.path.join(args.outroot, "immuannot_pilot", "analysis")
    os.makedirs(adir, exist_ok=True)
    fig_path = os.path.join(adir, "figures", "concordance_before_after_confidence.png")
    os.makedirs(os.path.dirname(fig_path), exist_ok=True)
    plot_before_after(before, after, fig_path)

    md = "\n".join(lines) + "\n"
    md_path = os.path.join(adir, "immuannot_specimmune_confidence_matched.md")
    with open(md_path, "w") as f:
        f.write(md)
    print(md)
    print(f"\n(written to {md_path}, figure at {fig_path} -- aggregate-only, keep raw calls on "
          f"the VM per the egress caveat; paste table/figure back)", file=sys.stderr)


if __name__ == "__main__":
    main()
