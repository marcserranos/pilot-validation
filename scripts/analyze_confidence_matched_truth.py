#!/usr/bin/env python3
"""Confidence-matched truth comparison -- Marc, 2026-07-24: "I want to be objective and find
what confidence interval we can tolerate, and then just cut everything that's below that
confidence... I don't want to treat different genes separately... If we choose whatever
threshold fits for SpecImmune and it leaves only half, we should adapt Immuannot so we only
get half too."

Design, exactly as specified:
  1. ONE global threshold per tool, applied identically across all 8 genes -- NOT tuned per
     gene (the per-gene percentile idea from the previous pass is explicitly rejected here).
  2. SpecImmune confident = NOT tied AND worst-haplotype si_identity >= --si-identity-threshold.
     Immuannot confident = NO real template_warning AND worst-haplotype template_distance
     <= --immuannot-distance-threshold. "Worst haplotype" = the genotype is only as confident
     as its less-confident haplotype (same logic as plot_confidence_distributions.py).
  3. SAMPLE-SIZE MATCHING: after applying each tool's own absolute threshold, if the two
     survivor counts differ, the LARGER set is rank-trimmed down to the smaller set's N --
     keeping only its own N most-confident calls (highest identity / lowest distance first),
     not a random subsample. This is the "adapt the other software so we only get half too"
     step. Both final truth sets end up the same size N.
  4. Per-gene retention is reported at 3 stages (all-resolved / after-threshold /
     after-N-matching) for BOTH tools -- diagnostic only. The threshold itself never changes
     per gene; this just shows which genes absorbed the cut, which is real information (e.g.
     DRB1 being hit hardest is itself a finding), not something to correct for.
  5. Fields 3-4 dropped from this analysis per Marc's own suggestion (Field 4 "doesn't give
     us nearly any information"; Field 3 carries the same DB-version-confound risk already
     logged in reports/immuannot_pilot/README.md). Field 1 is kept in the table for context;
     Field 2 is the plotted headline metric, same as everywhere else in this project.
  6. THE MERGED PLOT: one figure, Field 2 concordance by gene, AoU-native and SpecHLA each
     get their own subplot, and within each subplot every gene gets 4 grouped bars --
     SpecImmune-truth unfiltered, SpecImmune-truth confidence-matched, Immuannot-truth
     unfiltered, Immuannot-truth confidence-matched -- so both benchmarks' before/after are
     visible together instead of as two separate figures that are hard to eyeball against
     each other.

--si-identity-threshold and --immuannot-distance-threshold are REQUIRED, no default --
deliberately, so this never silently runs with a placeholder number. Fill them in from the
mathematically-grounded research pass (sequencing-error noise floor + IMGT inter-allele
spacing), not a guess.

Reads (no new VM run needed):
  ~/pipeline_outputs/experiment_d/cohort.tsv
  ~/pipeline_outputs/<pid>/comparison_log.csv   (aou_1/2, spechla_1/2, specimmune_1/2,
                                                  si_identity_1/2, si_tied_1/2)
  ~/pipeline_outputs/immuannot_calls.tsv        (immuannot_1/2)
  ~/pipeline_outputs/<pid>/immuannot_output/hap{1,2}.gtf.gz   (template_distance, template_warning)

Usage (via `pixi run -e spechla` for pandas+matplotlib):
  python3 scripts/analyze_confidence_matched_truth.py ~/pipeline_outputs/experiment_d/cohort.tsv \\
      --si-identity-threshold 0.999 --immuannot-distance-threshold 0
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
METHODS = ["aou", "sh"]
METHOD_LABEL = {"aou": "AoU-native", "sh": "SpecHLA"}
NULL = {"", "NA", "-", "nan", "None", ".", "na"}
ATTR_RE = re.compile(r'(\w+)\s+(?:"([^"]*)"|([^";]+))')

# Consistent with the project's existing palette (diagnose_immuannot_pilot.py): SpecImmune =
# orange family, Immuannot = blue family. Light shade = unfiltered, dark shade = confidence-
# matched, so "before vs after" reads as a shade change within the same color, not a new hue.
SI_UNFILTERED = "#F0B27A"
SI_MATCHED = "#B9540A"
IMM_UNFILTERED = "#A8C4E0"
IMM_MATCHED = "#2C5C8A"


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


def compare_allele(sr_fields, truth_fields):
    """Same 4-level cascade machinery as the other field-cascade scripts (kept so the best-
    pairing choice below still uses deep-field agreement to break ties) -- only levels 1-2
    are reported/plotted."""
    result, mismatched, unassessable = [], False, False
    for lvl in range(1, 5):
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
    return tuple((1 if res_a[i] is True else 0) + (1 if res_b[i] is True else 0) for i in range(4))


def compare_genotype(sr_pair, truth_pair):
    sr1, sr2 = parse_fields(sr_pair[0]), parse_fields(sr_pair[1])
    t1, t2 = parse_fields(truth_pair[0]), parse_fields(truth_pair[1])
    if sr1 is None or sr2 is None or t1 is None or t2 is None:
        return None
    pa = (compare_allele(sr1, t1), compare_allele(sr2, t2))
    pb = (compare_allele(sr1, t2), compare_allele(sr2, t1))
    return list(max([pa, pb], key=lambda p: pairing_score(p[0], p[1])))


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
    """{(person_id, gene): row-dict} with every field needed downstream: the 3 SR/consensus
    calls, SpecImmune's confidence pair, and Immuannot's confidence pair."""
    anc = dict(zip(cohort["person_id"].astype(str), cohort["ancestry"].astype(str)))

    # 1) comparison_log.csv -> aou/spechla/specimmune calls + SpecImmune confidence
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
                "ancestry": anc.get(pid),
                "aou_1": r.get("aou_1"), "aou_2": r.get("aou_2"),
                "spechla_1": r.get("spechla_1"), "spechla_2": r.get("spechla_2"),
                "specimmune_1": r.get("specimmune_1"), "specimmune_2": r.get("specimmune_2"),
                "si_worst_id": worst_id, "si_tied": tied,
                "immuannot_1": None, "immuannot_2": None,
                "imm_worst_td": None, "imm_warning": False,
            }

    # 2) immuannot_calls.tsv -> immuannot allele calls
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

    # 3) hap{1,2}.gtf.gz -> Immuannot confidence (template_distance / template_warning)
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


def resolved_keys(combined, truth_cols):
    """(person, gene) keys where the truth pair itself is a real, parseable genotype --
    the 'unfiltered' denominator, same discipline as every other script here."""
    out = []
    for key, row in combined.items():
        t1, t2 = row.get(truth_cols[0]), row.get(truth_cols[1])
        if parse_fields(t1) is not None and parse_fields(t2) is not None:
            out.append(key)
    return out


def confident_scores(combined, keys, tool, si_thresh, imm_thresh):
    """{(person,gene): score} for keys passing tool's absolute confidence bar. Score is
    'higher = more confident' for both tools so downstream ranking is uniform."""
    out = {}
    for key in keys:
        row = combined[key]
        if tool == "si":
            if row["si_tied"] or row["si_worst_id"] is None:
                continue
            if row["si_worst_id"] >= si_thresh:
                out[key] = row["si_worst_id"]
        else:
            if row["imm_warning"] or row["imm_worst_td"] is None:
                continue
            if row["imm_worst_td"] <= imm_thresh:
                out[key] = -row["imm_worst_td"]
    return out


def match_sample_size(si_scores, imm_scores):
    """Rank-trim the larger survivor set down to the smaller set's N, keeping its own most-
    confident calls (not a random subsample). Deterministic tie-break by key."""
    n = min(len(si_scores), len(imm_scores))
    si_sorted = sorted(si_scores.items(), key=lambda kv: (-kv[1], kv[0]))
    imm_sorted = sorted(imm_scores.items(), key=lambda kv: (-kv[1], kv[0]))
    return dict(si_sorted[:n]), dict(imm_sorted[:n]), n


def retention_report(combined, si_thresh, imm_thresh):
    si_all = resolved_keys(combined, ("specimmune_1", "specimmune_2"))
    imm_all = resolved_keys(combined, ("immuannot_1", "immuannot_2"))
    si_pass = confident_scores(combined, si_all, "si", si_thresh, imm_thresh)
    imm_pass = confident_scores(combined, imm_all, "imm", si_thresh, imm_thresh)
    si_matched, imm_matched, n = match_sample_size(si_pass, imm_pass)

    def per_gene_counts(keys):
        c = {g: 0 for g in GENES}
        for (_, gene) in keys:
            c[gene] += 1
        return c

    stages = {
        "si": {"all_resolved": per_gene_counts(si_all), "after_threshold": per_gene_counts(si_pass),
               "after_matching": per_gene_counts(si_matched)},
        "imm": {"all_resolved": per_gene_counts(imm_all), "after_threshold": per_gene_counts(imm_pass),
                "after_matching": per_gene_counts(imm_matched)},
    }

    lines = [f"## Retention by gene -- SpecImmune-truth vs Immuannot-truth "
             f"(matched to N={n} each)\n",
             f"Global thresholds used: SpecImmune si_identity >= {si_thresh} (and not tied); "
             f"Immuannot template_distance <= {imm_thresh} (and no template_warning). Same "
             f"threshold for every gene -- the per-gene counts below are a RESULT of that "
             f"single global bar, not a knob that was tuned per gene.\n",
             "| Gene | SpecImmune: all resolved | after threshold | after N-matching "
             "| Immuannot: all resolved | after threshold | after N-matching |",
             "|---|---|---|---|---|---|---|"]
    for g in GENES:
        lines.append(
            f"| {g} | {stages['si']['all_resolved'][g]} | {stages['si']['after_threshold'][g]} "
            f"| {stages['si']['after_matching'][g]} | {stages['imm']['all_resolved'][g]} "
            f"| {stages['imm']['after_threshold'][g]} | {stages['imm']['after_matching'][g]} |")

    si_total_all = sum(stages['si']['all_resolved'].values())
    si_total_thr = sum(stages['si']['after_threshold'].values())
    imm_total_all = sum(stages['imm']['all_resolved'].values())
    imm_total_thr = sum(stages['imm']['after_threshold'].values())
    si_share_thr = 100 * si_total_thr / si_total_all if si_total_all else float("nan")
    imm_share_thr = 100 * imm_total_thr / imm_total_all if imm_total_all else float("nan")
    si_share_final = 100 * n / si_total_all if si_total_all else float("nan")
    imm_share_final = 100 * n / imm_total_all if imm_total_all else float("nan")
    share_diff_final = abs(si_share_final - imm_share_final)

    lines.append(
        f"\n**Totals: SpecImmune {si_total_all} -> {si_total_thr} ({si_share_thr:.1f}% "
        f"retained by its own threshold) -> {n} ({si_share_final:.1f}% of its total pool); "
        f"Immuannot {imm_total_all} -> {imm_total_thr} ({imm_share_thr:.1f}% retained by its "
        f"own threshold) -> {n} ({imm_share_final:.1f}% of its total pool).**\n")
    lines.append(
        f"\n**SHARE-PARITY CHECK (Marc, 2026-07-24: 'make sure they are equal or nearly "
        f"equal, 80% vs 81% at most'): final retained share differs by {share_diff_final:.1f} "
        f"percentage points between the two truth sources.**\n")
    return ("\n".join(lines), stages, si_matched, imm_matched, n,
            si_share_thr, imm_share_thr, si_share_final, imm_share_final, share_diff_final)


def rate(n_true, n_false):
    d = n_true + n_false
    return (100 * n_true / d, d) if d else (float("nan"), 0)


def field12_rates(combined, truth_cols, allowed_keys):
    """{method: {gene: {field_level(1,2): (pct, n)}}} -- allowed_keys=None means no
    confidence restriction (the 'unfiltered' baseline); otherwise restrict to that key set."""
    out = {m: {g: {} for g in GENES} for m in METHODS}
    counts = {m: {g: {1: [0, 0], 2: [0, 0]} for g in GENES} for m in METHODS}  # [n_true, n_false]
    for key, row in combined.items():
        gene = key[1]
        if gene not in GENES:
            continue
        truth_pair = (row.get(truth_cols[0]), row.get(truth_cols[1]))
        if parse_fields(truth_pair[0]) is None or parse_fields(truth_pair[1]) is None:
            continue
        if allowed_keys is not None and key not in allowed_keys:
            continue
        for method, cols in (("aou", ("aou_1", "aou_2")), ("sh", ("spechla_1", "spechla_2"))):
            sr_pair = (row.get(cols[0]), row.get(cols[1]))
            res = compare_genotype(sr_pair, truth_pair)
            if res is None:
                continue
            for allele_res in res:
                for lvl in (1, 2):
                    v = allele_res[lvl - 1]
                    if v is True:
                        counts[method][gene][lvl][0] += 1
                    elif v is False:
                        counts[method][gene][lvl][1] += 1
    for m in METHODS:
        for g in GENES:
            for lvl in (1, 2):
                nt, nf = counts[m][g][lvl]
                out[m][g][lvl] = rate(nt, nf)
    return out


def plot_merged_field2(rates_by_scenario, path):
    """rates_by_scenario: {'si_unf':..., 'si_matched':..., 'imm_unf':..., 'imm_matched':...},
    each a field12_rates()-shaped dict. One figure, one subplot per SR method, 4 grouped bars
    per gene -- the 'overlaid, not side by side' merged view."""
    scenario_order = ["si_unf", "si_matched", "imm_unf", "imm_matched"]
    colors = {"si_unf": SI_UNFILTERED, "si_matched": SI_MATCHED,
              "imm_unf": IMM_UNFILTERED, "imm_matched": IMM_MATCHED}
    labels = {"si_unf": "vs SpecImmune (unfiltered)", "si_matched": "vs SpecImmune (confident, N-matched)",
              "imm_unf": "vs Immuannot (unfiltered)", "imm_matched": "vs Immuannot (confident, N-matched)"}

    fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
    width = 0.2
    x = range(len(GENES))
    for ax, method in zip(axes, METHODS):
        for i, scen in enumerate(scenario_order):
            heights, ns = [], []
            for gene in GENES:
                pct, n = rates_by_scenario[scen][method][gene][2]
                heights.append(0 if pct != pct else pct)
                ns.append(n)
            offs = [xi + (i - 1.5) * width for xi in x]
            bars = ax.bar(offs, heights, width=width, color=colors[scen],
                          label=labels[scen] if method == METHODS[0] else None)
            for xoff, h, n in zip(offs, heights, ns):
                ax.text(xoff, h + 1.5, f"{n}", ha="center", va="bottom", fontsize=6, color="#555",
                        rotation=90)
        ax.set_title(f"{METHOD_LABEL[method]} -- Field 2 (protein) concordance", fontsize=10)
        ax.set_ylim(0, 118)
        ax.set_yticks(range(0, 101, 20))
        ax.set_ylabel("Match rate (%)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xticks(list(x))
    axes[-1].set_xticklabels(GENES)
    fig.suptitle("Field 2 concordance vs SpecImmune-truth and Immuannot-truth, before/after a "
                 "single global confidence threshold + sample-size matching\n(small numbers above "
                 "each bar = n)", fontsize=10, y=0.995)
    handles, labs = axes[0].get_legend_handles_labels()
    fig.legend(handles, labs, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.94),
               fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_retention(stages, n_matched, path):
    """Per-gene retention diagnostic -- which genes absorbed the cut, for both tools, at each
    stage. Purely informational (per Marc: not used to adjust the threshold)."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)
    stage_order = ["all_resolved", "after_threshold", "after_matching"]
    stage_color = {"all_resolved": "#CCCCCC", "after_threshold": "#888888", "after_matching": None}
    tool_color = {"si": SI_MATCHED, "imm": IMM_MATCHED}
    tool_label = {"si": "SpecImmune", "imm": "Immuannot"}
    for ax, tool in zip(axes, ["si", "imm"]):
        width = 0.25
        x = range(len(GENES))
        for i, stage in enumerate(stage_order):
            heights = [stages[tool][stage][g] for g in GENES]
            color = stage_color[stage] if stage != "after_matching" else tool_color[tool]
            offs = [xi + (i - 1) * width for xi in x]
            ax.bar(offs, heights, width=width, color=color,
                   label=stage.replace("_", " ") if tool == "si" else None)
        ax.set_xticks(list(x))
        ax.set_xticklabels(GENES)
        ax.set_title(f"{tool_label[tool]} truth -- retention by gene (final N={n_matched} total)",
                     fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("n (person, gene) calls")
    handles, labs = axes[0].get_legend_handles_labels()
    fig.legend(handles, labs, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02),
               fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def field12_table(rates, label):
    lines = [f"### {label}\n", "| Gene | " + " | ".join(METHOD_LABEL[m] + " F1" for m in METHODS) +
             " | " + " | ".join(METHOD_LABEL[m] + " F2" for m in METHODS) + " |",
             "|---|" + "---|" * (2 * len(METHODS))]
    for g in GENES:
        cells = []
        for m in METHODS:
            pct, n = rates[m][g][1]
            cells.append(f"{pct:.0f}% (n={n})" if n else "— (no data)")
        for m in METHODS:
            pct, n = rates[m][g][2]
            cells.append(f"{pct:.0f}% (n={n})" if n else "— (no data)")
        lines.append(f"| {g} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cohort", help="cohort.tsv from build_experiment_d_cohort.py")
    ap.add_argument("--si-identity-threshold", type=float, required=True,
                    help="Global SpecImmune worst-haplotype si_identity cutoff (e.g. 0.999). "
                         "No default on purpose -- use the mathematically-grounded value.")
    ap.add_argument("--immuannot-distance-threshold", type=float, required=True,
                    help="Global Immuannot worst-haplotype template_distance cutoff (e.g. 0). "
                         "No default on purpose -- use the mathematically-grounded value.")
    ap.add_argument("--outroot", default=os.path.expanduser("~/pipeline_outputs"))
    ap.add_argument("--analysis-dir", default=None)
    ap.add_argument("--max-share-diff-pct", type=float, default=5.0,
                    help="If the two truth sources' final retained share (n_matched / its own "
                         "total resolved pool) differs by more than this many percentage "
                         "points, print a loud warning (does not abort -- you still get the "
                         "numbers and figures either way). Default 5.0.")
    args = ap.parse_args()

    cohort = pd.read_csv(args.cohort, sep="\t", dtype=str)
    combined = load_combined(cohort, args.outroot)

    (retention_md, stages, si_matched, imm_matched, n, si_share_thr, imm_share_thr,
     si_share_final, imm_share_final, share_diff_final) = retention_report(
        combined, args.si_identity_threshold, args.immuannot_distance_threshold)
    if n == 0:
        print("FATAL: at least one truth source has ZERO calls passing its threshold -- "
              "loosen the threshold(s) and rerun.", file=sys.stderr)
        sys.exit(1)

    print(f"\nSHARE CHECK: SpecImmune retains {si_share_final:.1f}% of its pool, Immuannot "
          f"retains {imm_share_final:.1f}% of its pool (post-matching, N={n} each). "
          f"Difference: {share_diff_final:.1f} points.", file=sys.stderr)
    if share_diff_final > args.max_share_diff_pct:
        print(f"*** WARNING: share difference ({share_diff_final:.1f} pts) exceeds "
              f"--max-share-diff-pct ({args.max_share_diff_pct}). One truth source is being "
              f"squeezed much harder than the other RELATIVE TO ITS OWN POOL, even though the "
              f"absolute N now matches -- that is exactly the bias Marc flagged. Do not treat "
              f"this run's comparison as apples-to-apples. Fix: tighten (never loosen) the "
              f"MORE LENIENT tool's own threshold (the one with the higher 'retained by its "
              f"own threshold' % above) and rerun -- do not adjust per gene. ***", file=sys.stderr)
    else:
        print(f"Share difference within tolerance ({args.max_share_diff_pct} pts) -- the two "
              f"truth sources are cutting a comparable fraction of their own pool.",
              file=sys.stderr)

    rates = {
        "si_unf": field12_rates(combined, ("specimmune_1", "specimmune_2"), None),
        "si_matched": field12_rates(combined, ("specimmune_1", "specimmune_2"), set(si_matched)),
        "imm_unf": field12_rates(combined, ("immuannot_1", "immuannot_2"), None),
        "imm_matched": field12_rates(combined, ("immuannot_1", "immuannot_2"), set(imm_matched)),
    }

    expdir = os.path.join(args.outroot, "experiment_d")
    adir = args.analysis_dir or os.path.join(expdir, "analysis_confidence_matched")
    os.makedirs(adir, exist_ok=True)

    plot_merged_field2(rates, os.path.join(adir, "field2_merged_confidence_matched.png"))
    plot_retention(stages, n, os.path.join(adir, "retention_by_gene.png"))

    md = ("# Confidence-matched truth comparison\n\n"
          f"Global thresholds: SpecImmune si_identity >= {args.si_identity_threshold} "
          f"(not tied); Immuannot template_distance <= {args.immuannot_distance_threshold} "
          f"(no template_warning). Both truth sets rank-trimmed to the same N={n}.\n\n"
          + retention_md + "\n"
          + field12_table(rates["si_unf"], "AoU/SpecHLA vs SpecImmune-truth -- UNFILTERED") + "\n\n"
          + field12_table(rates["si_matched"], f"AoU/SpecHLA vs SpecImmune-truth -- "
                          f"confident + N-matched (N={n})") + "\n\n"
          + field12_table(rates["imm_unf"], "AoU/SpecHLA vs Immuannot-truth -- UNFILTERED") + "\n\n"
          + field12_table(rates["imm_matched"], f"AoU/SpecHLA vs Immuannot-truth -- "
                          f"confident + N-matched (N={n})") + "\n")
    md_path = os.path.join(adir, "confidence_matched_truth.md")
    with open(md_path, "w") as f:
        f.write(md)

    print(md)
    print(f"\n(written to {md_path} + 2 PNGs in {adir} -- aggregate only, keep on the VM "
          f"per the egress caveat; paste table/describe PNGs back)", file=sys.stderr)


if __name__ == "__main__":
    main()
