#!/usr/bin/env python3
"""Capstone figure: the six independent, converging lines of evidence that DRB1 is this
project's hardest/least-reliable HLA locus -- this project's single most-repeated finding
(context/DECISIONS.md, context/STATUS.md), never presented as one unified figure until now.

**This script is DIFFERENT IN KIND from the other three in this directory.** It does not read
live production output -- it synthesizes six ALREADY-PUBLISHED findings from earlier, smaller
pilots (Experiment D n=60, the AoU callset validation report, Experiment F n=60), each with its
own metric, cohort, and date. Every number below is transcribed, not recomputed, with an exact
citation to its source file/line. If a cited finding is ever revised, this script's EVIDENCE list
must be updated by hand to match -- it will not pick up the change automatically.

The six lines (chronological):
  1. Experiment D field-cascade (2026-07-12, DECISIONS.md / EXPERIMENTS.md): SpecHLA's DRB1
     errors vs SpecImmune-LR truth are REAL wrong-allele-family errors (49% partial-miss + 6%
     full-miss = 55%), not near-misses -- unlike AoU-native's DQA1 errors, which are all near-misses.
  2. AoU callset validation, Hardy-Weinberg proxy (2026-07-19, reports/aou_callset_validation/
     README.md lines 67-69): DRB1 observed/expected homozygosity ratio is elevated everywhere --
     pooled 1.39, and within every single ancestry group (AMR 1.40, EAS 1.40, MID 1.70, SAS 1.30)
     -- vs. DPA1 (a locus with genuinely high but population-structure-explained homozygosity)
     staying ~1.0 in every group. A population-genetics signal, independent of any cross-tool
     comparison.
  3. Experiment F cross-method concordance (2026-07-22, reports/immuannot_pilot/README.md):
     DRB1 Field 2 concordance is 40%, "roughly 25 points below the next-worst gene (A, 65%)";
     other genes range 65-94%.
  4. SpecImmune's own confidence flags (2026-07-22, same report): worst of any gene -- half of
     its DRB1 calls are flagged `lower_identity`.
  5. Immuannot's own confidence signal (2026-07-22, same report): `template_distance` is highest
     of any gene at DRB1, 1.82.
  6. Confidence-matched-truth collapse (2026-07-24, reports/confidence_matched_truth/README.md +
     reports/immuannot_pilot/README.md section 4): under a strict, math-grounded confidence
     filter, DRB1's confident-call count collapses (SpecImmune-truth 60->4, Immuannot-truth
     51->3(matched to 4)), and the direct SpecImmune-vs-Immuannot confident-overlap join drops to
     ZERO DRB1 pairs -- explicitly logged as "a sixth independent line of evidence."

No single one of these would be conclusive alone (per reports/immuannot_pilot/README.md's own
words); together they are. This script exists to make that convergence visible in one figure.

Writes: a 6-panel PNG + markdown table under --out-dir (no pipeline data read, so no --outroot/
--calls/--cohort args -- this is pure presentation of prior results).

Usage:
  python3 scripts/production_analysis/summarize_drb1_evidence_capstone.py
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DRB1_COLOR = "#C44E52"
OTHER_COLOR = "#B0B0B0"
DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")

EVIDENCE = [
    {
        "title": "1. Experiment D field-cascade\n(2026-07-12)",
        "citation": "context/DECISIONS.md, EXPERIMENTS.md 2026-07-12",
        "kind": "bar_pair",
        "ylabel": "Real wrong-family error rate (%)",
        "bars": [("DRB1\n(SpecHLA)", 55, DRB1_COLOR), ("DQA1\n(AoU-native)", 0, OTHER_COLOR)],
        "note": "DRB1: 49% partial-miss + 6% full-miss vs SpecImmune-LR truth.\n"
                "DQA1's AoU errors are all near-misses (0% real miss) -- contrast, not a same-locus pair.",
    },
    {
        "title": "2. AoU callset validation --\nHardy-Weinberg proxy (2026-07-19)",
        "citation": "reports/aou_callset_validation/README.md L67-69",
        "kind": "bar_group",
        "ylabel": "Observed/expected homozygosity ratio",
        "groups": ["Pooled", "AMR", "EAS", "MID", "SAS"],
        "drb1_vals": [1.39, 1.40, 1.40, 1.70, 1.30],
        "ref_val": 1.0,
        "ref_label": "DPA1 (real pop. structure, stays ~1.0 in every group)",
        "note": "DRB1 elevated pooled AND within every single ancestry group -- not a pooling (Wahlund) artifact.",
    },
    {
        "title": "3. Experiment F concordance\n(2026-07-22)",
        "citation": "reports/immuannot_pilot/README.md",
        "kind": "bar_pair",
        "ylabel": "Field 2 concordance (%)",
        "bars": [("DRB1", 40, DRB1_COLOR), ("A (next-worst)", 65, OTHER_COLOR)],
        "note": "Other 6 genes range 65-94%. DRB1 sits alone, ~25 points below the next-worst gene.",
    },
    {
        "title": "4. SpecImmune's own\nconfidence flags (2026-07-22)",
        "citation": "reports/immuannot_pilot/README.md",
        "kind": "single_stat",
        "value": "50%",
        "note": "Half of SpecImmune's DRB1 calls flagged `lower_identity` -- worst of any gene.",
    },
    {
        "title": "5. Immuannot's own\ntemplate_distance (2026-07-22)",
        "citation": "reports/immuannot_pilot/README.md",
        "kind": "single_stat",
        "value": "1.82",
        "note": "Highest mean template_distance of any gene (lower = more confident everywhere else).",
    },
    {
        "title": "6. Confidence-matched-truth\ncollapse (2026-07-24)",
        "citation": "reports/confidence_matched_truth/README.md,\nreports/immuannot_pilot/README.md sec.4",
        "kind": "before_after",
        "ylabel": "n confident DRB1 calls",
        "pairs": [("SpecImmune-truth", 60, 4), ("Immuannot-truth", 51, 3)],
        "note": "Direct SpecImmune-vs-Immuannot confident-overlap join: DRB1 drops to ZERO pairs.",
    },
]


def draw_panel(ax, spec):
    kind = spec["kind"]
    if kind == "bar_pair":
        labels = [b[0] for b in spec["bars"]]
        vals = [b[1] for b in spec["bars"]]
        colors = [b[2] for b in spec["bars"]]
        bars = ax.bar(labels, vals, color=colors)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v, f"{v}%", ha="center", va="bottom", fontsize=9)
        ax.set_ylabel(spec["ylabel"], fontsize=9)
        ax.set_ylim(0, max(vals) * 1.3 + 5)
    elif kind == "bar_group":
        x = range(len(spec["groups"]))
        ax.bar(x, spec["drb1_vals"], color=DRB1_COLOR)
        ax.axhline(spec["ref_val"], color=OTHER_COLOR, linestyle="--", linewidth=1.5,
                   label=spec["ref_label"])
        for xi, v in zip(x, spec["drb1_vals"]):
            ax.text(xi, v, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
        ax.set_xticks(list(x))
        ax.set_xticklabels(spec["groups"], fontsize=8)
        ax.set_ylabel(spec["ylabel"], fontsize=9)
        ax.legend(fontsize=7, frameon=False, loc="upper right")
    elif kind == "single_stat":
        ax.text(0.5, 0.55, spec["value"], ha="center", va="center", fontsize=32,
                color=DRB1_COLOR, transform=ax.transAxes, weight="bold")
        ax.axis("off")
    elif kind == "before_after":
        n = len(spec["pairs"])
        x = range(n)
        before = [p[1] for p in spec["pairs"]]
        after = [p[2] for p in spec["pairs"]]
        labels = [p[0] for p in spec["pairs"]]
        width = 0.35
        ax.bar([xi - width / 2 for xi in x], before, width=width, color=OTHER_COLOR, label="before filter")
        ax.bar([xi + width / 2 for xi in x], after, width=width, color=DRB1_COLOR, label="after filter")
        for xi, b, a in zip(x, before, after):
            ax.text(xi - width / 2, b, str(b), ha="center", va="bottom", fontsize=8)
            ax.text(xi + width / 2, a, str(a), ha="center", va="bottom", fontsize=8)
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel(spec["ylabel"], fontsize=9)
        ax.legend(fontsize=7, frameon=False)

    ax.set_title(spec["title"], fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)


def plot_capstone(out_path):
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    for ax, spec in zip(axes.flat, EVIDENCE):
        draw_panel(ax, spec)
        ax.text(0.5, -0.32, spec["note"], transform=ax.transAxes, ha="center", va="top",
                fontsize=7.5, color="#444", wrap=True)
        ax.text(0.5, -0.5 if spec["kind"] != "single_stat" else -0.15, spec["citation"],
                transform=ax.transAxes, ha="center", va="top", fontsize=6.5, color="#888",
                style="italic")
    fig.suptitle("DRB1 is this project's hardest locus -- six independent, converging lines of evidence",
                 fontsize=15, y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=os.path.join(DEFAULT_OUTROOT, "production_analysis", "drb1_capstone"))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    fig_path = os.path.join(args.out_dir, "drb1_evidence_capstone.png")
    plot_capstone(fig_path)

    md = ["# DRB1 -- six independent, converging lines of evidence\n",
          "This script synthesizes prior findings; it does not recompute them. Every number here "
          "is transcribed with a citation, not derived from live production data.\n"]
    for e in EVIDENCE:
        md.append(f"\n## {e['title'].replace(chr(10), ' ')}\n")
        md.append(f"Source: `{e['citation']}`\n\n{e['note']}\n")
    md.append(f"\nFigure: `{fig_path}`\n")
    md_text = "\n".join(md)
    md_path = os.path.join(args.out_dir, "drb1_capstone_report.md")
    with open(md_path, "w") as f:
        f.write(md_text)
    print(md_text)
    print(f"\n(written to {md_path} + 1 PNG in {args.out_dir})")


if __name__ == "__main__":
    main()
