#!/usr/bin/env python3
"""Experiment A: AoU-native HLA callset distribution stats. See EXPERIMENTS.md
'Roadmap locked 2026-07-09' for design. Read-only analysis, no participant-level
output leaves this process (prints aggregate stats only)."""
import argparse
import sys

import pandas as pd

CLASSICAL_GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
HARDEST_LOCI = ["DRB1", "DPB1"]  # our own hardest loci, per 2026-07-09 bake-off


def field_count(allele):
    if pd.isna(allele) or str(allele).strip() == "":
        return None
    return len(str(allele).split(":"))


def per_locus_stats(df, gene):
    c1, c2 = f"{gene}_1", f"{gene}_2"
    if c1 not in df.columns or c2 not in df.columns:
        return {"gene": gene, "error": f"columns {c1}/{c2} not found"}

    n = len(df)
    missing_either = df[c1].isna() | df[c2].isna()
    missing_both = df[c1].isna() & df[c2].isna()

    resolved = df.loc[~missing_either]
    fields1 = resolved[c1].map(field_count)
    fields2 = resolved[c2].map(field_count)
    all_fields = pd.concat([fields1, fields2])
    field_dist = all_fields.value_counts(normalize=True).sort_index()

    homozygous = (resolved[c1] == resolved[c2]).sum()

    return {
        "gene": gene,
        "n_total": n,
        "missing_either_pct": round(100 * missing_either.sum() / n, 2),
        "missing_both_pct": round(100 * missing_both.sum() / n, 2),
        "n_resolved": len(resolved),
        "field_dist": {int(k): round(100 * v, 1) for k, v in field_dist.items()},
        "homozygosity_pct": round(100 * homozygous / len(resolved), 2) if len(resolved) else None,
    }


def print_report(df, stats):
    n = len(df)
    print(f"## Experiment A — AoU-native distribution stats\n")
    print(f"Total rows (research_ids) in TSV: {n}\n")

    print("| Gene | Missing (either allele) | Missing (both) | Homozygosity | Field dist (of resolved) |")
    print("|---|---|---|---|---|")
    for s in stats:
        if "error" in s:
            print(f"| {s['gene']} | ERROR: {s['error']} | | | |")
            continue
        fd = ", ".join(f"{k}-field={v}%" for k, v in sorted(s["field_dist"].items()))
        hz = f"{s['homozygosity_pct']}%" if s["homozygosity_pct"] is not None else "n/a"
        print(f"| {s['gene']} | {s['missing_either_pct']}% | {s['missing_both_pct']}% | {hz} | {fd} |")

    print()
    hardest = [s for s in stats if s.get("gene") in HARDEST_LOCI and "error" not in s]
    others = [s for s in stats if s.get("gene") not in HARDEST_LOCI and "error" not in s]
    if hardest and others:
        avg_hardest = sum(s["missing_either_pct"] for s in hardest) / len(hardest)
        avg_others = sum(s["missing_either_pct"] for s in others) / len(others)
        print(f"**Hardest-loci cross-check:** DRB1/DPB1 avg missingness {avg_hardest:.2f}% "
              f"vs. other 6 genes avg {avg_others:.2f}%. "
              f"{'Consistent with AoU also struggling there (corroboration, not a methodology flaw).' if avg_hardest > avg_others else 'AoU is cleaner at our hardest loci than we are — may point to a methodology-specific weakness on our side, worth flagging.'}\n")

    all_resolved_mask = pd.Series(True, index=df.index)
    for gene in CLASSICAL_GENES:
        c1, c2 = f"{gene}_1", f"{gene}_2"
        if c1 in df.columns and c2 in df.columns:
            all_resolved_mask &= df[c1].notna() & df[c2].notna()
    n_fully_resolved = all_resolved_mask.sum()
    print(f"**Rows with all 8 classical genes fully resolved (both alleles present):** "
          f"{n_fully_resolved} / {n} ({100*n_fully_resolved/n:.2f}%)\n")


def ancestry_preview(df, ancestry_tsv_path, id_col="research_id"):
    print("## Ancestry join preview (bonus)\n")
    try:
        anc = pd.read_csv(ancestry_tsv_path, sep=None, engine="python")
    except Exception as e:
        print(f"Could not read ancestry TSV at {ancestry_tsv_path}: {e}", file=sys.stderr)
        print(f"SKIPPED — could not read {ancestry_tsv_path} (see stderr).\n")
        return

    anc_id_col = next((c for c in anc.columns if c.lower() in ("research_id", "person_id", "s")), None)
    ancestry_col = next((c for c in anc.columns if "ancestry" in c.lower()), None)
    if anc_id_col is None or ancestry_col is None:
        print(f"Ancestry TSV columns found: {list(anc.columns)}")
        print("SKIPPED — could not auto-identify id/ancestry columns; "
              "update `ancestry_preview()` with the real column names once seen above.\n")
        return

    merged = df.merge(anc[[anc_id_col, ancestry_col]], left_on=id_col, right_on=anc_id_col, how="inner")
    print(f"Joined {len(merged)} / {len(df)} rows to an ancestry label.\n")

    gene = "DQA1"  # the one robust cross-technology open finding — look here first
    c1, c2 = f"{gene}_1", f"{gene}_2"
    if c1 not in merged.columns:
        print(f"{gene} columns not found in TSV — skipping locus preview.\n")
        return

    print(f"**{gene} missingness/homozygosity by ancestry** (the AoU-outlier locus from the 3-person bake-off):\n")
    print("| Ancestry | n | Missing either | Homozygosity |")
    print("|---|---|---|---|")
    for anc_label, sub in merged.groupby(ancestry_col):
        missing = (sub[c1].isna() | sub[c2].isna())
        resolved = sub.loc[~missing]
        hz = (resolved[c1] == resolved[c2]).sum() / len(resolved) * 100 if len(resolved) else float("nan")
        print(f"| {anc_label} | {len(sub)} | {100*missing.sum()/len(sub):.2f}% | {hz:.2f}% |")
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tsv", required=True, help="Path to hla_genotypes.tsv (mounted path)")
    ap.add_argument("--ancestry-tsv", default=None, help="Path to AoU genetic-ancestry TSV, once located")
    ap.add_argument("--id-col", default="research_id")
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")

    stats = [per_locus_stats(df, gene) for gene in CLASSICAL_GENES]
    print_report(df, stats)

    if args.ancestry_tsv:
        ancestry_preview(df, args.ancestry_tsv, id_col=args.id_col)
    else:
        print("## Ancestry join preview (bonus)\n\nSKIPPED — no --ancestry-tsv path given "
              "(locate it first via gsutil ls, see EXPERIMENTS.md Experiment A quirks).\n")


if __name__ == "__main__":
    main()
