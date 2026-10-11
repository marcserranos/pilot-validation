#!/usr/bin/env python3
"""Collects HLA calls (AoU-native / SpecHLA-SR / SpecImmune-LR) for one person into a
side-by-side matrix, for the 8 classical genes this pilot cares about.

This is a REPORTER, not a JUDGE: it does not compute concordance or declare a "winning"
method. It surfaces the raw calls plus SpecImmune's own confidence signals (identity,
read support, ambiguous-candidate-list size, whether One_guess was a tie) so a human (or
a later, separate analysis step) can weigh discordances with actual evidence instead of
a naive majority vote -- see DECISIONS.md for why vote-based "outlier" framing was wrong
(AoU-native and SpecHLA share the same short-read data and are not independent).

Usage:
  python3 compare_hla_results.py <person_id> [--run-label LABEL] \
      [--spechla-result PATH] [--specimmune-result PATH] [--aou-tsv PATH]

--run-label tags this run (e.g. "bwa_fullpad", "minimap2_tightpad") so results from
different SpecImmune configs (aligner / padding sweeps) don't overwrite each other and
can be aggregated later. Defaults to "default".

Defaults assume the slice_and_fastq.sh + SpecHLA/SpecImmune conventions used in this
project (see ENVIRONMENT.md / slice_and_fastq.sh):
  --spechla-result    ~/pipeline_outputs/<person_id>/spechla_output/<person_id>/hla.result.txt
  --specimmune-result ~/pipeline_outputs/<person_id>/specimmune_output/<person_id>/<person_id>.HLA.final.type.result.formatted.txt
  --aou-tsv           ~/mnt/aou-controlled/v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv

Writes two things per run:
  ~/pipeline_outputs/<person_id>/comparison_<run-label>.md   -- human-readable matrix
  ~/pipeline_outputs/<person_id>/comparison_log.csv           -- long-format row-per-gene,
    APPENDED (not overwritten) so a multi-config sweep accumulates one master log to
    analyze afterward. Never commit either file -- both contain real genotype data.
"""
import argparse
import csv
import os
import sys
from datetime import datetime, timezone

import pandas as pd

CLASSICAL_GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]

CSV_FIELDS = [
    "timestamp", "person_id", "run_label", "gene",
    "aou_1", "aou_2", "spechla_1", "spechla_2",
    "specimmune_1", "specimmune_2",
    "si_identity_1", "si_identity_2",
    "si_reads_1", "si_reads_2",
    "si_ambig_n_1", "si_ambig_n_2",
    "si_tied_1", "si_tied_2",
    "si_step1_1", "si_step1_2",
]


def load_aou_native(tsv_path: str, person_id: str) -> dict:
    if not os.path.exists(tsv_path):
        print(f"WARNING: AoU-native TSV not found at {tsv_path} (mount down?) "
              f"-- logging with NA for AoU columns", file=sys.stderr)
        return {g: ("NA", "NA") for g in CLASSICAL_GENES}
    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    except Exception as e:
        print(f"WARNING: failed to read AoU-native TSV ({e}) "
              f"-- logging with NA for AoU columns", file=sys.stderr)
        return {g: ("NA", "NA") for g in CLASSICAL_GENES}
    row = df[df["research_id"] == str(person_id)]
    if row.empty:
        print(f"WARNING: person {person_id} not found in {tsv_path}", file=sys.stderr)
        return {g: ("NA", "NA") for g in CLASSICAL_GENES}
    row = row.iloc[0]
    return {g: (row.get(f"{g}_1", "NA"), row.get(f"{g}_2", "NA")) for g in CLASSICAL_GENES}


def load_spechla(result_path: str) -> dict:
    if not os.path.exists(result_path):
        print(f"WARNING: SpecHLA result not found at {result_path}", file=sys.stderr)
        return {g: ("NA", "NA") for g in CLASSICAL_GENES}
    df = pd.read_csv(result_path, sep="\t", dtype=str, comment="#")
    row = df.iloc[0]
    return {g: (row.get(f"HLA_{g}_1", "NA"), row.get(f"HLA_{g}_2", "NA")) for g in CLASSICAL_GENES}


def _parse_match_info(match_info: str, allele: str):
    """Match_info is 'allele1|score1|identity1;allele2|score2|identity2;...'.
    Returns (identity, is_tied) for the given allele -- tied means another candidate
    shares its top score+identity (a One_guess pick among equals, not a clear winner)."""
    if not match_info or match_info in ("NA", "nan") or pd.isna(match_info):
        return None, False
    entries = []
    for part in str(match_info).split(";"):
        bits = part.split("|")
        if len(bits) == 3:
            a, score, ident = bits
            try:
                entries.append((a, float(score), float(ident)))
            except ValueError:
                continue
    if not entries:
        return None, False
    best_score = max(e[1] for e in entries)
    top_tier = [e for e in entries if e[1] == best_score]
    tied = len(top_tier) > 1
    for a, score, ident in entries:
        if a == allele:
            return ident, tied
    return None, tied


def load_specimmune(result_path: str) -> dict:
    """Returns {gene: {haplotype 1/2: {call, identity, reads, ambig_n, tied, step1}}}."""
    if not os.path.exists(result_path):
        print(f"WARNING: SpecImmune result not found at {result_path}", file=sys.stderr)
        return {g: {} for g in CLASSICAL_GENES}
    df = pd.read_csv(result_path, sep="\t", comment="#", dtype=str)
    out = {}
    for g in CLASSICAL_GENES:
        locus = f"HLA-{g}"
        sub = df[df["Locus"] == locus].sort_values("Chromosome")
        if len(sub) < 2:
            out[g] = {}
            continue
        haps = {}
        for hap_idx, (_, r) in enumerate(sub.iterrows(), start=1):
            call = r.get("One_guess", "NA")
            if call in (None, "NA", "nan") or pd.isna(call):
                call = str(r.get("Genotype", "NA")).split(";")[0]
            genotype_list = str(r.get("Genotype", "") or "")
            ambig_n = len([x for x in genotype_list.split(";") if x]) if genotype_list else 0
            identity, tied = _parse_match_info(r.get("Match_info", ""), call)
            step1 = r.get("Step1_type", "NA")
            haps[hap_idx] = {
                "call": call, "identity": identity, "reads": r.get("Reads_num", "NA"),
                "ambig_n": ambig_n, "tied": tied, "step1": step1,
            }
        out[g] = haps
    return out


def build_markdown(pid, run_label, aou, spechla, specimmune) -> str:
    lines = [
        f"### HLA call matrix — person {pid} — run: {run_label}",
        "",
        "Reporting only -- no verdict. SpecImmune columns include identity / read support / "
        "tie flag so discordances can be weighed with actual evidence (see DECISIONS.md).",
        "",
        "| Gene | AoU-native | SpecHLA (SR) | SpecImmune (LR) hap1 | hap2 |",
        "|---|---|---|---|---|",
    ]
    for g in CLASSICAL_GENES:
        a1, a2 = aou[g]
        s1, s2 = spechla[g]
        si = specimmune[g]
        h1 = si.get(1, {})
        h2 = si.get(2, {})

        def fmt_hap(h):
            if not h:
                return "NA"
            ident = f"{h['identity']:.4f}" if h["identity"] is not None else "NA"
            tie = " TIED" if h["tied"] else ""
            return f"{h['call']} (id={ident}, reads={h['reads']}, ambig={h['ambig_n']}{tie})"

        lines.append(f"| {g} | {a1} / {a2} | {s1} / {s2} | {fmt_hap(h1)} | {fmt_hap(h2)} |")
    return "\n".join(lines)


def append_csv_log(csv_path, pid, run_label, aou, spechla, specimmune):
    ts = datetime.now(timezone.utc).isoformat()
    is_new = not os.path.exists(csv_path)
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        if is_new:
            writer.writeheader()
        for g in CLASSICAL_GENES:
            a1, a2 = aou[g]
            s1, s2 = spechla[g]
            si = specimmune[g]
            h1 = si.get(1, {})
            h2 = si.get(2, {})
            writer.writerow({
                "timestamp": ts, "person_id": pid, "run_label": run_label, "gene": g,
                "aou_1": a1, "aou_2": a2, "spechla_1": s1, "spechla_2": s2,
                "specimmune_1": h1.get("call", "NA"), "specimmune_2": h2.get("call", "NA"),
                "si_identity_1": h1.get("identity"), "si_identity_2": h2.get("identity"),
                "si_reads_1": h1.get("reads"), "si_reads_2": h2.get("reads"),
                "si_ambig_n_1": h1.get("ambig_n"), "si_ambig_n_2": h2.get("ambig_n"),
                "si_tied_1": h1.get("tied"), "si_tied_2": h2.get("tied"),
                "si_step1_1": h1.get("step1"), "si_step1_2": h2.get("step1"),
            })


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("person_id")
    ap.add_argument("--run-label", default="default",
                     help="Tag for this SpecImmune config, e.g. 'bwa_fullpad', 'minimap2_tightpad'.")
    ap.add_argument("--spechla-result", default=None)
    ap.add_argument("--specimmune-result", default=None)
    ap.add_argument(
        "--aou-tsv",
        default=os.path.expanduser(
            "~/mnt/aou-controlled/v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv"
        ),
    )
    args = ap.parse_args()
    pid = args.person_id

    spechla_path = args.spechla_result or os.path.expanduser(
        f"~/pipeline_outputs/{pid}/spechla_output/{pid}/hla.result.txt"
    )
    specimmune_path = args.specimmune_result or os.path.expanduser(
        f"~/pipeline_outputs/{pid}/specimmune_output/{pid}/{pid}.HLA.final.type.result.formatted.txt"
    )

    aou = load_aou_native(args.aou_tsv, pid)
    spechla = load_spechla(spechla_path)
    specimmune = load_specimmune(specimmune_path)

    out_text = build_markdown(pid, args.run_label, aou, spechla, specimmune)
    print(out_text)

    outdir = os.path.expanduser(f"~/pipeline_outputs/{pid}")
    os.makedirs(outdir, exist_ok=True)
    md_path = os.path.join(outdir, f"comparison_{args.run_label}.md")
    with open(md_path, "w") as f:
        f.write(out_text + "\n")

    csv_path = os.path.join(outdir, "comparison_log.csv")
    append_csv_log(csv_path, pid, args.run_label, aou, spechla, specimmune)

    print(f"\n(written to {md_path} and appended to {csv_path} -- paste back to Claude, "
          f"never commit either, they contain real genotype data)", file=sys.stderr)


if __name__ == "__main__":
    main()
