#!/usr/bin/env python3
"""Builds the 3-way HLA comparison table (AoU-native / SpecHLA-SR / SpecImmune-LR)
for one person, for the 8 classical genes this pilot cares about.

Usage:
  python3 compare_hla_results.py <person_id> \
      [--spechla-result PATH] [--specimmune-result PATH] [--aou-tsv PATH]

Defaults assume the slice_and_fastq.sh + SpecHLA/SpecImmune conventions used in this
project (see ENVIRONMENT.md / slice_and_fastq.sh):
  --spechla-result    ~/pipeline_outputs/<person_id>/spechla_output/<person_id>/hla.result.txt
  --specimmune-result ~/pipeline_outputs/<person_id>/specimmune_output/<person_id>/<person_id>.HLA.final.type.result.formatted.txt
  --aou-tsv           ~/mnt/aou-controlled/v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv

Prints a markdown table to stdout and writes it to
~/pipeline_outputs/<person_id>/comparison.md -- paste that file's contents back into chat,
it gets folded into SMOKE_TEST_PICKS.local.md (gitignored -- never commit this file
or its contents, it's real participant genotype data).
"""
import argparse
import os
import re
import sys

import pandas as pd

CLASSICAL_GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]


def two_field(allele: str) -> str:
    """'HLA-A*11:01:01:01' or 'A*02:03:01' -> '11:01' -- for cross-method comparison,
    since different tools/databases report different field depths."""
    if not allele or allele in ("NA", "nan"):
        return "NA"
    m = re.search(r"\*([0-9]+):([0-9]+)", allele)
    return f"{m.group(1)}:{m.group(2)}" if m else allele


def load_aou_native(tsv_path: str, person_id: str) -> dict:
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
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


def load_specimmune(result_path: str) -> dict:
    if not os.path.exists(result_path):
        print(f"WARNING: SpecImmune result not found at {result_path}", file=sys.stderr)
        return {g: ("NA", "NA") for g in CLASSICAL_GENES}
    df = pd.read_csv(result_path, sep="\t", comment="#", dtype=str)
    out = {}
    for g in CLASSICAL_GENES:
        locus = f"HLA-{g}"
        sub = df[df["Locus"] == locus].sort_values("Chromosome")
        if len(sub) < 2:
            out[g] = ("NA", "NA")
            continue
        # Prefer One_guess; fall back to Genotype (may be a ;-separated ambiguous list)
        calls = []
        for _, r in sub.iterrows():
            val = r.get("One_guess", "NA")
            if val in (None, "NA", "nan") or pd.isna(val):
                val = str(r.get("Genotype", "NA")).split(";")[0]
            calls.append(val)
        out[g] = tuple(calls[:2])
    return out


def concordance_note(aou, spechla, specimmune) -> str:
    aou_set = {two_field(a) for a in aou}
    spechla_set = {two_field(a) for a in spechla}
    specimmune_set = {two_field(a) for a in specimmune}
    if aou_set == spechla_set == specimmune_set:
        return "3-way concordant"
    pairs_agree = []
    if aou_set == spechla_set:
        pairs_agree.append("AoU=SpecHLA")
    if aou_set == specimmune_set:
        pairs_agree.append("AoU=SpecImmune")
    if spechla_set == specimmune_set:
        pairs_agree.append("SpecHLA=SpecImmune")
    if pairs_agree:
        return f"Partial: {', '.join(pairs_agree)}, rest differ"
    return "All three differ"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("person_id")
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

    lines = [
        f"### Three-way comparison — person {pid}",
        "",
        "| Gene | AoU-native | SpecHLA (SR) | SpecImmune (LR) | Verdict |",
        "|---|---|---|---|---|",
    ]
    for g in CLASSICAL_GENES:
        a1, a2 = aou[g]
        s1, s2 = spechla[g]
        i1, i2 = specimmune[g]
        note = concordance_note((a1, a2), (s1, s2), (i1, i2))
        lines.append(f"| {g} | {a1} / {a2} | {s1} / {s2} | {i1} / {i2} | {note} |")

    out_text = "\n".join(lines)
    print(out_text)

    outdir = os.path.expanduser(f"~/pipeline_outputs/{pid}")
    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, "comparison.md")
    with open(out_path, "w") as f:
        f.write(out_text + "\n")
    print(f"\n(also written to {out_path} -- paste its contents back to Claude, "
          f"never commit this file, it's real genotype data)", file=sys.stderr)


if __name__ == "__main__":
    main()
