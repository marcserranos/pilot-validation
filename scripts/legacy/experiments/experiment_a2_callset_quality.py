#!/usr/bin/env python3
"""Experiment A2: AoU-native HLA callset QUALITY & DISTRIBUTION characterization.

Complements experiment_a_aou_stats.py (which covers missingness / resolution / homozygosity).
This one produces the population-genetic content a reliability report needs, WITHOUT any
per-sample ground truth, by leaning on properties real HLA data must satisfy:

  1. Allelic DIVERSITY per locus (distinct 2-field alleles) — the MHC is hyperpolymorphic;
     a suspiciously small allele set would signal calling collapse.
  2. Allele FREQUENCY spectra, overall and per genetic-ancestry group — the raw material for
     an EXTERNAL concordance check against published population frequencies (e.g. Allele
     Frequency Net Database): if AoU's EUR A*02:01 ~ the known ~28-30%, the callset is behaving.
  3. Hardy-Weinberg PROXY: observed vs. expected homozygosity per locus per ancestry, where
     expected = sum(freq_i^2) at 2-field. Observed >> expected => excess homozygosity, a classic
     signature of allele dropout / reference bias (or real consanguinity/selection — flagged, not
     auto-attributed). Observed ~ expected => consistent with a clean callset.
  4. Rare-allele burden — share of the allele set that is rare (<1%) / singleton, per locus.

All output is aggregate (allele frequencies, counts) — NO participant-level rows — so it is
egress-safe under the same caveat as experiment_a_aou_stats.py. Frequencies are computed at
2-FIELD (collapsing AoU's 2-3 field calls) for comparability with reference databases and to
avoid splitting counts across resolution depths.

Usage (via `pixi run -e spechla` for pandas):
  python3 experiment_a2_callset_quality.py \
      --tsv ~/mnt/aou-controlled/v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv \
      --ancestry-tsv ~/mnt/aou-controlled/v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv \
      --out-dir ~/pipeline_outputs/experiment_a2
"""
import argparse
import os
import sys

import pandas as pd

CLASSICAL_GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
TOP_N = 5          # top alleles to print per locus
RARE_THRESHOLD = 0.01  # <1% = "rare"


def to2field(allele):
    """Normalize any AoU call to a 2-field key: 'A*01:01:01'->'01:01', '01:01'->'01:01'.
    Returns None if uncallable or below 2 fields."""
    if pd.isna(allele):
        return None
    a = str(allele).strip()
    if a == "" or a.lower() in {"na", "nan", "none", "-", "."}:
        return None
    if "*" in a:
        a = a.split("*", 1)[1]
    fields = a.split(":")
    if len(fields) < 2:
        return None
    return f"{fields[0]}:{fields[1]}"


def locus_alleles(df, gene):
    """Return a Series of all 2-field alleles pooled from both haplotype columns, dropping
    individuals uncallable at either allele (so each contributes 0 or 2 alleles — HWE-clean)."""
    c1, c2 = f"{gene}_1", f"{gene}_2"
    a1 = df[c1].map(to2field)
    a2 = df[c2].map(to2field)
    both = a1.notna() & a2.notna()
    return a1[both], a2[both]


def hwe_proxy(a1, a2):
    """Observed vs expected 2-field homozygosity. Expected = sum(p_i^2) over 2-field freqs."""
    n = len(a1)
    if n == 0:
        return None
    obs_hom = (a1.values == a2.values).mean()
    pooled = pd.concat([a1, a2])
    freqs = pooled.value_counts(normalize=True)
    exp_hom = float((freqs.values ** 2).sum())
    return {
        "n_individuals": int(n),
        "obs_hom_pct": round(100 * obs_hom, 2),
        "exp_hom_pct": round(100 * exp_hom, 2),
        # >1 => excess homozygosity (dropout/ref-bias signature); ~1 => HWE-consistent
        "obs_over_exp": round(obs_hom / exp_hom, 2) if exp_hom else None,
        "n_distinct_alleles": int(pooled.nunique()),
    }


def print_diversity_and_freq(df):
    print("## Allelic diversity + top alleles per locus (2-field, whole cohort)\n")
    print("| Gene | Callable individuals | Distinct 2-field alleles | Rare (<1%) share | "
          "Top alleles (freq %) |")
    print("|---|---|---|---|---|")
    freq_tables = {}
    for gene in CLASSICAL_GENES:
        a1, a2 = locus_alleles(df, gene)
        pooled = pd.concat([a1, a2])
        if pooled.empty:
            print(f"| {gene} | 0 | — | — | NONE CALLABLE |")
            continue
        freqs = pooled.value_counts(normalize=True)
        freq_tables[gene] = freqs
        n_distinct = freqs.size
        rare_share = (freqs < RARE_THRESHOLD).sum() / n_distinct
        top = ", ".join(f"{al} ({100*f:.1f})" for al, f in freqs.head(TOP_N).items())
        print(f"| {gene} | {len(a1)} | {n_distinct} | {100*rare_share:.0f}% | {top} |")
    print()
    return freq_tables


def print_hwe(df, anc):
    print("## Hardy-Weinberg proxy — observed vs expected homozygosity (2-field)\n")
    print("Expected = sum(freq^2). **obs/exp near 1.0 = clean; >>1 = excess homozygosity "
          "(possible dropout/reference-bias, or real biology — flag, don't auto-attribute).**\n")
    groups = [("ALL", df)]
    if anc is not None:
        for label, sub in df.groupby("_ancestry"):
            groups.append((label, sub))
    print("| Gene | Group | n | Obs hom % | Exp hom % | obs/exp | Distinct alleles |")
    print("|---|---|---|---|---|---|---|")
    for gene in CLASSICAL_GENES:
        for label, sub in groups:
            a1, a2 = locus_alleles(sub, gene)
            h = hwe_proxy(a1, a2)
            if h is None:
                continue
            print(f"| {gene} | {label} | {h['n_individuals']} | {h['obs_hom_pct']} | "
                  f"{h['exp_hom_pct']} | {h['obs_over_exp']} | {h['n_distinct_alleles']} |")
    print()


def dump_freq_by_ancestry(df, out_dir):
    """Per-locus per-ancestry 2-field allele frequencies -> CSV, for the external
    (AFND / published) concordance check. Aggregate only."""
    rows = []
    for gene in CLASSICAL_GENES:
        for label, sub in df.groupby("_ancestry"):
            a1, a2 = locus_alleles(sub, gene)
            pooled = pd.concat([a1, a2])
            if pooled.empty:
                continue
            freqs = pooled.value_counts(normalize=True)
            n_ind = len(a1)
            for al, f in freqs.items():
                rows.append({"gene": gene, "ancestry": label, "allele": al,
                             "freq": round(float(f), 6), "n_individuals": n_ind})
    out = pd.DataFrame(rows)
    path = os.path.join(out_dir, "allele_freq_by_ancestry.csv")
    out.to_csv(path, index=False)
    print(f"(wrote per-ancestry 2-field allele frequencies to {path} — aggregate only, "
          f"{len(out)} allele-rows — use for the external reference-frequency comparison)\n",
          file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--ancestry-tsv", default=None)
    ap.add_argument("--id-col", default="research_id")
    ap.add_argument("--out-dir", default=os.path.expanduser("~/pipeline_outputs/experiment_a2"))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    df = pd.read_csv(args.tsv, sep="\t", dtype=str)
    print(f"## Experiment A2 — AoU-native callset quality/distribution\n")
    print(f"Loaded {len(df)} research_ids from {os.path.basename(args.tsv)}.\n")

    anc = None
    if args.ancestry_tsv:
        anc = pd.read_csv(args.ancestry_tsv, sep=None, engine="python", dtype=str)
        anc_id = next((c for c in anc.columns if c.lower() in (args.id_col, "research_id", "person_id", "s")), None)
        anc_col = next((c for c in anc.columns if "ancestry" in c.lower()), None)
        if anc_id and anc_col:
            df = df.merge(anc[[anc_id, anc_col]].rename(columns={anc_col: "_ancestry"}),
                          left_on=args.id_col, right_on=anc_id, how="left")
            dist = df["_ancestry"].value_counts(dropna=False)
            print("## Cohort composition by genetic ancestry\n")
            print("| Ancestry | n | % |")
            print("|---|---|---|")
            for label, cnt in dist.items():
                print(f"| {label} | {cnt} | {100*cnt/len(df):.1f}% |")
            print()
        else:
            print(f"WARNING: could not identify ancestry columns in {list(anc.columns)}; "
                  f"skipping ancestry breakdowns.\n", file=sys.stderr)
            anc = None

    print_diversity_and_freq(df)
    print_hwe(df, anc)
    if anc is not None:
        dump_freq_by_ancestry(df, args.out_dir)


if __name__ == "__main__":
    main()
