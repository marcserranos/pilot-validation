#!/usr/bin/env python3
"""Experiment A3: AoU-native callset vs. the FULL catalogued IPD-IMGT/HLA allele space.

Answers: of every allele IPD-IMGT/HLA has ever officially named for our 8 classical genes, how
many were actually observed in AoU's 535,658-person cohort -- and how does that coverage change
by resolution (2-field vs 3-field)? Also attempts a DB-VERSION LOWER BOUND: if AoU's ensemble
successfully called an allele-group that IMGT only added in release X, AoU's underlying reference
must be >= X (a lower bound, not a point estimate -- absence of a newer group doesn't prove an
older DB, since a rare group may just not appear in this cohort).

Reference data (public, NOT participant data -- download once, no need to commit to git, same
pattern as the ~/ref/ hg38 FASTA):
  mkdir -p ~/ref/imgt
  curl -sL -o ~/ref/imgt/Allelelist.txt \
    https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt
  curl -sL -o ~/ref/imgt/Allelelist_history.txt \
    https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist_history.txt   # ~10-20MB, optional

Allelelist.txt schema (confirmed by fetching it directly, 2026-07-18): 6 leading '#' comment
lines (file/date/version/origin/repository/author), then a header row `AlleleID,Allele`, then one
row per officially-named allele, e.g. `HLA00001,A*01:01:01:01`. Current release noted in the
comment header (e.g. "IPD-IMGT/HLA 3.65.0").

Allelelist_history.txt schema is NOT independently confirmed here (file exceeds this dev
environment's fetch size) -- documented by the IMGT/HLA repo as "the official name used in each
release" per allele, i.e. expected to be AlleleID + one column per historical release version.
The version-inference section below is written defensively: it prints exactly what columns it
detects and skips cleanly with a clear message if the assumed shape doesn't hold, rather than
silently producing a wrong answer. Treat any version-inference output as directional, not proof,
until that detection print is eyeballed.

All output is aggregate (allele counts/coverage, optionally per-ancestry frequency) -- no
participant-level rows -- egress-safe under this project's standing caveat.

Usage (via `pixi run -e spechla` for pandas):
  python3 experiment_a3_allele_space_coverage.py \
      --tsv ~/mnt/aou-controlled/v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv \
      --ancestry-tsv ~/mnt/aou-controlled/v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv \
      --allelelist ~/ref/imgt/Allelelist.txt \
      --allelelist-history ~/ref/imgt/Allelelist_history.txt \
      --out-dir ~/pipeline_outputs/experiment_a3
"""
import argparse
import os
import re
import sys

import pandas as pd

CLASSICAL_GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
NULL = {"", "NA", "-", "nan", "None", ".", "na"}
TOP_N_GROUPS = 15  # per-gene groups to print in the compact "tree-ready" table
VERSION_RE = re.compile(r"^\d+\.\d+\.\d+$")


def strip_star(a):
    return a.split("*", 1)[1] if "*" in a else a


def gene_of(allele_str):
    return allele_str.split("*", 1)[0] if "*" in allele_str else None


def trunc(allele, n):
    """Truncate a normalized (post strip_star) allele string to n dot-separated fields.
    Returns None if uncallable or has fewer than n fields."""
    if allele is None:
        return None
    a = str(allele).strip()
    if a in NULL:
        return None
    a = strip_star(a)
    fields = a.split(":")
    if len(fields) < n:
        return None
    return ":".join(fields[:n])


def load_catalogue(path):
    """Returns dict[gene] -> set of full-resolution allele strings (post strip_star)."""
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            if line.startswith("AlleleID"):
                continue  # header
            parts = line.strip().split(",")
            if len(parts) != 2:
                continue
            allele_id, allele = parts
            gene = gene_of(allele)
            if gene in CLASSICAL_GENES:
                rows.append((allele_id, gene, strip_star(allele)))
    df = pd.DataFrame(rows, columns=["AlleleID", "gene", "allele_full"])
    print(f"Loaded catalogue: {len(df)} named alleles across the 8 classical genes "
          f"(file may contain many more across all ~30 HLA/pseudogene loci -- filtered).",
          file=sys.stderr)
    return df


def observed_set(df, gene, n):
    c1, c2 = f"{gene}_1", f"{gene}_2"
    a1 = df[c1].map(lambda x: trunc(x, n))
    a2 = df[c2].map(lambda x: trunc(x, n))
    return set(a1.dropna()) | set(a2.dropna())


def coverage_report(aou_df, cat_df, out_dir):
    print("## Catalogued vs. observed allele-space coverage\n")
    print("Coverage = distinct alleles AoU observed in 535k people, as a fraction of every "
          "allele IPD-IMGT/HLA has ever officially named for that gene, at that resolution. "
          "'Observed, not catalogued' would flag a data-quality problem (a called allele that "
          "doesn't exist in the reference) -- expect this to be ~0.\n")
    print("| Gene | Catalogued (2f) | Observed (2f) | Coverage 2f | Catalogued (3f) | "
          "Observed (3f) | Coverage 3f | Observed-not-catalogued (2f) |")
    print("|---|---|---|---|---|---|---|---|")

    tree_rows = []
    for gene in CLASSICAL_GENES:
        cat_gene = cat_df[cat_df.gene == gene]["allele_full"]
        cat_2f = set(a for a in cat_gene.map(lambda x: trunc(x, 2)) if a)
        cat_3f = set(a for a in cat_gene.map(lambda x: trunc(x, 3)) if a)

        obs_2f = observed_set(aou_df, gene, 2)
        obs_3f = observed_set(aou_df, gene, 3)

        cov2 = len(obs_2f & cat_2f) / len(cat_2f) if cat_2f else None
        cov3 = len(obs_3f & cat_3f) / len(cat_3f) if cat_3f else None
        not_cat = obs_2f - cat_2f

        print(f"| {gene} | {len(cat_2f)} | {len(obs_2f)} | "
              f"{f'{100*cov2:.1f}%' if cov2 is not None else '—'} | "
              f"{len(cat_3f)} | {len(obs_3f)} | "
              f"{f'{100*cov3:.1f}%' if cov3 is not None else '—'} | "
              f"{len(not_cat)}{' <-- CHECK' if not_cat else ''} |")

        for g2 in obs_2f:
            cat_sub3 = {a for a in cat_3f if a.startswith(g2 + ":")}
            obs_sub3 = {a for a in obs_3f if a.startswith(g2 + ":")}
            tree_rows.append({"gene": gene, "group_2field": g2,
                              "catalogued_3field_variants": len(cat_sub3),
                              "observed_3field_variants": len(obs_sub3)})
    print()
    return pd.DataFrame(tree_rows)


def add_frequency_and_ancestry(tree_df, aou_df, anc_df, id_col):
    """Adds n_individuals (carriers of >=1 copy) and, if ancestry available, dominant_ancestry
    per 2-field group. Aggregate counts only."""
    merged = aou_df
    has_anc = anc_df is not None
    if has_anc:
        anc_id = next((c for c in anc_df.columns
                       if c.lower() in (id_col, "research_id", "person_id", "s")), None)
        anc_col = next((c for c in anc_df.columns if "ancestry" in c.lower()), None)
        if anc_id and anc_col:
            merged = aou_df.merge(anc_df[[anc_id, anc_col]].rename(columns={anc_col: "_ancestry"}),
                                  left_on=id_col, right_on=anc_id, how="left")
        else:
            has_anc = False

    counts = []
    for gene in CLASSICAL_GENES:
        c1, c2 = f"{gene}_1", f"{gene}_2"
        t1 = merged[c1].map(lambda x: trunc(x, 2))
        t2 = merged[c2].map(lambda x: trunc(x, 2))
        for g2 in pd.concat([t1, t2]).dropna().unique():
            carries = (t1 == g2) | (t2 == g2)
            row = {"gene": gene, "group_2field": g2, "n_individuals": int(carries.sum())}
            if has_anc:
                sub = merged.loc[carries, "_ancestry"]
                if len(sub):
                    vc = sub.value_counts(normalize=True)
                    row["dominant_ancestry"] = vc.index[0]
                    row["dominant_ancestry_share"] = round(100 * vc.iloc[0], 1)
            counts.append(row)
    freq_df = pd.DataFrame(counts)
    return tree_df.merge(freq_df, on=["gene", "group_2field"], how="left")


def print_tree_summary(tree_df, out_dir):
    tree_df = tree_df.sort_values(["gene", "n_individuals"], ascending=[True, False])
    path = os.path.join(out_dir, "allele_tree_data.csv")
    tree_df.to_csv(path, index=False)
    print(f"## Allele-space tree data (top {TOP_N_GROUPS} groups per gene by carrier count, "
          f"for a sunburst/treemap) -- full data in {path}\n")
    has_anc = "dominant_ancestry" in tree_df.columns
    hdr = "| Gene | 2-field group | Carriers (n) | Catalogued 3f variants | Observed 3f variants |"
    hdr += " Dominant ancestry (share %) |" if has_anc else ""
    print(hdr)
    print("|---|---|---|---|---|" + ("---|" if has_anc else ""))
    for gene, sub in tree_df.groupby("gene", sort=False):
        for _, r in sub.head(TOP_N_GROUPS).iterrows():
            line = (f"| {r['gene']} | {r['group_2field']} | {int(r['n_individuals'])} | "
                    f"{int(r['catalogued_3field_variants'])} | {int(r['observed_3field_variants'])} |")
            if has_anc:
                da = r.get("dominant_ancestry", "—")
                sh = r.get("dominant_ancestry_share", "")
                line += f" {da} ({sh}%) |"
            print(line)
    print()


def version_inference(history_path, cat_df, aou_df):
    print("## DB-version lower bound (exploratory -- see module docstring caveat)\n")
    try:
        with open(history_path) as f:
            header_line = None
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                header_line = line.strip()
                break
    except FileNotFoundError:
        print(f"SKIPPED -- {history_path} not found.\n")
        return

    cols = header_line.split(",")
    version_cols = [c for c in cols if VERSION_RE.match(c.strip())]
    print(f"Detected {len(cols)} columns; {len(version_cols)} look like release-version strings "
          f"(e.g. {version_cols[:3]}{'...' if len(version_cols) > 3 else ''}).\n")
    if len(version_cols) < 2:
        print("SKIPPED -- fewer than 2 version-looking columns detected; the assumed schema "
              "(one column per IMGT release) doesn't hold for this file as fetched. Needs a "
              "manual look at the real header before this section can run:\n")
        print(f"```\n{header_line[:500]}\n```\n")
        return

    id_col = cols[0].strip()
    print(f"Using '{id_col}' as the row key (assumed = AlleleID).\n")

    def vtuple(v):
        return tuple(int(x) for x in v.split("."))
    version_cols_sorted = sorted(version_cols, key=vtuple)

    hist = pd.read_csv(history_path, comment="#", dtype=str)
    if id_col not in hist.columns:
        print(f"SKIPPED -- expected row-key column '{id_col}' not found after full parse; "
              f"actual columns: {list(hist.columns)[:10]}...\n")
        return

    cat_ids = set(cat_df["AlleleID"])
    hist = hist[hist[id_col].isin(cat_ids)]
    print(f"Matched {len(hist)} / {len(cat_ids)} classical-gene AlleleIDs into the history file.\n")

    def earliest_release(row):
        for v in version_cols_sorted:
            val = row.get(v)
            if pd.notna(val) and str(val).strip() not in NULL:
                return v
        return None

    hist["_first_release"] = hist.apply(earliest_release, axis=1)
    id_to_release = dict(zip(hist[id_col], hist["_first_release"]))
    id_to_gene_allele = dict(zip(cat_df["AlleleID"], zip(cat_df["gene"], cat_df["allele_full"])))

    aou_observed_2f = set()
    for gene in CLASSICAL_GENES:
        aou_observed_2f |= {(gene, a) for a in observed_set(aou_df, gene, 2)}

    allele_to_release = {}
    for aid, rel in id_to_release.items():
        ga = id_to_gene_allele.get(aid)
        if ga is None or rel is None:
            continue
        gene, full = ga
        g2 = trunc(full, 2)
        if g2 and (gene, g2) not in allele_to_release:
            allele_to_release[(gene, g2)] = rel  # first AlleleID seen wins == earliest by construction

    newest_used = None
    for key in aou_observed_2f:
        rel = allele_to_release.get(key)
        if rel and (newest_used is None or vtuple(rel) > vtuple(newest_used)):
            newest_used = rel

    print(f"**Newest IMGT release that introduced a 2-field allele-group AoU actually called: "
          f"{newest_used or 'not determined'}**\n")
    print("This is a LOWER BOUND on AoU's underlying reference DB version -- it must be at "
          "least this new to have named this group. It is NOT an upper bound: AoU's DB could be "
          "newer and simply never have called an even-more-recent group in this cohort. Compare "
          "against this project's own tool DB versions for context: SpecHLA 3.38.0, SpecImmune "
          "3.64.0 (DECISIONS.md).\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--ancestry-tsv", default=None)
    ap.add_argument("--id-col", default="research_id")
    ap.add_argument("--allelelist", required=True)
    ap.add_argument("--allelelist-history", default=None)
    ap.add_argument("--out-dir", default=os.path.expanduser("~/pipeline_outputs/experiment_a3"))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("## Experiment A3 -- AoU-native callset vs. full IPD-IMGT/HLA allele space\n")

    aou_df = pd.read_csv(args.tsv, sep="\t", dtype=str)
    print(f"Loaded {len(aou_df)} research_ids.\n", file=sys.stderr)

    cat_df = load_catalogue(args.allelelist)

    anc_df = None
    if args.ancestry_tsv:
        anc_df = pd.read_csv(args.ancestry_tsv, sep=None, engine="python", dtype=str)

    tree_df = coverage_report(aou_df, cat_df, args.out_dir)
    tree_df = add_frequency_and_ancestry(tree_df, aou_df, anc_df, args.id_col)
    print_tree_summary(tree_df, args.out_dir)

    if args.allelelist_history:
        version_inference(args.allelelist_history, cat_df, aou_df)
    else:
        print("## DB-version lower bound\n\nSKIPPED -- no --allelelist-history path given.\n")


if __name__ == "__main__":
    main()
