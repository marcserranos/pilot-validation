#!/usr/bin/env python3
"""The main-text callout list: novel HLA alleles that are common in a population and absent from
IPD-IMGT/HLA.

Supervisor ask A6 (Cole, 2026-09-17 call):

    "I do also want to know if there's any crazy alleles that are, like, super common and they're
     in some ancestry, and we've just been missing them from the database. I want a really great
     list of those, because those are things we can literally call out in main text. Like, we found
     these seven examples of this."

This script runs **on committed aggregates only** -- `reports/hla_popgen/24_novelty_by_field/` and
`reports/hla_popgen/27_allele_space_coverage/` -- so it needs no VM and no participant data.

Selection, and why the disclosure rule and the science point the same way
------------------------------------------------------------------------
S01's cluster table writes every count below 20 as `<20` (AoU small-cell rule). That looks like a
constraint but here it is a feature: an allele we can quote a number for is, by construction, an
allele carried by at least 20 unrelated people. Those are exactly the alleles worth naming in main
text, and they are the only ones whose frequencies we could publish anyway.

So the tiers are:

  * **Tier 1 (main text).** Clean, recurrent in >= 20 unrelated people, count disclosable.
  * **Tier 2 (IMGT submission).** Clean and recurrent in >= 2 unrelated people, count censored.
    Too many to name, but they are the submission queue.
  * Everything else is reported only as a total.

`clean` means no member haplotype carried an artifact flag (`partial_CDS` / `inframe_stop`) and the
cluster is not homopolymer-indel-only -- S01's definition, unchanged.

Ancestry concentration is computed on the **strict** ancestry assignment (probability >= 0.9), and
carrier rates use the per-ancestry haplotype denominators from script 27. A rate is only printed
when both numerator and denominator are disclosable; otherwise the concentration share is given
without a rate.

Usage:

    python3 scripts/hla_popgen/32_novel_allele_callouts.py

Usage (local fixtures): scripts/hla_popgen/tests/test_novel_allele_callouts.py
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_THIS_DIR))
DEFAULT_REPORTS = os.path.join(_REPO_ROOT, "reports", "hla_popgen")
DEFAULT_CLUSTERS = os.path.join(DEFAULT_REPORTS, "24_novelty_by_field",
                                "protein_level_novel_clusters.tsv")
DEFAULT_COVERAGE = os.path.join(DEFAULT_REPORTS, "27_allele_space_coverage",
                                "coverage_by_gene_ancestry.tsv")
DEFAULT_NONCODING = os.path.join(DEFAULT_REPORTS, "25_noncoding_novelty_paf",
                                 "noncoding_signature_clusters.tsv")
DEFAULT_IMGT_TODAY = os.path.join(DEFAULT_REPORTS, "27_allele_space_coverage",
                                  "imgt_coverage_today.tsv")
DEFAULT_OUT_DIR = os.path.join(DEFAULT_REPORTS, "32_novel_callouts")

ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
ANCESTRY_COLORS = {"AFR": "#E69F00", "AMR": "#56B4E9", "EAS": "#009E73",
                   "EUR": "#0072B2", "MID": "#D55E00", "SAS": "#CC79A7"}
SUPPRESS_BELOW = 20
CENSORED = "<%d" % SUPPRESS_BELOW


def parse_count(v):
    """A suppressed cell -> NaN, not 0 and not 20.

    This is the one function in this script that can quietly corrupt everything. `<20` means
    "somewhere in 1..19, we are not allowed to say". Reading it as 0 would erase real carriers;
    reading it as 20 would invent them and let a censored allele outrank a disclosable one.
    NaN forces every downstream comparison to be explicit about censoring.
    """
    if v is None:
        return float("nan")
    s = str(v).strip()
    if s == "" or s.upper() in {"NA", "NAN", "NONE"}:
        return float("nan")
    if s.startswith("<"):
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def is_censored(v):
    return str(v).strip().startswith("<")


def ancestry_profile(row, prefix="n_persons_strict_"):
    """-> {ancestry: count} over disclosable cells, plus the set of censored ancestries."""
    counts, censored = {}, []
    for a in ANCESTRY_ORDER:
        raw = row.get(prefix + a)
        if is_censored(raw):
            censored.append(a)
            continue
        c = parse_count(raw)
        if not np.isnan(c) and c > 0:
            counts[a] = c
    return counts, censored


def concentration(counts):
    """Share of disclosable carriers in the single most common ancestry, and that ancestry.

    Reported instead of a formal test because the denominators differ hugely between ancestries
    and the censored cells are not missing at random -- a chi-square here would be theatre.
    """
    tot = sum(counts.values())
    if tot <= 0:
        return float("nan"), None
    top = max(counts.items(), key=lambda kv: kv[1])
    return top[1] / tot, top[0]


def load_denominators(path):
    """gene -> ancestry -> people typed (haplotypes / 2), from script 27's coverage table.

    Haplotypes are converted to people because the cluster table counts people. This assumes
    both haplotypes of a person were typed for that gene, which is not exactly true where an
    assembly is fragmented, so the resulting rates are slight UNDER-estimates. Stated, not hidden.
    """
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    if "resolution" in df.columns:
        # protein-level rows only; the CDS-level rows would double-count the same haplotypes
        res = df["resolution"].astype(str)
        df = df[res.str.contains("protein", case=False)] if res.str.contains(
            "protein", case=False).any() else df
    out = {}
    for _, r in df.iterrows():
        g = str(r.get("gene", "")).replace("HLA-", "")
        a = str(r.get("ancestry", "")).upper()
        n = parse_count(r.get("n_haplotypes"))
        if g and a in ANCESTRY_ORDER and not np.isnan(n):
            out.setdefault(g, {})[a] = n / 2.0
    return out


def catalogue_gap_by_gene(path):
    """gene -> % of this cohort's haplotypes whose protein is NOT yet in IPD-IMGT.

    This is the control the callout list needs. A gene can top the novel-allele table for two
    opposite reasons: it is genuinely diverse in under-sampled populations, or IPD-IMGT simply
    never catalogued it. Script 27 measured the second quantity directly, so we print it next to
    every callout instead of arguing about it. If the genes with the most novel alleles are also
    the genes with the largest catalogue gap, the two independent estimates agree -- which is
    evidence, not a confound.
    """
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    if "resolution" in df.columns:
        res = df["resolution"].astype(str)
        if res.str.contains("protein", case=False).any():
            df = df[res.str.contains("protein", case=False)]
    df["n"] = pd.to_numeric(df.get("n_haplotypes"), errors="coerce")
    df["f"] = pd.to_numeric(df.get("frac_in_imgt_today"), errors="coerce")
    df = df.dropna(subset=["n", "f"])
    out = {}
    for g, sub in df.groupby("gene"):
        tot = sub["n"].sum()
        if tot > 0:
            out[str(g)] = round(100.0 * (1.0 - (sub["f"] * sub["n"]).sum() / tot), 3)
    return out


def build_tiers(clusters, denominators, catalogue_gap=None):
    df = clusters.copy()
    df["gene_b"] = df["gene"].astype(str).str.replace("HLA-", "", regex=False)
    df["n_unrel_clean"] = df["n_persons_unrelated_clean"].map(parse_count)
    df["n_unrel_clean_censored"] = df["n_persons_unrelated_clean"].map(is_censored)
    df["recurrent"] = df["clean_recurrent_unrelated"].astype(str).str.lower().isin(
        {"true", "1", "yes"})

    rows = []
    for _, r in df.iterrows():
        counts, censored_anc = ancestry_profile(r)
        share, top_anc = concentration(counts)
        rate = float("nan")
        if top_anc is not None:
            den = denominators.get(r["gene_b"], {}).get(top_anc)
            if den and den > 0:
                rate = 1000.0 * counts[top_anc] / den
        rows.append({
            "cluster_id": r.get("cluster_id"),
            "cluster_type": r.get("cluster_type"),
            "gene": r.get("gene"),
            "gene_class": r.get("gene_class"),
            "nearest_known_allele": r.get("nearest_known_allele"),
            "nearest_known_protein": r.get("nearest_known_protein"),
            "n_aa_diffs": r.get("n_aa_diffs"),
            "aa_diffs": r.get("aa_diffs"),
            "in_groove": r.get("in_groove"),
            "n_groove_diffs": r.get("n_groove_diffs"),
            "n_persons_unrelated_clean": r.get("n_persons_unrelated_clean"),
            "n_unrel_clean": r["n_unrel_clean"],
            "recurrent_clean": r["recurrent"],
            "top_ancestry": top_anc,
            "top_ancestry_share": round(share, 3) if not np.isnan(share) else None,
            "carriers_per_1000_in_top_ancestry": round(rate, 2) if not np.isnan(rate) else None,
            "ancestries_disclosable": ",".join("%s:%d" % (a, int(c))
                                               for a, c in sorted(counts.items())) or "none",
            "ancestries_censored": ",".join(censored_anc) or "none",
            "pct_gene_haplotypes_not_in_imgt": (catalogue_gap or {}).get(str(r.get("gene"))),
        })
    out = pd.DataFrame(rows)

    tier1 = out[(out["n_unrel_clean"] >= SUPPRESS_BELOW) & out["recurrent_clean"]].copy()
    tier1 = tier1.sort_values("n_unrel_clean", ascending=False)
    tier2 = out[out["recurrent_clean"] & ~(out["n_unrel_clean"] >= SUPPRESS_BELOW)].copy()
    return out, tier1, tier2


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
def fig_callouts(tier1, path, top_n=16):
    """One row per named allele: carriers by ancestry, so the concentration is visible at a glance.

    This is the figure that makes Cole's sentence ('we found these N examples') checkable.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = tier1.head(top_n).copy()
    if d.empty:
        return
    d = d.iloc[::-1]
    fig, ax = plt.subplots(figsize=(9.5, 0.45 * len(d) + 2.0))
    labels = []
    for yi, (_, r) in enumerate(d.iterrows()):
        left = 0.0
        parts = dict(kv.split(":") for kv in r["ancestries_disclosable"].split(",")
                     if ":" in kv)
        for a in ANCESTRY_ORDER:
            if a not in parts:
                continue
            w = float(parts[a])
            ax.barh(yi, w, left=left, color=ANCESTRY_COLORS[a], height=0.68,
                    label=a if yi == len(d) - 1 else None)
            left += w
        if r["ancestries_censored"] != "none":
            ax.barh(yi, 0.0, left=left, color="none")
            ax.text(left + 0.6, yi, "+", va="center", fontsize=8, color="#666666")
        near = str(r["nearest_known_protein"])
        nd = r["n_aa_diffs"]
        labels.append("%s  (%s aa from %s)" % (r["gene"], nd, near))
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("unrelated carriers (strict-ancestry, disclosable counts only)")
    handles, lab = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, lab, frameon=False, ncol=len(lab), fontsize=8, loc="lower right")
    ax.set_title("Novel HLA proteins carried by >= %d unrelated people and absent from "
                 "IPD-IMGT/HLA\n'+' marks an allele with additional carriers in an ancestry whose "
                 "count is below the disclosure threshold" % SUPPRESS_BELOW, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def write_readme(path, all_df, tier1, tier2, noncoding_note):
    L = []
    L.append("# 32 — main-text callout list: common novel HLA alleles absent from IPD-IMGT/HLA\n")
    L.append("*Supervisor ask A6 (Cole, 2026-09-17): \"a really great list of those, because those "
             "are things we can literally call out in main text\".*\n")
    L.append("Built from committed aggregates only (scripts 24, 25, 27). No VM, no participant "
             "data.\n")

    L.append("\n## Tier 1 — nameable in main text (%d alleles)\n" % len(tier1))
    L.append("Clean (no artifact flag, not homopolymer-only), recurrent, and carried by at least "
             "%d unrelated people — which is also exactly the threshold above which we are "
             "permitted to publish the count.\n" % SUPPRESS_BELOW)
    if not tier1.empty:
        cols = ["gene", "nearest_known_protein", "n_aa_diffs", "in_groove",
                "n_persons_unrelated_clean", "top_ancestry", "top_ancestry_share",
                "carriers_per_1000_in_top_ancestry", "pct_gene_haplotypes_not_in_imgt",
                "ancestries_disclosable"]
        prot = tier1[tier1["cluster_type"] == "novel_protein"]
        syn = tier1[tier1["cluster_type"] != "novel_protein"]
        L.append("### 1a. Novel **proteins** (%d)\n" % len(prot))
        L.append(prot[cols].to_markdown(index=False) + "\n")
        if not syn.empty:
            L.append("\n### 1b. Novel synonymous CDS — same protein, new coding sequence (%d)\n"
                     % len(syn))
            L.append("These are *not* new proteins (`n_aa_diffs` is 0 by construction) and must "
                     "not be counted in a 'new HLA proteins' sentence. They are still real, "
                     "uncatalogued coding sequences.\n")
            L.append(syn[cols].to_markdown(index=False) + "\n")
        L.append("\n**The reference-depth control.** `pct_gene_haplotypes_not_in_imgt` is the "
                 "share of this cohort's haplotypes at that gene whose protein is already absent "
                 "from IPD-IMGT/HLA, measured independently in script 27. Read every callout "
                 "against it: a gene with a large gap was under-catalogued to begin with. The "
                 "classical genes sit at 0.1-0.2%; TAP1/TAP2 are the extreme at ~7-8%.\n")
        L.append("\n`in_groove` is whether any of the amino-acid differences fall in the "
                 "peptide-binding groove exons — a groove difference is the one most likely to "
                 "change which peptides the molecule presents, so those are the interesting ones "
                 "biologically, not just numerically.\n")

    L.append("\n## Tier 2 — IMGT submission queue (%d alleles)\n" % len(tier2))
    L.append("Clean and seen in at least two unrelated people, but with fewer than %d carriers, so "
             "their counts stay censored. Too many to name individually; this is the list to "
             "submit, not to quote.\n" % SUPPRESS_BELOW)
    if not tier2.empty:
        by_gene = (tier2.groupby("gene").size().reset_index(name="n_alleles")
                   .sort_values("n_alleles", ascending=False))
        L.append(by_gene.to_markdown(index=False) + "\n")

    L.append("\n## How to read the ancestry columns\n")
    L.append("`top_ancestry_share` is the fraction of **disclosable** carriers in the single most "
             "common ancestry. `ancestries_censored` lists ancestries where this allele has "
             "carriers but fewer than %d, so the true share is at least as concentrated as shown "
             "and possibly less. Rates per 1,000 use script 27's per-ancestry haplotype counts "
             "divided by two; where an assembly was fragmented the person was still typed on one "
             "haplotype, so these rates are slight under-estimates.\n" % SUPPRESS_BELOW)

    L.append("\n## What this list is not\n")
    L.append("- **Not validated orthogonally.** Every allele here rests on the long-read assembly "
             "alone. Cole's own suggestion — realigning the short reads, which every one of these "
             "people also has, to the long-read assembly — is the confirmation step, and he "
             "explicitly parked it for now.\n")
    L.append("- **Not an IMGT submission yet.** Submission needs full-length genomic sequence and "
             "a documented typing method; what we have is the candidate set.\n")
    L.append("- **Not a frequency estimate for the general population.** These are All of Us "
             "participants with long-read data, which is not a random sample of anything.\n")
    if noncoding_note:
        L.append("\n## Non-coding novelty\n")
        L.append(noncoding_note + "\n")

    L.append("\n## Files\n")
    L.append("- `main_text_candidates.tsv` — tier 1.\n")
    L.append("- `imgt_submission_candidates.tsv` — tier 2.\n")
    L.append("- `all_novel_clusters_annotated.tsv` — every cluster with the ancestry annotation.\n")
    L.append("- `fig_novel_callouts.png`.\n")
    with open(path, "w") as fh:
        fh.write("\n".join(L))


def run(args):
    os.makedirs(args.out_dir, exist_ok=True)
    if not os.path.exists(args.clusters):
        sys.exit("FATAL: cluster table not found: %s (run script 24 first)" % args.clusters)
    clusters = pd.read_csv(args.clusters, sep="\t", dtype=str)
    denominators = load_denominators(args.coverage)
    if not denominators:
        print("[32] WARNING: no per-ancestry denominators (%s); carrier rates will be blank"
              % args.coverage, flush=True)

    catalogue_gap = catalogue_gap_by_gene(args.imgt_today)
    if not catalogue_gap:
        print("[32] WARNING: no catalogue-gap table (%s); the reference-depth control will be "
              "missing from the callout list" % args.imgt_today, flush=True)
    all_df, tier1, tier2 = build_tiers(clusters, denominators, catalogue_gap)

    noncoding_note = ""
    if os.path.exists(args.noncoding):
        nc = pd.read_csv(args.noncoding, sep="\t", dtype=str)
        noncoding_note = ("Script 25 resolved the pooled non-coding clusters into %d distinct "
                          "genomic signatures. Those are catalogue gaps in IPD-IMGT's *genomic* "
                          "sequences, not new proteins, and they are reported separately in "
                          "`reports/hla_popgen/25_noncoding_novelty_paf/`." % len(nc))

    all_df.to_csv(os.path.join(args.out_dir, "all_novel_clusters_annotated.tsv"), sep="\t",
                  index=False, na_rep="NA")
    tier1.to_csv(os.path.join(args.out_dir, "main_text_candidates.tsv"), sep="\t", index=False,
                 na_rep="NA")
    tier2.to_csv(os.path.join(args.out_dir, "imgt_submission_candidates.tsv"), sep="\t",
                 index=False, na_rep="NA")
    fig_callouts(tier1, os.path.join(args.out_dir, "fig_novel_callouts.png"))
    write_readme(os.path.join(args.out_dir, "README.md"), all_df, tier1, tier2, noncoding_note)

    with open(os.path.join(args.out_dir, "summary.json"), "w") as fh:
        json.dump({"n_clusters_total": int(len(all_df)),
                   "n_tier1_main_text": int(len(tier1)),
                   "n_tier2_submission": int(len(tier2)),
                   "tier1_genes": sorted(set(tier1["gene"].astype(str)))}, fh, indent=2)
    print("[32] tier 1 (main text): %d; tier 2 (submission): %d -> %s"
          % (len(tier1), len(tier2), args.out_dir), flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clusters", default=DEFAULT_CLUSTERS)
    ap.add_argument("--coverage", default=DEFAULT_COVERAGE)
    ap.add_argument("--noncoding", default=DEFAULT_NONCODING)
    ap.add_argument("--imgt-today", default=DEFAULT_IMGT_TODAY)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    run(ap.parse_args(argv))


if __name__ == "__main__":
    main()
