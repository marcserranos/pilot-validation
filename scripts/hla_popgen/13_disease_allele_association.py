#!/usr/bin/env python3
"""Does actually carrying a disease-associated HLA allele predict actually being diagnosed with
that disease, in this cohort -- properly tested, not eyeballed off a UMAP.

## Why this exists (2026-09 chat discussion)

`08_embedding_compare.py --real-disease-labels` showed real EHR-confirmed diagnoses do NOT
visually cluster in allele-driven UMAP space, unlike the allele-carrier-proxy coloring (which
clusters close to tautologically, since it's colored by the same identity the embedding is built
from). A scatter plot isn't sensitive enough to detect a real but comparatively weak, allele-
specific association riding on top of a much stronger population-structure signal. This script
tests it properly instead: for each of a small set of allele<->diagnosis pairs where both a
carrier-status definition (08's DISEASE_ALLELE_GROUPS) and an EHR diagnosis definition
(12_disease_phenotypes.py's HLA_LINKED) already exist, build the 2x2 (diagnosed x carrier) table
and test association two ways:
  1. Raw (unstratified) Fisher's exact test + odds ratio.
  2. Ancestry-adjusted Cochran-Mantel-Haenszel test -- stratifies by the 6 continental ancestry
     groups before combining, since both allele frequency AND diagnosis/coding rates vary
     independently by ancestry (confirmed in this project's own ancestry gradient work), so a raw
     unstratified OR can be confounded in either direction (Simpson's paradox risk).

## Allele<->diagnosis pairing (8 of the original 12 allele groups; 4 have no diagnosis-side
## definition -- drug-hypersensitivity/rare phenotypes 12_disease_phenotypes.py's HLA_LINKED list
## doesn't cover, not tested here rather than paired with a guessed definition)

| allele proxy (08's group)          | EHR diagnosis (12's HLA_LINKED)             | note |
|-------------------------------------|----------------------------------------------|------|
| B*27 (spondyloarthritis)            | Ankylosing spondylitis (HLA-B27)              | exact |
| DQA1*05:01 / DQ2 (celiac)           | Celiac disease (HLA-DQ2/DQ8)                  | exact |
| DQB1*03:02 / DQ8 (T1D)              | Type 1 diabetes (HLA-DR3/DR4)                 | proxy is DQ8; diagnosis's classic association is DR3/DR4 (strong LD, not identical -- flagged) |
| C*06:02 (psoriasis)                 | Psoriasis (HLA-Cw6)                           | exact |
| DRB1*04 (rheumatoid arthritis)      | Rheumatoid arthritis (HLA shared epitope)     | DRB1*04 alleles ARE the shared-epitope carriers |
| DRB1*15:01 (multiple sclerosis)     | Multiple sclerosis (HLA-DRB1*15:01)           | exact |
| DQB1*06:02 (narcolepsy)             | Narcolepsy (HLA-DQB1*06:02)                   | exact |
| B*51 (Behcet's disease)             | Behcet disease (HLA-B51)                      | exact |

## Compliance

Everything read here (mat_full, ancestry) is already aggregate-derived or lives in the VM-local
embedding cache; the per-person diagnosis join happens VM-side and this script's OUTPUT is
counts/OR/p-values only (a 2x2-per-stratum table has 4x6=24 numbers per disease, all counts >=0,
no person_id) -- safe to commit, unlike either of the two per-person TSVs it reads.

Usage (VM): python3 scripts/hla_popgen/13_disease_allele_association.py
"""
import argparse
import os
import pickle
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, chi2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _viz_common as vc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

# (allele-group label from 08_embedding_compare.py's DISEASE_ALLELE_GROUPS,
#  diagnosis label from 12_disease_phenotypes.py's HLA_LINKED, short display name)
PAIRS = [
    ("B*27 (spondyloarthritis)", "Ankylosing spondylitis (HLA-B27)", "B*27 / ank. spondylitis"),
    ("DQA1*05:01 / DQ2 (celiac)", "Celiac disease (HLA-DQ2/DQ8)", "DQ2 / celiac"),
    ("DQB1*03:02 / DQ8 (T1D)", "Type 1 diabetes (HLA-DR3/DR4)", "DQ8 / T1D"),
    ("C*06:02 (psoriasis)", "Psoriasis (HLA-Cw6)", "Cw6 / psoriasis"),
    ("DRB1*04 (rheumatoid arthritis)", "Rheumatoid arthritis (HLA shared epitope)", "DRB1*04 / RA"),
    ("DRB1*15:01 (multiple sclerosis)", "Multiple sclerosis (HLA-DRB1*15:01)", "DRB1*15:01 / MS"),
    ("DQB1*06:02 (narcolepsy)", "Narcolepsy (HLA-DQB1*06:02)", "DQB1*06:02 / narcolepsy"),
    ("B*51 (Behcet's disease)", "Behcet disease (HLA-B51)", "B*51 / Behcet"),
]


def carrier_columns(mat, gene, allele_group):
    prefix = f"{vc.gene_display(gene)}*{allele_group.split('*', 1)[1]}"
    return [c for c in mat.columns if c.split(":", 1)[-1].startswith(prefix)]


def haldane_or_ci(a, b, c, d, alpha=0.05):
    """Odds ratio + Wald CI on log(OR), Haldane-Anscombe 0.5 correction applied uniformly when any
    cell is 0 (standard fix for the otherwise-undefined/infinite OR and CI)."""
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    orr = (a * d) / (b * c)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    z = 1.959963984540054
    lo, hi = np.exp(np.log(orr) - z * se), np.exp(np.log(orr) + z * se)
    return orr, lo, hi


def mantel_haenszel(tables):
    """tables: list of (a,b,c,d) 2x2 counts, one per stratum (ancestry). Returns (OR_MH, chi2_stat,
    p_value, n_strata_used) -- strata with n<2 or a zero margin are dropped (standard practice;
    contribute no information and can destabilize the variance term)."""
    num_or, den_or = 0.0, 0.0
    sum_a, sum_e, sum_v = 0.0, 0.0, 0.0
    n_used = 0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n < 2:
            continue
        row1, row2 = a + b, c + d
        col1, col2 = a + c, b + d
        if row1 == 0 or row2 == 0 or col1 == 0 or col2 == 0:
            continue
        n_used += 1
        num_or += (a * d) / n
        den_or += (b * c) / n
        sum_a += a
        sum_e += (row1 * col1) / n
        sum_v += (row1 * row2 * col1 * col2) / (n ** 2 * (n - 1)) if n > 1 else 0.0
    if den_or == 0 or sum_v == 0:
        return np.nan, np.nan, np.nan, n_used
    or_mh = num_or / den_or
    stat = (abs(sum_a - sum_e) - 0.5) ** 2 / sum_v  # continuity-corrected CMH chi-square, df=1
    p = float(1 - chi2.cdf(stat, df=1))
    return or_mh, stat, p, n_used


def build_2x2(carrier_bool, diag_bool):
    a = int((carrier_bool & diag_bool).sum())        # carrier, diagnosed
    b = int((carrier_bool & ~diag_bool).sum())        # carrier, not diagnosed
    c = int((~carrier_bool & diag_bool).sum())        # non-carrier, diagnosed
    d = int((~carrier_bool & ~diag_bool).sum())       # non-carrier, not diagnosed
    return a, b, c, d


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--embedding-cache",
                    default=os.path.expanduser("~/repos/pilot-validation/reports/hla_popgen/"
                                               "08_embedding_compare/lr/embedding_cache.pkl"))
    ap.add_argument("--wide-diagnosis-labels",
                    default=os.path.expanduser("~/pipeline_outputs/rnaseq/pheno/"
                                               "person_disease_labels_wide.tsv"))
    ap.add_argument("--out-dir", default=None,
                    help="Default: reports/hla_popgen/13_disease_allele_association/")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(os.path.dirname(here))
    out_dir = args.out_dir or os.path.join(repo_root, "reports", "hla_popgen",
                                           "13_disease_allele_association")
    os.makedirs(out_dir, exist_ok=True)

    print(f"Loading {args.embedding_cache!r} ...", file=sys.stderr)
    with open(os.path.expanduser(args.embedding_cache), "rb") as f:
        cache = pickle.load(f)
    mat_full = pd.DataFrame(cache["mat_full_values"], index=cache["mat_full_index"],
                            columns=cache["mat_full_columns"])
    ancestry_by_person = cache["ancestry_by_person"]

    print(f"Loading {args.wide_diagnosis_labels!r} ...", file=sys.stderr)
    diag = pd.read_csv(os.path.expanduser(args.wide_diagnosis_labels), sep="\t",
                       dtype={"person_id": str}).set_index("person_id")
    diag = diag.reindex(mat_full.index).fillna(False)

    ancestry = pd.Series({pid: ancestry_by_person.get(pid) for pid in mat_full.index})
    strata = [a for a in vc.ANCESTRY_ORDER]

    allele_groups = _import_disease_groups()
    rows = []
    for allele_label, disease_label, short_name in PAIRS:
        gene = next(g for lbl, g, _grp in allele_groups if lbl == allele_label)
        cols = carrier_columns(mat_full, gene, allele_label.split(" ")[0])
        carrier = (mat_full[cols].sum(axis=1) > 0) if cols else pd.Series(False, index=mat_full.index)
        diagnosed = diag[disease_label].astype(bool) if disease_label in diag.columns else \
            pd.Series(False, index=mat_full.index)

        a, b, c, d = build_2x2(carrier, diagnosed)
        or_raw, or_lo, or_hi = haldane_or_ci(a, b, c, d)
        _, p_fisher = fisher_exact([[a, b], [c, d]])

        strat_tables = []
        for anc in strata:
            mask = ancestry == anc
            if mask.sum() == 0:
                continue
            strat_tables.append(build_2x2(carrier[mask.index[mask]], diagnosed[mask.index[mask]]))
        or_mh, chi2_stat, p_cmh, n_strata = mantel_haenszel(strat_tables)

        rows.append({
            "pair": short_name, "allele_group": allele_label, "diagnosis": disease_label,
            "n_carriers": a + b, "n_diagnosed": a + c, "n_both": a, "n_cohort": a + b + c + d,
            "OR_raw": round(or_raw, 2), "OR_raw_lo": round(or_lo, 2), "OR_raw_hi": round(or_hi, 2),
            "p_fisher_raw": p_fisher, "OR_ancestry_adjusted": round(or_mh, 2) if or_mh == or_mh
                else None, "p_CMH_ancestry_adjusted": p_cmh, "n_ancestry_strata_used": n_strata,
        })

    df = pd.DataFrame(rows)
    n_tests = len(df)
    bonf = 0.05 / n_tests
    df["significant_raw_p05"] = df["p_fisher_raw"] < 0.05
    df["significant_bonferroni"] = df["p_fisher_raw"] < bonf
    df["significant_cmh_p05"] = df["p_CMH_ancestry_adjusted"] < 0.05
    df["significant_cmh_bonferroni"] = df["p_CMH_ancestry_adjusted"] < bonf

    csv_path = os.path.join(out_dir, "disease_allele_association.tsv")
    df.to_csv(csv_path, sep="\t", index=False)
    print(f"Wrote {csv_path!r}", file=sys.stderr)

    # Forest plot: log-scale OR with 95% CI, raw vs ancestry-adjusted side by side.
    fig, ax = plt.subplots(figsize=(9, 0.6 * len(df) + 1.5))
    y = np.arange(len(df))
    ax.errorbar(df["OR_raw"], y + 0.15,
               xerr=[df["OR_raw"] - df["OR_raw_lo"], df["OR_raw_hi"] - df["OR_raw"]],
               fmt="o", color="#999999", label="raw (unstratified)", capsize=3)
    ax.scatter(df["OR_ancestry_adjusted"], y - 0.15, color="#1f7a5c", marker="D", s=40,
              label="ancestry-adjusted (CMH)", zorder=3)
    ax.axvline(1, color="black", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(df["pair"])
    ax.set_xlabel("Odds ratio (log scale) -- carrier vs non-carrier, diagnosed vs not")
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    for yi, row in zip(y, df.itertuples()):
        sig = "**" if row.p_fisher_raw < bonf else ("*" if row.p_fisher_raw < 0.05 else "")
        ax.annotate(f"n={row.n_both}/{row.n_diagnosed} diagnosed carry it{sig}",
                   xy=(max(row.OR_raw_hi, row.OR_ancestry_adjusted or row.OR_raw_hi), yi),
                   xytext=(6, 0), textcoords="offset points", fontsize=7.5, va="center")
    fig.suptitle("Does carrying the risk allele predict the real diagnosis? "
                f"(* p<0.05, ** p<Bonferroni {bonf:.4f}, n={n_tests} tests)", fontsize=10)
    fig.tight_layout()
    fig_path = os.path.join(out_dir, "disease_allele_association_forest.png")
    fig.savefig(fig_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {fig_path!r}", file=sys.stderr)

    lines = ["# HLA-disease-allele association: real diagnosis vs real carrier status\n",
             "Methodology: see this script's module docstring for the pairing rationale and the "
             "Mantel-Haenszel ancestry-adjustment formula. All numbers below are counts/statistics "
             "only -- no person_id, safe to commit (unlike the two per-person files this script "
             f"reads).\n\n**{n_tests} tests, Bonferroni threshold = {bonf:.4f}.**\n"]
    lines.append("\n| pair | n_cohort | n_carriers | n_diagnosed | n_both | OR (raw, 95% CI) | "
                 "p (Fisher) | OR (ancestry-adj, CMH) | p (CMH) | sig? |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for r in df.itertuples():
        sig = "Bonferroni" if r.significant_bonferroni else ("p<0.05" if r.significant_raw_p05
                                                              else "no")
        lines.append(f"| {r.pair} | {r.n_cohort} | {r.n_carriers} | {r.n_diagnosed} | {r.n_both} "
                     f"| {r.OR_raw} ({r.OR_raw_lo}-{r.OR_raw_hi}) | {r.p_fisher_raw:.2e} | "
                     f"{r.OR_ancestry_adjusted} | {r.p_CMH_ancestry_adjusted:.2e} | {sig} |")
    lines.append(f"\nFigure: `{fig_path}`\n")
    md_path = os.path.join(out_dir, "disease_allele_association_report.md")
    vc.write_report(md_path, lines)
    print(f"Wrote {md_path!r}", file=sys.stderr)
    print(df[["pair", "OR_raw", "p_fisher_raw", "OR_ancestry_adjusted",
             "p_CMH_ancestry_adjusted"]].to_string(index=False), file=sys.stderr)


def _import_disease_groups():
    """Reuses 08_embedding_compare.py's DISEASE_ALLELE_GROUPS list (label, gene, allele_group) as
    a sibling module rather than duplicating it."""
    import importlib.util
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location(
        "hla_popgen_embedding_compare", os.path.join(here, "08_embedding_compare.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.DISEASE_ALLELE_GROUPS


if __name__ == "__main__":
    main()
