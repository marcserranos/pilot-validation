#!/usr/bin/env python3
"""Deep immune/tumor/HLA/autoimmune disease breakdown for the LR x RNA-seq cohort.

Runs entirely against the phenotype files query_overlap_phenotypes.py already cached at
~/pipeline_outputs/rnaseq/pheno/ (person.tsv, observation_period.tsv, conditions.tsv).
NO BigQuery call, NO mount needed -- just pandas over data already on the VM's persistent
disk. Cheap to re-run.

WHY THIS EXISTS: analyze_disease_burden.py's immune_feasibility table groups diseases into
10 broad FAMILIES using crude substring matching on the SNOMED name ("arthritis" catches
osteoarthritis, "hypersensitivity" catches trivial drug notes). That was the right first
pass -- fast, told us the families were viable at all -- but it cannot answer "how many
people have rheumatoid arthritis specifically" or "how many have CLL vs Hodgkin lymphoma".

This does the finer cut, on FOUR axes the user asked about directly:
  1. HLA-LINKED  -- diseases with an established, well-known HLA association (ankylosing
     spondylitis/B27, celiac/DQ2-DQ8, type 1 diabetes/DR3-DR4, narcolepsy/DQB1*06:02,
     psoriasis/Cw6, MS/DRB1*15:01...). These are the ones where a repertoire-HLA link is
     most directly testable.
  2. AUTOIMMUNE -- 20 specific diseases, not 5 families, each with its own count.
  3. IMMUNODEFICIENCY -- specific subtypes (CVID, IgA deficiency, HIV, neutropenia...)
     rather than one lumped "immunodeficiency" row.
  4. TUMORS / MALIGNANCY -- split into LYMPHOID (where the tumor cells themselves ARE
     clonal receptor sequences -- CLL, lymphoma, myeloma, MGUS -- directly relevant to
     TRUST4 output) versus SOLID (general immune-surveillance context, not repertoire-
     specific).

METHOD, and why it is more trustworthy than the family-level screen: each specific disease
is matched primarily by its 3-CHARACTER ICD-10 CODE (E10 = type 1 diabetes vs E11 = type 2,
a distinction pure name-matching cannot make reliably), with a SNOMED name-substring
fallback only where ICD alone doesn't disambiguate (e.g. graft-versus-host disease has no
clean dedicated ICD-10 category). Still a screen, not a validated phenotype -- flag any of
these before quoting them in a paper -- but a materially better one than family-level
substring matching.

Usage:
  pixi run python3 deep_immune_breakdown.py [--pheno ~/pipeline_outputs/rnaseq/pheno]
      [--cohort ~/pipeline_outputs/rnaseq/lr_rnaseq_overlap_cohort.tsv] [--min-cell 20]
      [--outdir ../results]
"""
import argparse
import os
import sys

import pandas as pd

# ---------------------------------------------------------------------------
# Each entry: (label, [icd3_codes], [name_substrings_OR_None])
# If icd3_codes is non-empty, a person matches if EITHER an icd3 code matches OR a name
# substring matches (fallback catches cases where the source vocabulary isn't ICD-10, or
# the mapping is imprecise). If icd3_codes is empty, name substring is the only route.
# ---------------------------------------------------------------------------

HLA_LINKED = [
    ("Ankylosing spondylitis (HLA-B27)",        ["M45"], ["ankylosing spondylitis"]),
    ("Celiac disease (HLA-DQ2/DQ8)",            [], ["celiac", "coeliac"]),  # K90 dropped: includes other malabsorption syndromes
    ("Type 1 diabetes (HLA-DR3/DR4)",           ["E10"], ["type 1 diabetes"]),
    ("Psoriasis (HLA-Cw6)",                     ["L40"], ["psoriasis"]),
    ("Psoriatic arthritis (HLA-Cw6/B27)",       [], ["psoriatic arthritis"]),  # M07 dropped: includes enteropathic arthropathy
    ("Multiple sclerosis (HLA-DRB1*15:01)",     ["G35"], ["multiple sclerosis"]),
    ("Rheumatoid arthritis (HLA shared epitope)",["M05","M06"], ["rheumatoid arthritis"]),
    ("Systemic lupus erythematosus (HLA-DR2/DR3)",["M32"], ["systemic lupus", "lupus erythematosus"]),
    ("Graves disease (HLA-DR3)",                [], ["graves"]),  # E05 dropped (toxic nodular goiter); single word safe against apostrophe
    ("Narcolepsy (HLA-DQB1*06:02)",             [], ["narcolepsy"]),  # G47 dropped: dominated by sleep apnea, not narcolepsy-specific
    ("Behcet disease (HLA-B51)",                [], ["behcet"]),  # M35 dropped: shared block with Sjogren/PMR/other connective tissue disease
]

AUTOIMMUNE_SPECIFIC = [
    ("Rheumatoid arthritis",        ["M05","M06"], ["rheumatoid arthritis"]),
    ("Systemic lupus erythematosus",["M32"], ["systemic lupus", "lupus erythematosus"]),
    ("Sjogren syndrome",            [], ["sjogren", "sjögren"]),  # M35 dropped: shared block with Behcet/PMR/other
    ("Systemic sclerosis",          ["M34"], ["systemic sclerosis", "scleroderma"]),
    ("Ankylosing spondylitis",      ["M45"], ["ankylosing spondylitis"]),
    ("Psoriatic arthritis",         [], ["psoriatic arthritis"]),  # M07 dropped: includes enteropathic arthropathy
    ("Vasculitis (systemic)",       ["M30","M31"], ["vasculitis", "polyarteritis"]),
    ("Type 1 diabetes",             ["E10"], ["type 1 diabetes"]),
    ("Hashimoto / autoimmune thyroiditis", [], ["hashimoto", "autoimmune thyroiditis"]),  # E06 dropped: includes non-autoimmune (viral) thyroiditis
    ("Graves disease",              [], ["graves"]),  # single word: safe against "Graves' disease" apostrophe
    ("Addison disease",             [], ["addison"]),  # E27 dropped: includes Cushing syndrome (opposite condition)
    ("Crohn disease",               ["K50"], ["crohn"]),
    ("Ulcerative colitis",          ["K51"], ["ulcerative colitis"]),
    ("Celiac disease",              [], ["celiac", "coeliac"]),  # K90 dropped: includes other malabsorption syndromes
    ("Autoimmune hepatitis",        [], ["autoimmune hepatitis"]),  # K75 dropped: includes other inflammatory liver disease
    ("Multiple sclerosis",          ["G35"], ["multiple sclerosis"]),
    ("Myasthenia gravis",           [], ["myasthenia gravis"]),  # G70 dropped: includes other myoneural disorders
    ("Guillain-Barre syndrome",     [], ["guillain"]),  # G61 dropped: includes other inflammatory polyneuropathy
    ("Psoriasis",                   ["L40"], ["psoriasis"]),
    ("Vitiligo",                    ["L80"], ["vitiligo"]),
    ("Pemphigus",                   ["L10"], ["pemphigus"]),
    ("Alopecia areata",             ["L63"], ["alopecia areata"]),
]

IMMUNODEFICIENCY_SPECIFIC = [
    ("HIV",                                    ["B20"], ["hiv", "human immunodeficiency virus"]),
    ("Common variable immunodeficiency (CVID)", ["D83"], ["common variable immunodeficiency", "cvid"]),
    ("Selective IgA / antibody deficiency",     ["D80"], ["iga deficiency", "agammaglobulin", "hypogammaglobulin"]),
    ("Other immunodeficiency (D84)",            ["D84"], ["immunodeficiency"]),
    ("Immune-mechanism disorder, other (D89)",  ["D89"], []),
    ("Neutropenia",                             ["D70"], ["neutropenia"]),
    ("Transplant status (any organ)",           ["Z94"], ["transplant status"]),
    ("Graft-versus-host disease",               [], ["graft-versus-host", "graft versus host", " gvhd"]),
    ("Immunosuppressed state (drug-induced)",   [], ["immunosuppress"]),
]

TUMOR_LYMPHOID = [
    ("Hodgkin lymphoma",                ["C81"], ["hodgkin lymphoma"]),
    ("Non-Hodgkin lymphoma",            ["C82","C83","C85"], ["non-hodgkin lymphoma", "non hodgkin lymphoma"]),
    ("Mycosis fungoides / cutaneous T-cell lymphoma", ["C84"], ["mycosis fungoides", "cutaneous t-cell lymphoma"]),
    ("Lymphoid leukemia (incl. CLL/ALL)",["C91"], ["lymphocytic leukemia", "lymphoblastic leukemia"]),
    ("Myeloid leukemia",                ["C92"], ["myeloid leukemia"]),
    ("Multiple myeloma",                ["C90"], ["multiple myeloma", "plasma cell myeloma"]),
    ("MGUS (monoclonal gammopathy)",    [], ["monoclonal gammopathy", "mgus"]),  # D47 dropped: includes other uncertain-behavior lymphoid neoplasms
    ("Myelodysplastic syndrome",        ["D46"], ["myelodysplastic"]),
]

TUMOR_SOLID = [
    ("Melanoma",        ["C43"], ["melanoma"]),
    ("Breast cancer",   ["C50"], ["breast cancer", "malignant neoplasm of breast"]),
    ("Lung cancer",     ["C34"], ["lung cancer", "malignant neoplasm of lung", "bronchus"]),
    ("Colorectal cancer",["C18","C19","C20"], ["colon cancer", "colorectal", "rectal cancer"]),
    ("Prostate cancer", ["C61"], ["prostate cancer", "malignant neoplasm of prostate"]),
    ("Bladder cancer",  ["C67"], ["bladder cancer"]),
    ("Kidney cancer",   ["C64"], ["kidney cancer", "renal cell"]),
]

# Added 2026-08-27 per Cole Shanks (Slack): "another disease we should actually look at is
# Alzheimers... Parkinsons disease, Atherosclerosis/CAD, Long Covid, which actually has an
# ICD code now (ICD-10 U09.9)". These don't share the HLA-linked/autoimmune/tumor mechanism
# of the sections above -- they're a feasibility check (do we even have enough cases?), not
# a repertoire-mechanism claim. Alzheimer's and Parkinson's both have reported (weaker,
# non-classical) HLA-region GWAS hits; atherosclerosis/CAD has an established inflammatory-
# immune component; Long Covid is post-viral immune dysregulation -- so a TCR/BCR repertoire
# link isn't absurd for any of them, just far less direct than e.g. T1D or ankylosing
# spondylitis.
OTHER_CONDITIONS_OF_INTEREST = [
    ("Alzheimer disease",                       ["G30"], ["alzheimer"]),
    ("Parkinson disease",                       ["G20"], ["parkinson"]),
    ("Atherosclerosis / coronary artery disease",["I70","I25"], ["atherosclero", "coronary artery disease"]),
    ("Long COVID (post COVID-19 condition)",    ["U09"], ["post covid", "post-covid", "long covid"]),
]

SECTIONS = [
    ("HLA-LINKED DISEASES", HLA_LINKED,
     "Diseases with an established genetic association to specific HLA alleles. These "
     "are the most direct candidates for a repertoire<->HLA link, since HLA type shapes "
     "which receptors get positively/negatively selected in the thymus."),
    ("AUTOIMMUNE -- specific diseases", AUTOIMMUNE_SPECIFIC,
     "The immune system attacking the body's own tissue. Specific diseases, not the 5 "
     "broad families from the first-pass screen."),
    ("IMMUNODEFICIENCY -- specific subtypes", IMMUNODEFICIENCY_SPECIFIC,
     "States where the immune system is weakened -- inherited, acquired, or induced by "
     "transplant drugs. Directly relevant to repertoire interpretation: an "
     "immunosuppressed person's repertoire reflects a suppressed immune system, not a "
     "normal one."),
    ("TUMORS -- lymphoid / hematologic (repertoire-relevant)", TUMOR_LYMPHOID,
     "*** Most relevant tumor category for this project ***. In these cancers the "
     "malignant cells ARE a single expanded T-cell or B-cell clone -- so TRUST4 output "
     "isn't just informative about these diseases, in principle it can directly detect "
     "them, the same way clinical MRD (measurable residual disease) testing works."),
    ("TUMORS -- solid (general immune-surveillance context)", TUMOR_SOLID,
     "Common solid cancers. Included for scale/context, not because repertoire data "
     "is expected to be as directly diagnostic as for the lymphoid tumors above."),
    ("OTHER CONDITIONS OF INTEREST (req. Cole Shanks, 2026-08-27)", OTHER_CONDITIONS_OF_INTEREST,
     "Alzheimer's, Parkinson's, atherosclerosis/CAD, and Long COVID. A feasibility check "
     "(do we have enough cases at all), not an established HLA/repertoire mechanism claim "
     "the way the sections above are. NOTE on Long COVID: U09 is a young code (introduced "
     "Oct 2021 in ICD-10-CM) -- anyone diagnosed with post-COVID sequelae before their EHR "
     "records started using it will be undercounted here."),
]


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def suppress(n, min_cell):
    return f"<{min_cell}" if 0 < n < min_cell else str(int(n))


def match_mask(cond, icd3_codes, name_subs):
    mask = pd.Series(False, index=cond.index)
    if icd3_codes:
        mask |= cond["icd3"].isin(icd3_codes)
    if name_subs:
        names = cond["condition_name"].astype(str).str.lower()
        for t in name_subs:
            mask |= names.str.contains(t, regex=False, na=False)
    return mask


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pheno", default=os.path.expanduser("~/pipeline_outputs/rnaseq/pheno"))
    ap.add_argument("--cohort", default=os.path.expanduser(
                        "~/pipeline_outputs/rnaseq/lr_rnaseq_overlap_cohort.tsv"))
    ap.add_argument("--min-cell", type=int, default=20)
    ap.add_argument("--outdir", default=os.path.join(
                        os.path.dirname(os.path.abspath(__file__)), "..", "results"))
    args = ap.parse_args()

    cond_path = os.path.join(os.path.expanduser(args.pheno), "conditions.tsv")
    obs_path = os.path.join(os.path.expanduser(args.pheno), "observation_period.tsv")
    for p in (cond_path, obs_path, os.path.expanduser(args.cohort)):
        if not os.path.exists(p):
            die(f"not found: {p}. Run query_overlap_phenotypes.py first (already cached "
                f"from the earlier run -- this script should not need BigQuery).")

    cohort = pd.read_csv(os.path.expanduser(args.cohort), sep="\t", dtype=str)
    obs = pd.read_csv(obs_path, sep="\t", dtype={"research_id": str})
    cond = pd.read_csv(cond_path, sep="\t", dtype={"research_id": str}, low_memory=False)

    ehr_ids = set(obs.loc[pd.to_numeric(obs["ehr_years"], errors="coerce") > 0,
                          "research_id"])
    cond = cond[cond["research_id"].isin(ehr_ids)].copy()
    denom = len(ehr_ids)
    print(f"Denominator (people with a non-zero EHR window): {denom:,}\n")

    cond["icd3"] = cond["source_code"].astype(str).str[:3]

    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    os.makedirs(outdir, exist_ok=True)
    all_rows = []

    for section_title, terms, blurb in SECTIONS:
        print("=" * 72)
        print(section_title)
        print("=" * 72)
        print(blurb)
        print()
        rows = []
        for label, icd3_codes, name_subs in terms:
            mask = match_mask(cond, icd3_codes, name_subs)
            n = cond.loc[mask, "research_id"].nunique()
            rows.append({"section": section_title, "condition": label, "n_people": n,
                        "pct_of_cohort": round(n / max(denom, 1) * 100, 2)})
        sec_df = pd.DataFrame(rows).sort_values("n_people", ascending=False)
        for _, r in sec_df.iterrows():
            print(f"  {r['condition']:<48} {suppress(r['n_people'], args.min_cell):>8}"
                  f"   ({r['pct_of_cohort']:.2f}%)")
        print()
        all_rows.append(sec_df)

    full = pd.concat(all_rows, ignore_index=True)
    full["n_people_reportable"] = full["n_people"].apply(
        lambda n: suppress(n, args.min_cell))
    full.loc[full["n_people"] < args.min_cell, "pct_of_cohort"] = pd.NA
    out = os.path.join(outdir, "lr_rnaseq_deep_immune_breakdown.csv")
    full[["section", "condition", "n_people_reportable", "pct_of_cohort"]].to_csv(
        out, index=False)
    print(f"De-identified detail (safe to commit): {out}")
    print("\nCAVEAT: ICD-10 3-char code matching is more precise than name substring alone "
          "but is still a screen, not a validated phenotype. Name-substring fallbacks "
          "(used where ICD alone underspecifies, e.g. GVHD) inherit the same over-matching "
          "risk as the first-pass family screen.")


if __name__ == "__main__":
    main()
