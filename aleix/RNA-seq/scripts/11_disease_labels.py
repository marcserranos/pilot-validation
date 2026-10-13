#!/usr/bin/env python3
"""11 -- W2: case / control labels for ~32 chronic diseases (report: reports/11_disease_labels/).

Cole's request (call 2026-10-06): about 30 chronic conditions, classic autoimmune plus
not-classically-immune ones (hypertension, CKD, heart disease, T2D), defined cleanly even at
the cost of data. Chronic only: blood-draw dates exist for few people, and a chronic
diagnosis is permanent, so its timing relative to the draw matters less.

DEFINITIONS (DISEASES below; written to disease_definitions.csv):
  case      >= 2 recorded occurrences of the disease's ICD-10-CM or ICD-9-CM codes (summed
            over codes). One occurrence is often a rule-out or a billing artefact.
  excluded  exactly 1 occurrence (ambiguous), or any code of the wider exclusion family
            (e.g. any diabetes code for a T2D control), or no EHR window
  control   no code of the disease's case or exclusion family, an EHR window >= 2 years,
            and at least one ICD-coded condition of any kind (shows the person's record
            carries ICD data at all; a person with no ICD rows cannot be called disease-free)

TIERS by cases in analysis_main: >= 150 classify (W5) and associate (W4); 50-149 associate
only; 20-49 described only; < 20 suppressed.

Population: analysis_main from report 09 (read floor, EHR window, adult, unrelated).

Outputs: VM-local analysis/disease_labels.pkl (person x disease: 1 case, 0 control, NaN
excluded). Committed: disease_definitions.csv, case_counts.csv (overall and by ancestry,
cells < 20 suppressed), summary.csv.

Usage (from aleix/RNA-seq/, after 09):
  pixi run python3 -u scripts/11_disease_labels.py
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import repfig as R  # noqa: E402

MIN_CASE_OCC = 2
MIN_CONTROL_EHR_YEARS = 2.0
TIERS = [(150, "classify + associate"), (50, "associate"), (20, "describe")]

# key, label, group, ICD-10 case prefixes, ICD-9 case prefixes (or regex), extra ICD-10
# exclusion prefixes, extra ICD-9 exclusion prefixes. Codes are written without dots.
DISEASES = [
    ("htn", "Hypertension", "cardiometabolic", "I10 I11 I12 I13 I15 I16", "401 402 403 404 405", "", ""),
    ("t2d", "Type 2 diabetes", "cardiometabolic", "E11", r"re:^250\d[02]", "E08 E09 E10 E13", "250"),
    ("t1d", "Type 1 diabetes", "autoimmune", "E10", r"re:^250\d[13]", "E08 E09 E11 E13", "250"),
    ("hld", "Hyperlipidaemia", "cardiometabolic", "E78", "272", "", ""),
    ("obesity", "Obesity", "cardiometabolic", "E66", "2780", "", ""),
    ("cad", "Coronary artery disease", "cardiometabolic", "I20 I21 I22 I23 I24 I25", "410 411 412 413 414", "", ""),
    ("hf", "Heart failure", "cardiometabolic", "I50 I110 I130 I132", "428", "I11 I13", "402 404"),
    ("af", "Atrial fibrillation", "cardiometabolic", "I48", "4273", "", ""),
    ("ckd", "Chronic kidney disease", "renal", "N18 I12 I13", "585 403 404", "N19 Z992 Z940", "586 V420 V451"),
    ("stroke", "Ischaemic stroke", "cardiometabolic", "I63", "433 434", "I64 I69 G45", "436 438"),
    ("pad", "Peripheral artery disease", "cardiometabolic", "I739 I702", "4439 4402", "I70 I73", "440 443"),
    ("asthma", "Asthma", "respiratory / allergic", "J45", "493", "", ""),
    ("copd", "COPD", "respiratory / allergic", "J44 J43", "491 492 496", "J41 J42", "490"),
    ("ar", "Allergic rhinitis", "respiratory / allergic", "J30", "477", "", ""),
    ("ad", "Atopic dermatitis", "respiratory / allergic", "L20", "6918", "", ""),
    ("urticaria", "Urticaria", "respiratory / allergic", "L50", "708", "", ""),
    ("ra", "Rheumatoid arthritis", "autoimmune", "M05 M06", "7140", "M08", "714"),
    ("psoriasis", "Psoriasis (incl. psoriatic arthritis)", "autoimmune", "L40", "6960 6961", "", ""),
    ("sle", "Systemic lupus erythematosus", "autoimmune", "M32", "7100", "", ""),
    ("sjogren", "Sjogren syndrome", "autoimmune", "M350", "7102", "", ""),
    ("hypothyroid", "Hypothyroidism / autoimmune thyroiditis", "autoimmune", "E03 E063", "244 2452", "E06 E02", "245"),
    ("thyrotox", "Thyrotoxicosis / Graves", "autoimmune", "E05", "242", "", ""),
    ("crohn", "Crohn disease", "autoimmune", "K50", "555", "K51", "556"),
    ("uc", "Ulcerative colitis", "autoimmune", "K51", "556", "K50", "555"),
    ("celiac", "Coeliac disease", "autoimmune", "K900", "5790", "", ""),
    ("ms", "Multiple sclerosis", "autoimmune", "G35", "340", "G36 G37", ""),
    ("as", "Ankylosing spondylitis", "autoimmune", "M45", "7200", "", ""),
    ("vitiligo", "Vitiligo", "autoimmune", "L80", "70901", "", ""),
    ("alopecia", "Alopecia areata", "autoimmune", "L63", "70401", "", ""),
    ("hiv", "HIV", "chronic infection", "B20 Z21", "042 V08", "", ""),
    ("hep", "Chronic viral hepatitis B / C", "chronic infection", "B18", "0702 0703 0704 0705 0707", "B15 B16 B17 B19", "070"),
    ("lymphoid", "Lymphoid / plasma-cell neoplasm", "haematologic", "C81 C82 C83 C84 C85 C86 C88 C90 C91 D472", "200 201 202 203 204 2731", "", ""),
    ("parkinson", "Parkinson disease", "neurologic", "G20", "332", "", ""),
]


def matcher(spec):
    """Space-separated code prefixes, or 're:<regex>'. Returns a function on code arrays."""
    spec = spec.strip()
    if not spec:
        return lambda codes: np.zeros(len(codes), bool)
    if spec.startswith("re:"):
        rx = spec[3:]
        return lambda codes: pd.Series(codes).str.match(rx).fillna(False).to_numpy()
    pre = tuple(spec.split())
    return lambda codes: pd.Series(codes).str.startswith(pre).fillna(False).to_numpy()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    base = "~/pipeline_outputs/rnaseq"
    ap.add_argument("--person-table", default=f"{base}/analysis/person_table.pkl")
    ap.add_argument("--conditions", default=f"{base}/pheno/conditions.tsv")
    ap.add_argument("--local-out", default=f"{base}/analysis")
    ap.add_argument("--outdir", default=f"{base}/reports/11_disease_labels")
    args = ap.parse_args()
    outdir, local = os.path.expanduser(args.outdir), os.path.expanduser(args.local_out)
    os.makedirs(outdir, exist_ok=True)

    pt = pd.read_pickle(os.path.expanduser(args.person_table))
    pt = pt[pt["analysis_main"]].set_index("research_id")
    cond = pd.read_csv(os.path.expanduser(args.conditions), sep="\t",
                       dtype={"research_id": str, "source_code": str, "source_vocabulary": str},
                       usecols=["research_id", "source_vocabulary", "source_code", "n_occurrences"])
    cond = cond[cond["research_id"].isin(pt.index)].copy()
    voc = cond["source_vocabulary"].fillna("")
    cond["v"] = np.where(voc.str.startswith("ICD10"), "10", np.where(voc.str.startswith("ICD9"), "9", ""))
    frac_icd = float((cond["v"] != "").mean())
    cond = cond[cond["v"] != ""].copy()
    cond["code"] = cond["source_code"].fillna("").str.replace(".", "", regex=False).str.upper()
    cond["n_occurrences"] = pd.to_numeric(cond["n_occurrences"], errors="coerce").fillna(1)
    has_icd = set(cond["research_id"])
    print(f"analysis_main: {len(pt):,} people; {100 * frac_icd:.1f}% of their condition rows are "
          f"ICD-coded; {len(has_icd):,} people have >= 1 ICD row", file=sys.stderr)

    c10, c9 = cond[cond["v"] == "10"], cond[cond["v"] == "9"]
    codes10, codes9 = c10["code"].to_numpy(), c9["code"].to_numpy()
    labels, defs, counts = {}, [], []
    for key, label, group, i10, i9, x10, x9 in DISEASES:
        m10, m9 = matcher(i10)(codes10), matcher(i9)(codes9)
        occ = pd.concat([c10.loc[m10, ["research_id", "n_occurrences"]],
                         c9.loc[m9, ["research_id", "n_occurrences"]]]).groupby("research_id")["n_occurrences"].sum()
        e10 = m10 | matcher(x10)(codes10)
        e9 = m9 | matcher(x9)(codes9)
        touched = set(c10.loc[e10, "research_id"]) | set(c9.loc[e9, "research_id"])
        lab = pd.Series(np.nan, index=pt.index)
        occ = occ.reindex(pt.index).fillna(0)
        lab[occ >= MIN_CASE_OCC] = 1.0
        ctrl = (~pt.index.isin(touched)) & (pt["ehr_years"] >= MIN_CONTROL_EHR_YEARS).to_numpy() \
            & pt.index.isin(has_icd)
        lab[ctrl] = 0.0
        labels[key] = lab
        n_case, n_ctrl = int((lab == 1).sum()), int((lab == 0).sum())
        tier = next((t for k, t in TIERS if n_case >= k), "suppressed (< 20 cases)")
        defs.append({"key": key, "disease": label, "group": group, "icd10_case": i10,
                     "icd9_case": i9, "icd10_extra_exclusion": x10, "icd9_extra_exclusion": x9,
                     "min_case_occurrences": MIN_CASE_OCC,
                     "min_control_ehr_years": MIN_CONTROL_EHR_YEARS})
        row = {"key": key, "disease": label, "group": group, "tier": tier,
               "cases": n_case if n_case >= R.MIN_PEOPLE or n_case == 0 else "<20",
               "controls": n_ctrl, "excluded": int(lab.isna().sum()),
               "case_prevalence_pct": round(100 * n_case / max(n_case + n_ctrl, 1), 2)
               if n_case >= R.MIN_PEOPLE else np.nan}
        for a in R.ANCESTRY_ORDER:
            k = int(((lab == 1) & (pt["ancestry"] == a)).sum())
            row[f"cases_{a}"] = k if k >= R.MIN_PEOPLE or k == 0 else "<20"
        counts.append(row)

    L = pd.DataFrame(labels)
    L.index.name = "research_id"
    L.to_pickle(os.path.join(local, "disease_labels.pkl"))
    pd.DataFrame(defs).to_csv(os.path.join(outdir, "disease_definitions.csv"), index=False)
    ct = pd.DataFrame(counts)
    ct.to_csv(os.path.join(outdir, "case_counts.csv"), index=False)
    pd.DataFrame([("people_analysis_main", len(pt)), ("people_with_icd_rows", len(has_icd)),
                  ("frac_condition_rows_icd", frac_icd), ("diseases", len(DISEASES))]
                 + [(f"n_{t}", int((ct["tier"] == t).sum())) for _, t in TIERS],
                 columns=["metric", "value"]).to_csv(os.path.join(outdir, "summary.csv"), index=False)
    print("\n" + ct[["disease", "tier", "cases", "controls", "excluded"]
                    + [f"cases_{a}" for a in R.ANCESTRY_ORDER]].to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
