#!/usr/bin/env python3
"""Unit tests for 38b_deletion_supplement_fig.py (S03 Task B).

What matters here is not the pixels -- it's the disclosure/censoring contract:

  1. `<20` in a count column must parse to NaN, never 0 and never 20 (feedback: "suppressed
     counts are not zero").
  2. A row whose bridged-haplotype denominator is itself censored must be marked
     `plottable=False` and excluded from the Wilson CI / bar, not silently coerced into one.
  3. The end-to-end script runs against a small synthetic copy of the 30_hla_sv tables and
     produces both a .pdf and a .png plus a README.md.

Run: python3 scripts/hla_popgen/tests/test_deletion_supplement_fig.py
"""
import importlib.util
import math
import os
import sys
import tempfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HLA_POPGEN_DIR = os.path.dirname(HERE)
if HLA_POPGEN_DIR not in sys.path:
    sys.path.insert(0, HLA_POPGEN_DIR)


def _load_module(filename, modname):
    path = os.path.join(HLA_POPGEN_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


m = _load_module("38b_deletion_supplement_fig.py", "hla_deletion_supp")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


def test_parse_count_censored_is_nan_not_zero_or_twenty():
    check("'<20' parses to NaN", math.isnan(m.parse_count("<20")))
    check("'0' parses to 0.0", m.parse_count("0") == 0.0)
    check("'137' parses to 137.0", m.parse_count("137") == 137.0)
    check("None parses to NaN", math.isnan(m.parse_count(None)))
    check("is_censored('<20') is True", m.is_censored("<20") is True)
    check("is_censored('137') is False", m.is_censored("137") is False)


def _fake_gene_rates_df():
    return pd.DataFrame([
        {"gene": "A", "pct_deleted_bridged": "1.5", "n_bridged_disp": "5000",
         "n_deleted_disp": "75", "gene_class": "classical_I"},
        {"gene": "DRB3", "pct_deleted_bridged": "49.0", "n_bridged_disp": "4000",
         "n_deleted_disp": "1960", "gene_class": "class_II_paralog"},
        # thin denominator -- must be dropped from the plot, not given a fabricated CI
        {"gene": "HFE", "pct_deleted_bridged": "0.0", "n_bridged_disp": "<20",
         "n_deleted_disp": "0", "gene_class": "other"},
    ])


def test_add_wilson_ci_drops_censored_denominator_rows():
    df = _fake_gene_rates_df()
    df["n_bridged"] = df["n_bridged_disp"].map(m.parse_count)
    df["n_bridged_censored"] = df["n_bridged_disp"].map(m.is_censored)
    df["n_deleted_disp_censored"] = df["n_deleted_disp"].map(m.is_censored)
    df["pct_deleted_bridged"] = df["pct_deleted_bridged"].astype(float)
    out = m.add_wilson_ci(df)

    row_a = out[out["gene"] == "A"].iloc[0]
    row_hfe = out[out["gene"] == "HFE"].iloc[0]
    check("a gene with a real denominator is plottable", bool(row_a["plottable"]))
    check("a gene with a censored (<20) denominator is NOT plottable",
          not bool(row_hfe["plottable"]))
    check("the plottable row gets a finite Wilson CI",
          not math.isnan(row_a["ci_lo"]) and not math.isnan(row_a["ci_hi"]))
    check("the dropped row's CI is NaN, not a fabricated interval",
          math.isnan(row_hfe["ci_lo"]) and math.isnan(row_hfe["ci_hi"]))
    check("the plottable row's CI brackets the point estimate",
          row_a["ci_lo"] <= row_a["pct_deleted_bridged"] <= row_a["ci_hi"])


def test_end_to_end_writes_pdf_png_and_readme():
    with tempfile.TemporaryDirectory() as tmp:
        sv_dir = os.path.join(tmp, "sv")
        out_dir = os.path.join(tmp, "out")
        os.makedirs(sv_dir, exist_ok=True)

        gene_rows = [
            {"gene": g, "pct_deleted_bridged": p, "n_bridged_disp": nb, "n_deleted_disp": nd,
             "gene_class": gc}
            for g, p, nb, nd, gc in [
                ("A", 1.8, 4000, 72, "classical_I"), ("B", 0.6, 4000, 24, "classical_I"),
                ("DRB3", 49.0, 3500, 1715, "class_II_paralog"),
                ("DRB4", 69.6, 3500, 2436, "class_II_paralog"),
                ("DRB5", 83.1, 3500, 2909, "class_II_paralog"),
                ("C4A", 11.0, 3600, 396, "complement"), ("C4B", 19.6, 3600, 706, "complement"),
                ("H", 12.2, 3900, 476, "pseudogene_I"),
            ]
        ]
        pd.DataFrame(gene_rows).to_csv(
            os.path.join(sv_dir, "deletion_rates_by_gene.tsv"), sep="\t", index=False)

        anc_rows = []
        for anc in ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]:
            for g, p, gc in [("DRB3", 45.0, "class_II_paralog"), ("C4A", 10.0, "complement")]:
                anc_rows.append({"ancestry": anc, "gene": g, "pct_deleted_bridged": p,
                                 "n_bridged_disp": 600, "n_deleted_disp": 60, "gene_class": gc})
        pd.DataFrame(anc_rows).to_csv(
            os.path.join(sv_dir, "deletion_rates_by_gene_ancestry.tsv"), sep="\t", index=False)

        old_argv = sys.argv
        sys.argv = ["38b_deletion_supplement_fig.py", "--sv-dir", sv_dir, "--out-dir", out_dir]
        try:
            m.main()
        finally:
            sys.argv = old_argv

        pdf_path = os.path.join(out_dir, "supp_deletions.pdf")
        png_path = os.path.join(out_dir, "supp_deletions.png")
        readme_path = os.path.join(out_dir, "README_38b_supplement.md")
        check("wrote supp_deletions.pdf",
              os.path.exists(pdf_path) and os.path.getsize(pdf_path) > 0)
        check("wrote supp_deletions.png",
              os.path.exists(png_path) and os.path.getsize(png_path) > 0)
        check("wrote README.md", os.path.exists(readme_path) and os.path.getsize(readme_path) > 0)
        with open(readme_path) as f:
            readme_text = f.read()
        check("README mentions the negative controls", "negative control" in readme_text.lower())
        check("README mentions the CNV genes", "DRB3" in readme_text)


def main():
    print("test_deletion_supplement_fig.py")
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
