#!/usr/bin/env python3
"""Unit tests for 37b_ld_supplement_table.py (S03 Task C).

The contract that matters: `n_hap_ij_disp` -> status must distinguish a genuine, disclosable
zero ("0" -> not_observed) from a censored 1-19 ("<20" -> suppressed_lt20, NaN n_hap, never read
as 0 or 20), and the heatmap must never fold a suppressed/not_observed cell into the D' color
scale.

Run: python3 scripts/hla_popgen/tests/test_ld_supplement_table.py
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


m = _load_module("37b_ld_supplement_table.py", "hla_ld_supp")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


def test_status_from_disp():
    check("'0' -> not_observed", m.status_from_disp("0") == m.STATUS_NOT_OBSERVED)
    check("'<20' -> suppressed_lt20", m.status_from_disp("<20") == m.STATUS_SUPPRESSED)
    check("'137' -> estimated", m.status_from_disp("137") == m.STATUS_ESTIMATED)
    check("'20' -> estimated (the floor itself is disclosable)",
          m.status_from_disp("20") == m.STATUS_ESTIMATED)


def test_n_hap_from_disp_never_guesses_a_suppressed_value():
    check("suppressed cell -> NaN n_hap",
          math.isnan(m.n_hap_from_disp("<20", m.STATUS_SUPPRESSED)))
    check("not_observed cell -> 0.0 n_hap",
          m.n_hap_from_disp("0", m.STATUS_NOT_OBSERVED) == 0.0)
    check("estimated cell -> real n_hap",
          m.n_hap_from_disp("137", m.STATUS_ESTIMATED) == 137.0)


def test_build_tidy_table_has_the_required_columns_and_no_fabricated_counts():
    raw = pd.DataFrame([
        {"pair": "DQA1~DQB1", "ancestry": "AFR", "allele_a": "DQA1*01:01",
         "allele_b": "DQB1*05:01", "freq_a": "0.3", "freq_b": "0.2", "r2": "0.5",
         "Dprime": "0.9", "n_hap_ij_disp": "301"},
        {"pair": "DQA1~DQB1", "ancestry": "AFR", "allele_a": "DQA1*01:01",
         "allele_b": "DQB1*06:03", "freq_a": "0.3", "freq_b": "0.01", "r2": "0.02",
         "Dprime": "0.1", "n_hap_ij_disp": "<20"},
        {"pair": "DQA1~DQB1", "ancestry": "AFR", "allele_a": "DQA1*01:01",
         "allele_b": "DQB1*02:01", "freq_a": "0.3", "freq_b": "0.07", "r2": "0.006",
         "Dprime": "-0.2", "n_hap_ij_disp": "0"},
    ])
    raw["status"] = raw["n_hap_ij_disp"].map(m.status_from_disp)
    raw["n_hap"] = raw.apply(lambda r: m.n_hap_from_disp(r["n_hap_ij_disp"], r["status"]), axis=1)
    tidy = m.build_tidy_table(raw)

    check("tidy table has exactly the required columns",
          list(tidy.columns) == ["gene_pair", "ancestry", "allele_a", "allele_b", "n_hap",
                                  "freq_a", "freq_b", "r2", "Dprime", "status"],
          str(list(tidy.columns)))
    check("gene_pair column is populated from 'pair'",
          set(tidy["gene_pair"]) == {"DQA1~DQB1"})
    row_sup = tidy[tidy["status"] == m.STATUS_SUPPRESSED].iloc[0]
    check("the suppressed row's n_hap is NaN in the tidy table", math.isnan(row_sup["n_hap"]))
    row_not_obs = tidy[tidy["status"] == m.STATUS_NOT_OBSERVED].iloc[0]
    check("the not_observed row's n_hap is 0 in the tidy table", row_not_obs["n_hap"] == 0.0)
    row_est = tidy[tidy["status"] == m.STATUS_ESTIMATED].iloc[0]
    check("the estimated row keeps its real n_hap", row_est["n_hap"] == 301.0)


def test_end_to_end_writes_csv_heatmap_and_readme():
    with tempfile.TemporaryDirectory() as tmp:
        ld_dir = os.path.join(tmp, "ld")
        out_dir = os.path.join(tmp, "out")
        os.makedirs(ld_dir, exist_ok=True)

        rows = []
        alleles_a = ["DQA1*01:01", "DQA1*01:02", "DQA1*02:01"]
        alleles_b = ["DQB1*02:01", "DQB1*05:01", "DQB1*06:03"]
        for anc in ["AFR", "EUR"]:
            for i, aa in enumerate(alleles_a):
                for j, bb in enumerate(alleles_b):
                    disp = "0" if (i + j) % 3 == 0 else ("<20" if (i + j) % 3 == 1 else "150")
                    rows.append({"allele_a": aa, "allele_b": bb, "freq_a": 0.2, "freq_b": 0.15,
                                "freq_hap": 0.03, "D": 0.01, "Dprime": 0.4, "r2": 0.05,
                                "chi2_1df": 1.0, "pair": "DQA1~DQB1", "ancestry": anc,
                                "n_hap_ij_disp": disp})
            # also a pair NOT in the heatmap scope, to check the CSV still covers it
            rows.append({"allele_a": "A*01:01", "allele_b": "B*07:02", "freq_a": 0.1,
                        "freq_b": 0.1, "freq_hap": 0.01, "D": 0.001, "Dprime": 0.2, "r2": 0.01,
                        "chi2_1df": 0.5, "pair": "A~B", "ancestry": anc, "n_hap_ij_disp": "45"})
        pd.DataFrame(rows).to_csv(os.path.join(ld_dir, "ld_pairwise.tsv"), sep="\t", index=False)

        old_argv = sys.argv
        sys.argv = ["37b_ld_supplement_table.py", "--ld-dir", ld_dir, "--out-dir", out_dir]
        try:
            m.main()
        finally:
            sys.argv = old_argv

        csv_path = os.path.join(out_dir, "supp_table_ld_pairs.csv")
        readme_path = os.path.join(out_dir, "README.md")
        check("wrote the CSV", os.path.exists(csv_path) and os.path.getsize(csv_path) > 0)
        check("wrote README.md", os.path.exists(readme_path) and os.path.getsize(readme_path) > 0)

        out_df = pd.read_csv(csv_path)
        check("CSV covers the non-heatmap pair too (A~B)", "A~B" in set(out_df["gene_pair"]))
        check("CSV has all three status values",
              {m.STATUS_NOT_OBSERVED, m.STATUS_SUPPRESSED, m.STATUS_ESTIMATED} <=
              set(out_df["status"]))

        pdf_path = os.path.join(out_dir, "supp_heatmap_dq_dp.pdf")
        png_path = os.path.join(out_dir, "supp_heatmap_dq_dp.png")
        check("wrote the heatmap .pdf",
              os.path.exists(pdf_path) and os.path.getsize(pdf_path) > 0)
        check("wrote the heatmap .png",
              os.path.exists(png_path) and os.path.getsize(png_path) > 0)


def main():
    print("test_ld_supplement_table.py")
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
