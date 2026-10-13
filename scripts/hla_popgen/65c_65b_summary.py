#!/usr/bin/env python3
"""65c: summary-only PUBLIC tables for the 65b reads-vs-assembly validation (TRB insertion, TRBV4-3/3-2).

Why a separate step: the 65 PUBLIC count tables (genotype_by_stratum, confusion cells) cover people who are
also in the committed 65 cuts tables, so publishing a second set of exact counts lets small counts be
recovered by subtraction (DECISIONS 2026-10-13, "Publishing the TRB mixture tables"; Marc: summary only).
This step reads the INTERNAL per-person tables on the VM and writes rates only; counts are shown only as
denominators rounded to the nearest 10 ('<20' for 1-19). Rates follow _disclosure.mask_rate (blank when the
numerator or denominator is 1-19).

Inputs (INTERNAL, VM only)
  --reads-gene    ~/s05/internal/65_reads/per_person_gene.tsv (whole cohort; the mixture is fitted on it)
  --asm-dir       ~/s05/internal/65_assembly (per_person_gene.tsv, per_haplotype.tsv, failures.tsv)
  --people-file   ~/s05/internal/65_people_strat2000.tsv (the stratified sample, unrelated by design)
Outputs (PUBLIC, --out-dir)
  65c_completion.tsv   per stratum: sampled (rounded), % genotyped, % haplotypes complete, % TRBV4-3 called
  65c_failures.tsv     per failure reason: people (rounded / '<20'), ever-failed and still-failed
  65c_concordance.tsv  per reads genotyper (stored, cuts, mixture) x feature x stratum: n both called
                       (rounded), % concordant + Wilson CI, and per copy number c: % of reads-c calls the
                       assembly also calls c (ppv_c) and % of assembly-c calls the reads also call c (sens_c)
  65c_hwe.tsv          per caller (assembly, reads_cuts, reads_mixture) on the same people (called by all
                       three) x feature x ancestry: n (rounded), insertion AF, % 0 copies, observed and
                       HWE-expected % heterozygotes, ratio, HWE exact p

VM:
    cd ~/repos/<checkout>/scripts/hla_popgen
    $PY 65c_65b_summary.py --out-dir ~/s05/results/65c
Python 3.8 compatible.
"""
import argparse
import importlib.util
import os
import sys

import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _THIS_DIR)
import _disclosure as _disc  # noqa: E402

_spec = importlib.util.spec_from_file_location("_tg65", os.path.join(_THIS_DIR, "65_trb_assembly_genotype.py"))
tg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tg)

FEATURES = ["TRBV4-3", "TRBV3-2"]


def rnd(n):
    """Denominator display: 0 exact, 1-19 '<20', else rounded to the nearest 10."""
    n = int(n)
    if n == 0:
        return "0"
    if n < _disc.SUPPRESS_BELOW:
        return f"<{_disc.SUPPRESS_BELOW}"
    return str(int(round(n / 10.0)) * 10)


def rate(k, n):
    r = _disc.mask_rate(k, n, _disc.PUBLIC, pct_decimals=1)
    return r["pct"], r["lo"], r["hi"]


def ok_copies(gene):
    d = gene[gene["status"] == "ok"][["person_id", "feature", "copies"]].copy()
    d["copies"] = pd.to_numeric(d["copies"], errors="coerce")
    return d[d["copies"].isin([0, 1, 2])]


def completion_table(sample, hap, asm_gene, people):
    d = sample[["person_id"]].merge(people, on="person_id", how="left")
    hap = hap.merge(people, on="person_id", how="left")
    hap["complete"] = tg._true(hap["complete"])
    called = set(ok_copies(asm_gene[asm_gene["feature"] == "TRBV4-3"])["person_id"])
    done = set(hap["person_id"])
    rows = []
    for kind, val, g in tg.strata(d):
        n = len(g)
        ng = int(g["person_id"].isin(done).sum())
        nc = int(g["person_id"].isin(called).sum())
        h = hap[hap["person_id"].isin(g["person_id"])]
        row = dict(stratum_type=kind, stratum=val, n_sampled=rnd(n), n_genotyped=rnd(ng))
        row["pct_genotyped"] = rate(ng, n)[0]
        row["pct_haplotypes_complete"] = rate(int(h["complete"].sum()), len(h))[0]
        row["pct_genotyped_trbv4_3_called"] = rate(nc, ng)[0]
        rows.append(row)
    return pd.DataFrame(rows)


def failures_table(fails, hap, sample):
    f = fails[fails["person_id"].isin(sample["person_id"])]
    ever = f.drop_duplicates("person_id", keep="last")["reason"].value_counts()
    still = f[~f["person_id"].isin(hap["person_id"])].drop_duplicates("person_id", keep="last")["reason"]
    still = still.value_counts()
    reasons = sorted(set(ever.index) | set(still.index))
    return pd.DataFrame([dict(reason=r, n_people_ever_failed=rnd(ever.get(r, 0)),
                              n_people_still_failed=rnd(still.get(r, 0))) for r in reasons],
                        columns=["reason", "n_people_ever_failed", "n_people_still_failed"])


def concordance_rows(reads_ok, asm_ok, people, genotyper):
    m = reads_ok.merge(asm_ok, on=["person_id", "feature"], suffixes=("_reads", "_asm"))
    m = m.merge(people, on="person_id", how="left")
    rows = []
    for feat in FEATURES:
        f = m[m["feature"] == feat]
        for kind, val, g in tg.strata(f):
            n = len(g)
            k = int((g["copies_reads"] == g["copies_asm"]).sum())
            pct, lo, hi = rate(k, n)
            row = dict(reads_genotyper=genotyper, feature=feat, stratum_type=kind, stratum=val,
                       n_both_called=rnd(n), pct_concordant=pct, lo=lo, hi=hi)
            for c in (0, 1, 2):
                rc, ac = g["copies_reads"] == c, g["copies_asm"] == c
                both = int((rc & ac).sum())
                row[f"n_reads{c}"] = rnd(int(rc.sum()))
                row[f"ppv{c}"] = rate(both, int(rc.sum()))[0]
                row[f"n_asm{c}"] = rnd(int(ac.sum()))
                row[f"sens{c}"] = rate(both, int(ac.sum()))[0]
            rows.append(row)
    return rows


def hwe_rows(calls, people, caller):
    d = calls.merge(people, on="person_id", how="left")
    rows = []
    for feat in FEATURES:
        f = d[d["feature"] == feat]
        for kind, val, g in tg.strata(f, platform=False):
            c0, c1, c2 = (int((g["copies"] == k).sum()) for k in (0, 1, 2))
            n = c0 + c1 + c2
            row = dict(caller=caller, feature=feat, ancestry=val, n_called=rnd(n))
            blank = n == 0 or any(0 < c < _disc.SUPPRESS_BELOW for c in (n, c0, c1, c2))
            if blank:
                row.update(insertion_af="", pct_0copies="", pct_het_obs="", pct_het_exp="",
                           het_obs_exp_ratio="", hwe_p="")
            else:
                p = (c1 + 2 * c2) / (2.0 * n)
                exp = 2 * p * (1 - p)
                row.update(insertion_af=round(p, 3), pct_0copies=round(100.0 * c0 / n, 1),
                           pct_het_obs=round(100.0 * c1 / n, 1), pct_het_exp=round(100 * exp, 1),
                           het_obs_exp_ratio=round(c1 / n / exp, 3) if exp > 0 else "",
                           hwe_p=float(f"{tg.hwe_exact_p(c0, c1, c2):.2g}"))
            rows.append(row)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reads-gene", default=os.path.join(tg.INTERNAL_ROOT, "65_reads", "per_person_gene.tsv"))
    ap.add_argument("--asm-dir", default=os.path.join(tg.INTERNAL_ROOT, "65_assembly"))
    ap.add_argument("--people-file", default=os.path.join(tg.INTERNAL_ROOT, "65_people_strat2000.tsv"))
    ap.add_argument("--cohort", default=tg.DEFAULT_COHORT)
    ap.add_argument("--min-reads", type=int, default=6)
    ap.add_argument("--af-cuts", default="0.2,0.8")
    ap.add_argument("--out-dir", default=os.path.join(tg.RESULTS_ROOT, "65c"))
    args = ap.parse_args(argv)

    people = tg._cohort_people(args.cohort)
    sample = pd.DataFrame(dict(person_id=sorted(tg.read_people_file(args.people_file))))
    ids = set(sample["person_id"])
    asm_gene = tg._read(os.path.join(args.asm_dir, "per_person_gene.tsv"), tg.GENE_COLS)
    asm_gene = asm_gene.drop_duplicates(["person_id", "feature"], keep="last")
    hap = tg._read(os.path.join(args.asm_dir, "per_haplotype.tsv"), tg.HAP_COLS)
    hap = hap.drop_duplicates(["person_id", "hap"], keep="last")
    hap = hap[hap["person_id"].isin(ids)]
    asm_gene = asm_gene[asm_gene["person_id"].isin(set(hap["person_id"]))]  # 'done' people only
    fails = tg._read(os.path.join(args.asm_dir, "failures.tsv"), tg.FAIL_COLS)

    reads = tg._read(args.reads_gene, tg.GENE_COLS).drop_duplicates(["person_id", "feature"], keep="last")
    cuts = tuple(float(x) for x in args.af_cuts.split(","))
    reads_by = {"stored": reads,
                "cuts": tg.rederive_reads_calls(reads, args.min_reads, cuts),
                "mixture": tg.mixture_reads_calls(reads, people, args.min_reads)[0]}
    asm_ok = ok_copies(asm_gene)
    asm_ok = asm_ok[asm_ok["feature"].isin(FEATURES)]

    os.makedirs(args.out_dir, exist_ok=True)
    w = lambda df, name: df.to_csv(os.path.join(args.out_dir, name), sep="\t", index=False)  # noqa: E731
    w(completion_table(sample, hap, asm_gene, people), "65c_completion.tsv")
    w(failures_table(fails, hap, sample), "65c_failures.tsv")
    conc = []
    for gname, r in reads_by.items():
        conc += concordance_rows(ok_copies(r[r["person_id"].isin(ids)]), asm_ok, people, gname)
    w(pd.DataFrame(conc), "65c_concordance.tsv")
    # HWE on one common set of people: called by the assembly and by both reads genotypers, per feature
    rc, rm = ok_copies(reads_by["cuts"]), ok_copies(reads_by["mixture"])
    key = ["person_id", "feature"]
    common = asm_ok[key].merge(rc[key], on=key).merge(rm[key], on=key)
    hw = []
    for caller, d in (("assembly", asm_ok), ("reads_cuts", rc), ("reads_mixture", rm)):
        hw += hwe_rows(d.merge(common, on=key), people, caller)
    w(pd.DataFrame(hw), "65c_hwe.tsv")
    tg.log(f"65c: {len(hap['person_id'].unique())} people done of {len(ids)}; tables in {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
