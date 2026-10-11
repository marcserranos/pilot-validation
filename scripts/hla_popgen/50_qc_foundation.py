#!/usr/bin/env python3
"""S05 axis-1/2 QC foundation for HLA: four aggregate tables that the P0/P1 QC asks need.

  a. waterfall.tsv          cohort exclusion waterfall (N per ancestry after each filter), with the
                            three ancestry-purity thresholds used across the repo (0.90/0.95/0.98)
                            side by side, so one exclusion list can be agreed.
  b. same_contig_matrix.tsv gene x gene "physical phasing level": among haplotypes where both genes
                            are called (copy 1), the fraction with both on ONE contig. By platform.
     mhc_contigs_per_hap.tsv number of distinct contigs carrying MHC gene calls per haplotype.
  c. novelty_by_platform.tsv per gene x ancestry x platform novelty rates (any / protein / syn /
                            non-coding), + a Mantel-Haenszel odds ratio revio vs sequel2e pooled over
                            ancestry strata (platform is confounded with ancestry and with time).
  d. protein_novel_waterfall.tsv per gene, distinct protein-novel (field-2 'new') sequences that
                            survive: >=2 unrelated carriers -> not homopolymer-indel-only -> no
                            template_warning -> seen on >=2 platforms -> seen in >=2 ancestries.
                            Distinct-sequence key = (gene, cds_mut): identical CDS sequences align
                            identically to the same refdata, so the key is a deterministic proxy for
                            nucleotide identity (it can only split, never merge, true recurrences).
  e. dq_ascertainment.tsv  DQA1/DQB1 contig-break rate by G1/G2 class (is the same-contig-only
                            "0 cross-group" result biased by which haplotypes break?), plus the
                            random-pairing test on split haplotypes.

Aggregate-only outputs. Participant counts 1-19 -> "<20"; distinct-allele counts exact (PUBLIC rule).

    python3 50_qc_foundation.py --out-dir ~/s05/results/50
    python3 50_qc_foundation.py --synthetic --out-dir /tmp/s05_50   # local self-test
"""
import argparse
import os
import re
from collections import defaultdict
from itertools import combinations

import numpy as np
import pandas as pd

import _s05_common as C

MHC_GENES = ["A", "C", "B", "E", "F", "G", "MICA", "MICB", "C4A", "C4B", "DRA", "DRB1", "DRB3",
             "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1", "TAP1", "TAP2", "DMA", "DMB", "DOA",
             "DOB"]
CLASSICAL = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1", "DRB3", "DRB4", "DRB5"]
PURITY = [0.90, 0.95, 0.98]
CS_TOKEN_RE = re.compile(r"(:\d+|\*[a-z]{2}|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})", re.I)


def is_homopolymer_indel_only(cds_mut):
    """Same rule as 03_novel_alleles.is_homopolymer_indel_only (copied to keep deploy flat)."""
    if not isinstance(cds_mut, str) or not cds_mut:
        return False
    any_indel = False
    for entry in cds_mut.split(","):
        parts = entry.split("|")
        if len(parts) < 2:
            continue
        for tok in CS_TOKEN_RE.findall(parts[1]):
            if tok.startswith("*"):
                return False
            if tok[0] in "+-":
                any_indel = True
                if len(set(tok[1:].upper())) != 1:
                    return False
    return any_indel


# --------------------------------------------------------------------------- a. waterfall
def waterfall(cohort, t1_people, unrelated, mode):
    rows = []
    c = cohort.copy()
    steps = [
        ("all_cohort_membership_rows", np.ones(len(c), bool)),
        ("in_lr", c["in_lr"].astype(str).str.lower().isin(["true", "1"])),
    ]
    m = steps[-1][1]
    m = m & c["platform"].isin(C.ASSEMBLY_PLATFORMS)
    steps.append(("assembly_platform", m))
    m = m & c["person_id"].isin(t1_people)
    steps.append(("has_hla_calls", m))
    m = m & c["person_id"].isin(unrelated)
    steps.append(("unrelated_kin_lt_0.0442", m))
    for thr in PURITY:
        steps.append((f"unrelated_and_ancestry_p_ge_{thr:.2f}", m & (c["p_max"] >= thr)))
    for name, mask in steps:
        sub = c[mask]
        row = {"step": name, "total": C.mask_n(len(sub), mode)}
        for a in C.ANCESTRY_ORDER:
            row[a] = C.mask_n(int((sub["ancestry_pred"] == a).sum()), mode)
        for p in C.ASSEMBLY_PLATFORMS:
            row[f"platform_{p}"] = C.mask_n(int((sub["platform"] == p).sum()), mode)
        rows.append(row)
    # where do in_lr people without HLA calls come from (platform / trim tier)?
    lost = c[steps[2][1] & ~c["person_id"].isin(t1_people)]
    for p in C.ASSEMBLY_PLATFORMS:
        rows.append({"step": f"lost_no_hla_calls_platform_{p}",
                     "total": C.mask_n(int((lost["platform"] == p).sum()), mode)})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- b. same-contig
def hap_gene_contigs(t1):
    """{(pid,hap): {gene: contig}} for copy_index 1 of MHC genes."""
    sub = t1[(t1["copy_index"] == 1) & t1["gene"].isin(MHC_GENES)]
    out = defaultdict(dict)
    for pid, hap, gene, contig in sub[["person_id", "hap", "gene", "contig"]].itertuples(index=False):
        out[(pid, hap)][gene] = contig
    return out


def same_contig_matrix(hgc, platform_of, mode):
    rows = []
    strata = {"all": None, "revio": "revio", "sequel2e": "sequel2e", "sequel2": "sequel2"}
    for sname, plat in strata.items():
        both = defaultdict(int)
        same = defaultdict(int)
        for (pid, hap), gc in hgc.items():
            if plat and platform_of.get(pid) != plat:
                continue
            genes = [g for g in MHC_GENES if g in gc]
            for g1, g2 in combinations(genes, 2):
                both[(g1, g2)] += 1
                same[(g1, g2)] += gc[g1] == gc[g2]
        for (g1, g2), d in both.items():
            r = C.rate_cells(same[(g1, g2)], d, mode, digits=1)
            rows.append({"stratum": sname, "gene_a": g1, "gene_b": g2, "n_haps_both_called": r["d"],
                         "n_same_contig": r["n"], "pct_same_contig": r["pct"],
                         "lo": r["lo"], "hi": r["hi"]})
    return pd.DataFrame(rows)


def contigs_per_hap(hgc, platform_of, ancestry_of, mode):
    rec = []
    for (pid, hap), gc in hgc.items():
        rec.append((platform_of.get(pid, "NA"), ancestry_of.get(pid, "NA"), len(set(gc.values()))))
    df = pd.DataFrame(rec, columns=["platform", "ancestry", "n_contigs"])
    df["n_contigs_bin"] = df["n_contigs"].clip(upper=5).astype(str).replace({"5": "5+"})
    rows = []
    for by in ["platform", "ancestry"]:
        for key, g in df.groupby(by):
            tot = len(g)
            for b, gb in g.groupby("n_contigs_bin"):
                r = C.rate_cells(len(gb), tot, mode, digits=1)
                rows.append({"by": by, "stratum": key, "n_contigs": b, "n_haps": r["n"],
                             "of": r["d"], "pct": r["pct"]})
            rows.append({"by": by, "stratum": key, "n_contigs": "mean",
                         "n_haps": "", "of": C.mask_n(tot, mode), "pct": f"{g['n_contigs'].mean():.2f}"})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- c. platform
def novelty_by_platform(t1, platform_of, ancestry_of, unrelated, mode):
    sub = t1[t1["gene"].isin(CLASSICAL) & t1["person_id"].isin(unrelated)].copy()
    sub = sub[sub["consensus"].notna() & (sub["consensus"] != "undetermined")]
    sub["platform"] = sub["person_id"].map(platform_of)
    sub["ancestry"] = sub["person_id"].map(ancestry_of)
    sub["nov_any"] = sub["novelty_depth"].notna() & (sub["novelty_depth"] >= 2)
    sub["nov_protein"] = sub["novelty_depth"] == 2
    sub["nov_syn"] = sub["novelty_depth"] == 3
    sub["nov_noncoding"] = sub["novelty_depth"] == 4
    rows, mh = [], []
    for gene, g in sub.groupby("gene"):
        for (anc, plat), gg in g.groupby(["ancestry", "platform"]):
            row = {"gene": gene, "ancestry": anc, "platform": plat}
            for k in ["nov_any", "nov_protein", "nov_syn", "nov_noncoding"]:
                r = C.rate_cells(int(gg[k].sum()), len(gg), mode)
                row[f"{k}_pct"] = r["pct"]
                row[f"{k}_ci"] = f"{r['lo']}-{r['hi']}" if r["pct"] else ""
            row["n_calls"] = C.mask_n(len(gg), mode)
            rows.append(row)
        # Mantel-Haenszel OR, revio (exposed) vs sequel2e, strata = ancestry
        for k in ["nov_any", "nov_protein", "nov_noncoding"]:
            num = den = 0.0
            n_str = 0
            for anc, ga in g.groupby("ancestry"):
                r = ga[ga["platform"] == "revio"][k]
                s = ga[ga["platform"] == "sequel2e"][k]
                if len(r) < 20 or len(s) < 20:
                    continue
                a_, b_ = r.sum(), len(r) - r.sum()
                c_, d_ = s.sum(), len(s) - s.sum()
                n = a_ + b_ + c_ + d_
                num += a_ * d_ / n
                den += b_ * c_ / n
                n_str += 1
            mh.append({"gene": gene, "outcome": k, "n_ancestry_strata": n_str,
                       "mh_or_revio_vs_sequel2e": f"{num/den:.3f}" if den > 0 else ""})
    return pd.DataFrame(rows), pd.DataFrame(mh)


# --------------------------------------------------------------------------- d. waterfall
def protein_novel_waterfall(t1, platform_of, ancestry_of, unrelated):
    sub = t1[(t1["novelty_depth"] == 2) & t1["person_id"].isin(unrelated)].copy()
    sub = sub[sub["cds_mut"].notna()]
    sub["key"] = sub["gene"] + "|" + sub["cds_mut"]
    rows = []
    for gene, g in sub.groupby("gene"):
        agg = g.groupby("key").agg(
            carriers=("person_id", "nunique"),
            platforms=("person_id", lambda s: len({platform_of.get(p) for p in s})),
            ancestries=("person_id", lambda s: len({ancestry_of.get(p) for p in s})),
            warn=("template_warning", lambda s: s.notna().any() and (s.astype(str) != "nan").any()),
            cds_mut=("cds_mut", "first"))
        agg["homopolymer"] = agg["cds_mut"].map(is_homopolymer_indel_only)
        s0 = len(agg)
        f1 = agg[agg["carriers"] >= 2]
        f2 = f1[~f1["homopolymer"]]
        f3 = f2[~f2["warn"]]
        f4 = f3[f3["platforms"] >= 2]
        f5 = f4[f4["ancestries"] >= 2]
        f6 = f3[f3["carriers"] >= 20]
        rows.append({"gene": gene, "distinct_protein_novel": s0, "ge2_unrelated_carriers": len(f1),
                     "not_homopolymer_only": len(f2), "no_template_warning": len(f3),
                     "ge2_platforms": len(f4), "ge2_ancestries": len(f5),
                     "ge20_carriers_clean": len(f6)})
    df = pd.DataFrame(rows)
    if len(df):
        tot = df.drop(columns="gene").sum()
        tot["gene"] = "ALL"
        df = pd.concat([df, pd.DataFrame([tot])], ignore_index=True)
    return df


# --------------------------------------------------------------------------- e. DQ ascertainment
def dq_group(gene, consensus):
    """Petersdorf 2022 rule, as in 37_dq_g1g2_signed_ld.dq_group: G1 = DQA1*02-06 / DQB1*02-04,
    G2 = DQA1*01 / DQB1*05-06. Uses field 1 only, so novel field-2+ calls still classify."""
    if not isinstance(consensus, str) or "*" not in consensus:
        return None
    f1 = C.norm_gene(consensus).split("*", 1)[1].split(":")[0]
    if not f1.isdigit():
        return None
    f1 = int(f1)
    if gene == "DQA1":
        return "G2" if f1 == 1 else ("G1" if 2 <= f1 <= 6 else None)
    if gene == "DQB1":
        return "G1" if 2 <= f1 <= 4 else ("G2" if f1 in (5, 6) else None)
    return None


def dq_ascertainment(t1, platform_of, ancestry_of, unrelated, mode):
    """Does the chance that DQA1 and DQB1 land on different contigs depend on their G-group?
    If cross-group haplotypes were real but tended to break contigs, the same-contig-only
    '0 cross-group' result would be biased. Also: on split haplotypes the hap label across
    contigs is arbitrary (non-trio hifiasm), so cross-group pairings there should look random."""
    sub = t1[(t1["copy_index"] == 1) & t1["gene"].isin(["DQA1", "DQB1"]) & t1["person_id"].isin(unrelated)]
    w = sub.pivot_table(index=["person_id", "hap"], columns="gene",
                        values=["contig", "consensus"], aggfunc="first")
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w = w.dropna(subset=["contig_DQA1", "contig_DQB1"]).reset_index()
    w["ga"] = [dq_group("DQA1", c) for c in w["consensus_DQA1"]]
    w["gb"] = [dq_group("DQB1", c) for c in w["consensus_DQB1"]]
    w = w.dropna(subset=["ga", "gb"])
    w["split"] = w["contig_DQA1"] != w["contig_DQB1"]
    w["ancestry"] = w["person_id"].map(ancestry_of)
    w["platform"] = w["person_id"].map(platform_of)
    rows = []
    for scope, d in [("ALL", w)] + [(f"anc_{a}", g) for a, g in w.groupby("ancestry")] + \
                    [(f"plat_{p}", g) for p, g in w.groupby("platform")]:
        for (ga, gb), g in d.groupby(["ga", "gb"]):
            r = C.rate_cells(int(g["split"].sum()), len(g), mode, 2)
            rows.append({"scope": scope, "DQA1_group": ga, "DQB1_group": gb, "n_haps": r["d"],
                         "n_split": r["n"], "pct_split": r["pct"], "lo": r["lo"], "hi": r["hi"]})
        # break rate by DQA1 group alone (ascertainment depends only on what is observable)
        for ga, g in d.groupby("ga"):
            r = C.rate_cells(int(g["split"].sum()), len(g), mode, 2)
            rows.append({"scope": scope, "DQA1_group": ga, "DQB1_group": "any", "n_haps": r["d"],
                         "n_split": r["n"], "pct_split": r["pct"], "lo": r["lo"], "hi": r["hi"]})
    # Random-pairing expectation on split haps: P(cross) = pA1*pB2 + pA2*pB1 within split haps.
    sp = w[w["split"]]
    if len(sp):
        pa1, pb1 = (sp["ga"] == "G1").mean(), (sp["gb"] == "G1").mean()
        exp_cross = pa1 * (1 - pb1) + (1 - pa1) * pb1
        obs_cross = int((sp["ga"] != sp["gb"]).sum())
        rows.append({"scope": "split_haps_cross_group", "DQA1_group": "obs_vs_random_pairing",
                     "DQB1_group": f"expected_pct={100*exp_cross:.1f}",
                     "n_haps": C.mask_n(len(sp), mode), "n_split": C.mask_n(obs_cross, mode),
                     "pct_split": C.rate_cells(obs_cross, len(sp), mode, 1)["pct"]})
        # Within-person relabel test: for people with BOTH haps split, does swapping DQB1 between
        # haps make both haplotypes group-concordant? (tests 'labels arbitrary across contigs')
        both = sp.groupby("person_id").filter(lambda g: len(g) == 2)
        fix = 0
        tot = 0
        for pid, g in both.groupby("person_id"):
            a = list(g["ga"])
            b = list(g["gb"])
            disc_now = (a[0] != b[0]) + (a[1] != b[1])
            if disc_now == 0:
                continue
            tot += 1
            fix += (a[0] == b[1]) and (a[1] == b[0])
        rows.append({"scope": "both_haps_split_people_with_cross", "DQA1_group": "swap_DQB1_resolves",
                     "DQB1_group": "", "n_haps": C.mask_n(tot, mode), "n_split": C.mask_n(fix, mode),
                     "pct_split": C.rate_cells(fix, tot, mode, 1)["pct"]})
    same = w[~w["split"]]
    rows.append({"scope": "same_contig_haps_cross_group", "DQA1_group": "observed", "DQB1_group": "",
                 "n_haps": C.mask_n(len(same), mode),
                 "n_split": C.mask_n(int((same["ga"] != same["gb"]).sum()), mode), "pct_split": ""})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- synthetic
def make_synthetic(n=600, seed=1):
    rng = np.random.default_rng(seed)
    pids = [f"P{i:05d}" for i in range(n)]
    anc = rng.choice(C.ANCESTRY_ORDER, n, p=[.2, .15, .1, .45, .03, .07])
    plat = rng.choice(["revio", "sequel2e", "sequel2"], n, p=[.8, .12, .08])
    coh = pd.DataFrame({"person_id": pids, "in_lr": True, "in_sr": True,
                        "ancestry_pred": [a.lower() for a in anc], "platform": plat,
                        "trim_tier": "paf_region"})
    for a in C.ANCESTRY_ORDER:
        coh[f"p_{a.lower()}"] = np.where(anc == a, rng.uniform(.85, 1, n), rng.uniform(0, .05, n))
    rows = []
    for p in pids[:-10]:
        for hap in ["hap1", "hap2"]:
            nctg = rng.choice([1, 2, 3], p=[.3, .5, .2])
            for i, gene in enumerate(MHC_GENES):
                ctg = f"{hap}_ctg{min(i * nctg // len(MHC_GENES), nctg - 1)}"
                nd = rng.choice([np.nan, 2, 3, 4], p=[.85, .03, .04, .08])
                cons = f"HLA-{gene}*01:{'new' if nd == 2 else '01'}"
                rows.append({"person_id": p, "hap": hap, "contig": ctg, "gene": f"HLA-{gene}",
                             "copy_index": 1, "consensus": cons, "novelty_depth": nd,
                             "template_distance": 0, "gene_start": i * 1000, "gene_end": i * 1000 + 900,
                             "cds_mut": f"{gene}*01:01|:10*ag:20|x{rng.integers(0, 30)}" if nd == 2 else np.nan,
                             "template_warning": np.nan})
    return coh, pd.DataFrame(rows), []


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=C.DEFAULT_TABLE1)
    ap.add_argument("--cohort-membership", default=C.DEFAULT_COHORT)
    ap.add_argument("--relatedness-table", default=C.DEFAULT_RELATEDNESS)
    ap.add_argument("--out-dir", default=os.path.expanduser("~/s05/results/50"))
    ap.add_argument("--disclosure", default="public", choices=["public", "internal"])
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--synthetic", action="store_true")
    args = ap.parse_args()
    mode = args.disclosure

    if args.synthetic:
        cohort, t1, pairs = make_synthetic()
        cohort["ancestry_pred"] = cohort["ancestry_pred"].str.upper()
        cohort["p_max"] = cohort[C.PROB_COLS].max(axis=1)
        t1["gene"] = t1["gene"].map(C.norm_gene)
    else:
        cohort = C.load_cohort(args.cohort_membership)
        t1 = C.load_table1(args.table1, args.limit)
        pairs = C.load_relatedness_pairs(args.relatedness_table)
    C.log(f"table1 rows={len(t1)} people={t1['person_id'].nunique()}")
    t1_people = set(t1["person_id"].unique())
    unrelated = C.greedy_unrelated(t1_people, pairs)
    C.log(f"unrelated among HLA-called: {len(unrelated)}")
    platform_of = dict(zip(cohort["person_id"], cohort["platform"]))
    ancestry_of = dict(zip(cohort["person_id"], cohort["ancestry_pred"]))

    out = os.path.expanduser(args.out_dir)
    C.write_tsv(waterfall(cohort, t1_people, unrelated, mode), f"{out}/waterfall.tsv")
    hgc = hap_gene_contigs(t1)
    C.write_tsv(same_contig_matrix(hgc, platform_of, mode), f"{out}/same_contig_matrix.tsv")
    C.write_tsv(contigs_per_hap(hgc, platform_of, ancestry_of, mode), f"{out}/mhc_contigs_per_hap.tsv")
    nbp, mh = novelty_by_platform(t1, platform_of, ancestry_of, unrelated, mode)
    C.write_tsv(nbp, f"{out}/novelty_by_platform.tsv")
    C.write_tsv(mh, f"{out}/novelty_platform_mh_or.tsv")
    C.write_tsv(protein_novel_waterfall(t1, platform_of, ancestry_of, unrelated),
                f"{out}/protein_novel_waterfall.tsv")
    C.write_tsv(dq_ascertainment(t1, platform_of, ancestry_of, unrelated, mode),
                f"{out}/dq_ascertainment.tsv")
    C.log("DONE 50")


if __name__ == "__main__":
    main()
