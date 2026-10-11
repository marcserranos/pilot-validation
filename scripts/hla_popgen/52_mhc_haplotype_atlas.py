#!/usr/bin/env python3
"""S05 axis-3 pilots on physically phased MHC haplotypes (existing outputs only).

  a. absent_cells.tsv / absent_pair_summary.tsv -- "forbidden-haplotype atlas". For gene pairs on
     ONE contig (copy 1, unrelated people), cross-tabulate allele lineages (field 1) per ancestry.
     A cell is "absent" when Observed = 0 and Expected (independence, n*pA*pB) is large enough that
     P(0 | Poisson(E)) = exp(-E) < 0.05 / n_cells_tested (Bonferroni over every cell tested in that
     pair x ancestry). LD alone makes many zero cells, so the gate is: absent in >=2 ancestries AND
     zero pooled. Positive controls that must come out: DQA1~DQB1 G1/G2 cross-group cells, and the
     DRB1 lineage -> DRB3/4/5 paralog rules (DRB1~DRBX pair, DRBX = which paralog sits on the hap).
  b. amr_source_decomposition.tsv / pooled_richness.tsv -- why AMR beats AFR at equal N. Equal-N
     samples of haplotypes from strict (p>=0.95) AFR, EUR, EAS, AMR; each AMR 2-field allele is
     classed by which reference pools also carry it. Then the richness of a synthetic 50/50
     EUR+AFR pool vs AFR alone at the same N: if mixing two sources beats the richer source,
     admixture-pooling is sufficient to explain the AMR excess.
  c. ase_power_census.tsv -- the non-coding "natural experiment" for allele-specific expression:
     per gene, unrelated people heterozygous for two PROTEIN-IDENTICAL alleles (same 2-field) that
     differ beyond the CDS / synonymously (field 3/4 differ, or non-coding/synonymous novel), vs
     identical (homozygous at full resolution) and protein-different; how many have RNA-seq.

Aggregate-only. Participant/haplotype counts 1-19 -> "<20"; allele-distinct counts exact.

    python3 52_mhc_haplotype_atlas.py --out-dir ~/s05/results/52
"""
import argparse
import os
from collections import Counter

import numpy as np
import pandas as pd

import _s05_common as C

ATLAS_PAIRS = [("A", "C"), ("A", "B"), ("B", "C"), ("C", "DRB1"), ("DRB1", "DQA1"),
               ("DRB1", "DQB1"), ("DQA1", "DQB1"), ("DPA1", "DPB1"), ("DQB1", "DPB1"),
               ("DRB1", "DRBX")]
CLASSICAL = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
DEFAULT_RNA_MANIFEST = os.path.expanduser(
    "~/workspace/vwb-aou-datasets-controlled-v9/v9/multiomics/rnaseq/manifest.tsv")


def field(consensus, k):
    """First k fields of a consensus allele, or None if any of them is 'new'/missing."""
    if not isinstance(consensus, str) or "*" not in consensus:
        return None
    gene, rest = C.norm_gene(consensus).split("*", 1)
    f = rest.split(":")
    if len(f) < k or any((x == "new" or x == "") for x in f[:k]):
        return None
    return f"{gene}*{':'.join(f[:k])}"


def hap_table(t1, unrelated):
    """One row per (person, hap): gene -> (contig, consensus) for copy 1, plus DRBX paralog."""
    sub = t1[(t1["copy_index"] == 1) & t1["person_id"].isin(unrelated)]
    recs = {}
    for pid, hap, gene, contig, cons in sub[["person_id", "hap", "gene", "contig", "consensus"]].itertuples(index=False):
        recs.setdefault((pid, hap), {})[gene] = (contig, cons)
    return recs


def drbx_state(gc):
    """Which DR paralog sits on the same contig as DRB1: DRB3/DRB4/DRB5/none/multi."""
    if "DRB1" not in gc:
        return None
    ctg = gc["DRB1"][0]
    present = [g for g in ("DRB3", "DRB4", "DRB5") if g in gc and gc[g][0] == ctg]
    # 'none' is only trustworthy if the contig also reaches DQA1 (spans the DR-DQ block)
    if not present:
        return "DRBX*none" if ("DQA1" in gc and gc["DQA1"][0] == ctg) else None
    return f"DRBX*{present[0][3:]}" if len(present) == 1 else "DRBX*multi"


def atlas(recs, ancestry_of, mode, level=1):
    rows, summ = [], []
    for ga, gb in ATLAS_PAIRS:
        obs = []
        for (pid, hap), gc in recs.items():
            if ga not in gc:
                continue
            if gb == "DRBX":
                b = drbx_state(gc)
            else:
                if gb not in gc or gc[gb][0] != gc[ga][0]:
                    continue
                b = field(gc[gb][1], level)
            a = field(gc[ga][1], level)
            if a and b:
                obs.append((ancestry_of.get(pid, "NA"), a, b))
        df = pd.DataFrame(obs, columns=["ancestry", "a", "b"])
        strata = [("POOLED", df)] + [(an, g) for an, g in df.groupby("ancestry") if an in C.ANCESTRY_ORDER]
        absent_by = {}
        for an, g in strata:
            n = len(g)
            if n < 200:
                continue
            pa = g["a"].value_counts() / n
            pb = g["b"].value_counts() / n
            oc = Counter(zip(g["a"], g["b"]))
            cells = [(x, y, n * pa[x] * pb[y]) for x in pa.index for y in pb.index]
            tested = [c for c in cells if c[2] >= 1]
            thr = np.log(max(len(tested), 1) / 0.05)  # exp(-E) < 0.05/m  <=>  E > ln(m/0.05)
            absent = [(x, y, e) for x, y, e in tested if oc.get((x, y), 0) == 0 and e > thr]
            summ.append({"pair": f"{ga}~{gb}", "level": f"field{level}", "ancestry": an,
                         "n_haps": C.mask_n(n, mode), "n_cells_tested_E_ge_1": len(tested),
                         "E_threshold": round(thr, 2), "n_absent_sig": len(absent),
                         "n_zero_cells_E_ge_5": sum(1 for x, y, e in tested if e >= 5 and oc.get((x, y), 0) == 0)})
            for x, y, e in absent:
                absent_by.setdefault((x, y), {})[an] = e
        for (x, y), d in absent_by.items():
            n_anc = sum(1 for k in d if k != "POOLED")
            if "POOLED" in d and n_anc >= 2:
                rows.append({"pair": f"{ga}~{gb}", "level": f"field{level}", "allele_a": x, "allele_b": y,
                             "E_pooled": round(d["POOLED"], 1), "n_ancestries_absent": n_anc,
                             "ancestries": ",".join(sorted(k for k in d if k != "POOLED")),
                             "observed_pooled": 0})
    return pd.DataFrame(rows), pd.DataFrame(summ)


def richness(alleles):
    return len(set(a for a in alleles if a))


def amr_decomposition(recs, cohort, n_boot=50, seed=7):
    strict = cohort[cohort["p_max"] >= 0.95]
    pools = {a: set(strict.loc[strict["ancestry_pred"] == a, "person_id"]) for a in ["AFR", "EUR", "EAS", "AMR"]}
    by_anc = {a: [gc for (pid, hap), gc in recs.items() if pid in s] for a, s in pools.items()}
    n_eq = min(len(v) for v in by_anc.values())
    rng = np.random.default_rng(seed)
    dec, rich = [], []
    for gene in CLASSICAL:
        acc = Counter()
        rr = {k: [] for k in ["AFR", "EUR", "EAS", "AMR", "MIX_EUR50_AFR50", "MIX_EUR45_AFR10_EAS45"]}
        for _ in range(n_boot):
            samp = {a: [field(by_anc[a][i].get(gene, (None, None))[1], 2)
                        for i in rng.choice(len(by_anc[a]), n_eq, replace=False)] for a in by_anc}
            ref = {a: set(samp[a]) - {None} for a in ["AFR", "EUR", "EAS"]}
            for al in set(samp["AMR"]) - {None}:
                hit = [a for a in ref if al in ref[a]]
                acc["private_to_AMR" if not hit else ("only_" + hit[0] if len(hit) == 1 else "shared_" + "+".join(sorted(hit)))] += 1
            for a in ["AFR", "EUR", "EAS", "AMR"]:
                rr[a].append(richness(samp[a]))
            half = n_eq // 2
            rr["MIX_EUR50_AFR50"].append(richness(samp["EUR"][:half] + samp["AFR"][:n_eq - half]))
            k1, k2 = int(.45 * n_eq), int(.10 * n_eq)
            rr["MIX_EUR45_AFR10_EAS45"].append(
                richness(samp["EUR"][:k1] + samp["AFR"][:k2] + samp["EAS"][:n_eq - k1 - k2]))
        tot = sum(acc.values())
        for k, v in sorted(acc.items()):
            dec.append({"gene": gene, "class": k, "pct_of_AMR_alleles": round(100 * v / tot, 1)})
        for k, v in rr.items():
            rich.append({"gene": gene, "pool": k, "n_haps_equal": n_eq, "richness_mean": round(np.mean(v), 1),
                         "richness_sd": round(np.std(v), 1)})
    return pd.DataFrame(dec), pd.DataFrame(rich)


def ase_census(t1, unrelated, rna_ids, mode):
    sub = t1[(t1["copy_index"] == 1) & t1["gene"].isin(CLASSICAL) & t1["person_id"].isin(unrelated)]
    rows = []
    for gene, g in sub.groupby("gene"):
        cnt = Counter()
        for pid, gg in g.groupby("person_id"):
            if len(gg) != 2:
                continue
            c1, c2 = list(gg["consensus"])
            d1, d2 = list(gg["novelty_depth"])
            if not (isinstance(c1, str) and isinstance(c2, str)):
                continue
            p1, p2 = field(c1, 2), field(c2, 2)
            if p1 is None or p2 is None or (d1 == 2) or (d2 == 2):
                cls = "protein_unresolved_or_novel"
            elif p1 != p2:
                cls = "protein_different"
            elif c1 == c2 and not (d1 >= 3 or d2 >= 3):
                cls = "identical_full_resolution"
            else:
                cls = "protein_identical_seq_different"
            cnt[(cls, pid in rna_ids)] += 1
        for cls in ["identical_full_resolution", "protein_identical_seq_different", "protein_different",
                    "protein_unresolved_or_novel"]:
            n_all = cnt[(cls, True)] + cnt[(cls, False)]
            rows.append({"gene": gene, "pair_class": cls, "n_people": C.mask_n(n_all, mode),
                         "n_with_rnaseq": C.mask_n(cnt[(cls, True)], mode)})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=C.DEFAULT_TABLE1)
    ap.add_argument("--cohort-membership", default=C.DEFAULT_COHORT)
    ap.add_argument("--relatedness-table", default=C.DEFAULT_RELATEDNESS)
    ap.add_argument("--rna-manifest", default=DEFAULT_RNA_MANIFEST)
    ap.add_argument("--out-dir", default=os.path.expanduser("~/s05/results/52"))
    ap.add_argument("--disclosure", default="public", choices=["public", "internal"])
    args = ap.parse_args()
    mode = args.disclosure
    cohort = C.load_cohort(args.cohort_membership)
    t1 = C.load_table1(args.table1)
    unrelated = C.greedy_unrelated(set(t1["person_id"]), C.load_relatedness_pairs(args.relatedness_table))
    ancestry_of = dict(zip(cohort["person_id"], cohort["ancestry_pred"]))
    recs = hap_table(t1, unrelated)
    C.log(f"haplotypes: {len(recs)}")
    out = os.path.expanduser(args.out_dir)
    for level in (1, 2):
        cells, summ = atlas(recs, ancestry_of, mode, level)
        C.write_tsv(cells, f"{out}/absent_cells_field{level}.tsv")
        C.write_tsv(summ, f"{out}/absent_pair_summary_field{level}.tsv")
    dec, rich = amr_decomposition(recs, cohort)
    C.write_tsv(dec, f"{out}/amr_source_decomposition.tsv")
    C.write_tsv(rich, f"{out}/pooled_richness.tsv")
    rna_ids = set()
    if os.path.exists(args.rna_manifest):
        rna_ids = set(pd.read_csv(args.rna_manifest, sep="\t", dtype=str)["research_id"])
    C.log(f"RNA-seq research_ids: {len(rna_ids)}")
    C.write_tsv(ase_census(t1, unrelated, rna_ids, mode), f"{out}/ase_power_census.tsv")
    C.log("DONE 52")


if __name__ == "__main__":
    main()
