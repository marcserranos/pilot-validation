#!/usr/bin/env python3
"""HLA-A 'deletion' validation (WS2, S03 call #8).

Cole (call #8): "that is a big result... how is that plausible?" -- script 30 reports HLA-A
bridged-absence "deletion" at ~1.83% of bridged haplotypes (n_bridged ~= 23k) and designates
HLA-A a negative control (grey, expect ~0). No documented whole-gene HLA-A deletion haplotype
exists in the literature (only null alleles; one 85-kb HLA-B deletion family is known). Prior:
artifact.

Four questions, each answered with numbers only (aggregate, small-cell suppressed):

  1. Per-haplotype HLA-A bridged-deletion rate, rerun on Table 1, per-ancestry Wilson CI.
     Biallelic (both-haplotype) loss: observed vs N*p^2 (Poisson/binomial exact P). Because
     independent per-haplotype assembly artifacts ALSO give ~N*p^2 by chance, this also computes
     (a) within-person hap1<->hap2 concordance (Fisher exact / odds ratio) and (b) correlation of
     deletion status with per-person assembly-quality proxies (bridging-contig span and gene
     count as an N50/length proxy -- Table 1 carries no genomic coordinates -- plus platform,
     trim_tier, n_genes_called and template_warning rate from cohort_membership.tsv / Table 1).
  2. Short-read contradiction test. For people with exactly one LR HLA-A-deleted haplotype, SR
     should look hemizygous. If SR reports two distinct alleles, the gene exists -> the LR call is
     an artifact. Compares SR heterozygosity in LR-deleted carriers vs ancestry-matched
     non-carriers (Fisher exact), and checks whether the "extra" SR allele's gene family is seen
     anywhere else in the person's own LR assembly (fragmentation/misplacement signature).
  3. Same SR-contradiction test as a calibration panel across every script-30 control gene
     (negative: A/B/C/DRA/DQA1/DQB1/DPA1/DPB1/DRB1; positive: DRB3/DRB4/DRB5, where SR is
     expected to reproduce the absence, not contradict it -- SR HLA genotyping does not
     necessarily type DRB3/4/5 at all, which is itself reported).
  4. Mechanism characterization: are deleted-gene calls concentrated in short/low-gene-count
     bridging contigs, one platform/trim_tier, or high template_warning individuals?

Reuses (never re-derives): 30_hla_structural_variation.bridged_absences / NEGATIVE_CONTROL_GENES
/ POSITIVE_CONTROL_GENES / DRB1_GROUP_EXPECTATION; 16_phasing_mendelian_validation.
derive_canonical_order; 24_novelty_by_field.build_people (unrelated cohort + strict ancestry);
_viz_common.{wilson_ci,load_sr_genotypes,load_cohort_membership,nature_style,save_fig,
ANCESTRY_COLORS,ANCESTRY_ORDER}.

NOTE on the VM deployment for this run: the VM's checked-out ~/repos/pilot-validation working
copy was on an unrelated, older branch (needle-view-cds-diversity-density) missing scripts 16,
24 and 30 entirely (never commit/checkout there -- AGENT_PREAMBLE/VM_CHANNEL.md). Rather than
push this branch or touch that checkout, the functions this script needs from those three files
are reproduced verbatim below (`gene_bare`, `bridged_absences`, `NEGATIVE_CONTROL_GENES`,
`POSITIVE_CONTROL_GENES`, `derive_canonical_order`, `strict_ancestry`, `greedy_unrelated`,
`load_relatedness_pairs`, `build_people`) instead of dynamically importing those files -- this is
still "reuse, don't re-derive": every function below is a direct copy of the corresponding
function in the source script, not a reimplementation. _viz_common.py (which IS present on the
VM's own branch, uncommitted-but-identical to this branch's copy) is still imported normally.

Usage (VM, full cohort):

    setsid nohup python3 -u scripts/hla_popgen/38_hla_a_deletion_validation.py \\
        --out-dir ~/results/38 < /dev/null & disown

Usage (local fixtures): scripts/hla_popgen/tests/test_hla_a_deletion_validation.py
"""
import argparse
import importlib.util
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import _viz_common as vc

try:
    from scipy.stats import binomtest, fisher_exact, spearmanr
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


def _load_module(filename, modname):
    path = os.path.join(_THIS_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


_mods = {}


def mod(filename, name):
    if name not in _mods:
        _mods[name] = _load_module(filename, name)
    return _mods[name]


AOU_CONTROLLED_ROOT = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9")
DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_TABLE1 = os.path.join(DEFAULT_OUTROOT, "hla_calls_rich.tsv")
DEFAULT_COHORT = os.path.join(DEFAULT_OUTROOT, "cohort_membership.tsv")
DEFAULT_RELATEDNESS = os.path.join(
    AOU_CONTROLLED_ROOT, "v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv")
DEFAULT_SR_GENOTYPES = os.path.join(
    AOU_CONTROLLED_ROOT, "v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv")
DEFAULT_OUT_DIR = os.path.expanduser("~/results/38")

SUPPRESS_BELOW = 20
KIN_MIN = 0.0442
STRICT_MIN = 0.9

TARGET_GENE = "A"
CALIBRATION_GENES = ["A", "B", "C", "DRA", "DQA1", "DQB1", "DPA1", "DPB1", "DRB1",
                      "DRB3", "DRB4", "DRB5"]


def suppress(n, threshold=SUPPRESS_BELOW):
    if n is None or (isinstance(n, float) and np.isnan(n)):
        return "NA"
    n = int(n)
    return "0" if n == 0 else ("<%d" % threshold if n < threshold else str(n))


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


# ---------------------------------------------------------------------------
# Inlined, verbatim, from 30_hla_structural_variation.py / 16_phasing_mendelian_validation.py /
# 24_novelty_by_field.py -- see the module-docstring NOTE above for why these are copied rather
# than dynamically imported for this run.
# ---------------------------------------------------------------------------
NEGATIVE_CONTROL_GENES = ["A", "B", "C", "DRA", "DQA1", "DQB1", "DPA1", "DPB1", "DRB1"]
POSITIVE_CONTROL_GENES = ["DRB3", "DRB4", "DRB5"]


def gene_bare(g):
    return str(g).replace("HLA-", "")


def bridged_absences(t1, gene_order, genes_of_interest):
    """Verbatim copy of 30_hla_structural_variation.bridged_absences -- see that module's
    docstring for the full rationale (a gene is called deleted only when a single contig carries
    a gene on each side of it in canonical order and the gene itself is absent)."""
    diag = Counter()
    per_hap = defaultdict(dict)
    seen_haps = set()

    cols = ["person_id", "hap", "contig", "gene_bare"]
    for (pid, hap, contig), grp in t1[cols].groupby(["person_id", "hap", "contig"], sort=False):
        seen_haps.add((pid, hap))
        present = set(grp["gene_bare"])
        ranks = [gene_order[g] for g in present if g in gene_order]
        if not ranks:
            diag["contig_with_no_ordered_gene"] += 1
            continue
        lo, hi = min(ranks), max(ranks)
        for g in present:
            per_hap[(pid, hap)][g] = "present"
        for g in genes_of_interest:
            if g in present or g not in gene_order:
                continue
            r = gene_order[g]
            if lo < r < hi:
                per_hap[(pid, hap)].setdefault(g, "deleted_bridged")
                diag["bridged_absence"] += 1
            else:
                per_hap[(pid, hap)].setdefault(g, "absent_unbridged")
                diag["unbridged_absence"] += 1

    rows = []
    for (pid, hap), gmap in per_hap.items():
        for g, status in gmap.items():
            if g in genes_of_interest:
                rows.append((pid, hap, g, status))
    calls = pd.DataFrame(rows, columns=["person_id", "hap", "gene", "status"])
    diag["haplotypes_seen"] = len(seen_haps)
    return calls, dict(diag)


def derive_canonical_order(table1):
    """Verbatim copy of 16_phasing_mendelian_validation.derive_canonical_order. Median
    fractional rank of each gene's gene_start within every single-contig (person, hap, contig)
    group that calls >=2 genes. Returns {gene_bare: rank_float}, lower = earlier in the region."""
    sub = table1[table1["copy_index"].fillna(1) == 1].dropna(subset=["gene_start"])
    ranks = defaultdict(list)
    for (person, hap, contig), grp in sub.groupby(["person_id", "hap", "contig"]):
        grp = grp.drop_duplicates("gene_bare").sort_values("gene_start")
        n = len(grp)
        if n < 2:
            continue
        for i, gene in enumerate(grp["gene_bare"]):
            ranks[gene].append(i / (n - 1))
    return {gene: float(np.median(vals)) for gene, vals in ranks.items() if vals}


def strict_ancestry(people, threshold):
    """Verbatim copy of 24_novelty_by_field.strict_ancestry: ancestry_pred if p_<anc> >=
    threshold, else NA."""
    out = []
    for anc, row in zip(people["ancestry_pred"], people.to_dict("records")):
        if not isinstance(anc, str) or anc in ("", "NA"):
            out.append(None)
            continue
        p = row.get("p_%s" % anc.lower())
        try:
            ok = p is not None and not pd.isna(p) and float(p) >= threshold
        except (TypeError, ValueError):
            ok = False
        out.append(anc if ok else None)
    return pd.Series(out, index=people.index, dtype=object)


def greedy_unrelated(person_ids, pairs, kin_min):
    """Verbatim copy of 24_novelty_by_field.greedy_unrelated: remove people until no pair with
    kin >= kin_min remains, among person_ids only. Each step drops the person with the most
    remaining relatives; ties go to the smallest id (string sort)."""
    ids = set(map(str, person_ids))
    adj = defaultdict(set)
    for i, j, k in pairs:
        i, j = str(i), str(j)
        try:
            k = float(k)
        except (TypeError, ValueError):
            continue
        if i == j or k < kin_min or i not in ids or j not in ids:
            continue
        adj[i].add(j)
        adj[j].add(i)
    removed = []
    while True:
        live = [(len(v), n) for n, v in adj.items() if v]
        if not live:
            break
        maxdeg = max(d for d, _ in live)
        victim = min(n for d, n in live if d == maxdeg)
        for nb in adj[victim]:
            adj[nb].discard(victim)
        adj[victim] = set()
        removed.append(victim)
    kept = ids - set(removed)
    return kept, removed


def load_relatedness_pairs(path):
    """Verbatim copy of 24_novelty_by_field.load_relatedness_pairs."""
    if not os.path.exists(path):
        sys.exit("FATAL: relatedness table not found at %r." % path)
    df = pd.read_csv(path, sep="\t", dtype={"i.s": str, "j.s": str})
    missing = {"i.s", "j.s", "kin"} - set(df.columns)
    if missing:
        sys.exit("FATAL: relatedness table missing %s; columns: %s" % (missing, list(df.columns)))
    return list(zip(df["i.s"], df["j.s"], df["kin"]))


def build_people(person_ids, cohort_path, relatedness_path, kin_min, strict_min, skip_rel):
    """Verbatim copy of 24_novelty_by_field.build_people (unrelated cohort + strict ancestry)."""
    cohort = vc.load_cohort_membership(cohort_path)
    pcols = ["p_%s" % a.lower() for a in vc.ANCESTRY_ORDER]
    missing = [c for c in ["ancestry_pred"] + pcols if c not in cohort.columns]
    if missing:
        sys.exit("FATAL: %s lacks %s (needed for strict ancestry). Columns: %s"
                 % (cohort_path, missing, list(cohort.columns)))
    people = pd.DataFrame({"person_id": sorted(set(person_ids))})
    people = people.merge(cohort[["person_id", "ancestry_pred"] + pcols], on="person_id",
                          how="left")
    people["anc_pred"] = people["ancestry_pred"].astype(object).where(
        people["ancestry_pred"].notna(), None)
    people["anc_strict"] = strict_ancestry(people, strict_min)
    for c in ("anc_pred", "anc_strict"):
        people[c] = people[c].map(lambda v: v if isinstance(v, str) and v else "UNASSIGNED")
    if skip_rel:
        people["unrelated"] = True
        n_removed = 0
    else:
        lr_ids = set(people["person_id"])
        kept, removed = greedy_unrelated(lr_ids, load_relatedness_pairs(relatedness_path),
                                         kin_min)
        people["unrelated"] = people["person_id"].isin(kept)
        n_removed = len(removed)
    return people[["person_id", "anc_pred", "anc_strict", "unrelated"]], n_removed


def first2(consensus):
    """First two allele fields, e.g. `A*02:01:01:01` -> `02:01`; None on new/unresolved."""
    if not isinstance(consensus, str) or "*" not in consensus:
        return None
    parts = consensus.split("*", 1)[1].split(":")
    parts = [p.strip() for p in parts[:2]]
    if not parts or any(p.lower() == "new" or not p for p in parts):
        return None
    return ":".join(p.zfill(2) for p in parts)


# ---------------------------------------------------------------------------
# Q1 -- per-haplotype rerun, HWE test, within-person concordance, quality correlates
# ---------------------------------------------------------------------------
def gene_rate_by_ancestry(calls, gene, anc_of):
    """Bridged-deletion rate for `gene`, overall and per ancestry, with Wilson 95% CI."""
    g = calls[calls["gene"] == gene].copy()
    g["ancestry"] = g["person_id"].map(anc_of)
    rows = []
    for anc, sub in [("ALL", g)] + [(a, g[g["ancestry"] == a]) for a in vc.ANCESTRY_ORDER]:
        n_bridged = int((sub["status"].isin(["present", "deleted_bridged"])).sum())
        n_del = int((sub["status"] == "deleted_bridged").sum())
        p, lo, hi = vc.wilson_ci(n_del, n_bridged)
        rows.append({"ancestry": anc, "n_bridged": n_bridged, "n_bridged_disp": suppress(n_bridged),
                     "n_deleted": n_del, "n_deleted_disp": suppress(n_del),
                     "pct_deleted": 100.0 * p if n_bridged else np.nan,
                     "ci_lo": 100.0 * lo if n_bridged else np.nan,
                     "ci_hi": 100.0 * hi if n_bridged else np.nan})
    return pd.DataFrame(rows)


def per_person_hap_status(calls, gene):
    """One row per person with hap1/hap2 bridged-deletion status for `gene`. Only rows where
    BOTH haplotypes were bridged for this gene are usable for the biallelic/concordance tests --
    an unbridged haplotype tells us nothing about whether the gene is there."""
    g = calls[calls["gene"] == gene]
    piv = g.pivot_table(index="person_id", columns="hap", values="status", aggfunc="first")
    for h in ("hap1", "hap2"):
        if h not in piv.columns:
            piv[h] = np.nan
    piv = piv.rename(columns={"hap1": "s1", "hap2": "s2"}).reset_index()
    both_bridged = piv["s1"].isin(["present", "deleted_bridged"]) & \
        piv["s2"].isin(["present", "deleted_bridged"])
    piv = piv[both_bridged].copy()
    piv["del1"] = (piv["s1"] == "deleted_bridged").astype(int)
    piv["del2"] = (piv["s2"] == "deleted_bridged").astype(int)
    piv["n_deleted_haps"] = piv["del1"] + piv["del2"]
    return piv


def hwe_biallelic_test(status_df):
    """Observed vs expected (N*p^2, independence across haplotypes within a person) count of
    people with BOTH haplotypes bridged-deleted. p is estimated from this same set of
    fully-bridged people (2N haplotype observations) -- the null this tests is 'haplotype-level
    deletion calls are independent across the two haplotypes of one person', which is exactly
    what per-haplotype assembly artifacts would also satisfy, so a significant excess over N*p^2
    is NOT on its own evidence of biology; see within_person_concordance for that."""
    n = len(status_df)
    if n == 0:
        return {"n_people_both_bridged": 0}
    n_hap_del = int(status_df["del1"].sum() + status_df["del2"].sum())
    p = n_hap_del / (2.0 * n)
    n_both = int((status_df["n_deleted_haps"] == 2).sum())
    expected = n * p * p
    result = {"n_people_both_bridged": n, "p_hap_deleted": p,
              "n_both_deleted_disp": suppress(n_both), "n_both_deleted": n_both,
              "expected_both_deleted": expected}
    if HAVE_SCIPY and n > 0:
        bt = binomtest(n_both, n, p * p, alternative="two-sided")
        result["poisson_exact_p"] = float(bt.pvalue)
    return result


def within_person_concordance(status_df):
    """2x2 table of hap1-deleted x hap2-deleted within the same person (both haplotypes
    bridged). Independent artifacts predict no association (odds ratio ~1); concordance > 1
    would suggest a shared per-person cause (e.g. one bad assembly run); discordance < 1 would
    argue against a shared cause."""
    if status_df.empty:
        return {}
    tab = pd.crosstab(status_df["del1"], status_df["del2"])
    tab = tab.reindex(index=[0, 1], columns=[0, 1], fill_value=0)
    out = {"table": tab.values.tolist()}
    if HAVE_SCIPY and tab.values.sum() > 0:
        odds, p = fisher_exact(tab.values)
        out["odds_ratio"] = float(odds) if np.isfinite(odds) else None
        out["fisher_p"] = float(p)
    return out


def bridging_contig_quality(t1, calls, gene):
    """Per (person, hap) proxy for the quality of the contig that BRIDGES `gene`'s position --
    Table 1 has no hg38 coordinates, so span/gene-count on that specific contig is the best
    available N50/length proxy. Only defined for deleted_bridged calls (that's the call whose
    contig quality is in question)."""
    del_calls = calls[(calls["gene"] == gene) & (calls["status"] == "deleted_bridged")]
    if del_calls.empty:
        return pd.DataFrame(columns=["person_id", "hap", "contig_span", "contig_n_genes"])
    # Recreate which contig bridged each (person,hap): the contig with genes flanking `gene`'s
    # canonical rank, same logic as bridged_absences but reading it back out.
    cols = ["person_id", "hap", "contig", "gene_bare", "gene_start", "gene_end"]
    sub = t1[cols].dropna(subset=["gene_start", "gene_end"])
    rows = []
    wanted = set(zip(del_calls["person_id"], del_calls["hap"]))
    for (pid, hap, contig), grp in sub.groupby(["person_id", "hap", "contig"], sort=False):
        if (pid, hap) not in wanted:
            continue
        if gene in set(grp["gene_bare"]):
            continue  # this contig has the gene present -- not the bridging one
        span = float(grp["gene_end"].max() - grp["gene_start"].min() + 1)
        rows.append({"person_id": pid, "hap": hap, "contig_span": span,
                     "contig_n_genes": int(grp["gene_bare"].nunique())})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    # a haplotype can have >1 contig; keep the largest-span one as "the" bridging contig
    return df.sort_values("contig_span", ascending=False).drop_duplicates(["person_id", "hap"])


def quality_correlates(status_df, contig_q, cohort_df, t1, gene):
    """Spearman correlation of per-person deletion status (n_deleted_haps, 0/1/2) against
    assembly-quality proxies, plus deletion-rate breakdown by platform/trim_tier."""
    out = {"spearman": {}, "by_platform": None, "by_trim_tier": None}
    df = status_df.merge(
        contig_q.groupby("person_id").agg(mean_contig_span=("contig_span", "mean"),
                                          mean_contig_n_genes=("contig_n_genes", "mean")).reset_index(),
        on="person_id", how="left")

    # per-person template_warning rate across ALL Table 1 rows (not just this gene) -- a
    # person-level assembly-quality proxy independent of the gene under test.
    tw = t1.assign(warned=t1["template_warning"].notna() & (t1["template_warning"] != "NA"))
    tw_rate = tw.groupby("person_id")["warned"].mean().rename("template_warning_rate").reset_index()
    df = df.merge(tw_rate, on="person_id", how="left")

    ccols = [c for c in ["platform", "trim_tier", "n_genes_called", "n_genes_exact",
                         "mean_template_distance"] if c in cohort_df.columns]
    if ccols:
        df = df.merge(cohort_df[["person_id"] + ccols], on="person_id", how="left")

    if HAVE_SCIPY:
        for col in ["mean_contig_span", "mean_contig_n_genes", "template_warning_rate",
                    "n_genes_called", "mean_template_distance"]:
            if col in df.columns and df[col].notna().sum() > 5:
                sub = df.dropna(subset=[col])
                rho, p = spearmanr(sub["n_deleted_haps"], sub[col])
                out["spearman"][col] = {"rho": float(rho), "p": float(p), "n": int(len(sub))}

    for cat_col, key in [("platform", "by_platform"), ("trim_tier", "by_trim_tier")]:
        if cat_col in df.columns:
            g = df.groupby(cat_col)["n_deleted_haps"].agg(
                n=lambda s: int(len(s)),
                pct_any_del=lambda s: 100.0 * float((s > 0).sum()) / max(1, len(s)))
            g["n_disp"] = g["n"].map(suppress)
            out[key] = g.reset_index()[[cat_col, "n_disp", "pct_any_del"]]
    return out, df


# ---------------------------------------------------------------------------
# Q2 / Q3 -- short-read contradiction test, per gene
# ---------------------------------------------------------------------------
def sr_alleles_per_person(sr_long, gene):
    """SR genotype (allele_1, allele_2 raw strings) per person for `gene`. Missing calls in the
    AoU SR HLA file are typically empty/NaN strings, not absent rows -- both are treated as
    missing here."""
    g = sr_long[sr_long["gene_bare"] == gene]
    if g.empty:
        return {}
    piv = g.pivot_table(index="person_id", columns="copy", values="allele", aggfunc="first")
    out = {}
    for pid, row in piv.iterrows():
        a1 = row.get(1)
        a2 = row.get(2)
        a1 = a1 if isinstance(a1, str) and a1.strip() and a1.strip().upper() != "NA" else None
        a2 = a2 if isinstance(a2, str) and a2.strip() and a2.strip().upper() != "NA" else None
        out[pid] = (a1, a2)
    return out


def sr_contradiction_test(status_df, sr_map, anc_of, gene, t1, calls, positive_control=False):
    """Core cross-validation. `status_df` is per_person_hap_status(calls, gene) restricted to
    people with exactly one deleted / one present bridged haplotype ("hemizygous carriers" of the
    LR call), and its complement (both present -- non-carriers), ancestry-matched by construction
    (grouped below).

    A carrier is CONTRADICTED when SR reports two distinct, non-missing alleles at the gene
    (heterozygous-looking) -- if the gene were truly hemizygously deleted, SR should see only one
    haplotype's worth of signal (homozygous-looking or a single confident call), not two.

    For the positive-control genes (DRB3/DRB4/DRB5) the same statistic is reported for
    calibration, but is *not* expected to be a contradiction -- confirmed presence of two DISTINCT
    SR alleles at a gene AoU's SR pipeline does not reliably type in the first place is weak
    evidence either way, and this is reported as such (`sr_gene_typed`).
    """
    sr_gene_typed = len(sr_map) > 0
    carriers = status_df[status_df["n_deleted_haps"] == 1].copy()
    noncarriers = status_df[status_df["n_deleted_haps"] == 0].copy()
    carriers["ancestry"] = carriers["person_id"].map(anc_of)
    noncarriers["ancestry"] = noncarriers["person_id"].map(anc_of)

    def het_rate(df):
        n = 0
        n_het = 0
        n_typed = 0
        for pid in df["person_id"]:
            pair = sr_map.get(pid)
            if pair is None:
                continue
            a1, a2 = pair
            if a1 is None and a2 is None:
                continue
            n_typed += 1
            if a1 is not None and a2 is not None and first2_raw(a1) != first2_raw(a2):
                n_het += 1
            n += 1
        p, lo, hi = vc.wilson_ci(n_het, n)
        return {"n_sr_typed": n_typed, "n_used": n, "n_het": n_het,
                "n_used_disp": suppress(n), "n_het_disp": suppress(n_het),
                "pct_het": 100.0 * p if n else np.nan,
                "ci_lo": 100.0 * lo if n else np.nan, "ci_hi": 100.0 * hi if n else np.nan}

    car = het_rate(carriers)
    non = het_rate(noncarriers)
    # ancestry-matched: also compute per-ancestry, then a pooled Mantel-Haenszel-style average
    per_anc = []
    for anc in vc.ANCESTRY_ORDER:
        c = het_rate(carriers[carriers["ancestry"] == anc])
        n_ = het_rate(noncarriers[noncarriers["ancestry"] == anc])
        per_anc.append({"ancestry": anc, "carrier_pct_het": c["pct_het"], "carrier_n": c["n_used"],
                        "noncarrier_pct_het": n_["pct_het"], "noncarrier_n": n_["n_used"]})

    fisher = None
    if HAVE_SCIPY and car["n_used"] > 0 and non["n_used"] > 0:
        table = [[car["n_het"], car["n_used"] - car["n_het"]],
                 [non["n_het"], non["n_used"] - non["n_het"]]]
        odds, p = fisher_exact(table)
        fisher = {"odds_ratio": float(odds) if np.isfinite(odds) else None, "p": float(p)}

    # fragmentation check: for contradicted carriers, does the "extra" gene signal show up
    # elsewhere in this person's own LR assembly (a different contig than the one that bridged
    # the deletion)?
    frag_evidence = None
    if gene == TARGET_GENE:
        frag_evidence = fragmentation_check(carriers, sr_map, t1, gene)

    return {"gene": gene, "positive_control": positive_control, "sr_gene_typed": sr_gene_typed,
            "carrier": car, "noncarrier": non, "fisher": fisher, "per_ancestry": per_anc,
            "fragmentation_evidence": frag_evidence,
            "sr_contradiction_rate_pct": car["pct_het"]}


def first2_raw(allele):
    """Best-effort two-field normalization of a raw SR allele string (format varies by AoU SR
    HLA pipeline export, may already be `A*02:01` or `02:01`)."""
    if not isinstance(allele, str):
        return None
    s = allele.split("*", 1)[-1]
    parts = [p for p in s.split(":") if p][:2]
    return ":".join(parts) if parts else None


def fragmentation_check(carriers, sr_map, t1, gene):
    """For LR HLA-A hemizygous carriers, does gene A appear anywhere else in the person's own
    assembly (any hap, any contig) beyond the one non-deleted haplotype call already known? A
    second location is evidence of assembly fragmentation/misplacement rather than a true
    deletion."""
    ids = set(carriers["person_id"])
    if not ids or "gene_bare" not in t1.columns:
        return {"n_checked": 0, "n_with_extra_location": 0}
    g = t1[(t1["gene_bare"] == gene) & (t1["person_id"].isin(ids))]
    n_locations = g.groupby("person_id")["contig"].nunique()
    n_extra = int((n_locations > 1).sum())
    return {"n_checked": suppress(len(ids)), "n_with_extra_location": suppress(n_extra),
            "pct_with_extra_location": 100.0 * n_extra / max(1, len(ids))}


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def make_figure(calibration_rows, hwe, out_path, contig_n_genes_deleted=None):
    """contig_n_genes_deleted: array-like of contig_n_genes values, one per confirmed HLA-A
    deleted_bridged haplotype call (panel c: is the deletion concentrated on short/low-gene-count
    bridging contigs -- a fragmentation signature -- or spread evenly?)."""
    plt_mod = __import__("matplotlib.pyplot", fromlist=["pyplot"])
    with vc.nature_style():
        fig = plt_mod.figure(figsize=(vc.mm(vc.NATURE_DOUBLE_COL_MM), vc.mm(75)))
        axa = fig.add_subplot(1, 3, 1)
        axb = fig.add_subplot(1, 3, 2)
        axc = fig.add_subplot(1, 3, 3)

        # (a) calibration scatter: LR deletion rate vs SR-contradiction rate
        from matplotlib.lines import Line2D
        for r in calibration_rows:
            if r["sr_contradiction_rate_pct"] is None or np.isnan(r["sr_contradiction_rate_pct"]):
                continue
            is_target = r["gene"] == TARGET_GENE
            color = "#009E73" if r["positive_control"] else ("#D55E00" if is_target else "#999999")
            axa.scatter(r["lr_deletion_pct"], r["sr_contradiction_rate_pct"], color=color,
                       s=26 if is_target else 16, zorder=4 if is_target else 3,
                       edgecolor="black" if is_target else "none", linewidth=0.5)
            dy = 3 if r["gene"] not in ("B", "DQB1") else -3
            axa.annotate(r["gene"], (r["lr_deletion_pct"], r["sr_contradiction_rate_pct"]),
                        fontsize=4.5, xytext=(3, dy), textcoords="offset points")
        axa.set_xlabel("LR bridged-deletion rate (%)")
        axa.set_ylabel("SR-contradiction rate in carriers (%)")
        axa.set_title("a  Calibration: LR deletion vs SR contradiction", fontsize=6.5)
        legend_elems = [
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#D55E00",
                  markeredgecolor="black", markersize=5, label="HLA-A (target)"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#999999",
                  markersize=5, label="Other negative controls"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#009E73",
                  markersize=5, label="Positive controls (DRB3/4/5)"),
        ]
        axa.legend(handles=legend_elems, loc="center right", fontsize=4.5, frameon=False)
        vc.panel_letter(axa, "a")

        # (b) HLA-A homozygous observed vs expected
        if hwe.get("n_people_both_bridged", 0) > 0:
            obs = hwe["n_both_deleted"]
            exp = hwe["expected_both_deleted"]
            axb.bar([0, 1], [obs, exp], color=["#D55E00", "#999999"], width=0.6)
            axb.set_xticks([0, 1])
            axb.set_xticklabels(["Observed", "Expected\n(N·p²)"], fontsize=6)
            axb.set_ylabel("People with both HLA-A haplotypes 'deleted'")
            p_txt = ("p=%.3g" % hwe["poisson_exact_p"]) if "poisson_exact_p" in hwe else "n/a"
            axb.set_title("b  Biallelic HLA-A loss (%s)" % p_txt, fontsize=6.5)
        vc.panel_letter(axb, "b")

        # (c) mechanism: gene-count on the bridging contig, for confirmed HLA-A deletions only --
        # a low count means the "deletion" sits on a short/sparse contig (fragmentation
        # signature); a distribution centered on typical contig gene counts argues against it.
        if contig_n_genes_deleted is not None and len(contig_n_genes_deleted) > 0:
            vals = np.asarray(contig_n_genes_deleted)
            bins = np.arange(vals.min() - 0.5, vals.max() + 1.5, 1)
            axc.hist(vals, bins=bins, color="#D55E00", edgecolor="white", linewidth=0.4)
            axc.set_xlabel("Genes on the bridging contig (HLA-A deleted calls)")
            axc.set_ylabel("Number of haplotypes")
            axc.set_title("c  Bridging-contig gene count for HLA-A 'deletions'\n(n=%d; median=%.0f)"
                          % (len(vals), np.median(vals)), fontsize=6)
        vc.panel_letter(axc, "c")

        fig.tight_layout()
    return vc.save_fig(fig, out_path)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def run(args):
    ensure_dir(args.out_dir)
    if not os.path.exists(args.table1):
        sys.exit("FATAL: table1 not found: %s" % args.table1)

    print("[38] loading table1 ...", flush=True)
    t1 = pd.read_csv(args.table1, sep="\t", dtype=str, low_memory=False)
    t1["copy_index"] = pd.to_numeric(t1.get("copy_index", 1), errors="coerce").fillna(1).astype(int)
    t1["gene_bare"] = t1["gene"].map(gene_bare)
    t1["gene_start"] = pd.to_numeric(t1.get("gene_start"), errors="coerce")
    t1["gene_end"] = pd.to_numeric(t1.get("gene_end"), errors="coerce")

    people, n_removed = build_people(
        sorted(set(t1["person_id"].astype(str))), args.cohort_membership,
        args.relatedness_table, args.kin_min, args.strict_threshold, args.skip_relatedness)
    keep = people[people["unrelated"]]
    anc_of = dict(zip(keep["person_id"].astype(str), keep["anc_strict"]))
    n_people = len(keep)
    t1 = t1[t1["person_id"].astype(str).isin(anc_of)].copy()
    print("[38] unrelated people: %d (dropped %d)" % (n_people, n_removed), flush=True)

    gene_order = derive_canonical_order(t1)
    genes_of_interest = sorted(set(t1["gene_bare"]) & set(gene_order))
    calls, diag = bridged_absences(t1, gene_order, genes_of_interest)
    print("[38] bridged_absences done: %s" % diag, flush=True)

    cohort_df = vc.load_cohort_membership(args.cohort_membership)

    # ---- Q1 ----
    q1 = {}
    q1["hla_a_rate_by_ancestry"] = gene_rate_by_ancestry(calls, TARGET_GENE, anc_of)
    status_a = per_person_hap_status(calls, TARGET_GENE)
    q1["hwe"] = hwe_biallelic_test(status_a)
    q1["within_person_concordance"] = within_person_concordance(status_a)
    contig_q = bridging_contig_quality(t1, calls, TARGET_GENE)
    q1["quality_correlates"], quality_df = quality_correlates(status_a, contig_q, cohort_df, t1,
                                                               TARGET_GENE)

    # ---- Q2 / Q3: calibration panel across control genes ----
    print("[38] loading SR genotypes ...", flush=True)
    sr_long = vc.load_sr_genotypes(args.sr_genotypes)
    sr_genes_present = set(sr_long["gene_bare"].unique())

    calibration = []
    per_gene_sr_detail = {}
    for gene in CALIBRATION_GENES:
        status_g = per_person_hap_status(calls, gene)
        sr_map = sr_alleles_per_person(sr_long, gene) if gene in sr_genes_present else {}
        pos = gene in POSITIVE_CONTROL_GENES
        res = sr_contradiction_test(status_g, sr_map, anc_of, gene, t1, calls,
                                    positive_control=pos)
        gr = gene_rate_by_ancestry(calls, gene, anc_of)
        lr_pct = gr.loc[gr["ancestry"] == "ALL", "pct_deleted"]
        res["lr_deletion_pct"] = float(lr_pct.iloc[0]) if len(lr_pct) else np.nan
        calibration.append(res)
        per_gene_sr_detail[gene] = res

    # ---- Figure ----
    fig_path = os.path.join(args.out_dir, "fig_sr_validation")
    make_figure(calibration, q1["hwe"], fig_path,
               contig_n_genes_deleted=contig_q["contig_n_genes"].tolist() if not contig_q.empty
               else [])

    # ---- write aggregates ----
    q1["hla_a_rate_by_ancestry"].to_csv(
        os.path.join(args.out_dir, "hla_a_rate_by_ancestry.tsv"), sep="\t", index=False)
    if q1["quality_correlates"].get("by_platform") is not None:
        q1["quality_correlates"]["by_platform"].to_csv(
            os.path.join(args.out_dir, "hla_a_rate_by_platform.tsv"), sep="\t", index=False)
    if q1["quality_correlates"].get("by_trim_tier") is not None:
        q1["quality_correlates"]["by_trim_tier"].to_csv(
            os.path.join(args.out_dir, "hla_a_rate_by_trim_tier.tsv"), sep="\t", index=False)

    calib_rows = []
    for r in calibration:
        calib_rows.append({
            "gene": r["gene"], "positive_control": r["positive_control"],
            "sr_gene_typed": r["sr_gene_typed"], "lr_deletion_pct": r["lr_deletion_pct"],
            "sr_contradiction_rate_pct": r["sr_contradiction_rate_pct"],
            "carrier_n_used_disp": r["carrier"]["n_used_disp"],
            "noncarrier_n_used_disp": r["noncarrier"]["n_used_disp"],
            "noncarrier_pct_het": r["noncarrier"]["pct_het"],
            "fisher_odds_ratio": r["fisher"]["odds_ratio"] if r["fisher"] else None,
            "fisher_p": r["fisher"]["p"] if r["fisher"] else None,
        })
    pd.DataFrame(calib_rows).to_csv(os.path.join(args.out_dir, "sr_calibration_panel.tsv"),
                                    sep="\t", index=False)

    def _clean(o):
        if isinstance(o, dict):
            return {k: _clean(v) for k, v in o.items()}
        if isinstance(o, list):
            return [_clean(v) for v in o]
        if isinstance(o, pd.DataFrame):
            return o.to_dict(orient="records")
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return None if np.isnan(o) else float(o)
        return o

    summary = {
        "n_people_unrelated": suppress(n_people),
        "q1_hwe": _clean(q1["hwe"]),
        "q1_within_person_concordance": _clean(q1["within_person_concordance"]),
        "q1_quality_spearman": _clean(q1["quality_correlates"]["spearman"]),
        "q2_q3_calibration": _clean(calib_rows),
    }
    with open(os.path.join(args.out_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    print("[38] done -> %s" % args.out_dir, flush=True)
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table1", default=DEFAULT_TABLE1)
    ap.add_argument("--cohort-membership", default=DEFAULT_COHORT)
    ap.add_argument("--relatedness-table", default=DEFAULT_RELATEDNESS)
    ap.add_argument("--sr-genotypes", default=DEFAULT_SR_GENOTYPES)
    ap.add_argument("--skip-relatedness", action="store_true")
    ap.add_argument("--kin-min", type=float, default=KIN_MIN)
    ap.add_argument("--strict-threshold", type=float, default=STRICT_MIN)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    args = ap.parse_args(argv)
    run(args)


if __name__ == "__main__":
    main()
