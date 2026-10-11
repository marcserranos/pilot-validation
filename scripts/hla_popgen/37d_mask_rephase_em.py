#!/usr/bin/env python3
"""WS1 deliverable 3 (S03 call #8, DQ G1/G2): the mask-and-rephase experiment.

Cole's own DQ G1/G2 figure (built from *statistically* phased data) shows light-blue cells in
the "predicted incompatible" cross-group quadrants. The 2-field committed pass
(`37c_dq_g1g2_from_committed.py`) and the 4-field VM pass (`37_vm_run.py` / `37_dq_g1g2_signed_ld.py`)
both show that *physically* phased (same-contig) DQA1~DQB1 haplotypes are completely purged of
cross-group pairs. This script turns "statistical phasing errors vs real rare recombinants" into a
measured number:

1. **Truth set** -- people with DQA1 and DQB1 physically phased (same contig) on BOTH hap1 and
   hap2 (`extract_haplotypes` from `37_vm_run.py`, reused verbatim, resolution=2 by default).
   Their true cis pairs are known exactly.
2. **Mask + re-phase.** Take only the truth set's unordered two-locus genotype (drop which allele
   pairs with which). Fit population haplotype frequencies by EM (Excoffier & Slatkin 1995) on the
   unphased genotypes of every unrelated person of the SAME ancestry who carries two calls at each
   locus (truth set + everyone else -- the realistic "what would a statistical phaser see" input).
   Re-phase each truth-set person with the most-likely diplotype under the converged frequencies.
   Count switch errors (most-likely resolution != true phase, doubly-heterozygous people only) and
   count how many EM-inferred cis haplotypes are G1/G2 cross-group ("predicted_incompatible") vs 0
   in the true phase.
3. **Real-world inferred-phase group.** People whose DQA1 and DQB1 calls sit on DIFFERENT contigs
   on a given hap (no physical cis evidence) -- what the pipeline's own hap1/hap2 assembly label
   would imply if naively used as the cis pairing. Count implied cross-group pairs, per ancestry.

Per ancestry (AFR expected worst, per the task) plus pooled. All counts are the raw *cell* counts
here (haplotype pairs / people), not participant-level rows -- aggregated before any TSV is
written, and any count 1-19 is suppressed to "<20" per the disclosure rule.

Usage (VM):
    setsid nohup python3 -u 37d_mask_rephase_em.py --out-dir ~/s03/results/37 \
        --resolution 2 < /dev/null > ~/s03/results/37/37d_run.log 2>&1 & disown
"""
import argparse
import importlib.util
import json
import os
import random
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import _viz_common as vc  # noqa: E402

SUPPRESS_BELOW = 20
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
EM_MAX_ITER = 200
EM_TOL = 1e-10


def _load_module(filename, modname):
    path = os.path.join(_THIS_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


def _sup(n, t=SUPPRESS_BELOW):
    if n == 0:
        return "0"
    if n < t:
        return "<%d" % t
    return str(n)


# ---------------------------------------------------------------------------
# EM haplotype-frequency estimator (Excoffier & Slatkin 1995), generic 2-locus,
# multiallelic, most-likely-pair (MLP) diplotype assignment.
# ---------------------------------------------------------------------------
def em_haplotype_freqs(genotypes, max_iter=EM_MAX_ITER, tol=EM_TOL, seed=0):
    """genotypes: list of ((a1, a2), (b1, b2)) unordered allele pairs at locus A / locus B
    (homozygotes as (a, a)). Returns dict {(a, b): freq} over all haplotypes with any support,
    converged by log-likelihood or max-iter, whichever first."""
    rng = random.Random(seed)
    a_alleles = sorted({a for (a1, a2), _ in genotypes for a in (a1, a2)})
    b_alleles = sorted({b for _, (b1, b2) in genotypes for b in (b1, b2)})
    haps = [(a, b) for a in a_alleles for b in b_alleles]
    idx = {h: i for i, h in enumerate(haps)}
    freq = np.full(len(haps), 1.0 / len(haps))

    # Pre-enumerate, per genotype, the compatible haplotype-pair resolutions.
    # Unambiguous (<=1 het locus): a single resolution. Both loci het: two resolutions
    # (cis/coupling vs trans/repulsion), weighted by current freq each E-step.
    resolved = []
    for (a1, a2), (b1, b2) in genotypes:
        a_het = a1 != a2
        b_het = b1 != b2
        if not a_het or not b_het:
            resolved.append([((a1, b1), (a2, b2))])
        else:
            resolved.append([((a1, b1), (a2, b2)), ((a1, b2), (a2, b1))])

    ll_prev = -np.inf
    for _ in range(max_iter):
        counts = np.zeros(len(haps))
        ll = 0.0
        for res in resolved:
            weights = []
            for h1, h2 in res:
                w = freq[idx[h1]] * freq[idx[h2]]
                if h1 == h2:
                    pass  # already the right diplotype probability (p^2), no factor of 2
                weights.append(w)
            tot = sum(weights)
            if tot <= 0:
                # Should not happen once every haplotype has nonzero support; guard anyway.
                weights = [1.0] * len(res)
                tot = float(len(res))
            ll += np.log(tot) if tot > 0 else 0.0
            for (h1, h2), w in zip(res, weights):
                post = w / tot
                counts[idx[h1]] += post
                counts[idx[h2]] += post
        new_freq = counts / counts.sum()
        if abs(ll - ll_prev) < tol:
            freq = new_freq
            break
        freq = new_freq
        ll_prev = ll
    return {h: float(freq[idx[h]]) for h in haps}


def most_likely_pair(a1, a2, b1, b2, freq):
    """Return the higher-posterior diplotype [(h1a,h1b), (h2a,h2b)] under converged `freq`,
    tie-broken deterministically (first enumeration) since ties are rare (freq is continuous)."""
    a_het = a1 != a2
    b_het = b1 != b2
    if not a_het or not b_het:
        return [(a1, b1), (a2, b2)]
    opt1 = [(a1, b1), (a2, b2)]
    opt2 = [(a1, b2), (a2, b1)]
    w1 = freq.get(opt1[0], 0.0) * freq.get(opt1[1], 0.0)
    w2 = freq.get(opt2[0], 0.0) * freq.get(opt2[1], 0.0)
    return opt1 if w1 >= w2 else opt2


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def run(args):
    m37 = _load_module("37_vm_run.py", "s03_37_vm_run")
    vc.ensure_dir(args.out_dir)

    if not os.path.exists(m37.TABLE1):
        sys.exit("FATAL: %s not found -- needs a live VM session." % m37.TABLE1)

    print("[37d] loading table1 ...", flush=True)
    t1 = pd.read_csv(m37.TABLE1, sep="\t", dtype=str, low_memory=False)
    unrelated_ids, _anc_of_buggy = m37.build_people(t1, m37.COHORT, m37.RELATEDNESS,
                                                      m37.KIN_MIN, m37.STRICT_MIN)
    # NOTE: 37_vm_run.py's own build_people() picks the ancestry column by substring match
    # ("anc" in colname), and "max_template_distance"/"mean_template_distance" both contain
    # "anc" (distANCe) and sort before "ancestry_pred" in cohort_membership.tsv's column
    # order -- so its returned anc_of silently maps everyone to a template-distance bucket,
    # not AFR/AMR/EAS/EUR/MID/SAS (verified: values were '1.0','2.0','5.0', nan, ... not
    # ancestry strings). unrelated_ids (relatedness-derived) is unaffected and reused as-is;
    # ancestry is instead re-derived directly here from cohort_membership.tsv's actual
    # 'ancestry_pred' column.
    cohort_df = pd.read_csv(m37.COHORT, sep="\t", dtype=str)
    anc_of = dict(zip(cohort_df["person_id"].astype(str), cohort_df["ancestry_pred"]))
    t1 = t1[t1["person_id"].astype(str).isin(unrelated_ids)].copy()
    print("[37d] unrelated people: %s" % _sup(len(unrelated_ids)), flush=True)

    res = args.resolution
    gene_a, gene_b = "DQA1", "DQB1"

    # --- physically-phased cis pairs, per (person, hap) -----------------------------------
    haps_df = m37.extract_haplotypes(t1, gene_a, gene_b, res)  # person_id, hap, allele_a, allele_b
    haps_df["ancestry"] = haps_df["person_id"].astype(str).map(anc_of)

    per_person_phys = defaultdict(dict)  # pid -> {hap: (a, b)}
    for _, r in haps_df.iterrows():
        per_person_phys[r["person_id"]][r["hap"]] = (r["allele_a"], r["allele_b"])

    truth_people = {pid: h for pid, h in per_person_phys.items()
                    if set(h.keys()) >= {"hap1", "hap2"}}
    print("[37d] truth-set people (both haps physically phased, res=%d): %s"
          % (res, _sup(len(truth_people))), flush=True)

    # --- assembly-labelled per-hap calls regardless of contig (for the "different contig"
    # real-world group) --------------------------------------------------------------------
    bare = t1["gene"].astype(str).str.replace("^HLA-", "", regex=True)
    sub = t1[bare.isin([gene_a, gene_b])].copy()
    sub["gene_bare"] = bare[sub.index]
    per_hap_calls = defaultdict(dict)  # (pid, hap) -> {gene: (allele, contig)}
    for (pid, hap), grp in sub.groupby(["person_id", "hap"], sort=False):
        for gene in (gene_a, gene_b):
            grows = grp[grp["gene_bare"] == gene]
            if len(grows) != 1:
                continue
            allele = m37.to_nfield(grows.iloc[0]["consensus"], res)
            if allele is None:
                continue
            per_hap_calls[(pid, hap)][gene] = (allele, grows.iloc[0]["contig"])

    diffcontig_rows = []  # pid, ancestry, hap, allele_a, allele_b, incompatible
    for (pid, hap), d in per_hap_calls.items():
        if gene_a not in d or gene_b not in d:
            continue
        (aa, ca), (bb, cb) = d[gene_a], d[gene_b]
        if ca == cb:
            continue  # physically phased -- already covered above
        cls = m37.classify_pair(aa, bb)
        diffcontig_rows.append({
            "person_id": pid, "ancestry": anc_of.get(pid), "hap": hap,
            "allele_a": aa, "allele_b": bb,
            "incompatible": cls == "predicted_incompatible",
        })
    diffcontig_df = pd.DataFrame(diffcontig_rows)

    # --- per-person unordered genotype (for EM fitting population), all unrelated people ---
    genos_by_person = {}
    geno_a = defaultdict(list)
    geno_b = defaultdict(list)
    for (pid, hap), d in per_hap_calls.items():
        if gene_a in d:
            geno_a[pid].append(d[gene_a][0])
        if gene_b in d:
            geno_b[pid].append(d[gene_b][0])
    for pid in list(geno_a.keys()):
        if len(geno_a[pid]) != 2 or pid not in geno_b or len(geno_b[pid]) != 2:
            continue
        genos_by_person[pid] = (tuple(sorted(geno_a[pid])), tuple(sorted(geno_b[pid])))

    print("[37d] people with 2/2 calls at both loci (EM population pool): %s"
          % _sup(len(genos_by_person)), flush=True)

    # ---------------------------------------------------------------------------------------
    # Per-ancestry (+ pooled "ALL") EM fit + truth-set re-phasing.
    # ---------------------------------------------------------------------------------------
    em_summary = []
    switch_rows = []
    spurious_hap_rows = []

    for anc in ANCESTRY_ORDER + ["ALL"]:
        pool_pids = [pid for pid in genos_by_person
                     if anc == "ALL" or anc_of.get(pid) == anc]
        genos = [genos_by_person[pid] for pid in pool_pids]
        n_pool = len(genos)
        if n_pool < SUPPRESS_BELOW:
            em_summary.append({"ancestry": anc, "n_em_pool": _sup(n_pool),
                                "n_truth": "0", "n_truth_doublehet": "0",
                                "switch_errors": "0", "switch_rate": None,
                                "spurious_incompatible_cis": "0",
                                "spurious_incompatible_rate": None,
                                "note": "pool below disclosure floor, EM not run"})
            continue

        freq = em_haplotype_freqs(genos, seed=hash(anc) % (2 ** 31))

        anc_truth_pids = [pid for pid in truth_people
                           if pid in genos_by_person and (anc == "ALL" or anc_of.get(pid) == anc)]
        n_truth = len(anc_truth_pids)
        n_doublehet = 0
        n_switch = 0
        n_spurious = 0
        for pid in anc_truth_pids:
            (a1, a2), (b1, b2) = genos_by_person[pid]
            true_h1, true_h2 = truth_people[pid]["hap1"], truth_people[pid]["hap2"]
            true_pair = sorted([true_h1, true_h2])
            a_het = a1 != a2
            b_het = b1 != b2
            inferred = most_likely_pair(a1, a2, b1, b2, freq)
            inferred_pair = sorted(inferred)
            if a_het and b_het:
                n_doublehet += 1
                if inferred_pair != true_pair:
                    n_switch += 1
                    switch_rows.append({"ancestry": anc, "person_id": "REDACTED",
                                         "true_h1": true_h1, "true_h2": true_h2,
                                         "em_h1": inferred[0], "em_h2": inferred[1]})
            for h in inferred:
                if m37.classify_pair(h[0], h[1]) == "predicted_incompatible":
                    n_spurious += 1
                    spurious_hap_rows.append({"ancestry": anc, "allele_a": h[0], "allele_b": h[1]})

        em_summary.append({
            "ancestry": anc,
            "n_em_pool": _sup(n_pool),
            "n_truth": _sup(n_truth),
            "n_truth_doublehet": _sup(n_doublehet),
            "switch_errors": _sup(n_switch),
            # Rate is only released when BOTH the numerator and the denominator clear the
            # disclosure floor -- a rate computed from a <20 numerator over a disclosed
            # denominator would let a reader back out the exact suppressed count (e.g.
            # 0.0026 x 1926 = 5), defeating the "<20" masking. See feedback: suppressed
            # counts are not zero, and must not be recoverable via a derived ratio either.
            "switch_rate": round(n_switch / n_doublehet, 4)
            if n_doublehet >= SUPPRESS_BELOW and n_switch >= SUPPRESS_BELOW else None,
            "spurious_incompatible_cis": _sup(n_spurious),
            "spurious_incompatible_rate": round(n_spurious / (2 * n_truth), 4)
            if n_truth >= SUPPRESS_BELOW and n_spurious >= SUPPRESS_BELOW else None,
            "note": "",
        })
        print("[37d] %s: pool=%d truth=%d doublehet=%d switch=%d spurious_incompatible=%d"
              % (anc, n_pool, n_truth, n_doublehet, n_switch, n_spurious), flush=True)

    em_summary_df = pd.DataFrame(em_summary)
    em_summary_df.to_csv(os.path.join(args.out_dir, "em_mask_rephase_summary.tsv"),
                          sep="\t", index=False)

    # -----------------------------------------------------------------------------------
    # Corrected per-ancestry (+ pooled) O/E, deliverable 1 -- 37_vm_run.py's own O/E table
    # only ever emitted "ALL" rows because its build_people() ancestry-column lookup picks
    # "max_template_distance" (contains the substring "anc") ahead of "ancestry_pred" in
    # cohort_membership.tsv's column order -- verified directly (its anc_of values were
    # '1.0','2.0','5.0', NaN, ... i.e. a template-distance bucket, not AFR/AMR/...). Redone
    # here with the corrected anc_of, at both resolutions, no separate re-extraction needed
    # for res=2 (haps_df already at args.resolution); res=4 is re-extracted if args.resolution
    # was 2, so both are always reported.
    # -----------------------------------------------------------------------------------
    def _group(gene, allele):
        # same G1/G2 rule as m37.dq_group, called directly to avoid depending on the exact
        # return-signature of m37.observed_expected_incompatible (not re-verified here).
        return m37.dq_group(gene, allele)

    def _oe_for(sub_df, gene_a_, gene_b_):
        nn = len(sub_df)
        if nn == 0:
            return None
        ga = sub_df["allele_a"].map(lambda a: _group(gene_a_, a))
        gb = sub_df["allele_b"].map(lambda b: _group(gene_b_, b))
        mask = ga.notna() & gb.notna()
        ga, gb = ga[mask], gb[mask]
        nnc = len(ga)
        if nnc == 0:
            return None
        p_g1a, p_g2a = (ga == "G1").mean(), (ga == "G2").mean()
        p_g1b, p_g2b = (gb == "G1").mean(), (gb == "G2").mean()
        observed = int((ga != gb).sum())
        expected = nnc * (p_g1a * p_g2b + p_g2a * p_g1b)
        return observed, expected, nnc

    oe_rows = []
    for oe_res, oe_haps in (
        (res, haps_df),
        (4 if res == 2 else 2,
         m37.extract_haplotypes(t1, gene_a, gene_b, 4 if res == 2 else 2)),
    ):
        oe_h = oe_haps.copy()
        oe_h["ancestry"] = oe_h["person_id"].astype(str).map(anc_of)
        for anc in ANCESTRY_ORDER + ["ALL"]:
            sub = oe_h if anc == "ALL" else oe_h[oe_h["ancestry"] == anc]
            row = {"ancestry": anc, "resolution": "%dfield" % oe_res,
                   "n_haplotypes": _sup(len(sub))}
            if len(sub) < SUPPRESS_BELOW:
                row.update({"observed_disp": None, "expected": None, "oe_ratio": None,
                            "note": "below floor"})
                oe_rows.append(row)
                continue
            oe = _oe_for(sub, gene_a, gene_b)
            if oe is None:
                continue
            observed, expected, nnc = oe
            # oe_ratio is only released when observed is exactly 0 (a real, disclosable
            # zero -- not suppressed) or >=20; a <20 nonzero observed with a disclosed
            # expected would let oe_ratio*expected reveal the suppressed count.
            releasable = observed == 0 or observed >= SUPPRESS_BELOW
            row.update({
                "n_classified": _sup(nnc),
                "observed_disp": _sup(observed),
                "expected": round(expected, 2),
                "oe_ratio": round(observed / expected, 4) if releasable and expected > 0 else None,
                "note": "",
            })
            oe_rows.append(row)
    pd.DataFrame(oe_rows).to_csv(
        os.path.join(args.out_dir, "oe_purge_table_by_ancestry_CORRECTED.tsv"),
        sep="\t", index=False)

    # spurious cis-haplotype pairs recurring >=20 times (aggregate only, no person IDs)
    if spurious_hap_rows:
        sdf = pd.DataFrame(spurious_hap_rows)
        agg = (sdf.groupby(["ancestry", "allele_a", "allele_b"]).size()
               .reset_index(name="n").sort_values("n", ascending=False))
        agg["n_disp"] = agg["n"].map(_sup)
        agg.drop(columns=["n"]).to_csv(
            os.path.join(args.out_dir, "em_spurious_incompatible_pairs.tsv"),
            sep="\t", index=False)

    # --- real-world different-contig implied-pairing summary -------------------------------
    dc_summary = []
    if not diffcontig_df.empty:
        for anc in ANCESTRY_ORDER + ["ALL"]:
            d = diffcontig_df if anc == "ALL" else diffcontig_df[diffcontig_df["ancestry"] == anc]
            n_total = len(d)
            n_incompat = int(d["incompatible"].sum())
            dc_summary.append({
                "ancestry": anc,
                "n_diffcontig_haps": _sup(n_total),
                "n_diffcontig_incompatible": _sup(n_incompat),
                # Same disclosure guard as switch_rate/spurious_incompatible_rate above:
                # only release the ratio when the numerator itself also clears the floor.
                "rate": round(n_incompat / n_total, 4)
                if n_total >= SUPPRESS_BELOW and n_incompat >= SUPPRESS_BELOW else None,
            })
    dc_df = pd.DataFrame(dc_summary)
    dc_df.to_csv(os.path.join(args.out_dir, "diffcontig_implied_pairing.tsv"),
                 sep="\t", index=False)

    # --- figure: switch-error rate + spurious-incompatible rate per ancestry ---------------
    fig_path = os.path.join(args.out_dir, "fig_mask_rephase_em")
    try:
        _fig_bar(em_summary_df, dc_df, fig_path)
    except Exception as e:  # noqa: BLE001 -- figure is best-effort, numbers already saved
        print("[37d] figure render failed (numbers still saved): %r" % e, flush=True)

    print("[37d] DONE. wrote em_mask_rephase_summary.tsv, "
          "em_spurious_incompatible_pairs.tsv, diffcontig_implied_pairing.tsv", flush=True)


def _fig_bar(em_df, dc_df, path_stem):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_df = em_df[em_df["switch_rate"].notna()]
    with vc.nature_style():
        fig, axes = plt.subplots(1, 2, figsize=(vc.mm(183), vc.mm(70)))

        ax = axes[0]
        ax.bar(plot_df["ancestry"], plot_df["switch_rate"].astype(float), color="#4C72B0")
        ax.set_ylabel("EM switch-error rate\n(doubly-het truth people)")
        ax.set_title("Mask-and-rephase: switch errors", fontsize=6)

        ax2 = axes[1]
        ax2.bar(plot_df["ancestry"], plot_df["spurious_incompatible_rate"].astype(float),
                color="#C44E52")
        ax2.set_ylabel("Spurious G1/G2-incompatible\ncis haplotypes / (2 x truth N)")
        ax2.set_title("Statistical phasing manufactures\nincompatible haplotypes (truth=0)",
                       fontsize=6)

        for a in axes:
            a.spines["top"].set_visible(False)
            a.spines["right"].set_visible(False)
        fig.tight_layout()
        vc.save_fig(fig, path_stem)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=os.path.expanduser("~/s03/results/37"))
    ap.add_argument("--resolution", type=int, default=2, choices=[2, 4])
    args = ap.parse_args()
    run(args)


if __name__ == "__main__":
    main()
