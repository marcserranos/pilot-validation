#!/usr/bin/env python3
"""S03 critic #1 fix (blocker #2 + major #3): rerun 37's O/E and mask-and-rephase EM restricted
to the CORRECT unrelated cohort (matching scripts 29/38's definition), and make the
mask-and-rephase EM experiment non-circular.

Two problems this fixes, both flagged by critic #1 (CRITIC_1.md):

1. **Wrong unrelated set.** `37_vm_run.py`'s own `build_people()` starts from ALL of
   `cohort_membership.tsv` (13,252 rows, relatives included) and greedily drops one member of
   each related PAIR in isolation -- not the maximal-independent-set algorithm
   (`24_novelty_by_field.greedy_unrelated`) every other S03 result (29, 38, 39) uses on the
   actual Table-1/LR-called cohort. This script re-derives the unrelated set the CORRECT way:
   `24_novelty_by_field.build_people(t1["person_id"], ...)`, which calls `greedy_unrelated` on
   Table-1's own person_ids against the real relatedness graph. Expect ~11,856 (38's number),
   not 13,252.

2. **Circularity in the original mask-and-rephase EM (37d).** The EM population pool for a given
   ancestry INCLUDED the truth-set people it was then evaluated against -- "0 spurious
   incompatible haplotypes" could just mean "the EM was fit on data that already has none",
   not "EM is safe in general". Fixed two ways, both non-circular:
   (a) **K-fold held-out EM**: split each ancestry's EM pool into K folds; for each fold, fit EM
       frequencies on the OTHER folds only, then re-phase truth-set people who fall in the
       held-out fold. Pooled across folds, no truth-set person is ever rephased using a model
       that saw their own data.
   (b) **Naive linkage-equilibrium (LE) baseline**: instead of the EM-fit joint frequencies, use
       independent per-locus marginal allele frequencies (p_i x q_j) learned from the training
       fold. Under strict independence, both cis/trans resolutions of a doubly-het genotype have
       EXACTLY equal posterior weight (p(a1)q(b1)p(a2)q(b2) == p(a1)q(b2)p(a2)q(b1)) -- a formal
       tie -- so this baseline breaks ties with a fixed-seed coin flip. This bounds how many
       incompatible haplotypes a phaser with ZERO real linkage information would manufacture,
       i.e. the "prior does all the work" end of the spectrum the EM sits somewhere below.

Usage (VM):
    setsid nohup python3 -u 37e_unrelated_fix_kfold_em.py --out-dir ~/s03/results/37 \
        --resolution 2 --k-folds 5 < /dev/null > ~/s03/results/37/37e_run.log 2>&1 & disown
"""
import argparse
import hashlib
import importlib.util
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


def _fold_of(pid, k_folds):
    """Deterministic fold assignment: stable across runs/processes (unlike hash(), which is
    salted per-process in Python by default) -- important since the same fold assignment must be
    reproducible if this script is ever rerun/audited."""
    h = int(hashlib.sha256(str(pid).encode()).hexdigest(), 16)
    return h % k_folds


# ---------------------------------------------------------------------------
# EM haplotype-frequency estimator (Excoffier & Slatkin 1995), reused verbatim from 37d.
# ---------------------------------------------------------------------------
def em_haplotype_freqs(genotypes, max_iter=EM_MAX_ITER, tol=EM_TOL, seed=0):
    rng = random.Random(seed)  # noqa: F841 -- kept for signature parity with 37d, unused (freq
    # init is uniform/deterministic, no randomness needed in the EM itself)
    a_alleles = sorted({a for (a1, a2), _ in genotypes for a in (a1, a2)})
    b_alleles = sorted({b for _, (b1, b2) in genotypes for b in (b1, b2)})
    haps = [(a, b) for a in a_alleles for b in b_alleles]
    idx = {h: i for i, h in enumerate(haps)}
    freq = np.full(len(haps), 1.0 / len(haps))

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
            weights = [freq[idx[h1]] * freq[idx[h2]] for h1, h2 in res]
            tot = sum(weights)
            if tot <= 0:
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
    a_het = a1 != a2
    b_het = b1 != b2
    if not a_het or not b_het:
        return [(a1, b1), (a2, b2)]
    opt1 = [(a1, b1), (a2, b2)]
    opt2 = [(a1, b2), (a2, b1)]
    w1 = freq.get(opt1[0], 0.0) * freq.get(opt1[1], 0.0)
    w2 = freq.get(opt2[0], 0.0) * freq.get(opt2[1], 0.0)
    return opt1 if w1 >= w2 else opt2


def marginal_freqs(genotypes):
    """Per-locus marginal allele frequencies from unphased genotypes (allele-copy counting).
    Returns (p_dict for locus A, q_dict for locus B)."""
    ca, cb = defaultdict(int), defaultdict(int)
    for (a1, a2), (b1, b2) in genotypes:
        ca[a1] += 1
        ca[a2] += 1
        cb[b1] += 1
        cb[b2] += 1
    na, nb = sum(ca.values()), sum(cb.values())
    p = {a: c / na for a, c in ca.items()} if na else {}
    q = {b: c / nb for b, c in cb.items()} if nb else {}
    return p, q


def le_naive_pair(a1, a2, b1, b2, p, q, pid, seed):
    """Naive linkage-equilibrium phaser: under strict independence, freq(a,b)=p(a)*q(b), so both
    cis/trans resolutions of a doubly-het genotype are an EXACT tie
    (p(a1)q(b1)p(a2)q(b2) == p(a1)q(b2)p(a2)q(b1)) -- there is no information to break the tie
    with, so this breaks it with a per-person, seeded coin flip (never a fixed default side,
    which would silently bias the "spurious incompatible" count high or low by construction)."""
    a_het = a1 != a2
    b_het = b1 != b2
    if not a_het or not b_het:
        return [(a1, b1), (a2, b2)]
    rng = random.Random("%r:%s" % (seed, pid))
    return [(a1, b1), (a2, b2)] if rng.random() < 0.5 else [(a1, b2), (a2, b1)]


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def run(args):
    m37 = _load_module("37_vm_run.py", "s03_37_vm_run")
    m24 = _load_module("24_novelty_by_field.py", "s03_24_novelty")
    vc.ensure_dir(args.out_dir)

    if not os.path.exists(m37.TABLE1):
        sys.exit("FATAL: %s not found -- needs a live VM session." % m37.TABLE1)

    print("[37e] loading table1 ...", flush=True)
    t1 = pd.read_csv(m37.TABLE1, sep="\t", dtype=str, low_memory=False)
    all_pids = t1["person_id"].astype(str).unique().tolist()

    # --- CORRECT unrelated set: same definition scripts 29/38/39 use ----------------------
    people_df, n_removed = m24.build_people(
        all_pids, m37.COHORT, m37.RELATEDNESS, m37.KIN_MIN, m37.STRICT_MIN, skip_rel=False)
    unrelated_ids = set(people_df.loc[people_df["unrelated"], "person_id"].astype(str))
    anc_of = dict(zip(people_df["person_id"].astype(str), people_df["anc_pred"]))
    print("[37e] Table-1 people: %s ; CORRECTED unrelated (greedy_unrelated, matches 29/38/39): "
          "%s ; removed for relatedness: %s"
          % (_sup(len(all_pids)), _sup(len(unrelated_ids)), _sup(n_removed)), flush=True)

    t1 = t1[t1["person_id"].astype(str).isin(unrelated_ids)].copy()

    res = args.resolution
    gene_a, gene_b = "DQA1", "DQB1"

    # --- physically-phased cis pairs (truth set) -------------------------------------------
    haps_df = m37.extract_haplotypes(t1, gene_a, gene_b, res)
    haps_df["ancestry"] = haps_df["person_id"].astype(str).map(anc_of)

    per_person_phys = defaultdict(dict)
    for _, r in haps_df.iterrows():
        per_person_phys[r["person_id"]][r["hap"]] = (r["allele_a"], r["allele_b"])
    truth_people = {pid: h for pid, h in per_person_phys.items()
                     if set(h.keys()) >= {"hap1", "hap2"}}
    print("[37e] truth-set people (both haps physically phased, res=%d): %s"
          % (res, _sup(len(truth_people))), flush=True)

    # --- assembly-labelled per-hap calls (for genotype pool + O/E) -------------------------
    bare = t1["gene"].astype(str).str.replace("^HLA-", "", regex=True)
    sub = t1[bare.isin([gene_a, gene_b])].copy()
    sub["gene_bare"] = bare[sub.index]
    per_hap_calls = defaultdict(dict)
    for (pid, hap), grp in sub.groupby(["person_id", "hap"], sort=False):
        for gene in (gene_a, gene_b):
            grows = grp[grp["gene_bare"] == gene]
            if len(grows) != 1:
                continue
            allele = m37.to_nfield(grows.iloc[0]["consensus"], res)
            if allele is None:
                continue
            per_hap_calls[(pid, hap)][gene] = (allele, grows.iloc[0]["contig"])

    genos_by_person = {}
    geno_a, geno_b = defaultdict(list), defaultdict(list)
    for (pid, hap), d in per_hap_calls.items():
        if gene_a in d:
            geno_a[pid].append(d[gene_a][0])
        if gene_b in d:
            geno_b[pid].append(d[gene_b][0])
    for pid in list(geno_a.keys()):
        if len(geno_a[pid]) != 2 or pid not in geno_b or len(geno_b[pid]) != 2:
            continue
        genos_by_person[pid] = (tuple(sorted(geno_a[pid])), tuple(sorted(geno_b[pid])))
    print("[37e] people with 2/2 calls at both loci (EM population pool), corrected unrelated "
          "set: %s" % _sup(len(genos_by_person)), flush=True)

    # =========================================================================================
    # Corrected O/E (deliverable 1), restricted to the CORRECT unrelated set.
    # =========================================================================================
    def _group(gene, allele):
        return m37.dq_group(gene, allele)

    def _oe_for(sub_df, gene_a_, gene_b_):
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
        (4 if res == 2 else 2, m37.extract_haplotypes(t1, gene_a, gene_b, 4 if res == 2 else 2)),
    ):
        oe_h = oe_haps.copy()
        oe_h["ancestry"] = oe_h["person_id"].astype(str).map(anc_of)
        for anc in ANCESTRY_ORDER + ["ALL"]:
            s = oe_h if anc == "ALL" else oe_h[oe_h["ancestry"] == anc]
            row = {"ancestry": anc, "resolution": "%dfield" % oe_res, "n_haplotypes": _sup(len(s))}
            if len(s) < SUPPRESS_BELOW:
                row.update({"observed_disp": None, "expected": None, "oe_ratio": None,
                            "note": "below floor"})
                oe_rows.append(row)
                continue
            oe = _oe_for(s, gene_a, gene_b)
            if oe is None:
                continue
            observed, expected, nnc = oe
            releasable = observed == 0 or observed >= SUPPRESS_BELOW
            row.update({
                "n_classified": _sup(nnc), "observed_disp": _sup(observed),
                "expected": round(expected, 2),
                "oe_ratio": round(observed / expected, 4) if releasable and expected > 0 else None,
                "note": "",
            })
            oe_rows.append(row)
    pd.DataFrame(oe_rows).to_csv(
        os.path.join(args.out_dir, "oe_purge_table_by_ancestry_UNRELATED_FIXED.tsv"),
        sep="\t", index=False)

    # =========================================================================================
    # K-fold held-out EM + naive-LE baseline mask-and-rephase, per ancestry + pooled ALL.
    # =========================================================================================
    K = args.k_folds
    kfold_rows = []
    le_rows = []
    spurious_hap_rows_kfold = []
    spurious_hap_rows_le = []

    for anc in ANCESTRY_ORDER + ["ALL"]:
        pool_pids = [pid for pid in genos_by_person if anc == "ALL" or anc_of.get(pid) == anc]
        n_pool = len(pool_pids)
        anc_truth_pids = [pid for pid in truth_people
                           if pid in genos_by_person and (anc == "ALL" or anc_of.get(pid) == anc)]
        n_truth = len(anc_truth_pids)

        if n_pool < SUPPRESS_BELOW or n_truth < SUPPRESS_BELOW:
            for rows, label in ((kfold_rows, "kfold"), (le_rows, "le_naive")):
                rows.append({"ancestry": anc, "scheme": label, "n_em_pool": _sup(n_pool),
                             "n_truth_evaluated": _sup(n_truth), "n_truth_doublehet": "0",
                             "switch_errors": "0", "switch_rate": None,
                             "spurious_incompatible_cis": "0", "spurious_incompatible_rate": None,
                             "note": "pool or truth below disclosure floor"})
            continue

        fold_of_pid = {pid: _fold_of(pid, K) for pid in pool_pids}

        for scheme in ("kfold", "le_naive"):
            n_doublehet = n_switch = n_spurious = n_eval = 0
            for f in range(K):
                held_out_truth = [pid for pid in anc_truth_pids if fold_of_pid.get(pid) == f]
                if not held_out_truth:
                    continue
                train_pids = [pid for pid in pool_pids if fold_of_pid[pid] != f]
                train_genos = [genos_by_person[pid] for pid in train_pids]
                if len(train_genos) < 2:
                    continue
                if scheme == "kfold":
                    freq = em_haplotype_freqs(train_genos, seed=hash((anc, f)) % (2 ** 31))
                else:
                    p, q = marginal_freqs(train_genos)

                for pid in held_out_truth:
                    (a1, a2), (b1, b2) = genos_by_person[pid]
                    true_h1, true_h2 = truth_people[pid]["hap1"], truth_people[pid]["hap2"]
                    true_pair = sorted([true_h1, true_h2])
                    a_het, b_het = a1 != a2, b1 != b2
                    if scheme == "kfold":
                        inferred = most_likely_pair(a1, a2, b1, b2, freq)
                    else:
                        inferred = le_naive_pair(a1, a2, b1, b2, p, q, pid, seed=args.seed)
                    inferred_pair = sorted(inferred)
                    n_eval += 1
                    if a_het and b_het:
                        n_doublehet += 1
                        if inferred_pair != true_pair:
                            n_switch += 1
                    for h in inferred:
                        if m37.classify_pair(h[0], h[1]) == "predicted_incompatible":
                            n_spurious += 1
                            (spurious_hap_rows_kfold if scheme == "kfold"
                             else spurious_hap_rows_le).append(
                                {"ancestry": anc, "allele_a": h[0], "allele_b": h[1]})

            row = {
                "ancestry": anc, "scheme": scheme, "n_em_pool": _sup(n_pool),
                "n_truth_evaluated": _sup(n_eval), "n_truth_doublehet": _sup(n_doublehet),
                "switch_errors": _sup(n_switch),
                "switch_rate": round(n_switch / n_doublehet, 4)
                if n_doublehet >= SUPPRESS_BELOW and n_switch >= SUPPRESS_BELOW else None,
                "spurious_incompatible_cis": _sup(n_spurious),
                "spurious_incompatible_rate": round(n_spurious / (2 * n_eval), 4)
                if n_eval >= SUPPRESS_BELOW and n_spurious >= SUPPRESS_BELOW else None,
                "note": "",
            }
            (kfold_rows if scheme == "kfold" else le_rows).append(row)
            print("[37e] %s/%s: pool=%d eval=%d doublehet=%d switch=%d spurious_incompatible=%d"
                  % (anc, scheme, n_pool, n_eval, n_doublehet, n_switch, n_spurious), flush=True)

    kfold_df = pd.DataFrame(kfold_rows)
    le_df = pd.DataFrame(le_rows)
    combined = pd.concat([kfold_df, le_df], ignore_index=True)
    combined.to_csv(os.path.join(args.out_dir, "vm_em_mask_rephase_kfold_and_le.tsv"),
                     sep="\t", index=False)

    for rows, name in ((spurious_hap_rows_kfold, "vm_spurious_pairs_kfold.tsv"),
                        (spurious_hap_rows_le, "vm_spurious_pairs_le_naive.tsv")):
        if rows:
            sdf = pd.DataFrame(rows)
            agg = (sdf.groupby(["ancestry", "allele_a", "allele_b"]).size()
                   .reset_index(name="n").sort_values("n", ascending=False))
            agg["n_disp"] = agg["n"].map(_sup)
            agg.drop(columns=["n"]).to_csv(os.path.join(args.out_dir, name), sep="\t", index=False)

    # --- figure: kfold vs LE-naive spurious-incompatible rate, per ancestry ----------------
    fig_path = os.path.join(args.out_dir, "vm_fig_kfold_vs_le")
    try:
        _fig_kfold_vs_le(kfold_df, le_df, fig_path)
    except Exception as e:  # noqa: BLE001
        print("[37e] figure render failed (numbers still saved): %r" % e, flush=True)

    print("[37e] DONE. wrote oe_purge_table_by_ancestry_UNRELATED_FIXED.tsv, "
          "vm_em_mask_rephase_kfold_and_le.tsv", flush=True)


def _fig_kfold_vs_le(kfold_df, le_df, path_stem):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    kf = kfold_df[kfold_df["spurious_incompatible_rate"].notna()]
    le = le_df[le_df["spurious_incompatible_rate"].notna()]
    ancs = [a for a in ANCESTRY_ORDER + ["ALL"]
            if a in set(kf["ancestry"]) | set(le["ancestry"])]
    with vc.nature_style():
        fig, ax = plt.subplots(figsize=(vc.mm(120), vc.mm(70)))
        x = np.arange(len(ancs))
        width = 0.35
        kf_vals = [float(kf.set_index("ancestry")["spurious_incompatible_rate"].get(a, np.nan))
                   for a in ancs]
        le_vals = [float(le.set_index("ancestry")["spurious_incompatible_rate"].get(a, np.nan))
                   for a in ancs]
        ax.bar(x - width / 2, kf_vals, width, label="K-fold held-out EM", color="#4C72B0")
        ax.bar(x + width / 2, le_vals, width, label="Naive linkage-equilibrium (no LD info)",
               color="#C44E52")
        ax.set_xticks(x)
        ax.set_xticklabels(ancs, fontsize=6)
        ax.set_ylabel("Spurious G1/G2-incompatible cis haplotypes\n/ (2 x held-out truth N)")
        ax.set_title("Non-circular mask-and-rephase: held-out EM vs. no-LD-information baseline",
                      fontsize=6)
        ax.legend(fontsize=5.5, frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        vc.save_fig(fig, path_stem)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=os.path.expanduser("~/s03/results/37"))
    ap.add_argument("--resolution", type=int, default=2, choices=[2, 4])
    ap.add_argument("--k-folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    run(args)


if __name__ == "__main__":
    main()
