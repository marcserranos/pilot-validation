#!/usr/bin/env python3
"""WS-A (Sprint S04) -- KIR recurrence + saturation, on the same footing as HLA (WS-B setup).

Repeats for KIR the recurrence-class and discovery/saturation analyses already done for HLA
(`03_novel_alleles.py`, `34_novel_protein_recurrence.py`, `24_novelty_by_field.py`,
`04_allele_saturation.py`, `18_discovery_rate.py`, and the corrected equal-N per-ancestry
methodology in `39_saturation_by_ancestry.py` -- its README documents the S03 fixes and is
followed here, not re-derived), and runs the IDENTICAL analysis on HLA at the same time, on the
same unrelated people (KIR ∩ HLA, then one shared `greedy_unrelated()` pass) and the same
ancestry scheme and subsampling seeds, so WS-B (script 46, KIR-vs-HLA catalogue coverage) can
compare the two regions on equal footing without re-deriving anything.

REUSE (imported via importlib, not re-derived):
  - `41_kir_pilot.py`: parse_hap_gtf, classify_novelty_tier, greedy_unrelated,
    load_relatedness_pairs, KIN_MIN, KIR_GENES, ANCESTRY_ORDER.
  - `43_kir_full_aggregate.py`: discover_people, parse_all/parse_person (multiprocessing GTF
    parse, returns {person_id, haps: [{genes_called, kir_rows}], n_invalid_haps}), suppressed(),
    rate_row(), wilson_ci(), load_ancestry_map().
  - `39_saturation_by_ancestry.py`: build_labeled_calls (the full HLA call-labeling pipeline --
    field_class/seq_class/prot_id/cds_id/keep_clean, ancestry + unrelated flags), its own imports
    of `24_novelty_by_field.py` (m24) and `04_allele_saturation.py` (m04: incidence_freqs, chao2,
    chao2_extrapolate, clench_curve/fit_clench_asymptote).

NOVELTY LEVELS (both species, both counted from the SAME per-call identity fields already
computed by the reused pipelines -- no new novelty classification logic here):
  - "all"           -- every clean call (known + novel), for a saturation baseline.
  - "any_novel"     -- any sequence difference from the reference catalogue at all.
                       KIR: `classify_novelty_tier(row) != "known"` (excludes "undetermined").
                       HLA: `field_class != "known"` (excludes "uncalled"), i.e. f2_protein /
                       f3_synonymous / f4_noncoding pooled -- this is 24_novelty_by_field.py's own
                       "any depth of novelty" grouping, not a new definition.
  - "protein_novel" -- a real amino-acid-level novel protein.
                       KIR: `novelty_tier == "novel_protein"` (Immuannot's own cds_mut-derived
                       call, see 41_kir_pilot.classify_novelty_tier docstring).
                       HLA: `field_class == "f2_protein" AND seq_class == "novel_protein"` --
                       exactly 39_saturation_by_ancestry's own "novel" category, reused verbatim.

ALLELE IDENTITY (documented granularity mismatch between species -- flagged, not hidden):
  - KIR: the Immuannot `consensus` string for a call, at ALL three levels. This matches
    43_kir_full_aggregate.py's own convention (`gene_known_allele_counts`/`gene_novel_protein_seqs`
    both key on `consensus`) -- i.e. genomic-consensus granularity (introns/UTR included), since
    Immuannot's KIR output does not carry a separate coarser field-level allele name the way HLA
    nomenclature does.
  - HLA "all"/"any_novel": `cds_id` (CDS-level identity; 24_novelty_by_field.allele_ids() already
    computes this per seq_class branch -- see that function's docstring). HLA "protein_novel":
    `prot_id` (2-field protein-level identity), matching 39's own "novel" category exactly.
  CAVEAT: KIR's genomic-consensus identity is a FINER granularity than HLA's cds_id (which
  collapses synonymous/intronic differences within a matched CDS core less aggressively than a
  literal genomic string would -- see 24's allele_ids() branches). This is the finest identity
  each project's existing pipeline already computes without building new machinery for either
  species; WS-B (script 46) must not read raw novelty-% or recurrence-class counts across species
  as perfectly like-for-like without repeating this caveat.

DISCLOSURE (hard rule, same as 39/43 -- see AGENT_PREAMBLE.md / ORCHESTRATOR_HANDOFF.md sect 3):
  - Every raw count 1-19 is written as the literal string "<20" (never 0, never omitted). A true
    0 stays "0". Applies to recurrence-class ALLELE counts too, even though those are allele
    counts rather than participant counts -- per this sprint's own explicit instruction ("still
    write 1-19 as '<20' to be safe"): an allele carried by exactly 1-19 people is disclosive about
    those specific people's rare genotype regardless of whether the number being reported is
    labeled "people" or "alleles".
  - No allele name/string (KIR consensus or HLA prot_id/cds_id) is ever written next to a
    carrier count < 20.
  - Any rate whose numerator or denominator is 1-19 is blanked (point estimate AND CI).
  - Good-Turing coverage / Chao2 are computed on-VM from f1/f2 (Q1/Q2), which can themselves be
    small; only the resulting coverage/Chao2/SE ESTIMATES are exported, never the raw f1/f2 pair,
    per this sprint's instruction ("compute Chao on VM, export only the estimate if f1/f2 are
    small"). f1/f2 are written to the export TSV only when both are >= SUPPRESS_BELOW.

RECURRENCE CLASSES (allele counts, not mutually exclusive by design -- documented, not a bug):
  eq1  -- alleles carried by EXACTLY 1 unrelated person.
  eq2  -- alleles carried by EXACTLY 2 unrelated people.
  gt2  -- alleles carried by >= 3 unrelated people (a superset that INCLUDES ge20).
  ge20 -- alleles carried by >= 20 unrelated people (a subset of gt2, called out separately
          because 20 is this project's own disclosure floor -- ORCHESTRATOR_HANDOFF.md WS-A).

STATISTICS:
  - Good-Turing sample coverage, incidence-data (Chao & Jost 2012, Ecol. 93:2533-2547) bias
    correction: C_hat = 1 - (Q1/n) * [(m-1)Q1 / ((m-1)Q1 + 2Q2)] (Q2>0 case), n = total number of
    incidences (sum of per-allele occurrence counts across sampling units), m = number of units.
  - Chao2 richness + its classic incidence-based SE (Chao 1987, Biometrics 43:783-791; formulas as
    reproduced in Colwell/EstimateS's user guide eq. 21b/23b) -- `chao2()` itself imported
    verbatim from 04_allele_saturation.py; the SE is new (04/39 never computed/exported it) and is
    unit-tested here for non-negativity and the Q2==0 bias-corrected branch.
  - Equal-N discovery slope: linear fit of the mean rarefaction curve over the window
    [0.9*N*, N*], exactly 39_saturation_by_ancestry's own "primary" cross-ancestry-comparable
    number (see that script's README) -- generalized here to ALSO run per recurrence class, not
    just the >=1-carrier curve.
  - Rarefaction/discovery curves: adapted from 39_saturation_by_ancestry.rarefaction_threshold_curves
    (itself already generic over an arbitrary `thresholds` list) to accept a PRECOMPUTED list of
    permutation orders instead of generating its own -- so KIR and HLA curves, within one ancestry,
    walk the EXACT SAME people in the EXACT SAME order across all `n_permutations` replicates (the
    brief's "same subsampling seeds" requirement for WS-B). thresholds=[1,2,3,20]; exact-class
    curves are derived by differencing (eq1=thr[1]-thr[2], eq2=thr[2]-thr[3], gt2=thr[3],
    ge20=thr[20]).
  - Optional iNEXT-style 2N extrapolation via 04's `chao2_extrapolate` (cheap -- same Q1/f0_hat
    already computed for Chao2; included when --extrapolate-2n is passed, default on).

OUTPUTS (--out-dir, default ~/s04/results/44/; only aggregates, safe to pull off the VM):
  recurrence_classes.tsv   gene, ancestry, level, species, eq1, eq2, gt2, ge20  (masked <20)
  saturation_curves.tsv    gene, ancestry, level, species, n, mean_distinct, lo2_5, hi97_5,
                           mean_eq1, mean_eq2, mean_gt2, mean_ge20 (rarefaction curve points;
                           n=1..N per (gene,ancestry,level,species), subsampled to
                           --curve-stride points to keep the file small)
  equal_n_slope.tsv        gene, ancestry, level, species, n_star, slope_per_1000,
                           mean_distinct_at_n_star
  coverage_chao2.tsv       gene, ancestry, level, species, n_people, s_obs, good_turing_coverage,
                           chao2, chao2_se, chao2_undetected_f0hat, chao_new_by_2n (if requested),
                           q1, q2 (blank unless both >= SUPPRESS_BELOW)
  STATUS.txt               aggregate-progress-only status file, rewritten as the run proceeds.

USAGE (VM; see the module's own --help for every flag):
  cd ~/s04 && PYTHONPATH=~/s04:~/s03:~/repos/pilot-validation/scripts/hla_popgen \\
    python3 -u 44_kir_recurrence_saturation.py \\
      --kir-outroot ~/pipeline_outputs_kir \\
      --hla-table ~/pipeline_outputs/hla_calls_rich.tsv \\
      --cohort-membership ~/pipeline_outputs/cohort_membership.tsv \\
      --relatedness-table ~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv \\
      --refdata ~/tools/Immuannot_refdata \\
      --out-dir ~/s04/results/44 --workers 4 2>&1 | tee -a ~/s04/results/44/run.log

Never writes anywhere under ~/pipeline_outputs* (read-only inputs; verify with a pre/post md5
check of the 14 HLA production tables per this project's standing rule, same as 43's run did).
"""
import argparse
import importlib.util
import math
import os
import sys
import time
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SUPPRESS_BELOW = 20
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
LEVELS = ["all", "any_novel", "protein_novel"]
SPECIES = ["kir", "hla"]
KIR_THRESHOLDS = [1, 2, 3, 20]
N_PERMUTATIONS_DEFAULT = 25       # Pakistan Fig 3e convention, matches 39's default
MIN_PEOPLE_PER_ANCESTRY = 100      # "well-powered" floor for N*, matches 39
CURVE_STRIDE_DEFAULT = 25          # export every Nth curve point (keeps saturation_curves.tsv small)


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def _load_module(filename, modname):
    path = os.path.join(_THIS_DIR, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


_m41 = _m43 = _m39 = _m04 = None


def m41():
    global _m41
    if _m41 is None:
        _m41 = _load_module("41_kir_pilot.py", "kir_pilot_44")
    return _m41


def m43():
    global _m43
    if _m43 is None:
        _m43 = _load_module("43_kir_full_aggregate.py", "kir_full_aggregate_44")
        _m43._kir = m41()
    return _m43


def m39():
    global _m39
    if _m39 is None:
        _m39 = _load_module("39_saturation_by_ancestry.py", "saturation_by_ancestry_44")
    return _m39


def m04():
    global _m04
    if _m04 is None:
        _m04 = _load_module("04_allele_saturation.py", "allele_saturation_44")
    return _m04


# ---------------------------------------------------------------------------
# Disclosure helpers (reused verbatim from 43, not re-derived).
# ---------------------------------------------------------------------------
def suppressed(n):
    return m43().suppressed(n)


def rate_row(n, d, pct_decimals=1):
    return m43().rate_row(n, d, pct_decimals)


def write_status(out_dir, msg):
    try:
        with open(os.path.join(out_dir, "STATUS.txt"), "a") as f:
            f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}\n")
    except OSError as e:
        log(f"WARNING: could not write STATUS.txt: {e}")


# ---------------------------------------------------------------------------
# Core, pure, unit-tested statistics -- shared by both species.
# ---------------------------------------------------------------------------
def carrier_counts(unit_sets):
    """unit_sets: dict person_id -> set(allele_id). Returns Counter(allele_id -> n_carriers)."""
    c = Counter()
    for s in unit_sets.values():
        for a in s:
            c[a] += 1
    return c


def classify_recurrence(counts):
    """counts: Counter/dict allele_id -> n_carriers (unrelated people). Returns raw (unmasked)
    dict eq1/eq2/gt2/ge20 -> n_alleles. gt2 and ge20 overlap by design -- see module docstring."""
    eq1 = eq2 = gt2 = ge20 = 0
    for n in counts.values():
        if n == 1:
            eq1 += 1
        elif n == 2:
            eq2 += 1
        else:
            gt2 += 1  # n >= 3, includes n >= 20
        if n >= 20:
            ge20 += 1
    return {"eq1": eq1, "eq2": eq2, "gt2": gt2, "ge20": ge20}


def incidence_freqs(unit_sets):
    """dict person_id -> set(allele_id) form of 04_allele_saturation.incidence_freqs (which takes
    a list of sets) -- thin adapter, not a re-derivation."""
    return m04().incidence_freqs(list(unit_sets.values()))


def good_turing_coverage(Q, m, n_incidences):
    """Chao & Jost (2012), Ecology 93:2533-2547 -- incidence-data sample coverage estimator.
    Q: Counter mapping occurrence-count -> number of alleles with that occurrence count (Q[1]=f1
    uniques, Q[2]=f2 duplicates). m: number of sampling units. n_incidences: total number of
    (unit, allele) incidences = sum(occ.values()) in 04's incidence_freqs() terms.
    C_hat = 1 - (Q1/n) * [(m-1)Q1 / ((m-1)Q1 + 2Q2)]           if Q2 > 0
          = 1 - (Q1/n) * [(m-1)(Q1-1) / ((m-1)(Q1-1) + 2)]     if Q2 == 0, Q1 > 1
          = 1                                                   if Q2 == 0, Q1 <= 1 (no unseen mass)
    """
    Q1, Q2 = Q.get(1, 0), Q.get(2, 0)
    if n_incidences <= 0:
        return float("nan")
    if Q2 > 0:
        corr = ((m - 1) * Q1) / ((m - 1) * Q1 + 2 * Q2)
    elif Q1 > 1:
        corr = ((m - 1) * (Q1 - 1)) / ((m - 1) * (Q1 - 1) + 2)
    else:
        corr = 0.0
    return 1.0 - (Q1 / n_incidences) * corr


def chao2_se(S_obs, Q, m):
    """Incidence-based Chao2 standard error (Chao 1987, Biometrics 43:783-791; formulas as
    reproduced in Colwell's EstimateS User's Guide sect. "Variance of Chao2"). Two branches,
    matching 04_allele_saturation.chao2()'s own Q2>0 / Q2==0 split so the SE is always consistent
    with whichever point-estimate branch chao2() used."""
    Q1, Q2 = Q.get(1, 0), Q.get(2, 0)
    if m <= 1:
        return 0.0
    A = (m - 1) / m
    if Q2 > 0:
        r = Q1 / Q2
        var = Q2 * (0.5 * A * r ** 2 + A ** 2 * r ** 3 + 0.25 * A ** 2 * r ** 4)
    else:
        s_c = S_obs + A * Q1 * (Q1 - 1) / 2.0
        var = (A * Q1 * (Q1 - 1) / 2.0
               + A ** 2 * Q1 * (2 * Q1 - 1) ** 2 / 4.0
               - (A ** 2 * Q1 ** 4 / (4.0 * s_c) if s_c > 0 else 0.0))
    return math.sqrt(max(var, 0.0))


def make_orders(n_units, n_permutations, seed):
    """Precomputed list of permutation index-arrays over range(n_units), one draw per replicate.
    Sharing ONE such list between the KIR and HLA extraction for a given ancestry (same n_units,
    same pids order feeding unit_sets_list) is what guarantees identical subsampling paths across
    species -- the brief's WS-B requirement."""
    rng = np.random.default_rng(seed)
    return [rng.permutation(n_units) for _ in range(n_permutations)]


def rarefaction_curves_fixed_orders(unit_sets_list, orders, thresholds=KIR_THRESHOLDS):
    """Adapted from 39_saturation_by_ancestry.rarefaction_threshold_curves: identical inner loop,
    generalized to take PRECOMPUTED permutation orders (see make_orders()) instead of drawing its
    own, so two calls with the same `orders` (one per species) walk exactly the same people in
    exactly the same sequence. unit_sets_list: list of sets, index-aligned with `orders`' values
    (i.e. unit_sets_list[i] is person i's allele set for whichever gene/level/species is being
    curved). Returns {'distinct': (n_perm, m) array, thresholds[k]: (n_perm, m) array}."""
    m = len(unit_sets_list)
    n_perm = len(orders)
    distinct = np.zeros((n_perm, m))
    thr = {k: np.zeros((n_perm, m)) for k in thresholds}
    for p, order in enumerate(orders):
        carrier_count = {}
        n_distinct = 0
        n_reached = {k: 0 for k in thresholds}
        for step, i in enumerate(order):
            for allele in unit_sets_list[i]:
                c = carrier_count.get(allele, 0)
                if c == 0:
                    n_distinct += 1
                c += 1
                carrier_count[allele] = c
                if c in n_reached:
                    n_reached[c] += 1
            distinct[p, step] = n_distinct
            for k in thresholds:
                thr[k][p, step] = n_reached[k]
    out = {"distinct": distinct}
    out.update(thr)
    return out


def summarize_curve(arr):
    """arr: (n_permutations, m). Returns steps(1..m), mean, lo2.5, hi97.5 -- identical convention
    to 39_saturation_by_ancestry.summarize_curve."""
    m = arr.shape[1]
    steps = np.arange(1, m + 1)
    return steps, arr.mean(axis=0), np.percentile(arr, 2.5, axis=0), np.percentile(arr, 97.5, axis=0)


def class_curves_from_thresholds(curves):
    """curves: output of rarefaction_curves_fixed_orders (KIR_THRESHOLDS=[1,2,3,20]). Returns
    {'eq1':arr,'eq2':arr,'gt2':arr,'ge20':arr} via differencing of the >=k "reached" curves --
    e.g. eq1(N) = (#alleles with >=1 carrier at N) - (#alleles with >=2 carriers at N)."""
    return {
        "eq1": curves[1] - curves[2],
        "eq2": curves[2] - curves[3],
        "gt2": curves[3],          # >=3 carriers, includes >=20
        "ge20": curves[20],
    }


def equal_n_slope(steps, mean, n_star, window_frac=0.9):
    """Linear fit of the mean rarefaction curve over [window_frac*n_star, n_star], reported as new
    alleles per next 1,000 people -- 39_saturation_by_ancestry's own "primary, cross-ancestry-
    comparable" equal-N slope, generalized here to any curve (overall or one recurrence class)."""
    if n_star is None or n_star < 2 or n_star > steps[-1]:
        return float("nan")
    lo = max(1, int(round(window_frac * n_star)))
    mask = (steps >= lo) & (steps <= n_star)
    if mask.sum() < 2:
        return float("nan")
    slope = np.polyfit(steps[mask], mean[mask], 1)[0]
    return float(slope * 1000)


def pick_n_star(n_people_by_ancestry, min_people=MIN_PEOPLE_PER_ANCESTRY, exclude=("MID",)):
    """39's own convention: the largest N every well-powered (>=min_people), non-excluded
    ancestry reaches -- i.e. min(count) over those ancestries. Returns None if none qualify."""
    powered = {a: n for a, n in n_people_by_ancestry.items()
               if a not in exclude and n >= min_people}
    return min(powered.values()) if powered else None


# ---------------------------------------------------------------------------
# KIR-side extraction: reuses 43's parse_all()/parse_person() (already returns, per person, per
# haplotype: genes_called + kir_rows with gene/consensus/novelty_tier) -- no GTF parsing here.
# ---------------------------------------------------------------------------
def kir_level_mask(tier, level):
    if level == "all":
        return tier != "undetermined"
    if level == "any_novel":
        return tier not in ("known", "undetermined")
    if level == "protein_novel":
        return tier == "novel_protein"
    raise ValueError(level)


def kir_person_gene_sets(results_by_pid, pids, level, kir):
    """Returns {gene: {person_id: set(consensus)}}, pooling both haplotypes per person (unit =
    person, matching HLA's own sampling unit in 39)."""
    out = defaultdict(lambda: defaultdict(set))
    for pid in pids:
        r = results_by_pid.get(pid)
        if not r:
            continue
        for hap in r["haps"]:
            for row in hap["kir_rows"]:
                if kir_level_mask(row["novelty_tier"], level):
                    out[row["gene"]][pid].add(row["consensus"])
    return out


# ---------------------------------------------------------------------------
# HLA-side extraction: reuses 39_saturation_by_ancestry.build_labeled_calls() (the full
# field_class/seq_class/prot_id/cds_id pipeline) -- no re-derivation of sequence matching.
# ---------------------------------------------------------------------------
def hla_level_mask(calls, level):
    if level == "all":
        return calls["keep_clean"]
    if level == "any_novel":
        return calls["keep_clean"] & (calls["field_class"] != "known")
    if level == "protein_novel":
        return (calls["keep_clean"] & (calls["field_class"] == "f2_protein")
                & (calls["seq_class"] == "novel_protein"))
    raise ValueError(level)


def hla_id_col(level):
    return "prot_id" if level == "protein_novel" else "cds_id"


def hla_person_gene_sets(calls, level, genes_bare):
    """Returns {gene: {person_id: set(allele_id)}}, restricted to `genes_bare` and `level`."""
    id_col = hla_id_col(level)
    mask = hla_level_mask(calls, level) & calls["gene_b"].isin(genes_bare) & calls[id_col].notna()
    sub = calls.loc[mask, ["person_id", "gene_b", id_col]]
    out = defaultdict(lambda: defaultdict(set))
    for pid, gene, aid in zip(sub["person_id"], sub["gene_b"], sub[id_col]):
        out[gene][pid].add(aid)
    return out


# ---------------------------------------------------------------------------
# Joint unrelated set (WS-B "same unrelated people who have both").
# ---------------------------------------------------------------------------
def joint_unrelated_set(kir_pids, hla_pids, relatedness_path, kin_min):
    """Intersect the two species' available-person pools FIRST, then run ONE greedy_unrelated()
    pass over that intersection -- so both species get the exact same unrelated set (running
    greedy_unrelated separately on each species' own full pool would not guarantee this, since the
    kinship graph is restricted to whichever id set is passed in)."""
    kir = m41()
    common = sorted(set(map(str, kir_pids)) & set(map(str, hla_pids)))
    pairs = kir.load_relatedness_pairs(relatedness_path)
    kept, removed = kir.greedy_unrelated(common, pairs, kin_min=kin_min)
    return sorted(kept), len(removed), len(common)


# ---------------------------------------------------------------------------
# Table builders (apply disclosure masking here, at the export boundary).
# ---------------------------------------------------------------------------
def build_recurrence_row(gene, ancestry, level, species, unit_sets):
    counts = carrier_counts(unit_sets)
    cls = classify_recurrence(counts)
    return {"gene": gene, "ancestry": ancestry, "level": level, "species": species,
            "n_people": len(unit_sets), "n_distinct_alleles": suppressed(len(counts)),
            "eq1": suppressed(cls["eq1"]), "eq2": suppressed(cls["eq2"]),
            "gt2": suppressed(cls["gt2"]), "ge20": suppressed(cls["ge20"])}


def build_coverage_chao2_row(gene, ancestry, level, species, unit_sets, extrapolate_2n=True):
    S_obs, m, Q = incidence_freqs(unit_sets)
    n_incidences = sum(k * v for k, v in Q.items())
    coverage = good_turing_coverage(Q, m, n_incidences)
    chao = m04().chao2(S_obs, Q, m)
    se = chao2_se(S_obs, Q, m)
    row = {"gene": gene, "ancestry": ancestry, "level": level, "species": species,
           "n_people": m, "s_obs": S_obs,
           "good_turing_coverage": round(coverage, 4) if coverage == coverage else "",
           "chao2": round(chao, 1), "chao2_se": round(se, 1),
           "chao2_undetected_f0hat": round(chao - S_obs, 1)}
    q1, q2 = Q.get(1, 0), Q.get(2, 0)
    row["q1"] = suppressed(q1) if q1 < SUPPRESS_BELOW else str(q1)
    row["q2"] = suppressed(q2) if q2 < SUPPRESS_BELOW else str(q2)
    if q1 < SUPPRESS_BELOW or q2 < SUPPRESS_BELOW:
        row["q1"] = row["q2"] = ""  # blank both together, never one alone (back-reveal risk)
    if extrapolate_2n:
        row["chao_new_by_2n"] = round(
            m04().chao2_extrapolate(S_obs, Q, m, 2 * m) - S_obs, 1) if m > 0 else ""
    return row


def build_curve_rows(gene, ancestry, level, species, unit_sets_list, orders, n_star, stride):
    curves = rarefaction_curves_fixed_orders(unit_sets_list, orders, KIR_THRESHOLDS)
    class_curves = class_curves_from_thresholds(curves)
    steps, mean_d, lo_d, hi_d = summarize_curve(curves["distinct"])
    means = {k: summarize_curve(v)[1] for k, v in class_curves.items()}
    rows = []
    idxs = list(range(0, len(steps), max(1, stride)))
    if idxs[-1] != len(steps) - 1:
        idxs.append(len(steps) - 1)
    for i in idxs:
        rows.append({
            "gene": gene, "ancestry": ancestry, "level": level, "species": species,
            "n": int(steps[i]), "mean_distinct": round(float(mean_d[i]), 2),
            "lo2_5": round(float(lo_d[i]), 2), "hi97_5": round(float(hi_d[i]), 2),
            "mean_eq1": round(float(means["eq1"][i]), 2),
            "mean_eq2": round(float(means["eq2"][i]), 2),
            "mean_gt2": round(float(means["gt2"][i]), 2),
            "mean_ge20": round(float(means["ge20"][i]), 2),
        })
    slope_rows = []
    if n_star is not None:
        slope_rows.append({
            "gene": gene, "ancestry": ancestry, "level": level, "species": species,
            "n_star": n_star, "curve": "distinct",
            "slope_per_1000": round(equal_n_slope(steps, mean_d, n_star), 2),
            "mean_at_n_star": round(float(np.interp(n_star, steps, mean_d)), 2),
        })
        for cname, arr in means.items():
            slope_rows.append({
                "gene": gene, "ancestry": ancestry, "level": level, "species": species,
                "n_star": n_star, "curve": cname,
                "slope_per_1000": round(equal_n_slope(steps, arr, n_star), 2),
                "mean_at_n_star": round(float(np.interp(n_star, steps, arr)), 2),
            })
    return rows, slope_rows


# ---------------------------------------------------------------------------
# Orchestration: for one gene x level, build unit sets for KIR + HLA at ALL + each ancestry,
# using a SHARED pid ordering/permutation-order set per ancestry (species-agnostic).
# ---------------------------------------------------------------------------
def run_gene_level(gene, level, species_unit_sets_all_ancestry, pids_by_ancestry, orders_by_ancestry,
                    n_star_by_ancestry, stride, extrapolate_2n):
    """species_unit_sets_all_ancestry: {species: {person_id: set(allele_id)}} already restricted
    to `gene`/`level` (ALL people, not yet sliced by ancestry). Returns the four output-row lists
    for this (gene, level) across species x {ALL + each ancestry}."""
    rec_rows, curve_rows, slope_rows, cov_rows = [], [], [], []
    for anc, pids in pids_by_ancestry.items():   # "ALL" is included as a pseudo-ancestry
        orders = orders_by_ancestry[anc]
        for species, full_sets in species_unit_sets_all_ancestry.items():
            unit_sets = {p: full_sets.get(p, set()) for p in pids}
            rec_rows.append(build_recurrence_row(gene, anc, level, species, unit_sets))
            cov_rows.append(build_coverage_chao2_row(gene, anc, level, species, unit_sets,
                                                      extrapolate_2n))
            unit_sets_list = [unit_sets[p] for p in pids]
            n_star = n_star_by_ancestry.get(anc)
            c_rows, s_rows = build_curve_rows(gene, anc, level, species, unit_sets_list, orders,
                                               n_star, stride)
            curve_rows.extend(c_rows)
            slope_rows.extend(s_rows)
    return rec_rows, curve_rows, slope_rows, cov_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kir-outroot", default=os.path.expanduser("~/pipeline_outputs_kir"))
    ap.add_argument("--hla-table", default=os.path.expanduser("~/pipeline_outputs/hla_calls_rich.tsv"))
    ap.add_argument("--cohort-membership", default=os.path.expanduser("~/pipeline_outputs/cohort_membership.tsv"))
    ap.add_argument("--relatedness-table", default=os.path.expanduser(
        "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/relatedness/"
        "samples_relatedness.tsv"))
    ap.add_argument("--refdata", default=os.path.expanduser("~/tools/Immuannot_refdata"))
    ap.add_argument("--out-dir", default=os.path.expanduser("~/s04/results/44"))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--threads", type=int, default=4, help="threads for HLA sequence matching (39's match_sequences)")
    ap.add_argument("--n-permutations", type=int, default=N_PERMUTATIONS_DEFAULT)
    ap.add_argument("--curve-stride", type=int, default=CURVE_STRIDE_DEFAULT)
    ap.add_argument("--seed", type=int, default=44)
    ap.add_argument("--kin-min", type=float, default=0.0442)
    ap.add_argument("--strict-threshold", type=float, default=0.9)
    ap.add_argument("--limit", type=int, default=None, help="smoke test: cap discovered KIR people")
    ap.add_argument("--no-extrapolate-2n", dest="extrapolate_2n", action="store_false")
    ap.add_argument("--skip-relatedness", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    write_status(args.out_dir, f"START outroot={args.kir_outroot} out_dir={args.out_dir}")
    t0 = time.perf_counter()

    kir = m41()
    kir_agg = m43()
    m39mod = m39()

    log("[44] discovering KIR people ...")
    kir_pids = kir_agg.discover_people(args.kir_outroot, args.limit)
    write_status(args.out_dir, f"discovered {len(kir_pids)} KIR people")

    log("[44] parsing KIR GTFs ...")
    kir_results = kir_agg.parse_all(kir_pids, args.kir_outroot, args.workers,
                                     os.path.join(_THIS_DIR, "41_kir_pilot.py"))
    kir_results_by_pid = {r["person_id"]: r for r in kir_results}
    write_status(args.out_dir, f"parsed {len(kir_results)} KIR people")

    log("[44] loading HLA labeled calls (39's pipeline) ...")
    hla_args = argparse.Namespace(
        table1=args.hla_table, limit=args.limit, cohort_membership=args.cohort_membership,
        relatedness_table=args.relatedness_table, kin_min=args.kin_min,
        strict_threshold=args.strict_threshold, skip_relatedness=args.skip_relatedness,
        refdata=args.refdata, outroot=os.path.expanduser("~/pipeline_outputs"), threads=args.threads)
    calls, ref_catalogue_size, mstats, n_removed = m39mod.build_labeled_calls(hla_args)
    hla_pids = sorted(calls["person_id"].unique())
    write_status(args.out_dir, f"HLA calls loaded: {len(hla_pids)} people")

    log("[44] computing joint unrelated set (KIR ∩ HLA, one greedy_unrelated pass) ...")
    unrelated_pids, n_removed_joint, n_common = joint_unrelated_set(
        kir_pids, hla_pids, args.relatedness_table, args.kin_min)
    log(f"[44] joint unrelated set: {len(unrelated_pids)}/{n_common} common people "
        f"({n_removed_joint} relatives dropped)")
    write_status(args.out_dir, f"joint unrelated set: {len(unrelated_pids)}")

    ancestry_of = kir_agg.load_ancestry_map(args.cohort_membership)
    pids_by_ancestry_base = defaultdict(list)
    for p in unrelated_pids:
        anc = ancestry_of.get(str(p), "UNASSIGNED")
        if anc in ANCESTRY_ORDER:
            pids_by_ancestry_base[anc].append(p)
    pids_by_ancestry = {"ALL": unrelated_pids, **pids_by_ancestry_base}
    n_people_by_ancestry = {a: len(ps) for a, ps in pids_by_ancestry_base.items()}
    n_star = pick_n_star(n_people_by_ancestry)
    n_star_by_ancestry = {a: (n_star if a != "ALL" else len(unrelated_pids)) for a in pids_by_ancestry}
    log(f"[44] N* (equal-N reference, excl. MID) = {n_star}")

    orders_by_ancestry = {a: make_orders(len(ps), args.n_permutations, args.seed + hash(a) % 10_000)
                          for a, ps in pids_by_ancestry.items()}

    kir_genes = kir.KIR_GENES
    hla_genes = m39mod.CLASSICAL_GENES_BARE

    all_rec, all_curve, all_slope, all_cov = [], [], [], []
    for level in LEVELS:
        log(f"[44] level={level}: building KIR per-gene sets ...")
        kir_sets_all = kir_person_gene_sets(kir_results_by_pid, unrelated_pids, level, kir)
        log(f"[44] level={level}: building HLA per-gene sets ...")
        hla_sets_all = hla_person_gene_sets(calls, level, hla_genes)

        for gene in kir_genes:
            species_sets = {"kir": kir_sets_all.get(gene, {})}
            rec, curve, slope, cov = run_gene_level(
                gene, level, species_sets, pids_by_ancestry, orders_by_ancestry,
                n_star_by_ancestry, args.curve_stride, args.extrapolate_2n)
            all_rec += rec; all_curve += curve; all_slope += slope; all_cov += cov
        for gene in hla_genes:
            species_sets = {"hla": hla_sets_all.get(gene, {})}
            rec, curve, slope, cov = run_gene_level(
                gene, level, species_sets, pids_by_ancestry, orders_by_ancestry,
                n_star_by_ancestry, args.curve_stride, args.extrapolate_2n)
            all_rec += rec; all_curve += curve; all_slope += slope; all_cov += cov
        write_status(args.out_dir, f"level={level} done ({len(all_rec)} recurrence rows so far)")

    pd.DataFrame(all_rec).to_csv(os.path.join(args.out_dir, "recurrence_classes.tsv"),
                                  sep="\t", index=False)
    pd.DataFrame(all_curve).to_csv(os.path.join(args.out_dir, "saturation_curves.tsv"),
                                    sep="\t", index=False)
    pd.DataFrame(all_slope).to_csv(os.path.join(args.out_dir, "equal_n_slope.tsv"),
                                    sep="\t", index=False)
    pd.DataFrame(all_cov).to_csv(os.path.join(args.out_dir, "coverage_chao2.tsv"),
                                  sep="\t", index=False)
    write_status(args.out_dir, f"DONE in {time.perf_counter()-t0:.0f}s -- 4 tables written")
    log(f"[44] done in {time.perf_counter()-t0:.0f}s")


if __name__ == "__main__":
    main()
