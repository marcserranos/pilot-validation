#!/usr/bin/env python3
"""Allele-space saturation analysis -- rarefaction/extrapolation curves and nonparametric richness
estimators (Chao2, ACE/ICE, jackknife), per gene per ancestry, for "all alleles", "novel alleles
only", and "known alleles only" -- per scripts/hla_popgen/research/NOVEL_LIT.md section 3.4's
concrete recommendation. This is the paper-worthy headline analysis: "% of allele space discovered"
per gene per ancestry, plus a projection of how many more alleles a 2x/5x/10x larger cohort would
find.

## Why Chao2 (incidence-based), not Chao1 (abundance-based) -- NOVEL_LIT.md section 3.4

The sampling unit here is the HAPLOTYPE (not the person, and not raw allele copies): what we
actually observe per allele is presence/absence across independent haplotypes, which is exactly the
"incidence" data structure Chao2/ACE/jackknife are built for, and is less sensitive than an
abundance-based estimator to a single person happening to carry an allele on both haplotypes
(homozygosity), which would otherwise inflate a naive Chao1 count. N ~= 2 x n_people (haplotypes),
not n_people.

## The balancing-selection problem -- addressed head-on, not buried in a caveat (NOVEL_LIT.md
## section 4, required by the task spec)

HLA is one of the best-documented violations of the neutral infinite-alleles model that the
classical Ewens sampling formula / Watterson's theta assume (Ewens-Watterson homozygosity tests
reject neutrality in the majority of tested population-locus combinations -- NOVEL_LIT.md section
4, citing the Genetics 2022 joint selection/mutation-rate paper and the classical Hedrick/Thomson-era
result). Balancing selection keeps more alleles at more even frequencies than neutral drift-mutation
predicts, which would make a neutral-model-based unseen-allele estimate (Ewens/Watterson) SYSTEMATICALLY
BIASED, most likely underestimating undiscovered diversity. Chao-family estimators are nonparametric:
they use only the counts of rare classes (uniques/duplicates across sampling units) and a general
inequality argument (Chao 1984; Chao & Colwell's Hill-number/iNEXT extension), with NO assumption
about the shape of the underlying frequency spectrum -- this is exactly why they remain valid under
selection-driven distributional shifts that break the neutral model, and exactly why this script uses
Chao2 as its headline estimator and computes the Ewens/Watterson neutral expectation ONLY as a
secondary, explicitly-labeled comparison to quantify the selection effect as a RESULT (how much the
neutral model would have underestimated undiscovered diversity), never as the headline number. See
`ewens_watterson_expected_k()` / the "neutral model comparison" report section.

Residual risk even with Chao2 (documented, not glossed over): Chao-family estimators still assume
each existing allele has a roughly homogeneous per-sampling-unit detection probability; ACE/ICE
partially corrects for heterogeneity in that probability among rare classes, which is why it's run
alongside Chao2 as a robustness check rather than treated as redundant. No existing richness
estimator was designed with balancing selection specifically in mind -- Chao2's numbers here should
be read as the best available nonparametric LOWER BOUND on undiscovered diversity, not a mechanistic
population-genetic model fit (NOVEL_LIT.md section 4, final bullet).

## No off-the-shelf Python iNEXT equivalent (NOVEL_LIT.md section 3.1)

The formulas below (Chao2, ACE/ICE, jackknife 1st/2nd order, Good-Turing-style extrapolation) are
implemented directly from their closed forms (Chao 1987; Chao & Lee 1992; Colwell et al. 2012's
unified rarefaction/extrapolation formula) rather than via `rpy2`+R's `iNEXT`, per NOVEL_LIT.md's own
finding that no mature Python equivalent exists. `tests/test_novel.py` unit-tests each formula
against a hand-computed toy example.

Usage:
    python3 scripts/hla_popgen/04_allele_saturation.py --outroot /tmp/hla_fixtures \\
        --table1 /tmp/hla_fixtures/hla_calls_rich.sample.tsv \\
        --table3 /tmp/hla_fixtures/reports/novel_alleles.sample.tsv \\
        --cohort-membership /tmp/hla_fixtures/cohort_membership.sample.tsv \\
        --out-dir /tmp/hla_fixtures/reports --sample --n-bootstrap 200 --n-permutations 200

    # Real run on the VM:
    python3 scripts/hla_popgen/04_allele_saturation.py
"""
import argparse
import importlib.util
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

DEFAULT_OUTROOT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_REPORTS_DIR_NAME = os.path.join("reports", "hla_popgen")
ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
CLASSICAL_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
                   "HLA-DRB1"]
CATEGORIES = ["all", "known", "novel"]
EXTRAPOLATION_FACTORS = [2, 5, 10]


def _load_sibling_module(filename, modname):
    """Load a digit-prefixed sibling module (can't `import 03_novel_alleles` directly)."""
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, filename)
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# Richness estimator formulas -- pure functions of incidence-frequency counts, unit-tested in
# tests/test_novel.py against hand-computed examples.
# ---------------------------------------------------------------------------
def incidence_freqs(unit_sets):
    """unit_sets: list of sets, one per sampling unit (haplotype), each holding the allele-identity
    strings observed in that unit. Returns (S_obs, m, Q: Counter mapping occurrence-count -> number
    of species with that occurrence count across units)."""
    m = len(unit_sets)
    occ = Counter()
    for s in unit_sets:
        for allele in s:
            occ[allele] += 1
    S_obs = len(occ)
    Q = Counter()
    for allele, n_units in occ.items():
        Q[n_units] += 1
    return S_obs, m, Q


def chao2(S_obs, Q, m):
    """Incidence-based Chao2 (Chao 1987). Q1 = uniques (seen in exactly 1 unit), Q2 = duplicates.
    Bias-corrected form used when Q2 == 0 (Chao & Colwell), to avoid a divide-by-zero blowup."""
    Q1, Q2 = Q.get(1, 0), Q.get(2, 0)
    if m <= 1:
        return float(S_obs)
    factor = (m - 1) / m
    if Q2 > 0:
        return S_obs + factor * (Q1 ** 2) / (2 * Q2)
    return S_obs + factor * Q1 * (Q1 - 1) / 2


def ace_incidence(S_obs, Q, m, kappa=10):
    """Incidence-based Coverage-corrected estimator (ICE; the incidence analogue of Chao & Lee
    1992's ACE, sharing the same formula with occurrence-counts across sampling units in place of
    abundance counts -- see NOVEL_LIT.md's ACE/ICE row). `kappa`=10 is the standard rare/abundant
    species cutoff from Chao & Lee 1992."""
    S_rare = sum(n for i, n in Q.items() if i <= kappa)
    S_abund = S_obs - S_rare
    N_rare = sum(i * n for i, n in Q.items() if i <= kappa)
    Q1 = Q.get(1, 0)
    if S_rare == 0 or N_rare == 0:
        return float(S_obs)
    C_ace = 1 - Q1 / N_rare
    if C_ace <= 0:
        # Coverage estimate degenerates (almost everything is a unique) -- fall back to Chao2,
        # which handles this regime more gracefully, rather than dividing by ~0.
        return chao2(S_obs, Q, m)
    if N_rare > 1:
        sum_term = sum(i * (i - 1) * n for i, n in Q.items() if i <= kappa)
        gamma2 = max((S_rare / C_ace) * sum_term / (N_rare * (N_rare - 1)) - 1, 0.0)
    else:
        gamma2 = 0.0
    return S_abund + S_rare / C_ace + (Q1 / C_ace) * gamma2


def jackknife1_incidence(S_obs, Q, m):
    """First-order incidence jackknife (Burnham & Overton 1978)."""
    if m <= 0:
        return float(S_obs)
    Q1 = Q.get(1, 0)
    return S_obs + Q1 * (m - 1) / m


def jackknife2_incidence(S_obs, Q, m):
    """Second-order incidence jackknife -- more conservative, less sensitive to the Q1/Q2 ratio
    alone (NOVEL_LIT.md's jackknife row)."""
    if m <= 1:
        return float(S_obs)
    Q1, Q2 = Q.get(1, 0), Q.get(2, 0)
    return S_obs + Q1 * (2 * m - 3) / m - Q2 * ((m - 2) ** 2) / (m * (m - 1))


def chao2_extrapolate(S_obs, Q, m, t):
    """Colwell et al. 2012's unified incidence extrapolation formula: expected richness at t total
    sampling units (t >= m), using the Chao2-estimated undetected-species count f0_hat. Falls back
    to S_obs (no extrapolation signal) when there's no unseen-mass evidence (Q1==0 or f0_hat<=0)."""
    if t <= m:
        return float(S_obs)
    Q1 = Q.get(1, 0)
    f0_hat = chao2(S_obs, Q, m) - S_obs
    if f0_hat <= 0 or Q1 == 0:
        return float(S_obs)
    return S_obs + f0_hat * (1 - (1 - Q1 / (m * f0_hat + Q1)) ** (t - m))


def ewens_watterson_expected_k(theta, n):
    """E[number of distinct alleles] in a neutral-infinite-alleles sample of size n, given theta
    (Ewens 1972 sampling formula's own expectation, closed-form -- Watterson's theta estimator
    inverts this). SECONDARY comparison only -- see module docstring's balancing-selection section;
    never the headline estimator in this script."""
    n = int(n)
    if n <= 0:
        return 0.0
    i = np.arange(n)
    return float(np.sum(theta / (theta + i)))


def solve_theta_for_k(k_obs, n, lo=1e-6, hi=1e6, iters=100):
    """Invert ewens_watterson_expected_k() for theta given an observed distinct-allele count k_obs
    at sample size n. E[K] is strictly increasing in theta, so bisection is exact and stable."""
    if k_obs <= 0 or n <= 0:
        return 0.0
    f = lambda th: ewens_watterson_expected_k(th, n) - k_obs
    if f(hi) < 0:
        hi *= 100  # k_obs implausibly close to n; widen the bracket rather than fail silently
    for _ in range(iters):
        mid = (lo + hi) / 2
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


# ---------------------------------------------------------------------------
# Bootstrap CI (resample sampling units with replacement -- standard nonparametric approach used by
# iNEXT itself for these estimators, per NOVEL_LIT.md's own recommendation to report CIs, not just
# point estimates, precisely because selection-driven distributional differences should widen them).
# ---------------------------------------------------------------------------
def bootstrap_ci(unit_sets, estimator_fn, n_bootstrap=500, seed=0, alpha=0.05):
    rng = np.random.default_rng(seed)
    m = len(unit_sets)
    if m == 0:
        return (np.nan, np.nan)
    idx = np.arange(m)
    vals = []
    for _ in range(n_bootstrap):
        sample_idx = rng.choice(idx, size=m, replace=True)
        boot_units = [unit_sets[i] for i in sample_idx]
        S_obs, mb, Q = incidence_freqs(boot_units)
        vals.append(estimator_fn(S_obs, Q, mb))
    lo, hi = np.percentile(vals, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi)


# ---------------------------------------------------------------------------
# Rarefaction curve (permutation-averaged accumulation, per NOVEL_LIT.md 3.1's own recommended
# method when a closed-form rarefaction isn't implemented).
# ---------------------------------------------------------------------------
def rarefaction_curve(unit_sets, n_permutations=200, seed=0):
    """Returns (steps: 1..m, mean distinct-count per step, lo95, hi95) averaged over random orderings
    of the sampling units."""
    m = len(unit_sets)
    if m == 0:
        return np.array([]), np.array([]), np.array([]), np.array([])
    rng = np.random.default_rng(seed)
    curves = np.zeros((n_permutations, m))
    idx = np.arange(m)
    for p in range(n_permutations):
        order = rng.permutation(idx)
        seen = set()
        for step, i in enumerate(order):
            seen |= unit_sets[i]
            curves[p, step] = len(seen)
    steps = np.arange(1, m + 1)
    mean = curves.mean(axis=0)
    lo = np.percentile(curves, 2.5, axis=0)
    hi = np.percentile(curves, 97.5, axis=0)
    return steps, mean, lo, hi


# ---------------------------------------------------------------------------
# Building the per-(gene, ancestry, category) incidence data from Table 1 + Table 3 + Table 4
# ---------------------------------------------------------------------------
def build_allele_identity_table(table1, matched_rows):
    """Per (person_id, hap, gene) row of Table 1, derive an allele-identity string suitable for
    incidence analysis: the `consensus` string for known (non-novel) rows, `novel_id` for novel rows
    that were unambiguously matched to a cds.fa.gz sequence (03_novel_alleles.py's matching logic,
    re-used here rather than re-derived), and NaN (excluded) for `undetermined` calls or novel rows
    that could not be resolved (ambiguous copy_index / missing cds.fa.gz record) -- SCHEMA.md
    discipline: never guess an identity for a row we can't actually resolve.
    """
    t = table1.copy()
    t["is_novel_bool"] = t["is_novel"].astype(str).str.lower().isin(["true", "1"])

    novel_id_lookup = {}
    for r in matched_rows:
        novel_id = f"{r['gene']}_nov_{r['cds_seq_sha1'][:8]}"
        novel_id_lookup[(r["person_id"], r["hap"], r["gene"])] = novel_id

    def identity(row):
        if row["is_novel_bool"]:
            return novel_id_lookup.get((row["person_id"], row["hap"], row["gene"]))
        c = row.get("consensus")
        if not isinstance(c, str) or c == "undetermined":
            return None
        return c

    t["allele_id"] = t.apply(identity, axis=1)
    return t


def build_unit_sets(table1_ident, ancestry_by_person, gene, ancestry, category):
    """Returns list of sets (one per haplotype sampling unit) of allele-identity strings observed
    for `gene`, restricted to `ancestry` (or None for pooled) and `category` in {"all","known",
    "novel"}. Haplotypes with zero resolved alleles for this gene are legitimately absent (not a
    zero-padded unit) -- they simply don't appear as sampling units for this gene, matching the
    real-world fact that gene detection can fail per-haplotype (SCHEMA.md: absence is absence, never
    NA-imputed)."""
    sub = table1_ident[table1_ident["gene"] == gene].copy()
    sub = sub[sub["allele_id"].notna()]
    if category == "novel":
        sub = sub[sub["is_novel_bool"]]
    elif category == "known":
        sub = sub[~sub["is_novel_bool"]]
    if ancestry is not None:
        sub["ancestry"] = sub["person_id"].map(ancestry_by_person)
        sub = sub[sub["ancestry"] == ancestry]
    by_unit = defaultdict(set)
    for _, r in sub.iterrows():
        by_unit[(r["person_id"], r["hap"])].add(r["allele_id"])
    return list(by_unit.values())


# ---------------------------------------------------------------------------
# Main analysis loop
# ---------------------------------------------------------------------------
def analyze(table1_ident, ancestry_by_person, genes, n_bootstrap, n_permutations):
    """Returns (richness_df, extrapolation_df, curves: {(gene,category,ancestry): (steps,mean,lo,hi)},
    neutral_df)."""
    richness_rows = []
    extrap_rows = []
    neutral_rows = []
    curves = {}

    groups = ANCESTRY_ORDER + [None]  # None = pooled
    for gene in genes:
        for category in CATEGORIES:
            for ancestry in groups:
                unit_sets = build_unit_sets(table1_ident, ancestry_by_person, gene, ancestry,
                                             category)
                m = len(unit_sets)
                anc_label = ancestry or "POOLED"
                if m == 0:
                    continue
                S_obs, m, Q = incidence_freqs(unit_sets)

                c2 = chao2(S_obs, Q, m)
                c2_lo, c2_hi = bootstrap_ci(unit_sets, chao2, n_bootstrap=n_bootstrap, seed=hash(
                    (gene, category, anc_label)) % (2 ** 31))
                ace = ace_incidence(S_obs, Q, m)
                jk1 = jackknife1_incidence(S_obs, Q, m)
                jk2 = jackknife2_incidence(S_obs, Q, m)

                # BUG FIX (found once real novel-allele recurrence entered the fixtures): a
                # percentile bootstrap CI is computed over BOOTSTRAP REPLICATES and is not
                # mathematically guaranteed to bracket the point estimate computed on the full
                # sample -- this bites hardest exactly when Q2_duplicates==0, where the
                # bias-corrected Chao2 form (S_obs + (m-1)/m * Q1(Q1-1)/2) blows up quadratically on
                # the full sample while most bootstrap resamples don't reproduce Q2==0 that cleanly
                # (a resample can easily pick up a duplicate by chance). Observed on
                # /tmp/hla_fx2 fixtures: HLA-B/novel/AFR had chao2=11.85 but a raw percentile CI of
                # [3.00, 10.85] -- the point estimate sat OUTSIDE its own interval. Rather than
                # silently emit an incoherent interval, widen it to guarantee it brackets the point
                # estimate and flag the row so a reviewer can see this happened rather than reading
                # a clean-looking number.
                chao2_ci_widened = bool(
                    c2_lo == c2_lo and c2_hi == c2_hi and not (c2_lo <= c2 <= c2_hi))
                if chao2_ci_widened:
                    c2_lo, c2_hi = min(c2_lo, c2), max(c2_hi, c2)

                # BUG FIX: "% of allele space discovered" cannot exceed 100% by definition --
                # Chao2 >= S_obs always, by construction of every formula above, so the POINT
                # estimate of pct_discovered was already <=100 mathematically, but the CI bounds
                # were not: pct = S_obs / chao2 is a DECREASING function of chao2, so pct's lower
                # bound must come from chao2's UPPER bound and pct's upper bound from chao2's LOWER
                # bound (this inversion direction was already correct in the code -- the apparent
                # "crossed" cases reported, e.g. HLA-B/novel/AFR pct=50.6 with ci_lo=55.3, were a
                # symptom of the point estimate sitting outside its own un-widened chao2 CI, now
                # fixed above). Separately, a bootstrap replicate can have fewer distinct alleles
                # than the true full-sample S_obs (it's a resample, not a superset), which can push
                # a chao2 CI bound below the true S_obs and make the naive pct bound exceed 100%
                # even with a correctly bracketed chao2 CI -- clamping is what actually enforces the
                # definitional bound. `np.clip` is monotonic non-decreasing, so applying it to all
                # three of (pct_lo, pct, pct_hi) cannot invert the lo<=point<=hi ordering already
                # established by the chao2 bracket fix above.
                def _pct(s_obs, c2_val):
                    if not (c2_val and c2_val > 0):
                        return np.nan
                    return float(np.clip(100 * s_obs / c2_val, 0.0, 100.0))

                pct_discovered = _pct(S_obs, c2)
                pct_lo = _pct(S_obs, c2_hi)
                pct_hi = _pct(S_obs, c2_lo)

                richness_rows.append({
                    "gene": gene, "category": category, "ancestry": anc_label,
                    "n_haplotypes": m, "S_obs": S_obs, "Q1_uniques": Q.get(1, 0),
                    "Q2_duplicates": Q.get(2, 0),
                    "chao2": round(c2, 2), "chao2_ci_lo": round(c2_lo, 2) if c2_lo == c2_lo else None,
                    "chao2_ci_hi": round(c2_hi, 2) if c2_hi == c2_hi else None,
                    "ace": round(ace, 2), "jackknife1": round(jk1, 2), "jackknife2": round(jk2, 2),
                    "pct_discovered": round(pct_discovered, 1),
                    "pct_discovered_ci_lo": round(pct_lo, 1) if pct_lo == pct_lo else None,
                    "pct_discovered_ci_hi": round(pct_hi, 1) if pct_hi == pct_hi else None,
                    # Q2==0 (no doubletons) is Chao2's known degenerate regime: the bias-corrected
                    # form (S_obs + (m-1)/m * Q1(Q1-1)/2) grows QUADRATICALLY in Q1 with nothing to
                    # anchor it, and can be wildly unstable/implausible -- flag rather than silently
                    # report a headline number iNEXT itself would caveat. Real data with genuine
                    # allele recurrence should rarely hit this; synthetic per-call-random fixtures
                    # will hit it constantly for the "novel" category by construction (see
                    # tests/test_novel.py and the run-time report to Marc for why).
                    "low_confidence_chao2": (Q.get(2, 0) == 0 and Q.get(1, 0) > 1),
                    # Set True whenever the raw bootstrap percentile CI failed to bracket the point
                    # estimate and had to be widened (see fix above) -- a reviewer-visible signal
                    # distinct from (though correlated with) low_confidence_chao2.
                    "chao2_ci_widened": chao2_ci_widened,
                })

                for factor in EXTRAPOLATION_FACTORS:
                    t_target = int(round(m * factor))
                    est = chao2_extrapolate(S_obs, Q, m, t_target)
                    extrap_rows.append({
                        "gene": gene, "category": category, "ancestry": anc_label,
                        "n_haplotypes_observed": m, "factor": factor,
                        "n_haplotypes_target": t_target, "S_obs": S_obs,
                        "extrapolated_richness": round(est, 2),
                        "additional_alleles_projected": round(est - S_obs, 2),
                    })

                if category == "all" and gene in CLASSICAL_GENES:
                    theta = solve_theta_for_k(S_obs, m)
                    ewens_at_m = ewens_watterson_expected_k(theta, m)
                    ewens_at_2m = ewens_watterson_expected_k(theta, 2 * m)
                    neutral_rows.append({
                        "gene": gene, "ancestry": anc_label, "n_haplotypes": m, "S_obs": S_obs,
                        "watterson_theta": round(theta, 3),
                        "neutral_expected_at_2x": round(ewens_at_2m, 2),
                        "neutral_additional_at_2x": round(ewens_at_2m - S_obs, 2),
                        "chao2_estimated_total": round(c2, 2),
                        "chao2_additional_unseen": round(max(c2 - S_obs, 0), 2),
                    })

                if category != "known" or gene in CLASSICAL_GENES:
                    # Rarefaction curves: keep this bounded (all genes x categories x pooled+6
                    # ancestries would be hundreds of curves) -- compute for every
                    # (gene, category, ancestry) pair but only PLOT classical genes downstream.
                    curves[(gene, category, anc_label)] = rarefaction_curve(
                        unit_sets, n_permutations=n_permutations,
                        seed=hash((gene, category, anc_label)) % (2 ** 31))

    return (pd.DataFrame(richness_rows), pd.DataFrame(extrap_rows), curves,
            pd.DataFrame(neutral_rows))


# ---------------------------------------------------------------------------
# Plotting -- headline figure: novel-allele discovery curves per classical gene, by ancestry.
# ---------------------------------------------------------------------------
ANCESTRY_COLOR = {
    "AFR": "#4C72B0", "AMR": "#DD8452", "EAS": "#55A868",
    "EUR": "#C44E52", "MID": "#8172B2", "SAS": "#937860", "POOLED": "#333333",
}


def plot_discovery_curves(curves, genes, category, out_path, title_suffix):
    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    for ax, gene in zip(axes.flat, genes):
        any_line = False
        for anc in ANCESTRY_ORDER + ["POOLED"]:
            key = (gene, category, anc)
            if key not in curves:
                continue
            steps, mean, lo, hi = curves[key]
            if len(steps) == 0:
                continue
            any_line = True
            color = ANCESTRY_COLOR[anc]
            lw = 2.2 if anc == "POOLED" else 1.3
            ax.plot(steps, mean, color=color, lw=lw, label=anc if gene == genes[0] else None)
            ax.fill_between(steps, lo, hi, color=color, alpha=0.12)
        ax.set_title(gene, fontsize=11)
        ax.set_xlabel("Haplotypes sampled", fontsize=8)
        ax.set_ylabel("Distinct alleles", fontsize=8)
        # "Distinct alleles" is a discrete count -- force integer y-ticks so a small-range panel
        # doesn't render fractional ticks (1.0, 1.5, 2.0, ...). Same bug class this repo already
        # hit once before (context/STATUS.md: "non-integer histogram bins on a discrete variable"
        # in the previous analysis round's fixture-caught bugs).
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.spines[["top", "right"]].set_visible(False)
        if not any_line:
            ax.axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=7, bbox_to_anchor=(0.5, 1.03),
               frameon=False)
    fig.suptitle(f"Rarefaction / accumulation curves -- {title_suffix}", y=1.07, fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def load_table1(path):
    if not os.path.exists(path):
        sys.exit(f"FATAL: Table 1 not found at {path!r}. Run 01_extract_rich.py first.")
    return pd.read_csv(path, sep="\t", dtype={"person_id": str})


def load_ancestry(path):
    if not os.path.exists(path):
        sys.exit(f"FATAL: Table 4 not found at {path!r}. Run 02_build_cohorts.py first.")
    df = pd.read_csv(path, sep="\t", dtype={"person_id": str})
    df["ancestry_pred"] = df["ancestry_pred"].astype(str).str.upper()
    df.loc[~df["ancestry_pred"].isin(ANCESTRY_ORDER), "ancestry_pred"] = None
    return dict(zip(df["person_id"], df["ancestry_pred"]))


def write_report(md_path, richness_df, extrap_df, neutral_df, fig_paths, is_sample=False):
    md = []
    md.append("# Allele-space saturation analysis (`04_allele_saturation.py`)\n")
    md.append(
        "Methodology: scripts/hla_popgen/research/NOVEL_LIT.md section 3.4. Headline estimator: "
        "**Chao2** (incidence-based, haplotypes as sampling units), with ACE and 1st/2nd-order "
        "jackknife as nonparametric robustness checks. HLA violates the neutral infinite-alleles "
        "model that Ewens/Watterson-style estimators assume (balancing selection -- see module "
        "docstring and the neutral-model comparison section below); Chao-family estimators make no "
        "assumption about the shape of the allele-frequency spectrum and are used as the headline "
        "number for exactly that reason.\n"
    )

    md.append("\n## Headline: % of allele space discovered, novel alleles only, classical genes, "
               "pooled cohort\n")
    headline = richness_df[(richness_df["category"] == "novel") &
                            (richness_df["gene"].isin(CLASSICAL_GENES)) &
                            (richness_df["ancestry"] == "POOLED")]
    md.append("| gene | n_haplotypes | S_obs | Chao2 (95% CI) | % discovered (95% CI) | |")
    md.append("|---|---|---|---|---|---|")
    any_flagged = False
    for _, r in headline.sort_values("gene").iterrows():
        flags = []
        if r.get("low_confidence_chao2"):
            flags.append("Q2=0")
        if r.get("chao2_ci_widened"):
            flags.append("CI widened")
        flag = f" LOW-CONFIDENCE ({', '.join(flags)})" if flags else ""
        any_flagged = any_flagged or bool(flags)
        md.append(f"| {r['gene']} | {r['n_haplotypes']} | {r['S_obs']} | "
                   f"{r['chao2']} ({r['chao2_ci_lo']}-{r['chao2_ci_hi']}) | "
                   f"{r['pct_discovered']}% ({r['pct_discovered_ci_lo']}-"
                   f"{r['pct_discovered_ci_hi']}%) |{flag} |")
    if any_flagged:
        md.append(
            "\n**LOW-CONFIDENCE rows** carry one or both of: `Q2=0` (zero doubletons -- Chao2's "
            "bias-corrected form is degenerate in this regime, growing quadratically in the "
            "singleton count with nothing to anchor it; ACE/jackknife columns in the full TSV are "
            "more conservative fallbacks) and `CI widened` (the raw bootstrap percentile interval "
            "did not bracket the point estimate on the full sample and was widened to guarantee it "
            "does -- see `04_allele_saturation.py`'s `analyze()` for why this happens and is "
            "expected precisely in the Q2=0 regime). Neither should be read as a reliable point "
            "estimate on its own. This mostly hits small, low-haplotype-count "
            "(gene, ancestry) strata even on the current recurrence-bearing fixtures -- the "
            "pooled/`known`-category rows, which have larger N and a finite underlying allele pool, "
            "saturate cleanly and rarely carry this flag.\n"
        )

    md.append("\n## Ancestry-stratified % discovered, novel alleles only, classical genes\n")
    strat = richness_df[(richness_df["category"] == "novel") &
                         (richness_df["gene"].isin(CLASSICAL_GENES)) &
                         (richness_df["ancestry"] != "POOLED")]
    if len(strat):
        pivot = strat.pivot_table(index="gene", columns="ancestry", values="pct_discovered")
        pivot = pivot.reindex(columns=[a for a in ANCESTRY_ORDER if a in pivot.columns])
        md.append("| gene | " + " | ".join(pivot.columns) + " |")
        md.append("|---|" + "---|" * len(pivot.columns))
        for gene, row in pivot.iterrows():
            cells = [f"{v:.1f}%" if v == v else "n/a" for v in row]
            md.append(f"| {gene} | " + " | ".join(cells) + " |")
        md.append(
            "\nIPD-IMGT/HLA's European bias predicts EUR should show the highest % discovered "
            "(closest to saturation) and AFR the lowest -- this table is where that claim either "
            "holds up or doesn't against the data actually loaded.\n"
        )

    md.append("\n## Extrapolation: projected additional novel alleles at 2x/5x/10x cohort size "
               "(pooled, classical genes)\n")
    ex = extrap_df[(extrap_df["category"] == "novel") & (extrap_df["gene"].isin(CLASSICAL_GENES)) &
                   (extrap_df["ancestry"] == "POOLED")]
    if len(ex):
        piv = ex.pivot_table(index="gene", columns="factor", values="additional_alleles_projected")
        md.append("| gene | " + " | ".join(f"{f}x" for f in piv.columns) + " |")
        md.append("|---|" + "---|" * len(piv.columns))
        for gene, row in piv.iterrows():
            md.append(f"| {gene} | " + " | ".join(f"{v:.1f}" for v in row) + " |")

    md.append(
        "\n## Neutral (Ewens/Watterson) comparison -- SECONDARY, quantifying the selection effect\n"
        "Watterson's theta fit to the observed distinct-allele count, then the neutral model's own "
        "expectation of additional alleles at 2x the current haplotype count -- compared against "
        "Chao2's nonparametric estimate of currently-unseen alleles. **Do not read the neutral "
        "column as a discovery projection** -- it is here specifically to show how much a "
        "neutral-model fit would have missed, per NOVEL_LIT.md section 4's instruction to turn the "
        "assumption violation into a result.\n"
    )
    if len(neutral_df):
        nd = neutral_df[neutral_df["ancestry"] == "POOLED"].sort_values("gene")
        md.append("| gene | S_obs | Watterson theta | neutral additional @2x | "
                   "Chao2 unseen (now) |")
        md.append("|---|---|---|---|---|")
        for _, r in nd.iterrows():
            md.append(f"| {r['gene']} | {r['S_obs']} | {r['watterson_theta']} | "
                       f"{r['neutral_additional_at_2x']} | {r['chao2_estimated_total'] - r['S_obs']:.1f} |")

    if is_sample:
        md.append(
            "\n## Note on the discovery-curve figures below (synthetic fixture run)\n"
            "**This run is against synthetic test fixtures (`--sample`), not real data.** The "
            "rarefaction/accumulation curves below saturate hard and quickly -- that is expected "
            "and is NOT a saturation finding: `make_fixtures.py` draws novel alleles from a small, "
            "FINITE per-gene synthetic pool by construction (needed so recurrence across people is "
            "testable at all), so the curves necessarily flatten once that small pool is exhausted. "
            "Real IPD-scale HLA data has an allele space many orders of magnitude larger and is not "
            "expected to saturate this way, especially for non-European ancestries against the "
            "European-biased IPD-IMGT/HLA reference -- do not read a fixture run's curve shape or "
            "'% discovered' numbers as a real saturation result.\n"
        )

    for p in fig_paths:
        md.append(f"\nFigure: `{p}`")
    md.append("")

    text = "\n".join(md)
    with open(md_path, "w") as f:
        f.write(text)
    return text


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT,
                    help="Root holding cds.fa.gz (needed to re-derive novel_id per haplotype -- "
                         "see module docstring on why this re-uses 03's matching logic).")
    ap.add_argument("--table1", default=None, help="Path to hla_calls_rich.tsv. Default: "
                                                     "<outroot>/hla_calls_rich.tsv")
    ap.add_argument("--table3", default=None,
                    help="Path to novel_alleles.tsv -- used only for a sanity cross-check against "
                         "the count of clusters this script independently re-derives.")
    ap.add_argument("--cohort-membership", default=None,
                    help="Path to cohort_membership.tsv (Table 4). Default: <outroot>/"
                         "cohort_membership.tsv")
    ap.add_argument("--out-dir", default=None,
                    help="Where to write TSVs/figures/report. Default: <repo_root>/reports/"
                         "hla_popgen")
    ap.add_argument("--sample", action="store_true")
    ap.add_argument("--genes", default=None,
                    help="Comma-separated gene list to restrict analysis to (default: every gene "
                         "present in Table 1 -- can be slow on the full real cohort; classical "
                         "genes are always included in the headline report sections).")
    ap.add_argument("--n-bootstrap", type=int, default=500)
    ap.add_argument("--n-permutations", type=int, default=200)
    args = ap.parse_args()

    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    table1_path = args.table1 or os.path.join(args.outroot, "hla_calls_rich.tsv")
    cohort_path = args.cohort_membership or os.path.join(args.outroot, "cohort_membership.tsv")
    out_dir = args.out_dir or os.path.join(repo_root, DEFAULT_REPORTS_DIR_NAME)
    os.makedirs(out_dir, exist_ok=True)
    suffix = ".sample" if args.sample else ""

    novel_mod = _load_sibling_module("03_novel_alleles.py", "hla_popgen_novel_alleles")

    print(f"Loading Table 1 from {table1_path!r} ...", file=sys.stderr)
    table1 = load_table1(table1_path)
    print(f"Loading Table 4 (ancestry) from {cohort_path!r} ...", file=sys.stderr)
    ancestry_by_person = load_ancestry(cohort_path)

    print("Re-deriving novel-allele haplotype matches (03_novel_alleles.py's matching logic) ...",
          file=sys.stderr)
    matched_rows, _seqs, match_stats = novel_mod.match_novel_rows(table1, args.outroot)
    print(f"  {match_stats}", file=sys.stderr)
    if args.table3 and os.path.exists(args.table3):
        t3 = pd.read_csv(args.table3, sep="\t")
        n_derived_clusters = len({(r["gene"], r["cds_seq_sha1"]) for r in matched_rows})
        if n_derived_clusters != len(t3):
            print(f"  WARNING: re-derived {n_derived_clusters} clusters vs {len(t3)} rows in "
                  f"{args.table3!r} -- check both scripts were run against the same Table 1 "
                  f"snapshot.", file=sys.stderr)

    table1_ident = build_allele_identity_table(table1, matched_rows)

    genes = (sorted(g.strip() for g in args.genes.split(",")) if args.genes
             else sorted(table1["gene"].unique()))
    # Classical genes must always be present in the headline sections even under --genes.
    for g in CLASSICAL_GENES:
        if g not in genes and g in set(table1["gene"].unique()):
            genes.append(g)

    print(f"Analyzing {len(genes)} genes x {len(CATEGORIES)} categories x "
          f"{len(ANCESTRY_ORDER) + 1} ancestry groups ...", file=sys.stderr)
    richness_df, extrap_df, curves, neutral_df = analyze(
        table1_ident, ancestry_by_person, genes, args.n_bootstrap, args.n_permutations)

    richness_path = os.path.join(out_dir, f"allele_richness{suffix}.tsv")
    extrap_path = os.path.join(out_dir, f"allele_extrapolation{suffix}.tsv")
    neutral_path = os.path.join(out_dir, f"neutral_model_comparison{suffix}.tsv")
    richness_df.to_csv(richness_path, sep="\t", index=False)
    extrap_df.to_csv(extrap_path, sep="\t", index=False)
    neutral_df.to_csv(neutral_path, sep="\t", index=False)

    fig_paths = []
    for category in ["novel", "all"]:
        fig_path = os.path.join(out_dir, f"discovery_curves_{category}{suffix}.png")
        plot_discovery_curves(curves, CLASSICAL_GENES, category, fig_path,
                               f"{category} alleles, classical genes")
        fig_paths.append(fig_path)

    md_path = os.path.join(out_dir, f"04_allele_saturation_report{suffix}.md")
    write_report(md_path, richness_df, extrap_df, neutral_df, fig_paths, is_sample=args.sample)

    print(f"\nWrote {len(richness_df)} richness rows to {richness_path!r}", file=sys.stderr)
    print(f"Wrote {len(extrap_df)} extrapolation rows to {extrap_path!r}", file=sys.stderr)
    print(f"Wrote {len(neutral_df)} neutral-comparison rows to {neutral_path!r}", file=sys.stderr)
    print(f"Wrote {len(fig_paths)} figures + report to {out_dir!r}", file=sys.stderr)


if __name__ == "__main__":
    main()
