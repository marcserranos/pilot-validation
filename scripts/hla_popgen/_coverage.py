#!/usr/bin/env python3
"""Sample-coverage (Good-Turing / Chao & Jost) estimators for HLA allele abundance data.

CONTEXT
-------
Table 1 gives, per gene, a set of alleles each observed on some number of *haplotypes*
(sampling units: this project's HLA calls are per-haplotype, one allele per gene per
haplotype, so counting haplotypes IS counting individuals of abundance data in the
ecological sense). If X_i is the number of haplotypes carrying allele i and
n = sum_i X_i (== m, the number of haplotypes typed for that gene), this is exactly the
"abundance data" setting of Chao & Jost (2012).

WHY COVERAGE, NOT RICHNESS
---------------------------
04_allele_saturation.py / 18_discovery_rate.py already estimate *richness* (how many
alleles exist) via Chao1/Chao2 -- known to be unstable for heavy-tailed allele-frequency
distributions (huge Q1 relative to Q2, "low_confidence_chao2" flag). Good's (1953) /
Good & Toulmin's (1956) insight, formalized by Chao & Jost (2012), is that *coverage* --
the fraction of the population's allele-frequency MASS already represented by catalogued
alleles, equivalently P(the next haplotype sampled carries an already-seen allele) -- is
far better behaved: it depends only on the shape of the observed frequency-of-frequency
spectrum near k=1,2, not on extrapolating an unbounded tail. This module implements
coverage (not richness) estimators for use alongside 18's richness/discovery-rate work,
and extends them to two questions specific to the AoU-validation study: (a) how much of a
NEW population's allele mass would be caught by a reference catalogue built from a
DIFFERENT population (cross_coverage), and (b) how a multi-population sampling design's
catalogue coverage behaves as population-specific subsample sizes vary
(pooled_design_coverage).

KEY FORMULAS AND CITATIONS
---------------------------
All hypergeometric-type ratios C(a,k)/C(b,k) (a,b,k possibly in the thousands) are
evaluated in log-space via scipy.special.gammaln, exactly as `_log_comb` /
`rarefied_richness` do in `_viz_common.py` and `exact_rate` does in
`18_discovery_rate.py` -- direct binomial-coefficient ratios overflow long before n
reaches AoU-scale haplotype counts.

1. Sample coverage (Chao & Jost 2012, eq. 1 as refined via Good-Turing with f1/f2
   bias correction; also Chao, Ma & Yang 1993):

       f1 = # singleton alleles (X_i == 1), f2 = # doubleton alleles (X_i == 2), n = sum(X)

       Chat = 1 - (f1/n) * A,   A = (n-1) f1 / ((n-1) f1 + 2 f2)                    (f2 > 0)
                                A = (n-1) f1 / ((n-1) f1 + 2)          (Chao & Jost's f2=0
                                                                          convention)
       Chat = 1.0 if f1 == 0 (no singletons -> nothing suggests undercatalogued mass)

2. Rarefied (interpolated) coverage for a hypothetical subsample of size m < n, drawn
   without replacement from the n haplotypes already typed (Chao & Jost 2012 eq. 4; the
   exact closed form used here matches the `Chat.Ind`-family estimator in the iNEXT
   reference implementation):

       Chat(m) = 1 - sum_i (X_i/n) * C(n - X_i, m) / C(n - 1, m)                (0 < m < n)

   IMPORTANT AMBIGUITY, RESOLVED AND VERIFIED (task instructions explicitly flagged this):
   the C(n-1, m) denominator above is *not* interchangeable with the more "obvious"
   C(n, m) denominator that would come from directly weighting "P(allele i's copies are
   present among the m sampled)" by allele i's population frequency X_i/n over ALL n
   haplotypes -- those are two different quantities. We verified algebraically (and by
   brute-force enumeration in tests/test_coverage.py::test_interpolation_brute_force)
   that the C(n-1,m) formula equals

       Chat(m) = E[ (1/(n-m)) * #{h not in subsample : allele(h) appears in subsample} ]

   i.e. the expected FRACTION OF THE HELD-OUT HAPLOTYPES whose allele was already seen in
   the m-subsample -- the direct finite-population generalization of Good-Turing "P(next
   draw is a repeat)". This is the quantity iNEXT/Chao & Jost actually publish as rarefied
   coverage, and it is what reduces, at m=1, to the hand-derivable
   sum_i (X_i/n)*(X_i-1)/(n-1) ("probability the single next haplotype after the one just
   drawn repeats its allele"). The alternative ("what fraction of the FULL n-haplotype
   sample's allele-frequency mass is covered by alleles present in an m-subsample", using
   C(n,m) in the denominator) is a genuinely different, also well-defined quantity, but it
   is NOT what coverage_at_size computes; we discard it here to match the published
   Chao & Jost / iNEXT estimator, and both are exercised/contrasted in the tests.

   At m == n the formula above is 0/0 (C(n-1,n) == 0); coverage_at_size special-cases
   m == n to return sample_coverage(X) exactly, which is the correct limit.

3. Extrapolated coverage beyond the observed sample, m = n + m* (m* > 0) (Chao & Jost
   2012, the abundance-data analogue of their incidence extrapolation; consistent with
   sample_coverage at m*=0):

       Chat(n + m*) = 1 - (f1/n) * A^(m* + 1)

   Chao & Jost note this extrapolation is only trustworthy out to about m <= 2n-3n; beyond
   that the geometric-decay assumption embedded in A is unverifiable from the data at hand.
   `size_for_coverage` returns a `reliable` flag for exactly this reason.

4. Expected richness at size m (Colwell et al. 2012 / iNEXT; the abundance-data twin of
   `rarefied_richness` in `_viz_common.py`, reused here to keep `richness_at_size` and
   `coverage_at_size` self-consistent within this module):

       interpolation (m <= n):   S(m) = sum_i [1 - C(n - X_i, m) / C(n, m)]
       extrapolation (m = n+m*): S(n+m*) = S_obs + f0_hat * (1 - (1 - f1/(n f0_hat + f1))^m*)

       Chao1 (Chao 1984, bias-corrected form Chao 1987):
           f0_hat = (n-1)/n * f1^2 / (2 f2)              if f2 > 0
                  = (n-1)/n * f1 (f1 - 1) / 2             if f2 == 0

   Unlike coverage_at_size, richness_at_size's interpolation formula needs NO special case
   at m == n: plugging m = n makes every term C(n-X_i, n)/C(n, n) = 0/1 = 0 for X_i >= 1,
   so S(n) = S_obs automatically.

5. `cross_coverage` / `pooled_design_coverage`: not from a single citation -- straight
   applications of hypergeometric "probability allele i is absent from an m-subsample of a
   population with X_i copies out of N total", P(absent) = C(N - X_i, m) / C(N, m),
   combined either as a single reference population (cross_coverage) or as a product
   over several independently-subsampled populations, "allele i counted as catalogued if
   observed in ANY population's subsample" (pooled_design_coverage). Both are documented
   in-line at the function.

ASSUMPTIONS
-----------
* X is abundance data over one gene at a time: one allele call per haplotype, so
  n = sum(X) is both the number of haplotypes and the number of "individuals" in the
  ecological sense. Mixing genes (different n) in one call is meaningless.
* All counts are treated as exact (no genotyping-error model); QC/warning filtering is
  the caller's responsibility (upstream Table 1 loaders already do this).
* No plotting here by design -- this module is pure numerics, consumed by figure scripts
  the way `_viz_common.py`'s statistics helpers are.
"""
import numpy as np
from scipy.special import gammaln

# ---------------------------------------------------------------------------
# Log-space combinatorics (same approach as _viz_common._log_comb /
# 18_discovery_rate.exact_rate: gammaln, never a direct binomial-coefficient ratio).
# ---------------------------------------------------------------------------


def _log_comb(n, k):
    """log C(n, k); -inf for an out-of-range (n, k) rather than raising, so callers can mask
    with np.exp(...) == 0 instead of branching everywhere."""
    n = np.asarray(n, dtype=float)
    k = np.asarray(k, dtype=float)
    out = np.full(np.broadcast(n, k).shape, -np.inf)
    ok = (k >= 0) & (k <= n) & (n >= 0)
    nb, kb = np.broadcast_arrays(n, k)
    out = np.where(ok, gammaln(nb + 1) - gammaln(kb + 1) - gammaln(nb - kb + 1), -np.inf)
    return out


def _comb_ratio(a, b, k):
    """exp(log C(a,k) - log C(b,k)), 0 wherever C(a,k) is out of range (a<k)."""
    log_a = _log_comb(a, k)
    log_b = _log_comb(b, k)
    with np.errstate(invalid="ignore"):
        return np.where(np.isfinite(log_a), np.exp(log_a - log_b), 0.0)


def _prep(X):
    X = np.asarray(X, dtype=float)
    X = X[X > 0]
    n = float(X.sum())
    f1 = float(np.sum(X == 1))
    f2 = float(np.sum(X == 2))
    return X, n, f1, f2


def _A_hat(n, f1, f2):
    """The geometric-decay ratio shared by sample_coverage's f1/f2 correction and the
    extrapolation formula. f2==0 uses Chao & Jost's own substitute (their eq. for f2=0)."""
    if f1 <= 0:
        return 0.0
    if f2 > 0:
        return ((n - 1) * f1) / ((n - 1) * f1 + 2 * f2)
    return ((n - 1) * f1) / ((n - 1) * f1 + 2)


# ---------------------------------------------------------------------------
# 1. Sample coverage
# ---------------------------------------------------------------------------


def sample_coverage(X):
    """Chao & Jost (2012) sample-coverage estimator for abundance data X (1-D array of
    positive haplotype counts per allele). Returns a float in [0, 1]."""
    X, n, f1, f2 = _prep(X)
    if n == 0:
        return float("nan")
    if f1 == 0:
        return 1.0
    A = _A_hat(n, f1, f2)
    return float(1.0 - (f1 / n) * A)


# ---------------------------------------------------------------------------
# 2 & 3. Coverage at an arbitrary sample size (interpolation + extrapolation)
# ---------------------------------------------------------------------------


def _interp_coverage_scalar(X, n, m):
    if m <= 0:
        return 0.0
    log_denom = _log_comb(n - 1, m)
    log_num = _log_comb(n - X, m)
    with np.errstate(invalid="ignore"):
        ratio = np.where(np.isfinite(log_num), np.exp(log_num - log_denom), 0.0)
    return float(1.0 - np.sum((X / n) * ratio))


def coverage_at_size(X, m):
    """Expected sample coverage of a (sub/super-)sample of size m, per the module
    docstring's formulas 2-3. `m` may be a scalar or an array; returns the same shape
    (scalar in -> scalar out, via a 0-d/1-elem array check)."""
    X, n, f1, f2 = _prep(X)
    if n == 0:
        return np.nan if np.isscalar(m) else np.full_like(np.asarray(m, dtype=float), np.nan)
    scalar_in = np.isscalar(m) or np.ndim(m) == 0
    m_arr = np.atleast_1d(np.asarray(m, dtype=float))
    C_hat = sample_coverage(X)
    A = _A_hat(n, f1, f2)
    out = np.empty_like(m_arr)
    for idx, mm in enumerate(m_arr):
        if mm == n:
            out[idx] = C_hat
        elif mm < n:
            out[idx] = _interp_coverage_scalar(X, n, mm)
        else:
            m_star = mm - n
            out[idx] = 1.0 - (f1 / n) * (A ** (m_star + 1)) if f1 > 0 else 1.0
    return float(out[0]) if scalar_in else out


# ---------------------------------------------------------------------------
# Smallest sample size reaching a target coverage
# ---------------------------------------------------------------------------


def size_for_coverage(X, target, max_reliable_factor=2.0):
    """Smallest integer sample size m (m may exceed n) with expected coverage >= target.

    Returns (m, reliable): m is np.inf when the extrapolation model cannot reach `target`
    (f1 > 0 and the geometric-decay ratio A >= 1, or target >= 1.0 exactly, which this
    model only approaches in the limit). `reliable` is False whenever m falls beyond
    `max_reliable_factor` * n -- Chao & Jost recommend not trusting extrapolation past
    about 2x-3x the observed sample size.
    """
    X, n, f1, f2 = _prep(X)
    if n == 0:
        return np.inf, False
    C_hat = sample_coverage(X)
    if target <= C_hat:
        lo, hi = 1, int(n)
        # coverage_at_size is monotone non-decreasing in m on [1, n] (verified in tests);
        # plain integer binary search for the smallest m clearing `target`.
        while lo < hi:
            mid = (lo + hi) // 2
            if coverage_at_size(X, mid) >= target:
                hi = mid
            else:
                lo = mid + 1
        return int(lo), True
    # target > C_hat: needs extrapolation, hence f1 > 0 (f1 == 0 => C_hat == 1.0 >= target).
    if f1 <= 0 or target >= 1.0:
        return np.inf, False
    A = _A_hat(n, f1, f2)
    if A <= 0 or A >= 1:
        return np.inf, False
    rhs = n * (1.0 - target) / f1
    if rhs <= 0:
        return np.inf, False
    m_star = np.log(rhs) / np.log(A) - 1.0
    m_star = max(m_star, 0.0)
    m = int(np.ceil(n + m_star))
    # Guard against floating-point rounding landing just short of target.
    while coverage_at_size(X, m) < target:
        m += 1
    reliable = m <= max_reliable_factor * n
    return m, reliable


# ---------------------------------------------------------------------------
# 4. Expected richness at an arbitrary sample size (iNEXT-style; abundance-data twin of
#    _viz_common.rarefied_richness / 18_discovery_rate's incidence extrapolation).
# ---------------------------------------------------------------------------


def _chao1_f0(n, f1, f2):
    if f1 <= 0:
        return 0.0
    if f2 > 0:
        return ((n - 1) / n) * (f1 ** 2) / (2 * f2)
    return ((n - 1) / n) * f1 * (f1 - 1) / 2.0


def _interp_richness_scalar(X, n, m):
    if m <= 0:
        return 0.0
    log_denom = _log_comb(n, m)
    log_num = _log_comb(n - X, m)
    with np.errstate(invalid="ignore"):
        ratio = np.where(np.isfinite(log_num), np.exp(log_num - log_denom), 0.0)
    return float(np.sum(1.0 - ratio))


def richness_at_size(X, m):
    """Expected number of distinct alleles in a (sub/super-)sample of size m, per the
    module docstring's formula 4. `m` may be scalar or array."""
    X, n, f1, f2 = _prep(X)
    S_obs = len(X)
    if n == 0:
        return 0.0 if np.isscalar(m) else np.zeros_like(np.asarray(m, dtype=float))
    scalar_in = np.isscalar(m) or np.ndim(m) == 0
    m_arr = np.atleast_1d(np.asarray(m, dtype=float))
    f0 = _chao1_f0(n, f1, f2)
    out = np.empty_like(m_arr)
    for idx, mm in enumerate(m_arr):
        if mm <= n:
            out[idx] = _interp_richness_scalar(X, n, mm)
        else:
            m_star = mm - n
            if f0 <= 0 or f1 <= 0:
                out[idx] = S_obs
            else:
                g = f1 / (n * f0 + f1)
                out[idx] = S_obs + f0 * (1.0 - (1.0 - g) ** m_star)
    return float(out[0]) if scalar_in else out


# ---------------------------------------------------------------------------
# 5. Cross-population coverage
# ---------------------------------------------------------------------------


def cross_coverage(counts_ref, counts_target, m_ref=None):
    """How much of population B's (`counts_target`) allele-frequency mass is expected to
    be caught by a catalogue built from a random subsample of `m_ref` haplotypes of
    population A (`counts_ref`)?

        cc = sum_i p_B(i) * P(allele i observed in an m_ref-subsample of A)
        p_B(i) = X_Bi / n_B
        P(observed) = 0                                   if X_Ai == 0
                    = 1 - C(n_A - X_Ai, m_ref) / C(n_A, m_ref)   otherwise

    Returns (cc, unseen_in_B) where unseen_in_B = 1 - sample_coverage(counts_target) is
    B's OWN estimated undetected allele-frequency mass (alleles B hasn't even sampled
    itself yet). Interpretation: cc is the fraction of B's *observed* frequency spectrum
    that A's catalogue would already cover; the true cross-coverage of B's full
    (including unobserved) population is bounded above by
    cc * (1 - unseen_in_B) + unseen_in_B (A can, at best, also catch some of B's unseen
    mass, but this module makes no attempt to estimate how much -- report both numbers and
    let the caller reason about the gap rather than over-fitting a correction).

    counts_ref / counts_target: dict allele_id -> positive haplotype count.
    m_ref: subsample size of A; default (None) is all of A (m_ref = n_A). Must be
    <= n_A (without-replacement subsample of an already-observed sample).
    """
    n_A = float(sum(counts_ref.values()))
    n_B = float(sum(counts_target.values()))
    if n_A == 0 or n_B == 0:
        return 0.0, 1.0
    if m_ref is None:
        m_ref = n_A
    if m_ref > n_A:
        raise ValueError(f"m_ref ({m_ref}) cannot exceed n_A ({n_A}): a subsample cannot "
                          f"exceed the observed reference sample.")
    cc = 0.0
    for allele_id, X_Bi in counts_target.items():
        p_Bi = X_Bi / n_B
        X_Ai = counts_ref.get(allele_id, 0)
        if X_Ai == 0:
            p_obs = 0.0
        elif m_ref == n_A:
            p_obs = 1.0
        else:
            p_obs = 1.0 - float(_comb_ratio(n_A - X_Ai, n_A, m_ref))
        cc += p_Bi * p_obs
    target_counts = np.array(list(counts_target.values()), dtype=float)
    unseen_in_B = 1.0 - sample_coverage(target_counts)
    return float(cc), float(unseen_in_B)


# ---------------------------------------------------------------------------
# 6. Pooled multi-population sampling design coverage
# ---------------------------------------------------------------------------


def pooled_design_coverage(counts_by_pop, sizes_by_pop, target_weights):
    """Expected coverage of a design that draws `sizes_by_pop[p]` haplotypes (without
    replacement) from each population p's observed sample `counts_by_pop[p]`, where an
    allele counts as "catalogued" if it is observed in ANY population's subsample:

        P(allele i observed) = 1 - prod_p [ C(n_p - X_pi, s_p) / C(n_p, s_p) ]

    Coverage of target population q's own allele-frequency mass:

        coverage_q = sum_i p_q(i) * P(allele i observed),   p_q(i) = X_qi / n_q

    `target_weights`: dict population -> weight (must sum to ~1) naming the populations to
    report coverage for and how to combine them into a single pooled/global figure.

    Returns (coverage_by_pop, global_coverage): coverage_by_pop has one entry per key of
    `target_weights`; global_coverage = sum_q target_weights[q] * coverage_by_pop[q].

    Raises ValueError if any sizes_by_pop[p] exceeds that population's observed n_p.
    """
    if abs(sum(target_weights.values()) - 1.0) > 1e-8:
        raise ValueError(f"target_weights must sum to 1, got {sum(target_weights.values())}")
    n_by_pop = {p: float(sum(c.values())) for p, c in counts_by_pop.items()}
    for p, s_p in sizes_by_pop.items():
        if s_p > n_by_pop.get(p, 0.0):
            raise ValueError(f"sizes_by_pop[{p!r}] ({s_p}) exceeds observed n_p "
                              f"({n_by_pop.get(p, 0.0)}); cannot subsample more than observed.")

    allele_ids = set()
    for c in counts_by_pop.values():
        allele_ids.update(c.keys())

    def p_not_observed_any(allele_id):
        prod = 1.0
        for p, n_p in n_by_pop.items():
            s_p = sizes_by_pop.get(p, 0.0)
            if n_p == 0 or s_p <= 0:
                continue
            X_pi = counts_by_pop[p].get(allele_id, 0)
            if X_pi == 0:
                continue
            prod *= float(_comb_ratio(n_p - X_pi, n_p, s_p))
        return prod

    p_observed = {a: 1.0 - p_not_observed_any(a) for a in allele_ids}

    coverage_by_pop = {}
    for q in target_weights:
        n_q = n_by_pop.get(q, 0.0)
        if n_q == 0:
            coverage_by_pop[q] = 0.0
            continue
        cov_q = 0.0
        for allele_id, X_qi in counts_by_pop[q].items():
            p_qi = X_qi / n_q
            cov_q += p_qi * p_observed.get(allele_id, 0.0)
        coverage_by_pop[q] = float(cov_q)

    global_coverage = sum(target_weights[q] * coverage_by_pop[q] for q in target_weights)
    return coverage_by_pop, float(global_coverage)


# ---------------------------------------------------------------------------
# 7. Multinomial bootstrap CI
# ---------------------------------------------------------------------------


def bootstrap_ci(X, fn, n_boot=200, seed=0):
    """Simple multinomial-resampling bootstrap CI for a statistic `fn` of abundance data X.

    Each replicate resamples n = sum(X) haplotypes with per-allele probability X_i/n
    (i.e. resamples FROM THE OBSERVED CATALOGUE ONLY), rebuilds a positive-count vector,
    and evaluates fn on it. Returns (lo, hi), the 2.5th/97.5th percentiles across
    n_boot replicates.

    LIMITATION (documented, not fixed here): because resampling draws only from alleles
    already observed, it can never generate a "new" unseen allele, so this systematically
    UNDERSTATES uncertainty for anything richness-related (an unseen-species correction
    would require a full Chao-style bootstrap that also injects a modeled unseen mass,
    which Chao & Jost themselves note is overkill for routine use). It is adequate for
    coverage-type statistics, whose main uncertainty driver (the observed f1/f2 counts) is
    exactly what this resampling reproduces.
    """
    X, n, _, _ = _prep(X)
    if n == 0 or len(X) == 0:
        return float("nan"), float("nan")
    p = X / n
    rng = np.random.default_rng(seed)
    n_int = int(round(n))
    vals = []
    for _ in range(n_boot):
        boot_counts = rng.multinomial(n_int, p)
        boot_counts = boot_counts[boot_counts > 0]
        if len(boot_counts) == 0:
            continue
        vals.append(fn(boot_counts))
    vals = np.asarray(vals, dtype=float)
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi)
