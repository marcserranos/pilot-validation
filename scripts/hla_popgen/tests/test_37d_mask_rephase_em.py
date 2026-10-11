#!/usr/bin/env python3
"""Unit tests for scripts/hla_popgen/37d_mask_rephase_em.py, on synthetic data.

Covers the two pieces of new statistical machinery in that script: the Excoffier-Slatkin (1995)
EM haplotype-frequency estimator (`em_haplotype_freqs`) and the most-likely-pair diplotype
assignment (`most_likely_pair`) built on top of it. Does not touch the VM-only `run()` (needs
~/pipeline_outputs and 37_vm_run.py's build_people/extract_haplotypes) -- that is exercised only
on the real cohort, per AGENT_PREAMBLE.md.

Run: python3 scripts/hla_popgen/tests/test_37d_mask_rephase_em.py
"""
import importlib.util
import os
import random
import sys

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


m = _load_module("37d_mask_rephase_em.py", "mask_rephase_em")

passed = 0
failed = 0


def check(name, cond):
    global passed, failed
    if cond:
        passed += 1
    else:
        failed += 1
        print("FAIL:", name)


def make_genotypes_from_true_haplotypes(true_pairs, rng):
    """true_pairs: list of ((a1,b1), (a2,b2)) cis haplotype pairs (one per synthetic person).
    Returns the corresponding unphased genotypes: ((a1,a2) sorted, (b1,b2) sorted)."""
    out = []
    for (a1, b1), (a2, b2) in true_pairs:
        out.append((tuple(sorted([a1, a2])), tuple(sorted([b1, b2]))))
    return out


# ---------------------------------------------------------------------------
# 1. Two haplotypes only, in strong (near-complete) coupling: A1-B1 and A2-B2 co-occur almost
#    always; A1-B2 / A2-B1 (the "repulsion" cis pairs) are rare. EM should recover frequencies
#    close to the true generating frequencies and, more importantly for this script's purpose,
#    most_likely_pair should resolve the doubly-heterozygous genotype to the coupling phase.
# ---------------------------------------------------------------------------
rng = random.Random(1)
n = 4000
true_pairs = []
for _ in range(n):
    r = rng.random()
    if r < 0.55:
        true_pairs.append((("A1", "B1"), ("A2", "B2")))  # coupling, doubly-het (common)
    elif r < 0.57:
        true_pairs.append((("A1", "B2"), ("A2", "B1")))  # repulsion, doubly-het (rare)
    elif r < 0.78:
        true_pairs.append((("A1", "B1"), ("A1", "B1")))  # homozygous both loci -- unambiguous,
        # anchors the EM fit and breaks the label symmetry a purely doubly-heterozygous
        # population would leave unresolved (a real EM degeneracy, not a bug: with every
        # genotype doubly heterozygous and two alleles per locus, the coupling/repulsion
        # resolutions are exactly symmetric and a uniform-start EM cannot break the tie).
        pass
    else:
        true_pairs.append((("A2", "B2"), ("A2", "B2")))  # homozygous both loci, other allele
genos = make_genotypes_from_true_haplotypes(true_pairs, rng)

freq = m.em_haplotype_freqs(genos, seed=1)
check("EM: the two coupling haplotypes (A1-B1, A2-B2) together dominate the two repulsion "
      "haplotypes (A1-B2, A2-B1)",
      (freq[("A1", "B1")] + freq[("A2", "B2")]) > 0.9 * (
          freq[("A1", "B1")] + freq[("A2", "B2")] + freq[("A1", "B2")] + freq[("A2", "B1")]))
check("EM: repulsion haplotypes are rare (<10%) under a population dominated by coupling",
      freq[("A1", "B2")] < 0.10 and freq[("A2", "B1")] < 0.10)
check("EM: frequencies sum to 1", abs(sum(freq.values()) - 1.0) < 1e-6)

# most_likely_pair should recover the coupling phase for a doubly-het genotype under these freqs.
inferred = m.most_likely_pair("A1", "A2", "B1", "B2", freq)
check("most_likely_pair resolves the doubly-het genotype to the dominant (coupling) phase",
      sorted(inferred) == sorted([("A1", "B1"), ("A2", "B2")]))

# ---------------------------------------------------------------------------
# 2. Homozygous-at-one-locus genotypes are unambiguous regardless of freq -- most_likely_pair
#    must return the forced resolution, not go looking at the other (irrelevant) phase option.
# ---------------------------------------------------------------------------
forced = m.most_likely_pair("A1", "A1", "B1", "B2", {("A1", "B1"): 0.01, ("A1", "B2"): 0.99})
check("most_likely_pair: homozygous locus A forces the only possible diplotype",
      sorted(forced) == sorted([("A1", "B1"), ("A1", "B2")]))

# ---------------------------------------------------------------------------
# 3. EM under true independence (no LD) should converge close to the product of marginals --
#    checks the estimator isn't biased toward one resolution by construction.
# ---------------------------------------------------------------------------
rng2 = random.Random(2)
alleles_a = ["X1", "X2"]
alleles_b = ["Y1", "Y2"]
indep_pairs = []
for _ in range(6000):
    a1, a2 = rng2.choice(alleles_a), rng2.choice(alleles_a)
    b1, b2 = rng2.choice(alleles_b), rng2.choice(alleles_b)
    indep_pairs.append(((a1, b1), (a2, b2)))
indep_genos = make_genotypes_from_true_haplotypes(indep_pairs, rng2)
freq_indep = m.em_haplotype_freqs(indep_genos, seed=2)
# Under independence with uniform marginals (~0.5/0.5 each locus), every haplotype should sit
# near 0.25; allow a generous tolerance since this EM only sees unphased genotypes, not the
# original haplotype pairs, and double-heterozygotes are genuinely ambiguous under independence.
check("EM under independence: all four haplotype frequencies within [0.15, 0.35]",
      all(0.15 <= f <= 0.35 for f in freq_indep.values()))

# ---------------------------------------------------------------------------
# 4. Sanity: em_haplotype_freqs handles a population with only unambiguous genotypes (no doubly
#    heterozygous individuals at all) without dividing by zero or crashing.
# ---------------------------------------------------------------------------
unambig_genos = [(("A1", "A1"), ("B1", "B2")), (("A1", "A2"), ("B1", "B1")),
                 (("A2", "A2"), ("B2", "B2"))] * 30
freq_unambig = m.em_haplotype_freqs(unambig_genos, seed=3)
check("EM on an all-unambiguous population still returns a normalized distribution",
      abs(sum(freq_unambig.values()) - 1.0) < 1e-6)

print("\n%d/%d checks passed" % (passed, passed + failed))
if failed:
    sys.exit(1)
