#!/usr/bin/env python3
"""Unit tests for scripts/hla_popgen/37e_unrelated_fix_kfold_em.py, on synthetic data.

Covers the new statistical machinery specific to 37e (not already covered by
test_37d_mask_rephase_em.py, which tests em_haplotype_freqs/most_likely_pair, reused verbatim
here): deterministic fold assignment (`_fold_of`), the marginal (linkage-equilibrium) frequency
estimator (`marginal_freqs`), and the naive-LE tie-breaking phaser (`le_naive_pair`) -- in
particular, the formal tie property the whole LE-baseline design depends on (both cis/trans
resolutions of a doubly-het genotype have EXACTLY equal weight under independence).

Does not touch the VM-only `run()` (needs ~/pipeline_outputs, 37_vm_run.py and
24_novelty_by_field.py's build_people) -- that is exercised only on the real cohort, per
AGENT_PREAMBLE.md.

Run: python3 scripts/hla_popgen/tests/test_37e_unrelated_fix_kfold_em.py
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


m = _load_module("37e_unrelated_fix_kfold_em.py", "unrelated_fix_kfold_em")

passed = 0
failed = 0


def check(name, cond):
    global passed, failed
    if cond:
        passed += 1
    else:
        failed += 1
        print("FAIL:", name)


# ---------------------------------------------------------------------------
# 1. _fold_of: deterministic (same pid -> same fold every call/process), and roughly balanced
#    across folds for a reasonably sized id set (not all piled into one fold).
# ---------------------------------------------------------------------------
ids = [str(i) for i in range(2000)]
K = 5
folds_run1 = [m._fold_of(pid, K) for pid in ids]
folds_run2 = [m._fold_of(pid, K) for pid in ids]
check("_fold_of is deterministic across repeated calls", folds_run1 == folds_run2)
check("_fold_of only returns values in [0, K)", all(0 <= f < K for f in folds_run1))
counts = [folds_run1.count(k) for k in range(K)]
check("_fold_of roughly balances 2000 ids across 5 folds (each fold 300-500)",
      all(300 <= c <= 500 for c in counts))

# ---------------------------------------------------------------------------
# 2. marginal_freqs: simple allele-copy-counting sanity check.
# ---------------------------------------------------------------------------
genos = [
    (("A1", "A1"), ("B1", "B2")),
    (("A1", "A2"), ("B1", "B1")),
    (("A2", "A2"), ("B2", "B2")),
]
p, q = m.marginal_freqs(genos)
# locus A copies: A1,A1,A1,A2,A2,A2 -> 3/6 each
check("marginal_freqs: locus A frequencies sum to 1", abs(sum(p.values()) - 1.0) < 1e-9)
check("marginal_freqs: locus A is 50/50 A1/A2 for this synthetic set",
      abs(p["A1"] - 0.5) < 1e-9 and abs(p["A2"] - 0.5) < 1e-9)
check("marginal_freqs: locus B frequencies sum to 1", abs(sum(q.values()) - 1.0) < 1e-9)

# ---------------------------------------------------------------------------
# 3. le_naive_pair: the core correctness property the whole LE baseline depends on is that a
#    doubly-het genotype is an EXACT tie under independence -- verify both resolutions are
#    picked roughly 50/50 across many distinct (seed, pid) draws, i.e. the tie-break is not
#    silently biased toward one side.
# ---------------------------------------------------------------------------
p_indep = {"A1": 0.3, "A2": 0.7}
q_indep = {"B1": 0.4, "B2": 0.6}
n_coupling = 0
n_trials = 4000
for i in range(n_trials):
    pair = m.le_naive_pair("A1", "A2", "B1", "B2", p_indep, q_indep, pid=f"person{i}", seed=0)
    if sorted(pair) == sorted([("A1", "B1"), ("A2", "B2")]):
        n_coupling += 1
frac = n_coupling / n_trials
check("le_naive_pair: tie-break lands close to 50/50 across many distinct people (%.3f)" % frac,
      0.44 <= frac <= 0.56)

# Determinism: the same (seed, pid) must always resolve the same way (so a rerun of the pipeline
# reproduces the same spurious/switch counts, even though the tie-break is itself "random").
r1 = m.le_naive_pair("A1", "A2", "B1", "B2", p_indep, q_indep, pid="fixed_person", seed=7)
r2 = m.le_naive_pair("A1", "A2", "B1", "B2", p_indep, q_indep, pid="fixed_person", seed=7)
check("le_naive_pair is deterministic for a fixed (seed, pid)", r1 == r2)

# Unambiguous genotypes (homozygous at one locus) must return the forced resolution regardless
# of the coin flip -- same contract as most_likely_pair in 37d.
forced = m.le_naive_pair("A1", "A1", "B1", "B2", p_indep, q_indep, pid="x", seed=0)
check("le_naive_pair: homozygous locus A forces the only possible diplotype",
      sorted(forced) == sorted([("A1", "B1"), ("A1", "B2")]))

print("\n%d/%d checks passed" % (passed, passed + failed))
if failed:
    sys.exit(1)
