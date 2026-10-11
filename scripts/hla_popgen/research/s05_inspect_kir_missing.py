#!/usr/bin/env python3
"""S05 follow-up to 51 (VM-only; run from ~/s05). On complete KIR haplotypes, tabulates the joint pattern of missing
framework genes (3DP1, 2DL4) with 3DL1/S1 presence and gene count, and the 3DL1+3DS1 co-call rate by ancestry.
Prints aggregate counts; the README reports only patterns/rates (cells <20 are not quoted)."""
import importlib.util
import os
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

import _s05_common as C

spec = importlib.util.spec_from_file_location('k', '51_kir_structure_qc.py')
k = importlib.util.module_from_spec(spec)
sys.modules['k'] = k
spec.loader.exec_module(k)

if __name__ == '__main__':
    root = C.DEFAULT_KIR_OUTROOT
    pids = sorted(d for d in os.listdir(root) if d.isdigit())
    haps = []
    with ProcessPoolExecutor(4) as ex:
        for r in ex.map(k.person_records, [(p, root) for p in pids], chunksize=64):
            haps += r
    coh = C.load_cohort(C.DEFAULT_COHORT)
    anc = dict(zip(coh.person_id, coh.ancestry_pred))
    miss, dup, comp_n = Counter(), Counter(), Counter()
    for pid, hap, calls in haps:
        if k.classify_hap(calls) != 'complete':
            continue
        a = anc.get(pid, 'NA')
        comp_n[a] += 1
        g = {x[1] for x in calls}
        if 'KIR2DL4' not in g or 'KIR3DP1' not in g:
            miss[('3DP1' if 'KIR3DP1' in g else '-3DP1', '2DL4' if 'KIR2DL4' in g else '-2DL4',
                  '3DL1S1' if ('KIR3DL1' in g or 'KIR3DS1' in g) else '-3DL1S1',
                  'n_genes=%d' % len(g) if len(g) <= 6 else 'n_genes>6')] += 1
        if k.pair_event(calls, 'KIR3DL1', 'KIR3DS1'):
            dup[a] += 1
    print('missing patterns:', sorted([(C.mask_n(v), kk) for kk, v in miss.items()], key=str, reverse=True)[:8])
    print('3DL1+3DS1 dup pct by anc:', {a: round(100 * dup[a] / comp_n[a], 2) for a in comp_n if comp_n[a] >= 200})
