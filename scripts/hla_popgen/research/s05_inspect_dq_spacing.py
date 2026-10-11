#!/usr/bin/env python3
"""S05 follow-up to 50e (VM-only; run from ~/s05 next to _s05_common.py). Prints aggregate checks only:
how many haplotypes carry >=2 DQB1 copies, and the distribution of same-contig DQA1-DQB1 spacing, including how
many haplotypes have spacing >= 100 kb (the README reports these as <20 / exact quantiles). An INTERNAL companion
(inspect_dq.py on the VM) printed the allele names of the same-contig cross-group haplotypes; it is not kept here."""
import numpy as np
import _s05_common as C

t = C.load_table1(C.DEFAULT_TABLE1)
d = t[t.gene.isin(['DQA1', 'DQB1'])]
cn = d.groupby(['person_id', 'hap', 'gene']).size().unstack(fill_value=0)
print('haps with DQB1 copies>=2:', int((cn.get('DQB1', 0) >= 2).sum()), 'of', len(cn))
w = d[d.copy_index == 1].pivot_table(index=['person_id', 'hap'], columns='gene',
                                      values=['contig', 'gene_start'], aggfunc='first').dropna()
w = w[w[('contig', 'DQA1')] == w[('contig', 'DQB1')]]
dist = (w[('gene_start', 'DQB1')] - w[('gene_start', 'DQA1')]).abs()
print('same-contig DQA1-DQB1 dist quantiles', np.percentile(dist, [1, 5, 50, 95, 99]).astype(int).tolist(),
      'n>=100kb:', C.mask_n(int((dist >= 100000).sum())), 'of', len(dist))
