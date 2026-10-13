# 57d: do the parents' long-read data drive trio error rates?

Script: `scripts/hla_popgen/57d_trio_parent_data.py` (VM, INTERNAL tables only). Asked by Marc on 2026-10-13: yak builds its
truth from the parents' reads, so the error measured in 57 can come from the parents' side as well as the child's.

- `57d_parent_depth.tsv`: parents by platform, absolute MAPQ-0 depth over the MHC window (from 58) and genome-wide.
- `57d_region_error_vs_parents.tsv`: per 57 target and metric, Spearman of per-trio error with the thinner parent's depth
  (all trios and within child platform), median error by parents' platform mix (groups of at least 20 trios), and OLS of
  log10(error + 0.01%) on child Sequel IIe, number of non-Revio parents, log2 thinner-parent depth and log2 child depth, with
  2,000-draw bootstrap 95% CIs.
- `57d_origin_vs_parents.tsv`: 57b gene-level origin calls; ambiguous share by thirds of thinner-parent depth, and the share of
  resolved calls assigned to the deeper parent where one parent is at least 25% deeper.

The trios-by-platform-mix table is not published: 57's public `trio_selection.tsv` gives the trios per child platform, so its
exact cells would recover the masked ones. Report text: S05 report, trio figures in Part 1–2.
