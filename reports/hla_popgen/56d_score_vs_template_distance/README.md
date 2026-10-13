# 56d: is the 7-point confidence score better than template distance alone?

Script: `scripts/hla_popgen/56d_score_vs_template_distance.py` (VM for the table, `--figure` here). Asked by Marc on
2026-10-13 on the typing-sweep figure.

- `score_vs_distance_sweep.tsv`: 56b's 6,455 MHC-sharing relative pair-genes; error (no shared 2-field allele) and share
  kept for `score_min` (both calls score at least t) and `template_distance_max` (both calls within t edits of the closest
  known allele; thresholds at the 5%–100% quantiles). Rows whose record or error step from the neighbouring kept row is
  1–19 are dropped, so cumulative rows cannot be subtracted into a small count (score 1 and distance 3 are absent).
- `fig_score_vs_distance.png`: error against % kept for both filters (score 7 left off the plot, 0.6% kept).

Reading: the score wins at light filtering (about 90% kept: 2.0% vs 3.0%), the two tie near 70% kept, and distance 0
(an exact known allele) keeps 46% at 1.0%, as clean as score 6 at 36% kept. A distance filter drops novel and rarely
catalogued alleles by construction. Report: S05 report, Part 3–4, next to the typing sweep.
