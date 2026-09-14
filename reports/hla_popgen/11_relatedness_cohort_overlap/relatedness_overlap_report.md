# Relatedness x long-read cohort overlap -- phasing-validation scoping

Long-read (phased) cohort size: **12233**

## Pair counts by bucket

| bucket | n_pairs |
|---|---|
| other | 52621 |
| lr_but_relative_uncalled | 2622 |
| both_lr | 639 |
| lr_sr | 25 |


## Kinship-degree distribution, usable buckets only

| bucket | kin_degree | n_pairs |
|---|---|---|
| both_lr | duplicate_or_MZ_twin | 5 |
| both_lr | first_degree_parentchild_or_sibling | 569 |
| both_lr | second_degree | 65 |
| lr_sr | first_degree_parentchild_or_sibling | 13 |
| lr_sr | second_degree | 12 |


## Coverage of the long-read cohort

- People with >=1 relative who ALSO has a long-read assembly (`both_lr`, the easy/direct-comparison case): **827**

- People with >=1 relative who only has short-read data (`lr_sr`, needs read-realignment): **25**

- People with >=1 usable relative of EITHER kind: **849** out of 12233 long-read-cohort people (6.9%)


**Known limitation:** the relatedness table has no IBD0 column, so `first_degree_parentchild_or_sibling` pairs cannot be split into parent-child vs. full-sibling from this data alone.

**Run provenance:** real VM run, 2026-09-08, against AoU v9 `samples_relatedness.tsv` (55,907 total AoU relative pairs) cross-referenced with `cohort_membership.tsv` (13,252 people; 12,233 in the long-read/Immuannot cohort, 13,228 flagged `in_sr` in this build of `cohort_membership.tsv` -- note this is NOT the full ~500K AoU-native short-read population, see caveat below).

**Important scoping caveat, not yet resolved:** `in_sr` in the current `cohort_membership.tsv` reflects only 13,228 people, far short of AoU's full short-read-WGS population. If that flag were joined against the FULL srWGS release instead, the `lr_sr` bucket (currently only 25 pairs) would very likely grow substantially, since AoU's genome-wide relatedness table itself is computed over the full ~245K+ short-read cohort, not just people who already have AoU-native HLA calls. This does not affect the headline `both_lr` number (827 people), which only depends on both cohort membership flags being complete for the long-read cohort itself.
