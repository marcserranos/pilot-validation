# 37 -- LD supplement table (r2 / D') and DQ/DP heatmap

Source: script `37b_ld_supplement_table.py`, re-rendered from the already-committed
`reports/hla_popgen/29_hla_ld/ld_pairwise.tsv`. No VM access, no re-computation from raw calls --
every number here was already public.

Cole (call #8): "this should be a table too -- a simple CSV with r2 and D-prime", with a distinct
colour for "not observed / below n=20".

## Data-gap check (per the task)

The committed `ld_pairwise.tsv` already carries the **full per-allele-pair signed table**: for
every allele pair (per gene pair, per ancestry) where BOTH alleles individually clear the
20-haplotype floor in their own marginal, it has `r2`, `D`, `Dprime`, `freq_a`, `freq_b`,
`freq_hap`, and the suppressed joint-haplotype count `n_hap_ij_disp`. This is **not** a top-12-only
or multiallelic-summary-only table (compare `ld_multiallelic.tsv`, which IS a whole-table summary
and was not sufficient on its own) -- so this script produces the complete CSV for **all five**
committed gene pairs (A~B, B~C, DPA1~DPB1, DQA1~DQB1, DRB1~DQB1), not only the two the heatmap is
scoped to. **No gap; nothing here is deferred to the VM version.**

## `supp_table_ld_pairs.csv`

Tidy, one row per (gene_pair, ancestry, allele_a, allele_b): `gene_pair`, `ancestry`, `allele_a`,
`allele_b`, `n_hap` (the joint haplotype count, blank/NaN where suppressed), `freq_a`, `freq_b`,
`r2`, `Dprime`, `status`.

`status` is one of:
- **`estimated`** -- the joint haplotype count is disclosable (>=20 or a genuine 0) and `r2`/
  `Dprime` are estimated from real data.
- **`not_observed`** -- the joint count is exactly 0 (this allele pair never co-occurred on a
  phased haplotype in this ancestry). Disclosable; not suppressed.
- **`suppressed_lt20`** -- the joint count is 1-19 under the AoU small-cell rule, written `<20` in
  the source table. **Parsed as censored throughout -- never read as 0 and never as 20.**

## Heatmap panel (DQA1~DQB1, DPA1~DPB1)

One small heatmap per (gene pair, ancestry): rows are the gene pair's first-locus alleles,
columns the second-locus alleles, colour is signed D' on a blue-white-red diverging scale
centered at 0 (`_viz_common.diverging_cmap()`/`diverging_norm()`). Cells with status
`not_observed` or `suppressed_lt20` are drawn as a distinct light-grey hatched cell
(`_viz_common.hatch_suppressed()` / `SUPPRESSED_COLOR`) and are **excluded from the colour scale
entirely** -- they never appear as if they were a real, low-magnitude D' value.

**How to read it:** a solid block of one colour within a gene pair/ancestry panel indicates
strong, consistently-signed linkage; hatched cells mean the underlying allele pair either never
co-occurred in this ancestry's phased haplotypes or its count fell under the 20-haplotype
disclosure floor -- in both cases, no D' is estimated or shown for that cell.

## Caveats

- This is a restyle + tabulation of already-committed, already-suppressed numbers -- no new
  statistic is computed here.
- r2/D' at low allele frequency are noisy; see `29_hla_ld/README.md`'s own caveats on rare-allele
  LD estimates (which this script does not repeat or override).
- `A~B`/`B~C`/`DRB1~DQB1` are in the CSV but not the heatmap panel, which is scoped to
  DQA1~DQB1/DPA1~DPB1 per the task.

---

# WS1 addendum — G1/G2 purge, signed-D' main figure, artifact-vs-recombinant checks

*Added by a separate S03 agent (WS1) on 2026-09-22, after the `37b_ld_supplement_table.py` work
above. Script: `scripts/hla_popgen/37_dq_g1g2_signed_ld.py`. Does NOT supersede
`supp_table_ld_pairs.csv`/`supp_heatmap_dq_dp.{png,pdf}` above — this session had no VM access to
produce a VM-derived version to supersede it with (see "Status" below), so 37b's committed,
already-public-data supplement stands as-is. This addendum instead covers WS1's own scope: the
G1/G2-classified main figure recreating Cole's target layout, the O/E purge statistic, and the
four artifact-vs-biology checks (a)-(d) on incompatible-pair carriers — none of which 37b's script
computes (37b re-renders 29's raw signed D' with no G1/G2 grouping, no crosshair/quadrant layout,
no O/E, no artifact checks).*

## Status: implementation complete, VM run NOT done

`workbench.verily.com` in this session's browser pane redirected straight to a sign-in page — no
authenticated tab was available (contrary to CLAUDE.md's assumption that Marc hands over a
logged-in tab), and entering credentials is outside what this agent does. **No number below is a
real result.** `37_dq_g1g2_signed_ld.py` exits with a clear FATAL message if
`~/pipeline_outputs/hla_calls_rich.tsv` is missing rather than fabricating output.

## What was built (all in `scripts/hla_popgen/37_dq_g1g2_signed_ld.py`)

- `signed_dprime_table()` — signed D' (Lewontin 1964, sign kept, range [-1,1]) per allele pair,
  from 29's `extract_haplotypes` reused verbatim, at 4-field (Cole's resolution) and 2-field, each
  allele needing >=20 haplotype carriers on its own margin, pooled and per ancestry.
- `dq_group()`/`classify_pair()` — Petersdorf 2022 G1 (DQA1*02/03/04/05/06 x DQB1*02/03/04) / G2
  (DQA1*01 x DQB1*05/06) / `predicted_incompatible` classification.
- `fig_g1g2_signed_ld()` — recreates the target PDF's layout exactly: DQA1*01 (G2) row block on
  top, DQB1*05/06 (G2) column block on the left, black crosshair, diverging [-1,1] colormap
  labelled "Signed phased D'", "Predicted incompatible" text in both off-diagonal quadrants,
  carrier-frequency marginal bars, suppressed/not-observed cells hatched grey
  (`_viz_common.hatch_suppressed`), Nature house style (`nature_style()`/`save_fig()`).
- `observed_expected_incompatible()` — O/E for the cross-group ("incompatible") haplotype count
  under a group-marginal independence null, with a bootstrap CI (resample haplotypes, recompute
  O and E each replicate). Pooled and per ancestry; counts <20 written `<20`.
- Four artifact-check functions for deliverable 4 — `check_a_swap_explicable` (VM),
  `build_phasing_confidence_table`/`check_b_phasing_confidence_crosstab` (per-person file stays
  on the VM at `~/s03/results/37/phasing_confidence_per_person.tsv`, never pulled back; only the
  aggregate cross-tab leaves), `check_c_mendelian_transmission` (16's related pairs),
  `check_d_short_read_concordance` (SR genotype table). `run_artifact_checks()` is the VM entry
  point (`--artifact-checks`).
- `cluster_dp_matrix()`/`fig_dp_clustered()` — hierarchical clustering (profile-correlation
  distance, average linkage) of the DPA1xDPB1 signed-D' matrix, looking for a data-driven
  "forbidden block" analogous to G1/G2 — an original angle the task calls out, unrun so whether DP
  shows one is still open.

## Unit test

`scripts/hla_popgen/tests/test_dq_g1g2_signed_ld.py` — 16/16 checks pass locally: signed D' hits
exactly -1 under synthetic perfect repulsion, +1 under perfect coupling, ~0 (|D'|<0.05, n=20k)
under independence; the 20-haplotype floor correctly drops a rare synthetic allele; the G1/G2
classifier is checked on four worked examples (both incompatible directions included); the O/E
purge statistic returns ~1 (CI bracketing 1) on a synthetic independent population and near-0 O/E
(observed=0) on a synthetic fully-purged population.

## Next step

Deploy `37_dq_g1g2_signed_ld.py` to the VM per `VM_CHANNEL.md`, run with `--artifact-checks`, pull
back everything under `~/s03/results/37/` except `phasing_confidence_per_person.tsv`, and fill in
this addendum's real O/E numbers, the artifact-check verdict, and the main/DP-clustered figures.

