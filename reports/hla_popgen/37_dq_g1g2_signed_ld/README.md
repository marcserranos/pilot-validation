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

---

# `--from-committed` mode (2026-09-23): real numbers, from committed 2-field aggregates only

*Added by a third S03 agent, local-only (no browser/VM in this session). Script:
`scripts/hla_popgen/37c_dq_g1g2_from_committed.py`. Answers deliverables 1-3 of the WS1 task using
ONLY the already-committed `reports/hla_popgen/29_hla_ld/ld_pairwise.tsv` -- no VM access, no
recomputation from raw haplotypes, no new number that wasn't already public. Deliverable 4
(artifact-vs-biology checks) still needs the VM (Table 1, relatedness, SR genotypes) and is
untouched by this addendum -- see the WS1 section above.*

## What this mode discovered about the source table

`ld_pairwise.tsv`'s `Dprime` column is **already signed** (verified: `D` and `Dprime` share sign
on every one of the 921 DQA1~DQB1 rows, zero mismatches) -- so no new D' computation was needed;
this script reuses the committed value directly. It also already enforces the exact 20-haplotype-
carrier floor the task asked for: `29_hla_ld_by_ancestry.py`'s `pairwise_ld()` drops any allele
whose own marginal (row or column sum) is below `MIN_ALLELE_HAPS = 20` **before** emitting any row
for it, so every allele appearing anywhere in `ld_pairwise.tsv` already clears the floor by
construction. **No pooled/ALL ancestry row set exists in the committed table** (only per-ancestry:
AFR/AMR/EAS/EUR/MID/SAS) -- confirmed by grep; a pooled panel needs the VM run.

## 1. Main figure + per-ancestry supplement

- `fig_dq_g1g2_committed_MAIN_AFR.png/.pdf` -- **AFR** (4,146 haplotypes -- the most of any
  ancestry for this gene pair, and one of the task's named candidates), Cole's exact layout:
  DQA1\*01 (G2) row block on top, DQB1\*05/06 (G2) column block on the left, black crosshair,
  diverging [-1, 1] "Signed phased D'" colorbar, "Predicted incompatible" text in both
  off-diagonal quadrants, carrier-frequency marginal bars, suppressed/not-estimable cells hatched
  grey.
- `fig_dq_g1g2_committed_supp_{AFR,AMR,EAS,EUR,SAS}.png/.pdf` -- the same layout per ancestry
  (MID excluded per the task's own list -- only 563 haplotypes, too few DQB1 alleles clear the
  floor for a readable panel; still included in the O/E table below with a caveat).
- N (haplotypes) per ancestry, recovered exactly from `n_hap_ij / freq_hap` on uncensored rows
  (internally consistent to floating-point precision -- a genuine identity check, not an
  approximation): **AFR 4146, AMR 2913, EAS 2377, EUR 2383, MID 563, SAS 2147**.
- **Renderer note:** `37_dq_g1g2_signed_ld.py`'s own `fig_g1g2_signed_ld()` has a latent gap for
  this input -- a pair present in the table with status `suppressed_lt20`/`not_observed` is
  neither given a color NOR hatched (it silently renders as an empty NaN cell), because its
  hatching only covers cells absent from the table entirely. `37c` uses its own renderer
  (`fig_g1g2_from_table`) that hatches every non-`estimated` cell, whatever the reason (censored,
  not observed, or simply below an allele's floor and absent from the table) -- this is what makes
  the "Predicted incompatible" quadrants render correctly (see finding below). Iterated once on a
  label-visibility bug (the "Predicted incompatible" text was being drawn under the hatch pattern,
  zorder 3 vs 4 -- fixed by bumping the text to zorder 10) after reading the first PNG.

### Key finding: complete purge at 2-field, no "light-blue" cells

At this resolution, **every cross-group (predicted-incompatible) DQA1xDQB1 allele pair has
`n_hap_ij_disp == "0"` in every ancestry -- not even a single ancestry has one below the 20-count
disclosure floor.** Both off-diagonal quadrants in every panel above are therefore rendered fully
hatched (grey), not the near-uniform deep blue Cole's own target figure shows. This is a *stronger*
purge signal than the target figure at 4-field/allele resolution, where some residual "light blue"
(non−-1) cross-group cells are visible -- plausible since 2-field pools more haplotypes per
cell (better power to detect a truly-zero cross-group count) while 4-field spreads the same
haplotypes across many more distinct alleles (a few of which may show 1 or 2 true co-occurrences
Cole flagged in his own figure). **This means Cole's phase-error-vs-real-recombinant question
(the light-blue cells) cannot be answered at 2-field -- it is specifically a 4-field-resolution
question**, which needs the VM run (deliverable 1's 4-field pass, already implemented in
`37_dq_g1g2_signed_ld.py`, not yet executed).

## 2. Purge statistic (O/E), censored as an interval -- `oe_purge_committed.tsv`

For each ancestry, summed over the cross-group ("predicted incompatible") cells only:
`observed_lower` treats every censored (`<20`) cell as 0 haplotypes, `observed_upper` treats every
one as 19 (task's exact bounds); `expected = N * sum(freq_a * freq_b)` over the same cells under
group-marginal independence. **No point estimate is ever built on a censored cell.**

| ancestry | N haplotypes | cross-group cells | estimated | not_observed | censored (`<20`) | observed | expected | O/E |
|---|---|---|---|---|---|---|---|---|
| AFR | 4146 | 77 | 0 | 77 | 0 | **0 (exact)** | 2018.8 | **0.0** |
| AMR | 2913 | 91 | 0 | 91 | 0 | **0 (exact)** | 1146.3 | **0.0** |
| EAS | 2377 | 107 | 0 | 107 | 0 | **0 (exact)** | 1126.9 | **0.0** |
| EUR | 2383 | 73 | 0 | 73 | 0 | **0 (exact)** | 1111.7 | **0.0** |
| MID | 563 | 35 | 0 | 35 | 0 | **0 (exact)** | 169.5 | **0.0** (small N -- caveat) |
| SAS | 2147 | 86 | 0 | 86 | 0 | **0 (exact)** | 1041.0 | **0.0** |

Every ancestry's cross-group cells are `not_observed` (disclosable exact zero, never censored) --
so `observed_lower == observed_upper == 0` in every row: this is a real, disclosable point value
(0 is not a suppressed count), not an interval collapsing by coincidence. Under independence,
1,000+ cross-group haplotypes were expected in every ancestry; **zero were observed anywhere** --
the strongest possible O/E result the disclosure rules allow reporting exactly.

## 2b. Recurrent cross-group pairs (Cole's "light blue" cells) -- `recurrent_cross_group_pairs.tsv`

**Empty (0 rows).** No cross-group DQA1xDQB1 pair clears the 20-haplotype floor in any ancestry at
2-field resolution -- consistent with the O/E finding above. This list will only be non-empty, if
at all, from the 4-field VM run.

## 3. Bimodality -- `fig_bimodality_committed.png/.pdf`, `bimodality_stats.json`

Pooled over AFR/AMR/EAS/EUR/SAS estimated cells (n=99 compatible pairs; 0 incompatible, per the
finding above, so the originally-planned two-histogram comparison collapsed to one distribution
with an explicit annotation rather than a misleading empty second series):
- **Compatible (G1/G1, G2/G2) signed D' is itself bimodal**: median 0.95, 63% of cells >= 0.90
  (near-perfect coupling), with a second, smaller cluster spread across [0, 0.6] and a handful of
  weakly negative outliers (min ~-0.4) -- i.e. *within* the "compatible" side of Cole's rule, most
  allele pairs are strongly coupled but a nontrivial minority are only weakly associated or even
  mildly repulsed, which is the bimodality he originally noticed in the unsigned/multiallelic
  summary (slides 12-13) showing up again here at the pairwise level.
- The incompatible side contributes no distribution at all (0 estimated cells, see above) -- so at
  2-field, the "bimodal" pattern is a property of the compatible quadrants alone, not a
  compatible-vs-incompatible split.

## Caveats (from-committed mode)

- **2-field (protein) resolution only.** Cole's target figure and the phase-error question are
  4-field; this mode cannot address them (see "Key finding" above). The 4-field pass, G1/G2
  labelling, O/E, and artifact-vs-biology checks are already implemented in
  `37_dq_g1g2_signed_ld.py` and need the VM run.
- **No pooled/ALL ancestry.** Main panel uses AFR (most haplotypes) as the best available
  substitute; a true pooled panel needs the VM run (which pools before the ancestry split).
- MID (563 haplotypes) is thin -- included in the O/E table for completeness but excluded from the
  supplement small-multiple per the task's own named ancestry list.
- The "zero cross-group haplotypes observed anywhere" finding is a real, disclosable result (exact
  0, not a suppressed small count) under the current genotyping/phasing pipeline and this cohort's
  ancestry-stratified sample sizes -- it does not by itself rule out rare recombinants existing
  below what ~2,000-4,000 haplotypes per ancestry can detect; a wider (pooled, VM) N or higher
  resolution could still surface some.
- Unit test: `scripts/hla_popgen/tests/test_37c_dq_g1g2_from_committed.py` (14/14 checks pass) --
  covers status parsing, N recovery, and the O/E censored-interval bounds on synthetic data with a
  hidden true count, including the specific guardrail that an all-censored ancestry must report a
  lower bound of exactly 0.

