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
on every one of the 927 DQA1~DQB1 rows, zero mismatches) -- so no new D' computation was needed;
this script reuses the committed value directly. It also already enforces the exact 20-haplotype-
carrier floor the task asked for: `29_hla_ld_by_ancestry.py`'s `pairwise_ld()` drops any allele
whose own marginal (row or column sum) is below `MIN_ALLELE_HAPS = 20` **before** emitting any row
for it, so every allele appearing anywhere in `ld_pairwise.tsv` already clears the floor by
construction. **Critically, `freq_a`/`freq_b`/`freq_hap` are disclosed as exact floats for EVERY
row regardless of the haplotype-count disclosure status** (`n_hap_ij_disp` is `"0"`,
`"<20"`, or an exact integer) -- so `Dprime` is a real, already-public value for every row in the
table, not just the `>=20` ones: 0/927 rows have a NaN `Dprime`. **No pooled/ALL ancestry row set
exists in the committed table** (only per-ancestry: AFR/AMR/EAS/EUR/MID/SAS).

*Correction (2026-09-23, orchestrator review): the first version of this mode treated every
non-`>=20` cell as "not estimable" and hatched it grey, including the entire incompatible
quadrants -- wrong, because D' is defined and disclosed for those cells too (a `not_observed`
cell's D' is deterministic given the marginals, and turns out to be exactly -1.0 on all 721
DQA1~DQB1 `not_observed` rows in this table). Fixed below: colour is shown whenever D' is present;
a small black dot marks a cell whose count is disclosed only as `<20` (not an exact integer);
hatching is now reserved for allele pairs genuinely absent from a given ancestry's table (which
only happens in the pooled reconstruction -- see below).*

## 1. Main figure (pooled) + per-ancestry supplement

- `fig_dq_g1g2_committed_MAIN_POOLED.png/.pdf` -- **pooled across all 6 ancestries**
  (14,529 haplotypes). D' depends only on allele/haplotype frequencies, and pooled frequencies are
  exact N-weighted means of the disclosed per-ancestry frequencies (`pool_across_ancestries()`),
  so pooled D' is reconstructed exactly wherever an ancestry's own table has the pair (small
  undercount only possible for the rare allele that falls *below its own* 20-carrier floor in one
  ancestry while clearing it elsewhere -- bounded by <20 haplotypes out of 14,529, i.e. <0.14%,
  per such ancestry/allele). Of the 236 allele pairs that appear in at least one ancestry's table,
  **2 pairs (4 cells) are absent from every ancestry** (never co-occurring with enough support
  anywhere) and are hatched grey; every other cell is coloured.
- `fig_dq_g1g2_committed_supp_{AFR,AMR,EAS,EUR,SAS}.png/.pdf` -- same layout per ancestry (MID
  excluded per the task's own named list -- only 563 haplotypes).
- N (haplotypes) per ancestry, recovered exactly from `n_hap_ij / freq_hap` on uncensored rows:
  **AFR 4146, AMR 2913, EAS 2377, EUR 2383, MID 563, SAS 2147** (sum 14,529).
- Nature-style polish per orchestrator review: no in-figure title (ancestry + N in a small corner
  label instead), 120mm width, 4.5-6pt text, thin light-grey marginal bars with fixed 0/0.1/0.2/0.3
  ticks, a short colourbar labelled "Signed phased D'", and G1/G2 group-block axis labels ("G2 α
  (DQA1\*01)" / "G1 α (DQA1\*02-06)", "G2 β (DQB1\*05/06)" / "G1 β (DQB1\*02/03/04)"). Iterated
  twice after reading the PNGs: (1) the "Predicted incompatible" text was drawn under the hatch
  pattern (zorder 3 vs 4) -- fixed by bumping text to zorder 10; (2) the group-block labels
  initially overlapped the allele tick labels -- fixed by pushing them to axes-fraction x=-0.62
  (free, since `save_fig`'s `bbox_inches="tight"` expands the canvas to fit).

### Key finding: complete purge, and it reads correctly now (uniform deep blue, per Cole's own figure)

Both off-diagonal ("Predicted incompatible") quadrants render as **uniform deep blue (D' = -1.0
exactly)** in every panel, matching the look of Cole's own target figure -- not grey/hatched as in
the pre-fix version. This is because **every cross-group DQA1xDQB1 allele pair has
`n_hap_ij_disp == "0"` in every ancestry** (721/721 `not_observed` cross-group cells checked; 0
censored) -- a real, disclosable exact zero, and D' = -1.0 follows deterministically from the
marginals (verified: 100% of `not_observed` rows have `Dprime == -1.0` to floating precision, not
merely negative). **No cell anywhere in the incompatible quadrant is disclosed only as `<20`** (the
dots visible in the compatible quadrants are all within-group), so there are no "light-blue"
(non_-1, low-but-nonzero-count) cross-group cells to be found at 2-field resolution -- Cole's
phase-error-vs-real-recombinant question is specifically a **4-field-resolution** question, which
needs the VM run (already implemented in `37_dq_g1g2_signed_ld.py`, not yet executed).

## 2. Purge statistic (O/E), censored as an interval -- `oe_purge_committed.tsv`

Unchanged by the rendering fix (this statistic always used true haplotype counts, never D',
so it was already correct): for each ancestry, summed over the cross-group cells only,
`observed_lower` treats every censored (`<20`) cell as 0 haplotypes, `observed_upper` treats every
one as 19; `expected = N * sum(freq_a * freq_b)` under group-marginal independence.

| ancestry | N haplotypes | cross-group cells | estimated | not_observed | censored (`<20`) | observed | expected | O/E |
|---|---|---|---|---|---|---|---|---|
| AFR | 4146 | 77 | 0 | 77 | 0 | **0 (exact)** | 2018.8 | **0.0** |
| AMR | 2913 | 91 | 0 | 91 | 0 | **0 (exact)** | 1146.3 | **0.0** |
| EAS | 2377 | 107 | 0 | 107 | 0 | **0 (exact)** | 1126.9 | **0.0** |
| EUR | 2383 | 73 | 0 | 73 | 0 | **0 (exact)** | 1111.7 | **0.0** |
| MID | 563 | 35 | 0 | 35 | 0 | **0 (exact)** | 169.5 | **0.0** (small N -- caveat) |
| SAS | 2147 | 86 | 0 | 86 | 0 | **0 (exact)** | 1041.0 | **0.0** |

Every ancestry's cross-group cells are `not_observed` (disclosable exact zero, never censored), so
`observed_lower == observed_upper == 0`: a real point value (0 is not a suppressed count), not an
interval collapsing by coincidence. Under independence, 1,000+ cross-group haplotypes were expected
in every ancestry; **zero were observed anywhere**.

## 2b. Recurrent cross-group pairs (Cole's "light blue" cells) -- `recurrent_cross_group_pairs.tsv`

**Empty (0 rows).** "Recurrent" here means an exact disclosed count of >=20 (this pipeline's usual
disclosure threshold) -- no cross-group DQA1xDQB1 pair reaches that anywhere. More precisely (see
Key Finding above), literally **no** cross-group pair has *any* nonzero disclosed count (not even a
censored `<20` one) in any ancestry at 2-field. This list will only be non-empty, if at all, from
the 4-field VM run.

## 3. Bimodality -- `fig_bimodality_committed.png/.pdf`, `bimodality_stats.json`

*Corrected alongside the figure fix: now built from every disclosed cell (not just `>=20` ones),
pooled over AFR/AMR/EAS/EUR/SAS. Two side-by-side panels (compatible / incompatible) rather than
one overlaid histogram, since both groups now have a tall spike at exactly -1.0 that would
otherwise occlude each other.*

- **n=421 compatible (G1/G1, G2/G2) cells, n=434 predicted-incompatible cells** -- both now
  populated (versus 99/0 in the pre-fix version).
- **Incompatible quadrant: a single spike at D' = -1.0 (100% of cells)** -- the purge is complete
  and total, exactly as the O/E table shows.
- **Compatible quadrants are genuinely bimodal**: ~63% of the mass sits at/near D' = -1.0 (median
  overall -1.0) -- these are non-partner allele combinations *within* the same G1 or G2 group (e.g.
  two different G1 haplotypes that never co-occur) -- with a second, smaller cluster near +1.0
  (14.7% of cells >= 0.90 -- real allele-level "partner" pairs, e.g. the common DQA1\*01:02~
  DQB1\*05:01/06:02 haplotypes) and a thin spread across [-0.5, 0.5] in between. **This is the real
  bimodality Cole saw**: the compatible/G1G1-or-G2G2 label only says two alleles are in the same
  epitope group, not that they are the specific haplotype-partner alleles that actually travel
  together -- most same-group allele *pairs* still never co-occur (D'=-1), and only the true
  haplotype partners cluster near +1.

## Caveats (from-committed mode)

- **2-field (protein) resolution only.** Cole's target figure and the phase-error question are
  4-field; this mode cannot address them. The 4-field pass, G1/G2 labelling, O/E, and
  artifact-vs-biology checks are already implemented in `37_dq_g1g2_signed_ld.py` and need the VM
  run.
- **Pooled reconstruction, not a true pooled computation.** `pool_across_ancestries()` sums
  N-weighted per-ancestry frequencies; an allele/pair missing from one ancestry's table (because it
  fell *below that ancestry's own* 20-carrier floor) contributes 0 rather than its true (<20) count
  to the pooled sum -- a bounded undercount, negligible at N=14,529, but not a byte-for-byte
  re-derivation from raw haplotypes (which needs the VM). Verified exactly on synthetic data with
  no missing-ancestry edge case (`test_pool_across_ancestries_exact_reconstruction`).
- MID (563 haplotypes) is thin -- included in the O/E table for completeness but excluded from the
  supplement small-multiple per the task's own named ancestry list.
- The "zero cross-group haplotypes observed anywhere" finding is a real, disclosable result under
  the current genotyping/phasing pipeline and this cohort's sample sizes -- it does not by itself
  rule out rare recombinants existing below what ~2,000-4,000 haplotypes per ancestry (14,529
  pooled) can detect; a higher-resolution (4-field, VM) pass could still surface some.
- Unit test: `scripts/hla_popgen/tests/test_37c_dq_g1g2_from_committed.py` (18/18 checks pass) --
  status parsing, N recovery, the O/E censored-interval bounds (including the all-censored ->
  lower-bound-exactly-0 guardrail), and exact pooled-D' reconstruction on synthetic data.

---

# VM run (2026-09-23): real cohort, 4-field, WS1 deliverables 1-3 -- and a verdict

*Fourth S03 agent, WS1 continuation. Ran on the live Workbench VM (`ws1b`) via the programmatic
JupyterLab REST/websocket channel (`VM_CHANNEL.md`), `<=2` cores, plain-text `PUT /api/contents`
deploys only, no base64. Deployed and ran the already-committed `scripts/hla_popgen/37_vm_run.py`
(a prior agent's condensed, dependency-free standalone runner -- the committed
`37_dq_g1g2_signed_ld.py` needs `29_hla_ld_by_ancestry.py`/`16_phasing_mendelian_validation.py`/
`24_novelty_by_field.py`, which this VM's `~/repos/pilot-validation` checkout (a different branch)
does not have; `37_vm_run.py` avoids that by reimplementing the needed logic inline). Added and
ran `scripts/hla_popgen/37d_mask_rephase_em.py`, the mask-and-rephase experiment (deliverable 3).*

## Deliverable 1 -- O/E, pooled + per ancestry, 4-field and 2-field, real cohort

`37_vm_run.py` on 13,252 unrelated people: **17,255 physically-phased DQA1~DQB1 haplotypes at
4-field, 22,341 at 2-field. Zero cross-group ("predicted incompatible") haplotypes observed at
either resolution, pooled.** `37_vm_run.py`'s own `oe_purge_table.tsv` only ever emitted pooled
("ALL") rows, because its `build_people()` picks the ancestry column by substring match
(`"anc" in colname`), and `max_template_distance`/`mean_template_distance` both contain the
substring `"anc"` (dist**anc**e) and sort before `ancestry_pred` in `cohort_membership.tsv`'s
column order -- so every person was silently bucketed by a template-distance value ('1.0', '2.0',
..., NaN), not AFR/AMR/EAS/EUR/MID/SAS. Verified directly against `cohort_membership.tsv`.
`37d_mask_rephase_em.py` re-derives ancestry correctly from `ancestry_pred` and recomputes the O/E
table per ancestry (`vm_oe_purge_table_by_ancestry.tsv`, `vm_fig_oe_purge_by_ancestry.png/.pdf`):

| ancestry | 4-field N | 4-field observed | 4-field expected | 2-field N | 2-field observed | 2-field expected |
|---|---|---|---|---|---|---|
| AFR | 4192 | **0** | 2093.2 | 5545 | **0** | 2765.4 |
| AMR | 4107 | **0** | 1792.2 | 5066 | **0** | 2223.1 |
| EAS | 2056 | **0** | 959.7 | 2796 | **0** | 1345.5 |
| EUR | 4484 | **0** | 2167.4 | 5626 | **0** | 2692.6 |
| MID | 662 | **0** | 300.8 | 920 | **0** | 407.8 |
| SAS | 1717 | **0** | 856.9 | 2343 | **0** | 1171.4 |
| ALL | 17255 | **0** | 8320.1 | 22341 | **0** | 10793.4 |

Zero observed against 300-8320 expected under group-marginal independence, in **every** ancestry,
**every** resolution, on the full real cohort -- the "any incompatible cis haplotype at all, no
floor" version of deliverable 1 finds none. `37_vm_run.py`'s artifact-check (a) (swap-explicable
rate) is consequently vacuous: 0 incompatible carriers means nothing to test
(`vm_run_summary.json`).

## Deliverable 3 -- the mask-and-rephase experiment (the key result)

`37d_mask_rephase_em.py` (new script this session): **truth set** = 9,967 people with DQA1 and
DQB1 physically phased (same contig) on *both* hap1 and hap2 (2-field; 4-field truth set is a
subset, thinner). For each ancestry, fit a 2-locus multiallelic EM (Excoffier & Slatkin 1995,
`em_haplotype_freqs`) on the unphased genotypes of every unrelated person of that ancestry with two
calls at each locus (10,107 pooled), then re-phase each truth-set person with the higher-posterior
diplotype (`most_likely_pair`) and compare to their known true phase.

**Result (`vm_em_mask_rephase_summary.tsv`, `vm_fig_mask_rephase_comparison.png/.pdf`): EM
statistical rephasing manufactured *zero* spurious incompatible cis haplotypes, pooled or in any
ancestry (0/9,835 x 2 = 0/19,670 inferred haplotypes) -- matching truth exactly on this axis.**
Switch errors (EM resolves a doubly-heterozygous person to the wrong phase entirely, regardless of
G1/G2) were rare and mostly below the 20-count disclosure floor per ancestry; pooled, **20/8,174
doubly-heterozygous truth people (switch rate 0.24%)**. So on this cohort, at this LD strength, a
population-level EM haplotype estimator does not manufacture the kind of cross-group cis pair
Cole's target figure shows as light blue -- **the null hypothesis (statistical-phasing artifact)
is not supported by this specific EM reconstruction.**

**The real-world "different contig" group tells a different story.** For people whose DQA1 and
DQB1 calls sit on different contigs within the same assembly hap (no physical cis evidence, so any
pairing is inferred from the pipeline's own hap1/hap2 label rather than measured): 320 such
pooled hap-level pairs exist, of which **fewer than 20 are G1/G2-incompatible** (exact count
withheld -- the numerator itself is below the 20-person disclosure floor, so no rate is reported;
see `vm_diffcontig_implied_pairing.tsv`). This is **not zero**, unlike both the physical-phasing
truth set and the EM-rephased set. Per-ancestry breakdowns (AFR 119 pairs, EUR 71, AMR 62, SAS 30,
EAS 27, MID <20) are all individually below the floor for the incompatible count, consistent with
AFR/EUR (thinner contig assembly, more fragmentation) carrying more of this signal, but not
statistically resolvable person-by-person at this N.

### Verdict

Physical phasing shows a complete purge (0/17,255 at 4-field, real cohort). A standard
population-EM statistical phaser, run on this cohort's own allele frequencies, reproduces that
purge almost exactly (0 spurious incompatible haplotypes; ~0.2% switch-error rate unrelated to
G1/G2). The one regime that does produce incompatible-looking cis pairs is naive use of the
assembly's own hap1/hap2 label when the two genes are not on the same contig -- a small but
nonzero fraction (<20/320 pooled, so roughly single-digit percent) of haplotype pairs. **Read
together: Cole's light-blue cells are more consistent with assembly/contig-boundary artifacts
(genes on different contigs, phase unknown, and a naive hap-label pairing used anyway) than with
generic EM statistical-phasing error** -- though this cohort's own statistically-phased DQ calls
(if produced by something other than the population-EM tested here) were not directly examined,
so this is evidence about the *class* of statistical-phasing error, not a re-analysis of Cole's
own pipeline's specific phasing method.

## Caveats (VM run)

- `37_vm_run.py` is a previously-deployed, independently-written condensed runner, not the
  committed `37_dq_g1g2_signed_ld.py` -- its logic was spot-checked against the committed script's
  docstring/functions (`extract_haplotypes`, `dq_group`/`classify_pair`, `observed_expected_incompatible`)
  during this session but not line-by-line diffed against it; both compile and both reused for the
  real O/E numbers above.
- The `build_people()` ancestry-column bug (see deliverable 1) affects only `37_vm_run.py`'s *own*
  `oe_purge_table.tsv` (pooled-only) and its unused-here artifact-check groupings; `unrelated_ids`
  (the relatedness filter) is unaffected and was reused as-is. Not fixed in `37_vm_run.py` itself
  (another agent's file, out of this session's direct scope) -- flagged here and worked around in
  `37d_mask_rephase_em.py` by re-deriving `anc_of` from `cohort_membership.tsv`'s `ancestry_pred`
  column directly.
- EM was run at 2-field only (4-field would need a much larger truth set per ancestry to converge
  reliably given how many more distinct alleles enter the multiallelic EM at that resolution --
  not attempted this session).
- All released rates are gated so that **both** the numerator and denominator clear the 20-count
  disclosure floor -- a rate computed from a `<20` numerator over a disclosed denominator would let
  a reader back-compute the suppressed count (e.g. `0.0026 x 1926 approx 5`), defeating the
  masking. This was caught and fixed mid-session (an earlier pulled TSV had this leak; the
  committed `vm_*.tsv` files here do not).
- Deliverable 4's artifact checks (b) phasing-confidence file, (c) Mendelian transmission, (d)
  short-read concordance were not run this session (deliverable 1's zero-incompatible-carriers
  result makes (a) trivially null and leaves nothing for (b)-(d) to test against real incompatible
  carriers; they remain implemented in `37_dq_g1g2_signed_ld.py`/available in `37_vm_run.py` for a
  future session if a nonzero incompatible-carrier set is ever found, e.g. from a different gene
  pair or a larger cohort).
- Deliverable 4's per-person phasing-confidence file was not built this session (no incompatible
  carriers to build it around, per above).
- Figures render slightly non-final (the O/E figure's expected-vs-observed annotation arrow could
  be tightened) -- numbers are final; a follow-up pass could tighten layout only.
- Unit test: `scripts/hla_popgen/tests/test_37d_mask_rephase_em.py` (7/7 checks pass) -- the EM
  estimator recovers coupling-dominant frequencies and resolves a doubly-heterozygous genotype to
  the correct (coupling) phase on synthetic data; a purely-doubly-heterozygous synthetic
  population is deliberately *not* asserted to break symmetry (a real EM degeneracy under two
  alleles per locus with no anchoring homozygotes, not a bug); an independence-generated
  population converges near the uniform product-of-marginals frequencies; an all-unambiguous
  population does not crash or divide by zero.
