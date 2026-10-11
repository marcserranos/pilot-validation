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
- **Layout redesigned S04 WS-C phase 2 (2026-09-25)**, replacing the version above after Marc
  rejected the first style pass ("I am not seeing nearly any difference"; overlapping text
  "clearly inadmissible" -- see `sprints/S04_kir_recurrence_style_share/CRITIC_WSC.md`). The old
  layout had independently-sized marginal axes that didn't line up with the heatmap's own
  columns/rows, a colourbar floating far to the right with large empty gaps, and G1/G2 labels
  pushed to a large negative-axes-fraction offset well outside the panel. Fix: one `GridSpec` with
  `sharex`/`sharey` between the heatmap and both marginal bar charts (alignment is now structural,
  guaranteed by matplotlib, not a matched-width/height coincidence to maintain by hand), a compact
  colourbar in its own thin column immediately beside the right marginal, and G1/G2 group brackets
  drawn as plain square brackets in dedicated thin axes directly adjacent to the row/column tick
  labels (`_bracket_v`/`_bracket_h` in `37c_dq_g1g2_from_committed.py`) instead of a large offset.
  183mm width (was 120mm -- needed for the wider bracket + colourbar columns), a short bold panel
  label ("DQ G1/G2, pooled" / ", AFR" / etc., N moved to this README/the corner-note stats below)
  replaces the old sentence-length corner label. New `_viz_common.check_layout()`/`mark_marginal()`
  linter (added the same session) runs on every `save_fig()` call by default (`strict=True`) and
  is the acceptance test for this redesign -- all pooled + per-ancestry DQ figures and the
  bimodality figure pass with zero layout violations.
- **Encoding change**: the two "predicted incompatible" quadrants are no longer painted on the
  same D' colour scale as the informative G1/G1, G2/G2 blocks (previously a solid dark-blue block
  duplicated across both quadrants, most of the figure's ink for zero information). They are now a
  single flat light-grey fill (`NEUTRAL_INCOMPAT_COLOR`, distinct from and unhatched vs. the
  disclosure-censored grey/hatch) with one small annotation giving the actual headline number
  directly, e.g. **"0/469 observed (all 6 ancestries)"** for the pooled panel, or "0/77 observed"
  per ancestry -- read straight from `oe_purge_committed.tsv`'s `n_cells_estimated`/
  `n_cross_group_cells`. A pair genuinely absent from an ancestry's table (never clearing the
  20-haplotype floor for either allele) still gets the usual hatched `SUPPRESSED_COLOR` treatment,
  drawn on top of the neutral fill, so "confirmed zero" and "never disclosable" stay visually
  distinct. Alleles are still ordered by descending carrier frequency within each G1/G2 group.

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

*S04 WS-C phase 2: title shortened from a sentence ("Bimodality within compatible vs incompatible
quadrants (pooled AFR, AMR, EAS, EUR, SAS)") to "DQ G1/G2 bimodality" (the pooled-ancestry list
lives here instead), and the x-axis switched from matplotlib's automatic ~7-tick scale (which
crowded into overlapping tick labels at this panel's 89mm single-column width) to a fixed
[-1, 0, 1] scale -- caught by the new layout linter's text_overlap check.*

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

## Redesign pass 2 (S04 WS-C, 2026-09-25) -- linter false negative + remaining layout faults

Marc rejected the pass-1 redesign above on sight: the rotated "DQA1" y-axis label ran straight
through the G1/G2 bracket's vertical stem, and the "DQB1" x-axis label collided with the bracket's
horizontal cap -- yet `check_layout(strict=True)` had passed every figure with zero violations.

**Root cause**: `check_layout()`'s overlap check only ever compared Text artists against other Text
artists. The G1/G2 brackets (`_bracket_v`/`_bracket_h`) are drawn as plain `ax.plot(...)` Line2D
segments in their own thin axes, never as Text -- so the linter had nothing in its `texts` list to
compare either collision against. It wasn't a tolerance or threshold bug; the check category simply
didn't exist.

**Fix** (`_viz_common.py`): a new `mark_decoration(artist)` helper lets a script opt a Line2D/Patch
into the linter's field of view as a structural "must never be covered by text" object (used here
for the three bracket segments per `_bracket_v`/`_bracket_h` call, and the `axhline`/`axvline`
G1/G2 divider lines). `check_layout()` gained two new error-severity checks: (a2) any Text
overlapping a `mark_decoration()`-registered artist, and (a3) any Text overlapping a visible spine
of an Axes it does not itself belong to (own-axes spines are exempt -- a tick label touching its
own axis is normal). Opt-in rather than blanket Line2D/Patch scanning was deliberate: most
Line2D/Patch objects in a figure ARE the data (bar rectangles, scatter markers) and text
legitimately sits near/on them with no fault; blanket-flagging would have drowned real faults in
false positives. See `_viz_common.mark_decoration()`'s docstring and
`scripts/hla_popgen/tests/test_layout_linter.py`'s `test_text_over_decoration_line_detected` /
`test_text_over_foreign_spine_detected` / `test_bracket_regression_*` for the regression coverage
(these tests replicate the exact fault using `_bracket_v`/`_bracket_h` from this script, not just a
synthetic stand-in).

**Once the linter could see it**, re-running this script surfaced a second, previously-invisible
real fault at the same class: the rotated DQB1 tick labels (e.g. "DQB1\*05:01") are long enough at
90 degrees that they reached down into the bracket axes below and clipped the bracket's own top
rail. Fixed by lowering the bracket's rail (`y=0.26` instead of the default `0.55`) to leave
clearance under the longest tick label.

**Other layout fixes, this pass** (Marc's redesign brief, independent of the linter bug):
- Dropped the redundant `ax_heat.set_xlabel("DQB1")`/`set_ylabel("DQA1")` axis labels entirely --
  every tick label already spells out the full allele name including gene (e.g. "DQA1\*01:02"), and
  the G1/G2 brackets already anchor which axis is which, so the axis label carried no information,
  only a collision surface right next to the bracket.
- Closed the dead band between the bold corner panel label and the top marginal bar chart
  (`fig.text` moved from figure-fraction (0.005, 0.995) to (0.01, 0.985), gridspec `top=` tightened
  0.95 -> 0.91, top row height ratio 1.9 -> 1.5).
- Marginal-axis ticks are now round numbers via `_nice_ceiling()` (0/0.1/0.2, 0/0.125/0.25, ...)
  instead of `max(observed)*1.05` rounded to 2dp (previously e.g. 0/0.155/0.31).

**Encoding change -- fading thin-evidence cells (item 4 of the brief)**: the compatible-quadrant
D'-scale heatmap is mostly dark blue (D'=-1) because most non-haplotype-partner allele *pairs*
genuinely never co-occur (see the bimodality section above) -- that is real signal, not an
artifact. But some of those D'=-1 (or D'=+1) cells are built from very few *expected* haplotypes
under independence (`E = freq_a * freq_b * N_haplotypes`), even though both alleles individually
cleared the 20-haplotype disclosure floor on their own marginals (29's per-allele floor, not a
per-cell one) -- a cell with `E < 5` could easily read `-1.0` from a handful of expected co-
occurrences rather than a well-powered exclusion. Rather than hatching these (hatching is reserved
for hard disclosure censorship, `vc.SUPPRESSED_COLOR`) or swapping the colormap (which would
recolor confidently-estimated cells too, misleadingly), thin cells (`EXPECTED_THIN_THRESHOLD = 5`)
get a translucent white overlay (alpha 0.55) on top of the same D'-scale color: nothing is hidden or
altered, a confidently-estimated D'=-1/+1 simply reads visually darker/more salient than a
thin-evidence one. This does not change any number in any TSV -- it is a rendering-only visual
weight, exactly like the pre-existing censored-cell dot marker.

## Redesign pass 3 (S04 WS-C, 2026-09-26) -- orchestrator full-size review

Three more things were fixed after a full-size (not thumbnail) review:
- **Marginal tick steps didn't match between the two axes**: `_nice_ceiling()` picked a step
  independently for each marginal, so a real max just over a round number (e.g. DQA1's pooled max
  0.193, `*1.05` padding = 0.203) rounded UP to the next step (0.25, ticks 0/0.125/0.25) while the
  other marginal's max rounded to a smaller step (0.2, ticks 0/0.1/0.2) -- two different-looking
  scales side by side on the same figure. Replaced by `_nice_axis_ticks()`, which picks the
  smallest step for which *2 steps* already cover the padded max (and uses a tighter 1.02 pad, not
  1.05), so both marginals land on 0/0.1/0.2 whenever their real maxima are this close.
- **The "faded: E<5" caption was dropped** from the panel entirely -- at the size it had to be to
  fit under the narrow colorbar column, it read as a tiny, orphaned fragment rather than a legend.
  The fading itself is unchanged; its explanation lives only in the "Encoding change" section above
  now, which is where a reader who notices the visual difference and wants to know why would look.
- **The dead band above the top marginal, still visible after pass 2's partial fix**, was closed
  further: gridspec top-row height ratio 1.5->1.05, `top=` 0.91->0.965, corner label pulled from
  y=0.985 to y=0.995.

Re-running after these fixes surfaced a genuinely new (not cosmetic) layout bug the linter had
never been able to see: `ax.axis("off")` on the G1/G2 bracket axes (`ax_brk_y`/`ax_brk_x`) does not
survive a SECOND `fig.canvas.draw()` in this matplotlib version -- `check_layout()`'s own draw (or
`save_fig()`'s) could regenerate that axis's default numeric tick labels ("1", "1.0", ...) from its
raw 0-1 range, rendered right at the seam between the bracket axes and the heatmap, colliding with
each other. Fixed with an explicit private `NullLocator` on both axes' otherwise-unused private
axis (their SHARED axis, with the heatmap, was already fixed the same way in pass 2's Ticker-
replacement code) -- see the script's own comments for the full writeup, since this bug class also
hit Figure 1 v5 independently (documented there too).

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

`37_vm_run.py` on **13,252 unrelated people (later found to be WRONG -- see "Unrelated-set fix"
below; the correct figure is 11,856, matching scripts 29/38/39)**: **17,255 physically-phased
DQA1~DQB1 haplotypes at 4-field, 22,341 at 2-field. Zero cross-group ("predicted incompatible")
haplotypes observed at either resolution, pooled.** `37_vm_run.py`'s own `oe_purge_table.tsv` only ever emitted pooled
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

---

# Critic #1 fix (2026-09-23): unrelated-set correction + non-circular mask-and-rephase

*Fifth S03 agent, fixing CRITIC_1.md blocker #2 and major #3. Ran
`scripts/hla_popgen/37e_unrelated_fix_kfold_em.py` (new script) on the live Workbench VM
(`fix1`), `<=2` cores, plain-text `PUT /api/contents` deploys only. Re-verified against real
cohort data, not a reanalysis of the numbers already above.*

## Unrelated-set fix (blocker)

**The "13,252 unrelated people" reported above was wrong.** `37_vm_run.py`'s own `build_people()`
starts from ALL of `cohort_membership.tsv` (relatives included) and greedily drops one member of
each related PAIR in isolation -- not the maximal-independent-set algorithm
(`24_novelty_by_field.greedy_unrelated`, iteratively removes the highest-remaining-degree person
until no related pair remains) every other S03 result (29, 38, 39) uses on the actual
Table-1/LR-called cohort. Re-deriving the unrelated set the correct way --
`24_novelty_by_field.build_people(t1["person_id"], ...)`, same call 38 makes -- on this same
Table-1 cohort gives **11,856 unrelated people (377 removed for relatedness out of 12,233
Table-1 people)**, exactly matching 38/39's number. This is now the authoritative unrelated count
for script 37; the "13,252" figure above is superseded.

**Deliverable 1 (O/E), restricted to the corrected 11,856-person unrelated set**
(`oe_purge_table_by_ancestry_UNRELATED_FIXED.tsv`):

| ancestry | 2-field N | 2-field observed | 2-field expected | 4-field N | 4-field observed | 4-field expected |
|---|---|---|---|---|---|---|
| AFR | 5383 | **0** | 2686.19 | 4071 | **0** | 2033.64 |
| AMR | 4844 | **0** | 2127.19 | 3924 | **0** | 1710.83 |
| EAS | 2729 | **0** | 1313.62 | 2007 | **0** | 938.96 |
| EUR | 5440 | **0** | 2599.29 | 4335 | **0** | 2093.13 |
| MID | 907 | **0** | 402.22 | 650 | **0** | 295.45 |
| SAS | 2304 | **0** | 1151.99 | 1689 | **0** | 843.25 |
| ALL | 21652 | **0** | 10468.97 | 16713 | **0** | 8066.21 |

**The headline result is unchanged by the fix**: zero observed cross-group cis haplotypes against
10,469 (2-field)/8,066 (4-field) expected under independence, pooled and in every ancestry, on
the correct unrelated cohort. The complete purge was never an artifact of the wrong denominator --
it holds on the right one too, at a slightly smaller but comparable N (21,652 vs. the original
22,341 2-field haplotypes; the difference is exactly the ~1,400-person gap between the two
unrelated-set definitions).

## Non-circular mask-and-rephase (major -- fixes the EM circularity)

**The problem, precisely:** the original mask-and-rephase EM (deliverable 3 above) fit its
population haplotype-frequency model on a pool that INCLUDED the truth-set people it then
rephased and scored. "0 spurious incompatible haplotypes" from an EM trained on a pool that
already contains zero cross-group haplotypes is closer to a tautology than an independent test --
the EM could simply be reproducing information it was given, not demonstrating anything about
real statistical phasers on unseen data.

**Fix: 5-fold held-out EM.** Each ancestry's EM population pool (9,793 people pooled, corrected
unrelated set) is split into 5 folds by a deterministic hash of `person_id`. For each fold, EM
haplotype frequencies are fit on the OTHER 4 folds only, then truth-set people who fall in the
held-out fold are re-phased with those frequencies and scored -- pooled across all 5 folds, no
truth-set person is ever rephased by a model that saw their own data.

**Second, harder stress test: naive linkage-equilibrium (LE) baseline.** Same 5-fold held-out
design, but instead of the EM-fit joint frequencies, phase with independent per-locus marginal
allele frequencies (p_i x q_j) learned from the training fold -- i.e. assume zero real linkage.
Under strict independence, both cis/trans resolutions of a doubly-heterozygous genotype are an
EXACT tie (`p(a1)q(b1)p(a2)q(b2) == p(a1)q(b2)p(a2)q(b1)`), so ties are broken with a seeded coin
flip. This bounds how many incompatible haplotypes a phaser with ZERO real linkage information
would manufacture -- the "prior does all the work" end of the spectrum, against which the
held-out EM's real performance can be judged.

**Results** (`vm_em_mask_rephase_kfold_and_le.tsv`, `vm_fig_kfold_vs_le.png/.pdf`):

| ancestry | scheme | truth evaluated | doubly-het | switch errors | switch rate | spurious incompatible | spurious rate |
|---|---|---|---|---|---|---|---|
| AFR | k-fold EM | 2303 | 1863 | 20 | 1.07% | 22 | 0.48% |
| AMR | k-fold EM | 2136 | 1779 | 20 | 1.12% | `<20` | -- |
| EAS | k-fold EM | 1230 | 1036 | `<20` | -- | `<20` | -- |
| EUR | k-fold EM | 2381 | 2011 | 22 | 1.09% | 28 | 0.59% |
| MID | k-fold EM | 406 | 322 | `<20` | -- | `<20` | -- |
| SAS | k-fold EM | 1052 | 888 | `<20` | -- | `<20` | -- |
| **ALL** | **k-fold EM** | **9530** | **7919** | **52** | **0.66%** | **46** | **0.24%** |
| AFR | LE-naive (no LD) | 2303 | 1863 | 964 | 51.7% | 1178 | 25.6% |
| AMR | LE-naive (no LD) | 2136 | 1779 | 891 | 50.1% | 922 | 21.6% |
| EAS | LE-naive (no LD) | 1230 | 1036 | 535 | 51.6% | 624 | 25.4% |
| EUR | LE-naive (no LD) | 2381 | 2011 | 1005 | 50.0% | 1144 | 24.0% |
| MID | LE-naive (no LD) | 406 | 322 | 168 | 52.2% | 184 | 22.7% |
| SAS | LE-naive (no LD) | 1052 | 888 | 461 | 51.9% | 554 | 26.3% |
| **ALL** | **LE-naive (no LD)** | **9530** | **7919** | **4033** | **50.9%** | **4612** | **24.2%** |

**This changes the verdict from deliverable 3.** The non-circular, held-out EM does NOT reproduce
a perfect 0-spurious-haplotype purge: pooled, it manufactures **46 spurious incompatible cis
haplotypes out of 9,530 x 2 = 19,060 held-out inferred haplotypes (0.24%)** -- small, but
genuinely nonzero, unlike the circular version's exact 0. The naive LE baseline (zero real linkage
information) manufactures spurious incompatible haplotypes at **~24%**, roughly 100x the held-out
EM's rate -- confirming the EM's real population-LD signal is doing substantial, real work (this
cohort's DQA1~DQB1 LD is strong enough that a phaser with no information at all is ~100x worse
than one using it), while also showing the earlier "0 spurious, exactly matching truth" claim
somewhat overstated how perfect a held-out EM actually is on unseen data.

**Revised verdict:** physical phasing still shows a complete purge (0/21,652 at 2-field, corrected
unrelated cohort). A non-circular, held-out population-EM statistical phaser comes very close but
not exact (0.24% spurious rate, pooled) -- closer to Cole's light-blue-cell observation than the
original circular test suggested, though still far short of explaining a large light-blue signal;
a phaser using zero real LD information would manufacture ~100x more spurious pairs (~24%), so the
"assembly/contig-boundary artifact" explanation from deliverable 3 (different-contig implied
pairing, <20/320 incompatible) remains the stronger single explanation for Cole's own figure, with
generic EM statistical-phasing error now a small-but-real secondary contributor rather than a
ruled-out one.

## Caveats (critic-fix run)

- K=5 folds chosen (not 2) to keep each fold's training pool larger for the smaller ancestries
  (MID, SAS); fold assignment is a deterministic SHA-256 hash of `person_id`, reproducible across
  reruns.
- Per-(ancestry, allele-pair) recurrence tables for the spurious cis haplotypes produced by both
  schemes (`vm_spurious_pairs_kfold.tsv`, `vm_spurious_pairs_le_naive.tsv`) stay VM-only -- every
  cell is already `<20` and pulling them back added no disclosure-safe information beyond the
  aggregate table above.
- The LE-naive tie-break's random-number seeding was fixed post-hoc for Python-version
  portability (`random.Random()` no longer accepts an arbitrary tuple as of Python 3.9+; changed
  to a string-formatted seed, behaviorally equivalent -- still an unbiased 50/50 coin flip per
  person, just a different draw of it). The numbers in the table above were produced by the VM's
  Python 3.7 runtime before this portability fix landed; a rerun with the fixed script would give
  a statistically equivalent but not bit-identical LE-naive draw (the k-fold EM numbers, which
  don't use this function, are unaffected).
- Unit tests: `scripts/hla_popgen/tests/test_37e_unrelated_fix_kfold_em.py` (9/9 checks pass) --
  deterministic and balanced fold assignment, marginal-frequency sanity, and the LE-naive tie's
  core correctness property (both resolutions of a doubly-het genotype are picked ~50/50 across
  many people, and the tie-break is itself deterministic per person for reproducibility).
