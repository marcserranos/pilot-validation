# 40 — Figure 1 v5

*WS4, S03 sprint. Implements Cole's 2026-09-22 panel feedback
(`sprints/CALL_SUMMARY_2026-09-22.md` Sec 6), then the orchestrator's 2026-09-23 review of the
first v5 draft. Built by `scripts/hla_popgen/40_figure1_v5.py` (+ `40a_admixture_bins.py` for
panel a).*

## 2026-09-23 follow-up: panel c is still mostly hatched, and why

The censoring fix above made panel c *correct* but revealed it is also mostly hatched: the
committed `24_novelty_by_field` table splits each gene x ancestry cell into several
independently-censored sub-cells (by field class and artifact label), so even a well-powered
combined total often has *some* sub-cell below 20. **`scripts/hla_popgen/40b_novelty_rate_export.py`**
fixes this the right way: it recomputes the combined totals directly from Table 1 on the VM and
applies the `<20` rule exactly once, to the number actually written out, with Wilson 95% CIs from
the true counts. `40_figure1_v5.py compose` prefers `panel_c_novelty_totals.tsv` when present
(falls back to the heavier-censored rendering otherwise) and auto-selects the most specific metric
(protein > CDS > any-field) that clears 20 in at least half of the 48 cells (`--c-metric` to
override). **2026-09-24: this VM export has now landed** (`panel_c_novelty_totals.tsv`, see
"Panels a and c: what they show now" below) — panel c is composed from it.

## 2026-09-24: panels a and c are no longer placeholders

Both VM exports panel a and c were waiting on now exist locally
(`panel_a_admixture_bins.tsv`, `panel_c_novelty_totals.tsv`), already aggregate-only and
disclosure-checked. `compose` was re-run against them (see "Panels a and c: what they show now"
below for exact thresholds and the panel-d legend overlap fix that came out of reviewing the
result). **All six panels (a-f) are now real, populated data** — no placeholder boxes remain.

## The 2026-09-22 correctness bug and how it's fixed

The orchestrator's review caught panel c rendering `MID x HLA-A` as a bare **"0"** — a point value
computed from `_bounds()`'s lower bound (censored `<20` cells counted as 0 in the numerator). This
is the exact fabricated-zero pattern that burned an earlier Figure 1 (`feedback_suppressed_counts_are_not_zero`
in project memory: a `<20` cell silently became 0 and produced a bogus 0.00% row).

**Fix:** every cell in panel c is now checked for a censored numerator contribution
(`n_censored_cells > 0` in `33_figure1_v3_compose.py::rate_by_gene_ancestry`'s own accounting).
If any of that cell's underlying field-class counts were written `<20`, the cell is drawn
**hatched** and labelled as an upper-bound interval, `"≤x"` — never a bare number. Only cells with
zero censored contribution show a plain point value (e.g. `HLA-B x AFR = 12`, exact). In practice
this means *most* cells in the classical-gene table are hatched, because the underlying per-field
counts (especially the protein-field slice) are frequently below 20 even when the gene x ancestry
total is well powered — that is itself an honest finding, not a defect: the true rate could be
anywhere up to the printed upper bound, and the figure now says so instead of guessing 0.
Denominator-too-small cells (`n_called_upper < 20`) are separately hatched and marked `n/a`.

Panel d's stacked counts (novel protein alleles by recurrence class) are **not** subject to this
same check: they are counts of distinct *alleles*, not of participants, per the established
convention already documented in `36_figure1_native.py`'s own panel-d docstring. A recurrence-class
legend entry is dropped rather than shown if it sums to zero across the 8 displayed genes (the
`≥20 unrelated people` class is in fact all-zero for these classical genes — the reportable novel
proteins are concentrated in TAP1/TAP2/MIC, per `36_figure1_native.py`'s own finding — so that
entry does not appear).

## What "absent from IPD-IMGT/HLA" means in panel c, and which metric was chosen

Three metrics exist in this project for "novelty":
1. **Any-field / genomic**: sequence differs from the catalogue at *any* nomenclature field
   (protein, synonymous CDS, or non-coding only) — what 33/36 plotted, and what Cole referred to
   as "the 800,000 one." Dominated by non-coding (intronic/UTR) differences: e.g. HLA-DRB1's
   ≈65% rate is mostly field-4 (non-coding-only) novelty, not new protein sequence.
2. **Synonymous-CDS-level**: coding sequence differs but translates to a known protein.
3. **Protein-level (field 2)**: the translated protein itself is not in IPD-IMGT/HLA — the
   scientifically strongest claim, and the one panel (d) counts.

**Panel c uses metric 1 (any-field/genomic)**, explicitly labelled `"% called haplotypes with
sequence (any field) absent from IPD-IMGT/HLA — mostly non-coding"` so the axis itself states the
granularity. This was chosen over protein-level for the **main panel** because at protein-level
resolution nearly every gene x ancestry cell's numerator is below 20 (novel protein alleles are
rare enough that few reach 20 carriers within one ancestry — see panel d, 21-60 total per gene
pooled across *all six* ancestries), so a protein-level heatmap is almost uniformly
hatched/uninformative at this granularity. The any-field rate is well powered (few fully-censored
cells) and is literally the panel Cole endorsed. Panel (d), immediately adjacent with the same gene
rows, carries the protein-level claim instead — so both metrics are on the figure, each where it is
legible. The full three-way field breakdown remains in `33_figure1_v3/` and `24_novelty_by_field/`.

## Panels a and c: what they show now, and the exact thresholds used

**Panel a** (`panel_a_admixture_bins.tsv`): 11,833 unrelated participants, binned into **589
consecutive bins of >=20 people each**, sorted by predicted genetic ancestry then by
dominant-ancestry-component probability (descending) within ancestry. Each bin's plotted value is
the **mean** of the six admixture proportions (`p_afr` ... `p_sas`) across its >=20 members —
disclosure-safe by construction, since no bin can be smaller than the AoU small-cell floor and no
individual's proportions are ever exported. This uses **predicted** ancestry (soft assignment,
continuous probability), distinct from the **strict** threshold panel c uses below.

**Panel c** (`panel_c_novelty_totals.tsv`): 48 gene x ancestry cells (8 classical genes x 6
ancestries), restricted to participants at **strict ancestry >=0.9** (a single dominant-ancestry
posterior probability threshold — looser than panel b's 0.98, see "Open issues" #2). Per cell,
`n_novel_*` counts 1-19 are written as the string `"<20"` and never resolved to a point value;
rates and Wilson 95% CIs are computed only where both the numerator and denominator clear 20. Of
the 48 cells, only **1 (HLA-A x MID, n_total=590) has a censored numerator** at the any-field
level, and it renders hatched with the label `<20` — never as `0`. At the protein-coding level
only **1/48 cells clears >=20 novel-protein carriers**, so `compose`'s auto metric-selection
(protein > CDS > any-field, first to clear >=50% of cells) correctly falls through to **any-field**
for the main heatmap — the same choice documented above under "which metric was chosen," now
confirmed by the real counts rather than assumed.

**Headline numbers from the real data:** HLA-DRB1 has by far the highest any-field novelty rate in
every ancestry (46-71%, driven by non-coding diversity in a gene that is itself hypervariable), and
the highest single cell is **HLA-DRB1 x EAS at 70.9%** (95% CI 69.1-72.7%, n=2,431). Among the
7 non-DRB1 genes, HLA-DPB1 is highest, peaking at **HLA-DPB1 x EAS = 27.2%** (n=2,438). The
lowest rates are at HLA-A (2.3-4.4% across the 5 uncensored ancestries).

**One rendering bug found and fixed this round**: panel d's recurrence-class legend
(`seen once` / `2-19 unrelated people`) was originally anchored inside the axes at
`loc="lower right"`, which — once real data replaced the placeholder — landed directly on top of
the HLA-B and HLA-A bars (both long, bottom two rows). Moved the legend fully below the x-axis
(`bbox_to_anchor=(0.5, -0.16)`, `loc="upper center"`); no more overlap with any bar or its label.

## Layout: fixed 183 x 150 mm canvas, orchestrator's row spec

- Row 1 (~25 mm): panel a, full width.
- Row 2 (~58-62 mm): b (ternary, ~55 mm) | c (heatmap) + d (marginal bar, same gene rows) | a
  small vertical colorbar (~3 x 30 mm, manually sized via `ax.set_position()` after gridspec
  layout — a plain gridspec cell would give it the full 62 mm row height).
- Row 3 (~55 mm): e (discovery curves) | f (DQ G1/G2 O/E purge).

Two variants were built (`layout_A`, `layout_B` — B gives the heatmap a touch more width and the
ternary a touch less). **A is chosen**: B's extra heatmap width added no legibility (8 rows already
read cleanly), and its narrower ternary crowded the tick labels.

No caption text is drawn inside the figure — the legend below is the only place it lives, per the
orchestrator's instruction.

## Panel f: DQ G1/G2 O/E purge, chosen over a KIR preview

`reports/hla_popgen/41_kir_scoping/` has **no aggregate data yet** — it is scoping-only, pending
the WS3 KIR rerun. `reports/hla_popgen/37_dq_g1g2_signed_ld/oe_purge_committed.tsv` is already
committed, disclosure-cleared, and well powered (563-4,146 DQA1~DQB1 cis haplotypes per ancestry).
It shows **zero observed G1xG2 cross-group haplotypes** against an expectation of 170-2,019 under
independence, in every one of the 6 ancestries — a clean, striking, already-public result. DQ O/E
was the clear choice on data availability alone.

## Draft figure legend (Nature style)

**Figure 1 | Long-read HLA typing across the All of Us cohort: cohort structure, reference
catalogue gaps, and allele-discovery saturation by ancestry.**
**(a)** Genetic-ancestry composition of unrelated long-read participants (n=11,833), binned into
589 consecutive groups of ≥20 people (sorted by predicted ancestry, then by dominant-ancestry
probability descending) — mean admixture proportion per bin.
**(b)** HLA-B alleles (≥20 carriers, strict ancestry probability ≥0.98) on the AFR/EUR/AMR simplex
at the mean renormalised ancestry composition of their carriers; point area/opacity ~ carrier
count (log). Triangles: alleles first observed in this cohort. **(c)** Percentage of called
haplotypes carrying sequence, at any IPD-IMGT/HLA nomenclature field, absent from the catalogue
(mostly non-coding), per classical gene (rows) x ancestry (columns), strict ancestry ≥0.9. Hatched
cells: the numerator is censored (<20), so only the "<20" label is shown, never a point value or a
bare "0" (1/48 cells: HLA-A x MID). **(d)** Distinct protein alleles not in
IPD-IMGT/HLA, same gene rows as c, ancestry-pooled, stacked by unrelated-carrier recurrence (seen
once / 2-19 / ≥20); ancestry-split version in `34_novel_recurrence/`. **(e)** Per-ancestry
allele-discovery curves (8 classical genes pooled), mean ± 95% band over 25 random orderings per
ancestry (`pred` ancestry scheme), direct end labels. **(f)** DQA1~DQB1 cis-haplotypes observed in
Petersdorf et al. (2022) G1xG2 cross-group cells, vs. expected under independence, by ancestry —
essentially zero in every group, consistent with background selection purging the
non-functional heterodimer.

## N per panel

- a: 589 bins of ≥20 unrelated people each, 11,833 people total, predicted-ancestry scheme.
- b: 53 HLA-B alleles, strict ancestry ≥0.98, ≥20 carriers.
- c/d: 8 classical genes x 6 ancestries, strict ancestry ≥0.9 for c; per-cell n_total 580-4,504.
  d's per-gene novel-protein-allele totals: HLA-A 34, HLA-B 60, HLA-C 51, HLA-DPA1 21, HLA-DPB1 32,
  HLA-DQA1 38, HLA-DQB1 30, HLA-DRB1 35.
- e: 6 ancestries, 25 permutations/point, N from 487 (MID) to 3,020 (AFR).
- f: 6 ancestries, 563-4,146 DQA1~DQB1 cis-haplotypes each.

## Distilled

- All six panels (a-f) are now real, aggregate, disclosure-checked data — no placeholders remain.
- Panel a: 589 bins (≥20 people each) x 11,833 unrelated people, predicted-ancestry scheme, mean
  admixture proportions only.
- Panel c: strict ancestry ≥0.9, any-field novelty rate (auto-selected because protein/CDS are
  almost entirely <20-censored); only 1/48 cells (HLA-A x MID) is itself censored, hatched and
  labelled `<20`, never drawn as `0`.
- Headline: HLA-DRB1 x EAS has the highest any-field novelty rate on the figure, 70.9%
  (n=2,431); HLA-DPB1 x EAS is highest among the 7 non-DRB1 genes, 27.2% (n=2,438).
- One bug fixed this round: panel d's legend overlapped the HLA-B/HLA-A bars once real data
  replaced the placeholder — moved below the axis.
- Open: panel b (0.98) and panel c/d (0.9) use different strict-ancestry thresholds — not unified,
  out of scope for this figure (would need a script-24 rerun at 0.98).

## S04 WS-C phase 2 redesign (2026-09-25) -- 5 layout-linter errors fixed

The S04 diagnostic `check_layout()` pass (`sprints/S04_kir_recurrence_style_share/FIGURES_INDEX.md`)
found 5 error-severity layout violations in both composed layouts (A and B) of this figure:
- `'100'`/`'0'` ternary tick numerals in panel b overlapping the bold `AFR`/`AMR` vertex labels
  (the tick sits exactly at the vertex, same place the vertex's own name is anchored).
- Panel b's marker legend (`bbox_to_anchor=(1.14, -0.08)`, placed outside the axes toward panel c)
  overlapping panel c's rotated two-line y-axis label, and separately the `EUR` vertex label also
  overlapping that same y-axis label.

**Fixes** (`draw_ternary_grid`, `draw_panel_b`, `_METRIC_LABEL` in `40_figure1_v5.py`):
- Ternary tick numerals: only the midpoint ("50") is labelled per edge now. The 0/100 endpoints are
  already unambiguous from the bold vertex name; the numeral added nothing but a collision surface.
- Panel c/d y-axis label shortened from a two-line sentence ("% called haplotypes with sequence
  (any field) absent from IPD-IMGT/HLA -- mostly non-coding") to one direct line ("any-field
  novelty (%)"). The second line was the actual cause of the cross-panel collision: a rotated
  (90°) matplotlib Text with an embedded newline lays its second line out WIDTH-wise (perpendicular
  to the vertical reading direction), so a 2-line rotated label is measurably wider, not taller --
  that extra width is what reached left into panel b's space. The full metric definition (why
  any-field vs CDS vs protein, what "mostly non-coding" means) stays in this README's "What
  'absent from IPD-IMGT/HLA' means" section above -- a direct label states the unit, the caveat
  belongs in prose, not on the axis (de-AI checklist item 7/13).
- Panel b's marker legend moved from outside the axes (`(1.14, -0.08)`) to below the triangle's own
  base edge (`loc="upper center", bbox_to_anchor=(0.5, -0.02)`, `ncol=2`). The padding band below
  the triangle's bottom edge is guaranteed free of both data points and other panels' artifacts by
  construction (ternary data only ever falls inside the triangle), so this placement can't collide
  with panel c regardless of exact figure sizing, and also stopped the legend from sitting on top
  of scattered data markers (a fault the mechanical linter can't see, but a reviewer would).

No panel content, statistic, or number changed -- every value in panels a-f still comes from the
same committed TSVs (`panel_a_admixture_bins.tsv`, `panel_b_ternary_alleles.tsv`,
`panel_c_novelty_totals.tsv`, `39_saturation_by_ancestry/curves.tsv`, `oe_purge_committed.tsv`);
`git diff` on `40_figure1_v5.py` for this pass touches only the label text and the ternary-grid
tick loop, none of the data-loading/aggregation functions.

`check_layout(strict=True)` (the `save_fig()` default) now passes both composed layouts (A and B)
with **0 errors, 0 warnings**.

## Open issues

1. Panel b (strict ancestry ≥0.98) and panel c/d (strict ancestry ≥0.9, from
   `panel_c_novelty_totals.tsv`) use different thresholds — not unified across the figure; would
   need a VM rerun of script 24 at 0.98, out of scope for this pass.
2. Panel e's Clench extrapolation is being corrected concurrently elsewhere; this script only
   plots raw curve points, which are unaffected, but a final diff against
   `39_saturation_by_ancestry/README.md` is worth doing before submission.
3. AFR/AMR end labels in panel e sit close together (both curves converge near N=3,000) — legible
   but tight; could use a leader line if a reviewer flags it.
4. A/B layout choice was a visual call at 100% zoom, not quantitative.

## Files

- `40_figure1_v5.py` — compose script (`vm-panels` / `compose` modes).
- `40a_admixture_bins.py` — VM-side panel (a) binning script.
- `panel_a_admixture_bins.tsv`, `panel_c_novelty_totals.tsv` — the two VM exports this round
  composed against (aggregate-only, disclosure-checked).
- `layout_A.png/.pdf`, `layout_B.png/.pdf` — both variants.
- `figure1_v5.png/.pdf` — chosen layout (A), all six panels populated.
- `compose_summary.json` — parameters used for this compose run.
