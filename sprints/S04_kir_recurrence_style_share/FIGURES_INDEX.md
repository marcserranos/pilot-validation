# S04 figure index

| Figure | Path | Caption |
|---|---|---|
| Figure 1 v5 (a-f) | `reports/hla_popgen/40_figure1_v5/figure1_v5.png` (+ `.pdf`) | Admixture, HLA-B ternary novelty, per-gene x ancestry novelty heatmap, novel-protein-allele counts, per-ancestry saturation curves, G1xG2 O/E purge. **S04 WS-C phase 2 redesign (2026-09-25)**: fixed all 5 `check_layout(strict=True)` errors the diagnostic pass found -- panel b's redundant 0/100 ternary tick numerals (which sat on top of the AFR/AMR vertex labels) dropped to a single midpoint "50" per edge; panel b's marker legend moved from outside the axes (colliding with panel c's y-label) to below the triangle's own base edge, where the padding band is data-free by construction; panel c's two-line sentence-length y-axis label ("% called haplotypes with sequence (any field) absent from IPD-IMGT/HLA -- mostly non-coding") shortened to a direct one-line label ("any-field novelty (%)", full metric definition kept in the report README) -- the second line was what gave the rotated label the extra width that reached into panel b. No panel-content/number changes; `check_layout(strict=True)` now passes with 0 errors for both composed layouts (A and B). |
| DQ G1/G2 pooled | `reports/hla_popgen/37_dq_g1g2_signed_ld/fig_dq_g1g2_committed_MAIN_POOLED.png` (+ `.pdf`) | Signed phased D' between DQA1 and DQB1 alleles, pooled across ancestries; G1/G2 haplotype-group brackets adjacent to tick labels. **S04 WS-C phase 2 redesign** (layout linter + redesign, 2026-09-25): one GridSpec with sharex/sharey ties the top/right marginal carrier-frequency bars to the heatmap's own columns/rows (previously misaligned); colourbar compact in its own column beside the right marginal (previously floating with large gaps); short bold panel label ("DQ G1/G2, pooled") replaces the sentence title; the two "predicted incompatible" quadrants are now a flat neutral-grey fill with one annotation ("0/469 observed") instead of a solid dark-blue block duplicated twice. **Redesign pass 2, same day**: Marc rejected pass 1 on sight -- the "DQA1"/"DQB1" axis labels visibly ran through the G1/G2 brackets, a fault `check_layout()` had missed entirely (it only ever compared Text against Text; the brackets are Line2D, invisible to it). Fixed the linter itself (new `mark_decoration()`-based checks, see below) AND the figure: dropped the now-redundant axis labels, closed the dead band between the corner label and the top marginal, switched marginal ticks to round numbers (`_nice_ceiling()`), and faded (translucent white overlay) compatible-quadrant cells whose independence-expectation count is <5 haplotypes so a confidently-estimated D'=-1/+1 reads more salient than a thin-evidence one (justified in the report README). Passes `_viz_common.check_layout()` with zero violations (`strict=True`, the `save_fig()` default) including the two new checks. |
| DQ G1/G2 per-ancestry supplements | `reports/hla_popgen/37_dq_g1g2_signed_ld/fig_dq_g1g2_committed_supp_{AFR,AMR,EAS,EUR,SAS}.png` | Same redesign (both passes) as the pooled panel, split by ancestry; each panel's incompatible-quadrant annotation gives that ancestry's own "n/N observed" count from `oe_purge_committed.tsv`. |
| DQ bimodality | `reports/hla_popgen/37_dq_g1g2_signed_ld/fig_bimodality_committed.png` | Histogram of signed D' within G1/G1 and G2/G2 compatible pairs vs predicted-incompatible pairs. **S04 WS-C phase 2**: sentence title shortened to "DQ G1/G2 bimodality" (ancestry list moved to the report README); x-axis ticks fixed to [-1, 0, 1] (was matplotlib's auto ~7-tick scale, which overlapped at this panel's narrow width -- caught by the layout linter). |
| KIR saturation, main panel (a-c) | `reports/hla_popgen/39_saturation_by_ancestry/fig1_main_saturation_panel.png` | Per-ancestry discovery curves: cumulative distinct protein alleles, alleles with >=2 carriers, novel-protein alleles vs cohort size. **S04 WS-C redesign (2026-09-25, this pass)**: fixed the diagnostic pass's 5 `check_layout(strict=True)` errors (panel letters a/b/c overlapping their nearest y-tick label -- a real, mechanically-clipped-in-the-unclipped-check phantom tick matplotlib generates just past each axis's view limit; see `_prune_offview_ticklabels()` in the script). Also: `render_only()` no longer regenerates `README.md` from a template (was silently dropping hand-written sections, caught twice by the critic -- it now writes figures only); axes made honest (start at 0, no negative padding on a count that can't be negative); added a dashed equal-N reference line + label (`N=487`) on panels a/b so the equal-N comparison point the README's headline claim rests on ("AFR has the steepest equal-N slope") is visible directly on the figure, not only in a table. |
| KIR saturation supplement grid | `reports/hla_popgen/39_saturation_by_ancestry/fig2_supplement_grid_by_gene.png` | Per-gene saturation curves, all ancestries. **S04 WS-C redesign**: fixed the diagnostic pass's 13 errors -- per-subplot phantom-tick clipping (same `_prune_offview_ticklabels()` fix) plus the `fig.supxlabel()`/`fig.legend()` pair sitting close enough to touch; replaced with explicit `fig.text()` y-positions stacking plots -> x-axis label -> legend with a fixed, clearly separated margin. |
| KIR saturation per-ancestry threshold panels | `reports/hla_popgen/39_saturation_by_ancestry/fig3_thresholds_{AFR,AMR,EAS,EUR,MID,SAS}.png` | Recurrence-class saturation (singleton / seen-twice / seen>=N) per ancestry. **S04 WS-C redesign**: fixed the diagnostic pass's 1 error (EAS panel's phantom y-tick clipping) plus widened the right margin so the in-line ">= k carriers" end labels never sit at/past the right spine. |
| KIR / IPD-KIR allele-space coverage | `reports/hla_popgen/39_saturation_by_ancestry/fig4_imgt_allele_space_explored.png` | Fraction of the IPD-KIR catalogue's alleles observed in this cohort. Re-rendered under the updated style only (was already clean, 0 errors). |
| KIR full-cohort, multi-panel (a-f) | `reports/hla_popgen/43_kir_full_cohort/fig_kir_full_cohort.png` (+ `.pdf`) | Gene presence, novelty composition, presence-deviation-from-pooled x ancestry, protein-novelty x ancestry, cA/cB content, QC framework-gene presence. WS-C fixes: (1) panel f converted from a truncated bar chart (baseline 86%, exaggerating differences by area) to a zero-free point + 95% CI plot (CRITIC_2 known debt, resolved); (2) panel b/c gridspec `wspace` widened 0.6->0.75 for margin. **S04 WS-C redesign pass 1**: fixed the diagnostic pass's 1 error (panel f's co-occurrence text clipped off the figure's right edge -- moved into an in-figure caption); removed the "†" marker for bold gene labels + existing darker-bar color coding; per-panel heatmap number-format precision documented. **Redesign pass 2 (orchestrator full-size review, same day)**: pass 1 was rejected on sight -- "not paper-ready yet". Fixed for real: (a) panel b's `n=` labels were STILL sitting on top of the bars (right-aligned-at-fixed-x means a long string like "n=23,142" grows leftward onto the data) -- rebuilt as a genuine left-aligned right-hand column past a thin rule, which can only grow away from the bars regardless of digit count; (b) panel c redesigned from a flat presence-per-ancestry heatmap (redundant with panel a) to each ancestry's DEVIATION from panel a's pooled value, diverging colormap centered at 0 (`_viz_common.diverging_cmap/_norm`) -- more informative and original, per review; (c) all heatmap cell text now gets its color from the CELL'S OWN rendered luminance (dark on light, white on dark), not a linear vmin/vmax-midpoint guess, fixing unreadable dark-on-dark numbers on viridis/OrRd's darkest cells (e.g. AFR/3DS1 "8", AFR/2DL5B "38.9"); panel d also switched to integer %; (d) the co-occurrence text and the long disclosure footnote were removed from the figure entirely (they now live only in the README/this caption -- an in-figure caption baked into the raster was itself flagged as against house style); (e) the novelty-tier legend now anchors directly under panel b's own axes (figure-fraction placement pinned to panel b's horizontal center) instead of `fig.legend()` centered across all four columns, closing the dead band between the rows; (f) panel e's unexplained dashed 50% reference line removed. `check_layout(strict=True)`: 0 errors both passes. |
| Deletion-rate supplement (a-b) | `reports/hla_popgen/38_hla_a_deletion_validation/supp_deletions.png` (+ `.pdf`) | Panel a: bridged-deletion % per gene (Wilson 95% CI), positive/negative controls marked by bracket + label, not colour. Panel b: CNV-gene (DRB3/4/5, C4A/C4B) deletion frequency by ancestry, dot-plot with Wilson 95% CI. **S04 WS-C redesign (this pass)**: fixed a `check_layout(strict=True)` failure on clipped tick labels -- same root cause as 39's fix, a phantom off-view tick (`'60'` on panel b's x-axis, `'0'`/`'>500'`-scale ticks elsewhere) that already renders invisible (clipped by the axes' own bbox) but that `check_layout()`'s unclipped `get_window_extent()` read as extending past the whole figure; fixed with the same `_prune_offview_ticklabels()` helper (careful with panel b's `invert_yaxis()` -- an inverted axis returns `get_xlim()`/`get_ylim()` reversed, `hi < lo`, which the first version of this helper didn't account for and briefly blanked panel b's own gene-name y-tick labels; fixed with `sorted()`). Also fixed a real, separate bug found while re-running the production path: `write_readme()` used to write straight to the SHARED `README.md` in this same directory (script 38's own hand-written cross-validation writeup, into which a human had previously hand-merged a condensed version of 38b's content) -- unconditionally overwriting it destroyed that hand-written file on every re-run. Now writes a standalone `README_38b_supplement.md` instead; `README.md` is left untouched (merging any updated numbers into it remains a deliberate human step). |

Before/after PNGs for all of the above: `sprints/S04_kir_recurrence_style_share/style_before/` and
`style_after/` (same filenames).

## Layout linter (S04 WS-C phase 2, 2026-09-25)

`scripts/hla_popgen/_viz_common.py` gained `check_layout(fig)`, `mark_marginal(ax_margin,
ax_main, axis=)`, and (2026-09-25, phase 2) `mark_decoration(artist)`. `check_layout` catches, on a
fully-drawn figure: (a) pairwise overlaps among all visible Text artists (titles, axis/tick labels,
legend text, annotations); (a2) any Text overlapping a Line2D/Patch explicitly registered via
`mark_decoration()` as a structural bracket/divider/connector -- added after (a) alone missed the
rotated "DQA1" ylabel drawn straight through the DQ G1/G2 bracket's Line2D stem (a Line2D was never
a Text artist, so it was structurally invisible to (a); see `mark_decoration()`'s docstring for the
full writeup); (a3) any Text overlapping a visible spine of an Axes it doesn't belong to (own-axes
spines exempt); (b) text or Axes extending far beyond the figure's own bbox (tolerant of the few-px
overflow `bbox_inches="tight"` normally absorbs, but not a large layout-breaking offset); (c) for
any axes registered via `mark_marginal()`, its data-to-display mapping must coincide with its main
axes' mapping within 0.5px -- catches a marginal whose bars don't actually line up with the main
panel's columns/rows even when nothing overlaps; (d) excess whitespace (warning only). `save_fig()`
calls `check_layout()` and raises `RuntimeError` on any error-severity violation by default
(`strict=True`); `strict=False` requires a mandatory `reason=` string. Unit tests:
`scripts/hla_popgen/tests/test_layout_linter.py` (14 tests: clean figure passes, each fault type
planted and caught including the two new (a2)/(a3) checks and two direct regressions using the
real `_bracket_v`/`_bracket_h` helpers from `37c_dq_g1g2_from_committed.py`, `save_fig`
strict/override behavior).

Report-only `check_layout()` pass (figures re-drawn in-process from their own committed input
TSVs via each script's own `compose`/`render_only` entry point, `save_fig()` monkey-patched to
call `check_layout()` and report instead of writing/raising -- current committed PNGs on disk are
NOT overwritten by this pass) against other key S04/S03 figures already in this index, to gauge
how much of the earlier "no overlap found on visual inspection" critic passes hold up against a
mechanical check:

| Figure | Errors | Warnings | Notes |
|---|---|---|---|
| `40_figure1_v5.py compose` layout A (= committed `figure1_v5.png`, `--pick A` default) | **0** (was 5, fixed 2026-09-25 -- see FIGURES_INDEX row above and `40_figure1_v5.py`) | 0 | Fixed: `100`/`0` tick labels overlapped `AFR`/`AMR` ancestry labels (dropped the redundant 0/100 numerals); `EUR` ancestry label and the "first observed here" annotation both overlapped panel c's two-line y-axis label (legend moved below the triangle, y-label shortened to one line). Re-run with `--panel-b-table reports/hla_popgen/36_figure1/panel_b_ternary_alleles.tsv` (40's own default path's file isn't committed locally; 36's is numerically the same panel-b data per 40's own docstring lineage) -- caveat: not verified byte-identical to whatever panel-b table produced the currently-committed PNG. |
| `40_figure1_v5.py compose` layout B (alternate, not picked) | **0** (was 5) | 0 | Same fault pattern as layout A, same fix. |
| `39_saturation_by_ancestry.py --render-only` fig1 (`fig1_main_saturation_panel.png`) | **0** (was 5, fixed 2026-09-25) | 0 | Was: panel letters a/b/c overlap the nearest y-axis tick label in each of the 3 stacked panels (e.g. `'a' overlaps '600'`). Root cause (confirmed by direct inspection, not guessed): matplotlib's tick Locator always generates one extra major tick just past each end of the actual view (e.g. ylim=(-16, 536) but a `'600'` y-tick Text object still exists, `clip_on=True` with the axes bbox as clip path) -- these already render invisible in the saved PNG/PDF, but `check_layout()`'s overlap check calls `Text.get_window_extent()` directly, which ignores clipping. Fixed in the script (not `_viz_common.py`, per this task's constraint) with a local `_prune_offview_ticklabels(ax)` that explicitly hides any tick Text whose data position is outside the axes' own final xlim/ylim -- a no-op for the rendered pixels, but now agrees with what `check_layout()` sees. CRITIC_WSC.md's "no overlap found" note for this figure was an eyeball miss at thumbnail size; the underlying fault was real (mis-clipped-tick bookkeeping) even though no reader would ever have seen it. |
| `39_saturation_by_ancestry.py --render-only` fig2 (`fig2_supplement_grid_by_gene.png`) | **0** (was 13, fixed 2026-09-25) | 0 | Was: repeated x-tick-label collisions (`'3000' overlaps '-1000'`, `'4000' overlaps '-5'/'-10'/'-20'`) across several of the 8 per-gene subplots (same phantom-tick cause as fig1, fixed the same way), plus `fig.supxlabel("cohort size (people)")` overlapping the ancestry legend below it -- this second fault WAS real and visible (CRITIC_WSC.md flagged it as "not re-fixed"); fixed by replacing the auto-placed `supxlabel`/`legend(bbox_to_anchor=negative)` pair with explicit `fig.text()` y-positions, stacking plots -> x-axis label -> legend with a fixed, clearly separated margin. |
| `39_saturation_by_ancestry.py --render-only` fig3 thresholds (AFR/AMR/EAS/EUR/MID/SAS) | 0/0/**0**/0/0/0 (EAS was 1, fixed 2026-09-25) | 0 | Was: EAS panel `'-200' overlaps '-50'` (same phantom off-view y-tick cause, same `_prune_offview_ticklabels()` fix); also widened the right margin so the in-line ">= k carriers" end labels never sit at/past the right spine (the actual EAS instance the critic's visual pass caught, distinct from the linter's own tick finding). |
| `39_saturation_by_ancestry.py --render-only` fig4 (`fig4_imgt_allele_space_explored.png`) | 0 | 0 | Clean; re-rendered with `_prune_offview_ticklabels()` + honest zero-based y-axis anyway for consistency with the other panels in this figure family. |
| `43b_kir_full_figure.py` (`fig_kir_full_cohort.png`) | **0** (was 1, fixed 2026-09-25) | 0 | Was: the panel f co-occurrence annotation text (`"2DL2 & 2DL3 co-occur..."` two-line block) extends past the figure's right edge -- a real, previously-unflagged clipping fault (CRITIC_WSC.md's note for 43b was "not re-verified at zoomed/native resolution"). Fixed by moving the co-occurrence numbers out of the axes-relative annotation (which had nothing to its right to expand into -- this is the rightmost panel) and into the figure-wide caption, centered and already sized to the full figure width. |
| `38b_deletion_supplement_fig.py` (`supp_deletions.png`) | **0** (was failing `test_deletion_supplement_fig.py`'s strict `save_fig()`, fixed 2026-09-25) | 0 | Was: a phantom off-view x-tick (`'60'`, panel b) read as extending past the figure's right edge, plus a phantom `'0'` off the left edge -- same root cause and fix (`_prune_offview_ticklabels()`) as 39/43b above, adapted for panel b's `invert_yaxis()` (which reverses `get_ylim()`'s order; the naive first version of the fix blanked panel b's own gene-name labels before this was caught by re-reading the rendered PNG, not just trusting the passing linter). |

Only 38b was explicitly in scope for a linter-failure fix per the task brief; 39/43b were listed as
diagnostic-only findings in the prior session but are now fixed too, since leaving 3 committed
figures at 0/13/5/1 known `check_layout(strict=True)` errors while claiming the redesign pass
"passes the strict linter" would not have been honest. WS-A/WS-D's own scripts (44/45/47/48) are
out of scope and untouched; `test_44_kir_recurrence_saturation.py`'s one failing test
(`45_kir_recurrence_figure.py`'s `fig_coverage_chao2`, 3 clipping errors) is a pre-existing failure
from the same linter-addition commit, confirmed via `git stash` to fail identically before this
session's changes -- flagged here, not fixed (not one of this session's three assigned figures).
Net finding: the mechanical linter surfaces several real layout faults that two rounds of
manual/critic visual inspection missed (fig1/fig2/43b/38b), and also correctly reproduces the one
fault CRITIC_WSC.md already knew about and left unfixed (fig2's supxlabel/legend gap) -- i.e. it
agrees with human review where human review was right, and catches what it missed. It also has a
real false-positive mode of its own (see "Linter false positives found" below), distinct from the
false negative the task brief already knew about (axis labels vs. tick labels/brackets).

## Linter false positives found (S04 WS-C, this pass, 2026-09-25)

`check_layout()`'s clip/overlap checks call `Text.get_window_extent()` directly, which returns a
Text artist's UNCLIPPED geometry. matplotlib's default tick Locator deliberately generates one
extra major tick just past each end of an axis's actual view (confirmed empirically, see the fig1
row above) -- those tick Text objects exist and have `get_visible() == True`, but are never
actually drawn on screen because they have `clip_on=True` with the axes' own bbox as their clip
path. `check_layout()` doesn't know about clipping, so it can flag two such phantom, nobody-ever-
sees-them ticks as colliding with a panel letter, with each other across adjacent subplots, or
(when the axis view is much narrower than one tick step, as in 38b's panel b) as extending far
past the whole figure bbox -- all without a single visible pixel being wrong. Worked around at the
script level in 39/43b/38b with a small `_prune_offview_ticklabels(ax)` helper (hides any tick
Text whose data position is outside the axes' own final xlim/ylim -- provably a no-op for the
rendered raster, confirmed by diffing rendered PNGs before/after). Not fixed in `_viz_common.py`
itself per this task's explicit instruction not to edit that file while another agent owns it;
flagged here and in `LOG.md` as a suggested enhancement for whoever next touches `check_layout()`
(e.g. clip each Text's extent to its own `get_clip_box()` before the overlap/bbox checks, or skip
tick Texts whose own data position already falls outside `ax.get_xlim()`/`get_ylim()`).

## Proposed new lint rule: text-vs-bar/patch overlap for labels (orchestrator, 2026-09-25)

Orchestrator's full-size review of `fig_kir_full_cohort.png` caught a fault none of the mechanical
checks above (nor two rounds of visual review) had: panel b's `n=23,142`-style labels sitting
directly ON TOP of the orange/blue bar ends -- text over data, "clearly inadmissible" (same
standard as the original overlapping-text rejection). `check_layout()`'s existing checks don't
cover this: check (a) is Text-vs-Text only; `mark_decoration()`/(a2) is deliberately for
non-data decorations (brackets/dividers) ONLY -- its own docstring explicitly excludes bar/patch
DATA, since a bar's own in-bar value label legitimately sits on top of its bar. What's missing is
a THIRD category: a text label that sits *near* data (not on top of its own bar, but drifted onto
a *different* bar/patch than the one it annotates) has no opt-in mechanism at all.

Proposed rule (not implemented -- `_viz_common.py` is owned by another agent this session, per
instruction): a `mark_label(text_artist)` opt-in (mirroring `mark_decoration()`'s pattern) for a
text artist that is a per-item VALUE LABEL for a specific bar/point/cell (e.g. panel b's `n=`
counts, panel f's would-be co-occurrence annotations) rather than an axis/tick/legend label. A new
check (a4) in `check_layout()` would test each `mark_label()`-registered Text's bbox against every
Patch/Rectangle in the same Axes EXCEPT the one item it's paired with (the pairing passed
explicitly, e.g. `mark_label(text, paired_patch=bar_patch)`), flagging overlap with any other bar
as an error. This is deliberately opt-in and pairs the label with its own patch, so it never
flags the common, legitimate case (a label inside/at the end of its OWN bar) while still catching
this exact fault (a label drifting onto a bar it does NOT belong to, i.e. text sitting on data it
isn't labeling). Root cause in this case: right-aligning a variable-length string at a fixed
anchor x lets its rendered width extend backward over the neighboring bar; the actual fix applied
here was structural (a left-aligned, fixed-start label column, so a label can only grow away from
the data), but the lint rule would have caught the original fault mechanically instead of needing
a human to view the PNG at full size.
