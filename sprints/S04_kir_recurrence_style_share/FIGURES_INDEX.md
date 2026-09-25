# S04 figure index

| Figure | Path | Caption |
|---|---|---|
| Figure 1 v5 (a-f) | `reports/hla_popgen/40_figure1_v5/figure1_v5.png` (+ `.pdf`) | Admixture, HLA-B ternary novelty, per-gene x ancestry novelty heatmap, novel-protein-allele counts, per-ancestry saturation curves, G1xG2 O/E purge. Re-rendered under WS-C style pass (cnsplots rcParams); no panel-content changes. Panel b (strict ancestry >=0.98) vs panel c/d (>=0.9) threshold mismatch is real, baked into VM-exported aggregates, documented as an open item in the report README -- not fixed here (would need a VM rerun). |
| DQ G1/G2 pooled | `reports/hla_popgen/37_dq_g1g2_signed_ld/fig_dq_g1g2_committed_MAIN_POOLED.png` (+ `.pdf`) | Signed phased D' between DQA1 and DQB1 alleles, pooled across ancestries; G1/G2 haplotype-group brackets adjacent to tick labels. **S04 WS-C phase 2 redesign** (layout linter + redesign, 2026-09-25): one GridSpec with sharex/sharey ties the top/right marginal carrier-frequency bars to the heatmap's own columns/rows (previously misaligned); colourbar compact in its own column beside the right marginal (previously floating with large gaps); short bold panel label ("DQ G1/G2, pooled") replaces the sentence title; the two "predicted incompatible" quadrants are now a flat neutral-grey fill with one annotation ("0/469 observed") instead of a solid dark-blue block duplicated twice. Passes `_viz_common.check_layout()` with zero violations (`strict=True`, the `save_fig()` default). |
| DQ G1/G2 per-ancestry supplements | `reports/hla_popgen/37_dq_g1g2_signed_ld/fig_dq_g1g2_committed_supp_{AFR,AMR,EAS,EUR,SAS}.png` | Same redesign as the pooled panel, split by ancestry; each panel's incompatible-quadrant annotation gives that ancestry's own "n/N observed" count from `oe_purge_committed.tsv`. |
| DQ bimodality | `reports/hla_popgen/37_dq_g1g2_signed_ld/fig_bimodality_committed.png` | Histogram of signed D' within G1/G1 and G2/G2 compatible pairs vs predicted-incompatible pairs. **S04 WS-C phase 2**: sentence title shortened to "DQ G1/G2 bimodality" (ancestry list moved to the report README); x-axis ticks fixed to [-1, 0, 1] (was matplotlib's auto ~7-tick scale, which overlapped at this panel's narrow width -- caught by the layout linter). |
| KIR saturation, main panel (a-c) | `reports/hla_popgen/39_saturation_by_ancestry/fig1_main_saturation_panel.png` | Per-ancestry discovery curves: cumulative distinct protein alleles, alleles with >=2 carriers, novel-protein alleles vs cohort size. Already compliant (direct line labels, no legend, honest axes); re-rendered under the updated style only. |
| KIR saturation supplement grid | `reports/hla_popgen/39_saturation_by_ancestry/fig2_supplement_grid_by_gene.png` | Per-gene saturation curves, all ancestries. Re-rendered under the updated style only. |
| KIR saturation per-ancestry threshold panels | `reports/hla_popgen/39_saturation_by_ancestry/fig3_thresholds_{AFR,AMR,EAS,EUR,MID,SAS}.png` | Recurrence-class saturation (singleton / seen-twice / seen>=N) per ancestry. Re-rendered under the updated style only. |
| KIR / IPD-KIR allele-space coverage | `reports/hla_popgen/39_saturation_by_ancestry/fig4_imgt_allele_space_explored.png` | Fraction of the IPD-KIR catalogue's alleles observed in this cohort. Re-rendered under the updated style only. |
| KIR full-cohort, multi-panel (a-f) | `reports/hla_popgen/43_kir_full_cohort/fig_kir_full_cohort.png` (+ `.pdf`) | Gene presence, novelty composition, presence x ancestry, protein-novelty x ancestry, cA/cB content, QC framework-gene presence + co-occurrence. WS-C fixes: (1) panel f converted from a truncated bar chart (baseline 86%, exaggerating differences by area) to a zero-free point + 95% CI plot -- point position carries no baseline claim, so it can honestly show the informative 93-97% range (CRITIC_2 known debt, resolved); (2) panel b's `n=` calls-count labels right-aligned + clipped instead of left-aligned overflow, fixing the "n=23,097" / "3DP1 †" text collision with panel c's row labels (CRITIC_2 known debt, resolved); (3) panel b/c gridspec `wspace` widened 0.6->0.75 for margin. |

Before/after PNGs for all of the above: `sprints/S04_kir_recurrence_style_share/style_before/` and
`style_after/` (same filenames).

## Layout linter (S04 WS-C phase 2, 2026-09-25)

`scripts/hla_popgen/_viz_common.py` gained `check_layout(fig)` and `mark_marginal(ax_margin,
ax_main, axis=)`. `check_layout` catches, on a fully-drawn figure: (a) pairwise overlaps among all
visible Text artists (titles, axis/tick labels, legend text, annotations); (b) text or Axes
extending far beyond the figure's own bbox (tolerant of the few-px overflow `bbox_inches="tight"`
normally absorbs, but not a large layout-breaking offset); (c) for any axes registered via
`mark_marginal()`, its data-to-display mapping must coincide with its main axes' mapping within
0.5px -- catches a marginal whose bars don't actually line up with the main panel's
columns/rows even when nothing overlaps; (d) excess whitespace (warning only). `save_fig()` now
calls `check_layout()` and raises `RuntimeError` on any error-severity violation by default
(`strict=True`); `strict=False` requires a mandatory `reason=` string. Unit tests:
`scripts/hla_popgen/tests/test_layout_linter.py` (8 tests: clean figure passes, each fault type
planted and caught, `save_fig` strict/override behavior).

Report-only `check_layout()` pass (figures re-drawn in-process from their own committed input
TSVs via each script's own `compose`/`render_only` entry point, `save_fig()` monkey-patched to
call `check_layout()` and report instead of writing/raising -- current committed PNGs on disk are
NOT overwritten by this pass) against other key S04/S03 figures already in this index, to gauge
how much of the earlier "no overlap found on visual inspection" critic passes hold up against a
mechanical check:

| Figure | Errors | Warnings | Notes |
|---|---|---|---|
| `40_figure1_v5.py compose` layout A (= committed `figure1_v5.png`, `--pick A` default) | **5** | 0 | `100`/`0` tick labels overlap `AFR`/`AMR` ancestry labels; `EUR` ancestry label and the "first observed here" annotation both overlap panel d's two-line y-axis label. Re-run with `--panel-b-table reports/hla_popgen/36_figure1/panel_b_ternary_alleles.tsv` (40's own default path's file isn't committed locally; 36's is numerically the same panel-b data per 40's own docstring lineage) -- caveat: not verified byte-identical to whatever panel-b table produced the currently-committed PNG. |
| `40_figure1_v5.py compose` layout B (alternate, not picked) | 5 | 0 | Same fault pattern as layout A. |
| `39_saturation_by_ancestry.py --render-only` fig1 (`fig1_main_saturation_panel.png`) | **5** | 0 | Panel letters a/b/c overlap the nearest y-axis tick label in each of the 3 stacked panels (e.g. `'a' overlaps '600'`) -- CRITIC_WSC.md's "no overlap found" note for this figure was an eyeball miss at thumbnail size, not a false negative in the fix. |
| `39_saturation_by_ancestry.py --render-only` fig2 (`fig2_supplement_grid_by_gene.png`) | **13** | 0 | Repeated x-tick-label collisions (`'3000' overlaps '-1000'`, `'4000' overlaps '-5'/'-10'/'-20'`) across several of the 8 per-gene subplots, plus `fig.supxlabel("cohort size (people)")` overlapping the EAS/EUR ancestry legend labels below it -- this is exactly the still-open gap CRITIC_WSC.md flagged as "not re-fixed" in the prior session; the linter confirms it mechanically rather than relying on another eyeball pass. |
| `39_saturation_by_ancestry.py --render-only` fig3 thresholds (AFR/AMR/EAS/EUR/MID/SAS) | 0/0/**1**/0/0/0 | 0 | EAS panel: `'-200' overlaps '-50'` (y-axis tick crowding); the other 5 ancestries are clean. |
| `39_saturation_by_ancestry.py --render-only` fig4 (`fig4_imgt_allele_space_explored.png`) | 0 | 0 | Clean. |
| `43b_kir_full_figure.py` (`fig_kir_full_cohort.png`) | **1** | 0 | The panel f co-occurrence annotation text (`"2DL2 & 2DL3 co-occur..."` two-line block) extends past the figure's right edge -- a real, previously-unflagged clipping fault (CRITIC_WSC.md's note for 43b was "not re-verified at zoomed/native resolution"; this is exactly the kind of fault that pass admitted it might miss). |

None of these were fixed in this session (out of scope -- WS-A/WS-D own scripts 44/45/47/48; this
pass is diagnostic only, per the task brief). Net finding: the mechanical linter surfaces several
real layout faults that two rounds of manual/critic visual inspection missed (fig1/fig2/43b), and
also correctly reproduces the one fault CRITIC_WSC.md already knew about and left unfixed (fig2's
supxlabel/legend gap) -- i.e. it agrees with human review where human review was right, and catches
what it missed.
