# Critic — WS-C figure style pass (S04, 2026-09-25)

Scope note: cut short by coordinator mid-review to hand a layout-fault list to a redesign
agent. This is not a complete pass over every checklist item in `reference/FIGURE_STYLE.md`
for every figure — treat the "not yet checked" note per figure as real gaps, not clean bills.

## Marc's verdict: style pass insufficient
Marc reviewed the re-rendered figures directly and said "I am not seeing nearly any
difference," and flagged overlapping text as clearly inadmissible. Confirmed on inspection:
the rcParam port (fonts/spines/legend frame) is real but subtle at thumbnail size, and it did
not touch panel/subplot *layout* (spacing, alignment, sizing) at all — the actual visible
problem. Full list of overlap/misalignment/layout faults found, by figure:

## reports/hla_popgen/37_dq_g1g2_signed_ld/fig_dq_g1g2_committed_*.png (panel c-style heatmap)
- **Top marginal carrier-frequency bar chart is not column-aligned with the heatmap below it.**
  The bars sit in their own `Axes` sized independently of the main heatmap's column count/width;
  any mismatch in `Axes` width or in the number of x-ticks vs. columns throws off alignment.
  Confirmed misaligned in `fig_dq_g1g2_committed_MAIN_POOLED.png` (and per-ancestry supplements,
  same shared plotting code) — bars visibly do not line up with their DQB1 columns.
- **Right marginal bar chart + colorbar float with large empty gaps** between the heatmap's right
  edge, the marginal-bar axes, and the colorbar axes — `gridspec`/`wspace` between these three
  axes is too large and/or their relative widths are not proportioned to the heatmap.
- Not yet fixed. Root cause is almost certainly in the shared plotting code for this figure family
  (need to locate the script — likely `scripts/hla_popgen/37_dq_g1g2_signed_ld.py`) where the
  marginal-axes `gridspec` is defined; the fix is a `gridspec_kw` width/height-ratio and
  `sharex`/`sharey` correction, not a style/rcParam change.

## reports/hla_popgen/39_saturation_by_ancestry/
- **fig2_supplement_grid_by_gene.png**: gene-name subplot titles (row 2: DQA1/DQB1/DPA1/DPB1)
  overlapped the x-axis tick labels of row 1 (e.g. "1000" collided with "DQA1"). **Fixed this
  session**: added `gridspec_kw={"hspace": 0.55, "wspace": 0.35}` and increased per-row height
  45mm->50mm in `fig_supplement_grid()` (`scripts/hla_popgen/39_saturation_by_ancestry.py`),
  re-rendered, re-inspected — row overlap gone. **New minor overlap introduced**: the
  `fig.supxlabel("cohort size (people)")` now sits very close to / touches the legend row
  immediately below it — not re-fixed, needs more bottom margin or moving the legend up.
- **fig1_main_saturation_panel.png, fig3_thresholds_*.png, fig4_imgt_allele_space_explored.png**:
  visually inspected, no text/label overlaps found. Not exhaustively checked against every
  checklist item (e.g. exact pt sizes not measured with a ruler, only eyeballed).

## reports/hla_popgen/40_figure1_v5/figure1_v5.png
- No overlap found in panels a-f on visual inspection (admixture bars, ternary, heatmap, bar
  chart, saturation curves, cA/cB bar chart). Panel c's colorbar and panel d's colorbar sit
  close together in the shared row but do not visibly overlap text. Not re-checked at full
  native resolution / zoomed crops — only whole-figure inspection.

## reports/hla_popgen/43_kir_full_cohort/fig_kir_full_cohort.png (43b)
- No new text overlap found on whole-figure inspection (panel b's n= labels were already fixed
  under CRITIC_2 in S03; panel f is now points+CI). **Not re-verified at zoomed/native
  resolution** — CRITIC_2 (S03) found overlaps only visible at native res that whole-figure
  inspection missed once before; the same risk applies here and this pass did not zoom in.
- Panel c integer formatting vs. panel d one-decimal formatting is an inconsistency (checklist
  item 11, "consistent number formats") — not a layout/overlap fault but flagged, not fixed.

## Not checked at all this session
- `reports/hla_popgen/37_dq_g1g2_signed_ld/fig_bimodality_committed.png` layout (title is
  sentence-length, checklist item 7 — flagged, not fixed).
- Exact point widths (89/183mm) not re-measured for any figure this session.

## Status
README.md for `39_saturation_by_ancestry` verified byte-identical to the pre-session
hand-written committed version after restoring content `render_only()` silently drops (see
commit). `format_p()` bug fixed and confirmed end-to-end. `pytest scripts/hla_popgen/tests -q`:
**444 passed, 9 errors** (errors are pre-existing missing `fixtures_dir` fixture / unrelated
setup issues in `test_allele_geometry.py`, `test_extraction.py`, `test_figures.py`,
`test_novel.py` — reproduced on a clean environment, not caused by this session's changes).

**Handoff to redesign agent**: fix 37c's marginal-axes alignment/gridspec (highest priority —
Marc called this out by name), then zoom-inspect 43b and 40 at native resolution for
sub-pixel-visible overlaps neither this pass nor CRITIC_2 fully ruled out.
