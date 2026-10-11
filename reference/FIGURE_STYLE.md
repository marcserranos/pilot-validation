# Figure style reference: cnsplots

Source: [github.com/faridrashidi/cnsplots](https://github.com/faridrashidi/cnsplots), commit
`634482c304cbbdc9e1b5a56e9acd1f64aeb1b4a6` (2026-09-13, tag/version `0.7.0`). License:
**BSD-3-Clause** (Farid Rashidi, 2023-2026) — permissive, attribution required for redistributed
code; exact rcParam values and hex codes are facts/data, not copyrightable expression, but we cite
the source below anyway. Install: `pip install cnsplots` (PyPI, Python >=3.10,<4.0; pulls in
matplotlib>=3.10, seaborn>=0.13.2, scanpy, gseapy, lifelines, statsmodels, scikit-learn,
pycomplexheatmap, and ~10 more — it is a large, opinionated dependency, not a thin style module).

Supervisor sent this as the guideline for why our figures look "AI-like." It is a matplotlib
styling + plotting library aimed at Cell/Nature/Science conventions; we already have our own
`_viz_common.py: nature_style()` doing the same job for our 05/06/07/36+ figures. The
recommendations below are extracted from its `_settings.py` (rcParam-equivalent defaults),
`_palettes.py` (journal color sets), and `_sizing.py` / docs (physical-unit figure sizing).

## Concrete rules (quoted defaults from `_settings.py`)

- **Font**: `font_family="sans-serif"`, `font_sans_serif=("Helvetica", "Helvetica Neue", "Arial",
  "Nimbus Sans", "Liberation Sans", "DejaVu Sans")`. `mathtext_fontset="custom"`.
- **Text sizes**: `title_fontsize=8` (titles + axis labels), `title_fontweight="bold"`,
  `legend_fontsize=7` (legend/tick/colorbar text), `pvalue_fontsize="small"`.
- **Spines/ticks**: `axes_linewidth=0.5`, `axes_spines_top=False`, `axes_spines_right=False`,
  `axes_edgecolor="black"`, `axes_grid=False` (no gridlines, ever). Tick marks: `xtick_major_size=2`,
  `xtick_major_width=0.6`, `xtick_major_pad=1` (same for y, `ytick_alignment="center_baseline"`).
- **Legend**: `legend_frameon=False` (no box), `legend_markerscale=0.5`, `legend_handlelength=0.7`,
  `legend_handletextpad=0.3` — small, unobtrusive legend glyphs; a `take_legend_out()` helper moves
  legends outside the axes rather than overlapping data.
- **Export**: `savefig_dpi=288` (72×4), `savefig_bbox="tight"`, `savefig_transparent=True`,
  `pdf_fonttype=42` (embedded TrueType, matches our `save_fig()`), `svg_fonttype="none"` (SVG text
  stays as editable text, not outlined paths — "Adobe Illustrator compatible" claim).
- **Figure sizing is unit-aware, not hardcoded**: `figure()`/`multipanel()` take `unit="pt"|"in"|"mm"`
  (default points, 1/72 in); `multipanel_max_width=540pt` (~190mm, close to our
  `NATURE_DOUBLE_COL_MM=183`). Panels are placed with explicit `width=`/`height=` per panel, not
  derived from data.
- **Panel labels**: `add_panel_label()`, bold, lowercase-by-convention like ours; font weight
  `panel_label_fontweight="bold"`.
- **Significance annotations**: `pvalue_format="star"` (stars, not raw p on the plot by default),
  `pvalue_loc="inside"`.
- **Color-by-luminance contrast**: `annotation_auto_contrast=True` — in-plot text (e.g. heatmap
  cell labels) auto-switches white/black by background luminance rather than a single fixed color.

## Palettes (hex, from `_palettes.py`)

| Palette | Hex (first colors) |
|---|---|
| `Nature` (npg) | `#E64B35 #4DBBD5 #00A087 #3C5488 #F39B7F #8491B4 #91D1C2 #DC0000 #7E6148 #B09C85` |
| `Science` (aaas) | `#3B4992 #EE0000 #008B45 #631879 #008280 #BB0021 #5F559B #A20056 #808180 #1B1919` |
| `Lancet` | `#00468B #ED0000 #42B540 #0099B4 #925E9F #FDAF91 #AD002A #ADB6B6 #1B1919` |
| `NEJM` | `#BC3C29 #0072B5 #E18727 #20854E #7876B1 #6F99AD #FFDC91 ...` |
| `Cell` | `#C84C3A #2F7E8F #E1A22E #4E5A8A #5F9862 #D07A6A #8B6FA8 #7B8C9E #B85F7A #6B6B6B` |
| `OkabeIto` | colorblind-safe qualitative (same family we already use for `ANCESTRY_COLORS`) |

Default qualitative palette is `Ecotyper1` (not `Nature`) — pick explicitly per figure, don't rely
on the library default. Note our `ANCESTRY_COLORS` (Okabe-Ito) is already the same colorblind-safe
choice cnsplots offers as `OkabeIto` — no change needed there.

## How to apply in this repo: **port, don't depend**

**Recommendation: port the rcParams + palette hex values into `_viz_common.py`, do not add
`cnsplots` as a dependency.** Reasons: (1) cnsplots pulls ~15 heavy transitive packages
(scanpy, gseapy, lifelines, pycomplexheatmap...) for what we need in one function,
`nature_style()`, which already exists and is 90% aligned; (2) `pixi.toml` is explicitly out of
scope to touch for this task; (3) our pipeline never plots via cnsplots' own plot functions
(`boxplot`, `heatmapplot`, etc.) — we'd only ever use its style/palette layer, so the dependency
weight buys us nothing; (4) BSD-3-Clause permits copying small config values with attribution,
which is what we're doing.

Concrete port into `scripts/hla_popgen/_viz_common.py`:
- `NATURE_RC`: add `"legend.frameon": False`, `"legend.fontsize": 7` (already close),
  `"axes.grid": False` (already implicit — make explicit), bump `xtick.major.size`/`ytick.major.size`
  to `2`, add `"mathtext.fontset": "custom"`, extend `font.sans-serif` with `"Helvetica Neue"`,
  `"Nimbus Sans"`, `"Liberation Sans"` for cross-platform fallback parity.
- Add a `JOURNAL_PALETTES` dict (hex lists above, attributed "via cnsplots / palettable, BSD-3")
  next to `ANCESTRY_COLORS`, for any figure needing an alternate qualitative palette (keep
  `ANCESTRY_COLORS` — Okabe-Ito — as the one true ancestry palette; don't reuse `Nature`/`Cell` for
  ancestry, since that would fight the existing convention).
- Keep `save_fig()`, `mm()`, `panel_letter()`, `SUPPRESSED_COLOR` / `hatch_suppressed()` as-is —
  they already do what cnsplots' `savefig`/`add_panel_label`/`multipanel` do, at a fraction of the
  dependency cost.

## De-AI checklist for reviewers

From cnsplots' encoded conventions:
1. No axes.grid — are gridlines off? *(cnsplots)*
2. Top/right spines removed? *(cnsplots: `axes_spines_top/right=False`)*
3. Legend has no frame/box? *(cnsplots: `legend_frameon=False`)*
4. Panel labels are bold, lowercase or per-house-convention, not sentence-case titles?
   *(cnsplots: `panel_label_fontweight="bold"`)*
5. Significance shown as stars/thresholds, not raw uninterpreted p on every panel?
   *(cnsplots: `pvalue_format="star"`)*
6. Fonts are Helvetica/Arial family throughout, one size scale (title vs tick vs legend), not
   mixed default-matplotlib sizes? *(cnsplots font stack + size hierarchy)*

General journal practice (not from cnsplots specifically):
7. No sentence-length panel titles ("Figure showing that X increases with Y") — a short label or
   none.
8. Direct labels on lines/bars where there are ≤5-6 series, instead of a legend forcing eye travel.
9. Restrained color: one accent hue per comparison, not a rainbow across unrelated categories.
10. Bar charts start at zero (no truncated y-axis exaggerating a difference) — this is our own
    known debt (43b panel f starts at 86%; flagged in `CRITIC_2.md`).
11. Consistent number formatting (decimal places, % vs proportion) across all panels of one figure.
12. No decorative boxes, drop shadows, 3D bar effects, or emoji/symbol clutter in labels
    (e.g. the "†" collision noted in 43b panel c).
13. Axis labels state the unit; tick labels aren't redundant with the axis label.
14. Suppressed/censored cells (`<20`) are visually distinct (hatched/grey), never rendered as zero
    or an ordinary data point — our own existing rule (`SUPPRESSED_COLOR`, `hatch_suppressed()`).

## Plot-type → our-figure mapping (cnsplots function names, vocabulary only — not calling these)

| Our figure need | cnsplots function |
|---|---|
| Bar with 95% CI (allele frequency by ancestry) | `barplot` |
| Heatmap (LD, cross-cohort concordance) | `heatmapplot` (returns a clustermap-like object, not a plain `Axes`) |
| Discovery/saturation curves | `lineplot` |
| Stacked composition (novelty tier by gene) | `stackplot` |
| Scatter (template_distance vs n_aa_changes) | `scatterplot` / `regplot` |
| ROC-style discovery curve compare | `rocplot` (bootstrap CI bands + paired DeLong test — a possible upgrade over our current CI method) |
| Forest plot (effect sizes, e.g. KIR-HLA ligand OR) | `forestplot` |
| Ternary/co-occurrence-like set overlap | `vennplot` / `upsetplot` |

## Adopted in `_viz_common.py` (S04 WS-C, 2026-09-25)

The port above landed in `scripts/hla_popgen/_viz_common.py`'s `NATURE_RC` / `nature_style()`.
Full rationale (including deliberate deviations) is inline in that file's own comment block just
above `NATURE_RC`; summary:

- **Adopted verbatim**: `legend.frameon=False`; `axes.grid=False` made explicit; `xtick.major.size`
  / `ytick.major.size` shortened 3.5 -> 2; `mathtext.fontset="custom"`; `font.sans-serif` extended
  with `"Helvetica Neue"`, `"Nimbus Sans"`, `"Liberation Sans"` fallbacks; a new `LEGEND_KW` dict
  (`markerscale=0.5, handlelength=0.7, handletextpad=0.3, frameon=False`) for the rare unavoidable
  legend; a new `JOURNAL_PALETTES` dict (nature/science/lancet/nejm/cell hex lists) plus
  `ACCENT_COLOR` for restrained one-hue-per-comparison figures, kept strictly separate from
  `ANCESTRY_COLORS` (Okabe-Ito stays the one ancestry palette).
- **Deliberately NOT adopted**: `savefig.dpi=288` (we keep 600 -- already our project-wide
  convention, matches Nature's own >=300dpi combination-art floor, and CRITIC_2 already reviewed
  figures assuming 600dpi); `savefig.transparent=True` (we keep opaque backgrounds -- these PNGs
  are shared flat in reports/Slack, where transparency on a dark viewer background is illegible);
  cnsplots' own plotting functions / pvalue-star / forest / venn helpers (out of scope -- "port,
  don't depend" covers rcParams + palette values only).
- **Backward compatible**: every existing name (`nature_style`, `save_fig`, `mm`, `panel_letter`,
  `SUPPRESSED_COLOR`, `hatch_suppressed`, `diverging_cmap`, `diverging_norm`, `NATURE_RC`,
  `NATURE_SINGLE_COL_MM`/`NATURE_DOUBLE_COL_MM`, `ANCESTRY_COLORS`) is unchanged; older scripts
  that import `_viz_common` need no changes. Covered by
  `scripts/hla_popgen/tests/test_viz_common_style.py`.
- **Figures re-rendered under the new style** (content unchanged, only rcParams + two known-debt
  fixes): Figure 1 v5, DQ G1/G2 (main + per-ancestry supplements + bimodality), KIR saturation
  (main/supplement/threshold/allele-space panels), KIR full-cohort multi-panel. See
  `sprints/S04_kir_recurrence_style_share/FIGURES_INDEX.md` for paths and per-figure notes, and
  `style_before/` / `style_after/` in that sprint folder for the PNG diffs.
