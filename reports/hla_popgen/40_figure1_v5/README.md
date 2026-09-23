# 40 — Figure 1 v5

*WS4, S03 sprint. Implements Cole's 2026-09-22 panel feedback
(`sprints/CALL_SUMMARY_2026-09-22.md` Sec 6), then the orchestrator's 2026-09-23 review of the
first v5 draft. Built by `scripts/hla_popgen/40_figure1_v5.py` (+ `40a_admixture_bins.py` for
panel a).*

## Status: panels b/c/d/e/f complete and real; panel a still pending a VM run

Panel **a** needs a fresh VM aggregation (`40a_admixture_bins.py`, written this session, not yet
run — see "Open issues"). Every other panel is built from genuine, already-committed,
disclosure-safe aggregate tables. The figure ships with an honest placeholder box in panel a's
position.

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
**(a)** Genetic-ancestry composition of unrelated long-read participants, binned into consecutive
groups of ≥20 people (sorted by predicted ancestry, then by dominant-ancestry probability
descending) — mean admixture proportion per bin *(pending VM run of `40a_admixture_bins.py`)*.
**(b)** HLA-B alleles (≥20 carriers, strict ancestry probability ≥0.98) on the AFR/EUR/AMR simplex
at the mean renormalised ancestry composition of their carriers; point area/opacity ~ carrier
count (log). Triangles: alleles first observed in this cohort. **(c)** Percentage of called
haplotypes carrying sequence, at any IPD-IMGT/HLA nomenclature field, absent from the catalogue
(mostly non-coding), per classical gene (rows) x ancestry (columns). Hatched, "≤x"-labelled cells:
the numerator has a censored (<20) contribution, so only an upper bound is shown, never a point
value; "n/a": denominator <20 called haplotypes. **(d)** Distinct protein alleles not in
IPD-IMGT/HLA, same gene rows as c, ancestry-pooled, stacked by unrelated-carrier recurrence (seen
once / 2-19 / ≥20); ancestry-split version in `34_novel_recurrence/`. **(e)** Per-ancestry
allele-discovery curves (8 classical genes pooled), mean ± 95% band over 25 random orderings per
ancestry (`pred` ancestry scheme), direct end labels. **(f)** DQA1~DQB1 cis-haplotypes observed in
Petersdorf et al. (2022) G1xG2 cross-group cells, vs. expected under independence, by ancestry —
essentially zero in every group, consistent with background selection purging the
non-functional heterodimer.

## N per panel

- a: pending VM run (bin size ≥20; ~500 bins over ~11,800 unrelated people expected).
- b: 53 HLA-B alleles, strict ancestry ≥0.98, ≥20 carriers.
- c/d: 8 classical genes x 6 ancestries; d's per-gene novel-protein-allele totals: HLA-A 34,
  HLA-B 60, HLA-C 51, HLA-DPA1 21, HLA-DPB1 32, HLA-DQA1 38, HLA-DQB1 30, HLA-DRB1 35.
- e: 6 ancestries, 25 permutations/point, N from 487 (MID) to 3,020 (AFR).
- f: 6 ancestries, 563-4,146 DQA1~DQB1 cis-haplotypes each.

## Open issues

1. **Panel a not rendered** — `40a_admixture_bins.py` is written (bins unrelated people into
   ≥20-person groups sorted by ancestry then dominant probability, exports per-bin means only —
   disclosure-safe by construction) but has not been run on the VM this session. VM Chrome tools
   were refused by the auto-mode permission classifier in the prior iteration of this task; this
   iteration was explicitly local-only per the orchestrator's instruction. Next step:
   ```
   python3 scripts/hla_popgen/40a_admixture_bins.py --out reports/hla_popgen/40_figure1_v5/panel_a_admixture_bins.tsv
   python3 scripts/hla_popgen/40_figure1_v5.py vm-panels --out-dir ~/s03/results/40   # panel b, if re-pulling
   python3 scripts/hla_popgen/40_figure1_v5.py compose --pick A
   ```
2. Panel c/d's ancestry scheme (24's own "strict" default threshold) is not unified with panels
   b's 0.98 — would need a VM rerun of script 24 at 0.98, out of scope for this figure.
3. Panel e's Clench extrapolation is being corrected concurrently elsewhere; this script only
   plots raw curve points, which are unaffected, but a final diff against
   `39_saturation_by_ancestry/README.md` is worth doing before submission.
4. AFR/AMR end labels in panel e sit close together (both curves converge near N=3,000) — legible
   but tight; could use a leader line if a reviewer flags it.
5. A/B layout choice was a visual call at 100% zoom, not quantitative.

## Files

- `40_figure1_v5.py` — compose script (`vm-panels` / `compose` modes).
- `40a_admixture_bins.py` — VM-side panel (a) binning script (not yet run).
- `layout_A.png/.pdf`, `layout_B.png/.pdf` — both variants.
- `figure1_v5.png/.pdf` — chosen layout (A), current state (panel a pending).
- `compose_summary.json` — parameters used for this compose run.
