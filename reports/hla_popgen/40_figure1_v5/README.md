# 40 — Figure 1 v5

*WS4, S03 sprint. Implements Cole's panel-by-panel feedback from the 2026-09-22 call
(`sprints/CALL_SUMMARY_2026-09-22.md` Sec 6) on top of `36_figure1_native.py` (panels a/b) and
`33_figure1_v3_compose.py` (panel c). Built by `scripts/hla_popgen/40_figure1_v5.py`.*

## Status: panels b/c/d/e complete and real; panel a blocked this session — see "Open issues"

`40_figure1_v5.py compose` was run against genuine, already-committed, disclosure-safe aggregate
tables (not synthetic/placeholder data): `36_figure1/panel_b_ternary_alleles.tsv` (per-allele
simplex centroids, strict ancestry ≥0.98, ≥20 carriers — exactly the thresholds Cole asked to
keep) for panel b, and `24_novelty_by_field` + `34_novel_protein_recurrence`'s committed tables for
the merged panel c/d. Panel e reads `39_saturation_by_ancestry/curves.tsv` directly, so it
automatically picks up the concurrent extrapolation-fit correction another agent is making (the
curve *points* it plots — `mean_distinct`/`lo2_5`/`hi97_5` per N — are independent of the Clench
fit that agent is touching).

Panel **a** could not be regenerated this session: it needs a fresh VM render (per-person
admixture must never leave the VM) with the '≥98%' annotation dropped, per Cole's ask. The VM
Chrome tools (`javascript_tool` via `browser_batch`) were refused twice by the Claude Code
auto-mode permission classifier (`[Auto-Mode Bypass]`) when opening the VM's terminal
websocket — this matches `AGENT_PREAMBLE.md`'s explicit instruction: *"if browser tools get
refused by a permission classifier, STOP and report"*. The figure below therefore ships with an
honest placeholder box in panel a's position rather than a fabricated or stale image.

**To finish:** on a session where the VM Chrome tools are not refused, run:
```
python3 scripts/hla_popgen/40_figure1_v5.py vm-panels --out-dir ~/s03/results/40
# pull back: panel_a.png, panel_a.pdf, panel_a_group_sizes.tsv, panel_b_ternary_alleles.tsv
python3 scripts/hla_popgen/40_figure1_v5.py compose \
    --panel-a-image reports/hla_popgen/40_figure1_v5/panel_a.png \
    --panel-a-table reports/hla_popgen/40_figure1_v5/panel_a_group_sizes.tsv \
    --panel-b-table reports/hla_popgen/40_figure1_v5/panel_b_ternary_alleles.tsv \
    --pick A
```

## Chosen layout: A

Two layouts were built (`layout_A.png/.pdf`, `layout_B.png/.pdf`) and inspected at 100% zoom.

- **A**: full-width admixture strip (a) on top, then b | c+d, then e | a reserved cell for the
  pending KIR panel (WS3).
- **B**: (a) top strip, then b and e stacked in the left column, c+d spanning the full right-hand
  height.

**A is chosen.** B's taller c+d heatmap does not add legibility (8 gene rows already read cleanly
at A's height) and B's stacked b/e left column forces both panels narrower than A's, making
panel b's ternary tick labels tighter. A also leaves a natural slot for the KIR panel the sprint
is adding in parallel, which B does not.

## Panel-by-panel, vs Cole's feedback

| Panel | Ask | Done |
|---|---|---|
| a | drop '≥98%' annotation | Annotation removed from the plot in `draw_panel_a_alone()`; the per-group % clearing the strict threshold moves to `panel_a_group_sizes.tsv` only, where it belongs as data, not a plot annotation. **Rendering pending VM access** (see above). |
| b | 0–100 tick marks on all 3 axes; strict ancestry 95–98% (kept 0.98, `36` already used it); ≥20-carrier floor | Done — `draw_ternary_grid()` adds tick labels 0/20/40/60/80/100 on all three edges. Reused `36`'s already-committed table: **N = 53 HLA-B alleles**, strict ancestry ≥0.98, ≥20 carriers, from the same VM run that produced `panel_a_group_sizes.tsv` (11,833 unrelated participants; per-ancestry strict-threshold retention 20.5%–92.6%, see table below). |
| c | replace reference-incompleteness bar chart with the ancestry × gene novelty-rate panel (33's, "the 800,000 one") | Done, as a heatmap (gene rows × ancestry columns, lower-bound %), reusing `33_figure1_v3_compose.py::rate_by_gene_ancestry` verbatim — no re-derivation. Hatched cells mark a denominator <20 haplotypes. |
| d | novel proteins by gene, labelled 'protein alleles not in IPD-IMGT/HLA' (not 'not in classical genes'); expand genes; ancestry split → supplement | Done as a marginal bar attached to panel c's gene rows, correctly labelled. Restricted to the 8 classical genes shown in c (for row alignment); the full, ancestry-split gene list stays in `34_novel_recurrence/` (referenced in the caption, not duplicated). |
| c+d merge | "consider... d as a marginal bar on c's heatmap" | Done, per above. |
| e | new panel: per-ancestry allele-discovery curves from 39, direct end labels | Done — pooled classical genes, `pred` ancestry scheme (39's primary), mean ± 95% permutation band, direct labels placed with a minimum-gap declutter pass so AFR/AMR and MID/SAS (which converge) stay legible. |

## N per panel

- **a** (pending render): 11,833 unrelated participants, predicted ancestry; per-group %
  clearing strict ≥0.98 — AFR 49.2%, AMR 45.7%, EAS 84.0%, EUR 32.4%, MID 20.5%, SAS 92.6%
  (`panel_a_group_sizes.tsv`, real data, already committed by `36`).
- **b**: 53 HLA-B alleles reaching ≥20 carriers among strict-ancestry (≥0.98) unrelated people.
- **c/d**: 8 classical genes × 6 ancestries; denominators per cell are exact haplotype counts from
  `24_novelty_by_field` (`ancestry_scheme=strict`, 24's own default threshold — not re-run at 0.98
  this session; noted as an open item below). Novel-protein counts (d): HLA-A 34, HLA-B 60,
  HLA-C 51, HLA-DPA1 21, HLA-DPB1 32, HLA-DQA1 38, HLA-DQB1 30, HLA-DRB1 35.
- **e**: 6 ancestries, 25 permutations/point, `pred` scheme; N at curve end ranges MID (n=487) to
  AFR (n=3,020) — see `39_saturation_by_ancestry/README.md` for the full per-ancestry N table.

## Draft figure legend (Nature style)

**Figure 1 | Long-read HLA typing across the All of Us cohort: cohort structure, reference
catalogue gaps, and allele-discovery saturation by ancestry.**
**(a)** Genetic-ancestry composition of 11,833 unrelated long-read participants, ordered within
each predicted-ancestry block from least to most admixed *(pending VM re-render this session)*.
**(b)** HLA-B alleles (n = 53, ≥20 carriers) placed on the AFR/EUR/AMR simplex at the mean
renormalised ancestry composition of their carriers, restricted to participants with strict
ancestry probability ≥0.98; point area and opacity scale with carrier count (log). Triangles are
alleles first observed in this cohort. **(c)** Percentage of called haplotypes carrying sequence
absent from IPD-IMGT/HLA, per classical gene (rows) and ancestry (columns); hatched cells mark a
called-haplotype denominator below 20 (All of Us small-cell rule). **(d)** Distinct protein alleles
not in IPD-IMGT/HLA, per gene (same rows as c, ancestry-pooled), stacked by how many unrelated
people carry them — seen once, 2–19, or ≥20 (the reportable threshold). **(e)** Per-ancestry
allele-discovery curves: distinct HLA protein alleles (8 classical genes, pooled) found as a
function of people sampled, mean ± 95% band over 25 random orderings per ancestry (Pakistan Genome
Resource, *Nature* 2026, Fig. 3e convention); direct end labels.

## Open issues

1. **Panel a not rendered this session** — VM Chrome tools refused by the auto-mode permission
   classifier (see Status above). This is the only blocking item; script, layout, and the other
   four panels are otherwise finished.
2. Panels c/d use `24_novelty_by_field`'s own "strict" ancestry scheme (its own default
   threshold), not the 0.98 used in a/b — unifying would need a VM rerun of script 24 at 0.98,
   out of scope for this figure alone. Flagging rather than silently mixing thresholds.
3. Panel e's Clench extrapolation is being corrected concurrently by another agent; this script
   only plots the raw curve points (`mean_distinct`/`lo2_5`/`hi97_5`), which are unaffected, but a
   final check against the corrected `39_saturation_by_ancestry/README.md` before submission is
   worthwhile.
4. Two-layout comparison (A chosen over B) was a visual judgement call at 100% zoom, not a
   quantitative one — worth a second pair of eyes.
5. KIR panel (WS3) is reserved as empty space in layout A, not yet populated.

## Files

- `40_figure1_v5.py` — the script (two modes: `vm-panels`, `compose`).
- `layout_A.png/.pdf`, `layout_B.png/.pdf` — both variants, for comparison.
- `figure1_v5.png/.pdf` — the chosen layout (A), current state (panel a pending).
- `compose_summary.json` — parameters used for this compose run.
