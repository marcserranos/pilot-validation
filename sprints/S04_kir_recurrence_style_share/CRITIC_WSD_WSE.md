# Critic — WS-D (48) / WS-E (49) / WS-D explainer, fresh-context review (2026-09-27)

Reviewed: `reports/hla_popgen/48_kir_hla_ligand_cooccurrence/` (README, all 4 TSVs,
`fig_kir_hla_ligand_forest.png` read at full size + cropped), `scripts/hla_popgen/
48_kir_hla_ligand_cooccurrence.py`/`48b_ligand_figure.py`; `reports/hla_popgen/
49_cole_share_package/README.md` and the README/SCHEMA templates in `scripts/hla_popgen/
49_cole_share_package.py`; `sprints/S04_kir_recurrence_style_share/WS-D_EXPLAINER_AND_IDEAS.md`.

## Disclosure

No violation found. Every 2x2 cell in `kir_hla_ligand_cooccurrence.tsv` in [1,19] renders `<20`
(never a bare number, never 0); the corresponding OR/CI/p columns are genuinely blank (empty
string), not `0` — an earlier `column -t` rendering of the raw TSV made a masked row look like
`odds_ratio=0`, but the raw tab-delimited fields are empty, confirmed with `grep -P` (verified
EAS/2DL1xC2, EAS/2DL3xC1, MID/2DL1xC2, MID/3DL1xBw4, SAS/2DL1xC2 — all 5 masked rows). No allele
name appears next to a <20 count anywhere (the crosscheck table is catalogue-level, two-field
group names only, never paired with a carrier count). No person IDs in any of the 4 TSVs or the
49 templates' printed output.

## Number traceability (5+ spot-checks)

- EAS C1/C2/Bw4 (1307/1460=89.5%, 555/1460=38.0%, 1102/1460=75.5%) — exact match, README table.
- Bonferroni threshold 0.05/30 = 0.0016667 ≈ "0.0017" as stated — correct, and both the README and
  figure caption consistently say neither nominal signal (EAS 2DL2xC1 p=0.011, EUR 2DL3xC1
  p=0.013) clears it — no overclaiming.
- Sequence-vs-lookup source mix (HLA-C 20,671/2,533; Bw 41,910/4,481) matches
  `ligand_lookup_qc.tsv` exactly.
- `n_two_field_groups_seq_vs_lookup_{compared,agree,disagree}` = 235/228/7 matches TSV and README.
- AFR 2DL1xC2 OR 1.504/CI 0.989-2.285/perm_p 0.098 — matches TSV row exactly (README rounds to
  1.50/0.99-2.29, consistent).

## Biology (ligand definitions + citations)

C1=Asn80(+Ser77)/C2=Lys80 and Bw4=Arg83+I80/T80 vs. Bw6=Ser77-Asn80-Leu81-Arg82-Gly83 are stated
correctly and consistently in the script docstring, README, and test file, matching the standard
KIR-ligand literature (Colonna 1993, Winter & Long 1997 for C1/C2; Gumperz 1995, Cella 1994, Parham
2005 review for Bw4/Bw6). Receptor-ligand pairings (2DL1-C2, 2DL2/3-C1, 3DL1-Bw4, 3DS1-Bw4-80I via
the Martin 2002 epidemiologic association, explicitly flagged as not a confirmed direct-binding
pair) are the standard set and correctly hedged. I do not have literature access to independently
verify the exact page/table citations, only that the biology stated is textbook-consistent.

## Figure

**[Fixed — major, legibility] Panel b dot occlusion.** At full resolution, panel b's compact
carrier-frequency dot plot placed all 6 ancestries' points at the literal y-value for each epitope
row; when two ancestries' `pct` values landed within a marker-width of each other (e.g. Bw4: AFR
74.0 vs EUR 74.6 vs EAS 75.5), a later `scatter()` call painted over an earlier point of a
different color. Cropped inspection of the original PNG showed only 4-5 of 6 ancestry dots visible
per epitope row (AFR's Bw4 point was fully hidden behind EUR/EAS). Not a disclosure issue (exact
values are still in `epitope_freq_by_ancestry.tsv` and the README table) but a real
readability/interpretability bug in a figure meant to stand on its own. **Fixed** in
`scripts/hla_popgen/48b_ligand_figure.py`: added a small, fixed, order-dependent per-ancestry
vertical jitter (±0.175 max) in panel b only, so all 6 dots stay visible regardless of how close
the underlying percentages are. Re-rendered `fig_kir_hla_ligand_forest.{png,pdf}` and re-inspected
crops of the Bw4/C2/C1 rows — all 6 ancestry colors now distinct in every row.

Panel a (forest plot) was legible at full resolution: no overlapping labels, hatched/censored rows
correctly distinct from real points, `*` annotations correctly placed on only the two nominal rows.

## READMEs — question/method/results/how-to-read/caveats/Distilled

**48 README**: has all six sections plus a plain-language paragraph for Marc. Complete.

**49 README (report-tier)**: has Question/Method/Results(TBD, correctly labeled since the VM run
hadn't happened at write time)/How to read/Caveats/Distilled. Complete, though "Results" is a
placeholder — acceptable since 49's real output (the package) is VM-only and never enters this
repo; the report-tier README's job is describing shape/provenance, which it does.

**[Fixed — minor, cosmetic] 49's package README/SCHEMA date bug.** The template stamped
`release {date_str}` using `time.strftime("%Y-%m-%d")` (wall-clock at build time), independent of
the bucket upload prefix chosen separately in the module's own docstring
(`release_2026-09-25/`). The already-uploaded package's docs read "release 2026-09-26" against a
`release_2026-09-25` bucket path — confirmed by reading `scripts/hla_popgen/
49_cole_share_package.py` lines ~107 (docstring, hardcoded `release_2026-09-25`) vs. the old
`date_str = time.strftime("%Y-%m-%d")` at the call site. Not a disclosure or data-correctness
issue — only the doc header string. **Fixed**: added `--release-tag` (defaults to today's date if
omitted, with a docstring warning against relying on the default); `date_str` now comes from
`args.release_tag`. Updated the module's own "Usage" block to show `--release-tag 2026-09-25` and
noted the bug in the report-tier README's Caveats so the next release remembers to pass it.
**Judgment call**: did not re-upload/regenerate the already-shared 2026-09-25 package to fix its
cosmetic header — re-touching an already-delivered in-perimeter share for a doc-string typo isn't
worth the VM/consent overhead; flagging here for Marc in case he disagrees.

## WS-D explainer

Part A (WS6 in plain language) is genuinely readable by a non-ML PI: it defines every term
inline (V/J usage, 4-mers, L1, CV, AUROC with concrete bands) and states the ancestry-confound
warning in causal, non-jargon terms. No overclaiming — "no disease could be tested yet" is stated
plainly as a power problem, not spun as a negative result.

The Part B rankings are sensible: cheapest-first (48, then its free byproduct epitope frequencies,
then the QC battery) ahead of the higher-value-but-expensive WS6 full run, which is correctly
flagged as deserving its own dedicated session rather than being squeezed in.

**[Fixed — major, staleness] "What's pending" section contradicted the "Results" section directly
above it.** The Results section (added 2026-09-26, dated "supersedes... where they conflict")
correctly states both scripts ran on the full 11,845-person cohort. But the very next section,
"What's pending," was never updated and still said "Both scripts (47, 48) are local-only and
synthetic-tested so far — not yet run on real AoU data," and separately claimed 48's epitope
lookup tables "need a VM-side cross-check ... before the co-occurrence numbers are treated as
final" — but 48 was rewritten (commit `d84348a`) to do sequence-derived classification as primary
with the lookup table as fallback+cross-check, and that cross-check is now a standing, already-run
output (`ligand_seq_vs_lookup_crosscheck.tsv`). A reader who read straight through would hit a
flat contradiction two sections apart. **Fixed**: struck through and annotated both stale bullets
in `WS-D_EXPLAINER_AND_IDEAS.md` with what actually happened, rather than deleting the history.

## Counts

- Blockers: 0
- Major: 2 (both fixed — panel b dot occlusion; WS-D explainer stale/contradictory "pending"
  section)
- Minor: 1 (fixed — 49 package README/SCHEMA release-date mismatch, template only; the
  already-uploaded package's cosmetic header was left as-is, a judgment call for Marc)

## Tests

`pytest scripts/hla_popgen/tests -q`: 560 passed, 9 pre-existing errors, all in
`test_allele_geometry.py`/`test_extraction.py`/`test_figures.py`/`test_novel.py` — unrelated to
44/45/46/47 or to this review's scope (missing `fixtures_dir`/other fixtures not touched here);
did not investigate further per the "don't touch 44/45/46/47" instruction, and this review made no
changes to any of those files or their tests.
