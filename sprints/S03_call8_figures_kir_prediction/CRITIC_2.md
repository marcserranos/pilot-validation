# Critic #2 — fresh-context review of 43_kir_full_cohort (S03, 2026-09-25)

Reviewed: `reports/hla_popgen/43_kir_full_cohort/` (README.md, all 6 TSVs, `fig_kir_full_cohort.png`
read directly) and `scripts/hla_popgen/43_kir_full_aggregate.py` / `43b_kir_full_figure.py`,
against `reports/hla_popgen/41_kir_scoping/README.md` (pilot/methods reference) and
`sprints/S03_call8_figures_kir_prediction/CRITIC_1.md` (format). No disclosure violation found
(no 1-19 count printed as an exact number anywhere in README prose or the TSVs; every `<20`-band
rate is blanked, not just the count; hatched cells verified to be the correct 7 gene×MID cells;
no person_id anywhere). Arithmetic and figure-correctness spot-checks did turn up real issues.

## Ranked findings

1. **[major] cA/cB overall percentage doesn't match its own numerator/denominator.**
   `reports/hla_popgen/43_kir_full_cohort/README.md:85` (before fix) stated "13,196 cA (55.7%) vs
   10,441 cB (44.3%)". `kir_run_summary.tsv`'s unrelated row gives exactly these two counts
   (13,196 / 10,441, total 23,637); 13,196/23,637 = 55.83% and 10,441/23,637 = 44.17%, i.e.
   **55.8%/44.2%, not 55.7%/44.3%** — off by 0.1 pp in a way plain rounding doesn't explain. Not
   disclosive, doesn't change any interpretation, but fails the "arithmetically right" check.
   **Fixed** (see below).

2. **[major] "104 people produced only one haplotype" is not reproducible from the committed
   TSVs.** `README.md:50` (original). The only failure-count derivable from `kir_run_summary.tsv`'s
   "all" row (12,261-people cohort) is **9** haplotype-level parse failures (24,522 attempted vs
   24,513 parsed valid), and `kir_qc.tsv`'s `people_missing_or_corrupt_gtf` row independently
   confirms 0 people have *both* haplotypes missing — so at most 9 people could have exactly one
   haplotype missing, not 104. A nearby but distinct quantity, zero-KIR-call haplotypes in that
   same "all" bundle, is 123 (24,513 valid − 24,390 with ≥1 call) — also not 104. Whatever "104"
   measures, it isn't computable from these six tables as written, and the prose conflates
   "produced only one haplotype" (implies a missing file) with "doesn't span KIR on one haplotype"
   (implies a zero-call haplotype that *did* parse) — two different failure modes. **Not fixed**
   (left as a documented open item — needs the source log/script that actually produced "104",
   which this reviewer doesn't have VM access to confirm).

3. **[major] Figure text below Nature's 5–7 pt floor and below this project's own type scale.**
   `scripts/hla_popgen/43b_kir_full_figure.py` set panel c/d heatmap cell numbers at **3.6 pt**
   (line 148, pre-fix), axis tick/colorbar labels at 4.2–4.6 pt (lines 151, 158, 174), panel b's
   `n=` annotations at 4.3 pt (line 121), panel f's co-occurrence text at 4.6 pt (line 232), and
   the legend/footnote at 4.6–4.8 pt (lines 266, 273) — all well under both the stated 5–7 pt
   floor and this project's own `_viz_common.nature_style()` scale (`font.size: 7`,
   `legend.fontsize: 6`). **Fixed**: bumped every one of these to 5.0–5.5 pt and re-rendered;
   re-inspected the new PNG at native resolution — no new overlaps introduced by the larger text
   (row labels in panel c/d still clear of the heatmap and of each other; panel a/b y-labels,
   panel e/f tick labels all legible with margin).

4. **[minor] Undetermined-call rate stated at higher precision than its source TSV.**
   `README.md:70` (original) said "0.02% undetermined" while `kir_run_summary.tsv`'s
   `pct_undetermined` column (1-decimal format) reads "0.0" for the same row — same underlying
   count (40/222,771 = 0.018%), different display precision, so a reader checking the TSV
   literally would see a mismatch. **Fixed**: added the raw count and a one-line note explaining
   the precision difference.

5. **[minor] "170-person pilot" comparison-table numbers aren't sourced in the cited report.**
   `README.md:94-96`'s comparison table cites 9.06 genes/haplotype and 61.0% novel for the
   170-person pilot, attributing it (via the section header) to `41_kir_scoping/README.md` — but
   that README (read in full for this review) only describes *launching* the 170-person run
   (§4a) and was never updated with its results. The numbers are correct — verified against
   `sprints/S03_call8_figures_kir_prediction/LOG.md`'s `[09-24 ~23:40]` entry and
   `reports/hla_popgen/41_kir_scoping/41_kir_pilot_summary.kir41.tsv`
   (`q_mean_kir_genes_per_hap=9.06`, `q_pct_novel=61.0`) — just untraceable from the one document
   a reader would naturally check. **Fixed**: added a footnote citing the actual sources.

6. **[minor] Ancestry people-counts in the "Run / cost / QC facts" section aren't in any TSV.**
   `README.md:57-58`, "AFR 3,025 / AMR 2,665 / EAS 1,463 / EUR 2,981 / MID 488 / SAS 1,237" (sum
   11,859) — these per-ancestry person counts exist only as a stderr log line in
   `43_kir_full_aggregate.py` (`log(f"[43]   {a}: {n} people (unrelated)")`, never written to a
   TSV). Cross-checked for internal plausibility (sum 11,859 + 23 unassigned-ancestry people not
   broken out ≈ 11,882 unrelated total; 23 people × 2 haplotypes ≈ the 46-haplotype gap between
   `kir_content_by_ancestry.tsv`'s ancestry-summed total of 23,591 and `kir_run_summary.tsv`'s
   23,637 unrelated-with-calls figure — consistent, not a red flag) but not independently
   verifiable against a committed table. **Not fixed** (would mean adding a new output column to
   `43_kir_full_aggregate.py`, a pipeline change rather than a text fix — flagged for the human).

7. **[minor, scientific-judgment — not fixed] Gene physical order (panels a–d) asserted against
   two citations this reviewer could not independently verify.** `43b_kir_full_figure.py:14-20`'s
   `GENE_ORDER` cites "Kulkarni et al 2008 Fig 1, Middleton & Gonzalez 2010 Table 1" for placing
   `KIR2DL5B/2DS3/2DS5` centromeric of `2DP1/2DL1` and `2DL5A/2DS1/2DS4` telomeric of `3DL1/3DS1`.
   This is plausible (2DL5B vs 2DL5A's centromeric/telomeric split in particular is a genuinely
   established fact) but the duplicated-locus genes (2DS3/2DS5 especially) are documented in the
   literature as variably placed depending on haplotype structure, and this reviewer has no
   literature access to confirm the exact cited figures place them where this script does. Not a
   confirmed error — flagging so a human with the actual references can do a 30-second check
   before this ships externally.

## Everything checked and found correct (no action needed)

- Every headline number in "Results — key numbers" and the QC section traced cleanly to a TSV cell
  and the arithmetic reproduced (9.23 genes/hap; 58.9% novel/222,771 calls; 15.3–96.5% presence
  range; 22.5% max novel-protein and ≤10% elsewhere; 4,301/693/52 framework-miss classification
  summing to 5,046; 178/22,409 and 390/21,890 co-occurrence figures; EAS/SAS cA extremes and CIs).
- Panel c/d heatmap cell values spot-checked against `kir_gene_by_ancestry.tsv` for KIR3DL3 (all 6
  ancestries) and KIR2DL5B (AFR) — exact match after rounding.
- The 7 hatched (censored) cells in panel d correspond exactly to the 7 gene×MID rows in
  `kir_gene_by_ancestry.tsv` where `n_novel_protein` is `<20` — never drawn as 0, matches the
  disclosure rule and the figure's own caption.
- Method statements checked against code: unrelatedness (`greedy_unrelated`) is computed only over
  `discover_people(outroot)` — i.e. only people this KIR run actually attempted, not the full AoU
  cohort, matching the README's claim. Ancestry is explicitly the `ancestry_pred` column
  (predicted, not self-report). Every CI in every table is `wilson_ci(z=1.96)`, consistent with
  the rest of this project.
- Scientific hedging (item 4 of the review brief): the 2DL2/2DL3 0.8% co-occurrence, the
  framework-gene-miss classification, the Tier-3/sequel2 exclusion, and the novelty-vs-IPD-KIR-size
  framing are all stated with explicit "not yet cross-validated" / "flagged, not resolved" /
  "cannot distinguish X from Y" language — no overclaim found in any of the four.

## Fixes applied (2026-09-25, this pass)

1. **[major]** cA/cB overall percentage corrected 55.7%/44.3% → **55.8%/44.2%**
   (`reports/hla_popgen/43_kir_full_cohort/README.md`).
2. **[major]** Figure text sizes raised from 3.6–4.8 pt to 5.0–5.5 pt across panels b/c/d/f, the
   shared legend, and the footnote (`scripts/hla_popgen/43b_kir_full_figure.py`); re-rendered
   `fig_kir_full_cohort.{png,pdf}` and re-inspected the new PNG at native resolution — no
   label/cell overlaps introduced.
3. **[minor]** Undetermined-rate precision mismatch — added raw count (40/222,771) and an
   explanatory note next to the 0.02% figure.
4. **[minor]** Added a source footnote for the 170-person pilot's comparison-table numbers
   (LOG.md + `41_kir_pilot_summary.kir41.tsv`), since `41_kir_scoping/README.md` itself doesn't
   contain them.
5. **[major]** Added an in-place "Open item" annotation on the "104 people produced only one
   haplotype" claim, spelling out the two closest-but-different derivable quantities (9 parse
   failures, 123 zero-call haplotypes) so a reader doesn't take 104 as TSV-verified.

## Left open (not fixed — scientific judgment or out-of-scope pipeline change)

- Finding 2 (104 vs 9/123) — needs the actual source of "104," not resolvable from committed
  aggregates alone.
- Finding 6 (ancestry people-counts not in any TSV) — would need a new `43_kir_full_aggregate.py`
  output column; a pipeline change, not a text fix.
- Finding 7 (gene physical order citations) — needs a human with literature access to confirm
  Kulkarni 2008 Fig 1 / Middleton & Gonzalez 2010 Table 1 actually place the duplicated loci
  (2DS3/2DS5 especially) where `GENE_ORDER` puts them.
