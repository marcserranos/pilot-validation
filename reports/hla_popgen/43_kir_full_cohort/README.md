# 43 — KIR full-cohort run: figures + writeup (S03 WS6 follow-on)

## Question

Does the 20-/170-person KIR pilot (`41_kir_scoping`, go/no-go: GO) hold up at full-cohort scale
(~12K people), and what does the full-cohort KIR gene-content, novelty, and ancestry picture look
like once every unrelated person with an AoU long-read assembly has been run through Immuannot's
chr19 KIR window? This report aggregates and visualizes that full run; it does not re-derive any
new statistic beyond what `43_kir_full_aggregate.py` already computed (this script only adds the
figure).

## Method

- **Pipeline**: same `run_immuannot_person.py` used for HLA, retargeted with `--region
  chr19:54,600,000-54,920,000 --pad 100,000` (the 320 kb window scoped and piloted in
  `reports/hla_popgen/41_kir_scoping/README.md`, anchored on LILRB1/FCAR outside the ~143 kb core
  array) and `--out-suffix .kir` so outputs never collide with the canonical HLA tables (verified:
  md5 of all 14 HLA tables identical before/after this run).
- **Cohort**: the paf_region (Tier-1/Tier-2) trim tier only. Tier-3/`self_align_needed` sequel2
  people (n=991) were excluded upstream of this run — the orchestrator's phase-2 self-align path
  would reuse a **chr6** reference cache for a chr19 target (a known bug documented, not fixed, in
  `41_kir_scoping/README.md` sect 2) — so running them here risked either a crash or a silently
  wrong chr6-anchored self-alignment. See Caveats.
- **Unrelatedness**: same KING-kinship greedy-unrelated-set definition as `24_novelty_by_field.py`/
  `29_hla_ld_by_ancestry.py` (`KIN_MIN = 0.0442`, third-degree-or-closer removed), reused via
  `41_kir_pilot.py`'s `greedy_unrelated()` — not re-derived.
- **Novelty tiers**: 4-tier decomposition (known / novel genomic-only / novel CDS-synonymous /
  novel protein) derived from Immuannot's own `cds_distance`/`cds_mut` GTF fields against IPD-KIR's
  `CDSseq/*.fa.gz` reference, identical logic to the pilot's `classify_novelty_tier()` (imported,
  not re-implemented, by `43_kir_full_aggregate.py`).
- **Statistics**: Wilson score interval (z=1.96) for every proportion (presence %, novelty %,
  co-occurrence %, cA/cB %) — preferred over Wald at these allele/gene frequencies, consistent with
  every other figure in this project. Framework-gene miss classification (`flanked` /
  `ambiguous_partial_flank` / `edge_fragmented`) reuses the pilot's `classify_framework_miss()`
  rule verbatim (checks whether genes centromeric/telomeric of a missing framework gene were
  themselves called on that haplotype).
- **Disclosure**: every raw count 1–19 is masked as the literal string `<20` (never 0); any rate
  whose numerator or denominator falls in that band is blanked entirely (both point estimate and
  CI) — applied uniformly in `43_kir_full_aggregate.py`, not re-derived here. No person_id was ever
  written to a file or printed. The figure hatches the 7 gene×MID cells this masking blanks in
  `kir_gene_by_ancestry.tsv`'s `novel_protein_pct` column, rather than dropping or zeroing them.

## Run / cost / QC facts

- **Run**: full-cohort KIR Immuannot pass on AoU long-read assemblies, window
  `chr19:54,600,000-54,920,000` (+100 kb Tier-1 `.paf` pad), on `n2-highcpu-80` (80 vCPU),
  concurrency 20 × 4 threads/person, wall time **31.1 h** (+0.58 h benchmark, +0.16 h retry pass),
  ~$2.98/h → **~$95 compute total**.
- **Attempted**: 12,261 people (paf_region tier). **991 Tier-3/self_align_needed sequel2 people
  excluded** upstream (see Caveats). 12,157 people (99.15%) produced both haplotypes; 104 people
  (0.85%) produced only one haplotype (deterministic — the retry pass reproduced the exact same
  104, plausibly a real assembly that doesn't span KIR on one haplotype, not a transient failure);
  **0 people produced no output at all**. **Open item (Critic #2, not reconciled here):** this
  "104" is not directly reproducible from the six committed TSVs -- the closest derivable
  quantities are 9 haplotype-level parse failures (`kir_run_summary.tsv`'s "all" row:
  24,522 attempted vs 24,513 parsed valid) and 123 zero-KIR-call haplotypes in that same "all"
  bundle (24,513 valid vs 24,390 with ≥1 call), neither of which is 104. If "104" comes from a
  different source (e.g. the orchestrator's own per-person completion log, external to this
  aggregation script), cite that source here; otherwise reconcile the number against the TSVs
  before treating it as a QC fact.
- **HLA outputs verified untouched**: md5 of all 14 HLA production tables identical before and
  after this run.
- **Unrelated set**: 11,882 people (379 dropped, KING kin ≥ 0.0442) → 23,764 haplotypes attempted,
  23,755 parsed valid, 23,637 (99.5%, CI 99.4–99.6%) with ≥1 KIR call.
- **Ancestry (predicted, unrelated set)**: AFR 3,025 / AMR 2,665 / EAS 1,463 / EUR 2,981 / MID 488 /
  SAS 1,237.
- **QC**: framework-gene presence 93.7–96.5% (KIR3DL3 96.5%, KIR3DP1 94.3%, KIR3DL2 94.1%, KIR2DL4
  93.7%); of the framework-gene misses, **4,301/5,046 (85%) are "flanked"** (genes on both sides
  called — contig clearly spans the gap, so not a truncated-assembly artifact), 693 (14%)
  ambiguous, only 52 (1%) at an array edge with no flanking gene at all. 118/23,755 haplotypes
  (0.5%) had zero KIR calls. 0/12,261 people had missing/corrupt GTF output.

## Results — key numbers

- **9.23 mean KIR genes called per haplotype** (unrelated set, 23,637 haplotypes with ≥1 call),
  **58.9% of the 222,771 total allele calls are novel** (no exact IPD-KIR match; Wilson CI
  58.7–59.1%), 0.02% undetermined (40/222,771 -- `kir_run_summary.tsv`'s own 1-decimal
  `pct_undetermined` column rounds this to "0.0"; both are the same underlying count, just
  different display precision).
- **Novelty is still overwhelmingly non-coding**: pooling across genes, the "novel, genomic-only"
  tier dominates every gene's novel calls; "novel protein" tops out at **22.5% (KIR2DL5B)** and is
  ≤10% for every other gene (see `kir_gene_summary.tsv`, Figure panel b).
- **Presence per gene** ranges from 15.3% (KIR2DS3, a "B-content" variable gene) to 96.5%
  (KIR3DL3, the centromeric-most framework gene) — panel a.
- **cA/cB haplotype content**: 13,196 cA (55.8%) vs 10,441 cB (44.2%) overall, with real ancestry
  structure: EAS highest cA (66.7%, CI 65.0–68.4%), SAS lowest (41.4%, CI 39.4–43.3%) — panel e.
- **QC sanity checks**: KIR2DL2/KIR2DL3 co-occurrence 0.8% (178/22,409) and KIR3DL1/KIR3DS1
  co-occurrence 1.8% (390/21,890) — both low as expected for allelic-alternative gene pairs, but
  **not exactly 0%** (see Caveats).
- **Allele diversity**: the most common KIR2DL1 allele (KIR2DL1*0030204) still accounts for only
  15.2% of known KIR2DL1 calls — up to 147 distinct known alleles observed per gene
  (`kir_allele_freq.tsv`).

## Comparison with the 20- and 170-person pilots

| Metric | 20-person pilot | 170-person pilot | **Full cohort (n=11,882 unrelated)** |
|---|---|---|---|
| Mean KIR genes/haplotype | 9.32 | 9.06 | **9.23** |
| % novel alleles | 60.8% | 61.0% | **58.9%** |

20-person figures per `41_kir_scoping/README.md`. The 170-person figures are **not** in that
README (it only describes launching the extended run, never got a results update after it
finished) -- they are traceable to `sprints/S03_call8_figures_kir_prediction/LOG.md`
(`[09-24 ~23:40]` entry: "170/170 people, 340/340 haps, 9.06 genes/hap, 61.0% novel") and to
`reports/hla_popgen/41_kir_scoping/41_kir_pilot_summary.kir41.tsv`
(`q_mean_kir_genes_per_hap`/`q_pct_novel` columns), which this table reproduces exactly.

The full-cohort numbers land squarely between the two pilots for gene count and just a couple of
points below both pilots for novelty — directionally consistent, not a surprise given the extra
statistical power (a 2–3 point drop at 100x the sample size is well within what a couple of
higher-novelty genes' pilot sampling variance could produce). The pilots' qualitative conclusions
— high but plausible novelty, framework genes near-universal, correct/near-correct mutual
exclusivity for the two allelic-alternative gene pairs — all replicate at full scale.

## How to read the figure (`fig_kir_full_cohort.png`/`.pdf`)

Genes in panels a–d are ordered centromeric → telomeric per the standard published KIR gene map
(framework genes marked with **†**); this ordering is for panel layout only (see the figure
script's docstring for why it differs cosmetically from `41_kir_pilot.py`'s internal
`FRAMEWORK_GENES` adjacency list, which is used for a different purpose and unaffected by the
choice).

- **a — Presence per gene**: % of haplotypes carrying each gene, Wilson 95% CI, unrelated set.
  Framework genes (dark blue) sit near 94–97%; "B-content" genes (KIR2DS*, KIR2DL2, KIR2DL5A/B,
  KIR3DS1) are the ones under ~30%, as expected from KIR's cA/cB haplotype-content biology.
- **b — Novelty per gene**: stacked bar, % of that gene's allele calls in each of the 4 novelty
  tiers (known / novel genomic-only / novel CDS-synonymous / novel protein); `n=` at each bar's
  end is the total call count for that gene.
- **c — Presence × ancestry**: heatmap, % of haplotypes carrying each gene, split by predicted
  ancestry (unrelated set). Row order matches panel a/b.
- **d — Protein novelty × ancestry**: heatmap, % of each gene's calls that are novel-protein-tier,
  by ancestry. **Hatched cells** are the 7 gene×MID cells where the underlying numerator is <20 —
  censored per this project's disclosure rule, never plotted as 0 or omitted silently.
- **e — cA vs cB content by ancestry**: % of haplotypes classified cA (fixed, inhibitory-dominated
  gene content) vs cB (carries ≥1 activating/variable gene), Wilson CI, dashed line at 50%.
- **f — QC strip**: framework-gene presence % (left) and the two mutual-exclusivity sanity checks
  (right, text annotation) — KIR2DL2/KIR2DL3 and KIR3DL1/KIR3DS1 are allelic alternatives at the
  same locus and are expected to co-occur on a haplotype only rarely/never.

## Caveats

- **Tier-3 exclusion (991 sequel2 self_align_needed people)**: excluded upstream of this run, not
  by this aggregation script, because of the chr6-reference-cache bug documented in
  `41_kir_scoping/README.md` sect 2 (never fixed, only worked around for the pilot's own much
  smaller Tier-3 test). This full run therefore has **zero Tier-3 coverage** — if sequel2 people are
  not uniformly distributed across ancestry groups (plausible; platform mix has skewed by ancestry
  elsewhere in this project), the ancestry breakdown here inherits whatever platform-ancestry skew
  exists in the sequel2-only population. Not quantified in this report — would need a
  platform × ancestry cross-tab of the excluded 991, which risks small-cell disclosure and is
  flagged rather than computed here.
- **Novelty is not yet cross-validated by recurrence.** As in the pilot, "novel" means "no exact
  match in IPD-KIR's 1,530-allele genomic reference panel" — a real signal given how much smaller
  IPD-KIR is than IPD-IMGT/HLA, but this run does not repeat script 03's HLA-side recurrence check
  (an allele seen in ≥2 unrelated people is far more likely to be real than a true assembly/calling
  artifact). That cross-check is future QC work, not done here.
- **Framework-gene misses are mostly "flanked" (85%)**, meaning contig-too-short/fragmentation is
  an unlikely explanation for most of them — genuine rare deletion haplotypes (documented in the
  literature for KIR3DP1/KIR2DL4 in particular) or a miscall-as-a-different-gene are more plausible,
  but this aggregate-only pass cannot distinguish those two explanations from each other.
- **KIR2DL2/KIR2DL3 co-occurrence is 0.8%, not the ~0% expected** for genes occupying the same
  locus as strict allelic alternatives (the pilot's 20-person run measured a clean 0.0%). At full
  scale, 178 of 22,409 haplotypes carrying either gene carry both. Possible explanations not
  distinguished by this aggregate pass: genuine rare structural duplication (both genes truly
  present on one haplotype — documented in the literature, though rare), a phasing error placing
  both genes' calls on the same haplotype when they're actually on different ones, or an
  Immuannot detection edge case at KIR's higher sequence divergence. Worth a targeted per-haplotype
  look before treating any individual "both present" call as ground truth.
- **Ancestry is a prediction**, not self-report (same caveat as everywhere else in this project).
- **Sample sizes still small enough to matter for MID**: MID is the smallest ancestry group here
  (488 people / 966 haplotypes with ≥1 call) — its per-gene protein-novelty% is hatched/censored
  in panel d for 7 of 17 genes, and its other point estimates (panel c, e) carry visibly wider CIs.

## Distilled

- Full-cohort KIR run: 11,882 unrelated people, 23,637 haplotypes with ≥1 KIR call (99.5%), ~31.1 h
  / ~$95 on `n2-highcpu-80`. HLA tables verified untouched (md5 match).
- **9.23 genes/haplotype, 58.9% novel alleles** — replicates the 20-/170-person pilots (9.32/9.06
  genes, 60.8%/61.0% novel) closely; novelty is still overwhelmingly non-coding (protein-novelty
  ≤22.5% per gene, usually single digits).
- Framework genes 93.7–96.5% present; 85% of the misses are "flanked" (not truncated assembly).
  cA/cB content varies meaningfully by ancestry (EAS 66.7% cA vs SAS 41.4% cA).
- **KIR2DL2/KIR2DL3 co-occurrence is 0.8%, above the expected ~0%** — flagged, not resolved; needs
  a targeted per-haplotype look, not treated as ground truth here.
- **991 Tier-3 (sequel2) people excluded** (known chr6-reference-cache bug, not fixed) — full run
  has zero Tier-3 coverage; possible undocumented ancestry skew from that exclusion.
- Novel-allele recurrence cross-validation (the equivalent of HLA's script-03 check) is still
  open — full-cohort scale is exactly what would make that check trustworthy, and it isn't done
  yet.
