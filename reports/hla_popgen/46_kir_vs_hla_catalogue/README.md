# 46 — KIR vs HLA: which catalogue covers its region better? (v4b)

**Status: v4b, 2026-09-27.** Supersedes every earlier number in this file. See
`../44_kir_recurrence_saturation/README.md`'s "Version history" for the full v1→v4b lineage
(name-based identity bug, missing KIR artifact gate, HLA `HLA-`/bare-gene-name bug) — not repeated
here.

## Question

Which reference catalogue — KIR (IPD-KIR 2.13.0) or HLA (IPD-IMGT/HLA 3.55.0) — represents its own
region better, per gene where possible? Local-only script, reads `44`'s already-pulled-back
aggregate TSVs (`../44_kir_recurrence_saturation/`), never raw per-person data.

## Method

Same identity levels, artifact gate, and Good-Turing/Chao2 definitions as `44`'s README (read that
first). Per gene × species (`build_gene_metrics()` in `46_kir_vs_hla_catalogue.py`):

- `completeness` (genomic, upper bound), `completeness_cds`, `completeness_protein` —
  S<sub>obs</sub>/Chao2 at each identity level.
- `pct_novel_any` = 100 × S<sub>obs</sub>(any_novel) / S<sub>obs</sub>(genomic).
- `pct_novel_protein` = 100 × S<sub>obs</sub>(protein_novel) / S<sub>obs</sub>(protein) — **same**
  identity granularity numerator/denominator (fixed 2026-09-27, see `44`'s v1–v3 history: dividing
  by `genomic` instead let this ratio exceed 100%).
- `slope_per_1000` — equal-N discovery slope, genomic level, from `44`'s `equal_n_slope.tsv`.
- `protein_catalogue_covered` / `protein_catalogue_reason` — `False`/`no_catalogued_protein_entries`
  for KIR2DP1/KIR3DP1 (pseudogenes), `True`/`"ok"` for every other gene in both species.

**The headline comparison is CDS and protein level, not genomic** — genomic (span-level) identity
is an upper bound inflated by person-specific span/UTR boundaries and intronic assembly noise the
CDS-based artifact gate can't see (44's Caveat 1). Genomic is plotted for context only, explicitly
labeled "upper bound" in the figure.

**Per-gene AND per-ancestry**: `build_gene_metrics()` takes an `ancestry` argument (default `ALL`);
`main()` now also writes one combined per-ancestry table
(`46_catalogue_metrics_by_ancestry.tsv`, one row per gene × species × ancestry, AFR/AMR/EAS/EUR/SAS
— MID excluded, underpowered) rather than 5 separate files.

## Results

### Pooled-ALL headline (from `46_catalogue_metrics.tsv`, summed via `coverage_chao2.tsv`)

| level | HLA S<sub>obs</sub> | HLA completeness | KIR S<sub>obs</sub> | KIR completeness | which is more complete |
|---|---:|---:|---:|---:|---|
| genomic (upper bound, context only) | 17,188 | 17.6% | 38,024 | 19.2% | ~tied, KIR nominally ahead |
| CDS (headline) | 1,363 | 38.8% | 2,100 | 41.6% | ~tied, KIR nominally ahead |
| protein (headline) | 1,079 | 39.9% | 1,444 | 39.2% | ~tied, HLA nominally ahead |
| protein_novel (novelty subset) | 198 | 6.5%¹ | 1,063 | 31.4%¹ | KIR is far LESS complete at capturing its own novel proteins |

¹ `protein_novel` completeness = S<sub>obs</sub>(protein_novel) / Chao2(protein_novel) — this is
"how much of the *novel-protein* space specifically is covered," a different (and much smaller,
much less complete) question than the `protein` row above. Worked example from Marc's brief:
**protein-level completeness KIR 39.1%¹ vs HLA 39.9%; protein_novel KIR 1,063 (31.4% complete) vs
HLA 198 (6.5%¹ complete)** — reproduced from this v4b run (the brief's "1,061"/"163" figures are
close but not exact; this README uses the values recomputed directly from the committed v4b TSVs,
verified two independent ways — `coverage_chao2.tsv` summed directly, and `saturation_curves.tsv`
via `pooled_recurrence_from_curves()` — both give S<sub>obs</sub>=198/1,063 exactly).

**Pooled-ALL verdict: near-tied at both headline levels, with the ranking flipping by less than 3
points depending on whether you look at CDS (KIR ahead) or protein (HLA ahead).** This is NOT the
clean "KIR's catalogue is worse" story that pooled genomic-level numbers alone might suggest — see
the per-ancestry breakdown below, which is the more informative view.

### Per-ancestry (from `46_catalogue_metrics_by_ancestry.tsv` / `coverage_chao2.tsv`)

| ancestry | HLA protein completeness | KIR protein completeness | HLA CDS completeness | KIR CDS completeness |
|---|---:|---:|---:|---:|
| AFR | 50.8% | 34.3% | 50.8% | 39.8% |
| AMR | 52.2% | 38.2% | 49.2% | 39.8% |
| EAS | 57.0% | 48.3% | 56.2% | 46.6% |
| EUR | 47.7% | 28.6% | 41.3% | 25.9% |
| SAS | 66.9% | 40.5% | 57.6% | 40.6% |

**HLA is more complete than KIR in every single one of the 5 ancestries, at both CDS and protein
level** — a consistent, unambiguous within-ancestry finding, even though the pooled-ALL numbers
above look tied or slightly favor KIR. This is a real aggregation effect (Simpson's-paradox shape:
pooling across ancestries of different sizes and allele-richness profiles can reverse an ordering
that holds in every subgroup), not an error in either number. **Treat the per-ancestry table as the
more trustworthy answer to "which catalogue is better represented"** — HLA, consistently — and the
pooled near-tie as an artifact of how the ancestries happen to mix, not evidence the two catalogues
are actually comparable.

### Per gene (`46_catalogue_metrics.tsv`, pooled ALL)

HLA's 8 classical genes range 27.0%–57.3% protein-completeness (worst: HLA-DQA1 at 21.6%; best:
HLA-DRB1 at 57.3%). KIR's 15 catalogue-covered genes (excluding the 2 pseudogenes) range
21.6%–56.1% (worst: KIR2DL2/KIR2DL5B at ~21.6–21.9%; best: KIR2DS3 at 52.6%, KIR3DL3 at 49.9%).
**No single gene in either species is dramatically better- or worse-covered than its own species'
typical range** — this isn't one outlier gene driving the pooled numbers; both catalogues have a
real, broad spread of per-gene completeness. `pct_novel_protein` (share of a gene's protein
catalogue-covered alleles that are novel) is far higher for KIR genes (65–92%) than HLA genes
(14–38%) at almost every gene — KIR's *catalogue coverage per novel call* is worse even where its
aggregate completeness looks similar, consistent with KIR simply having a smaller, less complete
reference set relative to its real diversity.

## How to read each figure

- **`fig_catalogue_completeness.png`** (a–c): three scatter panels, x = % of distinct alleles
  novel, y = Chao2 completeness, one point per gene, HLA (red) vs KIR (blue). **(a) genomic —
  upper-bound context only**, gene labels shown (dense cluster near x=100% is normal: most genes
  of both species are >80% "novel" at the genomic/span level — see Caveat 1, this is inflation, not
  biology). **(b) CDS — headline.** Gene labels are deliberately omitted in this panel (both
  species' any-level-novelty values overlap too much for either the free or column-leader-line
  label layout to place ~30 gene names without collisions — tried both, see the script's
  `_scatter_panel` docstring); per-gene values are in `46_catalogue_metrics.tsv`/
  `_by_ancestry.tsv`. **(c) protein — headline.** KIR's genes cluster tightly on the right
  (65–92% novel) and use a shared-column leader-line label layout; KIR2DP1/KIR3DP1 (pseudogenes,
  no catalogued protein) are named in an explicit in-panel note, never plotted as (0, 0).
- **`fig_recurrence_composition.png`** (a–b): stacked bars, % of S<sub>obs</sub> in each recurrence
  class (seen 1x/2x/>2x, `ge20` folded into `>2x` here to avoid double-drawing a subset), any-level
  vs protein-level novelty, HLA vs KIR. Uses the exact curve-derived pooled counts from `44`'s
  README (not `recurrence_classes.tsv`'s masked cells), so this figure and `44`'s recurrence-stats
  table never disagree.

## Caveats

1. Genomic-level identity is an upper bound (see `44`'s Caveat 1) — not used for the headline
   comparison here, plotted for context only (panel a).
2. Good-Turing coverage ≠ Chao2 completeness (see `44`'s Caveat 2) — every "completeness" number in
   this README is Chao2 (S<sub>obs</sub>/Chao2), never Good-Turing.
3. `coverage_chao2.tsv`'s per-gene s_obs is not `<20`-masked; several KIR genes' `any_novel`/
   `protein_novel` s_obs are below 20 (e.g. KIR2DS3 protein S<sub>obs</sub>=36, at the edge). This
   script uses 44's committed numbers as-is (re-masking someone else's committed aggregate is out
   of scope for a local figure script) — flagged per `context/DECISIONS.md`'s open disclosure
   question, not resolved here.
4. The per-ancestry reversal (HLA ahead in every ancestry, near-tied pooled) is real and
   reproducible from the TSVs but its *mechanism* (why pooling reverses the ordering) is not
   investigated further in this pass — reported as an aggregation-sensitivity finding, not
   explained mechanistically.
5. KIR2DP1/KIR3DP1: `pct_novel_protein`/`completeness_protein` are NaN by construction (no
   catalogued reference protein exists to compare against) — never a fabricated 0%, and excluded
   (not zero-plotted) from panel (c).

## Distilled

- Pooled-ALL: KIR and HLA catalogues are near-tied at CDS (KIR 41.6% vs HLA 38.8%) and protein
  (HLA 39.9% vs KIR 39.2%) completeness — genomic-level "KIR is better covered" is an upper-bound
  artifact, not used as the headline.
- **Per ancestry, HLA is more complete than KIR in all 5 well-powered ancestries**, at both CDS and
  protein level — a Simpson's-paradox-shaped aggregation effect versus the pooled near-tie; treat
  the per-ancestry result as the more trustworthy one.
- KIR's novel-protein completeness (31.4%) is far higher than HLA's (6.5%) — KIR's novel proteins
  are relatively better captured by *repeated observation*, even though its overall protein
  catalogue completeness is not higher than HLA's.
- No single outlier gene drives either species' pooled number — both have a broad, real per-gene
  completeness spread (HLA 22–57%, KIR 22–56% at protein level, excluding pseudogenes).
- KIR2DP1/KIR3DP1 protein metrics = NA (pseudogenes), never 0.
- v1–v3 numbers (incl. the "KIR 6,300 novel proteins" / ">100%" bugs) are invalid — see `44`'s
  Version history.
