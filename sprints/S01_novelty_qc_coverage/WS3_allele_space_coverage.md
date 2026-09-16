# WS3 — How much of the HLA allele space have we seen, and what would it take to see more? (L2 brief)

*Status: estimator library `_coverage.py` in progress → script 27 after WS1's script 24 produces
allele count tables.*

## Objective (plain language)
Marc's questions:
- "We have explored X% of the coding/protein space."
- "If we sampled this many people from Africa, this many from Europe and so on, how much would we
  have explored?"
- "How many more people does it take to reach X%?"

## Why the previous answer (04) can't be used as-is
- It counted flagged artifacts as novel alleles (81% of novel clusters) and pooled all non-coding
  novelty under CDS hashes. See WS1.
- It did not remove relatives.
- It estimated **richness** (the total number of alleles that exist). For a distribution with a
  huge tail of rare alleles, richness is notoriously unstable: every Chao2 row was
  LOW-CONFIDENCE.

## The better target: coverage
**Sample coverage** C is the fraction of the population's allele *frequency mass* belonging to
alleles already catalogued. Equivalently, it is the probability that the next person's haplotype
carries an allele we have already seen. It is estimated precisely from singleton and doubleton
counts (Good–Turing; Chao & Jost 2012), and it answers the practical question directly: "a
catalogue built from this sample would type X% of haplotypes in this population." Its complement,
1−C, is the "unseen mass". Richness stays a secondary, explicitly caveated number.

## Design
- **Resolutions:**
  - `protein`: 2-field for known calls; for novel alleles, the protein-sequence cluster of a clean
    `novel_protein` (script 24)
  - `cds`: 3-field / CDS cluster
  - `genomic`: added later from script 25 signatures if they are clean
- **Units:** haplotypes, **unrelated people only**. Flagged artifacts excluded (`clean_only`); the
  `all` variant is a sensitivity analysis.
- **Ancestry:** strict (admixture proportion ≥ 0.9) as primary; `ancestry_pred` as sensitivity.
- **Analyses (script 27 `27_allele_space_coverage.py`, VM, reads 24's local
  `allele_counts_by_resolution.tsv`):**
  1. Table per gene × resolution × ancestry: n haplotypes, S_obs, f1, f2, Ĉ with bootstrap CI,
     unseen mass, Chao1 (flagged secondary), and the size needed for 95/99/99.9% coverage, with an
     extrapolation-reliability flag (>2n).
  2. Coverage curves: rarefaction to n plus extrapolation to 2n, one line per ancestry, small
     multiples per gene, protein resolution as the main figure.
  3. **Cross-ancestry catalogue matrix:**
     - Build a catalogue from m haplotypes of ancestry A, with m = the smallest strict group's
       size, so the comparison is fair.
     - Measure the fraction of ancestry B's haplotypes it types.
     - Show it as a heatmap per gene and as the classical-gene mean.
     - The diagonal is "self" coverage. Off-diagonal loss shows how poorly a catalogue from one
       ancestry serves another.
  4. **Reference coverage today:** the fraction of each ancestry's haplotypes whose protein is
     already in IPD-IMGT/HLA (known 2-field or cds/protein_known). This is the reference-bias
     number at the level that matters clinically.
  5. **Sampling-design scenarios** at a fixed total N (the unrelated long-read cohort size):
     - the actual AoU LR composition
     - equal allocation across 6 ancestries
     - a EUR-dominated design (e.g. 94% EUR, the UK Biobank-like mix)
     - a greedy allocation that maximizes the population-weighted mean coverage

     Report per-target-ancestry and global coverage for each (subsampling without replacement,
     so N ≤ observed per group). Figure: grouped bars.
  6. **"How many more"** per ancestry and gene: the haplotypes needed to reach 99% and 99.9%
     protein-level coverage, against the current n. Extrapolation is only trusted up to about 2n;
     beyond that, give a lower bound and label it as such.
- **Outputs:** `reports/hla_popgen/27_allele_space_coverage/` holding the aggregate tables (no
  person counts <20 as bare numbers; allele-level counts are not written at all, only per-gene
  summaries), 4–5 figures, a plain-language README, and `summary.json`.

## Doubts to address in the report
- Coverage estimates assume haplotypes are independent draws. Relatives are removed, but
  population structure within an "ancestry" violates this mildly. Strict labels help.
- Novel-protein clusters depend on script 24's artifact labels. Report `clean_only` vs `all`.
- AoU is not a random sample of any population; state that coverage is "of the AoU-sampled
  population".

## Distilled (for SPRINT.md)
(pending results)

## Results — script 27, full cohort (2026-09-17)
Unrelated people, clean alleles, strict ancestry (admixture proportion ≥0.9).

**1. How much have we seen? Essentially all of it, at protein level.** Sample coverage per gene ×
ancestry is **99.0–99.8%** for the eight classical genes (mean by ancestry: AFR 98.9%, MID 99.2%,
SAS 99.2%, EAS 99.7%, AMR 99.7%, EUR 99.7%). CDS resolution is nearly identical. So a catalogue
built from this cohort already types about 99 of every 100 haplotypes drawn from the same
population — the classical HLA protein space is close to saturated at n≈12,000, and the older
"only 5–20% discovered" figure was a richness estimate inflated by artifacts and by pooled
non-coding clusters.

**2. The reference is not the problem at protein level.** 99.6–100% of haplotypes in every
ancestry already carry a protein present in IPD-IMGT/HLA. Reference bias in this data is a
*genomic-sequence* gap (field 4), not a protein gap.

**3. The real finding — catalogues do not transfer across ancestries.** Building a catalogue from
an equal number of haplotypes (~570, the smallest strict group) from one ancestry and using it to
type another, averaged over the 8 classical genes (protein resolution):

| catalogue from ↓ / types → | AFR | AMR | EAS | EUR | MID | SAS |
|---|---|---|---|---|---|---|
| AFR | 98.7 | 90.3 | 80.6 | **97.3** | 92.6 | 88.3 |
| AMR | 97.2 | 98.8 | 85.9 | 98.9 | 96.9 | 93.5 |
| EAS | 71.8 | 85.2 | 99.1 | 92.7 | 85.1 | 95.3 |
| EUR | **83.1** | 90.2 | 82.3 | 99.2 | 95.2 | 91.4 |
| MID | 95.4 | 91.2 | 79.7 | 98.1 | 100.0 | 93.8 |
| SAS | 77.1 | 87.1 | 91.0 | 95.7 | 92.1 | 99.2 |

The matrix is strongly **asymmetric**: an African catalogue types 97.3% of European haplotypes,
while a European catalogue of the same size types only 83.1% of African ones. African catalogues
are the best all-rounders; East Asian catalogues transfer worst to Africa (71.8%). This is a
direct, quantitative reference-panel design argument, and the cleanest answer to "how many people
from which ancestry would we need".

**4. Sampling designs at the full cohort size are uninformative** — every design (actual, equal,
94% EUR, greedy) reaches ~100% global coverage at n≈22,000 haplotypes, because the classical
protein space is saturated. The informative comparison is the equal-size matrix above; a budget
curve (coverage vs catalogue size per design) is the natural next figure.

**5. Where coverage is NOT saturated:** the non-classical genes — DPB2 90.3%, TAP2 99.0%,
MICA 99.1%, HFE 99.3%, DRB5 99.3%, TAP1 99.3%. Same conclusion as WS1: the undiscovered MHC is
outside the classical genes.

### Distilled (for SPRINT.md)
- Protein-level coverage of the classical genes is 99.0–99.8% per ancestry: near-saturation, not
  "5–20% discovered".
- 99.6–100% of haplotypes carry an IMGT-catalogued protein; the reference gap is genomic, not protein.
- Equal-size catalogues transfer asymmetrically: AFR→EUR 97.3% vs EUR→AFR 83.1%; EAS→AFR 71.8%.
- Non-classical genes (DPB2, TAP1/2, MIC, HFE, DRB5) are the unsaturated part.

## Sensitivity: does the artifact filter drive the near-saturation? (2026-09-17)
Re-ran script 27 with `--exclude-flagged all`, i.e. counting every flagged (mostly homopolymer)
allele as real. Mean protein coverage per ancestry falls from **99.3% to ~97.6%**
(AFR 97.6, AMR 97.8, EAS 97.8, EUR 97.4, MID 97.5, SAS 97.8); the weakest genes drop to ~96.4%.
Direction and conclusion are unchanged — coverage is high either way and far above the old
"5-20% discovered" framing — but the exact figure is filter-dependent, because flagged alleles are
overwhelmingly singletons and singletons are what coverage is most sensitive to. **Report both
numbers in the paper.** Files: `sensitivity_all_alleles_summary.json`,
`sensitivity_all_alleles_coverage_by_gene_ancestry.tsv`.


## Correction (2026-09-17): coverage range and catalogue size
- Protein coverage across the 8 classical genes spans **97.07-99.84%**, not "99.0-99.8% in every
  ancestry": MID falls below 99% for four of eight genes (HLA-B 97.07%, HLA-A 97.74%, DPB1 98.44%,
  DQB1 98.79%). Every other ancestry's minimum is >= 99.03%. MID is also the thinnest group, so this
  is partly sample size.
- The cross-ancestry catalogue size is **566-585 haplotypes, set per gene by the smallest group**,
  not a single 575. Describe it as "equal size per gene (~570)".
