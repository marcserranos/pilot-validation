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
