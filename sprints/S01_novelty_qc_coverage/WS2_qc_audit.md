# WS2 — Load-bearing QC: where does the 5% relative discordance come from? (L2 brief)

*Status: audit done → implementation (script 26).*

## Objective (plain language)
The paper will rest on "long-read calls are accurate." Our evidence was "545/545 relative pairs have
zero phase switches" and "~95% of genes agree between relatives." Marc's question: if phasing is
perfect, what is the 5%? We find out exactly, and replace the circular pieces with non-circular
tests.

## Audit findings (2026-09-16, `11/16/17` code + committed outputs)
1. **Pair selection is circular.** A pair is "high_sharing" if ≥80% of its genes share an allele
   (`16:185-190`), and the error rate is then the fraction of genes in those same pairs that don't
   share (`16:237-256`). Pairs are selected on the statistic being reported. At the pair level this
   mostly removes IBD0 siblings, so the bias is modest but unquantified.
2. **A mismatch is phase-agnostic set intersection at 2-field** (`16:161-165`). Genes missing in
   either person are **dropped from the denominator** instead of counted. Allele dropout and
   fragmentation are therefore invisible. `copy_index>1` rows are excluded.
3. **"0 switches" is structurally weak.**
   - A switch is only counted between consecutive *resolved* genes, meaning uniquely shared and
     unambiguous on both sides, and only on the same contig for both people (`16:200-224`).
   - Mismatched genes never get a state, so a localized error can never appear as a switch.
     The two statistics measure disjoint failure modes; they don't corroborate each other.
   - The number of testable transitions per pair is never reported, and 80% of haplotype files
     span more than one contig, so the test's power is unknown.
   - **Verdict: "0/545" is not yet evidence of correct phasing. It needs a denominator and an
     empirical power estimate.**
4. **Raw-sequence divergence is bimodal** (17, 19,700 comparisons, 1,020 divergent).
   - 39% are exactly 1 base and 57% are ≤3 bases (point-error-like).
   - 36% are ≥5 scattered blocks (a genuinely different allele).
   - Classical genes: 228/4,354 divergent (5.2%). Only 28% of those are 1-base; 47% are ≥5 blocks.
     The median number of differing bases runs from 2 (C) to 19 (DRB1).
   - The detail file has no pair id and no zygosity, so concentration can't be checked.
5. **Implied per-base error if 1-base divergences are assembly errors:** classical genes
   ≈1.8×10⁻⁵ (**QV ≈ 47**), all genes ≈2.5×10⁻⁵ (QV ≈ 46). This is inside the typical HiFi
   assembly range (Q40–50), so the point-error half of the discordance is fully explained by normal
   consensus error.
6. **The many-block half is the real concern.** Candidates:
   - allele dropout (a collapsed heterozygous haplotype makes a person look homozygous)
   - a wrong relationship or IBD state
   - gene conversion
   - paralog mis-assignment (MICA, HLA-H/K/U at 12–14%)

## Implementation spec — script 26 `26_qc_relatives_v2.py` (VM)
Inputs: Table 1, cohort_membership, per-haplotype `cds.fa.gz`, the relatedness TSV, the SR genotype
file (`--sr-genotypes`, as in 07), and an optional `--yob-tsv` (person_id, year_of_birth).
Aggregate outputs only.

1. **Pair universe without selection:** all long-read pairs with kin ≥ 0.177 (first degree), plus
   duplicates/MZ twins (kin > 0.354). Record the kin value.
2. **IBD state per pair, not selected on:**
   - Per gene, compute the phase-agnostic share of the observed CDS **sequence** multiset (exact
     sequence, not name).
   - Pair-level sharing fraction over classical + class-II-accessory genes.
   - Plot the histogram of sharing fraction for first-degree pairs. Expect modes near 0 (IBD0
     siblings) and near 1.
   - Fit a 2-component beta/binomial mixture, or use a sweep of thresholds (0.5, 0.8, 0.9, 0.95),
     and report the discordance rate at each. That shows the sensitivity instead of hiding it.
   - If `--yob-tsv` is given, flag age gap ≥15 y as likely parent-child and report its rates
     separately.
3. **Discordance decomposition** (IBD≥1 pairs, genes with copy_index==1). Every gene-comparison is
   exactly one of:
   - `concordant`
   - `missing_one` (gene called in one person only; count, don't drop)
   - `missing_both`
   - `point_diff` (≤3 bases, one block; best pairing as in 17)
   - `dropout_candidate` (not concordant, and at least one person is homozygous for the observed
     CDS at that gene, i.e. both haps carry identical sequence or only one hap has the gene)
   - `fragmentation_candidate` (the gene sits on a contig that carries no other classical gene, in
     either person)
   - `other_multiblock`
   
   Precedence is in that order after `concordant`. Report per gene and per gene class, with Wilson
   CIs.
4. **Technical replicates:** duplicate/MZ pairs get full zygosity-aware equality of both haplotype
   sequences, per gene, plus per-base differences, giving a direct QV.
5. **Dropout detector, cohort-wide:**
   - Per classical gene and strict ancestry: observed homozygosity (identical CDS on both haps)
     vs the HWE expectation from the same ancestry's CDS frequencies, as an obs/exp ratio.
   - The same for the **SR 2-field calls of the same people**, and for LR at 2-field, so the two
     compare like for like.
   - LR excess over SR suggests collapse.
   - Also compare homozygosity by fragmentation: the number of distinct contigs that carry the 8
     classical genes on each hap (1, 2 or ≥3).
6. **Switch-test power:**
   - Recompute states as in 16, but report the per-pair number of resolved genes and testable
     same-contig transitions (distribution and total).
   - **Empirical power:** for each IBD1 pair, simulate a switch by swapping hap labels of one
     person's genes downstream of a random breakpoint (within a contig), re-run the detector, and
     report the detection rate.
7. **Orthogonal SR concordance** for every person in both cohorts:
   - Per gene, LR 2-field vs SR 2-field genotype concordance (0, 1 or 2 alleles matching), by
     strict ancestry and by LR field class (known vs f2/f3/f4 novel).
   - LR f4 novelty should *not* reduce concordance; f2 novelty should.
8. **Outputs:**
   - `qc_relatives_v2_report.md` in plain language
   - `discordance_decomposition.tsv`
   - `sharing_histogram.png`
   - `decomposition_by_gene.png` (stacked bars)
   - `homozygosity_obs_exp.png` (LR vs SR, per gene × ancestry)
   - `switch_power.png`
   - `sr_lr_concordance.tsv` / `.png`
   - `summary.json`

   Never write person ids; counts 1–19 are written as `<20`.

## Distilled (for SPRINT.md)
- The "0/545 switches" result is not yet informative: errors are excluded before switches are
  counted, and there is no denominator. Script 26 adds the denominator and empirical power.
- Half of relative discordance is 1–3-base differences, consistent with HiFi consensus error
  (implied QV ≈ 46–47). The other half is multi-block; dropout, relationship and paralog
  explanations are under test.
