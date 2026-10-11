# 56 / 56b / 56c: phasing and typing quality benchmark

Scripts: `scripts/hla_popgen/56_phasing_typing_benchmark.py` (benchmark), `56b_relative_pair_refinement.py` (de-confounding), `56c_benchmark_figures.py` (figures, laptop-only). Tables: `*.tsv` in this folder (PUBLIC; counts 1-19 written `<20`).

## What and why
- Every claim (novel alleles, DQ purge, haplotypes, deletions) depends on typing accuracy and phasing accuracy.
- Earlier hints: platform shifts novelty; the MHC is split over about 4-5 contigs per haplotype; A~B is on one contig only 31% of the time.
- Goal: measure error with evidence independent of the caller, and build a per-call confidence score that lets each claim use only good-enough calls.

## How
**Truth sources available inside AoU**
- Twins / duplicates: same genome sequenced and assembled twice, so any difference is error.
- First-degree relatives: share at least one MHC haplotype unless they are IBD0 siblings. A gene with no shared allele is an error (or a rare recombination). A haplotype seen on one contig in one relative should appear intact in the other.
- Platform strata (Revio vs Sequel IIe) and same-contig flags.
- Not available: trio k-mer phasing (yak). No yak binary, the perimeter blocks downloads, and there is no IBS0 to identify trios.

**Confidence score (0-7)**: +1 each for same contig as phasing partner; at most 2 MHC contigs on that haplotype; template_distance at or below the 75th percentile; not homopolymer-indel-only; no template warning; no tie between candidate alleles; Revio. Pair-level score = minimum over the calls involved. Sweep = error among pairs/calls with score >= t.

**Two confounds in the first pass (56), fixed in 56b**
- Relatedness has only kinship (no IBS0), so parent-offspring and siblings cannot be separated, and IBD0 siblings inflate "no shared allele". Fix: keep only MHC-sharing pairs (at most 1 of 8 classical genes without a shared allele).
- The haplotype check penalised the relative's own fragmentation (A~B "77.5% inconsistent" was mostly A and B on different contigs in the relative). Fix: test haplotypes only where both relatives have the gene pair on one contig on both haplotypes.

## Results
Pairs: 569 first-degree (419 Revio-Revio, 38 mixed, 112 Sequel IIe-Sequel IIe), median kinship 0.251; twins/duplicates <20 pairs.

Pair sharing (genes with no shared allele): 0 -> 377 pairs; 1 -> 121; 2 -> 36; 3, 4, 5+ -> <20 each. The set is dominated by haplotype-sharing relatives (`fig_pair_sharing.png`).

**Typing sweep** (MHC-sharing pairs, all 17 genes; `fig_sweep_typing.png`)

| min score | pair-genes (n) | no shared allele % (95% CI) |
|---|---|---|
| 0 | 6,455 | 4.09 (3.63-4.60) |
| 2 | 6,389 | 3.79 |
| 3 | 6,188 | 2.96 |
| 4 | 5,884 | 2.02 |
| 5 | 4,517 | 1.48 |
| 6 | 2,338 (36% retained) | 1.15 (0.79-1.68) |

Monotone decrease, so the score captures error. At score 0, per gene: B 4.0%, DRB1 5.4%, TAP2 4.2%, MICA 12.3%, DRB3 13.3% (others masked, <20 events). At score >= 6, A, DQA1, DPA1, DQA2, TAP1, TAP2 have 0 events.

**Haplotype transmission** (both relatives same-contig; `fig_haplotype_transmission.png`)

| gene pair | n | not transmitted |
|---|---|---|
| all | 1,611 | 6.52% (5.41-7.83) |
| DQA1~DQB1 | 413 | <20 failures (<5%); 0 at score >= 6 (n=198) |
| DPA1~DPB1 | 435 | <20 failures; 0 at score >= 6 (n=205) |
| DRB1~DQB1 | 344 | <20 failures |
| B~C | 351 | 7.98% (5.58-11.29) |
| A~B | 68 | 54.4% (42.7-65.7) |

The haplotype sweep is non-monotone (6.5 -> 4.4 at >=4 -> 6.1 at >=6) because the A~B failures are not captured by the call-level score.

**Platform** (`fig_platform_quality.png`)

| | Revio-Revio | Sequel IIe-Sequel IIe |
|---|---|---|
| typing: no shared allele | 4.76% (n=4,688) | 2.07% (n=1,350) |
| haplotype not transmitted | 5.10% (n=1,059) | 10.34% (7.84-13.51; n=445) |

**Twins**: 96 novel calls in twin/duplicate pairs, <20 do not replicate in the co-twin, so at least 80% replicate. Too few pairs for gene-level twin rates.

## Conclusions
1. A per-call score from assembly features lowers relative-pair discordance about 3.5x (4.1% -> 1.15%) at 36% retention; score >= 4 halves error (2.0%) keeping about 91%, a candidate default.
2. Physical phasing is excellent within class II blocks (DQ, DP about 0 failures at high confidence), supporting the DQ G1/G2 result; acceptable for B~C (about 8%); unreliable for A~B (about 54%): chimeric or misjoined contigs across the roughly 1.3 Mb A-B interval, so physical A~B haplotypes should not be trusted.
3. Sequel IIe has about 2x the phasing failures; platform must be a stratum for haplotype claims.

## Caveats
- Selecting pairs with at most 1 failing classical gene conditions (bounds) the classical-gene rates; MICA, TAP, DRB3 are not conditioned.
- Rare MHC recombination (about 1-4% per meiosis over the class I-II span) sets a small floor on "error".
- Typing non-sharing is lower in Sequel IIe pairs, the opposite of the novelty odds ratio; likely ancestry/site confounding (Sequel IIe pairs come from a different collection). Not interpreted as better Sequel typing.
- Parent-offspring vs sibling cannot be split (no IBS0).
- Small cells are masked; the haplotype A~B estimate rests on n=68.

## Still missing
- Trio k-mer switch rate (yak): no binary, downloads blocked.
- Per-sample assembly QV.
- More twins/duplicates.
- A recombination-aware model.

## Files
`pairs_summary.tsv`, `twin_concordance_by_gene.tsv`, `novel_replication.tsv`, `transmission_consistency.tsv`, `threshold_sweep.tsv` (56, confounded first pass), `pair_sharing_distribution.tsv`, `refined_sweep.tsv` (56b, authoritative); `fig_sweep_typing.png`, `fig_haplotype_transmission.png`, `fig_pair_sharing.png`, `fig_platform_quality.png`.
