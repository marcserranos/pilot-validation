# 57 — Trio k-mer benchmark of phasing and typing (yak trioeval, parent-of-origin per gene)

Scripts: `scripts/hla_popgen/57_yak_trioeval_pilot.py` (regions), `57b_trio_gene_origin.py` (per gene),
`57c_trio_figures.py` (figures from the PUBLIC tables). VM 2026-10-06; yak 0.1-r69 from bioconda in `~/envs/yak`.
Aggregate PUBLIC output only; per-trio rows stay on the VM (`~/s05/internal/57*/`).

## Why
- 56/56b judged typing and phasing indirectly: relatives should share alleles and haplotypes.
- A trio is the direct truth. Parental k-mers say which parent each stretch of the child's haplotype contigs came from.
- yak trioeval is the standard metric of the HPRC and hifiasm papers.

## Trios
- Kinship only (no IBS0 in AoU): long-read child + two first-degree relatives (2^-2.5 ≤ kin ≤ 2^-1.5) that are not
  in the relatedness table as a pair (unrelated to each other).
- Any other configuration (parent + sib, parent + grandchild, two half-sibs) puts the pair in the table → these are
  parent–parent–child. Which parent is the father is irrelevant: switch and hamming are symmetric.
- `trio_selection.tsv`: 285 candidate trios → 181 with a child assembly and both parents' HiFi BAMs on the mount →
  179 distinct children (133 Revio, 46 Sequel IIe children).
- 57 used 25 per child platform (seed 57); 57b ran on all 179.
- Sanity: in child genes with both haplotypes resolved, the two copies almost always came from different parents
  (`gene_origin_summary.tsv`). A non-parent pair would give ~50%.

## How (57)
- Region design, to fit a 4-vCPU app:
  - k-mer DBs (`yak count -k31 -b30`) from the reads aligned to a padded window, for each parent and the child
    (MHC chr6:28.5–34.5 Mb, KIR chr19:54.5–55.1 Mb, all of chr20).
  - The child's hap1/hap2 hifiasm contigs are trimmed to the target with the AoU hifiasm→hg38 PAF.
  - `yak trioeval` is run on the trimmed pieces; `yak qv` uses the child's own window reads.
- Targets:
  - MHC 29.5–33.5 Mb
  - class I A..B 29.9–31.4 Mb
  - class II DR–DQ 32.4–32.7 Mb
  - KIR 54.7–54.9 Mb
  - chr20 5–60 Mb (genome background)
- Metrics (yak definitions):
  - **switch** = adjacent haplotype-specific k-mers on one contig that come from different parents;
  - **hamming** = haplotype-specific k-mers from the contig's minority parent.
  - Hamming grows with contig length when hifiasm (no trio, no Hi-C) joins phase blocks at random across
    homozygous stretches.
- `trioeval_summary.tsv`: per target × child platform, pooled % (all k-mers) and median/IQR across trios. QV is
  capped at 60. "% trios with zero error" was dropped: with 25–50 trios it encodes counts of 1–19.

## Results (57, 50 trios) — `trioeval_summary.tsv`, `fig_trioeval_regions.png`
| target | switch median % (IQR) | hamming median % (IQR) | Revio sw / hm | Sequel IIe sw / hm |
|---|---|---|---|---|
| class II DR–DQ | 0.19 (0.06–0.33) | 0.19 (0.08–0.66) | 0.31 / 0.33 | 0.13 / 0.14 |
| KIR | 0.34 (0–0.77) | 0.53 (0–2.46) | 0.39 / 1.18 | 0 / 0 |
| whole MHC | 0.46 (0.18–0.82) | 6.6 (2.6–12.5) | 0.68 / 5.1 | 0.17 / 9.7 |
| class I A..B | 0.42 (0.17–0.94) | 7.1 (0.9–16.4) | 0.83 / 9.5 | 0.17 / 3.5 |
| chr20 control | 1.56 (0.96–1.89) | 25.7 (16.4–32.2) | 1.75 / 17.8 | 0.96 / 31.6 |

- Class II is the best-phased region in the genome sample (switch and hamming ≈ 0.2%). This supports the
  DQ G1/G2 physical-phasing result.
- MHC switch error is below the chr20 background, as expected from high heterozygosity. MHC hamming error comes
  from the class I span (A..B hamming 7%, IQR up to 16%), i.e. contigs that join the parents' class I haplotypes
  → consistent with 56b's A~B failure.
- Platform:
  - Sequel IIe children have fewer, longer pieces (chr20 median 73 vs 183) → lower switch but higher hamming
    (chr20 32% vs 18%; MHC 9.7% vs 5.1%).
  - Longer contigs carry more long-range parental mixing. This is consistent with 56b's ~2× Sequel IIe haplotype
    failures.
  - Revio has more short-range flips.
- Median QV in the windows is 39–60. It is a window-read proxy, not a genome QV.

## How (57b)
- For every child HLA call (copy 1), take gene ±10 kb on its contig.
  - The call-table contig is an extraction range `ctg:S-E` and gene coordinates are relative to S; offset
    verified on 17,486 rows.
- Count paternal-only vs maternal-only k-mers (trioeval S lines, parental MHC DBs).
  - origin = pat/mat if one side has ≥5 k-mers and ≥3× the other, else ambiguous.
- `gene_origin_summary.tsv`: % resolved; % of children whose hap1 and hap2 copies came from different parents.
- `pair_phase_truth.tsv`: genes on one contig with different parents of origin (a physical-phasing error).
- `mendel_typing.tsv`: the parent the copy came from does not carry the child's allele (2-field, full), by 56's
  confidence score. This counts errors in the child's or the parent's call.

## Results (57b, all eligible trios: 179 children eligible; calls from ~172)
Figures: `fig_pair_phase_truth.png`, `fig_mendel_typing.png`.

- **Trio sanity** (`gene_origin_summary.tsv`):
  - 10–31% of calls per gene are ambiguous (too few parent-specific k-mers within ±10 kb; DQA2 31%, DPB1 25%).
  - Both copies from one parent: <20 children per gene out of 93–141 (0 for DPB1 and DQA2). These are real
    trios, and the rare same-parent cases are collapsed or duplicated haplotypes.
- **Physical phasing truth** (`pair_phase_truth.tsv`, "all" rows) — share of same-contig pairs whose genes come
  from different parents:

  | pair | different parents | 95% CI |
  |---|---|---|
  | A~B | 36/88 = 40.9% | 31.2–51.4 |
  | A~C | 32/81 = 39.5% | 29.6–50.4 |
  | DQB1~DPB1 | 37/174 = 21.3% | 15.8–27.9 |
  | B~C | <20/260 | < 7.7% |
  | DQA1~DQB1 | <20/261 | < 7.7% |
  | DRB1~DQA1 | <20/262 | < 7.7% |
  | DRB1~DQB1 | <20/253 | < 7.7% |
  | DPA1~DPB1 | <20/248 | < 7.7% |

  - The 50-trio subset had 0/69–76 for the class II pairs.
  - The A~B rate agrees with 56b's indirect 54.4% (42.7–65.7).
  - Platform rows are complementarily suppressed (`<20*`) so that "all" minus a platform cannot recover a masked
    count.
- **Mendelian typing** (`mendel_typing.tsv`) — child calls absent from the transmitting parent:

  | min score | 2-field, all 15 genes | full resolution | classical 8, 2-field |
  |---|---|---|---|
  | ≥0 | 4.15% (3.55–4.85; n = 3,611) | 12.6% | 4.09% |
  | ≥4 | 2.61% | 11.0% | 2.56% |
  | ≥6 | 2.34% | 9.3% | 2.45% |

  - The trio 2-field rate (4.15%) matches the relatives estimate (4.09%) almost exactly.
  - Full resolution at score ≥0, per gene: A 7.2% (lowest) to DQA1 and MICA 19.4%. 2-field per gene is masked
    except MICA 7.1%, DQA2 at score ≥6 (0/149) and DRB5 at score ≥4 (0/42).

### What this means
1. The DQ G1/G2 physical-phasing result rests on the best-phased pairs in the MHC: <7.7% trio mismatch at
   n ≈ 260, 0/69–76 in the first 50 trios, and ≈0.2% switch/hamming in the region.
2. Class I physical haplotypes spanning A (A~B, A~C) are wrong about 40% of the time. A~B and A~C haplotype
   claims need statistical phasing or trios; B~C is acceptable.
3. DQ–DP haplotypes are wrong about 21% of the time, so DQ–DP LD from physical phasing must be flagged.
4. 2-field typing error is about 4%, cut by about 40% at score ≥4. Third/fourth-field (non-coding) discordance is
   about 3× higher (12.6%), which matters for any non-coding-novelty claim. Part of it may be parent error rather
   than child error.

## Caveats
- Window-limited DBs: parental k-mers whose reads map outside the window (paralogs, KIR/MHC structural haplotypes
  absent from GRCh38) are missing. This adds false haplotype-specific k-mers in a symmetric way, so it inflates
  errors slightly. The class II ≈ 0.2% shows the floor is low.
- Pieces are PAF-trimmed alignment blocks, not contigs (pieces ≠ contig count).
- Platform is confounded with sequencing centre, hifiasm version and coverage.
- 50 trios; a few strata have wide IQRs (KIR).
- Trio = kinship inference, unverified by IBS0. The k-mer complementarity (both copies from one parent in <20 children per gene) confirms parentage.
- Masked mismatch counts are reported as an upper bound 19/n.
- Score thresholds are nested: threshold differences are call counts, not people.
- Mendelian mismatch mixes child and parent errors and needs the parent typed (HiFi parents only).
