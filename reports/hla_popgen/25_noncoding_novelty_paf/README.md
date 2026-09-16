# Depth-4 ('beyond CDS') novelty, split by the real non-coding sequence

## In plain words

A depth-4 call such as `HLA-A*01:01:01:new` means the person's coding sequence is already known, but something outside it (an intron or UTR) differs from every catalogued allele. The earlier novel-allele table grouped these calls only by their coding sequence, so one 'novel allele' could hold thousands of different non-coding sequences. This report compares each call's contig with its template allele and writes the exact differences as a signature. Calls with identical signatures are grouped together.

- Depth-4 haplotypes examined: **32207**; with a readable difference string: **31886**. Excluded before clustering: **0** split alignments, **300** below the --min-qcov floor (0.98).
- PAF selection: **9805** secondary (`tp:A:S`) records were discarded before selection; **0** calls needed multiple primary records chained into one alignment; **0** could not be chained (conflicting/overlapping primary records) and were excluded.
- Distinct signatures: **15112**. Former clusters split: **708** former clusters contain a median of **1.0** distinct signatures using each cluster's single most common (modal) template (max 287); using the raw, template-confounded count instead the median is **2.0** (max 788) -- the modal-template number is the honest headline (see Caveats).
- Fraction of depth-4 haplotypes whose only differences are homopolymer indels (the typical long-read error pattern): **0.2126**.
- Homopolymer context was checked against the real allele sequence for **1.0** of haplotypes. For the rest, the flag only means 'the inserted or deleted bases are all one letter', which also catches every 1-bp indel, so it over-counts.

Counts from 1 to 19 are written as `<20`. Signatures are listed only when they have at least 20 unrelated carriers. Any ratio whose denominator count is below 20 is blanked (`NA`) rather than shown, in every committed table, `summary.json`, and figure.

## How the difference string is read

minimap2 aligned every known genomic allele (query) to the person's contig (target). In its cs string, `*ab` means contig base a and allele base b, `+seq` means bases the contig lacks, and `-seq` means extra bases in the contig. Events are written on the known allele's own 1-based coordinates as `posA>G`, `posdelX` or `p_p+1insX`. For alignments on the reverse strand, positions are counted back from the alignment end and the bases are reverse-complemented, so the same difference gives the same signature on either strand.

## Per gene

| gene | n_depth4_haplotypes | n_parsed | n_no_row | n_cs_error | n_split_alignment | n_low_qcov | n_multi_record_chained | n_fallback_3field | n_context_sequence | n_minus_strand | n_signatures | n_signatures_ge20_unrelated | n_former_novel_ids | frac_homopolymer_only | median_qcov | n_intron_only | n_involves_utr | n_touches_cds | n_unannotated | n_other | n_no_difference |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HLA-A | 758 | 743 | 0 | 0 | 0 | <20 | 0 | 0 | 743 | 402 | 620 | <20 | 74 | 0.341 | 1.000 | 517 | 88 | 52 | 0 | 0 | 86 |
| HLA-B | 2344 | 2309 | <20 | 0 | 0 | 33 | 0 | <20 | 2309 | 1286 | 1283 | <20 | 134 | 0.216 | 1.000 | 524 | 949 | 363 | 0 | 0 | 473 |
| HLA-C | 1720 | 1679 | 0 | 0 | 0 | 41 | 0 | 0 | 1679 | 871 | 1136 | <20 | 79 | 0.232 | 1.000 | 585 | 560 | 231 | 0 | 0 | 303 |
| HLA-DPA1 | 2707 | 2679 | <20 | 0 | 0 | 25 | 0 | <20 | 2679 | 1427 | 1408 | <20 | 56 | 0.140 | 1.000 | 1261 | 960 | 288 | 0 | 0 | 170 |
| HLA-DPB1 | 4548 | 4460 | <20 | 0 | 0 | 78 | 0 | <20 | 4460 | 2001 | 2690 | <20 | 110 | 0.172 | 1.000 | 3677 | 527 | 243 | 0 | 0 | <20 |
| HLA-DQA1 | 3524 | 3491 | <20 | 0 | 0 | 30 | 0 | <20 | 3491 | 1640 | 1689 | <20 | 47 | 0.389 | 1.000 | 2304 | 742 | 422 | 0 | 0 | 23 |
| HLA-DQB1 | 2280 | 2253 | <20 | 0 | 0 | 25 | 0 | <20 | 2253 | 1172 | 1531 | <20 | 62 | 0.379 | 1.000 | 1753 | 396 | 94 | 0 | 0 | <20 |
| HLA-DRB1 | 14326 | 14272 | <20 | 0 | 0 | 53 | 0 | <20 | 14272 | 7634 | 4755 | 111 | 146 | 0.160 | 1.000 | 12745 | 411 | 1115 | 0 | 0 | <20 |

## Former clusters by number of distinct signatures (modal template, the headline)

| n_distinct_signatures_bin | n_former_clusters |
|---|---|
| 1 | 382 |
| 2 | 49 |
| 3-5 | 63 |
| 6-10 | 44 |
| 11-50 | 131 |
| 51-100 | 30 |
| 101-1000 | <20 |
| >1000 | 0 |

## `no_difference` calls by alignment-coverage decile

Motivates the `--min-qcov` floor: a call that only partly covers the allele can wrongly look identical to it. Computed before the floor excludes low-coverage calls from clustering.

| qcov_decile | n_calls | n_no_difference | frac_no_difference |
|---|---|---|---|
| (-0.001, 0.1] | <20 | 0 | nan |
| (0.1, 0.2] | 0 | 0 | nan |
| (0.2, 0.3] | 0 | 0 | nan |
| (0.3, 0.4] | 0 | 0 | nan |
| (0.4, 0.5] | 0 | 0 | nan |
| (0.5, 0.6] | 0 | 0 | nan |
| (0.6, 0.7] | 0 | 0 | nan |
| (0.7, 0.8] | 0 | 0 | nan |
| (0.8, 0.9] | 0 | 0 | nan |
| (0.9, 1.0] | 32185 | 1175 | 0.037 |

## Homopolymer-only fraction (all genes, strict ancestry, unrelated only)

| ancestry | n_haplotypes | n_homopolymer_only | n_snv_only | n_snv_plus_indel | n_indel_nonhomopolymer | frac_homopolymer_only |
|---|---|---|---|---|---|---|
| POOLED | 30962 | 6568 | 8483 | 4896 | 9980 | 0.212 |
| AFR | 6998 | 1345 | 2223 | 1259 | 1893 | 0.192 |
| AMR | 3771 | 743 | 983 | 661 | 1265 | 0.197 |
| EAS | 3966 | 1028 | 1238 | 646 | 954 | 0.259 |
| EUR | 2524 | 589 | 467 | 346 | 1003 | 0.233 |
| MID | 789 | 148 | 211 | 132 | 282 | 0.188 |
| SAS | 3385 | 622 | 1190 | 432 | 1051 | 0.184 |

## Where the differences are

Each event is placed on the template allele's own map from `alleles.csv.gz` (UTR5, exonN, intronN, UTR3). A depth-4 call has a known coding sequence, so its differences should sit in introns or UTRs. Haplotypes whose differences touch the coding sequence point either to a template that is not the allele the call matched, or to a pipeline inconsistency.

- Share of parsed haplotypes by region class: `{"intron_only": 0.7328, "involves_utr": 0.1453, "touches_cds": 0.0881, "no_difference": 0.0338}`
- Touching CDS: `{"n": 2808, "n_template_shares_call_3fields": 25, "n_template_differs_at_3fields": 2783, "by_selection_path": {"template": 2808}, "n_signatures": 1343}`

## Cohort coverage

- **0** haplotypes belong to a person absent from `cohort_membership.tsv`; their relatedness status is recorded as `NA` (a third state, not silently treated as related or unrelated).
- Cluster-level `signature_class`/`region_class` are taken by majority vote across member calls; disagreements are recorded per cluster (`n_signature_class_disagree`, `n_region_class_disagree`) -- they are expected to be zero, since both are pure functions of the signature.

## Caveats

- The signature is keyed on the template allele. Two identical sequences with different templates count as two raw signatures; the modal-template number used as the headline above avoids this confound by restricting to each cluster's single most common template.
- Differences outside the aligned part of the allele are not seen (see `median_qcov`); calls covering less than 0.98 of the allele are excluded from clustering entirely.
- Without the allele sequence, an indel inside a repeat may be placed differently on the two strands and split one signature into two.
- Former-cluster mapping skips genes with more than one copy on a contig (same rule as 03).
- Multi-record chaining assumes non-overlapping, order-consistent primary alignments; anything else is excluded as `split_alignment` rather than guessed at.

Run: genes=['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRB1'], limit=None, homopolymer_min_run=3, min_qcov=0.98, threads=3.
