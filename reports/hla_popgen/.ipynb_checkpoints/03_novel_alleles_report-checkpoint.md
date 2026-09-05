# Novel HLA allele discovery -- Table 3 (`novel_alleles.tsv`)

Methodology: scripts/hla_popgen/research/NOVEL_LIT.md. Novel-allele candidates = Immuannot's `new`-tagged consensus calls (Table 1), clustered by exact observed-CDS sequence identity (sha1). Recurrence in >=2 unrelated persons is the primary evidence a candidate is real biology rather than an assembly artifact (NOVEL_LIT.md section 1.4 item 5) -- this is the single hardest gate in `passes_qc`.


## Matching Table 1 -> `cds.fa.gz` sequences

| Metric | Count |
|---|---|
| n_novel_table1_rows | 285814 |
| n_matched | 285661 |
| n_ambiguous_copy | 152 |
| n_missing_cds | 1 |

Ambiguous-copy rows (`n_ambiguous_copy`) are genuinely unresolvable per reference/IMMUANNOT_GTF_SPEC.md part D (the detection-order join key `tmp.gene.csv` is deleted by Immuannot's own cleanup) -- excluded from Table 3 by design, not a bug.


## `template_warning` QC policy for this run

Disqualifying warning token(s) used in this run's `passes_qc` gate: **inframe_stop, partial_CDS**. Per SCHEMA.md's `template_warning` policy, mere presence of ANY warning is NOT a gate -- Immuannot writes the literal string `NA` to mean *no warning* on 57.4% of real transcript rows, so the true warning rate is ~38%, not ~95% (a real warning describes whether the template's CDS could be cleanly reconstructed, not whether the typing call is wrong). Only the token(s) listed above disqualify a candidate; every other token is treated as benign by this run. Override with `--disqualifying-warnings`.


Warning-token breakdown among the 285661 matched novel candidates (a candidate may carry more than one token, so counts need not sum to the total; 141684 carried at least one token):

| warning token | n candidates | disqualifying in this run? |
|---|---|---|
| partial_CDS | 107747 | YES |
| no-stop_codon | 30390 | no (benign by default) |
| no-start_codon | 29190 | no (benign by default) |

## Cluster-level summary

Headline counts below cover only clusters with a RESOLVED gene-level identity (`novelty_class != "undetermined"`). Depth-1 (`undetermined`) clusters -- even the gene-level field unresolved -- are reported separately immediately after, per SCHEMA.md's Fix 1 steer: an allele whose gene identity is itself unresolved is not a defensible "novel allele" claim.

| Metric | Count |
|---|---|
| Distinct novel-allele clusters, resolved identity (`novel_id`) | 21463 |
| Passing all QC gates (`passes_qc`), resolved identity | 1190 |
| Singleton (1 person only) -- excluded by the recurrence gate | 18760 |
| Flagged homopolymer-indel-only artifact | 13735 |

**`undetermined` (gene-level-unresolved) clusters: 75** -- included in the Table 3 TSV, EXCLUDED from every headline count above. Of those, 6 also pass the recurrence/warning/homopolymer gates (i.e. would look like real novel alleles by every gate except gene-level resolution) -- reported here explicitly rather than folded into either total.


(Total clusters across both categories: 21538.)


## Synonymous vs non-synonymous breakdown (`novelty_class`, resolved only)

| novelty_class | n clusters | % |
|---|---|---|
| protein_altering | 19486 | 90.8% |
| beyond_cds | 1300 | 6.1% |
| synonymous | 677 | 3.2% |

A strong excess of `protein_altering` (non-synonymous) clusters over what a uniform random-error process would produce is itself evidence of real biology under balancing selection (Zhou et al. 2024's own validation argument, NOVEL_LIT.md section 0) -- this table reports the raw breakdown; 04_allele_saturation.py's neutral-model comparison is where that claim gets a formal statistical treatment.


## By gene_class (resolved-identity clusters only)

| gene_class | n clusters | passing QC |
|---|---|---|
| pseudogene_I | 5093 | 158 |
| mic_tap | 4537 | 293 |
| classical_II | 2807 | 238 |
| classical_I | 2617 | 176 |
| class_II_paralog | 2158 | 113 |
| class_II_accessory | 2123 | 125 |
| nonclassical_I | 2090 | 76 |
| other | 38 | 11 |

## Ancestry gradient in raw novel-call rate (Table 1, key integration assertion)
Fraction of `(person, hap, gene)` calls tagged novel, by ancestry. IPD-IMGT/HLA's European bias predicts non-European ancestries should show a materially higher novel rate -- against synthetic fixtures this must recover the encoded gradient (afr 0.22 vs eur 0.05); against real data this is the paper's central empirical claim.

| ancestry | novel_rate | n_calls | n_novel |
|---|---|---|---|
| AFR | 0.380 | 226680 | 86119 |
| EAS | 0.345 | 108395 | 37408 |
| SAS | 0.316 | 92599 | 29257 |
| MID | 0.307 | 36129 | 11100 |
| AMR | 0.303 | 202820 | 61531 |
| EUR | 0.267 | 224541 | 59903 |

Outputs: `/home/jupyter/repos/pilot-validation/reports/hla_popgen/novel_alleles.tsv`, `/home/jupyter/pipeline_outputs/novel_alleles_seqs.fa` (VM-only, sequences never leave pipeline_outputs)
