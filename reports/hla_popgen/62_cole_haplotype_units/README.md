# 62 — Phased class II haplotype units: data dictionary

Script: `scripts/hla_popgen/62_cole_haplotype_units.py`. This file is the data dictionary
only; it contains no data. The per-person table (`haplotype_units.tsv.gz`) is INTERNAL:
it lives on the VM under `~/s05/internal/62/` and in the in-perimeter collaborator bucket,
never in this repo.

## Grain
One row per (`person_id`, `hap`, `unit`); `unit` is `DQA1~DQB1` or `DRB1~DQA1~DQB1`.
Every haplotype present in HLA Table 1 gets both units, including incomplete ones.
The file's first line is `# INTERNAL ...`; read with `comment.char = "#"` (R) or
`comment="#"` (pandas).

## Columns
| column | type | meaning |
|---|---|---|
| `person_id` | str | AoU research id. INTERNAL table only; never in PUBLIC outputs. |
| `hap` | str | hap1 / hap2: which assembled haplotype (hifiasm, no trio/Hi-C). |
| `unit` | str | DQA1~DQB1 or DRB1~DQA1~DQB1. Every (person, hap) present in Table 1 gets both rows. |
| `allele_<GENE>` | str | Full Immuannot consensus call of DRB1 / DQA1 / DQB1 on this hap (all copies, '/'-joined, copy order). Empty if no call or gene not in the unit. May contain 'new' (novel at that field). |
| `two_field_<GENE>` | str | Same call truncated to 2 fields (e.g. DQB1*06:02). Empty when field 1 or 2 is novel/undetermined. |
| `unit_full` | str | '~'-joined full calls of the unit's genes (copy 1). Empty if incomplete. |
| `unit_two_field` | str | '~'-joined 2-field calls (copy 1). Empty if incomplete or any gene unresolved at 2 fields. Typing claims are supported at 2 fields (57b trio: ~4% Mendelian error, 12.6% at full resolution). |
| `n_contigs` | int | Distinct assembly contigs carrying the unit's genes (copy 1) on this hap. |
| `same_contig` | bool | All unit genes (copy 1) on one contig. |
| `pair_evidence` | str | Per gene pair of the unit: 'GENE~GENE:<same_contig|split|missing>/<allowed|flagged|not_allowed>', ';'-joined. |
| `evidence_level` | str | physical | physical_flagged | not_allowed_by_rule | split_contigs | incomplete (rule in the script docstring). Use `physical` for haplotype-level claims. |
| `dq_g1g2` | str | Petersdorf 2022 DQ group of the copy-1 DQA1/DQB1 pair: G1 | G2 | predicted_incompatible | unclassified | not_called. |
| `dq_spacing_kb` | float | |gene_start(DQB1) - gene_start(DQA1)| in kb when both (copy 1) are on one contig; empty otherwise. |
| `dq_structure` | str | canonical | anomalous_spacing | anomalous_multicopy | anomalous_spacing;anomalous_multicopy | not_assessable. |
| `any_novel` | bool | Any call of the unit's genes contains 'new'. |
| `multi_copy` | bool | Any unit gene has copy_index > 1 on this hap. |
| `mhc_rel_depth` | float | Long-read depth over the MHC window / the person's genome-wide control depth (58, raw). Empty if 58's table was not read or the person is not in it. |
| `coverage_scale` | str | Scale of coverage_value and coverage_cutoff: norm (depth / cohort median of the MHC window on the person's platform; 1.0 = typical) or abs (raw rel_depth). |
| `coverage_value` | float | The person's MHC-window depth on coverage_scale, as compared with the cutoff. Empty if not measured. |
| `coverage_flag` | str | pass | low_coverage | not_measured. |
| `coverage_cutoff` | float | Cutoff used for coverage_flag (on coverage_scale). Empty when 59 proposed no MHC cutoff. |
| `coverage_cutoff_status` | str | 'set_by_59': cutoff taken from 59's exclude_people.tsv (its MHC-window proposal); 'provisional_until_59': fallback cutoff on 58's depth, a placeholder. |
| `platform` | str | revio / sequel2e (sequel2 people have no HLA calls yet). Platform affects novelty rates (50c); use as a covariate. |
| `ancestry_pred` | str | Predicted genetic ancestry (AoU), upper-cased. |
| `unrelated` | str | True / False from the greedy kin >= 0.0442 rule over HLA-called people; not_assessed when no relatedness table. |

## Rules and parameters
- Pair rules (`--pair-rules`): `DQA1~DQB1=allowed,DRB1~DQA1=allowed,DRB1~DQB1=allowed,DPA1~DPB1=allowed,B~C=allowed,DQB1~DPB1=flagged,DQA1~DPB1=flagged,A~B=not_allowed,A~C=not_allowed`; unlisted pairs: `not_allowed`. Source: 57b trio truth (same-contig pairs whose
  genes came from different parents: class II pairs and B~C <7.7%, DQB1~DPB1 21%,
  A~B / A~C ~40%). Proposed S05 rule, not yet ratified (context/DECISIONS.md).
- DQ spacing cut (`--dq-spacing-max-kb`): 100.0 kb, measured between
  `gene_start` coordinates on the shared contig (50e follow-up: 98% of same-contig
  haplotypes lie in 16.6–24.2 kb).
- Coverage, preferred source: 59's INTERNAL `exclude_people.tsv` (`--exclude-people`). Its
  `window_MHC` rows give each person's MHC depth on 59's scale and 59's proposed MHC cutoff;
  `coverage_flag = low_coverage` below it, status `set_by_59`. 59's `exclude` column (which
  also reflects the KIR window) is not used here.
- Coverage, fallback (no 59 file): 58's depth over `window_MHC` at MAPQ ≥ 0, scale `norm` (norm = rel_depth / the cohort
  median on the person's platform, platforms with < 20 people use the all-platform median),
  cutoff 0.6, status `provisional_until_59`.
- The status of the source actually used is in `coverage_cutoff_status` on every row.
- Unrelated: greedy removal at kin ≥ 0.0442 among HLA-called people.

## Recommended use
- Haplotype-level association: `evidence_level == "physical"`,
  `dq_structure == "canonical"`, `coverage_flag != "low_coverage"`, `unrelated == True`.
- Use `unit_two_field` as the default allele resolution; full resolution needs read support.
- Keep `platform` (and ancestry) as covariates.

## Known limits
- The 991 `sequel2` people (mostly AFR) have no HLA calls yet (Tier-3 trim; S05 WS-K).
- Without 59's file the coverage cutoff is a placeholder (`provisional_until_59`).
- `split_contigs` haplotypes keep their hap label, but the label is not reliable across a
  contig break (50e).

## Companion PUBLIC files (`~/s05/results/62/`)
- `unit_freq_by_ancestry.tsv`: 2-field unit frequency by ancestry among unrelated, physical,
  not-low-coverage haplotypes; units named only with ≥ 20 haplotypes (smaller ones pooled as
  `other`), counts 1–19 shown `<20`.
- `unit_qc_summary.tsv`: haplotype counts per evidence level, DQ class, DQ structure, coverage
  flag, platform and unrelated flag, per ancestry. Counts 1–19 shown `<20`; a cell masked only
  so that it cannot be recovered by subtraction from the total is shown `<20*`.

## Upload (run by a person, after review; not done by the script)
```
gsutil -m rsync -n -r ~/s05/internal/62/ \
    gs://hla-calls-share-wb-cordial-leechee-9743/release_<YYYY-MM-DD>/haplotype_units/
gsutil -m rsync -r ~/s05/internal/62/ \
    gs://hla-calls-share-wb-cordial-leechee-9743/release_<YYYY-MM-DD>/haplotype_units/
gsutil -m rsync -n -r -c ~/s05/internal/62/ \
    gs://hla-calls-share-wb-cordial-leechee-9743/release_<YYYY-MM-DD>/haplotype_units/
```
