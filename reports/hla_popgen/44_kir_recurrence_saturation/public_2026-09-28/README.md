# 44 PUBLIC tables from the 2026-09-28 run (disclosure-layer version)

Pulled 2026-10-11 from the VM, `~/s04/results/44_public/`, md5-verified against the VM copy.
They are the tables that were left on the VM when access was lost on 2026-09-28 (S04 LOG).

The parent folder holds the 2026-09-27 v4b tables (`~/s04/results/44_v4b/`), which the committed
44/46 figures and README numbers were built from. Of the 10 tables in `44_public`, 6 are
byte-identical to the parent folder (`artifact_qc`, `coverage_chao2`, `diagnostics_identity`,
`equal_n_slope`, `genomic_identity_qc`, `kir_protein_catalogue_qc`) and are not repeated here.
These 4 differ:

| Table | Difference from the parent-folder version |
|---|---|
| `recurrence_classes.tsv` | adds the `review_flag` column (strata with <100 carriers); same masking of `<20` cells |
| `saturation_curves.tsv` | full recompute (same columns); per-step means differ slightly from v4b (resampling) |
| `genomic_artifact_qc.tsv` | row order differs; includes a `KIR3DL3 partial_cds` row |
| `kir_cds_match_qc.tsv` | column order differs; `n_protein_gene_uncovered_calls` is an exact count (<20, a count of calls, not people) where the parent version shows `<20` |

Nothing here is per person: these are counts and means by gene × ancestry × level.
`genomic_artifact_qc.tsv` has some `n_calls` cells below 20 (calls per gene and label, not people).
Under DECISIONS 2026-10-11, the <20 rule applies at publication (figures, manuscript, slides).

Which version is canonical for the paper is still open. Until it is decided, the parent-folder
(v4b) tables remain the source of the committed 44/46 figures.
