# Status — live session state (RNA-seq / repertoire workstream)

> **Role:** where we are *right now*, plus the literal next commands. The only file fully
> rewritten each session. Anything durable graduates to ENVIRONMENT, DECISIONS or EXPERIMENTS.

## As of 2026-09-28 — report 01 done; reports 02–04 written and tested, not yet run on the VM

**Scope** (Aleix, 2026-09-27; widened 2026-09-28 to "Nature Comms-ready, with side
approaches that move the needle"): the TRB repertoire, its SCEPTR embedding, and validation
against known biology. HLA and disease association are still the next phase.

**Report 01** (embedding maps) is done and committed with its README.

**Reports 02–04** have been built and tested end to end on synthetic data with planted
effects. The effects were recovered, and every figure was inspected.
- **02 `repertoire_atlas`:** full-repertoire recovery, TRBV usage, CDR3 length, and the
  sharing spectrum. It also tests rarefied diversity and clonal expansion against age
  (Britanova 2014).
- **03 `publicness_and_robustness`:**
  - OLGA generation probability against publicness (Elhanati 2018).
  - Robustness of the person vector to the top-k cap and to read weighting.
  - Cross-validated ancestry signal, SCEPTR vectors against TRBV usage.
- **04 `antigen_specificity`:**
  - VDJdb epitope benchmark against 3-mer and V-only baselines.
  - Cohort matching to known CMV, EBV, flu and SARS-CoV-2 TCRs.
  - CMV carriage against age, with a low-Pgen sensitivity model.
  - CMV epitopes by restricting HLA allele, which is the bridge to the HLA phase.
- **Shared helpers:** `scripts/repfig.py`. **Runner:** `scripts/run_reports_02_04.sh`.

## Pick up here

1. On the main VM, run all three in order (~30–45 min):
   ```bash
   cd ~/repos/pilot-validation && git pull
   nohup bash aleix/RNA-seq/scripts/run_reports_02_04.sh > ~/reports_02_04.log 2>&1 &
   tail -5 ~/reports_02_04.log
   ```
2. When the log ends with `all reports done`, copy the de-identified outputs (figures and
   CSVs, no `_cache`) into the repo and push:
   ```bash
   cd ~/repos/pilot-validation && for r in 02_repertoire_atlas 03_publicness_and_robustness 04_antigen_specificity; do mkdir -p aleix/RNA-seq/reports/$r && cp ~/pipeline_outputs/rnaseq/reports/$r/*.{png,pdf,svg,csv} aleix/RNA-seq/reports/$r/; done && git add aleix/RNA-seq/reports && git commit -m "Reports 02-04: VM outputs" && git push
   ```
3. Locally: read the real figures and summaries, then write each report's README and
   `reports/COMPREHENSIVE_REPORT.md` (the narrative across 01–04). Append EXPERIMENTS.md.

**Backup still owed:** the embeddings live only on the main VM's disk.
```bash
gcloud storage cp ~/pipeline_outputs/rnaseq/embeddings/*cohort_full_vcdr3* gs://aleix-rnaseq-wb-cordial-leechee-9743/embeddings/
```

Related: [[../../context/STATUS.md]] (root, Marc's).
