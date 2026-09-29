# Status — live session state (RNA-seq / repertoire workstream)

> **Role:** where we are *right now*, plus the literal next commands. The only file fully
> rewritten each session. Anything durable graduates to ENVIRONMENT, DECISIONS or EXPERIMENTS.

## As of 2026-09-28 — reports 01–03 have real results; 04 rewritten around a null; 05 new

**01 embedding maps** — done, committed with README.

**02 atlas** — done on real data. **Immunosenescence replicated at n=7,105**: rarefied
diversity −3.3%/decade (95% CI −3.5 to −3.1, p≈2e-168), top-10 clonal share +1.91 pp/decade
(p≈4e-186). Depth-adjusted effect is unchanged (−3.2%/decade). Median 1,136 TRB clonotypes
per person; 7.2% of clonotypes shared by ≥2 people.

**03 publicness/robustness** — done on real data. **Publicness mechanism replicated**:
Spearman 0.66 between sharing level and log10 Pgen, spanning −9.3 (private) to −6.6 (≥100
people). **Negative result that changes the plan**: the mean-pooled person vector is not
stable — top-250 vs top-500 gives pairwise-distance ρ = 0.35 and 10-NN overlap 0.08;
read-weighting gives ρ = 0.18. Ancestry is predictable (macro AUROC 0.717) but plain TRBV
usage does as well (0.722), so the person vector's ancestry signal is germline V-gene
composition, not CDR3 chemistry.

**04 antigen specificity** — first real run exposed that naive VDJdb matching measures
chance: 88% of people "matched" CMV, but 69% matched HBV (US prevalence <1%), 27% HIV-1 and
71% human self-peptides, and carriage *fell* with age tracking repertoire size. Rewritten
(`68aaf18`): exact matching only, against a null of OLGA-generated decoy TCRs matched to
each VDJdb TCR on TRBV gene, CDR3 length and Pgen; reported as enrichment over that null,
with published seroprevalence and negative-control pathogens as references, and age tested
on excess (observed − expected). **Not yet re-run on the VM.**

**05 embedding metaclusters** — new, tested on synthetic data only. Partitions the 3.84M
clonotypes, describes what separates clusters, and tests whether known specificities
concentrate in particular clusters (Fisher, BH). Writes `clusters_3d.csv` for an interactive
3-d cluster map.

## Pick up here

1. Run 02→05 on the main VM (~30 min; 02/03 reuse caches):
   ```bash
   cd ~/repos/pilot-validation && git pull
   nohup bash aleix/RNA-seq/scripts/run_reports_02_05.sh > ~/reports.log 2>&1 &
   tail -5 ~/reports.log
   ```
2. When the log ends with `all reports done`, push the de-identified outputs:
   ```bash
   cd ~/repos/pilot-validation && for r in 02_repertoire_atlas 03_publicness_and_robustness 04_antigen_specificity 05_embedding_clusters; do mkdir -p aleix/RNA-seq/reports/$r && cp ~/pipeline_outputs/rnaseq/reports/$r/*.{png,pdf,svg,csv} aleix/RNA-seq/reports/$r/; done && git add aleix/RNA-seq/reports && git commit -m "Reports 02-05: VM outputs" && git push
   ```
3. Locally: per-report READMEs, the interactive 3-d cluster map from `clusters_3d.csv`, and
   `reports/COMPREHENSIVE_REPORT.md` (the narrative across 01–05). Append EXPERIMENTS.md.

**Owed regardless:** back up the embeddings, which exist only on the main VM's disk.
```bash
gcloud storage cp ~/pipeline_outputs/rnaseq/embeddings/*cohort_full_vcdr3* gs://aleix-rnaseq-wb-cordial-leechee-9743/embeddings/
```

## Decisions this opens (for Aleix / supervisors)

- **Person representation must change before the HLA phase** (03d). Mean pooling is unstable
  and mostly encodes V-gene usage. Candidates: V-usage profile as an explicit baseline to
  beat; per-cluster abundance profiles from 05; or supervised pooling trained against HLA.
- **Depth is the binding constraint** for anything clonotype-specific (median 2,076 TRB
  reads/person). It caps what 04 can ever detect; worth stating as a limit rather than
  fighting.
- **Small-cell disclosure** for figures with individual-level dots — still unresolved with
  Marc; all current figures avoid it by aggregating.

Related: [[../../context/STATUS.md]] (root, Marc's).
