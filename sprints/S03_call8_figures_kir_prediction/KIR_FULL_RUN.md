# KIR full-cohort run — operations record

## Original app config (restore after the run)
Recorded 2026-09-24 from Workbench → Full Cohort HLA Calling → Apps → big_run → Edit (not saved):
- App: `AoU_Jupyter_ComputeEngine_20260805_big_run`, instance `aoujupytercomputeengine20260805bigrun`, zone us-central1-a
- Compute: General, **machine type `n2-highmem-4` (4 vCPU, 32 GB)**
- Data disk: **2000 GB**
- Autostop: **on, 1 hour idle**
- Cost estimate shown: $0.37/h running; disk $81.60/month

## Target for the run
Same as the HLA production run: `n2-highcpu-96` (96 vCPU/96 GB), orchestrator `--concurrency 24 --threads-per-person 4`,
separate `--outroot ~/pipeline_outputs_kir` (never the HLA outroot).

## Log
- 09-24 23:45 Pre-run checksums of the 14 HLA tables: `~/s03/kir_full/hla_tsv_md5_pre.txt` (+ `~/s03/kir_full/PRE_RUN_STAMP`).
- 09-24 23:50 Orchestrator deployed as untracked `~/repos/pilot-validation/scripts/production_orchestrator/run_production_orchestrator_s03kir.py` (md5 fe06b79c…, = local HEAD).
- 09-24 23:55 Dry run: 2 people (random, non-self-align), `--outroot ~/pipeline_outputs_kir_test`, conc 2 × 2 threads, log `~/s03/kir_full/dryrun.log`.
- 09-25 00:05 Dry run PASSED (2/2, HLA md5 all OK, 0 newer files in HLA outroot). Next: stop app → Edit machine type n2-highcpu-96 → start.
- 09-25 00:15 Resized to n2-highcpu-96 (Update saved; $3.55/h live). 00:16 Start failed: GCP capacity stockout in us-central1-a; retrying.
- 09-25 01:40 Marc reprovisioned to n2-highcpu-80 ($2.98/h). Run: conc 20 × 4 threads.
- 09-25 02:15 Benchmark: ~390 people/h, 0.7% fail, 6.5 MB/person. Projection ~31 h / ~$93. Waiter armed to launch full run when bench exits.
- 09-25 02:30 Bench done 198/200 (0.58 h). FULL RUN started 06:57:30 UTC, log ~/s03/kir_full/full_run.log.
