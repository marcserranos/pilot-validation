# Environment — RNA-seq Verily Workbench VM(s)

> **Role:** ops reference, on-demand-consult tier but read in full at session start (small).
> **Edit:** append a quirk the moment you hit one; never rewrite an existing entry.
> Marc's `../../context/ENVIRONMENT.md` quirks are for his shared Workbench VM and mostly do
> NOT apply here — separate project, separate VM(s), no Chrome computer-use (driven by
> terminal directly).

## Layout

- Workspace GCP project: `wb-cordial-leechee-9743`. CDR (BigQuery): `wb-silky-artichoke-2408.C2025Q4R6`.
- Workspace has both long-read calls and RNA-seq data, so nothing controlled moves between
  workspaces.
- Repo: `~/repos/pilot-validation` on each VM, branch `aleix/hla-resolve-phase1` — **NOT
  `main`** (main is 134 commits behind and has no `aleix/` directory at all).
- Data under `~/mnt/aou-controlled/v9/` once gcsfuse is mounted: RNA-seq BAMs =
  `v9/multiomics/rnaseq/manifest.tsv` (8,980 people, STAR `.md.bam` ~5-20 GB each, physically
  at `pooled/multiomics/v9_base/rnaseq/bam/`); long-read = `v9/wgs/long_read/`.
- Two CONTROLLED buckets, both SHARED_ACCESS (Marc can see them):
  - `gs://aleix-disease-counts-wb-cordial-leechee-9743` — phenotype pulls + burden summaries
    (Cole's disease-counts deliverable). **Not** for repertoire data.
  - `gs://aleix-rnaseq-wb-cordial-leechee-9743` — LR calls + embeddings + repertoire results.
    `repertoire_results/<research_id>/` and `repertoire_shards/` live here.

## Quirks & fixes

1. **`git clone` with no branch arg defaults to `main`.** `main` doesn't have `aleix/RNA-seq/`
   at all (different top-level structure: `context/`, `reference/`, `scripts/`, `reports/`,
   `sprints/`, `pixi.toml` at repo root). Every fresh clone: `git clone <url> && cd
   pilot-validation && git checkout aleix/hla-resolve-phase1 && git pull`. Hit repeatedly
   across 3 of 4 fresh VMs during the 2026-09-19 sharded batch launch.

2. **`pixi` is not preinstalled on a fresh VM.** Install: `curl -fsSL
   https://pixi.sh/install.sh | bash` then `export PATH="$HOME/.pixi/bin:$PATH"` for the
   current shell (the installer updates `~/.bashrc` for future shells, but not the current
   one). `run_rnaseq_batch_local.sh` checks for it on PATH and fails loud (`FATAL:`) if
   missing, rather than failing confusingly mid-batch.

3. **gcsfuse mount does not survive a VM stop/restart.** Remount before trusting anything
   that reads AoU data:
   `gcsfuse --billing-project "$GOOGLE_PROJECT" --implicit-dirs vwb-aou-datasets-controlled ~/mnt/aou-controlled`
   Only needed for `build_rnaseq_cohort.py` (reads manifests over the mount); **not** needed
   for `run_rnaseq_batch_local.sh` itself, which reads BAMs via direct `gcloud storage cp`
   (see quirk 5 for why that split exists). Verify with `ls` on a known manifest path before
   trusting it.

4. **VMs auto-stop when idle.** This happens legitimately after a batch finishes, not only as
   a crash/preemption signal — check `run.log` for `==== batch finished ::` before assuming
   data was lost. (Misdiagnosed once, 2026-09-20 — see root DECISIONS-equivalent lesson.)

5. **gcsfuse streams a BAM at ~11 MB/s; native `gcloud storage cp` pulls the same BAM at
   ~670 MB/s** (measured 2026-09-10, 4.8 GB in 15s, even while gcsfuse was saturated — it
   bypasses the fuse daemon). This is *why* `run_rnaseq_batch_local.sh` copies each BAM to
   local disk before running TRUST4 rather than reading through the mount directly — TRUST4's
   two full sequential passes over a 5-20 GB BAM don't parallelize over one shared gcsfuse
   pipe (Experiment B: 9 people/hr at `--jobs 4`, ~$925/38 days for 8,327 people). Local-copy
   + delete makes extraction CPU-bound; `--jobs 8-12` scales. Disk-bounded: only `--jobs`
   BAMs on disk at once, worst case ~20 GB/BAM.

6. **Sharding a cohort across N VMs is collision-safe by construction, not by luck.** Shard
   files are built via `split -l` on a single `remaining.tsv` (itself = cohort minus anyone
   with an existing local `report.tsv`), so a `research_id` never appears in more than one
   shard, and `research_id` is globally unique. No coordination needed between VMs — each
   syncs to the same results bucket, per-person, safely.

7. **The runner skips already-done people without syncing them.** If per-person bucket sync
   is added to the script *after* some people already finished locally (exactly what happened
   2026-09-19 — the first 500 people predated the sync feature), those people exist nowhere
   but that VM's disk until a manual backfill loop runs. Re-running the batch script won't
   fix this — `run_one()` returns early on an existing `report.tsv` before reaching the sync
   step. Backfill is idempotent (safe to re-run): loop over local dirs with a report.tsv,
   `gcloud storage cp -r` each into `$RESULTS_BUCKET/`.

8. **`gcloud storage cp SRC1 SRC2 DST1 DST2` is not "copy each source to its paired dest."**
   With 2+ sources the last arg must be a single existing destination directory. Copy the BAM
   and its `.bai` as two separate single-file `cp` calls, not one 4-arg call.

## Machine spec (2026-09-19 full-batch run)

n1-highmem-16, `--jobs 8`, 4 VMs sharding ~1,856-1,857 people each. Chosen over a cheaper
spec because TRUST4 batches are CPU-bound once BAMs are staged locally (quirk 5) — a smaller
VM would directly cost wall-clock time against the 3-4 day deadline, not just be slower to
provision. See `context/STATUS.md` / `../../aleix/RNA-seq/README.md` machine history for the
smaller single-VM spec used for the earlier 100/500-person runs.

Related: [[../../context/ENVIRONMENT.md]] (Marc's, different VM/project),
[[../ENVIRONMENT_LOCAL.md]] (Aleix's WSL2 laptop — HLA-Resolve workstream, not this one).
