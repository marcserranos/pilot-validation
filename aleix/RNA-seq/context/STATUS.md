# Status — live session state (RNA-seq / repertoire workstream)

> **Role:** where we are *right now* + the literal next commands. The only file rewritten
> fully each session. Anything durable graduates to ENVIRONMENT/DECISIONS/EXPERIMENTS.

## As of 2026-09-27 — full-cohort repertoire calling done; moving into aggregation + embeddings

**The 7,922-person TRUST4 batch is done, zero failures, fully reconciled** (EXPERIMENTS.md
2026-09-19/20 entry). All 3 shard VMs verified fully-synced (no `.sync.err`, empty
`_staging/`) and deleted. **Main VM (`00eb81c5cc77`) is the only one still running** — it has
`cohort_full.tsv` and all 7,922 flattened `<research_id>_report.tsv` files locally under
`~/pipeline_outputs/rnaseq/`, needed for the next two steps below. Don't delete it before
those are done and their outputs are pulled back.

`context/` for this workstream was just created (2026-09-27), mirroring the root repo's
doc-tier convention (TASK_CONTEXT/ENVIRONMENT/STATUS/DECISIONS/EXPERIMENTS), seeded from
existing README status + prior session history rather than written fresh. `reports/` created
alongside the existing flat `results/` — new numbered deliverables go in
`reports/NN_description/` from here on; existing flat files in `results/` stay where they are.

## Pick up here

1. **Aggregate the full cohort:**
   ```bash
   python3 ~/repos/pilot-validation/aleix/RNA-seq/scripts/aggregate_rnaseq_results.py \
     ~/pipeline_outputs/rnaseq/cohort_full.tsv
   ```
   Produces the full-cohort per-ancestry CDR3 recovery table (de-identified, safe to commit).
   Compare against the 500-person `results/rnaseq_cohort_ancestry_summary.csv` — the open
   question this directly informs is whether the AFR/EAS recovery gap (DECISIONS.md) reaches
   significance at 79x the sample.

2. **Run SCEPTR embeddings on the full cohort** (`scripts/embed_cdr3s.py`, `--pheno-dir
   ~/pipeline_outputs/rnaseq` default already matches). ESMC not planned at this scale
   (DECISIONS.md — 16hr projected runtime for no shown benefit); may add a `--models sceptr`
   flag to `embed_cdr3s.py` first so it doesn't default to running both.

3. **Write both results into a numbered `reports/` pair** (first entry under the new
   convention): e.g. `scripts/01_full_cohort_aggregate_and_embed.md`-style report folder, or
   fold into the existing `aggregate_rnaseq_results.py`/`embed_cdr3s.py` scripts if no new
   script is actually needed — check the pipeline table convention before assigning a number
   Marc's side already uses (`scripts/hla_popgen/README.md`) doesn't apply here; this
   workstream doesn't have a numbered-script convention yet, so the first one sets it.

4. After both land, update `context/EXPERIMENTS.md` (append, don't rewrite) and this file.

Related: [[../../context/STATUS.md]] (root, Marc's — check before reading/merging his branch).
