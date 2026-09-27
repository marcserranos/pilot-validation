# Status — live session state (RNA-seq / repertoire workstream)

> **Role:** where we are *right now* + the literal next commands. The only file rewritten
> fully each session. Anything durable graduates to ENVIRONMENT/DECISIONS/EXPERIMENTS.

## As of 2026-09-27 — SCEPTR embedding of the full cohort running; viz script ready

**Scope set by Aleix (2026-09-27):** the current task is only this: ready TRB repertoires,
embed them with SCEPTR, visualize the embeddings. HLA prediction, disease association and
ancestry-confound experiments are deliberately out of scope for now. The only ancestry check
kept is the free one inside the embedding run (fraction of clonotypes excluded for lacking
a usable V gene, per ancestry).

**Done this session:**
- Supervisor decisions recorded in DECISIONS.md: TRB-first (Cole), cnsplots figures (David).
- `embed_cdr3s.py` fixed (`8d4a719`): clonotype dedupe + rank by reads, chain from V/J/C,
  canonical-junction filter, report.tsv fallback off by default, `b_sceptr` on TRBV+CDR3B,
  exact O(n·d) V-gene contrast (the old n×n needed ~60 TB), tagged outputs.
- Found: ~5,566 of 7,922 people had no `cdr3.out` on the main VM (only `report.tsv` was
  flattened back after the sharded batch). Pulled from the bucket before the embed run.
- `scripts/01_sceptr_embedding_viz.py` written and tested end-to-end on synthetic data
  locally (figures inspected: legends, palettes, layering, clipping).

**Running on the main VM (`00eb81c5cc77`):** `embed_cdr3s.py cohort_full.tsv --models sceptr`
under nohup, log `~/embed_full.log`. Confirm the log's `Input sources:` line says 7,922
`cdr3.out` and 0 `report.tsv`.

## Pick up here

1. Check the embedding finished:
   ```bash
   grep -E "Input sources|skipped|excluded|accepts|embedded|Summary|FATAL|failed" ~/embed_full.log | tail -20
   ```
2. Pull and run the visualization (UMAP on 200k clonotypes is single-threaded for
   reproducibility, ~10-20 min):
   ```bash
   cd ~/repos/pilot-validation && git pull
   cd aleix/RNA-seq && nohup pixi run python3 -u scripts/01_sceptr_embedding_viz.py > ~/viz01.log 2>&1 &
   ```
3. Figures + `summary.csv` land in `~/pipeline_outputs/rnaseq/reports/01_sceptr_embedding_viz/`
   (VM-local). View in JupyterLab's file browser. Figure dots are individual clonotypes /
   people — **don't commit figures to the public repo until the disclosure question in
   DECISIONS.md is settled**; `summary.csv` (aggregate) is fine.
4. Write `reports/01_sceptr_embedding_viz/README.md` from the real figures + summary, append
   EXPERIMENTS.md.

**Expect in Fig 2F:** people with few clonotypes sit at the map's edges — a mean over fewer
vectors is noisier (scales ~1/√n). Arithmetic, not biology; call it out in the report.

Related: [[../../context/STATUS.md]] (root, Marc's).
