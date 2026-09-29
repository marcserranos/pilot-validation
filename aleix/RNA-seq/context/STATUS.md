# Status — live session state (RNA-seq / repertoire workstream)

> **Role:** where we are *right now*, plus the literal next commands. This is the only file
> rewritten in full each session. Anything durable graduates to ENVIRONMENT, DECISIONS or
> EXPERIMENTS.

## As of 2026-09-28: report 01 done (full-cohort SCEPTR TRB embedding + figures)

**Scope set by Aleix (2026-09-27):** get the TRB repertoires ready, embed them with SCEPTR,
and visualize them. Joining to HLA, disease prediction and ancestry-confound experiments are
out of scope for now.

**Done:**
- 7,922 people and 3,836,906 TRB clonotypes, embedded with `b_sceptr` (TRBV + CDR3B).
- Filters: `cdr3.out` for every person, score ≥ 0.02, canonical junction, ≤ 30 aa, top 500
  by reads.
- Report `reports/01_sceptr_embedding_viz/`: README, the binned (≥ 20 people per cell)
  cnsplots figures, and `summary.csv`. Figures are committed; per the rule above they're
  safe to share.

**Where things live:**
- Embeddings and aligned pools (real research_ids, VM-local only):
  `~/pipeline_outputs/rnaseq/embeddings/{embeddings,pool}_sceptr_cohort_full_vcdr3.*` on the
  main VM (`00eb81c5cc77`).
- The UMAP cache is in `~/pipeline_outputs/rnaseq/reports/01_sceptr_embedding_viz/_cache/`.

**Not yet backed up:** the embeddings exist only on the main VM's disk. Before deleting or
resizing that VM, copy them to the controlled bucket:
```bash
gcloud storage cp ~/pipeline_outputs/rnaseq/embeddings/*cohort_full_vcdr3* gs://aleix-rnaseq-wb-cordial-leechee-9743/embeddings/
```

## Pick up here

**Open follow-ups from report 01**, none scheduled; pick with Aleix. The first two are cheap
and aggregate:
1. Quantify publicness against CDR3 length (Fig 1B vs 1D).
2. Decide whether between-V expansion differences (Fig 1C) are biology or TRUST4 recovery.
3. Robustness across the top-N cap (100 / 250 / 500 / all).
4. Next phase, when scoped: join to Marc's HLA calls. Ancestry and repertoire size are the
   confounders to adjust for (Fig 2).

Optional, deferred: VDJdb known-epitope overlay (panel F in the original plan).

Related: [[../../context/STATUS.md]] (root, Marc's).
