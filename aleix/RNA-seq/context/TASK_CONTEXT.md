# Task context — RNA-seq / repertoire workstream

> **Role:** the "why are we doing this at all" register — stable, rarely edited.
> **Read:** first, every session, before touching STATUS/DECISIONS/EXPERIMENTS.

## The direction, stated once

AoU RNA-seq BAMs contain TCR/BCR (T-cell / B-cell receptor) reads that STAR discards during
alignment, because V(D)J-recombined receptor sequences don't exist in any reference genome.
The immune repertoire is completely absent from every expression file AoU ships. **TRUST4**
recovers it directly from the BAM via de novo assembly.

HLA type is predictable from TCR repertoire alone (AUC 0.95, beta chain, published result),
so Marc's HLA calls (AoU-native + long-read, from the HLA-Resolve workstream) become the
**labels** for this direction. Full reasoning and the three supervisor papers are in
`slides/RNAseq_directions*.pptx`.

## Pipeline shape (end to end)

1. **Cohort selection** — `scripts/build_rnaseq_cohort.py`, ancestry-stratified, restricted
   to the LR×RNA-seq overlap (8,327 people) and a BAM-size cap (≤12 GB → 7,922 people used
   for the full batch).
2. **Repertoire calling** — `scripts/run_rnaseq_batch_local.sh`, copy-local TRUST4 (see
   `context/ENVIRONMENT.md` for why copy-local, not gcsfuse-direct).
3. **Aggregation** — `scripts/aggregate_rnaseq_results.py`, per-person detail (VM-local,
   has research_ids) + de-identified ancestry-group `.csv` (safe to commit).
4. **Embeddings** — `scripts/embed_cdr3s.py`, SCEPTR (finalized choice — see DECISIONS.md).
5. **Downstream joins** — disease/phenotype (Cole's disease-counts work) and HLA labels
   (Marc's workstream) against the same 7,922/8,327-person cohort.

## Join point with Marc's workstream

`NEXT_STEPS_AND_RESEARCH_MAP.md` §3 — the HLA×TCR/BCR integration plan (J1–J5) connecting
Marc's HLA calls to this repertoire data across the shared 8,327-person LR×RNA cohort. Read
before starting joint work; check `context/DECISIONS.md` on both sides for concurrent-editing
notes before reading heavily from or merging his branch.

## Compute environment

Runs on **Verily Workbench (VWB)**, not classic All of Us/Terra — a different project from
Marc's shared Workbench VM. Full detail in `context/ENVIRONMENT.md`.

Related: [[../../context/DECISIONS.md]] (root, Marc's), `reference/` (the four foundational
research docs: AOU_RNASEQ_DATA_REPORT, TRUST4_DEEP_DIVE, POST_TRUST4_OPTIONS,
LR_RNASEQ_DISEASE_STUDY).
