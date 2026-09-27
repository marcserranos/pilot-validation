# Decisions & Open Questions — RNA-seq / repertoire workstream

> **Role:** the "why" register. Open questions at the top; resolved decisions with rationale
> below. Move an item down when settled — append, don't delete history.
> **Read:** to understand why the pipeline is shaped this way, or before reopening a call.

## Open questions

- **Real research_ids committed to a public GitHub repo, bound to individual-level findings**
  — in at least 6 files (incl. `results/1000291_trust4_smoke_test.md`,
  `../../context/ENVIRONMENT.md`). Flagged 2026-08-17, escalated to Marc, still undecided.
  **Do NOT add git push credentials to the VM** until resolved — paste results back so a
  human reviews every byte leaving the controlled environment. Same class of issue as the
  root repo's "bare person_id in a public repo" open question — don't resolve independently,
  needs a joint call with Marc + supervisors.
- **Ancestry recovery gap (AFR high / EAS low) — real but unexplained at n=100.** Depth
  explains ~40%, remaining confounds (RQS, cell composition, 16 QC metrics) all ruled out.
  Kruskal-Wallis was not significant at n=100 (p=0.18) despite the pattern looking consistent
  across chain types. **Now testable at n=7,922** — 79x the sample. If it holds and reaches
  significance, this becomes a real methodological finding (reference-genome/database bias
  against non-European TCR/BCR sequences) worth its own writeup, not just a QC footnote.
- **Rarefaction vs. post-hoc depth normalization for the recovery metric.** Flagged in
  `reference/TRUST4_DEEP_DIVE.md` §8.5 as probably the most important methodological fix
  outstanding — post-hoc normalization only closed 40% of the ancestry gap, which itself
  argues for switching the metric rather than continuing to patch around it. Not yet done.
- **Whether disease burden predicts CDR3 recovery.** New lead (2026-08-17) — both variables
  now exist per-person for the same people. Cheap to check, not yet run.
- **HLA-label attachment timing** — attach Marc's HLA calls (AoU-native + long-read) once the
  ancestry-recovery confound question is far enough along that it won't contaminate the join.
  Not yet started; gated on the full-cohort ancestry-gap result above, not on anything
  external.

## Resolved decisions

- **SCEPTR adopted as the embedding model (2026-09-10, EXPERIMENTS.md).** Compared against
  ESMC-300M on the same 500-person CDR3 pool: SCEPTR gap (same- vs diff-V-gene cosine
  similarity) = 0.019 in 23s; ESMC-300M gap = 0.0006 (30x worse, embedding space nearly
  collapsed) in 3,666s (156x slower). Not a close call — SCEPTR wins on both signal and cost.
  ESMC is not planned to be re-run at full-cohort scale (would cost ~16 hours vs. SCEPTR's
  ~6 minutes for no discrimination benefit shown so far).
- **Copy-local BAM staging (not gcsfuse-direct reads) for the TRUST4 batch runner
  (2026-09-10).** gcsfuse streams a BAM at ~11 MB/s and doesn't parallelize (one shared fuse
  pipe); native `gcloud storage cp` pulls the same BAM at ~670 MB/s and bypasses the fuse
  daemon entirely. Local-copy-then-delete makes extraction CPU-bound, so `--jobs` actually
  scales. Full rationale in `scripts/run_rnaseq_batch_local.sh` header and
  `context/ENVIRONMENT.md` quirk 5.
- **Sharded 4-VM batch (n1-highmem-16, `--jobs 8` each) for the full 7,922-person cohort
  (2026-09-19), not a cheaper single-VM run.** Explicit hard constraint from the 3-4 day
  deadline: only accept a cost reduction with zero performance/timeline risk. Sharding turns
  a multi-day sequential job into a same-order-of-magnitude parallel one; collision-safety is
  structural (globally-unique research_id + disjoint `split -l` shards), not best-effort, so
  the parallelism carries no correctness risk.
- **Per-person incremental sync to GCS (not end-of-batch only), added 2026-09-19.** Protects
  against VM preemption/shutdown/loss during a multi-day unattended run — a killed VM loses
  at most the one person it was mid-run on, not the whole shard. Retrofitted after the
  original 500-person run had already completed without it, which is why that batch needed a
  manual backfill (see ENVIRONMENT.md quirk 7) — the gap this decision closes going forward.
- **Results bucket is `gs://aleix-rnaseq-wb-cordial-leechee-9743`, never the disease-counts
  bucket (corrected 2026-09-19, commit `b9f7e65`).** The two buckets have different owners:
  disease-counts = Cole's phenotype/burden deliverable; rnaseq = LR calls + embeddings +
  repertoire results. Caught mid-batch-launch by a direct question ("why in disease counts,
  aren't we storing the rnaseq repertoires in the rnaseq bucket") — the default had been set
  wrong in the script; fixed before most of the batch ran.

Related: [[../../context/DECISIONS.md]] (root, Marc's — the shared public-repo/research_id
disclosure question is tracked jointly there and here).
