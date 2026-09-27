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
- **Ancestry recovery gap: reframed by the full cohort (2026-09-27), not yet tested.**
  At n=100: AFR high / EAS low, 1.63x, "consistent across all 5 major chain types," KW p=0.18.
  At n=7,922 (group means only, `results/rnaseq_cohort_ancestry_summary.csv`) **both halves of
  that framing fail**: (1) EUR, not EAS, is now lowest (AFR/EUR 1.44x; EAS mid-pack); (2) the
  gap is **carried by the B-cell chains** (IGH/IGK/IGL AFR/EUR 1.68-1.75x) while **TRB is
  near-flat** (AFR/EUR 1.06x, all-group spread 1.19x, AFR not even highest). So "reference/
  database bias against non-European receptor sequences" is now the *less* likely reading —
  that would hit TRB too. Leading hypothesis to test, not assert: a biological Ig-expression
  difference (plasmablast/plasma-cell Ig mRNA load; higher serum Ig in African-ancestry
  populations is documented) rather than a TRUST4 artifact. **Not yet a result:** means only,
  no dispersion, no per-person test, no depth normalization at this n. Needs the per-person
  per-chain test (VM-local `batch_summary_detail.tsv` exists) before any claim.
- **Group min/max in committed summaries are single-person values.** `cdr3_min`/`cdr3_max` in
  `results/rnaseq_cohort_ancestry_summary.csv` are each one individual's count. Not an
  identifier, but an n=1 statistic in a public repo — same small-cell disclosure class as the
  root repo's open item. Flag for Marc/supervisors; don't decide here. Leaning: drop min/max
  (or report p5/p95) in future committed summaries.
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

- **TRB is the primary chain for every HLA-facing analysis (Cole Shanks, 2026-09-27).**
  Cole: what matters most are chains that undergo V(D)J recombination, and he's "most
  interested specifically in the TCRBs." TRG/TRD don't really interact with HLA (γδ T cells
  are largely not MHC-restricted) — ignore for now; possibly an interesting disease signal
  later, not an immediate focus. Consequences: `--chain TRB` (already `embed_cdr3s.py`'s
  default) is the headline for embeddings and HLA joins; every ancestry/recovery result is
  reported per chain with TRB first, never only as an all-chain total (the all-chain total
  is dominated by Ig chains and gave the wrong story at n=7,922 — see open question above).
  TRA stays a secondary (it pairs with TRB, but bulk RNA-seq can't pair them). IG chains are
  kept for the recovery-confound analysis, not the HLA aim.
- **All figures use `cnsplots`, journal (Cell/Nature/Science) style (David Bonet,
  2026-09-27).** Current figures read as AI-generated; David asked that every figure from now
  on follow `github.com/faridrashidi/cnsplots` (matplotlib-based, BSD-3, `pip install
  cnsplots`, bundled Claude skill via `cnsplots skill install --agent claude --scope
  project`). Conventions to follow: sizes in points (`cns.figure(width, height)`),
  `cns.multipanel(max_width=...)` with lettered panels, built-in stats annotation
  (`test=`, `p_adjust=`) rather than hand-drawn asterisks, export SVG/PDF with editable text
  (`cns.savefig`), PNG 300 dpi only as a preview. **Where figures get rendered:** figures
  drawn purely from committed group-level CSVs can be built locally; any figure that needs
  per-person values (box/violin/ECDF/scatter) is rendered **on the VM** with cnsplots
  installed there, and only the finished figure leaves — no individual points drawn for
  groups/bins under n=20, no per-person data file ever egresses. So cnsplots goes into
  `aleix/RNA-seq/pixi.toml` (or the VM's base python) as well as locally.

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
