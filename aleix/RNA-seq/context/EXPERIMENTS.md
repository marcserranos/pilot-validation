# Experiments Log — RNA-seq / repertoire workstream

> **Role:** append-only record of every pipeline run — what was run, the result, the runtime.
> **Edit:** append a dated entry per run; never rewrite past entries.
> **Read:** for what's been tried and how long it took, top to bottom.

## 2026-08-10 — feasibility proven, single person

TRUST4 recovers a real repertoire from an AoU whole-blood RNA-seq BAM: **2,013 CDR3s from
person 1000291**, all 7 chain types represented, ~9 min end-to-end. Required
`--abnormalUnmapFlag` (AoU's STAR run uses `--outSAMunmapped Within`, which TRUST4 doesn't
handle by default) — fixed in `scripts/run_trust4_sample.sh`. Writeup:
`results/1000291_trust4_smoke_test.md`.

## 2026-08-11 — first 100-person cohort

Ancestry-stratified (17/17/17/17/16/16). `--jobs 1` vs `--jobs 3` timed on 3 people: 2.55x
speedup, near-ideal, no CPU contention. `--jobs` efficiency sweep (4/6/8/12/16) overnight:
non-monotonic (population noise between test slices dominated the signal), no evidence of a
hard ceiling up to 16 — all configs beat sequential by 1.9-4.1x. Full 100-person cohort
completed at `--jobs 12`.

**Ancestry recovery gap found:** AFR highest (7,406 mean CDR3s), EAS lowest (4,538) — 1.63x
spread, consistent across all 5 major chain types. `results/rnaseq_100person_ancestry_writeup.md`.

**Depth confound checked:** explains ~40% of the gap (1.63x -> 1.38x normalized), not all of
it. AFR/EAS stay the extremes either way. `results/rnaseq_depth_confound_check.csv`.

**Remaining confounds checked:** RQS (r=0.061) and cell-composition proxy (flat across
ancestry) both ruled out. Kruskal-Wallis H=7.62, p=0.18 — not significant at n=100. The
ancestry pattern is real in the data but doesn't clear significance at this n.

## 2026-08-17 — QC confounds, foundation docs, overlap measurement

**16 AoU-computed QC metrics checked against recovery** — nothing explains the ancestry gap.
`Mean 3' bias` trips the two-part test but direction is inconsistent (AMR: highest bias,
second-highest recovery) and 4% of variance can't produce a 1.63x gap. `Genes Detected` is
the best recovery predictor found (r=0.608, 37% of variance) but flat across ancestry — a
covariate, not a confound. `Base Mismatch` varies by ancestry (p=0.028, reference-bias
signature) but doesn't predict recovery. `results/rnaseq_qc_confounds_check.csv`.

Read length verified constant at 146 bp across all 8,980 samples — TRUST4's k-mer threshold
(29 bp) is uniform, per-sample extraction stringency ruled out as a confound.

**N1 resolved:** 8,980 RNA-seq samples = 8,980 unique people, no technical replicates.

**LR x RNA-seq overlap measured: 8,327 people**, 34.2x above independent-sampling
expectation, ancestrally balanced (EUR 29.7% ... MID 5.1%), 5 of 6 groups >1,000 — this
becomes the base cohort for every joint analysis from here on.

**Disease study complete:** 12,172 conditions, 2,463 reportable, 7 of 10 immune-mediated
families viable for case/control. `results/lr_rnaseq_disease_study_writeup.md`.

## 2026-09-10 — copy-local staging designed and validated

Measured: gcsfuse streams a BAM at ~11 MB/s; native `gcloud storage cp` pulls the same BAM at
~670 MB/s (4.8 GB in 15s, even while gcsfuse was saturated). Redesigned the batch runner
around local-copy-then-delete (`run_rnaseq_batch_local.sh`) instead of reading through the
mount directly — see `context/ENVIRONMENT.md` quirk 5 for the full mechanism. Projected
~$30-70 on-demand / ~$10-20 spot for the full 8,327-person cohort, 1-3 days sequential (hours
if sharded).

## 2026-09-10 — 500-person repertoire recovery + embedding comparison

Real 500-person TRB repertoire pool built. SCEPTR vs ESMC-300M compared on the same CDR3
pool: SCEPTR separates same-/diff-V-gene pairs by gap 0.019 in 23s (dim 64); ESMC-300M's gap
is 0.0006 (30x worse discrimination, embedding space nearly collapsed, cos-sim ~0.98
everywhere) and took 3,666s (156x slower). **SCEPTR adopted as the embedding model going
forward** (see DECISIONS.md). `results/cdr3_embedding_comparison.csv`,
`results/rnaseq_cohort_ancestry_summary.csv` (500-person version, later superseded).

## 2026-09-19 to 2026-09-20 — full 7,922-person cohort, sharded 4-VM batch

Cohort built: `build_rnaseq_cohort.py --per-group 2500 --max-bam-gb 12` → 7,922 people
(8,326 LR-overlap-eligible minus 404 with BAM > 12 GB: 81 AFR / 92 AMR / 38 EAS / 124 EUR /
23 MID / 46 SAS, exact match). Remaining 500 already-done people excluded via `remaining.tsv`
(7,422 people), split 4 ways (`shard_0..3.tsv`, ~1,856 each) across 4 n1-highmem-16 VMs,
`--jobs 8` each.

**Launch cascade (all fixed live, see context/ENVIRONMENT.md for the durable lessons):**
wrong default results bucket (disease-counts bucket instead of rnaseq bucket) → fixed, commit
`b9f7e65`; one VM's clone predated that fix (`unknown arg: --results-bucket`) → `git pull` +
relaunch; 2 fresh VMs cloned `main` instead of `aleix/hla-resolve-phase1` → branch checkout;
2 fresh VMs missing `pixi` → installed; gcsfuse mount lost after a VM restart mid-session →
remounted. Original 500 people predated the per-person sync feature and existed nowhere but
the original VM's disk until a one-time backfill loop pushed them to the bucket.

**Result: zero failures.** All 4 shards finished cleanly (confirmed via `run.log` "batch
finished" markers, not inferred from VM auto-shutdown). Bucket-vs-cohort reconciliation
(`comm -23`/`comm -13` on sorted research_id lists): **0 people in cohort but missing from
bucket**; 36 harmless extras (old `<rid>_expB_j4/j8` Experiment-B timing dirs + 4 stray pilot
people not in the current cohort: 1000955, 1001121, 1002660, 1004700). Full local flatten
confirmed: **7,922 / 7,922** report.tsv files present on the main VM
(`00eb81c5cc77`, `~/pipeline_outputs/rnaseq/`).

**Full-cohort aggregation** (`aggregate_rnaseq_results.py cohort_full.tsv`, commit
`05e727f`, run from the VM): 7,922/7,922 people with a report. Mean CDR3s per person: AFR
6,331 · SAS 5,998 · EAS 5,406 · AMR 5,230 · MID 5,209 · EUR 4,392. Per chain, AFR/EUR:
TRB 1.06x · TRA 1.12x · IGH 1.68x · IGK 1.75x · IGL 1.75x. **Reverses two parts of the n=100
story** — EUR (not EAS) is lowest, and the gap is Ig-chain-driven with TRB near-flat (see
DECISIONS.md open question). Group means only; per-person test not yet run.

## 2026-09-27/28 — full-cohort SCEPTR TRB embedding + first visualization (script 01)

**Embedding** (`embed_cdr3s.py cohort_full.tsv --models sceptr`, commit `8d4a719`, main VM,
n1-highmem-16, CPU): input 7,922/7,922 people from `cdr3.out` (0 unfiltered fallback — after
pulling ~5,566 missing `cdr3.out` from the bucket). Pool: **3,882,613 TRB clonotypes**
(unique TRBV+CDR3aa per person, canonical junction, score ≥ 0.02, top 500 by reads),
59 distinct TRBV genes. `b_sceptr` accepts 48 of the TRBV symbols seen; **45,567 clonotypes
(1.2%) excluded** for no/unusable V — **flat across ancestry (AFR 1.1%, all others 1.2%)**, so
the V-usability filter introduces no ancestry bias. 3,837,046 embedded (64-dim) in ~57 min
(~1,100 clonotypes/s single process). Summary:
`results/cdr3_embedding_summary_cohort_full_vcdr3.csv` (V-gene gap is circular here — V is
an input — recorded, not interpreted).

**Visualization** (`01_sceptr_embedding_viz.py`, commit `39e9da3`): Fig 1 clonotype UMAP (200k
seeded subsample), Fig 2 person UMAP (7,922, mean vector), S1 PCA versions. Figures VM-local
(`~/pipeline_outputs/rnaseq/reports/01_sceptr_embedding_viz/`), not committed — disclosure
question. Aggregate numbers:
- Median 494 clonotypes/person: nearly everyone hit the 500 cap, so persons are compared on
  a near-fixed-size sample (depth largely decoupled from the person vector for most people).
- **14.3% of clonotypes are public** (identical TRBV+CDR3aa in ≥2 people). Higher than the
  few-% often quoted — plausibly because we keep each person's most expanded clones and
  match at amino-acid level across 7,922 people. Unverified; check vs CDR3 length.
- Clonotype PCA: PC1 4.5%, PC2 3.8% — variance spread over many dimensions (expected).
- Person PCA: PC1 18.5%, PC2 12.8%. Person PC1 vs repertoire size: Spearman ρ = −0.21
  (p ≈ 1e-77); PC2 ρ = +0.12. A real but modest size effect on the main person axis — the
  small-repertoire minority, consistent with noisier means over fewer vectors.

**Post-batch cleanup (2026-09-2x):** confirmed every shard VM had zero `.sync.err` files and
an empty `_staging/` before deleting its disk — nothing local was the only copy of anything.
3 shard VMs deleted; main VM kept running (has `cohort_full.tsv` + all flattened reports,
needed for aggregation/embeddings).
