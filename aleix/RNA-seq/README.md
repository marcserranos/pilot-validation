# aleix/RNA-seq/ — immune repertoire workstream

Isolated workspace for the RNA-seq / TCR-BCR repertoire direction. HLA calling is handed off
to Marc; this folder is the new work. Kept separate from `../` (the HLA-Resolve workstream)
and from `../../context/` (Marc's Experiments A–D) — same isolation rule as before.

> **Data rule unchanged from the root repo: no AoU participant data here.** Public reference
> material (TRUST4's own small db files, tool docs) is fine — it is not participant data.

Confidence key, used throughout this folder same as `../reference/AOU_DATA_ACCESS_NOTES.md`:
**[HIGH]** confirmed live or from primary-source docs. **[MED]** inferred from credible
sources, not yet verified live. **[LOW]** best guess, needs live verification.

---

## The direction, stated once

AoU RNA-seq BAMs contain TCR/BCR (T-cell / B-cell receptor) reads that **STAR discards**
during alignment, because V(D)J-recombined receptor sequences don't exist in any reference
genome. That means the immune repertoire is completely absent from every expression file
AoU ships — nobody has mined it. **TRUST4** recovers it directly from the BAM via de novo
assembly.

HLA type is predictable from TCR repertoire alone (AUC 0.95, beta chain, published result) —
so our own HLA calls (AoU-native + long-read, from the HLA-Resolve workstream) become the
**labels** for this direction, rather than being discarded. Full reasoning, the three
supervisor papers, and the two-directions comparison (expression vs. repertoire) are in the
`RNAseq_directions*.pptx` decks in `slides/`.

**RNA-seq manifest — found 2026-08-10:**
`gs://vwb-aou-datasets-controlled/v9/multiomics/rnaseq/manifest.tsv`
(columns: `sampleid`, `research_id`, `markduplicates_bam_file_path`,
`markduplicates_bam_index_path`). Sibling dir `v9/multiomics/` also holds `proteomics/`.
The `rnaseq/` dir itself has `manifest.tsv`, `rnaseq_metadata.tsv`, `eqtl/`, `rnaseqc2/`,
`rsem/`, `sqtl/` — matches the deliverables table in the primary-source PDF. BAMs are
mark-duplicates-applied (`.md.bam`), physically under
`pooled/multiomics/v9_base/rnaseq/bam/`.

**Feasibility — confirmed 2026-08-10.** TRUST4 recovers a real repertoire from an AoU
whole-blood RNA-seq BAM: 2,013 CDR3s from person 1000291, all 7 chain types represented,
~9 min end-to-end. Full writeup: `results/1000291_trust4_smoke_test.md`. Required a
`--abnormalUnmapFlag` fix (now in `scripts/run_trust4_sample.sh`) because AoU's STAR run
uses `--outSAMunmapped Within`, which TRUST4 doesn't handle by default.

## Batch runs (N-person cohorts)

Three scripts, mirroring `../../../scripts/build_experiment_d_cohort.py`'s established
pattern exactly:

1. `python3 scripts/build_rnaseq_cohort.py --total 100` — ancestry-stratified pick,
   ground-truth-verified against the mount. Writes `~/pipeline_outputs/rnaseq/cohort.tsv`
   (VM-local, has real research_ids — see privacy note below).
2. `bash scripts/run_rnaseq_batch.sh <cohort.tsv> --jobs N` — resumable (skips anyone
   already done), includes the `--abnormalUnmapFlag` fix by default.
3. `pixi run python3 scripts/aggregate_rnaseq_results.py <cohort.tsv>` — produces a
   VM-local per-person detail file AND a de-identified, ancestry-group-level `.csv` in
   `results/` (the only one of the three that's safe to commit).

**Test small before trusting a time estimate for the full run.** `head -4 cohort.tsv >
cohort.test3.tsv`, then run that slice at `--jobs 1` vs `--jobs 3` and compare wall-clock
time. Extraction streams the whole BAM over the network-mounted bucket — if that's the
real bottleneck (not CPU), more parallel jobs may not scale the way more vCPUs would
suggest. Untested as of 2026-08-10; don't assume either direction.

**Naive sequential estimate for 100 people:** ~9 min/person (the one proven data point) ×
100 ≈ 15 hours, though real depth varies across the cohort (QC histogram shows most
people 95–145M read pairs, tail past 350M), so more honestly 12–25 hours sequential.
Could be much less with working parallelism — that's exactly what the small test above
answers before committing to an unattended multi-hour run.

### Privacy convention for this batch tooling

`cohort.tsv`, per-person TRUST4 output, and `batch_summary_detail.tsv` all contain real
`research_id`s and **stay VM-local under `~/pipeline_outputs/`, never committed** — same
posture `build_experiment_d_cohort.py` already established for the HLA workstream. Only
`aggregate_rnaseq_results.py`'s **ancestry-group-level** summary (no individual IDs) is
safe to commit, and it's the only one of the three scripts that writes into `results/`.

## How we work together on this (the git loop)

1. Scripts and docs are authored **here**, in this repo, on the local machine (with Claude).
2. Commit + push from here.
3. On the Workbench VM: `git pull` inside `~/repos/pilot-validation` to get the latest.
4. Run scripts from `~/repos/pilot-validation/aleix/RNA-seq/scripts/` on the VM.
5. Anything the VM produces that's worth keeping (numbers, short summaries, not raw
   BAMs/FASTQs) gets written back into `results/*.md` and pushed from the VM, or pasted back
   here — either way it ends up in git, not just on a VM that can be deleted.

No more pasting full command blocks back and forth — the VM always has the current script via
`git pull`.

## Reference documents (written 2026-08-17, week of "foundation, not batch runs")

Four research documents, primary-sourced, written without VM access. Read in this order:

1. **`reference/AOU_RNASEQ_DATA_REPORT.md`** — what the 8,980-sample RNA-seq dataset
   actually is: generation, quality, every deliverable, what it can and cannot support.
   Ends with 6 numbered lookups (N1–N6) that need a VM.
2. **`reference/TRUST4_DEEP_DIVE.md`** — how TRUST4 works, read off its **source code**:
   the 4 stages, the 3 candidate-capture routes, the **two-pass BAM read** that explains
   every performance observation we've made, every parameter, every output column, and the
   honest limitation list. Ends with 7 concrete design recommendations.
3. **`reference/POST_TRUST4_OPTIONS.md`** — shallow map of what comes after: summary
   statistics → clustering/database lookup → embeddings (SCEPTR) → supervised
   classification. Organized around the constraint that bulk RNA-seq gives **unpaired**
   chains.
4. **`reference/LR_RNASEQ_DISEASE_STUDY.md`** — design for the long-read × RNA-seq overlap
   disease study, with the 3-script pipeline and what's needed from Aleix (T1–T6).

## Layout

- **`scripts/`** — everything runnable. Read a script before running it; nothing here is
  proven yet (see confidence tags inside each file).
- **`reference/`** — the four research docs above, plus TRUST4's small db files once
  fetched (`hg38_bcrtcr.fa`, `human_IMGT+C.fa` — a few MB each, not participant data, safe
  to vendor once we have them).
- **`results/`** — `results/*` is gitignored except `.md`/`.csv` — same policy as
  `../results/`. Bulky outputs (BAMs, assembled contigs) stay on the VM; written-up numbers
  come back here.
- **`slides/`** — the progress decks for this direction.

## Machine

New, separate VM from the HLA one (`HLAcalling_pilot_v0_m`) — own siloed app instance,
matching the existing two-collaborator pattern. TRUST4 is CPU/IO-bound, no GPU, far lighter
than the long-read HLA stack (no DeepVariant/sawfish/pbsv/sniffles). Recommended spec:
4-8 vCPU, 16-25 GB RAM, ~100 GB disk, Ubuntu 24.04 — see chat log for full reasoning.

## Status checklist

- [x] VM created (`n1-highmem-8`, 8 vCPU / 52 GB RAM, 100 GB disk)
- [x] repo cloned on VM, gcsfuse mounted (same recipe as the HLA workstream)
- [x] pixi env built (`pixi install` against this folder's `pixi.toml`)
- [x] RNA-seq manifest located — `v9/multiomics/rnaseq/manifest.tsv`
- [x] TRUST4 reference files fetched, vendored into `reference/`
- [x] one participant's RNA BAM path resolved (1000291)
- [x] TRUST4 run on that one sample — needed `--abnormalUnmapFlag`, now fixed in the script
- [x] CDR3 count recorded in `results/` — **2,013, feasibility confirmed**
- [x] real 100-person cohort built (17/17/17/17/16/16 across ancestry groups)
- [x] `--jobs 1` vs `--jobs 3` timed on 3 people — 2.55x speedup, near-ideal, no CPU contention
- [x] `--jobs` efficiency sweep (4/6/8/12/16), overnight — non-monotonic, population
  noise between the different test slices dominated the signal; no evidence of a hard
  ceiling up to 16. All configs beat sequential by 1.9–4.1x. See `run_rnaseq_jobs_sweep.sh`.
- [x] full 100-person cohort completed at `--jobs 12`
- [x] **ancestry recovery check — done, and it's uneven.** AFR highest (7,406 mean
  CDR3s), EAS lowest (4,538) — 1.63x spread, consistent across all 5 major chain types.
  Full writeup: `results/rnaseq_100person_ancestry_writeup.md`.
- [x] **depth confound checked** — explains ~40% of the gap (1.63x -> 1.38x
  normalized), not all of it. AFR/EAS stay the extremes either way. Not fully closed.
  `results/rnaseq_depth_confound_check.csv`.
- [x] **remaining confounds run** — RQS (r=0.061) and cell-composition proxy (flat across
  ancestry) both ruled out. **Kruskal-Wallis H=7.62, p=0.18 — not significant.** The
  ancestry pattern is real in the data but does not clear the bar at this n.
- [ ] **five further QC confounds, never checked** — library complexity, 3' bias, rRNA
  rate, globin rate, expression profiling efficiency. All already computed by AoU and
  sitting in `v9/multiomics/rnaseq/rnaseqc2/`. Zero compute. See data report item **N4** —
  highest value-per-minute check available.
- [ ] **switch recovery metric to rarefaction** rather than post-hoc depth normalization —
  see TRUST4 deep dive §8.5. Probably the most important methodological fix outstanding.
- [ ] **LR × RNA-seq overlap count** — `scripts/check_lr_rnaseq_overlap.py`, written and
  ready, ~1 min on the VM. **Gates the whole disease study.** Expected ~240–530 under
  independent sampling; likely more, and likely AFR-skewed.
- [ ] disease study steps 2–3 (notebook + analysis) — see `reference/LR_RNASEQ_DISEASE_STUDY.md`
- [ ] attach HLA labels (AoU-native + our own long-read calls) once confounds are resolved
- [ ] second 100-person cohort on a **different machine** (decided 2026-08-11) — real
  machine-to-machine comparison, and doubles ancestry-group sample size (16-17 -> ~33-34
  per group). Use `build_rnaseq_cohort.py --total 100 --skip 17` on the new machine so
  it's a genuinely different 100 people, not a repeat of the first cohort.
