# Current Sprint Context — Omni-HLA Genomics

*Companion to TASK_CONTEXT.md. That file is the permanent project framing; this one is the live snapshot of where we are right now. Rewritten 2026-07-07 to consolidate a long trial-and-error session into a compact reference — see git history if the blow-by-blow is ever needed.*

## Status (2026-07-07)

Controlled Tier access confirmed working. pixi-based environments for SpecHLA (short-read) and SpecImmune (long-read) are built and validated on the Workbench VM. SpecHLA is fully proven end-to-end against synthetic test data. SpecImmune's environment is installed; the tool itself is mid-setup (building its local HLA reference db) — first real run still pending. Real AoU data has not been touched yet. Decision on record: smoke-test on 2-3 real individuals before scaling to the full ~100-person pilot.

**Next up:** find the AoU WGS CRAM/BAM file-manifest access pattern and the chr6-slicing approach (a research task is running in parallel on this), then run the smoke test against real data for those 2-3 people.

## Repo & environment structure

- This repo (`pilot-validation`, public on GitHub: `marcserranos/pilot-validation`) holds only docs + `pixi.toml`. No participant data and no code that touches real data belongs here — see the data-egress caveat under Open Questions.
- Workbench app `HLAcalling_pilot_v0_m` (Standard VM, 4 vCPU / 25GB RAM / 88GB disk, Ubuntu 24.04, us-central1-a). On that VM:
  - `~/repos/pilot-validation` — this repo (manually `git clone`d, see quirk #3 below)
  - `~/tools/SpecHLA` — upstream SpecHLA, built via `bash index.sh` inside `pixi shell -e spechla`. Has one local source patch not tracked in git (see "Fixes already made" below).
  - `~/tools/SpecImmune` — upstream SpecImmune, HLA reference db build in progress.
  - Two pixi environments in `pixi.toml`: `spechla` and `specimmune`. Activate with `pixi shell -e <name>` — **must be run from inside `~/repos/pilot-validation`** (see quirk #1).

## Workbench quirks & gotchas — read before doing anything new

1. **`pixi shell` only works from the manifest's directory.** `cd ~/repos/pilot-validation` first, always. Run it elsewhere and it fails silently into a *plain, non-activated* shell (no error) — every SpecHLA/SpecImmune command then breaks with confusing "command not found" errors. Check for the `(omni-hla-pilot:<env>)` prompt prefix before trusting the shell.
2. **`pixi shell` sometimes misses its 3-second hook window.** It'll print a fallback line like `. /tmp/pixi_env_XXX.sh` — just run that line manually to finish activation.
3. **Workbench's "Add repository" auto-clone feature never actually worked** in this session — registering the repo, restarting the app, even making the repo public didn't get it cloned into `~/repos/`. Don't waste time debugging it; just `git clone` manually.
4. **`Failed to get AoU version: no AoU resources in workspace`** prints on every fresh shell — this is benign noise, unrelated to real data access (confirmed: a real BigQuery query against the CDR returned actual participant rows despite this message).
5. **The CDR is not exposed via shell env vars** in this Workbench flavor (`env | grep -i cdr` → nothing). It lives in a separate GCP project from the workspace's own (`wb-silky-artichoke-2408`, dataset `C2025Q4R6` as of this writing) and is accessed via BigQuery with an explicit fully-qualified table path — the Cohort/Dataset Builder auto-generates this in exported notebook snippets, that's the intended access pattern.
6. **No sudo** (`jupyter` user not in sudoers) — can't apt-install anything. Work around missing system tools rather than installing them (see the `less` fix below).
7. **samtools in the `spechla` pixi env (v1.21) has built-in Google Cloud Storage support** (`gs+http, gs+https, gs` URL handlers in its htslib build) — meaning it may be able to slice a `gs://`-hosted CRAM/BAM by region directly without downloading the whole file, given a reachable `.crai`/`.bai` index. Relevant for chr6-only slicing; not yet confirmed in practice.
8. **Two collaborators, two siloed apps, no conflict.** Aleix (`aleixruf@researchallofus.org`) works in parallel on his own app instance in the same workspace. Check his `hla-pilot-cohort-snapshot` resource before redoing work — the actual stratified ~100-person ancestry selection may already exist there.

## Fixes already made to `pixi.toml` — don't rediscover these

- `arpack` — required by SpecHLA's SpecHap build (cmake fails without it); missing from upstream's own `environment.yml`.
- `pandas` — required by SpecHLA's `g_group_annotation.py`.
- `less` could not be added via pixi (conda-forge's `less` needs a newer `libzlib` than `blast`/`bwa` allow — unresolvable solver conflict) and can't be apt-installed (no sudo). **Fixed by patching SpecHLA's source directly** — `less $fa | grep -v ">"` is just a `cat`-equivalent misuse, so: `sed -i 's/less \$fa/cat \$fa/g' ~/tools/SpecHLA/script/whole/annoHLA.pl`. **This patch lives only on the VM, not in git** — reapply if SpecHLA is ever re-cloned fresh.
- SpecImmune needs its HLA reference db built locally, it's not shipped in the repo: `python scripts/make_db.py -o ./db -i HLA`, run once from `~/tools/SpecImmune` inside `pixi shell -e specimmune`.

## Known unresolved issues (tracked, non-blocking)

- A second, unlocated `less` call fires during SpecHLA's long-read phasing path (seen during HLA_DRB1 phasing in the hybrid test) — didn't block the run, low priority to chase.
- Real segfaults hit SpecHLA's long-read linkage-extraction step for some genes (DPB1/DQA1/DQB1/DRB1) in hybrid short+long mode. The pipeline falls back to short-read-only phasing for those genes rather than crashing — the run still completes with plausible output, but long-read data may not actually be contributing phasing info for affected genes. **Needs real investigation before trusting hybrid-mode result quality at scale** — possibly thread-count related (`Requested more threads for alignment (5) than system-wide available (4)`).
- SpecImmune has not yet produced a successful run (db build in progress as of this writing).

## Package management

We manage SpecHLA's and SpecImmune's dependencies with [pixi](https://pixi.sh), not conda, as two separate pixi environments (`spechla`, `specimmune`) in one `pixi.toml` — mirroring how upstream ships them as two separate, never-designed-to-coexist conda envs. Deps are loosely pinned rather than a full port of upstream's exact-pin `environment.yml` files (SpecHLA's own `CONDA_PACKAGING_PROPOSAL.md` says its pins are outdated/overly strict; SpecImmune's has 400+ auto-exported pins not worth hand-porting). SpecHLA itself is not conda-packageable today — it vendors pre-compiled binaries (`bin/bcftools`, `bin/novoalign`, `bin/fermikit/`) called by hardcoded relative path, not `$PATH` — pixi supplies only the surrounding toolchain; the tool itself is cloned from source and built via `bash index.sh`.

## Where this sits in the bigger picture

TASK_CONTEXT.md describes the full project: direct HLA allele calling (SpecHLA/SpecImmune) across the All of Us cohort (short-read + long-read sub-cohort, 535,000+ and 14,000+ individuals respectively as of CDRv9), used for ancestry-stratified allele frequency studies and non-linear PRS/epistasis modeling on autoimmune disease phenotypes. The earlier plan to validate on a 1000 Genomes Project pilot before touching AoU was dropped — we're validating directly on AoU's own short-read/long-read overlap population instead.

## Immediate objectives (this sprint)

1. ~~Environment onboarding~~ — done: Workbench, Cohort Builder, BigQuery access, pixi/SpecHLA/SpecImmune setup all working.
2. **AoU small-n pilot (methods validation, in-cohort)** — in progress:
   - ~~Identify AoU individuals with both short-read and long-read WGS~~ — done: cohort `HLA_Pilot_Matched_SR_LR_100` built, 9,358 candidates (SR+LR matched, 6 races selected).
   - Smoke-test on 2-3 individuals first (decided 2026-07-07, given this is the first real run on this setup) before committing to the full ~100.
   - Run SpecHLA + SpecImmune on short-read and long-read data for the same individuals, sliced to chr6 HLA region only.
   - Compare allele calls across technologies: discordance rate per gene, per ancestry group, at 2-field vs 4-field resolution.
   - Stratify sample picks across ancestry superpopulations on purpose once past the smoke test.
3. **Metaparameter characterization** (not yet started): coverage, long-read vs short-read error profiles, region padding effect on accuracy/runtime, compute efficiency levers, workspace/cloud quirks (this doc now covers the last one).
4. **Basic population QC pass** at scale (not yet started).
5. **Cost estimation** for the full-cohort run (not yet started — needs #3 first).

## Open questions

- **Data egress / compliance**: HLA genotypes are highly polymorphic per person (closer to a fingerprint than most variants) — even small derived result files may count as data requiring All of Us egress/publication review before leaving the Workbench, regardless of file size. Not yet confirmed with whoever sponsored Controlled Tier access. Until confirmed, default to doing comparison/visualization *inside* the Workbench, not downloading results locally.
- Exact WGS CRAM/BAM manifest table/schema for chr6 slicing — researched, see `AOU_DATA_ACCESS_NOTES.md`. Headline risk to test first: AoU's own IGV docs say direct `gs://` streaming is blocked by "data exfiltration controls" (workaround: `gsutil cp` locally first) — unclear yet whether this also blocks `samtools view gs://...` run from inside the VM itself, or only affects external viewer tools. First live test should be exactly this, with the `gsutil cp`-then-slice-locally fallback ready.
- Whether 4-field allele resolution is actually required downstream, or 2-field is sufficient.
- Reference panel choice and its effect on miscalls in non-European ancestries.
- Exact compute backend for scaling (Workbench-native likely; GPU resources are probably a mismatch since HLA calling is CPU/IO-bound).

## How to use this doc

Hand this file plus TASK_CONTEXT.md to any new session or agent: TASK_CONTEXT.md for the scientific premise, this file for current status, repo layout, and — critically — the quirks section, so nobody has to rediscover them. Keep rewriting this file compactly as the sprint progresses rather than appending indefinitely.
