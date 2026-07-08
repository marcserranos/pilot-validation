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
9. **Always check you're reading the current CDR release, not an archived one.** AoU republishes core docs (Data Dictionaries, Genomic Quality Report, Controlled CDR Directory, "How data are organized") per CDR version and keeps old ones live, clearly marked "(Archived ...)". Confirmed 2026-07-08: **v9 (`C2025Q4R6`) is current** — everything cited so far is from the current version, not stale. This can change again on a future release — check the "Current CDRs" table on the Data Dictionaries page before trusting any cited fact from an older session.
10. **The whole platform recently migrated** from the original AoU Researcher Workbench (Terra/"Firecloud"-branded) to the newer **Verily Workbench** — visible in things we've already hit directly (the `vwb-` vs legacy `fc-` bucket-naming split, UI navigation differences). Any doc, tutorial, or past-session memory that uses old-platform terminology needs live re-verification against the current UI/infrastructure — don't assume continuity by default. When something doesn't match what's on screen, suspect staleness first, not user error.

## Fixes already made to `pixi.toml` — don't rediscover these

- `arpack` — required by SpecHLA's SpecHap build (cmake fails without it); missing from upstream's own `environment.yml`.
- `pandas` — required by SpecHLA's `g_group_annotation.py`.
- `less` could not be added via pixi (conda-forge's `less` needs a newer `libzlib` than `blast`/`bwa` allow — unresolvable solver conflict) and can't be apt-installed (no sudo). **Fixed by patching SpecHLA's source directly** — `less $fa | grep -v ">"` is just a `cat`-equivalent misuse, so: `sed -i 's/less \$fa/cat \$fa/g' ~/tools/SpecHLA/script/whole/annoHLA.pl`. **This patch lives only on the VM, not in git** — reapply if SpecHLA is ever re-cloned fresh.
- SpecImmune needs its HLA reference db built locally, it's not shipped in the repo: `python scripts/make_db.py -o ./db -i HLA`, run once from `~/tools/SpecImmune` inside `pixi shell -e specimmune`.

## Known unresolved issues (tracked, non-blocking)

- A second, unlocated `less` call fires during SpecHLA's long-read phasing path (seen during HLA_DRB1 phasing in the hybrid test) — didn't block the run, low priority to chase.
- Real segfaults hit SpecHLA's long-read linkage-extraction step for some genes (DPB1/DQA1/DQB1/DRB1) in hybrid short+long mode. The pipeline falls back to short-read-only phasing for those genes rather than crashing — the run still completes with plausible output, but long-read data may not actually be contributing phasing info for affected genes. **Needs real investigation before trusting hybrid-mode result quality at scale** — possibly thread-count related (`Requested more threads for alignment (5) than system-wide available (4)`).
- SpecImmune's `test_HLA.sh` (2026-07-07) completed its actual typing successfully — plausible allele calls produced for every gene tested (HLA-B, DQA1, H, DRA, W). Needed two fixes: `pypdf2` added to pixi (missing dep in the PDF-report step), and the SVG→PDF conversion itself still fails (`rsvg-convert`/`inkscape` not installed, no sudo to fix). Since visualization is a cosmetic reporting step irrelevant to our actual scientific output (allele calls, discordance rates), **decision: disable it going forward rather than chase installing a converter** — pass `--visualization ""` on every real `main.py` call. Note: upstream's `--visualization` flag is defined with `type=bool` in argparse, which is buggy (any non-empty string, including the literal text `"False"`, evaluates truthy) — empty string is the only value that actually disables it.

## Package management

We manage SpecHLA's and SpecImmune's dependencies with [pixi](https://pixi.sh), not conda, as two separate pixi environments (`spechla`, `specimmune`) in one `pixi.toml` — mirroring how upstream ships them as two separate, never-designed-to-coexist conda envs. Deps are loosely pinned rather than a full port of upstream's exact-pin `environment.yml` files (SpecHLA's own `CONDA_PACKAGING_PROPOSAL.md` says its pins are outdated/overly strict; SpecImmune's has 400+ auto-exported pins not worth hand-porting). SpecHLA itself is not conda-packageable today — it vendors pre-compiled binaries (`bin/bcftools`, `bin/novoalign`, `bin/fermikit/`) called by hardcoded relative path, not `$PATH` — pixi supplies only the surrounding toolchain; the tool itself is cloned from source and built via `bash index.sh`.

## Where this sits in the bigger picture

TASK_CONTEXT.md describes the full project: direct HLA allele calling (SpecHLA/SpecImmune) across the All of Us cohort (short-read + long-read sub-cohort, 535,000+ and 14,000+ individuals respectively as of CDRv9), used for ancestry-stratified allele frequency studies and non-linear PRS/epistasis modeling on autoimmune disease phenotypes. The earlier plan to validate on a 1000 Genomes Project pilot before touching AoU was dropped — we're validating directly on AoU's own short-read/long-read overlap population instead.

## Immediate objectives (this sprint)

1. ~~Environment onboarding~~ — done: Workbench, Cohort Builder, BigQuery access, pixi/SpecHLA/SpecImmune setup all working.
2. **Methods bake-off (Stage A) — decide which method(s)/technology(ies) to scale, before scaling anything.** Redesigned 2026-07-08 after realizing the original "just run SpecHLA on SR + SpecImmune on LR" framing skipped a cheaper, more rigorous option. The real, tool-constrained comparison matrix (confirmed via each tool's actual argparse/docs, not assumed):

   | | Short-read | Long-read |
   |---|---|---|
   | **AoU-native** (HLA-HD+Polysolver+OptiType ensemble, already computed by AoU) | ✅ zero extra compute — just need to locate the output file/table | ❌ doesn't exist — AoU has no official lrWGS HLA callset |
   | **SpecHLA** | ✅ proven working | technically has its own long-read mode, but **decision: skip it** — even SpecHLA's own README defers to SpecImmune for long-read, no reason to run a mode the authors themselves don't recommend |
   | **SpecImmune** | ❌ **not supported by the tool at all** (confirmed via argparse: read-type flag only accepts `nanopore`/`pacbio`/`pacbio-hifi`) | ✅ proven working |

   So it's a clean 3-way comparison: **AoU-native-SR vs. SpecHLA-SR vs. SpecImmune-LR**, not a 2×2 grid.
   - ~~Identify AoU individuals with both short-read and long-read WGS~~ — done: cohort `HLA_Pilot_Matched_SR_LR_100` built, 9,358 candidates (SR+LR matched, 6 races selected).
   - Smoke-test/bake-off on **2-3 individuals** (confirmed 2026-07-08, not increased despite AoU-native being free — SpecHLA/SpecImmune compute+debugging time is still the actual bottleneck, not sample count).
   - Run all three methods above on the same 2-3 individuals, sliced to chr6 HLA region only for the two compute-heavy methods.
   - Compare allele calls across all three: discordance rate per gene, per ancestry group, at 2-field vs 4-field resolution. **Terminology decision (2026-07-08): call SpecImmune-LR the "presumed most-accurate reference," not "ground truth"** — there's no independently-validated clinical-grade typing set in this pilot to check any method against; long reads are presumed better for this hyper-polymorphic region, not proven-correct in an absolute sense. Keep that distinction explicit in any write-up.
   - Once trade-offs (accuracy proxy, runtime, compute cost) are understood from the bake-off, decide which method(s) to carry to the larger n=25/100 stratified run.
   - Stratify sample picks across ancestry superpopulations on purpose for that larger run — and prefer AoU's own genetically-inferred ancestry TSV (PCA-based, AFR/AMR/EAS/EUR/MID/SAS/OTH) over the self-reported race field the original cohort query used (see Open Questions).
3. **Metaparameter characterization** (not yet started): coverage, long-read vs short-read error profiles, region padding effect on accuracy/runtime, compute efficiency levers, workspace/cloud quirks (this doc now covers the last one).
4. **Basic population QC pass** at scale (not yet started).
5. **Cost estimation** for the full-cohort run (not yet started — needs #3 first, and now needs a cost model per method in the bake-off, not just per technology).

## Open questions

- **Data egress / compliance**: HLA genotypes are highly polymorphic per person (closer to a fingerprint than most variants) — even small derived result files may count as data requiring All of Us egress/publication review before leaving the Workbench, regardless of file size. Not yet confirmed with whoever sponsored Controlled Tier access. Until confirmed, default to doing comparison/visualization *inside* the Workbench, not downloading results locally.
- Exact WGS CRAM/BAM manifest schema — **confirmed from AoU's own primary-source PDF** ("How the All of Us Genomic and multi-omic data are organized," v9 — Marc has the full 59-page copy in the repo root; `file`'s "8 pages" report on it is a bug in that command, not a real page count), fetched live via a connected browser session (2026-07-07/08). srWGS manifest CSV columns: `person_id,cram_uri,cram_index_uri`. lrWGS manifest columns include `research_id,center,platform,grch38_bam,grch38_bai` (use these, not the `chm13v2.0_*` T2T columns). Genomics bucket confirmed: `gs://vwb-aou-datasets-controlled/` (project `wb-silky-artichoke-2408`) — exact manifest sub-path not yet found, next step is `gsutil ls` on that bucket directly. Full details in `AOU_DATA_ACCESS_NOTES.md`.
- **New: where do AoU's own native srWGS HLA calls live?** Same shape of problem as the CRAM manifest — the org PDF confirms the callset exists and its schema (`sample_id, [Gene]_1, [Gene]_2`) but not a bucket path. Fold into the same `gsutil ls gs://vwb-aou-datasets-controlled/` exploration step.
- **Bucket-naming trap**: older AoU support articles (CDRv5-v8 era, pre-Verily-rebrand "Controlled CDR Directory" pages) reference a **different, legacy bucket**: `gs://fc-aou-datasets-controlled/...` (old Firecloud/Terra-branded platform naming). This is **not** our current bucket — ours is `vwb-aou-datasets-controlled`, confirmed from the current Data Dictionary page. Don't reuse `fc-` paths found in archived articles.
- Headline risk still to test first: AoU's own IGV docs say direct `gs://` streaming is blocked by "data exfiltration controls" (workaround: `gsutil cp` locally first) — unclear yet whether this also blocks `samtools view gs://...` run from inside the VM itself, or only affects external viewer tools. First live test should be exactly this, with the `gsutil cp`-then-slice-locally fallback ready.
- **Compatibility flag**: lrWGS BAMs are aligned to `grch38_noalt` (no ALT/HLA contigs), but SpecImmune's docs want hg38 *with* ALT contigs for best sensitivity on divergent alleles. No workaround available (AoU doesn't ship an ALT-contig lrWGS BAM) — just something to watch for in the smoke test (e.g., compare long-read call completeness/confidence against short-read for the same person).
- Cost note: raw CRAM/BAM usage incurs egress charges ("we do not charge egress for variant data... raw data will be more expensive to use") — relevant to the sprint's cost-estimation objective, and now needs to be broken down per method in the bake-off (AoU-native costs $0 in egress since it uses pre-computed variant-tier data; SpecHLA/SpecImmune both require pulling raw CRAM/BAM).
- **All of Us Genomic Quality Report** (per-CDR-version doc, has real QC/coverage/sensitivity-precision methodology from GiaB control samples) — directly relevant to the coverage/error-profile characterization goal (#3 above). Only found archived v6-v8 versions via search so far; worth checking for a current v9 version once back in the Workbench UI.
- Whether 4-field allele resolution is actually required downstream, or 2-field is sufficient.
- Reference panel choice and its effect on miscalls in non-European ancestries.
- Exact compute backend for scaling (Workbench-native likely; GPU resources are probably a mismatch since HLA calling is CPU/IO-bound).

## How to use this doc

Hand this file plus TASK_CONTEXT.md to any new session or agent: TASK_CONTEXT.md for the scientific premise, this file for current status, repo layout, and — critically — the quirks section, so nobody has to rediscover them. Keep rewriting this file compactly as the sprint progresses rather than appending indefinitely.
