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
11. **Don't use raw `gs://` URLs with samtools/htslib against the genomics bucket — use a `gcsfuse` mount instead.** htslib's built-in GCS support doesn't reliably attach billing info to requests against this requester-pays bucket (`GCS_OAUTH_TOKEN` alone → `Invalid argument`; adding `GCS_REQUESTER_PAYS_PROJECT` → different error, `Operation not permitted` — neither actually works). The confirmed working method (2026-07-08): mount with `gcsfuse` (already preinstalled, v3.10.0, no install needed) and treat the result as a normal local path:
    ```bash
    mkdir -p ~/mnt/aou-controlled
    gcsfuse --billing-project wb-glacial-potato-8710 --implicit-dirs vwb-aou-datasets-controlled ~/mnt/aou-controlled
    ```
    **`--implicit-dirs` is not optional** — without it, nested paths silently don't resolve (`ls` on the mount looks empty) because GCS has no real directories and this bucket has no directory-placeholder objects. This gcsfuse version daemonizes automatically (no second terminal needed). Once mounted, e.g. `~/mnt/aou-controlled/pooled/wgs/cram/v8_base/wgs_<person_id>.cram` works with plain `samtools view -H <path>`, confirmed against real data. **When debugging remote/mount access issues, always redirect output to a file and check `$?` — don't pipe through `tail`/`grep` while diagnosing, a filtered pipe silently hid a real error behind what looked like stale duplicate output during this investigation.** Public (non-requester-pays) buckets, like the Broad reference genome bucket, don't have this problem — plain `gs://` URLs work fine for those.

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

## Progress log

- **2026-07-08: first real chr6 region slices succeeded, both srWGS and lrWGS, for smoke-test person 1017156.**
  - srWGS CRAM → `chr6:29,500,000-33,500,000` via the gcsfuse mount (quirk #11) + a locally-downloaded copy of the public hg38 reference (`~/ref/Homo_sapiens_assembly38.fasta` — `-T` against the `gs://` reference path also failed with a different error than the CRAM's, same general class of remote-access flakiness; downloading the ~3GB static reference once was simpler than debugging further). Result: 79MB BAM, 1,219,708 reads.
  - lrWGS BAM (BI/Revio platform, `pooled/longreads/v8_delta/BI/revio/bam/1017156/GRCh38/1017156.bam`, 21.8GB full file) → same region, **no `-T` needed** since BAM isn't reference-compressed (unlike CRAM) — one less moving part. Result: 28.8MB BAM, 3,213 reads.
  - Confirms the full access chain (gcsfuse mount → region slice) works end-to-end for both data types.
  - FASTQ conversion done for both (paired short-read + singletons; single-file long-read).
  - **First real SpecHLA run completed on person 1017156's short-read data** — full pipeline works end-to-end against real AoU data (not just synthetic test data). Result: high concordance with AoU's own native calls across most of the 8 classical genes, with one real discordance found (2nd-field mismatch on one gene) and one minor 3rd-field difference on another. **Actual allele calls and comparison detail live in `SMOKE_TEST_PICKS.local.md` (gitignored) — not here**, per the participant-data policy above. Next: SpecImmune on the matching long-read FASTQ, then the full 3-way comparison for this person, then repeat for the other 2 smoke-test people.

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
   - **Smoke-test individuals selected (2026-07-08).** Picked via an ad-hoc BigQuery query (same `has_whole_genome_variant`/`has_lr_whole_genome_variant` flags as the named cohort below) against the full SR+LR-eligible pool (15,400 candidates), *not* drawn from the named `HLA_Pilot_Matched_SR_LR_100` cohort — decided this is fine for a 3-person smoke test since ancestry stratification only matters for the later n=25/100 run, which will use the named cohort. CRAM URIs and long-read `grch38_bam` URIs confirmed real for all 3 via the manifest join. **Actual `person_id`s and AoU-native HLA calls are kept in `SMOKE_TEST_PICKS.local.md` (gitignored, not committed) — see that file, not this one, for the real values.** This repo is public on GitHub; per the repo-scope rule above, participant-level data (even just IDs paired with genotype calls) does not belong in a tracked file.
   - Run all three methods above on the same 2-3 individuals, sliced to chr6 HLA region only for the two compute-heavy methods.
   - Compare allele calls across all three: discordance rate per gene, per ancestry group, at 2-field vs 4-field resolution. **Terminology decision (2026-07-08): call SpecImmune-LR the "presumed most-accurate reference," not "ground truth"** — there's no independently-validated clinical-grade typing set in this pilot to check any method against; long reads are presumed better for this hyper-polymorphic region, not proven-correct in an absolute sense. Keep that distinction explicit in any write-up.
   - Once trade-offs (accuracy proxy, runtime, compute cost) are understood from the bake-off, decide which method(s) to carry to the larger n=25/100 stratified run.
   - Stratify sample picks across ancestry superpopulations on purpose for that larger run — and prefer AoU's own genetically-inferred ancestry TSV (PCA-based, AFR/AMR/EAS/EUR/MID/SAS/OTH) over the self-reported race field the original cohort query used (see Open Questions).
3. **Metaparameter characterization** (not yet started): coverage, long-read vs short-read error profiles, region padding effect on accuracy/runtime, compute efficiency levers, workspace/cloud quirks (this doc now covers the last one).
4. **Basic population QC pass** at scale (not yet started).
5. **Cost estimation** for the full-cohort run (not yet started — needs #3 first, and now needs a cost model per method in the bake-off, not just per technology).

## Open questions

- **Data egress / compliance**: HLA genotypes are highly polymorphic per person (closer to a fingerprint than most variants) — even small derived result files may count as data requiring All of Us egress/publication review before leaving the Workbench, regardless of file size. Not yet confirmed with whoever sponsored Controlled Tier access. Until confirmed, default to doing comparison/visualization *inside* the Workbench, not downloading results locally.
- **RESOLVED (2026-07-08) — all three file-location questions confirmed live on the VM** via `gsutil -u wb-glacial-potato-8710 ls ...` (note: this bucket is **requester-pays**, always pass `-u wb-glacial-potato-8710` — your own workspace project bills, not the data-hosting project):
  - **srWGS CRAM manifest**: `gs://vwb-aou-datasets-controlled/v9/wgs/cram/manifest.csv` — confirmed columns `person_id,cram_uri,cram_index_uri`.
  - **lrWGS manifest**: `gs://vwb-aou-datasets-controlled/v9/wgs/long_read/manifest.tsv` — confirmed columns include `research_id,center,platform,grch38_bam,grch38_bai,...` (also `longreads.auxiliary_data.tsv` and `longreads.samples_flagged_by_qc.tsv` alongside it).
  - **AoU's own native srWGS HLA calls**: `gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/hla_variants/hla_genotypes.tsv`.
  - Bucket internal structure: `<bucket>/v9/wgs/{cram, short_read, long_read}/...` for the current-CDR-pinned view; a separate `<bucket>/pooled/wgs/cram/{v7_base,v8_base,v9_delta}/...` and `<bucket>/pooled/longreads/v9_delta/...` hold the actual files (incremental storage across releases) — always resolve paths via the `v9/` manifests, never guess into `pooled/` directly.
- **Bucket-naming trap**: older AoU support articles (CDRv5-v8 era, pre-Verily-rebrand "Controlled CDR Directory" pages) reference a **different, legacy bucket**: `gs://fc-aou-datasets-controlled/...` (old Firecloud/Terra-branded platform naming). This is **not** our current bucket — ours is `vwb-aou-datasets-controlled`, confirmed from the current Data Dictionary page. Don't reuse `fc-` paths found in archived articles.
- **RESOLVED (2026-07-08) — headline risk answered.** Direct data access from inside the VM is **not** blocked by AoU's exfiltration controls — that concern (from the IGV docs) turned out to be about external viewer tools, not VM-internal access. The real obstacle was narrower and more mundane: htslib's `gs://` support doesn't handle this requester-pays bucket cleanly. Fixed via a `gcsfuse` mount instead of raw `gs://` URLs — see quirk #11 in the quirks section above for the full recipe.
- **Compatibility flag**: lrWGS BAMs are aligned to `grch38_noalt` (no ALT/HLA contigs), but SpecImmune's docs want hg38 *with* ALT contigs for best sensitivity on divergent alleles. No workaround available (AoU doesn't ship an ALT-contig lrWGS BAM) — just something to watch for in the smoke test (e.g., compare long-read call completeness/confidence against short-read for the same person).
- Cost note: raw CRAM/BAM usage incurs egress charges ("we do not charge egress for variant data... raw data will be more expensive to use") — relevant to the sprint's cost-estimation objective, and now needs to be broken down per method in the bake-off (AoU-native costs $0 in egress since it uses pre-computed variant-tier data; SpecHLA/SpecImmune both require pulling raw CRAM/BAM).
- **All of Us Genomic Quality Report** (per-CDR-version doc, has real QC/coverage/sensitivity-precision methodology from GiaB control samples) — directly relevant to the coverage/error-profile characterization goal (#3 above). Only found archived v6-v8 versions via search so far; worth checking for a current v9 version once back in the Workbench UI.
- Whether 4-field allele resolution is actually required downstream, or 2-field is sufficient.
- Reference panel choice and its effect on miscalls in non-European ancestries.
- Exact compute backend for scaling (Workbench-native likely; GPU resources are probably a mismatch since HLA calling is CPU/IO-bound).

## How to use this doc

Hand this file plus TASK_CONTEXT.md to any new session or agent: TASK_CONTEXT.md for the scientific premise, this file for current status, repo layout, and — critically — the quirks section, so nobody has to rediscover them. Keep rewriting this file compactly as the sprint progresses rather than appending indefinitely.
