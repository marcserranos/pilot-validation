# Current Sprint Context — Omni-HLA Genomics

*Companion to TASK_CONTEXT.md. That file is the permanent project framing; this one is the live snapshot of where we are right now. Rewritten 2026-07-07 to consolidate a long trial-and-error session into a compact reference — see git history if the blow-by-blow is ever needed.*

## Status (2026-07-08, end of session)

Full pipeline (data access → chr6 region slice → FASTQ → typing → comparison) is proven end-to-end on real AoU data. Person 1017156's complete 3-way bake-off (AoU-native, SpecHLA-SR, SpecImmune-LR) is done: 6/8 classical genes concordant, 2 real discordances found — detail in `SMOKE_TEST_PICKS.local.md` (gitignored, real genotype data, never goes here). Two more smoke-test people (2522883, 1253627) have SpecImmune-LR done; SpecHLA-SR for both is the next step, queued but not yet run as of end of session.

**Pick up here next session** — literal next commands, nothing else needs rediscovering:
```bash
cd ~/repos/pilot-validation
pixi shell -e spechla
cd ~/tools/SpecHLA
{ time bash script/whole/SpecHLA.sh -n 2522883 -1 ~/pipeline_outputs/2522883/R1.fastq.gz -2 ~/pipeline_outputs/2522883/R2.fastq.gz -o ~/pipeline_outputs/2522883/spechla_output/ ; } 2> ~/pipeline_outputs/2522883/spechla_timing.txt
{ time bash script/whole/SpecHLA.sh -n 1253627 -1 ~/pipeline_outputs/1253627/R1.fastq.gz -2 ~/pipeline_outputs/1253627/R2.fastq.gz -o ~/pipeline_outputs/1253627/spechla_output/ ; } 2> ~/pipeline_outputs/1253627/spechla_timing.txt
python3 ~/repos/pilot-validation/compare_hla_results.py 2522883
python3 ~/repos/pilot-validation/compare_hla_results.py 1253627
```
Then: look at cross-person patterns across all 3 people's 3-way tables before drawing any conclusion about which method(s) to carry to n=25/100 — n=1 wasn't enough, n=3 might not be either, but it's the next real signal.

**One open thread, not urgent but not forgotten**: `--align_method_1 minimap2` gives SpecImmune a 3.25x speedup but caused a real DPB1 accuracy regression on person 1017156 (see `COMPUTE_OPTIMIZATION_LOG.md`) — unresolved whether that's a one-off or systematic. Don't adopt it as default without testing on another person first.

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
12. **Pipeline outputs live under `~/pipeline_outputs/<person_id>/`, not `/tmp`** — keep it that way. `/tmp` isn't visible in the Jupyter file browser (which is confined to `$HOME`) and reads as disposable/invisible scratch space; `~/pipeline_outputs/` is the same persistent VM disk, just somewhere you can actually see it.
13. **The `has_lr_whole_genome_variant` BigQuery flag does not guarantee a usable aligned BAM exists.** Hit this directly (2026-07-08): two smoke-test candidates (person 2037879, 1413564) both passed the flag and had a `grch38_bam` manifest entry, but the path didn't resolve to a real file. Root cause: their `platform` column says `revio`, but the manifest's `grch38_bam` path string contains `pacbio` instead — and that `pacbio` folder only contains de novo assembly + variant-call outputs (`assembly/`, `single_sample_vcf/`), no aligned BAM at all. **Fix: before picking smoke-test/scale-up candidates, filter to `grch38_bam` paths that literally contain `/revio/`** (confirmed to reliably have a real aligned BAM, e.g. `pooled/longreads/v8_delta/<center>/revio/bam/<id>/GRCh38/<id>.bam`) rather than trusting the eligibility flag alone. `df[df["grch38_bam"].str.contains("/revio/", na=False)]` in the notebook. This will matter more, not less, at n=25/100 scale — build this filter into candidate selection from the start rather than discovering broken paths one at a time.

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

- **2026-07-08 — full pipeline proven end-to-end on real AoU data, one person's 3-way bake-off complete, two more in progress.**
  - Built and confirmed the real access chain: gcsfuse mount (quirk #11) → chr6 region slice (`chr6:29,500,000-33,500,000`, CRAM needs `-T` against a locally-downloaded reference at `~/ref/Homo_sapiens_assembly38.fasta`, BAM doesn't) → FASTQ conversion. Wrapped into a reusable script (`slice_and_fastq.sh`) so this never needs re-deriving.
  - **Person 1017156: full 3-way bake-off complete** (AoU-native, SpecHLA-SR, SpecImmune-LR). 6/8 classical genes concordant across all three; 2 real discordances found, one suggesting AoU-native might be the outlier (2 independent methods agree against it), one suggesting SpecImmune-LR might be (it alone disagrees with both others). n=1 — not enough to conclude anything about method superiority yet. Full detail in `SMOKE_TEST_PICKS.local.md` (gitignored).
  - Tested a SpecImmune speed optimization (`--align_method_1 minimap2` vs default `bwa`): 3.25x faster, but caused a real per-locus accuracy regression (DPB1) on this person — logged in `COMPUTE_OPTIMIZATION_LOG.md`, not adopted as default pending more testing.
  - **Two originally-picked smoke-test candidates (2037879, 1413564) turned out to have no real aligned lrWGS BAM delivered** despite passing the BigQuery eligibility flag — see quirk #13. Replaced with 2522883 and 1253627, both confirmed to have real files before committing compute time to them.
  - Built `compare_hla_results.py` to auto-generate the 3-way comparison table per person (reads the AoU-native TSV + both tools' result files directly, no manual transcription) — removes a source of error now that we're repeating this across people.
  - **2522883 and 1253627: SpecImmune-LR done** (bwa config, 49m51s and 48m21s respectively — consistent with person 1017156's 52min, suggesting that baseline is typical, not a fluke). **SpecHLA-SR for both is the literal next step — see "Pick up here" above.**
  - Moved pipeline outputs from `/tmp` to `~/pipeline_outputs/<person_id>/` (quirk #12) so they're visible in the Jupyter file browser.

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
   - **Smoke-test individuals: final 3 as of 2026-07-08.** Picked via an ad-hoc BigQuery query (same `has_whole_genome_variant`/`has_lr_whole_genome_variant` flags as the named cohort below) against the full SR+LR-eligible pool (15,400 candidates), *not* drawn from the named `HLA_Pilot_Matched_SR_LR_100` cohort — decided this is fine for a 3-person smoke test since ancestry stratification only matters for the later n=25/100 run, which will use the named cohort. Two originally-picked candidates were swapped out for lacking a real delivered BAM (see quirk #13) — filter to `/revio/`-path manifest entries for any future picks. **Actual `person_id`s and AoU-native HLA calls are kept in `SMOKE_TEST_PICKS.local.md` (gitignored, not committed) — see that file, not this one, for the real values.** This repo is public on GitHub; per the repo-scope rule above, participant-level data (even just IDs paired with genotype calls) does not belong in a tracked file.
   - Run all three methods above on the same 2-3 individuals, sliced to chr6 HLA region only for the two compute-heavy methods.
   - Compare allele calls across all three: discordance rate per gene, per ancestry group, at 2-field vs 4-field resolution. **Terminology decision (2026-07-08): call SpecImmune-LR the "presumed most-accurate reference," not "ground truth"** — there's no independently-validated clinical-grade typing set in this pilot to check any method against; long reads are presumed better for this hyper-polymorphic region, not proven-correct in an absolute sense. Keep that distinction explicit in any write-up.
   - Once trade-offs (accuracy proxy, runtime, compute cost) are understood from the bake-off, decide which method(s) to carry to the larger n=25/100 stratified run.
   - Stratify sample picks across ancestry superpopulations on purpose for that larger run — and prefer AoU's own genetically-inferred ancestry TSV (PCA-based, AFR/AMR/EAS/EUR/MID/SAS/OTH) over the self-reported race field the original cohort query used (see Open Questions).
3. **Metaparameter characterization** (started 2026-07-08): coverage, long-read vs short-read error profiles, region padding effect on accuracy/runtime, compute efficiency levers, workspace/cloud quirks (this doc now covers the last one). **Runtime/CPU-time measurements now tracked in `COMPUTE_OPTIMIZATION_LOG.md`** — first real data point captured (SpecImmune default config, 52min wall-clock for one person's chr6-sliced HiFi data on the 4-vCPU test VM), plus an untested optimization hypothesis (swap `--align_method_1` to `minimap2` for the read-binning step, since it's currently defaulting to short-read-tuned `bwa` against ~17kb HiFi reads).
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
