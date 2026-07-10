# Environment & Operations

> **Role:** how the machine actually works — VM/repo layout, the pipeline runbook, confirmed data paths, and every Workbench quirk we've paid for in debugging time.
> **Edit:** append-mostly. Add a quirk or fix the moment you find it; don't rewrite existing ones. This is the crown-jewel operational memory — losing a quirk here costs real time to rediscover.
> **Read:** before running anything on the VM.

## Repo & VM layout

- **This repo** (`pilot-validation`, public on GitHub `marcserranos/pilot-validation`): docs + `pixi.toml` + two pipeline scripts. No participant data; no code touching raw data beyond the two slicing/compare helpers.
- **Workbench VM** `HLAcalling_pilot_v0_m` (Standard, 4 vCPU / 25 GB RAM / 88 GB disk, Ubuntu 24.04, us-central1-a):
  - `~/repos/pilot-validation` — this repo (manually `git clone`d; auto-clone never worked — quirk #3).
  - `~/tools/SpecHLA` — upstream, built via `bash index.sh` in `pixi shell -e spechla`. Has one local source patch not in git (see Fixes).
  - `~/tools/SpecImmune` — upstream; HLA reference db built locally (see Fixes).
  - `~/ref/Homo_sapiens_assembly38.fasta` (+ `.fai`, `.dict`) — hg38 reference, downloaded once (~3 GB) because remote `-T gs://…` doesn't work (quirk #11).
  - `~/mnt/aou-controlled` — gcsfuse mount of the genomics bucket (quirk #11).
  - `~/pipeline_outputs/<person_id>/` — all pipeline outputs (quirk #12).
  - Two pixi envs (`spechla`, `specimmune`) in one `pixi.toml`; activate with `pixi shell -e <name>` **from inside `~/repos/pilot-validation`** (quirk #1).

## Pipeline runbook (one person, end to end)

Prereqs: gcsfuse mount active (quirk #11), reference at `~/ref/`, both pixi envs built.

```bash
# 0. Mount the bucket if not already (quirk #11)
mkdir -p ~/mnt/aou-controlled
gcsfuse --billing-project wb-glacial-potato-8710 --implicit-dirs vwb-aou-datasets-controlled ~/mnt/aou-controlled

# 1. Slice both technologies to the chr6 HLA region + convert to FASTQ (samtools lives in the spechla env)
cd ~/repos/pilot-validation && pixi shell -e spechla
bash slice_and_fastq.sh <id> <cram_path_rel_to_mount> <lr_bam_path_rel_to_mount>
#   paths e.g.:  pooled/wgs/cram/v8_base/wgs_<id>.cram
#                pooled/longreads/v8_delta/BI/revio/bam/<id>/GRCh38/<id>.bam

# 2. SpecHLA (short-read) — from ~/tools/SpecHLA, spechla env
{ time bash script/whole/SpecHLA.sh -n <id> \
    -1 ~/pipeline_outputs/<id>/R1.fastq.gz -2 ~/pipeline_outputs/<id>/R2.fastq.gz \
    -o ~/pipeline_outputs/<id>/spechla_output/ ; } 2> ~/pipeline_outputs/<id>/spechla_timing.txt

# 3. SpecImmune (long-read) — from ~/tools/SpecImmune, specimmune env
#    bwa config (proven; see DECISIONS re: minimap2) + visualization disabled (cosmetic, missing converter)
{ time python3 scripts/main.py -n <id> -o ~/pipeline_outputs/<id>/specimmune_output/ \
    -j 4 -y pacbio-hifi -i HLA -r ~/pipeline_outputs/<id>/LR.fastq --db ./db --visualization "" ; } \
    2> ~/pipeline_outputs/<id>/specimmune_timing.txt

# 4. 3-way comparison table (writes ~/pipeline_outputs/<id>/comparison.md; paste that back to fold into .local.md)
python3 ~/repos/pilot-validation/compare_hla_results.py <id>
```

Region window: **`chr6:29,500,000-33,500,000`** (4 Mb — classical MHC + ~440 kb pad each side). See DECISIONS / EXPERIMENTS for the open padding-narrowing experiment.

> Note: `slice_and_fastq.sh`'s *printed* SpecImmune hint omits `--visualization ""`; use the runbook form above — disabling it is the standing decision (DECISIONS).

## Confirmed data locations (v9 CDR, `C2025Q4R6`)

Genomics bucket: **`gs://vwb-aou-datasets-controlled/`** — **requester-pays**: every `gsutil`/`gcloud storage` call needs `-u wb-glacial-potato-8710` (your workspace project bills), and gcsfuse needs `--billing-project wb-glacial-potato-8710`.

| What | Path (under the bucket) | Schema |
|---|---|---|
| srWGS CRAM manifest | `v9/wgs/cram/manifest.csv` | `person_id,cram_uri,cram_index_uri` |
| lrWGS manifest | `v9/wgs/long_read/manifest.tsv` | `research_id,center,platform,grch38_bam,grch38_bai,…` |
| AoU-native srWGS HLA calls | `v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv` | `research_id` + 30 gene-pair cols (`<gene>_1/_2`), 2–3 field |
| AoU genetic-ancestry predictions | `v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv` | `research_id`, `ancestry_pred` (AFR/AMR/EAS/EUR/MID/SAS, excl. "other"), `probabilities` (array, ordered AFR/AMR/EAS/EUR/MID/SAS), `pca_features` |

- Always resolve file paths via the `v9/` **manifests** — never hand-build into `pooled/` (physical files live under `pooled/wgs/cram/{v7,v8,v9}_base/` and `pooled/longreads/v9_delta/`, stored incrementally across releases).
- **Wrong-bucket trap:** legacy `gs://fc-aou-datasets-controlled/…` (old Firecloud naming) appears in archived docs — **not** ours.
- AoU-native HLA typing = **HLA-HD + Polysolver + OptiType** ensemble on srWGS CRAMs; **no lrWGS equivalent exists**. Restrict comparisons to the **8 classical genes** our tools call (A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1) — AoU types 30.
- lrWGS is delivered as **BAM** aligned to **`grch38_noalt`** (no ALT/HLA contigs — compatibility caveat in DECISIONS). BAM needs no `-T`; only srWGS CRAM does.
- **Path correction (2026-07-09):** the HLA calls live under `snpindel/aux/hla_variants/`, **not** directly under `snpindel/hla_variants/` as earlier notes (incl. `reference/AOU_DATA_ACCESS_NOTES.md`) recorded — that wrong path is why the comparison failed at session start until traced via `gsutil ls`. Table above is corrected; the reference doc is left as its dated historical snapshot.
- Full data-access provenance (web research, confidence markers, sources) archived in `reference/AOU_DATA_ACCESS_NOTES.md`.

## Workbench quirks & gotchas — read before doing anything new

1. **`pixi shell` only works from the manifest's directory.** `cd ~/repos/pilot-validation` first, always. Elsewhere it fails silently into a *plain, non-activated* shell (no error) and every tool command then breaks with confusing "command not found". Check for the `(omni-hla-pilot:<env>)` prompt prefix before trusting the shell.
2. **`pixi shell` sometimes misses its 3-second hook window.** It prints a fallback line like `. /tmp/pixi_env_XXX.sh` — run that line manually to finish activation. **Corollary:** any commands pasted in the *same block* as the `pixi shell` line get swallowed during this hiccup — paste run commands *separately*, only after the `(omni-hla-pilot:<env>)` prompt appears.
3. **Workbench's "Add repository" auto-clone never worked** this project — registering the repo, restarting the app, making it public: none cloned it into `~/repos/`. Don't debug it; just `git clone` manually.
4. **`Failed to get AoU version: no AoU resources in workspace`** prints on every fresh shell — benign noise, unrelated to data access (a real BigQuery query against the CDR returns real rows despite it).
5. **The CDR is not exposed via shell env vars** (`env | grep -i cdr` → nothing). It lives in a separate GCP project (`wb-silky-artichoke-2408`, dataset `C2025Q4R6`) and is accessed via BigQuery with a fully-qualified table path — the Cohort/Dataset Builder auto-generates this in exported notebook snippets (the intended access pattern).
6. **No sudo** (`jupyter` user not in sudoers) — can't apt-install anything. Work around missing system tools rather than installing (see the `less` fix).
7. **samtools in the `spechla` env (v1.21) has built-in GCS support** (`gs+http`/`gs+https`/`gs` handlers) — but it does **not** cleanly handle this requester-pays bucket (quirk #11); use the gcsfuse mount instead.
8. **Two collaborators, siloed apps, no conflict.** Aleix (`aleixruf@researchallofus.org`) works in his own app instance in the same workspace. Check his `hla-pilot-cohort-snapshot` resource before redoing the stratified ~100-person selection — it may already exist.
9. **Always confirm you're reading the current CDR release, not an archived one.** AoU republishes core docs per CDR version and keeps old ones live, marked "(Archived …)". Confirmed 2026-07-08: **v9 (`C2025Q4R6`) is current.** Re-check the "Current CDRs" table before trusting any cited fact from an older session.
10. **The platform migrated** from the original AoU Researcher Workbench (Terra/Firecloud) to **Verily Workbench** — visible in the `vwb-` vs legacy `fc-` bucket split and UI differences. Any doc/tutorial/memory using old-platform terminology needs live re-verification; when something doesn't match the screen, suspect staleness first, not user error.
11. **Don't use raw `gs://` URLs with samtools/htslib against the genomics bucket — use a `gcsfuse` mount.** htslib doesn't reliably attach billing to requests against this requester-pays bucket (`GCS_OAUTH_TOKEN` alone → `Invalid argument`; adding `GCS_REQUESTER_PAYS_PROJECT` → `Operation not permitted`; neither works). The confirmed working method:
    ```bash
    mkdir -p ~/mnt/aou-controlled
    gcsfuse --billing-project wb-glacial-potato-8710 --implicit-dirs vwb-aou-datasets-controlled ~/mnt/aou-controlled
    ```
    **`--implicit-dirs` is not optional** — without it nested paths silently don't resolve (`ls` looks empty) because GCS has no real directories and this bucket has no directory-placeholder objects. This gcsfuse (v3.10.0, preinstalled) daemonizes automatically — no second terminal. Then `~/mnt/aou-controlled/…/<file>.cram` works with plain `samtools view`. **When debugging remote/mount access, always redirect output to a file and check `$?` — never pipe through `tail`/`grep` while diagnosing; a filtered pipe once silently hid a real error behind stale-looking success output.** Public buckets (e.g. the Broad reference bucket) don't have this problem — plain `gs://` works for those.
12. **Pipeline outputs live under `~/pipeline_outputs/<person_id>/`, not `/tmp`.** `/tmp` isn't visible in the Jupyter file browser (confined to `$HOME`) and reads as disposable; `~/pipeline_outputs/` is the same persistent disk, just visible.
13. **The `has_lr_whole_genome_variant` BigQuery flag does not guarantee a usable aligned BAM.** Two candidates (2037879, 1413564) passed the flag and had a `grch38_bam` manifest entry, but the path didn't resolve — their `platform` said `revio` while the path string said `pacbio`, and that `pacbio` folder holds only assembly + variant-call outputs, no aligned BAM. **Fix: filter candidate `grch38_bam` paths to those literally containing `/revio/`** (`df[df["grch38_bam"].str.contains("/revio/", na=False)]`) — reliably a real BAM. Build this into candidate selection from the start; it matters more at n=25/100 scale.
14. **The VM auto-stops after ~1 hour of inactivity — assume it has been restarted at the start of every session.** Consequences, and a large fraction of session-start errors trace to them:
    - **The gcsfuse mount does not survive a stop/restart** — it's gone every new session. Remount (quirk #11) *before* anything that reads bucket data. A start-of-session `FileNotFoundError` on `~/mnt/aou-controlled/…` almost always means "the VM slept, remount," not a real path problem.
    - **After remounting, the FUSE layer needs a beat to become ready.** gcsfuse prints `File system has been successfully mounted` but a command fired in the *same paste* can still race ahead and get `FileNotFoundError`. Run the mount on its own, then **verify with `ls` of a known path** (which also warms the implicit-dir metadata cache) before running dependent commands — don't chain the mount and the consumer in one paste (same lesson as quirk #2).
    - **Running processes and `/tmp` don't survive; the persistent home disk does.** `~/pipeline_outputs/`, `~/repos/`, `~/tools/`, `~/ref/` all survive the restart fine — only the mount, background processes, and the activated pixi shell need redoing.
16. **Check `reference/` before doing any external research (web fetch, search, browser) for AoU documentation.** The 59-page "How the All of Us Genomic Data are Organized v9.pdf" is already saved in `reference/`, and `reference/AOU_DATA_ACCESS_NOTES.md` already has a running log of what's been extracted from it (bucket paths, schemas, confidence markers). **Confirmed 2026-07-10:** re-derived the genetic-ancestry TSV's location and schema entirely via live web search + browser automation (support.researchallofus.org search, a blocked 403 fetch, then manually paging through the same PDF via Chrome), when `AOU_DATA_ACCESS_NOTES.md` section 6 already had the schema written down from this exact PDF, sourced in a prior session — the only missing piece (the literal bucket path) took real time and tokens to re-find that a `grep reference/AOU_DATA_ACCESS_NOTES.md` and a `Read` of the local PDF would have gotten in seconds. **Fix, now standard practice: `grep` (or read) `reference/AOU_DATA_ACCESS_NOTES.md` and check for a matching PDF/doc in `reference/` first — only go external if the local notes don't cover it.**
17. **SpecImmune must always run inside the `specimmune` pixi env — never assume the ambient shell's env doesn't matter.** Its long-read pipeline shells out to `sniffles` (SV caller), which is only on `PATH` inside `specimmune`, not `spechla` (which has `samtools` but not SpecImmune's own deps). Running `python3 scripts/main.py` from a `spechla`-activated shell fails deep inside the pipeline (`sniffles: command not found`) on every gene, cascading into missing-VCF errors — **but `main.py` still exits 0**, so a naive `$?` check reports success while producing zero real output. Confirmed 2026-07-09: an entire unattended 9-config sweep silently produced only NA rows this way. **Fix, now standard practice: invoke SpecImmune via `pixi run --manifest-path ~/repos/pilot-validation/pixi.toml -e specimmune -- python3 scripts/main.py ...`** instead of relying on `cd`-and-hope — this pins the environment regardless of what shell launched the script. **Also: after any automated SpecImmune call, verify the actual expected output file exists — never trust exit code alone** (`run_aligner_pad_sweep.sh` does both now).

## Fixes already made — don't rediscover these

- `pixi.toml`: added `arpack` (SpecHLA's SpecHap cmake build needs it; missing from upstream's `environment.yml`), `pandas` (SpecHLA's `g_group_annotation.py`), `pypdf2` (SpecImmune PDF-report step).
- `less` couldn't be added via pixi (solver conflict: conda-forge `less` needs a newer `libzlib` than `blast`/`bwa` allow) and can't be apt-installed (no sudo). **Patched SpecHLA source instead** — `less $fa | grep -v ">"` is just a `cat` misuse: `sed -i 's/less \$fa/cat \$fa/g' ~/tools/SpecHLA/script/whole/annoHLA.pl`. **This patch lives only on the VM, not in git** — reapply if SpecHLA is re-cloned.
- SpecImmune needs its HLA reference db built locally (not shipped): `python scripts/make_db.py -o ./db -i HLA`, once from `~/tools/SpecImmune` in `pixi shell -e specimmune`.

## Known unresolved technical issues (tracked, non-blocking)

- A second, unlocated `less` call fires during SpecHLA's long-read phasing (seen in HLA_DRB1 phasing in the hybrid test) — didn't block, low priority.
- Real segfaults hit SpecHLA's long-read linkage-extraction for some genes (DPB1/DQA1/DQB1/DRB1) in hybrid short+long mode; the pipeline falls back to short-read-only phasing for those and still completes with plausible output — but long-read data may not be contributing phasing info there. Possibly thread-count related (`Requested more threads for alignment (5) than system-wide available (4)`). **Needs investigation before trusting hybrid-mode quality at scale.**
- SpecImmune's SVG→PDF visualization fails (`rsvg-convert`/`inkscape` missing, no sudo). Cosmetic, irrelevant to allele calls — **disabled going forward via `--visualization ""`** (empty string is the only value argparse's buggy `type=bool` actually treats as false).

## Package management

Deps for SpecHLA and SpecImmune are managed with [pixi](https://pixi.sh), not conda — two separate pixi envs (`spechla`, `specimmune`) in one `pixi.toml`, mirroring how upstream ships them as two never-designed-to-coexist conda envs. Deps are loosely pinned rather than a full port of upstream's exact-pin `environment.yml` (SpecHLA's own `CONDA_PACKAGING_PROPOSAL.md` calls its pins outdated/overly strict; SpecImmune's has 400+ auto-exported pins not worth hand-porting). SpecHLA itself isn't conda-packageable — it vendors pre-compiled binaries (`bin/bcftools`, `bin/novoalign`, `bin/fermikit/`) called by hardcoded relative path, so pixi supplies only the surrounding toolchain; the tool is cloned from source and built via `bash index.sh`.
