# Current Sprint Context — Omni-HLA Genomics

*Companion to TASK_CONTEXT.md. That file is the permanent project framing; this one is the live snapshot of where we are right now and what's immediately in front of us. Update/replace as the sprint moves.*

## Status
Access to the All of Us Researcher Workbench has just been granted (via UCSC-affiliated credentials, supervisor-sponsored). Registered Tier + Controlled Tier (genomic data) access is in progress. Data is read in-place in the cloud environment — no local download needed. This is a new compute/data environment for us; expect a learning-curve phase before the first real pipeline runs.

## Package management
We manage the SpecHLA environment and its dependencies with [pixi](https://pixi.sh), not conda. Upstream SpecHLA docs (see README.HLA.md) only describe conda-based setup (`conda create`/`conda env create -f environment.yml`) — when standing up the environment on the Workbench, translate that to a pixi-managed environment instead of following the conda steps literally. No `pixi.toml` exists in this repo yet; creating one is still open work.



## Where this sits in the bigger picture
TASK_CONTEXT.md describes the full project: direct HLA allele calling (SpecHLA/SpecImmune) across the All of Us cohort (short-read + long-read sub-cohort, now 535,000+ and 14,000+ individuals respectively as of CDRv9), used for ancestry-stratified allele frequency studies and non-linear PRS/epistasis modeling on autoimmune disease phenotypes. The earlier plan to validate the pipeline on a small 1000 Genomes Project pilot (public, matched short-/long-read samples) before touching AoU has been dropped — we are moving directly into the AoU dataset itself, using its own short-read/long-read overlap population for methods characterization instead of an external reference cohort.

## Immediate objectives (this sprint)

1. **Environment onboarding** — get comfortable in the Researcher Workbench: Cohort Builder, Dataset Builder, Jupyter/RStudio notebooks, how cloud compute is billed/allocated, how to read data in-place rather than downloading.

2. **AoU small-n pilot (methods validation, in-cohort)**
   - Identify AoU individuals with *both* short-read and long-read WGS via the Cohort Builder/Dataset Builder (the long-read sub-cohort, 14,000+ as of CDRv9). Done: cohort `HLA_Pilot_Matched_SR_LR_100` built in Cohort Builder, 9,358 candidates (SR+LR matched, 6 races selected).
   - Decision (2026-07-07): since this is the first time running SpecHLA/SpecImmune on this Workbench setup, smoke-test the full pipeline (chr6 slicing → SpecHLA → SpecImmune) on 2-3 individuals first before committing to the full ~100. Cheap to redo if the pixi env, ref genome build, or region coordinates are wrong; avoids burning compute finding that out at n=100.
   - Run SpecHLA on short-read and long-read data for the same ~100 individuals, sliced to the chr6 HLA region only (not whole-genome) for speed.
   - Compare allele calls across technologies: discordance rate per HLA gene, per ancestry group, at 2-field vs 4-field resolution.
   - Stratify sample picks across ancestry superpopulations on purpose (not just default EUR-heavy samples) to test where short-read underperforms most.

3. **Metaparameter characterization** — before scaling, understand and document:
   - **Coverage**: supervisors want as high as feasible, no need to deliberately oversample beyond what's available.
   - **Long-read vs short-read error profiles**: where each technology fails (e.g., Class II genes, paralog confusion, rare alleles, non-European haplotypes).
   - **Region padding**: how much flanking sequence around the core MHC region affects call accuracy and runtime.
   - **Compute efficiency levers**: chr6-slicing before calling, parallelization strategy, coverage subsampling tradeoffs.
   - **Workspace/cloud quirks**: how the Workbench's compute is provisioned and billed, data-reading patterns specific to this environment, any gotchas discovered during onboarding.

4. **Basic population QC pass** on the data we're about to retrieve at scale — sanity-check ancestry composition, coverage distributions, phenotype (EHR) data quality/completeness before committing to the full-cohort pull (535,000+ short-read, 14,000+ long-read as of CDRv9).

5. **Cost estimation** — once the above is characterized, produce a cost projection for the full-cohort (535,000+ short-read, 14,000+ long-read) HLA calling run on Workbench compute.

## Open questions / things still being figured out
- Exact compute backend to use for the pilot (Workbench-native compute is the likely answer; GPU resources like Hugging Face A100/B100 are probably a mismatch since HLA calling is CPU/IO-bound, not GPU-bound).
- Whether 4-field allele resolution is actually required downstream, or whether 2-field is sufficient for the PRS/epistasis modeling — affects compute/accuracy tradeoffs.
- Reference panel choice and its effect on miscalls in non-European ancestries.
- `pixi.toml` drafted (2026-07-07) with two separate pixi environments (`spechla`, `specimmune`), mirroring upstream's two separate conda envs rather than one merged environment — avoids forcing a single solve across two tool stacks that were never designed to coexist. Deps are loosely pinned (not a full port of the upstream exact-pin `environment.yml` files, which are outdated/overly strict per SpecHLA's own `CONDA_PACKAGING_PROPOSAL.md`).
- Both `pixi install -e spechla` and `pixi install -e specimmune` solved and installed cleanly on first try (2026-07-07, on the Workbench VM: 4 vCPU / 25GB RAM / 88GB free disk, Ubuntu 24.04). Sanity-checked spechla's binaries (samtools, minimap2, bedtools, blastn, blat, freebayes, bowtie2, pbmm2, pbsv, longshot) and Python imports (numpy, scipy, Bio, pysam, pulp) all resolve. specimmune throws a harmless warning about pip-installed `pysam` (pulled in transitively by `dysgu`) overwriting the conda-installed `pysam`'s files — cosmetic, not yet known to cause runtime issues.
- SpecHLA cloned from source into `~/tools/SpecHLA` on the VM and `bash index.sh` completed successfully (2026-07-07) inside the `pixi shell -e spechla` environment — bowtie2 indexes for the HLA reference built, SpecHap + ExtractHAIRs compiled via cmake. Needed one dependency fix: `arpack` was missing from the pixi manifest (required by SpecHap's build, per SpecHLA's own `CONDA_PACKAGING_PROPOSAL.md` host deps) — added to `pixi.toml`. `"The installation is finished! Please start use SpecHLA."` — ready for the built-in example/test suite next, before touching real AoU data.
- Workflow note: `pixi shell -e spechla` only works when run from inside `~/repos/pilot-validation` (where `pixi.toml` lives) — running it from another directory silently drops you into a *non-activated* plain shell instead of erroring loudly, causing confusing "command not found" errors downstream. Always confirm the `(omni-hla-pilot:spechla)` prompt prefix is present before running SpecHLA commands. Also: `pixi shell` sometimes fails to hook the current shell within its 3s timeout — when that happens it prints a fallback line like `. /tmp/pixi_env_XXX.sh`; just run that line manually to finish activation.
- SpecHLA is not actually conda/bioconda-installable today (confirmed via its `CONDA_PACKAGING_PROPOSAL.md`): it vendors pre-compiled binaries (`bin/bcftools`, `bin/novoalign`, `bin/fermikit/`) called by hardcoded relative path, not `$PATH`. pixi supplies the surrounding toolchain only; SpecHLA itself is still cloned from source and run via `bash index.sh` / `script/whole/SpecHLA.sh` per README.HLA.md.

## How to use this doc
This file plus TASK_CONTEXT.md should be handed together to any new LLM/agent session: TASK_CONTEXT.md for the full project's scientific premise, this file for what's actively being worked on right now. Replace this file's content as the sprint progresses rather than appending indefinitely.
