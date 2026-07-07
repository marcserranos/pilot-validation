# Current Sprint Context — Omni-HLA Genomics

*Companion to TASK_CONTEXT.md. That file is the permanent project framing; this one is the live snapshot of where we are right now and what's immediately in front of us. Update/replace as the sprint moves.*

## Status
Access to the All of Us Researcher Workbench has just been granted (via UCSC-affiliated credentials, supervisor-sponsored). Registered Tier + Controlled Tier (genomic data) access is in progress. Data is read in-place in the cloud environment — no local download needed. This is a new compute/data environment for us; expect a learning-curve phase before the first real pipeline runs.

## Where this sits in the bigger picture
TASK_CONTEXT.md describes the full project: direct HLA allele calling (SpecHLA/SpecImmune) across the All of Us cohort (short-read + long-read sub-cohort, now 535,000+ and 14,000+ individuals respectively as of CDRv9), used for ancestry-stratified allele frequency studies and non-linear PRS/epistasis modeling on autoimmune disease phenotypes. The earlier plan to validate the pipeline on a small 1000 Genomes Project pilot (public, matched short-/long-read samples) before touching AoU has been dropped — we are moving directly into the AoU dataset itself, using its own short-read/long-read overlap population for methods characterization instead of an external reference cohort.

## Immediate objectives (this sprint)

1. **Environment onboarding** — get comfortable in the Researcher Workbench: Cohort Builder, Dataset Builder, Jupyter/RStudio notebooks, how cloud compute is billed/allocated, how to read data in-place rather than downloading.

2. **AoU small-n pilot (methods validation, in-cohort)**
   - Identify AoU individuals with *both* short-read and long-read WGS via the Cohort Builder/Dataset Builder (the long-read sub-cohort, 14,000+ as of CDRv9).
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

## How to use this doc
This file plus TASK_CONTEXT.md should be handed together to any new LLM/agent session: TASK_CONTEXT.md for the full project's scientific premise, this file for what's actively being worked on right now. Replace this file's content as the sprint progresses rather than appending indefinitely.
