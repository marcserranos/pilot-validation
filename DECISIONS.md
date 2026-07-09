# Decisions & Open Questions

> **Role:** the "why" register. Open questions at the top (not yet settled); resolved decisions with rationale below (so settled calls aren't relitigated or rewritten away).
> **Edit:** move an item Open → Resolved when decided, keeping the rationale and a date. Append; don't delete history.
> **Read:** to understand why the project is shaped this way, or before reopening a settled call.

## Open questions

### Strategic
- **AoU-native callset: trust & sequencing.** AoU already ships a full-cohort srWGS HLA callset (HLA-HD + Polysolver + OptiType) — this reframes the "massive calling is our contribution" pitch. Fork: **(a)** if AoU-native holds up vs. long-read → pivot toward validate-then-build-downstream on the existing calls (faster path to the aims); **(b)** if it shows ancestry-correlated discordance → the calling work becomes *more* necessary and the bias finding is itself novel (every existing AoU-derived HLA study would inherit that bias). **Lean: validate against long-read first** — an ancestry-distribution study on AoU's calls is circular if we don't yet know whether those calls are ancestry-biased in exactly the dimension being studied. Questions for supervisors: were they aware AoU ships this? any known ancestry benchmarks of that specific ensemble? validate-first vs. study-first preference? which outcome do they expect (shapes whether the pilot needs to be stratified/larger to detect ancestry bias vs. just overall concordance)? (raised 2026-07-09).
- **Data egress / compliance.** HLA genotypes are near-fingerprint-identifying; even small derived files may need AoU egress/publication review before leaving the Workbench, regardless of size. Not confirmed with the Controlled-Tier sponsor. Until then: do comparison/analysis *inside* the Workbench, don't download results.
- **Bare `person_id`s in a public repo.** Operational commands (STATUS, EXPERIMENTS) carry bare research_ids. Genotypes are quarantined to the gitignored `.local.md`, but whether bare ids belong in a *public* repo at all is unresolved — flag to the sponsor alongside egress. (Current practice preserved, not expanded; the ids are already in git history.)

### Methodological
- **minimap2 vs bwa for SpecImmune read-binning.** `--align_method_1 minimap2` = 3.25× faster but caused a real DPB1 regression on person 1017156 (EXPERIMENTS). Unknown if one-off or systematic. **Not adopted; bwa stays default** until tested on ≥1 more person, checking DPB1 specifically.
- **4-field vs 2-field resolution** — is 4-field actually required downstream, or is 2-field enough for the frequency/PRS aims?
- **Reference panel choice** and its effect on miscalls in non-European ancestries.
- **`grch38_noalt` compatibility.** lrWGS BAMs are aligned to `grch38_noalt` (no ALT/HLA contigs); SpecImmune's docs want ALT contigs for divergent-allele sensitivity. No workaround (AoU ships no ALT-contig lrWGS). Watch long-read call completeness/confidence vs. short-read in the bake-off.
- **Ancestry variable for the scale-up.** Prefer AoU's genetic-ancestry TSV (PCA-based AFR/AMR/EAS/EUR/MID/SAS, keyed by research_id) over the self-reported race field the original cohort query used, for the stratified n=25/100 pick.
- **Compute backend for scaling** — Workbench-native likely; GPU is a mismatch (HLA calling is CPU/IO-bound). Needs a per-method cost model: AoU-native = $0 egress (pre-computed variant tier); SpecHLA/SpecImmune both pull raw CRAM/BAM = egress cost.
- **Genomic Quality Report (v9)** — has GiaB-based QC/coverage/sensitivity methodology relevant to the coverage/error-profile aim. Only archived v6–v8 versions found so far; look for v9 in the Workbench UI.

## Resolved decisions

- **3-way bake-off, not a 2×2 (2026-07-08).** SpecImmune can't do short-read (argparse `-y` only accepts `nanopore`/`pacbio`/`pacbio-hifi`); AoU has no native long-read caller. So the real comparison is AoU-native-SR vs. SpecHLA-SR vs. SpecImmune-LR — two of the four naive cells don't exist.
- **Skip SpecHLA's long-read mode (2026-07-08).** SpecHLA's own README defers to SpecImmune for long-read — no reason to test a mode its authors don't recommend.
- **"Presumed most-accurate reference," not "ground truth," for SpecImmune-LR (2026-07-08).** No independently-validated clinical typing set exists in this pilot; long reads are presumed better for this hyper-polymorphic region, not proven-correct in an absolute sense. Keep the distinction in any write-up.
- **Smoke test = 2–3 people (2026-07-08).** Not increased despite AoU-native being free — SpecHLA/SpecImmune compute + debugging time is the bottleneck, not sample count.
- **Candidate selection filters to `/revio/` BAM paths (2026-07-08).** The eligibility flag doesn't guarantee a real aligned BAM (ENVIRONMENT quirk #13). Build the filter in from the start.
- **gcsfuse mount, not raw `gs://` (2026-07-08).** htslib doesn't cleanly bill the requester-pays bucket. Recipe in ENVIRONMENT quirk #11.
- **Disable SpecImmune visualization, `--visualization ""` (2026-07-07).** SVG→PDF converter missing, no sudo; visualization is cosmetic, irrelevant to allele calls. Empty string is the only value argparse's buggy `type=bool` treats as false.
- **Pipeline outputs to `~/pipeline_outputs/`, not `/tmp` (2026-07-08).** `/tmp` isn't visible in the Jupyter file browser (quirk #12).
- **Dropped the 1000 Genomes pilot.** Validating directly on AoU's own SR/LR overlap population instead.
- **pixi, not conda** — see ENVIRONMENT / package management for the full rationale.
