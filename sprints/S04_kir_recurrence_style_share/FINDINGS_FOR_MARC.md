# S04 findings for Marc — KIR recurrence, KIR vs HLA catalogue, figure style, ligand co-occurrence, Cole package, cleanup

*Sprint `s04-kir-recurrence-style-share`, opened 2026-09-26 off S03. Written 2026-09-27, updated
2026-09-28 with the final WS-A/WS-B (v4b) numbers, by the WS-G report writer, from `LOG.md` (the
primary source) and the other sprint docs — not from raw session transcripts.*

## 1. Summary

This sprint answered three of your post-S03 asks (how often do novel KIR/HLA alleles recur across
unrelated people, does KIR's reference catalogue cover this cohort better or worse than HLA's, and
does a KIR receptor pair up with its HLA ligand more than chance by ancestry), delivered a plain-
language explainer + ranked idea list for the repertoire (WS6) work, redid the figure style pass
after you rejected the first attempt, built and uploaded the data-sharing package for Cole, and
proposed (and partly executed) a VM and repo cleanup. The single most important thing to know:
**the KIR-vs-HLA catalogue-coverage numbers went through three internal reversals** as real bugs
were found in how "novel" was being counted, before a fourth, checkpointed VM run (v4b) produced
the final, trustworthy numbers — see section 2 below and "WS-A/WS-B results — final (v4b)" for
what those numbers actually say (short version: neither region is fully catalogued, and once you
look within each ancestry rather than at the pooled cohort, HLA's catalogue is consistently more
complete than KIR's, the opposite of what the earlier, buggy runs suggested). Everything
genotype-derived stayed inside the Workbench; nothing with fewer than 20 people is reported
anywhere.

## 2. Per workstream

### WS-A — KIR recurrence & saturation (scripts 44/45)

**What ran and why:** you asked how often a novel KIR (and, for comparison, HLA) allele shows up
in more than one unrelated person — is KIR's ~59% "novel" rate mostly one-off sequencing/assembly
noise, or the same variant recurring because the reference catalogue is missing it systematically?
Script 44 (VM) exports, per person and gene, whether each called allele is "novel" at three levels
of strictness (genomic sequence, coding-sequence/CDS, and protein) and how many *unrelated* people
carry it, bucketed into recurrence classes (seen by 1 / 2 / >2 / >=20 unrelated people — the >=20
cutoff doubles as the disclosure floor, so that class can be described in more detail).

**Results (final, v4b, 11,856 unrelated people):** neither species' allele discovery curve
saturates in this cohort at any identity level — the "seen once" line dominates every panel and
keeps climbing, for both HLA and KIR, at every novelty definition, pooled and within every
ancestry. Good-Turing incidence coverage (chance a new person's allele has already been seen) is
high (92-99%) for both — a number that sounds like "done," but is a different, less informative
question than catalogue completeness (Chao2, 8-42% depending on species/level — see WS-B). Novelty
is mostly private (seen in only one unrelated person) at the coding-sequence and protein levels for
both species. Full numbers are in the "WS-A/WS-B results — final (v4b)" section below.

**Figures:** `reports/hla_popgen/44_kir_recurrence_saturation/fig_saturation_by_recurrence.png`,
`fig_saturation_by_recurrence_ancestry_{kir,hla}.png`, `fig_saturation_per_ancestry.png`,
`fig_coverage_completeness.png` (+ `_cds.png`, `_protein.png`, `_ancestry_protein.png`; all with
matching `.pdf`).

**Status: done.** v1-v3 were superseded by real bugs (name-based identity, the KIR artifact gap,
the HLA protein-path bug, an HLA gene-name prefix mismatch — see Problems 2/3/4 below); v4b is the
final, sequence-hash-identity, checkpointed run and its numbers are the ones to cite.

### WS-B — KIR vs HLA catalogue coverage (script 46)

**What ran and why:** given 44's per-allele novelty/recurrence data, which reference catalogue
(IPD-KIR for KIR, IPD-IMGT/HLA for HLA) better represents what's actually circulating in this
cohort? Measured with Chao2 (a standard ecology method: estimates how many *total* distinct
alleles likely exist from how many singletons vs. doubletons you've already seen, then compares
that to what you've actually observed).

**Figures:** `reports/hla_popgen/46_kir_vs_hla_catalogue/fig_catalogue_completeness.png` (currently
being redesigned for label-overlap fixes — same filename, will update in place),
`fig_recurrence_composition.png` (+ `.pdf`).

**Status: done.** The early conclusion ("KIR's catalogue is better represented than HLA's, pooled")
does not survive once the identity bug (Problem 2) is fixed and the comparison is done per
ancestry rather than pooled — see below. The corrected, final answer is the opposite in the
comparison that matters (per-ancestry): HLA's catalogue is the more complete one.

#### WS-A/WS-B results — final (v4b)

**Do not cite any v1-v3 recurrence, novelty-rate, Chao2-completeness, or "KIR vs HLA catalogue"
number from this sprint** — those were built on a broken identity definition and a missing
artifact filter (Problem 2/3) and are formally superseded. The numbers below are from the final
run ("v4b"), reconciled two independent ways against the committed tables and reviewed by a
fresh-context critic; use these.

**Headline distinct-allele counts (S_obs), pooled over all ancestries, 11,856 unrelated people:**

| identity level | HLA S_obs | KIR S_obs |
|---|---:|---:|
| genomic (upper bound — see caveat below) | 17,188 | 38,024 |
| CDS (coding sequence) | 1,363 | 2,100 |
| protein | 1,079 | 1,444 |
| protein, novel-only | 198 | 1,063 |

**Key findings:**
- **Neither region is saturated.** At every identity level and in every ancestry, the "seen in
  only one person" line keeps climbing as more people are added — this cohort has not come close
  to exhausting either species' real allele diversity.
- **Genomic-level identity is an upper bound, not a clean measurement**, and should not be used for
  the KIR-vs-HLA comparison. It counts ~79% (KIR) and ~78% (HLA) of alleles as "singletons" and
  produces implausibly large distinct-allele counts (38,024 for KIR, 17,188 for HLA) because it
  can't yet tell apart a genuinely different allele from the same allele assembled with a slightly
  different start/end boundary, or from leftover noise in the non-coding sequence around the gene.
  The coding-sequence (CDS) and protein levels, which don't have this problem, are the levels the
  real comparison is built on.
- **Pooled across all ancestries, HLA and KIR's catalogues look nearly tied** (CDS completeness:
  KIR 41.6% vs HLA 38.8%; protein: HLA 39.9% vs KIR 39.2%) — but **within every one of the 5
  well-powered ancestry groups (African, Admixed American, East Asian, European, South Asian),
  HLA's catalogue is more complete than KIR's**, at both levels, by a consistent 10-25 percentage
  points (e.g. protein-level completeness: South Asian 66.9% HLA vs 40.5% KIR; European 47.7% vs
  28.6%). **This was checked directly and is a genuine statistical effect of combining ancestries
  with different allele pools (the same shape as Simpson's paradox), not an artifact of some
  ancestries having more people than others** — weighting each ancestry's own completeness by its
  sample size reproduces the same ~17-point HLA lead as the raw per-ancestry table, not the pooled
  near-tie. **The per-ancestry result is the one to trust and cite: HLA's catalogue represents this
  cohort's real diversity better than KIR's does**, even though a naive pooled-cohort number would
  suggest they're about the same.
- **Novelty is mostly private (found in only one unrelated person) for both species**, at the
  coding-sequence level (50% of HLA's novel CDS alleles are singletons, 56% of KIR's) and,
  especially, at the protein level (97% of HLA's 198 novel proteins are singletons, vs 69% of
  KIR's 1,063) — i.e. most of what each catalogue is missing looks like ordinary rare variation,
  not a systematic gap the catalogue could easily close. KIR's novel *proteins* recur across
  unrelated people somewhat more often than HLA's (consistent with KIR's known gene-content/copy-
  number diversity independently producing the same derived protein on different haplotypes), but
  neither species' novelty is predominantly shared.
- **Version history:** v1/v2 used a name-based, not sequence-based, "genomic identity" (collapsing
  distinct novel alleles that shared a nearest-reference name) and had two further bugs — HLA's
  protein path pointed at the wrong folder (returning 0 novel proteins for every gene) and KIR's
  protein-catalogue coverage check was silently skipping some bundled reference files. v3 fixed
  those but kept the name-based genomic identity and still lacked KIR's artifact filter, producing
  an internally impossible number (more distinct "novel" KIR proteins than distinct KIR genomic
  alleles) that triggered the final investigation. v4 replaced identity with a true sequence hash
  at every level for both species and added the shared artifact gate; a subsequent full-cohort run
  ("v4b") was killed mid-run by an app stop, and was resolved by adding per-checkpoint saves and
  rerunning to completion (commit `22bf543`) rather than losing partial progress again — this is
  the version behind all numbers in this section.

### WS-C — Figure style pass (rcParams port + redesign + layout linter)

**What ran and why:** port the visual style used in the `cnsplots` reference gallery
(github.com/faridrashidi/cnsplots) into this project's shared plotting code, and re-render the key
figures (Figure 1, DQ G1/G2, KIR saturation, KIR full-cohort) to match.

**Results:** the first pass (fonts, spines, legend styling) was real but, in your words after
reviewing it, "not seeing nearly any difference" — it changed colors/fonts but not the actual
layout problems (misaligned panels, overlapping text, floating colorbars) that made the figures
look unfinished. A second, deeper pass followed: an automated "layout linter" was added to the
shared plotting code that checks a rendered figure for overlapping text, text sitting on top of
other elements, and marginal charts that don't line up column-for-column with the main panel it
sits beside, and the figures were structurally redesigned (not just restyled) against that check.
Three more rounds of full-size visual review (by the orchestrator, since the linter still misses
some real faults — see Problem 3/6) each found and fixed additional issues. Net effect: Figure 1,
the DQ G1/G2 heatmap, the KIR saturation panels, and the KIR full-cohort figure all now pass the
mechanical layout check with zero errors, and multiple real faults invisible to that check (a
reversed bar order, a legend sitting on a title, dark text on a dark cell) were caught by eye and
fixed. Full before/after detail per figure is in
`sprints/S04_kir_recurrence_style_share/FIGURES_INDEX.md`; before/after PNGs are in
`sprints/S04_kir_recurrence_style_share/style_before/` and `style_after/`.

**Status: done.**

### WS-D — WS6 explainer, ranked ideas, naive ML / ligand co-occurrence (scripts 47/48)

**What ran and why:** three things. (1) A plain-language explainer of the WS6 repertoire-prediction
work (`sprints/S04_kir_recurrence_style_share/WS-D_EXPLAINER_AND_IDEAS.md`) — what BenchRep-T does,
what the numbers mean, and why an ancestry-predictability check matters (any future
"repertoire predicts disease" result has to be checked against ancestry as a confound, since
ancestry alone predicts repertoire well). (2) A ranked list of follow-on experiment ideas. (3) Two
new scripts run on the full cohort: 47, a "naive ML" QC battery (can you predict ancestry,
sequencing platform, or KIR gene-content from HLA/KIR carriage alone — used as a batch-effect
check, not a scientific result) and 48, KIR receptor-HLA ligand co-occurrence by ancestry.

**Results:** 47 — ancestry is strongly predictable from HLA+KIR carriage (AUROC 0.81-0.98 across
the six ancestry groups, as expected — HLA is the textbook ancestry-informative locus). Sequencing
platform looked weakly predictable pooled across everyone (AUROC 0.574), but splitting by ancestry
group settles the question: within European-ancestry people it's 0.558 (barely above its own
chance-shuffled range), and within African-ancestry people it's 0.489 (no better than a coin flip),
and a side-by-side model shows ancestry alone already explains almost all of the pooled signal
(ancestry-adjusted platform effect only +0.004) — **the "guess the sequencing machine" result was
mostly "different ancestry groups happened to run on different machines more often," not real
sequencing-machine fingerprinting.** HLA does not predict KIR gene-content class (cA/cB) beyond
chance either (pooled AUROC 0.550), and accounting for ancestry makes that guess very slightly
*worse* (ancestry-adjusted effect −0.023) — about as clean a "no" as this kind of test can give,
and consistent with HLA (chromosome 6) and KIR (chromosome 19) not having to travel together
genetically. 48 — across 30 receptor-ligand-pair-by-ancestry tests, no combination
survived the multiple-testing-corrected threshold; two nominally interesting signals (an East Asian
and a European pairing) did not clear it once you account for testing 30 things at once, so they
are leads, not findings. A side product, HLA-C1/C2 and Bw4/Bw6 epitope frequencies by ancestry, is
a clean, reusable table (e.g., East Asian: 89.5% carry C1, 38.0% carry C2, 75.5% carry Bw4).

**Figures:** `reports/hla_popgen/47_naive_ml_tests/fig_naive_ml_auroc.png`,
`reports/hla_popgen/48_kir_hla_ligand_cooccurrence/fig_kir_hla_ligand_forest.png` (+ `.pdf` for
both).

**Status: done.**

### WS-E — Cole data-sharing package (script 49)

**What ran and why:** build the package Cole needs (HLA calls, KIR calls, a person table, a
README/schema, and checksums) and upload it to the shared bucket.

**Results:** built and uploaded. `gs://hla-calls-share-wb-cordial-leechee-9743/release_2026-09-25/`
now holds 6 files (~10.65 MB total: HLA calls, KIR calls, persons table, README, SCHEMA, MANIFEST
with checksums), verified byte-for-byte against the build. A pre-upload sanity check (numbers you
ran and pasted back) confirmed the package's own KIR novelty rate matches the already-established
S03 number and doesn't carry the same protein-level bug found in 44 (Problem 2) — the two scripts
compute novelty independently, so this is a real, independent cross-check, not a coincidence.

**Status: done** (one cosmetic issue — see Problem 8).

### WS-F — VM and repo cleanup, branch unification

**What ran and why:** you asked for the VM and repo decluttered, and eventually for the four
sprint branches (S03 and earlier work, all fully contained in this branch) unified into one.

**Results:** repo side is done — 37 pre-2026-convention one-off scripts, 3 loose sprint call notes,
and 2 pipeline briefs were moved into clearly-labeled subfolders (`scripts/legacy/experiments/`,
`sprints/_calls/`, `reports/hla_popgen/_briefs/`), with references updated everywhere except the
two append-only history files (`context/EXPERIMENTS.md`, `context/DECISIONS.md`), which got one
new pointer entry each instead, per this repo's own rule against rewriting history. Five files with
real production wiring on the VM were deliberately left in place after a dependency audit found
they're referenced by hardcoded paths (a `git mv` would have silently broken the production
pipeline the next time it runs). One judgment call was surfaced rather than made: the legacy
`outputs/` folder was left untouched instead of being moved and re-tracked in git, because it had
been deliberately `.gitignore`'d in an earlier commit as policy — that's flagged as an open item
for you below, not decided unilaterally.

VM side: a cleanup plan exists (`VM_CLEANUP_PLAN_v2.md`) but the actual deletions require you to
run the commands yourself — see section 6.

**Status: repo reorg done; VM cleanup proposed, not executed** (blocked on the same
permission-classifier issue as the Cole package upload — see Problem 5).

### WS-G — This report + tracking docs

LOG.md, FIGURES_INDEX.md, VM_RUNS.md, the two critic reviews, and this document. **Status: done**
(this file).

## 3. Problems and how they were resolved

Being direct about what went wrong, because several of these changed what's actually true about
the science, not just the presentation:

1. **Figure 1 panel d was silently showing bars against the wrong genes, inherited from S03.**
   The bug (in code from before this sprint) came from combining a heatmap that flips its row order
   with an unlabeled bar chart that flips its row order independently — the two flips didn't
   cancel out, so panel d's novel-protein-allele counts lined up with the wrong gene labels in
   every S03 and early-S04 render. Fixed by making the two panels share one axis so a single flip
   applies to both, with a new automated test that locks the value-to-gene binding going forward.
   **Any copy of Figure 1 v5 you have from before 2026-09-26 has this panel reversed — don't use
   it.**

2. **The name-based "genomic identity" bug — this invalidates the early KIR-vs-HLA catalogue
   conclusions.** This is the most consequential bug of the sprint. The pipeline's "genomic" novelty
   check was, underneath, comparing each call's nearest-reference-allele *name* (e.g., "closest
   match to KIR2DL1*003, flagged new") rather than the actual DNA sequence. Two different novel
   sequences that both happened to be closest to the same named reference allele got silently
   merged into one "recurrent" allele instead of being counted as two distinct novel ones — for
   both species, but with a bigger effect on KIR. Separately, KIR's protein-level check was missing
   an artifact filter that HLA's pipeline already had, so partial or garbled assemblies were being
   counted as "novel proteins" instead of being screened out first. Both bugs were caught the same
   way: an internal sanity check flagged 6,300 distinct "novel" KIR proteins against only 1,460
   distinct KIR genomic alleles — a mathematical impossibility, since a protein-level match can
   only ever be coarser (fewer distinct categories) than a full-sequence match, never finer. That
   contradiction is what triggered the investigation. The fix rebuilds identity as a true sequence
   hash at matched levels for both species and applies one shared artifact-filtering gate to both;
   see "WS-A/WS-B results — final (v4b)" above for the corrected numbers. One more incident on the
   way to the final run: the full-cohort v4b job was killed mid-run when the VM app stopped
   unexpectedly (~33 minutes in, cause not fully confirmed — possibly the account-wide rate limit
   and a Chrome-extension disconnect happening at the same time). Rather than just retrying and
   risking losing progress again, the script was given per-checkpoint saving (every 1,000 people,
   resumable) before the rerun — the successful rerun (commit `22bf543`) finished in one pass and
   didn't end up needing to resume from a checkpoint, but the safeguard is now standard for any
   future long VM job in this pipeline.

3. **The KIR artifact gap.** Related to #2: KIR calls were never being screened for the same
   assembly-quality problems (partial coding sequence, frameshifts, homopolymer-run indels) that
   HLA calls already are. Fixed as part of the same identity rework — same gate, both species.

4. **The HLA protein-path bug.** Independently of #2/#3, HLA's protein-level novelty count came out
   as exactly zero for every one of the 8 classical genes — obviously wrong, since S03 had already
   found real novel HLA proteins. Root cause: 44's export pointed at the wrong folder for HLA's
   per-person coding-sequence files (one path segment short), so the matching step silently found
   nothing to compare against and defaulted every record to "not novel." Fixed; the corrected
   version found 198 distinct novel HLA proteins in the same cohort. A parallel VM run afterward
   found and fixed a second HLA-gene-naming mismatch (files call the gene `HLA-A`, code expected
   bare `A`), which had been silently zeroing out HLA's true genomic-identity level specifically;
   both fixes are in the final v4b run cited above.

5. **Classifier refusals, and the resulting handoffs to you.** Several actions this sprint were
   refused outright by Claude Code's own safety classifier, categorized as things like "data
   exfiltration" or "auto-mode bypass" even though every one of them stayed entirely inside the
   Workbench perimeter: building the Cole data package, uploading it, and navigating to the
   JupyterLab tab to run a pre-upload check. Per this project's own rule, these were not worked
   around — they were handed to you as exact copy-pasteable terminal commands
   (`MARC_UPLOAD_COMMANDS.md`) instead, which you then ran yourself. The VM cleanup (section 6) is
   in the same position and needs the same treatment.

6. **Rate limits interrupted work mid-task twice.** An account-wide session limit killed three
   agents mid-work around midnight on 09-26 (one mid-VM-operation, two with uncommitted figure
   edits) and again on 09-27 before a VM deployment finished. Each time, the next session picked up
   from the interrupted agent's own saved transcript rather than restarting from scratch — no work
   was lost, but it added real wall-clock delay.

7. **The layout linter has real false negatives — it catches most overlap problems but not all.**
   The mechanical check added in WS-C (which flags overlapping text, text on top of chart brackets,
   and misaligned marginal charts) missed several real, visible problems across the sprint: text
   sitting on top of a data point rather than another piece of text (it only ever compared text to
   text), a legend overlapping a panel title (because an invisible axis's placeholder tick labels
   don't actually get drawn, but the check didn't know that), and one instance of a reversed value
   order that no layout check could catch because the layout was structurally fine — only the
   data-to-label pairing was wrong (this is bug #1 above). The lesson written into the sprint's
   process: an automated check plus one human visual pass is still not enough on its own; both a
   layout check *and* a full-resolution visual review *and*, where a figure encodes specific
   numbers, a direct check against the underlying data table are now the standard for any new
   figure.

8. **A cosmetic date mismatch in the uploaded Cole package.** The package's own README says
   "release 2026-09-26" (the date it happened to be built) while it actually lives at the
   `release_2026-09-25/` bucket path. Not a data or disclosure issue, just a documentation
   inconsistency. It was fixed going forward (a `--release-tag` flag now ties the two together) but
   the already-uploaded copy was not re-uploaded just to fix a string — flagged for you to decide
   (section 6).

## 4. Trade-offs and judgment calls

- **KIR's novel-allele protein identity, for calls that don't cleanly match a catalogued protein,
  falls back to a two-field lookup approximation rather than deriving the epitope directly from
  that specific call's own sequence.** Accepted because protein-changing novelty at the relevant
  positions is expected to be rare; documented rather than hidden.
- **The recommended VM cleanup organizes rather than deletes** in almost every case, per your own
  steer mid-sprint ("lower threshold for deletion but a higher [bar] for organization... if
  something is evidently trash or corrupted... that can be deleted"). Only ~0.4 GB of definitely-
  dead material (corrupted backups, one-off patch fragments, regenerable caches, a byte-identical
  duplicate of the already-uploaded Cole package) is proposed for deletion; everything else with any
  conceivable reuse value moves into a new `~/archive/` folder instead of being removed.
- **The repo reorg left `outputs/` alone rather than moving it**, because it had been deliberately
  untracked in an earlier commit — treated as "not this reorg's call," not decided for you.
  Content was checked and contains only aggregate frequency tables, no participant-level data, so
  there's no urgency either way.
- **The permutation-test p-values in script 47 hit a floor** (only 20 shuffles were run in the
  first pass, so the smallest possible p-value is about 0.05) — the report reads the observed
  AUROC against the shuffled-null range directly rather than trusting the p-value number, and a
  follow-up run increased this to 200 shuffles.
- **48's forest plot censors, rather than shows, any receptor-ligand-by-ancestry cell with fewer
  than 20 people** (a grey hatched bar, explicitly labeled, not a zero and not silently dropped) —
  consistent with this project's disclosure floor.

## 5. Cost and time (estimate)

VM: `AoU_Jupyter_ComputeEngine_..._big_run`, `n2-highmem-4`, no resizes this sprint, billed at
roughly $0.37/h while running (per `VM_RUNS.md`). **These are estimates, not a reconciled bill** —
no session logged a precise stop time, so wall-clock windows below are read off `LOG.md`'s
timestamped entries, not a billing record. Three broad sessions:

| Session | Window | Jobs | Est. duration | Est. cost |
|---|---|---|---|---|
| 1 | 09-25 ~21:10 -> 09-26 early morning (ended by an account rate-limit hit, ~00:00) | 47 v1, 48, 44 v1 | ~2h50m | ~$1.05 |
| 2 | 09-26 ~01:10 (restart) -> ~09-26 10:00 | 44 v1 pulls/debugging; a second rate-limit hit occurred within this window (~05:45), so the VM likely auto-stopped and was left stopped for part of it — this is an upper bound on wall-clock, not continuous billed time | up to ~8h50m | up to ~$3.27 |
| 3 | 09-27, multiple restarts across the day | 44 v2, 47 v2, 44 v3, 44 v4 smoke test, 44 v4b (one run killed by an app stop, then a successful ~46-minute rerun after adding checkpointing) | actual job time sums to ~1h35m (608s+1837s+292s+132s+2742s); wall-clock across the restarts and debugging gaps was longer, roughly ~3h | ~$0.59 (job time only) to ~$1.11 (incl. gaps) |

**Rough total: on the order of 13-15 VM-hours across the sprint, roughly $5-6 at $0.37/h** — this
is a coarse upper-bound estimate built from log timestamps and known auto-stop behavior, not a
reconciled bill; `VM_RUNS.md` now carries these rows with the same caveat, and should be replaced
with real start/stop timestamps in a future sprint if a firmer number is needed.

## 6. Open items — exact actions for you to take

1. **VM cleanup.** Run the commands in
   `sprints/S04_kir_recurrence_style_share/VM_CLEANUP_PLAN_v2.md` yourself (agent-side deletes are
   blocked by the permission classifier). Summary: archive `~/results`, `~/allele_geometry_scratch`,
   `~/pipeline_outputs_kir_test` into a new `~/archive/`; delete ~0.4 GB of confirmed-dead material
   (corrupted backups, patch fragments, `__pycache__`, the pip cache, and the already-uploaded
   duplicate `~/s04/share_release_2026-09/`); then drop the `~/VM_LAYOUT.md` map file included in
   that plan.

2. **Push and merge this branch, then delete the subsumed branches.**
   `s04-kir-recurrence-style-share` fully contains `main`, `fig1-drafts-and-research-map`,
   `needle-view-cds-diversity-density`, and `s03-call8-figures-kir-prediction` (confirmed via
   `git merge-base --is-ancestor`) — none of them have any commit this branch is missing.
   Recommended:
   ```
   git push -u origin s04-kir-recurrence-style-share
   git checkout main && git merge --ff-only s04-kir-recurrence-style-share   # falls back to a normal merge if ff fails
   git push origin main
   git branch -d fig1-drafts-and-research-map needle-view-cds-diversity-density s03-call8-figures-kir-prediction
   git push origin --delete fig1-drafts-and-research-map needle-view-cds-diversity-density s03-call8-figures-kir-prediction
   ```
   Do **not** touch `aleix/hla-resolve-phase1` — it's Aleix's own active branch, not part of this
   unification.

3. **Decide on `outputs/`** (legacy frequency/tree CSVs + plots): leave it gitignored where it is,
   fold it into `scripts/legacy/experiments_outputs/` and start tracking it in git (the original
   reorg proposal), or delete it as fully superseded by the numbered pipeline. No disclosure
   concern either way — it's aggregate-only.

4. **Small-cell / disclosure policy.** No violations were found this sprint (checked directly:
   every masked cell is a genuine blank field, not a printed 0, in all reviewed TSVs), but the
   general open question about exactly which columns must always be `<20`-masked (e.g., `46`'s
   `coverage_chao2.tsv` currently leaves per-gene `s_obs` unmasked even when it's under 20, matching
   an existing convention but never explicitly ratified) is still open in `context/DECISIONS.md`.
   Worth a deliberate policy pass rather than continuing to extend precedent implicitly.

5. **Cole package README date.** The already-uploaded package's header says "release 2026-09-26"
   against its actual `release_2026-09-25/` bucket path — cosmetic only, no data changed. Decide
   whether it's worth a re-upload just to fix the string, or leave it (the fix is already in place
   for the *next* release either way).

6. **Follow-ups noted but not done this sprint:**
   - A PAF-trimmed-core comparison for genomic-level identity (aligning each observed gene span to
     the reference and trimming to the aligned core before hashing) would let the genomic level
     move from "upper bound" to a real measurement by separating UTR/boundary differences from
     true novel sequence — not done this pass, flagged as the natural next step if the genomic
     level is ever needed as more than context.
   - The mechanism behind the pooled-vs-per-ancestry completeness reversal (Simpson's-paradox
     shape) was confirmed real but not decomposed further (e.g. how much of each species' pooled
     allele count is ancestry-private) — a candidate follow-up analysis, not done here.
   - A dedicated disclosure-scanning script does not exist yet; checks this sprint were manual
     `awk`/`grep` passes. Worth building `scripts/hla_popgen/check_disclosure.py` (flags any bare
     1-19 count in a `.tsv` outside an explicit allowlist for richness/`s_obs`-style columns that
     are allowed to stay unmasked).
   - Whether platform's residual ~0.55 within-EUR AUROC signal in script 47 is worth chasing
     further, or is noise at this sample size, is unresolved.
   - The WS6 full-cohort run (7,640 people + a stronger classifier) is still pending, deprioritized
     behind KIR twice now — recommend scheduling it as its own dedicated session rather than folding
     it into a future sprint's margins again.
   - Nearest-neighbor rare-allele sharing check (a relatedness-filter sanity check suggested in
     WS-D's idea list) was not implemented.


## Addendum 2026-09-28: Marc's decisions on the open items
- `outputs/`: deleted locally, stays gitignored, history not scrubbed (Marc's call; see DECISIONS).
- Package README date: left as is (the prefix is the version of record).
- Small-cell policy: Marc wants it looser for allele/richness counts. A research agent is checking AoU policy and precedent; the outcome is recorded in DECISIONS.md.
- Stray sprint downloads in ~/Downloads (25 aggregate TSVs, a STATUS.txt, 3 diff-chunk files) are deleted. Every one was either identical to the committed copy or a superseded version already in git history. Marc's own screenshots and notes are untouched.
