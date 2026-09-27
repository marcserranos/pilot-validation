# S04 findings for Marc — KIR recurrence, KIR vs HLA catalogue, figure style, ligand co-occurrence, Cole package, cleanup

*Sprint `s04-kir-recurrence-style-share`, opened 2026-09-26 off S03. Written 2026-09-27 by the
WS-G report writer, from `LOG.md` (the primary source) and the other sprint docs — not from raw
session transcripts.*

## 1. Summary

This sprint answered three of your post-S03 asks (how often do novel KIR/HLA alleles recur across
unrelated people, does KIR's reference catalogue cover this cohort better or worse than HLA's, and
does a KIR receptor pair up with its HLA ligand more than chance by ancestry), delivered a plain-
language explainer + ranked idea list for the repertoire (WS6) work, redid the figure style pass
after you rejected the first attempt, built and uploaded the data-sharing package for Cole, and
proposed (and partly executed) a VM and repo cleanup. The single most important thing to know:
**the KIR-vs-HLA catalogue-coverage numbers went through three internal reversals** as real bugs
were found in how "novel" was being counted, and the final, only-trustworthy version needs one
more VM run (44 v4) that is not done yet — see section 2 below for what is solid now (the
recurrence/singleton-sharing work, ligand co-occurrence, naive-ML QC) and section "WS-A/WS-B
results — pending 44 v4" for what to treat as provisional. Everything genotype-derived stayed
inside the Workbench; nothing with fewer than 20 people is reported anywhere.

## 2. Per workstream

### WS-A — KIR recurrence & saturation (scripts 44/45)

**What ran and why:** you asked how often a novel KIR (and, for comparison, HLA) allele shows up
in more than one unrelated person — is KIR's ~59% "novel" rate mostly one-off sequencing/assembly
noise, or the same variant recurring because the reference catalogue is missing it systematically?
Script 44 (VM) exports, per person and gene, whether each called allele is "novel" at three levels
of strictness (genomic sequence, coding-sequence/CDS, and protein) and how many *unrelated* people
carry it, bucketed into recurrence classes (seen by 1 / 2 / >2 / >=20 unrelated people — the >=20
cutoff doubles as the disclosure floor, so that class can be described in more detail).

**Results (as of the last completed full run, "v3"):** at the genomic level (all sequence
differences), pooled over ancestries: HLA had 2,823 distinct alleles (28.5% novel), KIR had 1,460
(51.6% novel). Good-Turing incidence coverage (chance a new person's allele has already been seen)
was >=99% for both — a number that sounds like "done," but is not the same question as catalogue
completeness (see the "coverage vs completeness" figure and the Chao2 numbers under WS-B).

**Figures:** `reports/hla_popgen/44_kir_recurrence_saturation/fig_saturation_by_recurrence.png`,
`fig_saturation_per_ancestry.png`, `fig_coverage_completeness.png` (+ matching `.pdf`).

**Status: numbers superseded, method sound, final run pending — see the "WS-A/WS-B results —
pending 44 v4" section.** The genomic/protein novelty *definitions* underneath these figures were
found to be broken partway through the sprint (see Problem 2 below) and are being redone; a fixed
version ("v4") is deployed on the VM but blocked on one more bug (Problem 4). Do not quote any
genomic-level "KIR vs HLA" comparison from before that run completes.

### WS-B — KIR vs HLA catalogue coverage (script 46)

**What ran and why:** given 44's per-allele novelty/recurrence data, which reference catalogue
(IPD-KIR for KIR, IPD-IMGT/HLA for HLA) better represents what's actually circulating in this
cohort? Measured with Chao2 (a standard ecology method: estimates how many *total* distinct
alleles likely exist from how many singletons vs. doubletons you've already seen, then compares
that to what you've actually observed).

**Figures:** `reports/hla_popgen/46_kir_vs_hla_catalogue/fig_catalogue_completeness.png`,
`fig_recurrence_composition.png` (+ `.pdf`).

**Status: same as WS-A — pending 44 v4.** The early conclusion ("KIR's catalogue is better
represented than HLA's") does not survive the bug found in Problem 2 and must not be repeated
until the fixed run is in.

#### WS-A/WS-B results — pending 44 v4

The method is right; the numbers are not final. **Do not cite any v1-v3 recurrence, novelty-rate,
Chao2-completeness, or "KIR vs HLA catalogue" number from this sprint as a finding** — three
different versions of this comparison were produced and each of the first three had a real,
confirmed bug (detailed in Problem 2 below): identity at the "genomic" level was actually
name-based rather than sequence-based, which artificially collapsed distinct novel alleles into
one bucket and inflated "how well-represented" KIR looked; and KIR protein translations lacked the
artifact-filtering (partial assemblies, homopolymer-indel garbage) that HLA already had, so some
"novel proteins" were sequencing junk rather than real biology. The fix (committed, not yet run at
full scale) makes identity a true sequence hash at three matched levels for both species (genomic,
CDS, protein), applies the same artifact filter to both, and adds an automatic internal check that
raises an error if protein-level novelty ever exceeds CDS-level, which would be a physical
impossibility (protein is coarser than the sequence that encodes it) and was in fact how the bug
was caught. A smoke test of this fixed version (300 people) completed cleanly except that it
exported "genomic identity = NA" for every HLA record, traced to yet another bug (the HLA
gene-annotation file labels genes as `HLA-A` while the code expected bare `A`) — also fixed but
**not yet re-run on the full cohort**. The next and (hopefully) final action is: run 44 v4 on the
VM at full scale, then regenerate 45/46's figures and numbers from that output. Until then, treat
the "which catalogue is better represented" question as open.

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
platform is only weakly predictable (AUROC 0.573 pooled), and a follow-up run that adjusted for
ancestry showed that weak signal is essentially an ancestry echo, not a real platform batch effect
(ancestry-adjusted platform effect ~0.004, within noise). HLA does not predict KIR gene-content
class beyond chance (~0.55). 48 — across 30 receptor-ligand-pair-by-ancestry tests, no combination
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
   see the pending-results section above for what's still outstanding.

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
   bare `A`) that is what's currently blocking the final full-scale run (see pending-results
   section).

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

VM: `AoU_Jupyter_ComputeEngine_..._big_run`, `n2-highmem-4`, billed at roughly $0.37/h while
running (per `VM_RUNS.md`). The VM was started and left running across several work sessions
spanning 09-25 evening through 09-27, with jobs for scripts 44 (multiple versions), 47, and 48. Based
on the run durations recorded in `LOG.md` — 44 v1 in 282s, 44 v2 in 608s, 44 v3 in 292s, 47 v1 a
few minutes, 47 v2 in 1,837s (~31 min), 48 a few minutes, plus the 44 v4 smoke test at 132s — the
actual compute-job time sums to well under an hour, but **the VM was left running idle between
jobs for review/debugging stretches across multiple sessions**, and it auto-stops after only 1 hour
idle, so the true billed running time is higher than the job time alone and was not tracked
end-to-end in `VM_RUNS.md` (only one session row was logged with a start time, no end time). **This
is a rough estimate, not a reconciled bill: likely on the order of a few VM-hours total (roughly
$1-3 at $0.37/h), but treat the actual number as unconfirmed** — `VM_RUNS.md` should be filled in
with real start/stop timestamps if you want a firmer number.

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
   - Run 44 v4 on the full cohort (VM is deployed and blocked only on the HLA gene-prefix fix,
     which is already committed) — this is the actual final deliverable for WS-A/WS-B.
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
