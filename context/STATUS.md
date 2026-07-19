# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-07-19 — AoU-native callset validation report complete; Marc ran additional experiments on the VM not yet documented here

**Experiments A, B, C, D are all done** (see EXPERIMENTS.md). Experiment D's n=60 headline: no
single method wins outright — AoU-native is safer than our own SpecHLA at DRB1; SpecHLA is
better at DQA1. Full resolution-tier writeup already in EXPERIMENTS.md/DECISIONS.md. The
2026-07-12 repo reorg (context/scripts/reference folders) has long since synced to the VM — not
a live concern anymore.

**New thread (supervisor-requested, started 2026-07-18): validating the AoU-native callset
itself** (separate from Experiment D's cross-tool comparison) — origin, exhaustiveness, and
quality/reliability of the 535,658-person pre-computed HLA table. **This is now a complete,
ready-to-send draft: `reports/aou_callset_validation.md`.** Built from:
- `experiment_a2_callset_quality.py` — allelic diversity, external frequency concordance
  (checked against 3 published anchors, all consistent), Hardy-Weinberg proxy (DRB1 and B flagged
  as least reliable even within single ancestry groups — independent, population-genetics-based
  corroboration of Experiment D's DRB1 finding, by a completely different method).
- `experiment_a3_allele_space_coverage.py` — catalogued-vs-observed allele space (99.96% of
  observed alleles match the IMGT catalogue), and an IMGT DB-version lower bound for AoU's
  ensemble (>= release 3.60.0, decoded from `Allelelist_history.txt`'s digit-coded release
  columns — decode rule cross-validated against SpecHLA's known 3.38.0 and SpecImmune's known
  3.64.0 DB versions, both landed correctly).
- The AoU v9 "Genomics & Multi-omics Quality Report" (public, no Workbench login) — found and
  folded into Section 1: strong general-genome GIAB accuracy, but **HLA is not locus-specifically
  benchmarked anywhere in it** — this project's own evidence remains the only HLA-specific
  reliability check that exists for this dataset.

Both new scripts + the report + doc updates are committed and pushed (`origin/main`, through
commit `d2c858b`) — confirmed in sync with the VM as of the last check this session.

## Pick up here

**Open, unresolved as of this session's end: Marc ran unspecified additional pipelines/experiments
on the VM in parallel with the callset-validation work above, and this has NOT yet been described
to or reconciled with this repo's context docs.** Do not assume what they were. Next session:
ask Marc directly what was run, check `ls ~/pipeline_outputs/` on the VM for anything unfamiliar,
and check `cd ~/repos/pilot-validation && git status --short && git log --oneline -3` on the VM
itself (not just the Mac/GitHub) in case anything was committed directly from there.

Once reconciled, remaining open items from the callset-validation report (not blocking, just not
yet done): a tighter EUR-only-vs-EUR-only frequency concordance check (current version compares
pooled-cohort vs. EUR-specific published figures — directionally solid, but not apples-to-apples);
and the still-open strategic fork decision (DECISIONS.md) on how much to build directly on
AoU-native vs. call ourselves, now better-informed by this report's per-locus reliability picture.

## Watch / blockers

- **VM sleeps after ~1h idle (quirk #14)** — every new session starts with the mount gone; remount (quirk #11) and `ls`-verify before running anything that reads bucket data.
- Paste the `gcsfuse` mount and the consumer command *separately* — chaining them in one paste races (quirks #2, #14).
- **`~/ref/`, `~/repos/`, `~/tools/`, `~/pipeline_outputs/` survive a VM restart; the mount, background processes, and activated pixi shell do not** — includes `~/ref/imgt/Allelelist.txt`/`Allelelist_history.txt` (no need to re-curl after a restart).
- **Any automated SpecImmune call must use `pixi run -e specimmune`** (quirk #15) and verify its output file exists — never trust exit code 0 alone. Same discipline applies to SpecHLA (quirk #18) and to any QC metric.
- **Bash `>` redirects need the target directory to exist first** — `mkdir -p <dir>` before any `... > <dir>/file` command; Python's own `os.makedirs()` inside the script runs too late for bash's own redirect setup.
- **pad100k (SpecImmune-LR) and pad2000-10000 (SpecHLA-SR) are two separate, tool-specific recommendations** — don't conflate them.
- **Gene-panel restriction is a closed question** (Experiment C) — don't re-attempt without a specific new reason.
- **The long-read BAM pool (2,763) is a floor, not a settled ceiling** — see DECISIONS.md open question before assuming that's the whole story for a future larger cohort.
- **DRB1's SpecHLA-pad10000 padding-degradation hypothesis is still open, untested** (DECISIONS.md, 2026-07-12 entry) — re-run SpecHLA at wider padding for already-processed people if this becomes a priority.
