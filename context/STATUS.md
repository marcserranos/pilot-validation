# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-07-19 — AoU-native callset validation delivered; repo housekeeping pass done

**Experiments A, B, C, D are done** (EXPERIMENTS.md). **The AoU-native callset validation
sub-thread is also done and delivered** — `reports/aou_callset_validation.md` is complete,
sent to supervisors, with 3 embedded figures (`reports/figures/`). Headline: the callset passes
every reliability check run; DRB1 is independently flagged as the weakest locus by two methods
that share no assumptions (Experiment D's cross-tool pilot, and this thread's Hardy-Weinberg
proxy). Full detail: EXPERIMENTS.md's 2026-07-18/19 entry (pointer only) → the report itself.

**Housekeeping pass (2026-07-19):** `.gitignore` now excludes `__pycache__/`, `.venv*/`,
`outputs/` (was causing repeated untracked-file noise in `git status`). `README.md`'s layout
tree now lists `reports/` as a fourth tier (consult on demand, not read at session start).
Two durable quirks graduated into ENVIRONMENT.md (#20 bash redirect needing the target dir to
pre-exist, #21 IMGT's digit-coded release-history columns). DECISIONS.md's "AoU-native callset:
trust" thread got one new dated sub-bullet tying this report back into that open question.

## Pick up here

**Two open items from this session, unresolved, need your input before they can be closed out:**

1. **What were the "new experiments" run in parallel on the VM** (mentioned mid-session, never
   described)? Not reconciled with any context doc. Check on the VM: `cd ~/repos/pilot-validation
   && git status --short && git log --oneline -3`, and `ls ~/pipeline_outputs/` for anything
   unfamiliar, then this can get folded in properly.
2. **`scripts/hla_disease_sanity_check.py` (Celiac/Narcolepsy positive-control check) references
   a companion file that does not exist.** Its own `--help` text points to
   `scripts/build_phenotype_csv_template.sql`; an earlier, differently-named draft
   (`build_phenotype_csv.md`, with manual Athena/BigQuery lookup instructions) was seen once
   this session but vanished from disk before it could be committed — cause unknown, not deleted
   by this assistant. Neither file exists now, in git history or on disk. The script itself is
   real, syntax-checked, and committed, but **not runnable end-to-end without one of these** —
   flagging rather than reconstructing from a partial memory of the vanished draft. Worth a
   decision: rebuild the companion doc, or is this thread paused for now.

Not blocking otherwise. Next natural step once the above is resolved: the strategic fork in
DECISIONS.md (build downstream directly on AoU-native vs. call independently), now better
informed by both the cross-tool pilot and the callset validation report.

## Watch / blockers

- **VM sleeps after ~1h idle (quirk #14)** — every new session starts with the mount gone; remount (quirk #11) and `ls`-verify before running anything that reads bucket data.
- Paste the `gcsfuse` mount and the consumer command *separately* — chaining them in one paste races (quirks #2, #14).
- **`~/ref/`, `~/repos/`, `~/tools/`, `~/pipeline_outputs/` survive a VM restart; the mount, background processes, and activated pixi shell do not** — includes `~/ref/imgt/` (no need to re-curl the IMGT reference files after a restart).
- **Any automated SpecImmune call must use `pixi run -e specimmune`** (quirk #15/#17) and verify its output file exists — never trust exit code 0 alone. Same discipline applies to SpecHLA (quirk #18) and to any QC metric.
- **`mkdir -p` the target directory before any `... > <dir>/file` redirect** (quirk #20).
- **pad100k (SpecImmune-LR) and pad2000-10000 (SpecHLA-SR) are two separate, tool-specific recommendations** — don't conflate them.
- **Gene-panel restriction is a closed question** (Experiment C) — don't re-attempt without a specific new reason.
- **The long-read BAM pool (2,763) is a floor, not a settled ceiling** — see DECISIONS.md open question before assuming that's the whole story for a future larger cohort.
- **DRB1's SpecHLA-pad10000 padding-degradation hypothesis is still open, untested** (DECISIONS.md, 2026-07-12 entry).
