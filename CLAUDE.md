# Agent operating guide — pilot-validation

This file is read automatically at the start of every Claude Code session in this repo.
It covers **how an agent operates here** — the human↔agent↔VM loop, when to stop and hand
off to a human, and the small set of quirks that bite immediately if skipped. It does not
duplicate `README.md` (the doc-tier map) or `context/ENVIRONMENT.md` (the full ops
reference) — read both of those next; this file just tells you the order and the rules of
engagement.

## Read this, then read in order

1. This file (rules of engagement).
2. **[README.md](README.md)** — the four-tier doc system (`context/` full-read,
   `reference/` on-demand, `scripts/` used-not-read, `reports/` on-demand).
3. **[context/TASK_CONTEXT.md](context/TASK_CONTEXT.md)**, **[context/ENVIRONMENT.md](context/ENVIRONMENT.md)**,
   **[context/STATUS.md](context/STATUS.md)**, **[context/DECISIONS.md](context/DECISIONS.md)**,
   **[context/EXPERIMENTS.md](context/EXPERIMENTS.md)** — in that order, per README.md's read order.

Do not skip straight to writing code because the request "sounds simple" — `context/STATUS.md`
tells you what's already in flight, and `context/DECISIONS.md` tells you what's already been
tried and rejected.

## The working loop

1. **Discuss intent with the human first.** Before planning experiments, confirm the
   direction, scope, and any open questions with whoever you're working with (Marc or
   Aleix) — don't infer a full research plan from a vague prompt. Check
   `context/DECISIONS.md` (open questions) and `reports/hla_popgen/NEXT_STEPS_AND_RESEARCH_MAP.md`
   if present for already-agreed directions before proposing new ones.
2. **Plan as a numbered script + report pair.** Convention:
   `scripts/hla_popgen/NN_description.py` → `reports/hla_popgen/NN_description/`. Numbers
   are assigned in build order, are never reused, and gaps are fine. Check the highest
   existing number (`scripts/hla_popgen/README.md`'s pipeline table) before assigning a new
   one.
3. **Run on the VM via Chrome computer-use**, not locally, for anything touching real AoU
   data — see "Driving the VM" below.
4. **Pull results back**, integrate figures/tables into the relevant report, and update
   `context/EXPERIMENTS.md` (append-only run log) and `context/STATUS.md` (rewritten,
   current state + literal next commands).
5. **Iterate** based on what the result shows — new question → back to step 1, informed by
   the result rather than the original plan.
6. **At the end of a cycle, distill.** Update `reports/hla_popgen/COMPREHENSIVE_REPORT.md`
   or the relevant per-experiment report so the narrative, rationale, and figures are
   traceable without reading the raw session transcript. A decision that will outlive this
   session belongs in `context/DECISIONS.md`, not just in the report prose.

## Driving the VM (Chrome computer-use)

Full detail: `context/ENVIRONMENT.md` quirks #31–34. The essentials:

- Marc hands over his own logged-in Chrome tab; you never see or need credentials.
- JupyterLab's Terminal is a `<canvas>` (xterm.js) — no DOM text. Work in the real terminal
  as normal, but start each terminal session with:
  ```bash
  exec > >(tee -a ~/.claude_session.log) 2>&1
  ```
  then read progress from a **Notebook/Console cell** via `get_page_text` (strip ANSI/OSC
  with the two regexes in quirk #31) — **never screenshot the terminal to read text.**
- **`type` and `Return` must always be two separate tool calls.** Batching them gets
  silently blocked by the permission classifier, regardless of content.
- The Workbench app embeds JupyterLab in a cross-origin iframe — navigate to the inner
  JupyterLab URL directly before using DOM/JS tools; the wrapper page returns nothing
  useful.
- Every VM/app instance is its own filesystem — no shared storage between instances in the
  same workspace. A fresh VM needs its own `git clone` (auto-clone doesn't work — quirk #3),
  its own pixi env build, its own `gcsfuse` mount.
- **Assume the VM has restarted at the start of every session** (auto-stops after ~1h idle,
  clock starts only once a running process exits — quirk #14). Remount the bucket, re-verify
  with `ls`, before trusting anything that reads AoU data.

## When to stop and hand off to a human

- **`git commit` (and sometimes `git push`) typed via computer-use gets blocked by the
  permission classifier**, inconsistently, even as a single standalone action (quirk #34).
  Don't burn more than 2–3 retries — hand the exact command to Marc or Aleix to run
  themselves; they aren't subject to the same classifier.
- **Before any destructive or shared-state VM operation** (killing a long-running process,
  force-restarting, rerunning something that writes to a path another invocation might also
  be writing to) — confirm with the human first. Concurrent orchestrator runs against the
  same per-person output paths have silently corrupted each other before (quirk #23).
- **Before merging, rebasing, or reading heavily from another collaborator's branch** — check
  `context/DECISIONS.md`/`STATUS.md` for concurrent-editing notes, and note in your own
  planning which branch/commit you read from, so the human can reconcile intent later.
- **Before publishing anything with small participant counts** (cell counts, per-group
  breakdowns) to a public report or GitHub — this repo is public; the disclosure question
  for small-cell counts is tracked as an open item in `context/DECISIONS.md`. Flag it rather
  than deciding it yourself.
- **Genuinely ambiguous research direction** — when the next step depends on a judgment call
  about research priority, not a technical detail, ask rather than picking one.

## Data-publication constraint

This repo is public on GitHub. It holds **docs + pipeline code only** — no participant data.
The one exception is `context/SMOKE_TEST_PICKS.local.md`, which is gitignored and must stay
that way. Real genotypes, controlled-tier identifiers, or anything derived from them at a
disclosive granularity do not belong in a committed file, a report, or a commit message.
When in doubt about whether a number is disclosive, treat it as withheld until the open
question in `context/DECISIONS.md` is resolved.

## Two collaborators, two branches

Marc and Aleix work in siloed Workbench app instances in the same workspace (no conflict by
construction — `context/ENVIRONMENT.md` quirk #8) and on separate git branches. Aleix's
active workstream lives under `aleix/` (its own `README.md`/`ENVIRONMENT_LOCAL.md`/`README.md`
mirroring this same doc-tier convention, scoped to his local WSL2 setup rather than the
shared Workbench VM). Before starting work that might overlap — a cohort selection, a shared
resource, a joint analysis referenced in `reports/hla_popgen/NEXT_STEPS_AND_RESEARCH_MAP.md` —
check whether the relevant resource or result already exists on the other branch rather than
redoing it.

## Report/traceability convention

- `context/EXPERIMENTS.md` — append-only, one entry per run (what, result, timing). This is
  the audit trail; never rewrite past entries.
- `context/DECISIONS.md` — open questions at top, resolved decisions with rationale below.
  Move an item down when settled; never delete history.
- `context/STATUS.md` — the only file fully rewritten each session. Compact, current-state
  only. Anything durable graduates out of it into one of the other three files above.
- `reports/hla_popgen/NN_description/` — the delivery-ready writeup for script `NN`, with
  figures. `reports/hla_popgen/COMPREHENSIVE_REPORT.md` is the standing narrative rollup
  across experiments — update it, don't fork a parallel summary doc.
