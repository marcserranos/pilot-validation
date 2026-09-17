# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-09-17 — Sprints S01 and S02 both complete, branch not pushed

**Most recent sprint board:** [`sprints/S02_supervisor_items/SPRINT.md`](../sprints/S02_supervisor_items/SPRINT.md),
with the plain-language writeup in that folder's `FINDINGS_FOR_MARC.md`. S01's board and findings
sit alongside it. Read those, not this file, for tasks and results. Branch:
`fig1-drafts-and-research-map` (off `main` @ `30b14ac`), **committed but never pushed** — pushing
is blocked for the agent, so Marc must do it.

S02 turned the 2026-09-17 supervisor call into scripts 29-33, all run on the full cohort. Open
items are listed in that board's §5: Figure 1 panels a/b need scripts 06 and 10 rerun (committing
the frequency tables, not just figures), the structural peptide-contact transfer function is
unfinished, and the VM is showing "Reboot is required" (quirk #41).

Context since the last STATUS (2026-09-05): scripts 10–22 ran on the real cohort and were merged to
main (`30b14ac` renumbered them); a supervisor meeting (~2026-09-14) framed the work as a paper with
a 4-panel Figure 1 and 8 action items; the plan for those items and future research lines is in
`reports/hla_popgen/NEXT_STEPS_AND_RESEARCH_MAP.md`; `23_fig1_draft_panels.py` drafted panels from
committed aggregates and showed novelty is mostly non-coding — S01 re-derives that foundation.

Open from 2026-09-05, still true: RUNBOOK.md does not document that 05/07 need the gcsfuse mount +
explicit `--sr-genotypes ~/mnt/aou-controlled/v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv`.
