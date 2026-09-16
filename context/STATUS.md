# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-09-16 — autonomous Sprint S01 in progress

**Active sprint board:** [`sprints/S01_novelty_qc_coverage/SPRINT.md`](../sprints/S01_novelty_qc_coverage/SPRINT.md)
(read that, not this, for tasks and findings). Branch: `fig1-drafts-and-research-map` (off `main` @ `30b14ac`).

Context since the last STATUS (2026-09-05): scripts 10–22 ran on the real cohort and were merged to
main (`30b14ac` renumbered them); a supervisor meeting (~2026-09-14) framed the work as a paper with
a 4-panel Figure 1 and 8 action items; the plan for those items and future research lines is in
`reports/hla_popgen/NEXT_STEPS_AND_RESEARCH_MAP.md`; `23_fig1_draft_panels.py` drafted panels from
committed aggregates and showed novelty is mostly non-coding — S01 re-derives that foundation.

Open from 2026-09-05, still true: RUNBOOK.md does not document that 05/07 need the gcsfuse mount +
explicit `--sr-genotypes ~/mnt/aou-controlled/v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv`.
