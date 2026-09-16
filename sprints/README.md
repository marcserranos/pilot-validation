# sprints/ — autonomous work sprints, with layered context

Each sprint is one folder. Read only as deep as your task needs:

| Level | File | Who reads it | Content |
|---|---|---|---|
| L0 | `context/STATUS.md` | everyone, every session | one paragraph + pointer to the active sprint |
| L1 | `sprints/<sprint>/SPRINT.md` | the orchestrator, Marc | goals, reflection, task board, headline findings (distilled, kept current) |
| L1 | `sprints/<sprint>/JOURNAL.md` | Marc, when tracing what happened | append-only timeline: what was done, what changed, what went wrong |
| L2 | `sprints/<sprint>/WS*_*.md` | the agent working that workstream | full method, decisions, quirks, numbers, figure explanations |
| L3 | `scripts/hla_popgen/NN_*.py`, `reports/hla_popgen/NN_*/README.md` | implementers | code, tests, and each result folder's own explainer |

Rules:
- Findings graduate upward in distilled form: an L2 brief ends with a "Distilled" section of ≤10
  lines, and only that section is copied into SPRINT.md.
- Durable facts still graduate into `context/` per the root README (quirks → ENVIRONMENT.md,
  calls → DECISIONS.md, results → EXPERIMENTS.md pointer entries).
- Aggregate-only in git. Participant counts 1–19 are written as `<20` in anything new.
