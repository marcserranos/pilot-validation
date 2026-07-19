# Omni-HLA Genomics — pilot-validation

Stage-A validation of direct HLA allele-calling (**SpecHLA** short-read, **SpecImmune** long-read) against All of Us biobank data, feeding a larger HLA allele-frequency + non-linear PRS project. This repo holds **docs + pipeline code only** — no participant data.

## Repository layout

```
pilot-validation/
├── README.md          this file — start here, always
├── pixi.toml           environment manifest (stays at root — `pixi shell`/`pixi run` need it there)
├── context/            MUST be fully read at the start of every new session — see "Read order" below
├── scripts/            all pipeline code (.py / .sh) — nothing to read up front, used as needed
├── reference/          upstream tool docs + AoU data-access research — consult, don't ingest wholesale
└── reports/            finished, delivery-ready writeups (e.g. for supervisors) — read on demand,
                        pointed to from EXPERIMENTS.md, never read wholesale at session start
```

Four tiers, by how much attention each deserves in a fresh session:
- **`context/` — read in full, every time.** Small, dense, and designed to be cheap to fully ingest. This is the project's working memory.
- **`reference/` — consult on demand.** SpecHLA/SpecImmune upstream READMEs, the AoU org PDF, and `AOU_DATA_ACCESS_NOTES.md`. Don't read these unprompted; grep/check them *before* researching something externally (a repeated, expensive lesson — see ENVIRONMENT.md quirk #16).
- **`scripts/` — used, not read.** Pipeline code. Read a specific script when you're about to run or modify it, not as part of session startup.
- **`reports/` — finished deliverables, consult on demand.** Full writeups of a completed sub-thread (with figures). EXPERIMENTS.md carries a short pointer entry, not the content — read the report itself only if you need the detail behind the pointer.

## How the context-update system works

`context/` is not one undifferentiated pile of notes — each file has a different *durability* and a different edit rule, and mixing them up is how a project's memory rots:

| File | Durability | Edit rule |
|---|---|---|
| `TASK_CONTEXT.md` | Permanent (the science) | Almost never — deliberate, structural changes only |
| `ENVIRONMENT.md` | Durable (ops knowledge) | **Append** a quirk/fix the moment you find one — never rewrite an existing entry |
| `DECISIONS.md` | Durable (the "why") | Move an item Open → Resolved when settled; append, don't delete history |
| `EXPERIMENTS.md` | Durable (the log) | **Append-only** — one entry per run, results and timing, never rewrite past entries |
| `STATUS.md` | **Ephemeral** (this session only) | **The only file fully rewritten each session** — compact, current-state-only |
| `SMOKE_TEST_PICKS.local.md` | Ephemeral, sensitive | Gitignored — real genotypes live only here, **never commit** |

**The one rule that makes this work:** only `STATUS.md` gets rewritten. Everything else is stable or append-only. When something discovered mid-session turns out to be durable, it *graduates out* of `STATUS.md` into whichever of the other four files actually owns that kind of fact — a quirk into `ENVIRONMENT.md`, a call into `DECISIONS.md`, a result into `EXPERIMENTS.md`. That's what keeps `STATUS.md` short and prevents a decision or a hard-won lesson from silently vanishing the next time `STATUS.md` gets overwritten. If you're ever unsure where a new fact belongs, ask: "will this still be true/relevant next session?" — if yes, it doesn't belong in `STATUS.md` alone.

## Read order for a new session / agent

1. **[context/TASK_CONTEXT.md](context/TASK_CONTEXT.md)** — the science: why this exists, the three aims, the competitive framing. *Stable — read once.*
2. **[context/ENVIRONMENT.md](context/ENVIRONMENT.md)** — how the machine works: VM/repo layout, the pipeline runbook, confirmed data locations, and the hard-won Workbench quirks. *Read before running anything.*
3. **[context/STATUS.md](context/STATUS.md)** — where we are right now + the literal next commands. *Read to pick up work.*
4. **[context/DECISIONS.md](context/DECISIONS.md)** — open questions (top) and resolved decisions with rationale (below). *Read to understand why, or before reopening a settled call.*
5. **[context/EXPERIMENTS.md](context/EXPERIMENTS.md)** — append-only log of every run, its result, and its runtime. *Read for what's been tried.*

## Data-handling rule

Repo is public on GitHub. Participant genotype calls live **only** in `context/SMOKE_TEST_PICKS.local.md` (gitignored). Bare `person_id`s appear in operational commands here out of pragmatic necessity — see the open compliance question in [context/DECISIONS.md](context/DECISIONS.md).
