# Omni-HLA Genomics — pilot-validation

Stage-A validation of direct HLA allele-calling (**SpecHLA** short-read, **SpecImmune** long-read) against All of Us biobank data, feeding a larger HLA allele-frequency + non-linear PRS project. This repo holds **docs + pipeline code only** — no participant data.

## Read order for a new session / agent

1. **[TASK_CONTEXT.md](TASK_CONTEXT.md)** — the science: why this exists, the three aims, the competitive framing. *Stable — read once.*
2. **[ENVIRONMENT.md](ENVIRONMENT.md)** — how the machine works: VM/repo layout, the pipeline runbook, confirmed data locations, and the 13 hard-won Workbench quirks. *Read before running anything.*
3. **[STATUS.md](STATUS.md)** — where we are right now + the literal next commands. *Read to pick up work.*
4. **[DECISIONS.md](DECISIONS.md)** — open questions (top) and resolved decisions with rationale (below). *Read to understand why, or before reopening a settled call.*
5. **[EXPERIMENTS.md](EXPERIMENTS.md)** — append-only log of every run, its result, and its runtime. *Read for what's been tried.*

## What each file is for — and when to edit it

| File | Job | You edit it… |
|---|---|---|
| `TASK_CONTEXT.md` | Scientific premise & aims | Almost never — deliberate only |
| `ENVIRONMENT.md` | Ops knowledge, runbook, quirks | Append a quirk/fix when you discover one |
| `STATUS.md` | Live session state, next command | **Every session — rewrite it compactly** |
| `DECISIONS.md` | Open questions + decision record | When a question opens or closes |
| `EXPERIMENTS.md` | Run / result / timing log | **Append** one entry per run — never rewrite past entries |
| `SMOKE_TEST_PICKS.local.md` | Real `person_id`s + genotype results | Gitignored — **never commit** |
| `reference/` | Upstream tool READMEs, AoU org PDF, data-access research | Never — immutable reference |
| `pixi.toml`, `slice_and_fastq.sh`, `compare_hla_results.py` | Env + pipeline code | As needed (usage in ENVIRONMENT runbook) |

**The one rule:** only `STATUS.md` is rewritten. Everything else is stable or append-only. A durable fact discovered mid-session *graduates out* of STATUS into ENVIRONMENT (a quirk), DECISIONS (a call), or EXPERIMENTS (a result) — that's how STATUS stays short and nothing gets lost on rewrite.

## Data-handling rule

Repo is public on GitHub. Participant genotype calls live **only** in `SMOKE_TEST_PICKS.local.md` (gitignored). Bare `person_id`s appear in operational commands here out of pragmatic necessity — see the open compliance question in [DECISIONS.md](DECISIONS.md).
