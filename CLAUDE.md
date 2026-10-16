# Agent operating guide — Aleix's branch (`aleix/hla-resolve-phase1`)

Read automatically at the start of every Claude Code session on this branch. This repo is the
working memory for Aleix's part of the Omni-HLA project (All of Us biobank). Sessions on a
new machine start with no other memory, so everything that matters must be in git.

## Which branch, which folder

- **This branch, not `main`.** `main` is Marc's line (HLA calling, KIR, TRB germline, sprints
  S01–S05) and has no `aleix/` folder. This branch has not merged `main` since ~August, so its
  root `context/`, `scripts/` and `reference/` are an old snapshot of Marc's; read Marc's
  current state from `origin/main` (`git show origin/main:context/STATUS.md`), don't edit the
  copies here.
- **Active work: `aleix/RNA-seq/`** — TCR/BCR repertoire from AoU RNA-seq (TRUST4, SCEPTR
  embeddings, ancestry, disease). `aleix/` itself (HLA-Resolve, Phase 1) was handed to Marc
  on 2026-08-03 and is historical.

## Read order at session start

1. `aleix/RNA-seq/context/TASK_CONTEXT.md` (why), `ENVIRONMENT.md` (VMs, buckets, quirks),
   `STATUS.md` (now + next steps), `DECISIONS.md` (open questions on top), and the latest
   entries of `EXPERIMENTS.md`.
2. `aleix/RNA-seq/README.md` for layout; `reference/` docs only on demand.
3. `git log -20` on this branch and `git log -10 origin/main` for what changed lately.
4. Cross-workstream items: `origin/main:context/STATUS.md` and
   `origin/main:reports/hla_popgen/NEXT_STEPS_AND_RESEARCH_MAP.md` §3 (HLA × TCR/BCR joins).

## How to keep it current

Same file rules as the root README:
- `STATUS.md` — the only file rewritten each session: current state + literal next steps.
- `EXPERIMENTS.md` — append-only, one dated `## YYYY-MM-DD — ...` entry per run or update.
- `DECISIONS.md` — open questions on top; move to Resolved with date and rationale.
- `ENVIRONMENT.md` — append a quirk when found; never rewrite old ones.
- New analysis = numbered script `scripts/…` + `reports/NN_description/`.

When Aleix shares news (a call, a decision, a result, an idea): add a dated entry to the
owning file (result → EXPERIMENTS, decision → DECISIONS, ops → ENVIRONMENT), rewrite STATUS if
the picture changed, then commit with a clear message and `git push -u origin
aleix/hla-resolve-phase1`. Ask before deleting files or reorganizing folders.

## Hard rules

- **Public repo. No participant data**: no research_ids, genotypes or per-person rows in
  committed files, commit messages or reports. Counts 1–19 written as `<20`. Per-person tables
  stay on the VM or in the controlled buckets (`ENVIRONMENT.md`).
- Call notes are private: never quote or attribute what someone said in a call; state the
  decision on its own merits.
- Don't merge `main` into this branch, or push to `main`, without asking Aleix (and Marc).
- Write in plain, short sentences, numbers first, as the existing logs do. Use they/them for
  people unless their pronouns are known.
