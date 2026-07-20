# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-07-20 — new thread: long-read data-format landscape (SpecImmune input + AoU LR manifest census)

**New thread opened this session, not yet run.** Marc's question: for the ~14,521 people with
an LR manifest row, only 2,763 have a confirmed real aligned BAM (ENVIRONMENT.md quirk #13) — why,
what input formats does SpecImmune actually accept, and can the other people's data (phased
assemblies, per DECISIONS.md) be converted to something usable? Answered the conceptual part from
existing docs (no VM needed): SpecImmune's real required input is FASTQ reads, not BAM — any
read-preserving container (BAM/CRAM/FASTQ/uBAM) works, but a phased assembly FASTA is an
already-collapsed consensus with the read-level evidence gone, so it isn't a format-conversion
problem — full reasoning now in DECISIONS.md's "Assembly-based HLA typing" bullet. Wrote
`scripts/experiment_e_lr_data_census.py` (syntax-checked) to answer the exhaustive part — full
per-person distribution of which file types actually exist, auto-detecting path-like columns
across the *entire* LR manifest rather than just the 2 known BAM columns, same
verify-by-real-existence discipline as quirk #13. **Not yet run — needs the VM's gcsfuse mount,
which this session (local Mac, no bucket access) can't reach.**

## Pick up here

**Three open items, in priority order:**

1. **Run `scripts/experiment_e_lr_data_census.py` on the VM** (new, this session). Suggested
   first step: `--sample-only 50` to sanity-check the auto-detected path-like-column list against
   the printed column report before committing to the full ~15k-row existence-check pass (a few
   minutes per column-thousand-rows). Paste the per-person profile distribution back — this is
   the direct answer to "what's the full distribution of long-read data formats."
2. **What were the "new experiments" run in parallel on the VM** (mentioned mid-session, never
   described, carried over from 2026-07-19)? Not reconciled with any context doc. Check on the VM:
   `cd ~/repos/pilot-validation && git status --short && git log --oneline -3`, and
   `ls ~/pipeline_outputs/` for anything unfamiliar, then this can get folded in properly.
3. **`scripts/hla_disease_sanity_check.py` (Celiac/Narcolepsy positive-control check) references
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
informed by both the cross-tool pilot and the callset validation report — plus, once item 1
lands, a clearer picture of how much the long-read validation cohort could grow.

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
