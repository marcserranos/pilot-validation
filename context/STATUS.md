# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-07-21 — long-read format census: real 5h data-loss incident, script rewritten + re-queued

**Long-read format thread:** Marc's question: for the ~14,521 people with an LR manifest row,
only 2,763 have a confirmed real aligned BAM (ENVIRONMENT.md quirk #13) — why, what input formats
does SpecImmune actually accept, and can the other people's data (phased assemblies) be converted
to something usable? Conceptual part answered from existing docs: SpecImmune's real required input
is FASTQ reads, not BAM — a phased assembly FASTA is an already-collapsed consensus with the
read-level evidence gone, so it isn't a format-conversion problem — see DECISIONS.md's
"Assembly-based HLA typing" bullet, now also carrying a subagent literature-research lead:
**Immuannot** (github.com/YingZhou001/Immuannot) looks like a ready-to-use tool that types HLA
directly from phased assembly contigs, reportedly already used by HPRC — not yet independently
verified or evaluated by us.

**Exhaustive per-person census — real incident, now fixed (ENVIRONMENT.md quirk #22, corrected
account).** `scripts/lr_manifest_format_census.py`'s first version ran ~5 hours on the VM (16
columns x ~15,424 rows, serial existence checks) and only ever printed its real results to
stderr — no file write existed until the very last lines of `main()`. **Corrected mechanism
(Marc, 2026-07-21): the VM does not stop mid-job (quirk #14) — the run most likely completed
successfully and printed everything, then the VM sat idle for an hour after it finished and
stopped, wiping that terminal's scrollback.** Separately, a later `--sample-only 50` invocation
(after the real result already existed) silently overwrote it, since the old script shared one
output path between test and real runs. **Net effect is the same either way: nothing was
recoverable** — no terminal capture (no tmux/screen/redirect was used) and no on-disk copy
survived. Rewritten and locally smoke-tested (fake mount + manifest, simulated crash-and-resume):
now checkpoints every 500 rows (fsync'd immediately), resumes automatically from the last
completed chunk on rerun, parallelizes checks (16-worker thread pool — I/O-bound over FUSE), and
only checks "primary" data columns (skips companion index files, ~halving the check count). A
`--sample-only` run now always uses its own separate output/checkpoint paths, so it can never again
silently clobber a real result (confirmed this is what happened once already — the census.tsv
Marc retrieved after "5 hours of runtime" turned out to be the leftover 50-row sanity-check output,
not the real result). **Queued to rerun on the VM — not yet done.**

**Both carried-over items from 2026-07-19 turned out to be already resolved, found via git log
during this session (see below) — not by direct investigation this session.**

**Important, unresolved: possible concurrent editing.** While this session was mid-edit on this
exact local repo path, a commit (`112372c`) landed that bundled this session's in-progress
uncommitted changes together with unrelated work neither this session nor its subagents made.
Content landed correctly (verified byte-for-byte), but **confirm whether another session/terminal
is/was active on this same working directory at the same time** — if so, treat concurrent edits
to the same files as a real risk going forward, not a one-off.

## Pick up here

**Three open items, in priority order:**

1. **Run `scripts/lr_manifest_format_census.py` on the VM** (new, this session; renamed from an
   earlier `experiment_e_*` name that collided with the unrelated, already-committed "Experiment E:
   Operation DQ" disease-sanity-check thread below). Suggested first step: `--sample-only 50` to
   sanity-check the auto-detected path-like-column list against the printed column report before
   committing to the full ~15k-row existence-check pass (a few minutes per column-thousand-rows).
   Paste the per-person profile distribution back.
2. **RESOLVED, not by this session — the "parallel VM experiments" mystery (carried over from
   2026-07-19) was Aleix's `aleix/` workstream** (3 commits, `e9acd8d`/`c148c8a`): evaluating
   HLA-Resolve as a 4th long-read caller, calibrated against real external ground truth (44
   GIAB/HPRC samples) rather than treating SpecImmune-LR as presumed-truth the way Experiment D
   did — directly answers a real methodological gap. Isolated under `aleix/`, doesn't touch shared
   files. Worth reading `aleix/README.md` directly for the current phase/status rather than
   re-deriving from this note.
3. **RESOLVED, not by this session — `scripts/hla_disease_sanity_check.py`'s missing companion
   file is moot; the thread finished with real results** (`reports/disease_sanity_check/README.md`
   + `cv_results.log`, Marc's own "Experiment E: Operation DQ" commits). Celiac shows real,
   modest, ancestry-uneven predictive signal (EUR AUROC 0.681 vs AFR 0.542); narcolepsy is
   chance-level (AUROC 0.592) despite the model correctly finding the true haplotype partner
   allele. Read the report directly for the full writeup — not summarized further here to avoid
   drifting from the source.

Not blocking otherwise. Next natural step: the strategic fork in DECISIONS.md (build downstream
directly on AoU-native vs. call independently), now informed by the cross-tool pilot, the callset
validation report, and soon Aleix's ground-truth-anchored comparison — plus, once item 1 lands, a
clearer picture of how much the long-read validation cohort could grow, and whether Immuannot is
worth evaluating for the assembly-only people.

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
