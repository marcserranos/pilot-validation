# Where the full-cohort production run's results live

> Deliverable #3 from `BRIEF.md`. Standing rule this whole project follows: real per-person allele
> calls are participant data -- stay on the VM, never committed to git or downloaded (same as
> `SMOKE_TEST_PICKS.local.md` / `immuannot_calls.tsv` throughout this project).

## Per-person raw output

`~/pipeline_outputs/<person_id>/immuannot_output/`, same layout `run_immuannot_person.py` already
uses for every prior run this project has done:

- `hap1.gtf.gz`, `hap2.gtf.gz` — Immuannot's real per-haplotype output. **This is what the
  orchestrator checks for real resumability** — a person only counts as done when both exist.
- `hap1.trimmed.fa`, `hap2.trimmed.fa` — the trimmed input actually fed to Immuannot (kept, not
  pruned — DECISIONS.md's "no per-person pruning" call, 2026-08-04).
- `hap{1,2}.immuannot.log` — full stdout/stderr of the immuannot.sh invocation, always written
  (ENVIRONMENT.md quirk #17/#18 discipline: never trust exit code alone).
- `hap{1,2}.self_align.paf` — only present for people processed via Tier 3 (self-align fallback,
  `--enable-self-align-fallback`), the synthetic .paf minimap2 produced against the chr6-only slice.
- `.orchestrator_attempts`, `.orchestrator_gave_up`, `.orchestrator_last_run.log` — orchestrator's
  own per-person bookkeeping (attempt count, give-up marker + reason, last subprocess output).

## Aggregate results

- `~/pipeline_outputs/immuannot_calls.tsv` — the real deliverable: one row per (person_id, gene),
  columns `immuannot_1`/`immuannot_2` (the two haplotype calls). Includes low-confidence and novel
  ("new"-tagged) calls — never filtered, per DECISIONS.md's explicit "call everyone, keep every
  call" instruction. Assembled by the orchestrator's `merge_fragments()` from each worker's
  isolated `immuannot_calls.<person_id>.tsv` fragment (required because concurrent worker processes
  can't safely share one file — see `run_immuannot_person.py`'s `--out-suffix` docstring). Refreshed
  every 500 completions and once more at the end of the run — safe to read mid-run for a partial
  snapshot, not just after full completion.
- `~/pipeline_outputs/immuannot_timing.tsv` — per-haplotype timing/size breakdown (trim_method,
  n_contigs, whole_contig_mb/padded_mb, per-stage seconds) — the scale-up cost/timing record, not
  itself participant genotype data, but still VM-local by convention (keeps one discipline, not two).
- `~/pipeline_outputs/immuannot_cohort_full.tsv` — the resolved cohort list this run was launched
  against (`person_id`, `platform`, `trim_tier`, `n_rows`), from `build_immuannot_cohort.py`.

## Progress / monitoring (aggregate-only, safe to view off-VM)

The Hetzner dashboard (`http://46.225.123.54:8943/`, password-gated) shows only counts/rates/ETA/
cost — never person_ids or alleles. This is the thing to actually watch during the unattended
52-58 hour run; the aggregate TSVs above are the thing to pull results from afterward, from inside
the Workbench.

## How to actually launch it (tmux + tee — not optional)

ENVIRONMENT.md quirk #14's standing practice for any unattended run over a few minutes: launch
inside `tmux` (or `nohup ... & disown`) **and** `tee` full stdout/stderr to a persistent log file.
Without this, a dropped SSH/browser session sends the foreground process a `SIGHUP` — killing the
whole orchestrator near-instantly — and that failure mode is indistinguishable from "the VM
restarted" by inspecting output files alone (this already happened once this project, cost a
misdiagnosed multi-hour incident). Do not launch this bare in a plain terminal:

```bash
tmux new -s immuannot_prod
cd ~/repos/pilot-validation && pixi shell -e specimmune
export MONITOR_AUTH_TOKEN=<from ~/hla-monitor/.env on the Hetzner box>
python3 scripts/production_orchestrator/run_production_orchestrator.py \
    --vm-rate <REAL Workbench-UI-confirmed USD/hour> \
    2>&1 | tee -a ~/pipeline_outputs/production_orchestrator.log
# Ctrl-b d to detach; `tmux attach -t immuannot_prod` to reattach later.
```

**This one command runs both phases automatically** (2026-08-05) — phase 1 (everyone except the
sequel2/`self_align_needed` group) then phase 2 (just that group, self-align fallback auto-enabled)
the moment phase 1 finishes, no second command needed. Watch the dashboard during phase 2
specifically (own isolated rate/ETA/cost, `<outroot>/monitor_state_phase2.json`) to decide whether
to let it finish or Ctrl-C it.

The persistent log is what actually proves how far the run got if you come back to a restarted VM
and find thin/empty results — never trust "looks idle, must be done" (quirk #14 addendum).

## Compliance note (carried from DECISIONS.md, still open)

Whether even bare `person_id`s belong in a public-repo-adjacent operational record is an open
compliance question (DECISIONS.md "Bare `person_id`s in a public repo"), not resolved by this
brief. Nothing above changes that: no `person_id` or allele value from the production run should be
pasted into this git repo (STATUS.md/EXPERIMENTS.md entries about the launch should reference
counts and timings, the same discipline already used for the 60-person Immuannot pilot).
