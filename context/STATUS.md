# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-07-12 — Experiment D running on the VM; repo reorganized locally, not yet synced

**Experiments A, B, C are done** (see EXPERIMENTS.md for results). **Experiment D (fused ancestry-stratified 3-way comparison) is live on the VM right now** — a 60-person cohort (10 per ancestry group: AFR/AMR/EAS/EUR/MID/SAS), launched via `nohup bash scripts/run_experiment_d.sh ...`, using **SpecImmune with minimap2** (switched from bwa — no evidence bwa was more accurate, minimap2 matches SpecImmune's own README examples and is ~30% faster; DECISIONS.md). Mid-run, the VM's disk filled at ~person 10 (each person leaves ~3GB of intermediates); fixed with an auto-prune step added to `run_experiment_d.sh` after every person completes, so disk now stays flat. **Check progress:** `grep "run complete" ~/pipeline_outputs/experiment_d/nohup.out` (final tally) or `ls ~/pipeline_outputs/*/expd.done 2>/dev/null | wc -l` (people done so far).

**Long-read BAM manifest investigation (2026-07-11–12): closed, with an honest floor, not a ceiling.** While building Experiment D's cohort, discovered the `v9/wgs/long_read/manifest.tsv` has far more usable people than the original `/revio/`-path filter found — went through three rounds of "found the real rule" before landing on direct file-existence verification (`build_experiment_d_cohort.py`'s `find_existing_bam()`) as the only reliable method, after each pattern-based rule missed a real data source. **Confirmed floor: 2,763 of 14,521 people** have a real, sliceable BAM (up from 1,031) — but this is explicitly *not* treated as the ceiling; AoU's own release notes claim "14,000+ long-read genomes," in tension with only ~19% being sliceable, and the open question (where are the rest?) is logged in DECISIONS.md with concrete next checks. Full history — including two self-caught mistakes along the way — in ENVIRONMENT.md quirk #13. **Does not affect the currently-running 60-person job**, which was drawn from the original pool before this was discovered; a corrected, much larger cohort (`cohort_v3_verified.tsv`, already built on the VM) is ready for a future iteration once the ancestry-bias question is worth scaling up.

**Repo reorganized into folders (2026-07-12), on the Mac only — NOT yet synced to the VM.** New layout: `context/` (the 5-6 files a new session must fully ingest — this file plus TASK_CONTEXT/ENVIRONMENT/DECISIONS/EXPERIMENTS.md and the gitignored SMOKE_TEST_PICKS.local.md), `scripts/` (all 18 pipeline `.py`/`.sh` files), `reference/` (unchanged — upstream tool docs + AoU PDF), `pixi.toml` stays at repo root. See root `README.md` for the full explanation of the new layout and how the context-update system works. All internal script path references (`compare_hla_results.py`, `spechla_pad_helpers.py`, `analyze_experiment_d.py`) were updated to their new `scripts/` location and syntax-checked; `pixi.toml` references were deliberately left pointing at the repo root.

## Pick up here

**Do NOT `git pull` this reorg onto the VM until the live Experiment D job has fully finished.** `run_experiment_d.sh` is an actively-running background process there; its already-loaded copy references the *old* file locations for `compare_hla_results.py`/`spechla_pad_helpers.py`. Overwriting those files' locations mid-run (via `git pull`) would break every subsequent person's processing with `FileNotFoundError`, and editing the running script file itself is separately unsafe (undefined bash behavior in an active loop). Sequence once the job is confirmed done:

```bash
# 1. Confirm the job is actually finished first:
grep "run complete" ~/pipeline_outputs/experiment_d/nohup.out

# 2. Only then, sync the reorg:
cd ~/repos/pilot-validation && git pull

# 3. Sanity-check the new layout resolved cleanly:
ls context/ scripts/ reference/

# 4. Analyze the finished cohort:
python3 scripts/analyze_experiment_d.py ~/pipeline_outputs/experiment_d/cohort.tsv
```

If the job is still running, just keep checking back later — nothing needs to happen until it's done. If any commands from before the reorg are still in your terminal history referencing old paths (e.g. `python3 build_experiment_d_cohort.py`), they'll need a `scripts/` prefix after the pull.

## Watch / blockers

- **VM sleeps after ~1h idle (quirk #14)** — every new session starts with the mount gone; remount (quirk #11) and `ls`-verify before running anything that reads bucket data. (Not a concern for the currently-running job — it stays busy and won't idle-sleep.)
- Paste the `gcsfuse` mount and the consumer command *separately* — chaining them in one paste races (quirks #2, #14).
- **Any automated SpecImmune call must use `pixi run -e specimmune`** (quirk #15) and verify its output file exists — never trust exit code 0 alone. Same discipline applies to SpecHLA (quirk #18) and to any QC metric.
- **Disk fills fast without pruning** (quirk #19) — `run_experiment_d.sh` now auto-prunes after each person; `scripts/prune_pipeline_outputs.sh --dry-run` cleans up old cruft on demand.
- **pad100k (SpecImmune-LR) and pad2000-10000 (SpecHLA-SR) are two separate, tool-specific recommendations** — don't conflate them.
- **Gene-panel restriction is a closed question** (Experiment C) — don't re-attempt without a specific new reason.
- **The long-read BAM pool (2,763) is a floor, not a settled ceiling** — see DECISIONS.md open question before assuming that's the whole story for a future larger cohort.
