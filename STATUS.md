# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-07-09

**Stage-A bake-off complete for all 3 smoke-test people** (1017156, 2522883, 1253627). Confidence-adjudicated (raw 54% concordance was an undercount — several "SpecImmune outliers" were `One_guess` tie-break artifacts, not real disagreements). Two real residual findings survive: **DQA1 = AoU-native systematic outlier 3/3** (robust — corroborated by both a short-read and a long-read method independently), and one confident SpecImmune DPA1 divergence. Full detail: EXPERIMENTS.md 2026-07-09 entries, raw genotypes in `SMOKE_TEST_PICKS.local.md`. **Comparison script rewritten** — no longer declares a verdict/winner, just reports raw calls + confidence.

**RUNNING NOW (unattended):** `run_aligner_pad_sweep.sh 2522883` — combined bwa-vs-minimap2 x 5-padding-level sweep for SpecImmune-LR, launched via `nohup`, ~5h expected. Full design in EXPERIMENTS.md's "combined aligner x padding sweep" entry. **If you're a fresh session reading this: check `~/pipeline_outputs/2522883/sweep/progress.log` first** to see whether it finished, is still running, or stalled (possible VM-sleep interruption, see Watch/blockers).

## Pick up here

**If the sweep is still running or just finished:** read `~/pipeline_outputs/2522883/sweep/progress.log` and `comparisons.md`, report status. If complete, analyze the 10-config matrix in `~/pipeline_outputs/2522883/comparison_log.csv` — does padding degrade accuracy before pad40k ("theoretical minimum")? Does minimap2's speedup hold with less noise at smaller windows? Fold real findings into EXPERIMENTS.md.

**If the sweep hasn't been started yet:** on the VM (after remount, see Watch/blockers), run:
```bash
cd ~/repos/pilot-validation && git pull   # picks up run_aligner_pad_sweep.sh + rewritten compare script
nohup bash run_aligner_pad_sweep.sh 2522883 > ~/pipeline_outputs/2522883/sweep/nohup.out 2>&1 &
disown
```
Then check in later with `tail -f ~/pipeline_outputs/2522883/sweep/progress.log`.

**Once the sweep is analyzed**, the earlier two analysis steps are still open:
1. Check SpecImmune per-locus confidence on any remaining unadjudicated outlier calls (the rewritten compare script now surfaces this directly, no manual digging needed).
2. **Supervisor conversation** (DECISIONS → "AoU-native trust" + "Is SpecImmune-LR trustworthy"): surface the DQA1-AoU-bias finding; get validate-first vs. study-first steer; then plan the **ancestry-stratified** n=25/100 pick (use AoU genetic-ancestry TSV, `/revio/` BAM filter).

## Watch / blockers

- **VM sleeps after ~1h idle (quirk #14)** — every session starts with the mount gone; remount (quirk #11) and `ls`-verify before running anything that reads bucket data. Most session-start `FileNotFoundError`s are just this.
- Paste the `gcsfuse` mount and the consumer command *separately* — chaining them in one paste races (quirks #2, #14).
- **Unresolved risk on the running sweep:** unknown whether the VM's idle-sleep timer counts a backgrounded compute job as "activity." If the sweep stalled partway, `progress.log`'s last timestamp shows where — resume by rerunning just the missing configs (edit the `for` loops in `run_aligner_pad_sweep.sh` or rerun in full; the CSV log appends, doesn't duplicate-corrupt, but reruns WILL duplicate rows for configs already logged — dedupe by `run_label` + latest `timestamp` if needed).

## Parallel thread (non-blocking)

Discussing with supervisors whether AoU's existing srWGS HLA callset reframes the project (validate-first vs. ancestry-study-first). Fully mapped in DECISIONS.md → open question **"AoU-native callset: trust & sequencing."**
