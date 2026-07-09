# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-07-09

**Stage-A bake-off complete for all 3 smoke-test people** (1017156, 2522883, 1253627). Confidence-adjudicated (raw 54% concordance was an undercount — several "SpecImmune outliers" were `One_guess` tie-break artifacts, not real disagreements). Two real residual findings survive: **DQA1 = AoU-native systematic outlier 3/3** (robust — corroborated by both a short-read and a long-read method independently), and one confident SpecImmune DPA1 divergence. Full detail: EXPERIMENTS.md 2026-07-09 entries, raw genotypes in `SMOKE_TEST_PICKS.local.md`. **Comparison script rewritten** — no longer declares a verdict/winner, just reports raw calls + confidence.

**Sweep first attempt FAILED silently (2026-07-09)** — `run_aligner_pad_sweep.sh 2522883` ran to completion (all 10 configs, ~2h) but 9/10 produced zero real output: it was launched from a `spechla`-activated shell, and SpecImmune's pipeline needs `sniffles`, which only exists in the `specimmune` env — every gene failed with `sniffles: command not found`, but `main.py` still exited 0, so the script logged false successes. **Fixed** (ENVIRONMENT.md quirk #15): script now always invokes SpecImmune via `pixi run -e specimmune` and verifies the output file exists before declaring success. Garbage sweep dir + CSV rows need clearing before rerun (see Pick up here).

## Pick up here

**Rerun the sweep with the fix** — first clear the garbage from the failed attempt, pull the fix, smoke-test the env-pinning, then relaunch:
```bash
# 1. clear garbage from the failed attempt
rm -rf ~/pipeline_outputs/2522883/sweep
rm -f ~/pipeline_outputs/2522883/comparison_log.csv ~/pipeline_outputs/2522883/comparison_*.md

# 2. pull the fix
cd ~/repos/pilot-validation && git pull

# 3. cheap smoke test BEFORE committing to another multi-hour run
pixi run --manifest-path ~/repos/pilot-validation/pixi.toml -e specimmune -- which sniffles
# ^ must print a path. If it doesn't, STOP and report back -- don't relaunch yet.

# 4. relaunch, detached
mkdir -p ~/pipeline_outputs/2522883/sweep
nohup bash run_aligner_pad_sweep.sh 2522883 > ~/pipeline_outputs/2522883/sweep/nohup.out 2>&1 &
disown

# 5. confirm alive
sleep 5 && tail -20 ~/pipeline_outputs/2522883/sweep/progress.log
```
**This time, also verify after the FIRST new config finishes** (fullpad_minimap2, expect ~15-20min) that its actual result file exists before trusting the rest to run unattended:
```bash
ls ~/pipeline_outputs/2522883/sweep/fullpad_minimap2/specimmune_output/2522883/2522883.HLA.final.type.result.formatted.txt
```
If that's missing again, stop the sweep (`pkill -f run_aligner_pad_sweep`) and report back rather than letting it burn hours on more broken runs.

Once genuinely running with real output, check in later with `tail -f ~/pipeline_outputs/2522883/sweep/progress.log`.

**Once the sweep is analyzed**, the earlier two analysis steps are still open:
1. Check SpecImmune per-locus confidence on any remaining unadjudicated outlier calls (the rewritten compare script now surfaces this directly, no manual digging needed).
2. **Supervisor conversation** (DECISIONS → "AoU-native trust" + "Is SpecImmune-LR trustworthy"): surface the DQA1-AoU-bias finding; get validate-first vs. study-first steer; then plan the **ancestry-stratified** n=25/100 pick (use AoU genetic-ancestry TSV, `/revio/` BAM filter).

## Watch / blockers

- **VM sleeps after ~1h idle (quirk #14)** — every session starts with the mount gone; remount (quirk #11) and `ls`-verify before running anything that reads bucket data. Most session-start `FileNotFoundError`s are just this.
- Paste the `gcsfuse` mount and the consumer command *separately* — chaining them in one paste races (quirks #2, #14).
- **Unresolved risk on the running sweep:** unknown whether the VM's idle-sleep timer counts a backgrounded compute job as "activity." If the sweep stalled partway, `progress.log`'s last timestamp shows where — resume by rerunning just the missing configs (edit the `for` loops in `run_aligner_pad_sweep.sh` or rerun in full; the CSV log appends, doesn't duplicate-corrupt, but reruns WILL duplicate rows for configs already logged — dedupe by `run_label` + latest `timestamp` if needed).

## Parallel thread (non-blocking)

Discussing with supervisors whether AoU's existing srWGS HLA callset reframes the project (validate-first vs. ancestry-study-first). Fully mapped in DECISIONS.md → open question **"AoU-native callset: trust & sequencing."**
