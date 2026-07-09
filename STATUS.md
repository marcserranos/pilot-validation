# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-07-09

Full pipeline proven end-to-end on real AoU data. **All three smoke-test people now have SpecHLA-SR + SpecImmune-LR done.** Bake-off state:
- **Person 1017156** — complete 3-way: 6/8 classical genes concordant, 2 discordances (allele detail in `SMOKE_TEST_PICKS.local.md`).
- **2522883, 1253627** — SpecHLA-SR done this session (timings logged to EXPERIMENTS.md: 21m8s / 19m40s); SpecImmune-LR done last session. **3-way comparison table not yet generated** — was blocked on the gcsfuse mount being gone after the VM slept (now remounted).

## Pick up here

Mount is warm again. Generate the two remaining comparison tables:
```bash
python3 ~/repos/pilot-validation/compare_hla_results.py 2522883
python3 ~/repos/pilot-validation/compare_hla_results.py 1253627
```
Then fold both `~/pipeline_outputs/<id>/comparison.md` into `SMOKE_TEST_PICKS.local.md`, and look at **cross-person patterns across all 3 people** before concluding anything about method choice for the n=25/100 scale-up. n=1 wasn't enough; n=3 might not be either, but it's the next real signal.

## Watch / blockers

- **VM sleeps after ~1h idle (quirk #14)** — every session starts with the mount gone; remount (quirk #11) and `ls`-verify before running anything that reads bucket data. Most session-start `FileNotFoundError`s are just this.
- Paste the `gcsfuse` mount and the consumer command *separately* — chaining them in one paste races (quirks #2, #14).

## Parallel thread (non-blocking)

Discussing with supervisors whether AoU's existing srWGS HLA callset reframes the project (validate-first vs. ancestry-study-first). Fully mapped in DECISIONS.md → open question **"AoU-native callset: trust & sequencing."**
