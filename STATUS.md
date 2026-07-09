# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-07-09

Full pipeline proven end-to-end on real AoU data. Bake-off state:
- **Person 1017156** — complete 3-way (AoU-native / SpecHLA-SR / SpecImmune-LR): 6/8 classical genes concordant, 2 discordances (allele detail in `SMOKE_TEST_PICKS.local.md`).
- **2522883, 1253627** — SpecImmune-LR done; **SpecHLA-SR running now** (launched this session).

## Pick up here

SpecHLA-SR for 2522883 and 1253627 is running (or just finished). When it's done:
```bash
python3 ~/repos/pilot-validation/compare_hla_results.py 2522883
python3 ~/repos/pilot-validation/compare_hla_results.py 1253627
```
Then paste both `~/pipeline_outputs/<id>/comparison.md` back to fold into `SMOKE_TEST_PICKS.local.md`, and look at **cross-person patterns across all 3 people** before concluding anything about method choice for the n=25/100 scale-up. n=1 wasn't enough; n=3 might not be either, but it's the next real signal.

If the SpecHLA runs weren't actually launched (pixi activation hiccup — quirk #2), the full command block is in ENVIRONMENT.md's runbook (steps 2 & 4); substitute each id.

## Watch / blockers

- pixi 3s-timeout activation hiccup (quirk #2) swallowed the run commands once this session — confirm the `(omni-hla-pilot:spechla)` prompt, and paste run commands *separately* from the `pixi shell` line.
- The two SpecHLA-SR runs give our **first SR timing baseline** — log both to EXPERIMENTS.md when they finish.

## Parallel thread (non-blocking)

Discussing with supervisors whether AoU's existing srWGS HLA callset reframes the project (validate-first vs. ancestry-study-first). Fully mapped in DECISIONS.md → open question **"AoU-native callset: trust & sequencing."**
