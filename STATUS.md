# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-07-09, end of session — signed off, roadmap locked for next session

**Stage-A bake-off complete for all 3 smoke-test people** (1017156, 2522883, 1253627), confidence-adjudicated (raw 54% concordance was an undercount — several "SpecImmune outliers" were `One_guess` tie-break artifacts). Two real residual findings survive: **DQA1 = AoU-native systematic outlier 3/3** (robust, cross-technology), and one confident SpecImmune DPA1 divergence. Comparison script rewritten — reports raw calls + confidence, no verdict.

**Aligner x padding sweep COMPLETE with real data.** (First attempt was invalid — wrong pixi env, silent failure, see ENVIRONMENT.md quirk #15 — fixed and rerun cleanly.) **pad100k looks like a safe new default** (3.9x smaller than current 4Mb, zero degradation across all 8 genes, real runtime savings). **minimap2's earlier DPB1 regression did NOT reproduce** in this 2nd person. Both n=1-2, directional. Full detail + a presentable visual card built from these results: EXPERIMENTS.md 2026-07-09 entries.

**A 4-experiment roadmap was designed and locked this session** (full design notes + quirks for each are in EXPERIMENTS.md's "Roadmap locked 2026-07-09" section — read that in full before starting any of them, don't redesign from scratch). Summary, in intended order:

- **Experiment A — AoU-native distribution stats.** Cheap, no new compute. Missingness/resolution-depth/coverage/homozygosity per locus on the already-known TSV, plus a bonus ancestry-TSV join for an Aim-1 preview.
- **Experiment B — SpecHLA (short-read) padding sweep.** Same idea as the LR sweep, but **do not reuse the LR padding levels** — short reads need their own measured insert-size distribution to set a real floor (~150bp reads vs. 15-20kb long reads). More intermediate levels + one sub-floor level, per Marc's explicit design ask. New QC check not needed for LR: mate-pair-dropout risk at narrow windows.
- **Experiment C — SpecImmune gene-panel restriction.** Investigate whether restricting SpecImmune's DB to the 8 classical genes (vs. the ~30 it currently types) is buildable — could be a bigger lever than padding. Previously ruled out as speculative; reconsidered given B's/the LR sweep's real payoff evidence. Read `make_db.py` directly first.
- **Experiment D — Fused ancestry-stratified 3-way comparison.** This **is** the n=25-100 scale-up (decided: not a separate step). Pilot on one ancestry group first, then scale to the rest. Candidate cohort already exists and is confirmed adequate (Marc, 2026-07-09) — no feasibility check needed. **Proceeding without waiting on the supervisor conversation** (decided 2026-07-09 — cost is low enough, ~$10/<40h, that the conversation doesn't need to gate it). Must use the confidence-weighted concordance metric already built into `compare_hla_results.py`, not a naive vote.

## Pick up here

**Start with Experiment A** (cheapest, no dependencies) while deciding how to sequence B and C. B and C are both aimed at cutting Experiment D's compute cost before its big spend, but **D is not blocked waiting indefinitely on them** — launch D's single-ancestry pilot once a reasonable attempt at B/C has been made, even if imperfect. Full per-experiment design notes, exact quirks, and rationale are in EXPERIMENTS.md — this file only summarizes; read that section fully first, it has everything a fresh session needs to start designing/executing without re-deriving anything.

The supervisor conversation (DECISIONS → "AoU-native trust") remains open and worth having whenever convenient — it's no longer a blocking prerequisite for anything above.

## Watch / blockers

- **VM sleeps after ~1h idle (quirk #14)** — every session starts with the mount gone; remount (quirk #11) and `ls`-verify before running anything that reads bucket data.
- Paste the `gcsfuse` mount and the consumer command *separately* — chaining them in one paste races (quirks #2, #14).
- **Any future automated SpecImmune call must use `pixi run -e specimmune`** (quirk #15) and verify its output file exists — never trust exit code 0 alone for this tool.
- **Experiment B explicitly must not reuse Experiment (LR sweep)'s padding levels** — see EXPERIMENTS.md roadmap, this is a repeated-mistake risk worth flagging twice.
