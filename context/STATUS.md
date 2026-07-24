# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-07-24 — Experiment F follow-on complete: Immuannot-as-truth cascade + confidence-matched comparison

Two analyses finished this session, both on Experiment F's existing 60-person cohort/calls (no
new VM runs): (1) re-scored the AoU/SpecHLA field cascade with Immuannot substituted as truth
alongside the original SpecImmune-truth version; (2) built a single, mathematically-grounded
(sequencing-error-rate + IMGT inter-allele-spacing derived, not per-gene-tuned) confidence
filter for both truth sources, sample-size-matched by count against a shared denominator.
Headline: DRB1's apparent improvement under Immuannot-truth is likely a selection artifact, not
real; DQA1's AoU-native weakness is now cross-validated under two independent truth sources and
survives confidence-filtering. Full detail: `context/EXPERIMENTS.md` (2026-07-24 entry),
`context/DECISIONS.md` (Resolved decisions, same date).

## Pick up here

1. **Fix a real broken reference before trusting the report as done:** `reports/confidence_matched_truth/figures/` is empty — the two PNGs (`field2_merged_confidence_matched.png`, `retention_by_gene.png`) were never copied from the VM (`~/pipeline_outputs/experiment_d/analysis_confidence_matched/`) into the repo, so the committed README's figure references currently dangle. Download both from the Jupyter file browser, move into that folder, commit.
2. **Parallelization/cost-scaling experiment not yet run.** `scripts/run_core_scaling_experiment.sh` is written and untested-on-real-data (was smoke-tested with synthetic data only) — sweeps 2/4/8-core configs for Immuannot. Needs the VM resized to 8vCPU for the 8-core rows to be real measurements, not oversubscribed noise.
3. **HLA-Resolve is Aleix's own workstream**, not this thread's — `aleix/hla-resolve-phase1` branch, Phase 1 (ground-truth calibration) still in progress as of his last commit (2026-07-21). Read `aleix/README.md` directly rather than re-deriving status; don't duplicate his infra.
4. Compute-cost-optimization research done but not acted on — ranked levers in `context/DECISIONS.md` ("Compute backend for scaling" entry). Nothing here is blocking.

Not otherwise blocking. Next natural step: the strategic fork in DECISIONS.md (build downstream
directly on AoU-native vs. call independently) — now additionally informed by DQA1's
cross-truth-source-validated AoU weakness.

## Watch / blockers

- **VM sleeps after ~1h idle** — every new session starts with the mount gone; remount and `ls`-verify before running anything that reads bucket data. (Not needed for the confidence-matched-truth scripts — they only read `~/pipeline_outputs/`, not the AoU bucket.)
- Paste the `gcsfuse` mount and the consumer command *separately* — chaining them in one paste races. Same discipline for `pixi shell -e <env>` — wait for the `(omni-hla-pilot:<env>)` prompt prefix before pasting the next command.
- **Two different terminals, easy to mix up:** git commit/push happens on the Mac's own Terminal app (paths like `/Users/marcserrano/...` don't exist on the VM); `pixi`/`python3`/VM data paths happen in the Jupyter terminal tab (`jupyter@...` prompt). A command pasted into the wrong one fails loudly (`No such file or directory`) — check the prompt before pasting.
- `~/ref/`, `~/repos/`, `~/tools/`, `~/pipeline_outputs/` survive a VM restart; the mount, background processes, and activated pixi shell do not.
- **Gene-panel restriction is a closed question** (Experiment C) — don't re-attempt without a specific new reason.
- **DRB1 is now confirmed the hardest locus by many independent, converging lines of evidence** (cross-tool disagreement, both tools' own confidence signals, confidence-filtering nearly eliminating it under both truth sources) — treat any new DRB1 result skeptically-but-not-surprised; this is a well-established, real property of the locus at this point, not a fluke.
- Marc and Aleix both work in this repo directly and concurrently — normal, not an anomaly to flag.
