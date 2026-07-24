# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-07-24 (end of session) — confidence-filtering thread closed out

This session's confidence-matched-truth work (AoU/SpecHLA vs. two truth sources) got one more
follow-on: the same global, mathematically-grounded thresholds applied **directly** to
SpecImmune-LR vs. Immuannot's own comparison (not via AoU/SpecHLA). Overall Field 2 concordance
rises from 77.6% to 92.0% once both callers are held to their confidence bar — the clearest
evidence yet that most of this cohort's raw cross-method disagreement is a confidence/noise
effect, not two genuinely different answers. DRB1 is wiped to zero confident-overlap pairs — a
sixth independent line of evidence it's the hardest locus. Full detail: `context/DECISIONS.md`
(Resolved, sub-bullet under "Confidence-matched truth comparison"), `context/EXPERIMENTS.md`
(2026-07-24 cont. entry), `reports/immuannot_pilot/README.md` (section 4).

**Also fixed:** the previously-flagged empty `reports/confidence_matched_truth/figures/` folder —
Marc downloaded both PNGs, but one landed as `field2_merged_confidence_matched-2.png` (browser
download-collision naming) while the README referenced the name without `-2`, so the link was
still dangling even after the download. Renamed to match; verified clean now.

## Pick up here

1. **Commit today's doc updates** — `context/DECISIONS.md`, `context/EXPERIMENTS.md`,
   `reports/immuannot_pilot/README.md`, the figure rename, and the new
   `reports/immuannot_pilot/figures/concordance_before_after_confidence.png`. This closes out the
   confidence-filtering thread cleanly; nothing else pending on it.
2. **Parallelization/cost-scaling experiment still not run on real data.**
   `scripts/run_core_scaling_experiment.sh` — sweeps 2/4/8-core configs for Immuannot, needs an
   8vCPU VM resize for the 8-core rows to be real measurements.
3. **HLA-Resolve is Aleix's own workstream** — `aleix/hla-resolve-phase1` branch. Read
   `aleix/README.md` directly for status rather than re-deriving it here.
4. Compute-cost-optimization research done but not acted on — ranked levers in
   `context/DECISIONS.md` ("Compute backend for scaling" entry). Not blocking.

**Next natural step (strategic, not urgent):** the fork in DECISIONS.md — build downstream
directly on AoU-native vs. call independently — now informed by both DQA1's cross-truth-source
AoU weakness *and* today's finding that SpecImmune/Immuannot converge strongly once both are
confident (i.e. the long-read "truth" itself is trustworthy when filtered).

## Watch / blockers

- **VM sleeps after ~1h idle** — every new session starts with the mount gone; remount and
  `ls`-verify before running anything that reads bucket data. Not needed for any of the
  confidence-filtering scripts — they only read `~/pipeline_outputs/`.
- Paste the `gcsfuse` mount and the consumer command *separately* — chaining them in one paste
  races. Same discipline for `pixi shell -e <env>` — wait for the `(omni-hla-pilot:<env>)` prompt
  prefix before pasting the next command.
- **Two different terminals, easy to mix up:** git commit/push happens on the Mac's own Terminal
  app (paths like `/Users/marcserrano/...` don't exist on the VM); `pixi`/`python3`/VM data paths
  happen in the Jupyter terminal tab (`jupyter@...` prompt).
- `~/ref/`, `~/repos/`, `~/tools/`, `~/pipeline_outputs/` survive a VM restart; the mount,
  background processes, and activated pixi shell do not.
- **Figures downloaded from the VM can silently get a browser-collision suffix** (`-2`, `-3`,
  ...) that won't match the README's reference — `ls` the figures folder and diff filenames
  against what the README actually links before considering a report "done."
- **Gene-panel restriction is a closed question** (Experiment C) — don't re-attempt without a
  specific new reason.
- **DRB1 is confirmed the hardest locus by six independent, converging lines of evidence now**
  (cross-tool disagreement, both tools' own confidence signals, confidence-filtering nearly
  eliminating it under every truth-source and direct-comparison variant tried) — treat any new
  DRB1 result skeptically-but-not-surprised; this is a well-established, real property of the
  locus, not a fluke.
- Marc and Aleix both work in this repo directly and concurrently — normal, not an anomaly.
