# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-08-04 — scaling diagnostic done with real data; production orchestrator is the one thing genuinely not built

Two things converged this session for the full-cohort Immuannot run (~13,000-14,521 people):

**1. Monitoring subsystem — done, deployed, working.** Built in a parallel session
(`scripts/monitoring/`), live on Hetzner (`hermes-agent`, 46.225.123.54:8943) as systemd unit
`hla-monitor.service`. Password-gated dashboard, `POST /heartbeat` receiver, ntfy.sh alerting on
anomaly/silence. Tested end-to-end with synthetic heartbeats. **One loose end: the Hetzner
firewall rule for port 8943** — needs Marc's own Hetzner console access, check whether it's been
added yet before assuming external reachability works. Full detail: `scripts/monitoring/README.md`.

**2. Core-scaling diagnostic — done, real data, a decision reached (not yet formally committed).**
`scripts/scaling_probe.py` (a from-scratch Python rewrite, replacing a bash orchestrator that
became a full day/night incident — see ENVIRONMENT.md quirks #23-27 before touching anything that
orchestrates multiple Immuannot invocations) produced real throughput/memory data at a fixed
32-core budget: **the concurrency/threads-per-person split doesn't affect speed (all within ~2%),
only memory** — so prefer a moderate split (concurrency ≈ cores/4, threads=4) over maxing
concurrency, same speed, ~half the memory risk. Full numbers and the extrapolation to a 96-core
production machine (~52-58h, ~$180-200): `context/EXPERIMENTS.md` (2026-08-03/04 entry),
`context/DECISIONS.md` (Resolved, "Concurrency vs. threads-per-person split").

**Genuinely not built: the production orchestrator itself.** This is the next concrete task, not
optional cleanup. It needs, at minimum (each of these is a real lesson paid for this session, not
speculative):
- The 8×4-style concurrency/threads config from the Resolved decision above.
- A single-instance lock (PID-file, refuse to start if another instance is live) — ENVIRONMENT.md quirk #23.
- A hard mount check at startup that refuses to run at all without a verified working mount — quirk #26.
- Real resumability (skip already-done people, safe to kill and relaunch) — the pattern
  `run_experiment_d.sh` already used successfully at 60-person scale.
- Per-person disk pruning as it goes — at ~81MB/person × 14,521 people ≈ **1.1-1.2 TB** if nothing
  is ever cleaned up, which would exceed a smaller data disk; the disk-size decision when creating
  the production VM needs this number, not a guess.
- Wiring into the already-deployed heartbeat monitor (`scripts/monitoring/README.md`'s "Wiring
  into the real orchestrator" section has the exact call shape) — call it every ~15 min with
  `people_done`/`people_failed`/`people_total`/`--vm-rate <the chosen machine's $/hr>`.

## Pick up here

1. **Confirm the Hetzner firewall rule for port 8943 is actually in place** (`scripts/monitoring/README.md` step 2) — quick to check, blocks knowing whether the dashboard is externally reachable.
2. **Decide: commit to `n1-highcpu-96` now, or spend one more small `scaling_probe.py` test at a higher core count first** to confirm the flat-throughput extrapolation holds past 32 cores. Marc's call — not yet made as of this session's end.
3. **Build the production orchestrator** — the bullet list above is the spec. Consider Python from the start (not bash) given this session's experience — `scaling_probe.py` is a reasonable starting skeleton, not a from-scratch design.
4. **Build the full-cohort person list.** `scripts/build_immuannot_cohort.py` exists but was built for a small verified test batch (`-n` up to a few dozen) — check whether it needs adjustting to comfortably enumerate the full ~13,000-14,521 assembly-holding population, not just verify a small sample.
5. Once the orchestrator exists and is smoke-tested small, launch the real run — cost/time budget from step 2's decision, monitored live via the already-deployed dashboard.

## Watch / blockers

- **The gcsfuse mount does not survive a VM restart or a stop — this bit this session hard, repeatedly, across many hours.** Any new script that reads the mount should hard-fail-fast if it's missing (quirk #26), not assume a human remembered to check.
- **A restarted VM proves whatever was running already exited — there is no "VM kills a mid-run job" mode.** If you come back to a restarted VM after an unattended run and find thin/empty results, diagnose it as "why did a completed/crashed run produce nothing," not "was it cut off partway" (quirk #14 addendum, 2026-08-03).
- **Never run two instances of the same per-person orchestrator concurrently** — even with output-file isolation via `--out-suffix`, the shared per-person *intermediate* directory (trimmed FASTA, `.gtf.gz`) is not isolated between invocations and will get silently corrupted (quirk #23).
- **This session worked in a second Verily Workbench workspace** (institutional/Stanford-pod billing, project `wb-cordial-leechee-9743`, distinct from the original personal workspace `wb-glacial-potato-8710`) to get a bigger test VM. Getting that workspace to actually access Controlled Tier data required a real VPC Service Controls / app-policy troubleshooting saga (workspace-level data linking + a group policy scoped to that specific workspace) — not logged as a durable quirk since it was a one-time institutional-workspace-setup issue, but worth remembering if a THIRD workspace is ever created: link Controlled Tier data before creating any compute environment in it.
- `~/ref/`, `~/repos/`, `~/tools/`, `~/pipeline_outputs/` survive a VM restart on any workspace; the mount, background processes, and activated pixi shell do not.
- **DRB1 is confirmed the hardest locus by six independent, converging lines of evidence** (cross-tool disagreement, both tools' own confidence signals, confidence-filtering nearly eliminating it under every truth-source and direct-comparison variant tried) — treat any new DRB1 result skeptically-but-not-surprised; this is a well-established, real property of the locus, not a fluke. Unrelated to this session's scaling work, but still standing.
- **Gene-panel restriction is a closed question** (Experiment C) — don't re-attempt without a specific new reason.
- Marc and Aleix both work in this repo directly and concurrently — normal, not an anomaly.
