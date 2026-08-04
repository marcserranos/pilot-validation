# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-08-04 (end of session) — all pre-launch decisions locked; only two build tasks remain

Everything needed to decide *how* to run the full-cohort Immuannot production job is now decided
and recorded in `context/DECISIONS.md`. Nothing here is provisional anymore:

- **VM: `n2-highcpu-96`, region `us-central1`.**
- **Concurrency: 24 people at once, 4 threads each (96 total cores).**
- **Disk: 2TB, no per-person pruning.**
- **Cost/time estimate: ~$160-200, ~52-58 hours wall-clock** (extrapolated from real 32-core data — see caveat below).
- **Low-confidence Immuannot calls: keep all of them, with confidence signals as metadata — never filtered.**
- **Monitoring: already built and deployed** (Hetzner, `hla-monitor.service`, port 8943) — just needs the orchestrator to call it.

Full rationale and the real scaling data behind these numbers: `context/DECISIONS.md` ("Concurrency
vs. threads-per-person split" + the VM/disk sub-bullets under it), `context/EXPERIMENTS.md`
(2026-08-03/04 entry).

**Caveat carried forward, not yet resolved: the 96-core extrapolation is linear from real 32-core
data, never confirmed at a higher core count.** Marc's call whether to launch on this basis or
spend one more small confirmation test first — not decided as of this session's end, and not
blocking (the orchestrator can be built and smoke-tested regardless of when that call is made).

## What's actually left — two things, both spec'd, neither built

1. **The production orchestrator.** Full implementation brief:
   `scripts/production_orchestrator/BRIEF.md` — locked-in config, required safety behavior (single-
   instance lock, hard mount check, real resumability, no `du`/`df` on broad scope), and why each
   requirement exists (every one is a real incident from this session, see ENVIRONMENT.md quirks
   #23-#27). Recommended starting point: extend `scripts/scaling_probe.py`'s already-working
   patterns (mount check, lock, `ThreadPoolExecutor` concurrency) rather than starting from zero.
2. **The full cohort list.** `scripts/build_immuannot_cohort.py` exists but was built for small
   verified test batches (`-n`, defaults to 40, stops early) — needs extending to enumerate the
   real full population (~13,000-14,521, exact number TBD), and has one real open scope question:
   whether to include `sequel2`'s 991 people (AFR-enriched, missing the alignment files the current
   trim step needs) via a new whole-contig fallback, or deliberately exclude them. Detailed in the
   orchestrator brief above — resolve explicitly with Marc, don't default silently either way.

Both are one self-contained brief (`scripts/production_orchestrator/BRIEF.md`) — a fresh session
can be pointed directly at that file plus this one and have everything needed to start building.

## Watch / blockers (carried from earlier this session — still true, will matter for the orchestrator)

- **The gcsfuse mount does not survive a VM restart or a stop.** Any new script that reads the
  mount must hard-fail-fast if it's missing (quirk #26) — do not rely on a human remembering to
  check first, this failed multiple times this session, for both Marc and the assistant giving
  instructions.
- **A restarted VM proves whatever was running already exited** — there is no "VM kills a mid-run
  job" mode. Thin/empty results after a restart means diagnose "why did it produce nothing," not
  "was it cut off partway" (quirk #14 addendum).
- **Never run two instances of the same per-person orchestrator concurrently** — the shared
  per-person *intermediate* directory is not isolated between invocations even when final output
  files are, and will get silently corrupted (quirk #23).
- **Confirm the Hetzner firewall rule for port 8943 is actually in place** before relying on the
  dashboard being externally reachable — `scripts/monitoring/README.md` step 2, needs Marc's own
  Hetzner console access, not something a Workbench-side session can do.
- This session worked in a second Verily Workbench workspace (`wb-cordial-leechee-9743`,
  Stanford-pod billing) to get a bigger test VM — getting Controlled Tier data linked to that
  workspace required a real VPC-SC/app-policy troubleshooting saga (workspace-level data linking +
  a group policy scoped to that workspace). Not logged as a durable quirk (one-time
  institutional-workspace-setup issue), but if a third workspace is ever created: **link Controlled
  Tier data before creating any compute environment in it.**
- `~/ref/`, `~/repos/`, `~/tools/`, `~/pipeline_outputs/` survive a VM restart on any workspace; the
  mount, background processes, and activated pixi shell do not.
- **DRB1 is confirmed the hardest locus by six independent, converging lines of evidence** —
  unrelated to this session's scaling work, but still standing; treat any new DRB1 result
  skeptically-but-not-surprised.
- **Gene-panel restriction is a closed question** (Experiment C) — don't re-attempt without a
  specific new reason.
- Marc and Aleix both work in this repo directly and concurrently — normal, not an anomaly.
