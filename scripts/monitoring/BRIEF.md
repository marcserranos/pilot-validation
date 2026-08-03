# Monitoring/reporting pipeline for the full-cohort Immuannot run — implementation brief

> **Status: not started.** This is a locked-in design + open-decisions brief, written 2026-08-03
> so a *separate* session/agent can build this without needing the parent conversation's history.
> Self-contained on purpose — read this file alone, nothing else required as prerequisite context.

## Why this exists

The production run (~13,000-14,500 people, Immuannot, on a large Workbench VM) is expected to cost
around $300 and run for multiple days. That combination — expensive, long, and only checkable by
Marc occasionally — is exactly the situation where "is it still working, and are the results any
good" needs to be answerable at a glance, without SSHing in and reading logs. This subsystem is
that answer. It must be genuinely minimal: a small periodic sample, not continuous telemetry, and
it must not add real load or complexity to the pipeline it's watching.

**Explicitly not in scope:** the production orchestrator itself (the script that will actually
process 13k+ people). That doesn't exist yet — it's the next deliverable after the current 32-person
core-scaling diagnostic finishes and Marc picks a VM size from the results. This monitoring
subsystem must be buildable and testable *independently* of that orchestrator (see "Interface
boundary" below) so the two workstreams don't block each other.

## Decided architecture

```
[Workbench pipeline VM]                    [Hetzner box, always-on]
  orchestrator (not yet built)
    -> periodically calls a small
       "send_heartbeat(...)" client   --POST-->  heartbeat receiver (new, isolated
       function/script (JSON, ~200B)              service -- see "Hetzner facts" below)
                                                     |
                                                     +-> appends to a flat JSONL log
                                                     +-> serves a single auto-refreshing
                                                     |   HTML dashboard (GET /)
                                                     +-> watchdog thread: if no heartbeat
                                                         received in > STALE_THRESHOLD,
                                                         fire an ntfy.sh push (this is the
                                                         real stall-detector -- it lives on
                                                         the receiver, not the sender, because
                                                         a crashed VM can't self-report)
                                                     +-> also fires ntfy.sh push immediately
                                                         if a heartbeat arrives with
                                                         anomaly=true (quality problem, not
                                                         just liveness)
```

Key decisions already made, do not re-litigate without a real reason:
- **Push (VM → Hetzner), not pull.** The Workbench VM is ephemeral and outbound-only is already
  confirmed working in this project (plain `curl` to external hosts works fine — see
  `context/ENVIRONMENT.md`); no new inbound port needed on the pipeline VM.
- **Stall detection lives on the receiver (Hetzner), via silence, not on the sender.** A genuinely
  crashed/hung VM cannot be relied on to report its own failure.
- **Cadence: heartbeat every ~15 minutes**, not per-person (at real production concurrency, many
  people finish close together — a per-person ping would be noisy, not minimal).
- **Storage: flat JSONL file, no database.** At 15-min cadence over a multi-day run this is at
  most a few hundred lines — trivially small, no need for sqlite/postgres/anything else.
- **Dashboard: single server-rendered HTML page, no JS framework, no build step.** Auto-refresh via
  `<meta http-equiv="refresh">` is sufficient; a simple inline-SVG polyline is enough for the
  completions-over-time chart — do not reach for a charting library for this.
- **Alerting channel: ntfy.sh**, not a new service. Already used in this exact codebase
  (`scripts/run_experiment_d.sh`'s `notify()` function, `NTFY_TOPIC` env var, pattern:
  `curl -s -d "$MESSAGE" "ntfy.sh/$NTFY_TOPIC"`) — reuse that pattern, don't reinvent it. That
  script's topic was never actually set (`NTFY_TOPIC=""` by default, example format in its comment:
  `"marc-hla-9f3k2x"`, an unguessable string) — **a fresh topic needs to be chosen for this**, not
  reused from anywhere, since nothing has ever actually used it yet.
- **Access: default to SSH-tunnel-only, no new firewall rule, unless Marc says otherwise.**
  (`ssh -L <port>:localhost:<port> root@46.225.123.54`, then browse `localhost:<port>`.) This is
  the zero-new-attack-surface option and matches the Hetzner box's existing `ssh-only` firewall
  (only inbound TCP 22 + ICMP allowed). Opening a real port + adding a shared-secret token is an
  explicit upgrade path if Marc wants phone-browser access without SSH — not the default.

## Hetzner box facts (embedded here so this brief is self-contained — this is a *different*
project's infrastructure, documented in a separate repo Marc may not have open)

- Server name `hermes-agent`, public IP **46.225.123.54**, Hetzner Cloud CX23, Ubuntu 26.04,
  Nuremberg. SSH key `~/.ssh/id_ed25519` on Marc's Mac.
- Firewall `ssh-only`: **only inbound TCP 22 + ICMP allowed**; everything else blocked inbound,
  outbound unrestricted. Any new listening service needs an explicit new firewall rule to be
  reachable from outside — do not assume a port is open just because the service binds to it.
- The box already runs an unrelated project (**Hermes Agent**, a personal AI assistant, Python-
  based) as systemd services, including a WhatsApp gateway on **port 3000**. **This monitoring
  service must be a completely separate, isolated systemd unit on its own port — do not touch,
  depend on, or share a process/venv with the Hermes Agent setup.** Avoid port 3000 and whatever
  else Hermes is using; check with `sudo ss -tlnp` before picking a port, don't assume it's free.
- **Marc runs all commands on this box himself, in his own SSH'd-in terminal.** Whoever builds
  this (human or agent) should produce files + exact commands for Marc to run, not attempt to
  operate the server directly — same working pattern as the rest of that project.
- Python3 is already present on this box (Hermes Agent is Python-based) — a stdlib-only receiver
  (`http.server`, `json`, `threading`, `urllib.request` for the ntfy POST) needs no new
  dependencies and is the right choice; don't add Flask/FastAPI/etc. for something this small.

## Heartbeat payload — proposed schema (open to adjustment, not sacred)

```json
{
  "ts": "2026-08-05T14:32:00Z",
  "people_total": 14521,
  "people_done": 3812,
  "people_failed": 22,
  "elapsed_hours": 18.4,
  "rate_per_hour": 207.2,
  "eta_hours_remaining": 51.8,
  "cost_so_far_usd": 21.3,
  "budget_usd": 300,
  "disk_used_pct": 12,
  "mem_avail_pct": 64,
  "rolling_failure_rate_pct": 1.5,
  "anomaly": false,
  "anomaly_reason": null
}
```
No person_ids, alleles, or any participant-level data in this payload or the dashboard, ever —
same aggregate-only discipline this whole project already follows for egress
(`context/DECISIONS.md`). This is pure operational metadata.

## Open decisions the new session must settle with Marc before/while building

1. **Anomaly threshold(s).** Proposed: rolling failure rate over the last ~50 people >15-20%
   triggers `anomaly: true`. Not confirmed. Also consider whether a `disk_used_pct` or
   `mem_avail_pct` crossing a danger line (e.g. disk >90%, mem available <10%) should also set
   `anomaly: true` — not yet decided whether resource danger and quality danger should be the same
   flag or two separate ones.
2. **Stale-heartbeat threshold** for the watchdog. Proposed: ~25 min (roughly 1.6x the 15-min
   cadence, tolerating one missed beat before alarming). Not confirmed.
3. **Port number** for the receiver on the Hetzner box. Not chosen — check what's free first
   (see Hetzner facts above), don't just default to something like 8080/8420 without checking.
4. **ntfy.sh topic name.** Needs a fresh, unguessable string chosen now (topics on ntfy.sh are
   public by anyone who knows the name — treat the topic name itself as a lightweight secret).
5. **Auth token** for the receiver's endpoints (a simple shared-header check is proposed — e.g.
   `X-Auth-Token`). Value needs generating; **never commit the actual token value to the public
   repo** — it should live only in an env var / untracked local file on both the sender and
   receiver sides.
6. **VM hourly rate for the cost_so_far_usd calculation** — this is a real external dependency,
   not a design gap: it depends on which machine type the *production* run actually uses, which
   depends on the still-running core-scaling diagnostic's results. Build the client/receiver to
   take this as a configurable parameter (e.g. an env var or CLI flag), not a hardcoded number —
   it will be filled in once Marc picks the production VM.
7. **Does Marc already have the ntfy.sh app installed / a channel he's watching?** `NTFY_TOPIC` in
   `run_experiment_d.sh` was never actually exercised (default empty, never set) — confirm this is
   genuinely a fresh setup, not something already wired up elsewhere.

## Interface boundary with the not-yet-built production orchestrator

Do not block on the orchestrator's design. The clean seam: the orchestrator's only obligation,
whenever it exists, is to call something like

```
send_heartbeat(people_done, people_failed, people_total, current_disk_pct, current_mem_avail_pct)
```

every ~15 minutes (or the orchestrator can just write its own small `progress.json` snapshot on
its own schedule, and a *separate* tiny loop/cron reads that file and POSTs it — either shape is
fine). **Build the receiver, the dashboard, and this client-side sender as a fully working,
independently-testable unit now** (fake/synthetic heartbeats are enough to prove it end-to-end).
Wiring it into the real orchestrator is a small, later step once that script exists — it should
not require redesigning anything built now.

## Expected deliverables from the new session

1. `scripts/monitoring/heartbeat_receiver.py` — stdlib-only, deployed on Hetzner as a systemd unit
   (write the unit file too), appends to JSONL, serves the dashboard, runs the watchdog thread,
   fires ntfy on anomaly/staleness.
2. `scripts/monitoring/heartbeat_client.py` (or `.sh`) — the sender, callable standalone for
   testing and importable/callable from the future orchestrator.
3. A short README in `scripts/monitoring/` covering: how to deploy on the Hetzner box (exact
   commands for Marc to run himself), how to view the dashboard (SSH tunnel command), how to test
   end-to-end with synthetic heartbeats before the real orchestrator exists, and where the open
   decisions above got settled (so this BRIEF.md's "open decisions" section can be marked resolved
   or superseded, not silently ignored).
