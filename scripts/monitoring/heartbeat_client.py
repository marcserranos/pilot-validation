#!/usr/bin/env python3
"""Heartbeat sender for the full-cohort Immuannot run.

Stdlib-only, runs on the Workbench pipeline VM (or anywhere with outbound
internet -- confirmed working from that VM in context/ENVIRONMENT.md).

Two ways to use it:
  1. Import send_heartbeat() from the not-yet-built production orchestrator --
     the clean interface boundary described in scripts/monitoring/BRIEF.md.
  2. Run it directly from the CLI for testing (a single synthetic beat, or a
     full --simulate run that proves the receiver/dashboard/alerting work
     end-to-end before the orchestrator exists).

The only per-call inputs the orchestrator needs to provide are the raw
counters (people_done, people_failed, people_total, disk_used_pct,
mem_avail_pct) -- this module derives rate/ETA/cost/rolling-failure-rate/
anomaly from those plus a small local state file, so the orchestrator itself
stays dumb and simple.
"""
import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.request

STATE_PATH_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "client_state.json")
FAILURE_RATE_ANOMALY_PCT = float(os.environ.get("MONITOR_FAILURE_RATE_ANOMALY_PCT", "15"))
DISK_DANGER_PCT = float(os.environ.get("MONITOR_DISK_DANGER_PCT", "90"))
MEM_AVAIL_DANGER_PCT = float(os.environ.get("MONITOR_MEM_AVAIL_DANGER_PCT", "10"))


def _load_state(state_path):
    if os.path.exists(state_path):
        with open(state_path) as f:
            return json.load(f)
    return None


def _save_state(state_path, state):
    tmp = state_path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(state, f)
    os.replace(tmp, state_path)


def build_payload(
    people_done,
    people_failed,
    people_total,
    disk_used_pct=None,
    mem_avail_pct=None,
    *,
    start_ts=None,
    vm_hourly_rate_usd=None,
    budget_usd=None,
    state_path=STATE_PATH_DEFAULT,
    now=None,
    run_state="running",
    phase=None,
):
    """Pure-ish function (does read/write a small local state file for the
    rolling-window failure rate) that turns raw counters into the full
    heartbeat schema from BRIEF.md.

    `run_state` (2026-08-05) is one of "starting" / "running" / "complete" / "failed".
    Without it, a run that finishes successfully is INDISTINGUISHABLE from one that
    crashed: heartbeats simply stop in both cases, and the receiver's stale watchdog
    fires the same "no heartbeat in N min" alarm either way -- at 3am, on a run that
    actually succeeded. The receiver uses this to push a completion notice instead and
    to stop arming the stale alarm."""
    now = now if now is not None else time.time()
    prev = _load_state(state_path) or {}
    first_ts = start_ts if start_ts is not None else prev.get("first_ts", now)

    elapsed_hours = max((now - first_ts) / 3600.0, 1e-9)
    rate_per_hour = people_done / elapsed_hours if elapsed_hours > 0 else 0.0
    remaining = max(people_total - people_done, 0)
    eta_hours_remaining = (remaining / rate_per_hour) if rate_per_hour > 0 else None

    cost_so_far_usd = elapsed_hours * vm_hourly_rate_usd if vm_hourly_rate_usd is not None else None

    # Rolling failure rate = failures among people completed since the last
    # heartbeat (a time-windowed approximation of "last ~50 people" -- see
    # scripts/monitoring/README.md for why, given the orchestrator only
    # reports cumulative counts, not a full per-person outcome log).
    prev_done = prev.get("people_done", 0)
    prev_failed = prev.get("people_failed", 0)
    delta_done = max(people_done - prev_done, 0)
    delta_failed = max(people_failed - prev_failed, 0)
    if delta_done > 0:
        rolling_failure_rate_pct = 100.0 * delta_failed / delta_done
    elif people_done > 0:
        rolling_failure_rate_pct = 100.0 * people_failed / people_done
    else:
        rolling_failure_rate_pct = 0.0

    reasons = []
    if rolling_failure_rate_pct > FAILURE_RATE_ANOMALY_PCT:
        reasons.append(f"rolling_failure_rate_pct={rolling_failure_rate_pct:.1f}% > {FAILURE_RATE_ANOMALY_PCT:.0f}%")
    if disk_used_pct is not None and disk_used_pct > DISK_DANGER_PCT:
        reasons.append(f"disk_used_pct={disk_used_pct:.1f}% > {DISK_DANGER_PCT:.0f}%")
    if mem_avail_pct is not None and mem_avail_pct < MEM_AVAIL_DANGER_PCT:
        reasons.append(f"mem_avail_pct={mem_avail_pct:.1f}% < {MEM_AVAIL_DANGER_PCT:.0f}%")
    anomaly = bool(reasons)

    payload = {
        "ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(now)),
        "run_state": run_state,
        "phase": phase,
        "people_total": people_total,
        "people_done": people_done,
        "people_failed": people_failed,
        "elapsed_hours": round(elapsed_hours, 3),
        "rate_per_hour": round(rate_per_hour, 3),
        "eta_hours_remaining": round(eta_hours_remaining, 3) if eta_hours_remaining is not None else None,
        "cost_so_far_usd": round(cost_so_far_usd, 2) if cost_so_far_usd is not None else None,
        "budget_usd": budget_usd,
        "disk_used_pct": disk_used_pct,
        "mem_avail_pct": mem_avail_pct,
        "rolling_failure_rate_pct": round(rolling_failure_rate_pct, 2),
        "anomaly": anomaly,
        "anomaly_reason": "; ".join(reasons) if reasons else None,
    }

    _save_state(state_path, {
        "first_ts": first_ts,
        "people_done": people_done,
        "people_failed": people_failed,
    })
    return payload


def post_heartbeat(payload, receiver_url, auth_token, timeout=10):
    req = urllib.request.Request(
        receiver_url.rstrip("/") + "/heartbeat",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json", "X-Auth-Token": auth_token},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.status, resp.read()


def send_heartbeat(
    people_done,
    people_failed,
    people_total,
    current_disk_pct=None,
    current_mem_avail_pct=None,
    *,
    receiver_url=None,
    auth_token=None,
    start_ts=None,
    vm_hourly_rate_usd=None,
    budget_usd=None,
    state_path=STATE_PATH_DEFAULT,
    run_state="running",
    phase=None,
):
    """The function the orchestrator calls every ~2 minutes (cadence revised
    2026-08-05, see monitoring/README.md "Cadence"). Reads MONITOR_URL /
    MONITOR_AUTH_TOKEN from the environment if not passed explicitly."""
    receiver_url = receiver_url or os.environ.get("MONITOR_URL")
    auth_token = auth_token or os.environ.get("MONITOR_AUTH_TOKEN")
    if not receiver_url or not auth_token:
        raise RuntimeError("receiver_url/MONITOR_URL and auth_token/MONITOR_AUTH_TOKEN are required")

    payload = build_payload(
        people_done, people_failed, people_total,
        current_disk_pct, current_mem_avail_pct,
        start_ts=start_ts, vm_hourly_rate_usd=vm_hourly_rate_usd,
        budget_usd=budget_usd, state_path=state_path,
        run_state=run_state, phase=phase,
    )
    try:
        status, body = post_heartbeat(payload, receiver_url, auth_token)
        return status, payload
    except urllib.error.URLError as e:
        # Never let a monitoring failure take down the pipeline it's watching.
        print(f"WARNING: heartbeat POST failed: {e}", file=sys.stderr)
        return None, payload


def _cli():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--url", default=os.environ.get("MONITOR_URL"), help="e.g. http://46.225.123.54:8943")
    ap.add_argument("--token", default=os.environ.get("MONITOR_AUTH_TOKEN"))
    ap.add_argument("--people-done", type=int)
    ap.add_argument("--people-failed", type=int, default=0)
    ap.add_argument("--people-total", type=int, default=14521)
    ap.add_argument("--disk-pct", type=float, default=None)
    ap.add_argument("--mem-avail-pct", type=float, default=None)
    ap.add_argument("--vm-rate", type=float, default=None, help="VM cost, USD/hour, for cost_so_far_usd")
    ap.add_argument("--budget", type=float, default=300.0)
    ap.add_argument("--state-file", default=STATE_PATH_DEFAULT)
    ap.add_argument(
        "--simulate", type=int, default=None,
        help="Ignore --people-done; send N synthetic heartbeats a few seconds apart, "
             "ramping people_done from 0 to people-total, to prove the whole system "
             "works end-to-end before the real orchestrator exists.",
    )
    ap.add_argument("--simulate-interval-sec", type=float, default=5.0)
    ap.add_argument("--simulate-fail-at", type=int, default=None,
                     help="If set, injects a burst of failures around this step to test anomaly alerting.")
    args = ap.parse_args()

    if not args.url or not args.token:
        sys.exit("--url/MONITOR_URL and --token/MONITOR_AUTH_TOKEN are required")

    if args.simulate:
        # Fresh state for a clean simulated run.
        if os.path.exists(args.state_file):
            os.remove(args.state_file)
        start_ts = time.time() - 3600  # pretend we're 1h into the run already, for a nonzero rate
        for step in range(1, args.simulate + 1):
            done = int(args.people_total * step / args.simulate)
            failed = 0
            if args.simulate_fail_at and step >= args.simulate_fail_at:
                failed = int(done * 0.3)  # inject a real failure spike to test anomaly alerting
            status, payload = send_heartbeat(
                done, failed, args.people_total,
                args.disk_pct, args.mem_avail_pct,
                receiver_url=args.url, auth_token=args.token,
                start_ts=start_ts, vm_hourly_rate_usd=args.vm_rate,
                budget_usd=args.budget, state_path=args.state_file,
            )
            print(f"[{step}/{args.simulate}] status={status} done={done} failed={failed} "
                  f"anomaly={payload['anomaly']}")
            if step < args.simulate:
                time.sleep(args.simulate_interval_sec)
        return

    if args.people_done is None:
        sys.exit("--people-done is required unless using --simulate")

    status, payload = send_heartbeat(
        args.people_done, args.people_failed, args.people_total,
        args.disk_pct, args.mem_avail_pct,
        receiver_url=args.url, auth_token=args.token,
        vm_hourly_rate_usd=args.vm_rate, budget_usd=args.budget,
        state_path=args.state_file,
    )
    print(f"status={status}")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    _cli()
