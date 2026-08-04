#!/usr/bin/env python3
"""Heartbeat receiver for the full-cohort Immuannot run.

Stdlib-only. Runs on the Hetzner box (isolated from the Hermes Agent project --
own port, own systemd unit, own process). Responsibilities:
  - POST /heartbeat  -- accept a JSON heartbeat from the Workbench pipeline VM,
                         append it to a flat JSONL log, update in-memory state.
  - GET  /           -- serve a single auto-refreshing HTML dashboard.
  - background watchdog thread -- fire an ntfy.sh push if no heartbeat has been
    received in STALE_MINUTES, or when a heartbeat carries anomaly=true.

Config is via environment variables (see scripts/monitoring/README.md for the
full list and deploy instructions). Nothing here is hardcoded to a secret
value -- the auth token and ntfy topic are required env vars, not defaults.
"""
import base64
import hmac
import html
import json
import os
import sys
import threading
import time
import urllib.parse
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

# ---------------------------------------------------------------------------
# Config (env-driven; fail fast on missing secrets rather than silently using
# a weak default)
# ---------------------------------------------------------------------------
PORT = int(os.environ.get("MONITOR_PORT", "8943"))
AUTH_TOKEN = os.environ.get("MONITOR_AUTH_TOKEN")
NTFY_TOPIC = os.environ.get("MONITOR_NTFY_TOPIC")
DASHBOARD_USER = os.environ.get("MONITOR_DASHBOARD_USER", "hla")
DASHBOARD_PASSWORD = os.environ.get("MONITOR_DASHBOARD_PASSWORD")
LOG_PATH = os.environ.get("MONITOR_LOG_PATH", os.path.expanduser("~/hla-monitor/heartbeats.jsonl"))
STALE_MINUTES = float(os.environ.get("MONITOR_STALE_MINUTES", "25"))
FAILURE_RATE_ANOMALY_PCT = float(os.environ.get("MONITOR_FAILURE_RATE_ANOMALY_PCT", "15"))
DISK_DANGER_PCT = float(os.environ.get("MONITOR_DISK_DANGER_PCT", "90"))
MEM_AVAIL_DANGER_PCT = float(os.environ.get("MONITOR_MEM_AVAIL_DANGER_PCT", "10"))
WATCHDOG_INTERVAL_SEC = 60
DASHBOARD_REFRESH_SEC = 30
CHART_MAX_POINTS = 500

# --- Cost model (2026-08-05, real quoted Workbench UI numbers for the n2-highcpu-96 production
# VM -- NOT the ~$3.03/hr research estimate DECISIONS.md carried). Computed SERVER-side from the
# first heartbeat of the current run session onward, rather than trusting the client's own
# cost_so_far_usd: the receiver is the thing that's actually up continuously, and this keeps the
# cost clock honest across an orchestrator restart (which resets the client's own elapsed timer).
# The client's cost_so_far_usd is still recorded in the JSONL, just not what the dashboard shows.
VM_HOURLY_USD = float(os.environ.get("MONITOR_VM_HOURLY_USD", "3.55"))
# Disk is billed continuously per month, not per running-hour -- converted at 730 h/month (the
# standard GCP convention) so it can be added to the same elapsed-hours math. At $81.60/mo this is
# ~$0.112/hr: real, but ~3% of the VM rate, so it moves the total by a few dollars over a ~2.5-day
# run rather than changing any decision.
DISK_MONTHLY_USD = float(os.environ.get("MONITOR_DISK_MONTHLY_USD", "81.60"))
HOURS_PER_MONTH = 730.0
DISK_HOURLY_USD = DISK_MONTHLY_USD / HOURS_PER_MONTH
COMBINED_HOURLY_USD = VM_HOURLY_USD + DISK_HOURLY_USD
BUDGET_USD = float(os.environ.get("MONITOR_BUDGET_USD", "300"))
# A gap longer than this starts a NEW run session for cost purposes (so an earlier smoke test's
# heartbeats don't inflate the real run's rolling cost, while a systemd restart or a brief network
# blip mid-run does NOT reset the clock). Deliberately well above both the ~5-min cadence and the
# stale threshold.
RUN_GAP_RESET_MINUTES = float(os.environ.get("MONITOR_RUN_GAP_RESET_MINUTES", "90"))

if not AUTH_TOKEN:
    sys.exit("FATAL: MONITOR_AUTH_TOKEN env var is required (never hardcode it in code or git).")
if not NTFY_TOPIC:
    sys.exit("FATAL: MONITOR_NTFY_TOPIC env var is required (topic name is a lightweight secret).")
if not DASHBOARD_PASSWORD:
    sys.exit("FATAL: MONITOR_DASHBOARD_PASSWORD env var is required -- the dashboard is reachable "
              "over the open internet (port open to anyone with the link), so it must be password-gated.")

os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)

# ---------------------------------------------------------------------------
# Shared state
# ---------------------------------------------------------------------------
_lock = threading.Lock()
_state = {
    "last_payload": None,       # most recent heartbeat dict, as received
    "last_received_ts": None,   # server-side time.time() of last heartbeat
    "history": [],              # bounded list of (received_ts, payload) for the chart
    "stale_active": False,
    "anomaly_active": False,
    "run_first_ts": None,       # server-side start of the CURRENT run session (cost clock)
    "run_finished": False,      # a heartbeat reported run_state=complete/failed -- stop the
                                # stale watchdog, so a SUCCESSFUL run doesn't page at 3am
}


def _load_history_tail():
    """Rebuild in-memory history from disk on startup (systemd restarts lose memory,
    the JSONL file is the source of truth). Also recovers the current run session's
    start timestamp, so a receiver restart mid-run does not reset the cost clock."""
    if not os.path.exists(LOG_PATH):
        return
    try:
        with open(LOG_PATH, "r") as f:
            all_ts = []
            records = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                ts = rec.get("received_ts")
                if isinstance(ts, (int, float)):
                    all_ts.append(ts)
                records.append((ts, rec.get("payload")))
        _state["history"] = records[-CHART_MAX_POINTS:]
        if _state["history"]:
            _state["last_payload"] = _state["history"][-1][1]
            _state["last_received_ts"] = _state["history"][-1][0]
        # Walk backwards to the start of the current contiguous session (see
        # RUN_GAP_RESET_MINUTES) so an earlier smoke test's heartbeats don't get
        # counted into this run's rolling cost.
        if all_ts:
            run_start = all_ts[-1]
            for earlier, later in zip(all_ts[-2::-1], all_ts[::-1]):
                if (later - earlier) / 60.0 > RUN_GAP_RESET_MINUTES:
                    break
                run_start = earlier
            _state["run_first_ts"] = run_start
    except (OSError, json.JSONDecodeError) as e:
        print(f"WARNING: could not rebuild history from {LOG_PATH}: {e}", file=sys.stderr)


def notify(message: str):
    """Fire an ntfy.sh push. Best-effort -- never let a notification failure take
    down the receiver."""
    try:
        req = urllib.request.Request(
            f"https://ntfy.sh/{NTFY_TOPIC}",
            data=message.encode("utf-8"),
            method="POST",
        )
        urllib.request.urlopen(req, timeout=10).read()
    except Exception as e:
        print(f"WARNING: ntfy push failed: {e}", file=sys.stderr)


def _resource_anomaly_reasons(payload: dict) -> list:
    """Receiver-side, point-in-time re-check of disk/mem danger -- independent of
    whatever the client computed, since these fields don't need history to validate.
    Failure-rate anomaly is trusted from the client's own 'anomaly' flag, since only
    the client has the state (previous cumulative counts) to compute a rolling rate."""
    reasons = []
    disk = payload.get("disk_used_pct")
    mem = payload.get("mem_avail_pct")
    if isinstance(disk, (int, float)) and disk > DISK_DANGER_PCT:
        reasons.append(f"disk_used_pct={disk:.1f}% > {DISK_DANGER_PCT:.0f}%")
    if isinstance(mem, (int, float)) and mem < MEM_AVAIL_DANGER_PCT:
        reasons.append(f"mem_avail_pct={mem:.1f}% < {MEM_AVAIL_DANGER_PCT:.0f}%")
    return reasons


def _handle_heartbeat(payload: dict):
    received_ts = time.time()
    reasons = _resource_anomaly_reasons(payload)
    client_anomaly = bool(payload.get("anomaly"))
    if client_anomaly and payload.get("anomaly_reason"):
        reasons.append(str(payload["anomaly_reason"]))
    is_anomaly = client_anomaly or bool(reasons)

    record = {"received_ts": received_ts, "payload": payload}
    with open(LOG_PATH, "a") as f:
        f.write(json.dumps(record) + "\n")
        f.flush()
        os.fsync(f.fileno())

    with _lock:
        # Start (or restart) the cost clock: a first-ever heartbeat, or the first one
        # after a gap long enough to mean "this is a different run, not a blip."
        prev_ts = _state["last_received_ts"]
        if (_state["run_first_ts"] is None or prev_ts is None
                or (received_ts - prev_ts) / 60.0 > RUN_GAP_RESET_MINUTES):
            _state["run_first_ts"] = received_ts
            print(f"cost clock started/reset at {received_ts} "
                  f"(gap > {RUN_GAP_RESET_MINUTES:.0f} min, or first heartbeat)", file=sys.stderr)

        _state["last_payload"] = payload
        _state["last_received_ts"] = received_ts
        _state["history"].append((received_ts, payload))
        if len(_state["history"]) > CHART_MAX_POINTS:
            _state["history"] = _state["history"][-CHART_MAX_POINTS:]

        # Terminal run states: announce once, and disarm the stale watchdog. Without this,
        # a run that FINISHES SUCCESSFULLY looks exactly like a crash to the watchdog --
        # heartbeats just stop -- and fires a 🚨 "no heartbeat" push in the middle of the night.
        run_state = str(payload.get("run_state") or "running")
        if run_state in ("complete", "failed") and not _state["run_finished"]:
            _state["run_finished"] = True
            done_n, total_n = payload.get("people_done"), payload.get("people_total")
            failed_n = payload.get("people_failed")
            if run_state == "complete":
                notify(f"✅ HLA pipeline FINISHED: {done_n}/{total_n} done, {failed_n} failed.")
            else:
                notify(f"🛑 HLA pipeline STOPPED (run_state=failed): {done_n}/{total_n} done, "
                       f"{failed_n} failed. Check the VM.")
        elif run_state not in ("complete", "failed"):
            _state["run_finished"] = False  # a new run reusing this receiver re-arms the watchdog

        # Stale state clears the moment any heartbeat arrives.
        if _state["stale_active"]:
            _state["stale_active"] = False
            notify(f"✅ HLA pipeline: heartbeat resumed after a gap.")

        # Debounce: only alert on the healthy->anomaly transition, not every
        # 15-min beat while the same issue persists.
        was_anomaly = _state["anomaly_active"]
        _state["anomaly_active"] = is_anomaly
        if is_anomaly and not was_anomaly:
            reason_text = "; ".join(reasons) if reasons else "flagged by sender"
            notify(f"⚠️ HLA pipeline anomaly: {reason_text}")
        elif was_anomaly and not is_anomaly:
            notify("✅ HLA pipeline: anomaly cleared.")


def _watchdog_loop():
    while True:
        time.sleep(WATCHDOG_INTERVAL_SEC)
        with _lock:
            last_ts = _state["last_received_ts"]
            already_stale = _state["stale_active"]
            if last_ts is None or _state["run_finished"]:
                continue
            gap_min = (time.time() - last_ts) / 60.0
            if gap_min > STALE_MINUTES and not already_stale:
                _state["stale_active"] = True
                notify(
                    f"\U0001f6a8 HLA pipeline: no heartbeat in {gap_min:.0f} min "
                    f"(threshold {STALE_MINUTES:.0f} min) -- VM may be stalled or crashed."
                )


# ---------------------------------------------------------------------------
# Dashboard rendering
# ---------------------------------------------------------------------------
def _fmt(v, suffix="", digits=1):
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.{digits}f}{suffix}"
    return f"{v}{suffix}"


def _compute_cost(run_first_ts, payload, now=None):
    """Server-side rolling cost, from the first heartbeat of the current run session
    (see RUN_GAP_RESET_MINUTES) counting upward, plus a projection to completion
    using the client's own ETA. Returns a dict of floats/None -- never raises.

    Deliberately computed here rather than trusted from the client's cost_so_far_usd:
    the receiver runs continuously, so its elapsed-time clock survives an orchestrator
    restart that would reset the client's. Both numbers land in the JSONL either way.
    """
    now = now if now is not None else time.time()
    out = {
        "elapsed_hours": None, "vm_cost": None, "disk_cost": None, "cost_so_far": None,
        "eta_hours": None, "projected_total": None, "projected_over_budget": False,
    }
    if not isinstance(run_first_ts, (int, float)):
        return out
    elapsed_h = max((now - run_first_ts) / 3600.0, 0.0)
    out["elapsed_hours"] = elapsed_h
    out["vm_cost"] = elapsed_h * VM_HOURLY_USD
    out["disk_cost"] = elapsed_h * DISK_HOURLY_USD
    out["cost_so_far"] = elapsed_h * COMBINED_HOURLY_USD

    eta = (payload or {}).get("eta_hours_remaining")
    if isinstance(eta, (int, float)) and eta >= 0:
        out["eta_hours"] = eta
        out["projected_total"] = (elapsed_h + eta) * COMBINED_HOURLY_USD
        out["projected_over_budget"] = out["projected_total"] > BUDGET_USD
    return out


def _render_chart_svg(history, key, label, stroke_var="--accent", fill=True):
    """Minimal inline SVG line chart of one numeric payload field over time. No charting
    library, no external requests (the box serves this over plain HTTP to a browser that
    may be on a phone -- keep it self-contained and tiny). Colors come from CSS custom
    properties so the chart follows the page's light/dark theme."""
    points = [(ts, p.get(key)) for ts, p in history
              if p and isinstance(p.get(key), (int, float)) and isinstance(ts, (int, float))]
    if len(points) < 2:
        return f"<p class='muted'>Not enough data yet for {html.escape(label)}.</p>"

    w, h = 720, 160
    pad_l, pad_r, pad_t, pad_b = 8, 8, 12, 18
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(0, min(ys)), max(ys)
    if y_max == y_min:
        y_max = y_min + 1
    x_span = (x_max - x_min) or 1

    def sx(x):
        return pad_l + (x - x_min) / x_span * (w - pad_l - pad_r)

    def sy(y):
        return h - pad_b - (y - y_min) / (y_max - y_min) * (h - pad_t - pad_b)

    poly = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in points)
    area = ""
    if fill:
        area = (f'<polygon points="{sx(xs[0]):.1f},{h - pad_b:.1f} {poly} '
                f'{sx(xs[-1]):.1f},{h - pad_b:.1f}" fill="var({stroke_var})" opacity="0.07" />')
    span_h = (x_max - x_min) / 3600.0
    return f"""
    <figure class="chart">
      <figcaption>{html.escape(label)} <span class="muted">· {len(points)} beats · {span_h:.1f} h</span></figcaption>
      <svg viewBox="0 0 {w} {h}" width="100%" height="{h}" preserveAspectRatio="none" role="img"
           aria-label="{html.escape(label)}">
        <line x1="{pad_l}" y1="{h - pad_b}" x2="{w - pad_r}" y2="{h - pad_b}" stroke="var(--rule)" stroke-width="1" />
        {area}
        <polyline points="{poly}" fill="none" stroke="var({stroke_var})" stroke-width="1.75"
                  stroke-linejoin="round" stroke-linecap="round" vector-effect="non-scaling-stroke" />
      </svg>
      <div class="chart-axis"><span>{_fmt(float(y_min), digits=0)}</span><span>{_fmt(float(y_max), digits=0)}</span></div>
    </figure>
    """


def _stat(label, value, sub=None, tone=None):
    tone_cls = f" tone-{tone}" if tone else ""
    sub_html = f"<div class='stat-sub'>{sub}</div>" if sub else ""
    return (f"<div class='stat{tone_cls}'><div class='stat-label'>{html.escape(label)}</div>"
            f"<div class='stat-value'>{value}</div>{sub_html}</div>")


def _render_dashboard():
    with _lock:
        payload = _state["last_payload"]
        last_ts = _state["last_received_ts"]
        history = list(_state["history"])
        stale = _state["stale_active"]
        anomaly = _state["anomaly_active"]
        run_first_ts = _state["run_first_ts"]

    if payload is None:
        return _page(
            "<div class='empty'><p>No heartbeats received yet.</p>"
            "<p class='muted'>Waiting for the pipeline to start. The cost clock begins at the "
            "first heartbeat.</p></div>", status_label="WAITING", tone="idle")

    age_min = (time.time() - last_ts) / 60.0 if last_ts else None
    run_state = str(payload.get("run_state") or "running")
    if run_state == "complete":
        status_label, tone = "COMPLETE", "ok"
    elif run_state == "failed":
        status_label, tone = "STOPPED", "bad"
    elif run_state == "starting":
        status_label, tone = "STARTING", "idle"
    elif stale:
        status_label, tone = "STALE", "bad"
    elif anomaly:
        status_label, tone = "ANOMALY", "warn"
    else:
        status_label, tone = "RUNNING", "ok"
    phase = payload.get("phase")
    if phase:
        status_label = f"{status_label} · {phase}"

    done = payload.get("people_done") or 0
    total = payload.get("people_total") or 0
    failed = payload.get("people_failed") or 0
    pct_done = (100.0 * done / total) if total else None
    cost = _compute_cost(run_first_ts, payload)

    # --- Hero stats -------------------------------------------------------
    eta_h = payload.get("eta_hours_remaining")
    eta_str = "—"
    if isinstance(eta_h, (int, float)):
        eta_str = f"{eta_h:.1f}<span class='unit'>h</span>" if eta_h < 48 else f"{eta_h / 24:.1f}<span class='unit'>d</span>"

    cost_tone = "bad" if cost["projected_over_budget"] else None
    projected = (f"proj. ${cost['projected_total']:.0f} / ${BUDGET_USD:.0f} budget"
                 if cost["projected_total"] is not None else f"budget ${BUDGET_USD:.0f}")

    stats = "".join([
        _stat("Progress", f"{done:,}<span class='unit'>/{total:,}</span>",
              f"{pct_done:.1f}% complete" if pct_done is not None else None),
        _stat("Rate", f"{_fmt(payload.get('rate_per_hour'), digits=0)}<span class='unit'>/h</span>",
              "people per hour"),
        _stat("ETA", eta_str, "remaining"),
        _stat("Cost so far",
              f"<span class='unit'>$</span>{cost['cost_so_far']:.2f}" if cost["cost_so_far"] is not None else "—",
              projected, tone=cost_tone),
    ])

    bar = ""
    if pct_done is not None:
        bar = (f"<div class='bar' role='progressbar' aria-valuenow='{pct_done:.0f}'>"
               f"<div class='bar-fill' style='width:{min(pct_done, 100):.2f}%'></div></div>")

    # --- Detail table -----------------------------------------------------
    rows = [
        ("Last heartbeat", f"{age_min:.1f} min ago" if age_min is not None else "—"),
        ("People failed", f"{failed:,}"),
        ("Rolling failure rate", _fmt(payload.get("rolling_failure_rate_pct"), "%")),
        ("Elapsed (run session)", _fmt(cost["elapsed_hours"], " h")),
        ("VM cost", f"${cost['vm_cost']:.2f}" if cost["vm_cost"] is not None else "—"),
        ("Disk cost", f"${cost['disk_cost']:.2f}" if cost["disk_cost"] is not None else "—"),
        ("Projected total", f"${cost['projected_total']:.2f}" if cost["projected_total"] is not None else "—"),
        ("Rate card", f"${VM_HOURLY_USD:.2f}/h VM + ${DISK_MONTHLY_USD:.2f}/mo disk"),
        ("Disk used", _fmt(payload.get("disk_used_pct"), "%")),
        ("Mem available", _fmt(payload.get("mem_avail_pct"), "%")),
        ("Sender timestamp", html.escape(str(payload.get("ts", "—")))),
    ]
    if payload.get("anomaly_reason"):
        rows.append(("Anomaly reason", html.escape(str(payload["anomaly_reason"]))))

    table_html = "\n".join(
        f"<tr><th scope='row'>{html.escape(k)}</th><td>{v}</td></tr>" for k, v in rows
    )

    charts = (_render_chart_svg(history, "people_done", "People completed")
              + _render_chart_svg(history, "rate_per_hour", "Throughput (people/hour)", "--accent2")
              + _render_chart_svg(history, "mem_avail_pct", "Memory available (%)", "--accent3"))

    body = f"""
    <div class="stats">{stats}</div>
    {bar}
    <section>{charts}</section>
    <section><table>{table_html}</table></section>
    """
    return _page(body, status_label=status_label, tone=tone)


def _page(body: str, status_label: str = "", tone: str = "idle") -> str:
    """Minimal, light-first shell. Deliberately restrained: near-white or near-black
    ground, one hairline rule weight, a single accent per chart, light type. Follows
    the viewer's system light/dark preference -- no toggle, no framework, no external
    fonts or scripts (the page must stay self-contained over plain HTTP)."""
    return f"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta http-equiv="refresh" content="{DASHBOARD_REFRESH_SEC}">
<meta name="color-scheme" content="light dark">
<title>Omni-HLA monitor</title>
<style>
  :root {{
    --bg:#ffffff; --fg:#111318; --muted:#6b7280; --rule:#e8eaed; --panel:#fafbfc;
    --accent:#1a7f4b; --accent2:#2563eb; --accent3:#7c3aed;
    --ok:#1a7f4b; --warn:#b45309; --bad:#b91c1c;
  }}
  @media (prefers-color-scheme: dark) {{
    :root {{
      --bg:#0c0d10; --fg:#e8eaed; --muted:#8b8f98; --rule:#1e2126; --panel:#131519;
      --accent:#4ade80; --accent2:#60a5fa; --accent3:#a78bfa;
      --ok:#4ade80; --warn:#fbbf24; --bad:#f87171;
    }}
  }}
  * {{ box-sizing:border-box; }}
  body {{
    margin:0; background:var(--bg); color:var(--fg);
    font-family:ui-sans-serif,-apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
    font-weight:300; -webkit-font-smoothing:antialiased; line-height:1.5;
    padding:32px 20px 56px;
  }}
  main {{ max-width:760px; margin:0 auto; }}
  header {{ display:flex; align-items:baseline; justify-content:space-between;
            gap:16px; flex-wrap:wrap; padding-bottom:16px; border-bottom:1px solid var(--rule); }}
  h1 {{ font-size:15px; font-weight:400; letter-spacing:.01em; margin:0; }}
  .status {{ font-size:11px; font-weight:500; letter-spacing:.09em; text-transform:uppercase;
             display:inline-flex; align-items:center; gap:7px; }}
  .status::before {{ content:""; width:7px; height:7px; border-radius:50%; background:currentColor; }}
  .tone-ok,.status.ok {{ color:var(--ok); }}
  .status.warn {{ color:var(--warn); }}
  .status.bad {{ color:var(--bad); }}
  .status.idle {{ color:var(--muted); }}
  .stats {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
            gap:1px; background:var(--rule); border:1px solid var(--rule);
            margin:28px 0 0; border-radius:3px; overflow:hidden; }}
  .stat {{ background:var(--bg); padding:16px 18px; }}
  .stat-label {{ font-size:10px; letter-spacing:.09em; text-transform:uppercase;
                 color:var(--muted); font-weight:500; }}
  .stat-value {{ font-size:27px; font-weight:250; letter-spacing:-.02em; margin-top:7px;
                 font-variant-numeric:tabular-nums; }}
  .stat-sub {{ font-size:11.5px; color:var(--muted); margin-top:3px; }}
  .stat.tone-bad .stat-value {{ color:var(--bad); }}
  .unit {{ font-size:.5em; color:var(--muted); font-weight:400; letter-spacing:0; }}
  .bar {{ height:2px; background:var(--rule); margin:20px 0 0; border-radius:2px; overflow:hidden; }}
  .bar-fill {{ height:100%; background:var(--accent); transition:width .4s ease; }}
  section {{ margin-top:36px; }}
  .chart {{ margin:0 0 26px; }}
  .chart figcaption {{ font-size:10px; letter-spacing:.09em; text-transform:uppercase;
                       color:var(--fg); font-weight:500; margin-bottom:8px; }}
  /* Explicit CSS height: with width="100%" + a numeric height attribute, browsers derive
     height from the viewBox aspect ratio instead of honouring the attribute, which made the
     charts render ~1.7x taller than intended and dominate the page. */
  .chart svg {{ display:block; overflow:visible; height:112px; width:100%; }}
  .chart-axis {{ display:flex; justify-content:space-between; font-size:10.5px;
                 color:var(--muted); font-variant-numeric:tabular-nums; margin-top:2px; }}
  table {{ border-collapse:collapse; width:100%; font-size:13px; }}
  th, td {{ padding:9px 0; border-bottom:1px solid var(--rule); text-align:left; }}
  th {{ font-weight:400; color:var(--muted); width:52%; }}
  td {{ font-weight:400; font-variant-numeric:tabular-nums; }}
  .muted {{ color:var(--muted); font-weight:300; }}
  .empty {{ padding:48px 0; }}
  footer {{ margin-top:40px; padding-top:16px; border-top:1px solid var(--rule);
            font-size:11.5px; color:var(--muted); }}
</style>
</head><body>
<main>
  <header>
    <h1>Omni-HLA &middot; full-cohort run</h1>
    <span class="status {tone}">{html.escape(status_label)}</span>
  </header>
  {body}
  <footer>
    Auto-refreshes every {DASHBOARD_REFRESH_SEC}s. Aggregate metadata only &mdash;
    no person_ids or alleles ever appear here.
  </footer>
</main>
</body></html>"""


# ---------------------------------------------------------------------------
# HTTP handling
# ---------------------------------------------------------------------------
class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        print(f"{self.address_string()} - {fmt % args}", file=sys.stderr)

    def _send(self, code, body: bytes, content_type="text/html; charset=utf-8"):
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _dashboard_auth_ok(self):
        """HTTP Basic Auth -- the browser's native login prompt blocks any page
        content from rendering until a correct password is entered, which is
        exactly the 'nothing loads until you put the password' behavior wanted
        now that the dashboard is reachable over the open internet, not just
        via SSH tunnel."""
        auth = self.headers.get("Authorization", "")
        expected = "Basic " + base64.b64encode(f"{DASHBOARD_USER}:{DASHBOARD_PASSWORD}".encode()).decode()
        return hmac.compare_digest(auth, expected)

    def _require_dashboard_auth(self):
        if self._dashboard_auth_ok():
            return True
        body = b"Authentication required."
        self.send_response(401)
        self.send_header("WWW-Authenticate", 'Basic realm="HLA pipeline monitor"')
        self.send_header("Content-Type", "text/plain")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)
        return False

    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        if parsed.path == "/":
            if not self._require_dashboard_auth():
                return
            self._send(200, _render_dashboard().encode("utf-8"))
        elif parsed.path == "/health":
            # Deliberately unauthenticated -- returns only a liveness bool, no
            # pipeline data, useful for a quick uptime check without a password.
            self._send(200, b'{"ok": true}', "application/json")
        else:
            self._send(404, b"not found")

    def do_POST(self):
        parsed = urllib.parse.urlparse(self.path)
        if parsed.path != "/heartbeat":
            self._send(404, b"not found")
            return

        token = self.headers.get("X-Auth-Token", "")
        if not hmac.compare_digest(token, AUTH_TOKEN):
            self._send(401, b'{"error": "bad or missing X-Auth-Token"}', "application/json")
            return

        length = int(self.headers.get("Content-Length", 0))
        raw = self.rfile.read(length) if length else b""
        try:
            payload = json.loads(raw.decode("utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError):
            self._send(400, b'{"error": "invalid JSON"}', "application/json")
            return

        _handle_heartbeat(payload)
        self._send(200, b'{"ok": true}', "application/json")


def main():
    _load_history_tail()
    threading.Thread(target=_watchdog_loop, daemon=True).start()
    server = ThreadingHTTPServer(("0.0.0.0", PORT), Handler)
    print(f"heartbeat_receiver listening on 0.0.0.0:{PORT}, log={LOG_PATH}", file=sys.stderr)
    server.serve_forever()


if __name__ == "__main__":
    main()
