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
}


def _load_history_tail():
    """Rebuild in-memory history from disk on startup (systemd restarts lose memory,
    the JSONL file is the source of truth)."""
    if not os.path.exists(LOG_PATH):
        return
    try:
        with open(LOG_PATH, "r") as f:
            lines = f.readlines()[-CHART_MAX_POINTS:]
        for line in lines:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            _state["history"].append((rec.get("received_ts"), rec.get("payload")))
        if _state["history"]:
            _state["last_payload"] = _state["history"][-1][1]
            _state["last_received_ts"] = _state["history"][-1][0]
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
        _state["last_payload"] = payload
        _state["last_received_ts"] = received_ts
        _state["history"].append((received_ts, payload))
        if len(_state["history"]) > CHART_MAX_POINTS:
            _state["history"] = _state["history"][-CHART_MAX_POINTS:]

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
            if last_ts is None:
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


def _render_chart_svg(history):
    """Minimal inline SVG polyline of people_done over time. No charting library."""
    points = [(ts, p.get("people_done")) for ts, p in history if p and p.get("people_done") is not None]
    if len(points) < 2:
        return "<p class='muted'>Not enough data yet for a chart.</p>"

    w, h, pad = 700, 200, 30
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = 0, max(ys) or 1
    x_span = (x_max - x_min) or 1

    def sx(x):
        return pad + (x - x_min) / x_span * (w - 2 * pad)

    def sy(y):
        return h - pad - (y - y_min) / (y_max - y_min) * (h - 2 * pad)

    poly = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in points)
    return f"""
    <svg viewBox="0 0 {w} {h}" width="100%" height="{h}" style="background:#0b0e14;border-radius:6px">
      <polyline points="{poly}" fill="none" stroke="#4ade80" stroke-width="2" />
      <text x="{pad}" y="16" fill="#94a3b8" font-size="11">people_done over time (last {len(points)} heartbeats)</text>
      <text x="{pad}" y="{h - 8}" fill="#94a3b8" font-size="11">0</text>
      <text x="{w - pad - 40}" y="{h - 8}" fill="#94a3b8" font-size="11">{y_max} people</text>
    </svg>
    """


def _render_dashboard():
    with _lock:
        payload = _state["last_payload"]
        last_ts = _state["last_received_ts"]
        history = list(_state["history"])
        stale = _state["stale_active"]
        anomaly = _state["anomaly_active"]

    if payload is None:
        body = "<p>No heartbeats received yet. Waiting for the pipeline to start.</p>"
        return _page(body)

    age_min = (time.time() - last_ts) / 60.0 if last_ts else None
    status_label = "STALE" if stale else ("ANOMALY" if anomaly else "OK")
    status_color = "#f87171" if stale else ("#fbbf24" if anomaly else "#4ade80")

    pct_done = None
    if payload.get("people_total"):
        pct_done = 100.0 * (payload.get("people_done") or 0) / payload["people_total"]

    rows = [
        ("Status", f"<span style='color:{status_color};font-weight:700'>{status_label}</span>"),
        ("Last heartbeat", f"{age_min:.1f} min ago" if age_min is not None else "—"),
        ("People done / total", f"{payload.get('people_done', '—')} / {payload.get('people_total', '—')}"
                                 + (f" ({pct_done:.1f}%)" if pct_done is not None else "")),
        ("People failed", _fmt(payload.get("people_failed"), digits=0)),
        ("Rolling failure rate", _fmt(payload.get("rolling_failure_rate_pct"), "%")),
        ("Rate", _fmt(payload.get("rate_per_hour"), " people/hr")),
        ("ETA remaining", _fmt(payload.get("eta_hours_remaining"), " hr")),
        ("Elapsed", _fmt(payload.get("elapsed_hours"), " hr")),
        ("Cost so far / budget", f"${_fmt(payload.get('cost_so_far_usd'), digits=2)} / ${_fmt(payload.get('budget_usd'), digits=0)}"),
        ("Disk used", _fmt(payload.get("disk_used_pct"), "%")),
        ("Mem available", _fmt(payload.get("mem_avail_pct"), "%")),
        ("Sender timestamp", html.escape(str(payload.get("ts", "—")))),
    ]
    if payload.get("anomaly_reason"):
        rows.append(("Anomaly reason", html.escape(str(payload["anomaly_reason"]))))

    table_html = "\n".join(
        f"<tr><td class='k'>{html.escape(k)}</td><td class='v'>{v}</td></tr>" for k, v in rows
    )
    chart = _render_chart_svg(history)
    body = f"""
    <table>{table_html}</table>
    <h2>Progress</h2>
    {chart}
    """
    return _page(body)


def _page(body: str) -> str:
    return f"""<!doctype html>
<html><head>
<meta charset="utf-8">
<meta http-equiv="refresh" content="{DASHBOARD_REFRESH_SEC}">
<title>HLA pipeline monitor</title>
<style>
  body {{ background:#0b0e14; color:#e2e8f0; font-family: -apple-system, sans-serif; padding: 24px; }}
  h1 {{ font-size: 20px; }}
  table {{ border-collapse: collapse; width: 100%; max-width: 640px; }}
  td {{ padding: 6px 10px; border-bottom: 1px solid #1e293b; }}
  td.k {{ color: #94a3b8; width: 45%; }}
  td.v {{ font-weight: 600; }}
  .muted {{ color: #64748b; }}
</style>
</head><body>
<h1>Omni-HLA full-cohort run &mdash; monitor</h1>
{body}
<p class="muted">Auto-refreshes every {DASHBOARD_REFRESH_SEC}s. Aggregate metadata only &mdash; no person_ids or alleles ever appear here.</p>
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
