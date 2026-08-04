# Pre-flight checklist — before turning on the expensive VM and walking away

> **Rolling checklist**, updated as items are cleared. Last updated 2026-08-05.
> Companion to [`BRIEF.md`](BRIEF.md) (the why) and [`RESULTS_LOCATION.md`](RESULTS_LOCATION.md)
> (where output lands + the launch command).
>
> The point of this file: the run costs ~$3.55/h and is meant to be left unattended overnight for
> ~2.5 days. Everything below is something that, if wrong, is either expensive, silent, or both.
> **Do not skip an item because it "should be fine" — that specific instinct is what
> ENVIRONMENT.md quirks #23-#27 are a record of.**

---

## ✅ Already verified (2026-08-05) — no action needed

| # | Item | Evidence |
|---|---|---|
| 1 | **Hetzner firewall, port 8943, reachable from the open internet** | `GET /health` → 200 and `GET /` → 401 from a machine outside the box. Both the heartbeat POST and the dashboard use this one port. |
| 2 | **`hla-monitor` service is live, enabled, and on the current code** | `systemctl is-active/is-enabled` → active/enabled; box repo at the latest commit. |
| 3 | **Heartbeat POST → receiver → JSONL → dashboard works end to end** | 4 synthetic beats over the public IP: all `status=200`, all persisted, dashboard rendered them. **This had never actually been tested against the deployed box before** — prior testing in `monitoring/README.md` was loopback-only. |
| 4 | **Every heartbeat is persisted forever** | `MONITOR_LOG_PATH=/root/hla-monitor/heartbeats.jsonl`, appended + `fsync`ed per beat. The dashboard chart only *displays* the last 500 points; the file keeps everything, for later graphing. |
| 5 | **Cost model live, with the real quoted rates** | `$3.55/h` VM + `$81.60/mo` disk (at 730 h/mo). Verified: a 41.1 h ETA renders `$150.53` projected, matching hand-calc `$150.5` — and consistent with DECISIONS.md's $160-200 estimate. Turns red over the `$300` budget. |
| 6 | **Cost clock survives a receiver restart, ignores old test beats** | Computed server-side from the first beat of the current run *session*; a gap > 90 min starts a new session. |
| 7 | **Stale-alert threshold corrected for the ~5 min cadence** | `MONITOR_STALE_MINUTES` 25 → 15 (`monitoring/README.md` flagged this as "not yet done"). |
| 8 | **Synthetic test beats archived, dashboard starts clean** | Moved to `heartbeats.synthetic-test-2026-08-05.jsonl`; dashboard shows `WAITING`. |
| 9 | **Inverted memory metric fixed** | `mem_avail_pct` was sending memory *used* while both consumers alert when it drops *below* 10%, reading it as *available*. Would have fired a spurious anomaly push early in the run and stayed silent during a real near-OOM. |

---

## ⬜ Still to do — blocking

### A. Prove the Workbench VM can actually reach the Hetzner box

**The single biggest untested assumption.** Item 3 above proves the receiver accepts beats from the
public internet — but it was sent *from the Hetzner box itself*. The production VM is on a
different network, and the workspace UI explicitly warns:

> *"Full Cohort HLA Calling is part of a perimeter that restricts network services. Access to
> services outside the perimeter may fail."*

That's **VPC Service Controls**, and egress to an arbitrary IP on a non-standard port is exactly
what it's designed to restrict. If this is blocked, the run still works — but you are flying blind
for 2.5 days, which is the whole failure mode the monitoring was built to prevent.

**Test it from the production VM the moment it exists, before launching:**
```bash
export MONITOR_AUTH_TOKEN=<from ~/hla-monitor/.env on the Hetzner box>
cd ~/repos/pilot-validation/scripts/monitoring
python3 heartbeat_client.py --url http://46.225.123.54:8943 --token "$MONITOR_AUTH_TOKEN" \
  --people-done 1 --people-total 13252 --vm-rate 3.55
```
`status=200` → good. Timeout/connection refused → VPC-SC is blocking egress; **tell me before
launching** and we'll decide (allowlist the IP in the perimeter, or fall back to a VM-local
progress file you check by hand).

### B. Subscribe to the ntfy topic on your phone

Alerts are useless if nobody receives them. Get the topic from `MONITOR_NTFY_TOPIC` in
`~/hla-monitor/.env` on the box, then install the ntfy app (iOS/Android) and subscribe, **or** open
`https://ntfy.sh/<topic>` in a browser. Confirm you actually get a push before relying on it.

### C. Production VM setup

Run `bootstrap_vm.sh` on the fresh VM (it's idempotent, handles pixi + repo + Immuannot + mount, and
as of today also fetches the hg38 reference):
```bash
curl -fsSL https://raw.githubusercontent.com/marcserranos/pilot-validation/main/scripts/bootstrap_vm.sh | bash
```
Then confirm each of these on the production VM specifically:
- [ ] **2TB data disk attached** (DECISIONS.md: ~1.15TB unpruned + margin).
- [ ] **96 vCPU confirmed** — the orchestrator refuses to start if `concurrency × threads > os.cpu_count()`, so a smaller-than-expected VM fails loudly, but check anyway.
- [ ] **gcsfuse mount up and verified** with a real `ls` (quirk #11/#14).
- [ ] **`~/ref/Homo_sapiens_assembly38.fasta` present** — needed only by phase 2's self-align. The orchestrator now warns loudly at launch if it's missing, and you have ~2 days of phase-1 headroom to fix it, but easier to just have it.
- [ ] **The cohort file is on *this* VM**: `~/pipeline_outputs/immuannot_cohort_full.tsv`. It was built on a different VM — **if the production VM is a fresh machine, this file does not exist there.** Either re-run `build_immuannot_cohort.py` (slow-ish, needs the mount) or copy it across. Confirm it has 13,252 rows and an `ancestry_pred` column.

### D. Smoke test — BRIEF.md's "don't skip this"

Never yet run. On the production VM once it's up, before the real launch:
```bash
python3 scripts/production_orchestrator/run_production_orchestrator.py \
    --single-phase --skip-trim-tier self_align_needed \
    --concurrency 6 --threads-per-person 4 --vm-rate 3.55 \
    2>&1 | tee -a ~/pipeline_outputs/smoke_test.log
```
Let ~24 people finish, then verify all four safety behaviors:
- [ ] **Resumability** — Ctrl-C mid-run, relaunch the same command, confirm it reports the finished people as already done and doesn't redo them.
- [ ] **Lock** — with it running, launch a second copy in another terminal; it must refuse loudly (quirk #23).
- [ ] **Mount check** — `fusermount -u ~/mnt/aou-controlled`, try to launch, confirm it dies immediately with the mount error; then remount (quirk #26).
- [ ] **Heartbeat** — confirm beats appear on the dashboard with real numbers (this also re-proves item A under real conditions).

Then clear the smoke-test people so they don't skew phase 1's timing stats — or just leave them,
they're legitimately done and resumability will skip them.

---

## ⬜ Still to do — non-blocking, but decide before you walk away

- [ ] **The 96-core extrapolation is still unconfirmed.** Every timing number (~52-58 h, ~$160-200) extrapolates linearly from real 32-core data and has never been tested above 32 cores (DECISIONS.md's own caveat). The smoke test above will give the first real signal — if throughput at 96 cores is meaningfully below ~249 people/h, the ETA and cost projections shift. Not a reason to delay, just don't be surprised.
- [ ] **Phase 2 (sequel2) remains genuinely untested.** By design — that's the whole point of running it second, watched, abortable. Expect to actually look at the dashboard when phase 1 ends.
- [ ] **`--vm-rate 3.55` must be passed** or cost tracking silently does nothing (the orchestrator warns, but the warning scrolls past).
- [ ] **`MONITOR_AUTH_TOKEN` must be exported** in the launching shell — same deal.
- [ ] **Launch inside `tmux` with `tee`** (quirk #14). Not optional: a dropped browser/SSH session otherwise `SIGHUP`s the whole run, and the aftermath is indistinguishable from a VM restart.

---

## Launch command (once everything above is green)

See [`RESULTS_LOCATION.md`](RESULTS_LOCATION.md) — one command runs both phases automatically.
