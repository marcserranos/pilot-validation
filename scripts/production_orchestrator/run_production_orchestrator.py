#!/usr/bin/env python3
"""Full-cohort Immuannot production orchestrator. See BRIEF.md in this directory for the full
design rationale -- every safety behavior here maps to a specific real incident earlier in this
project (ENVIRONMENT.md quirks #11/#14/#17/#18/#23-#27). This script extends scripts/scaling_probe.py's
proven patterns (mount check, PID lock, ThreadPoolExecutor concurrency, real exceptions stop
immediately -- no bash, no set -e ambiguity) rather than starting from zero.

What it adds beyond scaling_probe.py (which only ever probes ONE config, once):
  - Loops over the ENTIRE cohort (from build_immuannot_cohort.py's output), not one fixed batch.
  - Real resumability: before dispatching ANYONE, scans the filesystem for each candidate's actual
    expected output (hap1.gtf.gz / hap2.gtf.gz under ~/pipeline_outputs/<pid>/immuannot_output/) --
    not a log entry, not run_immuannot_person.py's own canonical-file check (which can't see
    concurrent workers' --out-suffix-isolated results -- see that script's own --out-suffix
    docstring). Safe to kill (Ctrl-C, SIGTERM, VM restart) and relaunch with the exact same command.
  - Calls heartbeat_client.send_heartbeat() every --heartbeat-interval-sec (~5 min in production)
    from IN-MEMORY counters only -- never by re-scanning the output directory (monitoring/README.md
    "Cadence" hard requirement; a broad scan every heartbeat would reintroduce the du/df hang
    failure mode, ENVIRONMENT.md quirk #25).
  - Periodically (and always at the end) merges each worker's --out-suffix-isolated
    immuannot_calls.<pid>.tsv / immuannot_timing.<pid>.tsv fragments into the canonical
    immuannot_calls.tsv / immuannot_timing.tsv files -- a bounded glob + concat, not a du/df walk.
  - Tracks per-person attempt counts; after --max-attempts failed tries (default 3) a person is
    marked given-up (not silently retried forever, not silently dropped either -- logged loudly and
    left for manual review) instead of consuming worker slots on a deterministically-failing person
    for the whole multi-day run.
  - Local budget check every heartbeat (cost_so_far_usd vs --budget) -- loud stderr warning, not an
    automatic kill (BRIEF.md: "flag it", the $300 cap is a human decision point, not an autopilot one).

Usage (from ~/repos/pilot-validation, inside `pixi shell -e specimmune` or `pixi run -e specimmune --`):
  python3 scripts/production_orchestrator/run_production_orchestrator.py \\
      --cohort ~/pipeline_outputs/immuannot_cohort_full.tsv \\
      --concurrency 24 --threads-per-person 4 \\
      --vm-rate <REAL Workbench-UI-confirmed USD/hour for the n2-highcpu-96 VM> \\
      --monitor-url http://46.225.123.54:8943

--vm-rate and MONITOR_AUTH_TOKEN (env var, never on the command line) are required for real
monitoring to work end-to-end; the orchestrator still runs and processes people without them, but
prints a loud one-time warning, since an unattended multi-day run with no working heartbeat is
exactly the failure mode monitoring/BRIEF.md was built to prevent.
"""
import argparse
import concurrent.futures
import glob
import gzip
import os
import re
import subprocess
import sys
import threading
import time
from collections import Counter

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "monitoring"))
from heartbeat_client import send_heartbeat  # noqa: E402

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
MOUNT_DEFAULT = os.path.expanduser("~/mnt/aou-controlled")
MOUNT_CHECK_REL = "v9/wgs/long_read/manifest.tsv"
OUTROOT_DEFAULT = os.path.expanduser("~/pipeline_outputs")
DEFAULT_COHORT = os.path.join(OUTROOT_DEFAULT, "immuannot_cohort_full.tsv")
DEFAULT_MONITOR_URL = "http://46.225.123.54:8943"
DEFAULT_REGION = "chr6:29500000-33500000"
DEFAULT_PAD = 100_000


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def warn(msg):
    print(f"WARNING: {msg}", file=sys.stderr)


def check_mount(mount):
    check_path = os.path.join(mount, MOUNT_CHECK_REL)
    if not os.path.isfile(check_path):
        die(f"{check_path} not found -- the gcsfuse mount is not up. This does NOT survive a VM "
            f"restart/stop (ENVIRONMENT.md quirk #11/#14/#26). Remount, then verify with "
            f"`ls {check_path}` yourself, before running this again. Refusing to start rather than "
            f"wander deep into a multi-day run with a dead mount and fail silently hours later.")


def check_and_take_lock(lock_path):
    if os.path.exists(lock_path):
        old_pid = open(lock_path).read().strip()
        try:
            os.kill(int(old_pid), 0)
            die(f"another instance is already running (PID {old_pid}, per {lock_path}). Two "
                f"concurrent orchestrators against overlapping people WILL silently corrupt shared "
                f"per-person intermediate files (ENVIRONMENT.md quirk #23 -- this happened, more "
                f"than once, in this project). Confirm with `ps -p {old_pid}` that it's really "
                f"gone before removing {lock_path} yourself and retrying -- do not rm reflexively.")
        except (ProcessLookupError, ValueError):
            warn(f"stale lock (PID {old_pid} not running) -- removing.")
            os.remove(lock_path)
    with open(lock_path, "w") as f:
        f.write(str(os.getpid()))


def mem_total_mb():
    for line in open("/proc/meminfo"):
        if line.startswith("MemTotal"):
            return int(line.split()[1]) / 1024
    return None


def mem_available_mb():
    for line in open("/proc/meminfo"):
        if line.startswith("MemAvailable"):
            return int(line.split()[1]) / 1024
    return None


def mem_avail_pct():
    total, avail = mem_total_mb(), mem_available_mb()
    if not total or avail is None:
        return None
    return 100.0 * (total - avail) / total  # reported as "used pct", matches heartbeat schema


def disk_used_pct(path, timeout=5):
    """A single statvfs(2) syscall on the filesystem containing `path` (local persistent disk,
    NOT the gcsfuse mount -- ~/pipeline_outputs is never gcsfuse-backed, ENVIRONMENT.md quirk #12).
    Deliberately NOT `du`/`df` as a subprocess: quirk #25's corrected finding is that even `df`
    can hang for 2+ hours on this VM when ANY mounted filesystem (including an unrelated FUSE
    mount elsewhere on the box) is slow to respond, because it can end up touching the whole mount
    table. os.statvfs() queries just the one target filesystem's superblock in-process. Still run
    inside a bounded-timeout thread as a hard safety net (same discipline as quirk #25's `timeout 5`
    fix), since this project has learned twice now not to fully trust a command's own reputation
    for being fast."""
    result = {}

    def _stat():
        try:
            st = os.statvfs(path)
            total = st.f_blocks * st.f_frsize
            free = st.f_bavail * st.f_frsize
            if total:
                result["pct"] = 100.0 * (total - free) / total
        except OSError:
            pass

    t = threading.Thread(target=_stat, daemon=True)
    t.start()
    t.join(timeout=timeout)
    return result.get("pct")


def person_output_paths(outroot, pid):
    person_dir = os.path.join(outroot, str(pid), "immuannot_output")
    return (os.path.join(person_dir, "hap1.gtf.gz"), os.path.join(person_dir, "hap2.gtf.gz"))


_GTF_CALL_RE = re.compile(r'(consensus|allele) "[^"]+"')


def gtf_has_real_call(gtf_gz_path):
    """File existence alone is NOT sufficient (ENVIRONMENT.md quirks #17/#18: both SpecImmune and
    SpecHLA have separately, actually shipped a real, existing output file while exiting 0 and
    having silently produced near-zero real content -- a coarse success signal like exit code or
    file existence was not enough for either tool, and there's no reason to assume Immuannot is
    immune to the same class of bug just because it hasn't been caught here yet). Cheap check:
    the gtf.gz decompresses and contains at least one real gene/allele attribute line -- mirrors
    run_immuannot_person.py's own parse_gtf() regex, just as a boolean instead of building the
    full calls dict. A truly empty or corrupt gtf.gz (0 real calls) is treated as NOT done, so a
    relaunch retries it instead of permanently trusting a hollow success."""
    try:
        with gzip.open(gtf_gz_path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if _GTF_CALL_RE.search(line):
                    return True
        return False
    except OSError:
        return False


def person_done(outroot, pid):
    """Real resumability check, per BRIEF.md: the expected OUTPUT file, not a log entry and not
    run_immuannot_person.py's own canonical-immuannot_calls.tsv check (which is blind to
    --out-suffix-isolated concurrent results -- see that script's --out-suffix docstring). Both
    haplotype .gtf.gz files must exist AND actually contain >=1 real call each -- see
    gtf_has_real_call()."""
    hap1_gtf, hap2_gtf = person_output_paths(outroot, pid)
    return (os.path.isfile(hap1_gtf) and os.path.isfile(hap2_gtf)
            and gtf_has_real_call(hap1_gtf) and gtf_has_real_call(hap2_gtf))


def gave_up_path(outroot, pid):
    return os.path.join(outroot, str(pid), "immuannot_output", ".orchestrator_gave_up")


def attempts_path(outroot, pid):
    return os.path.join(outroot, str(pid), "immuannot_output", ".orchestrator_attempts")


def person_gave_up(outroot, pid):
    return os.path.isfile(gave_up_path(outroot, pid))


def get_attempts(outroot, pid):
    p = attempts_path(outroot, pid)
    if not os.path.isfile(p):
        return 0
    try:
        return int(open(p).read().strip() or "0")
    except ValueError:
        return 0


def bump_attempts(outroot, pid):
    n = get_attempts(outroot, pid) + 1
    os.makedirs(os.path.dirname(attempts_path(outroot, pid)), exist_ok=True)
    with open(attempts_path(outroot, pid), "w") as f:
        f.write(str(n))
    return n


def load_cohort(cohort_path, skip_trim_tiers):
    """Returns (people, ancestry_map). ancestry_map is {person_id: ancestry_pred_or_'NA'} --
    build_immuannot_cohort.py's ancestry join (2026-08-04), used ONLY for local per-ancestry
    progress logging (see log_ancestry_progress()) -- never sent to the remote heartbeat dashboard,
    which stays aggregate-counts-only per monitoring/README.md. If the cohort file predates the
    ancestry join (no ancestry_pred column), every person maps to 'NA' -- degrades gracefully,
    doesn't block a launch."""
    if not os.path.isfile(cohort_path):
        die(f"cohort file not found: {cohort_path} -- run build_immuannot_cohort.py first.")
    df = pd.read_csv(cohort_path, sep="\t", dtype=str)
    for c in ["person_id", "trim_tier"]:
        if c not in df.columns:
            die(f"cohort file {cohort_path} missing expected column '{c}'.")
    skip = set(skip_trim_tiers)
    if skip:
        before = len(df)
        df = df[~df["trim_tier"].isin(skip)]
        print(f"--skip-trim-tier {sorted(skip)}: excluded {before - len(df)}/{before} people from "
              f"this launch (they remain in the cohort file for a later batch).", file=sys.stderr)
    if "ancestry_pred" in df.columns:
        ancestry_map = {pid: (a if pd.notna(a) else "NA")
                        for pid, a in zip(df["person_id"], df["ancestry_pred"])}
    else:
        ancestry_map = {pid: "NA" for pid in df["person_id"]}
    return list(df["person_id"]), ancestry_map


def scan_already_done(outroot, people, max_workers=32):
    """One bounded, finite-size check (2 os.path.isfile() calls per candidate, on local disk, not
    the gcsfuse mount) -- NOT a recursive tree walk. This is the startup resumability scan; after
    this, all progress tracking is in-memory only (monitoring/README.md's hard requirement)."""
    done, gave_up = set(), set()
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
        done_flags = list(pool.map(lambda pid: person_done(outroot, pid), people))
        gave_up_flags = list(pool.map(lambda pid: person_gave_up(outroot, pid), people))
    for pid, d, g in zip(people, done_flags, gave_up_flags):
        if d:
            done.add(pid)
        elif g:
            gave_up.add(pid)
    return done, gave_up


def run_one_person(pid, args):
    cmd = [
        "pixi", "run", "--manifest-path", os.path.join(REPO_ROOT, "pixi.toml"), "-e", "specimmune", "--",
        "python3", os.path.join(REPO_ROOT, "scripts", "run_immuannot_person.py"), str(pid),
        "--mount", args.mount, "--outroot", args.outroot,
        "--immuannot-dir", args.immuannot_dir, "--refdir", args.refdir,
        "--threads", str(args.threads_per_person),
        "--region", args.region, "--pad", str(args.pad),
        "--out-suffix", f".{pid}", "--force",
    ]
    if args.enable_self_align_fallback:
        cmd.append("--enable-self-align-fallback")
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - t0
    log_path = os.path.join(args.outroot, str(pid), "immuannot_output", ".orchestrator_last_run.log")
    try:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        with open(log_path, "w") as f:
            f.write(f"exit={proc.returncode} elapsed={elapsed:.1f}s\n--- stdout ---\n{proc.stdout}\n"
                     f"--- stderr ---\n{proc.stderr}\n")
    except OSError as e:
        warn(f"person {pid}: could not write {log_path}: {e}")
    return pid, proc.returncode, elapsed


def merge_fragments(outroot):
    """Merges every --out-suffix-isolated fragment into the canonical immuannot_calls.tsv /
    immuannot_timing.tsv. A bounded glob over already-produced small files, not a du/df tree walk --
    safe to call periodically. Single-writer (only the orchestrator's own main thread calls this),
    so no read-modify-write race despite write_incremental()'s own non-atomicity."""
    total_fragments = 0
    for base in ("immuannot_calls", "immuannot_timing"):
        fragments = sorted(glob.glob(os.path.join(outroot, f"{base}.*.tsv")))
        if not fragments:
            continue
        total_fragments += len(fragments)
        frames = []
        canonical = os.path.join(outroot, f"{base}.tsv")
        if os.path.exists(canonical):
            frames.append(pd.read_csv(canonical, sep="\t", dtype=str))
        for frag in fragments:
            try:
                frames.append(pd.read_csv(frag, sep="\t", dtype=str))
            except pd.errors.EmptyDataError:
                continue
        if not frames:
            continue
        merged = pd.concat(frames, ignore_index=True)
        key = "person_id"
        merged = merged.drop_duplicates(subset=[key] if base == "immuannot_calls" else [key, "hap"],
                                         keep="last")
        merged.to_csv(canonical, sep="\t", index=False)
        # Delete merged fragments once safely folded into the canonical file -- keeps every future
        # merge bounded to only NEW fragments since the last merge, instead of re-reading the whole
        # growing fragment set every time (would otherwise be O(n^2) I/O over a multi-day run).
        for frag in fragments:
            try:
                os.remove(frag)
            except OSError:
                pass
    if total_fragments:
        print(f"  [merge] canonical immuannot_calls.tsv / immuannot_timing.tsv refreshed from "
              f"{total_fragments} fragment(s).", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", default=DEFAULT_COHORT)
    ap.add_argument("--mount", default=MOUNT_DEFAULT)
    ap.add_argument("--outroot", default=OUTROOT_DEFAULT)
    ap.add_argument("--immuannot-dir", default=os.path.expanduser("~/tools/Immuannot"))
    ap.add_argument("--refdir", default=os.path.expanduser("~/tools/Immuannot_refdata"))
    ap.add_argument("--concurrency", type=int, default=24,
                    help="People processed at once (default 24 -- the committed production "
                         "config, DECISIONS.md 2026-08-04).")
    ap.add_argument("--threads-per-person", type=int, default=4,
                    help="Threads per Immuannot invocation (default 4, same committed config).")
    ap.add_argument("--region", default=DEFAULT_REGION)
    ap.add_argument("--pad", type=int, default=DEFAULT_PAD)
    ap.add_argument("--enable-self-align-fallback", action="store_true",
                    help="Pass through to every worker -- attempts Tier 3 (self-align) for people "
                         "whose only resolving trim_tier is self_align_needed. UNTESTED AT SCALE "
                         "as of 2026-08-04 -- do not pass this for the real production launch until "
                         "the sequel2 smoke test (BRIEF.md) has actually been run and reviewed.")
    ap.add_argument("--skip-trim-tier", action="append", default=[],
                    help="Exclude a trim_tier from this launch (repeatable). E.g. "
                         "--skip-trim-tier self_align_needed to hold sequel2-shaped people back "
                         "for a later batch. Excluded people stay in the cohort file, untouched.")
    ap.add_argument("--max-attempts", type=int, default=3,
                    help="After this many failed attempts, a person is marked given-up (not "
                         "retried on future relaunches) instead of consuming a worker slot on a "
                         "deterministically-failing person for the rest of a multi-day run.")
    ap.add_argument("--heartbeat-interval-sec", type=float, default=300,
                    help="Default 300s (~5 min) -- the committed cadence, monitoring/README.md "
                         "'Cadence'. Do not go below ~60s; heartbeat overhead is negligible but "
                         "there's no benefit past human reaction time.")
    ap.add_argument("--monitor-url", default=os.environ.get("MONITOR_URL", DEFAULT_MONITOR_URL))
    ap.add_argument("--vm-rate", type=float, default=None,
                    help="REAL Workbench-UI-confirmed USD/hour for the running VM -- required for "
                         "cost_so_far_usd / budget tracking to mean anything. DECISIONS.md's "
                         "~$3.03/hr is a research estimate, not a quoted price -- confirm live.")
    ap.add_argument("--budget", type=float, default=300.0)
    args = ap.parse_args()

    total_cores = args.concurrency * args.threads_per_person
    n_cpu = os.cpu_count() or 0
    if total_cores > n_cpu:
        die(f"concurrency({args.concurrency}) * threads_per_person({args.threads_per_person}) = "
            f"{total_cores} cores, but this machine only has {n_cpu}. Refusing to oversubscribe.")

    check_mount(args.mount)
    os.makedirs(args.outroot, exist_ok=True)
    lock_path = os.path.join(args.outroot, "production_orchestrator.lock")
    check_and_take_lock(lock_path)

    auth_token = os.environ.get("MONITOR_AUTH_TOKEN")
    if not auth_token:
        warn("MONITOR_AUTH_TOKEN not set -- heartbeats will fail every ~5 min for the whole run "
             "(non-fatal to the pipeline, but you will be flying blind on an unattended multi-day "
             "job). Set it in the environment before launching for real.")
    if args.vm_rate is None:
        warn("--vm-rate not passed -- cost_so_far_usd will be null in every heartbeat and the "
             "local budget check below is disabled. Confirm the real n2-highcpu-96 rate in the "
             "Workbench UI and pass it.")

    try:
        people, ancestry_map = load_cohort(args.cohort, args.skip_trim_tier)
        print(f"Cohort: {len(people)} people loaded from {args.cohort} (after any --skip-trim-tier "
              f"filtering).", file=sys.stderr)

        print("Scanning for already-completed people (real .gtf.gz existence + real-call check, "
              "not a log entry -- this is what makes killing and relaunching safe)...", file=sys.stderr)
        done_set, gave_up_set = scan_already_done(args.outroot, people)
        worklist = [p for p in people if p not in done_set and p not in gave_up_set]
        print(f"  {len(done_set)} already done (skipped), {len(gave_up_set)} previously given up "
              f"after {args.max_attempts}+ attempts (skipped -- see .orchestrator_gave_up markers "
              f"for review), {len(worklist)} to process this launch.", file=sys.stderr)

        # Per-ancestry progress (2026-08-04) -- LOCAL log only, never sent to the remote heartbeat
        # dashboard (monitoring/README.md keeps that aggregate-counts-only by design). In-memory
        # counters only, same as the aggregate done/failed counts -- never re-derived by scanning
        # the output directory. Motivation: sequel2 (the self_align_needed / --skip-trim-tier
        # candidate) is ~95% AFR, so a failure mode that disproportionately hits one ancestry group
        # should be visible live, not discovered after the fact in a downstream analysis.
        anc_total = Counter(ancestry_map.get(p, "NA") for p in people)
        anc_done = Counter(ancestry_map.get(p, "NA") for p in done_set)
        print("  Cohort ancestry distribution: " +
              ", ".join(f"{a}={n}" for a, n in anc_total.most_common()), file=sys.stderr)

        if not worklist:
            print("Nothing left to do -- every person is already done or given-up.", file=sys.stderr)
            merge_fragments(args.outroot)
            return

        state = {"done": len(done_set), "failed": len(gave_up_set), "lock": threading.Lock(),
                 "anc_done": anc_done}
        people_total = len(people)
        stop_event = threading.Event()

        def heartbeat_loop():
            while not stop_event.is_set():
                with state["lock"]:
                    d, f = state["done"], state["failed"]
                disk_pct = disk_used_pct(args.outroot)
                mem_pct = mem_avail_pct()
                status, payload = send_heartbeat(
                    d, f, people_total, disk_pct, mem_pct,
                    receiver_url=args.monitor_url, auth_token=auth_token,
                    vm_hourly_rate_usd=args.vm_rate, budget_usd=args.budget,
                )
                cost = payload.get("cost_so_far_usd")
                if cost is not None and args.budget and cost > args.budget:
                    warn(f"cost_so_far_usd={cost:.2f} has EXCEEDED --budget={args.budget:.2f} -- "
                         f"this is a flag, not an auto-stop (BRIEF.md). Decide whether to let it "
                         f"keep running.")
                if payload.get("anomaly"):
                    warn(f"heartbeat anomaly: {payload.get('anomaly_reason')}")
                print(f"  [heartbeat] status={status} done={d} failed={f}/{people_total} "
                      f"rate={payload.get('rate_per_hour')}/hr eta_hr={payload.get('eta_hours_remaining')} "
                      f"cost=${cost}", file=sys.stderr)
                with state["lock"]:
                    anc_done_snapshot = dict(state["anc_done"])
                print("  [heartbeat] per-ancestry done/total (local log only, not sent to the "
                      "remote dashboard): " +
                      ", ".join(f"{a}={anc_done_snapshot.get(a, 0)}/{n}"
                                for a, n in anc_total.most_common()), file=sys.stderr)
                stop_event.wait(args.heartbeat_interval_sec)

        hb_thread = threading.Thread(target=heartbeat_loop, daemon=True)
        hb_thread.start()

        print(f"=== Launching: concurrency={args.concurrency}, threads_per_person="
              f"{args.threads_per_person} (total_cores={total_cores}), {len(worklist)} people "
              f"===", file=sys.stderr)

        n_since_merge = 0
        t0 = time.time()
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.concurrency) as pool:
            futures = {pool.submit(run_one_person, pid, args): pid for pid in worklist}
            for i, fut in enumerate(concurrent.futures.as_completed(futures), 1):
                pid = futures[fut]
                try:
                    _, rc, elapsed = fut.result()
                except Exception as e:  # a worker itself crashing is a real bug -- surface it loudly
                    warn(f"person {pid}: worker raised {e!r} (not a subprocess failure -- a real "
                         f"orchestrator-side exception).")
                    rc, elapsed = -1, 0.0

                ok = person_done(args.outroot, pid)  # never trust exit code alone
                if ok:
                    with state["lock"]:
                        state["done"] += 1
                        state["anc_done"][ancestry_map.get(pid, "NA")] += 1
                    print(f"  [{i}/{len(worklist)}] {pid}: DONE ({elapsed:.0f}s)", file=sys.stderr)
                else:
                    n_attempts = bump_attempts(args.outroot, pid)
                    if n_attempts >= args.max_attempts:
                        os.makedirs(os.path.dirname(gave_up_path(args.outroot, pid)), exist_ok=True)
                        with open(gave_up_path(args.outroot, pid), "w") as f:
                            f.write(f"gave up after {n_attempts} attempts, exit={rc}\n")
                        with state["lock"]:
                            state["failed"] += 1
                        warn(f"person {pid}: GAVE UP after {n_attempts} attempts (exit={rc}, "
                             f"{elapsed:.0f}s) -- see .orchestrator_gave_up / .orchestrator_last_run.log")
                    else:
                        print(f"  [{i}/{len(worklist)}] {pid}: FAILED (exit={rc}, {elapsed:.0f}s, "
                              f"attempt {n_attempts}/{args.max_attempts} -- will retry on next "
                              f"relaunch)", file=sys.stderr)

                n_since_merge += 1
                if n_since_merge >= 500:
                    merge_fragments(args.outroot)
                    n_since_merge = 0

        wall = time.time() - t0
        stop_event.set()
        hb_thread.join(timeout=5)
        merge_fragments(args.outroot)

        with state["lock"]:
            d, f = state["done"], state["failed"]
        print(f"=== DONE: {d}/{people_total} done, {f}/{people_total} given-up, wall={wall/3600:.2f}h "
              f"===", file=sys.stderr)
    finally:
        if os.path.exists(lock_path):
            os.remove(lock_path)


if __name__ == "__main__":
    main()
