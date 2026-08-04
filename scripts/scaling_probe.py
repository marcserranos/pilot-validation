#!/usr/bin/env python3
"""Minimal, single-purpose horizontal-concurrency probe for Immuannot.

v2 replacement for run_core_scaling_experiment.sh (2026-08-04), after a full day/night of that
script's compounding bash-specific failure modes: a lock defeated by a routine `rm -f`, `du`/`df`
hanging for anywhere from 15 minutes to 2+ hours on this VM's gcsfuse mounts, a dead mount
silently ignored because bash had no `set -e` (so a FATAL Phase 0 didn't stop the script -- it
just wandered into Phase 1 and got stuck there with no visible connection to the real cause), and
background subshells that could be orphaned. Deliberately minimal: ONE config per invocation, ONE
result row, no phases, no nested bash loops, no background subshells, no bash arithmetic. Run it
once per concurrency level you want to test -- simpler to watch, simpler to know exactly what's
happening at any moment, and a crash in one invocation can't corrupt or block another.

Also deliberately scoped down: only tests "many people, 1 thread each" (the horizontal-
concurrency question that actually decides the production VM size). The thread-vs-concurrency
comparison from the old script's Phase 1 is dropped -- interesting, but never the number the VM
decision hinges on, and cutting it removes real complexity for shrinking value. We already have a
solid, repeatedly-confirmed single-person baseline (3.8-6.9 min/person, ~81MB/person disk) from
Phase 0's many successful runs earlier -- no need to re-measure that here.

Usage (from ~/repos/pilot-validation -- does NOT need an activated pixi shell; only the stdlib is
used at this level, `pixi run` is invoked internally per person):
  python3 scripts/scaling_probe.py <concurrency> <person_id_1> [<person_id_2> ...]

Requires exactly `concurrency` person_ids. Appends exactly one row to
~/pipeline_outputs/scaling_probe_results.tsv. Per-person full logs go to
~/pipeline_outputs/scaling_probe_<person_id>.log (overwritten each run, not accumulated).
"""
import concurrent.futures
import os
import subprocess
import sys
import threading
import time

MOUNT = os.path.expanduser("~/mnt/aou-controlled")
MOUNT_CHECK = os.path.join(MOUNT, "v9/wgs/long_read/manifest.tsv")
OUTROOT = os.path.expanduser("~/pipeline_outputs")
RESULTS = os.path.join(OUTROOT, "scaling_probe_results.tsv")
LOCK = os.path.join(OUTROOT, "scaling_probe.lock")


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def check_mount():
    if not os.path.isfile(MOUNT_CHECK):
        die(f"{MOUNT_CHECK} not found -- the gcsfuse mount is not up. This does NOT survive a VM "
            f"restart. Remount, then verify with `ls {MOUNT_CHECK}` yourself, before running this "
            f"again. (This check exists so a dead mount fails loudly here, immediately -- not "
            f"silently, hours later, somewhere else.)")


def check_and_take_lock():
    if os.path.exists(LOCK):
        old_pid = open(LOCK).read().strip()
        try:
            os.kill(int(old_pid), 0)
            die(f"another instance is already running (PID {old_pid}, per {LOCK}). Confirm with "
                f"`ps -p {old_pid}` that it's really gone before removing {LOCK} yourself and "
                f"retrying -- do not just rm this file reflexively.")
        except (ProcessLookupError, ValueError):
            print(f"NOTE: stale lock (PID {old_pid} not running) -- removing.", file=sys.stderr)
            os.remove(LOCK)
    with open(LOCK, "w") as f:
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


def run_one_person(pid):
    t0 = time.time()
    proc = subprocess.run(
        ["pixi", "run", "-e", "specimmune", "--",
         "python3", "scripts/run_immuannot_person.py", str(pid),
         "--threads", "1", "--out-suffix", f".probe_{pid}", "--force"],
        capture_output=True, text=True,
    )
    elapsed = time.time() - t0
    log_path = os.path.join(OUTROOT, f"scaling_probe_{pid}.log")
    with open(log_path, "w") as f:
        f.write(f"exit={proc.returncode}\n--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}\n")
    print(f"  person {pid}: exit={proc.returncode}, {elapsed:.1f}s (full log: {log_path})",
          file=sys.stderr)
    return pid, proc.returncode, elapsed


def main():
    if len(sys.argv) < 3:
        die("usage: scaling_probe.py <concurrency> <person_id_1> [<person_id_2> ...]")
    try:
        concurrency = int(sys.argv[1])
    except ValueError:
        die(f"first argument must be an integer concurrency level, got {sys.argv[1]!r}")
    people = sys.argv[2:]
    if len(people) != concurrency:
        die(f"need exactly {concurrency} person_ids to match concurrency={concurrency}, got "
            f"{len(people)}.")

    check_mount()
    os.makedirs(OUTROOT, exist_ok=True)
    check_and_take_lock()
    try:
        print(f"=== Probing concurrency={concurrency}, people={people} ===", file=sys.stderr)

        mem_readings = []
        stop_event = threading.Event()

        def sample_memory():
            while not stop_event.is_set():
                m = mem_available_mb()
                if m is not None:
                    mem_readings.append(m)
                stop_event.wait(2)

        sampler = threading.Thread(target=sample_memory, daemon=True)
        sampler.start()

        t0 = time.time()
        with concurrent.futures.ThreadPoolExecutor(max_workers=concurrency) as pool:
            results = list(pool.map(run_one_person, people))
        wall = time.time() - t0

        stop_event.set()
        sampler.join(timeout=5)

        n_failed = sum(1 for _, rc, _ in results if rc != 0)
        total_mb = mem_total_mb()
        mem_min = min(mem_readings) if mem_readings else None
        peak_pct = (100 * (total_mb - mem_min) / total_mb) if (mem_min and total_mb) else None

        new_file = not os.path.exists(RESULTS)
        with open(RESULTS, "a") as f:
            if new_file:
                f.write("timestamp\tconcurrency\tn_people\tn_failed\twall_clock_seconds\t"
                         "seconds_per_person\tmem_used_peak_pct\n")
            f.write(f"{time.strftime('%Y-%m-%dT%H:%M:%S')}\t{concurrency}\t{len(people)}\t"
                    f"{n_failed}\t{wall:.1f}\t{wall/len(people):.1f}\t"
                    f"{f'{peak_pct:.1f}' if peak_pct is not None else '-'}\n")

        print(f"=== DONE: concurrency={concurrency}, wall={wall:.1f}s, {n_failed}/{len(people)} "
              f"failed, peak mem {peak_pct:.1f}% -- appended to {RESULTS} ===", file=sys.stderr)
    finally:
        if os.path.exists(LOCK):
            os.remove(LOCK)


if __name__ == "__main__":
    main()
