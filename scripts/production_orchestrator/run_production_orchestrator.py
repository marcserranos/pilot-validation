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
  - AUTOMATIC TWO-PHASE RUN (default, 2026-08-05, Marc: "I don't want to intervene... make
    everything on the same run"): a single invocation runs phase 1 (everyone except the untested
    self_align_needed/sequel2 tier) to completion, then AUTOMATICALLY continues straight into phase
    2 (only self_align_needed, with the self-align fallback enabled) -- same process, same lock, no
    second command, no manual step. Each phase gets its own heartbeat state file so phase 2's
    rate/ETA/cost are never diluted by phase 1's already-elapsed hours -- watch the dashboard during
    phase 2 and Ctrl-C if it looks too slow/expensive (this waits for the current in-flight wave to
    finish, not an instant kill; a relaunch with the same command resumes correctly either way).
    Pass --single-phase to disable this and run one filtered pass instead (see --skip-trim-tier/
    --only-trim-tier), e.g. for the smoke test.

Usage (from ~/repos/pilot-validation, inside `pixi shell -e specimmune` or `pixi run -e specimmune --`):
  python3 scripts/production_orchestrator/run_production_orchestrator.py \\
      --cohort ~/pipeline_outputs/immuannot_cohort_full.tsv \\
      --concurrency 24 --threads-per-person 4 \\
      --vm-rate <REAL Workbench-UI-confirmed USD/hour for the n2-highcpu-96 VM> \\
      --monitor-url http://46.225.123.54:8943
That single command runs BOTH phases, in order, automatically -- nothing else to type.

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
    """Percent of memory AVAILABLE (not used).

    BUGFIX 2026-08-05: this previously returned (total-avail)/total -- i.e. memory USED --
    while both consumers of the `mem_avail_pct` heartbeat field (heartbeat_client's
    build_payload() and heartbeat_receiver's _resource_anomaly_reasons()) alert when the
    value drops BELOW MEM_AVAIL_DANGER_PCT (10), i.e. they correctly read the field by its
    name, as AVAILABLE. The inverted sense meant a healthy machine early in a run (say 8%
    used, 92% free) would trip the "< 10%" danger rule and fire a spurious anomaly ntfy
    push, while a genuine near-OOM (95% used, 5% free) would report 95 and never alert --
    exactly backwards, on an unattended multi-day run. Now matches the field name and both
    consumers."""
    total, avail = mem_total_mb(), mem_available_mb()
    if not total or avail is None:
        return None
    return 100.0 * avail / total


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


SEQUEL2_TIER = "self_align_needed"
# Must match run_immuannot_person.py's own defaults -- duplicated as plain constants rather than
# imported, since that module pulls in pandas and lives outside this package. Only used for the
# up-front phase-2 prerequisite check below.
HG38_REF_DEFAULT = os.path.expanduser("~/ref/Homo_sapiens_assembly38.fasta")
CHR6_REF_CACHE_DEFAULT = os.path.expanduser("~/ref/chr6.fasta")


def send_lifecycle_heartbeat(args, auth_token, run_state, people_done, people_failed,
                             people_total, phase=None, state_path=None):
    """One-off heartbeat outside the periodic loop, for a lifecycle transition.

    Two reasons this exists (2026-08-05):
      - STARTUP: if the orchestrator dies early (bad flag, missing cohort, crash on the first
        person), the periodic loop may never have run, so the dashboard sits on WAITING forever
        AND the receiver's stale watchdog never arms -- it skips while last_received_ts is None.
        Result: total silence on a run you thought you'd launched. One beat at startup fixes both,
        and doubles as instant confirmation the VM->Hetzner path actually works, while you're
        still watching.
      - COMPLETE/FAILED: heartbeats stopping is otherwise ambiguous between "finished fine" and
        "crashed", and the watchdog would page in the middle of the night either way.
    """
    kwargs = dict(receiver_url=args.monitor_url, auth_token=auth_token,
                  vm_hourly_rate_usd=args.vm_rate, budget_usd=args.budget,
                  run_state=run_state, phase=phase)
    if state_path:
        kwargs["state_path"] = state_path
    try:
        status, _ = send_heartbeat(people_done, people_failed, people_total,
                                    disk_used_pct(args.outroot), mem_avail_pct(), **kwargs)
        print(f"  [heartbeat:{run_state}] status={status}", file=sys.stderr)
        return status
    except Exception as e:
        # Monitoring must never take down the pipeline it watches.
        warn(f"lifecycle heartbeat ({run_state}) failed: {e!r}")
        return None


def check_phase2_prereqs():
    """Tier 3 (self-align) needs a chr6 reference slice, carved once from the full hg38 FASTA.
    bootstrap_vm.sh does NOT download that FASTA, so a freshly-bootstrapped production VM can
    reach phase 2 without it -- and the failure mode is quiet: ensure_chr6_ref() returns None and
    every sequel2 person is SKIPped one by one, ~40 hours into the run, long after anyone is
    watching. Checked ONCE at launch instead, in the same spirit as the hard mount check
    (BRIEF.md requirement #1): surface it while the human is still at the keyboard.

    Returns a problem string, or None if phase 2 has what it needs."""
    if os.path.isfile(CHR6_REF_CACHE_DEFAULT) and os.path.getsize(CHR6_REF_CACHE_DEFAULT) > 0:
        return None  # already carved by a previous run -- the full FASTA is no longer needed
    if os.path.isfile(HG38_REF_DEFAULT):
        return None
    return (
        f"phase 2 (the {SEQUEL2_TIER} / sequel2 group) needs a chr6 reference slice, carved from "
        f"{HG38_REF_DEFAULT} -- but NEITHER that file NOR a previously-carved {CHR6_REF_CACHE_DEFAULT} "
        f"exists on this VM. bootstrap_vm.sh does not fetch it. Phase 1 is unaffected and will run "
        f"fine, but every phase-2 person would be skipped. Fix it any time BEFORE phase 1 ends "
        f"(~2 days of headroom) with:\n"
        f"    mkdir -p ~/ref && gcloud storage cp \\\n"
        f"      gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta \\\n"
        f"      gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai \\\n"
        f"      ~/ref/\n"
        f"  (~3GB, public bucket -- plain gs:// works, no requester-pays flags needed, "
        f"ENVIRONMENT.md quirk #11.)"
    )


def load_cohort(cohort_path):
    """Returns (people, ancestry_map, trim_tier_map) for the FULL cohort file, unfiltered.
    ancestry_map is {person_id: ancestry_pred_or_'NA'} -- build_immuannot_cohort.py's ancestry join
    (2026-08-04), used ONLY for local per-ancestry progress logging, never sent to the remote
    heartbeat dashboard (monitoring/README.md keeps that aggregate-counts-only). If the cohort file
    predates the ancestry join (no ancestry_pred column), every person maps to 'NA' -- degrades
    gracefully, doesn't block a launch. Filtering by trim_tier is the CALLER's job (see
    filter_people()) -- kept separate from loading so the same load can serve both the default
    automatic two-phase split and a --single-phase manual filter."""
    if not os.path.isfile(cohort_path):
        die(f"cohort file not found: {cohort_path} -- run build_immuannot_cohort.py first.")
    df = pd.read_csv(cohort_path, sep="\t", dtype=str)
    for c in ["person_id", "trim_tier"]:
        if c not in df.columns:
            die(f"cohort file {cohort_path} missing expected column '{c}'.")
    if "ancestry_pred" in df.columns:
        ancestry_map = {pid: (a if pd.notna(a) else "NA")
                        for pid, a in zip(df["person_id"], df["ancestry_pred"])}
    else:
        ancestry_map = {pid: "NA" for pid in df["person_id"]}
    trim_tier_map = dict(zip(df["person_id"], df["trim_tier"]))
    return list(df["person_id"]), ancestry_map, trim_tier_map


def filter_people(people, trim_tier_map, skip_tiers=None, only_tiers=None):
    skip, only = set(skip_tiers or []), set(only_tiers or [])
    out = people
    if skip:
        out = [p for p in out if trim_tier_map.get(p) not in skip]
    if only:
        out = [p for p in out if trim_tier_map.get(p) in only]
    return out


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


def run_one_person(pid, args, enable_self_align):
    """`enable_self_align` is passed explicitly per call (not read from
    args.enable_self_align_fallback) so a phase can control it independently of the CLI flag --
    the default two-phase mode always passes False for phase 1 and True for phase 2, regardless of
    what --enable-self-align-fallback was set to (that flag only matters in --single-phase mode)."""
    cmd = [
        "pixi", "run", "--manifest-path", os.path.join(REPO_ROOT, "pixi.toml"), "-e", "specimmune", "--",
        "python3", os.path.join(REPO_ROOT, "scripts", "run_immuannot_person.py"), str(pid),
        "--mount", args.mount, "--outroot", args.outroot,
        "--immuannot-dir", args.immuannot_dir, "--refdir", args.refdir,
        "--threads", str(args.threads_per_person),
        "--region", args.region, "--pad", str(args.pad),
        "--out-suffix", f".{pid}", "--force",
    ]
    if enable_self_align:
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


def run_phase(phase_label, people, ancestry_map, args, monitor_state_file, enable_self_align,
              auth_token, run_start_ts):
    """Runs one phase (a bounded list of people) to completion: resumability scan, heartbeat loop
    (own state file, own people_total -- so rate/ETA are phase-local, not polluted by any other
    phase that ran before it in the same process), concurrent dispatch, periodic + final fragment
    merge. Returns (n_done, n_given_up, wall_seconds). Mount/lock are the CALLER's responsibility
    (held once for the whole multi-phase run, not per phase).

    `run_start_ts` (2026-08-05) is the timestamp the WHOLE run began (before phase 1), used ONLY to
    compute a COMBINED cost-vs-budget check across every phase run so far in this process -- kept
    deliberately separate from the phase-local heartbeat state (monitor_state_file), which is what
    you actually watch on the dashboard to judge "is THIS phase worth continuing." The real $300
    cap (BRIEF.md) is a total-spend number, not a per-phase one -- since it's the same VM running
    continuously across both phases, real dollars are proportional to total wall-clock time since
    the orchestrator started, not phase-local elapsed hours."""
    print(f"\n########## PHASE: {phase_label} ({len(people)} people) ##########", file=sys.stderr)
    print("Scanning for already-completed people (real .gtf.gz existence + real-call check, "
          "not a log entry -- this is what makes killing and relaunching safe)...", file=sys.stderr)
    done_set, gave_up_set = scan_already_done(args.outroot, people)
    worklist = [p for p in people if p not in done_set and p not in gave_up_set]
    print(f"  {len(done_set)} already done (skipped), {len(gave_up_set)} previously given up "
          f"after {args.max_attempts}+ attempts (skipped -- see .orchestrator_gave_up markers "
          f"for review), {len(worklist)} to process this phase.", file=sys.stderr)

    # Per-ancestry progress (2026-08-04) -- LOCAL log only, never sent to the remote heartbeat
    # dashboard (monitoring/README.md keeps that aggregate-counts-only by design). In-memory
    # counters only, never re-derived by scanning the output directory. Motivation: the
    # self_align_needed phase (sequel2) is ~95% AFR, so a failure mode that disproportionately
    # hits one ancestry group should be visible live, not discovered after the fact.
    anc_total = Counter(ancestry_map.get(p, "NA") for p in people)
    anc_done = Counter(ancestry_map.get(p, "NA") for p in done_set)
    print(f"  [{phase_label}] ancestry distribution: " +
          ", ".join(f"{a}={n}" for a, n in anc_total.most_common()), file=sys.stderr)

    if not worklist:
        print(f"  [{phase_label}] nothing left to do -- every person already done or given-up.",
              file=sys.stderr)
        merge_fragments(args.outroot)
        return len(done_set), len(gave_up_set), 0.0

    state = {"done": len(done_set), "failed": len(gave_up_set), "lock": threading.Lock(),
             "anc_done": anc_done}
    people_total = len(people)
    stop_event = threading.Event()

    def heartbeat_loop():
        while not stop_event.is_set():
          # Belt-and-braces: ANY unexpected error in here must not kill the thread. Losing one
          # beat is a blip; losing the thread means the dashboard goes dark for the rest of a
          # multi-day run while the pipeline carries on looking healthy locally (exactly what
          # happened on the 2026-08-05 smoke test). Monitoring must never be more fragile than
          # the thing it monitors.
          try:
            with state["lock"]:
                d, f = state["done"], state["failed"]
            disk_pct = disk_used_pct(args.outroot)
            mem_pct = mem_avail_pct()
            heartbeat_kwargs = dict(
                receiver_url=args.monitor_url, auth_token=auth_token,
                vm_hourly_rate_usd=args.vm_rate, budget_usd=args.budget,
                run_state="running", phase=phase_label,
            )
            # Only pass state_path when we actually have one -- passing None overrides
            # heartbeat_client's own STATE_PATH_DEFAULT and blows up in os.path.exists(None).
            # (Real bug, hit on the 2026-08-05 smoke test: it killed the heartbeat THREAD on its
            # first tick while the pipeline itself kept running, i.e. a silent loss of all live
            # monitoring for the rest of the run -- the exact failure mode the monitoring exists
            # to prevent. The startup/terminal beats survived because send_lifecycle_heartbeat
            # guards this correctly; only the periodic loop was wrong.)
            if monitor_state_file:
                heartbeat_kwargs["state_path"] = monitor_state_file
            status, payload = send_heartbeat(d, f, people_total, disk_pct, mem_pct,
                                              **heartbeat_kwargs)
            phase_cost = payload.get("cost_so_far_usd")  # phase-local -- for judging THIS phase

            # Combined cost across every phase run so far in THIS process -- the real number to
            # check against the true $300 total cap, independent of the phase-local dashboard state.
            global_elapsed_hours = (time.time() - run_start_ts) / 3600.0
            global_cost = (global_elapsed_hours * args.vm_rate) if args.vm_rate else None
            if global_cost is not None and args.budget and global_cost > args.budget:
                warn(f"[{phase_label}] COMBINED cost_so_far_usd={global_cost:.2f} (all phases, "
                     f"this run) has EXCEEDED --budget={args.budget:.2f} -- this is a flag, not an "
                     f"auto-stop. Decide whether to let it keep running.")
            if payload.get("anomaly"):
                warn(f"[{phase_label}] heartbeat anomaly: {payload.get('anomaly_reason')}")
            print(f"  [{phase_label}] [heartbeat] status={status} done={d} failed={f}/{people_total} "
                  f"rate={payload.get('rate_per_hour')}/hr eta_hr={payload.get('eta_hours_remaining')} "
                  f"phase_cost=${phase_cost} combined_cost=${global_cost}", file=sys.stderr)
            with state["lock"]:
                anc_done_snapshot = dict(state["anc_done"])
            print(f"  [{phase_label}] [heartbeat] per-ancestry done/total (local log only): " +
                  ", ".join(f"{a}={anc_done_snapshot.get(a, 0)}/{n}"
                            for a, n in anc_total.most_common()), file=sys.stderr)
          except Exception as e:
            warn(f"[{phase_label}] heartbeat tick failed ({e!r}) -- thread stays alive, will "
                 f"retry next interval.")
          stop_event.wait(args.heartbeat_interval_sec)

    hb_thread = threading.Thread(target=heartbeat_loop, daemon=True)
    hb_thread.start()

    total_cores = args.concurrency * args.threads_per_person
    print(f"=== [{phase_label}] Launching: concurrency={args.concurrency}, threads_per_person="
          f"{args.threads_per_person} (total_cores={total_cores}), {len(worklist)} people, "
          f"self_align={'ON' if enable_self_align else 'off'} ===", file=sys.stderr)

    n_since_merge = 0
    t0 = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.concurrency) as pool:
        futures = {pool.submit(run_one_person, pid, args, enable_self_align): pid
                   for pid in worklist}
        for i, fut in enumerate(concurrent.futures.as_completed(futures), 1):
            pid = futures[fut]
            try:
                _, rc, elapsed = fut.result()
            except Exception as e:  # a worker itself crashing is a real bug -- surface it loudly
                warn(f"[{phase_label}] person {pid}: worker raised {e!r} (not a subprocess "
                     f"failure -- a real orchestrator-side exception).")
                rc, elapsed = -1, 0.0

            ok = person_done(args.outroot, pid)  # never trust exit code alone
            if ok:
                with state["lock"]:
                    state["done"] += 1
                    state["anc_done"][ancestry_map.get(pid, "NA")] += 1
                print(f"  [{phase_label}] [{i}/{len(worklist)}] {pid}: DONE ({elapsed:.0f}s)",
                      file=sys.stderr)
            else:
                n_attempts = bump_attempts(args.outroot, pid)
                if n_attempts >= args.max_attempts:
                    os.makedirs(os.path.dirname(gave_up_path(args.outroot, pid)), exist_ok=True)
                    with open(gave_up_path(args.outroot, pid), "w") as f:
                        f.write(f"gave up after {n_attempts} attempts, exit={rc}\n")
                    with state["lock"]:
                        state["failed"] += 1
                    warn(f"[{phase_label}] person {pid}: GAVE UP after {n_attempts} attempts "
                         f"(exit={rc}, {elapsed:.0f}s) -- see .orchestrator_gave_up / "
                         f".orchestrator_last_run.log")
                else:
                    print(f"  [{phase_label}] [{i}/{len(worklist)}] {pid}: FAILED (exit={rc}, "
                          f"{elapsed:.0f}s, attempt {n_attempts}/{args.max_attempts} -- will "
                          f"retry on next relaunch)", file=sys.stderr)

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
    print(f"=== [{phase_label}] DONE: {d}/{people_total} done, {f}/{people_total} given-up, "
          f"wall={wall/3600:.2f}h ===", file=sys.stderr)
    return d, f, wall


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
    ap.add_argument("--single-phase", action="store_true",
                    help="Disable the default automatic two-phase behavior (see module docstring) "
                         "and just run ONE pass over --cohort, filtered by --skip-trim-tier/"
                         "--only-trim-tier if given, using --enable-self-align-fallback as passed. "
                         "For the real production launch, do NOT pass this -- the default "
                         "(everyone except self_align_needed, then automatically continue straight "
                         "into self_align_needed with the fallback enabled, same process, same "
                         "lock, no manual second command) is what Marc asked for 2026-08-05: launch "
                         "once, both phases run without any further intervention.")
    ap.add_argument("--enable-self-align-fallback", action="store_true",
                    help="--single-phase mode only: attempts Tier 3 (self-align) for people whose "
                         "only resolving trim_tier is self_align_needed. In default two-phase mode "
                         "this is ignored -- phase 2 always enables it automatically, phase 1 never "
                         "does (it never has self_align_needed people to begin with).")
    ap.add_argument("--skip-trim-tier", action="append", default=[],
                    help="--single-phase mode only: exclude a trim_tier (repeatable). Dies with a "
                         "clear message if passed without --single-phase, to avoid silently "
                         "reinterpreting the default two-phase split.")
    ap.add_argument("--only-trim-tier", action="append", default=[],
                    help="--single-phase mode only: restrict to JUST the given trim_tier(s) "
                         "(repeatable) -- the inverse of --skip-trim-tier.")
    ap.add_argument("--max-attempts", type=int, default=3,
                    help="After this many failed attempts, a person is marked given-up (not "
                         "retried on future relaunches) instead of consuming a worker slot on a "
                         "deterministically-failing person for the rest of a multi-day run.")
    ap.add_argument("--heartbeat-interval-sec", type=float, default=120,
                    help="Default 120s (2 min), revised 2026-08-05 from 5 min: overhead is "
                         "negligible at any interval (a few counters and one small POST), and a "
                         "tighter cadence is worth real money in the first watched minutes of a "
                         "~$3.55/h run. Do not go below ~60s -- past that, human reaction time is "
                         "the bottleneck, not detection time. Keep MONITOR_STALE_MINUTES on the "
                         "receiver at ~3x this (6 min) so a single dropped beat isn't an alarm.")
    ap.add_argument("--monitor-url", default=os.environ.get("MONITOR_URL", DEFAULT_MONITOR_URL))
    ap.add_argument("--monitor-state-file", default=None,
                    help="Local state file heartbeat_client uses for phase 1 (or the single run, "
                         "under --single-phase). Default: heartbeat_client's own client_state.json. "
                         "Phase 2 (in default two-phase mode) ALWAYS gets its own separate, "
                         "auto-derived state file (<outroot>/monitor_state_phase2.json) regardless "
                         "of this flag -- automatic, no manual step needed -- so its rate/ETA/cost "
                         "are never polluted by phase 1's already-elapsed hours.")
    ap.add_argument("--vm-rate", type=float, default=None,
                    help="REAL Workbench-UI-confirmed USD/hour for the running VM -- required for "
                         "cost_so_far_usd / budget tracking to mean anything. DECISIONS.md's "
                         "~$3.03/hr is a research estimate, not a quoted price -- confirm live.")
    ap.add_argument("--budget", type=float, default=300.0,
                    help="Total spend cap (BRIEF.md's real $300 cap) -- tracked COMBINED across "
                         "both phases (elapsed wall-clock time since this process started x "
                         "--vm-rate), not reset per phase. A separate phase-local cost is also "
                         "shown in each heartbeat log line (from the dashboard state), for judging "
                         "an individual phase's own cost -- but this flag's warning is the real "
                         "total-spend check.")
    args = ap.parse_args()

    total_cores = args.concurrency * args.threads_per_person
    n_cpu = os.cpu_count() or 0
    if total_cores > n_cpu:
        die(f"concurrency({args.concurrency}) * threads_per_person({args.threads_per_person}) = "
            f"{total_cores} cores, but this machine only has {n_cpu}. Refusing to oversubscribe.")

    if (args.skip_trim_tier or args.only_trim_tier) and not args.single_phase:
        die("--skip-trim-tier/--only-trim-tier only apply under --single-phase -- without it, "
            "the default automatic two-phase split (everyone except self_align_needed, then "
            "self_align_needed with the fallback enabled) would silently ignore them. Add "
            "--single-phase if you want manual tier filtering instead of the automatic two-phase "
            "run.")

    check_mount(args.mount)
    os.makedirs(args.outroot, exist_ok=True)
    lock_path = os.path.join(args.outroot, "production_orchestrator.lock")
    check_and_take_lock(lock_path)  # ONE lock for the whole run, held across both phases

    auth_token = os.environ.get("MONITOR_AUTH_TOKEN")
    if not auth_token:
        warn("MONITOR_AUTH_TOKEN not set -- heartbeats will fail every ~5 min for the whole run "
             "(non-fatal to the pipeline, but you will be flying blind on an unattended multi-day "
             "job). Set it in the environment before launching for real.")
    if args.vm_rate is None:
        warn("--vm-rate not passed -- cost_so_far_usd will be null in every heartbeat and the "
             "local budget check below is disabled. Confirm the real n2-highcpu-96 rate in the "
             "Workbench UI and pass it.")

    run_start_ts = time.time()  # combined cost-vs-budget tracking spans every phase from here

    run_ok = False
    people_total_all = 0
    totals = {"done": 0, "failed": 0}
    try:
        people, ancestry_map, trim_tier_map = load_cohort(args.cohort)
        people_total_all = len(people)
        print(f"Cohort: {len(people)} people loaded from {args.cohort}.", file=sys.stderr)

        # Immediate startup beat -- confirms the VM->Hetzner path and arms the stale watchdog
        # before any long-running work begins (see send_lifecycle_heartbeat's docstring).
        if auth_token:
            st = send_lifecycle_heartbeat(args, auth_token, "starting", 0, 0, len(people),
                                           phase="startup",
                                           state_path=args.monitor_state_file)
            if st == 200:
                print("  Monitoring path confirmed: the dashboard should now read STARTING. "
                      "If it does, this VM can reach the Hetzner box and you're safe to walk away "
                      "once the run settles.", file=sys.stderr)
            else:
                warn("startup heartbeat did NOT get a 200 -- the dashboard will not update and "
                     "you will have NO remote visibility for this run. Most likely cause on a "
                     "Workbench VM is the VPC-SC perimeter blocking egress to the Hetzner box "
                     "(see production_orchestrator/PREFLIGHT.md item A). The pipeline itself is "
                     "unaffected and will keep running -- but decide NOW whether to launch blind.")

        if args.single_phase:
            phase_people = filter_people(people, trim_tier_map, args.skip_trim_tier,
                                          args.only_trim_tier)
            d, f, _ = run_phase("single", phase_people, ancestry_map, args,
                                 args.monitor_state_file, args.enable_self_align_fallback,
                                 auth_token, run_start_ts)
            totals["done"], totals["failed"] = d, f
            people_total_all = len(phase_people)
            run_ok = True
            return

        # --- Default: automatic two-phase run, no manual intervention between phases ---
        phase1_people = filter_people(people, trim_tier_map, skip_tiers=[SEQUEL2_TIER])
        phase2_people = filter_people(people, trim_tier_map, only_tiers=[SEQUEL2_TIER])
        print(f"Automatic two-phase run: phase 1 = {len(phase1_people)} people (everyone except "
              f"{SEQUEL2_TIER}), phase 2 = {len(phase2_people)} people (only {SEQUEL2_TIER}, "
              f"self-align fallback auto-enabled). Phase 2 starts automatically the moment phase 1 "
              f"finishes -- no further command needed. Cost is tracked COMBINED across both phases "
              f"against --budget={args.budget} (not reset per phase).", file=sys.stderr)

        if phase2_people:
            problem = check_phase2_prereqs()
            if problem:
                warn(f"PHASE 2 PREREQUISITE MISSING -- {problem}")

        d1, f1, _ = run_phase("phase1_main_cohort", phase1_people, ancestry_map, args,
                               args.monitor_state_file, enable_self_align=False,
                               auth_token=auth_token, run_start_ts=run_start_ts)
        totals["done"], totals["failed"] = d1, f1

        if not phase2_people:
            print("No phase-2 (self_align_needed) people in this cohort -- run complete.",
                  file=sys.stderr)
            run_ok = True
            return

        # Re-verify the mount before phase 2 -- phase 1 can run for many hours, and while the
        # orchestrator process staying alive across that whole span means no VM restart happened
        # (quirk #14: a restart implies whatever was running already exited), a defensive re-check
        # here is cheap and catches the rarer case of the gcsfuse process itself dying independently
        # of the VM/orchestrator.
        check_mount(args.mount)
        phase2_state_file = os.path.join(args.outroot, "monitor_state_phase2.json")
        print(f"\nAuto-continuing into phase 2 ({len(phase2_people)} people, self-align fallback "
              f"ON, isolated monitor state at {phase2_state_file}) -- no manual step needed.",
              file=sys.stderr)
        d2, f2, _ = run_phase("phase2_self_align", phase2_people, ancestry_map, args,
                               phase2_state_file, enable_self_align=True, auth_token=auth_token,
                               run_start_ts=run_start_ts)
        totals["done"], totals["failed"] = d1 + d2, f1 + f2
        run_ok = True
    finally:
        # Terminal heartbeat on EVERY exit path, including Ctrl-C and an unhandled exception.
        # This is what turns "heartbeats stopped" from ambiguous into a definite answer, and it
        # disarms the receiver's stale watchdog so a finished run doesn't page overnight.
        if auth_token:
            send_lifecycle_heartbeat(
                args, auth_token, "complete" if run_ok else "failed",
                totals["done"], totals["failed"], people_total_all or 1,
                phase="finished" if run_ok else "stopped-early",
                state_path=args.monitor_state_file)
        if os.path.exists(lock_path):
            os.remove(lock_path)


if __name__ == "__main__":
    main()
