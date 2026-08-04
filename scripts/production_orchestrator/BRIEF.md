# Full-cohort production orchestrator — implementation brief

> **Status: not started.** This is a locked-in-decisions + open-decisions brief, written 2026-08-04
> so a *separate* session/agent can build this without needing this session's full history.
> Self-contained on purpose — read this file, `context/DECISIONS.md`, `context/EXPERIMENTS.md`,
> and `context/ENVIRONMENT.md` (especially quirks #11/#14 and #23-#27) before writing any code.
> Read `README.md` at the repo root first if you haven't — it explains how `context/` works and
> why this brief lives under `scripts/` instead.

## What this is for

Run Immuannot allele calling across the full long-read-eligible All of Us cohort
(~13,000-14,521 people, exact number TBD by the cohort builder — see below), producing a
structured per-person, per-gene allele table for downstream non-linear PRS modeling (Celiac,
T1D, Psoriasis, MS). This is the actual production run this entire project has been building
toward. It costs real money (~$160-200 estimated) and runs for ~2-2.5 days — get it right before
launching, not after.

## Locked-in decisions (do not re-litigate without a real reason — see `context/DECISIONS.md` for full rationale)

| Decision | Value | Source |
|---|---|---|
| VM machine type | `n2-highcpu-96` (96 vCPU, 96GB RAM) | DECISIONS.md, 2026-08-04 |
| Region | `us-central1` | Matches the Stanford-pod workspace's region policy; already proven this session |
| Concurrency | 24 people at once, 4 threads each (24×4=96 cores) | DECISIONS.md — speed is flat across concurrency/thread splits at a fixed core budget, so use the split with the best memory headroom + existing track record |
| Data disk | 2TB (2000GB), **no per-person pruning** | ~81MB/person × 14,521 ≈ 1.15TB + margin; disk cost is trivial at this run length, pruning bugs are a real historical risk in this project |
| Low-confidence calls | **Call everyone. Keep every call, including low-confidence ones.** Record `template_distance`, `template_warning`, and the novel-allele `"new"` tag as metadata columns — never filter or drop a call for being low-confidence. | Decided earlier this session, explicit user instruction — do not silently reintroduce a confidence filter |
| Monitoring | POST a heartbeat to the already-deployed Hetzner receiver every ~15 min | `scripts/monitoring/README.md` "Wiring into the real orchestrator" section has the exact call shape; use `heartbeat_client.py`'s `send_heartbeat()` |
| Cost cap | **$300 total** | User's explicit budget; the 96-core config estimates ~$160-200, leaving real margin — don't casually blow past this without flagging it |

**GPU: confirmed not applicable — do not add it.** minimap2/Immuannot have no GPU code path.

## Required behavior — every one of these is a real lesson this project already paid for

Do not treat any of these as optional hardening — each corresponds to an actual multi-hour
incident this session. Full accounts: `context/ENVIRONMENT.md` quirks #11, #14, #23-#27.

1. **Hard mount check at startup, refuse to run without it.** The gcsfuse mount does not survive a
   VM restart/stop, ever. Check `os.path.isfile(<mount>/v9/wgs/long_read/manifest.tsv)` (or
   equivalent) before doing anything else; die loudly and immediately if missing. Do not let a
   dead mount cause a silent, confusing partial failure hours later — this happened, repeatedly.
2. **Single-instance lock (PID-file).** Two concurrent orchestrator instances against overlapping
   people silently corrupt shared per-person intermediate files — this happened, more than once,
   in one session, and cost hours to even notice. `scaling_probe.py` has a working reference
   implementation of this pattern (`check_and_take_lock()`) — reuse it, don't reinvent it.
3. **Real resumability.** Must be safe to kill and relaunch with the exact same command. Skip
   people already marked done (check for the expected output file, e.g. their `.gtf.gz` results,
   not just a log entry). `run_experiment_d.sh` (an older bash script, but the *resumability
   design* is sound) and `run_immuannot_person.py`'s own `--force`/already-done-skip logic are the
   reference patterns.
4. **Prefer Python over bash for the orchestrator itself.** This session tried a bash orchestrator
   first; it became a full day/night incident (duplicate-launch races, `du`/`df` hanging 15 min to
   2+ hours on this VM's gcsfuse mounts, silent continuation past failures because bash had no
   `set -e`). The replacement, `scripts/scaling_probe.py`, is a working, much simpler reference —
   real exceptions stop execution immediately with a traceback, no `set -e`/`pipefail` ambiguity.
   **This is the recommended starting skeleton for the orchestrator**, not a from-scratch design —
   it already has the mount check, the lock, and a working concurrent-worker pattern
   (`ThreadPoolExecutor` calling `run_immuannot_person.py` via `subprocess.run`). It just needs to
   be extended to (a) loop over the *entire* cohort instead of one config, (b) call the heartbeat
   client periodically, and (c) not delete/skip already-done people on restart.
5. **No `du`/`df` on the whole output tree or whole filesystem, ever.** If any disk-usage
   reporting is added, scope it to specific small directories only (quirk #25) — this caused the
   worst, hardest-to-diagnose stalls of the whole session.
6. **Every stage that touches a subprocess should verify real output exists, never trust exit code
   alone** (`run_immuannot_person.py` already does this correctly for the `.gtf.gz` — follow the
   same discipline for anything new).

## The cohort — what "everyone we can call" actually means, and what's still open

`scripts/build_immuannot_cohort.py` exists and does the right *kind* of check (existence-verifies
`assembly_hap1_fa`/`assembly_hap2_fa`/`assembly_hap1_aln2_hg38_bam`/`assembly_hap2_aln2_hg38_bam`
against the mount, not just a platform label) — but it was built for small verified test batches
(`-n`, defaulting to 40, stops as soon as N are found) and filters to `platform == "revio"` only.
For the full production run:

- **Remove or raise the `-n` cap** so it enumerates every verified-eligible person, not just the
  first N. Given `revio` is ~11,070 manifest rows, a serial per-row existence check (as the script
  currently does) may be slow — consider parallelizing (the pattern in
  `scripts/lr_manifest_format_census.py`, `ThreadPoolExecutor` since this is FUSE-I/O-bound, not
  CPU-bound, already checkpointed/resumable) rather than assuming a serial loop is fast enough at
  this scale. Verify by timing a real run, don't assume.
- **Open, deliberate decision needed, not a silent default: what about `sequel2`'s 991 people?**
  Per `reports/lr_data_census/README.md`, `sequel2` (the original AoU African-American-enriched
  long-read pilot cohort, ~95% AFR) has the assembly FASTA files but genuinely lacks the
  `aln2_hg38_bam`/`.paf` files `run_immuannot_person.py` currently requires for its trim step —
  `build_immuannot_cohort.py`'s `platform == "revio"` filter excludes them entirely, and
  `run_immuannot_person.py`'s own `process_person()` also hard-requires the aln columns even
  though Immuannot itself only actually needs the raw assembly FASTA (trimming is this project's
  own runtime optimization, not a real Immuannot requirement). **Including these 991 people would
  need a third fallback trim tier: whole-untrimmed-contig extraction when neither `.paf` nor
  aln-to-hg38 BAM exists** (untrimmed = slower per-person, whole assembly vs. a ~30MB region, but
  still processes correctly). This has never been built or tested. Given AFR representation
  matters for this project's ancestry aims, raise this explicitly with Marc rather than silently
  excluding or silently building the fallback — it's a real scope/cost trade-off (991 more
  people, but each slower and untested), not a bug to just fix.
- Also check `sequel2e` (1,219 rows) — per the census, same full-assembly-suite shape as `revio`,
  should already pass the existing filter if the platform check is widened from
  `== "revio"` to `.isin(["revio", "sequel2e"])`. Verify this against the mount directly, don't
  assume from the census table alone (same existence-check discipline as everything else).

## Interface with the monitoring subsystem (already built, already deployed)

`scripts/monitoring/heartbeat_client.py` has a working `send_heartbeat()` function and a
standalone CLI. Read `scripts/monitoring/README.md`'s "Wiring into the real orchestrator" section
for the exact call shape and required fields (`people_done`, `people_failed`, `people_total`,
`--vm-rate` — pass **the real confirmed `n2-highcpu-96` hourly rate from the Workbench UI when the
VM is created**, not the ~$3.03/hr estimate in DECISIONS.md, which is a research estimate, not a
quoted price). Call it roughly every 15 minutes from within the main processing loop — not more
often (this was designed as minimal sampling, not continuous telemetry, per the original ask).

## Before the real launch — smoke test, don't skip this

1. Run the orchestrator against a small batch (e.g. 24-48 people, one full concurrency wave) on
   whatever VM you're testing on, and confirm: real calls produced, resumability actually works
   (kill it mid-run, relaunch, confirm it skips the already-done people and finishes the rest),
   the lock actually blocks a second launch, the mount check actually blocks a launch with a dead
   mount, and a heartbeat actually reaches the Hetzner dashboard.
2. Only then create the real `n2-highcpu-96` VM and launch the full cohort.
3. Confirm the Hetzner firewall rule for port 8943 is in place (`scripts/monitoring/README.md` step
   2) before walking away from an unattended multi-day run — check with Marc if unsure, this needs
   his own Hetzner console access.

## Expected deliverables

1. The orchestrator script itself (Python, per the "prefer Python" guidance above), extending
   `scaling_probe.py`'s proven patterns rather than starting from zero.
2. An updated/extended cohort builder producing the real full-population `person_id` list, with
   the `sequel2` question above explicitly resolved (one way or the other, documented) rather than
   silently defaulted.
3. A short results-location note (where the real per-person calls end up, e.g.
   `~/pipeline_outputs/<person_id>/immuannot_output/`, and where the aggregate progress/summary
   lives) — remember the project's standing rule: real per-person allele calls are participant
   data, stay on the VM, never committed or downloaded (same as `SMOKE_TEST_PICKS.local.md` /
   `immuannot_calls.tsv` throughout this project).
4. `context/EXPERIMENTS.md`/`context/STATUS.md` updated with the real launch, once it happens —
   same discipline this whole project has followed throughout.
