# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-08-04 (cont.) — orchestrator + cohort builder built, not yet run; sequel2 decision blocked on one VM test

Both build tasks from `scripts/production_orchestrator/BRIEF.md` are now written (this repo, this
session, not yet committed/pushed — confirm with Marc before pushing to the public GitHub repo):

- **`scripts/production_orchestrator/run_production_orchestrator.py`** — extends `scaling_probe.py`.
  Mount check, PID lock, real resumability (checks each person's actual `hap{1,2}.gtf.gz`, not a
  log entry), heartbeat every ~5 min from in-memory counters only, periodic fragment-merge into
  canonical `immuannot_calls.tsv`/`immuannot_timing.tsv`, per-person attempt cap (default 3) with a
  `.orchestrator_gave_up` marker instead of retrying a deterministically-failing person forever,
  local budget-vs-cost warning. **Not yet run on any VM — needs the pre-launch smoke test
  (BRIEF.md "Before the real launch") before the real n2-highcpu-96 launch.**
- **`scripts/build_immuannot_cohort.py`** — rewritten for the full cohort. No default `-n` cap,
  platform filter widened to revio+sequel2e+sequel2, row-level existence check + person-level union
  across a person's multiple manifest rows (quirk #13 discipline), parallelized+checkpointed
  (`lr_manifest_format_census.py`'s pattern). Outputs `person_id, platform, trim_tier, n_rows` —
  `trim_tier` is `paf_region` / `bam_whole_contig` / `self_align_needed` (the sequel2 gap) / person
  excluded entirely if no assembly FASTA at all. **Not yet run for real either — needs the mount.**
- **`scripts/production_orchestrator/RESULTS_LOCATION.md`** — deliverable #3, written.

### sequel2 (991 people, ~95% AFR): fallback tier BUILT, explicitly NOT validated — blocked on Marc

Raised explicitly with Marc per BRIEF.md's instruction not to silently default either way.
Marc's read: worried the "whole assembly" framing in BRIEF.md could mean ~100x-plus blowup: correct
— literally feeding Immuannot the untrimmed ~3.1Gb/hap assembly (no way to know which contig is
chr6 without an aln-to-hg38 file) is more like **~700-750x** the normal ~4.2Mb trimmed input, not
~100x. **Built a cheaper Tier 3 instead of the brief's literal "whole assembly" framing**: self-align
each hap FASTA against a chr6-ONLY reference slice (cached once from the full hg38 ref already on
the VM) with `minimap2`, producing our own synthetic `.paf` that the existing `regions_from_paf()`
parses unmodified — turns "whole assembly" into "one extra fast minimap2 pass + a normal-sized
trim." Implemented in `run_immuannot_person.py` (`--enable-self-align-fallback`, opt-in, off by
default) — **but genuinely untested**, no VM access from this session to actually time it.

**Decided with Marc (2026-08-04):** don't guess the number — smoke-test it for real on 1-2 actual
sequel2 people, on the bigger Stanford-pod test VM, before deciding whether sequel2 is in or out of
the production launch. **I don't have gcloud/SSH access to that VM from this session — Marc runs
the commands below himself and pastes the output back.**

```bash
# From ~/repos/pilot-validation on the Stanford-pod VM, mount up, inside pixi shell -e specimmune
# (or `pixi run -e specimmune --` prefix). Pick 2 real sequel2 person_ids first:
python3 -c "
import pandas as pd
lr = pd.read_csv('~/mnt/aou-controlled/v9/wgs/long_read/manifest.tsv', sep='\t', dtype=str)
cand = lr[(lr.platform == 'sequel2') & lr.assembly_hap1_fa.notna() & lr.assembly_hap2_fa.notna()]
print(cand['research_id'].drop_duplicates().head(2).tolist())
"
# Then, for each of those 2 ids (replace <PID>):
{ time python3 scripts/run_immuannot_person.py <PID> --enable-self-align-fallback \
    --out-suffix .selfalign_test --force ; } 2> ~/pipeline_outputs/selfalign_test_<PID>.log
```

Paste back: the `real`/`user`/`sys` time line, and the `self_align_seconds` /
`whole_contig_mb` / `padded_mb` / `trim_seconds` / `immuannot_seconds` columns from
`~/pipeline_outputs/immuannot_timing.selfalign_test.tsv` for both people. That's what turns "I
think it's fine" into a real number to decide on (BRIEF.md: "$10 vs $100 vs an extra day" — Marc's
framing).

**2026-08-05 — bench test aborted, replaced with a live two-phase launch — now fully automatic in
ONE command (Marc: "I don't want to intervene... make everything on the same run").** Didn't want
to spend unknown debugging time on an untested ~20min-2hr bench test with no time budget for it,
and didn't want to have to manually type a second command to start phase 2 either. **Current
behavior, default, no flags needed:** a single `run_production_orchestrator.py` invocation runs
phase 1 (everyone except `self_align_needed`, ~12,261 people, normal proven trim path) to
completion, then AUTOMATICALLY continues into phase 2 (only `self_align_needed`, the 991 sequel2
people, self-align fallback auto-enabled) — same process, same lock, no second command. Each phase
gets its own heartbeat state file (phase 2's is auto-derived at
`<outroot>/monitor_state_phase2.json`) so phase 2's rate/ETA/cost on the dashboard are never
diluted by phase 1's already-elapsed hours — that's the number to watch to decide whether to Ctrl-C
phase 2 (note: this waits for the current in-flight wave to finish, not an instant kill; a relaunch
with the same command resumes correctly regardless of where it was aborted).
`--single-phase` + `--skip-trim-tier`/`--only-trim-tier` preserves the old one-shot-filtered
behavior for the smoke test.

Exact launch command: see chat / `RESULTS_LOCATION.md` / the orchestrator's own `--help`, not
duplicated here (this file is for durable state, not copy-paste command logs).

## Pre-flight state as of 2026-08-05 — see `scripts/production_orchestrator/PREFLIGHT.md`

That file is the rolling launch checklist (what's verified, what's still blocking). Cleared this
session, on the Hetzner box directly: firewall/port 8943 confirmed reachable externally, receiver
redeployed with a server-side cost model using the **real quoted rates ($3.55/h VM + $81.60/mo
disk**, not DECISIONS.md's ~$3.03/hr research estimate), dashboard redesigned, stale threshold
25→15 min, and the **first-ever real heartbeat end-to-end test against the deployed box** (all
prior testing was loopback-only). Two real bugs found and fixed while doing it: `mem_avail_pct` was
inverted (sent *used*, consumers read it as *available* — would have fired a spurious anomaly push
early in the run and stayed silent during a real near-OOM), and the budget check was resetting
per-phase instead of tracking the real $300 total.

**The biggest remaining unknown is deliberately called out in PREFLIGHT.md item A: the Workbench
VM has never actually POSTed a heartbeat to the Hetzner box**, and that workspace sits behind a
VPC-SC perimeter that explicitly warns about egress to outside services. Proven-from-the-box is not
proven-from-the-VM. Test it before launching, not after.

## Carried forward from earlier 2026-08-04 — all pre-launch VM/config decisions locked

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

## What's actually left before the real launch

1. **Get the sequel2 self-align timing answer back from Marc** (commands above) and decide
   in/out/fast-follow for real, not by default.
2. **Commit + push this session's two new scripts + the `run_immuannot_person.py` Tier 3 addition**
   to the public GitHub repo (not done yet — confirm with Marc first, per this session's own
   git-safety discipline) so the VM can `git pull` them.
3. **Run the pre-launch smoke test** (BRIEF.md "Before the real launch, smoke test, don't skip
   this"): 24-48 people, one full concurrency wave, confirm resumability (kill mid-run, relaunch,
   confirm it skips done people), confirm the lock blocks a second launch, confirm the mount check
   blocks a launch with a dead mount, confirm a heartbeat actually reaches the Hetzner dashboard.
   **Not yet done — the orchestrator has never been run against real data, this session had no VM
   access to do it.**
4. **Confirm the live `n2-highcpu-96` hourly rate in the Workbench UI** when the VM is actually
   created, pass it via `--vm-rate` — DECISIONS.md's ~$3.03/hr is a research estimate.
5. **Confirm the Hetzner firewall rule for port 8943** (see below, carried from earlier).
6. Only then create the real VM and launch.

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
