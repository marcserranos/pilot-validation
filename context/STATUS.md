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

**2026-08-05 — bench test aborted, replaced with a live two-phase launch instead.** Marc: didn't
want to spend unknown debugging time on an untested ~20min-2hr bench test with no time budget for
it. **New plan, decided with Marc:** launch the real production run in two sequential phases on the
same VM, watched live on the dashboard, abortable if phase 2 looks too slow/expensive:
- **Phase 1:** everyone except `self_align_needed` (`--skip-trim-tier self_align_needed`, no
  `--enable-self-align-fallback`) — the ~12,261 people with a normal, already-proven trim path.
- **Phase 2:** only `self_align_needed` (`--only-trim-tier self_align_needed
  --enable-self-align-fallback`, new flag pair added 2026-08-05) — the 991 sequel2 people, run
  right after phase 1 finishes, same VM. Watch the dashboard's rate/ETA for phase 2 specifically;
  if it implies days instead of hours or a real budget blowout, abort (Ctrl-C — note this waits for
  the current in-flight wave, up to `--concurrency` people, to finish naturally before exiting, not
  an instant kill) and reconsider rather than let it run unattended.
- **Required for phase 2's numbers to mean anything:** pass a FRESH `--monitor-state-file` (e.g.
  `~/pipeline_outputs/monitor_state_phase2.json`) for phase 2 — otherwise its rate/ETA get computed
  against phase 1's already-elapsed hours plus a much smaller people_done, giving a misleadingly low
  rate right when an honest number is needed most.

Exact two-phase launch commands: see chat, not duplicated here (this file is for durable state, not
copy-paste command logs — the orchestrator's own `--help` and `RESULTS_LOCATION.md` are the source
of truth for flags).

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
