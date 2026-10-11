# S04 VM runs (config, times, cost)

App: AoU_Jupyter_ComputeEngine_20260805_big_run (workspace full-cohort-hla-calling, us-central1-a).
Config throughout S04: n2-highmem-4, 2000 GB disk, autostop 1 h, no resizes. Stopped rate shown
$0.11/h (disk); running $0.37/h.

**All durations/costs below are estimates reconstructed from `LOG.md`'s timestamped entries, not a
reconciled billing record** — no session logged a precise stop time at the point of stopping.

| Session | Start | End | Config | Jobs | Duration (est.) | Cost est. |
|---|---|---|---|---|---|---|
| 1 | 09-25 ~21:10 | 09-26 early morning (ended by an account rate-limit hit, ~00:00) | n2-highmem-4 | 47 v1, 48, 44 v1 | ~2h50m | ~$1.05 |
| 2 | 09-26 ~01:10 (restart) | ~09-26 10:00 | n2-highmem-4 | 44 v1 pulls/debugging; a second rate-limit hit landed within this window (~05:45 per LOG), so the VM likely auto-stopped for part of it — this row is an upper-bound wall-clock window, not continuous billed time | up to ~8h50m | up to ~$3.27 |
| 3 | 09-27 (multiple restarts across the day) | 09-27 (v4b complete) | n2-highmem-4 | 44 v2 (608s), 47 v2 (1,837s), 44 v3 (292s), 44 v4 smoke test (132s), 44 v4b (one run killed by an app stop mid-run, then a successful 2,742s / ~46min rerun after adding per-checkpoint saving) | job time sums to ~1h35m; wall-clock across restarts/debugging gaps est. ~3h | ~$0.59 (job time) to ~$1.11 (incl. gaps) |

**Rough total: ~13-15 VM-hours across the sprint, roughly $5-6 at $0.37/h.** Treat as a coarse
upper-bound estimate, not a bill — see `FINDINGS_FOR_MARC.md` §5 for the same caveat. A future
sprint should record real start/stop timestamps per session if a firmer number is ever needed.
