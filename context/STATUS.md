# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-07-21 (cont.) — long-read census complete: the 2,763-person floor is obsolete

**Headline result, delivered:** the exhaustive long-read data-format census
(`scripts/lr_manifest_format_census.py`, checkpointed/resumable rewrite after a real 5h
data-loss incident — ENVIRONMENT.md quirk #22) ran clean on the full 14,521-person cohort.
**All 14,521 people now have a real, existing aligned GRCh38 BAM** — AoU backfilled the
previously-missing `revio`-platform BAMs sometime between 2026-07-11 and 2026-07-21, directly
confirmed by re-checking the exact case ENVIRONMENT.md quirk #13 had named as dangling (person
1008366's `revio` row). Full writeup, per-platform data-shape table, and caveats:
`reports/lr_data_census/README.md`. DECISIONS.md's "Where are the aligned BAMs" open question is
now Resolved. **The long-read validation cohort is no longer capped near ~2,763 — this changes
the practical scope of Experiment D's successor and the MID/SAS ancestry bottleneck.**

Conceptual side-question (SpecImmune's real input requirement, why assembly data isn't a simple
BAM-conversion target) already answered in DECISIONS.md's "Assembly-based HLA typing" bullet,
which also now carries a literature-research lead: **Immuannot**
(github.com/YingZhou001/Immuannot) looks like a ready-to-use tool that types HLA directly from
phased assembly contigs, reportedly already used by HPRC — not yet independently verified or
evaluated by us. ~86% of the cohort (revio/sequel2e/sequel2) has assembly data available as a
*second*, cross-validating product now, not a fallback for otherwise-unreachable people.

**Two real incidents this session, both now fixed and logged:** (1) the census script's original
version had no incremental output — see ENVIRONMENT.md quirk #22 for the fix (checkpointed,
resumable, parallelized) and its corrected root-cause account (the VM doesn't stop mid-job,
quirk #14 — the terminal's scrollback was lost after the job *finished*, and a later
`--sample-only` rerun separately clobbered the on-disk result via a shared-output-path bug, also
fixed). (2) A commit (`112372c`) mid-session bundled this session's in-progress edits with
unrelated concurrent work from Marc and Aleix, who both regularly work in this same repo
directly — **confirmed by Marc as expected, ongoing, normal working style, not an anomaly to
keep flagging.**

**Both carried-over items from 2026-07-19 were already resolved independently (found via git log,
not by this session's own investigation):** Aleix's `aleix/hla-resolve-phase1` workstream (4th
long-read caller, calibrated against real external ground truth) explains the earlier
"unexplained parallel VM experiments" note; `scripts/hla_disease_sanity_check.py`'s
disease-sanity-check thread (Celiac/narcolepsy) finished with real results independently of the
missing-companion-file concern.

## Pick up here

1. **Evaluate Immuannot** (DECISIONS.md, "Assembly-based HLA typing" bullet) as a second
   long-read caller for the ~86% of the cohort with assembly data — not yet started, no design
   work done beyond the literature lead.
2. **Re-scope the long-read validation cohort** now that all 14,521 people (not ~2,763) are
   viable — Experiment D's 60-person cohort and the MID/SAS ancestry bottleneck should be
   revisited against this much larger pool before assuming the old constraints still apply.
3. Read `aleix/README.md` directly for HLA-Resolve's current phase/status rather than
   re-deriving it from this note — that workstream is moving independently of this thread.

Not blocking otherwise. Next natural step: the strategic fork in DECISIONS.md (build downstream
directly on AoU-native vs. call independently), now informed by the cross-tool pilot, the callset
validation report, this census, and soon Aleix's ground-truth-anchored comparison.

## Watch / blockers

- **VM sleeps after ~1h idle (quirk #14)** — every new session starts with the mount gone; remount (quirk #11) and `ls`-verify before running anything that reads bucket data.
- Paste the `gcsfuse` mount and the consumer command *separately* — chaining them in one paste races (quirks #2, #14).
- **`~/ref/`, `~/repos/`, `~/tools/`, `~/pipeline_outputs/` survive a VM restart; the mount, background processes, and activated pixi shell do not** — includes `~/ref/imgt/` (no need to re-curl the IMGT reference files after a restart).
- **Any automated SpecImmune call must use `pixi run -e specimmune`** (quirk #15/#17) and verify its output file exists — never trust exit code 0 alone. Same discipline applies to SpecHLA (quirk #18) and to any QC metric.
- **`mkdir -p` the target directory before any `... > <dir>/file` redirect** (quirk #20).
- **pad100k (SpecImmune-LR) and pad2000-10000 (SpecHLA-SR) are two separate, tool-specific recommendations** — don't conflate them.
- **Gene-panel restriction is a closed question** (Experiment C) — don't re-attempt without a specific new reason.
- **The long-read BAM pool is now all 14,521 people, not the old 2,763 floor** (DECISIONS.md, resolved 2026-07-21) — don't cite the old number from memory; re-check DECISIONS.md if a stale "2,763" surfaces anywhere.
- **DRB1's SpecHLA-pad10000 padding-degradation hypothesis is still open, untested** (DECISIONS.md, 2026-07-12 entry).
