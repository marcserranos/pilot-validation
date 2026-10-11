# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-10-11 — All of Us access restored; repo harmonized; planning before Sprint S05

**Access.** Marc lost Workbench access ~2026-09-28 (workspace showed Standard Tier only); UCSC
reinstated it, confirmed restored 2026-10-11. Assume the VM restarted and its disk state is
whatever it was on 2026-09-28 (quirk #14): remount, `ls`, verify before trusting anything.

**Branches.** PR #1 merged 2026-10-11: `main` now holds S03 + S04 and the harmonization pass and
is the canonical branch. `needle-view-cds-diversity-density`
and `s04-kir-recurrence-style-share` are fully contained in `main` and can be deleted.
Aleix's `aleix/hla-resolve-phase1` shares no history with `main` (separate root); leave it alone.

**Where things are.**
- Index of every result folder, current vs superseded: [`reports/hla_popgen/README.md`](../reports/hla_popgen/README.md).
- Last sprint (S04) plain-language findings: [`sprints/S04_kir_recurrence_style_share/FINDINGS_FOR_MARC.md`](../sprints/S04_kir_recurrence_style_share/FINDINGS_FOR_MARC.md); figure catalogue: `FIGURES_INDEX.md` in the same folder.
- Last call: [`sprints/_calls/2026-10-06_CALL_SUMMARY.md`](../sprints/_calls/2026-10-06_CALL_SUMMARY.md) (from Marc's pasted transcript). Marc's agreed focus: QC and phasing (incl. regional MHC/KIR coverage filter) and KIR; later, disease curation.
- Aleix's live state: `aleix/RNA-seq/context/STATUS.md` on his branch (W0 analysis table and W1 V/J usage done 10-08; TRBV4-3 germline-absence lead).

**Done since 2026-09-28 (S04 close-out, before access was lost):** 44 v4b full run (KIR/HLA
recurrence + catalogue completeness, valid numbers in 44/46 READMEs); 4 of 9 PUBLIC tables of 44 pulled.

**Open, needs Marc:**
1. Run `sprints/S04_kir_recurrence_style_share/VM_CLEANUP_PLAN_v2.md` (agent deletes are classifier-blocked).

**Literal next commands (VM, once Marc steers):**
1. Reattach the big_run app, remount gcsfuse, `ls` the bucket and `~/s04/results/44_public/`.
2. Pull the remaining 5 of 9 `44_public` tables, commit them under `reports/hla_popgen/44_kir_recurrence_saturation/`.
3. Then start the S05 board from the plan Marc approves.
