# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-09-27 — Sprint S04 nearly complete; one VM run pending, rest ready for Marc

Board: [`sprints/S04_kir_recurrence_style_share/SPRINT.md`](../sprints/S04_kir_recurrence_style_share/SPRINT.md);
PI report: [`FINDINGS_FOR_MARC.md`](../sprints/S04_kir_recurrence_style_share/FINDINGS_FOR_MARC.md);
log: `LOG.md`. Branch `s04-kir-recurrence-style-share` (off S03 @ 5d5ef99), **committed, not pushed**
— it is a strict superset of `main`, `fig1-drafts-and-research-map`,
`needle-view-cds-diversity-density`, and `s03-call8-figures-kir-prediction`.

**Done:** WS-C (figure style + layout linter + redesign of Figure 1/DQ/39/43b), WS-D (47/48 naive
ML + KIR-HLA ligand co-occurrence, full cohort), WS-E (Cole package built + uploaded + verified to
`gs://hla-calls-share-wb-cordial-leechee-9743/release_2026-09-25/`), WS-F (repo reorg via `git mv`,
37 legacy scripts + 3 call notes + 2 briefs relocated).

**Pending — WS-A/WS-B (scripts 44/45/46), KIR recurrence & KIR-vs-HLA catalogue coverage:** v1-v3
numbers are INVALID (name-based "genomic identity" + missing KIR artifact filter — see
`context/DECISIONS.md` 2026-09-27 entries). Fix is committed (`e4211ad`); `44_v4` is deployed on
the VM and smoke-tested (300 people, clean) but the full-cohort run has not been started yet.

**Literal next commands:**
1. Marc: run the VM cleanup in [`VM_CLEANUP_PLAN_v2.md`](../sprints/S04_kir_recurrence_style_share/VM_CLEANUP_PLAN_v2.md) (agent deletes are classifier-blocked).
2. Next session: reattach the VM (app `AoU_Jupyter_ComputeEngine_..._big_run`, inner lab URL in
   `VM_OPERATOR_HANDOFF.md`), remount/verify per quirk #14, then run
   `44_kir_recurrence_saturation_v4.py` at full cohort scale (no `--limit`). Pull results, rerun
   45/46 against the new export, update both READMEs' numbers, update `FINDINGS_FOR_MARC.md`'s
   "pending 44 v4" section with real numbers.
3. Marc: push and merge the branch per `FINDINGS_FOR_MARC.md` §6 (item 2), then delete the 3
   subsumed branches.
4. Marc: decide the 3 open items in `FINDINGS_FOR_MARC.md` §6 (`outputs/` tracking, small-cell
   `s_obs` masking policy, Cole package README date re-upload).

## As of 2026-09-26 — Sprint S04 opened (not started)

Next orchestrator: read [`sprints/S04_kir_recurrence_style_share/ORCHESTRATOR_HANDOFF.md`](../sprints/S04_kir_recurrence_style_share/ORCHESTRATOR_HANDOFF.md)
first. Branch `s04-kir-recurrence-style-share` (off S03 @ 5d5ef99), not pushed. S03 state below.

## As of 2026-09-25 — Sprint S03 (call #8) complete; KIR called on the full cohort

Board: [`sprints/S03_call8_figures_kir_prediction/SPRINT.md`](../sprints/S03_call8_figures_kir_prediction/SPRINT.md);
PI report `FINDINGS_FOR_MARC.md`; run record `KIR_FULL_RUN.md`. Branch `s03-call8-figures-kir-prediction`,
committed, not pushed. Scripts 37–43 added (43 = KIR full-cohort aggregate + figure).
- Full-cohort KIR Immuannot run done (12,261 people, n2-highcpu-80, 31.1 h, ~$95). The per-person
  outputs stay on the VM under `~/pipeline_outputs_kir/`; the aggregate is in `reports/hla_popgen/43_kir_full_cohort/`.
- The big_run app is back on n2-highmem-4 / 2000 GB / autostop 1 h, and **Stopped**.
- The VM repo has an untracked `scripts/production_orchestrator/run_production_orchestrator_s03kir.py`
  (= this branch's orchestrator). It can be deleted once this branch is merged and the VM repo is updated.
- VM automation from this laptop worked this session via the Chrome tools + JupyterLab REST/websocket (VM_CHANNEL.md).

Next: (1) WS6 full 7,640 + XGBoost (VM, ~hours on 4 vCPU or brief on a resized app); (2) KIR
novel-allele recurrence check + 2DL2/2DL3 co-occurrence look (run 43-style on the VM against
`~/pipeline_outputs_kir`); (3) fix the Tier-3 self-align chr6-cache path for non-chr6 regions, then run
the 991 sequel2 people for KIR; (4) Marc pushes the branch.

## As of 2026-09-24 — Sprint S03 (call #8) complete except VM-only steps

Board: [`sprints/S03_call8_figures_kir_prediction/SPRINT.md`](../sprints/S03_call8_figures_kir_prediction/SPRINT.md);
PI report: `FINDINGS_FOR_MARC.md`; pending VM commands: `VM_HANDOFF.md` (Figure 1 panels a/c
exports, KIR extended-pilot aggregation, KIR full run after 96-core resize, WS6 full-cohort run).
Branch `s03-call8-figures-kir-prediction`, committed, not pushed. Scripts 37–42 added
(37/37b–e DQ G1/G2, 38/38b deletions, 39 saturation, 40/40a/40b Figure 1 v5, 41 KIR, 42 repertoire).
Agent VM automation from this laptop is currently denied by the auto-mode permission classifier;
run VM steps by hand or add a permission rule.

## As of 2026-09-17 — Sprints S01 and S02 both complete, branch not pushed

**Most recent sprint board:** [`sprints/S02_supervisor_items/SPRINT.md`](../sprints/S02_supervisor_items/SPRINT.md),
with the plain-language writeup in that folder's `FINDINGS_FOR_MARC.md`. S01's board and findings
sit alongside it. Read those, not this file, for tasks and results. Branch:
`fig1-drafts-and-research-map` (off `main` @ `30b14ac`), **committed but never pushed** — pushing
is blocked for the agent, so Marc must do it.

S02 turned the 2026-09-17 supervisor call into scripts 29-33, all run on the full cohort. Open
items are listed in that board's §5: Figure 1 panels a/b need scripts 06 and 10 rerun (committing
the frequency tables, not just figures), the structural peptide-contact transfer function is
unfinished, and the VM is showing "Reboot is required" (quirk #41).

Context since the last STATUS (2026-09-05): scripts 10–22 ran on the real cohort and were merged to
main (`30b14ac` renumbered them); a supervisor meeting (~2026-09-14) framed the work as a paper with
a 4-panel Figure 1 and 8 action items; the plan for those items and future research lines is in
`reports/hla_popgen/NEXT_STEPS_AND_RESEARCH_MAP.md`; `23_fig1_draft_panels.py` drafted panels from
committed aggregates and showed novelty is mostly non-coding — S01 re-derives that foundation.

Open from 2026-09-05, still true: RUNBOOK.md does not document that 05/07 need the gcsfuse mount +
explicit `--sr-genotypes ~/mnt/aou-controlled/v9/wgs/short_read/snpindel/aux/hla_variants/hla_genotypes.tsv`.
