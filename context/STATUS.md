# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-10-27 — all sprint branches consolidated into `main`

`main` now holds S05 (incl. 65b/65c), S06 call quality (70-79), S08 case browser (90-96) and S09 discovery loop (100-114).
Merged: `s05-65b-results`, `s06-call-quality`, `s08-case-browser`, `s09-discovery-loop`, `s05-oct06-work` (only its call-note untracking;
main's S05 text kept). Not carried over: `s05-65b-results`' early 67/68/69 scripts (superseded by main's 67/68, as noted in the S05 STATUS).
Aleix's `aleix/hla-resolve-phase1` stays separate. Pipeline table: `scripts/hla_popgen/README.md`; sprint index: `sprints/README.md`.
Open (human): call-note files remain in git history (`sprints/_calls/`, commits before 13e25ad) — history rewrite is a Marc decision, see DECISIONS.md.
Script 80 (case browser, S06 copy) was folded into 90 (S08); the 79 final report builder is in `scripts/hla_popgen/79_s06_final_report.py`.

## Pending VM run (2026-10-29): novel-allele ladder

Cloud session built `118_novel_ladder.py` + `119_novel_ladder_report.py` (selector, saturation curves, IMGT section) but had no VM. Newest IPD = 3.66.0 (14 Oct 2026).
Next commands: `sprints/S08_case_browser/RUNBOOK_novel_ladder.md` (extract with a large `--exhaustive-snv` cap, `90 sr`, then 118, then 119). Nothing in the report has real cohort-wide numbers until that runs.
The WIP desktop commit 8e2febc (S09_REPORT.html edits + README rows pointing to nonexistent 116/117 scripts) was un-committed on purpose; its patch is recoverable from the object 8e2febc.

## Where to pick up (2026-10-27)

- **S06** (call quality, 70–79): essentially complete. Open: sequence-level CDS comparison for KIR/HLA labels; haplotype-level yak origin check. Log: `sprints/S06_call_quality_metrics/LOG.md`; final report is a local build (`python3 scripts/hla_popgen/79_s06_final_report.py --tables reports/hla_popgen/79_s06_final/aggregate_tables --out reports/hla_popgen/79_s06_final`).
- **S05** closed 2026-10-13; open decisions (ancestry-purity rule, protein-novel demotion, KIR complete-only, physical-haplotype rule, call-level depth flag, TRB mixture genotyper, Sequel II re-call, disclosure residue) are in `DECISIONS.md`. Sprint rule: no VM job longer than one day. VM state: assume restarted; remount the bucket (ENVIRONMENT quirk #14).
- **S08/S09**: logs in `sprints/S08_case_browser/` and `sprints/S09_discovery_loop/`; S09 red-team checklist (`sprints/S09_discovery_loop/README.md`) applies to every new experiment.
- **Housekeeping for a human** (classifier-blocked for Claude): delete merged remote branches (`git push origin --delete s04-kir-recurrence-style-share s05-58-results s05-59-65-results s05-61-62-results s05-64-results s05-vm-pulls s05-oct06-work s05-65b-results s06-call-quality s09-discovery-loop needle-view-cds-diversity-density claude/project-thread-a3l14r`); push or drop the unpushed walkthrough commit on local `s04-kir-recurrence-style-share`.
- Layout after the 2026-10-27 cleanup: `reports/pre_sprint/` (July–Aug reports), `reports/hla_popgen/` (numbered, with by-era table), `scripts/README.md`, `sprints/*/{critics,ops}/`.

