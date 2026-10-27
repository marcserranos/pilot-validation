# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-10-24 ~00:50 — S06 call-quality sprint essentially complete (branch s06-call-quality)

Final report (local, git-excluded): `reports/hla_popgen/79_s06_final/S06_FINAL_REPORT.html` (build: `python3 scripts/hla_popgen/79_s06_final_report.py --tables reports/hla_popgen/79_s06_final/aggregate_tables --out reports/hla_popgen/79_s06_final`).
Log: `sprints/S06_call_quality_metrics/LOG.md`. VM idle; aggregates in `~/s06/results/*`, person-level only in `~/s06/internal`. Small-cell masking pass still open before anything is committed.
Open: sequence-level CDS comparison for KIR/HLA labels; haplotype-level yak origin check; think-tank ideas appended at the report end.

## As of 2026-10-13 ~13:30 UTC — S05 closed; report extended with Marc's 10-13 asks

**Report:** `sprints/S05_qc_phasing_kir_trb/S05_REPORT.html` (built from `report/src/*.html` by `report/build_report.py`),
published as the S05 Artifact. 10-13 additions: review folded into the text, data deep dives (access, coverage, same contig,
RNA-seq composition from 68), KIR in context (Part 8b), Part 4b "where errors happen, and pruning" (67), 65b/65c assembly
check in Part 9. Board/log: `sprints/S05_qc_phasing_kir_trb/{SPRINT,LOG,METHODS_LEDGER}.md`.

**Sprint rule (Marc, 2026-10-12): no job longer than one day.** Longer jobs are planned and costed in the report, not run.

**VM (big_run, n2-highmem-4)** is idle. 65b stopped at 913/2,000 people (Marc); 65c, 67 and 68 tables are in the repo,
transferred by message because VM git push and laptop pulls are blocked. VM-only extras: `~/s05/results/65c/65c_concordance.unrounded_VMonly.tsv`,
old pre-renumbering outputs `~/s05/results/{67_prev_composition_0934,68_prev_trioerr_0935}`, a stash with old 67/68 edits in
`~/repos/pilot-validation-main`. Branch `s05-65b-results` on GitHub: its 65/65c code is merged; its early 67/68/69 scripts are
superseded by main's 67/68 (do not merge them). Remount the bucket after any restart (ENVIRONMENT quirk #14).

**Open, needs Marc** (all in `DECISIONS.md` open questions): ancestry-purity rule; demote protein-novel counts; KIR
complete-only; physical-haplotype rule; no coverage filter / call-level depth flag (67); TRB mixture genotyper replacing the
cuts tables; Sequel II re-call (resize + cost); disclosure residue in git history and allele-count masking. Script 60 waits on
Aleix copying his TRBV usage table to the bucket.

**Housekeeping left for a human** (blocked for Claude by the permission classifier): delete the fully merged remote
branches `s04-kir-recurrence-style-share s05-58-results s05-59-65-results s05-61-62-results s05-64-results s05-vm-pulls`
(`git push origin --delete ...`); `s05-oct06-work` and `s05-65b-results` are superseded too (content on main). On the VM,
`git -C ~/repos/pilot-validation-main stash drop stash@{0}` (old 67/68 edits).

**Literal next steps (VM, each under a day):**
1. 63: 20-person Revio forced-Tier-3 comparison (~7 h) with peak memory (`/usr/bin/time -v`) → numbers to Marc for the resize.
2. 65b: optionally finish the remaining ~1,090 people (~7 h; `65 --assembly ... --retry-failed`, then `65c_65b_summary.py`
   and round concordance to whole percent, blank where the implied discordant count is 1–19) to get AFR HWE under assembly; add a probe-free check (contigs vs the insertion sequence) since both routes share probes.
3. 68: rerun once to replace the dropped globin row (script now divides by Mapped Reads).
4. Re-derive PUBLIC tables with the cross-table masking (`--aggregate-only` for 58/59/62/65) and check against committed tables.
5. 64: rerun reads (≈1.5 h) for the fixed homopolymer readout; locate the indel for the 345 CDS disagreements; short-read check.
