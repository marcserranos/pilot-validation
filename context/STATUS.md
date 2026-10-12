# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-10-13 morning — S05 closed (report delivered); VM jobs pending a restart

**Report:** `sprints/S05_qc_phasing_kir_trb/S05_REPORT.html` (built from `report/src/*.html` by `report/build_report.py`),
published as the S05 Artifact. It has a background course, data trails per analysis, the review dispositions, and a costed
long-job plan (`#plans`). Board/log: `sprints/S05_qc_phasing_kir_trb/{SPRINT,LOG,METHODS_LEDGER}.md`.

**Sprint rule (Marc, 2026-10-12): no job longer than one day.** Longer jobs are planned and costed in the report, not run.

**VM (big_run, n2-highmem-4).** It went "Bad Gateway" (high memory) ~21:50 UTC 10-12 while 65b (assembly mode, stratified
2,000, 4 workers, started 17:21 UTC) ran. Marc OK'd a restart and resuming 65b from its checkpoint with 2 workers. Cap 65b at
24 h total. Remount the bucket after any restart (ENVIRONMENT quirk #14).

**Open, needs Marc** (all in `DECISIONS.md` open questions): ancestry-purity rule; demote protein-novel counts; KIR
complete-only; physical-haplotype rule; no coverage filter; TRB mixture genotyper; Sequel II re-call (resize + cost);
disclosure residue in git history and allele-count masking. Script 60 waits on Aleix copying his TRBV usage table to the bucket.

**Literal next steps (VM, each under a day):**
1. 63 second smoke (2 people, stage-1 identity filter) → record time/person and peak memory → 20-person Revio forced-Tier-3
   comparison → bring numbers to Marc for the resize decision.
2. 65: re-aggregate reads mode with `--genotyper mixture` (aggregate only); when 65b ends, run the concordance by ancestry.
3. Re-derive PUBLIC tables with the cross-table masking (`--aggregate-only` for 58/59/62/65) and check they match the
   locally re-masked committed tables.
4. 64: rerun reads (≈1.5 h) for the fixed homopolymer readout; locate the indel for the 345 CDS disagreements; short-read check.
5. Commit the VM inspection scripts (DQ ≥ 100 kb, contracted KIR).
