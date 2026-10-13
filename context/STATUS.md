# Status — live session state

> **Role:** where we are *right now* + the literal next commands. The only file that gets fully rewritten each session.
> **Edit:** rewrite compactly at each session end. Nothing here is durable — a fact that outlives this session graduates to ENVIRONMENT (a quirk/runbook change), DECISIONS (a call), or EXPERIMENTS (a result).
> **Read:** to pick up work.

## As of 2026-10-13 ~09:00 UTC — S05 closed (report delivered); 65b running until 11:00 UTC

**Report:** `sprints/S05_qc_phasing_kir_trb/S05_REPORT.html` (built from `report/src/*.html` by `report/build_report.py`),
published as the S05 Artifact. It has a background course, data trails per analysis, the review dispositions, and a costed
long-job plan (`#plans`). Board/log: `sprints/S05_qc_phasing_kir_trb/{SPRINT,LOG,METHODS_LEDGER}.md`.

**Sprint rule (Marc, 2026-10-12): no job longer than one day.** Longer jobs are planned and costed in the report, not run.

**VM (big_run, n2-highmem-4).** Restarted ~08:17 UTC 10-13 (after the 10-12 "Bad Gateway" and an overnight storage stall:
27 timeouts 01:29–07:48 UTC). 65b resumed 08:23 UTC with 2 workers and `--retry-failed`: 673/2,000 done then, ~1,600 projected.
**Marc (10-13, revised 08:34 UTC): stop 65b gracefully at 11:00 UTC (13:00 CEST) 10-13 even if unfinished, and aggregate whatever is done** (PUBLIC mode, report
completed count by ancestry + concordance). Remount the bucket after any restart (ENVIRONMENT quirk #14).

**Open, needs Marc** (all in `DECISIONS.md` open questions): ancestry-purity rule; demote protein-novel counts; KIR
complete-only; physical-haplotype rule; no coverage filter; TRB mixture genotyper; Sequel II re-call (resize + cost);
disclosure residue in git history and allele-count masking. Script 60 waits on Aleix copying his TRBV usage table to the bucket.

**Literal next steps (VM, each under a day):**
1. 63: second smoke PASSED 10-13 (~21 min/person; 991 people ≈ 350 h, see report `#plans`). Next: 20-person Revio
   forced-Tier-3 comparison (~7 h) with peak memory (`/usr/bin/time -v`) → numbers to Marc for the resize decision.
2. 65: mixture re-aggregation DONE on the VM (`~/s05/results/65_reads{,_unrelated}_mixture`, run with `--probes ~/s05/results/65_v3/probes.fa`;
   not pulled). At 11:00 UTC stop 65b, aggregate, run the concordance by ancestry.
3. Re-derive PUBLIC tables with the cross-table masking (`--aggregate-only` for 58/59/62/65) and check they match the
   locally re-masked committed tables.
4. 64: rerun reads (≈1.5 h) for the fixed homopolymer readout; locate the indel for the 345 CDS disagreements; short-read check.
5. Commit the VM inspection scripts (DQ ≥ 100 kb, contracted KIR).
