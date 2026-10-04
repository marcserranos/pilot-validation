# Status — live session state (RNA-seq / repertoire workstream)

> **Role:** where we are *right now*, plus the literal next commands. The only file fully
> rewritten each session. Anything durable graduates to ENVIRONMENT, DECISIONS or EXPERIMENTS.

## As of 2026-10-04 — methods week: review done, polish run ready for the VM

This week is methods consolidation, not new analysis. The methods artifact documents every
stage (https://claude.ai/artifact/DR6G9cPp4DKPeBSVyeZ1vi). Building it found nine code issues;
all are fixed in code and tested on synthetic data (DECISIONS.md, 2026-10-04). The VM run
that applies them is `scripts/run_polish.sh`.

Results through 2026-09-29 (EXPERIMENTS.md) stand except where the polish run changes them:
01, 03, 05, 06 depend on the top-500 pool (tie-break fix); 06 and 07 person-level numbers
change by design; 04/07 gain CIs. Report 02 (ageing, the lead result) is unaffected.

## Pick up here

1. On the VM:
   `cd ~/repos/pilot-validation && git pull && nohup bash aleix/RNA-seq/scripts/run_polish.sh > ~/polish.log 2>&1 &`
   Progress: `grep '^====' ~/polish.log`. Resumable: rerun the same command after a failure.
   Expected ~6-8 h (full-cache embedding ~2 h, pool re-embed ~50 min, 04 and 07 ~1 h each).
2. When it finishes (it commits and pushes): pull, read 08 first (how big was the tie skew?),
   then compare 01/03/05/06/07 with the v1 numbers, append results to EXPERIMENTS.md, update
   the methods artifact (software versions from `reports/08_methods_audit/environment.txt`,
   open issues 1-9 to resolved), and decide the identifiability headline (DECISIONS open item).
3. Then the ageing write-up (02), and the HLA/disease join when Aleix wants it.

Related: [[../../context/STATUS.md]] (root, Marc's).
