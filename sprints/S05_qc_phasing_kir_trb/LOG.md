# S05 log (append-only)

One line per event: `[MM-DD] KIND — what happened (commit / file / number)`. KIND is one of
PLAN, CODE, TEST, RUN, RESULT, FIX, DECISION, BLOCKER, NOTE. Never rewrite past lines.

- [10-11] PLAN — Sprint opened from the 2026-10-06 supervisor asks and the REPORT #11 agenda; plan agreed with Marc in the project thread ("complete everything", with an exhaustive report to follow).
- [10-11] DECISION — Marc: `reports/hla_popgen/novel_alleles.tsv` stays as is; `<20` masking applies to figures and claims (DECISIONS.md).
- [10-11] CODE — `50_regional_depth.py`: one `samtools bedcov` pass per person at MAPQ 0 and 20 over control windows, MHC window, HLA/KIR/TRB genes; relative depth vs control median; TRB copy-number calls with paralog-assessability gate; HWE exact test.
- [10-11] TEST — `tests/test_50_regional_depth.py`: 5/5 pass (planted TRBV4-3 deletion recovered >97% per person, allele frequency within 0.07; paralog gene correctly flagged not assessable; masking; failure accounting; resume).
- [10-11] BLOCKER — This cloud thread cannot drive Marc's browser, so VM steps go through a Claude Code session on Marc's laptop or Marc's own hands (`VM_COMMANDS.md`).
- [10-11] NOTE — The 10-06 S05 work (scripts 50–57c) was uncommitted on Marc's laptop; it was committed by Marc and merged (PR #4). Yak trio QC (WS-D) and the incompatibility map (WS-J, as the 52 atlas) were already done; canonical counts (WS-E) mostly done (50a).
- [10-11] FIX — Renumbered `50_regional_depth.py` → `58_regional_depth.py` (outputs `~/s05/results/58`, `~/s05/internal/58`) to avoid clashing with `50_qc_foundation` in the shared `~/s05` tree. Board re-planned (SPRINT.md).
- [10-11] RESULT — Reading done (`READING.md`): the TRBV4-3/3-2/6-2 deletion is known prior art; the AoU long-read flagship preprint has no HLA allele, DQ or KIR analysis; no prior report of the DQ G1/G2 physical-phasing result found.

