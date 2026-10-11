# sprints/_calls/

Point-in-time call artifacts (briefs, summaries, handoffs from a specific supervisor/collaborator
call), not tied to any one sprint folder. Filed here instead of loose at `sprints/` top level, with
a leading `YYYY-MM-DD_` prefix for sort order and so they read as dated records rather than living
docs. The underscore prefix on the directory keeps it sorted above the `S0N_*` sprint dirs,
signalling "meta," not a sprint board.

Moved here from a loose `sprints/` top level in the 2026-09-27 WS-F reorg (see
`context/EXPERIMENTS.md`'s matching entry for the old→new path table). These are historical
narrative — read them, don't edit them; if a call needs a follow-up note, add it to the relevant
sprint's own `LOG.md` instead.

Current contents:
- `2026-09-22_CALL_BRIEF.md`, `2026-09-22_CALL_SUMMARY.md` — the 2026-09-22 supervisor call (Marc,
  Cole): brief going in, summary coming out. Referenced from several `scripts/hla_popgen/NN_*.py`
  docstrings and `reports/hla_popgen/NN_*/README.md` files as the source of a design decision.
- `2026-09-22_ALEIX_HANDOFF.md` — the same-day handoff to/from Aleix's workstream.

Calls with no written notes in the repo (as of 2026-10-11). Decks are in Marc's Google Drive,
folder of the `REPORT #N - HLA_AoU - Marc Serrano and Aleix Ruf` slides:
- **Call #9, ~2026-09-29** (deck `REPORT #9`): RNA-seq repertoires and embeddings (Aleix), DQ G1/G2
  signed LD and phasing, HLA-A deletion validation, saturation by ancestry, Figure 1, KIR full cohort.
- **Call #10, ~2026-10-06** (deck `REPORT #10`): agenda "#1 What are the counts exactly?",
  "#2 Phasing QC". Aleix's branch records the outcome only as "we own everything RNA-seq/TRB; map
  disease and ancestry onto it" (his W0–W6 plan).
- **Call #11, upcoming** (deck `REPORT #11`, created 2026-10-07): agenda "#1 What are the counts
  exactly?", "#2 Phasing QC — Yak trioeval still running", "#3 Other branches — expression of
  immune genes (gene → expression), incompatibility map".
If a summary for #9 or #10 is written later, file it here as `YYYY-MM-DD_CALL_SUMMARY.md`.

- `2026-10-06_CALL_SUMMARY.md` — call #10 (Cole, David, Marc), from a partial transcript Marc
  pasted on 2026-10-11. PheWAS of repertoire
  features, clonotype × HLA and V/J-usage GWAS plans, regional coverage QC, KIR.
