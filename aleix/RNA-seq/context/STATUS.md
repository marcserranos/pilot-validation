# Status — live session state (RNA-seq / repertoire workstream)

> **Role:** where we are *right now*, plus the literal next commands. The only file fully
> rewritten each session. Anything durable graduates to ENVIRONMENT, DECISIONS or EXPERIMENTS.

## As of 2026-10-06 — methods week: polish run done, methods artifact updated

The polish run finished 2026-10-04 (commit 5b83486); results in EXPERIMENTS.md. Methods
artifact (v3): https://claude.ai/artifact/DR6G9cPp4DKPeBSVyeZ1vi

Numbers that changed and must be used from now on:
- 06 headline: SCEPTR mean + TRBV usage **27.0% z-scored** (29.4% raw), depth R2 0.17.
- 03 person-vector stability: 0.76 (not 0.35).
- 01: no repertoire-size effect on the person map (the old one was a tie-break artifact).
- 05: report crude and within-V enrichment side by side.
- 07: SCEPTR variants statistically indistinguishable; baselines worse.
Report 02 (ageing) unchanged.
One-pager: https://claude.ai/artifact/CyJ3ktfg3xuJeJWogCZRMv

## Pick up here

1. Decide the identifiability headline (DECISIONS open item; recommendation: z-scored).
2. Rewrite `reports/01_sceptr_embedding_viz/README.md` (still describes the old size region).
3. BAM-cap comparison: DONE 2026-10-06 (excluded = 2.5x deeper; see EXPERIMENTS). Optional: TRUST4 on the 404.
4. Consider a cap near 300 for the pool (~200 of 500 slots are random among 1-2-read clonotypes).
5. Then the official methods document, the ageing write-up, and the HLA/disease join.

Related: [[../../context/STATUS.md]] (root, Marc's).
