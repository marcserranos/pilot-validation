# Status — live session state (RNA-seq / repertoire workstream)

> **Role:** where we are *right now*, plus the literal next commands. The only file fully
> rewritten each session. Anything durable graduates to ENVIRONMENT, DECISIONS or EXPERIMENTS.

## As of 2026-10-16 — W0–W5 of the 10-08 plan done; W6 handover and write-ups left

Last work: 2026-10-13 (commit b47314e). The 10-08 plan (we own everything RNA-seq/TRB; map
disease and ancestry onto it) ran W0 → W5; results are in EXPERIMENTS.md (10-08, 10-12, 10-13).

Headline numbers to use (detail in EXPERIMENTS.md):
- **W0 (09)**: analysis_main 6,764 people (sex is a covariate, not a filter); pure90 4,271.
- **W1 (10)**: V+J predicts ancestry, AUROC 0.804 (covariates only 0.614). **TRBV4-3 never used
  by 27.3%** (SAS 43.5% … EUR 16.7%), a germline signature; TRBV12-4 15.6% (may be 12-3
  call assignment). PC1 is almost only TRBV4-3, not ancestry.
- **W4 (13)**: HIV is the positive control; **heart failure** shows clonal expansion
  (top-10 share +0.32 SD, Bonferroni, survives utilisation adjustment); hypertension and T2D
  weaker; no autoimmune disease significant. No CMV serostatus.
- **W5 (14)**: repertoire adds at most ~1 AUROC point over age/sex/PCs/depth; only diversity
  ever helps; SCEPTR, 3-mers and V/J add nothing.
- Earlier, still current: 06 identifiability 27.0% z-scored; 03 stability 0.76; cap 500.

Links: methods artifact v3 https://claude.ai/artifact/DR6G9cPp4DKPeBSVyeZ1vi ·
one-pager https://claude.ai/artifact/CyJ3ktfg3xuJeJWogCZRMv (both predate W0–W5).

## Pick up here

1. **W6 handover notes** — the last step of the 10-08 plan; not started.
2. **Marc is waiting on the TRBV usage table**: his script 60 (`60_trb_germline_usage.py`,
   main) joins per-person TRBV usage to his TRBV copy-number / TRBV4-3 insertion genotypes
   (script 65, 12,233 people). Copy the table to the bucket Marc reads (confirm path with him).
   Natural check: our 27.3% TRBV4-3 non-users vs his 65 0-copy rates by ancestry.
3. Rewrite `reports/01_sceptr_embedding_viz/README.md` (still describes the old size region).
   Reports 02–14 have no README; their write-up lives in EXPERIMENTS.md.
4. Methods document, the ageing write-up (02), and the HLA-label join (J1–J5 in main's
   `reports/hla_popgen/NEXT_STEPS_AND_RESEARCH_MAP.md` §3).
5. Open questions in DECISIONS.md still unworked: per-person per-chain ancestry recovery
   test (B-cell chains carry the gap), rarefaction metric, research_ids in the public repo.

Related: [[../../context/STATUS.md]] (root, Marc's; on `main` it is
`context/STATUS.md`, last 2026-10-13).
