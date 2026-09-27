# Sprint S04 — KIR recurrence/saturation, KIR vs HLA catalogue coverage, figure style, Cole data package, cleanup

*Opened 2026-09-26. Branch `s04-kir-recurrence-style-share` (off S03 @ 5d5ef99). Owner: Marc.
Orchestrator brief: [`ORCHESTRATOR_HANDOFF.md`](ORCHESTRATOR_HANDOFF.md). Log: [`LOG.md`](LOG.md).*

## Board

| WS | Script | Ask | Needs VM | Status |
|---|---|---|---|---|
| WS-A | 44 (VM export) / 45 (local figure) | KIR recurrence + saturation: overall + per ancestry (equal-N as 39), recurrence classes 1 / 2 / >2 / >=20 unrelated carriers, any-level vs protein-level novelty side by side, Good-Turing sample coverage + Chao2 | yes | **PENDING 44 v4 (09-27)**: v1-v3 numbers INVALID (name-based genomic identity + missing KIR artifact filter, see DECISIONS.md); fix committed (`e4211ad`), 44 v4 smoke-tested (300 people, clean) but full-cohort run not yet started. Do not cite v1-v3 numbers. |
| WS-B | 46 (local, uses 44 export incl. matched HLA) | KIR vs HLA catalogue coverage: % novel (any/protein), sample coverage, equal-N slope, Chao unseen, per gene | export from 44 | **PENDING 44 v4 (09-27)**: same invalidation as WS-A ("KIR catalogue better represented" does not survive the identity fix); rebuild after 44 v4 completes. |
| WS-C | _viz_common.py + re-renders | Port cnsplots rules (reference/FIGURE_STYLE.md) into nature_style()/palette/save_fig; re-render 40, 37/37c, 39, 43b (+ fix 43b f/c, Fig1 b/c threshold) ; before/after PNGs | no | **DONE 09-26/27**: layout linter added + Figure 1/DQ/39/43b redesigned to pass it, incl. a real inherited Fig1 panel-d reversal bug found and fixed |
| WS-D | 47 (naive ML, VM) / 48 (ligand co-occurrence, VM) | WS6 plain-language explainer; ranked experiment ideas; naive ML + simple tests with permutation baselines | yes (aggregates only) | **DONE 09-26/27**: 47/48 run on full cohort (11,845), figures + critic review complete |
| WS-E | — | Cole package: novel HLA calls, KIR calls, person table, README/SCHEMA, checksums -> gs://hla-calls-share-wb-cordial-leechee-9743/release_2026-09/ | yes (consent given 09-25; sizes to Marc before upload) | **DONE 09-26**: uploaded + md5-verified to release_2026-09-25/ (6 files, 10.65 MB; Marc ran build/upload after classifier refused agent); cosmetic README date mismatch open for Marc |
| WS-F | — | VM declutter (dry-run list -> Marc OK), repo reorganisation (git mv + refs + tests), single unified branch | yes (consent) | **repo DONE 09-27** (37 scripts + 3 call notes + 2 briefs moved, `outputs/` left as an open question); **VM cleanup PROPOSED, not executed** (`VM_CLEANUP_PLAN_v2.md`, Marc must run it) |
| WS-G | — | LOG, FIGURES_INDEX, VM_RUNS, critics, FINDINGS_FOR_MARC.md | no | **DONE 09-27**: `FINDINGS_FOR_MARC.md` written; context/STATUS,EXPERIMENTS,DECISIONS,ENVIRONMENT updated |

Authorization: Marc's S04 authorization block received in chat 2026-09-25 (items 1-5: app start/stop/resize with config restore + cost/ETA in chat for resized runs >2 h; deploy/run in ~/s04/; Cole package upload after listing files+sizes; VM scratch delete after upload confirm + dry-run OK; repo reorg + one unified branch). App config before: n2-highmem-4, 2000 GB, autostop 1 h.

## Distilled findings

- **WS-A/WS-B: PENDING 44 v4 — no v1-v3 number is a finding.** Three rounds of "KIR catalogue is
  better represented than HLA's" and similar claims were each built on a broken identity
  definition (genomic identity was name-based, not sequence-based; KIR protein calls lacked HLA's
  artifact filter) — caught via an impossible result (more distinct novel KIR proteins than
  distinct KIR genomic alleles). Fix committed (`e4211ad`), smoke-tested clean, full-cohort run
  ("44 v4") not yet done. Rerun 44 v4 -> 45/46 before quoting any KIR-vs-HLA catalogue number.
- Never write "99% of the allele space explored" — that's the Good-Turing INCIDENCE coverage
  number (>=99% both species, pending v4 recompute), a different and less informative metric than
  Chao2 RICHNESS completeness for the "how well represented is this catalogue" question.
- WS-D/WS-E/WS-C are solid: ancestry strongly predicts HLA+KIR carriage (AUROC 0.81-0.98); platform
  AUROC 0.573 is an ancestry echo (adjusted delta +0.004); no KIR-HLA ligand pair survives
  Bonferroni across 30 tests; Cole package uploaded + verified; Figure 1's inherited panel-d
  reversal (S03 bug) is fixed with a regression test.
- Reports: `reports/hla_popgen/44_kir_recurrence_saturation/README.md`,
  `reports/hla_popgen/46_kir_vs_hla_catalogue/README.md` — both marked pending 44 v4 update.
  Figures listed in `FIGURES_INDEX.md`. Full narrative: `FINDINGS_FOR_MARC.md`.
