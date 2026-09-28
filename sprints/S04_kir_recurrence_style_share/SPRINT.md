# Sprint S04 — KIR recurrence/saturation, KIR vs HLA catalogue coverage, figure style, Cole data package, cleanup

*Opened 2026-09-26. Branch `s04-kir-recurrence-style-share` (off S03 @ 5d5ef99). Owner: Marc.
Orchestrator brief: [`ORCHESTRATOR_HANDOFF.md`](ORCHESTRATOR_HANDOFF.md). Log: [`LOG.md`](LOG.md).*

## Board

| WS | Script | Ask | Needs VM | Status |
|---|---|---|---|---|
| WS-A | 44 (VM export) / 45 (local figure) | KIR recurrence + saturation: overall + per ancestry (equal-N as 39), recurrence classes 1 / 2 / >2 / >=20 unrelated carriers, any-level vs protein-level novelty side by side, Good-Turing sample coverage + Chao2 | yes | **FINAL 09-28 (v4b, critic `bd4610f`)**: v1-v3 numbers INVALID (name-based genomic identity + missing KIR artifact filter, see DECISIONS.md); v4b (sequence-hash identity, checkpointed, `22bf543`/`9ca7053`) is the committed final run, 11,856 unrelated. Neither species saturates; novelty mostly private at CDS/protein. |
| WS-B | 46 (local, uses 44 export incl. matched HLA) | KIR vs HLA catalogue coverage: % novel (any/protein), sample coverage, equal-N slope, Chao unseen, per gene | export from 44 | **FINAL 09-28 (v4b, critic `bd4610f`)**: pooled CDS/protein completeness near-tied, but HLA more complete than KIR in every well-powered ancestry (genuine Chao2-pooling/Simpson effect, confirmed not an unequal-N artifact) — reverses the early "KIR better represented" pooled reading. |
| WS-C | _viz_common.py + re-renders | Port cnsplots rules (reference/FIGURE_STYLE.md) into nature_style()/palette/save_fig; re-render 40, 37/37c, 39, 43b (+ fix 43b f/c, Fig1 b/c threshold) ; before/after PNGs | no | **DONE 09-26/27**: layout linter added + Figure 1/DQ/39/43b redesigned to pass it, incl. a real inherited Fig1 panel-d reversal bug found and fixed |
| WS-D | 47 (naive ML, VM) / 48 (ligand co-occurrence, VM) | WS6 plain-language explainer; ranked experiment ideas; naive ML + simple tests with permutation baselines | yes (aggregates only) | **DONE 09-26/27**: 47/48 run on full cohort (11,845), figures + critic review complete |
| WS-E | — | Cole package: novel HLA calls, KIR calls, person table, README/SCHEMA, checksums -> gs://hla-calls-share-wb-cordial-leechee-9743/release_2026-09/ | yes (consent given 09-25; sizes to Marc before upload) | **DONE 09-26**: uploaded + md5-verified to release_2026-09-25/ (6 files, 10.65 MB; Marc ran build/upload after classifier refused agent); cosmetic README date mismatch open for Marc |
| WS-F | — | VM declutter (dry-run list -> Marc OK), repo reorganisation (git mv + refs + tests), single unified branch | yes (consent) | **repo DONE 09-27** (37 scripts + 3 call notes + 2 briefs moved, `outputs/` left as an open question); **VM cleanup PROPOSED, not executed** (`VM_CLEANUP_PLAN_v2.md`, Marc must run it) |
| WS-G | — | LOG, FIGURES_INDEX, VM_RUNS, critics, FINDINGS_FOR_MARC.md | no | **DONE 09-27**: `FINDINGS_FOR_MARC.md` written; context/STATUS,EXPERIMENTS,DECISIONS,ENVIRONMENT updated |

Authorization: Marc's S04 authorization block received in chat 2026-09-25 (items 1-5: app start/stop/resize with config restore + cost/ETA in chat for resized runs >2 h; deploy/run in ~/s04/; Cole package upload after listing files+sizes; VM scratch delete after upload confirm + dry-run OK; repo reorg + one unified branch). App config before: n2-highmem-4, 2000 GB, autostop 1 h.

## Distilled findings

- **WS-A/WS-B: FINAL (v4b), 11,856 unrelated.** S_obs (HLA/KIR): genomic 17,188/38,024 (upper
  bound, inflated by span/UTR + intronic noise — not the headline), CDS 1,363/2,100, protein
  1,079/1,444, protein-novel 198/1,063. Neither species saturates at any level. Pooled CDS/protein
  completeness is near-tied, but **HLA is more complete than KIR in every well-powered ancestry** —
  confirmed a genuine Chao2-pooling (Simpson's-paradox-shaped) effect, not an unequal-N artifact
  (N-weighting the per-ancestry numbers reproduces the per-ancestry gap, not the pooled near-tie).
  Novelty is mostly private (singleton) for both species at CDS/protein level. v1-v3 numbers are
  formally invalid (name-based genomic identity + missing KIR artifact filter, see DECISIONS.md);
  v4b's app-stop mid-run was resolved by adding checkpointing (`22bf543`) and rerunning (`9ca7053`).
- Never write "99% of the allele space explored" — that's the Good-Turing INCIDENCE coverage
  number (92-99% both species), a different and less informative metric than Chao2 RICHNESS
  completeness (8-42%) for the "how well represented is this catalogue" question.
- WS-D: ancestry strongly predicts HLA+KIR carriage (AUROC 0.81-0.98); platform's pooled 0.574
  AUROC collapses within-ancestry (EUR 0.558, AFR 0.489) and is mostly an ancestry echo
  (ancestry-adjusted delta +0.004); HLA->KIR cA/cB has no signal beyond ancestry (adjusted delta
  -0.023); no KIR-HLA ligand pair survives Bonferroni across 30 tests.
- WS-C/WS-E: Cole package uploaded + verified; Figure 1's inherited panel-d reversal (S03 bug) is
  fixed with a regression test.
- Reports: `reports/hla_popgen/44_kir_recurrence_saturation/README.md`,
  `reports/hla_popgen/46_kir_vs_hla_catalogue/README.md` — both v4b-final. Figures listed in
  `FIGURES_INDEX.md` (`fig_catalogue_completeness.png` mid-redesign, same filename). Full
  narrative: `FINDINGS_FOR_MARC.md`.
