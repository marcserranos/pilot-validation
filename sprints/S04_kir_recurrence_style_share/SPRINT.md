# Sprint S04 — KIR recurrence/saturation, KIR vs HLA catalogue coverage, figure style, Cole data package, cleanup

*Opened 2026-09-26. Branch `s04-kir-recurrence-style-share` (off S03 @ 5d5ef99). Owner: Marc.
Orchestrator brief: [`ORCHESTRATOR_HANDOFF.md`](ORCHESTRATOR_HANDOFF.md). Log: [`LOG.md`](LOG.md).*

## Board

| WS | Script | Ask | Needs VM | Status |
|---|---|---|---|---|
| WS-A | 44 (VM export) / 45 (local figure) | KIR recurrence + saturation: overall + per ancestry (equal-N as 39), recurrence classes 1 / 2 / >2 / >=20 unrelated carriers, any-level vs protein-level novelty side by side, Good-Turing sample coverage + Chao2 | yes | **done 09-26**: 3 figures (`fig_saturation_by_recurrence`, `fig_saturation_per_ancestry`, `fig_coverage_completeness`) + `recurrence_stats.tsv` (singleton-share z-test) on the real 44 exports; 47 tests pass, `check_layout(strict=True)`-clean |
| WS-B | 46 (local, uses 44 export incl. matched HLA) | KIR vs HLA catalogue coverage: % novel (any/protein), sample coverage, equal-N slope, Chao unseen, per gene | export from 44 | **done 09-26**: `46_catalogue_metrics.tsv` + 2 figures (`fig_catalogue_completeness`, `fig_recurrence_composition`); KIR catalogue better represented (86.8% vs 66.8% mean per-gene Chao2 completeness); 9 tests incl. a label<->TSV binding test |
| WS-C | _viz_common.py + re-renders | Port cnsplots rules (reference/FIGURE_STYLE.md) into nature_style()/palette/save_fig; re-render 40, 37/37c, 39, 43b (+ fix 43b f/c, Fig1 b/c threshold) ; before/after PNGs | no | wave 1: running |
| WS-D | 47 (naive ML, VM) / 48 (ligand co-occurrence, VM) | WS6 plain-language explainer; ranked experiment ideas; naive ML + simple tests with permutation baselines | yes (aggregates only) | wave 1: ideation + ranking |
| WS-E | — | Cole package: novel HLA calls, KIR calls, person table, README/SCHEMA, checksums -> gs://hla-calls-share-wb-cordial-leechee-9743/release_2026-09/ | yes (consent given 09-25; sizes to Marc before upload) | DONE 09-26: uploaded + md5-verified to release_2026-09-25/ (6 files, 10.65 MB; Marc ran build/upload after classifier refused agent) |
| WS-F | — | VM declutter (dry-run list -> Marc OK), repo reorganisation (git mv + refs + tests), single unified branch | yes (consent) | last |
| WS-G | — | LOG, FIGURES_INDEX, VM_RUNS, critics, FINDINGS_FOR_MARC.md | no | ongoing |

Authorization: Marc's S04 authorization block received in chat 2026-09-25 (items 1-5: app start/stop/resize with config restore + cost/ETA in chat for resized runs >2 h; deploy/run in ~/s04/; Cole package upload after listing files+sizes; VM scratch delete after upload confirm + dry-run OK; repo reorg + one unified branch). App config before: n2-highmem-4, 2000 GB, autostop 1 h.

## Distilled findings

- **WS-A/WS-B answer**: KIR's reference catalogue (IPD-KIR 2.13.0) is better represented in this
  cohort than HLA's (IPD-IMGT/HLA 3.55.0), both pooled (Chao2 completeness 84.2% vs 63.7%) and
  per-gene (mean 86.8% vs 66.8%) — despite KIR's higher raw novelty rate (54.8% vs 29.7% mean
  per-gene any-level-novel) — because KIR's excess novelty is disproportionately SHARED across
  unrelated people (80.2% non-singleton) rather than private (59.1% non-singleton for HLA),
  confirmed by a two-proportion z-test (p≈6e-20). This flips at the protein level: KIR's own
  genuine amino-acid-changing novelty is 88.2% private/singleton, similar in character to ordinary
  rare variation, not a systematic catalogue gap.
- Never write "99% of the allele space explored" — that's the Good-Turing INCIDENCE coverage
  number (≥99% both species), a different and less informative metric than Chao2 RICHNESS
  completeness (64-84%) for the "how well represented is this catalogue" question.
- Two data-quality flags for Marc/Aleix (not fixed by WS-A/B, both documented in the 44/46
  READMEs): (1) HLA's `protein_novel` S_obs is 0 for every one of the 8 classical genes,
  unreconciled against S03's non-zero novel-protein counts; (2) KIR's `pct_novel_protein` exceeds
  100% for every gene (identity-granularity mismatch between the protein-hash count and the
  genomic-level denominator) — never cite it.
- Reports: `reports/hla_popgen/44_kir_recurrence_saturation/README.md` (WS-A, updated),
  `reports/hla_popgen/46_kir_vs_hla_catalogue/README.md` (WS-B, new). Figures listed in
  `FIGURES_INDEX.md`.
