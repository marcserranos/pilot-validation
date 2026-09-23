# Sprint S03 — call #8 (2026-09-22): publication figures, DQ G1/G2, KIR, prediction baseline

*Started 2026-09-22 night. Branch `s03-call8-figures-kir-prediction` (off `fig1-drafts-and-research-map` @ bff2de9).
Orchestrator: Claude Opus; implementers: Sonnet subagents. Owner: Marc. Mode: autonomous overnight.*

Inputs: [`../CALL_SUMMARY_2026-09-22.md`](../CALL_SUMMARY_2026-09-22.md) (root), S02 board, recon briefs (in LOG.md).
Orchestrator log (schematic, source for the PI report): [`LOG.md`](LOG.md). VM recipe: [`VM_CHANNEL.md`](VM_CHANNEL.md).

## Resume here
1. Check board below; any WS "running" whose agent died → re-read its WS brief, restart from its own "Resume" line.
2. VM: Verily Workbench, workspace "Full Cohort HLA Calling"; programmatic channel per VM_CHANNEL.md.

## Board

| WS | Script | Ask (call §) | Needs VM | Status |
|---|---|---|---|---|
| WS0 | — | VM bring-up + REST channel | yes | running |
| WS1 | 37_dq_g1g2_signed_ld | Recreate Cole's DQ figure (signed phased D′, G1/G2 2×2); phase-error vs rare-recombination test; phasing-confidence per-person file (VM-only); LD supplement CSV (§2, §3) | yes | todo |
| WS2 | 38_hla_a_deletion_validation | Is HLA-A deletion real? homozygote expectation test; short-read contradiction test; simplified deletion colours (§4) | yes | todo |
| WS3 | 39_saturation_by_ancestry | Per-ancestry equal-N discovery/saturation curves (Pakistan Fig 3e style); % allele space explored vs IPD-IMGT (§6) | yes (uncensored counts) | todo |
| WS4 | 40_figure1_v5 | Figure 1 per panel feedback (§6) — Nature-grade | yes (panels a/b) | todo |
| WS5 | 41_kir_pilot | KIR chr19 scoping: coordinates, Immuannot KIR support, pilot runtime/cost, full run if pilot clean (§7) | yes | done -- GO, conditional (orchestrator flags + framework-gene spot check); see LOG.md and reports/hla_popgen/41_kir_scoping/README.md |
| WS6 | 42_repertoire_baseline | BenchRep-T-style VJ-k-mer baseline (L1-LR, XGBoost; 3-fold stratified CV; AUROC/AUPRC) on autoimmune/chronic phenotypes, repertoire-only (§9–10) | yes | todo |
| WS7 | — | Critic passes (fresh-context) on each figure/claim | no | todo |

## Distilled findings
(filled as workstreams land)
