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
| WS0 | — | VM bring-up + REST channel | yes | done -- app up; VM automation via auto-mode Chrome tools ultimately blocked by the permission classifier ("Auto-Mode Bypass"), so remaining VM steps are handed to Marc (VM_HANDOFF.md) |
| WS1 | 37_dq_g1g2_signed_ld | Recreate Cole's DQ figure (signed phased D′, G1/G2 2×2); phase-error vs rare-recombination test; phasing-confidence per-person file (VM-only); LD supplement CSV (§2, §3) | yes | done -- complete purge (0/17,255 at 4-field, unrelated set 11,856); held-out 5-fold EM rerun done (37e, 59fef06): 0.24% spurious vs ~24% for a naive linkage-equilibrium phaser |
| WS2 | 38_hla_a_deletion_validation | Is HLA-A deletion real? homozygote expectation test; short-read contradiction test; simplified deletion colours (§4) | yes | done -- see reports/hla_popgen/38_hla_a_deletion_validation/README.md |
| WS3 | 39_saturation_by_ancestry | Per-ancestry equal-N discovery/saturation curves (Pakistan Fig 3e style); % allele space explored vs IPD-IMGT (§6) | yes (uncensored counts) | done -- full cohort (11,856 unrelated), corrected extrapolation methodology; see reports/hla_popgen/39_saturation_by_ancestry/README.md |
| WS4 | 40_figure1_v5 | Figure 1 per panel feedback (§6) — Nature-grade | yes (panels a/b) | done -- all 6 panels; 40a/40b exports run 09-24, composed (0ec28ff) |
| WS5 | 41_kir_pilot, 43_kir_full_cohort | KIR chr19 scoping, pilots, full-cohort run (§7) | yes | done -- 170-person pilot aggregated; FULL COHORT run 09-24/25: 12,261 people, 31.1 h on n2-highcpu-80, ~$95; 43 aggregate + figure (8f58dee), critic #2 (abdea67); app restored to n2-highmem-4 |
| WS6 | 42_repertoire_baseline | BenchRep-T-style VJ-k-mer baseline (§9–10) | yes | partial -- 1,500-person subsample (b9ff562): AFR AUROC 0.82, EUR 0.67, sex 0.55, no disease ≥100 cases; full 7,640 + XGBoost pending |
| WS7 | — | Critic passes (fresh-context) on each figure/claim | no | done -- CRITIC_1 (37/38/39/41), CRITIC_2 (43: 0 blockers, 3 major, 4 minor; fixes in abdea67, 104-person reconciliation by orchestrator) |

## Distilled findings
(numbers exact from each result's own README; provisional items marked)

- **WS1 (37, DQ G1/G2 signed LD):** Physically phased cis haplotypes show a complete purge: **0 of 17,255** DQA1~DQB1 pairs at 4-field (0/22,341 at 2-field) are G1/G2-incompatible, in every one of 6 ancestries, vs 300-8,320 expected under independence. A population-EM statistical rephaser (`37d`) also manufactured 0 spurious incompatible pairs -- but that test is **circular** (critic #1 major #3: EM was fit on a pool that already includes the truth set) and the held-out-fold fix (`37e`) is written but not yet run on the VM (provisional pending VM_HANDOFF.md (b)). The one regime that does produce incompatible-looking pairs: naive same-hap-label pairing when DQA1/DQB1 sit on different contigs (<20/320 pooled pairs) -- read as assembly/contig-boundary artifact, not recombination or generic EM error.
- **WS2 (38, HLA-A deletion validation):** The 1.83% raw HLA-A "deletion" rate is **82.8% contradicted** by AoU short-read genotypes (Fisher P=4e-8 vs non-carriers), matching every other implausible control gene (37-83%) and unlike true DRB3/4/5 deletions (<1% contradicted). Recommended paper framing: true HLA-A deletion rate < ~0.3-0.4% (upper bound); 1.83% is not usable as a biological estimate.
- **WS3 (39, saturation by ancestry):** No ancestry is saturated (all discovery curves still rising). Equal-N discovery slope at N*=1,236 (primary, cross-ancestry-comparable metric): **AFR highest** (101.6 new alleles/1,000), supporting the call's original expectation. At equal-N richness, AMR leads in raw distinct-allele count (plausibly admixture, since AMR still leads AFR under strict >=0.95 ancestry assignment). MID's Chao2 estimate is flagged low-power (smallest N=487).
- **WS4 (Figure 1 v5):** complete (0ec28ff). Panel a: 589 admixture bins of ≥20 unrelated people. Panel c: any-field novelty by gene × ancestry (strict ≥0.9); only HLA-A×MID censored (hatched). Top cell DRB1×EAS 70.9%. Open: panel b uses 0.98 vs c at 0.9 ancestry threshold.
- **WS5 (41/43, KIR):** Full cohort called (12,261 people, 31.1 h, ~$95; HLA tables md5-verified untouched). Unrelated 11,882: 99.5% of haplotypes with KIR calls, 9.23 genes/hap, 58.9% novel (mostly non-coding; protein novelty 4–9%/gene, KIR2DL5B 22.5%), framework genes 93.7–96.5%, cA share EAS 66.7% → SAS 41.4%. Replicates the 20/170-person pilots. Open: 2DL2/2DL3 co-occurrence 0.8% (expected ~0); 991 sequel2 people excluded (self-align chr6-cache bug); novelty recurrence check not done.
- **WS6 (42, repertoire baseline):** 1,500-person subsample (b9ff562; 3×3-fold CV, L1-LR): repertoire→ancestry AFR 0.82, EUR 0.67; sex 0.55; no disease ≥100 cases. Full 7,640 + XGBoost still pending.
- **WS7 (critics):** Round 1 (CRITIC_1.md) found 2 blockers -- both addressed (small counts redacted; unrelated-N mislabel understood as a `37_vm_run.py`-only bug, worked around in `37e`) -- and 1 still-open major (EM circularity, see WS1). No round 2 critic pass run yet on 40/41-extended/42.
