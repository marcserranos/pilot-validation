# WS6 — BenchRep-T VJ-kmer protocol (extracted for replication)

Source: Im, Cohen-Lavi, Buendia, Kundaje, Boyd, bioRxiv 10.64898/2026.06.09.727013v1
("the onshore paper"). PMC mirror (PMC13277846) carries only the abstract — the author's
license does not permit PMC archiving of the full text, so all section/figure refs below are
from the bioRxiv full-text HTML, fetched 2026-09-22. Not independently re-verified against a
PDF; treat section numbers as [MED] confidence, everything else quoted as [HIGH].

## 1. Feature construction — "VJ-kmer"

Two feature blocks, concatenated per sample, each independently normalized to relative
frequency:

- **V/J gene usage.** Usage counted and normalized to relative frequency **separately within
  V and within J** (not jointly) — i.e. two simplexes, not one V×J joint table. (This is the
  "VJ" half of "VJ-kmer".)
- **CDR3 k-mers.** **4-mers, plus the four single-position gapped variants of each 4-mer**
  (one wildcard position each — "X_YZ" pattern), counted over the CDR3 and normalized by
  total k-mer count per sample (Appendix §B.2.2). Ungapped 4-mers + gapped variants together
  form the k-mer feature block.
- **CDR3 trimming.** First and last residue of each CDR3 trimmed off before k-mer counting —
  removes the flanking conserved Cys/Phe, matching Mal-ID's convention (Appendix §B.1). Do
  this before k-mer extraction, not after.
- **Productive-only filter.** Applied where the annotation supports it (Mal-ID cohorts drop
  non-productive rearrangements per §B.1); immunoSEQ-style cohorts without a productivity
  flag skip the filter. TRUST4 does carry `CDR3_amino_acids` with `_`/`?` markers for
  stop/ambiguous — treat those as non-productive and drop, matching the spirit of the rule.
- **Weighting.** Base VJ-kmer model counts **unique sequences, not reads or clone sizes** —
  no abundance weighting in the feature construction itself (read_count/frequency are not
  folded into the k-mer or V/J tallies).
- **Per-sample normalization.** Both blocks are per-repertoire relative-frequency vectors —
  makes samples of different depth comparable in vector space, though see §6 below: this does
  NOT substitute for a depth-confound check.

## 2. Models and hyperparameter grids

- **L1-regularized logistic regression ("VJ-kmer Reg").** Inverse-regularization grid
  `C ∈ {1.0, 0.2, 0.1, 0.05, 0.03}`, selected by **5-fold stratified CV using AUROC**
  (nested inside each outer training fold — §B.2.2). An ensemble weight
  `α ∈ {0.0, 0.1, …, 1.0}` (V/J-usage model vs. k-mer model) is tuned on an 80/20 internal
  validation split.
- **XGBoost ("VJ-kmer XG").** Two-stage grid search:
  - Stage 1: `max_depth ∈ {3,4,5,6} × learning_rate ∈ {0.01,0.03,0.05,0.1}` (16 combos).
  - Stage 2 (conditioned on stage-1 winner): `subsample ∈ {0.7,0.8,1.0} ×
    colsample_bytree ∈ {0.7,0.8,1.0} × min_child_weight ∈ {1,3,5}` (27 combos).
  - Early stopping: 20 rounds without validation-metric improvement, cap 1000 rounds.

## 3. Cross-validation

- **3 disjoint folds per dataset**, stratified by disease label where the dataset isn't
  pre-split; **the same 3 folds are reused across every method being benchmarked** (§2.2) —
  this is the part to not skip: the fold assignment must be identical for repertoire-only,
  covariate-only, and repertoire+covariate models, and for every model class, so all AUROC/
  AUPRC numbers are comparable.
- Hyperparameter tuning happens strictly inside each fold's training partition (internal
  80/20 split or internal 5-fold CV, never touching the held-out fold) — i.e. nested CV for
  model selection, not for the reported metric itself (the reported metric is computed only
  from the pooled held-out predictions).

## 4. Metrics

- **AUROC and AUPRC computed on pooled out-of-fold predictions** across the 3 folds (§2.2) —
  not averaged per-fold-then-averaged; pool the predictions first, then score once.
- No explicit bootstrap CI in the main pipeline; error bars only appear in the
  sequencing-depth robustness experiment (Fig. A.3). WS6 adds bootstrap CIs on top (resample
  people with replacement from the pooled OOF predictions, e.g. 1000 reps) since our own
  n per phenotype is much smaller than BenchRep-T's cohorts and a point estimate alone would
  be misleading.

## 5. Class imbalance

- VJ-kmer models (Reg and XG) do **not** document explicit reweighting — imbalance is left to
  AUROC/AUPRC as metrics that are more informative than accuracy under imbalance, plus
  stratified folds keeping the case fraction stable across folds. Deep-learning baselines
  (DeepRC, DeepTCR) DO use `n_neg/n_pos` positive-class weighting (§B.2.4) — not applicable to
  the VJ-kmer track WS6 replicates, so WS6 does not add class weighting to LR/XGBoost either,
  to stay faithful.

## 6. Depth / downsampling control

- Downsampling indices are **pre-generated once per repertoire with a fixed seed**, and
  **nested across depths** — the D sequences at depth D are always a subset of the sequences
  used at any larger depth for the same repertoire (§B.4). Depths tested: {1k, 5k, 10k, 25k,
  50k, 75k}, 5 independent replicates each.
- WS6 applies the same nesting principle for its rarefied-repertoire variant: downsample every
  person's repertoire to one common target read count (fixed seed, nested subsets) before
  refitting features, to show that repertoire signal isn't just tracking depth.

## What WS6 replicates vs. adapts

Replicated faithfully: k-mer construction (4-mer + single-gap, CDR3-trimmed, per-sample
normalized), V/J usage normalization, the 3-fold-stratified-shared-folds CV discipline, pooled
OOF AUROC/AUPRC, the LR C-grid and XGBoost two-stage grid.

Adapted for AoU: TRUST4 report.tsv gives read_count/frequency per clone rather than
per-read sequences — WS6 treats each report.tsv row as one clone (unique CDR3+V+J+C) and
applies the k-mer/V-J tallies at clone-level (matching BenchRep-T's "unique sequences, not
reads" convention). Two chains modeled separately (TRB, IGH) and combined. Added: a
covariate-only baseline (age, sex, ancestry, depth) and repertoire+covariate combination,
which BenchRep-T doesn't need (its cohorts are pre-matched) but our EHR-derived phenotypes
require, plus bootstrap CIs given smaller n.
