# 07 — Which SCEPTR setting is best?

Script: `scripts/07_embedding_settings.py`. Table from the 2026-10-04 polish run (5b83486).
Figure redrawn from that table on 2026-10-20 (`--plot-only`, c79ba26) to fix the row order of
panels b and c. Numbers unchanged.

## Question

SCEPTR's beta-chain model on TRBV + CDR3 was adopted after a check that is circular once V is
an input. Is it the right setting?

## How

- **13 settings:**
  - input: b_sceptr TRBV + CDR3; cdr3_only; paired default fed beta only;
  - capacity: tiny 16-d, small 32, default 64, large 128;
  - training: mlm_only, synthetic_data (trained on OLGA sequences), shuffled_data;
  - architecture: average_pooling, blosum;
  - not learned: CDR3 3-mer SVD, TRBV one-hot.
- **Clonotype level:** VDJdb epitope pairwise AUROC and 5-NN accuracy (4,193 TCRs, 45
  epitopes), with epitope-bootstrap CIs and paired differences from b_sceptr.
- **Person level:** 1,000 people × 250 clonotypes; z-scored split-half identifiability and
  out-of-fold age R², for mean pooling and for mean + TRBV usage.

## Results

- **No SCEPTR variant differs from b_sceptr on epitope AUROC** (every paired CI includes 0).
  b_sceptr 0.579, blosum 0.578, default 0.575, synthetic_data 0.575, large 0.575, tiny 0.567.
- TRBV one-hot (0.552) and 3-mer (0.541) are significantly worse.
  cdr3_only (0.552) is borderline (paired CI −0.078 to +0.0004).
- **Trained only on synthetic OLGA sequences: 0.575**, within 0.004 of the real model. Much of
  the epitope-relevant structure is recombination statistics.
- Tiny (16-d) reaches 0.567 at a quarter of the dimensions.
- Person level is near the floor at 250 clonotypes per person:
  - mean pooling identifiability 0.3–2.0%;
  - mean + TRBV usage 1.0–3.1% (panel f). Pooling matters more than the embedding.
- Best per axis (`summary.csv`): epitope b_sceptr; identifiability mlm_only (2.0%); age R² TRBV
  one-hot (0.13). The two person-level winners are not the epitope winner, and the axes agree
  only moderately (Spearman 0.52, panel e).

## Figure

`fig_embedding_settings`:
- a: epitope AUROC with 95% CI.
- b: identifiability.
- c: age R².
- d: AUROC vs embedding dimension.
- e: epitope AUROC vs identifiability.
- f: mean pooling vs mean + TRBV.

Panels a–c share one row order (a's), and b and c carry no labels of their own.

Redraw off the VM: `python3 scripts/07_embedding_settings.py --plot-only --outdir reports/07_embedding_settings`.

## Files

`summary.csv`, `embedding_settings.csv` (every setting: dimension, time, epitope metrics and
CIs, person metrics).
