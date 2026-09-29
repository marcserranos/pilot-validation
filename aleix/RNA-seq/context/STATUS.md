# Status — live session state (RNA-seq / repertoire workstream)

> **Role:** where we are *right now*, plus the literal next commands. The only file fully
> rewritten each session. Anything durable graduates to ENVIRONMENT, DECISIONS or EXPERIMENTS.

## As of 2026-09-29 — reports 01-07 all complete on real data; deliverables built

All seven reports have run on the full cohort and their outputs are committed under
`aleix/RNA-seq/reports/`. Full results in EXPERIMENTS.md. Headlines:

- **02 replicates immunosenescence at n=7,105**: rarefied diversity -3.3%/decade
  (95% CI -3.5 to -3.1), top-10 clonal share +1.9 pp/decade. Strongest result we have.
- **03 replicates the publicness mechanism** (Spearman 0.66 sharing vs log10 Pgen) and
  exposed that mean-pooled person vectors were unstable, which motivated 06.
- **04 is a rigorous negative**: even against a Pgen-matched generative null, self-peptide
  and low-prevalence-pathogen controls enrich as much as real pathogens. TCR-database
  matching cannot infer exposure at this depth.
- **05**: 60 reproducible motif-organised metaclusters; specificity concentrates, but so do
  self-peptides. Stability 0.41.
- **06**: best person representation is SCEPTR mean + TRBV usage (19.6% split-half
  identifiability vs 0.025% chance); TRBV usage alone gets 16.0%. Read-weighting destroys it.
- **07**: all SCEPTR variants within 0.550-0.579 on epitope structure; the synthetic-data-
  trained variant matches the real one; 16-d reaches 96% of full.

**Deliverables (2026-09-29):** 15-slide .pptx deck (figures embedded) handed to Aleix for
Drive import — the Drive connector cannot take a 700KB file inline. Explainer artifact:
https://claude.ai/artifact/XZJaR66fV9bn5MVbRVnrAa (private; share from the page's Share menu).

**Backups:** embeddings are in `gs://aleix-rnaseq-wb-cordial-leechee-9743/embeddings/` as well
as on the main VM disk.

## Pick up here

1. **Write up the aging result (02)** as the lead finding. Sensitivity analyses are scripted;
   what is owed is the cap-robustness check across top-N (100/250/500/all) noted in
   DECISIONS.md.
2. **Improve the person representation.** 06 gives a benchmark and a baseline to beat
   (TRBV usage, 16.0%). Candidates: supervised pooling; features from clonal structure rather
   than averages; explicitly removing the depth component that tracks the signal.
3. **HLA and disease join**, when Aleix wants it. Constraints now known: adjust for ancestry
   (AUROC 0.72 from repertoire alone) and repertoire size; TRBV usage is the baseline to beat.

Related: [[../../context/STATUS.md]] (root, Marc's).
