# 05 — Metaclusters of the SCEPTR TRB space

Script: `scripts/05_embedding_clusters.py`. Current outputs from the 2026-10-04 polish run
(5b83486). Counts re-read from `vdjdb_cluster_enrichment.csv` on 2026-10-20 (see Correction).

## Question

If the 3.84 M embedded clonotypes are split into k = 60 metaclusters, what separates them?
And do VDJdb TCRs of known specificity concentrate in particular clusters?

## How

- **Clustering:** MiniBatchKMeans, k = 60. This is a partition used to describe the space, not
  a claim of discrete cell types. Stability is the mean best Jaccard against a second seed.
- **Enrichment:** for each cluster × pathogen pair, Fisher's exact test of exact VDJdb matches
  inside vs outside the cluster, with BH correction (364 tests).
- **Within-TRBV version:** Cochran–Mantel–Haenszel test stratified by TRBV gene, since
  clusters are largely V-gene islands.

## Results

- Clusters are V-gene islands: the dominant TRBV gene holds 63% of a cluster's clonotypes on
  average. Stability is moderate (Jaccard 0.47), so treat cluster identities as approximate.
- Between clusters, shorter mean CDR3 goes with higher Pgen and more sharing (panel b), the
  same pattern as report 03.
- **Enrichment:**
  - 5,526 clonotypes exactly match VDJdb.
  - Of 364 tests, 169 have BH q < 0.05: **81 enrichments and 88 depletions**.
  - 126 tests are significant within TRBV.
  - 89 are significant both ways, of which **50 are enriched both ways**.
- CMV and EBV enrichments hold within TRBV (top CMV cluster: crude OR 7.8, within-TRBV q 5e-151;
  EBV: 8.4, q 8e-108).
- The strongest influenza A hit (cluster 50, OR 29.8) and self-peptide hit (cluster 12, OR
  23.8) are not significant within TRBV (q 0.35 and 0.19). They are V-gene effects.
- Self-peptide TCRs also concentrate in some clusters. Clusters capture sequence motifs, not
  shared targets.

## Correction (2026-10-20)

- The earlier log said "169 significant (q < 0.05)" without noting that 88 of these are
  depletions. Enrichment claims should use the 81 (crude) or 50 (both) counts.
- Within-TRBV odds ratios: the Mantel–Haenszel estimate is degenerate when strata have empty
  cells. It is infinite for 65 of 364 tests and above 1,000 for 13 more (e.g. CMV cluster 42:
  26,497). Use `p_cmh_trbv` / `q_cmh_trbv`; do not report the MH odds ratio as an effect size.
  Panel d's diamonds inherit this, so the far-right diamond and the missing ones are artifacts.
- "Synthetic" is VDJdb's label for engineered TCRs, not a pathogen. Its 44 matches all fall in
  cluster 50, the influenza A cluster. It should be excluded or labelled as such in a redraw.

## Figure

`fig_embedding_clusters`:
- a: centroid map.
- b: cluster properties.
- c: volcano of all tests.
- d: strongest enrichments, crude (dots) and within TRBV (diamonds).

## Files

`summary.csv`, `vdjdb_cluster_enrichment.csv` (every test, crude and within-TRBV),
`clusters_3d.csv` (centroids in 3-d PCA with per-cluster statistics; feeds the interactive
map).
