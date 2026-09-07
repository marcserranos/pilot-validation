# Why the UMAP/PCA manifold shows ancestry but not disease

Open research question raised 2026-09-06: `08_embedding_compare.py`'s unsupervised embeddings show
clear ancestry structure (as documented in `06_figures_structure.py` and repeated here) but no
visible clustering when colored by real EHR-confirmed diagnosis
(`disease_real_diagnosis.png`) -- even though `13_disease_allele_association.py` finds a real,
strong, statistically robust association (B\*27 / ankylosing spondylitis: OR=9.5 ancestry-adjusted,
p=4.3e-10 CMH). This doc is the working answer to "why," plus the experiments it implies.
`14_manifold_structure.py` implements items 1-4 and 6; 5 is documented as future work.

## The mechanism (working hypothesis)

Unsupervised embeddings (PCA, UMAP) find directions of maximal **shared covariance** across the
whole feature matrix -- they have no notion of an external label at all. Ancestry and disease status
differ enormously in how much of the ~7,000-column allele-dosage matrix they touch:

- **Ancestry is a matrix-wide, rank-many signal.** HLA is the most differentiated region in the
  human genome (highest Fst loci known, maintained by balancing selection) and its classical genes
  sit in strong extended LD -- ancestry shifts correlated blocks across most of the 8 genes at once.
  That is exactly the kind of low-rank, high-magnitude, matrix-wide covariance PCA/UMAP are built to
  find.
- **A disease association like B\*27/AS is a rank-1 signal in a single column** (or a small LD-linked
  set), and the outcome itself (diagnosis) isn't a matrix column at all -- it's an external label
  with no reason to correlate with the *other* 7,399 columns for any given carrier. 423 B\*27
  carriers aren't similar to each other anywhere else in allele-space, so there's no reason for them
  to sit near each other in an embedding whose axes are defined by aggregate variance.
- **Standardization doesn't cause this asymmetry, it just reweights everything toward rarer
  variants** (see `08_embedding_compare.py`'s own docstring on why PC1/PC2 explain <1% var at this
  encoding). It doesn't change *which* alleles are ancestry-correlated vs. disease-correlated.

This is the manifold-space analog of the standard GWAS fact that genome-wide PCA is used to
*correct for* population stratification, not to *find* single-locus hits -- unsupervised embeddings
are the wrong instrument for sparse single-allele effects, independent of which 8-11 diseases we
picked. `13_disease_allele_association.py`'s Fisher/CMH test is the right instrument, and it already
found the signal.

## Experiments (this doc's checklist -- update status as they land)

1. **[DONE] Quantify the needle-in-haystack claim.** From the existing `embedding_cache.pkl`
   (no re-fit needed): raw and standardized variance fraction per disease-proxy allele column, its
   loading magnitude on PC1-PC3, and an ancestry-differentiation score (between/within-ancestry
   variance ratio) per column -- then correlate loading magnitude against the ancestry score across
   all columns. If disease-relevant columns have near-zero PC1/PC2 loading while
   ancestry-differentiated columns dominate, that's a direct confirmation of the mechanism above,
   not just a plausible story.
2. **[DONE] Residualize out ancestry, then re-embed.** Subtract each column's per-ancestry-group
   mean (removing between-ancestry variance only) and re-fit PCA + UMAP on the residual matrix,
   colored by real diagnosis / disease burden. Tests directly: is there a *second-order* manifold
   that only appears once the dominant ancestry axis is removed?
3. **[DONE] Supervised UMAP.** `umap-learn` supports `y=` supervision (`target_metric`). Fitting
   UMAP with disease-diagnosis (or burden score) as the target asks a sharper question than
   unsupervised UMAP ever can: does *any* combination of the allele-dosage space separate cases from
   controls, even one an unsupervised fit would never surface on its own.
4. **[DONE] Replace "no visible clustering" with a number.** A KNN-label-enrichment permutation
   test (fixed neighbor graph per embedding, labels permuted) run on all three embeddings x two
   label sets (ancestry as a positive control that should show strong enrichment; disease
   burden/diagnosis as the real test). Turns the visual impression into a p-value, and is reusable
   for #2 and #3's outputs.
5. **[NOT IMPLEMENTED -- future work] Haplotype-block features instead of per-gene dosage.** If
   extended ancestral haplotypes (e.g. the 8.1 AH, which itself carries strong autoimmune-disease
   risk) carry the ancestry signal *and* the disease risk together as a block, per-gene dosage
   columns may be diluting a haplotype-level signal that a haplotype-identity feature space would
   recover. Needs cis-phased multi-gene haplotype calling from Immuannot's per-contig output (see
   `scripts/hla_popgen/README.md`'s point 1 on `(person_id, hap, contig)` cis-phasing) -- real
   additional engineering, not a quick add. Worth revisiting if #2/#3 come back interesting.
6. **[DONE] Composite disease-burden score instead of single-disease binary flags.** Each individual
   disease has thin `n_diagnosed` (24-253 in the real run), which is weak power for "does this
   cluster" on its own. `14_manifold_structure.py` adds a continuous burden score (count of
   `HLA_LINKED` diagnoses matched, 0-11) as an additional coloring/label option, trading
   disease-specificity for sample size before concluding "no structure" disease-by-disease.

## Result (fill in after the VM run)

See `reports/hla_popgen/14_manifold_structure/manifold_structure_report.md` once run.
