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

## Result (real run, 9,355 people x 5,453 allele-identity columns, 853 [9.1%] with >=1 HLA-linked
diagnosis)

Full numbers: `reports/hla_popgen/14_manifold_structure/manifold_structure_report.md`,
`variance_loadings.tsv`, `knn_label_enrichment.tsv`. Figure:
`reports/hla_popgen/14_manifold_structure/manifold_structure_embeddings.png`.

**The mechanism is confirmed, and the open question has a clean answer.** KNN-label-enrichment
(k=15, 500 permutations, one-sided p):

| embedding | ancestry (positive control) | any HLA-linked diagnosis |
|---|---|---|
| locked_raw (unsupervised) | obs=0.529 vs null=0.205, **p=0.002** | obs=0.834 vs null=0.834, p=0.64 (n.s.) |
| ancestry_residualized | obs=0.495 vs null=0.205, **p=0.002** | obs=0.834 vs null=0.834, p=0.74 (n.s.) |
| supervised (y=any_diag) | obs=0.516 vs null=0.205, **p=0.002** | obs=0.970 vs null=0.834, **p=0.002** |

1. **Unsupervised embeddings never find disease structure, residualized or not** -- confirms the
   mechanism write-up above directly: ancestry is a matrix-wide covariance signal UMAP/PCA are built
   to find; disease status isn't in the matrix at all, so there's no reason for it to emerge on its
   own, at any resolution we tried.
2. **The information is there -- supervised UMAP finds it cleanly** (obs=0.97 vs a 0.83 null,
   visually two separated clusters in the figure's third panel, sized about right for a 9.1%-prevalence
   split). This resolves the original question: the allele-dosage space is not *blind* to disease
   status, unsupervised methods are just the wrong instrument for a signal this sparse relative to
   the dominant ancestry axis -- exactly the GWAS-PCA analogy in the mechanism section above, now with
   a positive result to back it up rather than just an absence of a negative one.
3. **Caveat, stated plainly:** `umap.UMAP(..., y=labels)` is explicitly optimized to pull same-label
   points together, so a high same-label KNN score on the *training* fit is expected almost by
   construction and is not by itself proof the signal generalizes -- it demonstrates the information
   is extractable, not that it is strong out-of-sample. See "New hypothesis" #1 below for the direct
   fix (held-out validation).
4. **Unplanned finding: mean-only ancestry residualization barely moves the ancestry KNN score**
   (0.529 raw -> 0.495 residualized, both far above the 0.205 null). Subtracting each column's
   per-ancestry-group mean removes *additive* ancestry signal but evidently leaves most of the
   ancestry-driven structure intact. The likely explanation: ancestry differentiates HLA allele
   *co-occurrence* (linkage disequilibrium), not just individual allele frequencies -- classical HLA
   genes sit in strong, population-differentiated LD, so the *correlation structure* between columns
   differs by ancestry group even after each column's own mean is matched. A first-order (mean)
   residualization can't remove a second-order (covariance) effect. Not something the original
   checklist anticipated; worth its own line of investigation (see below).
5. **Variance/loadings audit is consistent but non-obvious in its own right:** rank correlation
   between |PC1 loading| and per-column ancestry F-stat is 0.395 (positive, moderate, across all
   5,453 columns) -- loadings do track ancestry differentiation, supporting the mechanism. Disease-proxy
   columns (145 of them) actually have a *higher* mean ancestry F-stat (10.01) than the rest of the
   matrix (6.25), i.e. well-known disease alleles (B\*27, DQB1\*03:02, DRB1\*15:01, ...) are
   themselves fairly ancestry-differentiated -- unsurprising given HLA population genetics, but a
   reminder that "disease-associated" and "ancestry-neutral" are not the same axis, and that the
   B\*27/AS association in `13_disease_allele_association.py` needed ancestry-adjustment (CMH) for
   exactly this reason.

## New hypotheses / lines of investigation opened by this result

1. **[DONE, RESOLVED NEGATIVE -- 2026-09-08] Held-out validation of the supervised-UMAP signal.**
   `15_manifold_holdout_validation.py`: stratified 70/30 train/test split (6,548/2,807, prevalence
   preserved), UMAP fit (supervised AND an unsupervised baseline) on train only using train-only
   standardization stats, test split `.transform()`-projected out-of-sample, k-NN-in-embedding
   scoring against TRAIN neighbors only, evaluated with AUROC + bootstrap CI + a label-permutation
   p-value on the held-out test set. Full results:
   `reports/hla_popgen/15_manifold_holdout_validation/holdout_validation_report.md`, figure
   `holdout_train_test_embedding.png`.

   **Verdict: the caveat was right -- the original supervised-UMAP finding does NOT generalize.**
   Held-out AUROC: supervised **0.482** (bootstrap 95% CI 0.460-0.508, permutation p=0.93 -- not
   significant, indistinguishable from chance); unsupervised baseline **0.528** (CI 0.494-0.563,
   p=0.061 -- also not significant, though a touch closer to the edge). Supervised does not clear
   0.5 and does not beat the unsupervised baseline -- both conditions in the pre-registered verdict
   logic for "this is real" fail.

   The figure makes *why* visually obvious: supervised UMAP carved out a small, tight, isolated
   pocket containing almost exclusively TRAIN diagnosed people (bottom-right of the left panel) --
   but held-out diagnosed people (teal triangles) essentially never land in that pocket; they're
   scattered through the main blob indistinguishably from held-out undiagnosed people. This is
   textbook `umap.UMAP(y=...)` overfitting: given a label, it can trivially wall off the exact
   training points that carry it into their own bubble without learning any transferable direction
   in allele-dosage space -- a classic manifold-learning failure mode, not a data problem.

   **Updated conclusion for the whole doc:** `14_manifold_structure.py`'s "supervised UMAP finds
   real separation" result (step 3 above) is retracted as evidence of disease-relevant structure.
   The honest state of the evidence is back to "no disease structure found in the raw allele-dosage
   embedding space by any method tried here, supervised or not" -- `13_disease_allele_association.py`'s
   direct Fisher/CMH test on individual alleles remains the only instrument in this whole project
   that has found real, validated disease signal (B\*27/ankylosing spondylitis). That is not a
   failure of this investigation -- ruling out an overfit positive result via honest held-out
   testing is exactly what this line of work was for, and the mechanism argument at the top of this
   doc (why unsupervised embeddings can't find sparse single-locus effects) is untouched by this
   result.
2. **Test the LD-not-just-frequency explanation for finding #4 directly.** Residualize each column
   against ancestry using a *covariance-aware* method (e.g. per-ancestry-group full standardization,
   not just mean-centering; or fit ancestry-specific PCA and compare loadings) rather than a mean
   shift, and see whether the ancestry KNN score finally drops toward null. If it does, that
   confirms LD-structure-by-ancestry as the driver, not lingering mean effects; if it still doesn't
   drop, there is a real open question about what unsupervised PCA/UMAP is actually tracking here
   that residualization can't reach. Not yet implemented.
3. **Haplotype-block features (checklist item 5, unchanged from before).** Still the most direct test
   of whether disease risk travels with ancestry along shared extended haplotypes (e.g. the 8.1 AH) --
   now more motivated by finding #4 (ancestry-driven LD structure) than it was originally. Real
   engineering lift (cis-phasing), not yet started.
4. **Per-disease supervised UMAP instead of the pooled `any_diag` burden flag.** The current
   supervised result pools all 11 HLA-linked diagnoses into one binary target; worth checking whether
   the separation is driven by a few common diagnoses (T1D n=204, psoriasis n=223, RA n=253) or holds
   up per-disease even at low n (e.g. does B\*27/ankylosing-spondylitis alone, n=28, show any
   supervised separation, matching its very strong Fisher/CMH result from `13_disease_allele_association.py`?). Not yet implemented.
