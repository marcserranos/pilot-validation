# HLA Population-Genetics Visualization Literature Review & Figure Catalog

**Purpose**: drive the plotting scripts for the Omni-HLA project (AoU-native short-read calls, N≈500k, 30 genes, 2–3 field; long-read Immuannot calls, N≈12k, phased, up to 4-field, includes non-classical genes; high-confidence subset filtered on `template_distance`; ancestry as both discrete `ancestry_pred` and continuous 6-way `probabilities` + `pca_features`).

Every citation below is a real, checked source (title + URL). Where I could not verify a specific numeric or design claim from the fetched text, it is marked **UNVERIFIED**.

---

## Part 0 — Data dictionary assumed throughout

| Column | Description |
|---|---|
| `person_id` | participant |
| `cohort` | {aou_native, immuannot_all, immuannot_hiconf} |
| `gene` | HLA-A, -B, -C, -DRB1, ..., MICA, MICB, TAP1, TAP2, DRB3/4/5, DQA2, DQB2, HLA-E/F/G, pseudogenes |
| `allele_1field`, `allele_2field`, `allele_3field`, `allele_4field` | nested nomenclature fields |
| `haplotype_id` / `phase_set` | which of the 2 chromosomal copies (long-read only) |
| `template_distance` | edit distance to nearest IPD-IMGT/HLA reference allele (0 = exact match; long-read only) |
| `ancestry_pred` | AFR/AMR/EAS/EUR/MID/SAS hard label |
| `probabilities` | length-6 vector, same order, continuous admixture proportion per participant |
| `pca_features` | genotype PCs |

---

## Part 1 — Canonical figures in HLA population genetics

### 1.1 Grouped/stacked allele frequency bar chart

**What it shows**: for one gene, frequency of each allele group (usually 2-field) as bars, grouped by population/ancestry, often stacked to sum to 1 (since alleles are mutually exclusive per haplotype).

**When right**: small number of alleles (<15) per gene, few populations (<10) — the default first figure in almost every HLA typing/registry paper.

**Failure modes**: unreadable past ~8 populations x 15 alleles; legend explosion; can't show CIs cleanly on a stacked bar; rare alleles get lost visually.

**Implementation sketch** (matplotlib/seaborn):
```python
# long-format df: gene, allele_2field, ancestry_pred, freq, n
pivot = df.pivot(index='ancestry_pred', columns='allele_2field', values='freq')
pivot.plot(kind='bar', stacked=True, colormap='tab20')
```

**Precedent**:
- Population-specific classical HLA allele/haplotype tables, e.g. the 17th IHIW joint report on high-resolution allele/haplotype frequencies — https://pmc.ncbi.nlm.nih.gov/articles/PMC8315142/
- "Classical HLA Allele and Haplotype Frequency Estimates in US Populations" (bioRxiv) — grouped bars by 5 broad + 21 detailed US populations — https://www.biorxiv.org/content/10.64898/2026.04.09.717537v1.full

### 1.2 Allele × population frequency heatmap (double-clustered)

**What it shows**: rows = alleles, columns = populations, cell color = frequency; both axes hierarchically clustered so that similar populations and co-varying alleles group together.

**When right**: many populations (tens–hundreds) and/or many alleles — this is the workhorse figure for "HLA landscape" papers.

**Failure modes**: linear color scale hides rare-allele signal (use log or perceptual sqrt scale); double clustering can be unstable with sparse/rare-allele rows — consider clustering only common alleles (freq > 1%) and appending a separate "rare" block.

**Implementation sketch**: `seaborn.clustermap(freq_matrix, method='average', metric='euclidean', cmap='viridis', z_score=None, standard_scale=None)`; use `robust=True` to clip outliers.

**Precedent**:
- "An HLA map of the world: A comparison of HLA frequencies in 200 worldwide populations" (Frontiers in Genetics 2023) — explicitly uses "unsupervised double clustering-based heat maps" of allele-group frequencies for HLA-A/B/DRB1 across 200 populations, plus dendrograms from minimum-variance clustering of squared Euclidean distances between populations — https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2023.866407/full (PMC mirror: https://pmc.ncbi.nlm.nih.gov/articles/PMC10076764/)
- HLA Diversity in the 1000 Genomes Dataset (PLOS ONE 2014) — https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0097282

### 1.3 Rank-frequency / allele-frequency-spectrum curve

**What it shows**: alleles ranked by frequency (x, log or rank) vs frequency (y, log) per population/cohort — a "Zipf plot" for HLA, analogous to the site-frequency spectrum in population genetics.

**When right**: comparing overall diversity/evenness structure between cohorts irrespective of allele identity; great for showing the AoU-native (2-3 field, fewer distinct calls) vs Immuannot (4-field, long tail) resolution difference in one panel.

**Failure modes**: obscures which specific allele is where; ties in rank need jittering; log-log can visually flatten real differences — always overlay N to contextualize sampling noise on the tail.

**Implementation sketch**: sort per cohort, `plt.step`/`plt.loglog(rank, freq)`, one line per cohort per gene, small multiples across genes.

**Precedent**: motivated directly by the extreme skew documented in HLA allele-count studies — "the extremely large number of known HLA alleles ... creates a formidable challenge" — thousands of low-frequency haplotypes vs a handful of common ones, per Mapping Bias paper — https://academic.oup.com/g3journal/article/5/5/931/6025555 . General rank-frequency plotting is standard in population genetics SFS analysis (methodological precedent, not HLA-specific — **UNVERIFIED** as an HLA-specific figure in the literature searched).

### 1.4 Ternary (triangle) plot / De Finetti diagram

**What it shows**: for 3-population (or 3-ancestry-component) systems, each individual or allele-frequency vector as a point inside a triangle whose vertices are the pure ancestral states.

**When right**: exactly 3 ancestry components dominate variance (e.g., AFR/EUR/AMR admixture in a Latino cohort) or 2-allele/diploid genotype frequencies (classic De Finetti use).

**Failure modes**: doesn't generalize past 3 components (we have 6 ancestries) — must either pick a dominant 3-way subset per analysis (e.g., "AMR admixture continuum" = AFR/EUR/AMR only) or use a different projection for the full 6-way space (see 2.3 below).

**Implementation sketch**: `python-ternary` or manual barycentric transform + matplotlib; color points by carrier status of a target allele.

**Precedent**:
- De Finetti diagram, foundational population-genetics ternary plot — https://en.wikipedia.org/wiki/De_Finetti_diagram
- "The geometry of admixture in population genetics: the blessing of dimensionality" — formalizes that admixed individuals with k parental populations lie in a (k-1)-simplex — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11639143/
- mapmixture R package/webapp for spatial admixture visualization (pie/ternary-style summaries) — https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13943

### 1.5 ADMIXTURE/STRUCTURE-style stacked "barcode" plot

**What it shows**: one vertical stacked bar per individual, segments colored by ancestry-component proportion, individuals sorted/grouped by predicted ancestry — the ubiquitous population-genetics figure.

**When right**: showing that continuous ancestry proportions (our `probabilities` vector) are real admixture and not just label noise; good companion panel next to any ancestry-stratified HLA analysis to justify why continuous treatment matters (visually — many "EUR"-labeled participants will show non-trivial AFR/AMR mass).

**Failure modes**: sorting strategy dominates the story (sort within each predicted label by dominant proportion, not randomly); overplotting with 500k people — use a downsampled/binned representation or 2D density instead of literal per-person bars past ~5,000 individuals.

**Implementation sketch**: `matplotlib.bar(x=range(n), height=1, bottom=cumsum, color=ancestry_colors)`, sort by `ancestry_pred` then descending value of that component.

**Precedent**: AncestryPainter 2.0 "sectorplot"/composition visualization for ancestry matrices — https://academic.oup.com/gbe/article/16/11/evae249/7900898 (PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC11604083/); general admixture stacked-bar precedent from PCA/ADMIXTURE model review — https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_4

### 1.6 PCA / UMAP of HLA genotype or allele dosage

**What it shows**: dimensionality reduction of per-person HLA allele-dosage vectors (one-hot or count-coded across all observed alleles at a locus, or across all loci) to see whether HLA genotype alone recovers ancestry structure.

**When right**: demonstrating HLA-driven population structure independent of genome-wide PCs; a nice validation that our ancestry labels are consistent with HLA-only signal; also useful diagnostic for cohort-specific batch effects (AoU-native vs Immuannot calling pipeline artifacts).

**Failure modes**: PCA linear components often fail to separate >2-3 clusters cleanly in 2D for MHC SNP data; UMAP not deterministic without a fixed seed and can create spurious separations — always report parameters (n_neighbors, min_dist) and pair with a PCA panel for honesty.

**Implementation sketch**: build a person × allele indicator/count matrix (sparse), reduce via `sklearn.decomposition.PCA` and `umap.UMAP`, plot colored by continuous ancestry proportion (see 2.3).

**Precedent**:
- "UMAP reveals cryptic population structure and phenotype heterogeneity in large genomic cohorts" (PLOS Genetics) — explicitly notes UMAP on MHC/SNP data gave good ancestry separation "concordant with the frequency difference of HLA alleles between populations" where 2D PCA failed — https://journals.plos.org/plosgenetics/article?id=10.1371%2Fjournal.pgen.1008432
- "Quantitative evaluation of nonlinear methods for population structure visualization" (G3) — https://academic.oup.com/g3journal/article/12/9/jkac191/6651067

### 1.7 Fst / genetic-distance matrix and dendrogram between populations

**What it shows**: pairwise population differentiation (Fst, or Euclidean/chord distance on allele-frequency vectors) as a heatmap matrix and/or hierarchical dendrogram.

**When right**: summarizing whether HLA differentiation between AoU ancestry groups matches expectations from genome-wide Fst, and whether the AoU-native vs Immuannot cohorts give concordant between-ancestry distances (a QC/validation figure, not just descriptive).

**Failure modes**: Fst estimators are sensitive to sample-size imbalance (badly needed here: AFR/AMR/EAS/MID/SAS will be much smaller N than EUR in AoU) — use a bias-corrected estimator (Weir & Cockerham) and bootstrap CIs on branch lengths.

**Implementation sketch**: `scikit-allel` or manual Weir-Cockerham Fst per gene per population pair → `scipy.cluster.hierarchy.linkage` + `dendrogram`; matrix heatmap with `seaborn.heatmap(annot=True)`.

**Precedent**:
- "An HLA map of the world" — minimum-variance clustering dendrograms of squared Euclidean distances on HLA-A/B/DRB1 frequencies across 200 populations — https://pmc.ncbi.nlm.nih.gov/articles/PMC10076764/
- Higher HLA-locus Fst than genome-wide Fst between same-continent populations is a documented and citable pattern (from search synthesis of HLA diversity literature) — **UNVERIFIED exact source**; treat as a hypothesis to test with our own data rather than a citation-backed number until traced to a primary paper (candidates: Sanchez-Mazas group publications).

### 1.8 Haplotype / linkage (cis-pair) frequency plots between HLA loci

**What it shows**: two-locus haplotype frequencies (e.g., DRB1~DQB1) as a bipartite/alluvial diagram or a frequency table/heatmap of locus-A-allele × locus-B-allele combinations, sometimes annotated with D' or Δ linkage-disequilibrium statistics.

**When right**: essential once we have true haplotype (cis) phase from long-read data — most published HLA papers can only *estimate* haplotypes statistically (EM algorithm via tools like PyPop/Arlequin); we can report them directly, which is a strong point of novelty (see Part 3).

**Failure modes**: statistically-imputed two-locus haplotypes (used in nearly all comparator papers) carry phase uncertainty we don't have — don't present our phased vs literature's imputed haplotype frequencies on the same axis without flagging the methodological asymmetry.

**Implementation sketch**: alluvial/Sankey plot (`plotly.graph_objects.Sankey`) linking DRB1 allele nodes to DQB1 allele nodes, edge width = haplotype frequency, edges colored by ancestry.

**Precedent**: DRB1~DQB1 two-locus haplotype and LD reporting is standard across regional HLA frequency papers, e.g. Lebanese population report (D' = 0.80 for DRB1~DQB1) — https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08682-7 ; Bahraini population report of top two-locus haplotype frequency — https://www.sciencedirect.com/science/article/abs/pii/S0378111920300688 ; bimodal LD sign pattern (positive among frequent, negative among rare haplotypes) documented across populations — synthesized from search, primary source not individually confirmed, **UNVERIFIED** for exact paper attribution.

### 1.9 Sunburst / treemap of allele nomenclature hierarchy

**What it shows**: nested hierarchy 1-field → 2-field → 3-field → 4-field, area or arc length proportional to frequency at each level.

**When right**: showing how a single 2-field call (all most databases and the AoU-native pipeline resolve to) "fans out" into multiple 4-field subtypes once you have long-read resolution — an excellent, intuitive one-glance summary of what extra resolution buys you.

**Failure modes**: sunbursts get visually noisy beyond ~50 leaf segments — restrict to top-N alleles by frequency, bucket the remainder into "other."

**Implementation sketch**: `plotly.express.sunburst(df, path=['gene','field1','field2','field3','field4'], values='count')`.

**Precedent**: general applicability of sunburst/treemap for hierarchical categorical data is well established in visualization literature (e.g. Hi-d maps arXiv paper on multi-dimensional categorical hierarchies — https://arxiv.org/pdf/2507.07890), but no HLA-specific published sunburst figure was found in this search — this is a genuinely underused technique in the HLA literature (**UNVERIFIED as prior HLA-specific use** — treat as a novel/underused contribution, not as replicating precedent). The WHO 4-field nomenclature structure itself is well documented, e.g. https://www.researchgate.net/figure/HLA-nomenclature-with-allele-resolution-at-four-fields-eight-digits-Commonly-the-HLA_fig1_328519906 and https://cytologicsbio.com/how-to-decipher-human-leukocyte-antigen-hla-nomenclature/

### 1.10 Diversity indices by population (heterozygosity, Shannon entropy, allelic richness)

**What it shows**: bar or box plot of per-population, per-gene diversity summary statistics; box plot of "number of allele groups detected" per continental group with pairwise significance is a direct precedent.

**When right**: the standard population-genetics summary that lets readers compare our cohort's diversity against published worldwide figures without needing full frequency tables.

**Failure modes**: heterozygosity and allelic richness are both sample-size sensitive — rarefaction (subsampling to a common N) is required for fair comparison across ancestries/cohorts of very different sizes (AFR N will be much smaller than EUR in AoU).

**Implementation sketch**: compute expected heterozygosity He = 1 - Σp_i², Shannon H' = -Σ p_i ln(p_i), rarefied allelic richness via `scikit-bio` or manual rarefaction curves; boxplot per ancestry with N annotated.

**Precedent**: "An HLA map of the world" Figure 5 — box plots of number of HLA allele groups detected across continental groups with p-value matrices for pairwise comparison — https://pmc.ncbi.nlm.nih.gov/articles/PMC10076764/ ; progressive reduction in HLA diversity from African to Native American populations, a widely cited pattern in HLA diversity reviews (serial-founder-effect signature) — synthesized from search results, exact primary citation **UNVERIFIED** in this pass (likely traces to Solberg et al. 2008 or Sanchez-Mazas — recommend direct follow-up read of the "17th IHIW" report: https://pmc.ncbi.nlm.nih.gov/articles/PMC8315142/).

---

## Part 2 — Continuous-ancestry approaches (differentiator)

Almost every HLA population paper found in this search bins participants into discrete continental labels before computing frequencies. This is the norm even in multi-ancestry UK Biobank HLA work (Communications Biology 2023, https://www.nature.com/articles/s42003-023-05496-5 — analyzed by 5 discrete UKB ancestry groups). We have continuous 6-way admixture (`probabilities`) plus genotype PCs, which is a genuine point of departure.

### 2.1 Allele frequency as a smooth function of ancestry proportion

**Claim supported**: "Allele X's population frequency changes continuously and monotonically with AFR (or AMR, etc.) admixture fraction, rather than jumping discretely between labeled groups" — directly demonstrates that discrete-bin frequency estimates are a coarse approximation and that admixed individuals' expected allele frequency is a mixture, not a third category.

**Data needed**: `person_id`, `probabilities[ancestry]` (continuous, one component at a time), allele dosage (0/1/2 copies) of the target allele.

**Implementation sketch**: local regression / GAM of dosage/2 (or carrier indicator) on ancestry-component proportion: `statsmodels.nonparametric.lowess` or `pygam.LogisticGAM`; bin into deciles of ancestry proportion with Wilson CIs per bin as a robustness check against the smooth curve; small multiples per gene × top alleles.

**Precedent**: this is exactly the logic of ancestry-specific allele-frequency estimation methods built for admixed populations —
- AFA: "Ancestry-specific allele frequency estimation in admixed populations" (Hispanic Community Health Study/SOL) — models frequency as a function of continuous ancestry proportions via maximum likelihood — https://www.sciencedirect.com/science/article/pii/S2666247722000124
- LEI: allele-frequency-based feature selection for multi-ancestry admixed populations — https://www.nature.com/articles/s41598-019-47012-y

### 2.2 Regression of allele dosage on ancestry fraction ("ancestry-specific AF" model)

**Claim supported**: quantifies, per allele, the marginal AF attributable to each of the 6 ancestral components — lets us report "HLA-B*07:02 has estimated ancestral frequency 0.14 in the EUR component vs 0.01 in the AFR component" even for individuals who are admixed, which a discrete-label approach cannot do without throwing away partially-admixed people.

**Data needed**: full `probabilities` 6-vector per person (not just argmax label), allele dosage per gene/allele.

**Implementation sketch**: for a biallelic dosage y_i (0/1/2) and ancestry proportion vector q_i (sums to 1), fit constrained linear/logistic model y_i ~ Σ_k β_k q_ik (no intercept, β_k ∈ [0,1] interpretable as ancestral AF); report β_k with bootstrap CIs; this is literally the AFA/LAI-based ancestry-specific-AF approach adapted to a 6-way vector instead of 3-way.

**Precedent**:
- AFA (as above) — https://www.sciencedirect.com/science/article/pii/S2666247722000124, preprint https://www.biorxiv.org/content/10.1101/2021.08.06.455462.full.pdf
- Validation approach of comparing ADMIXTURE-derived vs RFMix local-ancestry-derived proportions (r=0.98 in HCHS/SOL) as a template for how to sanity-check our own continuous ancestry proportions before trusting them as regressors — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5695820/

### 2.3 Ancestry-continuum scatter/density plot (2D ancestry-probability projection colored by allele carrier status)

**Claim supported**: visually shows whether allele carriers cluster in ancestry-probability space even when their discrete `ancestry_pred` label is homogeneous — reveals cryptic admixture-driven signal invisible to a bar chart by label.

**Data needed**: 2 of the 6 `probabilities` components (or a PCA/UMAP of the full 6-vector), carrier status (0/1) or dosage of a specific allele, `ancestry_pred` as an outline/shape aesthetic.

**Implementation sketch**: 2D KDE contours of carriers vs non-carriers overlaid on a scatter of all points in ancestry-probability space (e.g. AFR-fraction vs EUR-fraction plane), or hexbin of carrier *rate* per ancestry-probability cell.

**Precedent**: analogous to plotting inferred local/global ancestry proportions from two methods against each other for validation (ADMIXTURE vs RFMix scatter, r=0.98) — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5695820/ ; conceptually the "ancestry continuum" idea is the visualization-side counterpart of admixture-mapping methodology reviewed for the Latino admixture map — https://www.sciencedirect.com/science/article/pii/S0002929707610222

### 2.4 Ternary plot colored by allele carrier status (3-component subset)

**Claim supported**: for cohorts/alleles where a 3-way ancestry decomposition dominates (e.g., AMR individuals decomposed into AFR/EUR/AMR-indigenous proportions), shows spatially within the simplex where carriers concentrate — a more geometrically faithful version of 2.3 for exactly-3-component subsets.

**Implementation sketch**: as in 1.4, but color/size-encode by carrier status or local allele dosage; consider a hexbin-in-simplex carrier-rate surface rather than raw scatter once N is large.

**Precedent**: De Finetti/ternary admixture geometry as above (1.4); genotype-frequency ternary diagrams are the direct ancestor of this idea in classical population genetics — https://en.wikipedia.org/wiki/De_Finetti_diagram

---

## Part 3 — Novel / out-of-the-box visualizations exploiting our unique cohort structure

These exploit things essentially no published HLA paper has simultaneously: three cohorts of the same participants (well, overlapping) at different resolution/technology, a template-distance accuracy axis, true haplotype phase, and continuous ancestry.

### 3.1 "Calling bias" scatter: same-allele frequency, cohort A (AoU-native) vs cohort B (Immuannot)

**Claim it supports**: short-read ensemble calling (HLA-HD+Polysolver+OptiType) systematically over- or under-calls specific alleles relative to long-read ground truth — a direct, visual accuracy/bias diagnostic, and would be the single most citation-worthy new figure in the paper.

**Data needed**: allele frequency computed independently in each cohort (restricted to the overlapping/shared participant set if available, or the same ancestry-matched subpopulations if not), joined on `gene` + `allele_2field` (coarsest common resolution).

**Implementation sketch**: scatter, x = freq in AoU-native, y = freq in Immuannot, log-log, 1:1 reference line, point size = combined allele count, color = gene or ancestry-majority; label systematic outliers (alleles that deviate >2x). Facet by gene.

**Precedent for the *type* of comparison** (not identical figure, but the same logic of comparing sequencing-technology-derived AF to a reference/high-confidence set): Mapping Bias paper on 1000 Genomes HLA overestimating reference allele frequency due to short-read mapping bias — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4426377/ (also https://academic.oup.com/g3journal/article/5/5/931/6025555) is exactly the phenomenon we'd be re-demonstrating with a much better long-read comparator.

### 3.2 template_distance distribution by ancestry — reference-bias evidence plot

**Claim it supports**: IPD-IMGT/HLA reference database itself is ancestry-biased (built predominantly from European-ancestry donors historically), so participants of non-European ancestry should show systematically higher `template_distance` (more novel/divergent alleles relative to reference) — this is a strong, novel, testable claim unique to having a phased long-read dataset with a distance-to-reference field.

**Data needed**: `template_distance`, `ancestry_pred` or continuous `probabilities`, `gene`.

**Implementation sketch**: violin or ridgeline plot of `template_distance` by ancestry per gene; complement with a continuous version — mean `template_distance` vs continuous ancestry-component fraction (loess), analogous to 2.1; also report % of alleles at exact distance-0 by ancestry, since your project notes already flag distance-0 threshold usability issues on real data (per your own recent commit "Fix filter-first experiment: exact-distance-0 threshold was unusable on real data") — worth a companion histogram of the raw distance distribution, not just a hard 0-cutoff.

**Precedent**: the underlying premise — reference/database ascertainment bias toward European-ancestry alleles — is a well-documented general genomics phenomenon; the closest direct HLA analogue found is the reference-bias mapping literature (1000G HLA overestimation of reference AF, above) and general reference-bias reviews — https://www.sevenbridges.com/reference-bias-challenges-and-solutions/. The specific figure (template_distance × ancestry) does not have a direct published precedent found in this search — **mark as a genuinely novel figure**, motivated by analogy rather than replicated from a prior paper.

### 3.3 Phased haplotype (cis-pair) frequency plot exploiting true phase

**Claim it supports**: we can report *directly observed* cis-haplotype frequencies (e.g., DQA1~DQB1, or extended DRB1~DQA1~DQB1~DPB1) without EM-based statistical phasing, and can quantify how much the EM-imputed haplotype frequencies from unphased short-read (AoU-native) data disagree with our ground-truth phase — a second bias/validation figure paralleling 3.1 but at the haplotype level.

**Data needed**: `haplotype_id`/`phase_set`, allele calls per gene per haplotype (long-read only) vs EM-estimated haplotype frequencies computed from the AoU-native unphased genotypes (would need to run an EM step, e.g. via PyPop-style expectation-maximization, as a comparator).

**Implementation sketch**: Sankey/alluvial diagram (as in 1.8) drawn twice side by side — "true phased" vs "EM-imputed from short-read" — with discordant edges highlighted; or a scatter of true vs EM-imputed two-locus haplotype frequency (same idea as 3.1, one level up).

**Precedent**: EM-based two-locus haplotype/LD estimation from unphased genotypes is standard (PyPop, Arlequin), documented across many regional HLA frequency papers reviewed above (Lebanese, Bahraini, Serbian reports); comparing EM-estimated to truly phased haplotypes is the natural validation study but was not found as a published HLA-specific figure in this search — **novel figure**, methodologically motivated by the general phasing-accuracy literature rather than a direct HLA precedent.

### 3.4 Resolution cascade: 2-field → 4-field fan-out by ancestry

**Claim it supports**: for a fixed 2-field allele group (the resolution ceiling of the AoU-native short-read pipeline and of most public frequency databases), the distribution of 4-field subtypes differs by ancestry — i.e., resolution loss doesn't just lose precision uniformly, it *specifically erases ancestry-informative variation*, an argument for why higher resolution matters for equity in HLA-based clinical algorithms (transplant matching, pharmacogenomic risk alleles like HLA-B*57:01, HLA-B*15:02).

**Data needed**: `allele_2field`, `allele_4field`, `ancestry_pred` (or continuous), from Immuannot cohort only (has 4-field resolution).

**Implementation sketch**: for each 2-field group, a small stacked-bar or sunburst (as in 1.9) faceted by ancestry, showing 4-field subtype proportions; or a single "fan-out entropy" scalar (Shannon entropy of the 4-field distribution conditional on 2-field group) plotted by ancestry — turns the qualitative claim into one comparable number per group per ancestry.

**Precedent**: motivated by the WHO nomenclature hierarchy structure (https://www.researchgate.net/figure/HLA-nomenclature-with-allele-resolution-at-four-fields-eight-digits-Commonly-the-HLA_fig1_328519906) combined with the documented resolution-dependent allele/haplotype divergence seen across high-resolution NGS frequency studies (17th IHIW report — https://pmc.ncbi.nlm.nih.gov/articles/PMC8315142/); the specific "conditional entropy by ancestry" framing is **novel**, not found directly in the literature searched.

### 3.5 Non-classical gene diversity map (MICA/MICB/TAP1/TAP2/DRB3-4-5/DQA2/DQB2/HLA-E/F/G)

**Claim it supports**: non-classical and framework genes, essentially invisible to short-read ensemble callers restricted to the 30-gene classical panel typically reported, show their own ancestry-structured diversity — first population-scale look at some of these loci in AoU-scale data.

**Data needed**: allele/genotype calls for the non-classical loci, present only in the Immuannot cohort; `ancestry_pred`/`probabilities`.

**Implementation sketch**: same heatmap/diversity-index toolkit as 1.2/1.10, applied specifically to this gene set, explicitly labeled as "genes unavailable in AoU-native calls" to make the resolution/coverage advantage visually explicit (e.g., a coverage-availability matrix: genes × cohorts, binary, as a simple companion figure).

**Precedent**: general recognition that non-classical HLA genes (E/F/G, MIC, TAP) are under-characterized at population scale relative to classical loci is implicit throughout the classical-locus-focused literature surveyed (nearly every population frequency paper found is restricted to A/B/C/DRB1/DQB1/DPB1) — **no direct population-frequency-map precedent for these specific genes found in this search; treat as a novel contribution**.

### 3.6 Confidence/accuracy-stratified frequency estimate ("does filtering on template_distance change your answer?")

**Claim it supports**: whether population-level conclusions (e.g., allele X is common in ancestry Y) are robust to restricting to the high-confidence (`template_distance`-filtered) subset — a direct sensitivity/robustness figure that most papers cannot even attempt since they lack a comparable accuracy covariate.

**Implementation sketch**: paired dumbbell/slope plot per allele per ancestry: point A = frequency in full Immuannot cohort, point B = frequency in high-confidence subset, connected by a line; large shifts flagged. Equivalent to a Bland-Altman plot (mean vs difference) across the two frequency estimates.

**Precedent**: general Bland-Altman-style agreement analysis is standard method-comparison practice (not HLA-specific — **UNVERIFIED as an HLA-specific precedent**, but methodologically well established broadly).

---

## Part 4 — Design & technical guidance

### 4.1 Ancestry color palette

No standard, documented gnomAD or All of Us hex-code ancestry color palette could be found/verified in this search — gnomAD's genetic-ancestry announcement page describes group definitions and sizes but the fetched text did not surface a color specification, and gnomAD-browser source files searched did not resolve (**UNVERIFIED — the color scheme likely lives in the gnomAD-browser frontend JS/TS bundle, e.g. a `constants.ts` or theme file, which was not directly accessible via search/fetch in this pass; recommend inspecting `https://github.com/broadinstitute/gnomad-browser` locally with `git clone` + `grep -ri` for population color constants if pixel-exact consistency with gnomAD figures is required**).

Given no verified external standard, recommend a **colorblind-safe, print-safe qualitative palette** for the 6 AoU ancestry groups (AFR/AMR/EAS/EUR/MID/SAS) drawn from a vetted source rather than an invented one — e.g., Okabe-Ito or ColorBrewer "Set2"/"Dark2" (both explicitly designed for color-vision deficiency and are the de facto standard in genomics figures), fixed to a stable ancestry→color mapping used consistently across every figure in the paper:
- AFR → e.g. Okabe-Ito vermillion
- AMR → Okabe-Ito orange/yellow
- EAS → Okabe-Ito green
- EUR → Okabe-Ito blue
- MID → Okabe-Ito reddish-purple
- SAS → Okabe-Ito sky-blue
(Exact hex assignment left to implementation — the key methodological point, not independently verifiable as "the gnomAD standard," is: pick once, document it, and never let a color mean a different ancestry group in two figures.)

### 4.2 Long tails of rare alleles

- Always report N per allele/population cell; suppress or grey out cells below a minimum count threshold (e.g. N<5 raw observations) rather than plotting a noisy point estimate.
- Use log or symlog frequency axes for any figure spanning >2 orders of magnitude (virtually all HLA frequency figures do, given the documented handful-of-common vs thousands-of-rare-allele structure — https://academic.oup.com/g3journal/article/5/5/931/6025555).
- Bucket rare alleles into an explicit "other/rare (<X%)" category in bar/stacked/sunburst charts rather than silently truncating the legend.

### 4.3 Confidence intervals for frequency estimates

- For any allele frequency, especially rare alleles or small-N ancestry groups (AFR/AMR/EAS/MID/SAS will all be much smaller N than EUR in AoU), use **Wilson score** or **Jeffreys** intervals rather than the naive normal-approximation (Wald) interval — both are specifically recommended over Wald for small n / rare-event proportions, with Jeffreys noted as "the most dependable choice for rare event settings where preserving nominal confidence levels is crucial" — https://davidzhao1015.github.io/blog/2025/benchmark-interval-prop/ ; see also the broader binomial-CI-for-rare-events literature — https://arxiv.org/pdf/2109.02516 and https://www.tandfonline.com/doi/full/10.1080/00031305.2024.2350445.
- `statsmodels.stats.proportion.proportion_confint(count, nobs, method='wilson')` or `method='jeffreys'` implements both directly.
- Always show error bars/bands on any cross-ancestry frequency comparison bar chart — a bare point estimate for an N=40 ancestry group next to an N=300,000 group without a CI is actively misleading.

### 4.4 Log vs linear scale

- Frequency axes: log by default for anything with rare alleles (see 4.2); linear only for restricted-to-common-allele summary panels aimed at a lay/clinical audience.
- Sample-size axes (N per ancestry group): log, since AoU ancestry group sizes will span 2+ orders of magnitude.
- template_distance: likely a small-integer count (0, 1, 2, ...) — treat as discrete/ordinal, bar or step histogram, not continuous log scale; a project note already indicates a hard "distance==0" threshold was found unusable on real data, implying the full discrete distribution (not just an exact-match binary) should be shown.

### 4.5 Genes with wildly different allele counts

- Never share a single y-axis scale for allele-count/diversity figures across genes (e.g. HLA-B has vastly more known alleles than HLA-DQA1); use faceted small multiples with independent y-axes, or normalize to "% of gene's total observed alleles" rather than raw counts when comparing shape across genes.
- The 17th IHIW/US-population frequency papers and "HLA map of the world" both restrict cross-population comparison figures to per-gene panels for exactly this reason (https://pmc.ncbi.nlm.nih.gov/articles/PMC8315142/, https://pmc.ncbi.nlm.nih.gov/articles/PMC10076764/) — follow that convention rather than pooling genes into one undifferentiated frequency axis.

### 4.6 Continuous ancestry data provenance (AoU-specific)

All of Us derives both the discrete `ancestry_pred` (a PCA-projection-based classifier trained on a 1000 Genomes + HGDP reference, using the same categorical labels as gnomAD) and continuous ancestry fractions via the **Rye** program applied to the first 25 genotype PCs — per "Genetic ancestry and population structure in the All of Us Research Program cohort" (Nature Communications 2025 / bioRxiv preprint) — https://www.nature.com/articles/s41467-025-59351-8, https://www.biorxiv.org/content/10.1101/2024.12.21.629909v1, PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC12049439/. This paper is the authoritative methods reference for our `ancestry_pred`/`probabilities`/`pca_features` columns and should be the primary citation anywhere we describe how AoU ancestry was derived — it also reports 83–98% correspondence between discrete label and dominant continuous-fraction group, which is itself worth reproducing as a validation figure (discrete label vs continuous fraction agreement, e.g. a confusion-matrix-style heatmap of predicted label vs argmax(probabilities)).

---

## Part 5 — Ranked recommendations

### If we only build 8 figures, build these

1. **Allele × ancestry heatmap, double-clustered, per gene** (§1.2) — the canonical, expected "first figure" reviewers will look for; cheap to produce from any cohort.
2. **Same-allele frequency scatter: AoU-native vs Immuannot** (§3.1) — the single most novel, citation-anchored result (reference/short-read bias), directly leverages our two-cohort design.
3. **Diversity indices (heterozygosity/Shannon/rarefied richness) by ancestry, per gene** (§1.10) — standard, expected, comparable to published worldwide numbers.
4. **Allele frequency vs continuous ancestry proportion (smooth curve + decile CIs)** (§2.1) — the headline "continuous ancestry" differentiator figure.
5. **template_distance distribution by ancestry** (§3.2) — novel reference-bias evidence unique to our long-read + distance-to-reference field.
6. **Resolution cascade: 2-field → 4-field fan-out by ancestry** (§3.4) — makes the "resolution matters and matters unevenly" argument concretely and visually.
7. **Phased DRB1~DQA1~DQB1 haplotype frequency (alluvial/Sankey), true phase vs EM-imputed comparator** (§1.8 + §3.3) — showcases the phasing advantage directly.
8. **Fst/genetic-distance dendrogram between AoU ancestry groups, both cohorts overlaid** (§1.7) — ties everything back to expected population-genetics structure as a sanity/validation anchor.

### Figures that specifically exploit our unique advantages

- Same-allele frequency scatter, cohort vs cohort (§3.1) — exploits having two independent calling technologies on (overlapping) samples.
- template_distance by ancestry (§3.2) — exploits the novel accuracy/distance-to-reference field; nobody else has this.
- Phased cis-haplotype frequency, true vs EM-imputed (§1.8/§3.3) — exploits true phase.
- Resolution cascade 2-field → 4-field by ancestry (§3.4) — exploits 4-field resolution.
- Non-classical gene diversity map (§3.5) — exploits the extended 30+ gene panel including MICA/MICB/TAP1/TAP2/DRB3-4-5/DQA2/DQB2/HLA-E/F/G.
- Continuous ancestry-proportion regression / smooth AF curve (§2.1, §2.2) and ancestry-continuum density plots (§2.3) — exploit the continuous 6-way `probabilities` vector instead of collapsing to discrete labels.
- High-confidence-subset sensitivity dumbbell/Bland-Altman plot (§3.6) — exploits `template_distance`-based filtering as an internal accuracy stratifier.

---

## Note on verification

Every URL above was retrieved via live web search/fetch in this session. Claims explicitly marked **UNVERIFIED** are either (a) synthesized from search-result summaries without being traced to a specific primary figure/page I could independently confirm, or (b) design recommendations not tied to any single citation. Treat UNVERIFIED items as directionally reliable but re-check the primary source before citing them in the actual manuscript.
