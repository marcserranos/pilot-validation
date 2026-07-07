# Project Omni-HLA Genomics: High-Resolution Allele Calling & Non-Linear PRS on the All of Us Biobank

## Executive Summary
This repository contains the pipeline and analytical frameworks for executing full-resolution HLA typing and advanced Polygenic Risk Score (PRS) modeling on the *All of Us* (AoU) Biobank[cite: 1, 2]. By moving away from traditional, inaccurate HLA imputation and instead utilizing direct, short-read and long-read allele calling via SpecHLA, this project aims to deliver state-of-the-art disease risk stratification across diverse ancestral populations[cite: 1].

---

## 1. Project Premise & Dataset
* **The Cohort:** The *All of Us* Biobank features a highly diverse multi-ancestry cohort of over 400,000 individuals sequenced using short-read Whole-Genome Sequencing (srWGS)[cite: 1]. 
* **Long-Read Sub-Cohort:** Approximately 10,000 individuals have been sequenced using highly accurate long-read WGS[cite: 1]. This subset will serve as a validation layer to ensure our HLA calling and downstream PRS models remain accurate and consistent across sequencing technologies[cite: 1].
* **Biological Target:** The Human Leukocyte Antigen (HLA) region on Chromosome 6[cite: 1]. This hyper-polymorphic locus contains critical genes responsible for antigen presentation, driving massive associations with autoimmune diseases, oncology, and infectious pathogen response[cite: 1].

---

## 2. Core Scientific Aims

### Aim 1: Frequency of HLA Alleles Across Diverse Ancestries
* **Methodology:** Leverage pre-existing local ancestry calls (which segment individual DNA tracks by population origin) to map HLA allele distribution[cite: 1].
* **Hypothesis:** Because the HLA region is subject to intense evolutionary selection pressure, we anticipate observing major, distinct differences in allele frequencies across Native American, European, East Asian, and African ancestries[cite: 1].

### Aim 2: Non-Linear Polygenic Risk Scores (PRS) Using HLA
* **The Problem with Standard PRS:** Traditional PRS architectures assume purely additive, linear effects, calculating a basic weighted sum of variants derived from univariate, single-variant GWAS models[cite: 1]. This completely ignores the dense, highly correlated haplotype blocks unique to the HLA region[cite: 1].
* **Our Approach:** Model non-linear interactions and epistasis within the HLA region[cite: 1]. The pipeline will execute:
  1. Fine-mapping to  isolate independent, non-linkage-driven variant signals[cite: 1].
  2. Training non-linear models (such as Elastic-Net, LASSO, custom tabular architectures, and Gradient-Boosted Trees like CatBoost) to learn complex variant-variant interactions[cite: 1].
* **Target Phenotypes:** Autoimmune diseases with large HLA effect profiles, including Celiac Disease, Type 1 Diabetes (T1D), Psoriasis, and Multiple Sclerosis (MS)[cite: 1].

### Aim 3: Uncovering the Genetic Architecture & Tails of Autoimmune Disease
* **Focus:** Evaluate how selection pressures and gene-by-gene epistasis dictate risk at the extreme tail ends of the PRS distribution[cite: 1].
* **Clinical Significance:** Recent literature highlights that interaction-aware HLA modeling drastically optimizes clinical stratification[cite: 1]. For instance, modeling epistasis can boost the positive predictive value (PPV) in high-risk Celiac settings from 17.5% to 27.1% while maintaining a negative predictive value (NPV) above 99%[cite: 1]. Our goal is to shift focus from generic AUC metric optimization to real-world tail enrichment and tail calibration[cite: 1].

---

## 3. The Competitive Edge: "Eating Their Lunch"
Existing corporate pipelines evaluating the *All of Us* cohort have critical flaws that this project directly resolves[cite: 1]:
> **The Competitor's Limitation:** Previous large-scale approaches (e.g., HLA-ARC frameworks) did not perform direct HLA allele calling[cite: 1]. Instead, they imputed classical alleles via `SNP2HLA` using localized SNP arrays[cite: 1]. 
> 
> While imputation achieves 95-97% accuracy for common European haplotypes, it heavily degrades for rarer alleles and across diverse, underrepresented ancestries[cite: 1]. This introduces non-differential misclassification, attenuating true HLA effect sizes[cite: 1]. 
> 
> **Our Solution:** Direct short-read allele calling via `SpecHLA` significantly outperforms imputation, maximizes discovery power for autoimmune traits, and delivers accurate predictions across non-European ancestries[cite: 1, 2].

---

## 4. Software Stack & Ecosystem
* **Primary Phase Caller:** [SpecHLA](https://github.com/deepomicslab/SpecHLA) (leveraging read binning and local diploid assembly for short and long reads)[cite: 1, 2].
* **Future Long-Read Frameworks:** [SpecImmune](https://github.com/deepomicslab/SpecImmune) or collaborator-developed [HLA-Resolve](https://www.medrxiv.org/content/10.64898/2026.03.27.26349549v3)[cite: 1].
* **Downstream Modeling:** CatBoost / Gradient-Boosted Decision Trees for non-linear epistatic interaction mapping[cite: 1].

---

## 5. Reference Compendium & Literature Review
For the AI agent and collaborators, prioritize reading the following literature resources to understand baseline configurations:

### HLA Biology & Nomenclature
* [Nature Reviews Immunology: HLA Review](https://www.nature.com/articles/nri.2017.143#Sec1)[cite: 1]
* [PubMed: HLA Nomenclature Standards](https://pubmed.ncbi.nlm.nih.gov/32307125/)[cite: 1]

### Polygenic Risk Scores & Epistasis
* [Nature Reviews Genetics: PRS Review (2025)](https://www.nature.com/articles/s41576-025-00900-8)[cite: 1]
* [Nature: Rare Variants and PRS Distribution Tails (2026)](https://www.nature.com/articles/s41586-026-10516-5)[cite: 1]
* [Nature Genetics: CatBoost Application for Type 1 Diabetes (2026)](https://www.nature.com/articles/s41588-026-02578-y)[cite: 1]

### Linkage & HLA Interaction Methods
* [BioRxiv: Identifying Independent Signals in Correlated Haplotypes](https://www.biorxiv.org/content/10.1101/2020.05.28.119669v2.full.pdf)[cite: 1]
* [PubMed Central: Historical Epistasis Models](https://pmc.ncbi.nlm.nih.gov/articles/PMC4552599/)[cite: 1]
* [MedRxiv: Imputation vs. Calling Benchmarks (HLA-ARC Insights)](https://www.medrxiv.org/content/10.1101/2025.10.02.25337098v3.full.pdf)[cite: 1]

---

## 6. Initial Sandbox & Prototyping Tasks (The Technical Spike)
Before scaling pipelines to the full AoU environment, the following engineering tasks must be completed[cite: 1]:
1. **Environment Setup:** Build the environment from source using the provided `environment.yml` and ensure your compiler meets the `GCC 9.4.0+` prerequisite[cite: 2].
2. **Database Indexing:** Run `bash index.sh` to construct local reference indexes[cite: 2].
3. **Execution Verification:** Navigate to the `example/` folder and execute `bash test_all.sh`[cite: 2].
4. **Output Validation:** Confirm successful generation of `hla.result.txt` and check that the resulting allele strings match known control samples[cite: 2].