# WS5 Literature Research — HLA Novelty / Coverage / Positioning

Status: IN PROGRESS. Written incrementally; sections filled as searches complete.
Scope: literature support for Omni-HLA paper (~12,000 AoU long-read HLA calls, Immuannot/IPD-IMGT/HLA, novel alleles mostly non-coding, ancestry gradient, unsaturated discovery curves).

---

## A. Prior large-scale long-read/assembly HLA typing & novel allele discovery (2022-2026)

1. **"A pan-MHC reference graph with 246 fully contiguous phased sequences"**, bioRxiv preprint, 2023 (posted 2023-09-01). URL: https://www.biorxiv.org/content/10.1101/2023.09.01.555813v1.full (PDF: .../555813.full.pdf). Author list UNVERIFIED (not confirmed from a fetched page — WebFetch was rate-limited/blocked; only search snippets available).
   - Relevant: assembled 246 fully contiguous, phased MHC haplotypes (mostly from HPRC long-read data, plus 4 cell lines with a novel targeted enrichment method); built a "pan-MHC" reference graph across 39 loci (class I + II HLA genes). Reported **1,246 alleles absent from IMGT/HLA** at high resolution — directly comparable novelty-rate precedent for our ~12,000-sample AoU work, though their N is far smaller (order ~120-250 haplotypes/individuals vs our ~12,000 diploid assemblies, i.e. ~24,000 haplotypes).
   - Also reported C4A/C4B copy-number variation and LD between C4 haplotypes and 14 MHC loci (less relevant to our novelty/coverage angle).
   - UNVERIFIED: whether they give an ancestry breakdown of the 1,246 novel alleles or any allele-space-remaining estimate — could not confirm from snippets alone; needs a direct fetch of the full text (not obtained this session due to 429).

2. **"Population-scale Long-read Sequencing in the All of Us Research Program"**, medRxiv preprint, posted 2025-10-05 (DOI 10.1101/2025.10.02.25336942), also PMC12622093, PubMed 41256123. URL: https://www.medrxiv.org/content/10.1101/2025.10.02.25336942v1
   - Relevant: first large-scale AoU long-read sequencing analysis paper — establishes the precedent/venue for AoU LRS population genomics papers. Cohort: **1,027 individuals self-identifying as Black or African American**, ~8x coverage PacBio HiFi, cloud-native pipelines; focus is structural variants (SVs) and 226 trait associations (191 SV-disease pairs), with ancestry-specific (African-similar ancestry) associations showing larger effect sizes / lower allele frequencies.
   - This is a direct AoU LRS precedent paper but is SV/GWAS-focused, not HLA-focused — good for "AoU LRS is a validated framework" framing, not a competing HLA novelty paper. Sample size (1,027, single ancestry group, ~8x coverage) is much smaller/shallower than our ~12,000 multi-ancestry HiFi assemblies — worth stating explicitly as a scale contrast.
   - UNVERIFIED: whether this paper does any HLA-specific analysis at all (search snippets only mention SVs/CYP2D6/repeat expansions elsewhere as separate tool papers, not confirmed inside this specific paper).

3. **FuFiHLA: a tool for full-field HLA typing from long-read data**, Bioinformatics, 2026 (Oxford Academic, vol 42 issue 5, btag231); also bioRxiv 2025.10.23.684216 and PMC13221240 / PMC12633391. URL: https://academic.oup.com/bioinformatics/article/42/5/btag231/8669788
   - Relevant: motivates our work's methods choice — notes existing HLA typing tools (mostly commercial/proprietary) "lose efficacy" on lower-coverage WGS reads typical of large population cohorts like AoU and the Singapore pangenome project. Supports 6 genes (HLA-A/B/C/DRB1/DQA1/DQB1) from long reads.
   - Useful as evidence that low-coverage population-scale HLA typing is an open/active tooling problem — but does NOT report novel-allele discovery rates or large-N novelty statistics itself (a tool paper, not a discovery paper). Ancestry-stratified coverage of IPD-IMGT/HLA: UNVERIFIED, not mentioned in snippet.

4. **HLA-Resolve: High-Resolution HLA Haplotyping Using Long-Read Hybrid Capture**, medRxiv, 2026 (posted ~2026-03-27; DOI 10.64898/2026.03.27.26349549v2). URL: https://www.medrxiv.org/content/10.64898/2026.03.27.26349549v2.full
   - Relevant: another 2026 long-read HLA haplotyping method paper (hybrid capture, not WGS/assembly) — a candidate "already done?" comparator; details of sample size/novelty NOT fetched this session (UNVERIFIED).

5. **Zhou Y, Song L, Li H. "Full-resolution HLA and KIR gene annotations for human genome assemblies."** Genome Research, 2024, 34(11):1931-. Also bioRxiv 2024.01.20.576452; PMC10849470; PubMed 38839374. URL: https://genome.cshlp.org/content/34/11/1931
   - This is the Immuannot paper itself (our project's annotation tool) — recording its own headline novelty numbers is directly relevant as the baseline our ~12,000-sample AoU result must be compared against and framed as a scale-up of. Applied to **56 regional + 212 whole-genome assemblies** (268 assemblies total, i.e. ~2 orders of magnitude smaller than our ~12,000 diploid samples), annotated **9,931 HLA/KIR genes**, found **4,068 (~41%) with novel sequence vs IPD**, represented by **2,664 distinct novel alleles**, of which only **92 had novel protein-level (nonsynonymous) changes** — i.e., the vast majority of novel alleles were non-coding/synonymous differences. This is essentially the same qualitative pattern we report (most novelty outside coding sequence) but at ~45x smaller assembly count — strong precedent to cite AND a key numeric contrast (92/2,664 ≈ 3.5% protein-changing novel alleles in Immuannot's own paper vs whatever fraction we find at 12,000-sample scale).
   - The paper explicitly frames itself as expected to "speed up discovery of new HLA/KIR alleles" — supports our framing that a large-scale application (ours) is a natural, anticipated follow-up rather than something already done.

6. **"Novel alleles in the era of next-generation sequencing-based HLA typing calls for standardization and policy"**, PMC10611506 (journal/year not yet confirmed from snippet — UNVERIFIED bibliographic details beyond PMCID).
   - Relevant to section B/E (submission/validation standards for novel alleles); flagged for direct fetch.

7. **Lion T. et al.(?) "Long-Read Next-Generation Sequencing Technologies Can Address Some Limitations of Short-Read Technologies in HLA Typing."** International Journal of Immunogenetics, 2026, DOI 10.1111/iji.70046. URL: https://onlinelibrary.wiley.com/doi/10.1111/iji.70046
   - Relevant: 2026 review-type paper on long-read HLA typing generally; candidate reviewer comparator for section F — needs a direct fetch to check if it makes ancestry-stratified novelty or saturation claims (UNVERIFIED, not yet fetched).

8. UK Biobank HLA work found so far is short-read/exome based (HLA-HD calling on WES, 129 novel autoimmune-disease associations, medRxiv 2023.01.15.23284570; and a British-Africans admixture/HLA imputation paper in Eur J Hum Genet 2025, https://www.nature.com/articles/s41431-025-01888-9) — NOT long-read allele-discovery work, so lower priority for section A but useful as a "no UKB long-read HLA discovery competitor found" negative result.

Still to check for section A: HPRC release 2 formal MHC follow-up (beyond 2023 pan-MHC preprint), explicit published estimate of total undiscovered HLA allele space / ancestry-stratified IPD-IMGT/HLA coverage (not yet found — may not exist, relevant to sections C/F).

## B. IPD-IMGT/HLA database: release, allele counts, full-length vs exon-only, ancestry bias, submission requirements

(pending)

## C. Coverage-based rarefaction/extrapolation (Chao & Jost 2012, iNEXT, Good-Turing) — formulas and genomics uses

(pending)

## D. Reference panel design / "how many genomes from which ancestry" arguments

(pending)

## E. Validation standards for novel HLA alleles from long reads; MHC assembly failure modes

(pending)

## F. Papers that would make a reviewer say "already done"

(pending)

---

## Distilled (≤10 lines)

(pending)
