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

1. **Current release / total allele count.** Per "Database Growth" / Statistics pages at EBI (https://www.ebi.ac.uk/ipd/imgt/hla/about/statistics/growth ; https://hla.alleles.org/pages/genes/statistics/) as summarized in search snippets: **release 3.63 (2026-01)** reports **44,876 distinct allelic variants** (79,871 component/sub-sequence entries). UNVERIFIED at the primary-page level this session (search-snippet only; recommend a direct check of https://www.ebi.ac.uk/ipd/imgt/hla/about/statistics/growth before using in the paper, since exact release numbers move monthly).

2. **Full-length genomic vs exon-only coverage.** "IPD-IMGT/HLA database: recent developments in sequence submission", Nucleic Acids Research, 54(D1):D1152, 2026 (Oxford Academic: https://academic.oup.com/nar/article/54/D1/D1152/8326453). Relevant: as of ~July 2025 the database has "over 43,000 unique sequences" spanning exons/introns/UTRs (3-17 kb); since the prior (2023) NAR update, **8,148 alleles were added (+23%)**, and the **proportion of alleles with full genomic sequence rose from 51.2% to 57.6%** — i.e., even in the most recent release, ~42% of named alleles still lack full-length genomic sequence (many labs submit only exons 2/3, the peptide-binding-groove exons). This is a strong, directly citable number for our "reference completeness" framing.
   - Coverage is worse for HLA class II, "particularly HLA-DRB1," reflecting sequencing difficulty of those genes (same NAR paper family / EBI genomics help page: https://www.ebi.ac.uk/ipd/imgt/hla/help/genomics.html).
   - UNVERIFIED / GAP: no source found this session that explicitly quantifies **ancestry-stratified bias** in full-genomic-sequence coverage (i.e., "European-population alleles are more likely to have full-length sequence than African/AMR/EAS alleles"). This appears to be an actual gap in the literature — the closest statements are about exon-only submission being a general historical practice, not an ancestry-specific quantification. This is a genuine opportunity for our paper to make a novel, citable claim if we can show it empirically from IPD-IMGT/HLA allele-frequency-by-population data.

3. **Submission requirements for novel alleles.**
   - "IPD-IMGT/HLA database: recent developments in sequence submission" (NAR 2026, D1152, as above) and prior editions ("IPD-IMGT/HLA Database", NAR 51(D1):D1053, 2023, https://academic.oup.com/nar/article/51/D1/D1053/6814448; PMC9825470) describe manual + automated curation with "strict acceptance criteria."
   - Confirmation practice: independent confirmation is achieved by repeating PCR and sequencing ("dual confirmation"); it is "highly desirable that novel alleles are characterized and submitted in full length," and known partial alleles are encouraged to be extended to full gene length. Source: search snippet from NAR pages plus "Dual redundant sequencing strategy: Full-length gene characterisation of 1056 novel and confirmatory HLA alleles," PMC6084308 (journal/year UNVERIFIED beyond PMCID — appears to be an older paper, likely Human Immunology or Tissue Antigens; needs direct fetch to confirm).
   - Submission tooling: **TypeLoader / TypeLoader2** — "Automated submission of novel HLA and killer-cell immunoglobulin-like receptor alleles in full length," PMC6594033 (journal/year UNVERIFIED, appears ~2019 Human Immunology based on PMCID range).
   - Also relevant to validation standards (cross-ref section E): "Resolving unknown nucleotides in the IPD-IMGT/HLA database by extended and full-length sequencing of HLA class I and II alleles," Immunogenetics (Springer), 2024, https://link.springer.com/article/10.1007/s00251-024-01333-z (PMC10944811) — describes ongoing effort to close partial/unknown-nucleotide gaps in existing database entries via full-length resequencing; relevant precedent for "many database alleles are themselves incomplete," reinforcing our exon-vs-genomic novelty distinction.
   - "Novel alleles in the era of next-generation sequencing-based HLA typing calls for standardization and policy" (PMC10611506) is directly on-topic for validation/submission policy discussion — flagged for a follow-up fetch (not done this session; UNVERIFIED beyond PMCID).

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
