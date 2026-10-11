# READING.md — TRBV usage / MHC genetics / long-read HLA-KIR literature review

**Compiled:** 2026-10-11 (literature search via web; no AoU data touched)
**Scope:** (A) Sharon et al. 2016 in depth; (B) TRBV germline deletion/CNV and V-usage genetics;
(C) scoop check against the AoU long-read flagship and 2024–2026 long-read HLA/KIR work;
(D) implications for the four planned paper claims.

**Verification convention.** Everything below was read from the source indicated. Statements I
could not confirm from a primary source are marked **[unverified]**. Where a publisher page was
paywalled or robot-blocked I say so explicitly rather than filling the gap from memory.

Sources that could **not** be read in this session: PubMed Central HTML (repeated reCAPTCHA
interstitials on `pmc.ncbi.nlm.nih.gov`), eLife article HTML (JS challenge),
`link.springer.com/article/10.1007/s00251-004-0684-z` and `jimmunol.org` (rate-limited / JS-only),
Europe PMC REST API (blocked by egress policy). Where a paywalled version blocked me I used the
author-hosted PDF or the bioRxiv/medRxiv preprint and I flag which version the numbers come from.

---

## A. Sharon et al. 2016 — the closest prior art to the V-usage × MHC question

### A.1 Citation (verified)

> Sharon E, Sibener LV, Battle A, Fraser HB, Garcia KC, Pritchard JK.
> **Genetic variation in MHC proteins is associated with T cell receptor expression biases.**
> *Nature Genetics* **48**(9): 995–1002 (2016). doi:10.1038/ng.3625. PMCID: PMC5010864.

Verified against two independent sources: the Nature article page (volume 48, pages 995–1002,
published online 1 Aug 2016; received 16 Jan 2016, accepted 23 Jun 2016, doi 10.1038/ng.3625) and
the Pritchard lab publication list, which renders it as "Nat Genet, 2016, 48(9):995-1002".
Author order and spelling as above are from the Nature page / lab list. PMID **[unverified]**
(PubMed was not reachable; PMCID PMC5010864 is confirmed by search-result metadata and by the
PMC URL resolving to this article).

Methods/results numbers in A.2–A.6 are taken from the **author-hosted accepted PDF**
(`https://web.stanford.edu/group/pritchardlab/publications/pdfs/Sharon16.pdf`), i.e. the full text,
not the abstract.

### A.2 Cohort and data

- **n = 922** individuals of **European ancestry** from the **Depression Genes and Networks (DGN)**
  project (the Battle et al. 2014 whole-blood RNA-seq resource); **895** retained after QC filters.
- **Whole-blood bulk RNA-seq**, ~70 M **51-bp** reads per individual. No sorted T cells, no
  targeted TCR-seq (this matters a lot — see A.7).
- Final phenotype set: **44 Vα, 40 Vβ, 11 Vγ, 3 Vδ** TCR V genes, plus Ig V genes as negative
  controls (149 Ig V genes analysed for the GCTA control).

### A.3 How V usage was quantified

- Bowtie2 + TopHat alignment (default parameters); a **modified HTSeq** counting uniquely-mapped
  reads to each V / J / C gene, explicitly **allowing a read to map to more than one V gene**
  (i.e. multi-mapping V genes are not resolved to a single gene).
- Sample filter: individuals with <80,000 TCR/Ig reads dropped. Gene filters: mean >1 read per
  individual, ≤100 zero-count samples, ≥20 reads in ≥1% of individuals.
- Counts log-transformed (+0.1 pseudocount), regressed on technical/biological confounders **and on
  total chain read count**, then quantile-normalised to a normal distribution.
- Note what this measures: **V-gene expression in bulk blood**, not clonotype-level V usage from a
  repertoire assay. It conflates per-cell transcription, clonal expansion, and T-cell subset
  composition.
- Six Vα, one Vβ and three Vγ genes were **excluded from cis analysis** because of
  **haplotype-dependent mappability** differences — an acknowledged short-read artifact at the
  TR loci.

### A.4 HLA typing / imputation

- **SNP2HLA** imputation of 2- and 4-digit alleles for **HLA-A, -B, -C, -DRB1, -DQA1, -DQB1, -DPA1,
  -DPB1**, plus imputed **amino-acid** and SNP genotypes.
- Reference panel: **T1DGC, 5,225 European-ancestry individuals**.
- Imputation quality (Beagle r²): >0.95 for 97% and >0.99 for 80% of amino acids (excluding
  MAF < 2.5%).
- Consequences: European-only panel; allele resolution capped at 4-digit; **no novel alleles**, and
  rare alleles and non-European haplotypes are structurally invisible.

### A.5 Association model

- **cis**: Pearson correlation, variants within 1 Mb; empirical FDR from 10,000 permutations (5%).
- **trans**: Pearson correlation between V-gene phenotypes ≥1 Mb from the MHC and
  extended-MHC variants (hg19 chr6:25,912,984–33,290,793); genome-wide threshold P < 5 × 10⁻⁸;
  5% empirical FDR for MHC-restricted tests.
- **Conditional analysis**: forward stepwise linear regression with F-tests, adding positions while
  P < 0.05/n (≈1 × 10⁻⁵).
- **Fine-mapping**: Bayesian spike-and-slab model with MCMC, prior 0.1% per residue.
- **Heritability**: GCTA on MHC amino-acid variation.
- The accepted PDF does not spell out the PC/covariate set inside the conditional model
  **[unverified]**.

### A.6 Main findings

Scale of the MHC *trans* effect:

| Quantity | Value |
|---|---|
| Vα genes with genome-wide (P<5e-8) MHC association | **47.7%** |
| Vβ genes with genome-wide MHC association | **22.5%** |
| Vα genes associated in extended MHC at 5% FDR | **66%** |
| Vβ genes associated in extended MHC at 5% FDR | **35%** |
| Independent associations (conditional) | **66**, covering **28 Vα + 15 Vβ** genes |
| Variance explained per independent association | **2.2% – 37%** |
| Of the 66: in HLA genes / class II / HLA-DRB1 alone | **42 / 37 / 25** |
| Haplotype-level associations | **92**, of which 75 class II, 32 HLA-DRB1 |
| GCTA: Vα / Vβ genes with significant MHC-AA heritability | **33/44** and **16/40** |
| Share of that variance from MHC *protein* components | **92% (Vα), 88% (Vβ)** |
| Per-gene variance explained range (GCTA) | **5% – 75%** |
| Ig V-gene negative control significant | **1 / 149** |

Allele- and residue-level results:

- **HLA-DRB1\*03:01** associated with **12** V genes.
- Top Vα-influencing residues include **HLA-DRβ1 position 71** (third-ranked), **HLA-DRβ1 86**,
  **HLA-DRβ1 77** (posterior 0.71), and **HLA-DQβ1 57**; **HLA-DQβ1 185** is a spatially distal
  outlier. DRβ1 71/86 and DQβ1 57 are the classic RA / T1D / MS risk positions.
- **Vβ signals were markedly fewer and weaker than Vα** (detail relegated to their Supplementary
  Table 10) — important for us, because our partner's phenotype is **TRBV** usage.
- Structure: analysis of IMGT-curated TCR–pMHC structures (PDB as of 15 Jun 2015; 16 complexes in
  Fig. 5e). Residue 77 contacts germline-encoded TCR residues across class II complexes.
  Association probability vs distance to nearest TCR residue: Pearson **r = −0.28, P = 0.022**;
  Spearman **r = −0.35, P = 0.0044** — i.e. a real but modest structural correlation.

### A.7 cis / germline TRB findings — and what they did *not* do

- **cis-eQTLs for 78%** of TCR V genes and 46% of Ig V genes at 5% FDR. So local germline effects on
  V usage are large and pervasive in their own data.
- **The only genome-wide significant non-MHC trans hit for any TCR V gene was TRBV24-1**, near
  **ZNF443** (chr19). This is the same signal Russell et al. 2022 later recover (see B.3) — worth
  noting as an independently replicated locus.
- **No deletion or copy-number analysis at TRA/TRB.** Confirmed absent from the full text: germline
  TRBV CNV is neither genotyped nor used as a covariate. Their only nod in that direction is
  excluding genes with haplotype-dependent mappability.
- **No independent replication cohort.** Robustness came from internal checks: haplotype-based
  re-analysis, joint multi-phenotype regression, the Ig negative control, and conditioning on
  CD4:CD8-associated SNPs.

### A.8 Limitations the authors state themselves

1. Whole blood averages over CD4⁺ and CD8⁺ subsets; **CD4:CD8 ratio can confound**, especially at
   HLA-B.
2. The class II excess may partly reflect **greater power from CD4⁺ abundance**, not biology.
3. **LD prevents causal residue assignment**; they explicitly write that fine-mapping resolution
   should improve with **longer-read sequencing of sorted cells**.
4. Poorly imputed causal sites may have been missed.
5. Measurement noise from clonal expansion and mapping differences reduces power.
6. Structural correlation is modest; distal residues may act indirectly.

### A.9 What our project can add beyond Sharon 2016

| Axis | Sharon 2016 | Us |
|---|---|---|
| HLA determination | **Imputed** (SNP2HLA), 8 genes, 4-digit, European panel | **Called from phased long-read assemblies**, 65 genes via Immuannot, full-resolution incl. non-classical / pseudogenes / C4 |
| Novel alleles | Impossible by construction | Directly discoverable; discovery curves by ancestry (claim 2) |
| Phasing | **Statistical**, panel-based | **Physical**, per haplotype — enables DQA1~DQB1 haplotype-level analysis (claim 1) and residue–residue cis configuration |
| Amino-acid resolution | Imputed AA dosages, LD-limited | Exact AA sequence per haplotype, incl. positions absent from the imputation panel |
| Ancestry | European only (n=922) | ~12k AoU participants, deliberately diverse (the partner's SAS 43.5% vs EUR 16.7% TRBV4-3 contrast is exactly the regime they could not enter) |
| Germline TRB | Not modelled; TRBV24-1 trans hit unexplained | **TRB germline deletions genotyped from long-read depth** (claim 4) → usable as a covariate / as the exposure, separating "cannot use this V gene" from "chooses not to" |
| TCR phenotype | Bulk whole-blood V-gene read counts | Partner's RNA-seq TRBV usage; per-person zero-usage calls (TRBV4-3 unused in ~27%) |

The single sharpest methodological addition: **Sharon's phenotype cannot distinguish a deleted gene
from an unexpressed one.** A person homozygous for the ~21 kb TRBV4-3/TRBV3-2/TRBV6-2 deletion has
structurally zero usage, and in a European cohort with a quantile-normalised bulk-expression
phenotype that is just a low outlier. Our long-read CN genotype makes that explicit, and the
ancestry gradient in deletion frequency (B.5) is a direct candidate explanation for the
SAS-vs-EUR difference in never-used TRBV4-3 — which also means **we must test it before claiming an
immunological interpretation**.

---

## B. TRBV germline deletions / CNV and the genetics of V usage

### B.1 Omer et al. 2022 — TRB germline variability inferred from repertoires

> Omer A, et al. **T cell receptor beta germline variability is revealed by inference from
> repertoire data.** *Genome Medicine* **14**:2 (2022). doi:10.1186/s13073-021-01008-4.
> (Journal, volume and article number from search metadata and the journal URL; the author list
> beyond the first author is **[unverified]** — both the BMC page and PMC were unreachable this
> session. Numbers below are from the bioRxiv version, doi:10.1101/2021.05.17.444409.)

Four TCRβ repertoire datasets (DS1 n=28 and DS2 n=25 full-length; DS3 n=313 BIOMED-2; DS4 n≈786
Adaptive short reads — the preprint is internally inconsistent, giving both 768 and 786) were
processed with IgBLAST/IgDiscover, TIgGER genotyping, a binomial test for deletion calls and RAbHIT
for haplotype inference, validated against **35 diploid long-read assemblies**. They report the
**common ~21 kb deletion spanning TRBV4-3, TRBV3-2 and TRBV6-2**: in DS1, TRBV4-3 and TRBV3-2 were
absent in 8/28 individuals and were always inferred deleted together; in DS2, TRBV4-3 absent in 8/25.
Lower-frequency single-gene deletions were inferred for TRBV11-1, TRBV7-3, TRBV4-1, TRBV4-2, TRBV28
and TRBV30; no TRBD/TRBJ deletions. Haplotype anchoring on TRBJ1-6 showed **TRBV7-2\*02 marks the
TRBV4-3-deleted chromosome** (7/7 heterozygotes in DS1; 11/12 high-confidence DS3 individuals), and
linked the TRBJ1-6\*02 chromosome to TRBV24-1\*02, a TRBV29 deletion and TRBV20-1\*05. They catalogue
**38 undocumented TRBV alleles/SNPs** (15 corroborated in the long-read assemblies) and 31 novel
5′UTR variants. Crucially, their deletion calls are **absence from the expressed repertoire, not
proven absence from the genome** — the exact inferential gap our depth-based genotyping closes. No
population breakdown of deletion frequency is reported.

### B.2 Rodriguez / Watson et al. — targeted long-read phased assembly of TRA/D and TRB

> Preprint: **Targeted long-read sequencing facilitates phased diploid assembly and genotyping of
> the human T cell receptor alpha, delta and beta loci.** bioRxiv doi:10.1101/2022.05.24.493244.
> Published version indexed at Mount Sinai as a journal article (likely *Cell Genomics* 2023 —
> journal, volume and final DOI **[unverified]**, the publisher record was not reachable). Author
> list **[unverified]** (the bioRxiv PDF's author block did not render); commonly cited as
> Rodriguez OL, Watson CT and colleagues.

Long-fragment (4–6 kb) capture libraries across TRA/D and TRB in **11 individuals of the 1000
Genomes Project** (two trios plus five unrelated; Yoruba, Kinh Vietnamese, Japanese, Colombian,
Puerto Rican, European), PacBio Sequel IIe HiFi, producing phased diploid assemblies. In TRB they
find **three SVs >20 kb, two of them polymorphic among the nine unrelated samples**, one of which
**deletes TRBV6-2, TRBV3-2 and TRBV4-3** — i.e. independent assembly-level confirmation of the
Omer-inferred 21 kb event. Two trio probands (HG02059, NA18507) carry **20.2 kb and 20.4 kb TRB
insertions**. Smaller events: largest intergenic deletion 1,044 bp, largest insertion 665 bp (a 37 bp
tandem-repeat expansion), a 331 bp Alu insertion in 6 individuals, an 800 bp deletion in 5, a 592 bp
insertion in 6. They report **79 putative novel alleles (42 TRA, 35 TRB, 2 TRD)**, 13 supported by
>10 HiFi reads in ≥2 individuals and 66 singletons, a potential **+23% to IMGT**; 7 of the 35 TRB
alleles matched an AIRR-seq inference study; **samples of African ancestry carried the most novel
alleles**. No TRBV24-1 SV. No per-population SV frequencies (Supplementary Table 4 holds the
genotypes). This is the methodological template closest to our claim 4 — but **n=11, targeted
capture, no population genetics**.

### B.3 Russell et al. 2022 — GWAS of V(D)J recombination parameters and V usage

> Preprint read: **Combining genotypes and T cell receptor distributions to infer genetic loci
> determining V(D)J recombination probabilities.** bioRxiv doi:10.1101/2021.09.17.460747.
> Published in *eLife* **11**:e73475 (2022), doi:10.7554/eLife.73475 (eLife article number verified
> from the eLife URL `elifesciences.org/articles/73475`; the eLife HTML itself was JS-blocked, so
> final page/author list from the published version is **[unverified]**; commonly cited as
> Russell ML, Bradley P, Matsen FA IV et al.). Note: this is *eLife*, **not** *Sci Immunol* — the
> "Russell 2022 Sci Immunol" formulation in the task prompt appears to be a conflation and should
> not be cited as such.

Discovery: **398** individuals with paired TCRβ immunosequencing (Adaptive; from the 666 healthy
bone-marrow-donor cohort) and SNP arrays (dbGaP phs001918), clustered by PCA into
Caucasian-associated 322, Asian 23, Hispanic 30, African 8, Native American 10, Middle Eastern 5 —
i.e. **effectively a European-ancestry study**. Validation: 94 individuals from Managua, Nicaragua
with paired TCRα/β. Findings: **DNTT (TdT)** variants (rs2273892, rs12569756) associated with
N-insertion counts; variation near **DCLRE1C (Artemis)** with V- and J-gene trimming; in cis at TRB,
the strongest signal is **TRBV28 expression**, top SNP **rs17213** in its recombination signal
sequence; **MHC variation (top hits in HLA-DRB1) associated with productive TRBV10-3 usage**, with
**12.1% of productive V genes showing MHC associations** and **no MHC association for non-productive
usage** (an elegant internal control for thymic selection vs recombination); and a chr19 region near
**ZNF443 / ZNF709** associated with **TRBV24-1 and TRBV24/OR9-2** usage, mostly in non-productive
TCRs — replicating Sharon's lone non-MHC trans hit. They **do not mention TRBV4-3, TRBV3-2 or any
TRBV deletion/SV**; their TRB-locus controls are a TRBD2 genotype correction term and restriction to
TRBJ1 TCRs. The paper also does not discuss Sharon et al. by name in the preprint text
**[unverified]** whether the published version does.

### B.4 Ishigaki et al. 2022 — HLA restricts the CDR3, not V usage

> Ishigaki K, Lagattuta KA, Luo Y, James EA, Buckner JH, Raychaudhuri S. **HLA autoimmune risk
> alleles restrict the hypervariable region of T cell receptors.** *Nature Genetics* **54**(4):
> 393–402 (2022). PMID 35332318. DOI **[unverified]** (the lab page omits it; likely
> 10.1038/s41588-022-01032-z — **do not cite that DOI without checking**).
> Numbers below from the medRxiv version, doi:10.1101/2020.11.08.20227983.

"cdr3-QTL" analysis: **628** healthy individuals (TCRβ CDR3s from blood, ~242,461 unique CDR3s each),
replication in **169** individuals with RNA-seq of sorted naïve CD4⁺ T cells. Strongest signal:
**HLA-DRB1 site 13 × CDR3 position 109**, MANOVA **P = 4.1 × 10⁻¹³⁸**, ~**9.4%** of variance in
amino-acid usage at that position. Conditional analysis within HLA-DRB1 implicates sites
**13, 71, 32, 74, 86, 30** (13/71/74 line the P4 pocket), with independent signals at HLA-B,
-DQB1, -DPB1, -DQA1; **81.1%** of significant associations are class II, concentrated at CDR3
positions 108–112. Higher CDR3 hydrophobicity at position 109 tracks HLA risk scores for RA, T1D and
CD. For us the key methodological point is the inverse of ours: they treat **V/J usage as a
confounder** (V gene as a fixed effect in the LMM) and note that cis variants in the TCR locus
"exclusively affect flanking positions correlated with J genes". **TRBV germline deletions/CNV are
not mentioned.** So the CDR3-level MHC story is taken; the **V-usage-level** story with germline CN
control is not.

### B.5 Is the TRBV4-3 deletion already known, and at what frequency by population?

**Yes — the deletion polymorphism is long-established and is not novel.** Its lineage:

- A classical insertion/deletion polymorphism in TCRB containing functional V segments was described
  decades ago (J Exp Med; PMC2191692) and a high-throughput assay for a "frequent
  insertion–deletion polymorphism in the human T-cell receptor beta locus" was published in
  *Immunogenetics* (2004), doi:10.1007/s00251-004-0684-z, PMID 15160257 — **full citation details and
  the frequencies reported there are [unverified]**, both the Springer page and PubMed were
  unreachable this session.
- Its functional consequences were examined in **"The impact of a large and frequent deletion in the
  human TCR β locus on antiviral immunity", *Journal of Immunology* 188(6):2742 (2012)** — volume,
  issue and page verified from the jimmunol.org URL; **authors, DOI and the reported frequencies are
  [unverified]** (the page is JS-only and ResearchGate blocked).
- Omer 2022 (B.1): co-deletion of TRBV4-3 + TRBV3-2 in 8/28 and 8/25 individuals in two small
  full-length cohorts (≈29–32% of individuals carry at least one deleted haplotype) — repertoire
  inference, no ancestry breakdown.
- Rodriguez/Watson (B.2): assembly-level confirmation of the TRBV6-2/3-2/4-3 deletion as polymorphic
  in 9 unrelated 1000G individuals.
- Luo S, Yu JA, Li H, Song YS. **Worldwide genetic variation of the IGHV and TRBV immune receptor
  gene families in humans.** *Life Science Alliance* **2**(2):e201800221 (2019),
  doi:10.26508/lsa.201800221, PMID 30808649 — 286 Simons Genome Diversity Project individuals,
  >100 populations. They treat TRBV4-3 only as part of the **{TRBV4-2, TRBV4-3} / {TRBV6-2, TRBV6-3}
  insertion**, calling it "a frequent polymorphism" **without giving a frequency in the text**; the
  per-region frequencies live only in Figure S8 (an image). Their quantified result is a different
  event: a TRBV5-8/7-8/6-9 three-gene deletion that "is the major variant among genomes sampled from
  the Americas but appears in only 5% of our sample globally". SGDP TRBV sample sizes: 73 West
  Eurasian, 27 Central Asian–Siberian, 45 East Asian, 49 South Asian, 22 Oceanian.

**The only per-population numbers I could verify** come from a very recent preprint:

> **gAIRR-wgs: An Algorithm to Genotype T Cell Receptor Alleles Using Whole Genome Sequencing Data.**
> bioRxiv doi:10.1101/2025.10.11.681280 (authors **[unverified]** — the PDF's author block did not
> render; the gAIRR-Suite lineage is Taiwanese, Li/Chen group **[unverified]**).

- Taiwan Biobank, n = 1,492: TRBV3-2 missing in **37.13%**, TRBV4-3 in **38.50%** of individuals;
  **38.20%** carry the homozygous deletion haplotype and **57.91%** the insertion haplotype.
- HPRC release 2, n = 232 assemblies: **34.48%** carry the homozygous deletion.
- Per-population homozygous-deletion frequencies, as stated: **AFR 12.93% (n=70), SAS 9.48% (n=54),
  EAS 3.45% (n=33)**. (These labels are the preprint's own — AFR/NAM/SAM/EAS/SAS/EUP; the EAS figure
  sits oddly against the 38% Taiwan Biobank number, so **treat these as provisional** and re-derive
  from our own data rather than citing them as the reference frequencies.)

**Bottom line for claim 4 and for the partner's TRBV4-3 result.** The deletion itself is known since
the 1990s–2000s, confirmed by long-read assembly in 2022, and now genotyped at population scale from
short-read WGS (2025 preprint). **What is not established is a well-powered, ancestry-stratified
allele-frequency estimate from phased long reads, and nobody has used it as the explanatory variable
for ancestry differences in observed TRBV usage.** Our novelty has to be framed there — "first
population-scale physically-phased genotyping of TRB germline deletions in a diverse cohort, linked
to repertoire usage" — not "we discovered the deletion". And the partner's SAS 43.5% vs EUR 16.7%
never-used TRBV4-3 figure sits in the same numeric range as published deletion carrier frequencies,
so the **deletion is the leading null hypothesis for that observation.**

### B.6 Adjacent work worth knowing

- **Gao K, Chen L, Zhang Y, et al. Germline-encoded TCR–MHC contacts promote TCR V gene bias in
  umbilical cord blood T cell repertoire.** *Frontiers in Immunology* **10**:2064 (2019),
  doi:10.3389/fimmu.2019.02064 (with a 2020 corrigendum, doi:10.3389/fimmu.2020.01752). Cord-blood
  repertoires; reports **HLA-DRB1 A71R associated with TRBV13**; finds HLA-A/-B (not HLA-C) the most
  influential class I loci, contrasting with the adult European picture; notes **TRBV20-1 and TRBV9**
  were previously described in European adults by Sharon et al. Undetected TRBV genes were simply
  **set to zero** for QTL analysis and germline deletions were not considered — the same confound.
- **"HLA alleles shape distinct biases in the usage preferences of TCR Vβ segments"**, PMC13261014 —
  appeared in search results as a recent paper directly on HLA × TRBV usage. **Title/PMCID only;
  journal, authors, year, DOI and content all [unverified]** (PMC blocked). **Follow this one up — it
  is the nearest potential competitor to the HLA × TRBV-usage framing.**
- **A template for multitudes: Germline immune polymorphism of the T cell receptor loci** (PMC9903699)
  — review of TR germline polymorphism, useful for framing; citation **[unverified]**.
- **Haplotype analysis of the TCRB locus by long-amplicon TCRB repertoire sequencing**,
  *Innovations in Pharmacology & Oncology* 2(4):137 — exists; content **[unverified]**.

---

## C. Scoop check

### C.1 The AoU long-read flagship preprint — verified, and reassuring

> **Population-scale Long-read Sequencing in the All of Us Research Program.**
> Garimella KV, Li Q, Wertz J, Lee SK, Cunial F, Huang Y, Mostovoy Y, Lorig-Roach R, English A, Su H,
> Levy S, Muzny DM, Berngruber C, Danzi MC, Harvey WT, LaPlante EL, Patterson K, Rozanski AN,
> Schwartz S, Shifaw B, Wang Y, Wong I, Xu IRL, Zaheri S, Zuchner S, Zheng X, Dugan-Perez S,
> Izydorczyk M, Mehta H, Gibbs RA, Lichtenstein L, Gupta N, Lennon N, Gabriel S,
> *All of Us Research Program Long Read Working Group*, Timp W, Doheny KF, Dutka T, Musick A, Wei C-L,
> Sedlazeck FJ, Schatz MC, Talkowski ME, Eichler EE.
> medRxiv **2025.10.02.25336942** (doi:10.1101/2025.10.02.25336942), posted ~2 Oct 2025.
> Title and author list verified from the medRxiv full-text page.

What it actually analyses:

- **1,027 self-identified Black or African American participants** with PacBio HiFi at ~8× (mid-pass);
  50 of these also ONT ~35×; **10,000** short-read WGS participants with EHRs for the PheWAS; a
  reference panel of **1,074 AoU+HPRC** samples (the HPRC count is given inconsistently as 44 and 47).
- **hifiasm v0.16.1** diploid assemblies, HiPhase read-based physical phasing, SHAPEIT4 statistical
  phasing, KAGE/GLIMPSE for panel-based genotyping; DeepVariant, PAV, PBSV, Sniffles2 for variants.
- Immune loci: the abstract lists HLA among "haplotypes relevant to disease risk (HLA) and
  pharmacogenomics (CYP2D6)", and the result is SV-level: **"723 LRS-unique SVs were identified
  within the major histocompatibility complex, including an 816 bp deletion within HLA-A"**. The
  deep-dive locus is **CYP2D6/CYP2D7**, not HLA.
- **Not present in the readable text: HLA allele calling or typing, novel HLA alleles, DQA1/DQB1
  phasing, any KIR analysis, any TRB/TRA/IG locus analysis, immune-gene copy number.** The fetched
  text was truncated (~106 kB, last ~6 kB unread) and supplements were not retrieved, so absence is
  **strongly indicated but not proven** — **re-check the supplement before submission.**

**Overlap with our claims:** claim 1 **none** (they do no HLA allele/haplotype calling; they use the
same HiPhase/SHAPEIT4 toolchain we are critiquing, which if anything sets up our result); claim 2
**none**; claim 3 **none**; claim 4 **none** (no TR loci). The genuine risk is not this paper but the
**AoU Long Read Working Group** itself: an immune-loci companion paper from the same consortium is the
obvious next output, and its author list is the list above.

### C.2 Immuannot — our own tool's paper, and the novel-allele baseline

> Zhou Y, Song L, Li H. **Full-resolution HLA and KIR gene annotations for human genome assemblies.**
> *Genome Research* **34**(11): 1931–… (2024). doi:10.1101/gr.278985.124 (the DOI is from the Genome
> Research early-access URL `content/early/2024/06/05/gr.278985.124`; **page range end [unverified]**).
> Preprint: bioRxiv doi:10.1101/2024.01.20.576452.

Types **36 HLA genes** (IPD-IMGT/HLA v3.52.0) and **17 KIR genes** (IPD-KIR v2.12), plus C4.
Applied to **212 haploid assemblies**: HPRC n=94, CPC n=114, CHM13.0, GRCh38.p14, CN-T2T n=2, plus
50 KIR reference contigs and 6 HLA reference assemblies. Results: **9,931 HLA+KIR genes**, of which
**4,068 (41.0%)** carried sequences novel relative to IPD; HLA: 7,519 genes / 2,673 unique alleles /
2,696 novel; KIR: 2,412 genes / 1,503 unique alleles / 1,372 novel; **2,664 distinct novel alleles**
with **92 novel protein sequences**. Novel-CDS mutation rate **5.1 × 10⁻⁵ /bp**, >10× the HPRC contig
base-error rate — their argument that most novel calls are real, though they still want experimental
validation. Population observations: HLA-DQB1 most diverse (mean pairwise divergence approaching
**5%**); HLA-A/-B/-C/-DRB1 above **2.5%**; African HPRC and Chinese CPC samples similarly diverse for
most HLA genes, attributed to balancing selection; HLA-Y deletion in **71.3%** of HPRC vs **87.3%** of
CPC; HTKU-del **9.6%** vs **16.4%**; MICA deletion **3.2%** HPRC, **1.9%** CPC, **11.9%** in 1000G MXL.

**Overlap:** claim 1 **none** (no DQA1~DQB1 haplotype analysis). Claim 2 **partial, and this is the
most important one to handle** — they already quantify novel-allele burden from assemblies and already
make an ancestry comparison (HPRC-African vs CPC-Chinese), at n=212 haplotypes. Our 17,255 haplotypes
across many ancestries is an order-of-magnitude extension with proper discovery curves, but **the
qualitative claim "assemblies reveal many novel HLA/KIR alleles, more in African ancestries" is
already in print.** Claim 3 **partial** (they report KIR allele diversity, not gene-content haplotype
frequencies by ancestry). Claim 4 **none**.

### C.3 Chen et al.(?) 2024 — KIR in pangenome assemblies

> **Genetic complexity of killer-cell immunoglobulin-like receptor genes in human pangenome
> assemblies.** *Genome Research* **34**(8): 1211–1223 (2024). doi:10.1101/gr.278358.123.
> Author byline **[unverified]** (the page's contributions section gives initials only: P-L.C.,
> J.S.H., Y-C.Y., Y-C.L., T-K.H., W-C.L., H-W.C., H-Y.L., C-L.H., C-Y.C.).

47 HPRC Year-1 phased diploid assemblies = **94 haplotypes**. Haplotype groups: **cA-tA 67.02%
(63/94)**, cA-tB 7, cB-tA 14, cB-tB 10. **15/94 (15.96%)** haplotypes carry KIR structural variation:
8 (8.51%) with ≥1 gene duplication (6 involving the KIR3DP1–KIR2DL4–KIR3DL1/S1 triad), 6 (6.38%)
lacking both KIR3DP1 and KIR2DL4, one lacking KIR3DL2; a **novel KIR2DS4/3DL1 fusion** in HG02630-M
and NA19240-P and a known KIR2DL3/2DP1 fusion in HG02723-P. **570 novel alleles**; 15 novel
nonsynonymous alleles and 8 SVs independently validated; **75/94 (79.8%)** haplotypes carry ≥5 novel
alleles and all 94 carry ≥2. They explicitly **do not** break KIR gene content down by ancestry and
write that "expanding the analysis of KIR gene diversity to include a broader range of populations
could help".

**Overlap:** claim 3 **partial-to-direct in method, none in result** — same tooling and same
assembly-based KIR gene-content framing, but n=94 haplotypes, no ancestry stratification, and they
name our gap as their own limitation. Claim 2 **partial** (KIR novel alleles). Claims 1, 4 **none**.

### C.4 SpecImmune 2025 — the most under-appreciated scoop risk

> **SpecImmune accurately genotypes diverse immune-related gene families using long-read data.**
> bioRxiv doi:10.1101/2025.02.04.636381 (posted Feb 2025). Title from search metadata; **authors and
> journal status [unverified]** (the PDF's author block did not render).

Genotypes **HLA (35 loci after LD filtering, 39 MHC-region genes), KIR (8), IG (67), TCR (61) and CYP
(13)** from PacBio and ONT. Cohorts: **1000 Genomes ONT n = 1,019** (26 populations, 5
superpopulations), HPRC 47 HiFi + 37 ONT, HGSVC 42 CLR + 9 HiFi, plus targeted TCR and CYP2D6
amplicons. Phasing via WhatsHap/LongPhase; in 9 HGSVC HiFi samples phase blocks averaged 4.5
heterozygous loci (max 30) and trio consistency was **84.6% (115/136)**; cross-locus haplotypes were
inferred **only for IG/TCR**. Population results: **African populations show higher IG/TCR
heterozygosity, with 91.2% of 80 significant loci enriched there**; East Asian populations lower
CYP2D6 activity; geographic patterning of HLA, KIR and CYP allele frequencies; strong LD among HLA,
KIR and CYP loci. Novel-allele calling is not a clear feature (novel variants are only visualisable).

**Overlap:** claim 2 **partial** (population-scale long-read HLA/KIR allele frequencies across 26
populations, n≈1,019 — but from reads, not assemblies, and without a novel-allele discovery curve).
Claim 3 **partial** (KIR across 5 superpopulations, 8 loci only). Claim 1 **none**, though note their
honest admission that cross-locus haplotypes were resolved only for IG/TCR — DQA1~DQB1 cross-gene
phasing is explicitly *not* what they did. Claim 4 **partial for TR genotyping in principle, none for
deletion genotyping from depth.**

### C.5 Other 2024–2026 items checked

- **High-resolution MHC haplotyping with long-read hybrid capture.** Glasenapp MR, Yee M-C,
  Symons AE, Garcia OA, Cornejo OE. medRxiv doi:10.64898/2026.03.27.26349549, v3 posted 13 May 2026.
  33 cell lines (32 analysed after HG01891 failed), PacBio Revio + ONT, custom Twist panel, tool
  **HLA-Resolve**: **97% two-field concordance (215/222)**, 97% three-field (193/198), **81%
  four-field**; four-field 131/162 vs 121/162 for StarPhase and SpecImmune. Class II is their weak
  spot — only **51%** zero-edit-distance to HPRC assemblies, "driven largely by HLA-DQA1 (2/13) and
  HLA-DRB1 (1/13)", attributed to "complex structural variants in introns and UTRs that are difficult
  for HLA-Resolve to genotype and phase"; they also flag a DQB1 exon-5-skipping splice variant and
  revise HG002 to DQB1\*05:01/\*03:02. No novel-allele reporting, **no DQA1~DQB1 pairing accuracy**,
  no direct comparison to statistical phasing (they only note that "imputation accuracy drops
  substantially in non-European populations"), and **no KIR** (a KIR probe panel is "under
  development"). **Overlap: claim 1 partial-adjacent (they establish that class II phasing is hard,
  which supports our premise rather than pre-empting it); claims 2–4 none.** Their stated KIR roadmap
  is a future-competition signal.
- **Zhou C, Gong T, Li S, Jin L, Fan S. Long-read sequencing reveals novel genetic polymorphisms in
  the major histocompatibility complex region and their impacts on the Han Chinese population.**
  *Science China Life Sciences* **68**(5):1400–1409 (2025). doi:10.1007/s11427-024-2742-y.
  100 Han Chinese, MHC by LRS and SRS at ~30×; **78,249 SNPs, 26.0% undetectable by short reads**;
  an LRS-based MHC reference panel with ~2.6× more variants than the SRS Han-MHC panel; PheWAS over
  26,024 phenotypes → 530 independent signals across 409 phenotypes, 28 previously unreported from
  the eight classical HLA genes. Abstract-only (paywalled): novel HLA **alleles** not specifically
  reported, no haplotype/phasing results visible, **no KIR**. **Overlap: claim 1 none; claim 2
  partial-adjacent (novel MHC variation, not novel HLA alleles); claims 3–4 none.**
- **Complex genetic variation in nearly complete human genomes** (HGSVC), *Nature* (2025),
  doi:10.1038/s41586-025-09140-6 — exists and is relevant background for assembly-based SV at immune
  loci; its specific HLA/KIR/TRB content was **not read** this session **[unverified]**.
- **MHConstructor: a high-throughput, haplotype-informed solution to the MHC assembly challenge**
  (PMC11142050) — short-read MHC assembly; citation and content **[unverified]**.
- No paper or preprint was found that reports **DQA1~DQB1 haplotype-group (G1/G2) mixing rates from
  physically phased population data, or that demonstrates statistical phasing manufacturing such
  haplotypes.** Multiple searches on this specific framing returned nothing. **Claim 1 appears
  genuinely unscooped** — the strongest of the four.

### C.6 Overlap matrix

| Work | Claim 1 (DQ G1/G2 phasing) | Claim 2 (novel-allele curves) | Claim 3 (KIR gene content × ancestry) | Claim 4 (TRB deletions from depth) |
|---|---|---|---|---|
| AoU LRS flagship (medRxiv 2025.10.02.25336942) | none | none | none | none |
| Immuannot / Zhou, Song & Li 2024 | none | **partial (closest)** | partial | none |
| KIR in pangenome assemblies 2024 (GR 34:1211) | none | partial | **partial, method-direct** | none |
| SpecImmune 2025 | none | partial | partial | partial |
| HLA-Resolve / Glasenapp 2026 | partial-adjacent (supports us) | none | none | none |
| Zhou et al. 2025 Han Chinese MHC LRS | none | partial-adjacent | none | none |
| Omer 2022 Genome Med | none | none | none | **partial (inference, not genome)** |
| Rodriguez/Watson TRA-TRB long-read | none | none | none | **partial (n=11, assembly)** |
| gAIRR-wgs 2025 preprint | none | none | none | **direct-adjacent (short-read, frequencies by population)** |
| Sharon 2016 / Russell 2022 / Ishigaki 2022 | none | none | none | none (all three ignore TRBV CNV) |

---

## D. Implications for us

**Risks, in order of seriousness.**

1. **Claim 4 is the most exposed.** The TRBV4-3/3-2/6-2 deletion has been known for ~25 years,
   assembly-confirmed in 2022, and a 2025 preprint (gAIRR-wgs) already publishes per-population
   homozygous-deletion frequencies from short-read WGS at n=1,492 + n=232. Frame our contribution as
   **(i)** depth-based CN genotyping validated against phased assemblies at ~12k scale,
   **(ii)** ancestry-resolved frequencies in a US-diverse cohort, and above all **(iii)** the link to
   repertoire phenotype — never as discovery of the deletion.
2. **The partner's TRBV4-3 ancestry result may be the deletion.** 27% never-used overall, SAS 43.5%
   vs EUR 16.7% — these sit squarely in the published deletion-frequency range. **Genotype the
   deletion first, then interpret.** If the deletion explains it, that is still a good result (and a
   cautionary one for every V-usage GWAS); if it does not, the residual
   is the interesting finding. Either way this ordering protects us from a referee doing it for us.
3. **Claim 2 is partly pre-empted by Immuannot's own paper.** Zhou/Song/Li already report 41% of
   HLA/KIR genes carrying novel sequence and already note African-ancestry diversity. Differentiate
   on **scale (17,255 haplotypes vs 212), proper rarefaction/discovery curves, per-ancestry
   saturation behaviour, and the operational consequence** (what fraction of participants carry an
   allele absent from the imputation reference — this is the clinically legible version and nobody
   has it).
4. **Claim 3's method is already demonstrated** at n=94 haplotypes. Our edge is exactly the gap that
   paper names as its own limitation: ancestry-stratified gene-content and structural-haplotype
   frequencies. Keep the framing on frequencies and ancestry, not on "KIR assemblies are complex".
5. **Consortium risk.** The AoU Long Read Working Group is the obvious author of an immune-loci
   companion paper. Worth a direct question up the chain about whether such a paper is in flight
   before we invest further in claims 2–3.
6. **Small-cell disclosure.** Novel-allele counts and KIR gene-content by ancestry both tend toward
   small per-group counts. This is the open item in `context/DECISIONS.md`; flag any table before it
   leaves the repo.
7. **Unverified-citation risk.** Several DOIs/author lists above are marked **[unverified]** because
   PMC, PubMed, Europe PMC, Springer, eLife and jimmunol were unreachable from this container. **Every
   one must be checked against the publisher record before it enters a manuscript.** Specifically:
   the Ishigaki DOI, the Omer author list, the Rodriguez/Watson published venue, the 2004
   *Immunogenetics* and 2012 *J Immunol* TRB-deletion papers, the gAIRR-wgs authors, and
   PMC13261014.
8. **One prompt-level correction to carry forward:** Russell et al. 2022 is **eLife 11:e73475**, not
   *Science Immunology*. Citing it as the latter would be an immediate credibility hit.

**What to cite, and for what.**

- **Sharon 2016 (Nat Genet 48:995–1002)** — the founding V-usage × MHC result, and our explicit point
  of departure. Cite their own stated limitations (imputed European HLA; bulk blood; LD-limited
  fine-mapping; their call for "longer-read sequencing of sorted cells") as the motivation for our
  design. Also cite their **TRBV24-1 / ZNF443** hit.
- **Russell 2022 (eLife 11:e73475)** — productive-vs-non-productive MHC association as the
  thymic-selection control; TRBV28 RSS cis effect; independent TRBV24-1/ZNF443 replication; and as
  the second V-usage GWAS whose cohort was effectively European.
- **Ishigaki 2022 (Nat Genet 54:393–402)** — the CDR3-level HLA restriction story; cite to delimit
  scope (they condition V usage away; we study it) and for the DRB1 13/71/74 P4-pocket positions.
- **Omer 2022 (Genome Med 14:2)** and **Rodriguez/Watson (bioRxiv 2022.05.24.493244)** — the
  established existence, extent and haplotype context of TRB germline deletions, including the
  TRBV7-2\*02 tag; cite for prior art and for the repertoire-inference-vs-genome distinction.
- **gAIRR-wgs (bioRxiv 2025.10.11.681280)** — must be cited and compared to directly; it is the
  nearest thing to a population-frequency competitor for claim 4.
- **Immuannot (Genome Res 34:1931)** — our method, and the novel-allele baseline we must exceed.
- **KIR pangenome (Genome Res 34:1211–1223)** — the KIR baseline; quote their ancestry limitation.
- **SpecImmune (bioRxiv 2025.02.04.636381)** — population-scale long-read immune genotyping;
  acknowledge, and lean on their "cross-locus haplotypes only for IG/TCR" admission.
- **HLA-Resolve (medRxiv 2026.03.27.26349549)** — independent evidence that class II (DQA1, DRB1)
  genotyping and phasing is the hard part; good support for claim 1's premise.
- **AoU LRS flagship (medRxiv 2025.10.02.25336942)** — the mandatory cohort/resource citation.
- **Luo 2019 (LSA 2:e201800221)** — worldwide TRBV CN variation; cite for the global-diversity frame
  and for the TRBV5-8/7-8/6-9 Americas deletion.

**What to differentiate on, in one line each.**

- **Claim 1** — lead with it. No prior work reports DQA1~DQB1 haplotype-group mixing from physically
  phased population data, and two independent 2026 sources confirm class II phasing is where
  short-read and capture methods fail. This is our cleanest novelty.
- **Claim 2** — differentiate on scale, discovery-curve methodology, and the imputation-reference
  consequence, not on the existence of novel alleles.
- **Claim 3** — differentiate on ancestry-stratified frequencies, filling the gap the 2024 KIR
  pangenome paper names itself.
- **Claim 4** — differentiate on phased-assembly-validated CN genotyping at scale **and** on the
  repertoire link; concede the deletion's prior discovery up front, and use it to explain, not merely
  to describe, the ancestry differences in TRBV usage.
