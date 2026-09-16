# Literature review: HLA diversity, selection statistics, and ARS residue definitions

Scope: Cole asked us to find citable methodology for (a) nucleotide/AA diversity contrasting
peptide-binding-residue (PBR/ARS) vs non-ARS codons, (b) balancing selection at class I vs
class II, (c) whether/why class II is more diverse and more ancestry-differentiated than
class I, and (d) the standard statistics (dN/dS on ARS, Tajima's D, Ewens–Watterson/Slatkin,
Fst on alleles, HKA) and their pitfalls at a multiallelic locus. This is a literature survey
only — no new analysis was run. Everything below comes from a fetched abstract/full text
unless flagged **(unverified)**.

---

## 1. Canonical citations

### ARS vs non-ARS dN/dS

- **Hughes AL, Nei M (1988).** "Pattern of nucleotide substitution at major histocompatibility
  complex class I loci reveals overdominant selection." *Nature* 335:167–170.
  doi:10.1038/335167a0. **The** founding paper: dN significantly exceeds dS at the
  antigen-recognition-site (ARS) codons of classical class I genes (HLA-A/B/C), while dN≈dS≈~0
  outside the ARS. Interpreted as evidence for balancing (overdominant) selection maintaining
  diversity specifically at peptide-contact residues. The ARS codon set they used was taken
  from the Bjorkman/Saper crystal structure (see §2).
- **Hughes AL, Nei M (1989).** "Nucleotide substitution at major histocompatibility complex
  class II loci: evidence for overdominant selection." *PNAS* 86:958–962. The class-II analog
  of the 1988 paper — dN > dS at the class II (DRB1/DQ) ARS, dN≈dS outside it, again read as
  overdominant selection. Cited together with the 1988 paper, this pair is the standard
  citation for "ARS dN/dS shows selection at HLA."
- **Hughes AL, Hughes MK, Watkins DI (1990).** "Positive Darwinian selection promotes charge
  profile diversity in the antigen-binding cleft of class I MHC molecules." *Mol Biol Evol*
  7:515–524. Successor paper refining which ARS sub-region (binding cleft vs TCR-facing vs
  outward-facing residues) drives the signal for class I. **(unverified — found via search
  summary, not fetched in full)**.
- **Successor / modern reworking:** a bioRxiv preprint, "Heterogeneity of dN/dS ratios at the
  classical HLA class I genes over divergence time and across the allelic phylogeny"
  (bioRxiv 10.1101/008342) revisits the Hughes & Nei approach with a phylogenetic framework
  and reportedly finds the ARS dN/dS signal is not uniform across lineages/time. I could not
  fetch the full text (rate-limited) to confirm authors/year/journal of record —
  **(unverified, needs a direct re-check before citing)**.

### Ewens–Watterson / Slatkin exact test (allele-frequency homozygosity)

- **Watterson GA (1978).** "The homozygosity test of neutrality." *Genetics* 88:405–417.
  Establishes the observed-homozygosity-vs-Ewens-sampling-formula test; shown to be powerful
  against symmetric overdominance (i.e., balancing selection with heterozygote advantage) —
  exactly the alternative hypothesis of interest for HLA. **(unverified — not directly
  fetched, cited from the description in the Solberg et al. 2008 review below, which is
  itself confirmed)**.
- **Slatkin M (1994).** "An exact test for neutrality based on the Ewens sampling
  distribution." *Genetics Research* 64:71–74; correction **Slatkin M (1996)**, *Genetics
  Research* 68:259–260 (the 1994 test used the ordered rather than unordered allele
  configuration probability; the correction fixes this — this is the version actually
  implemented in most HLA software, so cite both together).
- **Solberg OD, Mack SJ, Lancaster AK, Single RM, Tsai Y, Sanchez-Mazas A, Thomson G (2008).**
  "Balancing selection and heterogeneity across the classical human leukocyte antigen loci: a
  meta-analytic review of 497 population studies." *Human Immunology* 69(7):443–464.
  **(verified, fetched).** Applies the Ewens–Watterson homozygosity test and a normalized
  deviate (Fnd, comparable across loci/sample sizes) to ~497 population samples (~66,800
  people) across A, C, B, DRB1, DQA1, DQB1, DPA1, DPB1. This is the single best "we used
  standard methodology, here's who else did" citation for an EW/Fnd-based homozygosity scan
  across loci and populations — exactly Cole's ask.

### Tajima's D at HLA

- **Tajima F (1989).** "Statistical method for testing the neutral mutation hypothesis by DNA
  polymorphism." *Genetics* 123:585–595. The origin of the statistic itself.
  **(unverified — not fetched this session, standard citation)**.
- **Meyer D, Single RM, Mack SJ, Erlich HA, Thomson G (2006).** "Signatures of Demographic
  History and Natural Selection in the Human Major Histocompatibility Complex Loci."
  *Genetics* 173(4):2121–2142. **(verified, fetched).** Applies a battery of methods
  (including allele-frequency-based neutrality tests) to HLA-A, -B, -C, -DRB1, -DQA1, -DQB1
  across populations, explicitly contrasting HLA against non-HLA genomic background loci to
  separate demographic history from selection; finds 4 of 6 loci deviate from neutral
  expectations in a way not explained by demography alone. Good citation for "how do you
  separate demography from selection at HLA," a pitfall Cole will want addressed.
  (I could not get the fetch to enumerate which named test — Tajima's D vs EW vs Fst — maps to
  which locus in this paper; treat the general framing as verified, the test-by-test
  breakdown as **(unverified)**.)

### Fst / allele-level differentiation

- **Brandt DYC, César J, Goudet J, Meyer D (2018).** "The Effect of Balancing Selection on
  Population Differentiation: A Study with HLA Genes." *G3: Genes|Genomes|Genetics*
  8(8):2805–2815. **(verified, fetched).** Directly addresses the "does balancing selection
  reduce or increase Fst" question. Method: Weir & Cockerham (1984) Fst estimator, "ratio of
  averages" (not "average of ratios") to avoid MAF-driven bias, MAF-matched bins for
  comparison against genome background. Finding: HLA class I and class II SNPs show **reduced**
  Fst vs matched genome-wide controls once MAF is properly controlled — except HLA-DPA1/DPB1,
  which show *elevated* Fst and are flagged as under directional rather than balancing
  selection. This is the citation to use for "why don't we just compute vanilla Fst on HLA
  alleles" and for the DP-locus caveat.
- **Weir BS, Cockerham CC (1984).** "Estimating F-statistics for the analysis of population
  structure." *Evolution* 38:1358–1370. The underlying Fst estimator used above.
  **(unverified — standard citation, not independently fetched)**.
- **Arrieta-Bolaños E, Hernández-Zaragoza DI, Barquera R (2023).** "An HLA map of the world: A
  comparison of HLA frequencies in 200 worldwide populations reveals diverse patterns for
  class I and class II." *Frontiers in Genetics* 14:866407. **(verified, fetched).** Notably
  they *avoid* Fst and instead use PCoA on Euclidean allele-frequency distances plus Barrier
  v2.2 genetic-discontinuity analysis on 200 populations — worth citing both for the
  class I (geography-tracking) vs class II (weaker geography, different clustering) contrast
  and as a worked example of the "don't just use Fst on multiallelic loci" alternative.

### HKA test

- **Hudson RR, Kreitman M, Aguadé M (1987).** "A test of neutral molecular evolution based on
  nucleotide data." *Genetics* 116:153–159. Origin of the polymorphism-vs-divergence
  neutrality test; per search summaries it has been applied to HLA with the strongest signal
  in/around DQA1/DQB1 exon 2. **(unverified — origin paper not independently fetched this
  session; the "applied to HLA, strongest in DQA1/DQB1" claim is from a secondary summary,
  not a primary source I read directly — treat as (unverified) until we pull the actual
  application paper)**.

### Long-read / large-cohort HLA typing

- **Robinson J et al. (2020).** "IPD-IMGT/HLA Database." *Nucleic Acids Research*
  48(D1):D948–D955. **(verified via search, standard reference database citation)** — cite
  for allele nomenclature/4-field naming conventions and as the sequence resource against
  which novel alleles are called.
- PacBio HiFi HLA-typing methodology and multi-cohort benchmarking: search results surfaced
  **SpecImmune** (benchmarked across 1000 Genomes ONT, HPRC PacBio/ONT samples) and
  **HLA-Resolve** (HiFi-optimized, benchmarked against GIAB/HPRC/IHWG references) as the
  current long-read HLA-typing tool papers, both on bioRxiv/medRxiv. I did not fetch either
  in full, so treat exact author lists/journal-of-record status as **(unverified)** — these
  are the right search terms to chase down for the "prior long-read HLA cohort work" citation,
  but neither is yet at the scale of our ~12k-assembly cohort based on what I could confirm
  (SpecImmune's PacBio arm was 47 samples; HLA-Resolve's validation set was 32 samples).
  **This is a real gap**: I did not find a published long-read HLA study at anywhere near our
  N, which supports point 4 below (scale is a genuine differentiator for us).

---

## 2. ARS residue definitions (the most load-bearing part of this review)

**Class I (HLA-A/B/C, exons 2–3 → α1/α2 domains).**
The residue-level source is the crystallographic identification of the peptide-binding groove:

- **Bjorkman PJ, Saper MA, Samraoui B, Bennett WS, Strominger JL, Wiley DC (1987).** "The
  foreign antigen binding site and T cell recognition regions of class I histocompatibility
  antigens." *Nature* 329:512–518 (companion structure paper: "Structure of the human class I
  histocompatibility antigen, HLA-A2," *Nature* 329:506–512). **(verified — paper located and
  described via search; groove-forming residue table itself not independently re-extracted
  from the primary PDF, since Nature blocked automated fetch this session)**. This structure
  is where the physical groove residues were first mapped; Hughes & Nei (1988) then translated
  this into the ARS codon list used for their dN/dS analysis.
- **Parham P, Lomen CE, Lawlor DA, Ways JP, Holmes N, Coppin HL, Salter RD, Wan AM,
  Ennis PD (1988).** "Nature of polymorphism in HLA-A, -B, and -C molecules." *PNAS*
  85:4005–4009. Classifies the 39 known class I sequences at the time into **20 highly
  variable amino acid positions** (clustered in the peptide-binding groove) versus 71 largely
  invariant positions. I could not re-fetch the PDF (403) to pull the literal position numbers
  this session, so the residue list itself is **(unverified)** — but the paper and its
  20-vs-71 classification is confirmed via search/abstract, and it's the standard companion
  citation to Bjorkman & Saper for "which class I residues are actually polymorphic /
  peptide-contacting." **Action item:** get the Parham 1988 Table (or its restatement in any
  IPD-IMGT/HLA or Marsh "HLA Facts Book" appendix) directly — this is the single most useful
  document to physically read for a defensible class-I ARS/PBR residue list, and it's not
  paywalled everywhere (try ResearchGate/PMC mirror).
- Practically, most downstream HLA dN/dS and diversity papers cite the **Hughes & Nei (1988)**
  ARS list itself (their Table 1) rather than re-deriving it from Bjorkman/Parham — i.e., cite
  Hughes & Nei for "the codon list used," and Bjorkman/Saper + Parham for "where that list
  structurally comes from."

**Class II (exon 2 of DRB1, and of DQA1/DQB1, DPA1/DPB1 → β1/α1 domains).**

- **Brown JH, Jardetzky TS, Gorga JC, Stern LJ, Urban RG, Strominger JL, Wiley DC (1993).**
  "Three-dimensional structure of the human class II histocompatibility antigen HLA-DR1."
  *Nature* 364:33–39, and the companion peptide-bound structure **Stern LJ, Brown JH,
  Jardetzky TS, Gorga JC, Urban RG, Strominger JL, Wiley DC (1994).** "Crystal structure of the
  human class II MHC protein HLA-DR1 complexed with an influenza virus peptide." *Nature*
  368:215–221. These define the class II groove: floor formed by a β-sheet from both α1 and
  β1 domains, flanked by two α-helices (one contributed by each chain), groove open at both
  ends (unlike class I). **(verified — papers located and described via search; exact
  pocket-lining residue numbers not independently re-extracted this session)**.
- **Bondinas GP, Moustakas AK, Papadopoulos GK (2007).** "The spectrum of HLA-DQ and HLA-DR
  alleles, 2006: a listing correlating sequence and structure with function."
  *Immunogenetics* 59:539–553. **This is the paper to cite for an actual, structurally
  cross-referenced, per-pocket residue list for DRB1/DQA1/DQB1** — it proposes a common
  numbering scheme across the class II loci so that, e.g., "pocket 1," "pocket 4," "pocket 6,"
  "pocket 9" residues are directly comparable between DR and DQ molecules. I confirmed the
  paper exists and its stated purpose via search but did **not** get the full pocket-residue
  table into this session (paywalled/redirect-blocked) — **(unverified pending direct read)**.
  This is the second concrete action item: get full-text access to Bondinas 2007 (library
  proxy or a review that reproduces its table, e.g. the Journal of Immunology 2009 "pocket 6
  and pocket 9" paper which explicitly cites and uses Bondinas numbering) before finalizing a
  class-II ARS residue set.
- A commonly used *shortcut* residue list independent of Bondinas: the classical
  "shared epitope" and individual pocket-residue literature (e.g., DRB1 positions 11, 13, 71,
  74 for pocket 4; position 57 of DQB1 for pocket 9 in T1D) gives well-validated single
  positions, but these are disease-association residues, not a systematic ARS/PBR list —
  don't substitute them for a real pocket-by-pocket table.

**Bottom line for computation:** we should not compute ARS/non-ARS diversity from a
Wikipedia-style guess. The defensible chain of citation is Bjorkman & Saper (1987) / Brown et
al. (1993) for the structural definition → Hughes & Nei (1988, 1989) for the codon list as
actually operationalized for dN/dS → Bondinas et al. (2007) for a modern, cross-locus-comparable
class II pocket numbering. **Before running any ARS-vs-non-ARS analysis, get the literal
residue-number tables from Hughes & Nei (1988) Table 1, Parham et al. (1988), and Bondinas et
al. (2007)** — I flagged these as unverified because I could not get full-text access this
session (Nature/PNAS/Springer all blocked the automated fetch), not because the lists don't
exist. This is worth a follow-up WebFetch pass or a manual PDF pull.

---

## 3. Pitfalls of these statistics on allele-level, biobank-scale HLA data

| Issue | What goes wrong | What the literature does about it |
|---|---|---|
| **Ascertainment / reference bias** | Short-read/array-based HLA calls are biased toward reference-like alleles; allele frequencies (and hence homozygosity/Fst) get skewed. Confirmed source: **Bettinotti/Diskin-style finding** reported for 1000 Genomes: mapping bias overestimates reference allele frequencies at HLA in 1000G Phase I data — ~18.6% of SNP genotype calls wrong, ±0.1 allele-frequency error at ~25% of HLA SNPs. **(source: G3 2015 paper "Mapping Bias Overestimates Reference Allele Frequencies at the HLA Genes in the 1000 Genomes Project Phase I Data," located via search, not independently fetched — (unverified) on the exact numbers, but the phenomenon and paper are real and citable)**. | Use assembly/long-read calls (which is exactly our advantage) rather than reference-guided short-read pipelines; population-specific reference panels for anywhere imputation is still used (e.g. FinnGen panel, NAR Genomics & Bioinformatics 2020). |
| **Multi-allelic vs biallelic statistic design** | Tajima's D, HKA, and most Fst estimators were derived for biallelic SNP data; applying them naively to a locus with hundreds of alleles per gene either requires translating to per-site nucleotide diversity (losing the allele-as-unit information Cole actually cares about) or requires allele-frequency-specific tools (Ewens–Watterson/Slatkin, which were purpose-built for exactly this). | Ewens–Watterson/Slatkin homozygosity tests are the standard workaround because they operate on the allele-frequency spectrum directly rather than assuming biallelic sites; Brandt et al. (2018) show that even Fst needs MAF-matched control to be interpretable at a multiallelic, high-diversity locus like HLA. |
| **Recombination / gene conversion within the locus** | HLA (especially DRB1, and class I between allelic lineages) undergoes intragenic and intergenic gene conversion, violating the single-genealogy assumption behind Tajima's D, HKA, and most coalescent-based selection tests; it also generates new "mosaic" alleles that inflate apparent diversity without new mutation. | Sequence-feature studies of DRB1 exon 2 (elevated CpG content, conserved flanking motifs facilitating conversion) are cited as the mechanistic explanation; papers applying HKA/Tajima's D to HLA generally flag recombination as a caveat on the interpretation of significant results as "selection" rather than "recombination/conversion echo." |
| **Admixture / population structure vs selection** | A locus with strongly non-neutral allele-frequency distributions in an admixed cohort (like ours) can show apparent departures from neutrality (positive Tajima's D, EW homozygosity deficit) driven by mixing multiple differently-selected or differently-drifted subpopulations rather than by selection acting within any one population. | Meyer et al. (2006, *Genetics*) explicitly contrast HLA loci against non-HLA background loci genotyped in the same individuals to separate demographic signal from locus-specific selection — this "HLA vs matched genomic background, same samples" design is the standard control and is directly usable with our AoU genome-wide data. |
| **Relatedness / cryptic family structure in biobanks** | Allele-frequency-based neutrality tests (EW, Tajima's D analogs) assume independent draws from a population; cryptic relatedness (common in large biobanks) inflates apparent homozygosity/rare-allele sharing and can bias both the observed statistic and its null distribution. | This is a general biobank/GWAS confound (extensively documented for association testing, e.g. work on cryptic relatedness and population stratification) rather than an HLA-specific literature; **(unverified whether anyone has applied this specifically to HLA EW/Tajima's-D testing — worth checking, but the general mitigation — kinship filtering / one-per-family subsampling before computing allele-frequency statistics — transfers directly)**. |
| **DP loci behaving differently** | DPA1/DPB1 are sometimes lumped in with "HLA class II is under balancing selection," but Brandt et al. (2018) and Solberg et al. (2008) both separately show DP loci deviate — DPB1 compatible with near-neutral expectations (Solberg), and DPA1/DPB1 showing *elevated*, not reduced, Fst consistent with directional rather than balancing selection (Brandt). | Don't report "class II" as a monolith; stratify DRB1/DQA1/DQB1 vs DPA1/DPB1 in any comparison. |

---

## 4. What is genuinely new for us vs. re-demonstration

Blunt assessment:

- **Re-demonstration, not new:** "Class I and class II show elevated diversity at ARS vs
  non-ARS codons" (Hughes & Nei, 1988/1989, and dozens of follow-ups over 35+ years) and
  "class II shows stronger/more heterogeneous balancing-selection signal than class I" (Solberg
  et al. 2008, Meyer et al. 2006) are both extremely well-established. Replicating them in our
  cohort is a **sanity check / methods-validation step**, not a finding — useful to show the
  method works on our novel data type, not publishable on its own as biology.
- **Re-demonstration but with a genuinely better instrument:** every ARS dN/dS and HKA/Tajima's
  D study on HLA to date that I could locate is built on short-read genotyping or Sanger/
  amplicon-based allele calls, largely European/East Asian-heavy reference panels, and sample
  sizes rarely exceeding a few thousand for full-CDS-resolution data (Solberg et al.'s 66,800
  individuals is allele-frequency data pooled across many small studies, not one cohort with
  full CDS + phasing). Running the *same* classical statistics but on ~12,000 phased,
  4-field, full-CDS, long-read assemblies with AoU's ancestry breadth is a real
  precision/resolution advance even though the qualitative conclusion (yes, ARS is more
  diverse; yes, class II differs from class I) is expected. This is worth including as a
  "biggest cohort of this kind analyzed this way" methods point, not as a "we discovered new
  biology" point.
- **Actually new, if it works:** (a) whether the ~1,400 novel protein alleles we've found are
  enriched or depleted at ARS vs non-ARS positions relative to known alleles — nobody has this
  data at our scale, and the direction isn't obvious a priori (novel *because* rare mutation
  at a hypervariable site, or novel *because* sequencing/assembly resolution catches variants
  that were always there but uncallable); (b) directly measuring Fst/differentiation *at the
  novel-allele layer* across AoU's self-reported ancestry groups, which by construction cannot
  exist in any legacy dataset built on known-allele imputation panels; (c) testing whether
  Brandt et al.'s DP-locus exception (directional rather than balancing selection) replicates
  or looks different when the allele set includes rare/novel DP alleles not previously
  characterized.
- I could not find a directly comparable long-read biobank-scale HLA selection-statistics paper
  in this search pass (SpecImmune/HLA-Resolve papers benchmark accuracy on tens of samples, not
  thousands, and don't appear to run population-genetics selection tests) — so the "n=12,000,
  fully phased, 4-field, long-read" combination applied to any of these classical tests is
  plausibly a first, but I'd want a more targeted search (e.g., "All of Us HLA long read",
  "biobank-scale HLA selection Tajima") before claiming that in the paper.

## 5. Recommended plan (ranked)

1. **ARS vs non-ARS nucleotide/AA diversity (π, Watterson's θ, or simple heterozygosity) split
   by codon class, per gene, per ancestry group.** Lowest lift given we already have full CDS
   and 4-field calls; directly answers Cole's stated interest; requires nailing down the
   residue lists first (§2 action items).
2. **Ewens–Watterson / Slatkin homozygosity test per classical locus (A, B, C, DRB1, DQA1,
   DQB1, DPA1, DPB1), overall and per ancestry group**, replicating and extending Solberg et
   al. (2008)'s design but on our single large, consistently-typed cohort instead of pooled
   heterogeneous studies — directly comparable to a paper we can cite.
3. **ARS-restricted vs whole-CDS dN/dS comparison per gene**, reproducing the Hughes & Nei
   design as a validation step, then testing whether novel alleles specifically elevate or
   dampen the ARS signal.
4. **Fst (Weir & Cockerham, MAF-matched) at the allele level across AoU ancestry groups**,
   following the Brandt et al. (2018) "ratio of averages, MAF-binned" protocol, stratifying
   DR/DQ vs DP explicitly rather than pooling all class II — tests whether their DP-exception
   finding replicates with novel alleles included.
5. **Novel-allele enrichment analysis at ARS vs non-ARS codons** — genuinely new, cohort- and
   technology-specific; frame as the paper's actual selling point rather than a repeat of
   classical HLA population genetics.
6. **(Lower priority / exploratory) HKA test contrasting HLA loci against matched non-HLA
   genomic background** in the same individuals, following the Meyer et al. (2006) design —
   valuable for the demography-vs-selection pitfall but needs a well-matched neutral
   comparison locus set and is more work to set up than 1–5.

---

*Caveat on this document*: several primary sources (Nature, PNAS, Springer) blocked automated
full-text fetch this session (403/redirect-to-login), so a handful of claims above — mainly the
literal ARS/pocket residue numbers and the exact Hudson–Kreitman–Aguadé HLA application paper —
are based on search-result summaries rather than a directly read primary source, and are marked
**(unverified)** accordingly. Before we compute anything, at minimum pull full text (via
Stanford library proxy) for: Hughes & Nei (1988) Table 1, Parham et al. (1988) PNAS, and
Bondinas et al. (2007) Immunogenetics.
