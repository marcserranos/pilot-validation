# Call summary — 2026-10-06, Progress update #10

*Written 2026-10-11 from a partial transcript Marc pasted into the project thread (it starts
mid-call, during Cole's PheWAS slides). Present: Cole Shanks, David Bonet Sole, Marc, and one more
person who left early (unnamed in the transcript). Date confirmed by Marc on 2026-10-11; deck
`REPORT #10`.*

## Transcription corrections

| In the transcript | Read as | Why |
|---|---|---|
| "cure genes", "curogenes" | **KIR genes** | Same mis-hearing as call #8; context is Marc's KIR work and the four framework genes |
| "Vcodex", "Vcodes" | **PhecodeX** (phecodes) | Hierarchical phenotype coding Cole switched to from raw ICD codes |
| "VWAS" | **PheWAS** | Phenome-wide scan of repertoire features |
| "VCRs" | **BCRs** | Paired with TCRs throughout |
| "KMRs" | **k-mers** | As in call #8 (BenchRep-T features) |
| "LEI" | **LAI** (local ancestry inference) | David's local-ancestry point |
| "Bonn-Faroni" | **Bonferroni** | |
| "Jonathan Pritchard 2016 Nature Genetics paper" | Very likely **Sharon et al. 2016, Nat Genet**, "Genetic variation in MHC proteins is associated with T cell receptor expression biases" (Pritchard lab) | <1,000 people, RNA-seq + imputed HLA, GWAS of V/J usage. **(verify)** |
| Cole: KIR on "chromosome 14" | KIR is on **chr19** | Same slip as call #8; chr14 holds TRA/IGH |

## What Cole presented (his stat-gen and repertoire work)

1. **PheWAS of repertoire features (TRB and BCR), now on PhecodeX.** Striking excess of renal and
   cardiovascular phenotypes for TRB diversity (hypertension in many forms, acute and chronic
   kidney failure, end-stage renal disease, heart failure, angina, glomerular disease). Also eye
   disorders, rheumatoid arthritis *without* rheumatoid factor (a possible endotype), leukemia /
   mature B-cell lymphoma, and the expected CMV and HIV signals. Pregnancy dropped out after the
   switch from ICD codes to PhecodeX, probably a coding artefact.
   - Cole: "If we can make even a small sliver of a story out of this, it would go big." Framing:
     can the repertoire classify hypertension, heart failure or chronic kidney disease? No
     causality claim needed.
   - David: **collapse the phecode hierarchy** (hypertension / hypertensive crisis / hypertensive
     heart disease are one test, not three), fewer correlated tests, cleaner figure. Cole agreed:
     collapsed version with Bonferroni for a main figure, full scan as a "how we got here"
     supplement.
   - Open: hypertension, kidney and heart disease are causally entangled (CKD usually follows
     hypertension); worth poking at, not solving.
2. **Public clonotype × HLA associations.** Common clonotypes (CDR3 + V/J), presence/absence
   regressed on each HLA allele with covariates. Hits look EBV-specific (e.g. an HLA-B\*08
   clonotype with a very large odds ratio). Cole: **treat DQ as one phased unit (DQA1~DQB1
   haplotype), not two independent genes.** Next: conditional fine-mapping on the top ~10 alleles.
3. **GWAS of V/J usage (replicating Sharon et al. 2016 at scale).** Expected hits at the TRB locus
   and the MHC; an unexplained hit on chr4 (3–4 candidate genes; may not survive correcting for
   the ~60–80 usage traits). Next: replace SNPs with HLA alleles, then fine-map which HLA alleles
   drive V/J usage, and connect to ancestry: V/J usage predicts ancestry well, and so does HLA.
   Suggested figure: PCA of per-person V/J frequencies coloured by dominant ancestry.
4. **Ancestry handling.** David raised admixture at the MHC (e.g. excess European HLA ancestry in
   Indigenous American–European admixed people, from post-contact selection; the opposite,
   African HLA, in African Americans). Cole: real, but **no local ancestry inference for now**;
   use ancestry-purity filtering (90–95%) as Marc has been doing, revisit only if it becomes
   apparent.

## Marc's questions and the answers

- Can a disease leave strong clonal patterns without changing diversity/entropy? **Yes.**
  Disease-driving clonotypes in most autoimmune disease (and allergy, celiac) are very rare and
  mostly unobservable in bulk RNA-seq. Entropy effects flag diseases worth a closer look, but
  signal may sit in ratios/expansions or in **motifs shared by similar CDR3s** (possibly found
  with XGBoost on k-mers).
- KIR structural/framework genes occasionally missing, disappearing at stricter phasing → Cole:
  **misassemblies**, probably low coverage. Marc proposed **per-region coverage**, not just
  per-person: Cole agreed, e.g. drop someone at 17–20× genome-wide but 8–10× over the MHC.

## Decisions

1. Phenotypes are collapsed along the PhecodeX hierarchy before testing (David, agreed by Cole).
2. DQ is represented as a phased DQA1~DQB1 unit in association tests.
3. No local ancestry inference for now; ancestry-purity filtering instead.
4. No single story locked yet; the renal/cardiac repertoire signal is the leading candidate hook.

## Action items

**Marc** (his stated focus, agreed by Cole):
- **QC and phasing**, including a **regional coverage filter** (read depth over the MHC and the
  KIR region vs the person's genome-wide depth) and checking whether it removes the missing
  framework-gene calls.
- **KIR**: understand the biology better; KIR gene and allele frequencies by ancestry as a likely
  Figure 1 addition.
- Later: **help curate diseases for the classification work.**

**Cole**: collapse phecodes and redo the PheWAS figure; check why pregnancy drops out under
PhecodeX; fine-map clonotype × HLA and HLA × V/J usage by conditional regression; look at the
chr4 V/J-usage hit; link the Pritchard-lab paper.

**David**: advise on the tabular models once they start.
