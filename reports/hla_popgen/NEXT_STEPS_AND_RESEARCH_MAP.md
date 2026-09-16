# Next steps and research map — HLA long-read paper + HLA×TCR/BCR integration

*Written 2026-09-15 on branch `fig1-drafts-and-research-map` (cut from `origin/main` @ `30b14ac`).
Inputs: the full `context/` set, `COMPREHENSIVE_REPORT.md`, every `reports/hla_popgen/*` report
on main, Aleix's `origin/aleix/hla-resolve-phase1:aleix/RNA-seq/` (through `b9f7e65`), and the
~2026-09-14 supervisor meeting notes. Citations marked **(verify)** come from memory, not from a
fetched source. Check them before they go in a manuscript.*

---

## 0. Where the project stands, in five lines

1. **HLA calling is done at scale.** Immuannot on ~12,233 AoU long-read assemblies gives 65 genes,
   4-field resolution, phased and contig-aware. Pipeline, QC and 20+ analyses are on main.
2. **The validated findings:**
   - Phasing: 545/545 relative pairs show zero phase switches, and 94.8% of relative gene
     comparisons are byte-identical.
   - Selection: CDS π exceeds non-coding π for A/B/C/DRB1/DPB1 but not for the DRA/F controls.
   - Disease signal: B\*27/AS replicates (OR 9.5, CMH p=4×10⁻¹⁰).
   - Reference bias: a template-distance gradient runs from EUR (84.6% exact) down to AFR/EAS
     (~78%).
3. **The negatives, established honestly:**
   - Unsupervised HLA embeddings don't separate disease.
   - Supervised UMAP failed held-out testing (AUROC 0.48).
   - Chao2 saturation estimates are low-confidence.
4. **Aleix's side:** TRUST4 works on AoU whole-blood RNA-seq, and 500 people are done (49,903
   TRB CDR3s). SCEPTR beats ESM-C. The LR∩RNA cohort is 8,327 people, and seven immune disease
   families are powered. Nothing is joined to the HLA calls yet.
5. **The supervisors now frame this as a paper.** Figure 1 = admixture, class I ternary, novel
   rate by ancestry × gene, and SFS novel vs known.

---

## 1. A finding from today that changes how Figure 1 must be framed

`scripts/hla_popgen/23_fig1_draft_panels.py` (new, laptop-only, 6 fixture tests) splits novelty
into three buckets. The split matters more than expected.

| gene | novel haplotypes: CDS-changing | non-coding only | flagged artifact | non-coding share of clean |
|---|---|---|---|---|
| A | 57 | 743 | 657 | 92.9% |
| B | 91 | 2,333 | 729 | 96.2% |
| C | 80 | 1,707 | 826 | 95.5% |
| DPA1 | 40 | 2,698 | 407 | 98.5% |
| DPB1 | 45 | 4,487 | 627 | 99.0% |
| DQA1 | 55 | 3,509 | 395 | 98.5% |
| DQB1 | 50 | 2,214 | 603 | 97.8% |
| DRB1 | 39 | 14,268 | 577 | **99.7%** |

What this means:

- **The "~3,000 novel alleles" headline is mostly non-coding.** After removing flagged artifacts,
  93–99.7% of novel classical haplotypes carry a CDS that already exists in IPD-IMGT/HLA. They
  differ only in introns/UTRs. IPD-IMGT's full-genomic sequences are much sparser than its exon
  sequences, so a large part of this is *catalogue incompleteness*, not new protein diversity.
  Reviewers will find this in minutes if we don't state it first.
- **Recurrent CDS-changing novel alleles are rare.** In the 8 classical genes, **only 7** are seen
  in ≥2 unrelated people, and every one has <20 carriers. Clean CDS-changing singletons are more
  common (~430), but they are unconfirmed by construction.
- **The ancestry gradient survives the split and lives in the non-coding bucket.** Non-coding
  novelty is lowest in EUR for 6 of 8 genes. Examples: DRB1 50.6 per 100 haplotypes in EUR vs
  71.1 in EAS; DPB1 11.7 vs 26.3. That is still reference bias, but it is bias in *genomic*
  sequence coverage of IMGT. It is not automatically evidence of undiscovered protein alleles.
- **Built-in negative control:** the flagged-artifact rate is flat across ancestries (≈2–3.5 per
  100 haplotypes in every gene). Assembly and annotation artifacts are not ancestry-structured, so
  the gradient in the other buckets is not an artifact of uneven assembly quality. Put this in the
  figure.
- **The SFS panel shows the same thing from another angle.** The CDS-changing curve collapses after
  k=2. Non-coding novel alleles have a *known-like* spectrum that reaches hundreds to >1,000
  haplotypes, which is what common, old alleles look like. The top one is a DPA1\*01:03:01
  genomic variant carried by 1,531 people across all 6 ancestries.

**Recommended framing:** "long reads reveal that a large fraction of common HLA haplotypes carry
full-gene sequences absent from IPD-IMGT/HLA, with coverage gaps largest outside Europe; new
*protein* alleles are rare and mostly private." That is more defensible than "3,000 novel alleles".
It is also more useful to the field, because the fix is a concrete IMGT submission (§4, idea W1).

Draft outputs, all aggregate, with participant counts 1–19 written as `<20`:
`reports/hla_popgen/23_fig1_draft_panels/`, containing `panel_novel_rate_by_gene_ancestry.png`,
`panel_sfs_known_vs_novel.png`, `recurrent_novel_cds_changing.tsv`, `recurrent_novel_noncoding.tsv`
(230 alleles in ≥11 people), and `fig1_draft_panels_report.md`.

---

## 2. Supervisor action items — status and concrete plan

Effort: S = hours, M = 1–2 days, L = a week or more. Where: 💻 laptop from committed aggregates;
🖥️ needs the Workbench VM.

### Item 1 — Phasing QC with relatives · **mostly done** · remaining S · 🖥️

Done: scripts 11, 16 and 17 (827 people, 545 high-sharing pairs, 0 switches, 94.8% byte-identical,
named-allele mismatch 5.65% in the classical genes).

Still missing, and cheap:

- **a. The 5 duplicate/MZ-twin pairs, reported separately.** They are the only true technical
  replicates, so their discordance *is* the end-to-end assembly+calling error rate. For relatives,
  "~95% concordance" is a transmission statistic, not an error rate. Supervisors and reviewers will
  read it as the error rate unless we give them the real one.
- **b. Split parent-child from sibling without IBD0.** Two independent handles:
  - Year of birth from the OMOP `person` table: an age gap of about ≥15 years means parent-child.
  - The MHC-block sharing state, legitimate because there are 0 switches. Each pair shares 0, 1 or
    2 whole haplotypes. Parent-child pairs never share 0. Full sibs share 0/1/2 at ≈25/50/25%. The
    29 "mixed" pairs are probably IBD0 siblings or 2nd-degree pairs. Say so explicitly.
- **c. Stratify concordance** by (i) resolution (2-field, 4-field, raw sequence), (ii) whether the
  gene sits on a contig that spans its neighbours, and (iii) ancestry. That tells us whether
  discordance concentrates in fragmented assemblies or in particular ancestries.
- **Figure (supplement):** per-gene concordance with CIs, MZ-pair points overlaid.

### Item 2 — LD for DPA1–DPB1 and DQA1–DQB1 by ancestry · **not started** · M · 🖥️

- **Input:** Table 2 cis-heterodimers, keyed by contig. DQ pairing succeeds for 97.4–99.9% of
  haplotypes. We have *observed* haplotypes, so there is **no EM phasing step**. That is itself a
  methods point: show the phased result next to an EM estimate from the unphased genotypes of the
  same people. Frequency-level agreement will be good. The difference will sit in rare pairs.
- **Metrics:** per allele pair, D′ and r² with counts. Per locus pair, **Cramér's V (Wn)** and
  **asymmetric LD (ALD; Thomson & Single 2014, Genetics — verify)**. ALD is the HLA-standard
  measure when two loci have very different allele counts (DPB1 ≫ DPA1), and plain r² is
  misleading there.
- **Resolution:** 2-field and protein/P-group. 4-field inflates allele counts and deflates r².
- **Ancestry:** apply the stricter admixture threshold the supervisors asked for (e.g. ≥0.9
  single-ancestry proportion), and report n per stratum. MID (~500 people) will be thin, so
  bootstrap the CIs.
- **Positive control:** DRB1–DQB1 (expected very strong). Expected result: DQ strong,
  **DP markedly weaker** (the DP region has higher recombination — verify the magnitude). Any
  ancestry difference in DP LD is interesting in itself.
- **Figure:** ALD heatmap (locus pair × ancestry) plus top-haplotype frequency bars per ancestry.

### Item 3 — Novel allele rate by ancestry × gene · **draft done 💻** · remaining S · 🖥️

- The draft panel exists (§1). It uses carrier-person counts per cluster as the numerator, a
  slight undercount for homozygotes.
- To do on the VM: an exact haplotype-level version (one groupby over Table 1), the stricter
  admixture filter, and the 57 non-classical genes for the supplement.
- **Decide with supervisors:** show all three buckets stacked (current draft), or CDS-changing and
  non-coding as two separate rows of panels. I recommend two rows. The CDS bucket is invisible at
  the shared scale.

### Item 4 — SFS, novel vs known overlaid · **draft done 💻** · remaining S–M

- The draft is pooled across ancestries, with log₂ bins normalized by bin width.
- Remaining: **per-ancestry SFS projected to a common sample size** (hypergeometric down-projection
  to the smallest group's haplotype count). Raw per-ancestry spectra aren't comparable because N
  ranges from 499 to 4,053.
- **Clarify with supervisors which SFS they mean:**
  - (a) the *allele* frequency spectrum, which is what the draft shows; or
  - (b) a true *site* frequency spectrum (derived-allele counts per nucleotide). We can build (b)
    from the PAF projection behind `21_hla_manhattan.py`, split into CDS vs non-coding and
    PBR vs non-PBR. A site SFS would plug straight into the selection analyses in item 7, e.g.
    Tajima's D per segment.

### Item 5 — High-frequency novel alleles missing from IMGT · **list done 💻** · remaining M

- Lists exist: 7 CDS-changing and 230 non-coding alleles in ≥11 people.
- **Before any callout:**
  1. **Re-check against the latest IPD-IMGT/HLA release.** Immuannot ships a fixed database
     snapshot, and some "novel" genomic sequences may have been added since. The release FASTA is
     public, so this can run locally against the CDS/genomic hashes, with no participant data
     leaving the Workbench.
  2. **Characterize the non-coding differences** for the top ~20: position (intron/UTR), SNV vs
     indel, homopolymer context. Check whether all carriers share the *identical* difference,
     which by construction they do (clusters are keyed on sequence hash).
  3. **RNA-seq support for the CDS-changing ones.** For carriers who are also in the 8,327 RNA-seq
     cohort, count RNA reads supporting the novel codon. Transcribed variant reads give orthogonal,
     cheap validation. See the §3 join.
- **Callout wording:** the common non-coding ones are the story, e.g. "a DPA1\*01:03:01 genomic
  haplotype carried by 1,531 participants (13%) across all six ancestries is absent from
  IPD-IMGT/HLA".

### Item 6 — KIR calls + HLA deletions/duplications · **not started** · chr6 part M, KIR part L · 🖥️

- **chr6 structural variation (cheap, existing Table 1):**
  - gene presence/absence per haplotype
  - DRB3/4/5 content, i.e. DR51/DR52/DR53 haplotype groups (expected near-perfect linkage to DRB1
    lineage, so this doubles as a QC)
  - HLA-Y and HLA-H presence
  - MICA deletion
  - C4A/C4B copy number (C4 is already annotated)
  - `copy_index>1` duplications
- **Guard against fragmentation:** call a deletion only when the flanking genes lie on the *same*
  contig, so absence is observed across continuous sequence, not inferred from a contig gap.
- **Validate against Immuannot-paper frequencies** (Zhou et al. 2024): HLA-Y deletion 71.3% HPRC /
  87.3% CPC; MICA deletion 11.9% in MXL.
- **KIR:** chr19 is outside the production trim window, and `01_extract_rich.py` deliberately
  asserts KIR is absent, so this needs a new extraction path. Immuannot annotates KIR natively.
  Reuse the Tier-1 `paf_region` trim with the LRC region (chr19 ~54.7–55.0 Mb, GRCh38 — verify
  coordinates) and re-run Immuannot on the trimmed contigs.
  - Cost: the chr6 production run was ~$200. The LRC window is smaller, but KIR haplotype
    assembly quality is lower.
  - **Pilot 50 people first**: completeness, gene-content haplotypes (cA01/tB01 frequencies vs
    literature), novel-allele rate.
- **Why KIR is worth the L effort:** HLA-C1/C2 and Bw4 ligands × KIR receptor genotype is a
  classic, underpowered-everywhere interaction (§4, idea R4).

### Item 7 — Literature review: HLA selection/diversity methods · **not started** · M · 💻

One document with a pitfall column, organized by what each method assumes:

| Method family | Examples (verify each) | HLA-specific pitfall |
|---|---|---|
| Frequency-spectrum neutrality | Ewens–Watterson homozygosity (Slatkin exact test); Tajima's D | Admixed strata violate panmixia, so run within strict-ancestry strata. Assembly singletons inflate the rare tail; use tiered QC. |
| Rate-based | dN/dS in PBR vs non-PBR codons (Hughes & Nei 1988) | Gene conversion and recombination break codon-model trees. Allele *names* ≠ sequences, so use CDS sequences. |
| Long-term balancing | β statistics (BetaScan, Siewert & Voight 2017); NCD (Bitarello 2018); trans-species polymorphism | Reference bias in the ancestral-allele call; strong LD hitchhiking across the MHC. |
| Diversity-within-individual | HLA evolutionary divergence (HED; Pierini & Lenz 2018; Chowell 2019) | Depends on the distance metric (Grantham) and the exon set. Novel alleles need the observed CDS, which we have. |
| Richness / discovery | Chao1/2, ACE, Good–Turing, iNEXT | No estimator models balancing selection. Already caveated (NOVEL_LIT.md §4). |
| Local-ancestry selection | MHC ancestry deviation in admixed genomes | Needs local ancestry at chr6p. Check AoU availability for the LR cohort. |

Deliverable: `research/SELECTION_METHODS_LIT.md`, then a short proposal of 2 analyses (not 6).
Suggested pair: the site-SFS Tajima's D by segment × ancestry, and HED by ancestry.

### Item 8 — Draft Figure 1 + supplement dump · **partly assemblable today** · S

| Panel | Source | State |
|---|---|---|
| a Admixture barcode | `06_figures_structure/lr/admixture_barcode.png` | exists |
| b Class I ternary | `10_allele_ancestry_geometry/multi_gene/ternary_gene_grid.png` | exists; crop to A/B/C, triangles for novel |
| c Novel rate by ancestry × gene | `23_fig1_draft_panels/panel_novel_rate_by_gene_ancestry.png` | draft |
| d SFS novel vs known | `23_fig1_draft_panels/panel_sfs_known_vs_novel.png` | draft |

Supplement dump: Manhattan panel (21/22), phasing (16/17), SR vs LR (07), saturation (04/18/19),
B\*27 forest (13), frequency heatmaps (05). All already on main.

**Next action:** a `24_fig1_compose.py` that places the 4 panels on one canvas with shared
ancestry colours and a/b/c/d labels. S, 💻.

### Suggested order for the week

1. **Mon–Tue:**
   - 💻 finalize panels c/d with the split decision
   - 💻 compose Figure 1
   - 💻 the IMGT latest-release check for item 5
   - 🖥️ exact novel rate + MZ-twin concordance
2. **Wed:** 🖥️ LD script (item 2); 🖥️ chr6 structural variation (item 6a).
3. **Thu:** 💻 selection-methods lit doc (item 7); 🖥️ KIR 50-person pilot launch (item 6b).
4. **Fri:** supplement dump; update `context/STATUS.md` (stale since 2026-09-05) and
   `EXPERIMENTS.md`.

---

## 3. HLA × TCR/BCR — the integration with Aleix's workstream

**The join is small, and nobody has done it yet.** 8,327 people have LR+RNA. Immuannot calls exist
for ~12k LR people, so the overlap with calls is probably ~7–8k. Measure it on the VM first (S).
That is the denominator for everything below.

### Near-term (weeks), ordered by risk

- **J1. Sanity replication: HLA-associated TCRβ in bulk RNA-seq depth.**
  - Aleix's plan step 2. Published precedent: Emerson et al. 2017 Nat Genet; DeWitt et al. 2018
    eLife (verify).
  - Realistic expectation: those used ~10⁵ unique TRB per person. TRUST4 on whole-blood RNA-seq
    gives ~10³ TRB, so per-person HLA prediction will be far weaker than AUC 0.95.
  - Test instead at the **population level**: for each HLA allele, which public CDR3 clusters
    (SCEPTR/tcrdist neighbourhoods) are enriched in carriers vs non-carriers, via Fisher/CMH within
    ancestry.
  - If A\*02:01 or DRB1\*15:01 associations reappear, the data can carry HLA signal. If they don't,
    it is a measured depth limit, which is also a result.
- **J2. HLA → V-gene usage and CDR3 composition ("CDR3-QTL").**
  - Ishigaki et al. 2022 Nat Genet (verify): HLA autoimmune-risk amino-acid positions (e.g.
    DRB1 pos 13/71, DQB1 57) shape position-specific CDR3 amino-acid usage. Sharon et al. 2016 Nat
    Genet (verify): MHC variants bias TRBV usage.
  - **Our edge:** multi-ancestry, phased 4-field long-read HLA, and amino-acid-position encoding
    from the *observed* CDS, so novel alleles aren't dropped.
  - V-usage frequencies are depth-robust (Aleix's own note), so this is the most powered HLA×TCR
    test available.
- **J3. HLA diversity → repertoire breadth.**
  - Does higher HLA evolutionary divergence (HED) or class I/II heterozygosity give a more diverse
    TCR repertoire (rarefied richness/clonality)?
  - Covariates: age, sex, depth, RIN, cell composition (CD3E/MS4A1 TPM), and a CMV-exposure proxy
    (VDJdb CMV-specific clonotype hits).
  - Clean, one model, n≈7k. Plausible either way, so both outcomes are publishable as a
    well-powered test.
- **J4. RNA read support for novel HLA alleles** (also feeds item 5). Personalized references from
  each carrier's own assembly, then count RNA reads on the novel codon. It also gives
  **allele-specific HLA expression** from phased assemblies for free: known HLA-C expression
  differences by allele (verify), now measurable per ancestry.
- **J5. HLA-restricted antigen-specific clones.**
  - VDJdb entries carry a restricting HLA. Is a person's hit rate for epitope-specific clonotypes
    higher when they *carry* the restricting allele?
  - It's a biological positive control, and it turns the exposure readout (CMV/EBV) into an
    HLA-aware feature for disease models.
  - Watch for the known ancestry skew of VDJdb.

### The standing confounds

State these on every HLA×repertoire result:

- **Healthcare-access gradient.** EUR records 41 median conditions vs SAS 16, so any
  disease×ancestry contrast mostly measures access. Use within-ancestry models, EHR-years
  normalization, and negative-control phenotypes.
- **Depth and chain recovery by ancestry.** Kruskal–Wallis p=0.18 at n=100; re-test at n=500.
  Rarefy.
- **Population structure.** HLA alleles are ancestry-differentiated, so HLA↔repertoire
  associations can be ancestry proxies. Adjust for genome-wide PCs *and* stratify.

---

## 4. Future lines of research

### Rational, high-probability

- **R1. Graph-genotype the full AoU short-read cohort with our pan-MHC.**
  - Build a pangenome/allele panel from ~24k phased long-read haplotypes (HPRC's pan-MHC used ~246
    haplotigs — verify).
  - Genotype the ~400k+ short-read WGS with a haplotype-aware locus genotyper (e.g. Locityper —
    verify tool and fit).
  - Use the ~12k LR∩SR overlap as truth, stratified by ancestry.
  - This closes the "AoU-native vs long-read" story from `DECISIONS.md` and gives LR-quality HLA,
    novel alleles included, to the whole biobank. It is the most *impactful* extension, because
    every AoU HLA study would use it.
- **R2. Multi-ancestry HLA PheWAS at amino-acid resolution.**
  - Phecode-based (not substring) phenotypes. Positions from the observed CDS. Conditional analysis
    to separate independent signals (the fine-mapping aim in `TASK_CONTEXT.md`).
  - Once R1 exists, run at ~400k rather than 12k, which is where the original non-linear PRS aim
    becomes powered.
- **R3. Structural haplotypes of the MHC by ancestry.**
  - DR haplotype groups, C4 copy number and C4A/C4B, HLA-Y/H/MICA deletions, DRB paralog content.
  - C4 copy number → SLE/schizophrenia literature (Sekar 2016; Kamitaki 2020 — verify) is a
    natural phenotype link.
- **R4. KIR–HLA ligand compatibility** (after item 6). KIR3DL1×Bw4 and KIR2DL/S×C1/C2 genotype
  combinations vs infection/autoimmune families in the powered disease list. Very few cohorts have
  sequence-resolved KIR *and* HLA at this scale, in multiple ancestries.

### Innovative, medium risk

- **I1. Immune germline atlas beyond HLA: IGH and TRB loci from the same assemblies.**
  - The IG/TR germline loci (IGH chr14q32; TRB chr7q34; IGK/IGL; TRA) are structurally complex and
    poorly catalogued outside Europeans, the same story as HLA.
  - Tooling exists for assembly-based IG/TR germline annotation (IGenotyper; digger/OGRDB-style
    annotation — verify current best).
  - **Then join to Aleix's expressed repertoire:** does a person's germline TRBV/IGHV allele
    content or deletion predict their expressed V usage? Which novel germline alleles are actually
    used?
  - One cohort gives genotype (long-read assembly), expressed repertoire (RNA) and phenotype (EHR)
    for the same ~8k people. That is the unifying "Omni-immune" paper and the natural merge of
    Marc's and Aleix's tracks.
- **I2. Functional diversity, not sequence diversity.**
  - Predict each allele's peptide-binding repertoire (NetMHCpan/NetMHCIIpan — verify licence for
    novel sequences). Compute per-person functional breadth (fraction of a pathogen proteome
    presentable).
  - Test whether the ancestry diversity gradient holds functionally, and whether functional breadth
    predicts TCR repertoire breadth (J3) and infection burden.
- **I3. Where the 7 recurrent CDS-changing novel alleles land.** Map their changes onto groove
  structures (PBR pocket residues), predict binding shifts, and check RNA expression (J4). A small,
  concrete vignette for a main-text panel.
- **I4. Repertoire-inferred HLA as a sample-identity and call-QC tool.** Predict HLA from TCR (J1)
  and flag people whose long-read call disagrees with their repertoire. Useful for detecting sample
  swaps between the LR and RNA assays at biobank scale (a real, underappreciated multi-omics
  problem).

### Out-of-the-box, high risk

- **W1. "The IMGT that AoU would build."**
  - Submit the recurrent full-genomic sequences (the 230 non-coding novel alleles in ≥11 people) to
    IPD-IMGT/HLA as a community resource.
  - Use the discovery-rate curves (18/19) to model *how many more people, from which ancestries*,
    it takes to reach X% of common-haplotype genomic coverage: a design argument for future
    reference panels.
  - Needs an AoU data-policy check before any sequence submission (DECISIONS.md egress question).
- **W2. Post-admixture selection at the MHC.**
  - African-American and Latino genomes show local-ancestry deviations at 6p (verify current
    evidence). With allele-level phased haplotypes, ask *which* HLA alleles ride the over-represented
    ancestry, turning a blurry admixture-mapping signal into an allele-level hypothesis.
  - Requires local ancestry on the LR cohort; check availability.
- **W3. Mosaic chr6p loss in blood.**
  - Copy-neutral LOH of 6p (HLA loss) is a known immune-escape mechanism in aplastic anaemia and
    post-transplant relapse (verify). In ~12k blood long-read samples, allelic imbalance between
    phased HLA haplotypes (read support per haplotype) could detect clonal 6p LOH.
  - Associate with age, CHIP-like phenotypes, infection and cancer. Speculative, but cheap to
    screen once per-haplotype read depth is extracted.
- **W4. Three-way co-evolution: HLA × KIR × TCR within individuals.** Do people whose HLA presents a
  narrower functional peptide space (I2) compensate with NK-receptor diversity (KIR) or repertoire
  breadth? A systems-level "immune portfolio" hypothesis. It only becomes testable once items 6
  and J3 exist.
- **W5. Allele age from haplotype background.** Novel alleles on long, unbroken extended haplotypes
  (shared across unrelated carriers) are young; those on diverse backgrounds are old or recurrent.
  In a 24k-haplotype multi-ancestry sample this is an allele-dating tool for HLA. It also refines
  the novel-allele QC: an "allele" found on many unrelated backgrounds with an identical rare
  change suggests a systematic artifact.

---

## 5. A possible paper architecture

**Paper 1 — long-read HLA atlas (Marc; supervisors' current focus)**

1. Cohort, reference gaps, novel alleles (three-bucket framing), SFS
2. Validation: MZ twins, relatives, SR vs LR, RNA support
3. Phased structure: LD and haplotypes by ancestry, ternaries
4. Selection: CDS vs non-coding π, PBR, site SFS, HED
5. Structural variation: gene content, C4, (KIR)
6. Translational hook: B\*27/AS replication, and/or R1 short-read graph genotyping

**Paper 2 — HLA × repertoire (Aleix + Marc):** J1 → J2 → J3 as the core, J4/J5 as validations, I1
as the ambitious extension.

---

## 6. Things to flag now

- **Disclosure policy on the public repo.** `reports/hla_popgen/novel_alleles.tsv` (committed)
  lists per-cluster `n_persons` and per-ancestry `ancestry_counts` down to 1. Several other
  aggregate TSVs carry small cells too. The AoU Data and Statistics Dissemination Policy forbids
  publishing participant counts of 1–19 (verify the current wording). The new `23_*` outputs
  suppress these. **The existing committed files were not changed**, because rewriting published
  history is a decision for Marc and the sponsor, not something to do quietly.
- `context/STATUS.md` is dated 2026-09-05 and doesn't mention scripts 10–22, the merge, or the
  supervisor meeting.
- Local branch `needle-view-cds-diversity-density` is fully merged into main. This work lives on
  `fig1-drafts-and-research-map`, which is not pushed.
