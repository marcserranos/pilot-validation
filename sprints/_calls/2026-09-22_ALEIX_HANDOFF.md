# Handoff for Aleix — supervisor call of 2026-09-22 (Progress update #8) and project state

*Written for Aleix to read himself or to give to his agents. Self-contained: no need to read the
call transcript. Aleix was not on the call. Source material: the call transcript, the slide deck
(`REPORT #8`), Cole's Slack messages that evening, and the repo. The fuller Marc-side writeup is
[`CALL_SUMMARY_2026-09-22.md`](CALL_SUMMARY_2026-09-22.md); this file is the Aleix-relevant
distillation plus the project context he needs around it.*

**Confidence flags used below:** **[said]** = Cole/David said this on the call or in Slack.
**[repo]** = read from the repo. **[inferred]** = my inference, check it. **[stale?]** = comes from
Marc's notes on Aleix's branch as of 2026-09-12, which I have not re-verified; Aleix knows his own
current state better than this file does.

---

## 0. The short version

1. The supervisors (Cole Shanks, David Bonet) are treating this as **one paper aimed at a big
   journal**, with a Figure 1 already drafted. They are worried about being scooped and want
   **focus, not ten side-projects**. **[said]**
2. The project is now formally **two halves**: **stat-gen / population genetics** (Cole leads it
   himself) and **prediction / ML** (Marc + Aleix). Aleix's half is the **immune-repertoire
   (TCR/BCR) → disease prediction** work, with HLA as a secondary, later-added feature. **[said]**
3. The instruction for the prediction half is explicit: **start simple, repertoire features only,
   no HLA at first, baseline regression then XGBoost, copy the methodology of one specific
   benchmark paper (BenchRep-T), and care about interpretability more than accuracy.** **[said]**
4. **Phenotypes: autoimmune and chronic conditions (HIV, B-cell disorders). Avoid cancer.** **[said]**
5. Nothing is joined between Marc's HLA calls and Aleix's repertoire data yet. The join is the
   next real integration step. **[repo]**

---

## 1. Who is who

| Person | Role in this project |
|---|---|
| **Cole Shanks** | Lead supervisor. Drives the science and the statistical-genetics half. Available on Slack at any hour and explicitly invites questions ("don't hesitate to bother me... even if it's a weird time"). **[said]** |
| **David Bonet** | Supervisor, described by Cole as "one of the world's leading experts in tabular models." Joined the call for the prediction discussion. Will give modeling advice. **[said]** |
| **Marc Serrano** | Owns the HLA long-read calling and everything HLA/population-genetic that was presented. Also on the prediction half. |
| **Aleix Ruf** | Owns RNA-seq immune-repertoire extraction (TRUST4), repertoire embeddings, and the disease-cohort definition. Prediction half. |

Aleix handed HLA calling to Marc on 2026-08-03 and moved to the RNA-seq repertoire work. **[repo]**
His work lives on branch `origin/aleix/hla-resolve-phase1` under `aleix/RNA-seq/`, not on `main`.

---

## 2. What the supervisors want from Aleix's half (the actionable part)

### 2.1 Scope: prediction, and why

Cole: "absolutely we should do prediction." Two reasons he gave: the team is good at it, and
**prediction is how you show there is signal at all** ("prediction tells you that there's something
going on"). The stat-gen half then explains it. He said the ML side is "the fun machine learning
stuff" and assigned it to Marc and Aleix, including "all the cool stuff with the embedding models."
He is taking the stat-gen half himself because "there's a lot there to screw up."

His view on where signal is: **HLA alone is low-signal. The repertoire is where the signal is, and
HLA structures the repertoire.** **[said]**

### 2.2 The order of work

1. **Repertoire features only, no HLA.** Cole's reason: HLA and TCR features are "so correlated to
   each other" that adding HLA early "doesn't add a ton of signal," only dimensionality, which
   "makes learning harder." Add HLA later, and possibly only for particular cases, because "we're
   not going to be observing the entirety of the immune repertoire, so maybe having the HLA is
   really good." **[said]**
   - **[inferred]** This is consistent with Aleix's own plan (his `POST_TRUST4_OPTIONS.md`, Tier 1
     step: "validate that the repertoire predicts HLA type, published AUC ~0.95"). If the
     repertoire predicts HLA that well, HLA is largely redundant as a *feature* and far more useful
     as a *conditioning variable or stratifier*. That is a framing worth stating to Cole.
2. **Baselines before anything fancy.** Regression on V/J-gene usage plus CDR3 short-motif
   **k-mer** features, then XGBoost. Cole: "You have to do better than a simple base model. Don't
   get cute." He would start with V/J + k-mers, then logistic regression, then XGBoost.
   - *Transcript note:* the raw transcript wrote "VJ CAMERS." That is a transcription error for
     **k-mers**. See §7.
3. **SCEPTR is allowed but is not the starting point.** "You can do the SCEPTR, because I think
   that might be really good. But definitely regression over the k-mers is the baseline for sure."
   Aleix's own SCEPTR-vs-ESM-C result is already in hand **[stale?]** and stays valid as an
   embedding option.
4. **Only then** consider fancier tabular models. David is the person to ask for that.

### 2.3 The reference paper to copy: BenchRep-T

**Im C, Cohen-Lavi L, Buendia A, Kundaje A, Boyd SD. "BenchRep-T: A Systematic Evaluation of
T-Cell Repertoire-Based Disease Diagnostics."** bioRxiv, 2026.
`https://www.biorxiv.org/content/10.64898/2026.06.09.727013v1` (Cole sent it at 8:04 PM).

- Benchmarks nine methods, statistical through deep learning, on disease classification from blood
  TCR sequences, plus limited-data performance and antigen-specific sequence detection.
- Headline finding: **"simple baselines prove competitive, with tree-based models trained on V- and
  J-gene usage and short sequence motifs approaching the classification performance of more complex
  methods."** This is published backing for Cole's "start simple."
- Cole says it has "really nice technical appendices, telling you how they're doing absolutely
  everything," a GitHub repo, and dataset availability. **[said]**
- **What Cole asked for specifically:** copy their methodology, "exactly the same but on All of Us."
  He stressed the *evaluation design*, not just the features: how they split cases and controls,
  whether they use **outer splits and inner splits for cross-validation, with the inner tuning
  separated from the outer test**, and "all these questions... that is important here." **[said]**
- Cole's caveat about our data vs theirs: their dataset is smaller but uses a better assay; ours is
  bigger but "we're also only going to get the high clonal counts... it's noisy." **[said]** That is
  a real difference to document: TRUST4 on bulk whole-blood RNA-seq recovers only the abundant
  clones.

**Concrete tasks:** read it and the appendix; extract their split protocol, features, baselines,
metrics, and hyperparameter search; reproduce that protocol on the AoU repertoire data; list every
place our data forces a deviation.

### 2.4 Interpretability over accuracy

Cole was emphatic. The point of a good predictor is that "if something is a good predictor, then
there's probably really something interesting to look at" with conventional statistical approaches.
If a feature is strong, do model interpretation, and if it interacts with another feature, that
becomes "a really concrete question that might actually show some pretty new biology." He warns:
"we're not just trying to absolutely maximize the accuracy, which might be pretty misleading...
we don't want to do anything silly." **[said]** So: report feature importances and interactions, and
do not chase a leaderboard number.

### 2.5 Phenotypes

- **Do not use cancers.** Cole thinks cleaning them is "maybe an impossible task": how many are in
  remission, how differently they are treated, patients on "crazy immunotherapies that we might not
  even be able to pick out from the EHR." **[said]**
- **Use autoimmune and chronic conditions.** He named **HIV** ("HIV is really good") and **any
  B-cell-related disorders** ("very good"), plus autoimmune disease generally. **[said]**
- He flagged that phenotype cleaning is the thing to be "super careful about." **[said]**
- **[repo/stale?]** Marc's notes on Aleix's cohort work: 8,327 people have both long reads and
  RNA-seq; 7,726 have a non-zero EHR window; seven of ten immune-disease families are powered
  (allergy 2,297; rheumatologic autoimmune 444; endocrine autoimmune 327; GI autoimmune only 183);
  there is a strong EHR healthcare-access gradient by ancestry.
- **[inferred] Things to reconcile with Cole's guidance:** (a) check whether any of the ten families
  are cancer-derived and drop or demote them; (b) check whether **HIV** and **B-cell disorders**
  are among the ten families or need adding, and how many people each has in the 8,327; (c) the
  ancestry gradient in EHR coverage is a confounder that Cole will care about.

### 2.5b Confounders and traps worth stating up front **[inferred]**

Not said on the call but implied by it and by Aleix's own pilot: sequencing depth and RNA quality
vary by person (Aleix's pilot showed mean CDR3 counts per person from 7,406 in AFR to 4,538 in EAS,
partly depth, Kruskal-Wallis p = 0.18 so not significant at n=100); ancestry structures both HLA and
repertoire; EHR follow-up length differs by ancestry; relatives in the cohort break independence,
so **splits must be grouped by kinship** (Marc's pipeline already has a relatedness table from the
AoU KING-style output). Ask Marc for the relatedness list before building splits.

### 2.6 The ideas that were floated but not assigned

- **Peptide-groove ↔ repertoire link (Cole, exploratory).** The peptide-binding groove of HLA is
  what presents antigen, so "it might be really interesting to connect to the TCRs." Picture: take
  ~100 TCR motifs associated with DRB1\*05 and another set with DRB1\*04, and ask whether
  differences in the groove correspond to distance in repertoire features. Marc floated an
  embedding of the groove mapped to an embedding of the TCRs/BCRs. Cole guessed **SCEPTR
  "probably mostly picks up on the peptide groove, the way it's trained,"** and suggested checking
  the SCEPTR paper for whether they discuss it. He said "I'll think about this." **No task
  assigned.** But Aleix is the SCEPTR person, so reading the SCEPTR paper with this question in mind
  is a cheap, useful thing to do.
- **AlphaMissense-style collapsing of rare protein-altering variants** into one burden feature.
  Parked. Not for this stage. **[said]**

---

## 3. What Marc is providing to Aleix's half

- **HLA calls:** common **two-field** HLA alleles only, already pushed to the bucket (Cole confirmed
  he used them). Cole's rule: "the novel ones are interesting, but they're kind of dangerous, and we
  can't say much with them because we don't have a lot of them." So **do not use novel alleles as
  features.** **[said]**
- **KIR calls** will be added after Marc reruns the pipeline on the KIR region (see §5.2).
- **Ancestry labels/proportions** and **relatedness** are computed on Marc's side and can be shared.
- **A per-individual phasing-confidence flag** (same-contig or not, for DQ/DP pairs) is being built
  because Cole asked for it. It matters for anything haplotype-level.
- **Coordination:** Marc said on the call that he would "talk with Aleix" because he is "not really
  familiar with the shape of the TCR and BCR data," and that he would start thinking about the
  prediction side. Expect Marc to ask Aleix for the repertoire data layout.

---

## 4. State of the whole project, so Aleix has the full picture

### 4.1 What the project is

Direct HLA calling from **long-read assemblies** (Immuannot, on ~12,000 AoU long-read participants;
65 genes, 4-field resolution, phased and contig-aware), used to characterize HLA diversity across
ancestries, and then to connect HLA to the immune repertoire and disease. **[repo]**

### 4.2 Headline results so far (Marc's side)

**Corrected novelty numbers (an earlier headline was wrong).** The old "~3,000 novel alleles" was
wrong and is superseded. The corrected picture **[repo]**:

- **1,404 distinct novel protein alleles**; **231** seen in ≥2 unrelated people; **29 (slide 4) to 38
  (callout list including synonymous)** reach ≥20 unrelated carriers.
- **None of the common novel alleles is in the classical genes.** They are **TAP1, TAP2, MICB, MICA,
  HLA-G, DQB2, and others**. The classical HLA-A/B/C/DR/DQ/DP genes are essentially well catalogued
  at the protein level. The 12 named "≥20 carriers" alleles are overwhelmingly **African-ancestry**.
- Most novelty in the classical genes is **non-coding** (93–99.7% of novel haplotypes carry a CDS
  already in IPD-IMGT/HLA and differ in introns/UTRs). That is largely *catalogue incompleteness*
  (IMGT has far fewer full-genomic than exon sequences), not new protein diversity.
- The reference-incompleteness gradient by ancestry is real, and an artifact-rate negative control
  stays flat (~1–3%) across ancestries, so it is not uneven assembly quality.

**Quality control.** Technical duplicates/twins (5 pairs) differ at **121 of 306,033 bases ≈ 1 error
per 2,500 bases** (Q34-ish). Phasing: relatives show **zero phase switches** across 3,021 testable
transitions with 100% power on injected synthetic switches. Physical phasing yield (both genes on
one contig): **DPA1~DPB1 96.96%, DQA1~DQB1 96.31%, DRB1~DQB1 91.72%, B~C 86.88%, A~B 30.24%.**

**Linkage disequilibrium.** Class II linkage is weaker in African-ancestry haplotypes than
European at equal sample size (DQA1~DQB1 D′ 0.908 vs 0.974). DPA1~DPB1 does not follow that pattern.
LD is **bimodal** (near-perfect or near-zero). Cole explained why: the **DQ G1/G2 rule** (§4.4).

**Deletions / structural variation.** A gene is called deleted only from a "bridged absence" on one
contig. Positive control passes: DRB3/4/5 presence matches the textbook DR51/52/53 expectation
(92–98%). Rates: DRB5 ~83%, DRB4 ~70%, DRB3 ~49%, C4B ~20%, C4A ~11%, MICA ~4%. Deletions of HLA-A,
DQB1, DRB1 appear at low rates and Cole treats **HLA-A deletion as a "big result if real"** that
needs short-read validation before anyone believes it.

**Selection.** Amino-acid diversity is concentrated in the peptide-groove exons (2.2×–8.0× enriched
in eight genes) with the DRA and HLA-F controls showing nothing. Diversity ranking recovers known
functional residues (DRB1 β11/13/71/74, DQB1 β57, HLA-B 77/80) without being told about disease.

### 4.3 The negatives, stated honestly **[repo]**

Unsupervised HLA embeddings did not separate disease. A supervised UMAP failed held-out testing
(AUROC 0.48). Chao2 saturation estimates are low-confidence. Cole's reaction on the call is
consistent with this: HLA alone is low-signal.

### 4.4 The DQ G1/G2 rule (Cole's main biological insight of the call)

Source: Petersdorf EW et al., *Blood* 2022;139(20):3009–3017, "HLA-DQ heterodimers in hematopoietic
cell transplantation" (PMC9121842). **Group 1 (G1):** DQA1\*02/03/04/05/06 paired with DQB1\*02/03/04.
**Group 2 (G2):** DQA1\*01 paired with DQB1\*05/06. A G1 alpha with a G2 beta (or the reverse) on the
same haplotype forms a non-functional heterodimer, and background selection removes it, which is why
LD is bimodal. Clinically, G1G2/G2G2 genotypes carry higher post-transplant relapse risk. DP has no
equivalent rule. Cole made a target figure (`reference/cole_dq_g1g2_target_figure.pdf`: pairwise
**signed phased D′**, DQA1 × DQB1, split into a 2×2 grid at the G1/G2 boundary, cross-group
quadrants labeled "Predicted incompatible") and asked Marc to recreate it. He also asked whether the
lighter-blue cells inside the "incompatible" quadrants are **phase errors or rare real
recombination**, to be tested against per-individual phasing confidence.

**Why this matters to Aleix:** it is the clearest example of HLA structure the repertoire work can be
conditioned on. Functional DQ heterodimer group is a defined, biologically meaningful covariate.

---

## 5. Things on the call that affect Aleix indirectly

### 5.1 Figure 1 (the paper's introductory figure)

Four panels: **(a)** admixture barplot of the long-read cohort (n = 11,833 unrelated); **(b)** HLA-B
ternary plot of allele ancestry composition; **(c)** reference-catalogue incompleteness by gene and
ancestry; **(d)** novel proteins per gene by recurrence class. Feedback: add 0–100 ticks to the
ternary, filter it to stricter ancestry (95–98%) because admixed people drag the cloud to the center,
drop the "≥98%" text on panel (a), and probably merge or replace (c)/(d), with an ancestry × gene
novelty-rate version preferred over (c). Per-ancestry **saturation/discovery curves** were also
requested, inspired by the **Pakistan Genome Resource** paper (*Nature* 2026, 173,303 exomes and
genomes), where Cole noted the curve fails to saturate in a strongly endogamous, structured
population. He expects African ancestry to be the least saturated in ours. This is Marc's work.

### 5.2 KIR genes are being added

Cole: the pipeline only calls the HLA region, so **KIR genes were never called**, and he thinks they
are "going to be very important," because KIR receptors on NK cells "have really particular
interactions with the HLA," and he expects more germline diversity and novel alleles there than in
HLA. He wants the long-read pipeline rerun on the **KIR region** (chr19), possibly with a panel in
Figure 1. Cost is not a concern to him. The KIR region is small (~150 kb) but structurally complex
(duplications, deletions). Marc will scope it and report next week. **For Aleix:** KIR calls will
become another HLA-side feature, and KIR–HLA ligand pairs are a classic interaction term for any
later interpretability work.

### 5.3 Deletions and copy number

Cole linked the DRB3/4/5 variability to a recent argument that selection favors **higher HLA copy
number** ("from an antigen-presenting point of view there's really no downside to another gene").
He wants DRB3/4/5 copy number by ancestry as a panel. Only tangential to Aleix, but it means HLA
copy-number variants exist as a possible covariate.

---

## 6. Decisions the supervisors made or implied on this call **[said]**

1. Short reads are for **validating** long-read calls, not discovering alleles. Stay in this cohort.
2. **KIR gets added.**
3. The project splits into **prediction (Marc + Aleix)** and **stat-gen (Cole)**, joined by Figure 1.
4. Prediction starts from **repertoire features alone**, common two-field HLA only, baselines first.
5. Phenotypes: **autoimmune/chronic (HIV, B-cell disorders), not cancer.**
6. Reject scope creep: "we shouldn't shove in 10 projects into one." Things can be spun out later.
7. Cole is worried about **being scooped** by another group doing this and wants a big journal.
   Consequence: **be careful and correct rather than ad hoc.** "It's very easy to do an ad hoc analysis
   and be like, 'let's throw it in.' It's usually not worth it."

---

## 7. Glossary and transcript corrections (so agents do not get misled)

The call was auto-transcribed and several terms are wrong in the raw text. Corrected here:

| Raw transcript | Means |
|---|---|
| "current genes" / "curve genes" | **KIR** genes (killer-cell immunoglobulin-like receptors) |
| "chromosome 14" | KIR is on **chromosome 19** in every repo reference |
| "VJ CAMERS" / "CAMERS" | **V/J-gene usage + CDR3 k-mers** |
| "the onshore paper" | **BenchRep-T** (Im et al., bioRxiv 2026) |
| "deep prime" (LD) | **D′** (D-prime), the LD statistic |
| "king sheep related score" | the **KING** relatedness coefficient (φ) used to find duplicates/relatives |
| "IgT" / "IMGT" | **IPD-IMGT/HLA**, the reference allele database |
| "scepter" | **SCEPTR**, the TCR embedding model Aleix already benchmarked |
| "cis" vs "trans" | on the same haplotype vs on different haplotypes |
| "G1 / G2" | the DQ heterodimer groups of §4.4 |

Other terms: **CDR3**, the hypervariable receptor region; **MIL**, multiple-instance learning
(Aleix's own option list); **AoU**, All of Us; **LR / SR**, long-read / short-read; **ARS**, antigen
recognition (peptide-binding) site; **D′ / r²**, LD statistics.

---

## 8. Constraints Aleix's agents must respect

- **The repo is public.** Docs and pipeline code only. No participant-level data, no controlled-tier
  identifiers, and **no small counts.** Counts of 1–19 are written as `<20` in anything committed.
  A past bug silently turned suppressed `<20` cells into zeros and produced a falsely clean figure,
  so **never coerce a suppressed count to 0.** **[repo]**
- **Do not commit** anything from the controlled workbench without checking with Marc. Whether
  small-cell counts may be published at all is an open item in `context/DECISIONS.md`.
- **Do not touch Marc's branches or VM concurrently** without coordinating. Concurrent runs against
  the same output paths have corrupted results before. Marc and Aleix work in separate Workbench
  app instances and on separate git branches, so a resource that exists on one side may not exist on
  the other; check before redoing it.
- The repo's working convention: number scripts `NN_description.py` with a matching
  `reports/hla_popgen/NN_description/`, and update `context/EXPERIMENTS.md` (append-only) and
  `context/STATUS.md` (rewritten each session).

---

## 9. Suggested work order for Aleix's agents **[inferred from the above; adjust freely]**

1. **Read BenchRep-T and its appendix.** Write a one-page protocol note: features, splits (outer and
   inner), baselines, metrics, tuning, and every deviation our data forces. *Deliverable before any
   modeling.*
2. **Reconcile the cohort with Cole's phenotype guidance** (§2.5): drop cancer families, check HIV and
   B-cell disorders, report n per phenotype per ancestry, and flag ancestry-linked EHR coverage bias.
3. **Build the feature table** from the TRUST4 output for the 8,327 people: V/J-gene usage, CDR3
   k-mer counts, simple diversity summaries. Keep per-person depth and read counts as covariates.
4. **Get the kinship list from Marc** and build **grouped** outer/inner splits (relatives never
   straddle train and test).
5. **Baselines:** logistic regression, then XGBoost, following BenchRep-T's protocol. Report against
   a trivial baseline (covariates only: age, sex, ancestry, depth).
6. **Interpretation:** feature importances, and interaction checks, especially with ancestry and HLA.
7. **Only then** add SCEPTR embeddings, then HLA (common two-field alleles, DQ G1/G2 group, later KIR)
   as extra features or as a stratifier, and test whether they add anything over the repertoire alone.
8. **Optional, low cost:** read the SCEPTR paper for whether it captures peptide-groove signal.

---

## 10. Open questions

- Which of the ten disease families are usable under Cole's rule, and is there enough HIV and
  B-cell-disorder signal in the 8,327? **[open]**
- How much of the AoU repertoire is high-clonal-count only, and how does that limit power vs
  BenchRep-T's data? **[open]**
- Does adding HLA add anything over the repertoire once repertoire features are in? Cole expects
  "not much, at first." **[open, his hypothesis]**
- Whether the peptide-groove ↔ repertoire idea is testable in this cohort. **[open, no owner]**
- Whether short-read validation of the 38 named novel alleles restarts. Marc wants it; Cole said to
  focus on validating long reads first. **[open]**

---

## 11. Who to ask, and how

Cole: Slack, any hour, invites questions and thoughts. David: modeling and tabular-method advice.
Marc: HLA calls, ancestry, relatedness, KIR timing, the repo and VM.

*Provenance: some statements about Aleix's own data (cohort counts, depth results) come from Marc's
notes dated 2026-09-12 and were not re-verified against Aleix's branch. Aleix should correct them.*
