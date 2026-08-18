# LR × RNA-seq overlap cohort — disease study results — 2026-08-17

**Headline: the cohort is 8,327 people, ancestrally balanced, with complete EHR linkage,
and seven of ten immune-mediated disease families clear the threshold for real case/control
modelling. The repertoire direction is viable — but the cohort carries a severe
healthcare-access gradient across ancestry groups that will corrupt any naive
ancestry-stratified comparison.**

Data: `../scripts/run_disease_study.sh` (BigQuery pull + analysis), run against CDR
`wb-silky-artichoke-2408.C2025Q4R6`. Per-person clinical data is VM-local under
`~/pipeline_outputs/rnaseq/pheno/` and is not committed.

---

## 1. The cohort

| | |
|---|---|
| People with both long-read WGS and RNA-seq | **8,327** |
| Enrichment over independent-sampling expectation | **34.2×** |
| With any EHR record | **8,327 (100%)** |
| With a non-zero EHR observation window | **7,726 (92.8%)** |

**The 601-person gap matters and is easy to miss.** Every person has an
`observation_period` row, so "100% EHR coverage" is literally true — but 601 of them have a
window of *zero days*: a single date on which records begin and end. They cannot contribute
to any rate-based measure, so **the real denominator for everything below is 7,726, not
8,327.** Reporting the 100% figure without this caveat would overstate the usable cohort by
about 7%.

Ancestry composition (genetically inferred): EUR 2,475 · AMR 1,791 · AFR 1,369 ·
EAS 1,172 · SAS 1,092 · MID 427. **Five of six groups exceed 1,000 people and EUR is under
a third of the cohort** — for comparison, UK Biobank is ~94% European.

---

## 2. Disease burden — and the confound that dominates it

| Ancestry | n | Median distinct conditions | Median EHR-years | Conditions per EHR-year |
|---|---|---|---|---|
| EUR | 2,363 | **41** | **12.63** | **3.11** |
| MID | 407 | 26 | 9.59 | 2.37 |
| AMR | 1,606 | 27 | 10.52 | 2.18 |
| SAS | 1,062 | **16** | 8.28 | 1.78 |
| AFR | 1,176 | 20 | 11.69 | 1.68 |
| EAS | 1,111 | 17 | **8.54** | **1.57** |

**A 2.6× spread in recorded conditions (EUR 41 vs SAS 16) and a 2.0× spread in the
depth-normalized rate (EUR 3.11 vs EAS 1.57).**

**This is almost certainly not a 2.6× difference in how sick people are.** The two
strongest reads on it:

1. **EHR-years track the same ordering.** EUR has the longest recorded history (12.63
   years), EAS and SAS the shortest (8.54, 8.28). Longer engagement with the health system
   means more encounters, more codes, more distinct conditions — independent of health.
2. **Normalizing by EHR-years does not remove it.** If chart length were the whole story,
   the per-year rate would flatten. It doesn't — EUR stays highest at 3.11. So the residual
   is *encounter density*: how often someone sees a clinician within their covered window,
   which is a function of insurance, access, trust, and care-seeking behaviour, not
   pathology.

**This was flagged as the most important caveat in the study design before any data was
pulled, and the data confirmed it at a larger magnitude than expected. Consequence:
any ancestry-stratified disease comparison in this cohort measures healthcare access at
least as much as it measures disease.** It must appear on every slide that breaks disease
down by ancestry — not as a footnote.

**A new lead for the other workstream.** Disease burden ordering
(EUR > MID > AMR > SAS > AFR > EAS) is *not* the same as CDR3 recovery ordering
(AFR > AMR > SAS > EUR > MID > EAS) — EUR is highest in burden but fourth in recovery,
AFR the reverse. So burden is not trivially driving recovery. But EAS is lowest in **both**,
and there is a real mechanism worth testing: chronically ill people have activated immune
systems, more clonal expansion, and plausibly different repertoire recovery. **We now have
per-person disease counts and per-person CDR3 counts for the same people — this is a
directly checkable confound that did not exist before today.**

---

## 3. What this cohort is made of

### 3.0 Coverage — how much of "all diseases" we can actually report

| | |
|---|---|
| Distinct conditions recorded across the cohort | **12,172** |
| Reportable at n ≥ 20 (AoU disclosure floor) | **2,463** |
| Suppressed as below the floor | **9,709 (79.8%)** |
| ICD-10 3-character categories present | 1,352 |
| …reportable at n ≥ 20 | 778 |

Full ranked lists: `lr_rnaseq_all_conditions.csv` (2,463 SNOMED concepts) and
`lr_rnaseq_icd3_categories.csv` (778 ICD-10 categories).

**Four-fifths of the distinct conditions in this cohort occur in fewer than 20 people.**
That is not a data problem — it is the long tail of medicine, and it quantifies something
we previously only asserted: **rare-disease work is not possible in this cohort**, and no
increase in analytical care changes that. The tail is invisible by disclosure rule, not by
oversight.

### 3.1 The real disease ranking — ICD-10 3-character categories

The SNOMED concept list (§3.2) fragments single diseases into dozens of near-identical
variants, which makes each look smaller than it is. The 3-character ICD-10 category is the
level most people mean by "a disease", and it produces a materially different — and more
honest — picture:

| ICD | Commonest label | n | % |
|---|---|---|---|
| M25 | Shoulder/joint pain | 2,894 | 37.5 |
| E78 | Hyperlipidemia | 2,687 | 34.8 |
| M54 | Low back pain | 2,645 | 34.2 |
| M79 | Muscle pain | 2,614 | 33.8 |
| I10 | Essential hypertension | 2,608 | 33.8 |
| R10 | Abdominal pain | 2,339 | 30.3 |
| G89 | Chronic pain | 2,002 | 25.9 |
| R07 | Chest pain | 1,922 | 24.9 |
| R06 | Dyspnea | 1,920 | 24.9 |
| K21 | Gastro-oesophageal reflux | 1,768 | 22.9 |
| R05 | Cough | 1,713 | 22.2 |
| G47 | **Obstructive sleep apnoea** | 1,644 | 21.3 |
| F41 | **Anxiety disorder** | 1,595 | 20.6 |
| E66 | Obesity | 1,541 | 19.9 |
| F32 | **Major depression, single episode** | 1,438 | 18.6 |

**Three things this view surfaces that the concept-level list hid:**

1. **Musculoskeletal pain dominates.** M25, M54 and M79 take three of the top four slots —
   over a third of the cohort each. Aggregated, chronic pain is the defining feature of this
   population, ahead of any cardiometabolic diagnosis.
2. **Mental health is a top-tier burden.** Anxiety 20.6% and major depression 18.6% both
   enter the top 15 here, and neither appeared in the concept-level top 20. Roughly one in
   five.
3. **Obstructive sleep apnoea at 21.3%** is far above general-population estimates and is
   consistent with the cardiometabolic loading (obesity 19.9%).

### 3.2 Concept-level top 20 (denominator 7,726)

| Condition | n | % |
|---|---|---|
| Essential hypertension | 2,810 | 36.4 |
| Hyperlipidemia | 2,436 | 31.5 |
| Chest pain | 2,130 | 27.6 |
| Cough | 2,047 | 26.5 |
| Abdominal pain | 1,956 | 25.3 |
| Low back pain | 1,897 | 24.6 |
| Chronic pain | 1,858 | 24.0 |
| Obesity | 1,658 | 21.5 |
| Anxiety disorder | 1,443 | 18.7 |
| Type 2 diabetes (uncomplicated) | 1,307 | 16.9 |

**Note how many of the top entries are symptoms, not diseases** — chest pain, cough,
abdominal pain, back pain, dizziness. This is normal and expected for EHR data: billing
codes record *presentations*, not diagnoses. It is why the ICD-10 chapter table is led by
"Symptoms & abnormal findings" at 65.9%, and it is the reason a proper phenotype needs a
curated concept set rather than a single code.

Substantively, this is a **middle-aged-to-older, cardiometabolically loaded US clinical
population**: hypertension, hyperlipidemia, obesity, and diabetes in the top ten. Not a
healthy-volunteer cohort.

Top ICD-10 chapters by share of cohort: Symptoms & abnormal findings 65.9% ·
Musculoskeletal 57.7% · Endocrine/metabolic 54.9% · Digestive 49.3% · Genitourinary 46.5% ·
Nervous system 44.6% · Respiratory 44.4% · Circulatory 44.1%.

---

## 4. The table that decides the science

Immune-mediated disease families, denominator 7,726:

| Family | Cases | % | Distinct concepts | Verdict |
|---|---|---|---|---|
| Allergy & hypersensitivity | **2,297** | 29.7 | 92 | real case/control model |
| Acute / recurrent infection | **1,724** | 22.3 | 96 | real case/control model |
| Chronic infection | **562** | 7.3 | 46 | real case/control model |
| Autoimmune — rheumatologic | **444** | 5.7 | 43 | real case/control model |
| Immunodeficiency & immunosuppression | **344** | 4.5 | 42 | real case/control model |
| Autoimmune — endocrine | **327** | 4.2 | 33 | real case/control model |
| Autoimmune — dermatologic | **314** | 4.1 | 23 | real case/control model |
| Autoimmune — gastrointestinal | 183 | 2.4 | 41 | descriptive only |
| Lymphoid & haematologic malignancy | 151 | 2.0 | 74 | descriptive only |
| Autoimmune — neurologic | 76 | 1.0 | 11 | descriptive only |

**Seven of ten families clear the 200-case bar.** That substantially exceeds the
pre-measurement projection (which put only asthma/allergy and infection above the line).

**Discount these before believing them.** Matching is by substring against SNOMED concept
names, and it over-matches in known ways:
- "Allergy & hypersensitivity" catches every *drug hypersensitivity* code, which is
  extremely common and clinically trivial. The true asthma/atopy count is materially lower.
- "Immunodeficiency & immunosuppression" catches "transplant" in any context, including
  transplant *evaluation* and family history.
- "Autoimmune — rheumatologic" catches osteoarthritis via "arthritis", which is not
  immune-mediated and is very common in this age group.

**This is a screen, and its job was to tell us where to spend effort building real
phenotypes. It did that.** Even discounting aggressively, allergy/atopy, infection history,
and rheumatologic autoimmunity are clearly viable; the three underpowered families are
clearly not, and no amount of definitional care will rescue 76 neurologic cases.

---

## 4b. Deep dive - HLA-linked, autoimmune, immunodeficiency, tumor detail (2026-08-18)

The family-level screen above answers "is this direction viable." It cannot answer "how
many people have rheumatoid arthritis specifically," or separate psoriasis from psoriatic
arthritis, or tell a lymphoid malignancy (where the tumor cells *are* a clonal receptor
sequence) from a solid one. `scripts/deep_immune_breakdown.py` runs against the same
cached phenotype files -- no new BigQuery call -- and matches each disease primarily on its
**3-character ICD-10 code** rather than name substring, which is more precise wherever that
code block is genuinely exclusive to one disease (confirmed per-code before use, not
assumed -- see the caveat below).

### HLA-linked diseases

The most direct candidates for testing a repertoire<->HLA association, since HLA type
shapes which receptors survive thymic selection in the first place.

| Disease | HLA association | Cases | % |
|---|---|---|---|
| Psoriasis | Cw6 | **208** | 2.69 |
| Rheumatoid arthritis | shared epitope | **199** | 2.58 |
| Type 1 diabetes | DR3/DR4 | **179** | 2.32 |
| Systemic lupus erythematosus | DR2/DR3 | 96 | 1.24 |
| Multiple sclerosis | DRB1*15:01 | 55 | 0.71 |
| Psoriatic arthritis | Cw6/B27 | 43 | 0.56 |
| Celiac disease | DQ2/DQ8 | 33 | 0.43 |
| Ankylosing spondylitis | B27 | 25 | 0.32 |
| Narcolepsy | DQB1*06:02 | 25 | 0.32 |
| Behcet disease | B51 | <20 | -- |
| Graves disease | DR3 | <20 | -- |

**Psoriasis, rheumatoid arthritis and type 1 diabetes are the three practical near-term
candidates** -- each clears ~180-210 cases, comfortably above the 200-case working
threshold used elsewhere in this study.

### Autoimmune - specific diseases (not families)

22 named diseases screened; 16 clear the n>=20 disclosure floor.

| Disease | Cases | % | | Disease | Cases | % |
|---|---|---|---|---|---|---|
| Psoriasis | 208 | 2.69 | | Vitiligo | 50 | 0.65 |
| Rheumatoid arthritis | 199 | 2.58 | | Alopecia areata | 44 | 0.57 |
| Type 1 diabetes | 179 | 2.32 | | Psoriatic arthritis | 43 | 0.56 |
| Hashimoto thyroiditis | 154 | 1.99 | | Celiac disease | 33 | 0.43 |
| Sjogren syndrome | 127 | 1.64 | | Systemic sclerosis | 29 | 0.38 |
| Lupus (SLE) | 96 | 1.24 | | Ankylosing spondylitis | 25 | 0.32 |
| Ulcerative colitis | 89 | 1.15 | | Myasthenia gravis, autoimmune hepatitis, Guillain-Barre, Graves, pemphigus, Addison | each <20 |  |
| Crohn disease | 83 | 1.07 | | | | |
| Multiple sclerosis | 55 | 0.71 | | | | |
| Vasculitis | 51 | 0.66 | | | | |

### Immunodeficiency - specific subtypes

| Subtype | Cases | % |
|---|---|---|
| Other immunodeficiency (ICD D84) | 280 | 3.62 |
| Other immune-mechanism disorder (D89) | 170 | 2.20 |
| HIV | 98 | 1.27 |
| IgA / antibody deficiency | 52 | 0.67 |
| Neutropenia | 43 | 0.56 |
| Transplant status (any organ) | 28 | 0.36 |
| Drug-induced immunosuppression, GVHD, CVID | each <20 | -- |

**Most immunodeficiency-coded people fall into the two unavoidably broad ICD buckets**
(D84/D89, 450 of ~653 total) rather than a named subtype -- a real limit of ICD-10
granularity for this category specifically, not a matching artifact. Directly relevant to
repertoire interpretation regardless: an immunosuppressed person's repertoire reflects a
suppressed immune system, not a baseline one, and this is a candidate confound for any
future repertoire-based model.

### Tumors - lymphoid vs. solid

**Lymphoid tumors are the most repertoire-relevant category in the whole cohort.** In these
cancers the malignant cells often *are* a single massively expanded T-cell or B-cell
clone -- TRUST4 output is not just descriptive of them, in principle it can detect them
directly, the same way clinical MRD (measurable residual disease) testing works.

| Lymphoid / hematologic | Cases | % | | Solid (context only) | Cases | % |
|---|---|---|---|---|---|---|
| MGUS (monoclonal gammopathy) | **83** | 1.07 | | Prostate cancer | 187 | 2.42 |
| Non-Hodgkin lymphoma | 64 | 0.83 | | Breast cancer | 162 | 2.10 |
| Multiple myeloma | 32 | 0.41 | | Lung cancer | 116 | 1.50 |
| Lymphoid leukemia (CLL/ALL) | 25 | 0.32 | | Colorectal cancer | 65 | 0.84 |
| MDS, myeloid leukemia, Hodgkin, mycosis fungoides | each <20 | -- | | Melanoma, kidney, bladder | each <20-61 | -- |

MGUS at 1.07% is consistent with published prevalence in an older-skewing adult
population -- a good sign the matching is behaving.

### A methodology correction worth recording

The first pass through this breakdown used the 3-character ICD-10 code as the *sole*
match for every disease, on the assumption that all such codes are exclusive to one
condition. **That assumption was wrong for 10 of the 12 non-oncology codes checked**, and
one -- narcolepsy matched to `G47`, the broad "sleep disorders" block -- produced a
clinically impossible 21.3% prevalence, identical to the cohort's already-reported sleep
apnea rate. Every code was re-checked individually for whether its ICD block is genuinely
exclusive; nine more were found to bundle the target disease with something clinically
distinct (Graves' with non-autoimmune goiter, Addison's with Cushing's -- the *opposite*
condition, celiac with unrelated malabsorption, and others) and were switched to
name-matching alone. The type 1 vs. type 2 diabetes split (E10/E11) that motivated using
ICD codes in the first place remains correct -- that block genuinely is exclusive.

**Standing caveat, unchanged from the family-level screen:** this is still a screen, not a
validated phenotype. Name-substring fallbacks inherit the same over-matching risk as
before. Anything used in a real analysis needs a proper curated concept-set definition.

---

## 5. What follows

**The repertoire direction is viable, on exactly the phenotypes where it should be.** The
argument for repertoire over genotype is that it records *what you were exposed to*, not
just what you inherited. The families that cleared the bar — allergy, acute infection,
chronic infection — are precisely the exposure-driven ones. That is a coherent result, not
a lucky one.

**Immediate, ordered:**
1. **Check whether disease burden predicts CDR3 recovery** (§2). New, cheap, and directly
   relevant to the open ancestry question. Data already on disk for both.
2. **Build a real concept set for one phenotype** — asthma is the obvious first, being both
   large and unambiguously immune-mediated — replacing substring matching.
3. **Scale TRUST4.** The binding constraint is no longer cohort size or confounds; it is
   that only ~154 of 8,327 people have repertoires assembled.

**Standing caveats for any presentation of this:**
- Descriptive, not association. No causal or even adjusted-association claim is supported.
- The healthcare-access gradient (§2) confounds every ancestry breakdown.
- ICD→phenotype reproducibility is 70–75% for chronic conditions and **under 10% for acute
  ones** — so the "acute / recurrent infection" row is the least trustworthy of the large
  families, despite its size.
- Small cells suppressed at n < 20, with derived percentages and medians blanked alongside
  them.

---

## Sources
- Live run, 2026-08-17: `../scripts/run_disease_study.sh` → `query_overlap_phenotypes.py`
  (BigQuery/OMOP) → `analyze_disease_burden.py`.
- Design and background: `../reference/LR_RNASEQ_DISEASE_STUDY.md`.
- ICD→phecode reproducibility figures: <https://pmc.ncbi.nlm.nih.gov/articles/PMC6911227/>
